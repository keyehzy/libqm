#include "Basis.h"
#include "Commutator.h"
#include "Expression.h"
#include "MatrixElements.h"
#include "NormalOrderer.h"
#include "Term.h"
#include "Logger.h"

#include <armadillo>
#include <set>
#include <random>

using namespace libqm;

template <typename F, typename... Args>
using is_operator_callable = std::is_invocable<F, Operator, Args...>;

template <typename F, typename... Args>
auto transform_operator(F&& f, Operator op, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value,
                        decltype(f(op, std::forward<Args>(args)...))> {
  return f(op, std::forward<Args>(args)...);
}

template <typename F, typename... Args>
auto transform_term(F&& f, const Term& term, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value, Expression> {
  Expression result(term.c);
  for (const auto& op : term.operators) {
    result *= transform_operator(std::forward<F>(f), op, std::forward<Args>(args)...);
  }
  return result;
}

template <typename F, typename... Args>
auto transform_expression(F&& f, const Expression& expr, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value, Expression> {
  Expression result;
  for (const auto& [ops, coeff] : expr.hashmap) {
    result += transform_term(std::forward<F>(f), Term(coeff, ops), std::forward<Args>(args)...);
  }
  return result;
}

static Term sigma_operator(const Term::container_type& element, size_t size) {
  Term::container_type ops;
  for (Operator o : element) {
    Operator transformed_op(Operator(o.type(), o.spin(), (size - o.value()) % size));
    ops.push_back(transformed_op);
  }
  return Term(ops);
}

static Expression fourier_transform_operator(Operator op, size_t size) {
  Expression result;
  const float type_sign = (op.type() == Operator::Type::Annihilation) ? -1.0f : 1.0f;
  const float factor = (2.0f * std::numbers::pi_v<float>) / static_cast<float>(size);

  for (size_t k = 0; k < size; ++k) {
    std::complex<float> coefficient(0.0f,  -type_sign * factor * static_cast<float>(k * op.value()));
    coefficient = std::exp(coefficient) / std::sqrt(static_cast<float>(size));

    Operator transformed_op(op.type(), op.spin(), k);
    result += Term(coefficient, {transformed_op});
  }

  return result;
}

std::pair<size_t, size_t> from_index(size_t index, size_t width, size_t height) {
  size_t x = index % width;
  size_t y = index / width;
  return {index % width, index / width};
}

size_t to_index(size_t x, size_t y, size_t width, size_t height) {
  return y * width + x;
}

Expression fourier_transform_operator_2d(Operator op, size_t width, size_t height) {
  Expression result;
  const float type_sign =
      (op.type() == Operator::Type::Annihilation) ? -1.0f : 1.0f;

  constexpr float pi = std::numbers::pi_v<float>;
  const float factor_x = 2.0f * pi / static_cast<float>(width);
  const float factor_y = 2.0f * pi / static_cast<float>(height);

  auto [x, y] = from_index(op.value(), width, height);

  for (size_t ky = 0; ky < height; ++ky) {
    for (size_t kx = 0; kx < width; ++kx) {
      std::complex<float> coefficient(
          0.0f, -type_sign * (factor_x * static_cast<float>(kx * x) +
                              factor_y * static_cast<float>(ky * y)));
      coefficient =
          std::exp(coefficient) /
          std::sqrt(static_cast<float>(width * height));

      Operator transformed_op(op.type(), op.spin(), to_index(kx, ky, width, height));
      result += Term(coefficient, {transformed_op});
    }
  }
  return result;
}

struct HubbardModel {
  HubbardModel(float t, float u, size_t size) : t(t), u(u), size(size) {}

  Expression ke() const {
    Expression result;
    for (size_t i = 0; i < size; ++i) {
      result += t * hopping(i, (i + 1)%size, Operator::Spin::Up);
      result += t * hopping(i, (i + 1)%size, Operator::Spin::Down);
    }
    return result;
  }

  Expression inter() const {
    Expression result;
    for (size_t i = 0; i < size; ++i) {
      result += u * density_density(Operator::Spin::Up, i, Operator::Spin::Down, i);
    }
    return result;
  }

  Expression hamiltonian() const {
    Expression result;
    for (size_t i = 0; i < size; ++i) {
      result += t * hopping(i, (i + 1)%size, Operator::Spin::Up);
      result += t * hopping(i, (i + 1)%size, Operator::Spin::Down);
    }
    for (size_t i = 0; i < size; ++i) {
      result += u * density_density(Operator::Spin::Up, i, Operator::Spin::Down, i);
    }
    return result;
  }

  float t;
  float u;
  size_t size;
};

struct HubbardModel2D {
  HubbardModel2D(float t, float u, size_t width, size_t height) : t(t), u(u), width(width), height(height) {}

  Expression hamiltonian() const {
    Expression result;
    for (size_t i = 0; i < width; ++i) {
      for (size_t j = 0; j < height; ++j) {
        size_t index = to_index(i, j, width, height);
        size_t index_forward = to_index((i+1)%width, j, width, height);
        size_t index_up = to_index(i, (j+1)%height, width, height);

        result += t * hopping(index, index_forward, Operator::Spin::Up);
        result += t * hopping(index, index_forward, Operator::Spin::Down);

        result += t * hopping(index, index_up, Operator::Spin::Up);
        result += t * hopping(index, index_up, Operator::Spin::Down);
      }
    }
    for (size_t i = 0; i < width*height; ++i) {
      result += u * density_density(Operator::Spin::Up, i, Operator::Spin::Down, i);
    }
    return result;
  }

  float t;
  float u;
  size_t width;
  size_t height;
};

template <typename T>
int sign_of(T val) {
  return (T(0) < val) - (val < T(0));
}

float find_sign(const Expression& H, const Expression& A, const Term& target, NormalOrderer& orderer, float tol) noexcept {
  auto comm1 = commutator(A, H, orderer);
  auto comm2 = commutator(A, comm1, orderer);
  std::complex<float> alpha = comm1.hashmap[target.operators];
  std::complex<float> beta = comm2.hashmap[target.operators];
  std::complex<float> gamma = target.c;

  if (std::norm(beta) < tol*tol) {
    return sign_of((-gamma / alpha).real());
  }

  std::complex<float> alpha_beta = alpha / beta;
  std::complex<float> gamma_beta = gamma / beta;

  std::complex<float> r1 =
      -alpha_beta + std::sqrt(alpha_beta * alpha_beta - 2.0f * gamma_beta);
  std::complex<float> r2 =
      -alpha_beta - std::sqrt(alpha_beta * alpha_beta - 2.0f * gamma_beta);

  if (std::norm(r1) < std::norm(r2)) {
    return sign_of(r1.real());
  } else {
    return sign_of(r2.real());
  }
}

/*
Expression compute_A(const Expression& H,  NormalOrderer& orderer, float tol) {
  Expression A;
  for (const auto& [ops, coeff] : H.hashmap) {
    if (!is_diagonal(ops)) {
      Term term(coeff, ops);
      float sign = find_sign(H, term, orderer, tol);
      A += sign * (term - term.adjoint());
    }
  }
  std::cout << "A:" << A.to_string() << "\n";
  return A;
}

Expression jacobi_diagonalization(const Expression& H1, NormalOrderer& orderer, float tol = 1e-7f) noexcept {
  // 1. Compute K1
  float dl = 0.01f;
  Expression A1 = compute_A(H1, orderer, tol);
  Expression K1 = dl * commutator(A1, H1, orderer);
  return H1 + K1;
  // // 2. Compute K2
  // Expression H2 = H1 + 0.5f * K1;
  // Expression A2 = compute_A(H2, orderer, tol);
  // Expression K2 = dl * commutator(A2, H2, orderer);

  // // 3. Compute K3
  // Expression H3 = H1 + 0.5f * K2;
  // Expression A3 = compute_A(H3, orderer, tol);
  // Expression K3 = dl * commutator(A3, H3, orderer);

  // // 4. Compute K4
  // Expression H4 = H1 + K3;
  // Expression A4 = compute_A(H4, orderer, tol);
  // Expression K4 = dl * commutator(A4, H4, orderer);

  // return H1 + (1.0f / 6.0f) * (K1 + 2.0f * K2 + 2.0f * K3 + K4);
}
*/

void purge(Expression& expr, float tol) {
  for (auto it = expr.hashmap.begin(); it != expr.hashmap.end();) {
    if (it->first.size() > 4/* || std::norm(it->second) < tol*/) {
      it = expr.hashmap.erase(it);
    } else {
      it++;
    }
  }
}

float compute_off_diagonal_norm(const Expression& e) {
  float result = 0;
  for (const auto& [ops, c] : e.hashmap) {
    if (!is_diagonal(ops)) {
      result += std::norm(c);
    }
  }
  return result;
}

// static std::random_device rd;
static std::mt19937 gen(42);
static std::uniform_real_distribution<> dis(0.0, 1.0);

Expression compute_A(const Expression& H, NormalOrderer& orderer, float tol = 1e-7f) {
  std::vector<Term> terms;

  for (const auto& [ops, c] : H.hashmap) {
    if (!is_diagonal(ops)) {
      terms.emplace_back(c, ops);
    }
  }

  Expression A;
  for (const auto& target : terms) {
    Expression a = target.adjoint() - target;
    float sign = find_sign(H, a, target, orderer, tol);
    A += dis(gen) * sign * a;
  }
  return A;
}

static bool finished = false;
Expression jacobi_diagonalization(const Expression& H1,
                                  NormalOrderer& orderer,
                                  float tol = 1e-4f) {
  // 1. Compute K1
  // float dl = 0.01f;
  // Expression A1 = compute_A(H1, orderer, tol, true);
  // Expression K1 = dl * commutator(A1, H1, orderer);
  // purge(K1, tol);

  // // 2. Compute K2
  // Expression H2 = H1 + 0.5f * K1;
  // Expression A2 = compute_A(H2, orderer, tol);
  // Expression K2 = dl * commutator(A2, H2, orderer);
  // purge(K2, tol);

  // // 3. Compute K3
  // Expression H3 = H1 + 0.5f * K2;
  // Expression A3 = compute_A(H3, orderer, tol);
  // Expression K3 = dl * commutator(A3, H3, orderer);
  // purge(K3, tol);

  // // 4. Compute K4
  // Expression H4 = H1 + K3;
  // Expression A4 = compute_A(H4, orderer, tol);
  // Expression K4 = dl * commutator(A4, H4, orderer);
  // purge(K4, tol);

  // Expression result = H1 + (1.0f / 6.0f) * (K1 + 2.0f * K2 + 2.0f * K3 + K4);

  Expression A = compute_A(H1, orderer, tol);

  /*
  bool accepted = false;
  float dl = 1.0f;
  Expression result2, result4;
  float y1, y2, error;

  do {
    result2 = BCH(A, H1, dl, orderer, 2);
    purge(result2, tol);

    result4 = BCH(A, H1, dl, orderer, 4);
    purge(result4, tol);

    y1 = compute_off_diagonal_norm(result2);
    y2 = compute_off_diagonal_norm(result4);

    error = std::abs(y2 - y1);

    if (error < tol) {
      break;
    }

    dl = std::max(dl * 0.1f, std::min(dl * 5.0f, dl * 0.95f * std::sqrt(tol / error)));

  } while(!accepted);

  std::cout << "size:" << result4.hashmap.size() << "\n";
  std::cout << "off-diagonal:" << y2 << "\n";
  std::cout << "error:" << error << "\n";
  std::cout << "dl:" << dl << "\n";

  if (y2 < tol) {
    finished = true;
  }
  */
  Expression result = BCH(A, H1, tol, orderer, 2);
  purge(result, tol);
  float y2 = compute_off_diagonal_norm(result);
  Logger::the().info_fmt("size: {}, off-diagonal: {}", result.hashmap.size(), y2);
  if (y2 < tol) {
    finished = true;
  }

  return result;

  /*
  std::vector<Term> terms;
  std::set<size_t> indices;

  auto collect_indices = [](const Term::container_type& container) {
    std::set<size_t> idxs;
    for (size_t i = 0; i < container.size(); i++) {
      idxs.insert(container[i].key());
    }
    return idxs;
  };

  std::vector<Term> Hterms;
  Hterms.reserve(H.hashmap.size());
  for (const auto& [ops, coeff] : H.hashmap) {
    Hterms.emplace_back(coeff, ops);
  };
  std::sort(Hterms.begin(), Hterms.end(), [](const Term& lhs, const Term& rhs) {
    return std::norm(lhs.c) < std::norm(rhs.c);
  });

  float off_diagonal_norm = 0.0f;

  for (const auto& term : Hterms) {
    if (!is_diagonal(term.operators)) {
      off_diagonal_norm += std::norm(term.c);
      auto term_indices = collect_indices(term.operators);
      bool has_no_intersection = std::none_of(
          term_indices.begin(), term_indices.end(), [&indices](size_t idx) {
            return indices.find(idx) != indices.end();
          });

      if (has_no_intersection) {
        terms.push_back(term);
        indices.insert(term_indices.begin(), term_indices.end());
      }
    }
  }

  std::cout << off_diagonal_norm << "\n";

  if (indices.empty()) {
    std::cout << "No off diagonal element found.\n";
    return H;
  }

  Expression A;
  for (const auto& target : terms) {
    float sign = find_sign(H, target, orderer, tol);
    A += 0.01f * sign * (target.adjoint() - target);
  }

  Expression h = BCH(A, H, 1.0f, orderer, 2);
  purge(h, tol);
  return h;
  */
}

int main2() {
  const float t = -1.0f;
  const float u = 100.0f;
  const size_t num_sites = 4;
  const size_t num_particles = 4;
  Basis basis(num_sites, num_particles, Basis::Strategy::Restrict);
  HubbardModel model(t, u, num_sites);
  NormalOrderer orderer;

  {
    Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));
    auto hamilt_matrix_cx = compute_matrix_elements<arma::cx_fmat>(basis, hamilt);
    arma::fvec vals;
    arma::eig_sym(vals, hamilt_matrix_cx);
    std::cout << vals << "\n";
  }


  {
    std::cout << "Basis size:" << basis.set.size() << "\n";

    Expression total_S;
    for (size_t i = 0; i < num_particles; i++) {
      total_S += spin_x(i) + spin_y(i) + spin_z(i);
    }

    Expression total_P;
    for (size_t i = 0; i < num_particles; i++) {
      float k = static_cast<float>(i) * (2.0f * std::numbers::pi_v<float> + std::numbers::pi_v<float>) / static_cast<float>(num_sites);
      total_P += k * density(Operator::Spin::Up, i);
      total_P += k * density(Operator::Spin::Down, i);
    }

    for (float target_S : {0, 1, 2}) {
    Basis S_basis = basis.refine([&] (const auto& elm) {
      Expression result = orderer.normal_order(Term(elm).adjoint() * total_S * Term(elm));
      // std::cout << "S:" << result.hashmap[{}] << "\n";
      return std::abs(std::abs(result.hashmap[{}]) - target_S) < 1e-4f;
    });


      for (float target_Sz : {-2, -1, 0, 1, 2}) {
       Basis Sz_basis = S_basis.refine([&] (const auto& elm) {
         int Sz = 0;
         for (Operator op : elm) {
           Sz += (op.spin() == Operator::Spin::Up ? +0.5f : -0.5f);
         }
         // std::cout << "Sz:" << Sz << "\n";
         return std::abs(Sz - target_Sz) < 1e-4f;
       });

        // for (float target_P : {0, 1, 2, 3}) {
        // Basis P_basis = Sz_basis.refine([&] (const auto& elm) {
        //   Expression result = orderer.normal_order(Term(elm).adjoint() * total_P * Term(elm));
        //   float P = std::fmod(result.hashmap[{}].real(), 2.0f * std::numbers::pi_v<float>);
        //   std::cout << "P:" << P/(std::numbers::pi_v<float>/4.f) << "\n";
        //   return std::abs(P/(std::numbers::pi_v<float>/4.f) - target_P) < 1e-4f;
        // });
        // 
        // arma::cx_fmat sigma(P_basis.set.size(), P_basis.set.size());
        // for (size_t j = 0; j < P_basis.set.size(); ++j) {
        //   Expression transf = orderer.normal_order(sigma_operator(P_basis.set[j], num_sites));
        //   for (const auto& term : transf.hashmap) {
        //     if (P_basis.set.contains(term.first)) {
        //       size_t i = P_basis.set.index_of(term.first);
        //       sigma(i, j) = term.second;
        //     }
        //   }
        // }
        // std::cout << "sigma:\n" << sigma;
        // 
        // arma::fvec sigma_vals;
        // arma::cx_fmat sigma_vecs;
        // {
        //   arma::eig_sym(sigma_vals, sigma_vecs, sigma);
        //   // std::cout << "sigma:\n" << sigma_vals << "\n";
        // }

        std::cout << "Reduced size:" << Sz_basis.set.size() << "\n";
        Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));
        auto hamilt_matrix_cx = compute_matrix_elements<arma::cx_fmat>(Sz_basis, hamilt);

        arma::cx_fmat sigma_transformed_h = hamilt_matrix_cx;
        std::cout << target_S << " " << target_Sz << "\n";
        std::cout << arma::conv_to<arma::fmat>::from(sigma_transformed_h.clean(0.01f)) << "\n";


        // for (const auto& elm : p1.set) {
        //   std::cout << Term(elm).to_string() << "\n";
        // }

        // std::cout << arma::conv_to<arma::fmat>::from(hamilt_matrix_cx.clean(0.01f)) << "\n";

        // LIBQM_ASSERT(arma::abs(arma::imag(hamilt_matrix_cx)).max() < 1e-5f);
        // auto hamilt_matrix = arma::conv_to<arma::fmat>::from(hamilt_matrix_cx);
        // hamilt_matrix.save("./raw4f.bin", arma::raw_binary);



        // std::cout << hamilt_matrix << "\n";
        // std::cout << basis.set.size() << "\n";
        // size_t idx = 1;
        // for (const auto& elm : basis.set) {
        //   std::cout << idx++ << " " << Term(elm).to_string() << "\n";
        // }
        // 
        // 
        // 
        arma::fvec vals;
        arma::eig_sym(vals, sigma_transformed_h);
        std::cout << vals << "\n";
        //}
        }
        }
  }
}

int main() {
  const float t = -1.0f;
  const float u = 100.0f;
  const size_t num_sites = 4;
  const size_t num_particles = 4;
  HubbardModel model(t, u, num_sites);
  NormalOrderer orderer;

  Logger::the().info_fmt("parameters: t={}, u={}, sites={}, particles={}", t, u, num_sites, num_particles);

  Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));
  Logger::the().info("hamiltonian");

  size_t iters = 0;
  while(!finished) {
    Logger::the().info_fmt("Iter: {}", iters++);
    hamilt = jacobi_diagonalization(hamilt, orderer);
  }
  std::cout << hamilt.to_string() << "\n\n";
  return 0;
}

int main4() {
  const float t = -1.0f;
  const float u = 100.0f;
  const size_t num_sites = 14;
  const size_t num_particles = 14;
  HubbardModel model(t, u, num_sites);
  NormalOrderer orderer;

  Logger::the().info_fmt("parameters: t={}, u={}, sites={}, particles={}", t, u, num_sites, num_particles);

  //{
  //  Basis basis(num_sites, num_particles, Basis::Strategy::Restrict);
  //  std::cout << "Basis size:" << basis.set.size() << "\n";
  //  Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));
  //  auto hamilt_matrix_cx = compute_matrix_elements<arma::cx_fmat>(basis, hamilt);
  //  arma::fvec vals;
  //  arma::eig_sym(vals, hamilt_matrix_cx);
  //  std::cout << vals << "\n";
  //}


  Basis reduced_basis(num_sites, num_particles, Basis::Strategy::Restrict_Sz_P);
  Logger::the().info_fmt("Reduced basis size: {}", reduced_basis.set.size());

  Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));
  auto imag_hamilt_matrix_cx = compute_matrix_elements<arma::sp_cx_fmat>(reduced_basis, hamilt);
  Logger::the().info_stream(arma::abs(arma::imag(imag_hamilt_matrix_cx)).max());
  LIBQM_ASSERT(arma::abs(arma::imag(imag_hamilt_matrix_cx)).max() < 1e-4f);
  auto hamilt_matrix = arma::conv_to<arma::sp_fmat>::from(imag_hamilt_matrix_cx);
  Logger::the().info("hamiltonian matrix");
  // std::cout << arma::conv_to<arma::fmat>::from(hamilt_matrix_cx.clean(0.01f)) << "\n";

  arma::fvec vals;
  arma::fmat vecs;
  arma::eigs_sym(vals, vecs, hamilt_matrix, 1, "sa");
  Logger::the().info_stream("eigenvalue:", vals(0));
  // Logger::the().info_stream("eigenvector:", vecs.col(0).t());

  Expression docc;
  for (size_t i = 0; i < num_sites; i++) {
    docc +=
      Expression({
          density_density(Operator::Spin::Up, i, Operator::Spin::Down, i) / static_cast<float>(num_sites),
        });
  }
  auto imag_docc_matrix_cx = compute_matrix_elements<arma::sp_cx_fmat>(reduced_basis, docc);
  LIBQM_ASSERT(arma::abs(arma::imag(imag_docc_matrix_cx)).max() < 1e-5f);
  auto docc_matrix = arma::conv_to<arma::sp_fmat>::from(imag_docc_matrix_cx);
  Logger::the().info("double occ matrix");

  arma::fmat a1 = (vecs.col(0).t() * docc_matrix * vecs.col(0));
  Logger::the().info_stream("double occ exp val:", num_sites, " ", a1);
}


int main3() {
  const float t = -1.0f;
  const float u = 100.f;

  const size_t num_sites = 6;
  const size_t num_particles = 6;

  for (size_t num_sites = 3; num_sites <= 10; num_sites++) {
  size_t num_particles = num_sites;
  Basis basis(num_sites, num_particles, Basis::Strategy::Restrict);
  // std::cout << "Size:" << basis.set.size() << "\n";

  // for (const auto& el : basis.set) {
  //   std::cout << Term(el).to_string() << "\n";
  // }
  // std::cout << "\n";

  HubbardModel model(t, u, num_sites);
  NormalOrderer orderer;

  Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));

  // Expression ke = orderer.normal_order(model.ke());
  // Expression inter = orderer.normal_order(model.inter());
  // Expression hamilt = orderer.normal_order(ke+inter);

  // Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));
  // Expression hamilt = orderer.normal_order(model.hamiltonian());


  arma::cx_fvec vals;
  arma::cx_fmat vecs;
  {
    auto hamilt_matrix_cx = compute_matrix_elements<arma::sp_cx_fmat>(basis, hamilt);
    // auto hamilt_matrix = arma::conv_to<arma::cx_fmat>::from(hamilt_matrix_cx);
    // LIBQM_ASSERT(arma::abs(arma::imag(hamilt_matrix_cx)).max() < 1e-5f);
    // auto hamilt_matrix = arma::conv_to<arma::sp_fmat>::from(hamilt_matrix_cx);
    arma::eigs_gen(vals, vecs, hamilt_matrix_cx, 1, "sr");
  }
  // std::cout << vals.size() << "\n";
  // std::cout << vals << "\n";

  // Expression sz;
  // for (size_t i = 0; i < num_sites; i++) {
  //   sz += spin_z(i)*spin_z(i);
  // }
  // auto sz_matrix = compute_matrix_elements<arma::cx_fmat>(basis, sz);

  Expression docc;
  for (size_t i = 0; i < num_sites; i++) {
    docc +=
      Expression({
          density_density(Operator::Spin::Up, i, Operator::Spin::Down, i) / static_cast<float>(num_sites),
        });
  }

  arma::cx_fmat a1;
  {
    auto docc_matrix_cx = compute_matrix_elements<arma::sp_cx_fmat>(basis, docc);
    // auto docc_matrix = arma::conv_to<arma::cx_fmat>::from(docc_matrix_cx);
    // LIBQM_ASSERT(arma::abs(arma::imag(docc_matrix_cx)).max() < 1e-5f);
    // auto docc_matrix = arma::conv_to<arma::sp_fmat>::from(docc_matrix_cx);
    // std::cout << angle << " " << vecs.col(0).t() * docc_matrix_cx * vecs.col(0);
    // auto proj = vecs(arma::span::all, arma::span(0, 1/*(1<<num_sites*num_sites) - 1*/));
    a1 = (vecs.col(0).t() * docc_matrix_cx * vecs.col(0));
  }
  std::cout << num_sites << " " << arma::real(a1);

  //{
  //  arma::fvec vs;
  //  arma::eig_sym(vs, a1);
  //  // std::cout << u << " " << vs.t();
  //  std::cout << num_sites << " "  << arma::sum(vs) << "\n";
  //}

  }

  // size_t idx;
  // for (size_t i = 0; i < basis.set.size(); i++) {
  //   auto x = vecs.col(0)[i];
  //   if (std::norm(x) < 1e-6f) continue;
  //   std::cout << Term(x, basis.set.at(i)).to_string() << " " << std::norm(x) << "\n";
  // }
  // std::cout << "\n\n";

  // size_t iters = 0;
  // while(!finished) {
  //      std::cout << "Iter:" << iters++ << "\n";
  //   h2 = jacobi_diagonalization(h2, orderer);
  // }
  // std::cout << h2.to_string() << "\n\n";

  return 0;
}



/*
int main() {
  const float t = -1.0f;
  const float u = 2.0f;

  const size_t orbitals = 4;
  const size_t particles = 3;
  HubbardModel model(t, u, orbitals);

  // Diagonalize one body using fourier transform
  Expression h2 = transform_expression(fourier_transform_operator, model.hamiltonian(), model.size);
  // std::cout << h2.to_string() << "\n\n";

  Basis basis(orbitals, particles, Basis::Strategy::All);
  NormalOrderer orderer;

  for (size_t iters = 0; iters < 1000; ++iters) {
    std::cout << "# Iter " << iters << "\n";

    std::rngstor<float> E(basis.set.size());
    for (size_t i = 0; i < basis.set.size(); ++i) {
      Expression b(basis.set[i]);
      E[i] = orderer.normal_order(b.adjoint() * h2 * b).hashmap[{}].real();
    }

    Expression result;
    for (size_t i = 0; i < basis.set.size(); ++i) {
      for (size_t j = 0; j < basis.set.size(); ++j) {
        Term b1(basis.set[i]);
        Term b2(basis.set[j]);
        Term b1b2 = b1 * b2.adjoint();
        std::sort(b1b2.operators.begin(), b1b2.operators.end());
        std::complex<float> alpha12 = h2.hashmap[b1b2.operators];
        if (std::norm(alpha12) < 1e-6f) {
          continue;
        }
        for (size_t k = 0; k < basis.set.size(); ++k) {
          for (size_t l = 0; l < basis.set.size(); ++l) {
            if (std::abs(E[k] - E[l]) < 1e-6f) {
              continue;
            }
            Term b3(basis.set[k]);
            Term b4(basis.set[l]);
            Term b3b4 = b3 * b4.adjoint();
            std::sort(b3b4.operators.begin(), b3b4.operators.end());
            std::complex<float> alpha34 = h2.hashmap[b3b4.operators];
            if (std::norm(alpha34) < 1e-6f) {
              continue;
            }
            // std::cout << "foo\n";
            Expression comm = commutator(b1b2, b3b4, orderer);
            result += alpha12 * alpha34 * comm / (E[k] - E[l]);
          }
        }
      }
    }
    h2 = orderer.normal_order(h2 + 3e-3f * result);
    Expression new_h;
    for (const auto& term : h2.hashmap) {
      if (term.first.size() <= 4) {
        new_h += Term(term.second, term.first);
      }
    }
    h2 = std::move(new_h);
    for (const auto& term : h2.hashmap) {
      std::cout << std::norm(term.second) << "\n";
    }
    std::cout << "\n";
  }

  return 0;
}
*/
