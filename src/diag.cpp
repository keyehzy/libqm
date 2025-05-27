#include "Basis.h"
#include "Commutator.h"
#include "Expression.h"
#include "MatrixElements.h"
#include "NormalOrderer.h"
#include "Term.h"
#include "Logger.h"

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

void purge(Expression& expr, size_t count, float tol) {
  for (auto it = expr.hashmap.begin(); it != expr.hashmap.end();) {
    if (it->first.size() > count || std::norm(it->second) < tol) {
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
      result = std::max(result, std::norm(c));
    }
  }
  return result;
}

// static std::random_device rd;
static std::mt19937 gen(42);
static std::uniform_real_distribution<> dis(0.0, 1.0);
static std::unordered_map<Term::container_type, float> momentum;

Expression compute_A(const Expression& H, NormalOrderer& orderer, float u, float tol = 1e-6f) {
  std::vector<Term> terms;

  for (const auto& [ops, c] : H.hashmap) {
    if (!is_diagonal(ops)) {
      terms.emplace_back(c, ops);
    }
  }

  Expression A;
  float alpha = 100.0f;
  float beta = 0.99f;
  for (const auto& target : terms) {
    Expression a = target.adjoint() - target;
    float sign = find_sign(H, a, target, orderer, tol);
    float grad = sign * dis(gen);
    float z = beta * momentum[target.operators] + grad;
    momentum[target.operators] = z;
    A += alpha * z * a;
  }
  return A;
}

static bool finished = false;
Expression jacobi_diagonalization(const Expression& H1, NormalOrderer& orderer, std::vector<Expression>& exprs, float u, float tol = 1e-6f) {
  Expression A = compute_A(H1, orderer, u, tol);
  Expression result = BCH(A, H1, tol, orderer, 2);
  purge(result, 6, tol);
  // for (Expression &e : exprs) {
  //   e = BCH(A, e, tol, orderer, 2);
  //   purge(e, 4, tol);
  // }
  float y2 = compute_off_diagonal_norm(result);
  Logger::the().info_fmt("size: {}, off-diagonal: {}", result.hashmap.size(), y2);
  if (y2 < tol) {
    finished = true;
  }

  return result;
}

void print_expression(const Expression& expr, const std::vector<Term::container_type>& basis, NormalOrderer& orderer, float factor) {
  Expression normed_expr = orderer.normal_order(expr);
  std::vector<Term> terms;
  terms.reserve(basis.size());
  for (const auto& e : basis) {
    if (auto it = normed_expr.hashmap.find(e); it != normed_expr.hashmap.end()) {
      terms.emplace_back(factor*it->second, e);
    } else {
      terms.emplace_back(0, e);
    }
  }
  /*
  std::sort(terms.begin(), terms.end(), [] (const auto& l, const auto& r) {
    if (l.size() != r.size()) {
      return l.size() < r.size();
    } else {
      return l.operators < r.operators;
    }
  });
    auto special_less = [] (const auto& l, const auto& r) {
    if (l.type() != r.type()) {
      return l.type() < r.type();
    } else if (l.spin() != r.spin()) {
      if (l.type() == Operator::Type::Creation && r.type() == Operator::Type::Annihilation) {
        return l.spin() < r.spin();
      } else {
        return l.spin() > r.spin();
      }
    } else {
      if (l.type() == Operator::Type::Creation && r.type() == Operator::Type::Annihilation) {
        return l.value() < r.value();
      } else {
        return l.value() > r.value();
      }
    }
  };
  for (auto& term : terms) {
    for (size_t i = 1; i < term.operators.size(); ++i) {
      size_t j = i;
      while (j > 0 && special_less(term.operators[j], term.operators[j - 1])) {
        if (term.operators[j].commutes(term.operators[j - 1])) {
          std::swap(term.operators[j], term.operators[j - 1]);
          term.c = -term.c;
          --j;
        }
      }
    }
  }
  */
  for (const auto& t : terms) {
    std::string ts = t.to_string();
    if (!ts.empty()) {
      Logger::the().info_stream(ts);
    }
  }
}

std::vector<Term::container_type> construct_basis_index(size_t num_sites) {
  //Basis one_body(num_sites, 1, Basis::Strategy::Restrict);
  // Basis two_body(num_sites, 2, Basis::Strategy::Restrict);
  Basis three_body(num_sites, 3, Basis::Strategy::Restrict);
  // Basis four_body(num_sites, 4, Basis::Strategy::Restrict);

  std::vector<Term::container_type> elements;

  // for (const auto& e1 : one_body.set) {
  //   for (const auto& e2 : one_body.set) {
  //     Term::container_type e2_adj(e2.rbegin(), e2.rend());
  //     std::transform(e2_adj.begin(), e2_adj.end(), e2_adj.begin(), [](auto& o){ return o.adjoint(); });
  //     auto e = merge(e1, e2_adj);
  //     std::sort(e.begin(), e.end());
  //     elements.push_back(std::move(e));
  //   }
  // }

  for (const auto& e1 : three_body.set) {
    for (const auto& e2 : three_body.set) {
      Term::container_type e2_adj(e2.rbegin(), e2.rend());
      std::transform(e2_adj.begin(), e2_adj.end(), e2_adj.begin(), [](auto& o){ return o.adjoint(); });
      auto e = merge(e1, e2_adj);
      std::sort(e.begin(), e.end());
      elements.push_back(std::move(e));
    }
  }
  return elements;
}

int main() {
  const float t = -1.0f;
  // const float u = 10.0f;
  const size_t num_sites = 5;
  const size_t num_particles = 5;
  NormalOrderer orderer;

  std::vector<Term::container_type> elements = construct_basis_index(num_sites);

  for (float u = 100.0f; u < 100.01f; u += 0.2f) {
    float tt = t / (0.5f*(-t+u));
    float uu = u / (0.5f*(-t+u));
    Logger::the().info_fmt("Effective: t={}, u={}", tt, uu);
    HubbardModel model(tt, uu, num_sites);

    Logger::the().info_fmt("parameters: t={}, u={}, sites={}, particles={}", t, u, num_sites, num_particles);

    Expression hamilt = orderer.normal_order(transform_expression(fourier_transform_operator, model.hamiltonian(), model.size));
    Logger::the().info("hamiltonian");

    std::vector<Expression> exprs = {
      Expression(Operator::creation(Operator::Spin::Up, 0)),
      // Expression(Term{1.0f, {Operator::creation(Operator::Spin::Up, 0),Operator::creation(Operator::Spin::Down, 0)}}),
    };

    size_t iters = 0;
    while(!finished) {
      // Logger::the().info_fmt("Iter: {}", iters++);
      hamilt = jacobi_diagonalization(hamilt, orderer, exprs, u, 1e-5f);
    }

    std::cout << "Expression:" << "\n";
    print_expression(hamilt, elements, orderer, (0.5f*(-t+u)));

    // for (const Expression& e : exprs) {
    //   print_expression(e, elements, orderer, (0.5f*(-t+u)));
    // }
    finished = false;
  }

  return 0;
}
