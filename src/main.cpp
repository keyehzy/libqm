#include "Basis.h"
#include "Commutator.h"
#include "Expression.h"
#include "MatrixElements.h"
#include "NormalOrderer.h"
#include "Term.h"
#include "eigen3/Eigen/Dense"

using namespace libqm;

template <typename F, typename... Args>
using is_operator_callable = std::is_invocable<F, const Operator&, Args...>;

template <typename F, typename... Args>
auto transform_operator(F&& f, const Operator& op, Args&&... args)
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

static Expression fourier_transform_operator(const Operator& op, size_t size) {
  Expression result;
  const float type_sign = (op.type() == Operator::Type::Annihilation) ? -1.0f : 1.0f;
  const float factor = 2.0f * std::numbers::pi_v<float> / static_cast<float>(size);

  for (size_t k = 0; k < size; ++k) {
    std::complex<float> coefficient(0.0f, -type_sign * factor * static_cast<float>(k * op.value()));
    coefficient = std::exp(coefficient) / std::sqrt(static_cast<float>(size));

    Operator transformed_op(op.type(), op.spin(), k);
    result += Term(coefficient, {transformed_op});
  }

  return result;
}

struct HubbardModel {
  HubbardModel(float t, float u, size_t size) : t(t), u(u), size(size) {}

  Expression hamiltonian() const {
    Expression result;
    for (size_t i = 0; i < size; ++i) {
      result += t * hopping(i, (i + 1) % size, Operator::Spin::Up);
      result += t * hopping(i, (i + 1) % size, Operator::Spin::Down);
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

    std::vector<float> E(basis.set.size());
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
