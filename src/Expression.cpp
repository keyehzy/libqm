#include "Expression.h"

#include "Assert.h"

namespace libqm {
Expression::Expression(complex_type c) { hashmap.emplace(container_type(), c); }

Expression::Expression(Operator op) { hashmap.emplace(container_type({op}), 1.0f); }

Expression::Expression(const Term& term) { hashmap.emplace(term.operators, term.c); }

Expression::Expression(Term&& term) { hashmap.emplace(std::move(term.operators), term.c); }

Expression::Expression(const container_type& container) { hashmap.emplace(container, 1.0f); }

Expression::Expression(container_type&& container) { hashmap.emplace(std::move(container), 1.0f); }

Expression::Expression(complex_type c, const container_type& container) {
  hashmap.emplace(container, c);
}

Expression::Expression(complex_type c, container_type&& container) {
  hashmap.emplace(std::move(container), c);
}

Expression::Expression(std::initializer_list<Term> lst) {
  hashmap.reserve(lst.size());
  for (const auto& term : lst) {
    hashmap[term.operators] += term.c;
  }
}

Expression Expression::adjoint() const {
  Expression result;
  result.hashmap.reserve(hashmap.size());
  for (const auto& [ops, c] : hashmap) {
    container_type adjoint_ops;
    for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
      adjoint_ops.push_back(it->adjoint());
    }
    result.hashmap.emplace(std::move(adjoint_ops), std::conj(c));
  }
  return result;
}

std::string Expression::to_string() const {
  if (hashmap.empty()) {
    return std::string();
  }
  
  std::string result;

  for (const auto& [ops, c] : hashmap) {
    std::string term_string = Term(c, ops).to_string();
    if (!term_string.empty()) {
      result += term_string + "\n";
    }
  }

  return result;
}

Expression& Expression::operator+=(complex_type other) {
  hashmap[{}] += other;
  return *this;
}

Expression& Expression::operator-=(complex_type other) {
  hashmap[{}] -= other;
  return *this;
}

Expression& Expression::operator*=(complex_type other) {
  for (auto& [ops, c] : hashmap) {
    c *= other;
  }
  return *this;
}

Expression& Expression::operator/=(complex_type other) {
  for (auto& [ops, c] : hashmap) {
    c /= other;
  }
  return *this;
}

Expression& Expression::operator+=(const Expression& other) {
  for (const auto& [ops, c] : other.hashmap) {
    hashmap[ops] += c;
  }
  return *this;
}

Expression& Expression::operator-=(const Expression& other) {
  for (const auto& [ops, c] : other.hashmap) {
    hashmap[ops] -= c;
  }
  return *this;
}

Expression& Expression::operator*=(const Expression& other) {
  map_type result;
  result.reserve(hashmap.size() * other.hashmap.size());
  for (const auto& [ops1, c1] : hashmap) {
    for (const auto& [ops2, c2] : other.hashmap) {
      auto combined_ops = merge(ops1, ops2);
      result[combined_ops] += c1 * c2;
    }
  }
  hashmap = std::move(result);
  return *this;
}

Expression& Expression::operator+=(const Term& other) {
  hashmap[other.operators] += other.c;
  return *this;
}

Expression& Expression::operator-=(const Term& other) {
  hashmap[other.operators] -= other.c;
  return *this;
}

Expression& Expression::operator*=(const Term& other) {
  map_type result;
  result.reserve(hashmap.size());
  for (const auto& [ops1, c1] : hashmap) {
    auto combined_ops = merge(ops1, other.operators);
    result[combined_ops] += c1 * other.c;
  }
  hashmap = std::move(result);
  return *this;
}

Expression operator+(const Expression& a, const Expression& b) {
  Expression result(a);
  result += b;
  return result;
}

Expression operator-(const Expression& a, const Expression& b) {
  Expression result(a);
  result -= b;
  return result;
}

Expression operator*(const Expression& a, const Expression& b) {
  Expression result(a);
  result *= b;
  return result;
}

Expression operator+(const Expression& a, const Term& b) {
  Expression result(a);
  result += b;
  return result;
}

Expression operator-(const Expression& a, const Term& b) {
  Expression result(a);
  result -= b;
  return result;
}

Expression operator*(const Expression& a, const Term& b) {
  Expression result(a);
  result *= b;
  return result;
}

Expression operator+(const Expression& a, Expression::complex_type b) {
  Expression result(a);
  result += b;
  return result;
}

Expression operator-(const Expression& a, Expression::complex_type b) {
  Expression result(a);
  result -= b;
  return result;
}

Expression operator*(const Expression& a, Expression::complex_type b) {
  Expression result(a);
  result *= b;
  return result;
}

Expression operator/(const Expression& a, Expression::complex_type b) {
  Expression result(a);
  result /= b;
  return result;
}

Expression operator+(Expression::complex_type a, const Expression& b) {
  Expression result(b);
  result += a;
  return result;
}

Expression operator-(Expression::complex_type a, const Expression& b) {
  Expression result(a);
  result -= b;
  return result;
}

Expression operator*(Expression::complex_type a, const Expression& b) {
  Expression result(b);
  result *= a;
  return result;
}

Expression operator+(const Term& a, const Expression& b) {
  Expression result(a);
  result += b;
  return result;
}

Expression operator-(const Term& a, const Expression& b) {
  Expression result(a);
  result -= b;
  return result;
}

Expression operator*(const Term& a, const Expression& b) {
  Expression result(a);
  result *= b;
  return result;
}

Expression operator+(const Term& a, const Term& b) {
  Expression result(a);
  result += b;
  return result;
}

Expression operator-(const Term& a, const Term& b) {
  Expression result(a);
  result -= b;
  return result;
}

Expression hopping(size_t f, size_t t, Operator::Spin s) {
  LIBQM_ASSERT(f != t);
  return one_body(s, t, s, f) + one_body(s, f, s, t);
}

Expression hopping(Expression::complex_type c, size_t f, size_t t, Operator::Spin s) {
  LIBQM_ASSERT(f != t);
  return c * one_body(s, t, s, f) + std::conj(c) * one_body(s, f, s, t);
}

Expression spin_x(size_t i) {
  return 0.5f * one_body(Operator::Spin::Up, i, Operator::Spin::Down, i) +
         0.5f * one_body(Operator::Spin::Down, i, Operator::Spin::Up, i);
}

Expression spin_y(size_t i) {
  return Expression::complex_type(0.0f, 0.5f) *
             one_body(Operator::Spin::Up, i, Operator::Spin::Down, i) -
         Expression::complex_type(0.0f, 0.5f) *
             one_body(Operator::Spin::Down, i, Operator::Spin::Up, i);
}

Expression spin_z(size_t i) {
  return 0.5f * one_body(Operator::Spin::Up, i, Operator::Spin::Up, i) -
         0.5f * one_body(Operator::Spin::Down, i, Operator::Spin::Down, i);
}

Expression spin_plus(size_t i) {
  return Expression(one_body(Operator::Spin::Up, i, Operator::Spin::Down, i));
}

Expression spin_minus(size_t i) {
  return Expression(one_body(Operator::Spin::Down, i, Operator::Spin::Up, i));
}

Expression spin_dot_product(size_t i, size_t j) {
  return spin_x(i) * spin_x(j) + spin_y(i) * spin_y(j) + spin_z(i) * spin_z(j);
}
}  // namespace libqm
