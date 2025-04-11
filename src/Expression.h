#pragma once

#include <unordered_map>

#include "Term.h"

namespace libqm {
struct Expression {
  using complex_type = Term::complex_type;
  using container_type = Term::container_type;
  using map_type = std::unordered_map<container_type, complex_type>;

  map_type hashmap;

  Expression() = default;
  ~Expression() = default;

  Expression(const Expression&) = default;
  Expression& operator=(const Expression&) = default;
  Expression(Expression&&) = default;
  Expression& operator=(Expression&&) = default;

  explicit Expression(complex_type c);
  explicit Expression(Operator op);
  explicit Expression(const Term& term);
  explicit Expression(Term&& term);
  explicit Expression(const container_type& container);
  explicit Expression(container_type&& container);
  explicit Expression(complex_type c, const container_type& container);
  explicit Expression(complex_type c, container_type&& container);
  explicit Expression(std::initializer_list<Term> lst);

  Expression adjoint() const;

  size_t size() const { return hashmap.size(); }

  std::string to_string() const;

  Expression& operator+=(complex_type value);
  Expression& operator-=(complex_type value);
  Expression& operator*=(complex_type value);
  Expression& operator/=(complex_type value);

  Expression& operator+=(const Expression& value);
  Expression& operator-=(const Expression& value);
  Expression& operator*=(const Expression& value);

  Expression& operator+=(const Term& value);
  Expression& operator-=(const Term& value);
  Expression& operator*=(const Term& value);
};

Expression operator+(const Expression& a, const Expression& b);
Expression operator-(const Expression& a, const Expression& b);
Expression operator*(const Expression& a, const Expression& b);

Expression operator+(const Expression& a, const Term& b);
Expression operator-(const Expression& a, const Term& b);
Expression operator*(const Expression& a, const Term& b);

Expression operator+(const Expression& a, Expression::complex_type b);
Expression operator-(const Expression& a, Expression::complex_type b);
Expression operator*(const Expression& a, Expression::complex_type b);
Expression operator/(const Expression& a, Expression::complex_type b);

Expression operator+(Expression::complex_type a, const Expression& b);
Expression operator-(Expression::complex_type a, const Expression& b);
Expression operator*(Expression::complex_type a, const Expression& b);

Expression operator+(const Term& a, const Expression& b);
Expression operator-(const Term& a, const Expression& b);
Expression operator*(const Term& a, const Expression& b);

Expression operator+(const Term& a, const Term& b);
Expression operator-(const Term& a, const Term& b);

Expression hopping(size_t from_orbital, size_t to_orbital, Operator::Spin spin);
Expression hopping(Expression::complex_type c, size_t from_orbital, size_t to_orbital,
                   Operator::Spin spin);

Expression spin_x(size_t i);
Expression spin_y(size_t i);
Expression spin_z(size_t i);
Expression spin_plus(size_t site);
Expression spin_minus(size_t site);
Expression spin_dot_product(size_t i, size_t j);
}  // namespace libqm
