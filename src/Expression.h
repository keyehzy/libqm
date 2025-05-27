#pragma once

#include <boost/unordered/unordered_flat_map.hpp>

#include "Term.h"

namespace libqm {
struct Expression {
  using complex_type = Term::complex_type;
  using container_type = Term::container_type;
  using map_type = boost::unordered_flat_map<container_type, complex_type>;

  map_type hashmap;

  Expression() = default;
  ~Expression() = default;

  Expression(const Expression&) = default;
  Expression& operator=(const Expression&) = default;
  Expression(Expression&&) noexcept = default;
  Expression& operator=(Expression&&) noexcept = default;

  explicit Expression(complex_type c);
  explicit Expression(Operator op);
  explicit Expression(const Term& term);
  explicit Expression(Term&& term);
  explicit Expression(const container_type& container);
  explicit Expression(container_type&& container);
  explicit Expression(std::initializer_list<Term> lst);

  template<typename Container>
  Expression(complex_type c, Container&& ops) {
    hashmap.emplace(std::forward<Container>(ops), c);
  }

  Expression adjoint() const;

  size_t size() const { return hashmap.size(); }

  std::string to_string() const;

  void normalize();

  Expression& operator+=(const complex_type& value);
  Expression& operator-=(const complex_type& value);
  Expression& operator*=(const complex_type& value);
  Expression& operator/=(const complex_type& value);

  Expression& operator+=(const Expression& value);
  Expression& operator-=(const Expression& value);
  Expression& operator*=(const Expression& value);

  Expression& operator+=(const Term& value);
  Expression& operator-=(const Term& value);
  Expression& operator*=(const Term& value);
};

Expression operator+(Expression a, const Expression& b);
Expression operator-(Expression a, const Expression& b);
Expression operator*(Expression a, const Expression& b);

Expression operator+(Expression a, const Term& b);
Expression operator-(Expression a, const Term& b);
Expression operator*(Expression a, const Term& b);

Expression operator+(Expression a, const Expression::complex_type& b);
Expression operator-(Expression a, const Expression::complex_type& b);
Expression operator*(Expression a, const Expression::complex_type& b);
Expression operator/(Expression a, const Expression::complex_type& b);

Expression operator+(const Expression::complex_type& a, Expression b);
Expression operator-(const Expression::complex_type& a, Expression b);
Expression operator*(const Expression::complex_type& a, Expression b);

Expression operator+(const Term& a, Expression b);
Expression operator-(const Term& a, Expression b);
Expression operator*(const Term& a, Expression b);

Expression operator+(const Term& a, const Term& b);
Expression operator-(const Term& a, const Term& b);

Expression hopping(size_t from_orbital, size_t to_orbital, Operator::Spin spin);
Expression hopping(const Expression::complex_type& c, size_t from_orbital, size_t to_orbital,
                   Operator::Spin spin);

Expression spin_x(size_t i);
Expression spin_y(size_t i);
Expression spin_z(size_t i);
Expression spin_plus(size_t site);
Expression spin_minus(size_t site);
Expression spin_dot_product(size_t i, size_t j);
}  // namespace libqm
