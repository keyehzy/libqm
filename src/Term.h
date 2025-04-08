#pragma once

#include <complex>

#include "Operator.h"
#include "StaticVector.h"

namespace libqm {
constexpr size_t cache_line_size = 64;  // NOTE: might be different

struct __attribute__((aligned(cache_line_size))) Term {
  using complex_type = std::complex<float>;

  using container_type =
      StaticVector<Operator, (cache_line_size - sizeof(complex_type)) / sizeof(Operator) - 1, uint8_t>;

  complex_type c{1.0f, 0.0f};
  container_type operators;

  Term() noexcept = default;
  ~Term() noexcept = default;

  Term(const Term& other) noexcept = default;
  Term& operator=(const Term& other) noexcept = default;
  Term(Term&& other) noexcept = default;
  Term& operator=(Term&& other) noexcept = default;

  explicit Term(const Operator& x) noexcept : operators({x}) {}
  explicit Term(const complex_type& x) noexcept : c(x) {}
  explicit Term(const container_type& ops) noexcept : operators(ops) {}
  explicit Term(container_type&& ops) noexcept : operators(std::move(ops)) {}
  explicit Term(const complex_type& x, const container_type& ops) noexcept : c(x), operators(ops) {}
  explicit Term(const complex_type& x, container_type&& ops) noexcept
      : c(x), operators(std::move(ops)) {}
  explicit Term(const std::initializer_list<Operator>& ops) noexcept : operators(ops) {}
  explicit Term(const complex_type& x, const std::initializer_list<Operator>& ops) noexcept
      : c(x), operators(ops) {}

  size_t size() const noexcept { return operators.size(); }

  std::string to_string() const;

  bool operator==(const Term& other) const {
    return c == other.c && operators == other.operators;
  }

  Term adjoint() const noexcept;

  Term& operator*=(const complex_type& value) noexcept;
  Term& operator*=(const Term& value) noexcept;
  Term& operator*=(const Operator& value) noexcept;
  Term& operator/=(const complex_type& value) noexcept;
};
static_assert(sizeof(Term::container_type) == 56);
static_assert(sizeof(Term) == 64, "Should fit a single cache line");

Term operator*(const Term& a, const Term& b) noexcept;
Term operator*(const Term& a, const Operator& b) noexcept;
Term operator*(const Term& a, const Term::complex_type& b) noexcept;
Term operator/(const Term& a, const Term::complex_type& b) noexcept;
Term operator*(const Term::complex_type& a, const Term& b) noexcept;
Term operator*(const Operator& a, const Term& b) noexcept;

bool is_diagonal(const Term::container_type& operators);

Term creation(Operator::Spin spin, size_t orbital) noexcept;
Term annihilation(Operator::Spin spin, size_t orbital) noexcept;
Term one_body(Operator::Spin spin1, size_t orbital1, Operator::Spin spin2,
              size_t orbital2) noexcept;
Term density(Operator::Spin spin, size_t orbital) noexcept;
Term density_density(Operator::Spin spin1, size_t i, Operator::Spin spin2, size_t j) noexcept;
}  // namespace libqm
