#pragma once

#include <complex>

#include "Operator.h"
#include "StaticVector.h"

namespace libqm {
constexpr size_t cache_line_size = 64;  // NOTE: might be different

struct __attribute__((aligned(cache_line_size))) Term {
  using complex_type = std::complex<float>;

  static constexpr size_t static_vector_size =
      (cache_line_size - sizeof(complex_type)) / sizeof(Operator) - 1;

  using container_type = StaticVector<Operator, static_vector_size, uint8_t>;

  complex_type c{1.0f, 0.0f};
  container_type operators{};

  constexpr Term() noexcept = default;
  constexpr ~Term() noexcept = default;

  constexpr Term(const Term& other) noexcept = default;
  constexpr Term& operator=(const Term& other) noexcept = default;
  constexpr Term(Term&& other) noexcept = default;
  constexpr Term& operator=(Term&& other) noexcept = default;

  explicit constexpr Term(const Operator& x) noexcept : operators({x}) {}
  explicit constexpr Term(const complex_type& x) noexcept : c(x) {}
  explicit constexpr Term(const container_type& ops) noexcept : operators(ops) {}
  explicit constexpr Term(container_type&& ops) noexcept : operators(std::move(ops)) {}
  explicit constexpr Term(const complex_type& x, const container_type& ops) noexcept
      : c(x), operators(ops) {}
  explicit constexpr Term(const complex_type& x, container_type&& ops) noexcept
      : c(x), operators(std::move(ops)) {}
  explicit constexpr Term(std::initializer_list<Operator> init) noexcept : operators(init) {}
  explicit constexpr Term(const complex_type& x, std::initializer_list<Operator> init) noexcept
      : c(x), operators(init) {}

  constexpr size_t size() const noexcept { return operators.size(); }

  std::string to_string() const;

  constexpr bool operator==(const Term& other) const {
    return c == other.c && operators == other.operators;
  }

  constexpr Term adjoint() const noexcept {
    Term result(std::conj(c));
    for (size_t i = 0; i < size(); i++) {
      result.operators.push_back(operators[size() - i - 1].adjoint());
    }
    return result;
  }

  constexpr Term& operator*=(const Term& value) noexcept {
    c *= value.c;
    operators.append_range(value.operators.begin(), value.operators.end());
    return *this;
  }

  constexpr Term& operator*=(const Operator& value) noexcept {
    operators.push_back(value);
    return *this;
  }

  constexpr Term& operator*=(const complex_type& value) noexcept {
    c *= value;
    return *this;
  }

  constexpr Term& operator/=(const complex_type& value) noexcept {
    c /= value;
    return *this;
  }
};
static_assert(sizeof(Term) == 64, "Should fit a single cache line");

constexpr Term operator*(const Term& a, const Term& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

constexpr Term operator*(const Term& a, const Operator& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

constexpr Term operator*(const Term& a, const Term::complex_type& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

constexpr Term operator/(const Term& a, const Term::complex_type& b) noexcept {
  Term result(a);
  result /= b;
  return result;
}

constexpr Term operator*(const Term::complex_type& a, const Term& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

constexpr Term operator*(const Operator& a, const Term& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

bool is_diagonal(const Term::container_type& operators);

constexpr Term creation(Operator::Spin spin, size_t orbital) noexcept {
  return Term({Operator::creation(spin, orbital)});
}

constexpr Term annihilation(Operator::Spin spin, size_t orbital) noexcept {
  return Term({Operator::annihilation(spin, orbital)});
}

constexpr Term one_body(Operator::Spin s1, size_t o1, Operator::Spin s2, size_t o2) noexcept {
  return Term({Operator::creation(s1, o1), Operator::annihilation(s2, o2)});
}

constexpr Term density(Operator::Spin s, size_t o) noexcept {
  return Term({Operator::creation(s, o), Operator::annihilation(s, o)});
}

constexpr Term density_density(Operator::Spin s1, size_t i, Operator::Spin s2, size_t j) noexcept {
  return density(s1, i) * density(s2, j);
}
}  // namespace libqm
