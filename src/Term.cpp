#include "Term.h"

#include <format>
#include <string>

namespace libqm {
std::string Term::to_string() const {
  if (std::norm(c) < 1e-12) {
    return std::string();
  }

  std::string result = std::format("({:.5f},{:.5f})", c.real(), c.imag());
  for (const auto& op : operators) {
    result += op.to_string();
  }

  return result;
}

Term Term::adjoint() const noexcept {
  Term result(std::conj(c));
  for (size_t i = 0; i < size(); i++) {
    result.operators.push_back(operators[size() - i - 1].adjoint());
  }
  return result;
}

Term& Term::operator*=(const Term& value) noexcept {
  c *= value.c;
  operators.append_range(value.operators.begin(), value.operators.end());
  return *this;
}

Term& Term::operator*=(const Operator& value) noexcept {
  operators.push_back(value);
  return *this;
}

Term& Term::operator*=(const complex_type& value) noexcept {
  c *= value;
  return *this;
}

Term& Term::operator/=(const complex_type& value) noexcept {
  c /= value;
  return *this;
}

Term operator*(const Term& a, const Term& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

Term operator*(const Term& a, const Operator& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

Term operator*(const Term& a, const Term::complex_type& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

Term operator/(const Term& a, const Term::complex_type& b) noexcept {
  Term result(a);
  result /= b;
  return result;
}

Term operator*(const Term::complex_type& a, const Term& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

Term operator*(const Operator& a, const Term& b) noexcept {
  Term result(a);
  result *= b;
  return result;
}

bool is_diagonal(const Term::container_type& operators) {
  std::unordered_map<size_t, int> counts;

  for (const auto& op : operators) {
    auto key = op.data & ~Operator::kFermionTypeTagMask;
    counts[key] += (op.type() == Operator::Type::Creation) ? 1 : -1;
  }

  return std::all_of(counts.begin(), counts.end(),
                     [](const auto& element) { return element.second == 0; });
}

Term creation(Operator::Spin spin, size_t orbital) noexcept {
  return Term({Operator::creation(spin, orbital)});
}

Term annihilation(Operator::Spin spin, size_t orbital) noexcept {
  return Term({Operator::annihilation(spin, orbital)});
}

Term one_body(Operator::Spin spin1, size_t orbital1, Operator::Spin spin2,
              size_t orbital2) noexcept {
  return Term({Operator::creation(spin1, orbital1), Operator::annihilation(spin2, orbital2)});
}

Term density(Operator::Spin spin, size_t orbital) noexcept {
  return Term({Operator::creation(spin, orbital), Operator::annihilation(spin, orbital)});
}

Term density_density(Operator::Spin spin1, size_t i, Operator::Spin spin2, size_t j) noexcept {
  return density(spin1, i) * density(spin2, j);
}
}  // namespace libqm
