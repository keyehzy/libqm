#include "Basis.h"

#include <vector>

#include "Assert.h"

namespace libqm {
constexpr uint64_t choose(uint64_t n, uint64_t m) {
  if (m > n) return 0;
  if (m == 0 || m == n) return 1;

  if (m > n - m) m = n - m;

  uint64_t result = 1;
  for (uint64_t i = 0; i < m; ++i) {
    result *= (n - i);
    result /= (i + 1);
  }

  return result;
}

constexpr uint64_t compute_basis_size(uint64_t orbitals, uint64_t particles) {
  return choose(2 * orbitals, particles);
}

Basis::Basis(size_t orbitals, size_t particles, Strategy strategy)
    : orbitals(orbitals), particles(particles) {
  LIBQM_ASSERT(orbitals <= Operator::max_index());
  LIBQM_ASSERT(particles <= 2 * orbitals);

  switch (strategy) {
    case Strategy::Restrict: {
      size_t basis_size = compute_basis_size(orbitals, particles);
      std::vector<key_type> acc;
      acc.reserve(basis_size);
      generate_restrict_combinations({}, 0, acc);
      LIBQM_ASSERT(acc.size() == basis_size);
      set = IndexedHashSet(std::move(acc));
      break;
    }

    case Strategy::All: {
      size_t basis_size = 0;
      for (size_t i = 1; i <= particles; ++i) {
        basis_size += compute_basis_size(orbitals, i);
      }
      std::vector<key_type> acc;
      acc.reserve(basis_size);
      generate_all_combinations({}, 0, acc);
      LIBQM_ASSERT(acc.size() == basis_size);
      std::sort(acc.begin(), acc.end(),
                [](const auto& a, const auto& b) { return a.size() < b.size(); });
      set = IndexedHashSet(std::move(acc));
      break;
    }

    default:
      break;
  }
}

// NOTE: we can almost generate them in the correct order, but we still need to sort them
void Basis::generate_all_combinations(key_type current, size_t first_orbital,
                                      std::vector<key_type>& acc) const {
  if (current.size() > 0 && current.size() <= particles) {
    std::sort(current.begin(), current.end(), [](const auto& a, const auto& b) { return a < b; });
    acc.emplace_back(current);
  }
  for (int spin_index = 0; spin_index < 2; ++spin_index) {
    for (size_t i = first_orbital; i < orbitals; i++) {
      Operator::Spin spin = static_cast<Operator::Spin>(spin_index);
      bool should_iterate =
          current.empty() || (current.back().value() < i ||
                              (current.back().value() == i && spin > current.back().spin()));
      if (should_iterate) {
        current.push_back(Operator::creation(spin, i));
        generate_all_combinations(current, i, acc);
        current.pop_back();
      }
    }
  }
}

void Basis::generate_restrict_combinations(key_type current, size_t first_orbital,
                                           std::vector<key_type>& acc) const {
  if (current.size() == particles) {
    std::sort(current.begin(), current.end(), [](const auto& a, const auto& b) { return a < b; });
    acc.emplace_back(current);
    return;
  }
  for (int spin_index = 0; spin_index < 2; ++spin_index) {
    for (size_t i = first_orbital; i < orbitals; i++) {
      Operator::Spin spin = static_cast<Operator::Spin>(spin_index);
      bool should_iterate =
          current.empty() || (current.back().value() < i ||
                              (current.back().value() == i && spin > current.back().spin()));
      if (should_iterate) {
        current.push_back(Operator::creation(spin, i));
        generate_restrict_combinations(current, i, acc);
        current.pop_back();
      }
    }
  }
}
}  // namespace libqm
