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

Basis::Basis(size_t orbitals, size_t particles) : orbitals(orbitals), particles(particles) {
  LIBQM_ASSERT(orbitals <= Operator::max_index());
  LIBQM_ASSERT(particles <= 2 * orbitals);

  size_t basis_size = compute_basis_size(orbitals, particles);

  std::vector<key_type> acc;
  acc.reserve(basis_size);
  generate_basis(acc);
  LIBQM_ASSERT(acc.size() == basis_size);
  basis_set = IndexedHashSet(std::move(acc));
}

void Basis::generate_basis(std::vector<key_type>& acc) const {
  key_type current;
  generate_combinations(current, 0, 0, acc);
}

void Basis::generate_combinations(key_type& current, size_t first_orbital, size_t depth,
                                  std::vector<key_type>& acc) const {
  if (depth == particles) {
    acc.emplace_back(current);
    return;
  }

  for (size_t i = first_orbital; i < orbitals; i++) {
    for (int spin_index = 0; spin_index < 2; ++spin_index) {
      Operator::Spin spin = static_cast<Operator::Spin>(spin_index);
      bool should_iterate =
          current.empty() || (current.back().value() < i ||
                              (current.back().value() == i && spin > current.back().spin()));
      if (should_iterate) {
        current.push_back(Operator::creation(spin, i));
        generate_combinations(current, i, depth + 1, acc);
        current.pop_back();
      }
    }
  }
}
}  // namespace libqm
