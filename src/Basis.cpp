#include "Basis.h"

#include <vector>

#include "Assert.h"

namespace libqm {
constexpr uint64_t choose(uint64_t n, uint64_t m) noexcept {
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

constexpr uint64_t compute_basis_size(uint64_t orbitals, uint64_t particles) noexcept {
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

    case Strategy::Restrict_Sz: {
      int required_sz = 0; // @@@
      int total = static_cast<int>(particles);
      LIBQM_ASSERT((total + required_sz) % 2 == 0 && (total - required_sz) % 2 == 0);

      size_t required_up = static_cast<size_t>((total + required_sz) / 2);
      size_t required_down = static_cast<size_t>((total - required_sz) / 2);
      LIBQM_ASSERT(required_up + required_down == particles);
      LIBQM_ASSERT(std::abs(required_sz) <= static_cast<int>(particles));

      size_t basis_size = choose(orbitals, required_up) * choose(orbitals, required_down);
      std::vector<key_type> acc;
      acc.reserve(basis_size);
      generate_restrict_sz_combinations({}, 0, required_up, required_down, acc);
      LIBQM_ASSERT(acc.size() == basis_size);
      set = IndexedHashSet(std::move(acc));
      break;
    }

    case Strategy::Restrict_Sz_P: {
      int required_sz = 0; // @@@
      int required_total_p = 0; // @@@
      int L = 14; // @@@
      int total = static_cast<int>(particles);
      LIBQM_ASSERT((total + required_sz) % 2 == 0 && (total - required_sz) % 2 == 0);

      size_t required_up = static_cast<size_t>((total + required_sz) / 2);
      size_t required_down = static_cast<size_t>((total - required_sz) / 2);
      LIBQM_ASSERT(required_up + required_down == particles);
      LIBQM_ASSERT(std::abs(required_sz) <= static_cast<int>(particles));

      // size_t basis_size = choose(orbitals, required_up) * choose(orbitals, required_down);
      std::vector<key_type> acc;
      // acc.reserve(basis_size);
      generate_restrict_sz_p_combinations({}, 0, required_up, required_down, required_total_p, 0, L, acc);
      // LIBQM_ASSERT(acc.size() == basis_size);
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
      // LIBQM_ASSERT(acc.size() == basis_size);
      std::sort(acc.begin(), acc.end(), [](const auto& a, const auto& b) { return a.size() < b.size(); });
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

void Basis::generate_restrict_sz_combinations(key_type current,
                                              size_t first_orbital,
                                              size_t remaining_up,
                                              size_t remaining_down,
                                              std::vector<key_type>& acc) {
  if (remaining_up == 0 && remaining_down == 0) {
    std::sort(current.begin(), current.end(), [](const auto& a, const auto& b) { return a < b; });
    acc.emplace_back(current);
    return;
  }

  for (size_t i = first_orbital; i < orbitals; i++) {
    if (remaining_up > 0) {
      Operator::Spin spin = Operator::Spin::Up;
      if (current.empty() || current.back().value() < i ||
          (current.back().value() == i && spin > current.back().spin())) {
        current.push_back(Operator::creation(spin, i));
        generate_restrict_sz_combinations(current, i, remaining_up - 1, remaining_down, acc);
        current.pop_back();
      }
    }
    if (remaining_down > 0) {
      Operator::Spin spin = Operator::Spin::Down;
      if (current.empty() || current.back().value() < i ||
          (current.back().value() == i && spin > current.back().spin())) {
        current.push_back(Operator::creation(spin, i));
        generate_restrict_sz_combinations(current, i, remaining_up, remaining_down - 1, acc);
        current.pop_back();
      }
    }
  }
}

void Basis::generate_restrict_sz_p_combinations(
    key_type current,
    size_t first_orbital,
    size_t remaining_up,
    size_t remaining_down,
    int target_K,
    int current_K_sum,
    int L,
    std::vector<key_type>& acc) {
  if (remaining_up == 0 && remaining_down == 0) {
    if (current_K_sum % L == target_K % L) {
      std::sort(current.begin(), current.end(),
                [](const auto& a, const auto& b) { return a < b; });
      acc.emplace_back(current);
    }
    return;
  }

  for (size_t i = first_orbital; i < orbitals; ++i) {
    if (remaining_up > 0) {
      Operator::Spin spin = Operator::Spin::Up;
      if (current.empty() || current.back().value() < i ||
          (current.back().value() == i && spin > current.back().spin())) {
        current.push_back(Operator::creation(spin, i));
        generate_restrict_sz_p_combinations(
            current,
            i,
            remaining_up - 1,
            remaining_down,
            target_K,
            current_K_sum + i,
            L,
            acc);
        current.pop_back();
      }
    }

    if (remaining_down > 0) {
      Operator::Spin spin = Operator::Spin::Down;
      if (current.empty() || current.back().value() < i ||
          (current.back().value() == i && spin > current.back().spin())) {
        current.push_back(Operator::creation(spin, i));
        generate_restrict_sz_p_combinations(
            current,
            i,
            remaining_up,
            remaining_down - 1,
            target_K,
            current_K_sum + i,
            L,
            acc);
        current.pop_back();
      }
    }
  }
}


Basis generate_product_basis(const Basis& basis1, const Basis& basis2) {
  Basis basis;

  std::vector<Basis::key_type> result;
  result.reserve(basis1.set.size() * basis2.set.size());
  for (const auto& b1 : basis1.set) {
    for (const auto& b2 : basis2.set) {
      result.emplace_back(b1);
      result.back().append_range_reverse(b2, [](const auto& a) { return a.adjoint(); });
    }
  }
  std::sort(result.begin(), result.end(),
            [](const auto& a, const auto& b) { return a.size() < b.size(); });
  basis.set = IndexedHashSet(std::move(result));
  return basis;
}

}  // namespace libqm
