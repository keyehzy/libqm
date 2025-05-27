#pragma once

#include "IndexedHashSet.h"
#include "Term.h"

namespace libqm {
struct Basis {
  using key_type = Term::container_type;
  using set_type = IndexedHashSet<key_type>;

  enum class Strategy {
    All,
    Restrict,
    Restrict_Sz,
    Restrict_Sz_P,
  };

  Basis() = default;

  Basis(size_t orbitals, size_t particles, Strategy strategy = Strategy::Restrict);

  void generate_all_combinations(key_type current, size_t first_orbital,
                                 std::vector<key_type>&) const;

  void generate_restrict_combinations(key_type current, size_t first_orbital,
                                      std::vector<key_type>&) const;

  void generate_restrict_sz_combinations(key_type current,
                                         size_t first_orbital,
                                         size_t remaining_up,
                                         size_t remaining_down,
                                         std::vector<key_type>&);

  void generate_restrict_sz_p_combinations(key_type current,
                                           size_t first_orbital,
                                           size_t remaining_up,
                                           size_t remaining_down,
                                           int target_K,
                                           int current_K_sum,
                                           int L,
                                           std::vector<key_type>&);

  template <typename Func>
  void sort(Func&& fn) {
    std::vector<key_type> elements = set.elements();
    std::sort(elements.begin(), elements.end(), fn);
    set = IndexedHashSet(std::move(elements));
  }

  template <typename Func>
  Basis refine(Func&& fn) {
    Basis basis;
    basis.orbitals = orbitals;
    basis.particles = particles;
    std::vector<key_type> elements;
    elements.reserve(set.elements().size());
    std::copy_if(set.elements().begin(), set.elements().end(), std::back_inserter(elements), fn);
    basis.set = IndexedHashSet(std::move(elements));
    return basis;
  }

  set_type set;
  size_t orbitals;
  size_t particles;
};

Basis generate_product_basis(const Basis& basis1, const Basis& basis2);
}  // namespace libqm
