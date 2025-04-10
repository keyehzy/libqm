#pragma once

#include "IndexedHashSet.h"
#include "Term.h"

namespace libqm {
struct Basis {
 public:
  using key_type = Term::container_type;

  Basis(size_t orbitals, size_t particles);

  void generate_basis(std::vector<key_type>&) const;
  void generate_combinations(key_type& current, size_t first_orbital, size_t depth,
                             std::vector<key_type>&) const;

 private:
  IndexedHashSet<key_type> basis_set;
  size_t orbitals;
  size_t particles;
};
}  // namespace libqm
