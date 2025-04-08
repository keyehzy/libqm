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

bool is_diagonal(const Term::container_type& operators) {
  std::unordered_map<size_t, int> counts;

  for (const auto& op : operators) {
    auto key = op.data & ~Operator::kFermionTypeTagMask;
    counts[key] += (op.type() == Operator::Type::Creation) ? 1 : -1;
  }

  return std::all_of(counts.begin(), counts.end(),
                     [](const auto& element) { return element.second == 0; });
}
}  // namespace libqm
