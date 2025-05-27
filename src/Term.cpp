#include "Term.h"

#include <boost/unordered/unordered_flat_map.hpp>

#include <format>
#include <string>

namespace libqm {
std::string Term::to_string() const {
  // if (std::norm(c) < 1e-12) {
  //   return std::string();
  // }

  std::string result;
  for (const auto& op : operators) {
    result += op.to_string();
  }
  result += std::format(" ({:.5f},{:.5f})", c.real(), c.imag());

  return result;
}

bool is_diagonal(const Term::container_type& operators) {
  if (operators.empty()) {
    return true;
  }

  boost::unordered_flat_map<uint64_t, int> counts;
  for (const auto& op : operators) {
    counts[op.key()] += (op.type() == Operator::Type::Creation) ? 1 : -1;
  }

  for (const auto& pair : counts) {
    if (pair.second != 0) {
      return false;
    }
  }

  return true;
}

}  // namespace libqm
