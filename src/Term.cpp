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
}  // namespace libqm
