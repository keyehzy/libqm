#include "Operator.h"

#include <format>

namespace libqm {
std::string Operator::to_string() const {
  return std::format("c{}({},{})", (type() == Type::Creation) ? "+" : "",
                     (spin() == Spin::Up) ? "↑" : "↓", value());
}
}  // namespace libqm
