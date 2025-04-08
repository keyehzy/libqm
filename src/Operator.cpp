#include "Operator.h"

namespace libqm {
std::string Operator::to_string() const {
  const Type type_val = type();
  const Spin spin_val = spin();
  const ubyte value_val = value();

  std::string result;
  result.reserve(4096);

  result += 'c';
  if (type_val == Type::Creation) {
    result += '+';
  }
  result += '(';
  result += (spin_val == Spin::Up) ? "↑" : "↓";
  result += ',';

  char num_buffer[256];
  int len = std::sprintf(num_buffer, "%u", value_val);
  result.append(num_buffer, len);

  result += ')';

  return result;
}
}  // namespace libqm
