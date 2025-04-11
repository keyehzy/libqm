#include "Commutator.h"

#include "NormalOrderer.h"

namespace libqm {
Expression commutator(const Term& A, const Term& B, NormalOrderer& orderer) {
  return orderer.normal_order(A * B) - orderer.normal_order(B * A);
}

Expression commutator(const Expression& A, const Expression& B, NormalOrderer& orderer) {
  return orderer.normal_order(A * B) - orderer.normal_order(B * A);
}

Expression anticommutator(const Term& A, const Term& B, NormalOrderer& orderer) {
  return orderer.normal_order(A * B) + orderer.normal_order(B * A);
}

Expression anticommutator(const Expression& A, const Expression& B, NormalOrderer& orderer) {
  return orderer.normal_order(A * B) + orderer.normal_order(B * A);
}
}  // namespace libqm
