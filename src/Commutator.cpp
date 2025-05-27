#include "Commutator.h"
#include "Logger.h"

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

  Expression BCH(const Expression& A, const Expression& B, float lambda, NormalOrderer& orderer, size_t order)
  {
    Expression current = B;
    float coeff = 1.0f;
    Expression result = current * coeff;

    for (size_t n = 1; n <= order; ++n) {
      current = commutator(A, current, orderer);
      coeff *= (lambda / static_cast<float>(n));
      result += current * coeff;
    }

    return orderer.normal_order(result);
  }

}  // namespace libqm
