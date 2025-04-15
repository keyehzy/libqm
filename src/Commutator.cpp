#include "Commutator.h"

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

Expression BCH(const Expression& A, const Expression& B, float lambda, NormalOrderer& orderer, size_t order) {
  std::vector<Expression> terms(order + 1);

  terms[0] = B;
  terms[1] = commutator(A, B, orderer);

  for (size_t n = 2; n <= order; ++n) {
    terms[n] = commutator(A, terms[n - 1], orderer);
  }

  Expression result;
  uint64_t factorial = 1;

  for (size_t n = 0; n <= order; ++n) {
    float coeff = std::powf(lambda, n) / factorial;
    result += terms[n] * coeff;
    factorial *= n + 1;
  }

  return orderer.normal_order(result);
}

}  // namespace libqm
