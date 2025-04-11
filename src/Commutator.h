#pragma once

#include "Expression.h"
#include "NormalOrderer.h"

namespace libqm {
Expression commutator(const Term& A, const Term& B, NormalOrderer& orderer);
Expression commutator(const Expression& A, const Expression& B, NormalOrderer& orderer);
Expression anticommutator(const Term& A, const Term& B, NormalOrderer& orderer);
Expression anticommutator(const Expression& A, const Expression& B, NormalOrderer& orderer);
}  // namespace libqm
