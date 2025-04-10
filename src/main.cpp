#include "Basis.h"
#include "Expression.h"
#include "NormalOrderer.h"
#include "MatrixElements.h"

#include "eigen3/Eigen/Dense"

using namespace libqm;

int main() {
  Basis basis(2, 1);

  for (const auto& term : basis.set) {
    std::cout << Term(term).to_string() << std::endl;
  }

  Expression hop = hopping(1.0f, 0, 1, Operator::Spin::Up);
  NormalOrderer orderer;
  auto matrix = compute_matrix_elements_serial<Eigen::MatrixXcd>(basis, hop, orderer);
  std::cout << matrix << std::endl;
  return 0;
}
