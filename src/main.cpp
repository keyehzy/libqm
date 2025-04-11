#include "Basis.h"
#include "Expression.h"
#include "MatrixElements.h"
#include "NormalOrderer.h"
#include "eigen3/Eigen/Dense"

using namespace libqm;

int main() {
  Basis basis(10, 10, Basis::Strategy::All);
  std::cout << basis.set.size() << std::endl;
  for (const auto& term : basis.set) {
    std::cout << Term(term).to_string() << std::endl;
  }
  return 0;
}
