#include "Expression.h"
#include "NormalOrderer.h"

using namespace libqm;

int main() {
  Term term(
      {Operator::annihilation(Operator::Spin::Up, 0), Operator::creation(Operator::Spin::Up, 0)});
  auto orderer = NormalOrderer();
  Expression e = orderer.normal_order(term);
  std::cout << e.to_string() << std::endl;
  orderer.print_cache_stats();
  return 0;
}