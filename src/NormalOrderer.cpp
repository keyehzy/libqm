#include "NormalOrderer.h"

#include <limits>
#include <thread>

#include "Expression.h"

namespace libqm {
constexpr float tolerance = std::numeric_limits<float>::epsilon();

Expression NormalOrderer::normal_order(const complex_type& c, const container_type& ops) {
  if (std::norm(c) < tolerance * tolerance) {
    return {};
  }
  return c * normal_order_recursive(ops);
}

Expression NormalOrderer::normal_order(const Term& term) {
  return normal_order(term.c, term.operators);
}

Expression NormalOrderer::normal_order(const Expression& expr) {
  Expression result;
  for (const auto& [ops, c] : expr.hashmap) {
    result += normal_order(c, ops);
  }
  return result;
}

Expression NormalOrderer::normal_order_recursive(container_type ops) {
  if (ops.size() < 2) {
    return Expression(ops);
  }

  if (has_consecutive_elements(ops)) {
    return {};
  }

  Expression result;
  float phase = 1.0f;
  bool early_exit = false;
  size_t non_commuting_index = 0;

  for (size_t i = 1; i < ops.size(); ++i) {
    size_t j = i;
    while (j > 0 && ops[j] < ops[j - 1]) {
      if (ops[j].commutes(ops[j - 1])) {
        std::swap(ops[j], ops[j - 1]);
        phase = -phase;
        --j;
      } else {
        return phase * handle_non_commuting(std::move(ops), j-1);
      }
    }
  }

  if (has_consecutive_elements(ops)) {
    return {};
  }

  return Expression(phase, std::move(ops));
}

Expression NormalOrderer::handle_non_commuting(container_type ops, size_t index) {
  container_type contracted;
  contracted.append_range(ops.begin(), ops.begin() + index);
  contracted.append_range(ops.begin() + index + 2, ops.end());
  std::swap(ops[index], ops[index + 1]);
  return normal_order_recursive(std::move(contracted)) - normal_order_recursive(std::move(ops));
}

void NormalOrderer::print_cache_stats() const {
  std::cout << "Total entries: " << cache.size() << std::endl;
  std::cout << "Cache hits: " << cache_hits << std::endl;
  std::cout << "Cache misses: " << cache_misses << std::endl;
  double hit_ratio =
      static_cast<double>(cache_hits) / static_cast<double>(cache_hits + cache_misses);
  std::cout << "Hit ratio: " << hit_ratio << std::endl;
}
}  // namespace libqm
