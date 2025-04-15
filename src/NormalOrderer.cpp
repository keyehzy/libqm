#include "NormalOrderer.h"

#include <limits>
#include <thread>

#include "Expression.h"

namespace libqm {
constexpr float tolerance = std::numeric_limits<float>::epsilon();

Expression NormalOrderer::normal_order(const complex_type& c, const container_type& ops) {
  if (std::norm(c) < tolerance * tolerance) {
    return Expression();
  }
  auto result = c * normal_order_recursive(ops);
  return result;
}

Expression NormalOrderer::normal_order(const Term& term) {
  Expression result = normal_order(term.c, term.operators);
  return result;
}

Expression NormalOrderer::normal_order(const Expression& expr) {
  size_t num_tasks = expr.hashmap.size();

  if (num_tasks == 0) {
    return Expression();
  }

  std::vector<Expression> results(num_tasks);
  std::vector<std::thread> threads;
  threads.reserve(num_tasks);

  size_t current_task_index = 0;
  for (const auto& [ops_key, c_value] : expr.hashmap) {
    Expression& result_slot = results[current_task_index];
    threads.emplace_back([this, ops = ops_key, c = c_value, &result_slot]() {
      result_slot = this->normal_order(c, ops);
    });
    current_task_index++;
  }

  for (std::thread& t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  Expression result;
  for (const auto& individual_result : results) {
    result += individual_result;
  }

  return result;
}

Expression NormalOrderer::normal_order_recursive(const container_type& ops) {
  if (ops.size() < 2) {
    return Expression(ops);
  }

  if (has_consecutive_elements(ops)) {
    return Expression();
  }

  return normal_order_recursive(ops, std::hash<container_type>()(ops));
}

Expression NormalOrderer::normal_order_recursive(container_type ops, size_t hash) {
  {
    std::shared_lock reader_lock(cache_mutex);
    auto it = cache.find(hash);
    if (it != cache.end()) {
      cache_hits++;
      return it->second;
    }
    cache_misses++;
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
        early_exit = true;
        non_commuting_index = j - 1;
        goto handle_computation;
      }
    }
  }

handle_computation:
  if (early_exit) {
    result = static_cast<complex_type>(phase) *
             handle_non_commuting(std::move(ops), non_commuting_index);
  } else {
    if (has_consecutive_elements(ops)) {
      result = Expression();
    } else {
      result = Expression(phase, std::move(ops));
    }
  }

  std::unique_lock writer_lock(cache_mutex);
  cache.emplace(hash, result);
  return result;
}

Expression NormalOrderer::handle_non_commuting(container_type ops, size_t index) {
  container_type contracted;
  contracted.append_range(ops.begin(), ops.begin() + index);
  contracted.append_range(ops.begin() + index + 2, ops.end());
  std::swap(ops[index], ops[index + 1]);
  Expression lhs = normal_order_recursive(std::move(contracted));
  Expression rhs = normal_order_recursive(std::move(ops));
  return lhs - rhs;
}

void NormalOrderer::print_cache_stats() const {
  std::shared_lock reader_lock(cache_mutex);
  std::cout << "Total entries: " << cache.size() << std::endl;
  std::cout << "Cache hits: " << cache_hits << std::endl;
  std::cout << "Cache misses: " << cache_misses << std::endl;
  double hit_ratio =
      static_cast<double>(cache_hits) / static_cast<double>(cache_hits + cache_misses);
  std::cout << "Hit ratio: " << hit_ratio << std::endl;
}
}  // namespace libqm
