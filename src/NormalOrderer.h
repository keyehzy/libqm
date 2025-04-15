#pragma once

#include <atomic>
#include <mutex>
#include <shared_mutex>

#include "Expression.h"

namespace libqm {
struct NormalOrderer {
  using complex_type = Term::complex_type;
  using container_type = Term::container_type;
  using cache_map_type = std::unordered_map<size_t, Expression>;

  NormalOrderer() = default;

  NormalOrderer(const NormalOrderer&) = delete;
  NormalOrderer& operator=(const NormalOrderer&) = delete;
  NormalOrderer(NormalOrderer&&) = delete;
  NormalOrderer& operator=(NormalOrderer&&) = delete;

  Expression normal_order(const complex_type& c, const container_type& ops);
  Expression normal_order(const Term& term);
  Expression normal_order(const Expression& expr);

  Expression normal_order_recursive(const container_type& ops);
  Expression normal_order_recursive(container_type ops, std::size_t ops_hash);
  Expression handle_non_commuting(container_type ops, size_t index);

  cache_map_type cache;

  void print_cache_stats() const;
  std::atomic<size_t> cache_hits{0};
  std::atomic<size_t> cache_misses{0};
  mutable std::shared_mutex cache_mutex;
};
}  // namespace libqm
