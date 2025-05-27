#pragma once

#include <atomic>
#include <mutex>

#include "Expression.h"
#include "LruCache.h"
#include <boost/unordered/unordered_flat_map.hpp>

template <typename K, typename V, typename... Args>
class boost_unordered_flat_map_wrapper : public boost::unordered_flat_map<K, V, Args...> {};

namespace libqm {
struct NormalOrderer {
  using complex_type = Term::complex_type;
  using container_type = Term::container_type;
  using cache_map_type = LruCache<size_t, Expression, boost_unordered_flat_map_wrapper>;

  NormalOrderer() = default;

  Expression normal_order(const complex_type& c, const container_type& ops);
  Expression normal_order(const Term& term);
  Expression normal_order(const Expression& expr);

  Expression normal_order_recursive(container_type ops);
  Expression normal_order_recursive(container_type ops, std::size_t ops_hash);
  Expression handle_non_commuting(container_type ops, size_t index);

  cache_map_type cache = cache_map_type(1'000'000);

  void print_cache_stats() const;
  size_t cache_hits{0};
  size_t cache_misses{0};
};
}  // namespace libqm
