#pragma once

#include <cstddef>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Assert.h"

namespace libqm {
struct Graph {
  std::vector<std::vector<size_t>> adj;

  Graph(size_t n) : adj(n) {}

  void add_edge(size_t u, size_t v) {
    LIBQM_ASSERT(u < adj.size());
    LIBQM_ASSERT(v < adj.size());
    adj[u].push_back(v);
    adj[v].push_back(u);
  }

  template <typename Visitor>
  void bfs(size_t start_node, Visitor&& visitor) const {
    LIBQM_ASSERT(start_node < adj.size());
    LIBQM_ASSERT(!adj.empty());
    std::queue<std::pair<size_t, int>> q;
    std::unordered_set<size_t> visited;
    std::vector<int> distances(adj.size(), -1);

    q.push({start_node, 0});
    visited.insert(start_node);
    distances[start_node] = 0;

    while (!q.empty()) {
      auto [current_node, distance] = q.front();
      q.pop();

      visitor(current_node, distance);

      for (size_t neighbor : adj[current_node]) {
        if (!visited.count(neighbor)) {
          visited.insert(neighbor);
          q.push({neighbor, distances[current_node] + 1});
          distances[neighbor] = distances[current_node] + 1;
        }
      }
    }
  }

  std::unordered_map<size_t, int> neighbors(size_t start_node) const;
  std::unordered_set<size_t> neighbors(size_t start_node, size_t level) const;
};

}  // namespace libqm