#include "Graph.h"

namespace libqm {
std::unordered_map<size_t, int> Graph::neighbors(size_t start_node) const {
  std::unordered_map<size_t, int> result;
  bfs(start_node, [&](size_t node, int distance) { result.insert({node, distance}); });
  return result;
}

std::unordered_set<size_t> Graph::neighbors(size_t start_node, size_t level) const {
  std::unordered_set<size_t> result;
  bfs(start_node, [&](size_t node, int distance) {
    if (static_cast<size_t>(distance) == level) {
      result.insert(node);
    }
  });
  return result;
}
}  // namespace libqm