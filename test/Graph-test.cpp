#include "../src/Graph.h"

#include <gmock/gmock.h>  // For advanced matchers like UnorderedElementsAre
#include <gtest/gtest.h>

#include <vector>

namespace libqm::testing {

using ::testing::IsEmpty;
using ::testing::UnorderedElementsAre;

// --- Test Fixture for Graph Tests ---
class GraphTest : public ::testing::Test {
 protected:
  void SetUp() override {}
  void TearDown() override {}
};

// --- Constructor Tests ---

TEST_F(GraphTest, ConstructEmptyGraph) {
  Graph g(0);
  EXPECT_EQ(g.adj.size(), 0);
  EXPECT_TRUE(g.adj.empty());
}

TEST_F(GraphTest, ConstructSingleNodeGraph) {
  Graph g(1);
  ASSERT_EQ(g.adj.size(), 1);
  EXPECT_TRUE(g.adj[0].empty());
}

TEST_F(GraphTest, ConstructMultipleNodesGraph) {
  Graph g(5);
  ASSERT_EQ(g.adj.size(), 5);
  for (size_t i = 0; i < 5; ++i) {
    EXPECT_TRUE(g.adj[i].empty());
  }
}

TEST_F(GraphTest, AddEdgeBasic) {
  Graph g(3);
  g.add_edge(0, 1);

  ASSERT_EQ(g.adj.size(), 3);
  // Check neighbors of 0
  ASSERT_EQ(g.adj[0].size(), 1);
  EXPECT_EQ(g.adj[0][0], 1);
  // Check neighbors of 1 (undirected)
  ASSERT_EQ(g.adj[1].size(), 1);
  EXPECT_EQ(g.adj[1][0], 0);
  // Check neighbors of 2 (should be empty)
  EXPECT_TRUE(g.adj[2].empty());
}

TEST_F(GraphTest, AddMultipleEdges) {
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(0, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 3);

  // Use UnorderedElementsAre because the order within adjacency lists isn't guaranteed
  EXPECT_THAT(g.adj[0], UnorderedElementsAre(1, 2));
  EXPECT_THAT(g.adj[1], UnorderedElementsAre(0, 3));
  EXPECT_THAT(g.adj[2], UnorderedElementsAre(0, 3));
  EXPECT_THAT(g.adj[3], UnorderedElementsAre(1, 2));
}

TEST_F(GraphTest, AddSelfLoop) {
  Graph g(3);
  g.add_edge(1, 1);

  ASSERT_EQ(g.adj[1].size(), 2);  // Adds edge from 1 to 1 twice due to undirected nature
  EXPECT_THAT(g.adj[1], UnorderedElementsAre(1, 1));
  EXPECT_THAT(g.adj[0], IsEmpty());
  EXPECT_THAT(g.adj[2], IsEmpty());
}

TEST_F(GraphTest, AddDuplicateEdge) {
  Graph g(3);
  g.add_edge(0, 1);
  g.add_edge(1, 0);
  g.add_edge(0, 1);

  // Expect duplicates in the adjacency list based on current implementation
  EXPECT_THAT(g.adj[0], UnorderedElementsAre(1, 1, 1));
  EXPECT_THAT(g.adj[1], UnorderedElementsAre(0, 0, 0));
  EXPECT_THAT(g.adj[2], IsEmpty());
}

// --- BFS Tests ---
TEST_F(GraphTest, BFSOnSingleNodeGraph) {
  Graph g(1);
  auto visited = g.neighbors(0);

  ASSERT_EQ(visited.size(), 1);
  EXPECT_EQ(visited.count(0), 1);
}

TEST_F(GraphTest, BFSOnLineGraph) {
  // 0 -- 1 -- 2 -- 3
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  auto visited = g.neighbors(0);

  // Expected distances from node 0
  std::unordered_map<size_t, int> expected = {{0, 0}, {1, 1}, {2, 2}, {3, 3}};
  EXPECT_EQ(visited, expected);

  // Start from a middle node
  visited = g.neighbors(2);
  expected = {{2, 0}, {1, 1}, {3, 1}, {0, 2}};
  EXPECT_EQ(visited, expected);
}

TEST_F(GraphTest, BFSOnStarGraph) {
  // Center 0 connected to 1, 2, 3
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(0, 2);
  g.add_edge(0, 3);

  // Start from center
  auto visited = g.neighbors(0);
  std::unordered_map<size_t, int> expected_center = {{0, 0}, {1, 1}, {2, 1}, {3, 1}};
  EXPECT_EQ(visited, expected_center);

  // Start from a leaf
  visited = g.neighbors(1);
  std::unordered_map<size_t, int> expected_leaf = {{1, 0}, {0, 1}, {2, 2}, {3, 2}};
  EXPECT_EQ(visited, expected_leaf);
}

TEST_F(GraphTest, BFSOnCompleteGraphK3) {
  // 0 -- 1 -- 2 -- 0
  Graph g(3);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 0);

  auto visited = g.neighbors(0);
  std::unordered_map<size_t, int> expected = {{0, 0}, {1, 1}, {2, 1}};
  EXPECT_EQ(visited, expected);
}

TEST_F(GraphTest, BFSOnDisconnectedGraph) {
  // 0 -- 1   2 -- 3
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(2, 3);

  // Start in first component
  auto visited = g.neighbors(0);
  std::unordered_map<size_t, int> expected1 = {{0, 0}, {1, 1}};
  EXPECT_EQ(visited, expected1);  // Should only visit nodes 0 and 1

  // Start in second component
  visited = g.neighbors(3);
  std::unordered_map<size_t, int> expected2 = {{3, 0}, {2, 1}};
  EXPECT_EQ(visited, expected2);  // Should only visit nodes 2 and 3
}

TEST_F(GraphTest, BFSWithCycleAndMultiplePaths) {
  //   0 -- 1 -- 2
  //   |    |    |
  //   3 -- 4 -- 5
  Graph g(6);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(0, 3);
  g.add_edge(1, 4);
  g.add_edge(2, 5);
  g.add_edge(3, 4);
  g.add_edge(4, 5);

  auto visited = g.neighbors(0);
  std::unordered_map<size_t, int> expected = {{0, 0}, {1, 1}, {3, 1}, {2, 2}, {4, 2}, {5, 3}};
  // Node 4 can be reached via 0->1->4 (dist 2) or 0->3->4 (dist 2) - BFS finds shortest
  // Node 5 can be reached via 0->1->2->5 (dist 3) or 0->1->4->5 (dist 3) or 0->3->4->5 (dist 3)
  EXPECT_EQ(visited, expected);
}

// --- Neighbors Tests ---
TEST_F(GraphTest, NeighborsLevel0) {
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(1, 2);

  EXPECT_THAT(g.neighbors(0, 0), UnorderedElementsAre(0));
  EXPECT_THAT(g.neighbors(1, 0), UnorderedElementsAre(1));
  EXPECT_THAT(g.neighbors(3, 0), UnorderedElementsAre(3));  // Isolated node
}

TEST_F(GraphTest, NeighborsLevel1LineGraph) {
  // 0 -- 1 -- 2 -- 3
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  EXPECT_THAT(g.neighbors(0, 1), UnorderedElementsAre(1));
  EXPECT_THAT(g.neighbors(1, 1), UnorderedElementsAre(0, 2));
  EXPECT_THAT(g.neighbors(2, 1), UnorderedElementsAre(1, 3));
  EXPECT_THAT(g.neighbors(3, 1), UnorderedElementsAre(2));
}

TEST_F(GraphTest, NeighborsHigherLevelLineGraph) {
  // 0 -- 1 -- 2 -- 3
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  EXPECT_THAT(g.neighbors(0, 2), UnorderedElementsAre(2));
  EXPECT_THAT(g.neighbors(1, 2),
              UnorderedElementsAre(3));  // Only neighbor at distance 2 from 1 is 3
  EXPECT_THAT(g.neighbors(0, 3), UnorderedElementsAre(3));
}

TEST_F(GraphTest, NeighborsBeyondReach) {
  // 0 -- 1 -- 2
  Graph g(3);
  g.add_edge(0, 1);
  g.add_edge(1, 2);

  EXPECT_THAT(g.neighbors(0, 3), IsEmpty());
  EXPECT_THAT(g.neighbors(1, 2),
              UnorderedElementsAre());  // Correction: Node 1 has no neighbors at level 2
  EXPECT_THAT(g.neighbors(2, 5), IsEmpty());
}

TEST_F(GraphTest, NeighborsDisconnectedGraph) {
  // 0 -- 1   2 -- 3
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(2, 3);

  EXPECT_THAT(g.neighbors(0, 1), UnorderedElementsAre(1));
  EXPECT_THAT(g.neighbors(0, 2), IsEmpty());  // Cannot reach level 2 from 0
  EXPECT_THAT(g.neighbors(1, 1), UnorderedElementsAre(0));
  EXPECT_THAT(g.neighbors(2, 1), UnorderedElementsAre(3));
  EXPECT_THAT(g.neighbors(3, 2), IsEmpty());
}

TEST_F(GraphTest, NeighborsStarGraph) {
  // Center 0 connected to 1, 2, 3
  Graph g(4);
  g.add_edge(0, 1);
  g.add_edge(0, 2);
  g.add_edge(0, 3);

  // From Center
  EXPECT_THAT(g.neighbors(0, 0), UnorderedElementsAre(0));
  EXPECT_THAT(g.neighbors(0, 1), UnorderedElementsAre(1, 2, 3));
  EXPECT_THAT(g.neighbors(0, 2), IsEmpty());

  // From Leaf (e.g., node 1)
  EXPECT_THAT(g.neighbors(1, 0), UnorderedElementsAre(1));
  EXPECT_THAT(g.neighbors(1, 1), UnorderedElementsAre(0));
  EXPECT_THAT(g.neighbors(1, 2),
              UnorderedElementsAre(2, 3));  // Reach 2 via 1->0->2, reach 3 via 1->0->3
  EXPECT_THAT(g.neighbors(1, 3), IsEmpty());
}

TEST_F(GraphTest, NeighborsCompleteGraphK3) {
  // 0 -- 1 -- 2 -- 0
  Graph g(3);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 0);

  EXPECT_THAT(g.neighbors(0, 0), UnorderedElementsAre(0));
  EXPECT_THAT(g.neighbors(0, 1), UnorderedElementsAre(1, 2));
  EXPECT_THAT(g.neighbors(0, 2), IsEmpty());  // No nodes exactly 2 steps away in K3
}

TEST_F(GraphTest, NeighborsWithSelfLoop) {
  // 0 -- 1 -- 1 (self loop)
  Graph g(2);
  g.add_edge(0, 1);
  g.add_edge(1, 1);  // Adds edge 1->1 twice

  // BFS from 0: dist(0)=0, dist(1)=1
  // BFS from 1: dist(1)=0, dist(0)=1
  EXPECT_THAT(g.neighbors(0, 0), UnorderedElementsAre(0));
  EXPECT_THAT(g.neighbors(0, 1), UnorderedElementsAre(1));
  EXPECT_THAT(g.neighbors(0, 2), IsEmpty());

  EXPECT_THAT(g.neighbors(1, 0), UnorderedElementsAre(1));
  // Neighbors of 1 are 0 and 1 itself. Both are distance 1 according to BFS logic.
  // The self-loop doesn't change the distance calculation in this BFS implementation.
  EXPECT_THAT(
      g.neighbors(1, 1),
      UnorderedElementsAre(0));  // Only 0 is distance 1 from 1. Node 1 itself is distance 0.
  EXPECT_THAT(g.neighbors(1, 2), IsEmpty());
}

}  // namespace libqm::testing