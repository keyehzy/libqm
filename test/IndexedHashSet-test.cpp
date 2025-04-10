#include "../src/IndexedHashSet.h"

#include <iterator>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"

template <typename Set>
void ExpectSetContents(const Set& s, const std::vector<typename Set::key_type>& expected) {
  EXPECT_EQ(s.size(), expected.size());
  EXPECT_EQ(s.empty(), expected.empty());

  std::vector<typename Set::key_type> actual_forward(s.begin(), s.end());
  EXPECT_EQ(actual_forward, expected);

  std::vector<typename Set::key_type> actual_const_forward(s.cbegin(), s.cend());
  EXPECT_EQ(actual_const_forward, expected);

  std::vector<typename Set::key_type> actual_reverse(s.crbegin(), s.crend());
  std::vector<typename Set::key_type> expected_reverse = expected;
  std::reverse(expected_reverse.begin(), expected_reverse.end());
  EXPECT_EQ(actual_reverse, expected_reverse);

  for (size_t i = 0; i < expected.size(); ++i) {
    EXPECT_TRUE(s.contains(expected[i]));
    EXPECT_EQ(s.at(i), expected[i]);
    EXPECT_EQ(s[i], expected[i]);
    EXPECT_EQ(s.index_of(expected[i]), i);
  }

  if constexpr (std::is_integral_v<typename Set::key_type>) {
    if (!expected.empty()) {
      typename Set::key_type non_existent = expected.back() + 1;

      while (s.contains(non_existent) && non_existent != expected.back()) {
        ++non_existent;
      }
      if (!s.contains(non_existent)) {
        EXPECT_FALSE(s.contains(non_existent));
        EXPECT_THROW(s.index_of(non_existent), std::out_of_range);
      }
    } else {
      EXPECT_FALSE(s.contains(0));
      EXPECT_THROW(s.index_of(0), std::out_of_range);
    }
  } else if constexpr (std::is_same_v<typename Set::key_type, std::string>) {
    std::string non_existent = "a value highly unlikely to be present";
    EXPECT_FALSE(s.contains(non_existent));
    EXPECT_THROW(s.index_of(non_existent), std::out_of_range);
  }

  EXPECT_THROW(s.at(s.size()), std::out_of_range);
}

class IndexedHashSetTest : public ::testing::Test {
 protected:
};

TEST_F(IndexedHashSetTest, DefaultConstructor) {
  IndexedHashSet<int> s;
  EXPECT_TRUE(s.empty());
  EXPECT_EQ(s.size(), 0);
  EXPECT_EQ(s.begin(), s.end());
  EXPECT_EQ(s.cbegin(), s.cend());
  EXPECT_EQ(s.crbegin(), s.crend());
  ExpectSetContents(s, {});
}

TEST_F(IndexedHashSetTest, InitializerListConstructorEmpty) {
  IndexedHashSet<std::string> s = {};
  EXPECT_TRUE(s.empty());
  EXPECT_EQ(s.size(), 0);
  ExpectSetContents(s, {});
}

TEST_F(IndexedHashSetTest, InitializerListConstructorUnique) {
  IndexedHashSet<int> s({10, 20, 5, 30});
  std::vector<int> expected = {10, 20, 5, 30};
  ExpectSetContents(s, expected);
}

TEST_F(IndexedHashSetTest, InitializerListConstructorStrings) {
  IndexedHashSet<std::string> s({"apple", "banana", "cherry"});
  std::vector<std::string> expected = {"apple", "banana", "cherry"};
  ExpectSetContents(s, expected);
}

TEST_F(IndexedHashSetTest, IteratorRangeConstructorVector) {
  std::vector<int> data = {5, 15, 25};
  IndexedHashSet<int> s(data.begin(), data.end());
  ExpectSetContents(s, data);
}

TEST_F(IndexedHashSetTest, IteratorRangeConstructorList) {
  std::list<std::string> data = {"one", "two", "three"};
  IndexedHashSet<std::string> s(data.begin(), data.end());
  std::vector<std::string> expected = {"one", "two", "three"};
  ExpectSetContents(s, expected);
}

TEST_F(IndexedHashSetTest, IteratorRangeConstructorEmptyRange) {
  std::vector<int> data = {};
  IndexedHashSet<int> s(data.begin(), data.end());
  ExpectSetContents(s, {});
}

TEST_F(IndexedHashSetTest, IteratorRangeConstructorInputIterator) {
  std::string data_str = "10 20 30 40";
  std::istringstream iss(data_str);
  std::istream_iterator<int> begin(iss);
  std::istream_iterator<int> end;  // Default constructor creates end iterator

  // Note: Input iterators can only be traversed once.
  // The IndexedHashSet constructor reads them.
  IndexedHashSet<int> s(begin, end);

  std::vector<int> expected = {10, 20, 30, 40};
  ExpectSetContents(s, expected);
}

TEST_F(IndexedHashSetTest, VectorCopyConstructor) {
  std::vector<std::string> data = {"x", "y", "z"};
  std::vector<std::string> data_copy = data;  // Keep original safe
  IndexedHashSet<std::string> s(data);
  ExpectSetContents(s, data_copy);
  EXPECT_EQ(data, data_copy);  // Ensure original vector wasn't modified
}

TEST_F(IndexedHashSetTest, VectorCopyConstructorEmpty) {
  std::vector<int> data = {};
  IndexedHashSet<int> s(data);
  ExpectSetContents(s, {});
}

TEST_F(IndexedHashSetTest, VectorMoveConstructor) {
  std::vector<std::string> data = {"move", "construct", "test"};
  std::vector<std::string> expected = data;
  IndexedHashSet<std::string> s(std::move(data));  // data is now in a valid but unspecified state
  ExpectSetContents(s, expected);
  // Cannot reliably test the state of 'data' after move as per standard
}

TEST_F(IndexedHashSetTest, VectorMoveConstructorEmpty) {
  std::vector<int> data = {};
  std::vector<int> expected = {};
  IndexedHashSet<int> s(std::move(data));
  ExpectSetContents(s, expected);
}

TEST_F(IndexedHashSetTest, CopyConstructor) {
  IndexedHashSet<int> s1({1, 2, 3});
  IndexedHashSet<int> s2 = s1;  // Copy construct

  std::vector<int> expected = {1, 2, 3};
  ExpectSetContents(s1, expected);  // Original should be unchanged
  ExpectSetContents(s2, expected);  // Copy should match
  EXPECT_TRUE(s1 == s2);
}

TEST_F(IndexedHashSetTest, CopyConstructorEmpty) {
  IndexedHashSet<int> s1;
  IndexedHashSet<int> s2 = s1;  // Copy construct

  ExpectSetContents(s1, {});
  ExpectSetContents(s2, {});
  EXPECT_TRUE(s1 == s2);
}

TEST_F(IndexedHashSetTest, MoveConstructor) {
  IndexedHashSet<int> s1({10, 20, 30});
  std::vector<int> expected = {10, 20, 30};
  IndexedHashSet<int> s2 = std::move(s1);  // Move construct

  ExpectSetContents(s2, expected);
  // s1 is now in a valid but unspecified state (likely empty)
  EXPECT_TRUE(s1.empty());  // Check common post-move state
  EXPECT_EQ(s1.size(), 0);
}

TEST_F(IndexedHashSetTest, MoveConstructorEmpty) {
  IndexedHashSet<int> s1;
  IndexedHashSet<int> s2 = std::move(s1);  // Move construct empty

  ExpectSetContents(s2, {});
  EXPECT_TRUE(s1.empty());
  EXPECT_EQ(s1.size(), 0);
}

TEST_F(IndexedHashSetTest, CopyAssignment) {
  IndexedHashSet<int> s1({1, 2, 3});
  IndexedHashSet<int> s2({10, 20});
  std::vector<int> expected1 = {1, 2, 3};
  std::vector<int> expected2 = {10, 20};

  ExpectSetContents(s1, expected1);
  ExpectSetContents(s2, expected2);

  s2 = s1;  // Copy assign

  ExpectSetContents(s1, expected1);  // Source unchanged
  ExpectSetContents(s2, expected1);  // Destination updated
  EXPECT_TRUE(s1 == s2);
}

TEST_F(IndexedHashSetTest, CopyAssignmentToEmpty) {
  IndexedHashSet<int> s1({1, 2, 3});
  IndexedHashSet<int> s2;
  std::vector<int> expected1 = {1, 2, 3};

  s2 = s1;

  ExpectSetContents(s1, expected1);
  ExpectSetContents(s2, expected1);
  EXPECT_TRUE(s1 == s2);
}

TEST_F(IndexedHashSetTest, CopyAssignmentFromEmpty) {
  IndexedHashSet<int> s1;
  IndexedHashSet<int> s2({1, 2, 3});
  std::vector<int> expected2 = {1, 2, 3};

  ExpectSetContents(s1, {});
  ExpectSetContents(s2, expected2);

  s2 = s1;

  ExpectSetContents(s1, {});
  ExpectSetContents(s2, {});
  EXPECT_TRUE(s1 == s2);
}

TEST_F(IndexedHashSetTest, CopyAssignmentSelf) {
  IndexedHashSet<int> s1({1, 2, 3});
  std::vector<int> expected1 = {1, 2, 3};
  s1 = s1;  // Self assignment
  ExpectSetContents(s1, expected1);
}

TEST_F(IndexedHashSetTest, MoveAssignment) {
  IndexedHashSet<std::string> s1({"a", "b"});
  IndexedHashSet<std::string> s2({"x", "y", "z"});
  std::vector<std::string> expected1 = {"a", "b"};
  std::vector<std::string> expected2 = {"x", "y", "z"};

  ExpectSetContents(s1, expected1);
  ExpectSetContents(s2, expected2);

  s2 = std::move(s1);  // Move assign

  ExpectSetContents(s2, expected1);  // Destination has source's content
  // s1 is now in a valid but unspecified state (likely empty)
  EXPECT_TRUE(s1.empty());  // Check common post-move state
  EXPECT_EQ(s1.size(), 0);
}

TEST_F(IndexedHashSetTest, MoveAssignmentToEmpty) {
  IndexedHashSet<int> s1({1, 2, 3});
  IndexedHashSet<int> s2;
  std::vector<int> expected1 = {1, 2, 3};

  s2 = std::move(s1);

  ExpectSetContents(s2, expected1);
  EXPECT_TRUE(s1.empty());
  EXPECT_EQ(s1.size(), 0);
}

TEST_F(IndexedHashSetTest, MoveAssignmentFromEmpty) {
  IndexedHashSet<int> s1;
  IndexedHashSet<int> s2({1, 2, 3});
  std::vector<int> expected2 = {1, 2, 3};

  ExpectSetContents(s1, {});
  ExpectSetContents(s2, expected2);

  s2 = std::move(s1);

  ExpectSetContents(s2, {});
  EXPECT_TRUE(s1.empty());  // Source remains empty
  EXPECT_EQ(s1.size(), 0);
}

TEST_F(IndexedHashSetTest, MoveAssignmentSelf) {
  IndexedHashSet<int> s1({1, 2, 3});
  std::vector<int> expected1 = {1, 2, 3};
  // Move self-assignment is tricky and potentially dangerous,
  // but the default implementation should handle it correctly (might be no-op or clear).
  // We test that it doesn't crash and the object remains valid.
  s1 = std::move(s1);
  // The state after move self-assignment is often unspecified, but it shouldn't crash.
  // It *might* still hold the original elements or be empty.
  // Let's check if it matches the original state OR is empty, as both are plausible outcomes
  // depending on the exact std::vector/std::unordered_map move assignment behavior.
  bool is_original = (s1.size() == 3 && s1[0] == 1 && s1[1] == 2 && s1[2] == 3);
  bool is_empty = s1.empty();
  EXPECT_TRUE(is_original || is_empty);  // Should be in one of these valid states
}

TEST_F(IndexedHashSetTest, EmptyAndSize) {
  IndexedHashSet<int> s_empty;
  EXPECT_TRUE(s_empty.empty());
  EXPECT_EQ(s_empty.size(), 0);

  IndexedHashSet<int> s_one({42});
  EXPECT_FALSE(s_one.empty());
  EXPECT_EQ(s_one.size(), 1);

  IndexedHashSet<std::string> s_many({"a", "b", "c", "d"});
  EXPECT_FALSE(s_many.empty());
  EXPECT_EQ(s_many.size(), 4);
}

TEST_F(IndexedHashSetTest, MaxSize) {
  IndexedHashSet<int> s;
  // max_size depends on system constraints, just check it's positive
  EXPECT_GT(s.max_size(), 0);

  // Check that it's related to vector and map max_sizes
  std::vector<int> v;
  std::unordered_map<int, size_t> m;
  EXPECT_LE(s.max_size(), v.max_size());
  EXPECT_LE(s.max_size(), m.max_size());
}

TEST_F(IndexedHashSetTest, Contains) {
  IndexedHashSet<std::string> s({"cat", "dog", "bird"});
  EXPECT_TRUE(s.contains("cat"));
  EXPECT_TRUE(s.contains("dog"));
  EXPECT_TRUE(s.contains("bird"));
  EXPECT_FALSE(s.contains("fish"));
  EXPECT_FALSE(s.contains(""));     // Empty string
  EXPECT_FALSE(s.contains("CAT"));  // Case-sensitive
}

TEST_F(IndexedHashSetTest, ContainsEmptySet) {
  IndexedHashSet<int> s;
  EXPECT_FALSE(s.contains(0));
  EXPECT_FALSE(s.contains(1));
}

TEST_F(IndexedHashSetTest, At) {
  IndexedHashSet<int> s({100, 200, 300});
  EXPECT_EQ(s.at(0), 100);
  EXPECT_EQ(s.at(1), 200);
  EXPECT_EQ(s.at(2), 300);

  const auto& cs = s;
  EXPECT_EQ(cs.at(0), 100);
  EXPECT_EQ(cs.at(1), 200);
  EXPECT_EQ(cs.at(2), 300);
}

TEST_F(IndexedHashSetTest, AtThrowsOutOfRange) {
  IndexedHashSet<int> s({100, 200});
  EXPECT_THROW(s.at(2), std::out_of_range);
  EXPECT_THROW(s.at(100), std::out_of_range);

  const auto& cs = s;
  EXPECT_THROW(cs.at(2), std::out_of_range);
}

TEST_F(IndexedHashSetTest, AtEmptySetThrows) {
  IndexedHashSet<int> s;
  EXPECT_THROW(s.at(0), std::out_of_range);
  const auto& cs = s;
  EXPECT_THROW(cs.at(0), std::out_of_range);
}

TEST_F(IndexedHashSetTest, OperatorSquareBrackets) {
  IndexedHashSet<int> s({100, 200, 300});
  EXPECT_EQ(s[0], 100);
  EXPECT_EQ(s[1], 200);
  EXPECT_EQ(s[2], 300);

  const auto& cs = s;
  EXPECT_EQ(cs[0], 100);
  EXPECT_EQ(cs[1], 200);
  EXPECT_EQ(cs[2], 300);
}

TEST_F(IndexedHashSetTest, IndexOf) {
  IndexedHashSet<std::string> s({"first", "second", "third"});
  EXPECT_EQ(s.index_of("first"), 0);
  EXPECT_EQ(s.index_of("second"), 1);
  EXPECT_EQ(s.index_of("third"), 2);

  const auto& cs = s;
  EXPECT_EQ(cs.index_of("first"), 0);
  EXPECT_EQ(cs.index_of("second"), 1);
  EXPECT_EQ(cs.index_of("third"), 2);
}

TEST_F(IndexedHashSetTest, IndexOfThrowsOutOfRange) {
  IndexedHashSet<std::string> s({"first", "second"});
  EXPECT_THROW(s.index_of("third"), std::out_of_range);
  EXPECT_THROW(s.index_of(""), std::out_of_range);

  const auto& cs = s;
  EXPECT_THROW(cs.index_of("third"), std::out_of_range);
}

TEST_F(IndexedHashSetTest, IndexOfEmptySetThrows) {
  IndexedHashSet<std::string> s;
  EXPECT_THROW(s.index_of("anything"), std::out_of_range);
  const auto& cs = s;
  EXPECT_THROW(cs.index_of("anything"), std::out_of_range);
}

TEST_F(IndexedHashSetTest, IteratorsForward) {
  IndexedHashSet<int> s({1, 2, 3});
  std::vector<int> result;
  for (auto it = s.begin(); it != s.end(); ++it) {
    result.push_back(*it);
  }
  EXPECT_EQ(result, std::vector<int>({1, 2, 3}));

  result.clear();
  for (auto it = s.cbegin(); it != s.cend(); ++it) {
    result.push_back(*it);
  }
  EXPECT_EQ(result, std::vector<int>({1, 2, 3}));

  result.clear();
  for (const auto& val : s) {
    result.push_back(val);
  }
  EXPECT_EQ(result, std::vector<int>({1, 2, 3}));
}

TEST_F(IndexedHashSetTest, IteratorsReverse) {
  IndexedHashSet<int> s({1, 2, 3});
  std::vector<int> result;
  // Note: The class provides const_iterator for rbegin/rend via vector's crbegin/crend
  for (auto it = s.crbegin(); it != s.crend(); ++it) {
    result.push_back(*it);
  }
  EXPECT_EQ(result, std::vector<int>({3, 2, 1}));
}

TEST_F(IndexedHashSetTest, IteratorsEmptySet) {
  IndexedHashSet<int> s;
  EXPECT_EQ(s.begin(), s.end());
  EXPECT_EQ(s.cbegin(), s.cend());
  EXPECT_EQ(s.crbegin(), s.crend());
}

// --- Comparison Operator Tests ---

TEST_F(IndexedHashSetTest, EqualityOperator) {
  IndexedHashSet<int> s1({1, 2, 3});
  IndexedHashSet<int> s2({1, 2, 3});
  IndexedHashSet<int> s3({1, 3, 2});
  IndexedHashSet<int> s4({1, 2});
  IndexedHashSet<int> s5({1, 2, 4});
  IndexedHashSet<int> e1;
  IndexedHashSet<int> e2;

  EXPECT_TRUE(s1 == s2);
  EXPECT_FALSE(s1 == s3);
  EXPECT_FALSE(s1 == s4);
  EXPECT_FALSE(s1 == s5);
  EXPECT_FALSE(s1 == e1);
  EXPECT_TRUE(e1 == e2);
  EXPECT_FALSE(e1 == s1);
}

TEST_F(IndexedHashSetTest, InequalityOperator) {
  IndexedHashSet<int> s1({1, 2, 3});
  IndexedHashSet<int> s2({1, 2, 3});
  IndexedHashSet<int> s3({1, 3, 2});
  IndexedHashSet<int> s4({1, 2});
  IndexedHashSet<int> s5({1, 2, 4});
  IndexedHashSet<int> e1;
  IndexedHashSet<int> e2;

  EXPECT_FALSE(s1 != s2);
  EXPECT_TRUE(s1 != s3);
  EXPECT_TRUE(s1 != s5);
  EXPECT_TRUE(s1 != e1);
  EXPECT_FALSE(e1 != e2);
  EXPECT_TRUE(e1 != s1);

  EXPECT_FALSE(s1 != s1);
  EXPECT_FALSE(e1 != e1);
}

struct Point {
  int x, y;

  Point() : x(0), y(0) {}
  Point(int x_val, int y_val) : x(x_val), y(y_val) {}

  bool operator==(const Point& other) const { return x == other.x && y == other.y; }
};

struct PointHash {
  std::size_t operator()(const Point& p) const {
    std::size_t h1 = std::hash<int>{}(p.x);
    std::size_t h2 = std::hash<int>{}(p.y);
    return h1 ^ (h2 << 1);
  }
};

struct PointEqual {
  bool operator()(const Point& p1, const Point& p2) const { return p1.x == p2.x && p1.y == p2.y; }
};

TEST_F(IndexedHashSetTest, CustomTypeWithHashAndEqual) {
  using PointSet = IndexedHashSet<Point, PointHash, PointEqual>;

  Point p1{1, 2}, p2{3, 4}, p3{1, 2};

  PointSet s({p1, p2});
  std::vector<Point> expected = {p1, p2};

  EXPECT_EQ(s.size(), 2);
  EXPECT_TRUE(s.contains(p1));
  EXPECT_TRUE(s.contains(Point{1, 2}));
  EXPECT_TRUE(s.contains(p2));
  EXPECT_FALSE(s.contains({5, 6}));

  EXPECT_EQ(s.index_of(p1), 0);
  EXPECT_EQ(s.index_of(p2), 1);
  EXPECT_EQ(s[0], p1);
  EXPECT_EQ(s[1], p2);

  std::vector<Point> actual(s.begin(), s.end());
  EXPECT_EQ(actual.size(), 2);
  EXPECT_EQ(actual[0], p1);
  EXPECT_EQ(actual[1], p2);
}

TEST_F(IndexedHashSetTest, DeductionGuideInitializerList) {
  IndexedHashSet s({1, 2, 3});
  static_assert(std::is_same_v<decltype(s)::key_type, int>);
  EXPECT_EQ(s.size(), 3);
  EXPECT_EQ(s[0], 1);
}

TEST_F(IndexedHashSetTest, DeductionGuideIterators) {
  std::vector<std::string> v = {"a", "b"};
  IndexedHashSet s(v.begin(), v.end());
  static_assert(std::is_same_v<decltype(s)::key_type, std::string>);
  EXPECT_EQ(s.size(), 2);
  EXPECT_EQ(s[0], "a");
}

TEST_F(IndexedHashSetTest, DeductionGuideVectorCopy) {
  std::vector<double> v = {1.1, 2.2};
  IndexedHashSet s(v);
  static_assert(std::is_same_v<decltype(s)::key_type, double>);
  EXPECT_EQ(s.size(), 2);
  EXPECT_EQ(s[0], 1.1);
}

TEST_F(IndexedHashSetTest, DeductionGuideVectorMove) {
  std::vector<char> v = {'x', 'y'};
  IndexedHashSet s(std::move(v));
  static_assert(std::is_same_v<decltype(s)::key_type, char>);
  EXPECT_EQ(s.size(), 2);
  EXPECT_EQ(s[0], 'x');
}
