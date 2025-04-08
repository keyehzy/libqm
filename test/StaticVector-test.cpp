#include "../src/StaticVector.h"

#include <array>
#include <compare>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::ContainerEq;
using ::testing::ElementsAre;
using ::testing::ElementsAreArray;

class StaticVectorTest : public ::testing::Test {
 protected:
  using VecInt5 = StaticVector<int, 5, uint8_t>;
  using VecDouble3 = StaticVector<double, 3, uint8_t>;
  using VecChar10 = StaticVector<char, 10, uint8_t>;
  using VecUInt16_20 = StaticVector<unsigned, 20, uint16_t>;

  template <typename T, typename U>
  struct Pair {
    T first;
    U second;

    constexpr Pair() = default;
    constexpr Pair(T first, U second) : first(first), second(second) {}

    constexpr bool operator==(const Pair& other) const = default;
  };

  using VecPoint4 = StaticVector<Pair<int, int>, 4>;
};

TEST_F(StaticVectorTest, DefaultConstructor) {
  constexpr VecInt5 vec{};

  static_assert(vec.empty());
  static_assert(vec.size() == 0);
  static_assert(vec.capacity() == 5);

  EXPECT_TRUE(vec.empty());
  EXPECT_EQ(vec.size(), 0);
  EXPECT_EQ(vec.capacity(), 5);
  EXPECT_EQ(vec.max_size(), 5);
  EXPECT_EQ(vec.remaining_capacity(), 5);
  EXPECT_FALSE(vec.full());
}

TEST_F(StaticVectorTest, FillConstructor) {
  constexpr VecInt5 vec(3, 42);
  static_assert(vec.size() == 3);
  static_assert(!vec.empty());
  static_assert(vec[0] == 42 && vec[1] == 42 && vec[2] == 42);

  EXPECT_THAT(vec, ElementsAre(42, 42, 42));
  EXPECT_EQ(vec.size(), 3);
  EXPECT_FALSE(vec.empty());
  EXPECT_FALSE(vec.full());
}

TEST_F(StaticVectorTest, FillConstructorFull) {
  constexpr VecInt5 vec(5, -1);
  static_assert(vec.size() == 5);
  static_assert(vec.full());

  EXPECT_THAT(vec, ElementsAre(-1, -1, -1, -1, -1));
  EXPECT_EQ(vec.size(), 5);
  EXPECT_TRUE(vec.full());
}

TEST_F(StaticVectorTest, FillConstructorEmpty) {
  constexpr VecInt5 vec(0, 99);
  static_assert(vec.empty());

  EXPECT_TRUE(vec.empty());
  EXPECT_EQ(vec.size(), 0);
}

TEST_F(StaticVectorTest, InitializerListConstructor) {
  constexpr VecInt5 vec = {10, 20, 30};
  static_assert(vec.size() == 3);
  static_assert(vec[0] == 10 && vec[1] == 20 && vec[2] == 30);

  EXPECT_THAT(vec, ElementsAre(10, 20, 30));
  EXPECT_EQ(vec.size(), 3);
}

TEST_F(StaticVectorTest, InitializerListConstructorFull) {
  constexpr VecInt5 vec = {1, 2, 3, 4, 5};
  static_assert(vec.size() == 5);
  static_assert(vec.full());

  EXPECT_THAT(vec, ElementsAre(1, 2, 3, 4, 5));
  EXPECT_EQ(vec.size(), 5);
  EXPECT_TRUE(vec.full());
}

TEST_F(StaticVectorTest, InitializerListConstructorEmpty) {
  constexpr VecInt5 vec = {};
  static_assert(vec.empty());

  EXPECT_TRUE(vec.empty());
  EXPECT_EQ(vec.size(), 0);
}

TEST_F(StaticVectorTest, CopyConstructor) {
  VecInt5 vec1 = {1, 2, 3};
  VecInt5 vec2 = vec1;  // Copy construct

  EXPECT_THAT(vec2, ContainerEq(vec1));
  EXPECT_THAT(vec2, ElementsAre(1, 2, 3));

  EXPECT_THAT(vec1, ElementsAre(1, 2, 3));
}

TEST_F(StaticVectorTest, CopyConstructorEmpty) {
  VecInt5 vec1;
  VecInt5 vec2 = vec1;  // Copy construct empty

  EXPECT_TRUE(vec1.empty());
  EXPECT_TRUE(vec2.empty());
  EXPECT_THAT(vec2, ContainerEq(vec1));
}

TEST_F(StaticVectorTest, MoveConstructor) {
  VecInt5 vec1 = {1, 2, 3};
  const VecInt5 vec_copy = vec1;   // Keep a copy for comparison
  VecInt5 vec2 = std::move(vec1);  // Move construct

  // Check moved-to vector
  EXPECT_THAT(vec2, ContainerEq(vec_copy));
  EXPECT_THAT(vec2, ElementsAre(1, 2, 3));

  // Check moved-from vector (should be empty as per implementation)
  EXPECT_TRUE(vec1.empty());
  EXPECT_EQ(vec1.size(), 0);
}

TEST_F(StaticVectorTest, MoveConstructorEmpty) {
  VecInt5 vec1;
  VecInt5 vec2 = std::move(vec1);  // Move construct empty

  EXPECT_TRUE(vec1.empty());  // Original state after move (empty)
  EXPECT_TRUE(vec2.empty());  // New vector is empty
  EXPECT_THAT(vec1, ContainerEq(vec2));
}

TEST_F(StaticVectorTest, CopyAssignment) {
  VecInt5 vec1 = {1, 2, 3};
  VecInt5 vec2 = {99, 88};

  vec2 = vec1;  // Copy assign

  EXPECT_THAT(vec2, ContainerEq(vec1));
  EXPECT_THAT(vec2, ElementsAre(1, 2, 3));
  EXPECT_THAT(vec1, ElementsAre(1, 2, 3));  // Original unchanged

  // Assign empty
  VecInt5 vec3 = {1, 2, 3, 4, 5};
  const VecInt5 vec_empty;
  vec3 = vec_empty;
  EXPECT_TRUE(vec3.empty());
  EXPECT_THAT(vec3, ContainerEq(vec_empty));

  // Self-assignment
  VecInt5 vec4 = {10, 20};
  const VecInt5 vec4_copy = vec4;
  vec4 = vec4;  // Should have no effect due to self-assignment check
  EXPECT_THAT(vec4, ContainerEq(vec4_copy));
  EXPECT_THAT(vec4, ElementsAre(10, 20));
}

TEST_F(StaticVectorTest, MoveAssignment) {
  VecInt5 vec1 = {1, 2, 3};
  VecInt5 vec2 = {99, 88};
  const VecInt5 vec_copy = vec1;  // Keep a copy for comparison

  vec2 = std::move(vec1);  // Move assign

  // Check moved-to vector
  EXPECT_THAT(vec2, ContainerEq(vec_copy));
  EXPECT_THAT(vec2, ElementsAre(1, 2, 3));

  // Check moved-from vector (should be empty)
  EXPECT_TRUE(vec1.empty());

  // Move assign empty
  VecInt5 vec3 = {1, 2, 3, 4, 5};
  VecInt5 vec_empty;
  vec3 = std::move(vec_empty);
  EXPECT_TRUE(vec3.empty());
  EXPECT_TRUE(vec_empty.empty());  // Source is also empty

  // Self-assignment (move)
  VecInt5 vec4 = {10, 20};
  const VecInt5 vec4_copy = vec4;
  VecInt5* p_vec4 = &vec4;
  vec4 = std::move(*p_vec4);  // Move self-assignment
  // Implementation checks for self-assignment, so it should be a no-op
  EXPECT_THAT(vec4, ContainerEq(vec4_copy));
  EXPECT_THAT(vec4, ElementsAre(10, 20));
}

TEST_F(StaticVectorTest, InitializerListAssignment) {
  VecInt5 vec = {1, 2, 3, 4, 5};

  vec = {10, 20};  // Assign shorter list
  EXPECT_THAT(vec, ElementsAre(10, 20));
  EXPECT_EQ(vec.size(), 2);

  vec = {100, 200, 300, 400, 500};  // Assign full list
  EXPECT_THAT(vec, ElementsAre(100, 200, 300, 400, 500));
  EXPECT_TRUE(vec.full());

  vec = {};  // Assign empty list
  EXPECT_TRUE(vec.empty());
}

TEST_F(StaticVectorTest, AccessOperatorBracket) {
  constexpr VecInt5 vec = {10, 20, 30};
  static_assert(vec[0] == 10);
  static_assert(vec[1] == 20);
  static_assert(vec[2] == 30);

  EXPECT_EQ(vec[0], 10);
  EXPECT_EQ(vec[1], 20);
  EXPECT_EQ(vec[2], 30);

  VecInt5 vec_mut = {1, 2};
  vec_mut[0] = 11;
  vec_mut[1] = 22;
  EXPECT_EQ(vec_mut[0], 11);
  EXPECT_EQ(vec_mut[1], 22);
}

TEST_F(StaticVectorTest, AccessAt) {
  constexpr VecInt5 vec = {10, 20, 30};
  static_assert(vec.at(0) == 10);

  EXPECT_EQ(vec.at(0), 10);
  EXPECT_EQ(vec.at(1), 20);
  EXPECT_EQ(vec.at(2), 30);

  VecInt5 vec_mut = {1, 2};
  vec_mut.at(0) = 11;
  vec_mut.at(1) = 22;
  EXPECT_EQ(vec_mut.at(0), 11);
  EXPECT_EQ(vec_mut.at(1), 22);

  const VecInt5& cvec = vec_mut;
  EXPECT_EQ(cvec.at(0), 11);
  EXPECT_EQ(cvec.at(1), 22);
}

TEST_F(StaticVectorTest, Front) {
  constexpr VecInt5 vec = {10, 20};
  static_assert(vec.front() == 10);
  EXPECT_EQ(vec.front(), 10);

  VecInt5 vec_mut = {1, 2};
  vec_mut.front() = 11;
  EXPECT_EQ(vec_mut.front(), 11);
  EXPECT_EQ(vec_mut[0], 11);

  const VecInt5& cvec = vec_mut;
  EXPECT_EQ(cvec.front(), 11);
}

TEST_F(StaticVectorTest, Back) {
  constexpr VecInt5 vec = {10, 20, 30};
  static_assert(vec.back() == 30);
  EXPECT_EQ(vec.back(), 30);

  VecInt5 vec_mut = {1, 2};
  vec_mut.back() = 22;
  EXPECT_EQ(vec_mut.back(), 22);
  EXPECT_EQ(vec_mut[1], 22);

  vec_mut.push_back(33);
  EXPECT_EQ(vec_mut.back(), 33);

  const VecInt5& cvec = vec_mut;
  EXPECT_EQ(cvec.back(), 33);
}

TEST_F(StaticVectorTest, Data) {
  VecInt5 vec = {10, 20, 30};
  int* ptr = vec.data();
  EXPECT_EQ(*ptr, 10);
  EXPECT_EQ(ptr[1], 20);
  ptr[0] = 11;
  EXPECT_EQ(vec[0], 11);

  const VecInt5& cvec = vec;
  const int* cptr = cvec.data();
  EXPECT_EQ(*cptr, 11);
  EXPECT_EQ(cptr[1], 20);
}

TEST_F(StaticVectorTest, GetSpan) {
  VecInt5 vec = {10, 20, 30};
  std::span<int> s = vec.get_span();
  // Verify span contents
  EXPECT_EQ(s.size(), 3);
  EXPECT_EQ(s[0], 10);
  EXPECT_EQ(s[1], 20);
  EXPECT_EQ(s[2], 30);
  s[0] = 11;
  EXPECT_EQ(vec[0], 11);

  const VecInt5& cvec = vec;
  std::span<const int> cs = cvec.get_span();
  EXPECT_EQ(cs.size(), 3);
  EXPECT_EQ(cs[0], 11);
  EXPECT_EQ(cs[1], 20);
  EXPECT_EQ(cs[2], 30);
}

TEST_F(StaticVectorTest, IteratorsBeginEnd) {
  VecInt5 vec = {10, 20, 30};
  const VecInt5& cvec = vec;

  EXPECT_EQ(*vec.begin(), 10);
  EXPECT_EQ(*(vec.end() - 1), 30);
  EXPECT_EQ(std::distance(vec.begin(), vec.end()), 3);
  *vec.begin() = 11;
  EXPECT_EQ(vec[0], 11);
  EXPECT_EQ(vec[1], 20);
  EXPECT_EQ(vec[2], 30);

  EXPECT_EQ(*cvec.begin(), 11);
  EXPECT_EQ(*cvec.cbegin(), 11);
  EXPECT_EQ(std::distance(cvec.cbegin(), cvec.cend()), 3);
}

TEST_F(StaticVectorTest, IteratorsEmpty) {
  VecInt5 vec;
  const VecInt5& cvec = vec;

  EXPECT_EQ(vec.begin(), vec.end());
  EXPECT_EQ(cvec.begin(), cvec.end());
  EXPECT_EQ(cvec.cbegin(), cvec.cend());
  EXPECT_EQ(std::distance(vec.begin(), vec.end()), 0);
}

TEST_F(StaticVectorTest, IteratorLoop) {
  VecInt5 vec = {1, 2, 3, 4, 5};
  int sum = 0;
  for (int x : vec) {
    sum += x;
  }
  EXPECT_EQ(sum, 15);

  for (int& x : vec) {
    x *= 2;
  }
  // Keep ElementsAre for concise verification of all elements
  EXPECT_THAT(vec, ElementsAre(2, 4, 6, 8, 10));

  const VecInt5& cvec = vec;
  sum = 0;
  for (const int& x : cvec) {
    sum += x;
  }
  EXPECT_EQ(sum, 30);
}

TEST_F(StaticVectorTest, StdAlgorithmsWithIterators) {
  VecInt5 vec = {5, 1, 4, 2, 3};

  std::sort(vec.begin(), vec.end());
  // Verify sorted state concisely
  EXPECT_THAT(vec, ElementsAre(1, 2, 3, 4, 5));

  int sum = std::accumulate(vec.cbegin(), vec.cend(), 0);
  EXPECT_EQ(sum, 15);

  auto it = std::find(vec.begin(), vec.end(), 4);
  ASSERT_NE(it, vec.end());
  EXPECT_EQ(*it, 4);
  EXPECT_EQ(std::distance(vec.begin(), it), 3);

  it = std::find(vec.begin(), vec.end(), 99);
  EXPECT_EQ(it, vec.end());
}

TEST_F(StaticVectorTest, CapacityMethods) {
  VecInt5 vec;
  EXPECT_TRUE(vec.empty());
  EXPECT_EQ(vec.size(), 0);
  EXPECT_EQ(vec.capacity(), 5);
  EXPECT_EQ(vec.remaining_capacity(), 5);
  EXPECT_FALSE(vec.full());

  vec.push_back(1);
  vec.push_back(2);
  vec.push_back(3);
  EXPECT_FALSE(vec.empty());
  EXPECT_EQ(vec.size(), 3);
  EXPECT_EQ(vec.remaining_capacity(), 2);
  EXPECT_FALSE(vec.full());

  vec.push_back(4);
  vec.push_back(5);
  EXPECT_EQ(vec.size(), 5);
  EXPECT_EQ(vec.remaining_capacity(), 0);
  EXPECT_TRUE(vec.full());
}

TEST_F(StaticVectorTest, Clear) {
  VecInt5 vec = {1, 2, 3};
  EXPECT_FALSE(vec.empty());

  vec.clear();
  EXPECT_TRUE(vec.empty());
  EXPECT_EQ(vec.size(), 0);
  EXPECT_EQ(vec.capacity(), 5);

  vec.clear();  // Clear already empty
  EXPECT_TRUE(vec.empty());
}

TEST_F(StaticVectorTest, PushBack) {
  VecInt5 vec;
  int v1 = 10;

  vec.push_back(v1);
  EXPECT_THAT(vec, ElementsAre(10));
  EXPECT_EQ(vec.size(), 1);

  vec.push_back(20);
  EXPECT_THAT(vec, ElementsAre(10, 20));
  EXPECT_EQ(vec.size(), 2);

  vec.push_back(30);
  vec.push_back(40);
  vec.push_back(50);
  EXPECT_THAT(vec, ElementsAre(10, 20, 30, 40, 50));
  EXPECT_EQ(vec.size(), 5);
  EXPECT_TRUE(vec.full());
}

TEST_F(StaticVectorTest, EmplaceBack) {
  VecPoint4 vec;

  auto& p1_ref = vec.emplace_back(1, 10);
  EXPECT_THAT(vec, ElementsAre(Pair{1, 10}));
  EXPECT_EQ(vec.size(), 1);
  EXPECT_EQ(&p1_ref, &vec.back());

  Pair<int, int> p_src{2, 20};
  auto& p2_ref = vec.emplace_back(p_src);
  EXPECT_THAT(vec, ElementsAre(Pair{1, 10}, Pair{2, 20}));
  EXPECT_EQ(vec.size(), 2);
  EXPECT_EQ(&p2_ref, &vec.back());

  auto& p3_ref = vec.emplace_back(3, 30);
  auto& p4_ref = vec.emplace_back(4, 40);
  EXPECT_THAT(vec, ElementsAre(Pair{1, 10}, Pair{2, 20}, Pair{3, 30}, Pair{4, 40}));
  EXPECT_EQ(vec.size(), 4);
  EXPECT_TRUE(vec.full());
  EXPECT_EQ(&p4_ref, &vec.back());
}

TEST_F(StaticVectorTest, PopBack) {
  VecInt5 vec = {1, 2, 3, 4, 5};

  vec.pop_back();
  EXPECT_THAT(vec, ElementsAre(1, 2, 3, 4));
  EXPECT_EQ(vec.size(), 4);

  vec.pop_back();
  vec.pop_back();
  EXPECT_THAT(vec, ElementsAre(1, 2));
  EXPECT_EQ(vec.size(), 2);

  vec.pop_back();
  vec.pop_back();
  EXPECT_TRUE(vec.empty());
  EXPECT_EQ(vec.size(), 0);
}

TEST_F(StaticVectorTest, ResizeCount) {
  VecInt5 vec = {1, 2, 3};

  vec.resize(2);
  EXPECT_THAT(vec, ElementsAre(1, 2));
  EXPECT_EQ(vec.size(), 2);

  vec.resize(4);
  EXPECT_EQ(vec.size(), 4);
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[1], 2);

  vec.resize(1);
  EXPECT_THAT(vec, ElementsAre(1));

  vec.resize(0);
  EXPECT_TRUE(vec.empty());

  vec.resize(3);
  EXPECT_EQ(vec.size(), 3);
}

TEST_F(StaticVectorTest, ResizeCountValue) {
  VecInt5 vec = {1, 2, 3};

  vec.resize(2, 99);
  EXPECT_THAT(vec, ElementsAre(1, 2));

  vec.resize(4, 42);
  EXPECT_THAT(vec, ElementsAre(1, 2, 42, 42));

  vec.resize(5, 55);
}

using TestVector = StaticVector<int, 10, uint8_t>;
using ConstTestVector = const StaticVector<int, 10, uint8_t>;

class StaticVectorReverseIteratorTest : public ::testing::Test {
 protected:
  static void SetUpTestSuite() {}

  static void TearDownTestSuite() {}

  void SetUp() override {
    sv_multi_ = {10, 20, 30, 40, 50};
    expected_multi_reverse_ = {50, 40, 30, 20, 10};
  }

  void TearDown() override {}

  TestVector sv_multi_;
  std::vector<int> expected_multi_reverse_;
};

TEST_F(StaticVectorReverseIteratorTest, EmptyVector) {
  TestVector sv_empty;
  ConstTestVector& csv_empty = sv_empty;

  EXPECT_EQ(sv_empty.rbegin(), sv_empty.rend());
  EXPECT_EQ(sv_empty.crbegin(), sv_empty.crend());
  EXPECT_EQ(csv_empty.rbegin(), csv_empty.rend());
  EXPECT_EQ(csv_empty.crbegin(), csv_empty.crend());

  EXPECT_EQ(sv_empty.rbegin().base(), sv_empty.end());
  EXPECT_EQ(sv_empty.rend().base(), sv_empty.begin());
  EXPECT_EQ(csv_empty.rbegin().base(), csv_empty.end());
  EXPECT_EQ(csv_empty.rend().base(), csv_empty.begin());
  EXPECT_EQ(sv_empty.crbegin().base(), sv_empty.cend());
  EXPECT_EQ(sv_empty.crend().base(), sv_empty.cbegin());
}

TEST_F(StaticVectorReverseIteratorTest, SingleElementVector) {
  TestVector sv_one = {42};
  ConstTestVector& csv_one = sv_one;

  EXPECT_NE(sv_one.rbegin(), sv_one.rend());
  EXPECT_NE(sv_one.crbegin(), sv_one.crend());
  EXPECT_NE(csv_one.rbegin(), csv_one.rend());
  EXPECT_NE(csv_one.crbegin(), csv_one.crend());

  EXPECT_EQ(*sv_one.rbegin(), 42);
  EXPECT_EQ(*csv_one.rbegin(), 42);
  EXPECT_EQ(*sv_one.crbegin(), 42);
  EXPECT_EQ(*csv_one.crbegin(), 42);

  auto rit = sv_one.rbegin();
  ASSERT_NE(rit, sv_one.rend());
  EXPECT_EQ(*rit, 42);
  ++rit;
  EXPECT_EQ(rit, sv_one.rend());

  auto crit = csv_one.crbegin();
  ASSERT_NE(crit, csv_one.crend());
  EXPECT_EQ(*crit, 42);
  ++crit;
  EXPECT_EQ(crit, csv_one.crend());
}

TEST_F(StaticVectorReverseIteratorTest, MultiElementIteration) {
  ConstTestVector& csv_multi = sv_multi_;

  std::vector<int> actual_reverse;

  actual_reverse.clear();
  std::copy(sv_multi_.rbegin(), sv_multi_.rend(), std::back_inserter(actual_reverse));
  EXPECT_EQ(actual_reverse, expected_multi_reverse_);
  EXPECT_THAT(actual_reverse, ::testing::ElementsAre(50, 40, 30, 20, 10));

  actual_reverse.clear();
  std::copy(csv_multi.rbegin(), csv_multi.rend(), std::back_inserter(actual_reverse));
  EXPECT_EQ(actual_reverse, expected_multi_reverse_);
  actual_reverse.clear();
  std::copy(csv_multi.rbegin(), csv_multi.rend(), std::back_inserter(actual_reverse));
  EXPECT_EQ(actual_reverse, expected_multi_reverse_);
  EXPECT_THAT(actual_reverse, ::testing::ElementsAre(50, 40, 30, 20, 10));

  actual_reverse.clear();
  std::copy(sv_multi_.crbegin(), sv_multi_.crend(), std::back_inserter(actual_reverse));
  EXPECT_EQ(actual_reverse, expected_multi_reverse_);
  EXPECT_THAT(actual_reverse, ::testing::ElementsAre(50, 40, 30, 20, 10));

  actual_reverse.clear();
  std::copy(csv_multi.crbegin(), csv_multi.crend(), std::back_inserter(actual_reverse));
  EXPECT_EQ(actual_reverse, expected_multi_reverse_);
  EXPECT_THAT(actual_reverse, ::testing::ElementsAre(50, 40, 30, 20, 10));
}

TEST_F(StaticVectorReverseIteratorTest, MultiElementModification) {
  ASSERT_FALSE(sv_multi_.empty());

  *sv_multi_.rbegin() = 55;
  EXPECT_EQ(sv_multi_.back(), 55);
  EXPECT_EQ(sv_multi_[4], 55);

  ASSERT_NE(sv_multi_.rbegin(), sv_multi_.rend());
  auto it = sv_multi_.rbegin();
  ++it;
  ASSERT_NE(it, sv_multi_.rend());
  *it = 44;
  EXPECT_EQ(sv_multi_[3], 44);

  EXPECT_THAT(sv_multi_, ::testing::ElementsAre(10, 20, 30, 44, 55));
}

TEST_F(StaticVectorReverseIteratorTest, MultiElementAlgorithms) {
  std::vector<int> result;
  result.resize(sv_multi_.size());

  std::reverse_copy(sv_multi_.begin(), sv_multi_.end(), result.begin());
  EXPECT_EQ(result, expected_multi_reverse_);
  EXPECT_THAT(result, ::testing::ElementsAre(50, 40, 30, 20, 10));

  std::fill(result.begin(), result.end(), 0);
  std::copy(sv_multi_.rbegin(), sv_multi_.rend(), result.begin());
  EXPECT_EQ(result, expected_multi_reverse_);
  EXPECT_THAT(result, ::testing::ElementsAre(50, 40, 30, 20, 10));
}

TEST_F(StaticVectorReverseIteratorTest, MultiElementBaseAndArithmetic) {
  ConstTestVector& csv_multi = sv_multi_;

  ASSERT_FALSE(sv_multi_.empty());
  EXPECT_EQ(*(sv_multi_.rbegin()), sv_multi_.back());
  EXPECT_EQ(*(csv_multi.rbegin()), csv_multi.back());
  EXPECT_EQ(*(sv_multi_.crbegin()), sv_multi_.back());

  EXPECT_EQ(sv_multi_.rbegin().base(), sv_multi_.end());
  EXPECT_EQ(csv_multi.rbegin().base(), csv_multi.end());
  EXPECT_EQ(sv_multi_.crbegin().base(), sv_multi_.cend());

  EXPECT_EQ(sv_multi_.rend().base(), sv_multi_.begin());
  EXPECT_EQ(csv_multi.rend().base(), csv_multi.begin());
  EXPECT_EQ(sv_multi_.crend().base(), sv_multi_.cbegin());

  EXPECT_EQ(*(sv_multi_.rbegin().base() - 1), *sv_multi_.rbegin());
  EXPECT_EQ(*(csv_multi.rbegin().base() - 1), *csv_multi.rbegin());
  EXPECT_EQ(*(sv_multi_.crbegin().base() - 1), *sv_multi_.crbegin());
  auto rit1 = sv_multi_.rbegin();
  auto crit1 = csv_multi.crbegin();

  auto rit2 = rit1 + 2;
  auto crit2 = crit1 + 2;
  EXPECT_EQ(*rit2, 30);
  EXPECT_EQ(*crit2, 30);

  EXPECT_EQ(rit2 - rit1, 2);
  EXPECT_EQ(crit2 - crit1, 2);
  EXPECT_EQ(rit1 - rit2, -2);
  EXPECT_EQ(crit1 - crit2, -2);

  auto rit3 = rit2 - 1;
  auto crit3 = crit2 - 1;
  EXPECT_EQ(*rit3, 40);
  EXPECT_EQ(*crit3, 40);

  rit1 += 3;
  EXPECT_EQ(*rit1, 20);

  rit1 -= 2;
  EXPECT_EQ(*rit1, 40);

  EXPECT_LT(sv_multi_.rbegin(), sv_multi_.rend());
  EXPECT_LE(sv_multi_.rbegin(), sv_multi_.rend());
  EXPECT_LE(sv_multi_.rbegin(), sv_multi_.rbegin());
  EXPECT_GT(sv_multi_.rend(), sv_multi_.rbegin());
  EXPECT_GE(sv_multi_.rend(), sv_multi_.rbegin());
  EXPECT_GE(sv_multi_.rbegin(), sv_multi_.rbegin());

  EXPECT_EQ(sv_multi_.rbegin() + sv_multi_.size(), sv_multi_.rend());
  EXPECT_EQ(csv_multi.crbegin() + csv_multi.size(), csv_multi.crend());
}

TEST_F(StaticVectorReverseIteratorTest, FullVector) {
  StaticVector<int, 5, uint8_t> sv_full;
  sv_full.push_back(1);
  sv_full.push_back(2);
  sv_full.push_back(3);
  sv_full.push_back(4);
  sv_full.push_back(5);

  ASSERT_TRUE(sv_full.full());
  ASSERT_EQ(sv_full.size(), 5);

  std::vector<int> expected_reverse = {5, 4, 3, 2, 1};
  std::vector<int> actual_reverse;

  std::copy(sv_full.rbegin(), sv_full.rend(), std::back_inserter(actual_reverse));
  EXPECT_EQ(actual_reverse, expected_reverse);
  EXPECT_THAT(actual_reverse, ::testing::ElementsAre(5, 4, 3, 2, 1));

  const auto& csv_full = sv_full;
  actual_reverse.clear();
  std::copy(csv_full.crbegin(), csv_full.crend(), std::back_inserter(actual_reverse));
  EXPECT_EQ(actual_reverse, expected_reverse);
  EXPECT_THAT(actual_reverse, ::testing::ElementsAre(5, 4, 3, 2, 1));
}