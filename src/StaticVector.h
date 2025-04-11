#pragma once

#include <algorithm>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>

#include "Assert.h"

template <typename T>
concept Trivial = std::is_trivially_copyable_v<T> && std::is_trivially_destructible_v<T>;

template <typename S>
concept UnsignedIntegralSizeType = std::unsigned_integral<S>;

template <size_t MaxCount, typename SizeType>
concept ValidStaticVectorCapacity = MaxCount > 0 &&
                                    MaxCount <= std::numeric_limits<SizeType>::max();

namespace detail {
inline constexpr void hash_combine(std::size_t& seed, std::size_t value) noexcept {
  constexpr std::size_t prime = 0x00000100000001b3ULL;
  seed ^= value;
  seed *= prime;
}
}  // namespace detail

/**
 * @brief A C++20 fixed-capacity container optimized for TRIVIALLY COPYABLE types. Similar to
 * boost::container::static_vector.
 */
template <Trivial T, size_t MaxCount, UnsignedIntegralSizeType SizeType = std::uint8_t>
  requires ValidStaticVectorCapacity<MaxCount, SizeType>
struct StaticVector {
  using value_type = T;
  using size_type = SizeType;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  static constexpr size_t static_capacity = MaxCount;

  constexpr StaticVector() noexcept : size_(0) {}
  ~StaticVector() noexcept = default;

  constexpr StaticVector(size_t count, const T& value) noexcept
      : size_(static_cast<size_type>(count)) {
    LIBQM_ASSERT(count <= static_capacity);
    std::fill_n(data_, count, value);
  }

  constexpr StaticVector(std::initializer_list<T> init) noexcept
      : size_(static_cast<size_type>(init.size())) {
    LIBQM_ASSERT(init.size() <= static_capacity);
    std::copy(init.begin(), init.end(), data_);
  }

  constexpr StaticVector(const StaticVector& other) noexcept : size_(other.size_) {
    std::copy_n(other.data_, size_, data_);
  }

  constexpr StaticVector(StaticVector&& other) noexcept : size_(other.size_) {
    std::copy_n(other.data_, size_, data_);
    other.size_ = 0;
  }
  constexpr StaticVector& operator=(const StaticVector& other) & noexcept {
    if (this != &other) {
      size_ = other.size_;
      std::copy_n(other.data_, size_, data_);
    }
    return *this;
  }

  constexpr StaticVector& operator=(StaticVector&& other) & noexcept {
    if (this != &other) {
      size_ = other.size_;
      std::copy_n(other.data_, size_, data_);
      other.size_ = 0;
    }
    return *this;
  }

  constexpr StaticVector& operator=(std::initializer_list<T> ilist) & noexcept {
    size_ = static_cast<size_type>(ilist.size());
    LIBQM_ASSERT(size_ <= static_capacity);
    std::copy(ilist.begin(), ilist.end(), data_);
    return *this;
  }

  constexpr reference at(size_t pos) noexcept {
    LIBQM_ASSERT(pos < size_);
    return data_[pos];
  }
  constexpr const_reference at(size_t pos) const noexcept {
    LIBQM_ASSERT(pos < size_);
    return data_[pos];
  }
  constexpr reference operator[](size_t pos) noexcept { return data_[pos]; }
  constexpr const_reference operator[](size_t pos) const noexcept { return data_[pos]; }
  constexpr reference front() noexcept {
    LIBQM_ASSERT(!empty());
    return data_[0];
  }
  constexpr const_reference front() const noexcept {
    LIBQM_ASSERT(!empty());
    return data_[0];
  }
  constexpr reference back() noexcept {
    LIBQM_ASSERT(!empty());
    return data_[size_ - 1];
  }
  constexpr const_reference back() const noexcept {
    LIBQM_ASSERT(!empty());
    return data_[size_ - 1];
  }
  constexpr T* data() noexcept { return data_; }
  constexpr const T* data() const noexcept { return data_; }

  constexpr std::span<T> get_span() noexcept { return std::span(data_, size_); }
  constexpr std::span<const T> get_span() const noexcept { return std::span(data_, size_); }

  constexpr iterator begin() noexcept { return data_; }
  constexpr const_iterator begin() const noexcept { return data_; }
  constexpr const_iterator cbegin() const noexcept { return begin(); }
  constexpr iterator end() noexcept { return data_ + size_; }
  constexpr const_iterator end() const noexcept { return data_ + size_; }
  constexpr const_iterator cend() const noexcept { return end(); }

  constexpr reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
  constexpr const_reverse_iterator rbegin() const noexcept { return const_reverse_iterator(end()); }
  constexpr const_reverse_iterator crbegin() const noexcept {
    return const_reverse_iterator(cend());
  }

  constexpr reverse_iterator rend() noexcept { return reverse_iterator(begin()); }
  constexpr const_reverse_iterator rend() const noexcept { return const_reverse_iterator(begin()); }
  constexpr const_reverse_iterator crend() const noexcept {
    return const_reverse_iterator(cbegin());
  }

  [[nodiscard]] constexpr bool empty() const noexcept { return size_ == 0; }
  [[nodiscard]] constexpr size_t size() const noexcept { return size_; }
  [[nodiscard]] constexpr size_t max_size() const noexcept { return static_capacity; }
  [[nodiscard]] constexpr size_t capacity() const noexcept { return static_capacity; }
  [[nodiscard]] constexpr size_t remaining_capacity() const noexcept {
    return static_capacity - size_;
  }
  [[nodiscard]] constexpr bool full() const noexcept { return size_ == static_capacity; }

  constexpr void clear() noexcept { size_ = 0; }

  constexpr void push_back(const T& value) noexcept {
    LIBQM_ASSERT(!full());
    data_[size_] = value;
    ++size_;
  }

  constexpr void push_back(T&& value) noexcept {
    LIBQM_ASSERT(!full());
    data_[size_] = std::move(value);
    ++size_;
  }

  template <typename... Args>
    requires std::constructible_from<T, Args...>
  constexpr reference emplace_back(Args&&... args) noexcept {
    LIBQM_ASSERT(!full());
    data_[size_] = T(std::forward<Args>(args)...);
    ++size_;
    return data_[size_ - 1];
  }

  constexpr void pop_back() noexcept {
    LIBQM_ASSERT(!empty());
    --size_;
  }

  constexpr void resize(size_t count) noexcept {
    LIBQM_ASSERT(count <= static_capacity);
    size_ = static_cast<size_type>(count);
  }

  constexpr void resize(size_t count, const value_type& value) noexcept {
    LIBQM_ASSERT(count <= static_capacity);
    if (count > size_) {
      std::fill_n(data_ + size_, count - size_, value);
    }
    size_ = static_cast<size_type>(count);
  }

  template <std::input_iterator InputIt, std::sentinel_for<InputIt> Sentinel>
    requires std::assignable_from<T&, std::iter_reference_t<InputIt>>
  constexpr void append_range(InputIt first, Sentinel last) noexcept {
    pointer current_dest = data_ + size_;
    size_type elements_available = static_cast<size_type>(static_capacity - size_);
    size_type elements_added = 0;

    for (; first != last && elements_added < elements_available;
         ++first, ++current_dest, ++elements_added) {
      *current_dest = *first;
    }
    size_ += elements_added;
    LIBQM_ASSERT(first == last);
  }

  template <size_t OtherMaxCount, typename OtherSizeType>
  constexpr void append_range(const StaticVector<T, OtherMaxCount, OtherSizeType>& other) noexcept {
    LIBQM_ASSERT(remaining_capacity() >= other.size());
    std::copy_n(other.begin(), other.size(), data_ + size_);
    size_ += other.size();
  }

  template <size_t OtherMaxCount, typename OtherSizeType>
  constexpr void append_range_reverse(
      const StaticVector<T, OtherMaxCount, OtherSizeType>& other) noexcept {
    LIBQM_ASSERT(remaining_capacity() >= other.size());
    std::copy_n(other.rbegin(), other.size(), data_ + size_);
    size_ += other.size();
  }

  template <size_t OtherMaxCount, typename OtherSizeType, typename Callback>
  constexpr void append_range(const StaticVector<T, OtherMaxCount, OtherSizeType>& other,
                              Callback&& callback) noexcept {
    LIBQM_ASSERT(remaining_capacity() >= other.size());
    for (const auto& item : other) {
      push_back(callback(item));
    }
  }

  template <size_t OtherMaxCount, typename OtherSizeType, typename Callback>
  constexpr void append_range_reverse(const StaticVector<T, OtherMaxCount, OtherSizeType>& other,
                                      Callback&& callback) noexcept {
    LIBQM_ASSERT(remaining_capacity() >= other.size());
    for (auto it = other.rbegin(); it != other.rend(); ++it) {
      push_back(callback(*it));
    }
  }

  constexpr void swap(StaticVector& other) noexcept {
    if (this != &other) {
      std::swap(size_, other.size_);
      for (size_t i = 0; i < static_capacity; ++i) {
        std::swap(data_[i], other.data_[i]);
      }
    }
  }

  [[nodiscard]] constexpr auto operator<=>(const StaticVector& other) const noexcept
    requires requires(const T& a, const T& b) {
      { a <=> b } -> std::convertible_to<std::strong_ordering>;
    }
  {
    return std::lexicographical_compare_three_way(begin(), end(), other.begin(), other.end(),
                                                  std::compare_three_way{});
  }

  [[nodiscard]] constexpr bool operator==(const StaticVector& other) const noexcept
    requires std::equality_comparable<T>
  {
    if (size() != other.size()) {
      return false;
    }
    return std::equal(begin(), end(), other.begin());
  }

 private:
  alignas(T) T data_[static_capacity]{};
  size_type size_;
};

template <Trivial T, size_t MaxCount, UnsignedIntegralSizeType SizeType>
  requires ValidStaticVectorCapacity<MaxCount, SizeType>
constexpr void swap(StaticVector<T, MaxCount, SizeType>& lhs,
                    StaticVector<T, MaxCount, SizeType>& rhs) noexcept {
  lhs.swap(rhs);
}

template <Trivial T, size_t MaxCount, UnsignedIntegralSizeType SizeType>
  requires ValidStaticVectorCapacity<MaxCount, SizeType>
struct std::hash<StaticVector<T, MaxCount, SizeType>> {
  [[nodiscard]] constexpr std::size_t operator()(
      const StaticVector<T, MaxCount, SizeType>& container) const noexcept {
    std::size_t seed = 0xcbf29ce484222325ULL;
    for (auto it = container.begin(); it != container.end(); ++it) {
      detail::hash_combine(seed, std::hash<T>{}(*it));
    }
    return seed;
  }
};

template <Trivial T, size_t MaxCount, UnsignedIntegralSizeType SizeType>
  requires ValidStaticVectorCapacity<MaxCount, SizeType>
constexpr bool has_consecutive_elements(const StaticVector<T, MaxCount, SizeType>& container) {
  return std::adjacent_find(container.begin(), container.end()) != container.end();
}

template <Trivial T, size_t MaxCount, UnsignedIntegralSizeType SizeType>
  requires ValidStaticVectorCapacity<MaxCount, SizeType>
constexpr StaticVector<T, MaxCount, SizeType> merge(
    const StaticVector<T, MaxCount, SizeType>& container1,
    const StaticVector<T, MaxCount, SizeType>& container2) {
  LIBQM_ASSERT(container1.remaining_capacity() >= container2.size());
  StaticVector<T, MaxCount, SizeType> result(container1);
  result.append_range(container2.begin(), container2.end());
  return result;
}