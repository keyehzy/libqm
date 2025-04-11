#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "Assert.h"

namespace libqm {
// Fermion operator:
// High Low
// XXXXXXYZ
// XXXXXXYZ
// Z: creation/annihilation bits
// Y: spin bits
// X: orbital bits

struct Operator {
  using ubyte = uint8_t;

  static constexpr size_t max_index() { return 1 << (8 * sizeof(ubyte) - 2); }
  static constexpr size_t max_unique_keys() { return 1 << (8 * sizeof(ubyte) - 1); }

  static constexpr ubyte kFermionTypeTagMask = 0x1;
  static constexpr ubyte kFermionSpinTagMask = 0x2;
  static constexpr ubyte kFermionBitsMask = 0xfc;
  static constexpr ubyte kFermionBitsShift = 2;

  enum class Type : ubyte { Creation = 0, Annihilation = 1 };
  enum class Spin : ubyte { Up = 0, Down = 1 };

  constexpr Operator() noexcept = default;

  constexpr explicit Operator(Type type, Spin spin, size_t value) noexcept
      : data(static_cast<ubyte>(value << kFermionBitsShift) | (static_cast<ubyte>(spin) << 1) |
             (static_cast<ubyte>(type) << 0)) {
    LIBQM_ASSERT(value < max_index());
  }

  constexpr explicit Operator(ubyte x) noexcept : data(x) {}

  constexpr Type type() const noexcept {
    return static_cast<Type>((data & kFermionTypeTagMask) >> 0);
  }

  constexpr Spin spin() const noexcept {
    return static_cast<Spin>((data & kFermionSpinTagMask) >> 1);
  }

  constexpr size_t value() const noexcept { return static_cast<size_t>(data >> kFermionBitsShift); }

  constexpr size_t key() const noexcept { return data & ~kFermionTypeTagMask; }

  constexpr Operator adjoint() const noexcept { return Operator(data ^ kFermionTypeTagMask); }

  std::string to_string() const;

  constexpr bool operator<(Operator other) const noexcept {
    if (type() != other.type()) {
      return type() < other.type();
    } else if (spin() != other.spin()) {
      return spin() < other.spin();
    } else {
      return value() < other.value();
    }
  }

  constexpr bool operator==(Operator other) const noexcept { return data == other.data; }

  constexpr bool commutes(Operator other) const noexcept {
    return (data ^ other.data) != kFermionTypeTagMask;
  }

  constexpr static Operator creation(Spin spin, size_t value) noexcept {
    return Operator(Type::Creation, spin, value);
  }

  constexpr static Operator annihilation(Spin spin, size_t value) noexcept {
    return Operator(Type::Annihilation, spin, value);
  }

  ubyte data{};
};
static_assert(sizeof(Operator) == 1);

}  // namespace libqm

template <>
struct std::hash<libqm::Operator> {
  [[nodiscard]] constexpr std::size_t operator()(libqm::Operator op) const noexcept {
    return op.data;
  }
};
