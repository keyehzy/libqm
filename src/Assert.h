#pragma once

#include <cstdlib>
#include <iostream>

namespace libqm {

inline void libqm_assert_fail(const char* expression, const char* file, int line,
                              const char* function) {
  std::cerr << file << ":" << line << ": " << function << ": Assertion `" << expression
            << "' failed." << std::endl;
  __builtin_trap();
}

#define LIBQM_ASSERT(condition) \
  ((condition) ? ((void)0) : ::libqm::libqm_assert_fail(#condition, __FILE__, __LINE__, __func__))
#define LIBQM_UNREACHABLE() LIBQM_ASSERT(false)

}  // namespace libqm
