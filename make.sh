#!/usr/bin/env bash

set -e
set -x
set -o pipefail

CXX="g++"
CXXFLAGS="-O2 -Wall -Wformat -Wformat=2 -Wconversion -Wimplicit-fallthrough -Werror=format-security -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -D_GLIBCXX_ASSERTIONS -fstrict-flex-arrays=3 -fstack-clash-protection -fstack-protector-strong -Wl,-z,nodlopen -Wl,-z,noexecstack -Wl,-z,relro -Wl,-z,now -Wl,--as-needed -Wl,--no-copy-dt-needed-entries"

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="$CXXFLAGS"
cmake --build build
./build/test_runner
