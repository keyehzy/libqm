#!/usr/bin/env bash

set -e
set -x
set -o pipefail

CXXFLAGS="-Wall -Wextra -O2 -march=native"

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="$CXXFLAGS"
cmake --build build
./build/test_runner
