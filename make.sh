#!/usr/bin/env bash

set -e
set -x
set -o pipefail

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
./build/test_runner
