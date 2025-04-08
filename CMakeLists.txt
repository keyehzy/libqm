cmake_minimum_required(VERSION 3.14)
project(libqm VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(GTest REQUIRED)

add_executable(
  test_runner
  test/Operator-test.cpp
  src/Operator.cpp
  test/StaticVector-test.cpp

  src/Term.cpp
  test/Term-test.cpp
)

target_link_libraries(
  test_runner
  GTest::GTest
  GTest::Main gmock
)

add_test(NAME OperatorTests COMMAND test_runner)
enable_testing()
