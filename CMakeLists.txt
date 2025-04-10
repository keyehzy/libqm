cmake_minimum_required(VERSION 3.14)
project(libqm VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(GTest REQUIRED)

add_executable(
  main
  src/main.cpp
)

target_link_libraries(
  main
  qm
)

add_library(
  qm
  src/Operator.h
  src/Operator.cpp

  src/StaticVector.h
  src/StaticVector.cpp

  src/Term.h
  src/Term.cpp

  src/Expression.h
  src/Expression.cpp

  src/IndexedHashSet.h
  src/IndexedHashSet.cpp

  src/Basis.h
  src/Basis.cpp
)

add_executable(
  test_runner
  test/Operator-test.cpp
  test/StaticVector-test.cpp
  test/Term-test.cpp
  test/Expression-test.cpp
  test/IndexedHashSet-test.cpp
)

target_link_libraries(
  test_runner
  qm
  GTest::GTest
  GTest::Main gmock
)

add_test(NAME OperatorTests COMMAND test_runner)
enable_testing()
