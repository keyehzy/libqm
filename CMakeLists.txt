cmake_minimum_required(VERSION 3.14)
project(libqm VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)
# find_package(Armadillo REQUIRED)

add_executable(
  main
  src/diag.cpp
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

  src/NormalOrderer.h
  src/NormalOrderer.cpp

  src/MatrixElements.h
  src/MatrixElements.cpp

  src/Commutator.h
  src/Commutator.cpp

  src/xxhash.h
  src/xxhash.c

  src/Logger.h
  src/Logger.cpp
)

target_include_directories(
  qm
  PUBLIC
  OpenMP::OpenMP_CXX
  # ${ARMADILLO_INCLUDE_DIRS}
)

target_link_libraries(
  qm
  PUBLIC
  OpenMP::OpenMP_CXX
  # ${ARMADILLO_LIBRARIES}
  # -DARMA_DONT_USE_WRAPPER
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
  PRIVATE
  qm
  OpenMP::OpenMP_CXX
  GTest::GTest
  GTest::Main gmock
)

add_test(NAME OperatorTests COMMAND test_runner)
enable_testing()
