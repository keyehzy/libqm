#include "../src/Term.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <complex>
#include <limits>
#include <vector>

using ::testing::ContainerEq;

namespace libqm::testing {
const Operator cU0 = Operator::creation(Operator::Spin::Up, 0);
const Operator aU0 = Operator::annihilation(Operator::Spin::Up, 0);
const Operator cD1 = Operator::creation(Operator::Spin::Down, 1);
const Operator aD1 = Operator::annihilation(Operator::Spin::Down, 1);
const Operator cU2 = Operator::creation(Operator::Spin::Up, 2);
const Operator aU2 = Operator::annihilation(Operator::Spin::Up, 2);
const Operator aD2 = Operator::annihilation(Operator::Spin::Down, 2);

const Term::complex_type complex_one{1.0f, 0.0f};
const Term::complex_type complex_zero{0.0f, 0.0f};
const Term::complex_type complex_i{0.0f, 1.0f};
const Term::complex_type complex_custom{2.5f, -1.5f};

class TermTest : public ::testing::Test {
 protected:
  Term::container_type empty_ops;
  Term::container_type single_op_vec{cU0};
  Term::container_type multi_op_vec{cU0, aD1, cU2};
};

TEST_F(TermTest, DefaultConstructor) {
  Term t;
  EXPECT_EQ(t.c, complex_one);
  EXPECT_TRUE(t.operators.empty());
}

TEST_F(TermTest, OperatorConstructor) {
  Term t(cU0);
  EXPECT_EQ(t.c, complex_one);
  EXPECT_EQ(t.operators, Term::container_type{cU0});
}

TEST_F(TermTest, ComplexConstructor) {
  Term t(complex_custom);
  EXPECT_EQ(t.c, complex_custom);
  EXPECT_TRUE(t.operators.empty());
}

TEST_F(TermTest, ContainerConstructorCopy) {
  Term t(multi_op_vec);
  EXPECT_EQ(t.c, complex_one);
  EXPECT_EQ(t.operators, multi_op_vec);
}

TEST_F(TermTest, ContainerConstructorMove) {
  Term::container_type ops = multi_op_vec;
  Term t(std::move(ops));
  EXPECT_EQ(t.c, complex_one);
  EXPECT_EQ(t.operators, multi_op_vec);
  EXPECT_TRUE(ops.empty());
}

TEST_F(TermTest, ComplexContainerConstructorCopy) {
  Term t(complex_custom, multi_op_vec);
  EXPECT_EQ(t.c, complex_custom);
  EXPECT_EQ(t.operators, multi_op_vec);
}

TEST_F(TermTest, ComplexContainerConstructorMove) {
  Term::container_type ops = multi_op_vec;
  Term t(complex_custom, std::move(ops));
  EXPECT_EQ(t.c, complex_custom);
  EXPECT_EQ(t.operators, multi_op_vec);
  EXPECT_TRUE(ops.empty());
}

TEST_F(TermTest, InitializerListConstructor) {
  Term t({cU0, aD1});
  EXPECT_EQ(t.c, complex_one);
  EXPECT_EQ(t.operators, (Term::container_type{cU0, aD1}));
}

TEST_F(TermTest, ComplexInitializerListConstructor) {
  Term t(complex_i, {cU0, aD1});
  EXPECT_EQ(t.c, complex_i);
  EXPECT_EQ(t.operators, (Term::container_type{cU0, aD1}));
}

TEST_F(TermTest, InitializerListConstructorEmpty) {
  Term t({});
  EXPECT_EQ(t.c, complex_one);
  EXPECT_TRUE(t.operators.empty());
}

TEST_F(TermTest, ZeroCoefficient) {
  Term t1(complex_custom, {cU0, aD1});
  Term::container_type t1_orig_ops(t1.operators);

  // Multiply by zero
  t1 *= complex_zero;
  EXPECT_EQ(t1, Term(complex_zero, t1_orig_ops));

  // Create a term that is already zero
  Term t2(complex_zero, {aD2, cU2});
  EXPECT_EQ(t2, Term(complex_zero, {aD2, cU2}));

  // Divide zero by non-zero
  t2 /= complex_i;
  EXPECT_EQ(t2, Term(complex_zero, {aD2, cU2}));

  // Multiply non-zero term by zero term
  Term t3(complex_i, {cU0});
  Term t4(complex_zero, {aD1});
  Term t5 = t3 * t4;
  EXPECT_EQ(t5, Term(complex_zero, {cU0, aD1}));

  // Multiply zero term by non-zero term
  Term t6 = t4 * t3;
  EXPECT_EQ(t6, Term(complex_zero, {aD1, cU0}));
}

TEST_F(TermTest, ComplexInitializerListConstructorEmpty) {
  Term t(complex_custom, {});
  EXPECT_EQ(t.c, complex_custom);
  EXPECT_TRUE(t.operators.empty());
}

TEST_F(TermTest, CopyConstructor) {
  Term t1(complex_custom, {cU0, aD1});
  Term t2(t1);

  EXPECT_EQ(t2, t1);

  t2 *= complex_i;
  t2 *= cU2;
  EXPECT_EQ(t2, Term(complex_custom * complex_i, {cU0, aD1, cU2}));
}

TEST_F(TermTest, CopyAssignment) {
  Term t1(complex_custom, {cU0, aD1});
  Term t2;
  Term t3(complex_i, {aD2});

  t2 = t1;
  EXPECT_EQ(t2, t1);

  t2 *= complex_i;
  EXPECT_EQ(t2, Term(complex_custom * complex_i, {cU0, aD1}));

  t1 = t1;
  EXPECT_EQ(t1, Term(complex_custom, {cU0, aD1}));
}

TEST_F(TermTest, MoveConstructor) {
  Term t1(complex_custom, {cU0, aD1});
  Term t1_copy = t1;

  Term t2(std::move(t1));

  EXPECT_EQ(t2, Term(complex_custom, {cU0, aD1}));
  EXPECT_TRUE(t1.operators.empty());
}

TEST_F(TermTest, MoveAssignment) {
  Term t1(complex_custom, {cU0, aD1});
  Term t1_copy = t1;
  Term t2(complex_i, {aD2});

  t2 = std::move(t1);

  EXPECT_EQ(t2, Term(complex_custom, {cU0, aD1}));

  // Check the state of the moved-from object t1
  EXPECT_TRUE(t1.operators.empty());
}

TEST_F(TermTest, Size) {
  Term t0;
  Term t1(cU0);
  Term t3(complex_custom, {cU0, aD1, cU2});
  EXPECT_EQ(t0.size(), 0);
  EXPECT_EQ(t1.size(), 1);
  EXPECT_EQ(t3.size(), 3);
}

TEST_F(TermTest, EqualityOperators) {
  Term t1(complex_custom, {cU0, aD1});
  Term t2(complex_custom, {cU0, aD1});
  Term t3(complex_i, {cU0, aD1});
  Term t4(complex_custom, {cU0, aD2});
  Term t5(complex_custom, {aD1, cU0});
  Term t6(complex_custom, {cU0});
  Term t7;
  Term t8;

  EXPECT_TRUE(t1 == t2);
  EXPECT_FALSE(t1 != t2);
  EXPECT_EQ(t1, t2);

  EXPECT_TRUE(t1 == t1);
  EXPECT_FALSE(t1 != t1);
  EXPECT_EQ(t1, t1);

  EXPECT_TRUE(t1 != t3);
  EXPECT_FALSE(t1 == t3);
  EXPECT_NE(t1, t3);

  EXPECT_TRUE(t1 != t4);
  EXPECT_FALSE(t1 == t4);
  EXPECT_NE(t1, t4);

  EXPECT_TRUE(t1 != t5);
  EXPECT_FALSE(t1 == t5);
  EXPECT_NE(t1, t5);

  EXPECT_TRUE(t1 != t6);
  EXPECT_FALSE(t1 == t6);
  EXPECT_NE(t1, t6);

  EXPECT_TRUE(t7 == t8);
  EXPECT_FALSE(t7 != t8);
  EXPECT_EQ(t7, t8);

  EXPECT_TRUE(t1 != t7);
  EXPECT_FALSE(t1 == t7);
  EXPECT_NE(t1, t7);
}

TEST_F(TermTest, Adjoint) {
  Term t1(complex_custom, {cU0, aD1, cU2});
  Term adj_t1 = t1.adjoint();

  Term::complex_type expected_c = std::conj(complex_custom);
  Term::container_type expected_ops = {aU2, cD1, aU0};

  EXPECT_EQ(adj_t1.c, expected_c);
  EXPECT_EQ(adj_t1.operators, expected_ops);
  EXPECT_EQ(adj_t1.size(), 3);

  Term t_empty(complex_i);
  Term adj_t_empty = t_empty.adjoint();
  EXPECT_EQ(adj_t_empty.c, std::conj(complex_i));
  EXPECT_TRUE(adj_t_empty.operators.empty());

  Term t_single(complex_one, {cD1});
  Term adj_t_single = t_single.adjoint();
  EXPECT_EQ(adj_t_single.c, complex_one);
  EXPECT_EQ(adj_t_single.operators, Term::container_type{aD1});
}

TEST_F(TermTest, MultiplyAssignComplex) {
  Term t(complex_custom, {cU0});
  t *= complex_i;
  Term::complex_type expected_c = complex_custom * complex_i;
  EXPECT_EQ(t.c, expected_c);
  EXPECT_EQ(t.operators, Term::container_type{cU0});
}

TEST_F(TermTest, MultiplyAssignOperator) {
  Term t(complex_custom, {cU0});
  t *= aD1;
  EXPECT_EQ(t.c, complex_custom);
  EXPECT_EQ(t.operators, (Term::container_type{cU0, aD1}));
}

TEST_F(TermTest, MultiplyAssignTerm) {
  Term t1(complex_custom, {cU0});
  Term t2(complex_i, {aD1, cU2});

  t1 *= t2;

  Term::complex_type expected_c = complex_custom * complex_i;
  Term::container_type expected_ops = {cU0, aD1, cU2};

  EXPECT_EQ(t1.c, expected_c);
  EXPECT_EQ(t1.operators, expected_ops);
}

TEST_F(TermTest, MultiplyAssignTermEmptyOps) {
  Term t1(complex_custom, {cU0});
  Term t2(complex_i, {});

  t1 *= t2;

  Term::complex_type expected_c = complex_custom * complex_i;
  Term::container_type expected_ops = {cU0};

  EXPECT_EQ(t1.c, expected_c);
  EXPECT_EQ(t1.operators, expected_ops);

  Term t3(complex_i, {});
  Term t4(complex_custom, {cU0});
  t3 *= t4;
  EXPECT_EQ(t3.c, expected_c);
  EXPECT_EQ(t3.operators, expected_ops);
}

TEST_F(TermTest, DivideAssignComplex) {
  Term t(complex_custom, {cU0});
  t /= complex_i;
  Term::complex_type expected_c = complex_custom / complex_i;
  EXPECT_EQ(t.c, expected_c);
  EXPECT_EQ(t.operators, Term::container_type{cU0});
}

TEST_F(TermTest, DivideAssignComplexByZero) {
  Term t(complex_custom, {cU0});
  Term::complex_type zero{0.0f, 0.0f};
  t /= zero;

  EXPECT_TRUE(std::isinf(t.c.real()) || std::isnan(t.c.real()));
  EXPECT_TRUE(std::isinf(t.c.imag()) || std::isnan(t.c.imag()));
  EXPECT_EQ(t.operators, Term::container_type{cU0});
}

TEST_F(TermTest, MultiplyTermTerm) {
  Term t1(complex_custom, {cU0});
  Term t2(complex_i, {aD1});
  Term t1_orig = t1;
  Term t2_orig = t2;

  Term result = t1 * t2;

  Term expected_result(complex_custom * complex_i, {cU0, aD1});
  EXPECT_EQ(result, expected_result);

  EXPECT_EQ(t1, Term(complex_custom, {cU0}));
  EXPECT_EQ(t2, Term(complex_i, {aD1}));
}

TEST_F(TermTest, MultiplyTermOperator) {
  Term t1(complex_custom, {cU0});
  Operator op = aD1;
  Term t1_orig = t1;

  Term result = t1 * op;

  Term expected_result(complex_custom, {cU0, aD1});
  EXPECT_EQ(result, expected_result);
  EXPECT_EQ(t1, Term(complex_custom, {cU0}));
}

TEST_F(TermTest, MultiplyTermComplex) {
  Term t1(complex_custom, {cU0});
  Term::complex_type val = complex_i;
  Term t1_orig = t1;

  Term result = t1 * val;

  Term expected_result(complex_custom * complex_i, {cU0});
  EXPECT_EQ(result, expected_result);
  EXPECT_EQ(t1, Term(complex_custom, {cU0}));
}

TEST_F(TermTest, DivideTermComplex) {
  Term t1(complex_custom, {cU0});
  Term::complex_type val = complex_i;
  Term t1_orig = t1;

  Term result = t1 / val;

  Term expected_result(complex_custom / complex_i, {cU0});
  EXPECT_EQ(result, expected_result);
  EXPECT_EQ(t1, t1_orig);
}

TEST_F(TermTest, MultiplyComplexTerm) {
  Term t1(complex_custom, {cU0});
  Term::complex_type val = complex_i;
  Term t1_orig = t1;

  Term result = val * t1;

  Term expected_result(complex_i * complex_custom, {cU0});
  EXPECT_EQ(result, expected_result);
  EXPECT_EQ(t1, Term(complex_custom, {cU0}));
}

TEST_F(TermTest, MultiplyOperatorTerm) {
  Term t1(complex_custom, {aD1});
  Operator op = cU0;
  Term t1_orig = t1;

  Term result = op * t1;

  Term expected_result(complex_custom, {cU0, aD1});
  EXPECT_EQ(result, expected_result);
  EXPECT_EQ(t1, Term(complex_custom, {aD1}));
}

TEST(TermHelperTest, IsDiagonal) {
  Term ops_empty({});
  Term ops_diag1({cU0, aU0});
  Term ops_diag2({cU0, cD1, aU0, aD1});
  Term ops_diag3({cD1, aD1});
  Term ops_non_diag1({cU0});
  Term ops_non_diag2({cU0, aD1});
  Term ops_non_diag3({cU0, cU0, aU0});

  EXPECT_TRUE(is_diagonal(ops_empty));
  EXPECT_TRUE(is_diagonal(ops_diag1));
  EXPECT_TRUE(is_diagonal(ops_diag2));
  EXPECT_TRUE(is_diagonal(ops_diag3));
  EXPECT_FALSE(is_diagonal(ops_non_diag1));
  EXPECT_FALSE(is_diagonal(ops_non_diag2));
  EXPECT_FALSE(is_diagonal(ops_non_diag3));
}

TEST(TermFactoryTest, Creation) {
  Term t = creation(Operator::Spin::Down, 5);
  Term expected_t(complex_one, {Operator::creation(Operator::Spin::Down, 5)});
  EXPECT_EQ(t, expected_t);
}

TEST(TermFactoryTest, Annihilation) {
  Term t = annihilation(Operator::Spin::Up, 3);
  Term expected_t(complex_one, {Operator::annihilation(Operator::Spin::Up, 3)});
  EXPECT_EQ(t, expected_t);
}

TEST(TermFactoryTest, OneBody) {
  Term t = one_body(Operator::Spin::Up, 1, Operator::Spin::Down, 2);
  Term expected_t(complex_one, {Operator::creation(Operator::Spin::Up, 1),
                                Operator::annihilation(Operator::Spin::Down, 2)});
  EXPECT_EQ(t, expected_t);
}

TEST(TermFactoryTest, Density) {
  Term t = density(Operator::Spin::Up, 4);
  Term expected_t(complex_one, {Operator::creation(Operator::Spin::Up, 4),
                                Operator::annihilation(Operator::Spin::Up, 4)});
  EXPECT_EQ(t, expected_t);
}

TEST(TermFactoryTest, DensityDensity) {
  Term t = density_density(Operator::Spin::Up, 1, Operator::Spin::Down, 2);
  Term expected_t(complex_one, {Operator::creation(Operator::Spin::Up, 1),
                                Operator::annihilation(Operator::Spin::Up, 1),
                                Operator::creation(Operator::Spin::Down, 2),
                                Operator::annihilation(Operator::Spin::Down, 2)});
  EXPECT_EQ(t, expected_t);
}
}  // namespace libqm::testing
