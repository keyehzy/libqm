
#include "../src/Expression.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <initializer_list>
#include <limits>

#include "../src/Operator.h"
#include "../src/Term.h"

// Use the namespace of the code being tested
using namespace libqm;

// Helper function to compare complex numbers with tolerance
// GTest might have built-in ways, but this is explicit
::testing::AssertionResult AssertComplexNear(const std::complex<float>& expected,
                                             const std::complex<float>& actual,
                                             float tolerance = 1e-6f) {
  if (std::abs(expected.real() - actual.real()) <= tolerance &&
      std::abs(expected.imag() - actual.imag()) <= tolerance) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure()
           << "Value (" << actual.real() << ", " << actual.imag() << ")"
           << " is not close to expected (" << expected.real() << ", " << expected.imag() << ")"
           << ", tolerance=" << tolerance;
  }
}

// Helper to check if an expression contains a specific term (operators + coefficient)
::testing::AssertionResult AssertTermExists(const Expression& expr, const Term::container_type& ops,
                                            const Expression::complex_type& expected_coeff,
                                            float tolerance = 1e-6f) {
  auto it = expr.hashmap.find(ops);
  if (it == expr.hashmap.end()) {
    return ::testing::AssertionFailure() << "Term with specified operators not found.";
  }
  return AssertComplexNear(expected_coeff, it->second, tolerance);
}

// Helper to check if an expression *does not* contain a specific term
::testing::AssertionResult AssertTermNotExists(const Expression& expr,
                                               const Term::container_type& ops) {
  if (expr.hashmap.count(ops)) {
    // Check if coefficient is non-zero before failing
    if (std::abs(expr.hashmap.at(ops).real()) > std::numeric_limits<float>::epsilon() ||
        std::abs(expr.hashmap.at(ops).imag()) > std::numeric_limits<float>::epsilon()) {
      return ::testing::AssertionFailure()
             << "Term with specified operators exists unexpectedly with coefficient "
             << expr.hashmap.at(ops) << ".";
    }
  }
  // Either not present or coefficient is effectively zero
  return ::testing::AssertionSuccess();
}

// Define some common operators and terms for convenience
const Operator op_c0u = Operator::creation(Operator::Spin::Up, 0);
const Operator op_a0u = Operator::annihilation(Operator::Spin::Up, 0);
const Operator op_c1d = Operator::creation(Operator::Spin::Down, 1);
const Operator op_a1d = Operator::annihilation(Operator::Spin::Down, 1);
const Operator op_c2u = Operator::creation(Operator::Spin::Up, 2);
const Operator op_a2u = Operator::annihilation(Operator::Spin::Up, 2);

const Term term_n0u({op_c0u, op_a0u});      // Number operator n_{0,up}
const Term term_n1d({op_c1d, op_a1d});      // Number operator n_{1,down}
const Term term_hop_02u({op_c0u, op_a2u});  // Hopping 2u -> 0u
const Term term_hop_20u({op_c2u, op_a0u});  // Hopping 0u -> 2u

const Expression::complex_type complex_one(1.0f, 0.0f);
const Expression::complex_type complex_zero(0.0f, 0.0f);
const Expression::complex_type complex_i(0.0f, 1.0f);
const Expression::complex_type complex_two(2.0f, 0.0f);
const Expression::complex_type complex_half(0.5f, 0.0f);
const Expression::complex_type complex_c1(2.0f, 3.0f);
const Expression::complex_type complex_c2(-1.0f, 1.5f);

// Test Fixture for common setup if needed (optional here, but good practice)
class ExpressionTest : public ::testing::Test {
 protected:
  // Per-test setup can go here if needed
  void SetUp() override {}

  // Per-test teardown can go here if needed
  void TearDown() override {}

  // Helper to create a container_type easily
  Term::container_type make_ops(std::initializer_list<Operator> ops) {
    return Term::container_type(ops);
  }
};

// --- Test Cases ---

TEST_F(ExpressionTest, DefaultConstructor) {
  Expression expr;
  EXPECT_EQ(expr.size(), 0);
  EXPECT_TRUE(expr.hashmap.empty());
}

TEST_F(ExpressionTest, ConstructFromComplex) {
  Expression expr(complex_c1);
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({}), complex_c1));
}

TEST_F(ExpressionTest, ConstructFromOperator) {
  Expression expr(op_c0u);
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_one));
}

TEST_F(ExpressionTest, ConstructFromTermCopy) {
  Term t1(complex_c1, {op_c0u, op_a1d});
  Expression expr(t1);
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u, op_a1d}), complex_c1));
  // Ensure original term is unchanged
  EXPECT_EQ(t1.c, complex_c1);
  EXPECT_EQ(t1.operators, make_ops({op_c0u, op_a1d}));
}

TEST_F(ExpressionTest, ConstructFromTermMove) {
  Term t1(complex_c1, {op_c0u, op_a1d});
  Term::container_type original_ops = t1.operators;  // Copy ops for checking later
  Expression expr(std::move(t1));

  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, original_ops, complex_c1));
  // Check moved-from state (operators should be empty or valid unspecified)
  // According to StaticVector move ctor, it should be empty.
  EXPECT_TRUE(t1.operators.empty());
  // Coefficient is trivially copyable, so it remains
  EXPECT_EQ(t1.c, complex_c1);
}

TEST_F(ExpressionTest, ConstructFromContainerCopy) {
  auto ops = make_ops({op_c0u, op_a1d});
  Expression expr(ops);
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, ops, complex_one));
  // Ensure original container is unchanged (it's a copy)
  ASSERT_EQ(ops.size(), 2);
  ASSERT_EQ(ops[0], op_c0u);
  ASSERT_EQ(ops[1], op_a1d);
}

TEST_F(ExpressionTest, ConstructFromContainerMove) {
  auto ops_orig = make_ops({op_c0u, op_a1d});
  auto ops_copy = ops_orig;  // Keep a copy for checking
  Expression expr(std::move(ops_orig));

  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, ops_copy, complex_one));
  // Check moved-from state (should be empty for StaticVector)
  EXPECT_TRUE(ops_orig.empty());
}

TEST_F(ExpressionTest, ConstructFromComplexAndContainerCopy) {
  auto ops = make_ops({op_c0u, op_a1d});
  Expression expr(complex_c1, ops);
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, ops, complex_c1));
  // Ensure original container is unchanged
  ASSERT_EQ(ops.size(), 2);
}

TEST_F(ExpressionTest, ConstructFromComplexAndContainerMove) {
  auto ops_orig = make_ops({op_c0u, op_a1d});
  auto ops_copy = ops_orig;
  Expression expr(complex_c1, std::move(ops_orig));

  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, ops_copy, complex_c1));
  EXPECT_TRUE(ops_orig.empty());
}

TEST_F(ExpressionTest, ConstructFromInitializerList) {
  Expression expr({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  EXPECT_EQ(expr.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2));
}

TEST_F(ExpressionTest, ConstructFromInitializerListCombineTerms) {
  Expression expr({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d}), Term(complex_one, {op_c0u})});
  EXPECT_EQ(expr.size(), 2);  // op_c0u terms should combine
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1 + complex_one));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2));
}

TEST_F(ExpressionTest, CopyConstructor) {
  Expression expr1({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Expression expr2(expr1);

  EXPECT_EQ(expr1.size(), 2);
  EXPECT_EQ(expr2.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_c0u}), complex_c1));
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_a1d}), complex_c2));

  // Modify copy, check original is unchanged
  expr2 += Term(complex_one, {op_c1d});
  EXPECT_EQ(expr1.size(), 2);
  EXPECT_EQ(expr2.size(), 3);
  ASSERT_TRUE(AssertTermNotExists(expr1, make_ops({op_c1d})));
}

TEST_F(ExpressionTest, MoveConstructor) {
  Expression expr1({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Expression expr2(std::move(expr1));

  EXPECT_EQ(expr2.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_c0u}), complex_c1));
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_a1d}), complex_c2));

  // Check moved-from state (hashmap should be empty or valid unspecified)
  // Standard map move leaves source empty.
  EXPECT_TRUE(expr1.hashmap.empty());
  EXPECT_EQ(expr1.size(), 0);
}

TEST_F(ExpressionTest, CopyAssignment) {
  Expression expr1({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Expression expr2({Term(complex_one, {op_c1d})});  // Different initial content
  expr2 = expr1;

  EXPECT_EQ(expr1.size(), 2);
  EXPECT_EQ(expr2.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_c0u}), complex_c1));
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_a1d}), complex_c2));
  ASSERT_TRUE(AssertTermNotExists(expr2, make_ops({op_c1d})));  // Original content gone

  // Modify copy, check original is unchanged
  expr2 += Term(complex_one, {op_c1d});
  EXPECT_EQ(expr1.size(), 2);
  EXPECT_EQ(expr2.size(), 3);
  ASSERT_TRUE(AssertTermNotExists(expr1, make_ops({op_c1d})));
}

TEST_F(ExpressionTest, MoveAssignment) {
  Expression expr1({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Expression expr2({Term(complex_one, {op_c1d})});  // Different initial content
  expr2 = std::move(expr1);

  EXPECT_EQ(expr2.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_c0u}), complex_c1));
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_a1d}), complex_c2));
  ASSERT_TRUE(AssertTermNotExists(expr2, make_ops({op_c1d})));  // Original content gone

  // Check moved-from state
  EXPECT_TRUE(expr1.hashmap.empty());
  EXPECT_EQ(expr1.size(), 0);
}

TEST_F(ExpressionTest, AdjointSimpleTerm) {
  Expression expr(Term(complex_c1, {op_c0u, op_a1d}));
  Expression adj = expr.adjoint();

  EXPECT_EQ(adj.size(), 1);
  // Adjoint: conj(c) * op_a1d.adjoint() * op_c0u.adjoint()
  // = conj(c) * op_c1d * op_a0u
  ASSERT_TRUE(AssertTermExists(adj, make_ops({op_c1d, op_a0u}), std::conj(complex_c1)));
}

TEST_F(ExpressionTest, AdjointIdentityTerm) {
  Expression expr(complex_c1);  // Identity term
  Expression adj = expr.adjoint();

  EXPECT_EQ(adj.size(), 1);
  ASSERT_TRUE(AssertTermExists(adj, make_ops({}), std::conj(complex_c1)));
}

TEST_F(ExpressionTest, AdjointMultipleTerms) {
  Expression expr({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_c1d, op_a2u})});
  Expression adj = expr.adjoint();

  EXPECT_EQ(adj.size(), 2);
  // Adjoint Term 1: conj(c1) * op_a0u
  ASSERT_TRUE(AssertTermExists(adj, make_ops({op_a0u}), std::conj(complex_c1)));
  // Adjoint Term 2: conj(c2) * op_a2u.adjoint() * op_c1d.adjoint()
  // = conj(c2) * op_c2u * op_a1d
  ASSERT_TRUE(AssertTermExists(adj, make_ops({op_c2u, op_a1d}), std::conj(complex_c2)));
}

TEST_F(ExpressionTest, AdjointIsSelfInverseForHermitian) {
  // Create a Hermitian expression, e.g., n_0u = c_0u^+ a_0u
  Expression expr(term_n0u);
  Expression adj = expr.adjoint();

  // n_0u.adjoint() = (c_0u^+ a_0u)^+ = a_0u^+ (c_0u^+)^+ = c_0u^+ a_0u = n_0u
  EXPECT_EQ(adj.size(), 1);
  ASSERT_TRUE(AssertTermExists(adj, make_ops({op_c0u, op_a0u}), complex_one));
  EXPECT_EQ(expr.hashmap, adj.hashmap);  // Check for equality

  // Check (adj).adjoint() == expr
  Expression adj_adj = adj.adjoint();
  EXPECT_EQ(adj_adj.size(), 1);
  ASSERT_TRUE(AssertTermExists(adj_adj, make_ops({op_c0u, op_a0u}), complex_one));
  EXPECT_EQ(expr.hashmap, adj_adj.hashmap);
}

TEST_F(ExpressionTest, OperatorPlusEqualsComplex) {
  Expression expr;
  expr += complex_c1;
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({}), complex_c1));

  expr += complex_c2;
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({}), complex_c1 + complex_c2));

  // Add to existing non-identity expression
  Expression expr2(op_c0u);
  expr2 += complex_c1;
  EXPECT_EQ(expr2.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_c0u}), complex_one));
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({}), complex_c1));
}

TEST_F(ExpressionTest, OperatorMinusEqualsComplex) {
  Expression expr;
  expr -= complex_c1;
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({}), -complex_c1));

  expr -= complex_c2;
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({}), -complex_c1 - complex_c2));

  // Subtract from existing non-identity expression
  Expression expr2(op_c0u);
  expr2 -= complex_c1;
  EXPECT_EQ(expr2.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_c0u}), complex_one));
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({}), -complex_c1));
}

TEST_F(ExpressionTest, OperatorMultiplyEqualsComplex) {
  Expression expr({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  expr *= complex_two;
  EXPECT_EQ(expr.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1 * complex_two));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2 * complex_two));
}

TEST_F(ExpressionTest, OperatorDivideEqualsComplex) {
  Expression expr({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  expr /= complex_two;
  EXPECT_EQ(expr.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1 / complex_two));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2 / complex_two));
}

TEST_F(ExpressionTest, OperatorPlusEqualsTerm) {
  Expression expr;
  Term t1(complex_c1, {op_c0u});
  expr += t1;
  EXPECT_EQ(expr.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1));

  Term t2(complex_c2, {op_a1d});
  expr += t2;
  EXPECT_EQ(expr.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2));

  // Add term with same operators
  Term t3(complex_one, {op_c0u});
  expr += t3;
  EXPECT_EQ(expr.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1 + complex_one));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2));
}

TEST_F(ExpressionTest, OperatorMinusEqualsTerm) {
  Expression expr({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Term t1(complex_one, {op_c0u});
  expr -= t1;
  EXPECT_EQ(expr.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1 - complex_one));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2));

  Term t2(complex_two, {op_c1d});  // New term
  expr -= t2;
  EXPECT_EQ(expr.size(), 3);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u}), complex_c1 - complex_one));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d}), complex_c2));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c1d}), -complex_two));
}

TEST_F(ExpressionTest, OperatorPlusEqualsExpression) {
  Expression expr1({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Expression expr2({Term(complex_one, {op_c0u}), Term(complex_two, {op_c1d})});
  expr1 += expr2;

  EXPECT_EQ(expr1.size(), 3);
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_c0u}), complex_c1 + complex_one));  // Combined
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_a1d}), complex_c2));                // From expr1
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_c1d}), complex_two));               // From expr2
}

TEST_F(ExpressionTest, OperatorMinusEqualsExpression) {
  Expression expr1({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Expression expr2({Term(complex_one, {op_c0u}), Term(complex_two, {op_c1d})});
  expr1 -= expr2;

  EXPECT_EQ(expr1.size(), 3);
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_c0u}), complex_c1 - complex_one));  // Combined
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_a1d}), complex_c2));                // From expr1
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_c1d}), -complex_two));  // From expr2 (negated)
}

TEST_F(ExpressionTest, OperatorMultiplyEqualsTerm) {
  Expression expr({Term(complex_c1, {op_c0u}), Term(complex_c2, {op_a1d})});
  Term term_mul(complex_two, {op_c1d});

  expr *= term_mul;  // expr = (c1*c0u + c2*a1d) * (2*c1d)
                     //      = 2*c1*(c0u c1d) + 2*c2*(a1d c1d)

  EXPECT_EQ(expr.size(), 2);
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_c0u, op_c1d}), complex_c1 * complex_two));
  ASSERT_TRUE(AssertTermExists(expr, make_ops({op_a1d, op_c1d}), complex_c2 * complex_two));
}

TEST_F(ExpressionTest, OperatorMultiplyEqualsExpression) {
  Expression expr1({Term(complex_c1, {op_c0u}), Term(complex_one, {})});        // c1*c0u + 1
  Expression expr2({Term(complex_c2, {op_a1d}), Term(complex_two, {op_c1d})});  // c2*a1d + 2*c1d

  // expr1 * expr2 = (c1*c0u + 1) * (c2*a1d + 2*c1d)
  // = c1*c2*(c0u a1d) + c1*2*(c0u c1d) + 1*c2*(a1d) + 1*2*(c1d)
  expr1 *= expr2;

  EXPECT_EQ(expr1.size(), 4);
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_c0u, op_a1d}), complex_c1 * complex_c2));
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_c0u, op_c1d}), complex_c1 * complex_two));
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_a1d}), complex_c2));
  ASSERT_TRUE(AssertTermExists(expr1, make_ops({op_c1d}), complex_two));
}

TEST_F(ExpressionTest, OperatorMultiplyEqualsExpressionEmpty) {
  Expression expr1({Term(complex_c1, {op_c0u})});
  Expression expr_empty;  // Default constructed

  expr1 *= expr_empty;
  // Multiplying by an empty expression (zero) should result in zero (empty expression)
  EXPECT_EQ(expr1.size(), 0);
  EXPECT_TRUE(expr1.hashmap.empty());

  Expression expr2({Term(complex_c1, {op_c0u})});
  Expression expr_ident(complex_one);  // Identity expression
  expr2 *= expr_ident;
  // Multiplying by identity should not change the expression
  EXPECT_EQ(expr2.size(), 1);
  ASSERT_TRUE(AssertTermExists(expr2, make_ops({op_c0u}), complex_c1));
}

TEST_F(ExpressionTest, NonMemberOperators) {
  Expression expr_a({Term(complex_c1, {op_c0u})});
  Expression expr_b({Term(complex_c2, {op_a1d})});
  Term term_t(complex_two, {op_c1d});
  Expression::complex_type complex_val = complex_i;

  // Expression + Expression
  Expression sum_expr = expr_a + expr_b;
  Expression sum_expr_manual = expr_a;
  sum_expr_manual += expr_b;
  EXPECT_EQ(sum_expr.hashmap, sum_expr_manual.hashmap);

  // Expression + Term
  Expression sum_term = expr_a + term_t;
  Expression sum_term_manual = expr_a;
  sum_term_manual += term_t;
  EXPECT_EQ(sum_term.hashmap, sum_term_manual.hashmap);

  // Term + Expression
  Expression sum_term_rev = term_t + expr_a;
  Expression sum_term_rev_manual(term_t);
  sum_term_rev_manual += expr_a;
  EXPECT_EQ(sum_term_rev.hashmap, sum_term_rev_manual.hashmap);

  // Expression + Complex
  Expression sum_complex = expr_a + complex_val;
  Expression sum_complex_manual = expr_a;
  sum_complex_manual += complex_val;
  EXPECT_EQ(sum_complex.hashmap, sum_complex_manual.hashmap);

  // Complex + Expression
  Expression sum_complex_rev = complex_val + expr_a;
  Expression sum_complex_rev_manual = expr_a;
  sum_complex_rev_manual += complex_val;  // Commutative
  EXPECT_EQ(sum_complex_rev.hashmap, sum_complex_rev_manual.hashmap);

  // Term + Term -> Expression
  Expression sum_term_term = term_n0u + term_n1d;
  EXPECT_EQ(sum_term_term.size(), 2);
  ASSERT_TRUE(AssertTermExists(sum_term_term, term_n0u.operators, term_n0u.c));
  ASSERT_TRUE(AssertTermExists(sum_term_term, term_n1d.operators, term_n1d.c));

  // Subtraction (similar checks)
  Expression diff_expr = expr_a - expr_b;
  Expression diff_expr_manual = expr_a;
  diff_expr_manual -= expr_b;
  EXPECT_EQ(diff_expr.hashmap, diff_expr_manual.hashmap);

  // Multiplication (similar checks)
  Expression prod_expr = expr_a * expr_b;
  Expression prod_expr_manual = expr_a;
  prod_expr_manual *= expr_b;
  EXPECT_EQ(prod_expr.hashmap, prod_expr_manual.hashmap);

  // Division (Expression / Complex)
  Expression div_complex = expr_a / complex_val;
  Expression div_complex_manual = expr_a;
  div_complex_manual /= complex_val;
  EXPECT_EQ(div_complex.hashmap, div_complex_manual.hashmap);
}
