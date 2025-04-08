#include "../src/Operator.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace libqm::testing {

using ::testing::StrEq;

using Type = Operator::Type;
using Spin = Operator::Spin;

class OperatorInterfaceTest : public ::testing::Test {
 protected:
  void SetUp() override {
    op_c_up_0_ = Operator::creation(Spin::Up, 0);
    op_c_up_5_ = Operator::creation(Spin::Up, 5);
    op_c_dn_5_ = Operator::creation(Spin::Down, 5);
    op_a_up_5_ = Operator::annihilation(Spin::Up, 5);
    op_a_dn_5_ = Operator::annihilation(Spin::Down, 5);
    op_a_dn_0_ = Operator::annihilation(Spin::Down, 0);
    op_c_up_max_ = Operator::creation(Spin::Up, Operator::max_index() - 1);
    op_a_dn_max_ = Operator::annihilation(Spin::Down, Operator::max_index() - 1);
  }

  Operator op_c_up_0_;
  Operator op_c_up_5_;
  Operator op_c_dn_5_;
  Operator op_a_up_5_;
  Operator op_a_dn_5_;
  Operator op_a_dn_0_;
  Operator op_c_up_max_;
  Operator op_a_dn_max_;
};

TEST_F(OperatorInterfaceTest, ParameterizedConstructorSetsProperties) {
  Operator op1(Type::Creation, Spin::Up, 0);
  EXPECT_EQ(Type::Creation, op1.type());
  EXPECT_EQ(Spin::Up, op1.spin());
  EXPECT_EQ(0, op1.value());

  Operator op2(Type::Annihilation, Spin::Down, 5);
  EXPECT_EQ(Type::Annihilation, op2.type());
  EXPECT_EQ(Spin::Down, op2.spin());
  EXPECT_EQ(5, op2.value());

  Operator op3(Type::Creation, Spin::Down, Operator::max_index() - 1);
  EXPECT_EQ(Type::Creation, op3.type());
  EXPECT_EQ(Spin::Down, op3.spin());
  EXPECT_EQ(Operator::max_index() - 1, op3.value());
}

TEST_F(OperatorInterfaceTest, StaticFactoryCreationSetsProperties) {
  Operator op1 = Operator::creation(Spin::Up, 5);
  EXPECT_EQ(Type::Creation, op1.type());
  EXPECT_EQ(Spin::Up, op1.spin());
  EXPECT_EQ(5, op1.value());
  EXPECT_EQ(op_c_up_5_, op1);

  Operator op2 = Operator::creation(Spin::Down, Operator::max_index() - 1);
  EXPECT_EQ(Type::Creation, op2.type());
  EXPECT_EQ(Spin::Down, op2.spin());
  EXPECT_EQ(Operator::max_index() - 1, op2.value());
}

TEST_F(OperatorInterfaceTest, StaticFactoryAnnihilationSetsProperties) {
  Operator op1 = Operator::annihilation(Spin::Down, 0);
  EXPECT_EQ(Type::Annihilation, op1.type());
  EXPECT_EQ(Spin::Down, op1.spin());
  EXPECT_EQ(0, op1.value());
  EXPECT_EQ(op_a_dn_0_, op1);

  Operator op2 = Operator::annihilation(Spin::Up, 10);
  EXPECT_EQ(Type::Annihilation, op2.type());
  EXPECT_EQ(Spin::Up, op2.spin());
  EXPECT_EQ(10, op2.value());
}

TEST_F(OperatorInterfaceTest, AdjointFlipsTypePreservesSpinAndValue) {
  Operator adjoint_c = op_c_up_5_.adjoint();
  EXPECT_EQ(Type::Annihilation, adjoint_c.type());
  EXPECT_EQ(op_c_up_5_.spin(), adjoint_c.spin());
  EXPECT_EQ(op_c_up_5_.value(), adjoint_c.value());
  EXPECT_EQ(op_a_up_5_, adjoint_c);

  Operator adjoint_a = op_a_dn_0_.adjoint();
  EXPECT_EQ(Type::Creation, adjoint_a.type());
  EXPECT_EQ(op_a_dn_0_.spin(), adjoint_a.spin());
  EXPECT_EQ(op_a_dn_0_.value(), adjoint_a.value());
  EXPECT_EQ(Operator::creation(Spin::Down, 0), adjoint_a);
}

TEST_F(OperatorInterfaceTest, AdjointIsSelfInverse) {
  EXPECT_EQ(op_c_up_5_, op_c_up_5_.adjoint().adjoint());
  EXPECT_EQ(op_a_dn_0_, op_a_dn_0_.adjoint().adjoint());
  EXPECT_EQ(op_a_dn_max_, op_a_dn_max_.adjoint().adjoint());
}

TEST_F(OperatorInterfaceTest, CommutationBehaviorFollowsFermionicRules) {
  // Rule 1: An operator commutes with itself.
  EXPECT_TRUE(op_c_up_5_.commutes(op_c_up_5_));
  EXPECT_TRUE(op_a_dn_0_.commutes(op_a_dn_0_));

  // Rule 2: An operator anticommutes (does NOT commute) with its exact adjoint.
  EXPECT_FALSE(op_c_up_5_.commutes(op_c_up_5_.adjoint()));
  EXPECT_FALSE(op_a_dn_0_.commutes(op_a_dn_0_.adjoint()));
  EXPECT_FALSE(op_c_up_5_.commutes(op_a_up_5_));
  EXPECT_FALSE(op_a_up_5_.commutes(op_c_up_5_));

  // Rule 3: Operators commute if they differ by more than just the type bit
  // (i.e., different spin, different value, or both).

  // Case 3a: Same type, different spin or value -> Commute
  EXPECT_TRUE(op_c_up_5_.commutes(op_c_dn_5_));
  EXPECT_TRUE(op_c_up_5_.commutes(op_c_up_0_));
  EXPECT_TRUE(op_a_dn_0_.commutes(op_a_dn_max_));
  EXPECT_TRUE(op_a_dn_0_.commutes(op_a_up_5_));

  // Case 3b: Different type, AND different spin or value -> Commute
  EXPECT_TRUE(op_c_up_5_.commutes(op_a_dn_5_));
  EXPECT_TRUE(op_c_up_5_.commutes(op_a_dn_0_));
  EXPECT_TRUE(op_c_up_0_.commutes(op_a_dn_max_));

  // Verify symmetry for commuting pairs
  EXPECT_TRUE(op_c_dn_5_.commutes(op_c_up_5_));
  EXPECT_TRUE(op_a_dn_5_.commutes(op_c_up_5_));
}

TEST_F(OperatorInterfaceTest, EqualityOperatorBehavesCorrectly) {
  EXPECT_EQ(op_c_up_5_, op_c_up_5_);

  Operator op_c_up_5_copy = Operator::creation(Spin::Up, 5);
  EXPECT_EQ(op_c_up_5_, op_c_up_5_copy);
  Operator op_c_up_5_copy2(Type::Creation, Spin::Up, 5);
  EXPECT_EQ(op_c_up_5_, op_c_up_5_copy2);

  EXPECT_NE(op_c_up_5_, op_a_up_5_);
  EXPECT_NE(op_c_up_5_, op_c_dn_5_);
  EXPECT_NE(op_c_up_5_, op_c_up_0_);
  EXPECT_NE(op_c_up_5_, op_a_dn_0_);
  EXPECT_NE(op_c_up_0_, op_a_dn_max_);
}

TEST_F(OperatorInterfaceTest, LessThanOperatorProvidesCorrectOrdering) {
  // Rule 1: Type (Creation < Annihilation)
  EXPECT_LT(op_c_up_0_, op_a_up_5_);
  EXPECT_LT(op_c_up_0_, op_a_dn_max_);
  EXPECT_TRUE(op_c_up_0_ < op_a_dn_max_);

  // Rule 2: Spin (Up < Down) when Type is the same
  EXPECT_LT(op_c_up_max_, op_c_dn_5_);
  EXPECT_LT(op_c_up_5_, op_c_dn_5_);
  EXPECT_FALSE(op_c_dn_5_ < op_c_up_5_);

  EXPECT_LT(op_a_up_5_, op_a_dn_0_);
  EXPECT_LT(op_a_up_5_, op_a_dn_max_);
  EXPECT_FALSE(op_a_dn_0_ < op_a_up_5_);

  // Rule 3: Value (ascending) when Type and Spin are the same
  EXPECT_LT(op_c_up_0_, op_c_up_5_);
  EXPECT_LT(op_c_up_5_, op_c_up_max_);
  EXPECT_FALSE(op_c_up_5_ < op_c_up_0_);

  EXPECT_LT(op_a_dn_0_, op_a_dn_max_);
  EXPECT_FALSE(op_a_dn_max_ < op_a_dn_0_);

  // Reflexivity: !(a < a)
  EXPECT_FALSE(op_c_up_5_ < op_c_up_5_);
  EXPECT_FALSE(op_a_dn_0_ < op_a_dn_0_);
}

TEST_F(OperatorInterfaceTest, ToStringFormatIsCorrect) {
  EXPECT_THAT(op_c_up_0_.to_string(), StrEq("c+(↑,0)"));
  EXPECT_THAT(op_c_up_5_.to_string(), StrEq("c+(↑,5)"));
  EXPECT_THAT(op_c_dn_5_.to_string(), StrEq("c+(↓,5)"));
  EXPECT_THAT(op_a_up_5_.to_string(), StrEq("c(↑,5)"));
  EXPECT_THAT(op_a_dn_0_.to_string(), StrEq("c(↓,0)"));
}
}  // namespace libqm::testing
