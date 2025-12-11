#include <gtest/gtest.h>
#include <rational_linalg/fast_rational.hpp>
#include <limits>
#include <sstream>
#include <unordered_set>

// Test Construction
TEST(Rational64Test, DefaultConstruction) {
    rational64 r;
    EXPECT_EQ(r.numerator(), 0);
    EXPECT_EQ(r.denominator(), 1);
}

TEST(Rational64Test, IntegerConstruction) {
    rational64 r(42);
    EXPECT_EQ(r.numerator(), 42);
    EXPECT_EQ(r.denominator(), 1);
}

TEST(Rational64Test, FractionConstruction) {
    rational64 r(3, 4);
    EXPECT_EQ(r.numerator(), 3);
    EXPECT_EQ(r.denominator(), 4);
}

TEST(Rational64Test, NormalizedConstruction) {
    rational64 r(6, 8);
    EXPECT_EQ(r.numerator(), 3);
    EXPECT_EQ(r.denominator(), 4);
}

TEST(Rational64Test, NegativeDenominatorNormalization) {
    rational64 r(3, -4);
    EXPECT_EQ(r.numerator(), -3);
    EXPECT_EQ(r.denominator(), 4);
}

TEST(Rational64Test, ZeroDenominatorThrows) {
    EXPECT_THROW(rational64(1, 0), rational_overflow);
}

TEST(Rational64Test, ZeroNumerator) {
    rational64 r(0, 5);
    EXPECT_EQ(r.numerator(), 0);
    EXPECT_EQ(r.denominator(), 1);
}

// Test Arithmetic Operations
TEST(Rational64Test, Addition) {
    rational64 a(1, 2);
    rational64 b(1, 3);
    rational64 result = a + b;
    EXPECT_EQ(result.numerator(), 5);
    EXPECT_EQ(result.denominator(), 6);
}

TEST(Rational64Test, Subtraction) {
    rational64 a(1, 2);
    rational64 b(1, 3);
    rational64 result = a - b;
    EXPECT_EQ(result.numerator(), 1);
    EXPECT_EQ(result.denominator(), 6);
}

TEST(Rational64Test, Multiplication) {
    rational64 a(2, 3);
    rational64 b(3, 4);
    rational64 result = a * b;
    EXPECT_EQ(result.numerator(), 1);
    EXPECT_EQ(result.denominator(), 2);
}

TEST(Rational64Test, Division) {
    rational64 a(2, 3);
    rational64 b(4, 5);
    rational64 result = a / b;
    EXPECT_EQ(result.numerator(), 5);
    EXPECT_EQ(result.denominator(), 6);
}

TEST(Rational64Test, DivisionByZeroThrows) {
    rational64 a(1, 2);
    rational64 b(0, 1);
    EXPECT_THROW(a / b, rational_overflow);
}

// Test Integer Arithmetic
TEST(Rational64Test, AdditionWithInteger) {
    rational64 a(1, 2);
    rational64 result = a + 3;
    EXPECT_EQ(result.numerator(), 7);
    EXPECT_EQ(result.denominator(), 2);
}

TEST(Rational64Test, IntegerAdditionWithRational) {
    rational64 a(1, 2);
    rational64 result = 3 + a;
    EXPECT_EQ(result.numerator(), 7);
    EXPECT_EQ(result.denominator(), 2);
}

TEST(Rational64Test, SubtractionWithInteger) {
    rational64 a(7, 2);
    rational64 result = a - 3;
    EXPECT_EQ(result.numerator(), 1);
    EXPECT_EQ(result.denominator(), 2);
}

TEST(Rational64Test, MultiplicationWithInteger) {
    rational64 a(2, 3);
    rational64 result = a * 3;
    EXPECT_EQ(result.numerator(), 2);
    EXPECT_EQ(result.denominator(), 1);
}

TEST(Rational64Test, DivisionWithInteger) {
    rational64 a(6, 1);
    rational64 result = a / 3;
    EXPECT_EQ(result.numerator(), 2);
    EXPECT_EQ(result.denominator(), 1);
}

// Test Compound Assignment
TEST(Rational64Test, CompoundAssignmentAddition) {
    rational64 a(1, 2);
    a += rational64(1, 3);
    EXPECT_EQ(a.numerator(), 5);
    EXPECT_EQ(a.denominator(), 6);
}

TEST(Rational64Test, CompoundAssignmentSubtraction) {
    rational64 a(1, 2);
    a -= rational64(1, 3);
    EXPECT_EQ(a.numerator(), 1);
    EXPECT_EQ(a.denominator(), 6);
}

TEST(Rational64Test, CompoundAssignmentMultiplication) {
    rational64 a(2, 3);
    a *= rational64(3, 4);
    EXPECT_EQ(a.numerator(), 1);
    EXPECT_EQ(a.denominator(), 2);
}

TEST(Rational64Test, CompoundAssignmentDivision) {
    rational64 a(2, 3);
    a /= rational64(4, 5);
    EXPECT_EQ(a.numerator(), 5);
    EXPECT_EQ(a.denominator(), 6);
}

// Test Comparisons
TEST(Rational64Test, Equality) {
    rational64 a(1, 2);
    rational64 b(2, 4);
    EXPECT_EQ(a, b);
}

TEST(Rational64Test, Inequality) {
    rational64 a(1, 2);
    rational64 b(1, 3);
    EXPECT_NE(a, b);
}

TEST(Rational64Test, LessThan) {
    rational64 a(1, 3);
    rational64 b(1, 2);
    EXPECT_LT(a, b);
}

TEST(Rational64Test, GreaterThan) {
    rational64 a(1, 2);
    rational64 b(1, 3);
    EXPECT_GT(a, b);
}

TEST(Rational64Test, LessThanOrEqual) {
    rational64 a(1, 2);
    rational64 b(2, 4);
    rational64 c(1, 3);
    EXPECT_LE(a, b);
    EXPECT_LE(c, a);
}

TEST(Rational64Test, GreaterThanOrEqual) {
    rational64 a(1, 2);
    rational64 b(2, 4);
    rational64 c(1, 3);
    EXPECT_GE(a, b);
    EXPECT_GE(a, c);
}

TEST(Rational64Test, NegativeComparison) {
    rational64 a(-1, 2);
    rational64 b(1, 2);
    EXPECT_LT(a, b);
    EXPECT_LT(a, rational64(0));
}

// Test Integer Comparisons
TEST(Rational64Test, EqualityWithInteger) {
    rational64 a(5, 1);
    EXPECT_EQ(a, 5);
    EXPECT_EQ(5, a);
}

TEST(Rational64Test, InequalityWithInteger) {
    rational64 a(5, 2);
    EXPECT_NE(a, 5);
    EXPECT_NE(5, a);
}

TEST(Rational64Test, LessThanWithInteger) {
    rational64 a(4, 1);
    EXPECT_LT(a, 5);
    EXPECT_LT(3, a);
}

// Test Unary Operators
TEST(Rational64Test, UnaryPlus) {
    rational64 a(1, 2);
    rational64 result = +a;
    EXPECT_EQ(result.numerator(), 1);
    EXPECT_EQ(result.denominator(), 2);
}

TEST(Rational64Test, UnaryMinus) {
    rational64 a(1, 2);
    rational64 result = -a;
    EXPECT_EQ(result.numerator(), -1);
    EXPECT_EQ(result.denominator(), 2);
}

TEST(Rational64Test, UnaryMinusNegative) {
    rational64 a(-1, 2);
    rational64 result = -a;
    EXPECT_EQ(result.numerator(), 1);
    EXPECT_EQ(result.denominator(), 2);
}

// Test Utility Functions
TEST(Rational64Test, Abs) {
    rational64 a(-3, 4);
    rational64 result = a.abs();
    EXPECT_EQ(result.numerator(), 3);
    EXPECT_EQ(result.denominator(), 4);
}

TEST(Rational64Test, AbsPositive) {
    rational64 a(3, 4);
    rational64 result = a.abs();
    EXPECT_EQ(result.numerator(), 3);
    EXPECT_EQ(result.denominator(), 4);
}

TEST(Rational64Test, Inverse) {
    rational64 a(3, 4);
    rational64 result = a.inverse();
    EXPECT_EQ(result.numerator(), 4);
    EXPECT_EQ(result.denominator(), 3);
}

TEST(Rational64Test, InverseZeroThrows) {
    rational64 a(0, 1);
    EXPECT_THROW(a.inverse(), rational_overflow);
}

TEST(Rational64Test, IsZero) {
    rational64 a(0, 1);
    rational64 b(0, 5);
    rational64 c(1, 2);
    EXPECT_TRUE(a.is_zero());
    EXPECT_TRUE(b.is_zero());
    EXPECT_FALSE(c.is_zero());
}

TEST(Rational64Test, IsPositive) {
    rational64 a(1, 2);
    rational64 b(-1, 2);
    rational64 c(0, 1);
    EXPECT_TRUE(a.is_positive());
    EXPECT_FALSE(b.is_positive());
    EXPECT_FALSE(c.is_positive());
}

TEST(Rational64Test, IsNegative) {
    rational64 a(1, 2);
    rational64 b(-1, 2);
    rational64 c(0, 1);
    EXPECT_FALSE(a.is_negative());
    EXPECT_TRUE(b.is_negative());
    EXPECT_FALSE(c.is_negative());
}

TEST(Rational64Test, IsInteger) {
    rational64 a(5, 1);
    rational64 b(5, 2);
    EXPECT_TRUE(a.is_integer());
    EXPECT_FALSE(b.is_integer());
}

TEST(Rational64Test, ToDouble) {
    rational64 a(1, 2);
    EXPECT_DOUBLE_EQ(a.to_double(), 0.5);
}

TEST(Rational64Test, ToDoubleLarge) {
    rational64 a(22, 7);
    EXPECT_NEAR(a.to_double(), 3.142857142857, 1e-10);
}

// Test String Conversion
TEST(Rational64Test, ToStringInteger) {
    rational64 a(5, 1);
    EXPECT_EQ(to_string(a), "5");
}

TEST(Rational64Test, ToStringFraction) {
    rational64 a(3, 4);
    EXPECT_EQ(to_string(a), "3/4");
}

// Test Stream Operators
TEST(Rational64Test, StreamOutput) {
    rational64 a(3, 4);
    std::ostringstream oss;
    oss << a;
    EXPECT_EQ(oss.str(), "3/4");
}

TEST(Rational64Test, StreamOutputInteger) {
    rational64 a(5, 1);
    std::ostringstream oss;
    oss << a;
    EXPECT_EQ(oss.str(), "5");
}

TEST(Rational64Test, StreamInput) {
    rational64 r;
    std::istringstream iss("3/4");
    iss >> r;
    EXPECT_EQ(r.numerator(), 3);
    EXPECT_EQ(r.denominator(), 4);
}

TEST(Rational64Test, StreamInputInteger) {
    rational64 r;
    std::istringstream iss("5");
    iss >> r;
    EXPECT_EQ(r.numerator(), 5);
    EXPECT_EQ(r.denominator(), 1);
}

TEST(Rational64Test, StreamInputZeroDenominator) {
    rational64 r;
    std::istringstream iss("1/0");
    iss >> r;
    EXPECT_TRUE(iss.fail());
}

// Test Hash Function
TEST(Rational64Test, HashConsistency) {
    rational64 a(3, 4);
    rational64 b(3, 4);
    std::hash<rational64> hasher;
    EXPECT_EQ(hasher(a), hasher(b));
}

TEST(Rational64Test, HashUnorderedSet) {
    std::unordered_set<rational64> s;
    s.insert(rational64(1, 2));
    s.insert(rational64(3, 4));
    s.insert(rational64(1, 2)); // Duplicate
    EXPECT_EQ(s.size(), 2);
}

// Test Static Constants
TEST(Rational64Test, StaticConstants) {
    EXPECT_EQ(rational64::ZERO.numerator(), 0);
    EXPECT_EQ(rational64::ZERO.denominator(), 1);
    EXPECT_EQ(rational64::ONE.numerator(), 1);
    EXPECT_EQ(rational64::ONE.denominator(), 1);
    EXPECT_EQ(rational64::MINUS_ONE.numerator(), -1);
    EXPECT_EQ(rational64::MINUS_ONE.denominator(), 1);
}

// Test Overflow Cases (if possible without actually overflowing)
TEST(Rational64Test, LargeValues) {
    rational64 a(INT64_MAX / 2, 1);
    rational64 b(1, 2);
    rational64 result = a * b;
    // Should not throw if within range
    EXPECT_GT(result.numerator(), 0);
}

// Test Edge Cases
TEST(Rational64Test, OneOverOne) {
    rational64 a(1, 1);
    EXPECT_EQ(a.numerator(), 1);
    EXPECT_EQ(a.denominator(), 1);
}

TEST(Rational64Test, NegativeOneOverOne) {
    rational64 a(-1, 1);
    EXPECT_EQ(a.numerator(), -1);
    EXPECT_EQ(a.denominator(), 1);
}

TEST(Rational64Test, LargeDenominator) {
    rational64 a(1, INT64_MAX);
    EXPECT_EQ(a.numerator(), 1);
    EXPECT_EQ(a.denominator(), INT64_MAX);
}

