#ifndef EAMXX_EXPRESSION_HELPERS_HPP
#define EAMXX_EXPRESSION_HELPERS_HPP

#include "share/expressions/base.hpp"
#include "share/expressions/field.hpp"
#include "share/expressions/scalar.hpp"
#include "share/expressions/binary_op.hpp"

#include "share/field/field.hpp"

namespace scream {

// Support the case where LHS and/or RHS of a binary op is a Field or scalar
// NOTE: we ONLY support the case of real-valued field here.
// If you have an int-valued field, you have to manually
// create the FieldEvaluate expression

inline void check_real (const Field& f)
{
  EKAT_REQUIRE_MSG (f.data_type()==DataType::RealType,
      "Error! Overload of math ops to create expressions from fields is only offered for Real valued fields.\n"
      "To create an expression from non-real valued fields, create the FieldExpressionBase manually first.\n"
      " - field name: " + f.name() << "\n"
      " - data type : " + e2str(f.data_type()) + "\n");
}

// ---------------- Plus ----------------- //

template<typename ERight>
BinaryExpression<RealFieldExpression,ERight,BinOp::Plus>
operator+ (const Field& left, const Expression<ERight>& right)
{
  check_real(left);
  return RealFieldExpression(left)+right;
}

template<typename ELeft>
BinaryExpression<ELeft,RealFieldExpression,BinOp::Plus>
operator+ (const Expression<ELeft>& left, const Field& right)
{
  check_real(right);
  return left+RealFieldExpression(right);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Plus>
operator+ (const Field& left, const Field& right)
{
  check_real(right);
  check_real(left);
  return RealFieldExpression(left)+RealFieldExpression(right);
}

// ---------------- Minus ----------------- //

template<typename ERight>
BinaryExpression<RealFieldExpression,ERight,BinOp::Minus>
operator- (const Field& left, const Expression<ERight>& right)
{
  check_real(left);
  return RealFieldExpression(left)-right;
}

template<typename ELeft>
BinaryExpression<ELeft,RealFieldExpression,BinOp::Minus>
operator- (const Expression<ELeft>& left, const Field& right)
{
  check_real(right);
  return left-RealFieldExpression(right);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Minus>
operator- (const Field& left, const Field& right)
{
  check_real(right);
  check_real(left);
  return RealFieldExpression(left)-RealFieldExpression(right);
}

// ---------------- Mult ----------------- //

template<typename ERight>
BinaryExpression<RealFieldExpression,ERight,BinOp::Mult>
operator* (const Field& left, const Expression<ERight>& right)
{
  check_real(left);
  return RealFieldExpression(left)*right;
}

template<typename ELeft>
BinaryExpression<ELeft,RealFieldExpression,BinOp::Mult>
operator* (const Expression<ELeft>& left, const Field& right)
{
  check_real(right);
  return left*RealFieldExpression(right);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Mult>
operator* (const Field& left, const Field& right)
{
  check_real(right);
  check_real(left);
  return RealFieldExpression(left)*RealFieldExpression(right);
}

// ---------------- Div ----------------- //

template<typename ERight>
BinaryExpression<RealFieldExpression,ERight,BinOp::Div>
operator/ (const Field& left, const Expression<ERight>& right)
{
  check_real(left);
  return RealFieldExpression(left)/right;
}

template<typename ELeft>
BinaryExpression<ELeft,RealFieldExpression,BinOp::Div>
operator/ (const Expression<ELeft>& left, const Field& right)
{
  check_real(right);
  return left/RealFieldExpression(right);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Div>
operator/ (const Field& left, const Field& right)
{
  check_real(right);
  check_real(left);
  return RealFieldExpression(left)/RealFieldExpression(right);
}

} // namespace scream

#endif // EAMXX_EXPRESSION_HELPERS_HPP
