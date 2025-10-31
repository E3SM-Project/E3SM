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

template<typename ERight, typename RetR>
BinaryExpression<RealFieldExpression,ERight,BinOp::Plus>
operator+ (const Field& fl, const Expression<ERight,RetR>& r)
{
  check_real(fl);
  return RealFieldExpression(fl)+r;
}

template<typename ELeft, typename RetL>
BinaryExpression<ERight,RealFieldExpression,BinOp::Plus>
operator+ (const Expression<ELeft,RetL>& l, const Field& fr)
{
  check_real(fr);
  return l+RealFieldExpression(fr);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Plus>
operator+ (const Field& fl, const Field& fr)
{
  check_real(fr);
  check_real(fl);
  return RealFieldExpression(fl)+RealFieldExpression(fr);
}

// ---------------- Minus ----------------- //

template<typename ERight, typename RetR>
BinaryExpression<RealFieldExpression,ERight,BinOp::Minus>
operator- (const Field& fl, const Expression<ERight,RetR>& r)
{
  check_real(fl);
  return RealFieldExpression(fl)+r;
}

template<typename ELeft, typename RetL>
BinaryExpression<ERight,RealFieldExpression,BinOp::Minus>
operator- (const Expression<ELeft,RetL>& l, const Field& fr)
{
  check_real(fr);
  return l+RealFieldExpression(fr);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Minus>
operator- (const Field& fl, const Field& fr)
{
  check_real(fr);
  check_real(fl);
  return RealFieldExpression(fl)+RealFieldExpression(fr);
}

// ---------------- Mult ----------------- //

template<typename ERight, typename RetR>
BinaryExpression<RealFieldExpression,ERight,BinOp::Mult>
operator* (const Field& fl, const Expression<ERight,RetR>& r)
{
  check_real(fl);
  return RealFieldExpression(fl)+r;
}

template<typename ELeft, typename RetL>
BinaryExpression<ERight,RealFieldExpression,BinOp::Mult>
operator* (const Expression<ELeft,RetL>& l, const Field& fr)
{
  check_real(fr);
  return l+RealFieldExpression(fr);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Mult>
operator* (const Field& fl, const Field& fr)
{
  check_real(fr);
  check_real(fl);
  return RealFieldExpression(fl)+RealFieldExpression(fr);
}

// ---------------- Div ----------------- //

template<typename ERight, typename RetR>
BinaryExpression<RealFieldExpression,ERight,BinOp::Div>
operator/ (const Field& fl, const Expression<ERight,RetR>& r)
{
  check_real(fl);
  return RealFieldExpression(fl)+r;
}

template<typename ELeft, typename RetL>
BinaryExpression<ERight,RealFieldExpression,BinOp::Div>
operator/ (const Expression<ELeft,RetL>& l, const Field& fr)
{
  check_real(fr);
  return l+RealFieldExpression(fr);
}

BinaryExpression<RealFieldExpression,RealFieldExpression,BinOp::Div>
operator/ (const Field& fl, const Field& fr)
{
  check_real(fr);
  check_real(fl);
  return RealFieldExpression(fl)+RealFieldExpression(fr);
}

} // namespace scream

#endif // EAMXX_EXPRESSION_HELPERS_HPP
