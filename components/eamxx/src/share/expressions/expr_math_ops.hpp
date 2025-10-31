#ifndef EAMXX_EXPR_MATH_OPS_HPP
#define EAMXX_EXPR_MATH_OPS_HPP

#include "share/expressions/base_expr.hpp"
#include "share/expressions/field_expr.hpp"
#include "share/expressions/scalar_expr.hpp"

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

// ---------------- Operator + ----------------- //

template<typename ERight, typename RetR>
SumExpression<RealFieldExpression,ERight>
operator+ (const Field& fl, const Expression<ERight,RetR>& r)
{
  check_real(fl);
  return RealFieldExpression(fl)+r;
}

template<typename ELeft, typename RetL>
SumExpression<ELeft,RealFieldExpression>
operator+ (const Expression<ELeft,RetL>& l, const Field& fr)
{
  check_real(fr);
  return l+RealFieldExpression(fr);
}

SumExpression<RealFieldExpression,RealFieldExpression>
operator+ (const Field& fl, const Field& fr)
{
  check_real(fr);
  check_real(fl);
  return RealFieldExpression(fl)+RealFieldExpression(fr);
}

} // namespace scream

#endif // EAMXX_EXPR_MATH_OPS_HPP
