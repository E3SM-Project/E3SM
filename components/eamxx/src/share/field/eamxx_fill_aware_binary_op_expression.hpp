#ifndef EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP
#define EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_expression_binary_op.hpp>
#include <ekat_expression_meta.hpp>

namespace ekat {

template<typename ELeft, typename ERight, BinOp OP>
class FillAwareBinaryExpression : public BinaryExpression<ELeft,ERight,OP> {
public:
  using base_t = BinaryExpression<ELeft,ERight,OP>;

  static constexpr bool expr_l = is_expr_v<ELeft>;
  static constexpr bool expr_r = is_expr_v<ERight>;

  using eval_left_t  = typename base_t::eval_left_t;
  using eval_right_t = typename base_t::eval_right_t;
  using eval_t       = typename base_t::eval_t;

  FillAwareBinaryExpression (const ELeft& left, const ERight& right)
    : base_t(left,right)
    , m_fv_l (scream::constants::fill_value<eval_left_t>)
    , m_fv_r (scream::constants::fill_value<eval_right_t>)
    , m_fv   (scream::constants::fill_value<eval_t>)
  {
    // Nothing to do here
  }

  template<typename... Args>
  KOKKOS_INLINE_FUNCTION
  auto eval(Args... args) const {
    if constexpr (not expr_l) {
      return eval_impl(this->m_left,this->m_right.eval(args...));
    } else if constexpr (not expr_r) {
      return eval_impl(this->m_left.eval(args...),this->m_right);
    } else {
      return eval_impl(this->m_left.eval(args...),this->m_right.eval(args...));
    }
  }

protected:

  KOKKOS_INLINE_FUNCTION
  eval_t eval_impl (const eval_left_t& l, const eval_right_t& r) const {
    if (l==m_fv_l or r==m_fv_l)
      return m_fv;
    else 
      if constexpr (OP==BinOp::Plus) {
        return l+r;
      } else if constexpr (OP==BinOp::Minus) {
        return l-r;
      } else if constexpr (OP==BinOp::Mult) {
        return l*r;
      } else {
        return l/r;
      }
  }

  eval_left_t  m_fv_l;
  eval_right_t m_fv_r;
  eval_t       m_fv;
};

template<typename ELeft, typename ERight, BinOp OP>
struct is_expr<FillAwareBinaryExpression<ELeft,ERight,OP>> : std::true_type {};
template<typename ELeft, typename ERight, BinOp OP>
struct eval_return<FillAwareBinaryExpression<ELeft,ERight,OP>> : eval_return<BinaryExpression<ELeft,ERight,OP>> {};

// Unary minus implemented as -1*expr
template<typename ERight>
std::enable_if_t<is_expr_v<ERight>,FillAwareBinaryExpression<int,ERight,BinOp::Mult>>
fa_neg (const ERight& r)
{
  return FillAwareBinaryExpression<int,ERight,BinOp::Mult>(-1,r);
}

// Overload arithmetic operators
template<typename ELeft, typename ERight>
std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Plus>>
fa_plus (const ELeft& l, const ERight& r)
{
  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Plus>(l,r);
}

template<typename ELeft, typename ERight>
std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Minus>>
fa_minus (const ELeft& l, const ERight& r)
{
  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Minus>(l,r);
}

template<typename ELeft, typename ERight>
std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Mult>>
fa_mult (const ELeft& l, const ERight& r)
{
  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Mult>(l,r);
}

template<typename ELeft, typename ERight>
std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Div>>
fa_div (const ELeft& l, const ERight& r)
{
  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Div>(l,r);
}

//// Overload max/min functions
//template<typename ELeft, typename ERight>
//std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Max>>
//fa_max (const ELeft& l, const ERight& r)
//{
//  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Max>(l,r);
//}

//template<typename ELeft, typename ERight>
//std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Min>>
//fa_min (const ELeft& l, const ERight& r)
//{
//  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Min>(l,r);
//}

} // namespace ekat

#endif // EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP
