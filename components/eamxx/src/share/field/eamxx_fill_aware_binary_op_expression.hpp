#ifndef EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP
#define EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_expression_binary_op.hpp>
#include <ekat_expression_base.hpp>

namespace ekat {

template<typename ELeft, typename ERight, BinOp OP>
class FillAwareBinaryExpression : public Expression<FillAwareBinaryExpression<ELeft,ERight,OP>>{
public:
  static constexpr bool expr_l = is_expr_v<ELeft>;
  static constexpr bool expr_r = is_expr_v<ERight>;

  using ret_left  = eval_t<ELeft>;
  using ret_right = eval_t<ERight>;

  static_assert (expr_l or expr_r,
    "[CmpExpression] At least one between ELeft and ERight must be an Expression type.\n");

  FillAwareBinaryExpression (const ELeft& left, const ERight& right)
    : m_left(left)
    , m_right(right)
    , m_fv_l (scream::constants::fill_value<ret_left>)
    , m_fv_r (scream::constants::fill_value<ret_right>)
  {
    // Nothing to do here
  }

  static constexpr int rank () { return BinaryExpression<ELeft,ERight,OP>::rank(); }
  int extent (int i) const {
    return BinaryExpression<ELeft,ERight,OP>(m_left,m_right).extent(i);
  }

  template<typename... Args>
  KOKKOS_INLINE_FUNCTION
  auto eval(Args... args) const {
    if constexpr (not expr_l) {
      return eval_impl(m_left,m_right.eval(args...));
    } else if constexpr (not expr_r) {
      return eval_impl(m_left.eval(args...),m_right);
    } else {
      return eval_impl(m_left.eval(args...),m_right.eval(args...));
    }
  }

  static auto ret_type () {
    return std::common_type_t<ret_left,ret_right>(0);
  }
protected:

  KOKKOS_INLINE_FUNCTION
  std::common_type_t<ret_left,ret_right>
  eval_impl (const ret_left& l, const ret_right& r) const {
    using ret_t = std::common_type_t<ret_left,ret_right>;
    if (l==m_fv_l)
      return static_cast<ret_t>(r);
    else if (r==m_fv_r)
      return static_cast<ret_t>(l);
    else
      if constexpr (OP==BinOp::Plus) {
        return static_cast<const ret_t&>(l)+static_cast<const ret_t&>(r);
      } else if constexpr (OP==BinOp::Minus) {
        return static_cast<const ret_t&>(l)-static_cast<const ret_t&>(r);
      } else if constexpr (OP==BinOp::Mult) {
        return static_cast<const ret_t&>(l)*static_cast<const ret_t&>(r);
      } else if constexpr (OP==BinOp::Div) {
        return static_cast<const ret_t&>(l)/static_cast<const ret_t&>(r);
      } else if constexpr (OP==BinOp::Max) {
        return Kokkos::max(static_cast<const ret_t&>(l),static_cast<const ret_t&>(r));
      } else if constexpr (OP==BinOp::Min) {
        return Kokkos::min(static_cast<const ret_t&>(l),static_cast<const ret_t&>(r));
      }
  }

  ELeft    m_left;
  ERight   m_right;

  ret_left  m_fv_l;
  ret_right m_fv_r;
};

template<typename ELeft, typename ERight, BinOp OP>
struct is_expr<FillAwareBinaryExpression<ELeft,ERight,OP>> : std::true_type {};

// Unary minus implemented as -1*expr
template<typename ERight>
FillAwareBinaryExpression<int,ERight,BinOp::Mult>
fa_neg (const Expression<ERight>& r)
{
  return FillAwareBinaryExpression<int,ERight,BinOp::Mult>(-1,r.cast());
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

// Overload max/min functions
template<typename ELeft, typename ERight>
std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Max>>
fa_max (const ELeft& l, const ERight& r)
{
  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Max>(l,r);
}

template<typename ELeft, typename ERight>
std::enable_if_t<is_expr_v<ELeft> or is_expr_v<ERight>,FillAwareBinaryExpression<ELeft,ERight,BinOp::Min>>
fa_min (const ELeft& l, const ERight& r)
{
  return FillAwareBinaryExpression<ELeft,ERight,BinOp::Min>(l,r);
}

} // namespace ekat

#endif // EAMXX_FILL_AWARE_BINARY_OP_EXPRESSION_HPP
