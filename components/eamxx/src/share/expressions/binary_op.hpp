#ifndef EAMXX_BIN_OP_EXPRESSION_HPP
#define EAMXX_BIN_OP_EXPRESSION_HPP

#include "share/expressions/base_expr.hpp"

namespace scream {

namespace impl {

template<typename ELeft, typename ERight>
using bin_op_ret_t = decltype(std::declval<typename ELeft::ret_t>() + 
                              std::declval<typename ERight::ret_t>());
}

enum class BinOp {
  Plus,
  Minus,
  Mult,
  Div,
  Max,
  Min
};

template<typename ELeft, typename ERight, BinOp OP>
class BinaryExpression : public Expression<BinaryExpression<ELeft,ERight,OP>,impl::bin_op_ret_t<ELeft,ERight>>{
public:
  using ret_left  = typename ELeft::ret_t;
  using ret_right = typename ERight::ret_t;
  using ret_t     = impl::bin_op_ret_t<ELeft,ERight>;
  static constexpr bool is_leaf = false;

  BinaryExpression (const ELeft& left, const ERight& right)
    : m_left(left)
    , m_right(right)
  {
    // Nothing to do here
  }

  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i) const {
    if constexpr (OP==BinOp::Plus) {
      return static_cast<ret_t>(m_left(i) + m_right(i));
    } else if constexpr (OP==BinOp::Minus) {
      return static_cast<ret_t>(m_left(i) - m_right(i));
    } else if constexpr (OP==BinOp::Mult) {
      return static_cast<ret_t>(m_left(i) * m_right(i));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<ret_t>(m_left(i) / m_right(i));
    }
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i, int j) const {
    if constexpr (OP==BinOp::Plus) {
      return static_cast<ret_t>(m_left(i,j) + m_right(i,j));
    } else if constexpr (OP==BinOp::Minus) {
      return static_cast<ret_t>(m_left(i,j) - m_right(i,j));
    } else if constexpr (OP==BinOp::Mult) {
      return static_cast<ret_t>(m_left(i,j) * m_right(i,j));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<ret_t>(m_left(i,j) / m_right(i,j));
    }
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i, int j, int k) const {
    if constexpr (OP==BinOp::Plus) {
      return static_cast<ret_t>(m_left(i,j,k) + m_right(i,j,k));
    } else if constexpr (OP==BinOp::Minus) {
      return static_cast<ret_t>(m_left(i,j,k) - m_right(i,j,k));
    } else if constexpr (OP==BinOp::Mult) {
      return static_cast<ret_t>(m_left(i,j,k) * m_right(i,j,k));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<ret_t>(m_left(i,j,k) / m_right(i,j,k));
    }
  }

protected:

  const ELeft    m_left;
  const ERight   m_right;
};

template<typename ELeft, typename RetL, typename ERight, typename RetR>
BinaryExpression<ELeft,ERight,BinOp::Plus>
operator+ (const Expression<ELeft,RetL>& l, const Expression<ERight,RetR>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Plus>(l.cast(),r.cast());
}

template<typename ELeft, typename RetL, typename ERight, typename RetR>
BinaryExpression<ELeft,ERight,BinOp::Minus>
operator- (const Expression<ELeft,RetL>& l, const Expression<ERight,RetR>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Minus>(l.cast(),r.cast());
}

template<typename ELeft, typename RetL, typename ERight, typename RetR>
BinaryExpression<ELeft,ERight,BinOp::Mult>
operator* (const Expression<ELeft,RetL>& l, const Expression<ERight,RetR>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Mult>(l.cast(),r.cast());
}

template<typename ELeft, typename RetL, typename ERight, typename RetR>
BinaryExpression<ELeft,ERight,BinOp::Div>
operator/ (const Expression<ELeft,RetL>& l, const Expression<ERight,RetR>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Div>(l.cast(),r.cast());
}

} // namespace scream

#endif // EAMXX_BIN_OP_EXPRESSION_HPP
