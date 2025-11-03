#ifndef EAMXX_BIN_OP_EXPRESSION_HPP
#define EAMXX_BIN_OP_EXPRESSION_HPP

#include "share/expressions/base.hpp"

namespace scream {

enum class BinOp {
  Plus,
  Minus,
  Mult,
  Div,
  Max,
  Min
};

template<typename ELeft, typename ERight, BinOp OP>
class BinaryExpression : public Expression<BinaryExpression<ELeft,ERight,OP>>{
public:
  static constexpr bool is_leaf = false;

  BinaryExpression (const ELeft& left, const ERight& right)
    : m_left(left)
    , m_right(right)
  {
    // Nothing to do here
  }

  int num_indices () const { return std::max(m_left.num_indices(),m_right.num_indices()); }

  // void set_eval_layout (const FieldLayout& fl) {}

  KOKKOS_INLINE_FUNCTION
  Real eval (int i) const {
    if constexpr (OP==BinOp::Plus) {
      return m_left.eval(i) + m_right.eval(i);
    } else if constexpr (OP==BinOp::Minus) {
      return m_left.eval(i) - m_right.eval(i);
    } else if constexpr (OP==BinOp::Mult) {
      return m_left.eval(i) * m_right.eval(i);
    } else if constexpr (OP==BinOp::Div) {
      return m_left.eval(i) / m_right.eval(i);
    } else if constexpr (OP==BinOp::Max) {
      return Kokkos::max(m_left.eval(i),m_right.eval(i));
    } else if constexpr (OP==BinOp::Div) {
      return Kokkos::min(m_left.eval(i),m_right.eval(i));
    }
  }
  KOKKOS_INLINE_FUNCTION
  Real eval (int i, int j) const {
    if constexpr (OP==BinOp::Plus) {
      return m_left.eval(i,j) + m_right.eval(i,j);
    } else if constexpr (OP==BinOp::Minus) {
      return m_left.eval(i,j) - m_right.eval(i,j);
    } else if constexpr (OP==BinOp::Mult) {
      return m_left.eval(i,j) * m_right.eval(i,j);
    } else if constexpr (OP==BinOp::Div) {
      return m_left.eval(i,j) / m_right.eval(i,j);
    } else if constexpr (OP==BinOp::Max) {
      return Kokkos::max(m_left.eval(i,j),m_right.eval(i,j));
    } else if constexpr (OP==BinOp::Div) {
      return Kokkos::min(m_left.eval(i,j),m_right.eval(i,j));
    }
  }
  KOKKOS_INLINE_FUNCTION
  Real eval (int i, int j, int k) const {
    if constexpr (OP==BinOp::Plus) {
      return m_left.eval(i,j,k) + m_right.eval(i,j,k);
    } else if constexpr (OP==BinOp::Minus) {
      return m_left.eval(i,j,k) - m_right.eval(i,j,k);
    } else if constexpr (OP==BinOp::Mult) {
      return m_left.eval(i,j,k) * m_right.eval(i,j,k);
    } else if constexpr (OP==BinOp::Div) {
      return m_left.eval(i,j,k) / m_right.eval(i,j,k);
    } else if constexpr (OP==BinOp::Max) {
      return Kokkos::max(m_left.eval(i,j,k),m_right.eval(i,j,k));
    } else if constexpr (OP==BinOp::Div) {
      return Kokkos::min(m_left.eval(i,j,k),m_right.eval(i,j,k));
    }
  }

protected:

  const ELeft    m_left;
  const ERight   m_right;
};

template<typename ELeft, typename ERight>
BinaryExpression<ELeft,ERight,BinOp::Plus>
operator+ (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Plus>(l.cast(),r.cast());
}

template<typename ELeft, typename ERight>
BinaryExpression<ELeft,ERight,BinOp::Minus>
operator- (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Minus>(l.cast(),r.cast());
}

template<typename ELeft, typename ERight>
BinaryExpression<ELeft,ERight,BinOp::Mult>
operator* (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Mult>(l.cast(),r.cast());
}

template<typename ELeft, typename ERight>
BinaryExpression<ELeft,ERight,BinOp::Div>
operator/ (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return BinaryExpression<ELeft,ERight,BinOp::Div>(l.cast(),r.cast());
}

} // namespace scream

#endif // EAMXX_BIN_OP_EXPRESSION_HPP
