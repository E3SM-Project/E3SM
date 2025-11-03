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

  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval (int i) const {
    if constexpr (OP==BinOp::Plus) {
      return static_cast<T>(m_left.template eval<T>(i) + m_right.template eval<T>(i));
    } else if constexpr (OP==BinOp::Minus) {
      return static_cast<T>(m_left.template eval<T>(i) - m_right.template eval<T>(i));
    } else if constexpr (OP==BinOp::Mult) {
      return static_cast<T>(m_left.template eval<T>(i) * m_right.template eval<T>(i));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<T>(m_left.template eval<T>(i) / m_right.template eval<T>(i));
    } else if constexpr (OP==BinOp::Max) {
      return static_cast<T>(Kokkos::max(m_left.template eval<T>(i),m_right.template eval<T>(i)));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<T>(Kokkos::min(m_left.template eval<T>(i),m_right.template eval<T>(i)));
    }
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval (int i, int j) const {
    if constexpr (OP==BinOp::Plus) {
      return static_cast<T>(m_left.template eval<T>(i,j) + m_right.template eval<T>(i,j));
    } else if constexpr (OP==BinOp::Minus) {
      return static_cast<T>(m_left.template eval<T>(i,j) - m_right.template eval<T>(i,j));
    } else if constexpr (OP==BinOp::Mult) {
      return static_cast<T>(m_left.template eval<T>(i,j) * m_right.template eval<T>(i,j));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<T>(m_left.template eval<T>(i,j) / m_right.template eval<T>(i,j));
    } else if constexpr (OP==BinOp::Max) {
      return static_cast<T>(Kokkos::max(m_left.template eval<T>(i,j),m_right.template eval<T>(i,j)));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<T>(Kokkos::min(m_left.template eval<T>(i,j),m_right.template eval<T>(i,j)));
    }
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval (int i, int j, int k) const {
    if constexpr (OP==BinOp::Plus) {
      return static_cast<T>(m_left.template eval<T>(i,j,k) + m_right.template eval<T>(i,j,k));
    } else if constexpr (OP==BinOp::Minus) {
      return static_cast<T>(m_left.template eval<T>(i,j,k) - m_right.template eval<T>(i,j,k));
    } else if constexpr (OP==BinOp::Mult) {
      return static_cast<T>(m_left.template eval<T>(i,j,k) * m_right.template eval<T>(i,j,k));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<T>(m_left.template eval<T>(i,j,k) / m_right.template eval<T>(i,j,k));
    } else if constexpr (OP==BinOp::Max) {
      return static_cast<T>(Kokkos::max(m_left.template eval<T>(i,j,k),m_right.template eval<T>(i,j,k)));
    } else if constexpr (OP==BinOp::Div) {
      return static_cast<T>(Kokkos::min(m_left.template eval<T>(i,j,k),m_right.template eval<T>(i,j,k)));
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
