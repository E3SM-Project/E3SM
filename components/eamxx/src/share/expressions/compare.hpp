#ifndef EAMXX_CMP_EXPRESSION_HPP
#define EAMXX_CMP_EXPRESSION_HPP

#include "share/expressions/base.hpp"

#include "share/core/eamxx_types.hpp" // For Comparison

#include <ekat_kernel_assert.hpp>

namespace scream {

template<typename ELeft, typename ERight>
class CmpExpression : public Expression<CmpExpression<ELeft,ERight>> {
public:
  static constexpr bool is_leaf = false;

  CmpExpression (const ELeft& left, const ERight& right, Comparison CMP)
    : m_left(left)
    , m_right(right)
    , m_cmp(CMP)
  {
    // Nothing to do here
  }

  template<typename T>
  KOKKOS_INLINE_FUNCTION
  int eval (int i) const {
    switch (m_cmp) {
      case Comparison::EQ: return m_left.template eval<T>(i) == m_right.template eval<T>(i);
      case Comparison::NE: return m_left.template eval<T>(i) != m_right.template eval<T>(i);
      case Comparison::GT: return m_left.template eval<T>(i) >  m_right.template eval<T>(i);
      case Comparison::GE: return m_left.template eval<T>(i) >= m_right.template eval<T>(i);
      case Comparison::LT: return m_left.template eval<T>(i) <  m_right.template eval<T>(i);
      case Comparison::LE: return m_left.template eval<T>(i) <= m_right.template eval<T>(i);
      default:
        EKAT_KERNEL_ERROR_MSG ("Unsupported cmp operator.\n");
    }
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  int eval (int i, int j) const {
    switch (m_cmp) {
      case Comparison::EQ: return m_left.template eval<T>(i,j) == m_right.template eval<T>(i,j);
      case Comparison::NE: return m_left.template eval<T>(i,j) != m_right.template eval<T>(i,j);
      case Comparison::GT: return m_left.template eval<T>(i,j) >  m_right.template eval<T>(i,j);
      case Comparison::GE: return m_left.template eval<T>(i,j) >= m_right.template eval<T>(i,j);
      case Comparison::LT: return m_left.template eval<T>(i,j) <  m_right.template eval<T>(i,j);
      case Comparison::LE: return m_left.template eval<T>(i,j) <= m_right.template eval<T>(i,j);
      default:
        EKAT_KERNEL_ERROR_MSG ("Unsupported cmp operator.\n");
    }
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  int eval (int i, int j, int k) const {
    switch (m_cmp) {
      case Comparison::EQ: return m_left.template eval<T>(i,j,k) == m_right.template eval<T>(i,j,k);
      case Comparison::NE: return m_left.template eval<T>(i,j,k) != m_right.template eval<T>(i,j,k);
      case Comparison::GT: return m_left.template eval<T>(i,j,k) >  m_right.template eval<T>(i,j,k);
      case Comparison::GE: return m_left.template eval<T>(i,j,k) >= m_right.template eval<T>(i,j,k);
      case Comparison::LT: return m_left.template eval<T>(i,j,k) <  m_right.template eval<T>(i,j,k);
      case Comparison::LE: return m_left.template eval<T>(i,j,k) <= m_right.template eval<T>(i,j,k);
      default:
        EKAT_KERNEL_ERROR_MSG ("Unsupported cmp operator.\n");
    }
  }

protected:

  const ELeft    m_left;
  const ERight   m_right;

  Comparison m_cmp;
};

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator== (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::EQ);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator!= (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::NE);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator> (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::GT);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator>= (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::GE);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator< (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::LT);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator<= (const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::LE);
}

} // namespace scream

#endif // EAMXX_CMP_EXPRESSION_HPP
