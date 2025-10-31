#ifndef EAMXX_CMP_EXPRESSION_HPP
#define EAMXX_CMP_EXPRESSION_HPP

#include "share/field_math/expression.hpp"

#include "share/core/eamxx_types.hpp" // For Comparison

#include <ekat_kernel_assert.hpp>

namespace scream {

template<typename ELeft, typename ERight>
class CmpExpression : public Expression<CmpExpression<ELeft,ERight>,int> {
public:
  static_assert (std::is_constructible<typename ELeft::ret_t,typename ERight::ret_t>::value,
      "Incompatible left and right expressions");
  static constexpr bool is_leaf = false;
  using ret_t = int;

  CmpExpression (const ELeft& left, const ERight& right, Comparison CMP)
    : m_left(left)
    , m_right(right)
    , m_cmp(CMP)
  {
    // Nothing to do here
  }

  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i) const {
    switch (m_cmp) {
      case Comparison::EQ: return m_left(i) == m_right(i);
      case Comparison::NE: return m_left(i) != m_right(i);
      case Comparison::GT: return m_left(i) >  m_right(i);
      case Comparison::GE: return m_left(i) >= m_right(i);
      case Comparison::LT: return m_left(i) <  m_right(i);
      case Comparison::LE: return m_left(i) <= m_right(i);
      default:
        EKAT_KERNEL_ERROR_MSG ("Unsupported cmp operator.\n");
    }
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i, int j) const {
    switch (m_cmp) {
      case Comparison::EQ: return m_left(i,j) == m_right(i,j);
      case Comparison::NE: return m_left(i,j) != m_right(i,j);
      case Comparison::GT: return m_left(i,j) >  m_right(i,j);
      case Comparison::GE: return m_left(i,j) >= m_right(i,j);
      case Comparison::LT: return m_left(i,j) <  m_right(i,j);
      case Comparison::LE: return m_left(i,j) <= m_right(i,j);
      default:
        EKAT_KERNEL_ERROR_MSG ("Unsupported cmp operator.\n");
    }
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i, int j, int k) const {
    switch (m_cmp) {
      case Comparison::EQ: return m_left(i,j,k) == m_right(i,j,k);
      case Comparison::NE: return m_left(i,j,k) != m_right(i,j,k);
      case Comparison::GT: return m_left(i,j,k) >  m_right(i,j,k);
      case Comparison::GE: return m_left(i,j,k) >= m_right(i,j,k);
      case Comparison::LT: return m_left(i,j,k) <  m_right(i,j,k);
      case Comparison::LE: return m_left(i,j,k) <= m_right(i,j,k);
      default:
        EKAT_KERNEL_ERROR_MSG ("Unsupported cmp operator.\n");
    }
  }

protected:

  std::conditional_t<ELeft::is_leaf, const ELeft&, const ELeft>   m_left;
  std::conditional_t<ERight::is_leaf,const ERight&,const ERight>  m_right;

  Comparison m_cmp;
};

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator== (const Expression<ELeft, typename ELeft::ret_t>& l, const Expression<ERight, typename ELeft::ret_t>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::EQ);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator!= (const Expression<ELeft, typename ELeft::ret_t>& l, const Expression<ERight, typename ELeft::ret_t>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::NE);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator> (const Expression<ELeft, typename ELeft::ret_t>& l, const Expression<ERight, typename ELeft::ret_t>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::GT);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator>= (const Expression<ELeft, typename ELeft::ret_t>& l, const Expression<ERight, typename ELeft::ret_t>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::GE);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator< (const Expression<ELeft, typename ELeft::ret_t>& l, const Expression<ERight, typename ELeft::ret_t>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::LT);
}

template<typename ELeft, typename ERight>
CmpExpression<ELeft,ERight> operator<= (const Expression<ELeft, typename ELeft::ret_t>& l, const Expression<ERight, typename ELeft::ret_t>& r)
{
  return CmpExpression<ELeft,ERight>(l.cast(),r.cast(),Comparison::LE);
}

} // namespace scream

#endif // EAMXX_CMP_EXPRESSION_HPP
