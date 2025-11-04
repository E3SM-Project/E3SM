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

  int num_indices () const { return std::max(m_left.num_indices(),m_right.num_indices()); }

  void set_eval_layout (const FieldLayout& fl) {
    m_left.set_eval_layout(fl);
    m_right.set_eval_layout(fl);
  }

  KOKKOS_INLINE_FUNCTION
  void set_eval_data (const EvalData& data) const {
    m_left.set_eval_data(data);
    m_right.set_eval_data(data);
  }


  KOKKOS_INLINE_FUNCTION
  Real eval () const {
    switch (m_cmp) {
      case Comparison::EQ: return m_left.eval() == m_right.eval();
      case Comparison::NE: return m_left.eval() != m_right.eval();
      case Comparison::GT: return m_left.eval() >  m_right.eval();
      case Comparison::GE: return m_left.eval() >= m_right.eval();
      case Comparison::LT: return m_left.eval() <  m_right.eval();
      case Comparison::LE: return m_left.eval() <= m_right.eval();
      default:
        EKAT_KERNEL_ERROR_MSG ("Unsupported cmp operator.\n");
    }
  }
protected:

  ELeft    m_left;
  ERight   m_right;

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
