#ifndef EAMXX_CONDITIONAL_EXPRESSION_HPP
#define EAMXX_CONDITIONAL_EXPRESSION_HPP

#include "share/expressions/base.hpp"

#include "share/core/eamxx_types.hpp" // For Comparison

#include <ekat_kernel_assert.hpp>

namespace scream {

template<typename ECond, typename ELeft, typename ERight>
class ConditionalExpression : public Expression<ConditionalExpression<ECond,ELeft,ERight>> {
public:
  static constexpr bool is_leaf = false;

  ConditionalExpression (const ECond& cmp, const ELeft& left, const ERight& right)
    : m_cmp(cmp)
    , m_left(left)
    , m_right(right)
  {
    // Nothing to do here
  }

  int num_indices () const {
    return std::max(m_cmp.num_indices(),std::max(m_left.num_indices(),m_right.num_indices()));
  }

  void set_eval_layout (const FieldLayout& fl) {
    m_cmp.set_eval_layout(fl);
    m_left.set_eval_layout(fl);
    m_right.set_eval_layout(fl);
  }

  KOKKOS_INLINE_FUNCTION
  void set_eval_data (const EvalData& data) const {
    m_cmp.set_eval_data(data);
    m_left.set_eval_data(data);
    m_right.set_eval_data(data);
  }


  KOKKOS_INLINE_FUNCTION
  Real eval () const {
    if (m_cmp.eval())
      return m_left.eval();
    else
      return m_right.eval();
  }
protected:

  ECond    m_cmp;
  ELeft    m_left;
  ERight   m_right;
};

template<typename ECond, typename ELeft, typename ERight>
ConditionalExpression<ECond,ELeft,ERight>
conditional(const Expression<ECond>& c, const Expression<ELeft>& l, const Expression<ERight>& r)
{
  return ConditionalExpression<ECond,ELeft,ERight>(c.cast(),l.cast(),r.cast());
}

} // namespace scream

#endif // EAMXX_CONDITIONAL_EXPRESSION_HPP
