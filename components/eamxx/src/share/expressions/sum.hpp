#ifndef EAMXX_SUM_EXPRESSION_HPP
#define EAMXX_SUM_EXPRESSION_HPP

#include "share/expressions/base.hpp"

namespace scream {

namespace impl {
template<typename ELeft, typename ERight>
using sum_ret_t = decltype(std::declval<typename ELeft::ret_t>() + 
                           std::declval<typename ERight::ret_t>());
}

template<typename ELeft, typename ERight>
class SumExpression : public Expression<SumExpression<ELeft,ERight>,impl::sum_ret_t<ELeft,ERight>>{
public:
  using ret_left  = typename ELeft::ret_t;
  using ret_right = typename ERight::ret_t;
  using ret_t     = impl::sum_ret_t<ELeft,ERight>;
  static constexpr bool is_leaf = false;

  SumExpression (const ELeft& left, const ERight& right)
    : m_left(left)
    , m_right(right)
  {
    // Nothing to do here
  }

  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int i) const {
    return static_cast<ret_t>(m_left(i) + m_right(i));
  }
  KOKKOS_INLINE_FUNCTION
  ret_t evaluate (int i, int j) const {
    return static_cast<ret_t>(m_left(i,j) + m_right(i,j));
  }
  KOKKOS_INLINE_FUNCTION
  ret_t evaluate (int i, int j, int k) const {
    return static_cast<ret_t>(m_left(i,j,k) + m_right(i,j,k));
  }

protected:

  const ELeft    m_left;
  const ERight   m_right;
};

template<typename ELeft, typename RetL, typename ERight, typename RetR>
SumExpression<ELeft,ERight>
operator+ (const Expression<ELeft,RetL>& l, const Expression<ERight,RetR>& r)
{
  return SumExpression<ELeft,ERight>(l.cast(),r.cast());
}

} // namespace scream

#endif // EAMXX_SUM_EXPRESSION_HPP
