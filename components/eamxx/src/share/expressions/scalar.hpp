#ifndef EAMXX_FIELD_EVALUATE_HPP
#define EAMXX_FIELD_EVALUATE_HPP

#include "share/expressions/base_expr.hpp"

namespace scream {

// TODO: support 4+ dim. Also 0d?
template<typename ST>
class ScalarExpression : public Expression<ScalarExpression<ST>,ST> {
public:
  static constexpr bool is_leaf = true;

  using ret_t = ST;

  ScalarExpression (const ST value)
  {
    m_value = value;
  }

  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int) const {
    return m_value;
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int, int) const {
    return m_value;
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator() (int, int, int) const {
    return m_value;
  }

protected:

  ST m_value;
};

} // namespace scream

#endif // EAMXX_FIELD_EVALUATE_HPP
