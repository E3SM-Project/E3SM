#ifndef EAMXX_SCALAR_EXPRESSION_HPP
#define EAMXX_SCALAR_EXPRESSION_HPP

#include "share/expressions/base.hpp"

namespace scream {

// TODO: support 4+ dim. Also 0d?
template<typename ST>
class ScalarExpression : public Expression<ScalarExpression<ST>> {
public:
  static constexpr bool is_leaf = true;

  using ret_t = ST;

  ScalarExpression (const ST value)
  {
    m_value = value;
  }

  int num_indices () override { return 0; }

  KOKKOS_INLINE_FUNCTION
  Real eval(int) const {
    return m_value;
  }
  KOKKOS_INLINE_FUNCTION
  Real eval(int,int) const {
    return m_value;
  }
  KOKKOS_INLINE_FUNCTION
  Real eval(int,int,int) const {
    return m_value;
  }

  // void set_eval_layout (const FieldLayout& fl) {}

protected:

  ST m_value;
};

} // namespace scream

#endif // EAMXX_SCALAR_EXPRESSION_HPP
