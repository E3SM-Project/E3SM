#ifndef EAMXX_EXPRESSION_HPP
#define EAMXX_EXPRESSION_HPP

#include <Kokkos_Core.hpp>

#include <share/core/eamxx_types.hpp>

namespace scream {

template<typename Derived>
class Expression {
public:
  static constexpr bool is_leaf = false;

  int num_indices () const { return cast().num_indices(); }

  KOKKOS_INLINE_FUNCTION
  Real eval(int i) const {
    return cast().eval(i);
  }
  KOKKOS_INLINE_FUNCTION
  Real eval(int i, int j) const {
    return cast().eval(i,j);
  }
  KOKKOS_INLINE_FUNCTION
  Real eval(int i,int j,int k) const {
    return cast().eval(i,j,k);
  }

  KOKKOS_INLINE_FUNCTION
  const Derived& cast () const { return static_cast<const Derived&>(*this); }

  // void set_eval_layout (const FieldLayout& fl) { cast().set_eval_layout(fl); }
};

} // namespace scream

#endif // EAMXX_EXPRESSION_HPP
