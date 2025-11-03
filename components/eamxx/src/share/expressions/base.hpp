#ifndef EAMXX_EXPRESSION_HPP
#define EAMXX_EXPRESSION_HPP

#include <Kokkos_Core.hpp>

namespace scream {

template<typename Derived>
class Expression {
public:
  static constexpr bool is_leaf = false;

  int num_indices () const { return cast().num_indices(); }

  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval(int i) const {
    cast()(i);
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval(int i, int j) const {
    cast()(i,j);
  }
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  T eval(int i,int j,int k) const {
    cast()(i,j,k);
  }

  KOKKOS_INLINE_FUNCTION
  const Derived& cast () const { return static_cast<const Derived&>(*this); }

  // void set_eval_layout (const FieldLayout& fl) { cast().set_eval_layout(fl); }
};

} // namespace scream

#endif // EAMXX_EXPRESSION_HPP
