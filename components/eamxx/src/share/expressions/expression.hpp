#ifndef EAMXX_EXPRESSION_HPP
#define EAMXX_EXPRESSION_HPP

#include <Kokkos_Core.hpp>

namespace scream {

template<typename Derived, typename RT>
class Expression {
public:
  static constexpr bool is_leaf = false;
  using ret_t = RT;

  KOKKOS_INLINE_FUNCTION
  ret_t operator()(int i) const {
    cast()(i);
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator()(int i, int j) const {
    cast()(i,j);
  }
  KOKKOS_INLINE_FUNCTION
  ret_t operator()(int i,int j,int k) const {
    cast()(i,j,k);
  }

  KOKKOS_INLINE_FUNCTION
  const Derived& cast () const { return static_cast<const Derived&>(*this); }
};

} // namespace scream

#endif // EAMXX_EXPRESSION_HPP
