#ifndef SCREAM_BOUND_FUNCTOR_HPP
#define SCREAM_BOUND_FUNCTOR_HPP

#include "share/util/eamxx_tuple.hpp"

#include <type_traits>

namespace scream
{

// ===========================================================================
// A functor that prepends some stored arguments to the ones it is called with
// ===========================================================================
//
// This is a device-friendly version of std::bind_front (std::bind_front is not
// guaranteed to be usable in device code). Like std::bind_front, and unlike std::bind,
// there are no placeholders: the stored arguments are always the first ones.
//
//   auto g = bind_functor(f,a,b);
//   g(x,y);                              // calls f(a,b,x,y)
//
// This is useful for Kokkos kernels. Kokkos calls the functor of a parallel_for or
// parallel_reduce with ONLY the indices (and, for reductions, the result argument):
//
//   parallel_for:     functor(i,j,...)
//   parallel_reduce:  functor(i,j,...,result)
//
// But sometimes we want to write a kernel that also needs some other arguments, whose
// type we only know in a templated context (e.g., views of runtime-dependent type):
//
//   func(a,b,i,j,...)
//
// A BoundFunctor bridges the two: it stores func and (a,b), and its operator() takes
// whatever Kokkos gives it, and calls func with the stored arguments first:
//
//   auto g = bind_functor(func,a,b);
//   g(i,j,...);                          // calls func(a,b,i,j,...)
//   g(i,j,...,result);                   // calls func(a,b,i,j,...,result)
//
// The call arguments are perfectly forwarded, so func can modify the reduction result.
// Since the bound functor is copied to device, func and the stored arguments must be
// device-copyable (e.g., func is a KOKKOS_LAMBDA or a functor struct, and the stored arguments
// are Kokkos views). Note that, if func needs to be generic (e.g., because the type of the stored
// arguments is not known where func is written), it cannot be a KOKKOS_LAMBDA, since nvcc does not
// allow generic extended lambdas: use a functor struct with a templated operator() instead.
// ===========================================================================

template<typename Func, typename... StoredArgs>
struct BoundFunctor {
  static_assert (sizeof...(StoredArgs)>0,
      "[BoundFunctor] Error! At least one argument must be bound (otherwise, just use the functor).\n");

  Func func;
  details::Tuple<StoredArgs...> stored;

  template<typename... Args>
  KOKKOS_INLINE_FUNCTION
  decltype(auto) operator() (Args&&... args) const {
    return details::apply(func,stored,static_cast<Args&&>(args)...);
  }
};

template<typename Func, typename... StoredArgs>
BoundFunctor<std::decay_t<Func>,std::decay_t<StoredArgs>...>
bind_functor (const Func& func, const StoredArgs&... stored)
{
  return {func, details::make_tuple(stored...)};
}

} // namespace scream

#endif // SCREAM_BOUND_FUNCTOR_HPP
