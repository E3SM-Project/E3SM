#ifndef SCREAM_TUPLE_HPP
#define SCREAM_TUPLE_HPP

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <utility>

namespace scream
{

namespace details {

// ===========================================================================
// A minimal tuple that can be stored in, and used from, device functors
// ===========================================================================
//
// std::tuple (and std::apply) are not guaranteed to be usable in device code:
// they are not marked as KOKKOS_INLINE_FUNCTION, and rely on compiler flags
// (like nvcc's --expt-relaxed-constexpr) to work in kernels. This is a small
// replacement, which provides only what we need:
//
//   auto t = make_tuple(a,b,c);      // Tuple<A,B,C>, with copies of a,b,c
//   get<1>(t);                       // b
//   apply(f,t,x,y);                  // f(a,b,c,x,y)
//
// NOTE: apply is NOT std::apply: it allows (and forwards) extra arguments,
// that are appended after the tuple entries.
//
// These live in scream::details (not scream), so that our make_tuple/get/apply don't
// compete with the std ones when called unqualified from inside namespace scream.
//
// The tuple is a "head and tail" chain: Tuple<A,B,C> holds an A and a Tuple<B,C>
// (and Tuple<C> holds just a C, so that the tuple has no size overhead).
// A tuple must have at least one entry.
// ===========================================================================

// NOTE: there is no empty tuple. Tuple<> is ill-formed, on purpose: we would not be able to store it
//       without paying for it (an empty member still has sizeof>=1), and we never need it.
template<typename T, typename... Ts>
struct Tuple {
  T head;
  Tuple<Ts...> tail;
};

// The last entry has no tail (so that the tuple is as big as its entries)
template<typename T>
struct Tuple<T> {
  T head;
};

template<typename T, typename... Ts>
KOKKOS_INLINE_FUNCTION
Tuple<T,Ts...> make_tuple (const T& t, const Ts&... ts)
{
  if constexpr (sizeof...(Ts)==0) {
    return {t};
  } else {
    return {t, details::make_tuple(ts...)};
  }
}

// Access the I-th entry (walks down the tail I times)
template<std::size_t I, typename T, typename... Ts>
KOKKOS_INLINE_FUNCTION
const auto& get (const Tuple<T,Ts...>& t)
{
  if constexpr (I==0) {
    return t.head;
  } else {
    return details::get<I-1>(t.tail);
  }
}

namespace impl {

template<typename F, typename Tup, std::size_t... Is, typename... Extra>
KOKKOS_INLINE_FUNCTION
decltype(auto) apply_impl (F&& f, const Tup& t, std::index_sequence<Is...>, Extra&&... extra)
{
  // NOTE: static_cast instead of std::forward, to stay clear of std in device code
  return f(details::get<Is>(t)..., static_cast<Extra&&>(extra)...);
}

} // namespace impl

// Call f(t_0,...,t_{N-1},extra...), where t_i are the entries of t.
// The extra arguments are perfectly forwarded (so they can be non-const references,
// e.g., the result argument of a parallel_reduce kernel).
template<typename F, typename... Ts, typename... Extra>
KOKKOS_INLINE_FUNCTION
decltype(auto) apply (F&& f, const Tuple<Ts...>& t, Extra&&... extra)
{
  return impl::apply_impl(static_cast<F&&>(f), t,
                          std::index_sequence_for<Ts...>(),
                          static_cast<Extra&&>(extra)...);
}

} // namespace details

} // namespace scream

#endif // SCREAM_TUPLE_HPP
