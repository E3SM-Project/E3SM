#ifndef INCLUDE_SCREAM_PACK_KOKKOS
#define INCLUDE_SCREAM_PACK_KOKKOS

#include "scream_pack.hpp"
#include "scream_kokkos_meta.hpp"

namespace scream {
namespace pack {

/* These functions combine Pack, Mask, and Kokkos::Views.
 */

// Index a scalar array with Pack indices, returning a compatible Pack of array
// values.
template<typename Array1, typename IdxPack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<IdxPack, Pack<typename Array1::value_type, IdxPack::n> >
index (const Array1& a, const IdxPack& i0,
       typename std::enable_if<Array1::Rank == 1>::type* = nullptr) {
  Pack<typename Array1::non_const_value_type, IdxPack::n> p;
  vector_simd for (int i = 0; i < IdxPack::n; ++i)
    p[i] = a(i0[i]);
  return p;
}
template<typename Array2, typename IdxPack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<IdxPack, Pack<typename Array2::value_type, IdxPack::n> >
index (const Array2& a, const IdxPack& i0, const IdxPack& i1,
       typename std::enable_if<Array2::Rank == 2>::type* = nullptr) {
  Pack<typename Array2::non_const_value_type, IdxPack::n> p;
  vector_simd for (int i = 0; i < IdxPack::n; ++i)
    p[i] = a(i0[i], i1[i]);
  return p;
}

template <typename T, typename ...Parms, int pack_size> KOKKOS_FORCEINLINE_FUNCTION
typename ko::Unmanaged<Kokkos::View<T**, Parms...> >::type
scalarize (const Kokkos::View<scream::pack::Pack<T, pack_size>**, Parms...>& vp) {
return typename ko::Unmanaged<Kokkos::View<T**, Parms...> >::type(
    reinterpret_cast<T*>(vp.data()), vp.extent_int(0), pack_size * vp.extent_int(1));
}

} // namespace pack
} // namespace scream

#endif // INCLUDE_SCREAM_PACK_KOKKOS
