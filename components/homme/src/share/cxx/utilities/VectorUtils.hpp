/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_VECTOR_UTILS_HPP
#define HOMMEXX_VECTOR_UTILS_HPP

#include "vector/vector_pragmas.hpp"
#include "vector/KokkosKernels_Vector.hpp"
#include "utilities/MathUtils.hpp"

namespace KokkosKernels {
namespace Batched {
namespace Experimental {

template <typename T, typename Space, int l>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<T, Space>, l> >
max (const Vector<VectorTag<SIMD<T, Space>, l> >& a,
     const Vector<VectorTag<SIMD<T, Space>, l> >& b)
{
  Vector<VectorTag<SIMD<T, Space>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<T, Space>, l>>::vector_length; i++) {
    r_val[i] = Homme::max(a[i],b[i]);
  }

  return r_val;
}

template <typename T, typename Space, int l>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<T, Space>, l> >
min (const Vector<VectorTag<SIMD<T, Space>, l> >& a,
     const Vector<VectorTag<SIMD<T, Space>, l> >& b)
{
  Vector<VectorTag<SIMD<T, Space>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<T, Space>, l>>::vector_length; i++) {
    r_val[i] = Homme::min(a[i],b[i]);
  }

  return r_val;
}

template <typename T, typename Space, int l, typename ExpType>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<T, Space>, l> >
pow (const Vector<VectorTag<SIMD<T,Space>,l>>& v, const ExpType p)
{
  using VectorType = Vector<VectorTag<SIMD<T,Space>,l>>;
  VectorType vp;
VECTOR_SIMD_LOOP
  for (int i = 0; i < VectorType::vector_length; ++i) {
    vp[i] = std::pow(v[i],p);
  }

  return vp;
}

template <typename T, typename Space, int l>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<T, Space>, l> >
log (const Vector<VectorTag<SIMD<T,Space>,l>>& v)
{
  using VectorType = Vector<VectorTag<SIMD<T,Space>,l>>;
  VectorType vp;
VECTOR_SIMD_LOOP
  for (int i = 0; i < VectorType::vector_length; ++i) {
    vp[i] = std::log(v[i]);
  }

  return vp;
}

} // namespace KokkosKernels
} // namespace Batched
} // namespace Experimental

#endif // HOMMEXX_VECTOR_UTILS_HPP
