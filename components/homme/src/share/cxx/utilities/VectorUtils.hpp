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

template <typename SpT, int l>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<double, SpT>, l> >
max (const Vector<VectorTag<SIMD<double, SpT>, l> >& a,
     const Vector<VectorTag<SIMD<double, SpT>, l> >& b)
{
  Vector<VectorTag<SIMD<double, SpT>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<double, SpT>, l>>::vector_length; i++) {
    r_val[i] = Homme::max(a[i],b[i]);
  }

  return r_val;
}

template <typename SpT, int l>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<double, SpT>, l> >
min (const Vector<VectorTag<SIMD<double, SpT>, l> >& a,
     const Vector<VectorTag<SIMD<double, SpT>, l> >& b)
{
  Vector<VectorTag<SIMD<double, SpT>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<double, SpT>, l>>::vector_length; i++) {
    r_val[i] = Homme::min(a[i],b[i]);
  }

  return r_val;
}

template <typename SpT, int l, typename ExpType>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<double, SpT>, l> >
pow (const Vector<VectorTag<SIMD<double,SpT>,l>>& v, const ExpType p)
{
  using VectorType = Vector<VectorTag<SIMD<double,SpT>,l>>;
  VectorType vp;
VECTOR_SIMD_LOOP
  for (int i = 0; i < VectorType::vector_length; ++i) {
    vp[i] = std::pow(v[i],p);
  }

  return vp;
}

} // namespace KokkosKernels
} // namespace Batched
} // namespace Experimental

#endif // HOMMEXX_VECTOR_UTILS_HPP
