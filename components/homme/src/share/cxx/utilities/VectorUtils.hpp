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

template <typename SpT>
inline
Vector<VectorTag<AVX<double, SpT>, 4> >
max (const Vector<VectorTag<AVX<double, SpT>, 4> >& a,
     const Vector<VectorTag<AVX<double, SpT>, 4> >& b)
{
  return _mm256_max_pd (a, b);
}

template <typename SpT>
inline
Vector<VectorTag<AVX<double, SpT>, 4> >
min (const Vector<VectorTag<AVX<double, SpT>, 4> >& a,
     const Vector<VectorTag<AVX<double, SpT>, 4> >& b)
{
  return _mm256_min_pd (a, b);
}

template <typename SpT>
inline
Vector<VectorTag<AVX<double, SpT>, 8> >
max (const Vector<VectorTag<AVX<double, SpT>, 8> >& a,
     const Vector<VectorTag<AVX<double, SpT>, 8> >& b)
{
  return _mm512_max_pd (a, b);
}

template <typename SpT>
inline
Vector<VectorTag<AVX<double, SpT>, 8> >
min (const Vector<VectorTag<AVX<double, SpT>, 8> >& a,
     const Vector<VectorTag<AVX<double, SpT>, 8> >& b)
{
  return _mm512_min_pd (a, b);
}

template <typename SpT, int l>
KOKKOS_INLINE_FUNCTION
Vector<VectorTag<SIMD<double, SpT>, l> >
max (const Vector<VectorTag<SIMD<double, SpT>, l> >& a,
     const Vector<VectorTag<SIMD<double, SpT>, l> >& b)
{
  Vector<VectorTag<SIMD<double, SpT>, l> > r_val;
VECTOR_IVDEP_LOOP
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
VECTOR_IVDEP_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<double, SpT>, l>>::vector_length; i++) {
    r_val[i] = Homme::min(a[i],b[i]);
  }

  return r_val;
}

// For pow, we use a standard SIMD loop, regardless of what Scalar is.
// This is because pow implementation may be available only in _some_
// avx 512 instruction sets (or have a slightly different name)
template<typename VectorType, typename ExpType>
KOKKOS_INLINE_FUNCTION
VectorType
pow (const VectorType& v, const ExpType p)
{
  VectorType vp;
VECTOR_IVDEP_LOOP
  for (int i = 0; i < VectorType::vector_length; ++i) {
    vp[i] = std::pow(v[i],p);
  }

  return vp;
}

} // namespace KokkosKernels
} // namespace Batched
} // namespace Experimental

#endif // HOMMEXX_VECTOR_UTILS_HPP
