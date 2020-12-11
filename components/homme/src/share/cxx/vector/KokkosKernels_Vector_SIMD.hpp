//@HEADER
// ************************************************************************
//
//               KokkosKernels: Linear Algebra and Graph Kernels
//                 Copyright 2016 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef __KOKKOSKERNELS_VECTOR_SIMD_HPP__
#define __KOKKOSKERNELS_VECTOR_SIMD_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <assert.h>
#include "vector_pragmas.hpp"

namespace KokkosKernels {
namespace Batched {
namespace Experimental {

template <typename T, typename SpT, int l>
class Vector<VectorTag<SIMD<T, SpT>, l> > {
public:
  using tag_type = VectorTag<SIMD<T, SpT>, l>;

  using type = Vector<tag_type>;
  using value_type = typename tag_type::value_type;
  using member_type = typename tag_type::member_type;
  using data_type = value_type[l];

  static const int vector_length = tag_type::length;

  KOKKOS_INLINE_FUNCTION
  static const char *label() { return "SIMD"; }

private:
  mutable data_type _data;

public:
  KOKKOS_FORCEINLINE_FUNCTION Vector() {
VECTOR_SIMD_LOOP
    for (int i = 0; i < vector_length; i++) {
      _data[i] = 0;
    }
  }
  template <typename ArgValueType>
  KOKKOS_FORCEINLINE_FUNCTION Vector(const ArgValueType val) {
VECTOR_SIMD_LOOP
    for (int i = 0; i < vector_length; i++) {
      _data[i] = val;
    }
  }
  KOKKOS_FORCEINLINE_FUNCTION Vector(const type &b) {
VECTOR_SIMD_LOOP
    for (int i = 0; i < vector_length; i++) {
      _data[i] = b._data[i];
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  type &loadAligned(value_type const *p) {
VECTOR_SIMD_LOOP
    for (int i = 0; i < vector_length; i++) {
      _data[i] = p[i];
    }
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  type &loadUnaligned(value_type const *p) {
VECTOR_SIMD_LOOP
    for (int i = 0; i < vector_length; i++) {
      _data[i] = p[i];
    }
    return *this;
  }

  // AVX has aligned version and unaligned version;
  // aligned load store are recommended if memory is aligned
  // in this version, there is no difference between aligned and unaligned

  KOKKOS_FORCEINLINE_FUNCTION
  void storeAligned(value_type *p) const {
VECTOR_SIMD_LOOP
    for (int i = 0; i < vector_length; i++) {
      p[i] = _data[i];
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void storeUnaligned(value_type *p) const {
VECTOR_SIMD_LOOP
    for (int i = 0; i < vector_length; i++) {
      p[i] = _data[i];
    }
  }

  KOKKOS_FORCEINLINE_FUNCTION
  value_type &operator[](const int i) const { return _data[i]; }

#ifdef NDEBUG
  // left, right specify the closed range of indices to set to quiet NaNs
  KOKKOS_FORCEINLINE_FUNCTION
  void debug_set_invalid(int left, int right) {}
#else
  KOKKOS_FORCEINLINE_FUNCTION
  void debug_set_invalid(int left, int right) {
    for(int i = left; i <= right; i++) {
      _data[i] = 0.0 / 0.0;
    }
  }
#endif

  KOKKOS_FORCEINLINE_FUNCTION
  void shift_left(int num_values) {
    assert(num_values > 0);
    for(int i = 0; i < vector_length - num_values; i++) {
      _data[i] = _data[i + num_values];
    }
    debug_set_invalid(vector_length - num_values, vector_length - 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void shift_right(int num_values) {
    assert(num_values > 0);
    for(int i = vector_length - 1; i >= num_values; i--) {
      _data[i] = _data[i - num_values];
    }
    debug_set_invalid(0, num_values - 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  value_type reduce_add() {
    value_type sum = 0.0;
    for (int i=0; i<vector_length; ++i) {
      sum += _data[i];
    }
    return sum;
  }
};

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator+(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  Vector<VectorTag<SIMD<T, SpT>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<T, SpT>, l> >::vector_length; i++) {
    r_val[i] = a[i] + b[i];
  }
  return r_val;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator+(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  return a + Vector<VectorTag<SIMD<T, SpT>, l> >(b);
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator+(const typename VectorTag<SIMD<T, SpT>, l>::value_type a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  return Vector<VectorTag<SIMD<T, SpT>, l> >(a) + b;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator+=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  a = a + b;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator+=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  a = a + b;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator++(Vector<VectorTag<SIMD<T, SpT>, l> > &a, int) {
  Vector<VectorTag<SIMD<T, SpT>, l> > a0 = a;
  a = a + 1.0;
  return a0;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator++(Vector<VectorTag<SIMD<T, SpT>, l> > &a) {
  a = a + 1.0;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator-(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  Vector<VectorTag<SIMD<T, SpT>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<T, SpT>, l> >::vector_length; i++) {
    r_val[i] = a[i] - b[i];
  }
  return r_val;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator-(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  return a - Vector<VectorTag<SIMD<T, SpT>, l> >(b);
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator-(const typename VectorTag<SIMD<T, SpT>, l>::value_type a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  return Vector<VectorTag<SIMD<T, SpT>, l> >(a) - b;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator-(Vector<VectorTag<SIMD<T, SpT>, l> > a) {
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<T, SpT>, l> >::vector_length; i++) {
    a[i] = -a[i];
  }
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator-=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  a = a - b;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator-=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  a = a - b;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator--(Vector<VectorTag<SIMD<T, SpT>, l> > &a, int) {
  Vector<VectorTag<SIMD<T, SpT>, l> > a0 = a;
  a = a - 1.0;
  return a0;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator--(Vector<VectorTag<SIMD<T, SpT>, l> > &a) {
  a = a - 1.0;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator*(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  Vector<VectorTag<SIMD<T, SpT>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<T, SpT>, l> >::vector_length; i++) {
    r_val[i] = a[i] * b[i];
  }
  return r_val;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator*(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  return a * Vector<VectorTag<SIMD<T, SpT>, l> >(b);
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator*(const typename VectorTag<SIMD<T, SpT>, l>::value_type a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  return Vector<VectorTag<SIMD<T, SpT>, l> >(a) * b;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator*=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  a = a * b;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator*=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  a = a * b;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator/(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  Vector<VectorTag<SIMD<T, SpT>, l> > r_val;
VECTOR_SIMD_LOOP
  for (int i = 0; i < Vector<VectorTag<SIMD<T, SpT>, l> >::vector_length; i++) {
    r_val[i] = a[i] / b[i];
  }
  return r_val;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator/(Vector<VectorTag<SIMD<T, SpT>, l> > const &a,
          const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  return a / Vector<VectorTag<SIMD<T, SpT>, l> >(b);
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> >
operator/(const typename VectorTag<SIMD<T, SpT>, l>::value_type a,
          Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  return Vector<VectorTag<SIMD<T, SpT>, l> >(a) / b;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator/=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           Vector<VectorTag<SIMD<T, SpT>, l> > const &b) {
  a = a / b;
  return a;
}

template <typename T, typename SpT, int l>
KOKKOS_FORCEINLINE_FUNCTION
static Vector<VectorTag<SIMD<T, SpT>, l> > &
operator/=(Vector<VectorTag<SIMD<T, SpT>, l> > &a,
           const typename VectorTag<SIMD<T, SpT>, l>::value_type b) {
  a = a / b;
  return a;
}

} // Experimental
} // Batched
} // KokkosKernels

#endif
