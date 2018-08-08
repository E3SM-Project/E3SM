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

#ifndef __KOKKOSKERNELS_VECTOR_AVX512D_HPP__
#define __KOKKOSKERNELS_VECTOR_AVX512D_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)
#if defined(__AVX512F__)

#include <immintrin.h>

namespace KokkosKernels {
namespace Batched {
namespace Experimental {
///
/// AVX512D double
///

template <typename SpT> class Vector<VectorTag<AVX<double, SpT>, 8> > {
public:
  using type = Vector<VectorTag<AVX<double, SpT>, 8> >;
  using value_type = double;
  using real_type = double;

  static const int vector_length = 8;

  union data_type {
    __m512d v;
    double d[8];
  };

  KOKKOS_INLINE_FUNCTION
  static const char *label() { return "AVX512"; }

private:
  mutable data_type _data;

public:
  inline Vector() { _data.v = _mm512_setzero_pd(); }
  inline Vector(const value_type val) { _data.v = _mm512_set1_pd(val); }
  inline Vector(const type &b) { _data.v = b._data.v; }
  inline Vector(__m512d const &val) { _data.v = val; }

  inline type &operator=(__m512d const &val) {
    _data.v = val;
    return *this;
  }

  inline operator __m512d() const { return _data.v; }

  inline type &loadAligned(value_type const *p) {
    _data.v = _mm512_load_pd(p);
    return *this;
  }

  inline type &loadUnaligned(value_type const *p) {
    _data.v = _mm512_loadu_pd(p);
    return *this;
  }

  inline void storeAligned(value_type *p) const { _mm512_store_pd(p, _data.v); }

  inline void storeUnaligned(value_type *p) const {
    _mm512_storeu_pd(p, _data.v);
  }

#ifdef NDEBUG
  // Does nothing in non-debug mode
  KOKKOS_INLINE_FUNCTION
  void debug_set_invalid(int left, int right) {}
#else
  // left, right specify the closed range of indices to set to quiet NaNs
  KOKKOS_INLINE_FUNCTION
  void debug_set_invalid(int left, int right) {
    for(int i = left; i <= right; i++) {
      _data.d[i] = 0.0 / 0.0;
    }
  }
#endif // NDEBUG

  KOKKOS_INLINE_FUNCTION
  void shift_left(int num_values) {
    assert(num_values > 0);
    for(int i = 0; i < vector_length - num_values; i++) {
      _data.d[i] = _data.d[i + num_values];
    }
    debug_set_invalid(vector_length - num_values, vector_length - 1);
  }

  KOKKOS_INLINE_FUNCTION
  void shift_right(int num_values) {
    assert(num_values > 0);
    for(int i = vector_length - 1; i >= num_values; i--) {
      _data.d[i] = _data.d[i - num_values];
    }
    debug_set_invalid(0, num_values - 1);
  }

  inline value_type &operator[](int i) const { return _data.d[i]; }
};

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator+(Vector<VectorTag<AVX<double, SpT>, 8> > const &a,
          Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return _mm512_add_pd(a, b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator+(Vector<VectorTag<AVX<double, SpT>, 8> > const &a, const double b) {
  return a + Vector<VectorTag<AVX<double, SpT>, 8> >(b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator+(const double a, Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return Vector<VectorTag<AVX<double, SpT>, 8> >(a) + b;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator+=(Vector<VectorTag<AVX<double, SpT>, 8> > &a,
           Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  a = a + b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator+=(Vector<VectorTag<AVX<double, SpT>, 8> > &a, const double b) {
  a = a + b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator++(Vector<VectorTag<AVX<double, SpT>, 8> > &a, int) {
  Vector<VectorTag<AVX<double, SpT>, 8> > a0 = a;
  a = a + 1.0;
  return a0;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator++(Vector<VectorTag<AVX<double, SpT>, 8> > &a) {
  a = a + 1.0;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator-(Vector<VectorTag<AVX<double, SpT>, 8> > const &a,
          Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return _mm512_sub_pd(a, b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator-(Vector<VectorTag<AVX<double, SpT>, 8> > const &a, const double b) {
  return a - Vector<VectorTag<AVX<double, SpT>, 8> >(b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator-(const double a, Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return Vector<VectorTag<AVX<double, SpT>, 8> >(a) - b;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator-=(Vector<VectorTag<AVX<double, SpT>, 8> > &a,
           Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  a = a - b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator-=(Vector<VectorTag<AVX<double, SpT>, 8> > &a, const double b) {
  a = a - b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator--(Vector<VectorTag<AVX<double, SpT>, 8> > &a, int) {
  Vector<VectorTag<AVX<double, SpT>, 8> > a0 = a;
  a = a - 1.0;
  return a0;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator--(Vector<VectorTag<AVX<double, SpT>, 8> > &a) {
  a = a - 1.0;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator*(Vector<VectorTag<AVX<double, SpT>, 8> > const &a,
          Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return _mm512_mul_pd(a, b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator*(Vector<VectorTag<AVX<double, SpT>, 8> > const &a, const double b) {
  return a * Vector<VectorTag<AVX<double, SpT>, 8> >(b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator*(const double a, Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return Vector<VectorTag<AVX<double, SpT>, 8> >(a) * b;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator*=(Vector<VectorTag<AVX<double, SpT>, 8> > &a,
           Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  a = a * b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator*=(Vector<VectorTag<AVX<double, SpT>, 8> > &a, const double b) {
  a = a * b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator/(Vector<VectorTag<AVX<double, SpT>, 8> > const &a,
          Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return _mm512_div_pd(a, b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator/(Vector<VectorTag<AVX<double, SpT>, 8> > const &a, const double b) {
  return a / Vector<VectorTag<AVX<double, SpT>, 8> >(b);
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator/(const double a, Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  return Vector<VectorTag<AVX<double, SpT>, 8> >(a) / b;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator/=(Vector<VectorTag<AVX<double, SpT>, 8> > &a,
           Vector<VectorTag<AVX<double, SpT>, 8> > const &b) {
  a = a / b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> > &
operator/=(Vector<VectorTag<AVX<double, SpT>, 8> > &a, const double b) {
  a = a / b;
  return a;
}

template <typename SpT>
inline static Vector<VectorTag<AVX<double, SpT>, 8> >
operator-(Vector<VectorTag<AVX<double, SpT>, 8> > const &a) {
  return -1 * a;
}

} // Experimental
} // Batched
} // KokkosKernels

#endif
#endif
