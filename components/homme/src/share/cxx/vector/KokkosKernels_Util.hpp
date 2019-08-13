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

#ifndef __KOKKOSKERNELS_UTIL_HPP__
#define __KOKKOSKERNELS_UTIL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>
#include <string>

#include <cmath>
#include <ctime>
#include <limits>

namespace KokkosKernels {
namespace Batched {
namespace Experimental {

// view manipulation
template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::Unmanaged |
                                          MemoryTraitsType::RandomAccess |
                                          //  MemoryTraitsType::Atomic |
                                          flag>;

template <typename ViewType>
using UnmanagedViewType = Kokkos::View<
    typename ViewType::data_type, typename ViewType::array_layout,
    typename ViewType::device_type,
    MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged>>;

template <typename ViewType>
using ConstViewType = Kokkos::View<
    typename ViewType::const_data_type, typename ViewType::array_layout,
    typename ViewType::device_type, typename ViewType::memory_traits>;
template <typename ViewType>
using ConstUnmanagedViewType = ConstViewType<UnmanagedViewType<ViewType>>;

template <typename ViewType>
using ScratchViewType = Kokkos::View<
    typename ViewType::data_type, typename ViewType::array_layout,
    typename ViewType::execution_space::scratch_memory_space,
    MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged>>;

// helper for vector type
template <typename T>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<std::is_fundamental<T>::value, size_t>::type
    adjustDimension(const size_t &m) {
  return m;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!std::is_fundamental<T>::value, size_t>::type
    adjustDimension(const size_t &m) {
  return (m / T::vector_length + (m % T::vector_length > 0));
}

template <size_t BufSize, typename SpaceType = Kokkos::DefaultExecutionSpace>
struct Flush {
  typedef double value_type;

  // flush a large host buffer
  Kokkos::View<value_type *, SpaceType> _buf;
  Flush() : _buf("Flush::buf", BufSize / sizeof(double)) {
    Kokkos::deep_copy(_buf, 1);
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &update) { update = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type &update, const volatile value_type &input) {
    update += input;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &update) const { update += _buf[i]; }

  void run() {
    double sum = 0;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<SpaceType>(0, BufSize / sizeof(double)), *this,
        sum);
    SpaceType::fence();
    FILE *fp = fopen("/dev/null", "w");
    fprintf(fp, "%f\n", sum);
    fclose(fp);
  }
};

template <typename ValueType> struct Random;

template <> struct Random<double> {
  Random(const unsigned int seed = 0) { srand(seed); }
  double value() { return rand() / ((double)RAND_MAX + 1.0); }
};

template <typename T, typename SpT = Kokkos::DefaultHostExecutionSpace>
struct ImplicitVector {
  using value_type = T;
  using exec_space = SpT;
};

// Requested compiler vectorization
template <typename T, typename SpT = Kokkos::DefaultHostExecutionSpace>
struct SIMD {
  static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                "KokkosKernels:: Invalid SIMD<> type.");

  using value_type = T;
  using exec_space = SpT;
};

// Intel AVX instruction device (explicit vectorization)
template <typename T, typename SpT = Kokkos::DefaultHostExecutionSpace>
struct AVX {
  static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                "KokkosKernels:: Invalid AVX<> type.");

  using value_type = T;
  using exec_space = SpT;
};

template <class T, int l> struct VectorTag {
  using value_type = typename T::value_type;
  using exec_space = typename T::exec_space;
  using member_type =
      typename Kokkos::Impl::TeamPolicyInternal<exec_space>::member_type;

  static_assert(
      std::is_same<T, SIMD<value_type, exec_space>>::value ||  // host compiler
                                                               // vectorization
          std::is_same<T, AVX<value_type, exec_space>>::value, // || // host AVX
                                                               // vectorization
      // std::is_same<T,SIMT<value_type,exec_space> >::value,   // cuda thread
      // vectorization
      "KokkosKernels:: Invalid VectorUnitTag<> type.");

  using type = VectorTag;
  static const int length = l;
};

// Tags for BLAS
struct Trans {
  struct Transpose {};
  struct NoTranspose {};
  struct ConjTranspose {};
};

struct Side {
  struct Left {};
  struct Right {};
};

struct Uplo {
  struct Upper {};
  struct Lower {};
};

struct Diag {
  struct Unit {
    enum : bool { use_unit_diag = true };
  };
  struct NonUnit {
    enum : bool { use_unit_diag = false };
  };
};

struct Algo {
  struct Level3 {
    struct Unblocked {
      static const char *name() { return "Unblocked"; }
    };
    struct Blocked {
      static const char *name() { return "Blocked"; }
// TODO:: for now harwire the blocksizes; this should reflect
// regieter blocking (not about team parallelism).
// this mb should vary according to
// - team policy (smaller) or range policy (bigger)
// - space (cuda vs host)
// - blocksize input (blk <= 4 mb = 2, otherwise mb = 4), etc.
#if defined(KOKKOS_ENABLE_CUDA)
      template <typename ActiveMemorySpaceType>
      KOKKOS_INLINE_FUNCTION static constexpr typename std::enable_if<
          std::is_same<ActiveMemorySpaceType, Kokkos::CudaSpace>::value,
          int>::type
      mb() {
        return 2;
      }
#endif
      template <typename ActiveMemorySpaceType>
      KOKKOS_INLINE_FUNCTION static constexpr typename std::enable_if<
          std::is_same<ActiveMemorySpaceType, Kokkos::HostSpace>::value,
          int>::type
      mb() {
        return 4;
      }
    };
    struct MKL {
      static const char *name() { return "MKL"; }
    };
    struct CompactMKL {
      static const char *name() { return "CompactMKL"; }
    };
  };

  using Gemm = Level3;
  using Trsm = Level3;
  using LU = Level3;

  struct Level2 {
    struct Unblocked {};
    struct Blocked {
// TODO:: for now harwire the blocksizes; this should reflect
// regieter blocking (not about team parallelism).
// this mb should vary according to
// - team policy (smaller) or range policy (bigger)
// - space (cuda vs host)
// - blocksize input (blk <= 4 mb = 2, otherwise mb = 4), etc.
#if defined(KOKKOS_ENABLE_CUDA)
      template <typename ActiveMemorySpaceType>
      KOKKOS_INLINE_FUNCTION static constexpr typename std::enable_if<
          std::is_same<ActiveMemorySpaceType, Kokkos::CudaSpace>::value,
          int>::type
      mb() {
        return 2;
      }
#endif
      template <typename ActiveMemorySpaceType>
      KOKKOS_INLINE_FUNCTION static constexpr typename std::enable_if<
          std::is_same<ActiveMemorySpaceType, Kokkos::HostSpace>::value,
          int>::type
      mb() {
        return 4;
      }
    };
    struct MKL {};
    struct CompactMKL {};
  };

  using Gemv = Level2;
  using Trsv = Level2;

  //         struct Level1 {
  //           struct Unblocked {};
  //           struct Blocked {
  //             // TODO:: for now harwire the blocksizes; this should reflect
  //             // regieter blocking (not about team parallelism).
  //             // this mb should vary according to
  //             // - team policy (smaller) or range policy (bigger)
  //             // - space (cuda vs host)
  //             // - blocksize input (blk <= 4 mb = 2, otherwise mb = 4), etc.
  // #if defined(KOKKOS_ENABLE_CUDA)
  //             template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION
  //             static constexpr
  //             typename
  //             std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::CudaSpace>::value,int>
  //             ::type mb() { return 4; }
  // #endif
  //             template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION
  //             static constexpr
  //             typename
  //             std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::HostSpace>::value,int>
  //             ::type mb() { return 4; }
  //           };
  //           //struct MKL {};
  //           //struct CompactMKL {};
  //         };
};

struct Util {

  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static void
  packColMajor(ValueType *__restrict__ A, const int m, const int n,
               const ValueType *__restrict__ B, const int bs0, const int bs1) {
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i)
        A[i + j * m] = B[i * bs0 + j * bs1];
  }

  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static void
  packRowMajor(ValueType *__restrict__ A, const int m, const int n,
               const ValueType *__restrict__ B, const int bs0, const int bs1) {
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
        A[i * n + j] = B[i * bs0 + j * bs1];
  }
};
}
}
}
#endif
