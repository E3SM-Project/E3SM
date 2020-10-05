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

// NOTE (LB, Jan 2020): Stripped of most content, for the needs of HOMME

namespace KokkosKernels {
namespace Batched {
namespace Experimental {

// Requested compiler vectorization
template <typename T, typename SpT = Kokkos::DefaultHostExecutionSpace>
struct SIMD {
  static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                "KokkosKernels:: Invalid SIMD<> type.");

  using value_type = T;
  using exec_space = SpT;
};

template <class T, int l> struct VectorTag {
  using value_type = typename T::value_type;
  using exec_space = typename T::exec_space;
  using member_type =
      typename Kokkos::Impl::TeamPolicyInternal<exec_space>::member_type;

  static_assert(
      std::is_same<T, SIMD<value_type, exec_space>>::value, // host compiler vectorization
      "KokkosKernels:: Invalid VectorUnitTag<> type.");

  using type = VectorTag;
  static const int length = l;
};

} // namespace KokkosKernels
} // namespace Batched
} // namespace Experimental
#endif
