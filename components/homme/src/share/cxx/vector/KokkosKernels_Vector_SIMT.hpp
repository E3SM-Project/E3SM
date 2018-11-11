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

#ifndef __KOKKOSKERNELS_VECTOR_SIMT_HPP__
#define __KOKKOSKERNELS_VECTOR_SIMT_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosKernels {
  namespace Batched {
    namespace Experimental {

      /// I am not sure if this is ever useful

      template<typename T, typename SpT, int l>
      class Vector<VectorTag<SIMT<T,SpT>,l> > {
      public:
        using tag_type = VectorTag<SIMT<T,SpT>,l>;

        using type = Vector<tag_type>;
        using value_type = typename tag_type::value_type;
        using member_type = typename tag_type::member_type;

        enum : int { vector_length = tag_type::length };

        typedef value_type data_type[vector_length];
    
      private:
        mutable data_type _data;
    
      public:
        KOKKOS_INLINE_FUNCTION Vector() { 
          Kokkos::parallel_for
            (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length), 
             [&](const int &i) {
              _data[i] = 0; 
            });
        }
        KOKKOS_INLINE_FUNCTION Vector(const value_type val) { 
          Kokkos::parallel_for
            (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length), 
             [&](const int &i) {
              _data[i] = val; 
            });
        }      
        KOKKOS_INLINE_FUNCTION Vector(const type &b) { 
          Kokkos::parallel_for
            (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length), 
             [&](const int &i) {
              _data[i] = b._data[i];
            });
        }      
    
        KOKKOS_INLINE_FUNCTION 
        type& loadAligned(value_type const *p) {
          Kokkos::parallel_for
            (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length), 
             [&](const int &i) {
              _data[i] = p[i];
            });
          return *this;
        }
    
        KOKKOS_INLINE_FUNCTION 
        type& loadUnaligned(value_type const *p) {
          Kokkos::parallel_for
            (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length), 
             [&](const int &i) {
              _data[i] = p[i];
            });
          return *this;
        }

        KOKKOS_INLINE_FUNCTION 
        void storeAligned(value_type *p) const {
          Kokkos::parallel_for
            (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length), 
             [&](const int &i) {
              p[i] = _data[i];
            });
        }
    
        KOKKOS_INLINE_FUNCTION 
        void storeUnaligned(value_type *p) const {
          Kokkos::parallel_for
            (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length), 
             [&](const int &i) {
              p[i] = _data[i];
            });
        }

        KOKKOS_INLINE_FUNCTION 
        value_type& operator[](const int i) const {
          return _data[i];
        }
    
      };

    }
  }
}
#endif
