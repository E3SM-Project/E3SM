
#pragma once

#include <type_traits>

#include "samxx_const.h"
#include "vars.h"
// #include "ekat/kokkos/ekat_kokkos_utils.hpp"
// #include "ekat/ekat_pack_kokkos.hpp"
// #include "ekat/ekat_assert.hpp"
// #include "scream_session.hpp"

// for kokkos debuge only
#if defined(DEBUG)
#include "Cuda/Kokkos_Cuda_Instance.hpp"
#endif

extern "C" {
 void scream_session_init();
 void scream_session_finalize();
}

YAKL_INLINE int is_same_str(const char *str_a, const char *str_b)
{
  int match = 0;
  unsigned i = 0;
  unsigned done = 0;
  while ((match == 0) && !done)
  {
    if ((str_a[i] == 0) || (str_b[i] == 0)) {
      done = 1;
    } else if (str_a[i] != str_b[i]) {
      match = i+1;
      if (((int)str_a[i] - (int)str_b[i]) < 0) match = 0 - (i + 1);
    }
    i++;
  }
  return match;
}

// // convert YAKL multidimensional array to Kokkos view on GPU
// template <typename SizeT, typename ViewT>
// void array_to_view(const typename ViewT::value_type::scalar* const data,
//                    const SizeT& size,
//                    ViewT& view)
// {
//   using PackT = typename ViewT::value_type;
//   EKAT_ASSERT(PackT::n >= 1);

//   const size_t npack = (size + PackT::n-1)/PackT::n;

// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {npack, PackT::n}), KOKKOS_LAMBDA(int k, int s) {
//     const size_t scalar_offset = k*PackT::n;
//     if (scalar_offset+s < size) view(k)[s] = data[scalar_offset+s];
//   });
// }

// // 2D YAKL array to Kokkos view
// template <typename SizeT, typename ViewT>
// void array_to_view(const typename ViewT::value_type::scalar* const data,
//                    const SizeT& dim1_size,
//                    const SizeT& dim2_size,
//                    ViewT& view)
// {
//   using PackT = typename ViewT::value_type;
//   EKAT_ASSERT(PackT::n >= 1);

//   const int pack_size = static_cast<int>(PackT::n);
//   const int npack     = (dim2_size+pack_size-1)/pack_size;

// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
//      const int num_scalars = k*pack_size;
//      const int scalar_offset = i*dim2_size+num_scalars;
//      if (num_scalars+s<dim2_size)  view(i,k)[s] = data[scalar_offset+s];
//   });
// }

// // 3D YAKL to Kokkos views
// template <typename SizeT, typename ViewT>
// void array_to_view(const typename ViewT::value_type::scalar* const data,
//                    const SizeT& dim1_size,
//                    const SizeT& dim2_size,
//                    const SizeT& dim3_size,
//                    ViewT& view)
// {
//   using PackT = typename ViewT::value_type;
//   EKAT_ASSERT(PackT::n >= 1);

//   const int pack_size = static_cast<int>(PackT::n);
//   const int npack     = (dim3_size+pack_size-1)/pack_size;

// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
//      const int num_scalars = k*pack_size;
//      const int scalar_offset = (i*dim2_size+j)*dim3_size+num_scalars;
//      if (num_scalars+s<dim3_size) view(i,j,k)[s] = data[scalar_offset+s];
//   });
// }

// // 4D YAKL to Kokkos views
// template <typename SizeT, typename ViewT>
// void array_to_view(const typename ViewT::value_type::scalar* const data,
//                    const SizeT& dim1_size,
//                    const SizeT& dim2_size,
//                    const SizeT& dim3_size,
//                    const SizeT& dim4_size,
//                    ViewT& view)
// {
//   using PackT = typename ViewT::value_type;
//   EKAT_ASSERT(PackT::n >= 1);

//   const int pack_size = static_cast<int>(PackT::n);
//   const int npack     = (dim3_size+pack_size-1)/pack_size;

// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<5>>({0, 0, 0, 0, 0}, {dim1_size, dim2_size, dim3_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int l, int s) {
//      const int num_scalars = k*pack_size;
//      const int scalar_offset = ((i*dim2_size+j)*dim3_size+l)*dim4_size+num_scalars;
//      if (num_scalars+s<dim4_size) view(i,j,k,l)[s] = data[scalar_offset+s];
//   });
// }

// // 5D YAKL to Kokkos views
// template <typename SizeT, typename ViewT>
// void array_to_view(const typename ViewT::value_type::scalar* const data,
//                    const SizeT& dim1_size,
//                    const SizeT& dim2_size,
//                    const SizeT& dim3_size,
//                    const SizeT& dim4_size,
//                    const SizeT& dim5_size,
//                    ViewT& view)
// {
//   using PackT = typename ViewT::value_type;
//   EKAT_ASSERT(PackT::n >= 1);

//   const int pack_size = static_cast<int>(PackT::n);
//   const int npack     = (dim3_size+pack_size-1)/pack_size;

// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<6>>({0, 0, 0, 0, 0, 0}, {dim1_size, dim2_size, dim3_size, dim4_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int l, int m, int s) {
//      const int num_scalars = k*pack_size;
//      const int scalar_offset = (((i*dim2_size+j)*dim3_size+l)*dim4_size+m)*dim5_size+num_scalars;
//      if (num_scalars+s<dim5_size) view(i,j,k,l,m)[s] = data[scalar_offset+s];
//   });
// }

// // 1D Kokkos view to YAKL array
// template <typename SizeT, typename ViewT, typename ArrayT>
// void view_to_array(const ViewT& view,
//                    const SizeT& size,
//                    ArrayT& array)
// {
//   using PackT      = typename ViewT::value_type;
//   using scalarType = typename ViewT::value_type::scalar;
//   using arrayType  = typename ArrayT::type;

//   auto is_same_type = std::is_same<scalarType, arrayType>::value;

//   EKAT_ASSERT(is_same_type);
//   EKAT_ASSERT(PackT::n >= 1);

//   const size_t npack = (size + PackT::n-1)/PackT::n;
    
// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {npack, PackT::n}), KOKKOS_LAMBDA(int k, int s) {
//       const size_t scalar_offset = k*PackT::n;
//       if (scalar_offset+s < size) array.myData[scalar_offset+s] = view(k)[s];
//   });
// }

// // 2D Kokkos view to YAKL array
// template <typename SizeT, typename ViewT, typename ArrayT>
// void view_to_array(const ViewT& view,
//                    const SizeT& dim1_size,
//                    const SizeT& dim2_size,
//                    ArrayT& array)
// {
//   using PackT      = typename ViewT::value_type;
//   using scalarType = typename ViewT::value_type::scalar;
//   using arrayType  = typename ArrayT::type;
 
//   auto is_same_type = std::is_same<scalarType, arrayType>::value;

//   EKAT_ASSERT(is_same_type);
//   EKAT_ASSERT(PackT::n >= 1);

//   const int pack_size = static_cast<int>(PackT::n);
//   const int npack     = (dim2_size+pack_size-1)/pack_size;

// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
//     const int num_scalars = k*pack_size;
//     const int scalar_offset = i*dim2_size + num_scalars;
//     if (num_scalars+s<dim2_size)  array.myData[scalar_offset+s] = view(i,k)[s];
//   });
// }

// // 3D Kokkos view to YAKL array
// template <typename SizeT, typename ViewT, typename ArrayT>
// void view_to_array(const ViewT& view,
//                    const SizeT& dim1_size,
//                    const SizeT& dim2_size,
//                    const SizeT& dim3_size,
//                    ArrayT& array)
// {
//   using PackT      = typename ViewT::value_type;
//   using scalarType = typename ViewT::value_type::scalar;
//   using arrayType  = typename ArrayT::type;

//   auto is_same_type = std::is_same<scalarType, arrayType>::value;

//   EKAT_ASSERT(is_same_type);
//   EKAT_ASSERT(PackT::n >= 1);

//   const int pack_size = static_cast<int>(PackT::n);
//   const int npack     = (dim3_size+pack_size-1)/pack_size;

// #if defined(DEBUG)
//   kokkos_impl_cuda_set_serial_execution(true);
// #endif
//   Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {dim1_size, dim2_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int j, int k, int s) {
//     const int num_scalars = k*pack_size;
//     const int scalar_offset = (i*dim2_size+j)*dim3_size+num_scalars;
//     if (num_scalars+s<dim3_size)  array.myData[scalar_offset+s] = view(i,j,k)[s];
//   });
// }

// // validation code 
// template <typename SizeT, typename ViewT>
// void array_to_view_2d(typename ViewT::value_type::scalar* data,
//                       const SizeT& dim1_sizes,
//                       const SizeT& dim2_sizes,
//                       ViewT& view)
// {
//   using PackT = typename ViewT::value_type;
//   EKAT_ASSERT(PackT::n >= 1);

//     const int dim1_size = static_cast<int>(dim1_sizes);
//     const int dim2_size = static_cast<int>(dim2_sizes);
//     const int pack_size = static_cast<int>(PackT::n);
//     const int npack     = (dim2_size+pack_size-1)/pack_size;

// #if defined(DEBUG)
//     kokkos_impl_cuda_set_serial_execution(true);
// #endif
//     Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {dim1_size, npack, pack_size}), KOKKOS_LAMBDA(int i, int k, int s) {
//        const int num_scalars = k*pack_size;
//        const int scalar_offset = i*dim2_size + num_scalars;
//        if (num_scalars+s<dim2_size) view(i,k)[s] = data[scalar_offset+s];
//    });
// }

