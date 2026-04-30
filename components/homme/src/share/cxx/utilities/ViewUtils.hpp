/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_VIEW_UTILS_HPP
#define HOMMEXX_VIEW_UTILS_HPP

#include "Types.hpp"
#include "ExecSpaceDefs.hpp"

#include <ekat_pack_kokkos.hpp>

namespace Homme {

// Structure to define the type of a view that has a const data type,
// given the type of an input view
template<typename ViewT>
struct ViewConst{};

// Note:: using ViewType::const_type may add explicit template arguments to ViewType.
//        Instead, we want the same 'template signature', simply with a const data type.
template<typename DataType, typename...Props>
struct ViewConst<ViewType<DataType,Props...>> {
  using type = ViewType<const DataType,Props...>;
};

// This is ugly, but prevents unnecessary copies
template<typename ViewT>
KOKKOS_INLINE_FUNCTION
typename ViewConst<ViewT>::type
viewConst(const ViewT& v) {
  return reinterpret_cast<const typename ViewConst<ViewT>::type&>(v);
}

template<int NUM_LEVELS>
KOKKOS_INLINE_FUNCTION
void print_col (const char* prefix,
                const ExecViewUnmanaged<const Scalar[ColInfo<NUM_LEVELS>::NumPacks]>& v) {
  printf("%s:",prefix);
  for (int k=0; k<NUM_LEVELS; ++k) {
    const int ilev = k / VECTOR_SIZE;
    const int ivec = k % VECTOR_SIZE;
    printf(" %3.15f",v(ilev)[ivec]);
  }
  printf("\n");
}

template<int NUM_LEVELS>
KOKKOS_INLINE_FUNCTION
void print_col (const char* prefix, const ExecViewUnmanaged<const Real[NUM_LEVELS]>& v) {
  printf("%s:",prefix);
  for (int k=0; k<NUM_LEVELS; ++k) {
    printf(" %3.15f",v(k));
  }
  printf("\n");
}

} // namespace Homme

#endif // HOMMEXX_VIEW_UTILS_HPP
