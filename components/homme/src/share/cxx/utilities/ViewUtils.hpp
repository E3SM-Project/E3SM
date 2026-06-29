/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_VIEW_UTILS_HPP
#define HOMMEXX_VIEW_UTILS_HPP

#include "Types.hpp"
#include "ExecSpaceDefs.hpp"

namespace Homme {

// This helper struct (and its shorter alias) simply provide
// the map 'Scalar->Real' and 'const Scalar->const Real'
namespace Impl {
template<typename In, typename Out>
struct ConstIfConst {
  static constexpr bool is_const = std::is_const<In>::value;
  using type = typename std::conditional<is_const,
                  typename std::add_const<Out>::type,
                  typename std::remove_const<Out>::type
               >::type;
};
} // namespace Impl

template<typename ScalarType>
using RealType = typename Impl::ConstIfConst<ScalarType,Real>::type;

// ================ Reinterpret a view of Scalar as a view of Real ======================= //
// Note: we template on ScalarType to allow both const and non-const, but the underlying
//       type (the one you get with std::remove_const) *must* be Scalar (as defined in Types.hpp).
template <typename ScalarType, int DIM1, typename... Properties>
KOKKOS_INLINE_FUNCTION
typename
std::enable_if<std::is_same<typename std::remove_const<ScalarType>::type,Scalar>::value,
               Unmanaged<ViewType<RealType<ScalarType>[DIM1*VECTOR_SIZE],Properties...>>
              >::type
viewAsReal(ViewType<ScalarType [DIM1], Properties...> v_in) {
  using ReturnST = RealType<ScalarType>;
  using ReturnView = Unmanaged<ViewType<RealType<ScalarType>[DIM1*VECTOR_SIZE],Properties...>>;
  return ReturnView(reinterpret_cast<ReturnST*>(v_in.data()));
}

template <typename ScalarType, int DIM1, int DIM2, typename... Properties>
KOKKOS_INLINE_FUNCTION
typename
std::enable_if<std::is_same<typename std::remove_const<ScalarType>::type,Scalar>::value,
               Unmanaged<ViewType<RealType<ScalarType>[DIM1][DIM2*VECTOR_SIZE],Properties...>>
              >::type
viewAsReal(ViewType<ScalarType [DIM1][DIM2], Properties...> v_in) {
  using ReturnST = RealType<ScalarType>;
  using ReturnView = Unmanaged<ViewType<RealType<ScalarType>[DIM1][DIM2*VECTOR_SIZE],Properties...>>;
  return ReturnView(reinterpret_cast<ReturnST*>(v_in.data()));
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, typename... Properties>
KOKKOS_INLINE_FUNCTION
typename
std::enable_if<std::is_same<typename std::remove_const<ScalarType>::type,Scalar>::value,
               Unmanaged<ViewType<RealType<ScalarType>[DIM1][DIM2][DIM3*VECTOR_SIZE],Properties...>>
              >::type
viewAsReal(ViewType<ScalarType [DIM1][DIM2][DIM3], Properties...> v_in) {
  using ReturnST = RealType<ScalarType>;
  using ReturnView = Unmanaged<ViewType<RealType<ScalarType>[DIM1][DIM2][DIM3*VECTOR_SIZE],Properties...>>;
  return ReturnView(reinterpret_cast<ReturnST*>(v_in.data()));
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4, typename... Properties>
KOKKOS_INLINE_FUNCTION
typename
std::enable_if<std::is_same<typename std::remove_const<ScalarType>::type,Scalar>::value,
               Unmanaged<ViewType<RealType<ScalarType>[DIM1][DIM2][DIM3][DIM4*VECTOR_SIZE],Properties...>>
              >::type
viewAsReal(ViewType<ScalarType [DIM1][DIM2][DIM3][DIM4], Properties...> v_in) {
  using ReturnST = RealType<ScalarType>;
  using ReturnView = Unmanaged<ViewType<RealType<ScalarType>[DIM1][DIM2][DIM3][DIM4*VECTOR_SIZE],Properties...>>;
  return ReturnView(reinterpret_cast<ReturnST*>(v_in.data()));
}

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
