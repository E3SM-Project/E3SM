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

template<typename ScalarType>
struct ScalarTypeHelper
{
  static_assert(std::is_same<typename std::remove_const<ScalarType>::type,Scalar>::value,
                "Error! The ScalarTypeHelper struct can only be templated on 'Scalar' and 'const Scalar'.\n");

  static constexpr bool is_const = std::is_const<ScalarType>::value;
  using type = typename std::conditional<is_const,const Real, Real>::type;
};

// ================ Reinterpret a view of Scalar as a view of Real ======================= //
// Note: we template on ScalarType to allow both const and non-const, but the underlying
//       type (the one you get with std::remove_const) *must* be Scalar (as defined in Types.hpp).
template <typename ScalarType, int DIM1, typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION
ViewUnmanaged<typename ScalarTypeHelper<ScalarType>::type[DIM1*VECTOR_SIZE], MemSpace>
viewAsReal(ViewType<ScalarType [DIM1], MemSpace, Properties...> v_in) {
  using ReturnST = typename ScalarTypeHelper<ScalarType>::type;
  return ViewUnmanaged<ReturnST [DIM1*VECTOR_SIZE], MemSpace>(reinterpret_cast<ReturnST*>(v_in.data()));
}

template <typename ScalarType, int DIM1, int DIM2, typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION
ViewUnmanaged<typename ScalarTypeHelper<ScalarType>::type[DIM1][DIM2*VECTOR_SIZE], MemSpace>
viewAsReal(ViewType<ScalarType [DIM1][DIM2], MemSpace, Properties...> v_in) {
  using ReturnST = typename ScalarTypeHelper<ScalarType>::type;
  return ViewUnmanaged<ReturnST [DIM1][DIM2*VECTOR_SIZE], MemSpace>(reinterpret_cast<ReturnST*>(v_in.data()));
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION
ViewUnmanaged<typename ScalarTypeHelper<ScalarType>::type[DIM1][DIM2][DIM3*VECTOR_SIZE], MemSpace>
viewAsReal(ViewType<ScalarType [DIM1][DIM2][DIM3], MemSpace, Properties...> v_in) {
  using ReturnST = typename ScalarTypeHelper<ScalarType>::type;
  return ViewUnmanaged<ReturnST [DIM1][DIM2][DIM3*VECTOR_SIZE], MemSpace>(reinterpret_cast<ReturnST*>(v_in.data()));
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4, typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION
ViewUnmanaged<typename ScalarTypeHelper<ScalarType>::type[DIM1][DIM2][DIM3][DIM4*VECTOR_SIZE], MemSpace>
viewAsReal(ViewType<ScalarType [DIM1][DIM2][DIM3][DIM4], MemSpace, Properties...> v_in) {
  using ReturnST = typename ScalarTypeHelper<ScalarType>::type;
  return ViewUnmanaged<ReturnST [DIM1][DIM2][DIM3][DIM4*VECTOR_SIZE], MemSpace>(reinterpret_cast<ReturnST*>(v_in.data()));
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

} // namespace Homme

#endif // HOMMEXX_VIEW_UTILS_HPP
