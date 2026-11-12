/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_SUBVIEW_UTILS_HPP
#define HOMMEXX_SUBVIEW_UTILS_HPP

#include "Types.hpp"
#include "ExecSpaceDefs.hpp"
#include "ViewUtils.hpp"

#include <functional>

namespace Homme {

// ================ Subviews of several ranks views ======================= //
// Note: we template on ScalarType to allow both Real and Scalar case, and
//       also to allow const/non-const versions.
// Note: we assume to have exactly one runtime dimension.
template <typename ScalarType, int DIM1, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM1], MemSpace>
subview(ViewType<ScalarType * [DIM1], MemSpace, Properties...> v_in, int ie) {
  return Kokkos::subview(v_in, ie, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM2], MemSpace>
subview(ViewType<ScalarType[DIM1][DIM2], MemSpace, Properties...> v_in,
        int idx_1) {
  return Kokkos::subview(v_in, idx_1, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM1][DIM2], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2], MemSpace, Properties...> v_in,
        int ie) {
  return Kokkos::subview(v_in, ie, Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM2], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2], MemSpace, Properties...> v_in,
        int idx_0, const int idx_1) {
  return Kokkos::subview(v_in, idx_0, idx_1, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM3], MemSpace>
subview(ViewType<ScalarType[DIM1][DIM2][DIM3], MemSpace, Properties...> v_in,
        int idx_1, int idx_2) {
  return Kokkos::subview(v_in, idx_1, idx_2, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM1][DIM2][DIM3], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3], MemSpace, Properties...> v_in,
        int ie) {
  return Kokkos::subview(v_in, ie, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM2][DIM3], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3], MemSpace, Properties...> v_in,
        int ie, int idim1) {
  return Kokkos::subview(v_in, ie, idim1, Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, typename MemSpace,
          typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM3], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3], MemSpace, Properties...> v_in,
        int ie, int idim1, int idim2) {
  return Kokkos::subview(v_in, ie, idim1, idim2, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION
ViewUnmanaged<ScalarType[DIM1][DIM2][DIM3][DIM4], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4], MemSpace, Properties...>
            v_in, int ie) {
  return Kokkos::subview(v_in, ie, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM2][DIM3][DIM4], MemSpace>
subview(
    ViewType<ScalarType[DIM1][DIM2][DIM3][DIM4], MemSpace, Properties...> v_in,
    int var) {
  return Kokkos::subview(v_in, var, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM4], MemSpace>
subview(ViewType<ScalarType [DIM1][DIM2][DIM3][DIM4], MemSpace, Properties...>
            v_in, int idim1, int idim2, int idim3) {
  return Kokkos::subview(v_in, idim1, idim2, idim3, Kokkos::ALL);
}


template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM4], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4], MemSpace, Properties...>
            v_in, int ie, int tl, int igp, int jgp) {
  return Kokkos::subview(v_in, ie, tl, igp, jgp, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM2][DIM3][DIM4], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4], MemSpace, Properties...>
            v_in, int ie, int idim1) {
  return Kokkos::subview(v_in, ie, idim1, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM3][DIM4], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4], MemSpace, Properties...>
            v_in, int ie, int idim1, int idim2) {
  return Kokkos::subview(v_in, ie, idim1, idim2, Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION
ViewUnmanaged<ScalarType[DIM1][DIM2][DIM3][DIM4][DIM5], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4][DIM5], MemSpace,
                 Properties...> v_in,
        int ie) {
  return Kokkos::subview(v_in, ie, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION
ViewUnmanaged<ScalarType[DIM2][DIM3][DIM4][DIM5], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4][DIM5], MemSpace,
                 Properties...> v_in,
        int ie, int idim1) {
  return Kokkos::subview(v_in, ie, idim1, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM3][DIM4][DIM5], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4][DIM5], MemSpace,
                 Properties...> v_in,
        int ie, int idim1, int idim2) {
  return Kokkos::subview(v_in, ie, idim1, idim2, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM5], MemSpace>
subview(ViewType<ScalarType * [DIM1][DIM2][DIM3][DIM4][DIM5], MemSpace,
                 Properties...> v_in,
        int ie, int idim1, int idim2, int idim3, int idim4) {
  return Kokkos::subview(v_in, ie, idim1, idim2, idim3, idim4, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType * [DIM1][DIM2][DIM3], MemSpace>
subview(ViewType<ScalarType ** [DIM1][DIM2][DIM3], MemSpace,
                 Properties...> v_in,
        int ie) {
  return Kokkos::subview(v_in, ie, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType * [DIM1][DIM2][DIM3][DIM4], MemSpace>
subview(ViewType<ScalarType ** [DIM1][DIM2][DIM3][DIM4], MemSpace,
                 Properties...> v_in,
        int ie) {
  return Kokkos::subview(v_in, ie, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM1][DIM2][DIM3], MemSpace>
subview(ViewType<ScalarType ** [DIM1][DIM2][DIM3], MemSpace,
                 Properties...> v_in,
        int ie, int remap_idx) {
  return Kokkos::subview(v_in, ie, remap_idx, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM3], MemSpace>
subview(ViewType<ScalarType ** [DIM1][DIM2][DIM3], MemSpace,
                 Properties...> v_in,
        int ie, int remap_idx, int idim1, int idim2) {
  return Kokkos::subview(v_in, ie, remap_idx, idim1, idim2, Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM3][DIM4], MemSpace>
subview(ViewType<ScalarType ** [DIM1][DIM2][DIM3][DIM4], MemSpace,
                 Properties...> v_in,
        int ie, int remap_idx, int idim1, int idim2) {
  return Kokkos::subview(v_in, ie, remap_idx, idim1, idim2, Kokkos::ALL,
                         Kokkos::ALL);
}

template <typename ScalarType, int DIM1, int DIM2, int DIM3, int DIM4,
          typename MemSpace, typename... Properties>
KOKKOS_INLINE_FUNCTION ViewUnmanaged<ScalarType[DIM1][DIM2][DIM3][DIM4], MemSpace>
subview(ViewType<ScalarType ** [DIM1][DIM2][DIM3][DIM4], MemSpace,
                 Properties...> v_in,
        int ie, int idim0) {
  return Kokkos::subview(v_in, ie, idim0, Kokkos::ALL, Kokkos::ALL,
                         Kokkos::ALL, Kokkos::ALL);
}

// Force a subview to be const
template<typename View, typename... Ints>
KOKKOS_INLINE_FUNCTION
auto
subviewConst(const View& v, const Ints... idx) ->
  // The compile will match the correct method above.
 typename ViewConst<decltype(Homme::subview(v,idx...))>::type
{
  // Upon return, the view will be copied into one with const scalar type
  return Homme::subview(v,idx...);
}

} // namespace Homme

#endif // HOMMEXX_SUBVIEW_UTILS_HPP
