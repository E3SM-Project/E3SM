#ifndef SCREAM_FIELD_PFOR_HPP
#define SCREAM_FIELD_PFOR_HPP

#include "share/field/field_visit.hpp"

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace scream
{

namespace details {

// ===========================================================================
// Rank-agnostic parallel_for over the entries of one or more fields
// ===========================================================================
//
//   field_pfor(func, view_as<T_1>(f_1), ..., view_as<T_K>(f_K))
//
// runs, for each multi-index (i_0,...,i_{N-1}) in the layout of the fields,
//
//   func(v_1,...,v_K, i_0,...,i_{N-1})
//
// where v_k is the rank-N view of f_k with scalar type T_k (see visit_field_views),
// and N is the (runtime) rank of the fields. So, a kernel that is valid for any rank is
//
//   field_pfor(KOKKOS_LAMBDA (const auto& y, const auto& x, auto... i) {
//                y(i...) += x(i...);
//              },
//              view_as<Real>(y_field), view_as<const Real>(x_field));
//
// Notes:
//  - Kokkos only ever calls a parallel_for functor with the indices (and nothing else).
//    So we do not hand func to Kokkos; we hand it an "adapter" (PForAdapter) that stores
//    the visited views and func, and whose operator() takes only the indices and forwards
//    (views..., indices...) to func. The adapter is the same for all ranks.
//  - For rank 0 Kokkos has no 0-dim policy, so we run a 1-entry RangePolicy, and the adapter
//    swallows the dummy index. Hence func is called with no indices, and y(i...) is y().
//  - All fields must have the same rank and extents (tags are not checked). The kernel is instantiated
//    for every rank in [0,MaxRank] and every layout combination (see field_visit.hpp)
//    so it should be small. If the kernel needs to treat ranks differently, it can use
//    'if constexpr (sizeof...(i)==2)' or the likes.
//  - func is copied to device, so it must be a KOKKOS_LAMBDA (or a device-callable functor)
//    and capture by value. Do not capture host references.
//  - This call is blocking (it fences before returning).
// ===========================================================================

// Run functor(i_0,...,i_{N-1}) for all indices in the range dims[0]x...xdims[N-1].
// For N=0, functor is called once, with a single dummy index (Kokkos has no 0-dim policy).
template<int N, typename ExecSpace, typename Functor>
void parallel_for_nd (const std::vector<int>& dims, const Functor& functor)
{
  static_assert (N>=0 && N<=Field::MaxRank,
      "[parallel_for_nd] Error! Unsupported rank (must be in [0,Field::MaxRank]).\n");
  EKAT_REQUIRE_MSG (static_cast<int>(dims.size())==N,
      "[parallel_for_nd] Error! Number of dimensions does not match the rank.\n"
      " - rank: " + std::to_string(N) + "\n"
      " - dims.size(): " + std::to_string(dims.size()) + "\n");

  if constexpr (N==0) {
    Kokkos::RangePolicy<ExecSpace> policy(0,1);
    Kokkos::parallel_for(policy,functor);
  } else if constexpr (N==1) {
    Kokkos::RangePolicy<ExecSpace> policy(0,dims[0]);
    Kokkos::parallel_for(policy,functor);
  } else {
    using MDRange = Kokkos::MDRangePolicy<
                      ExecSpace,
                      Kokkos::Rank<N,Kokkos::Iterate::Right,Kokkos::Iterate::Right>
                    >;
    Kokkos::Array<int64_t,N> lower, upper;
    for (int i=0; i<N; ++i) {
      lower[i] = 0;
      upper[i] = dims[i];
    }
    Kokkos::parallel_for(MDRange(lower,upper),functor);
  }
}

// Holds one view; the index makes the type unique, so that we can inherit from several
// of them even if the views have the same type. Inheritance (rather than std::tuple)
// keeps the adapter trivially usable on device.
template<std::size_t I, typename View>
struct PForViewHolder {
  View view;
};

template<typename IndexSeq, typename Func, typename... Views>
struct PForAdapter;

// Functor that Kokkos can call with the indices only. It calls
// func(views..., indices...), with exactly N indices (none, if N=0)
template<std::size_t... Is, typename Func, typename... Views>
struct PForAdapter<std::index_sequence<Is...>,Func,Views...>
 : PForViewHolder<Is,Views>...
{
  static constexpr int N = std::tuple_element_t<0,std::tuple<Views...>>::rank();
  static_assert (((static_cast<int>(Views::rank())==N) && ...),
      "[PForAdapter] Error! All views must have the same rank.\n");

  PForAdapter (const Func& f, const Views&... v)
   : PForViewHolder<Is,Views>{v}...
   , func(f)
  {}

  template<typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const Idx&... idx) const {
    if constexpr (N==0) {
      // Swallow the dummy index
      func(static_cast<const PForViewHolder<Is,Views>&>(*this).view...);
    } else {
      static_assert (sizeof...(Idx)==N,
          "[PForAdapter] Error! Wrong number of indices.\n");
      func(static_cast<const PForViewHolder<Is,Views>&>(*this).view...,idx...);
    }
  }

  Func func;
};

template<typename Func, typename... Views>
using PForAdapterFor = PForAdapter<std::index_sequence_for<Views...>,Func,Views...>;

template<typename Func, typename... Ts>
void field_pfor (const Func& func, const ViewAs<Ts>&... fields)
{
  static_assert (sizeof...(Ts)>0,
      "[field_pfor] Error! At least one field must be provided.\n");

  const Field* const all_fields[] = { &fields.field... };
  const auto& layout = all_fields[0]->get_header().get_identifier().get_layout();
  for (const auto* f : all_fields) {
    const auto& fl = f->get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG (fl.dims()==layout.dims(),
        "[field_pfor] Error! All fields must have the same extents.\n"
        " - field 0: " + all_fields[0]->name() + " (" + layout.to_string() + ")\n"
        " - field  : " + f->name() + " (" + fl.to_string() + ")\n");
  }

  const auto& dims = layout.dims();
  visit_field_views(layout.rank(),
    [&](const auto&... views) {
      using adapter_t = PForAdapterFor<Func,std::decay_t<decltype(views)>...>;
      using exec_space = typename std::tuple_element_t<0,std::tuple<std::decay_t<decltype(views)>...>>::traits::execution_space;
      parallel_for_nd<adapter_t::N,exec_space>(dims,adapter_t(func,views...));
    },
    fields...);
  Kokkos::fence();
}

} // namespace details

} // namespace scream

#endif // SCREAM_FIELD_PFOR_HPP
