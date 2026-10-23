#ifndef SCREAM_FIELD_PFOR_HPP
#define SCREAM_FIELD_PFOR_HPP

#include "share/field/field_visit.hpp"
#include "share/util/eamxx_bound_functor.hpp"
#include "share/util/eamxx_parallel_for_nd.hpp"

#include <tuple>
#include <type_traits>

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
// and N is the (runtime) rank of the fields. The views have types that depend on the rank
// and layout of the fields, so func must be a functor with a templated operator(). A kernel
// that is valid for any rank is
//
//   struct AddKernel {
//     template<typename YView, typename XView, typename... Idx>
//     KOKKOS_INLINE_FUNCTION
//     void operator() (const YView& y, const XView& x, const Idx&... i) const {
//       y(i...) += x(i...);
//     }
//   };
//   ...
//   field_pfor(AddKernel(), view_as<Real>(y_field), view_as<const Real>(x_field));
//
// Notes:
//  - func CANNOT be a generic KOKKOS_LAMBDA: nvcc does not allow extended __host__ __device__
//    lambdas to be generic (nor to first-capture a variable inside an 'if constexpr').
//    A functor struct with a templated operator() has neither restriction.
//  - Kokkos only calls a parallel_for functor with the indices. So we don't hand func to
//    Kokkos: we hand it a BoundFunctor, that stores func and the visited views, and
//    forwards (views..., indices...) to func.
//  - For rank 0, func is called with no indices (so y(i...) is y()).
//  - All fields must have the same rank and extents (tags are not checked). The kernel is
//    instantiated for every rank in [0,MaxRank] and every layout combination (see
//    field_visit.hpp), so it should be small. If the kernel needs to treat ranks differently,
//    it can use 'if constexpr (sizeof...(i)==2)' or the likes.
//  - func is copied to device, so it must be device-callable and device-copyable.
//  - Like Kokkos::parallel_for, this call is asynchronous: it does not fence. Kernels launched on
//    the same execution space run in order, so consecutive calls can be chained without fencing.
//    Fence (or sync the fields to host) before using the results on host.
// ===========================================================================

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
    [&](const auto& first_view, const auto&... other_views) {
      constexpr int N = std::decay_t<decltype(first_view)>::rank();
      using exec_space = typename std::decay_t<decltype(first_view)>::traits::execution_space;

      parallel_for_nd<N,exec_space>(dims,bind_functor(func,first_view,other_views...));
    },
    fields...);
}

} // namespace details

} // namespace scream

#endif // SCREAM_FIELD_PFOR_HPP
