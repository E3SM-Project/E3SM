#ifndef SCREAM_FIELD_VISIT_HPP
#define SCREAM_FIELD_VISIT_HPP

#include "share/field/field.hpp"

#include <type_traits>

namespace scream
{

namespace details {

// ===========================================================================
// Visiting the views of one or more fields
// ===========================================================================
//
// THE PROBLEM
// To access a field's data in a kernel, we need a Kokkos view whose *type* depends
// on things we only know at runtime:
//   - the rank N of the field,
//   - the scalar type T we want to see (double, const float, ...),
//   - the layout: a plain LayoutRight view if the field's allocation allows it, or a
//     LayoutStride view otherwise (e.g., a subfield that is not contiguous in its
//     parent, like parent.subfield(1,k) of a rank-2 parent).
// Kernel code (usually a generic lambda that calls a templated helper) needs these
// as compile-time types. So we "visit": we branch on the runtime properties, and
// in each branch we call the user's callable with views of the right concrete types.
//
// WHAT THE USER GETS
//   visit_field_views(rank, func, view_as<T_1>(f_1), ..., view_as<T_K>(f_K))
// calls func(v_1,...,v_K), where v_i is a rank-'rank' view of f_i with scalar type T_i,
// in the same order as the view_as tags. Each vi can independently be LayoutRight or
// LayoutStride, so func is instantiated for ALL 2^K layout combinations (the one
// matching the actual fields is the only one that runs). This is also true for the
// ranks: the runtime-rank overload instantiates func for each rank in [0,MaxRank].
// Since this multiplies compile time and code size, func should be small: have it
// forward to a templated helper (like cvh/cvmh) rather than hold a big body.
//
// THE FOUR PIECES
//   dispatch_on_rank(rank, func)
//       Runtime int -> compile-time constant. Calls func(std::integral_constant<int,N>).
//       It knows nothing about fields or views.
//   visit_field_view<N,T>(field, func)                               [overload 1]
//       ONE field. It is the only place with the strided/non-strided if: it checks
//       allows_layout_right() and calls func(view) with the matching view type.
//   visit_field_views<N>(func, first, rest...)                       [overload 3]
//       Peels the first field: delegates it to visit_field_view, and recurses on rest.
//   visit_field_views<N>(func)                                       [overload 2]
//       No fields left: the end of the recursion. Just calls func().
//   (plus the runtime-rank entry point, which calls dispatch_on_rank and then
//    overload 3 with the compile-time N.)
//
// HOW THE RECURSION WORKS
// visit_field_view is not called "at the end"; it is called by overload 3, once per
// field, at every level. Each level picks the layout for its field, and then, from
// inside the callback that visit_field_view invokes, recurses on the remaining fields.
// To make sure the final call sees ALL the views, each level wraps 'func' in a new
// lambda that prepends its own view:
//
//     [&](const auto&... vs) { func(v, vs...); }
//
// ('v' is captured by reference). So 'func' at depth d is not the user's callable but
// the lambda from depth d-1, and the real callable is reached only when the chain
// unwinds. For two fields (y,x), with Ty/Tx their requested types:
//
//   visit_field_views<N>(func, y, x)                                    [overload 3]
//   |- visit_field_view<N,Ty>(y, A)       picks y's view type, calls A(yv), where
//   |    A(yv) = visit_field_views<N>(B, x)                             [overload 3]
//   |            with B(vs...) = func(yv, vs...)
//   |- visit_field_view<N,Tx>(x, C)       picks x's view type, calls C(xv), where
//   |    C(xv) = visit_field_views<N>(D)                                [overload 2]
//   |            with D(vs...) = B(xv, vs...)
//   '- D()  ->  B(xv)  ->  func(yv, xv)   <- base case reached: user's callable
//                                            finally runs, with views in order.
//
// NOTES
//  - Lambdas capture by reference, which is safe: every call is synchronous, and
//    each view lives in the stack frame of visit_field_view until func returns.
//  - No heap allocation, no std::function: these are plain lambda objects that the
//    compiler can inline. The cost is paid once per call, not per element.
//  - The lambdas are host code only. Device kernels are launched inside the
//    helper that func calls, which takes the views by value; nothing captured by
//    reference ever reaches the device.
//  - The overload-2 base case could be merged into overload 3 via
//    'if constexpr (sizeof...(rest)==0)'. We keep them separate for readability.
// ===========================================================================

// Tag associating a field with the scalar type its view should have (e.g., 'const double').
template<typename T>
struct ViewAs {
  const Field& field;
};

template<typename T>
ViewAs<T> view_as (const Field& f) { return ViewAs<T>{f}; }

// Call func(rank_c), with rank_c a std::integral_constant<int,N>, where N is the runtime
// value 'rank' (which must be in [0,Field::MaxRank]). This turns the runtime rank
// into a compile-time constant.
template<typename Func>
void dispatch_on_rank (const int rank, Func&& func)
{
  switch (rank) {
    case 0: func(std::integral_constant<int,0>()); break;
    case 1: func(std::integral_constant<int,1>()); break;
    case 2: func(std::integral_constant<int,2>()); break;
    case 3: func(std::integral_constant<int,3>()); break;
    case 4: func(std::integral_constant<int,4>()); break;
    case 5: func(std::integral_constant<int,5>()); break;
    case 6: func(std::integral_constant<int,6>()); break;
    default:
      EKAT_ERROR_MSG ("Error! Field rank not supported (must be in [0,6]).\n"
          " - rank: " + std::to_string(rank) + "\n");
  }
}

// Call func(view), with view being the rank-N view of type T of the field.
// The view is a strided one if the field does not allow layout right
// (e.g., it is a subfield of another field), and a LayoutRight one otherwise.
template<int N, typename T, typename Func>
void visit_field_view (const Field& field, Func&& func)
{
  using DT = Field::data_nd_t<T,N>;
  if (field.get_header().get_alloc_properties().allows_layout_right())
    func(field.get_view<DT>());
  else
    func(field.get_strided_view<DT>());
}

template<int N, typename Func, typename T, typename... Rest>
void visit_field_views (Func&& func, const ViewAs<T>& first, const Rest&... rest)
{
  visit_field_view<N,T>(first.field, [&](const auto& v) {
    visit_field_views<N>([&](const auto&... vs) { func(v,vs...); }, rest...);
  });
}

template<int N, typename Func>
void visit_field_views (Func&& func)
{
  // When we reach this instantation, `func` is a lambda that already captured all
  // the fields' views, and does not need to take any additional arg. Just run it.
  func();
}

// Runtime-rank version of visit_field_views
template<typename Func, typename... Ts>
void visit_field_views (const int rank, Func&& func, const ViewAs<Ts>&... fields)
{
  static_assert (sizeof...(Ts)>0,
      "[visit_field_views] Error! At least one field must be provided.\n");
  static_assert (sizeof...(Ts)<=3,
      "[visit_field_views] Error! We only support up to 3 fields.\n");

  dispatch_on_rank(rank,[&](auto rank_c) {
    visit_field_views<decltype(rank_c)::value>(func,fields...);
  });
}

} // namespace details

} // namespace scream

#endif // SCREAM_FIELD_VISIT_HPP
