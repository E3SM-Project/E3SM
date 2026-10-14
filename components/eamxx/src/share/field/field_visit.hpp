#ifndef SCREAM_FIELD_VISIT_HPP
#define SCREAM_FIELD_VISIT_HPP

#include "share/field/field.hpp"

#include <type_traits>
#include <utility>
#include <variant>

namespace scream
{

namespace details {

// ===========================================================================
// Visiting the views of one or more fields
// ===========================================================================
//
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
//   visit_field_views(rank, func, view_as<T_1>(f_1), ..., view_as<T_K>(f_K))
//
// calls func(v_1,...,v_K), where v_i is a rank-'rank' view of f_i with scalar type T_i,
// in the same order as the view_as tags. Each v_i can independently be LayoutRight or
// LayoutStride, so func is instantiated for ALL 2^K layout combinations (the one
// matching the actual fields is the only one that runs). This is also true for the
// ranks: the runtime-rank overload instantiates func for each rank in [0,MaxRank].
// Since this multiplies compile time and code size, func should be small: have it
// forward to a templated helper (like cvh/cvmh) rather than hold a big body.
//
// How it works: a function can't return a LayoutRight view for some fields and a
// LayoutStride view for others, depending on a runtime check. But it can return a
// std::variant of the two view types (see get_view_variant). Once each field is
// turned into such a variant, std::visit(func, var_1, ..., var_K) does the rest: it
// calls func with the active alternative of each variant, instantiating func for the
// cross product of all alternatives. The runtime cost is one jump through std::visit's
// dispatch table per call (not per element); the views are copied into the variants
// (cheap: Kokkos views are reference counted handles), and only the host sees them:
// device kernels are launched inside the helper that func calls.
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

// Return the rank-N view of type T of the field, wrapped in a variant that holds
// a LayoutRight view if the field allows it, and a LayoutStride view otherwise
// (e.g., if the field is a subfield of another field).
template<int N, typename T>
auto get_view_variant (const Field& field)
{
  using DT = Field::data_nd_t<T,N>;
  using lr_view_t = std::decay_t<decltype(field.get_view<DT>())>;
  using st_view_t = std::decay_t<decltype(field.get_strided_view<DT>())>;
  using variant_t = std::variant<lr_view_t,st_view_t>;

  // Construct by index, so that we do not rely on the two view types being different
  if (field.get_header().get_alloc_properties().allows_layout_right())
    return variant_t(std::in_place_index<0>,field.get_view<DT>());
  else
    return variant_t(std::in_place_index<1>,field.get_strided_view<DT>());
}

// Call func(v_1,...,v_K), with v_i the rank-N view of the i-th field. See above for details.
template<int N, typename Func, typename... Ts>
void visit_field_views (Func&& func, const ViewAs<Ts>&... fields)
{
  std::visit(std::forward<Func>(func),get_view_variant<N,Ts>(fields.field)...);
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
