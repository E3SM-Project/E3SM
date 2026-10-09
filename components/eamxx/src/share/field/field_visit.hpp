#ifndef SCREAM_FIELD_VISIT_HPP
#define SCREAM_FIELD_VISIT_HPP

#include "share/field/field.hpp"

#include <type_traits>

namespace scream
{

namespace details {

// Tag associating a field with the scalar type its view should have (e.g., 'const double').
template<typename T>
struct ViewAs {
  const Field& field;
};

template<typename T>
ViewAs<T> view_as (const Field& f) { return ViewAs<T>{f}; }

// Call f(rank_c), with rank_c a std::integral_constant<int,N>, where N is the runtime
// value 'rank' (which must be in [0,Field::MaxRank]). This turns the runtime rank
// into a compile-time constant.
template<typename F>
void dispatch_on_rank (const int rank, F&& f)
{
  switch (rank) {
    case 0: f(std::integral_constant<int,0>()); break;
    case 1: f(std::integral_constant<int,1>()); break;
    case 2: f(std::integral_constant<int,2>()); break;
    case 3: f(std::integral_constant<int,3>()); break;
    case 4: f(std::integral_constant<int,4>()); break;
    case 5: f(std::integral_constant<int,5>()); break;
    case 6: f(std::integral_constant<int,6>()); break;
    default:
      EKAT_ERROR_MSG ("Error! Field rank not supported (must be in [0,6]).\n"
          " - rank: " + std::to_string(rank) + "\n");
  }
}

// Call f(view), with view being the rank-N view of type T of the field.
// The view is a strided one if the field does not allow layout right
// (e.g., it is a subfield of another field), and a LayoutRight one otherwise.
template<int N, typename T, typename F>
void visit_field_view (const Field& field, F&& f)
{
  using DT = Field::data_nd_t<T,N>;
  if (field.get_header().get_alloc_properties().allows_layout_right())
    f(field.get_view<DT>());
  else
    f(field.get_strided_view<DT>());
}

// Call f(v1,v2,...), with vi the rank-N view of the i-th field (strided or not, see
// visit_field_view), and with the scalar type specified by the i-th ViewAs tag.
// This instantiates f for all combinations of strided/non-strided views.
template<int N, typename F>
void visit_field_views (F&& f)
{
  f();
}

template<int N, typename F, typename T, typename... Rest>
void visit_field_views (F&& f, const ViewAs<T>& first, const Rest&... rest)
{
  visit_field_view<N,T>(first.field, [&](const auto& v) {
    visit_field_views<N>([&](const auto&... vs) { f(v,vs...); }, rest...);
  });
}

// Runtime-rank version of the above. Example:
//
//   visit_field_views(rank,
//                     [&](const auto& y_view, const auto& x_view) { ... },
//                     view_as<double>(y), view_as<const float>(x));
//
// Here, y_view/x_view are rank-'rank' views, with value types double/const float
template<typename F, typename... Ts>
void visit_field_views (const int rank, F&& f, const ViewAs<Ts>&... fields)
{
  static_assert (sizeof...(Ts)>0,
      "[visit_field_views] Error! At least one field must be provided.\n");

  dispatch_on_rank(rank,[&](auto rank_c) {
    visit_field_views<decltype(rank_c)::value>(f,fields...);
  });
}

} // namespace details

} // namespace scream

#endif // SCREAM_FIELD_VISIT_HPP
