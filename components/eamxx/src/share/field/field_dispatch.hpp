#ifndef SCREAM_FIELD_DISPATCH_HPP
#define SCREAM_FIELD_DISPATCH_HPP

#include "share/field/field.hpp"

#include <type_traits>

namespace scream
{

namespace details {

// Tag associating a field with the scalar type its view should have (e.g., 'const double').
template<typename T>
struct FieldAs {
  const Field& field;
};

template<typename T>
FieldAs<T> as (const Field& f) { return FieldAs<T>{f}; }

// Call f(rank_c), with rank_c a std::integral_constant<int,N>, where N is the runtime
// value 'rank' (which must be in [0,Field::MaxRank]). This turns the runtime rank
// into a compile-time constant.
template<typename F>
void dispatch_rank (const int rank, F&& f)
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
void dispatch_view (const Field& field, F&& f)
{
  using DT = Field::data_nd_t<T,N>;
  if (field.get_header().get_alloc_properties().allows_layout_right())
    f(field.get_view<DT>());
  else
    f(field.get_strided_view<DT>());
}

// Call f(v1,v2,...), with vi the rank-N view of the i-th field (strided or not, see
// dispatch_view), and with the scalar type specified by the i-th FieldAs tag.
// This instantiates f for all combinations of strided/non-strided views.
template<int N, typename F>
void dispatch_views (F&& f)
{
  f();
}

template<int N, typename F, typename T, typename... Rest>
void dispatch_views (F&& f, const FieldAs<T>& first, const Rest&... rest)
{
  dispatch_view<N,T>(first.field, [&](const auto& v) {
    dispatch_views<N>([&](const auto&... vs) { f(v,vs...); }, rest...);
  });
}

// Runtime-rank version of the above. Example:
//
//   dispatch_views(rank,
//                  [&](const auto& y_view, const auto& x_view) { ... },
//                  as<double>(y), as<const float>(x));
//
// Here, y_view/x_view are rank-'rank' views, with value types double/const float
template<typename F, typename... Fields>
void dispatch_views (const int rank, F&& f, const FieldAs<Fields>&... fields)
{
  dispatch_rank(rank,[&](auto rank_c) {
    dispatch_views<decltype(rank_c)::value>(f,fields...);
  });
}

} // namespace details

} // namespace scream

#endif // SCREAM_FIELD_DISPATCH_HPP
