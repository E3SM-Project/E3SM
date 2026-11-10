#ifndef SCREAM_FIELD_UPDATE_IMPL_HPP
#define SCREAM_FIELD_UPDATE_IMPL_HPP

#include "share/field/field.hpp"
#include "share/field/field_pfor.hpp"

namespace scream
{

namespace details {

// NOTE: kernels that go through field_pfor must be functors with a templated operator()
//       (not generic lambdas), since extended lambdas cannot be generic.
template<CombineMode CM, typename ST>
struct UpdateKernel {
  ST alpha;
  ST beta;
  ST gamma;

  template<typename YView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const XView& x, const Idx&... i) const {
    auto& y_val = y(i...);
    combine<CM>(x(i...),y_val,alpha,beta);
    if constexpr (CM==CombineMode::Update)
      y_val += gamma;
  }
};

} // namespace details

template<CombineMode CM, typename ST, typename XST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta, const ST gamma) const
{
  // Must handle the case where one of the two views is strided (or both)
  details::field_pfor(details::UpdateKernel<CM,ST>{alpha,beta,gamma},
                      details::view_as<ST>(*this), details::view_as<const XST>(x));

  // field_pfor is asynchronous, but the Field update methods are blocking
  Kokkos::fence();
}

template<typename ST>
void Field::deep_copy_impl (const ST value) const
{
  details::visit_field_views(get_header().get_identifier().get_layout().rank(),
    [&](const auto& view) {
      Kokkos::deep_copy(view,value);
    },
    details::view_as<ST>(*this));
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_IMPL_HPP
