#ifndef SCREAM_FIELD_UPDATE_FILL_AWARE_HPP
#define SCREAM_FIELD_UPDATE_FILL_AWARE_HPP

#include "share/field/field.hpp"
#include "share/field/field_pfor.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream
{

namespace details {

// NOTE: kernels that go through field_pfor must be functors with a templated operator()
//       (not generic lambdas), since extended lambdas cannot be generic.
template<CombineMode CM, typename ST>
struct UpdateFillAwareKernel {
  ST alpha;
  ST beta;
  ST gamma;

  template<typename YView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const XView& x, const Idx&... i) const {
    auto& y_val = y(i...);
    combine_fill_aware<CM>(x(i...),y_val,alpha,beta);
    if constexpr (CM==CombineMode::Update)
      if (gamma!=constants::fill_value<ST>)
        y_val += gamma;
  }
};

} // namespace details

template<CombineMode CM, typename ST, typename XST>
void Field::
update_fill_aware (const Field& x, const ST alpha, const ST beta, const ST gamma) const
{
  // Must handle the case where one of the two views is strided (or both)
  details::field_pfor(details::UpdateFillAwareKernel<CM,ST>{alpha,beta,gamma},
                      details::view_as<ST>(*this), details::view_as<const XST>(x));

  // field_pfor is asynchronous, but the Field update methods are blocking
  Kokkos::fence();
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_FILL_AWARE_HPP
