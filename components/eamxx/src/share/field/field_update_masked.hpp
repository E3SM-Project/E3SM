#ifndef SCREAM_FIELD_UPDATE_MASKED_HPP
#define SCREAM_FIELD_UPDATE_MASKED_HPP

#include "share/field/field.hpp"
#include "share/field/field_pfor.hpp"

namespace scream
{

namespace details {

// NOTE: kernels that go through field_pfor must be functors with a templated operator()
//       (not generic lambdas), since extended lambdas cannot be generic.
template<CombineMode CM, typename ST>
struct UpdateMaskedKernel {
  ST alpha;
  ST beta;
  ST gamma;

  template<typename YView, typename XView, typename MaskView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const XView& x, const MaskView& mask, const Idx&... i) const {
    if (mask(i...)) {
      auto& y_val = y(i...);
      combine<CM>(x(i...),y_val,alpha,beta);
      if constexpr (CM==CombineMode::Update)
        y_val += gamma;
    }
  }
};

template<bool negate_mask, typename ST>
struct SetValueMaskedKernel {
  ST value;

  template<typename YView, typename MaskView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const MaskView& mask, const Idx&... i) const {
    // Set the value where mask=1 (or, if negate_mask=true, where mask=0)
    const bool selected = negate_mask ? (mask(i...)==0) : (mask(i...)!=0);
    if (selected)
      y(i...) = value;
  }
};

} // namespace details

template<CombineMode CM, typename ST, typename XST>
void Field::
update_masked (const Field& x, const ST alpha, const ST beta, const ST gamma, const Field& mask) const
{
  // Must handle the case where any of the views is strided
  details::field_pfor(details::UpdateMaskedKernel<CM,ST>{alpha,beta,gamma},
                      details::view_as<ST>(*this), details::view_as<const XST>(x),
                      details::view_as<const int>(mask));

  // field_pfor is asynchronous, but the Field update methods are blocking
  Kokkos::fence();
}

template<bool negate_mask, typename ST>
void Field::deep_copy_masked (const ST value, const Field& mask) const
{
  EKAT_REQUIRE_MSG (mask.is_allocated(),
      "Error! Calling Field::deep_copy_masked with a mask field that is not allocated.\n"
      " - mask field name: " + mask.name() + "\n");

  details::field_pfor(details::SetValueMaskedKernel<negate_mask,ST>{value},
                      details::view_as<ST>(*this), details::view_as<const int>(mask));

  // field_pfor is asynchronous, but the Field update methods are blocking
  Kokkos::fence();
}

} // namespace scream

#endif // SCREAM_FIELD_UPDATE_MASKED_HPP
