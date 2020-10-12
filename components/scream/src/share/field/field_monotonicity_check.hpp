#ifndef SCREAM_FIELD_MONOTONICITY_CHECK_HPP
#define SCREAM_FIELD_MONOTONICITY_CHECK_HPP

#include "share/field/field.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// monotonically increasing or decreasing, and false if not. Currently, no
// repairs to nonmonotonic fields can be performed. There are ways of repairing
// non-monotonic fields given certain assumptions, but we do not make such
// assumptions here.
template<typename ScalarType, typename Device>
class FieldMonotonicityCheck: public FieldPropertyCheck<ScalarType, Device> {
public:

  // Default constructor.
  FieldMonotonicityCheck () {}

  // Overrides.

  bool check(const Field<ScalarType, Device>& field) const override {
    auto host_view = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_view, field.get_view());
    ScalarType sign;
    Kokkos::parallel_reduce(host_view.extent(0), KOKKOS_LAMBDA(Int i, ScalarType& s) {
      if ((i > 0) && (i < host_view.extent(0)-1)) {
        auto diff1 = host_view(i) - host_view(i-1);
        auto diff2 = host_view(i+1) - host_view(i);
        s *= (diff1 * diff2 > 0) ? (diff1 * diff2) : 0;
      } else {
        s *= 1;
      }
    }, Kokkos::Prod<ScalarType>(sign));
    return (sign > 0);
  }

  bool can_repair() const override {
    return false;
  }

  void repair(Field<ScalarType, Device>& field) const override {
    EKAT_REQUIRE_MSG(false, "Cannot repair a non-monotonic field!");
  }
};

} // namespace scream

#endif // SCREAM_FIELD_MONOTONICITY_CHECK_HPP
