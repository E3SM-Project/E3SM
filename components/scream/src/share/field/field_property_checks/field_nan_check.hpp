#ifndef SCREAM_FIELD_NAN_CHECK_HPP
#define SCREAM_FIELD_NAN_CHECK_HPP

#include "share/field/field.hpp"
#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// not a NaN.
// Default to no repair allowed, we don't really want to be repairing NaN's
// in our fields.
template<typename RealType>
class FieldNaNCheck: public FieldPropertyCheck<RealType> {
public:
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;

  // Default constructor -- cannot repair fields that fail the check.
  FieldNaNCheck () = default; 

  // Overrides.

  // The name of the field check
  std::string name () const override { return "NaN Field Check"; }

  bool check(const Field<const_RT>& field) const override {
    auto view = field.get_flattened_view();
    Int num_nans = 0; 
    Kokkos::parallel_reduce(view.extent(0), KOKKOS_LAMBDA(Int i, Int& m) {
      if (isnan(view(i))) {
        m += 1;
      }
    }, num_nans);
    return (num_nans == 0);
  }

  bool can_repair() const override {
    return false;
  }

  void repair(Field<non_const_RT>& /* field */) const override {
    // Do Nothing
  }

protected:

};

} // namespace scream

#endif // SCREAM_FIELD_NEN_CHECK_HPP
