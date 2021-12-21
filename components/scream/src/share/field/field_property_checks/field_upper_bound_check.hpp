#ifndef SCREAM_FIELD_UPPER_BOUND_CHECK_HPP
#define SCREAM_FIELD_UPPER_BOUND_CHECK_HPP

#include "share/field/field_property_checks/field_within_interval_check.hpp"

#include <limits>

namespace scream
{

// This field property check returns true if all data in a given field are
// bounded above by the given upper bound, and false if not.
// It can repair a field that fails the check by clipping all unbounded values
// to the upper bound.
// Note: The upper bound check is a sub-class of the within interval check with
// the upper bound set to the numeric limit.
template<typename RealType>
class FieldUpperBoundCheck: public FieldWithinIntervalCheck<RealType> {
public:
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;

  // No default constructor -- we need an upper bound.
  FieldUpperBoundCheck () = delete;

  // Constructor with upper bound. By default, this property check
  // can repair fields that fail the check by overwriting nonpositive values
  // with the given upper bound. If can_repair is false, the check cannot
  // apply repairs to the field.
  explicit FieldUpperBoundCheck (const_RT upper_bound,
                                 bool can_repair = true) : 
    FieldWithinIntervalCheck<RealType>(-std::numeric_limits<RealType>::max(), upper_bound, can_repair)
    {
      // Do Nothing
    }

  // Overrides.

  // The name of the field check
  std::string name () const override { return "Upper Bound Check of " + std::to_string(this->m_upper_bound); }

};

} // namespace scream

#endif // SCREAM_FIELD_UPPER_BOUND_CHECK_HPP
