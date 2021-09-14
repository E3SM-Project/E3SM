#ifndef SCREAM_FIELD_LOWER_BOUND_CHECK_HPP
#define SCREAM_FIELD_LOWER_BOUND_CHECK_HPP

#include "share/field/field_property_checks/field_within_interval_check.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// bounded below by the given lower bounds, and false if not.
// It can repair a field that fails the check by clipping all unbounded values
// to the lower bound.
// Note: The lower bound check is a sub-class of the within interval check with
// the upper bound set to the numeric limit.
template<typename RealType>
class FieldLowerBoundCheck: public FieldWithinIntervalCheck<RealType> {
public:
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;

  // No default constructor -- we need a lower bound.
  FieldLowerBoundCheck () = delete;

  // Constructor with lower and bound. By default, this property check
  // can repair fields that fail the check by overwriting nonpositive values
  // with the given lower bound. If can_repair is false, the check cannot
  // apply repairs to the field.
  explicit FieldLowerBoundCheck (const_RT lower_bound,
                                 bool can_repair = true) : 
    FieldWithinIntervalCheck<RealType>(lower_bound, std::numeric_limits<RealType>::max(), can_repair),
    m_lower_bound(lower_bound)
    {
      // Do Nothing
    }

  // Overrides.

  // The name of the field check
  std::string name () const override { return "Lower Bound Check of " + std::to_string(m_lower_bound); }

  protected:

  // Lower bound
  RealType m_lower_bound;

};

} // namespace scream

#endif // SCREAM_FIELD_LOWER_BOUND_CHECK_HPP
