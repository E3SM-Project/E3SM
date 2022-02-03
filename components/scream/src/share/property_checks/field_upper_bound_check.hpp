#ifndef SCREAM_FIELD_UPPER_BOUND_CHECK_HPP
#define SCREAM_FIELD_UPPER_BOUND_CHECK_HPP

#include "share/property_checks/field_within_interval_check.hpp"

#include <limits>

namespace scream
{

// Convenience implementation of check for interval (-\infty,U]. The class
// inherits from FieldWithinIntervalCheck, and sets lower bound to -infinity
class FieldUpperBoundCheck: public FieldWithinIntervalCheck {
public:
  // Constructor with upper bound. By default, this property check
  // can repair fields that fail the check by overwriting nonpositive values
  // with the given upper bound. If can_repair is false, the check cannot
  // apply repairs to the field.
  FieldUpperBoundCheck (const Field& field,
                        const double upper_bound,
                        const bool can_repair = true)
   : FieldWithinIntervalCheck(field,
                              -std::numeric_limits<double>::max(),
                              upper_bound,
                              can_repair)
  {
    // Do Nothing
  }

  // The name of the field check
  std::string name () const override {
    // NOTE: std::to_string does not do a good job with small numbers (like 1e-9).
    std::stringstream ss;
    ss << fields().front().name() << " upper bound check: " << this->m_upper_bound;
    return ss.str();
  }
};

} // namespace scream

#endif // SCREAM_FIELD_UPPER_BOUND_CHECK_HPP
