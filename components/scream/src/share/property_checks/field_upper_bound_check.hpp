#ifndef SCREAM_FIELD_UPPER_BOUND_CHECK_HPP
#define SCREAM_FIELD_UPPER_BOUND_CHECK_HPP

#include "share/property_checks/field_within_interval_check.hpp"

#include <limits>

namespace scream
{

// Convenience implementation of check for interval (-\infty,U]. The class
// inherits from FieldWithinIntervalCheck, and sets lower bound to -DBL_MAX
class FieldUpperBoundCheck: public FieldWithinIntervalCheck {
public:
  // Constructor with lower bound. By default, this property check
  // can *NOT* repair fields that fail the check. If can_repair is true,
  // this class will overwrite values out of bounds with the stored upper bound
  FieldUpperBoundCheck (const Field& field,
                        const std::shared_ptr<const AbstractGrid>& grid,
                        const double upper_bound,
                        const bool can_repair = false,
                        const double ub_repairable = s_max)
   : FieldWithinIntervalCheck (field, grid,
                              -s_max,
                               upper_bound,
                               can_repair,
                              -s_max,
                               ub_repairable)
  {
    // Do Nothing
  }

  // The name of the field check
  std::string name () const override {
    // NOTE: std::to_string does not do a good job with small numbers (like 1e-9).
    std::stringstream ss;
    ss << fields().front().name() << " upper bound check: " << this->m_ub;
    return ss.str();
  }
};

} // namespace scream

#endif // SCREAM_FIELD_UPPER_BOUND_CHECK_HPP
