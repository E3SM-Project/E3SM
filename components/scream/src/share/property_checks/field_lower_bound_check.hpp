#ifndef SCREAM_FIELD_LOWER_BOUND_CHECK_HPP
#define SCREAM_FIELD_LOWER_BOUND_CHECK_HPP

#include "share/property_checks/field_within_interval_check.hpp"

#include <limits>

namespace scream
{

// Convenience implementation of check for interval [U,-\infty). The class
// inherits from FieldWithinIntervalCheck, and sets lower bound to -infinity
class FieldLowerBoundCheck: public FieldWithinIntervalCheck {
public:
  // Constructor with lower bound. By default, this property check
  // can *NOT* repair fields that fail the check. If can_repair is true,
  // this class will overwrite values out of bounds with the stored lower bound
  FieldLowerBoundCheck (const Field& field,
                        const std::shared_ptr<const AbstractGrid>& grid,
                        const double lower_bound,
                        const bool can_repair = false)
   : FieldWithinIntervalCheck(field, grid,
                              lower_bound,
                              std::numeric_limits<double>::max(),
                              can_repair)
  {
    // Do Nothing
  }

  // The name of the field check
  std::string name () const override {
    // NOTE: std::to_string does not do a good job with small numbers (like 1e-9).
    std::stringstream ss;
    ss << fields().front().name() << " lower bound check: " << this->m_lower_bound;
    return ss.str();
  }
};

} // namespace scream

#endif // SCREAM_FIELD_LOWER_BOUND_CHECK_HPP
