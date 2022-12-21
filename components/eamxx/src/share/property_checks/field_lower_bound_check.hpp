#ifndef SCREAM_FIELD_LOWER_BOUND_CHECK_HPP
#define SCREAM_FIELD_LOWER_BOUND_CHECK_HPP

#include "share/property_checks/field_within_interval_check.hpp"

#include <limits>

namespace scream
{

// Convenience implementation of check for interval [L,-\infty). The class
// inherits from FieldWithinIntervalCheck, and sets upper bound to DBL_MAX
class FieldLowerBoundCheck: public FieldWithinIntervalCheck {
public:
  // Constructor with lower bound. By default, this property check
  // can *NOT* repair fields that fail the check. If can_repair is true,
  // this class will overwrite values out of bounds with the stored lower bound
  FieldLowerBoundCheck (const Field& field,
                        const std::shared_ptr<const AbstractGrid>& grid,
                        const double lower_bound,
                        const bool can_repair = false,
                        const double lb_repairable = -s_max)
   : FieldWithinIntervalCheck(field, grid,
                              lower_bound,
                              s_max,
                              can_repair,
                              lb_repairable,
                              s_max)
  {
    // Do Nothing
  }

  // The name of the field check
  std::string name () const override {
    // NOTE: std::to_string does not do a good job with small numbers (like 1e-9).
    std::stringstream ss;
    ss << fields().front().name() << " lower bound check: " << this->m_lb;
    return ss.str();
  }
};

} // namespace scream

#endif // SCREAM_FIELD_LOWER_BOUND_CHECK_HPP
