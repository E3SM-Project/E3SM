#ifndef SCREAM_FIELD_POSITIVITY_CHECK_HPP
#define SCREAM_FIELD_POSITIVITY_CHECK_HPP

#include "share/property_checks/field_lower_bound_check.hpp"

namespace scream
{

// Convenience implementation of lower bound check for L=0. The class
// inherits from FieldLowerBound, and sets L to machine epsilon

class FieldPositivityCheck: public FieldLowerBoundCheck {
public:

  // Default constructor -- cannot repair fields that fail the check.
  FieldPositivityCheck (const Field& f,
                        const std::shared_ptr<const AbstractGrid>& grid,
                        const bool can_repair = false)
   : FieldLowerBoundCheck (f,grid,
                           std::numeric_limits<double>::epsilon(),
                           can_repair)
  {
    // Nothing to do here
  }

  // The name of the field check
  std::string name () const override {
    return "Positivity check for field " + fields().front().name();
  }
};

} // namespace scream

#endif // SCREAM_FIELD_POSITIVITY_CHECK_HPP
