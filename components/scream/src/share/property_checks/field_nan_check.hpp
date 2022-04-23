#ifndef SCREAM_FIELD_NAN_CHECK_HPP
#define SCREAM_FIELD_NAN_CHECK_HPP

#include "share/property_checks/property_check.hpp"

namespace scream
{

// Inspect whether a field contains any NaN values.
// No repair allowed. If we find NaN's, we should crash.
class FieldNaNCheck: public PropertyCheck {
public:
  FieldNaNCheck (const Field& f);

  // The name of the field check
  std::string name () const override {
    return "NaN check for field " + fields().front().name();
  }

  CheckResult check() const override;

// CUDA requires the parent fcn of a KOKKOS_LAMBDA to have public access
#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif
  template<typename ST>
  CheckResult check_impl() const;
};

} // namespace scream

#endif // SCREAM_FIELD_NAN_CHECK_HPP
