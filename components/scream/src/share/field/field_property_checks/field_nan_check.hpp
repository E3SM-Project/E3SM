#ifndef SCREAM_FIELD_NAN_CHECK_HPP
#define SCREAM_FIELD_NAN_CHECK_HPP

#include "share/field/field_property_check.hpp"

namespace scream
{

// Inspect whether a field contains any NaN values.
// NOTE: we do NOT allow to repair the field.

class FieldNaNCheck: public FieldPropertyCheck {
public:
  // Default constructor -- cannot repair fields that fail the check.
  FieldNaNCheck () = default; 

  // The name of the field check
  std::string name () const override { return "NaN Field Check"; }

  bool check(const Field& field) const override;

  bool can_repair() const override {
    return false;
  }

  void repair(Field& /* field */) const override {
    EKAT_ERROR_MSG ("Error! FieldNaNCheck cannot repair the field.\n");
  }

protected:
  template<typename ST>
  bool check_impl(const Field& field) const;
};

} // namespace scream

#endif // SCREAM_FIELD_NAN_CHECK_HPP
