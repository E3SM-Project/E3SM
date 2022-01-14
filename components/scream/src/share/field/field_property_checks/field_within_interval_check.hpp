#ifndef SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP
#define SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP

#include <type_traits>
#include "share/field/field_property_check.hpp"
#include "share/field/field.hpp"
#include "share/util/scream_view_utils.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// bounded by the given upper and lower bounds (inclusively), and false if not.
// It can repair a field that fails the check by clipping all unbounded values
// to their closest bound.
class FieldWithinIntervalCheck: public FieldPropertyCheck {
public:
  // No default constructor -- we need lower and upper bounds.
  FieldWithinIntervalCheck () = delete;

  // Constructor with lower and upper bounds. By default, this property check
  // can repair fields that fail the check by overwriting nonpositive values
  // with the given lower bound. If can_repair is false, the check cannot
  // apply repairs to the field.
  FieldWithinIntervalCheck (const double lower_bound,
                            const double upper_bound,
                            bool can_repair = true) :
    m_lower_bound(lower_bound),
    m_upper_bound(upper_bound),
    m_can_repair(can_repair) {
    EKAT_ASSERT_MSG(lower_bound <= upper_bound,
                    "lower_bound must be less than or equal to upper_bound.");
  }

  // Overrides.

  // The name of the field check
  std::string name () const override {
    // NOTE: std::to_string does not do a good job with small numbers (like 1e-9).
    std::stringstream ss;
    ss << "Within Interval [" << m_lower_bound << ", " << m_upper_bound << "] Check";
    return ss.str();
  }

  bool check(const Field& field) const override;

  bool can_repair() const override {
    return m_can_repair;
  }

  void repair(Field& field) const override;

protected:

  template<typename ST>
  bool check_impl (const Field& field) const;

  template<typename ST>
  void repair_impl(Field& field) const;

  // Lower and upper bounds.
  double m_lower_bound, m_upper_bound;

  // Can we repair a field?
  bool m_can_repair;
};

} // namespace scream

#endif // SCREAM_FIELD_BOUNDEDNESS_CHECK_HPP
