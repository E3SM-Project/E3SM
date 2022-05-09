#ifndef SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP
#define SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP

#include "share/property_checks/property_check.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/field/field.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// bounded by the given upper and lower bounds (inclusively), and false if not.
// It can repair a field that fails the check by clipping all unbounded values
// to their closest bound.

class FieldWithinIntervalCheck: public PropertyCheck {
public:

  // Constructor with lower and upper bounds. By default, this property check
  // can *NOT* repair fields that fail the check. If can_repair is true,
  // this class will overwrite values out of bounds with the proper bound
  // (upper if v>upper_bound and lower if v<lower_bound).
  FieldWithinIntervalCheck (const Field& field,
                            const std::shared_ptr<const AbstractGrid>& grid,
                            const double lower_bound,
                            const double upper_bound,
                            const bool can_repair = false);

  // The name of the property check
  std::string name () const override;

  CheckResult check() const override;

// CUDA requires the parent fcn of a KOKKOS_LAMBDA to have public access
#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  template<typename ST>
  CheckResult check_impl () const;

  template<typename ST>
  void repair_impl() const;

protected:

  void repair_impl() const override;

  // Lower and upper bounds.
  double m_lower_bound, m_upper_bound;

  std::shared_ptr<const AbstractGrid>   m_grid;
};

} // namespace scream

#endif // SCREAM_FIELD_BOUNDEDNESS_CHECK_HPP
