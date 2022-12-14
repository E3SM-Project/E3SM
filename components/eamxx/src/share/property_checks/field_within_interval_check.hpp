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
  static constexpr double s_max = std::numeric_limits<double>::max();

  // Constructor with lower and upper bounds. By default, this property check
  // can *NOT* repair fields that fail the check. If can_repair is true,
  // this class will overwrite values out of bounds with the proper bound
  // (upper if v>upper_bound and lower if v<lower_bound).
  // The [lb|ub]_repairable bounds are no tighter than [lb,ub] (if not set,
  // they are set equal to lb/ub). If field is outside [lb,ub], but inside
  // [lb_rep,ub_rep], we return a "Repairable" check result, rather than a Fail.
  FieldWithinIntervalCheck (const Field& field,
                            const std::shared_ptr<const AbstractGrid>& grid,
                            const double lower_bound,
                            const double upper_bound,
                            const bool can_repair = false,
                            const double lb_repairable = -s_max,
                            const double ub_repairable =  s_max);

  // The name of the property check
  std::string name () const override;

  PropertyType type () const override { return PropertyType::PointWise; }

  ResultAndMsg check() const override;

// CUDA requires the parent fcn of a KOKKOS_LAMBDA to have public access
#ifndef EAMXX_ENABLE_GPU
protected:
#endif

  template<typename ST>
  ResultAndMsg check_impl () const;

  template<typename ST>
  void repair_impl() const;

  bool same_as (const PropertyCheck& pc) const override;

protected:

  void repair_impl() const override;

  // Lower and upper bounds.
  double m_lb, m_ub;

  // (Potentially) Less tight bounds
  double m_lb_repairable, m_ub_repairable;

  std::shared_ptr<const AbstractGrid>   m_grid;
};

} // namespace scream

#endif // SCREAM_FIELD_BOUNDEDNESS_CHECK_HPP
