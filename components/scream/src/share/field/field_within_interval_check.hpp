#ifndef SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP
#define SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP

#include "share/field/field.hpp"
#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// bounded by the given upper and lower bounds (inclusively), and false if not.
// It can repair a field that fails the check by clipping all unbounded values
// to their closest bound.
template<typename RealType>
class FieldWithinIntervalCheck: public FieldPropertyCheck<RealType> {
public:

  // No default constructor -- we need lower and upper bounds.
  FieldWithinIntervalCheck () = delete;

  // Constructor with lower and upper bounds. By default, this property check
  // can repair fields that fail the check by overwriting nonpositive values
  // with the given lower bound. If can_repair is false, the check cannot
  // apply repairs to the field.
  explicit FieldWithinIntervalCheck (RealType lower_bound,
                                     RealType upper_bound,
                                     bool can_repair = true) :
    m_lower_bound(lower_bound),
    m_upper_bound(upper_bound),
    m_can_repair(can_repair) {
    EKAT_ASSERT_MSG(lower_bound <= upper_bound,
                    "lower_bound must be less than or equal to upper_bound.");
  }

  // Overrides.

  bool check(const Field<RealType>& field) const override {
    auto view = field.get_view();
    typename Kokkos::MinMax<RealType>::value_type minmax;
    Kokkos::parallel_reduce(view.extent(0), KOKKOS_LAMBDA(Int i,
        typename Kokkos::MinMax<RealType>::value_type& mm) {
      if (i == 0) {
        mm.min_val = mm.max_val = view(0);
      } else {
        mm.min_val = ekat::impl::min(mm.min_val, view(i));
        mm.max_val = ekat::impl::max(mm.max_val, view(i));
      }
    }, Kokkos::MinMax<RealType>(minmax));
    return ((minmax.min_val >= m_lower_bound) && (minmax.max_val <= m_upper_bound));
  }

  bool can_repair() const override {
    return m_can_repair;
  }

  void repair(Field<RealType>& field) const override {
    if (m_can_repair) {
      auto view = field.get_view();
      RealType lower_bound = m_lower_bound;
      RealType upper_bound = m_upper_bound;
      Kokkos::parallel_for(view.extent(0), KOKKOS_LAMBDA(Int i) {
        auto fi = view(i);
        if (fi < lower_bound) {
          view(i) = lower_bound;
        } else if (fi > upper_bound) {
          view(i) = upper_bound;
        }
      });
      Kokkos::fence();
    } else {
      EKAT_REQUIRE_MSG(false, "Cannot repair the field!");
    }
  }

protected:

  // Lower and upper bounds.
  RealType m_lower_bound, m_upper_bound;

  // Can we repair a field?
  bool m_can_repair;

};

} // namespace scream

#endif // SCREAM_FIELD_BOUNDEDNESS_CHECK_HPP
