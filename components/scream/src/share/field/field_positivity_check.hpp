#ifndef SCREAM_FIELD_POSITIVITY_CHECK_HPP
#define SCREAM_FIELD_POSITIVITY_CHECK_HPP

#include "share/field/field.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// positive, and false if not. If constructed with a specific lower bound, it
// can repair a field that fails the check by setting all nonpositive values
// to that lower bound. If no specific lower bound is given (i.e. if the
// default constructor is used), no repairs can be made.
template<typename ScalarType, typename Device>
class FieldPositivityCheck: public FieldPropertyCheck<ScalarType, Device> {
public:

  // Default constructor -- cannot repair fields that fail the check.
  FieldPositivityCheck () : m_lower_bound(0) {}

  // Constructor with lower bound -- can repair fields that fail the check
  // by overwriting nonpositive values with the given lower bound.
  explicit FieldPositivityCheck (ScalarType lower_bound) :
    m_lower_bound(lower_bound) {
    EKAT_ASSERT_MSG(lower_bound > 0, "lower_bound must be positive.");
  }

  // Overrides.

  bool check(const Field<ScalarType, Device>& field) const override {
    auto f_view = field.get_view();
    auto host_view = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_view, field.get_view());
    ScalarType min_val = -90000000;//host_view(0);
    Kokkos::parallel_reduce(host_view.extent(0), KOKKOS_LAMBDA(Int i, ScalarType& m) {
      m = std::min(m, host_view(i));
    }, min_val);
    return (min_val > 0);
  }

  bool can_repair() const override {
    return (m_lower_bound > 0);
  }

  void repair(Field<ScalarType, Device>& field) const override {
    if (can_repair()) {
      auto host_view = Kokkos::create_mirror_view(field.get_view());
      Kokkos::deep_copy(host_view, field.get_view());
      Kokkos::parallel_for(host_view.extent(0), KOKKOS_LAMBDA(Int i) {
        if (host_view(i) < m_lower_bound) {
          host_view(i) = m_lower_bound;
        }
      });
      Kokkos::deep_copy(field.get_view(), host_view);
    }
  }

protected:

  // The given lower bound (0 if not supplied).
  ScalarType m_lower_bound;

};

} // namespace scream

#endif // SCREAM_FIELD_POSITIVITY_CHECK_HPP
