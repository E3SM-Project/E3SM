#ifndef SCREAM_FIELD_BOUNDEDNESS_CHECK_HPP
#define SCREAM_FIELD_BOUNDEDNESS_CHECK_HPP

#include "share/field/field.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// bounded by the given upper and lower bounds (inclusively), and false if not.
// It can repair a field that fails the check by clipping all unbounded values
// to their closest bound.
template<typename ScalarType, typename Device>
class FieldBoundednessCheck: public FieldPropertyCheck<ScalarType, Device> {
public:

  // No default constructor -- we need lower and upper bounds.
  FieldBoundednessCheck () = delete;

  // Constructor with lower and upper bounds -- can repair fields that fail the check
  // by overwriting nonpositive values with the given lower bound.
  explicit FieldBoundednessCheck (ScalarType lower_bound,
                                  ScalarType upper_bound) :
    m_lower_bound(lower_bound),
    m_upper_bound(upper_bound) {
    EKAT_ASSERT_MSG(lower_bound < upper_bound,
                    "lower_bound must be less than upper_bound.");
  }

  // Overrides.

  bool check(const Field<ScalarType, Device>& field) const override {
    auto host_view = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_view, field.get_view());
    bool bounded = true;
    Kokkos::parallel_reduce(host_view.extent(0), KOKKOS_LAMBDA(Int i, bool& b) {
      b = b && ((host_view(i) >= m_lower_bound) && ((host_view(i) <= m_upper_bound)));
    }, bounded);
    return bounded;
  }

  bool can_repair() const override {
    return true;
  }

  void repair(Field<ScalarType, Device>& field) const override {
    auto host_view = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_view, field.get_view());
    Kokkos::parallel_for(host_view.extent(0), KOKKOS_LAMBDA(Int i) {
      auto fi = host_view(i);
      if (fi < m_lower_bound) {
        host_view(i) = m_lower_bound;
      } else if (fi > m_upper_bound) {
        host_view(i) = m_upper_bound;
      }
    });
    Kokkos::deep_copy(field.get_view(), host_view);
  }

protected:

  // Lower and upper bounds.
  ScalarType m_lower_bound, m_upper_bound;

};

} // namespace scream

#endif // SCREAM_FIELD_BOUNDEDNESS_CHECK_HPP
