#ifndef SCREAM_FIELD_POSITIVITY_CHECK_HPP
#define SCREAM_FIELD_POSITIVITY_CHECK_HPP

#include "share/field/field.hpp"
#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// positive, and false if not. If constructed with a specific lower bound, it
// can repair a field that fails the check by setting all nonpositive values
// to that lower bound. If no specific lower bound is given (i.e. if the
// default constructor is used), no repairs can be made.
template<typename Realtype>
class FieldPositivityCheck: public FieldPropertyCheck<Realtype> {
public:

  // Default constructor -- cannot repair fields that fail the check.
  FieldPositivityCheck () : m_lower_bound(0) {}

  // Constructor with lower bound -- can repair fields that fail the check
  // by overwriting nonpositive values with the given lower bound.
  explicit FieldPositivityCheck (Realtype lower_bound) :
    m_lower_bound(lower_bound) {
    EKAT_ASSERT_MSG(lower_bound > 0, "lower_bound must be positive.");
  }

  // Overrides.

  bool check(const Field<Realtype>& field) const override {
    auto view = field.get_view();
    Realtype min_val;
    Kokkos::parallel_reduce(view.extent(0), KOKKOS_LAMBDA(Int i, Realtype& m) {
      if (i == 0) {
        m = view(0);
      } else {
        m = ekat::impl::min(m, view(i));
      }
    }, Kokkos::Min<Realtype>(min_val));
    return (min_val > 0);
  }

  bool can_repair() const override {
    return (m_lower_bound > 0);
  }

  void repair(Field<Realtype>& field) const override {
    if (can_repair()) {
      auto view = field.get_view();
      Kokkos::parallel_for(view.extent(0), KOKKOS_LAMBDA(Int i) {
        if (view(i) < 0) {
          view(i) = m_lower_bound;
        }
      });
    }
  }

protected:

  // The given lower bound (0 if not supplied).
  Realtype m_lower_bound;

};

} // namespace scream

#endif // SCREAM_FIELD_POSITIVITY_CHECK_HPP
