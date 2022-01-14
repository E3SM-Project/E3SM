#ifndef SCREAM_FIELD_POSITIVITY_CHECK_HPP
#define SCREAM_FIELD_POSITIVITY_CHECK_HPP

#include "share/field/field_property_check.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// positive, and false if not. If constructed with a specific lower bound, it
// can repair a field that fails the check by setting all nonpositive values
// to that lower bound. If no specific lower bound is given (i.e. if the
// default constructor is used), no repairs can be made.
template<typename RealType>
class FieldPositivityCheck: public FieldPropertyCheck<RealType> {
public:
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;
  using kt = typename Field<RealType>::kt;
  using exe_space = typename kt::ExeSpace;
  using RangePolicy = typename kt::RangePolicy;

  // Default constructor -- cannot repair fields that fail the check.
  FieldPositivityCheck () : m_lower_bound(-1) {}

  // Constructor with lower bound -- can repair fields that fail the check
  // by overwriting nonpositive values with the given lower bound.
  explicit FieldPositivityCheck (const_RT lower_bound) :
    m_lower_bound(lower_bound) {
    EKAT_ASSERT_MSG(lower_bound >= 0, "lower_bound must be positive.");
  }

  // Overrides.

  // The name of the field check
  std::string name () const override { return "Positivity Field Check"; }

  bool check(const Field<const_RT>& field) const override {
    using RT = non_const_RT;
    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto extents = layout.extents();
    const auto size = layout.size();

    RT min_val = std::numeric_limits<RT>::max();
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<const_RT*>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, RT& result) {
            result = ekat::impl::min(result, v(i));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 2:
        {
          auto v = field.template get_view<const_RT**>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, RT& result) {
            int i,j;
            unflatten_idx(idx,extents,i,j);
            result = ekat::impl::min(result, v(i,j));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 3:
        {
          auto v = field.template get_view<const_RT***>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, RT& result) {
            int i,j,k;
            unflatten_idx(idx,extents,i,j,k);
            result = ekat::impl::min(result, v(i,j,k));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 4:
        {
          auto v = field.template get_view<const_RT****>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, RT& result) {
            int i,j,k,l;
            unflatten_idx(idx,extents,i,j,k,l);
            result = ekat::impl::min(result, v(i,j,k,l));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 5:
        {
          auto v = field.template get_view<const_RT*****>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, RT& result) {
            int i,j,k,l,m;
            unflatten_idx(idx,extents,i,j,k,l,m);
            result = ekat::impl::min(result, v(i,j,k,l,m));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 6:
        {
          auto v = field.template get_view<const_RT******>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, RT& result) {
            int i,j,k,l,m,n;
            unflatten_idx(idx,extents,i,j,k,l,m,n);
            result = ekat::impl::min(result, v(i,j,k,l,m,n));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
    }
    return min_val>=0;
  }

  bool can_repair() const override {
    return (m_lower_bound >= 0);
  }

  void repair(Field<non_const_RT>& field) const override {
    EKAT_REQUIRE_MSG (can_repair(),
        "Error! Cannot repair check '" + name() + "', for field '" + field.get_header().get_identifier().name() + "'.\n");
    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto extents = layout.extents();
    const auto size = layout.size();

    const auto lb = m_lower_bound;
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<non_const_RT*>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int i) {
            v(i) = ekat::impl::max(lb, v(i));
          });
        }
        break;
      case 2:
        {
          auto v = field.template get_view<non_const_RT**>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j;
            unflatten_idx(idx,extents,i,j);
            v(i,j) = ekat::impl::max(lb, v(i,j));
          });
        }
        break;
      case 3:
        {
          auto v = field.template get_view<non_const_RT***>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k;
            unflatten_idx(idx,extents,i,j,k);
            v(i,j,k) = ekat::impl::max(lb, v(i,j,k));
          });
        }
        break;
      case 4:
        {
          auto v = field.template get_view<non_const_RT****>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l;
            unflatten_idx(idx,extents,i,j,k,l);
            v(i,j,k,l) = ekat::impl::max(lb, v(i,j,k,l));
          });
        }
        break;
      case 5:
        {
          auto v = field.template get_view<non_const_RT*****>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m;
            unflatten_idx(idx,extents,i,j,k,l,m);
            v(i,j,k,l,m) = ekat::impl::max(lb, v(i,j,k,l,m));
          });
        }
        break;
      case 6:
        {
          auto v = field.template get_view<non_const_RT******>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m,n;
            unflatten_idx(idx,extents,i,j,k,l,m,n);
            v(i,j,k,l,m,n) = ekat::impl::max(lb, v(i,j,k,l,m,n));
          });
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
    }
  }

protected:

  // The given lower bound (-1 if not supplied, and if not supplied won't be used for repair).
  RealType m_lower_bound;
};

} // namespace scream

#endif // SCREAM_FIELD_POSITIVITY_CHECK_HPP
