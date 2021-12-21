#ifndef SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP
#define SCREAM_FIELD_WITHIN_INTERVAL_CHECK_HPP

#include "share/field/field_property_check.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

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
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;

  // No default constructor -- we need lower and upper bounds.
  FieldWithinIntervalCheck () = delete;

  // Constructor with lower and upper bounds. By default, this property check
  // can repair fields that fail the check by overwriting nonpositive values
  // with the given lower bound. If can_repair is false, the check cannot
  // apply repairs to the field.
  FieldWithinIntervalCheck (const_RT lower_bound,
                            const_RT upper_bound,
                            bool can_repair = true) :
    m_lower_bound(lower_bound),
    m_upper_bound(upper_bound),
    m_can_repair(can_repair) {
    EKAT_ASSERT_MSG(lower_bound <= upper_bound,
                    "lower_bound must be less than or equal to upper_bound.");
  }

  // Overrides.

  // The name of the field check
  std::string name () const override { return "Within Interval [" + std::to_string(m_lower_bound) + "," + std::to_string(m_upper_bound) + "] Field Check"; }

  bool check(const Field<const_RT>& field) const override {
    using RT = non_const_RT;
    using minmax_t = typename Kokkos::MinMax<non_const_RT>::value_type;

    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto extents = layout.extents();
    const auto size = layout.size();

    minmax_t minmax;
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<const_RT*>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, minmax_t& result) {
            result.min_val = ekat::impl::min(result.min_val, v(i));
            result.max_val = ekat::impl::max(result.max_val, v(i));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 2:
        {
          auto v = field.template get_view<const_RT**>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            int i,j;
            unflatten_idx(idx,extents,i,j);
            result.min_val = ekat::impl::min(result.min_val, v(i,j));
            result.max_val = ekat::impl::max(result.max_val, v(i,j));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 3:
        {
          auto v = field.template get_view<const_RT***>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            int i,j,k;
            unflatten_idx(idx,extents,i,j,k);
            result.min_val = ekat::impl::min(result.min_val, v(i,j,k));
            result.max_val = ekat::impl::max(result.max_val, v(i,j,k));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 4:
        {
          auto v = field.template get_view<const_RT****>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            int i,j,k,l;
            unflatten_idx(idx,extents,i,j,k,l);
            result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l));
            result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 5:
        {
          auto v = field.template get_view<const_RT*****>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            int i,j,k,l,m;
            unflatten_idx(idx,extents,i,j,k,l,m);
            result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l,m));
            result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l,m));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 6:
        {
          auto v = field.template get_view<const_RT******>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            int i,j,k,l,m,n;
            unflatten_idx(idx,extents,i,j,k,l,m,n);
            result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l,m,n));
            result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l,m,n));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
    }
    return minmax.min_val>=m_lower_bound && minmax.max_val<=m_upper_bound;
  }

  bool can_repair() const override {
    return m_can_repair;
  }

  void repair(Field<non_const_RT>& field) const override {
    EKAT_REQUIRE_MSG (can_repair(),
        "Error! Cannot repair check '" + name() + "', for field '" + field.get_header().get_identifier().name() + "'.\n");

    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto extents = layout.extents();
    const auto size = layout.size();

    auto lb = m_lower_bound;
    auto ub = m_upper_bound;
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<non_const_RT*>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int i) {
            auto& ref = v(i);
            ref = ekat::impl::min(ub, ref);
            ref = ekat::impl::max(lb, ref);
          });
        }
        break;
      case 2:
        {
          auto v = field.template get_view<non_const_RT**>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j;
            unflatten_idx(idx,extents,i,j);

            auto& ref = v(i,j);
            ref = ekat::impl::min(ub, ref);
            ref = ekat::impl::max(lb, ref);
          });
        }
        break;
      case 3:
        {
          auto v = field.template get_view<non_const_RT***>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k;
            unflatten_idx(idx,extents,i,j,k);

            auto& ref = v(i,j,k);
            ref = ekat::impl::min(ub, ref);
            ref = ekat::impl::max(lb, ref);
          });
        }
        break;
      case 4:
        {
          auto v = field.template get_view<non_const_RT****>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l;
            unflatten_idx(idx,extents,i,j,k,l);

            auto& ref = v(i,j,k,l);
            ref = ekat::impl::min(ub, ref);
            ref = ekat::impl::max(lb, ref);
          });
        }
        break;
      case 5:
        {
          auto v = field.template get_view<non_const_RT*****>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m;
            unflatten_idx(idx,extents,i,j,k,l,m);

            auto& ref = v(i,j,k,l,m);
            ref = ekat::impl::min(ub, ref);
            ref = ekat::impl::max(lb, ref);
          });
        }
        break;
      case 6:
        {
          auto v = field.template get_view<non_const_RT******>();
          Kokkos::parallel_for(size, KOKKOS_LAMBDA(int idx) {
            int i,j,k,l,m,n;
            unflatten_idx(idx,extents,i,j,k,l,m,n);

            auto& ref = v(i,j,k,l,m,n);
            ref = ekat::impl::min(ub, ref);
            ref = ekat::impl::max(lb, ref);
          });
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
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
