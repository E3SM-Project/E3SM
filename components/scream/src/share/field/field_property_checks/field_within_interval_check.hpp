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
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;

  // No default constructor -- we need lower and upper bounds.
  FieldWithinIntervalCheck () = delete;

  // Constructor with lower and upper bounds. By default, this property check
  // can repair fields that fail the check by overwriting nonpositive values
  // with the given lower bound. If can_repair is false, the check cannot
  // apply repairs to the field.
  explicit FieldWithinIntervalCheck (const_RT lower_bound,
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
  std::string name () const override { return "Within Interval Field Check"; }

  bool check(const Field<const_RT>& field) const override {
    using RT = non_const_RT;
    using minmax_t = typename Kokkos::MinMax<non_const_RT>::value_type;

    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto& dims = layout.dims();
    const int dim0 = dims[0];

    minmax_t minmax;
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<const_RT*>();
          Kokkos::parallel_reduce(dim0, KOKKOS_LAMBDA(int i, minmax_t& result) {
            result.min_val = ekat::impl::min(result.min_val, v(i));
            result.max_val = ekat::impl::max(result.max_val, v(i));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 2:
        {
          auto v = field.template get_view<const_RT**>();
          const int dim1 = dims[1];
          Kokkos::parallel_reduce(dim0*dim1, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            const int i = idx / dim1;
            const int j = idx % dim1;
            result.min_val = ekat::impl::min(result.min_val, v(i,j));
            result.max_val = ekat::impl::max(result.max_val, v(i,j));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 3:
        {
          auto v = field.template get_view<const_RT***>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          Kokkos::parallel_reduce(dim0*dim1*dim2, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            const int i = (idx / dim2) / dim1;
            const int j = (idx / dim2) % dim1;
            const int k =  idx % dim2;
            result.min_val = ekat::impl::min(result.min_val, v(i,j,k));
            result.max_val = ekat::impl::max(result.max_val, v(i,j,k));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 4:
        {
          auto v = field.template get_view<const_RT****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          Kokkos::parallel_reduce(dim0*dim1*dim2*dim3, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            const int i = ((idx / dim3) / dim2) / dim1;
            const int j = ((idx / dim3) / dim2) % dim1;
            const int k =  (idx / dim3) % dim2;
            const int l =   idx % dim3;
            result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l));
            result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l));
          }, Kokkos::MinMax<RT>(minmax));
        }
        break;
      case 5:
        {
          auto v = field.template get_view<const_RT*****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          const int dim4 = dims[4];
          Kokkos::parallel_reduce(dim0*dim1*dim2*dim3*dim4, KOKKOS_LAMBDA(int idx, minmax_t& result) {
            const int i = (((idx / dim4) / dim3) / dim2) / dim1;
            const int j = (((idx / dim4) / dim3) / dim2) % dim1;
            const int k =  ((idx / dim4) / dim3) % dim2;
            const int l =   (idx / dim4) % dim3;
            const int m =    idx % dim4;
            result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l,m));
            result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l,m));
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
    const auto& dims = layout.dims();
    const int dim0 = dims[0];

    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<non_const_RT*>();
          Kokkos::parallel_for(dim0, KOKKOS_LAMBDA(int i) {
            auto& ref = v(i);
            ref = ekat::impl::min(m_upper_bound, ref);
            ref = ekat::impl::max(m_lower_bound, ref);
          });
        }
        break;
      case 2:
        {
          auto v = field.template get_view<non_const_RT**>();
          const int dim1 = dims[1];
          Kokkos::parallel_for(dim0*dim1, KOKKOS_LAMBDA(int idx) {
            const int i = idx / dim1;
            const int j = idx % dim1;
            auto& ref = v(i,j);
            ref = ekat::impl::min(m_upper_bound, ref);
            ref = ekat::impl::max(m_lower_bound, ref);
          });
        }
        break;
      case 3:
        {
          auto v = field.template get_view<non_const_RT***>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          Kokkos::parallel_for(dim0*dim1*dim2, KOKKOS_LAMBDA(int idx) {
            const int i = (idx / dim2) / dim1;
            const int j = (idx / dim2) % dim1;
            const int k =  idx % dim2;
            auto& ref = v(i,j,k);
            ref = ekat::impl::min(m_upper_bound, ref);
            ref = ekat::impl::max(m_lower_bound, ref);
          });
        }
        break;
      case 4:
        {
          auto v = field.template get_view<non_const_RT****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          Kokkos::parallel_for(dim0*dim1*dim2*dim3, KOKKOS_LAMBDA(int idx) {
            const int i = ((idx / dim3) / dim2) / dim1;
            const int j = ((idx / dim3) / dim2) % dim1;
            const int k =  (idx / dim3) % dim2;
            const int l =   idx % dim3;
            auto& ref = v(i,j,k,l);
            ref = ekat::impl::min(m_upper_bound, ref);
            ref = ekat::impl::max(m_lower_bound, ref);
          });
        }
        break;
      case 5:
        {
          auto v = field.template get_view<non_const_RT*****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          const int dim4 = dims[4];
          Kokkos::parallel_for(dim0*dim1*dim2*dim3*dim4, KOKKOS_LAMBDA(int idx) {
            const int i = (((idx / dim4) / dim3) / dim2) / dim1;
            const int j = (((idx / dim4) / dim3) / dim2) % dim1;
            const int k =  ((idx / dim4) / dim3) % dim2;
            const int l =   (idx / dim4) % dim3;
            const int m =    idx % dim4;
            auto& ref = v(i,j,k,l,m);
            ref = ekat::impl::min(m_upper_bound, ref);
            ref = ekat::impl::max(m_lower_bound, ref);
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
