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
template<typename RealType>
class FieldPositivityCheck: public FieldPropertyCheck<RealType> {
public:
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;
  using kt = typename Field<RealType>::kt;
  using exe_space = typename kt::ExeSpace;
  using RangePolicy = typename kt::RangePolicy;

  // Default constructor -- cannot repair fields that fail the check.
  FieldPositivityCheck () : m_lower_bound(0) {}

  // Constructor with lower bound -- can repair fields that fail the check
  // by overwriting nonpositive values with the given lower bound.
  explicit FieldPositivityCheck (const_RT lower_bound) :
    m_lower_bound(lower_bound) {
    EKAT_ASSERT_MSG(lower_bound > 0, "lower_bound must be positive.");
  }

  // Overrides.

  // The name of the field check
  std::string name () const override { return "Positivity Field Check"; }

  bool check(const Field<const_RT>& field) const override {
    using RT = non_const_RT;
    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto& dims = layout.dims();
    const int dim0 = dims[0];

    RT min_val = std::numeric_limits<RT>::max();
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<const_RT*>();
          Kokkos::parallel_reduce(dim0, KOKKOS_LAMBDA(int i, RT& result) {
            result = ekat::impl::min(result, v(i));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 2:
        {
          auto v = field.template get_view<const_RT**>();
          const int dim1 = dims[1];
          Kokkos::parallel_reduce(dim0*dim1, KOKKOS_LAMBDA(int idx, RT& result) {
            const int i = idx / dim1;
            const int j = idx % dim1;
            result = ekat::impl::min(result, v(i,j));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 3:
        {
          auto v = field.template get_view<const_RT***>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          Kokkos::parallel_reduce(dim0*dim1*dim2, KOKKOS_LAMBDA(int idx, RT& result) {
            const int i = (idx / dim2) / dim1;
            const int j = (idx / dim2) % dim1;
            const int k =  idx % dim2;
            result = ekat::impl::min(result, v(i,j,k));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 4:
        {
          auto v = field.template get_view<const_RT****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          Kokkos::parallel_reduce(dim0*dim1*dim2*dim3, KOKKOS_LAMBDA(int idx, RT& result) {
            const int i = ((idx / dim3) / dim2) / dim1;
            const int j = ((idx / dim3) / dim2) % dim1;
            const int k =  (idx / dim3) % dim2;
            const int l =   idx % dim3;
            result = ekat::impl::min(result, v(i,j,k,l));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      case 5:
        {
          auto v = field.template get_view<const_RT*****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          const int dim4 = dims[4];
          Kokkos::parallel_reduce(dim0*dim1*dim2*dim3*dim4, KOKKOS_LAMBDA(int idx, RT& result) {
            const int i = (((idx / dim4) / dim3) / dim2) / dim1;
            const int j = (((idx / dim4) / dim3) / dim2) % dim1;
            const int k =  ((idx / dim4) / dim3) % dim2;
            const int l =   (idx / dim4) % dim3;
            const int m =    idx % dim4;
            result = ekat::impl::min(result, v(i,j,k,l,m));
          }, Kokkos::Min<RT>(min_val));
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
    }
    return min_val>0;
  }

  bool can_repair() const override {
    return (m_lower_bound > 0);
  }

  void repair(Field<non_const_RT>& field) const override {
    EKAT_REQUIRE_MSG (can_repair(),
        "Error! Cannot repair check '" + name() + "', for field '" + field.get_header().get_identifier().name() + "'.\n");
    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto& dims = layout.dims();
    const int dim0 = dims[0];

    const auto lb = m_lower_bound;
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<non_const_RT*>();
          Kokkos::parallel_for(dim0, KOKKOS_LAMBDA(int i) {
            v(i) = ekat::impl::max(lb, v(i));
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
            v(i,j) = ekat::impl::max(lb, v(i,j));
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
            v(i,j,k) = ekat::impl::max(lb, v(i,j,k));
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
            v(i,j,k,l) = ekat::impl::max(lb, v(i,j,k,l));
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
            v(i,j,k,l,m) = ekat::impl::max(lb, v(i,j,k,l,m));
          });
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
    }
  }

protected:

  // The given lower bound (0 if not supplied).
  RealType m_lower_bound;
};

} // namespace scream

#endif // SCREAM_FIELD_POSITIVITY_CHECK_HPP
