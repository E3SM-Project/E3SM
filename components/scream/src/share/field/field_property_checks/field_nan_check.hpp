#ifndef SCREAM_FIELD_NAN_CHECK_HPP
#define SCREAM_FIELD_NAN_CHECK_HPP

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// not a NaN.
// Default to no repair allowed, we don't really want to be repairing NaN's
// in our fields.
template<typename RealType>
class FieldNaNCheck: public FieldPropertyCheck<RealType> {
public:
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;

  // Default constructor -- cannot repair fields that fail the check.
  FieldNaNCheck () = default; 

  // Overrides.

  // The name of the field check
  std::string name () const override { return "NaN Field Check"; }

  bool check(const Field<const_RT>& field) const override {
    const auto& layout = field.get_header().get_identifier().get_layout();
    const auto extents = layout.extents();
    const auto size = layout.size();

    int num_nans = 0;
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<const_RT*>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, int& result) {
            if (std::isnan(v(i))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 2:
        {
          auto v = field.template get_view<const_RT**>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
            int i,j;
            unflatten_idx(idx,extents,i,j);
            if (std::isnan(v(i,j))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 3:
        {
          auto v = field.template get_view<const_RT***>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
            int i,j,k;
            unflatten_idx(idx,extents,i,j,k);
            if (std::isnan(v(i,j,k))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 4:
        {
          auto v = field.template get_view<const_RT****>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
            int i,j,k,l;
            unflatten_idx(idx,extents,i,j,k,l);
            if (std::isnan(v(i,j,k,l))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 5:
        {
          auto v = field.template get_view<const_RT*****>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
            int i,j,k,l,m;
            unflatten_idx(idx,extents,i,j,k,l,m);
            if (std::isnan(v(i,j,k,l,m))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 6:
        {
          auto v = field.template get_view<const_RT******>();
          Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
            int i,j,k,l,m,n;
            unflatten_idx(idx,extents,i,j,k,l,m,n);
            if (isnan(v(i,j,k,l,m,n))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
    }
    return num_nans==0;
  }

  bool can_repair() const override {
    return false;
  }

  void repair(Field<non_const_RT>& /* field */) const override {
    // Do Nothing
  }

protected:

};

} // namespace scream

#endif // SCREAM_FIELD_NEN_CHECK_HPP
