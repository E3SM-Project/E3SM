#ifndef SCREAM_FIELD_NAN_CHECK_HPP
#define SCREAM_FIELD_NAN_CHECK_HPP

#include "share/field/field.hpp"

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
    const auto& dims = layout.dims();
    const int dim0 = dims[0];

    int num_nans = 0;
    switch (layout.rank()) {
      case 1:
        {
          auto v = field.template get_view<const_RT*>();
          Kokkos::parallel_reduce(dim0, KOKKOS_LAMBDA(int i, int& result) {
            if (isnan(v(i))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 2:
        {
          auto v = field.template get_view<const_RT**>();
          const int dim1 = dims[1];
          Kokkos::parallel_reduce(dim0*dim1, KOKKOS_LAMBDA(int idx, int& result) {
            const int i = idx / dim1;
            const int j = idx % dim1;
            if (isnan(v(i,j))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 3:
        {
          auto v = field.template get_view<const_RT***>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          Kokkos::parallel_reduce(dim0*dim1*dim2, KOKKOS_LAMBDA(int idx, int& result) {
            const int i = (idx / dim2) / dim1;
            const int j = (idx / dim2) % dim1;
            const int k =  idx % dim2;
            if (isnan(v(i,j,k))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 4:
        {
          auto v = field.template get_view<const_RT****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          Kokkos::parallel_reduce(dim0*dim1*dim2*dim3, KOKKOS_LAMBDA(int idx, int& result) {
            const int i = ((idx / dim3) / dim2) / dim1;
            const int j = ((idx / dim3) / dim2) % dim1;
            const int k =  (idx / dim3) % dim2;
            const int l =   idx % dim3;
            if (isnan(v(i,j,k,l))) {
              ++result;
            }
          }, Kokkos::Sum<int>(num_nans));
        }
        break;
      case 5:
        {
          auto v = field.template get_view<const_RT*****>();
          const int dim1 = dims[1];
          const int dim2 = dims[2];
          const int dim3 = dims[3];
          const int dim4 = dims[4];
          Kokkos::parallel_reduce(dim0*dim1*dim2*dim3*dim4, KOKKOS_LAMBDA(int idx, int& result) {
            const int i = (((idx / dim4) / dim3) / dim2) / dim1;
            const int j = (((idx / dim4) / dim3) / dim2) % dim1;
            const int k =  ((idx / dim4) / dim3) % dim2;
            const int l =   (idx / dim4) % dim3;
            const int m =    idx % dim4;
            if (isnan(v(i,j,k,l,m))) {
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
