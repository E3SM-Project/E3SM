#include "share/field/field_property_checks/field_nan_check.hpp"
#include "share/util//scream_view_utils.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

template<typename ST>
bool FieldNaNCheck::check_impl(const Field& field) const {
  using const_ST    = typename std::add_const<ST>::type;

  const auto& layout = field.get_header().get_identifier().get_layout();
  const auto extents = layout.extents();
  const auto size = layout.size();

  int num_nans = 0;
  switch (layout.rank()) {
    case 1:
      {
        auto v = field.template get_view<const_ST*>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, int& result) {
          if (std::isnan(v(i))) {
            ++result;
          }
        }, Kokkos::Sum<int>(num_nans));
      }
      break;
    case 2:
      {
        auto v = field.template get_view<const_ST**>();
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
        auto v = field.template get_view<const_ST***>();
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
        auto v = field.template get_view<const_ST****>();
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
        auto v = field.template get_view<const_ST*****>();
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
        auto v = field.template get_view<const_ST******>();
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

bool FieldNaNCheck::check(const Field& field) const {
  const auto& dt = field.get_header().get_identifier().data_type();

  bool check;
  if (dt=="int") {
    check = check_impl<int>(field);
  } else if (dt=="float") {
    check = check_impl<float>(field);
  } else if (dt=="double") {
    check = check_impl<double>(field);
  } else {
    EKAT_ERROR_MSG ("Error! Field data type not supported.\n");
  }

  return check;
}

} // namespace scream
