#include "share/property_checks/field_nan_check.hpp"
#include "share/util//scream_view_utils.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

FieldNaNCheck::
FieldNaNCheck (const Field& f)
{
  // We can't repair NaN's.
  set_fields ({f},{false});
}

template<typename ST>
PropertyCheck::CheckResult FieldNaNCheck::check_impl() const {
  using const_ST    = typename std::add_const<ST>::type;

  const auto& f = fields().front();

  const auto& layout = f.get_header().get_identifier().get_layout();
  const auto extents = layout.extents();
  const auto size = layout.size();

  int num_invalid = 0;
  switch (layout.rank()) {
    case 1:
      {
        auto v = f.template get_view<const_ST*>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, int& result) {
          if (ekat::is_invalid(v(i))) {
            ++result;
          }
        }, Kokkos::Sum<int>(num_invalid));
      }
      break;
    case 2:
      {
        auto v = f.template get_view<const_ST**>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j;
          unflatten_idx(idx,extents,i,j);
          if (ekat::is_invalid(v(i,j))) {
            ++result;
          }
        }, Kokkos::Sum<int>(num_invalid));
      }
      break;
    case 3:
      {
        auto v = f.template get_view<const_ST***>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k;
          unflatten_idx(idx,extents,i,j,k);
          if (ekat::is_invalid(v(i,j,k))) {
            ++result;
          }
        }, Kokkos::Sum<int>(num_invalid));
      }
      break;
    case 4:
      {
        auto v = f.template get_view<const_ST****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k,l;
          unflatten_idx(idx,extents,i,j,k,l);
          if (ekat::is_invalid(v(i,j,k,l))) {
            ++result;
          }
        }, Kokkos::Sum<int>(num_invalid));
      }
      break;
    case 5:
      {
        auto v = f.template get_view<const_ST*****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k,l,m;
          unflatten_idx(idx,extents,i,j,k,l,m);
          if (ekat::is_invalid(v(i,j,k,l,m))) {
            ++result;
          }
        }, Kokkos::Sum<int>(num_invalid));
      }
      break;
    case 6:
      {
        auto v = f.template get_view<const_ST******>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k,l,m,n;
          unflatten_idx(idx,extents,i,j,k,l,m,n);
          if (ekat::is_invalid(v(i,j,k,l,m,n))) {
            ++result;
          }
        }, Kokkos::Sum<int>(num_invalid));
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported f rank.\n");
  }

  return num_invalid==0 ? CheckResult::Pass : CheckResult::Fail;
}

PropertyCheck::CheckResult FieldNaNCheck::check() const {
  const auto& f = fields().front();

  switch (f.data_type()) {
    case DataType::IntType:
      return check_impl<int>();
    case DataType::FloatType:
      return check_impl<float>();
    case DataType::DoubleType:
      return check_impl<double>();
    default:
      EKAT_ERROR_MSG ("Error! Field data type not supported.\n");
  }
}

} // namespace scream
