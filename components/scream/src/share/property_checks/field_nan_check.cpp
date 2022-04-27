#include "share/property_checks/field_nan_check.hpp"
#include "share/util//scream_array_utils.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

FieldNaNCheck::
FieldNaNCheck (const Field& f)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (f.rank()<=6,
      "Error in FieldNaNCheck constructor: unsupported field rank.\n"
      "  - Field name: " + f.name() << "\n"
      "  - Field rank: " + std::to_string(f.rank()) + "\n");
  EKAT_REQUIRE_MSG (field_valid_data_types().has_v(f.data_type()),
      "Error in FieldNaNCheck constructor: field data type not supported.\n"
      "  - Field name: " + f.name() << "\n"
      "  - Field rank: " + std::to_string(f.rank()) + "\n");

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
      EKAT_ERROR_MSG (
          "Internal error in FieldNaNCheck: unsupported field rank.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }

  PropertyCheck::CheckResult check_result;
  check_result.pass = (num_invalid == 0);
  check_result.msg = "";
  if (not check_result.pass) {
    check_result.msg = std::string("FieldNaNCheck failed; ") + std::to_string(num_invalid) + " invalid values found\n";
  }
  return check_result;
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
      EKAT_ERROR_MSG (
          "Internal error in FieldNaNCheck: unsupported field data type.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }
}

} // namespace scream
