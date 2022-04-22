#include "share/property_checks/field_within_interval_check.hpp"
#include "share/util/scream_view_utils.hpp"

#include <ekat/util/ekat_math_utils.hpp>

namespace scream
{

FieldWithinIntervalCheck::
FieldWithinIntervalCheck (const Field& f,
                          const double lower_bound,
                          const double upper_bound,
                          const bool can_repair)
 : m_lower_bound(lower_bound)
 , m_upper_bound(upper_bound)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (f.rank()<=6,
      "Error in FieldWithinIntervalCheck constructor: unsupported field rank.\n"
      "  - Field name: " + f.name() << "\n"
      "  - Field rank: " + std::to_string(f.rank()) + "\n");
  EKAT_REQUIRE_MSG (field_valid_data_types().has_v(f.data_type()),
      "Error in FieldWithinIntervalCheck constructor: field data type not supported.\n"
      "  - Field name: " + f.name() << "\n"
      "  - Field rank: " + std::to_string(f.rank()) + "\n");
  EKAT_ASSERT_MSG(lower_bound <= upper_bound,
                  "lower_bound must be less than or equal to upper_bound.");

  set_fields ({f},{can_repair});
}

template<typename ST>
PropertyCheck::CheckResult FieldWithinIntervalCheck::check_impl () const
{
  using const_ST    = typename std::add_const<ST>::type;
  using nonconst_ST = typename std::remove_const<ST>::type;

  const auto& f = fields().front();

  const auto& layout = f.get_header().get_identifier().get_layout();
  const auto extents = layout.extents();
  const auto size = layout.size();

  using minmax_t = typename Kokkos::MinMax<nonconst_ST>::value_type;
  minmax_t minmax;
  switch (layout.rank()) {
    case 1:
      {
        auto v = f.template get_view<const_ST*>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, minmax_t& result) {
          result.min_val = ekat::impl::min(result.min_val, v(i));
          result.max_val = ekat::impl::max(result.max_val, v(i));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    case 2:
      {
        auto v = f.template get_view<const_ST**>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
          int i,j;
          unflatten_idx(idx,extents,i,j);
          result.min_val = ekat::impl::min(result.min_val, v(i,j));
          result.max_val = ekat::impl::max(result.max_val, v(i,j));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    case 3:
      {
        auto v = f.template get_view<const_ST***>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
          int i,j,k;
          unflatten_idx(idx,extents,i,j,k);
          result.min_val = ekat::impl::min(result.min_val, v(i,j,k));
          result.max_val = ekat::impl::max(result.max_val, v(i,j,k));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    case 4:
      {
        auto v = f.template get_view<const_ST****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
          int i,j,k,l;
          unflatten_idx(idx,extents,i,j,k,l);
          result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l));
          result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    case 5:
      {
        auto v = f.template get_view<const_ST*****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
          int i,j,k,l,m;
          unflatten_idx(idx,extents,i,j,k,l,m);
          result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l,m));
          result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l,m));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    case 6:
      {
        auto v = f.template get_view<const_ST******>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
          int i,j,k,l,m,n;
          unflatten_idx(idx,extents,i,j,k,l,m,n);
          result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l,m,n));
          result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l,m,n));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    default:
      EKAT_ERROR_MSG (
          "Internal error in FieldWithinIntervalCheck: unsupported field rank.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }
  PropertyCheck::CheckResult check_result;
  check_result.pass = minmax.min_val>=m_lower_bound && minmax.max_val<=m_upper_bound;
  check_result.msg = "";
  if (not check_result.pass) {
    check_result.msg = std::string("FieldWithinIntervalCheck failed")
                     + "; min = " + std::to_string(minmax.min_val)
                     + "; max = " + std::to_string(minmax.max_val) + "\n";
  }
  return check_result;
}

PropertyCheck::CheckResult FieldWithinIntervalCheck::check() const {
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
          "Internal error in FieldWithinIntervalCheck: unsupported field data type.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }
}

template<typename ST>
void FieldWithinIntervalCheck::repair_impl() const
{
  using nonconst_ST = typename std::remove_const<ST>::type;

  const auto& f = fields().front();

  const auto& layout = f.get_header().get_identifier().get_layout();
  const auto extents = layout.extents();
  const auto size = layout.size();

  ST lb = m_lower_bound;
  ST ub = m_upper_bound;
  switch (layout.rank()) {
    case 1:
      {
        auto v = f.template get_view<nonconst_ST*>();
        Kokkos::parallel_for(size, KOKKOS_LAMBDA(int i) {
          auto& ref = v(i);
          ref = ekat::impl::min(ub, ref);
          ref = ekat::impl::max(lb, ref);
        });
      }
      break;
    case 2:
      {
        auto v = f.template get_view<nonconst_ST**>();
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
        auto v = f.template get_view<nonconst_ST***>();
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
        auto v = f.template get_view<nonconst_ST****>();
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
        auto v = f.template get_view<nonconst_ST*****>();
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
        auto v = f.template get_view<nonconst_ST******>();
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
      EKAT_ERROR_MSG (
          "Internal error in FieldWithinIntervalCheck: unsupported field rank.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }
}

void FieldWithinIntervalCheck::repair_impl() const {

  const auto& f = fields().front();

  switch (f.data_type()) {
    case DataType::IntType:
      repair_impl<int>();
      break;
    case DataType::FloatType:
      repair_impl<float>();
      break;
    case DataType::DoubleType:
      repair_impl<double>();
      break;
    default:
      EKAT_ERROR_MSG (
          "Internal error in FieldWithinIntervalCheck: unsupported field data type.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }
}

} // namespace scream
