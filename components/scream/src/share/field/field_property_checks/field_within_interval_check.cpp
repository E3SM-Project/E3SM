#include "share/field/field_property_checks/field_within_interval_check.hpp"

#include <ekat/util/ekat_math_utils.hpp>

namespace scream
{

template<typename ST>
bool FieldWithinIntervalCheck::check_impl (const Field& field) const
{
  using const_ST    = typename std::add_const<ST>::type;
  using nonconst_ST = typename std::remove_const<ST>::type;

  const auto& layout = field.get_header().get_identifier().get_layout();
  const auto extents = layout.extents();
  const auto size = layout.size();

  using minmax_t = typename Kokkos::MinMax<nonconst_ST>::value_type;
  minmax_t minmax;
  switch (layout.rank()) {
    case 1:
      {
        auto v = field.template get_view<const_ST*>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, minmax_t& result) {
          result.min_val = ekat::impl::min(result.min_val, v(i));
          result.max_val = ekat::impl::max(result.max_val, v(i));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    case 2:
      {
        auto v = field.template get_view<const_ST**>();
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
        auto v = field.template get_view<const_ST***>();
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
        auto v = field.template get_view<const_ST****>();
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
        auto v = field.template get_view<const_ST*****>();
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
        auto v = field.template get_view<const_ST******>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmax_t& result) {
          int i,j,k,l,m,n;
          unflatten_idx(idx,extents,i,j,k,l,m,n);
          result.min_val = ekat::impl::min(result.min_val, v(i,j,k,l,m,n));
          result.max_val = ekat::impl::max(result.max_val, v(i,j,k,l,m,n));
        }, Kokkos::MinMax<ST>(minmax));
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }
  return minmax.min_val>=m_lower_bound && minmax.max_val<=m_upper_bound;
}

bool FieldWithinIntervalCheck::check(const Field& field) const {
  switch (field.data_type()) {
    case DataType::IntType:
      return check_impl<int>(field);
    case DataType::FloatType:
      return check_impl<float>(field);
    case DataType::DoubleType:
      return check_impl<double>(field);
    default:
      EKAT_ERROR_MSG ("Error! Field data type not supported.\n");
  }
}

template<typename ST>
void FieldWithinIntervalCheck::repair_impl(Field& field) const
{
  using nonconst_ST = typename std::remove_const<ST>::type;

  const auto& layout = field.get_header().get_identifier().get_layout();
  const auto extents = layout.extents();
  const auto size = layout.size();

  ST lb = m_lower_bound;
  ST ub = m_upper_bound;
  switch (layout.rank()) {
    case 1:
      {
        auto v = field.template get_view<nonconst_ST*>();
        Kokkos::parallel_for(size, KOKKOS_LAMBDA(int i) {
          auto& ref = v(i);
          ref = ekat::impl::min(ub, ref);
          ref = ekat::impl::max(lb, ref);
        });
      }
      break;
    case 2:
      {
        auto v = field.template get_view<nonconst_ST**>();
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
        auto v = field.template get_view<nonconst_ST***>();
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
        auto v = field.template get_view<nonconst_ST****>();
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
        auto v = field.template get_view<nonconst_ST*****>();
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
        auto v = field.template get_view<nonconst_ST******>();
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

void FieldWithinIntervalCheck::repair(Field& field) const {

  EKAT_REQUIRE_MSG (can_repair(),
      "Error! Field property check misses repair capability.\b"
      "  - Property check: " + name() + "\n"
      "  - field name    : " + field.name() + "\n");

  switch (field.data_type()) {
    case DataType::IntType:
      repair_impl<int>(field);
      break;
    case DataType::FloatType:
      repair_impl<float>(field);
      break;
    case DataType::DoubleType:
      repair_impl<double>(field);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Field data type not supported.\n");
  }
}

} // namespace scream
