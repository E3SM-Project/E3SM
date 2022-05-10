#include "share/property_checks/field_within_interval_check.hpp"
#include "share/util/scream_array_utils.hpp"

#include <ekat/util/ekat_math_utils.hpp>

#include <sstream>

namespace scream
{

FieldWithinIntervalCheck::
FieldWithinIntervalCheck (const Field& f,
                          const std::shared_ptr<const AbstractGrid>& grid,
                          const double lower_bound,
                          const double upper_bound,
                          const bool can_repair)
 : m_lower_bound(lower_bound)
 , m_upper_bound(upper_bound)
 , m_grid (grid)
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

  EKAT_REQUIRE_MSG (grid==nullptr || f.get_header().get_identifier().get_grid_name()==grid->name(),
      "Error! The name of the input grid does not match the grid name stored in the field identifier.\n"
      "  - Field name: " + f.name() + "\n"
      "  - Field grid name: " + f.get_header().get_identifier().get_grid_name() + "\n"
      "  - Input grid name: " + grid->name() + "\n");

  set_fields ({f},{can_repair});
}

std::string FieldWithinIntervalCheck::name () const {
  // NOTE: std::to_string does not do a good job with small numbers (like 1e-9).
  std::stringstream ss;
  ss << fields().front().name()
    << " within interval [" << m_lower_bound << ", " << m_upper_bound << "]";
  return ss.str();
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

  using space_t = typename Field::device_t::execution_space;
  using minmaxloc_t = Kokkos::MinMaxLoc<nonconst_ST,int,space_t>;
  using minmaxloc_value_t = typename minmaxloc_t::value_type;
  minmaxloc_value_t minmaxloc;
  switch (layout.rank()) {
    case 1:
      {
        auto v = f.template get_view<const_ST*>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, minmaxloc_value_t& result) {
          if (v(i)<result.min_val) {
            result.min_val = v(i);
            result.min_loc = i;
          }
          if (v(i)>result.max_val) {
            result.max_val = v(i);
            result.max_loc = i;
          }
        }, minmaxloc_t(minmaxloc));
      }
      break;
    case 2:
      {
        auto v = f.template get_view<const_ST**>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmaxloc_value_t& result) {
          int i,j;
          unflatten_idx(idx,extents,i,j);
          if (v(i,j)<result.min_val) {
            result.min_val = v(i,j);
            result.min_loc = idx;
          }
          if (v(i,j)>result.max_val) {
            result.max_val = v(i,j);
            result.max_loc = idx;
          }
        }, minmaxloc_t(minmaxloc));
      }
      break;
    case 3:
      {
        auto v = f.template get_view<const_ST***>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmaxloc_value_t& result) {
          int i,j,k;
          unflatten_idx(idx,extents,i,j,k);
          if (v(i,j,k)<result.min_val) {
            result.min_val = v(i,j,k);
            result.min_loc = idx;
          }
          if (v(i,j,k)>result.max_val) {
            result.max_val = v(i,j,k);
            result.max_loc = idx;
          }
        }, minmaxloc_t(minmaxloc));
      }
      break;
    case 4:
      {
        auto v = f.template get_view<const_ST****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmaxloc_value_t& result) {
          int i,j,k,l;
          unflatten_idx(idx,extents,i,j,k,l);
          if (v(i,j,k,l)<result.min_val) {
            result.min_val = v(i,j,k,l);
            result.min_loc = idx;
          }
          if (v(i,j,k,l)>result.max_val) {
            result.max_val = v(i,j,k,l);
            result.max_loc = idx;
          }
        }, minmaxloc_t(minmaxloc));
      }
      break;
    case 5:
      {
        auto v = f.template get_view<const_ST*****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmaxloc_value_t& result) {
          int i,j,k,l,m;
          unflatten_idx(idx,extents,i,j,k,l,m);
          if (v(i,j,k,l,m)<result.min_val) {
            result.min_val = v(i,j,k,l,m);
            result.min_loc = idx;
          }
          if (v(i,j,k,l,m)>result.max_val) {
            result.max_val = v(i,j,k,l,m);
            result.max_loc = idx;
          }
        }, minmaxloc_t(minmaxloc));
      }
      break;
    case 6:
      {
        auto v = f.template get_view<const_ST******>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, minmaxloc_value_t& result) {
          int i,j,k,l,m,n;
          unflatten_idx(idx,extents,i,j,k,l,m,n);
          if (v(i,j,k,l,m,n)<result.min_val) {
            result.min_val = v(i,j,k,l,m,n);
            result.min_loc = idx;
          }
          if (v(i,j,k,l,m,n)>result.max_val) {
            result.max_val = v(i,j,k,l,m,n);
            result.max_loc = idx;
          }
        }, minmaxloc_t(minmaxloc));
      }
      break;
    default:
      EKAT_ERROR_MSG (
          "Internal error in FieldWithinIntervalCheck: unsupported field rank.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }
  PropertyCheck::CheckResult check_result;
  check_result.pass = minmaxloc.min_val>=m_lower_bound && minmaxloc.max_val<=m_upper_bound;
  if (not check_result.pass) {
    check_result.msg  = "Check failed.\n";
    check_result.msg += "  - check name: " + this->name() + "\n";
    check_result.msg += "  - field id: " + f.get_header().get_identifier().get_id_string() + "\n";
  } else {
    check_result.msg  = "Check passed.\n";
    check_result.msg += "  - check name:" + this->name() + "\n";
    check_result.msg += "  - field id: " + f.get_header().get_identifier().get_id_string() + "\n";
  }

  auto idx_min = unflatten_idx(layout.dims(),minmaxloc.min_loc);
  auto idx_max = unflatten_idx(layout.dims(),minmaxloc.max_loc);

  using namespace ShortFieldTagsNames;

  int min_col_lid, max_col_lid;
  AbstractGrid::dofs_list_h_type gids;
  AbstractGrid::geo_view_h_type lat, lon;
  bool has_latlon;
  bool has_col_info = m_grid and layout.tag(0)==COL;

  if (has_col_info) {
    // We are storing grid info, and the field is over columns. Get col id and coords.
    min_col_lid = idx_min[0];
    max_col_lid = idx_max[0];
    gids = m_grid->get_dofs_gids_host();
    has_latlon = m_grid->has_geometry_data("lat") && m_grid->has_geometry_data("lon");
    if (has_latlon) {
      lat = m_grid->get_geometry_data_host("lat");
      lon = m_grid->get_geometry_data_host("lon");
    }
  }

  std::stringstream msg;
  msg << "  - minimum:\n";
  msg << "    - value: " << minmaxloc.min_val << "\n";
  if (has_col_info) {
    msg << "    - entry: (" << gids(min_col_lid);
    for (size_t i=1; i<idx_min.size(); ++i) {
      msg << "," << i;
    }
    msg << ")\n";
    if (has_latlon) {
      msg << "    - lat/lon: (" << lat(min_col_lid) << ", " << lon(min_col_lid) << ")\n";
    }
  }

  msg << "  - maximum:\n";
  msg << "    - value: " << minmaxloc.max_val << "\n";
  if (has_col_info) {
    msg << "    - entry: (" << gids(max_col_lid);
    for (size_t i=1; i<idx_max.size(); ++i) {
      msg << "," << i;
    }
    msg << ")\n";
    if (has_latlon) {
      msg << "    - lat/lon: (" << lat(max_col_lid) << ", " << lon(max_col_lid) << ")\n";
    }
  }

  check_result.msg += msg.str();

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
