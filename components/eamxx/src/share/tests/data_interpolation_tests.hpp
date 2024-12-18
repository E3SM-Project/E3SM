#ifndef EAMXX_DATA_INTERPOLATION_TESTS_HPP
#define EAMXX_DATA_INTERPOLATION_TESTS_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/field/field.hpp"
#include "share/util/scream_universal_constants.hpp"

#include <vector>

namespace scream
{

constexpr int ncmps = 2;
constexpr auto spd = constants::seconds_per_day;
constexpr int data_ngcols = 12;
constexpr int fine_ngcols = 2*data_ngcols-1; // stick one dof between each data dofs
constexpr int data_nlevs = 32;
constexpr int fine_nlevs = 64;
const std::string map_file_name = "map_file_for_data_interp.nc";

// At each month in the input data, we are adding a delta to the "base" value of the fields.
constexpr double delta_data[12] = {0, 30, 60, 90, 120, 150, 180, 150, 120, 90, 60, 30};

inline util::TimeStamp get_t_ref () {
  return util::TimeStamp ({2010,1,1},{0,0,0});
}

// Slices are at midnight between 15th and 16th of each month
// First slice is Jul 15th
inline util::TimeStamp get_first_slice_time () {
  // 15 days after the reference time
  auto t = get_t_ref() + spd*15;      // Mid Jan
  for (int mm=0; mm<6; ++mm) {
    t += spd*t.days_in_curr_month();
  }
  return t;
}

inline util::TimeStamp get_last_slice_time () {
  // 11 months after the 1st slice
  auto t = get_first_slice_time();
  for (int mm=0; mm<11; ++mm) {
    t += spd*t.days_in_curr_month();
  }
  return t;
}

std::vector<Field>
create_fields (const std::shared_ptr<const AbstractGrid>& grid,
               const bool init_values,
               const bool with_pressure = false)
{
  constexpr auto m  = ekat::units::m;
  constexpr auto Pa = ekat::units::Pa;
  const auto& gn = grid->name();

  auto int_same_as_mid = with_pressure;

  // Create fields

  auto layout_s2d   = grid->get_2d_scalar_layout();
  auto layout_v2d   = grid->get_2d_vector_layout(ncmps);
  auto layout_s3d_m = grid->get_3d_scalar_layout(true);
  auto layout_v3d_m = grid->get_3d_vector_layout(true,ncmps);
  auto layout_s3d_i = grid->get_3d_scalar_layout(int_same_as_mid);
  auto layout_v3d_i = grid->get_3d_vector_layout(int_same_as_mid,ncmps);

  Field s2d  (FieldIdentifier("s2d",   layout_s2d,   m, gn));
  Field v2d  (FieldIdentifier("v2d",   layout_v2d,   m, gn));
  Field s3d_m(FieldIdentifier("s3d_m", layout_s3d_m, m, gn));
  Field v3d_m(FieldIdentifier("v3d_m", layout_v3d_m, m, gn));
  Field s3d_i(FieldIdentifier("s3d_i", layout_s3d_i, m, gn));
  Field v3d_i(FieldIdentifier("v3d_i", layout_v3d_i, m, gn));

  s3d_m.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  v3d_m.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  s3d_i.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  v3d_i.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);

  s2d.allocate_view();
  v2d.allocate_view();
  s3d_m.allocate_view();
  v3d_m.allocate_view();
  s3d_i.allocate_view();
  v3d_i.allocate_view();

  if (init_values) {
    // We set horiz/vert values based on the position of the dof, assuming that
    //  - we have a 1d horiz grid
    //  - leftmost h dof gets a value of 0, rightmost a value of 1.0
    //  - bottom interface gets a value of 0, top interface get a value of 1.0
    //  - midpoints are the avg of interfaces
    //  - we do hvalue*prod_other_dims + icmp*num_v_levs + v_value
    int ncols = grid->get_num_local_dofs();
    int nlevs = grid->get_num_vertical_levels();
    int nlevsp1 = nlevs+1;
    int ngcols = grid->get_num_global_dofs();
    int nh_intervals = ngcols - 1;
    double h_value, v_value;
    double h_max = 1.0;
    double v_max = 1.0;
    double dh = h_max / nh_intervals;
    double dv = v_max / nlevs;
    auto gids = grid->get_dofs_gids().get_view<const int*,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      auto gid = gids[icol];
      h_value = gid*dh;
      // 3D quantities
      for (int ilev=0; ilev<nlevs; ++ilev) {
        v_value = ilev*dv;
        for (int icmp=0; icmp<ncmps; ++icmp) {
          v3d_m.get_view<Real***,Host>()(icol,icmp,ilev) = h_value*ncmps*nlevs   + icmp*nlevs   + v_value + dv/2;
          v3d_i.get_view<Real***,Host>()(icol,icmp,ilev) = h_value*ncmps*nlevsp1 + icmp*nlevsp1 + v_value;
        }
        s3d_m.get_view<Real**,Host>()(icol,ilev) = h_value*nlevs   + v_value;
        s3d_i.get_view<Real**,Host>()(icol,ilev) = h_value*nlevsp1 + v_value;
      }
      // Last interface (if mid!=int)
      if (not int_same_as_mid) {
        s3d_i.get_view<Real**,Host>()(icol,nlevs) = h_value*nlevsp1 + v_max;
        for (int icmp=0; icmp<ncmps; ++icmp) {
          v3d_i.get_view<Real***,Host>()(icol,icmp,nlevs) = h_value*ncmps*nlevsp1 + icmp*nlevsp1 + v_max;
        }
      }

      // 2D quantities
      for (int icmp=0; icmp<ncmps; ++icmp) {
        v2d.get_view<Real**,Host>()(icol,icmp) = h_value*ncmps + icmp;
      }
      s2d.get_view<Real*,Host>()(icol) = h_value;
    }

    s2d.sync_to_dev();
    v2d.sync_to_dev();
    s3d_m.sync_to_dev();
    v3d_m.sync_to_dev();
    s3d_i.sync_to_dev();
    v3d_i.sync_to_dev();
  }

  if (with_pressure) {
    // Don't just clone s3d_m, so we actually get "Pa" in the nc file units
    Field p3d(FieldIdentifier("p3d",layout_s3d_m,Pa,gn));
    p3d.allocate_view();
    p3d.deep_copy(s3d_m);

    return {s2d, v2d, s3d_m, v3d_m, s3d_i, v3d_i, p3d};
  } else {
    return {s2d, v2d, s3d_m, v3d_m, s3d_i, v3d_i};
  }
}

} // namespace scream

#endif // EAMXX_DATA_INTERPOLATION_TESTS_HPP
