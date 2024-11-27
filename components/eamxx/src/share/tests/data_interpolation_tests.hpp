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

// At each month in the input data, we are adding a delta to the "base" value of the fields.
constexpr double delta_data[12] = {0, 30, 60, 90, 120, 150, 180, 150, 120, 90, 60, 30};

inline util::TimeStamp get_t_ref () {
  return util::TimeStamp ({2010,1,1},{0,0,0});
}

inline util::TimeStamp get_first_slice_time () {
  // 15 days after the reference time
  return get_t_ref() + 86400*15;
}

inline util::TimeStamp get_last_slice_time () {
  // 11 months after the 1st slice
  auto t = get_first_slice_time();

  for (int mm=0; mm<23; ++mm) {
    t += constants::seconds_per_day*t.days_in_curr_month();
  }
  return t;
}

std::vector<Field>
create_fields (const std::shared_ptr<const AbstractGrid>& grid,
               const bool init_values)
{
  constexpr auto m = ekat::units::m;
  const auto& gn = grid->name();

  // Create fields and initialize their data as
  //   set data as f(icol,icmp,ilev) = icol*ncmp*nlev + icmp*nlev + ilev

  auto layout_s2d   = grid->get_2d_scalar_layout();
  auto layout_v2d   = grid->get_2d_vector_layout(ncmps);
  auto layout_s3d_m = grid->get_3d_scalar_layout(true);
  auto layout_v3d_m = grid->get_3d_vector_layout(true,ncmps);
  auto layout_s3d_i = grid->get_3d_scalar_layout(false);
  auto layout_v3d_i = grid->get_3d_vector_layout(false,ncmps);

  Field s2d  (FieldIdentifier("s2d",   layout_s2d,   m, gn));
  Field v2d  (FieldIdentifier("v2d",   layout_v2d,   m, gn));
  Field s3d_m(FieldIdentifier("s3d_m", layout_s3d_m, m, gn));
  Field v3d_m(FieldIdentifier("v3d_m", layout_v3d_m, m, gn));
  Field s3d_i(FieldIdentifier("s3d_i", layout_s3d_i, m, gn));
  Field v3d_i(FieldIdentifier("v3d_i", layout_v3d_i, m, gn));

  s2d.allocate_view();
  v2d.allocate_view();
  s3d_m.allocate_view();
  v3d_m.allocate_view();
  s3d_i.allocate_view();
  v3d_i.allocate_view();

  if (init_values) {
    int ncols = grid->get_num_local_dofs();
    int nlevs = grid->get_num_vertical_levels();
    int nlevsp1 = nlevs+1;
    for (int icol=0; icol<ncols; ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        for (int icmp=0; icmp<ncmps; ++icmp) {
          v3d_m.get_view<Real***,Host>()(icol,icmp,ilev) = icol*ncmps*nlevs   + icmp*nlevs   + ilev;
          v3d_i.get_view<Real***,Host>()(icol,icmp,ilev) = icol*ncmps*nlevsp1 + icmp*nlevsp1 + ilev;
        }
        s3d_m.get_view<Real**,Host>()(icol,ilev) = icol*nlevs   + ilev;
        s3d_i.get_view<Real**,Host>()(icol,ilev) = icol*nlevsp1 + ilev;
      }
      s3d_i.get_view<Real**,Host>()(icol,nlevs) = icol*nlevsp1 + nlevs;
      for (int icmp=0; icmp<ncmps; ++icmp) {
        v3d_i.get_view<Real***,Host>()(icol,icmp,nlevs) = icol*ncmps*nlevsp1 + icmp*nlevsp1 + nlevs;
        v2d.get_view<Real**,Host>()(icol,icmp) = icol*ncmps + icmp;
      }
      s2d.get_view<Real*,Host>()(icol) = icol;
    }

    s2d.sync_to_dev();
    v2d.sync_to_dev();
    s3d_m.sync_to_dev();
    v3d_m.sync_to_dev();
    s3d_i.sync_to_dev();
    v3d_i.sync_to_dev();
  }

  return {s2d, v2d, s3d_m, v3d_m, s3d_i, v3d_i};
}

} // namespace scream

#endif // EAMXX_DATA_INTERPOLATION_TESTS_HPP
