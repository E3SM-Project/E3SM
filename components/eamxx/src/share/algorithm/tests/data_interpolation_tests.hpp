#ifndef EAMXX_DATA_INTERPOLATION_TESTS_HPP
#define EAMXX_DATA_INTERPOLATION_TESTS_HPP

#include "share/algorithm/eamxx_data_interpolation.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <catch2/catch.hpp>

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
// The delta is almost periodic, but the 2nd year it adds 1 to the 1st year values
constexpr int num_data_months = 24;
constexpr double delta_data[num_data_months] = {
  0, 30, 60, 90, 120, 150, 180, 150, 120, 90, 60, 30,
  1, 31, 61, 91, 121, 151, 181, 151, 121, 91, 61, 31
};

// Give ourselves some room for roundoff errors, since our manual
// evaluation may be different (in finite prec) from the one in the class.
constexpr auto tol = std::numeric_limits<Real>::epsilon()*10;

constexpr auto NOP = DataInterpolation::None;
constexpr auto P1D = DataInterpolation::Static1D;
constexpr auto P2D = DataInterpolation::Dynamic3DRef;
constexpr auto P3D = DataInterpolation::Dynamic3D;

using strvec_t = std::vector<std::string>;

inline util::TimeStamp get_t_ref () {
  return util::TimeStamp ({0001,1,1},{0,0,0});
}

// Slices are at midnight between 15th and 16th of each month
// First slice is Jun 15th
inline util::TimeStamp get_first_slice_time () {
  return get_t_ref() + spd*15; // Mid of the month
}

inline util::TimeStamp get_last_slice_time () {
  // 11 months after the 1st slice
  auto t = get_first_slice_time();
  for (int mm=0; mm<(num_data_months-1); ++mm) {
    t += spd*t.days_in_curr_month();
  }
  return t;
}

util::TimeStamp reset_year (const util::TimeStamp& t_in, int yy)
{
  auto date = t_in.get_date();
  auto time = t_in.get_time();
  date[0] = yy;
  return util::TimeStamp(date,time);
}

std::vector<Field>
create_fields (const std::shared_ptr<const AbstractGrid>& grid,
               const bool init_values,
               const bool pad_for_packing = true)
{
  using namespace ShortFieldTagsNames;
  constexpr auto m  = ekat::units::m;
  const auto& gn = grid->name();

  int ncols = grid->get_num_local_dofs();
  int nlevs = grid->get_num_vertical_levels();
  bool p_grid = grid->get_vkind()==AbstractGrid::VKind::Pressure;

  // Create fields
  auto tag_mid = p_grid ? LEVP : LEV;
  auto tag_int = p_grid ? LEVP : ILEV;
  auto layout_s2d   = grid->get_2d_scalar_layout();
  auto layout_v2d   = grid->get_2d_vector_layout(ncmps);
  auto layout_s3d_m = grid->get_3d_scalar_layout(tag_mid);
  auto layout_v3d_m = grid->get_3d_vector_layout(tag_mid, ncmps);
  auto layout_s3d_i = grid->get_3d_scalar_layout(tag_int);
  auto layout_v3d_i = grid->get_3d_vector_layout(tag_int, ncmps);

  Field s2d  (FieldIdentifier("s2d",   layout_s2d,   m, gn));
  Field v2d  (FieldIdentifier("v2d",   layout_v2d,   m, gn));
  Field s3d_m(FieldIdentifier("s3d_m", layout_s3d_m, m, gn));
  Field v3d_m(FieldIdentifier("v3d_m", layout_v3d_m, m, gn));
  Field s3d_i(FieldIdentifier("s3d_i", layout_s3d_i, m, gn));
  Field v3d_i(FieldIdentifier("v3d_i", layout_v3d_i, m, gn));

  if (pad_for_packing) {
    s3d_m.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
    v3d_m.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
    s3d_i.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
    v3d_i.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  }

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
    //  - we do h_value*v_value + icmp
    int ngcols = grid->get_num_global_dofs();
    int nh_intervals = ngcols - 1;
    double h_value, v_value;
    double h_max = 1.0;
    double v_max = 1.0;
    double dh = nh_intervals>0 ? h_max / nh_intervals : 0;
    double dv = v_max / nlevs;
    auto gids = grid->get_dofs_gids().get_view<const int*,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      auto gid_0based = gids[icol]-1;
      h_value = gid_0based*dh + 1.0; // +1.0 to avoid a whole column full of zeros
      // 3D quantities
      for (int ilev=0; ilev<nlevs; ++ilev) {
        v_value = ilev*dv;
        for (int icmp=0; icmp<ncmps; ++icmp) {
          v3d_m.get_view<Real***,Host>()(icol,icmp,ilev) = h_value*(v_value + dv/2) + icmp;
          if (p_grid) {
            v3d_i.get_view<Real***,Host>()(icol,icmp,ilev) = h_value*(v_value + dv/2) + icmp;
          } else {
            v3d_i.get_view<Real***,Host>()(icol,icmp,ilev) = h_value*(v_value) + icmp;
          }
        }
        s3d_m.get_view<Real**,Host>()(icol,ilev) = h_value*(v_value + dv/2);
        if (p_grid) {
          s3d_i.get_view<Real**,Host>()(icol,ilev) = h_value*(v_value + dv/2);
        } else {
          s3d_i.get_view<Real**,Host>()(icol,ilev) = h_value*(v_value);
        }
      }
      // Last interface (if mid!=int), where v_value=1
      if (not p_grid) {
        s3d_i.get_view<Real**,Host>()(icol,nlevs) = h_value;
        for (int icmp=0; icmp<ncmps; ++icmp) {
          v3d_i.get_view<Real***,Host>()(icol,icmp,nlevs) = h_value + icmp;
        }
      }

      // 2D quantities
      for (int icmp=0; icmp<ncmps; ++icmp) {
        v2d.get_view<Real**,Host>()(icol,icmp) = h_value + icmp;
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

  auto p1d = s3d_m.subfield(0,0).clone("p1d", CloneFlags::All);
  auto comm = grid->get_comm();
  comm.broadcast(p1d.get_internal_view_data<Real,Host>(),nlevs,0);
  p1d.sync_to_dev();

  return {s2d, v2d, s3d_m, v3d_m, s3d_i, v3d_i, p1d};
}

std::shared_ptr<DataInterpolation>
create_interp (const std::shared_ptr<const AbstractGrid>& grid,
               const std::vector<Field>& fields)
{
  return std::make_shared<DataInterpolation>(grid,fields);
}

void root_print (const ekat::Comm& comm,
                 const std::string& msg)
{
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

// Run the data interpolation to the input grid, and check against expected values.
// t_beg is the start of the first data interval (used to compute expected deltas).
void run_tests (const std::shared_ptr<const AbstractGrid>& grid,
                const strvec_t& input_files, util::TimeStamp t_beg,
                const util::TimeLine timeline,
                const DataInterpolation::TimeInterpType time_interp_type,
                const DataInterpolation::VRemapType vr_type = DataInterpolation::None)
{
  using namespace ShortFieldTagsNames;

  auto t_end = t_beg + t_beg.days_in_curr_month()*spd;
  auto t0 = t_beg + (t_end-t_beg)/2;

  auto ncols = grid->get_num_local_dofs();
  auto nlevs = grid->get_num_vertical_levels();

  auto vcoarse_grid = grid->clone("vcoarse",true);
  auto vcoarse_vkind = vr_type==NOP ? AbstractGrid::VKind::Model
                                    : AbstractGrid::VKind::Pressure;

  vcoarse_grid->reset_vertical_configuration(data_nlevs, vcoarse_vkind);

  // These are the fields we will compute
  auto fields = create_fields(grid,false);
  fields.pop_back(); // We don't interpolate p1d...

  std::string map_file = grid->get_num_global_dofs()==data_ngcols ? "" : map_file_name;

  // These are used to check the answer
  auto base     = create_fields(grid,true);
  auto ones     = create_fields(grid,false);
  auto diff     = create_fields(grid,false);
  auto expected = create_fields(grid,false);
  for (auto& f : ones) {
    f.deep_copy(1);
  }

  auto model_pmid = base[2].clone("pmid",CloneFlags::All);
  auto model_pint = base[4].clone("pint",CloneFlags::All);
  if (vr_type==P1D) {
    // It's complicated to test the static profile, since we'd have to really do
    // a manual interpolation. But setting all model pressure equal to the 1st col
    // of the pressure makes things doable
    auto comm = grid->get_comm();
    for (const Field& p : {model_pmid, model_pint}) {
      auto col_0 = p.subfield(0,0);
      auto len = col_0.get_header().get_identifier().get_layout().size();
      comm.broadcast(col_0.get_internal_view_data<Real,Host>(),len,0);
      col_0.sync_to_dev();

      for (int icol=0; icol<ncols; ++icol) {
        auto col_i = p.subfield(0,icol);
        col_i.deep_copy(col_0);
      }
    }
  }

  std::vector<Field> expected_vcoarse, base_vcoarse, ones_vcoarse;
  if (vr_type!=NOP) {
    // If we do remap, there is some P0 extrapolation,
    // for which we need to know the data at the top/bot
    // NOTE: the data has NO ilev coord in this case
    expected_vcoarse = create_fields(vcoarse_grid,false);
    base_vcoarse = create_fields(vcoarse_grid,true);
    ones_vcoarse = create_fields(vcoarse_grid,false);
    for (auto& f : ones_vcoarse) {
      f.deep_copy(1);
    }
  }

  DataInterpolation::VertRemapData vremap_data;
  vremap_data.vr_type = vr_type;
  vremap_data.pname =  // if vr_type==None, it's not used anyways
    vr_type==P1D ? "p1d" :
                   (vr_type==P3D ? "p3d" : "p2d");
  vremap_data.pmid = model_pmid;
  vremap_data.pint = model_pint;

  int nfields = fields.size();
  auto interp = create_interp(grid,fields);
  if (timeline==util::TimeLine::Linear)
    interp->setup_linear_time_database(input_files);
  else
    interp->setup_periodic_time_database(input_files);
  interp->create_horiz_remappers (map_file);
  interp->create_vert_remapper (vremap_data);
  interp->init_time_interpolation(t0,time_interp_type);

  // We jump ahead by 2 months, but the shift interval logic cannot keep up with
  // a dt that long, so we should get an error due to the interpolation param being
  // outside the [0,1] interval.
  REQUIRE_THROWS (interp->run(t0+60*spd));

  // Loop for two year at a 20 day increment
  bool periodic = timeline==util::TimeLine::YearlyPeriodic;
  int dt = 20*spd;
  for (auto time = t0+dt; time.days_from(t0)<365; time+=dt) {
    if (t_end<time) {
      // update t_beg/t_end
      t_beg = t_end;
      t_end += t_end.days_in_curr_month()*spd;
    }

    // Compute the delta to add to field base value for mm_beg and mm_end
    int mm_beg = t_beg.get_month()-1;
    int mm_end = t_end.get_month()-1;

    // Since input data is at the 15th of the month, we need to compute
    // the distance between current time and the beg slice, then use it
    // to do a convex interpolation between f(t=t_beg) and f(t=t_end).
    util::TimeInterval time_from_beg(t_beg,time,timeline);
    double alpha = time_from_beg.length / t_beg.days_in_curr_month();

    if (time_interp_type==DataInterpolation::Nearest) {
      alpha = alpha>0.5 ? 1 : 0;
    }
    int data_mm_beg = mm_beg + (periodic ? 0 : 12*(t_beg.get_year()-get_t_ref().get_year()));
    int data_mm_end = mm_end + (periodic ? 0 : 12*(t_end.get_year()-get_t_ref().get_year()));
    double delta = delta_data[data_mm_beg]*(1-alpha) + delta_data[data_mm_end]*alpha;

    // Compute expected difference from base value
    interp->run(time);
    for (int i=0; i<nfields; ++i) {
      const auto& layout = fields[i].get_header().get_identifier().get_layout();
      bool midpoints  = layout.has_tag(LEV);
      bool interfaces = layout.has_tag(ILEV);
      bool has_levels = midpoints or interfaces;

      // Compute expected field, possibly clipping the top/bot value, if vr_type!=None
      // Since we run with 2x the number of levels as in the data, we get ONE level
      // out-of-bounds at both top and bot. The easiest way to convince yourself is to
      // make a quick diagram. Recall that in case of vremap, the input data is at MIDPOINTS
      // for both lev/ilev fields. The vfine grid has interfaces corresponding to the union
      // of midpoints and interfaces on the data grid.
      expected[i].update(base[i],1,0);
      expected[i].update(ones[i],delta,1.0);
      if ( vr_type!=NOP and has_levels) {
        auto vlen = midpoints ? nlevs : nlevs+1;
        auto scalar = not layout.is_vector_layout();

        expected_vcoarse[i].update(base_vcoarse[i],1,0);
        expected_vcoarse[i].update(ones_vcoarse[i],delta,1.0);

        expected[i].sync_to_host();
        expected_vcoarse[i].sync_to_host();
        if (scalar) {
          auto e = expected[i].get_view<Real**,Host>();
          auto e_vcoarse = expected_vcoarse[i].get_view<const Real**,Host>();

          for (int icol=0; icol<ncols; ++icol) {
            e(icol,0) = e_vcoarse(icol,0);
            e(icol,vlen-1) = e_vcoarse(icol,data_nlevs-1);
          }
        } else {
          auto e = expected[i].get_view<Real***,Host>();
          auto e_vcoarse = expected_vcoarse[i].get_view<const Real***,Host>();

          for (int icol=0; icol<ncols; ++icol) {
            for (int idim=0; idim<layout.get_vector_dim(); ++idim) {
              e(icol,idim,0) = e_vcoarse(icol,idim,0);
              e(icol,idim,vlen-1) = e_vcoarse(icol,idim,data_nlevs-1);
            }
          }
        }
        expected[i].sync_to_dev();
      }

      // Compute the rel error
      diff[i].update(fields[i],1,1);
      diff[i].update(expected[i],-1,1);
      diff[i].scale(1.0 / frobenius_norm(expected[i]).as<Real>());
      REQUIRE (frobenius_norm(diff[i]).as<Real>()<tol);
    }
  }
}

// Run tests for static (time-independent) data.
// expected_delta: the delta that was added to base_fields when writing the file
//                 (0.0 for a no-time-dep file; delta_data[mm_idx] for a timed snapshot).
void run_static_tests (const std::shared_ptr<const AbstractGrid>& grid,
                       const std::string& input_file,
                       int time_index,
                       double expected_delta,
                       const DataInterpolation::VRemapType vr_type = DataInterpolation::None)
{
  using namespace ShortFieldTagsNames;

  auto ncols = grid->get_num_local_dofs();
  auto nlevs = grid->get_num_vertical_levels();

  auto vcoarse_grid = grid->clone("vcoarse_s",true);
  auto vcoarse_vkind = vr_type==NOP ? AbstractGrid::VKind::Model
                                    : AbstractGrid::VKind::Pressure;
  vcoarse_grid->reset_vertical_configuration(data_nlevs, vcoarse_vkind);

  auto fields = create_fields(grid, false);
  fields.pop_back(); // Don't interpolate p1d

  std::string map_file = grid->get_num_global_dofs()==data_ngcols ? "" : map_file_name;

  auto base     = create_fields(grid, true);
  auto ones     = create_fields(grid, false);
  auto diff     = create_fields(grid, false);
  auto expected = create_fields(grid, false);
  for (auto& f : ones) { f.deep_copy(1); }

  auto model_pmid = base[2].clone("pmid",CloneFlags::CopyData | CloneFlags::MatchPacking);
  auto model_pint = base[4].clone("pint",CloneFlags::CopyData | CloneFlags::MatchPacking);
  if (vr_type==P1D) {
    auto comm = grid->get_comm();
    for (const Field& p : {model_pmid, model_pint}) {
      auto col_0 = p.subfield(0,0);
      auto len = col_0.get_header().get_identifier().get_layout().size();
      comm.broadcast(col_0.get_internal_view_data<Real,Host>(), len, 0);
      col_0.sync_to_dev();
      for (int icol=0; icol<ncols; ++icol) {
        p.subfield(0,icol).deep_copy(col_0);
      }
    }
  }

  std::vector<Field> expected_vcoarse, base_vcoarse, ones_vcoarse;
  if (vr_type!=NOP) {
    expected_vcoarse = create_fields(vcoarse_grid, false);
    base_vcoarse     = create_fields(vcoarse_grid, true);
    ones_vcoarse     = create_fields(vcoarse_grid, false);
    for (auto& f : ones_vcoarse) { f.deep_copy(1); }
  }

  DataInterpolation::VertRemapData vremap_data;
  vremap_data.vr_type = vr_type;
  vremap_data.pname   = vr_type==P1D ? "p1d" : (vr_type==P3D ? "p3d" : "p2d");
  vremap_data.pmid    = model_pmid;
  vremap_data.pint    = model_pint;

  int nfields = fields.size();
  auto interp = create_interp(grid, fields);
  interp->setup_static_database({input_file}, time_index);
  interp->create_horiz_remappers(map_file);
  interp->create_vert_remapper(vremap_data);

  // Run multiple times: output should always be the same
  for (int irun=0; irun<3; ++irun) {
    interp->run();

    for (int i=0; i<nfields; ++i) {
      const auto& layout = fields[i].get_header().get_identifier().get_layout();
      bool midpoints  = layout.has_tag(LEV);
      bool interfaces = layout.has_tag(ILEV);
      bool has_levels = midpoints or interfaces;

      expected[i].update(base[i],  1, 0);
      expected[i].update(ones[i], expected_delta, 1.0);

      if (vr_type!=NOP and has_levels) {
        auto vlen   = midpoints ? nlevs : nlevs+1;
        auto scalar = not layout.is_vector_layout();

        expected_vcoarse[i].update(base_vcoarse[i],  1, 0);
        expected_vcoarse[i].update(ones_vcoarse[i], expected_delta, 1.0);

        expected[i].sync_to_host();
        expected_vcoarse[i].sync_to_host();
        if (scalar) {
          auto e        = expected[i].get_view<Real**,Host>();
          auto e_coarse = expected_vcoarse[i].get_view<const Real**,Host>();
          for (int icol=0; icol<ncols; ++icol) {
            e(icol,0)      = e_coarse(icol,0);
            e(icol,vlen-1) = e_coarse(icol,data_nlevs-1);
          }
        } else {
          auto e        = expected[i].get_view<Real***,Host>();
          auto e_coarse = expected_vcoarse[i].get_view<const Real***,Host>();
          for (int icol=0; icol<ncols; ++icol) {
            for (int idim=0; idim<layout.get_vector_dim(); ++idim) {
              e(icol,idim,0)      = e_coarse(icol,idim,0);
              e(icol,idim,vlen-1) = e_coarse(icol,idim,data_nlevs-1);
            }
          }
        }
        expected[i].sync_to_dev();
      }

      diff[i].update(fields[i],   1, 1);
      diff[i].update(expected[i], -1, 1);
      diff[i].scale(1.0 / frobenius_norm(expected[i]).as<Real>());
      REQUIRE (frobenius_norm(diff[i]).as<Real>() < tol);
    }
  }
}

} // namespace scream

#endif // EAMXX_DATA_INTERPOLATION_TESTS_HPP
