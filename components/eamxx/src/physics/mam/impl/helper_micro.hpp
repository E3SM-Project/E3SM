#ifndef EAMXX_MAM_HELPER_MICRO
#define EAMXX_MAM_HELPER_MICRO

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/util/ekat_lin_interp.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace scream::mam_coupling {

using namespace ShortFieldTagsNames;
using view_1d_host = typename KT::view_1d<Real>::HostMirror;

using ExeSpace = typename KT::ExeSpace;
using ESU      = ekat::ExeSpaceUtils<ExeSpace>;
using C        = scream::physics::Constants<Real>;
using LIV      = ekat::LinInterp<Real, 1>;

enum TracerFileType {
  // file with PS ncol, lev, and time
  FORMULA_PS,
  // nc zonal file from ncremap
  ZONAL,
  // vertical emission files
  VERT_EMISSION,
};

/* Maximum number of tracers (or fields) that the tracer reader can handle.
 Note: We are not allocating memory for MAX_NVARS_TRACER tracers.
 Therefore, if a file contains more than this number, it is acceptable to
 increase this limit. Currently, Linoz files have 8 fields. */
constexpr int MAX_NVARS_TRACER = 10;
constexpr int MAX_NUM_VERT_EMISSION_FIELDS = 25;

// Linoz structures to help manage all of the variables:
struct TracerTimeState {
  // Whether the timestate has been initialized.
  // The current month
  int current_month = -1;
  // Julian Date for the beginning of the month, as defined in
  //           /src/share/util/scream_time_stamp.hpp
  // See this file for definition of Julian Date.
  Real t_beg_month;
  // Current simulation Julian Date
  Real t_now;
  // Number of days in the current month, cast as a Real
  Real days_this_month;
};  // TricerTimeState

struct TracerData {
  TracerData() = default;
  TracerData(const int ncol, const int nlev, const int nvars) {
    init(ncol, nlev, nvars);
  }
  void init(const int ncol, const int nlev, const int nvars) {
    ncol_  = ncol;
    nlev_  = nlev;
    nvars_ = nvars;
    EKAT_REQUIRE_MSG(
        nvars_ <= int(MAX_NVARS_TRACER),
        "Error! Number of variables is bigger than NVARS_MAXTRACER. \n");
  }

  int ncol_{-1};
  int nlev_{-1};
  int nvars_{-1};
  // We cannot use a std::vector<view_2d>
  // because we need to access these views from device.
  view_2d data[MAX_NVARS_TRACER];
  view_1d ps;
  const_view_1d hyam;
  const_view_1d hybm;

  TracerFileType file_type;

  void allocate_data_views() {
    EKAT_REQUIRE_MSG(ncol_ != int(-1), "Error! ncols has not been set. \n");
    EKAT_REQUIRE_MSG(nlev_ != int(-1), "Error! nlevs has not been set. \n");

    for(int ivar = 0; ivar < nvars_; ++ivar) {
      data[ivar] = view_2d("linoz_1", ncol_, nlev_);
    }
  }  // allocate_data_views

  void set_file_type(const TracerFileType file_type_in) {
    file_type = file_type_in;
  }

  void allocate_ps() {
    EKAT_REQUIRE_MSG(ncol_ != int(-1), "Error! ncols has not been set. \n");
    EKAT_REQUIRE_MSG(file_type == FORMULA_PS,
                     "Error! file does have the PS variable. \n");

    ps = view_1d("ps", ncol_);
  }

  void set_data_views(view_2d list_of_views[]) {
    for(int ivar = 0; ivar < nvars_; ++ivar) {
      EKAT_REQUIRE_MSG(list_of_views[ivar].data() != 0,
                       "Error! Insufficient memory  size.\n");
      data[ivar] = list_of_views[ivar];
    }
  }

  void set_data_ps(const view_1d &ps_in) {
    EKAT_REQUIRE_MSG(file_type == FORMULA_PS,
                     "Error! file does have the PS variable. \n");
    ps = ps_in;
  }

  void set_hyam_n_hybm(const std::shared_ptr<AbstractRemapper> &horiz_remapper,
                       const std::string &tracer_file_name) {
    EKAT_REQUIRE_MSG(file_type == FORMULA_PS,
                     "Error! file does have the PS variable. \n");

    // Read in hyam/hybm in start/end data
    auto nondim        = ekat::units::Units::nondimensional();
    const auto io_grid = horiz_remapper->get_src_grid();
    Field hyam_f(FieldIdentifier("hyam", io_grid->get_vertical_layout(true),
                                 nondim, io_grid->name()));
    Field hybm_f(FieldIdentifier("hybm", io_grid->get_vertical_layout(true),
                                 nondim, io_grid->name()));
    hyam_f.allocate_view();
    hybm_f.allocate_view();
    AtmosphereInput hvcoord_reader(tracer_file_name, io_grid, {hyam_f, hybm_f},
                                   true);
    hvcoord_reader.read_variables();
    hvcoord_reader.finalize();
    hyam = hyam_f.get_view<const Real *>();
    hybm = hyam_f.get_view<const Real *>();
  }
};

inline const_view_1d get_altitude_int(
    const std::shared_ptr<AbstractRemapper> &horiz_remapper,
    const std::string &tracer_file_name) {
  // Read in hyam/hybm in start/end data
  auto nondim        = ekat::units::Units::nondimensional();
  const auto io_grid = horiz_remapper->get_src_grid();
  Field altitude_int_f(FieldIdentifier("altitude_int",
                                       io_grid->get_vertical_layout(false),
                                       nondim, io_grid->name()));
  altitude_int_f.allocate_view();
  AtmosphereInput hvcoord_reader(tracer_file_name, io_grid, {altitude_int_f},
                                 true);
  hvcoord_reader.read_variables();
  hvcoord_reader.finalize();
  return altitude_int_f.get_view<const Real *>();
}  // set_altitude_int

// Direct port of components/eam/src/chemistry/utils/tracer_data.F90/vert_interp
// FIXME: I need to convert for loops to Kokkos loops.
KOKKOS_INLINE_FUNCTION
void vert_interp(int ncol, int levsiz, int pver, const view_2d &pin,
                 const const_view_2d &pmid, const view_2d &datain,
                 const view_2d &dataout,
                 // work array
                 const view_int_1d &kupper) {
  const int one = 1;
  // Initialize index array
  for(int i = 0; i < ncol; ++i) {
    kupper(i) = one;
  }  // ncol

  for(int k = 0; k < pver; ++k) {
    // Top level we need to start looking is the top level for the previous k
    // for all column points
    int kkstart = levsiz;
    for(int i = 0; i < ncol; ++i) {
      kkstart = haero::min(kkstart, kupper(i));
    }

    // Store level indices for interpolation
    for(int kk = kkstart - 1; kk < levsiz - 1; ++kk) {
      for(int i = 0; i < ncol; ++i) {
        if(pin(i, kk) < pmid(i, k) && pmid(i, k) <= pin(i, kk + 1)) {
          kupper(i) = kk;
        }  // end if
      }    // end for
    }      // end kk
    // Interpolate or extrapolate...
    for(int i = 0; i < ncol; ++i) {
      if(pmid(i, k) < pin(i, 0)) {
        dataout(i, k) = datain(i, 0) * pmid(i, k) / pin(i, 0);
      } else if(pmid(i, k) > pin(i, levsiz - 1)) {
        dataout(i, k) = datain(i, levsiz - 1);
      } else {
        Real dpu = pmid(i, k) - pin(i, kupper(i));
        Real dpl = pin(i, kupper(i) + 1) - pmid(i, k);
        dataout(i, k) =
            (datain(i, kupper(i)) * dpl + datain(i, kupper(i) + 1) * dpu) /
            (dpl + dpu);
      }  // end if
    }    // end col
  }      // end k

}  // vert_interp

KOKKOS_INLINE_FUNCTION
Real linear_interp(const Real &x0, const Real &x1, const Real &t) {
  return (1 - t) * x0 + t * x1;
}  // linear_interp

// time[3]={year,month, day}
inline util::TimeStamp convert_date(const int date) {
  constexpr int ten_thousand = 10000;
  constexpr int one_hundred  = 100;

  int year  = date / ten_thousand;
  int month = (date - year * ten_thousand) / one_hundred;
  int day   = date - year * ten_thousand - month * one_hundred;
  return util::TimeStamp(year, month, day, 0, 0, 0);
}
// FIXME: check if this function is implemented in eamxx
// Assumes 365 days/year, 30 days/month
inline int compute_days(const util::TimeStamp &ts) {
  return ts.get_year() * 365 + ts.get_month() * 30 + ts.get_day();
}

inline void create_linoz_chlorine_reader(
    const std::string &linoz_chlorine_file, const util::TimeStamp &model_time,
    const int chlorine_loading_ymd,  // in format YYYYMMDD
    std::vector<Real> &values, std::vector<int> &time_secs) {
  auto time_stamp_beg = convert_date(chlorine_loading_ymd);

  const int offset_time =
      compute_days(time_stamp_beg) - compute_days(model_time);
  scorpio::register_file(linoz_chlorine_file, scorpio::Read);
  const int nlevs_time = scorpio::get_time_len(linoz_chlorine_file);
  for(int itime = 0; itime < nlevs_time; ++itime) {
    int date;
    scorpio::read_var(linoz_chlorine_file, "date", &date, itime);
    if(date >= chlorine_loading_ymd) {
      Real value;
      scorpio::read_var(linoz_chlorine_file, "chlorine_loading", &value, itime);
      values.push_back(value);
      auto time_stamp = convert_date(date);
      time_secs.push_back(compute_days(time_stamp) - offset_time);
    }
  }  // end itime
  scorpio::release_file(linoz_chlorine_file);
}

inline Real chlorine_loading_advance(const util::TimeStamp &ts,
                                     std::vector<Real> &values,
                                     std::vector<int> &time_secs) {
  const int current_time = compute_days(ts);
  int index              = 0;
  // update index
  for(int i = 0; i < int(values.size()); i++) {
    if(current_time > time_secs[i]) {
      index = i;
      break;
    }
  }  //

  const Real delt = (current_time - time_secs[index]) /
                    (time_secs[index + 1] - time_secs[index]);
  return values[index] + delt * (values[index + 1] - values[index]);
}

inline std::shared_ptr<AbstractRemapper> create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &trace_data_file, const std::string &map_file,
    const std::vector<std::string> &var_names,
    TracerFileType &tracer_file_type) {
  using namespace ShortFieldTagsNames;

  scorpio::register_file(trace_data_file, scorpio::Read);

  // by default, I am assuming a zonal file.
  tracer_file_type = ZONAL;

  int nlevs_data = -1;
  if(scorpio::has_var(trace_data_file, "lev")) {
    nlevs_data = scorpio::get_dimlen(trace_data_file, "lev");
  }
  const bool has_altitude = scorpio::has_var(trace_data_file, "altitude");

  // This type of files use altitude (zi) for vertical interpolation
  if(has_altitude) {
    nlevs_data       = scorpio::get_dimlen(trace_data_file, "altitude");
    tracer_file_type = VERT_EMISSION;
  }

  EKAT_REQUIRE_MSG(
      nlevs_data != -1,
      "Error: The file does not contain either lev or altitude.   \n");

  const int ncols_data = scorpio::get_dimlen(trace_data_file, "ncol");

  // This type of files use model pressure (pmid) for vertical interpolation
  if(scorpio::has_var(trace_data_file, "PS")) {
    tracer_file_type = FORMULA_PS;
  }

  scorpio::release_file(trace_data_file);

  // We could use model_grid directly if using same num levels,
  // but since shallow clones are cheap, we may as well do it (less lines of
  // code)
  auto horiz_interp_tgt_grid =
      model_grid->clone("tracer_horiz_interp_tgt_grid", true);
  horiz_interp_tgt_grid->reset_num_vertical_lev(nlevs_data);

  if(has_altitude) {
    horiz_interp_tgt_grid->reset_field_tag_name(LEV, "altitude");
    horiz_interp_tgt_grid->reset_field_tag_name(ILEV, "altitude_int");
  }

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  if(ncols_data == ncols_model) {
    remapper = std::make_shared<IdentityRemapper>(
        horiz_interp_tgt_grid, IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG(ncols_data <= ncols_model,
                     "Error! We do not allow to coarsen spa data to fit the "
                     "model. We only allow\n"
                     "       spa data to be at the same or coarser resolution "
                     "as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG(
        map_file != "",
        "ERROR: Spa data is on a different grid than the model one,\n"
        "       but spa_remap_file is missing from SPA parameter list.");

    remapper =
        std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid, map_file);
  }

  remapper->registration_begins();
  const auto tgt_grid  = remapper->get_tgt_grid();
  const auto layout_2d = tgt_grid->get_2d_scalar_layout();

  const auto layout_3d_mid = tgt_grid->get_3d_scalar_layout(true);
  // FieldLayout  layout_3d_mid;
  // if ( has_altitude ) {

  //   auto make_layout = [](const std::vector<int>& extents,
  //                       const std::vector<std::string>& names)
  //   {
  //     std::vector<FieldTag> tags(extents.size(),CMP);
  //     return FieldLayout(tags,extents,names);
  //   };
  //   layout_3d_mid = make_layout({ncols_model, nlevs_data},
  //                               {"ncol","altitude"});
  //   // FieldLayout({FieldTag::Column,CMP},{ncols_model,nlevs_data});
  // } else {
  //   layout_3d_mid = tgt_grid->get_3d_scalar_layout(true);
  // }

  const auto nondim = ekat::units::Units::nondimensional();

  for(auto var_name : var_names) {
    Field ifield(
        FieldIdentifier(var_name, layout_3d_mid, nondim, tgt_grid->name()));
    ifield.allocate_view();
    remapper->register_field_from_tgt(ifield);
  }
  // zonal files do not have the PS variable.
  if(tracer_file_type == FORMULA_PS) {
    Field ps(FieldIdentifier("PS", layout_2d, nondim, tgt_grid->name()));
    ps.allocate_view();
    remapper->register_field_from_tgt(ps);
  }
  remapper->registration_ends();
  return remapper;

}  // create_horiz_remapper

inline std::shared_ptr<AtmosphereInput> create_tracer_data_reader(
    const std::shared_ptr<AbstractRemapper> &horiz_remapper,
    const std::string &tracer_data_file) {
  std::vector<Field> io_fields;
  for(int i = 0; i < horiz_remapper->get_num_fields(); ++i) {
    io_fields.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(tracer_data_file, io_grid, io_fields,
                                           true);
}  // create_tracer_data_reader

inline void update_tracer_data_from_file(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    const int time_index,  // zero-based
    AbstractRemapper &tracer_horiz_interp, TracerData &tracer_data) {
  // 1. read from field
  scorpio_reader->read_variables(time_index);
  // 2. Run the horiz remapper (it is a do-nothing op if spa data is on same
  // grid as model)
  tracer_horiz_interp.remap(/*forward = */ true);
  //
  const int nvars = tracer_data.nvars_;
  //
  for(int i = 0; i < nvars; ++i) {
    tracer_data.data[i] =
        tracer_horiz_interp.get_tgt_field(i).get_view<Real **>();
  }

  if(tracer_data.file_type == FORMULA_PS) {
    // Recall, the fields are registered in the order: tracers, ps
    // 3. Copy from the tgt field of the remapper into the spa_data
    tracer_data.ps =
        tracer_horiz_interp.get_tgt_field(nvars).get_view<Real *>();
  }

}  // update_tracer_data_from_file
inline void update_tracer_timestate(
    std::shared_ptr<AtmosphereInput> &scorpio_reader, const util::TimeStamp &ts,
    AbstractRemapper &tracer_horiz_interp, TracerTimeState &time_state,
    TracerData &data_tracer_beg, TracerData &data_tracer_end) {
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that SPA assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1;  // Make it 0-based
  if(month != time_state.current_month) {
    //
    const auto tracer_beg = data_tracer_beg.data;
    const auto tracer_end = data_tracer_end.data;
    const int nvars       = data_tracer_end.nvars_;

    // Update the SPA time state information
    time_state.current_month = month;
    time_state.t_beg_month =
        util::TimeStamp({ts.get_year(), month + 1, 1}, {0, 0, 0})
            .frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(), month + 1);

    // Copy spa_end'data into spa_beg'data, and read in the new spa_end
    for(int ivar = 0; ivar < nvars; ++ivar) {
      Kokkos::deep_copy(tracer_beg[ivar], tracer_end[ivar]);
    }
    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong
    // values
    //       to be assigned.  A timestep greater than a month is very unlikely
    //       so we will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    update_tracer_data_from_file(scorpio_reader, ts, next_month,
                                 tracer_horiz_interp, data_tracer_end);
  }

}  // END updata_spa_timestate

// This function is based on the SPA::perform_time_interpolation function.
inline void perform_time_interpolation(const TracerTimeState &time_state,
                                       const TracerData &data_tracer_beg,
                                       const TracerData &data_tracer_end,
                                       const TracerData &data_tracer_out) {
  // NOTE: we *assume* data_beg and data_end have the *same* hybrid v coords.
  //       IF this ever ceases to be the case, you can interp those too.
  // Gather time stamp info
  auto &t_now   = time_state.t_now;
  auto &t_beg   = time_state.t_beg_month;
  auto &delta_t = time_state.days_this_month;

  // We can ||ize over columns as well as over variables and bands
  const auto data_beg = data_tracer_beg.data;
  const auto data_end = data_tracer_end.data;
  const auto data_out = data_tracer_out.data;

  const auto file_type = data_tracer_out.file_type;

  const auto ps_beg = data_tracer_beg.ps;
  const auto ps_end = data_tracer_end.ps;
  const auto ps_out = data_tracer_out.ps;

  const int num_vars = data_tracer_end.nvars_;

  const int ncol     = data_tracer_beg.ncol_;
  const int num_vert = data_tracer_beg.nlev_;

  const int outer_iters = ncol * num_vars;

  const auto policy = ESU::get_default_team_policy(outer_iters, num_vert);

  auto delta_t_fraction = (t_now - t_beg) / delta_t;

  EKAT_REQUIRE_MSG(
      delta_t_fraction >= 0 && delta_t_fraction <= 1,
      "Error! Convex interpolation with coefficient out of [0,1].\n"
      "  t_now  : " +
          std::to_string(t_now) +
          "\n"
          "  t_beg  : " +
          std::to_string(t_beg) +
          "\n"
          "  delta_t: " +
          std::to_string(delta_t) + "\n");

  Kokkos::parallel_for(
      "linoz_time_interp_loop", policy, KOKKOS_LAMBDA(const Team &team) {
        // The policy is over ncols*num_vars, so retrieve icol/ivar
        const int icol = team.league_rank() / num_vars;
        const int ivar = team.league_rank() % num_vars;

        // Get column of beg/end/out variable
        auto var_beg = ekat::subview(data_beg[ivar], icol);
        auto var_end = ekat::subview(data_end[ivar], icol);
        auto var_out = ekat::subview(data_out[ivar], icol);

        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, num_vert), [&](const int &k) {
              var_out(k) =
                  linear_interp(var_beg(k), var_end(k), delta_t_fraction);
            });
        // linoz files do not have ps variables.
        if(ivar == 1 && file_type == FORMULA_PS) {
          ps_out(icol) =
              linear_interp(ps_beg(icol), ps_end(icol), delta_t_fraction);
        }
      });
  Kokkos::fence();
}  // perform_time_interpolation

inline void compute_source_pressure_levels(const view_1d &ps_src,
                                           const view_2d &p_src,
                                           const const_view_1d &hyam,
                                           const const_view_1d &hybm) {
  constexpr auto P0        = C::P0;
  const int ncols          = ps_src.extent(0);
  const int num_vert_packs = p_src.extent(1);
  const auto policy = ESU::get_default_team_policy(ncols, num_vert_packs);

  Kokkos::parallel_for(
      "tracer_compute_p_src_loop", policy, KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, num_vert_packs), [&](const int k) {
              p_src(icol, k) = ps_src(icol) * hybm(k) + P0 * hyam(k);
            });
      });
}  // compute_source_pressure_levels

// Linoz NetCDF files use levs instead of formula_terms.
// This function allocates a view, so we need to do it during initialization.
// Thus, we assume that source pressure is independent of time,
// which is the case for Linoz files (zonal file).
inline void compute_p_src_zonal_files(const std::string &tracer_file_name,
                                      const view_2d &p_src) {
  EKAT_REQUIRE_MSG(p_src.data() != 0,
                   "Error: p_src has not been allocated. \n");
  // Read in levs in start/end data
  // FIXME: units are mbar; how can I get units using scorpio interface
  auto nondim = ekat::units::Units::nondimensional();
  scorpio::register_file(tracer_file_name, scorpio::Read);
  const int nlevs_data = scorpio::get_dimlen(tracer_file_name, "lev");
  view_1d_host levs_h("levs_h", nlevs_data);
  scorpio::read_var(tracer_file_name, "lev", levs_h.data());
  scorpio::release_file(tracer_file_name);
  view_1d levs("levs", nlevs_data);
  Kokkos::deep_copy(levs, levs_h);

  const int ncol = p_src.extent(0);
  EKAT_REQUIRE_MSG(
      p_src.extent(1) == nlevs_data,
      "Error: p_src has a different number of levels than the source data. \n");

  const auto policy_pressure = ESU::get_default_team_policy(ncol, nlevs_data);
  const int pi               = haero::Constants::pi;
  Kokkos::parallel_for(
      "pressure_computation", policy_pressure, KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs_data),
                             [&](const Int &kk) {
                               // mbar->pascals
                               // FIXME: Does EAMxx have a better method to
                               // convert units?"
                               p_src(icol, kk) = levs(kk) * 100;
                             });
      });
  Kokkos::fence();
}

inline void perform_vertical_interpolation(const view_2d &p_src_c,
                                           const const_view_2d &p_tgt_c,
                                           const TracerData &input,
                                           const view_2d output[]) {
  // At this stage, begin/end must have the same horiz dimensions
  EKAT_REQUIRE(input.ncol_ == output[0].extent(0));

#if 1
  // FIXME: I was encountering a compilation error when using const_view_2d.
  // The issue is fixed by https://github.com/E3SM-Project/EKAT/pull/346.
  // I will keep this code until this PR is merged into the EKAT master branch
  // and we update the EKAT version in our code. I am converting const_view_2d
  // to view_2d.
  auto p_src_ptr = (Real *)p_src_c.data();
  view_2d p_src(p_src_ptr, p_src_c.extent(0), p_src_c.extent(1));
  auto p_tgt_ptr = (Real *)p_tgt_c.data();
  view_2d p_tgt(p_tgt_ptr, p_tgt_c.extent(0), p_tgt_c.extent(1));
#else
  const auto p_src = p_src_c;
  const auto p_tgt = p_tgt_c;
#endif

  const int ncols = input.ncol_;
  // FIXME: I am getting FPEs if I do not subtract 1 from nlevs_src.
  const int nlevs_src = input.nlev_ - 1;
  const int nlevs_tgt = output[0].extent(1);

  LIV vert_interp(ncols, nlevs_src, nlevs_tgt);

  // We can ||ize over columns as well as over variables and bands
  const int num_vars       = input.nvars_;
  const int num_vert_packs = nlevs_tgt;
  const auto policy_setup = ESU::get_default_team_policy(ncols, num_vert_packs);

  // Setup the linear interpolation object
  Kokkos::parallel_for(
      "tracer_vert_interp_setup_loop", policy_setup,
      KOKKOS_LAMBDA(typename LIV::MemberType const &team) {
        const int icol = team.league_rank();
        // Setup
        vert_interp.setup(team, ekat::subview(p_src, icol),
                          ekat::subview(p_tgt, icol));
      });
  Kokkos::fence();

  // Now use the interpolation object in || over all variables.
  const int outer_iters = ncols * num_vars;
  const auto policy_interp =
      ESU::get_default_team_policy(outer_iters, num_vert_packs);
  Kokkos::parallel_for(
      "tracer_vert_interp_loop", policy_interp,
      KOKKOS_LAMBDA(typename LIV::MemberType const &team) {
        const int icol = team.league_rank() / num_vars;
        const int ivar = team.league_rank() % num_vars;

        const auto x1 = ekat::subview(p_src, icol);
        const auto x2 = ekat::subview(p_tgt, icol);

        const auto y1 = ekat::subview(input.data[ivar], icol);
        const auto y2 = ekat::subview(output[ivar], icol);

        vert_interp.lin_interp(team, x1, x2, y1, y2, icol);
      });
  Kokkos::fence();
}

// rebin is a port from:
// https://github.com/eagles-project/e3sm_mam4_refactor/blob/ee556e13762e41a82cb70a240c54dc1b1e313621/components/eam/src/chemistry/utils/mo_util.F90#L12
inline void rebin(int nsrc, int ntrg, const const_view_1d &src_x,
                  const Real trg_x[], const view_1d &src, const view_1d &trg) {
  for(int i = 0; i < ntrg; ++i) {
    Real tl = trg_x[i];
    if(tl < src_x(nsrc)) {
      int sil = 0;
      for(; sil <= nsrc; ++sil) {
        if(tl <= src_x(sil)) {
          break;
        }
      }
      Real tu = trg_x[i + 1];
      int siu = 0;
      for(; siu <= nsrc; ++siu) {
        if(tu <= src_x(siu)) {
          break;
        }
      }
      Real y = 0.0;
      sil    = haero::max(sil, 1);
      siu    = haero::min(siu, nsrc);
      for(int si = sil; si <= siu; ++si) {
        int si1 = si - 1;
        Real sl = haero::max(tl, src_x(si1));
        Real su = haero::min(tu, src_x(si));
        y += (su - sl) * src(si1);
      }
      trg(i) = y / (trg_x[i + 1] - trg_x[i]);
    } else {
      trg(i) = 0.0;
    }
  }
}  // rebin

inline void perform_vertical_interpolation(const const_view_1d &altitude_int,
                                           const const_view_2d &zi,
                                           const TracerData &input,
                                           const view_2d output[]) {
  EKAT_REQUIRE_MSG(
      input.file_type == VERT_EMISSION,
      "Error! vertical interpolation only with altitude variable. \n");
  const int ncols          = input.ncol_;
  const int num_vars       = input.nvars_;
  const int ntrg           = output[0].extent(1);
  const int num_vert_packs = ntrg;
  const int outer_iters    = ncols * num_vars;
  const auto policy_interp =
      ESU::get_default_team_policy(outer_iters, num_vert_packs);
  // FIXME: Get m2km from emaxx.
  const Real m2km    = 1e-3;
  const auto &src_x  = altitude_int;
  const int nsrc     = input.nlev_;
  constexpr int pver = mam4::nlev;
  const int pverp    = pver + 1;

  Kokkos::parallel_for(
      "tracer_vert_interp_loop", policy_interp,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank() / num_vars;
        const int ivar = team.league_rank() % num_vars;

        const auto src = ekat::subview(input.data[ivar], icol);
        const auto trg = ekat::subview(output[ivar], icol);

        // trg_x
        Real trg_x[pver + 1];
        // I am trying to do this:
        // model_z(1:pverp) = m2km * state(c)%zi(i,pverp:1:-1)
        for(int i = 0; i < pverp; ++i) {
          trg_x[pverp - i - 1] = m2km * zi(icol, i);
        }
        team.team_barrier();

        rebin(nsrc, ntrg, src_x, trg_x, src, trg);
      });
}

inline void advance_tracer_data(
    std::shared_ptr<AtmosphereInput> &scorpio_reader,
    AbstractRemapper &tracer_horiz_interp, const util::TimeStamp &ts,
    TracerTimeState &time_state, TracerData &data_tracer_beg,
    TracerData &data_tracer_end, TracerData &data_tracer_out,
    const view_2d &p_src, const const_view_2d &p_tgt,
    const const_view_1d &zi_src, const const_view_2d &zi_tgt,
    const view_2d output[]) {
  /* Update the TracerTimeState to reflect the current time, note the addition
   * of dt */
  time_state.t_now = ts.frac_of_year_in_days();
  /* Update time state and if the month has changed, update the data.*/
  update_tracer_timestate(scorpio_reader, ts, tracer_horiz_interp, time_state,
                          data_tracer_beg, data_tracer_end);
  // Step 1. Perform time interpolation
  perform_time_interpolation(time_state, data_tracer_beg, data_tracer_end,
                             data_tracer_out);

  if(data_tracer_out.file_type == FORMULA_PS) {
    // Step 2. Compute source pressure levels
    compute_source_pressure_levels(data_tracer_out.ps, p_src,
                                   data_tracer_out.hyam, data_tracer_out.hybm);
  }

  // Step 3. Perform vertical interpolation
  if(data_tracer_out.file_type == FORMULA_PS ||
     data_tracer_out.file_type == ZONAL) {
    perform_vertical_interpolation(p_src, p_tgt, data_tracer_out, output);
  } else if(data_tracer_out.file_type == VERT_EMISSION) {
    perform_vertical_interpolation(zi_src, zi_tgt, data_tracer_out, output);
  }

}  // advance_tracer_data

}  // namespace scream::mam_coupling
#endif  // EAMXX_MAM_HELPER_MICRO
