#ifndef EAMXX_MAM_TRACER_READER_UTILS
#define EAMXX_MAM_TRACER_READER_UTILS

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
using view_2d_host = typename KT::view_2d<Real>::HostMirror;

using ExeSpace = typename KT::ExeSpace;
using ESU      = ekat::ExeSpaceUtils<ExeSpace>;
using C        = scream::physics::Constants<Real>;
using LIV      = ekat::LinInterp<Real, 1>;

// Linoz NetCDF files use levs instead of formula_terms.
// This function allocates a view, so we need to do it during initialization.
// Thus, we assume that source pressure is independent of time,
// which is the case for Linoz files (zonal file).

inline void compute_p_src_zonal_files(const view_1d &levs,
                                      const view_2d &p_src) {
  EKAT_REQUIRE_MSG(p_src.data() != 0,
                   "Error: p_src has not been allocated. \n");
  EKAT_REQUIRE_MSG(levs.data() != 0, "Error: levs has not been allocated. \n");
  const int ncol       = p_src.extent(0);
  const int nlevs_data = levs.extent(0);
  EKAT_REQUIRE_MSG(
      p_src.extent_int(1) == nlevs_data,
      "Error: p_src has a different number of levels than the source data. \n");
  Kokkos::parallel_for(
      "pressure_computation",
      Kokkos::MDRangePolicy<Kokkos::Rank<2> >({0, 0}, {ncol, nlevs_data}),
      KOKKOS_LAMBDA(const int icol, const int kk) {
        // mbar->pascals
        // FIXME: Does EAMxx have a better method to
        // convert units?"
        p_src(icol, kk) = levs(kk) * 100;
      });
  Kokkos::fence();
}

// We have a similar version in MAM4xx.
// This version was created because the data view cannot be modified
// inside the parallel_for.
// This struct will be used in init while reading nc files.
// The MAM4xx version will be used instead of parallel_for that loops over cols.
struct ForcingHelper {
  // This index is in Fortran format. i.e. starts in 1
  int frc_ndx;
  // does this file have altitude?
  bool file_alt_data;
  // number of sectors per forcing
  int nsectors;
  // offset in output vector from reader
  int offset;
};

enum TracerFileType {
  // file with PS ncol, lev, and time
  FORMULA_PS,
  // nc zonal file from ncremap
  ZONAL,
  // vertical emission files
  VERT_EMISSION,
};

enum TracerDataIndex { BEG = 0, END = 1, OUT = 2 };

/* Maximum number of tracers (or fields) that the tracer reader can handle.
 Note: We are not allocating memory for MAX_NVARS_TRACER tracers.
 Therefore, if a file contains more than this number, it is acceptable to
 increase this limit. Currently, Linoz files have 8 fields. */
constexpr int MAX_NVARS_TRACER             = 10;
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
  int nlevs_data;
  int ncols_data;

  //
  int offset_time_index_{0};

  // We cannot use a std::vector<view_2d>
  // because we need to access these views from device.
  // 0: beg 1: end 3: out
  view_2d data[3][MAX_NVARS_TRACER];
  // type of file
  TracerFileType file_type;

  // These views are employed in files with the PS variable
  // 0: beg 1: end 3: out
  view_1d ps[3];
  const_view_1d hyam;
  const_view_1d hybm;
  view_int_1d work_vert_inter[MAX_NVARS_TRACER];

  // External forcing file (vertical emission)
  // Uses altitude instead of pressure to interpolate data
  view_1d altitude_int_;
  //
  bool has_altitude_{false};
  // pressure source for ZONAL or FORMULA_PS files
  view_2d p_src_;

  // only for zonal files
  view_1d zonal_levs_;

  void allocate_temporal_views() {
    // BEG and OUT data views.
    EKAT_REQUIRE_MSG(ncol_ != int(-1), "Error! ncols has not been set. \n");
    EKAT_REQUIRE_MSG(nlev_ != int(-1), "Error! nlevs has not been set. \n");
    EKAT_REQUIRE_MSG(nvars_ != int(-1), "Error! nvars has not been set. \n");

    for(int ivar = 0; ivar < nvars_; ++ivar) {
      data[TracerDataIndex::OUT][ivar] = view_2d("linoz_1_out", ncol_, nlev_);
      data[TracerDataIndex::BEG][ivar] = view_2d("linoz_1_out", ncol_, nlev_);
    }

    // for vertical interpolation using rebin routine
    if(file_type == FORMULA_PS || file_type == ZONAL) {
      // we only need work array for FORMULA_PS or ZONAL
      for(int ivar = 0; ivar < nvars_; ++ivar) {
        work_vert_inter[ivar] =
            view_int_1d("allocate_work_vertical_interpolation", ncol_);
      }
      // we use ncremap and python scripts to convert zonal files to ne4pn4
      // grids.
      p_src_ = view_2d("pressure_src_invariant", ncol_, nlev_);
    }

    if(file_type == TracerFileType::FORMULA_PS) {
      ps[TracerDataIndex::OUT] = view_1d("ps", ncol_);
      ps[TracerDataIndex::BEG] = view_1d("ps", ncol_);
    }

    if(file_type == TracerFileType::ZONAL) {
      // we use ncremap and python scripts to convert zonal files to ne4pn4
      // grids.
      compute_p_src_zonal_files(zonal_levs_, p_src_);
    }
  }
};

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
// FIXME: This function is not implemented in eamxx.
// FIXME: Assumes 365 days/year, 30 days/month;
// NOTE: that this assumption is mainly used for plotting.
// NOTE: This is not a direct port from EAM.
//We only use this routine for chlorine.
inline int compute_number_days_from_zero(const util::TimeStamp &ts) {
  return ts.get_year() * 365 + ts.get_month() * 30 + ts.get_day();
}

inline void create_linoz_chlorine_reader(
    const std::string &linoz_chlorine_file, const util::TimeStamp &model_time,
    const int chlorine_loading_ymd,  // in format YYYYMMDD
    std::vector<Real> &values, std::vector<int> &time_secs) {
  auto time_stamp_beg = convert_date(chlorine_loading_ymd);

  const int offset_time =
      compute_number_days_from_zero(time_stamp_beg) - compute_number_days_from_zero(model_time);
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
      time_secs.push_back(compute_number_days_from_zero(time_stamp) - offset_time);
    }
  }  // end itime
  scorpio::release_file(linoz_chlorine_file);
}

// Gets the times from the NC file
// Given a date in the format YYYYMMDD, returns its index in the time dimension.
inline void get_time_from_ncfile(const std::string &file_name,
                                 const int cyclical_ymd,  // in format YYYYMMDD
                                 int &cyclical_ymd_index,
                                 std::vector<int> &dates) {
  // in file_name: name of the NC file
  // in cyclical_ymd: date in the format YYYYMMDD
  // out cyclical_ymd_index: time index for cyclical_ymd
  // out dates: date in YYYYMMDD format
  scorpio::register_file(file_name, scorpio::Read);
  const int nlevs_time = scorpio::get_time_len(file_name);
  cyclical_ymd_index   = -1;
  for(int itime = 0; itime < nlevs_time; ++itime) {
    int date;
    scorpio::read_var(file_name, "date", &date, itime);
    // std::cout << itime << " date: " << date << "\n";
    if(date >= cyclical_ymd && cyclical_ymd_index == -1) {
      cyclical_ymd_index = itime;
    }
    dates.push_back(date);
  }  // end itime

  EKAT_REQUIRE_MSG(cyclical_ymd_index >= 0,
                   "Error! Current model time (" +
                       std::to_string(cyclical_ymd) + ") is not within " +
                       "Tracer time period: [" + std::to_string(dates[0]) +
                       ", " + "(" + std::to_string(dates[nlevs_time - 1]) +
                       ").\n");
  scorpio::release_file(file_name);
}

inline Real chlorine_loading_advance(const util::TimeStamp &ts,
                                     std::vector<Real> &values,
                                     std::vector<int> &time_secs) {
  const int current_time = compute_number_days_from_zero(ts);
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

// It reads variables that are not time-dependent and independent of columns (no
// MPI involved here). We also obtain the offset_time_index using a date
// (cyclical_ymd) as input. We initialize a few members of tracer_data.
inline void setup_tracer_data(TracerData &tracer_data,             // out
                              const std::string &trace_data_file,  // in
                              const int cyclical_ymd)              // in
{
  scorpio::register_file(trace_data_file, scorpio::Read);
  // by default, I am assuming a zonal file.
  TracerFileType tracer_file_type = ZONAL;

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
    view_1d_host hyam_h("hyam_h", nlevs_data);
    view_1d_host hybm_h("hybm_h", nlevs_data);

    scorpio::read_var(trace_data_file, "hyam", hyam_h.data());
    scorpio::read_var(trace_data_file, "hybm", hybm_h.data());
    view_1d hyam("hyam", nlevs_data);
    view_1d hybm("hybm", nlevs_data);
    Kokkos::deep_copy(hyam, hyam_h);
    Kokkos::deep_copy(hybm, hybm_h);
    tracer_data.hyam = hyam;
    tracer_data.hybm = hybm;
  }

  if(tracer_file_type == ZONAL) {
    view_1d_host levs_h("levs_h", nlevs_data);
    view_1d levs("levs", nlevs_data);
    scorpio::read_var(trace_data_file, "lev", levs_h.data());
    Kokkos::deep_copy(levs, levs_h);
    tracer_data.zonal_levs_ = levs;
  }

  if(tracer_file_type == VERT_EMISSION) {
    const int nilevs_data =
        scorpio::get_dimlen(trace_data_file, "altitude_int");
    view_1d_host altitude_int_host("altitude_int_host", nilevs_data);
    view_1d altitude_int = view_1d("altitude_int", nilevs_data);
    scorpio::read_var(trace_data_file, "altitude_int",
                      altitude_int_host.data());
    Kokkos::deep_copy(altitude_int, altitude_int_host);
    tracer_data.altitude_int_ = altitude_int;
  }
  // time index
  {
    const int nlevs_time   = scorpio::get_dimlen(trace_data_file, "time");
    int cyclical_ymd_index = -1;
    for(int itime = 0; itime < nlevs_time; ++itime) {
      int date;
      scorpio::read_var(trace_data_file, "date", &date, itime);
      if(date >= cyclical_ymd) {
        cyclical_ymd_index = itime;
        break;
      }
    }  // end itime

    EKAT_REQUIRE_MSG(cyclical_ymd_index >= 0, "Error! Current model time (" +
                                                  std::to_string(cyclical_ymd) +
                                                  ") is not within " +
                                                  "Tracer time period.\n");

    tracer_data.offset_time_index_ = cyclical_ymd_index;
  }

  scorpio::release_file(trace_data_file);
  tracer_data.file_type     = tracer_file_type;
  tracer_data.nlevs_data    = nlevs_data;
  tracer_data.ncols_data    = ncols_data;
  tracer_data.has_altitude_ = has_altitude;
}
inline std::shared_ptr<AbstractRemapper> create_horiz_remapper(
    const std::shared_ptr<const AbstractGrid> &model_grid,
    const std::string &trace_data_file, const std::string &map_file,
    const std::vector<std::string> &var_names, TracerData &tracer_data) {
  using namespace ShortFieldTagsNames;
  // We could use model_grid directly if using same num levels,
  // but since shallow clones are cheap, we may as well do it (less lines of
  // code)
  auto horiz_interp_tgt_grid =
      model_grid->clone("tracer_horiz_interp_tgt_grid", true);
  horiz_interp_tgt_grid->reset_num_vertical_lev(tracer_data.nlevs_data);

  if(tracer_data.file_type == VERT_EMISSION) {
    horiz_interp_tgt_grid->reset_field_tag_name(LEV, "altitude");
    horiz_interp_tgt_grid->reset_field_tag_name(ILEV, "altitude_int");
  }

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  if(tracer_data.ncols_data == ncols_model) {
    remapper = std::make_shared<IdentityRemapper>(
        horiz_interp_tgt_grid, IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG(tracer_data.ncols_data <= ncols_model,
                     "Error! We do not allow to coarsen tracer external "
                     "forcing data to fit the "
                     "model. We only allow\n"
                     "       tracer external forcing data to be at the same or "
                     "coarser resolution "
                     "as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG(
        map_file != "",
        "ERROR: tracer external forcing data is on a different grid than the "
        "model one,\n"
        "       but tracer external forcing data remap file is missing from "
        "tracer external forcing data parameter list.");

    remapper =
        std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid, map_file);
  }

  remapper->registration_begins();
  const auto tgt_grid  = remapper->get_tgt_grid();
  const auto layout_2d = tgt_grid->get_2d_scalar_layout();

  const auto layout_3d_mid = tgt_grid->get_3d_scalar_layout(true);

  const auto nondim = ekat::units::Units::nondimensional();

  for(auto var_name : var_names) {
    Field ifield(
        FieldIdentifier(var_name, layout_3d_mid, nondim, tgt_grid->name()));
    ifield.allocate_view();
    remapper->register_field_from_tgt(ifield);
  }
  // zonal files do not have the PS variable.
  if(tracer_data.file_type == FORMULA_PS) {
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
    const std::shared_ptr<AtmosphereInput> &scorpio_reader,
    const int time_index,  // zero-based
    AbstractRemapper &tracer_horiz_interp, TracerData &tracer_data) {
  // 1. read from field
  scorpio_reader->read_variables(time_index);
  // 2. Run the horiz remapper (it is a do-nothing op if tracer external forcing
  // data is on same grid as model)
  tracer_horiz_interp.remap(/*forward = */ true);
  //
  const int nvars = tracer_data.nvars_;
  //
  for(int i = 0; i < nvars; ++i) {
    tracer_data.data[TracerDataIndex::END][i] =
        tracer_horiz_interp.get_tgt_field(i).get_view<Real **>();
  }

  if(tracer_data.file_type == FORMULA_PS) {
    // Recall, the fields are registered in the order: tracers, ps
    // 3. Copy from the tgt field of the remapper into the spa_data
    tracer_data.ps[TracerDataIndex::END] =
        tracer_horiz_interp.get_tgt_field(nvars).get_view<Real *>();
  }

}  // update_tracer_data_from_file
inline void update_tracer_timestate(
    const std::shared_ptr<AtmosphereInput> &scorpio_reader,
    const util::TimeStamp &ts, AbstractRemapper &tracer_horiz_interp,
    TracerTimeState &time_state, TracerData &data_tracer) {
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that tracer external forcing assumes monthly data to
  // update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1;  // Make it 0-based
  if(month != time_state.current_month) {
    const auto tracer_data = data_tracer.data;
    const int nvars        = data_tracer.nvars_;
    const auto ps          = data_tracer.ps;

    // Update the tracer external forcing time state information
    time_state.current_month = month;
    time_state.t_beg_month =
        util::TimeStamp({ts.get_year(), month + 1, 1}, {0, 0, 0})
            .frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(), month + 1);

    // Copy spa_end'data into spa_beg'data, and read in the new spa_end
    for(int ivar = 0; ivar < nvars; ++ivar) {
      Kokkos::deep_copy(tracer_data[TracerDataIndex::BEG][ivar],
                        tracer_data[TracerDataIndex::END][ivar]);
    }

    if(data_tracer.file_type == FORMULA_PS) {
      Kokkos::deep_copy(ps[TracerDataIndex::BEG], ps[TracerDataIndex::END]);
    }

    // Following SPA to time-interpolate data in MAM4xx
    // Assume the data is saved monthly and cycles in one year
    // Add offset_time_index to support cases where data is saved
    // from other periods of time.
    // Update the tracer external forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong
    // values
    //       to be assigned.  A timestep greater than a month is very unlikely
    //       so we will proceed.
    int next_month =
        data_tracer.offset_time_index_ + (time_state.current_month + 1) % 12;
    update_tracer_data_from_file(scorpio_reader, next_month,
                                 tracer_horiz_interp, data_tracer);
  }

}  // END update_tracer_timestate

// This function is based on the SPA::perform_time_interpolation function.
inline void perform_time_interpolation(const TracerTimeState &time_state,
                                       const TracerData &data_tracer) {
  // NOTE: we *assume* data_beg and data_end have the *same* hybrid v coords.
  //       IF this ever ceases to be the case, you can interp those too.
  // Gather time stamp info
  auto &t_now   = time_state.t_now;
  auto &t_beg   = time_state.t_beg_month;
  auto &delta_t = time_state.days_this_month;

  // We can ||ize over columns as well as over variables and bands
  const auto& data      = data_tracer.data;
  const auto& file_type = data_tracer.file_type;

  const auto& ps      = data_tracer.ps;
  const int num_vars = data_tracer.nvars_;

  const int ncol     = data_tracer.ncol_;
  const int num_vert = data_tracer.nlev_;

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
        auto var_beg = ekat::subview(data[TracerDataIndex::BEG][ivar], icol);
        auto var_end = ekat::subview(data[TracerDataIndex::END][ivar], icol);
        auto var_out = ekat::subview(data[TracerDataIndex::OUT][ivar], icol);

        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, num_vert), [&](const int &k) {
              var_out(k) =
                  linear_interp(var_beg(k), var_end(k), delta_t_fraction);
            });
        // linoz files do not have ps variables.
        if(ivar == 1 && file_type == FORMULA_PS) {
          ps[TracerDataIndex::OUT](icol) =
              linear_interp(ps[TracerDataIndex::BEG](icol),
                            ps[TracerDataIndex::END](icol), delta_t_fraction);
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

inline void perform_vertical_interpolation(const view_2d &p_src_c,
                                           const const_view_2d &p_tgt_c,
                                           const TracerData &input,
                                           const view_2d output[]) {
  const int ncol   = input.ncol_;
  const int levsiz = input.nlev_;
  const int pver   = mam4::nlev;

  const int num_vars = input.nvars_;
  // make a local copy of output
  view_2d output_local[MAX_NVARS_TRACER];
  EKAT_REQUIRE_MSG(
        num_vars <= int(MAX_NVARS_TRACER),
        "Error! Number of variables is bigger than NVARS_MAXTRACER. \n");
  for (int ivar = 0; ivar < num_vars; ++ivar)
  {
    // At this stage, begin/end must have the same horiz dimensions
    EKAT_REQUIRE(input.ncol_ == output[ivar].extent_int(0));
    output_local[ivar] = output[ivar];
  }
  const int outer_iters   = ncol * num_vars;
  const auto policy_setup = ESU::get_default_team_policy(outer_iters, pver);
  const auto& data         = input.data;

  Kokkos::parallel_for(
      "vert_interp", policy_setup,
      KOKKOS_LAMBDA(typename LIV::MemberType const &team) {
        // The policy is over ncols*num_vars, so retrieve icol/ivar
        const int icol             = team.league_rank() / num_vars;
        const int ivar             = team.league_rank() % num_vars;
        const auto pin_at_icol     = ekat::subview(p_src_c, icol);
        const auto pmid_at_icol    = ekat::subview(p_tgt_c, icol);
        const auto& datain          = data[TracerDataIndex::OUT][ivar];
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout = output_local[ivar];
        const auto dataout_at_icol = ekat::subview(dataout, icol);

        mam4::vertical_interpolation::vert_interp(
            team, levsiz, pver, pin_at_icol, pmid_at_icol, datain_at_icol,
            dataout_at_icol);
      });
}

inline void perform_vertical_interpolation(const const_view_1d &altitude_int,
                                           const const_view_2d &zi,
                                           const TracerData &input,
                                           const view_2d output[]) {
  EKAT_REQUIRE_MSG(
      input.file_type == VERT_EMISSION,
      "Error! vertical interpolation only with altitude variable. \n");
  const int ncols                   = input.ncol_;
  const int num_vars                = input.nvars_;
  const int num_vertical_lev_target = output[0].extent(1);
  const int num_vert_packs          = num_vertical_lev_target;
  const int outer_iters             = ncols * num_vars;
  const auto policy_interp =
      ESU::get_default_team_policy(outer_iters, num_vert_packs);
  // FIXME: Get m2km from emaxx.
  const Real m2km    = 1e-3;
  const auto &src_x  = altitude_int;
  const int nsrc     = input.nlev_;
  constexpr int pver = mam4::nlev;
  const int pverp    = pver + 1;
  const auto& data    = input.data;

  // make a local copy of output
  view_2d output_local[MAX_NVARS_TRACER];
  EKAT_REQUIRE_MSG(
        num_vars <= int(MAX_NVARS_TRACER),
        "Error! Number of variables is bigger than NVARS_MAXTRACER. \n");
  for (int ivar = 0; ivar < num_vars; ++ivar)
  {
    // At this stage, begin/end must have the same horiz dimensions
    EKAT_REQUIRE(input.ncol_ == output[ivar].extent_int(0));
    output_local[ivar] = output[ivar];
  }
  Kokkos::parallel_for(
      "tracer_vert_interp_loop", policy_interp,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank() / num_vars;
        const int ivar = team.league_rank() % num_vars;

        const auto src = ekat::subview(data[TracerDataIndex::OUT][ivar], icol);
        const auto trg = ekat::subview(output_local[ivar], icol);
        // FIXME: Try to avoid copy of trg_x by modifying rebin
        // trg_x
        Real trg_x[pver + 1];
        // I am trying to do this:
        // model_z(1:pverp) = m2km * state(c)%zi(i,pverp:1:-1)
        for(int i = 0; i < pverp; ++i) {
          trg_x[pverp - i - 1] = m2km * zi(icol, i);
        }
        team.team_barrier();
        mam4::vertical_interpolation::rebin(team, nsrc, num_vertical_lev_target,
                                            src_x, trg_x, src, trg);
      });
}

inline void advance_tracer_data(
    const std::shared_ptr<AtmosphereInput> &scorpio_reader,   // in
    AbstractRemapper &tracer_horiz_interp,                    // out
    const util::TimeStamp &ts,                                // in
    TracerTimeState &time_state, TracerData &data_tracer,     // out
    const const_view_2d &p_tgt, const const_view_2d &zi_tgt,  // in
    const view_2d output[]) {                                 // out
  /* Update the TracerTimeState to reflect the current time, note the addition
   * of dt */
  time_state.t_now = ts.frac_of_year_in_days();
  /* Update time state and if the month has changed, update the data.*/
  update_tracer_timestate(scorpio_reader, ts, tracer_horiz_interp, time_state,
                          data_tracer);
  // Step 1. Perform time interpolation
  perform_time_interpolation(time_state, data_tracer);

  if(data_tracer.file_type == FORMULA_PS) {
    // Step 2. Compute source pressure levels
    const auto ps = data_tracer.ps[TracerDataIndex::OUT];
    compute_source_pressure_levels(ps, data_tracer.p_src_, data_tracer.hyam,
                                   data_tracer.hybm);
  }

  // Step 3. Perform vertical interpolation
  if(data_tracer.file_type == FORMULA_PS || data_tracer.file_type == ZONAL) {
    perform_vertical_interpolation(data_tracer.p_src_, p_tgt, data_tracer,
                                   output);
  } else if(data_tracer.file_type == VERT_EMISSION) {
    perform_vertical_interpolation(data_tracer.altitude_int_, zi_tgt,
                                   data_tracer, output);
  }

}  // advance_tracer_data

}  // namespace scream::mam_coupling
#endif  // EAMXX_MAM_HELPER_MICRO
