#ifndef EAMXX_MAM_HELPER_MICRO
#define EAMXX_MAM_HELPER_MICRO

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/remap/coarsening_remapper.hpp"
#include "share/grid/remap/refining_remapper_p2p.hpp"
#include "share/grid/remap/identity_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include <ekat/util/ekat_lin_interp.hpp>
#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream::mam_coupling {

  using namespace ShortFieldTagsNames;

  // using npack equal to 1.
  using LIV = ekat::LinInterp<Real,1>;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;

  using view_1d_host    = typename KT::view_1d<Real>::HostMirror;

  constexpr int NVARS_LINOZ=8;
  constexpr int MAX_NVARS_TRACER=5;
  const std::vector<std::string>
  linoz_var_names={"o3_clim", "o3col_clim",
                   "t_clim", "PmL_clim", "dPmL_dO3",
                   "dPmL_dT", "dPmL_dO3col","cariolle_pscs"};

  struct LinozReaderParams {
  int nlevs{-1};
  int nlat{-1};
  // latitude array in linoz data.
  view_1d latitudes;
  // hybrid level pressure at interfaces (1000*(A+B))
  view_1d levs;

  // non_interpolated data from linoz files.
  view_2d data_orig[NVARS_LINOZ];

  // data arrays after horizontal interpolation.
  view_2d data_horiz[NVARS_LINOZ];

  // work arrays
  view_int_1d kupper;
  //
  view_2d pin;
  view_1d col_latitudes;
  };

  // Linoz structures to help manage all of the variables:
  struct LinozTimeState {
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
  }; // LinozTimeState

  struct TracerData{
    TracerData() = default;
    TracerData(const int ncol, const int nlev, const int nvars)
    {
      init (ncol,nlev, nvars);
    }
    void init (const int ncol,
               const int nlev,
               const int nvars){
     ncol_=ncol;
     nlev_=nlev;
     nvars_=nvars;
     EKAT_REQUIRE_MSG (nvars_ <= int(MAX_NVARS_TRACER),
      "Error! Number of variables is bigger than NVARS_MAXTRACER. \n");
    }

    int ncol_{-1};
    int nlev_{-1};
    int nvars_{-1};
    view_2d data[MAX_NVARS_TRACER];
    view_1d ps;
    const_view_1d hyam;
    const_view_1d hybm;

    void allocate_data_views()
    {
      EKAT_REQUIRE_MSG (ncol_ != int(-1),
      "Error! ncols has not been set. \n");
      EKAT_REQUIRE_MSG (nlev_ !=int(-1),
      "Error! nlevs has not been set. \n");

      for (int ivar = 0; ivar< nvars_; ++ivar) {
        data[ivar] = view_2d("linoz_1",ncol_,nlev_);
      }
    } //allocate_data_views

    void allocate_ps()
    {
      EKAT_REQUIRE_MSG (ncol_ != int(-1),
      "Error! ncols has not been set. \n");
      ps = view_1d("ps",ncol_);
    }

    void set_data_views(view_2d list_of_views[])
    {
     for (int ivar = 0; ivar< nvars_; ++ivar) {
      EKAT_REQUIRE_MSG(list_of_views[ivar].data() != 0,
                   "Error! Insufficient memory  size.\n");
      data[ivar] =list_of_views[ivar];
      }
    }

    void set_data_ps(const view_1d& ps_in)
    {
      ps = ps_in;
    }

    void set_hyam_n_hybm(const std::shared_ptr<AbstractRemapper>& horiz_remapper,
                              const std::string& tracer_file_name)
    {

      // Read in hyam/hybm in start/end data, and pad them
      auto nondim = ekat::units::Units::nondimensional();
      const auto io_grid = horiz_remapper->get_src_grid();
      Field hyam_f(FieldIdentifier("hyam",io_grid->get_vertical_layout(true),nondim,io_grid->name()));
      Field hybm_f(FieldIdentifier("hybm",io_grid->get_vertical_layout(true),nondim,io_grid->name()));
      hyam_f.allocate_view();
      hybm_f.allocate_view();
      AtmosphereInput hvcoord_reader(tracer_file_name,io_grid,{hyam_f,hybm_f},true);
      hvcoord_reader.read_variables();
      hvcoord_reader.finalize();
      hyam = hyam_f.get_view<const Real*>();
      hybm = hyam_f.get_view<const Real*>();
    }
  };

  struct LinozData {

    LinozData() = default;
    LinozData(const int ncol, const int nlev)
    {
      init (ncol,nlev);
    }
    void init (const int ncol, const int nlev){
     ncol_=ncol;
     nlev_=nlev;
    }
    int ncol_{-1};
    int nlev_{-1};
    int nvars_{NVARS_LINOZ};
    view_2d data[NVARS_LINOZ];

    void allocate_data_views()
    {
      EKAT_REQUIRE_MSG (ncol_ != int(-1),
      "Error! ncols has not been set. \n");
      EKAT_REQUIRE_MSG (nlev_ !=int(-1),
      "Error! nlevs has not been set. \n");

      for (int ivar = 0; ivar< nvars_; ++ivar) {
        data[ivar] = view_2d("data_tracer",ncol_,nlev_);
      }

    } //allocate_data_views

    void set_data_views(std::vector<view_2d>& list_of_views)
    {
     for (int ivar = 0; ivar< nvars_; ++ivar) {
      EKAT_REQUIRE_MSG(list_of_views[ivar].data() != 0,
                   "Error! Insufficient memory  size.\n");
      data[ivar] =list_of_views[ivar];
      }
    }

    void set_data_views(const view_2d& linoz_o3_clim,
                        const view_2d& linoz_o3col_clim,
                        const view_2d& linoz_t_clim,
                        const view_2d& linoz_PmL_clim,
                        const view_2d& linoz_dPmL_dO3,
                        const view_2d& linoz_dPmL_dT,
                        const view_2d& linoz_dPmL_dO3col,
                        const view_2d& linoz_cariolle_pscs)
    {
      data[0] = linoz_o3_clim;
      data[1] = linoz_o3col_clim;
      data[2] = linoz_t_clim;
      data[3] = linoz_PmL_clim;
      data[4] = linoz_dPmL_dO3;
      data[5] = linoz_dPmL_dT;
      data[6] = linoz_dPmL_dO3col;
      data[7] = linoz_cariolle_pscs;
    }

  };

  // define the different field layouts that will be used for this process
  inline std::shared_ptr<AtmosphereInput>
  create_linoz_data_reader (
      const std::string& linoz_data_file,
      LinozReaderParams& linoz_params,
      const int ncol,
      const const_view_1d const_col_latitudes,
      const ekat::Comm& comm)
  {

     auto make_layout = [](const std::vector<int>& extents,
                        const std::vector<std::string>& names)
  {
    std::vector<FieldTag> tags(extents.size(),CMP);
    return FieldLayout(tags,extents,names);
  };

  // query the file for its resolution
  scorpio::register_file(linoz_data_file,scorpio::Read);
  const int nlevs_data = scorpio::get_dimlen(linoz_data_file,"lev");
  const int nlat_data = scorpio::get_dimlen(linoz_data_file,"lat");
  scorpio::release_file(linoz_data_file);
  linoz_params.nlevs=nlevs_data;
  linoz_params.nlat=nlat_data;

  // create an IO grid, with that number of cols
  // linoz files do not have number of cols,
  // I will use nlat_data instead.

  std::string name="linoz_grid";
  const auto io_grid = std::make_shared<PointGrid>(name,nlat_data,nlevs_data,comm);

  const auto nondim = ekat::units::Units::nondimensional();

  auto scalar2d_layout_linoz = make_layout({nlevs_data, nlat_data},
                                             {"lev","lat"});
  auto scalar1d_lat_layout_linoz = make_layout({nlat_data},
                                             {"lat"});
  auto scalar1d_lev_layout_linoz = make_layout({nlevs_data},
                                             {"lev"});

  const int nvars=NVARS_LINOZ;

  Field lat (FieldIdentifier("lat",  scalar1d_lat_layout_linoz,  nondim,io_grid->name()));
  Field lev (FieldIdentifier("lev",  scalar1d_lev_layout_linoz,  nondim,io_grid->name()));
  lat.allocate_view();
  lev.allocate_view();

  std::vector<Field> io_fields;
  io_fields.push_back(lat);
  io_fields.push_back(lev);

  // FIXME: units are wrong.
  for (int ivar = 0; ivar < nvars; ++ivar) {
    auto var_name = linoz_var_names[ivar];
    // set and allocate fields
    Field f(FieldIdentifier(var_name,  scalar2d_layout_linoz,  nondim,io_grid->name()));
    f.allocate_view();
    io_fields.push_back(f);
    // get views
    linoz_params.data_orig[ivar] = io_fields[ivar+2].get_view<Real**>();
    // allocate views to store data after horizontal interpolation.
    linoz_params.data_horiz[ivar]=view_2d(var_name+"_test", nlevs_data, ncol);
  }

  const auto& latitudes_degree = io_fields[0].get_view<Real*>();
  linoz_params.levs = io_fields[1].get_view< Real*>();

  // make a copy of col_latitudes without const Real
  view_1d col_latitudes("col",ncol);
  Kokkos::deep_copy(col_latitudes, const_col_latitudes);
  linoz_params.col_latitudes = col_latitudes;

  // allocate temp views
  linoz_params.kupper = view_int_1d("kupper",ncol);
  linoz_params.pin = view_2d("pin", ncol,nlevs_data);

  const auto& pin = linoz_params.pin;
  view_1d latitudes("lat",nlevs_data );

  const auto& levs = linoz_params.levs;
  const auto policy_interp = ESU::get_default_team_policy(ncol, nlevs_data);
  const int pi =haero::Constants::pi;
  Kokkos::parallel_for("unit_convertion", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs_data), [&] (const Int& kk) {
      // mbar->pascals
      pin(icol, kk) = levs(kk)*100;
      // radians to degrees
      if (icol==0){
        latitudes(kk) = latitudes_degree (kk)* pi/180.;
      }
      });
    });
  Kokkos::fence();
  linoz_params.latitudes=latitudes;

  return std::make_shared<AtmosphereInput>(linoz_data_file, io_grid,io_fields,true);
  }


 inline void perform_horizontal_interpolation( const              LinozReaderParams& linoz_params, LinozData& linoz_data_out)
 {
  // FIXME: get this inputs from eamxx interface.
  const auto col_latitudes = linoz_params.col_latitudes;
  const int ncol = linoz_data_out.ncol_;
  const int linoz_data_nlev = linoz_data_out.nlev_;
  const int nvars = linoz_data_out.nvars_;

  // We can ||ize over columns as well as over variables and bands
  LIV horiz_interp(linoz_data_nlev, linoz_params.nlat, ncol);
  const auto policy_setup = ESU::get_default_team_policy(1, ncol);
  auto lat = linoz_params.latitudes;

  Kokkos::parallel_for("vert_interp_setup_loop", policy_setup,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    // Setup
    horiz_interp.setup(team, lat, col_latitudes);
  });
  Kokkos::fence();

  // Now use the interpolation object in || over all variables.
  const int outer_iters = linoz_params.nlevs;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, ncol);
  Kokkos::parallel_for("linoz_horizongal_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int kk = team.league_rank();
    const auto x1 = lat;
    const auto x2 = col_latitudes;
    for (int ivar = 0; ivar < nvars; ++ivar) {
        const auto var_org = linoz_params.data_orig[ivar];
        const auto y1 = ekat::subview(var_org,kk);
        const auto y2 = ekat::subview(linoz_params.data_horiz[ivar],kk);
       horiz_interp.lin_interp(team, x1, x2, y1, y2, kk);
    }
  });
  Kokkos::fence();


  // FIXME: Does ekat or kokkos have a call to transpose a view?
  const auto policy_transpose = ESU::get_default_team_policy(ncol*linoz_data_nlev,1);
  Kokkos::parallel_for("transpose_horization_data", policy_transpose,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int icol = team.league_rank() / linoz_data_nlev;
    const int ilev = team.league_rank() % linoz_data_nlev;
    for (int ivar = 0; ivar < nvars; ++ivar) {
      const auto input = linoz_params.data_horiz[ivar];
      const auto output = linoz_data_out.data[ivar];
      output(icol, ilev) = input(ilev, icol);
    }// ivar
  });
  Kokkos::fence();

 }
// Direct port of components/eam/src/chemistry/utils/tracer_data.F90/vert_interp
KOKKOS_INLINE_FUNCTION
void vert_interp(int ncol,
                 int levsiz,
                 int pver,
                 const view_2d&  pin,
                 const const_view_2d&  pmid,
                 const view_2d&  datain,
                 const view_2d&  dataout,
                 //work array
                 const view_int_1d& kupper
                 )
{
    const int one = 1;
    // Initialize index array
    for (int i = 0; i < ncol; ++i) {
      kupper(i)= one;
    } // ncol

    for (int k = 0; k < pver; ++k) {
        // Top level we need to start looking is the top level for the previous k for all column points
        int kkstart = levsiz;
        for (int i = 0; i < ncol; ++i) {
            kkstart = haero::min(kkstart, kupper(i));
        }

        // Store level indices for interpolation
        for (int kk = kkstart - 1; kk < levsiz - 1; ++kk) {
            for (int i = 0; i < ncol; ++i) {
                if (pin(i, kk) < pmid(i, k) && pmid(i, k) <= pin(i, kk + 1)) {
                    kupper(i) = kk;
                }// end if
            } // end for
        } // end kk
        // Interpolate or extrapolate...
        for (int i = 0; i < ncol; ++i) {
            if (pmid(i, k) < pin(i, 0)) {
                dataout(i, k) = datain(i, 0) * pmid(i, k) / pin(i, 0);
            } else if (pmid(i, k) > pin(i, levsiz - 1)) {
                dataout(i, k) = datain(i, levsiz - 1);
            } else {
                Real dpu = pmid(i, k) - pin(i, kupper(i));
                Real dpl = pin(i, kupper(i) + 1) - pmid(i, k);
                dataout(i, k) = (datain(i, kupper(i)) * dpl + datain(i, kupper(i) + 1) * dpu) / (dpl + dpu);
            }// end if
        } // end col
    } // end k

} // vert_interp

inline void perform_vertical_interpolation(const LinozReaderParams& linoz_params,
                                           const const_view_2d& p_mid,
                                           LinozData& non_interpolated_linoz,
                                           LinozData& interpolated_linoz)
{
  const int ncol = interpolated_linoz.ncol_;
  const int nlev = interpolated_linoz.nlev_;

  const int nvars = non_interpolated_linoz.nvars_;
  const int nlevs_linoz = non_interpolated_linoz.nlev_;

  const auto kupper = linoz_params.kupper;
  const auto levs = linoz_params.levs;
  const auto pin = linoz_params.pin;


  const auto policy_interp = ESU::get_default_team_policy(nvars, 1);
  Kokkos::parallel_for("vertical_interpolation_linoz", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
  const int ivar = team.league_rank();
  const auto var_non_inter = non_interpolated_linoz.data[ivar];
  const auto var_inter = interpolated_linoz.data[ivar];
  vert_interp(ncol,
              nlevs_linoz,
              nlev,
              pin,
              p_mid,
              var_non_inter,
              var_inter,
              //work array
              kupper);
    });
    Kokkos::fence();

}//perform_vertical_interpolation

// This function is based on update_spa_timestate
inline void update_linoz_timestate(const util::TimeStamp& ts,
                                   LinozTimeState& time_state,
                                   std::shared_ptr<AtmosphereInput> linoz_reader,
                                   const LinozReaderParams& linoz_params,
                                   LinozData&  data_beg,
                                   LinozData&  data_end)
{
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that SPA assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1; // Make it 0-based
  if (month != time_state.current_month) {
    // Update the SPA time state information
    time_state.current_month = month;
    time_state.t_beg_month = util::TimeStamp({ts.get_year(),month+1,1}, {0,0,0}).frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(),month+1);

    // // Copy spa_end'data into spa_beg'data, and read in the new spa_end
    for (int ivar = 0; ivar < data_beg.nvars_; ++ivar)
    {
      Kokkos::deep_copy(data_beg.data[ivar],data_end.data[ivar]);
    }

    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong values
    //       to be assigned.  A timestep greater than a month is very unlikely so we
    //       will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    linoz_reader->read_variables(next_month);
    perform_horizontal_interpolation(linoz_params, data_end);
    //
  } // end if
} // update_linoz_timestate

KOKKOS_INLINE_FUNCTION
Real linear_interp(const Real& x0, const Real& x1, const Real& t)
{
  return (1 - t)*x0 + t*x1;
} // linear_interp

KOKKOS_INLINE_FUNCTION
view_1d get_var_column (const LinozData& data,
                               const int icol,
                               const int ivar)
{
    return ekat::subview(data.data[ivar],icol);
} // get_var_column
// This function is based on the SPA::perform_time_interpolation function.
 inline void perform_time_interpolation(
  const LinozTimeState& time_state,
  const LinozData&  data_beg,
  const LinozData&  data_end,
  const LinozData&  data_out)
{
  // NOTE: we *assume* data_beg and data_end have the *same* hybrid v coords.
  //       IF this ever ceases to be the case, you can interp those too.
  // Gather time stamp info
  auto& t_now = time_state.t_now;
  auto& t_beg = time_state.t_beg_month;
  auto& delta_t = time_state.days_this_month;

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = data_beg.nvars_;
  const int outer_iters = data_beg.ncol_*num_vars;
  const int num_vert = data_beg.nlev_;
  const auto policy = ESU::get_default_team_policy(outer_iters, num_vert);

  auto delta_t_fraction = (t_now-t_beg) / delta_t;

  EKAT_REQUIRE_MSG (delta_t_fraction>=0 && delta_t_fraction<=1,
      "Error! Convex interpolation with coefficient out of [0,1].\n"
      "  t_now  : " + std::to_string(t_now) + "\n"
      "  t_beg  : " + std::to_string(t_beg) + "\n"
      "  delta_t: " + std::to_string(delta_t) + "\n");

  Kokkos::parallel_for("linoz_time_interp_loop", policy,
    KOKKOS_LAMBDA(const Team& team) {

    // The policy is over ncols*num_vars, so retrieve icol/ivar
    const int icol = team.league_rank() / num_vars;
    const int ivar = team.league_rank() % num_vars;

    // Get column of beg/end/out variable
    auto var_beg = get_var_column (data_beg,icol,ivar);
    auto var_end = get_var_column (data_end,icol,ivar);
    auto var_out = get_var_column (data_out,icol,ivar);

    Kokkos::parallel_for (Kokkos::TeamVectorRange(team,num_vert),
                          [&] (const int& k) {
      var_out(k) = linear_interp(var_beg(k),var_end(k),delta_t_fraction);
    });
  });
  Kokkos::fence();
} // perform_time_interpolation

  // time[3]={year,month, day}
  inline util::TimeStamp convert_date(const int date)
  {
   constexpr int ten_thousand = 10000;
   constexpr int one_hundred = 100;

   int year  = date / ten_thousand;
   int month = (date-year*ten_thousand)/one_hundred;
   int day   = date-year*ten_thousand-month*one_hundred;
   return util::TimeStamp(year,month,day,0,0,0);
  }
  // FIXME: check if this function is implemented in eamxx
  // Assumes 365 days/year, 30 days/month
  inline int compute_days(const util::TimeStamp& ts )
  {
    return ts.get_year()*365 + ts.get_month()*30 + ts.get_day();
  }

  inline void create_linoz_chlorine_reader(
      const std::string& linoz_chlorine_file,
      const util::TimeStamp& model_time,
      const int chlorine_loading_ymd, // in format YYYYMMDD
      std::vector<Real>& values,
      std::vector<int>& time_secs
      )
  {
  auto time_stamp_beg = convert_date(chlorine_loading_ymd);

  const int offset_time = compute_days(time_stamp_beg) - compute_days(model_time);
  scorpio::register_file(linoz_chlorine_file,scorpio::Read);
  const int nlevs_time = scorpio::get_time_len(linoz_chlorine_file);
  for (int itime = 0; itime < nlevs_time; ++itime)
  {
    int date;
    scorpio::read_var(linoz_chlorine_file,"date",&date,itime);
    if (date >= chlorine_loading_ymd ) {
      Real value;
      scorpio::read_var(linoz_chlorine_file,"chlorine_loading",&value, itime);
      values.push_back(value);
      auto time_stamp = convert_date(date);
      time_secs.push_back(compute_days(time_stamp)-offset_time);
    }
  }// end itime
  scorpio::release_file(linoz_chlorine_file);
  }

  inline Real chlorine_loading_advance(const util::TimeStamp& ts,
                                       std::vector<Real>& values,
                                       std::vector<int>& time_secs)
  {

  const int current_time = compute_days(ts);
  int index=0;
  // update index
  for(int i=0; i < int(values.size()); i++){
  if (current_time > time_secs[i] ) {
    index =i;
    break;
    }
  }//

  const Real delt = ( current_time - time_secs[index] ) / ( time_secs[index+1] - time_secs[index] );
  return values[index] + delt*( values[index+1] - values[index] );
  }

inline
std::shared_ptr<AbstractRemapper>
create_horiz_remapper (
    const std::shared_ptr<const AbstractGrid>& model_grid,
    const std::string& trace_data_file,
    const std::string& map_file,
    std::vector<std::string>& var_names
    )
{
  using namespace ShortFieldTagsNames;

  scorpio::register_file(trace_data_file,scorpio::Read);
  const int nlevs_data = scorpio::get_dimlen(trace_data_file,"lev");
  const int ncols_data = scorpio::get_dimlen(trace_data_file,"ncol");
  scorpio::release_file(trace_data_file);

  // We could use model_grid directly if using same num levels,
  // but since shallow clones are cheap, we may as well do it (less lines of code)
  auto horiz_interp_tgt_grid = model_grid->clone("tracer_horiz_interp_tgt_grid",true);
  horiz_interp_tgt_grid->reset_num_vertical_lev(nlevs_data);

  const int ncols_model = model_grid->get_num_global_dofs();
  std::shared_ptr<AbstractRemapper> remapper;
  if (ncols_data==ncols_model) {
    remapper = std::make_shared<IdentityRemapper>(horiz_interp_tgt_grid,IdentityRemapper::SrcAliasTgt);
  } else {
    EKAT_REQUIRE_MSG (ncols_data<=ncols_model,
      "Error! We do not allow to coarsen spa data to fit the model. We only allow\n"
      "       spa data to be at the same or coarser resolution as the model.\n");
    // We must have a valid map file
    EKAT_REQUIRE_MSG (map_file!="",
        "ERROR: Spa data is on a different grid than the model one,\n"
        "       but spa_remap_file is missing from SPA parameter list.");

    remapper = std::make_shared<RefiningRemapperP2P>(horiz_interp_tgt_grid,map_file);
  }

  remapper->registration_begins();
  const auto tgt_grid = remapper->get_tgt_grid();
  const auto layout_2d   = tgt_grid->get_2d_scalar_layout();
  const auto layout_3d_mid = tgt_grid->get_3d_scalar_layout(true);
  const auto nondim = ekat::units::Units::nondimensional();


  for(auto  var_name : var_names){
    Field ifield(FieldIdentifier(var_name,  layout_3d_mid,  nondim,tgt_grid->name()));
    ifield.allocate_view();
    remapper->register_field_from_tgt (ifield);
  }

  Field ps (FieldIdentifier("PS",        layout_2d,  nondim,tgt_grid->name()));
  ps.allocate_view();
  remapper->register_field_from_tgt (ps);
  remapper->registration_ends();

  return remapper;

} // create_horiz_remapper

inline
std::shared_ptr<AtmosphereInput>
create_tracer_data_reader
(
    const std::shared_ptr<AbstractRemapper>& horiz_remapper,
    const std::string& tracer_data_file)
{
  std::vector<Field> io_fields;
  for (int i=0; i<horiz_remapper->get_num_fields(); ++i) {
    io_fields.push_back(horiz_remapper->get_src_field(i));
  }
  const auto io_grid = horiz_remapper->get_src_grid();
  return std::make_shared<AtmosphereInput>(tracer_data_file,io_grid,io_fields,true);
} // create_tracer_data_reader

inline
void
update_tracer_data_from_file(
    std::shared_ptr<AtmosphereInput>& scorpio_reader,
    const util::TimeStamp&            ts,
    const int                         time_index, // zero-based
    AbstractRemapper&                 tracer_horiz_interp,
    TracerData& tracer_data)
{
 // 1. read from field
 scorpio_reader->read_variables(time_index);
 // 2. Run the horiz remapper (it is a do-nothing op if spa data is on same grid as model)
 tracer_horiz_interp.remap(/*forward = */ true);
 //
 const int nvars =tracer_data.nvars_;

  //
 for (int i = 0; i < nvars; ++i) {
  tracer_data.data[i] = tracer_horiz_interp.get_tgt_field (i).get_view< Real**>();
 }

  // Recall, the fields are registered in the order: tracers, ps
 // 3. Copy from the tgt field of the remapper into the spa_data
 tracer_data.ps = tracer_horiz_interp.get_tgt_field(nvars).get_view< Real*>();

} // update_tracer_data_from_file
inline void
update_tracer_timestate(
    std::shared_ptr<AtmosphereInput>& scorpio_reader,
    const util::TimeStamp&            ts,
    AbstractRemapper&                 tracer_horiz_interp,
    LinozTimeState&                   time_state,
    TracerData&  data_tracer_beg,
    TracerData&  data_tracer_end)
{
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that SPA assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1; // Make it 0-based
  if (month != time_state.current_month) {
    //
    const auto tracer_beg = data_tracer_beg.data;
    const auto tracer_end = data_tracer_end.data;
    const int nvars=data_tracer_end.nvars_;

    // Update the SPA time state information
    time_state.current_month = month;
    time_state.t_beg_month = util::TimeStamp({ts.get_year(),month+1,1}, {0,0,0}).frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(),month+1);

    // Copy spa_end'data into spa_beg'data, and read in the new spa_end
    for (int ivar = 0; ivar < nvars ; ++ivar)
    {
      Kokkos::deep_copy(tracer_beg[ivar],  tracer_end[ivar]);
    }
    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong values
    //       to be assigned.  A timestep greater than a month is very unlikely so we
    //       will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    update_tracer_data_from_file(scorpio_reader, ts,
                                 next_month,
                                 tracer_horiz_interp,
                                 data_tracer_end);
  }

} // END updata_spa_timestate

// This function is based on the SPA::perform_time_interpolation function.
 inline void perform_time_interpolation(
  const LinozTimeState& time_state,
  const TracerData& data_tracer_beg,
  const TracerData& data_tracer_end,
  const TracerData& data_tracer_out)
{
  // NOTE: we *assume* data_beg and data_end have the *same* hybrid v coords.
  //       IF this ever ceases to be the case, you can interp those too.
  // Gather time stamp info
  auto& t_now = time_state.t_now;
  auto& t_beg = time_state.t_beg_month;
  auto& delta_t = time_state.days_this_month;

  // We can ||ize over columns as well as over variables and bands
  const auto data_beg = data_tracer_beg.data;
  const auto data_end = data_tracer_end.data;
  const auto data_out = data_tracer_out.data;


  const auto ps_beg = data_tracer_beg.ps;
  const auto ps_end = data_tracer_end.ps;
  const auto ps_out = data_tracer_out.ps;

  const int num_vars = data_tracer_end.nvars_;

  const int ncol = data_tracer_beg.ncol_;
  const int num_vert = data_tracer_beg.nlev_;

  const int outer_iters = ncol*num_vars;

  const auto policy = ESU::get_default_team_policy(outer_iters, num_vert);

  auto delta_t_fraction = (t_now-t_beg) / delta_t;

  EKAT_REQUIRE_MSG (delta_t_fraction>=0 && delta_t_fraction<=1,
      "Error! Convex interpolation with coefficient out of [0,1].\n"
      "  t_now  : " + std::to_string(t_now) + "\n"
      "  t_beg  : " + std::to_string(t_beg) + "\n"
      "  delta_t: " + std::to_string(delta_t) + "\n");

  Kokkos::parallel_for("linoz_time_interp_loop", policy,
    KOKKOS_LAMBDA(const Team& team) {

    // The policy is over ncols*num_vars, so retrieve icol/ivar
    const int icol = team.league_rank() / num_vars;
    const int ivar = team.league_rank() % num_vars;

    // Get column of beg/end/out variable
      auto var_beg = ekat::subview(data_beg[ivar],icol);
      auto var_end = ekat::subview(data_end[ivar],icol);
      auto var_out = ekat::subview(data_out[ivar],icol);

    Kokkos::parallel_for (Kokkos::TeamVectorRange(team,num_vert),
                          [&] (const int& k) {
      var_out(k) = linear_interp(var_beg(k),var_end(k),delta_t_fraction);
    });

    if(ivar==1){
       ps_out(icol) = linear_interp(ps_beg(icol), ps_end(icol),delta_t_fraction);
    }

  });
  Kokkos::fence();
} // perform_time_interpolation

inline void
compute_source_pressure_levels(
  const view_1d& ps_src,
  const view_2d& p_src,
  const const_view_1d& hyam,
  const const_view_1d& hybm)
{
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using C = scream::physics::Constants<Real>;

  constexpr auto P0 = C::P0;
  const int ncols = ps_src.extent(0);
  const int num_vert_packs = p_src.extent(1);
  const auto policy = ESU::get_default_team_policy(ncols, num_vert_packs);

  Kokkos::parallel_for("tracer_compute_p_src_loop", policy,
    KOKKOS_LAMBDA (const Team& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_vert_packs),
                         [&](const int k) {
      p_src(icol,k) = ps_src(icol) * hybm(k)  + P0 * hyam(k);
    });
  });
} // compute_source_pressure_levels


inline void
perform_vertical_interpolation(
  const view_2d& p_src_c,
  const const_view_2d& p_tgt_c,
  const TracerData& input,
  const view_2d output[])
{
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;
  using LIV = ekat::LinInterp<Real,1>;

  // At this stage, begin/end must have the same horiz dimensions
  EKAT_REQUIRE(input.ncol_==output[0].extent(0));
#if 1
  // FIXME: I was encountering a compilation error when using const_view_2d.
  // The issue is fixed by https://github.com/E3SM-Project/EKAT/pull/346.
  // I will keep this code until this PR is merged into the EKAT master branch and
  // we update the EKAT version in our code.
  // I am converting const_view_2d to view_2d.
  auto p_src_ptr = (Real *)p_src_c.data();
  view_2d p_src(p_src_ptr,p_src_c.extent(0),p_src_c.extent(1));
  auto p_tgt_ptr = (Real *)p_tgt_c.data();
  view_2d p_tgt(p_tgt_ptr,p_tgt_c.extent(0),p_tgt_c.extent(1));
#else
  const auto p_src = p_src_c;
  const auto p_tgt = p_tgt_c;
#endif

  const int ncols     = input.ncol_;
  // FIXME: I am getting FPEs if I do not subtract 1 from nlevs_src.
  const int nlevs_src = input.nlev_-1;
  const int nlevs_tgt = output[0].extent(1);

  LIV vert_interp(ncols,nlevs_src,nlevs_tgt);

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = input.nvars_;
  const int num_vert_packs = nlevs_tgt;
  const auto policy_setup = ESU::get_default_team_policy(ncols, num_vert_packs);

  // Setup the linear interpolation object
  Kokkos::parallel_for("tracer_vert_interp_setup_loop", policy_setup,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank();
    // Setup
    vert_interp.setup(team, ekat::subview(p_src,icol),
                            ekat::subview(p_tgt,icol));
  });
  Kokkos::fence();

  // Now use the interpolation object in || over all variables.
  const int outer_iters = ncols*num_vars;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, num_vert_packs);
  Kokkos::parallel_for("tracer_vert_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank() / num_vars;
    const int ivar = team.league_rank() % num_vars;

    const auto x1 = ekat::subview(p_src,icol);
    const auto x2 = ekat::subview(p_tgt,icol);

    const auto y1 = ekat::subview(input.data[ivar],icol);
    const auto y2 = ekat::subview(output[ivar],icol);

    vert_interp.lin_interp(team, x1, x2, y1, y2, icol);
  });
  Kokkos::fence();
}

inline void
advance_tracer_data(std::shared_ptr<AtmosphereInput>& scorpio_reader,
                    AbstractRemapper& tracer_horiz_interp,
                    const util::TimeStamp& ts,
                    LinozTimeState& time_state,
                    TracerData& data_tracer_beg,
                    TracerData& data_tracer_end,
                    TracerData& data_tracer_out,
                    const view_2d& p_src,
                    const const_view_2d& p_tgt,
                    const view_2d output[])
{

  /* Update the TracerTimeState to reflect the current time, note the addition of dt */
  time_state.t_now = ts.frac_of_year_in_days();
  /* Update time state and if the month has changed, update the data.*/
  update_tracer_timestate(
    scorpio_reader,
    ts,
    tracer_horiz_interp,
    time_state,
    data_tracer_beg,
    data_tracer_end);
  // Step 1. Perform time interpolation
  perform_time_interpolation(
  time_state,
  data_tracer_beg,
  data_tracer_end,
  data_tracer_out);
  // Step 2. Compute source pressure levels
  compute_source_pressure_levels(
    data_tracer_out.ps,
    p_src,
    data_tracer_out.hyam,
    data_tracer_out.hybm);

  // Step 3. Perform vertical interpolation

  perform_vertical_interpolation(
  p_src,
  p_tgt,
  data_tracer_out,
  output);

}// advance_tracer_data


} // namespace scream::mam_coupling
#endif //EAMXX_MAM_HELPER_MICRO
