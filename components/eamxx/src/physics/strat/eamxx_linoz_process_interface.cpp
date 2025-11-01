#include <physics/strat/eamxx_linoz_process_interface.hpp>
#include <ekat_team_policy_utils.hpp>
#include "share/algorithm/eamxx_data_interpolation.hpp"
#include "physics/mam/readfiles/vertical_remapper_mam4.hpp"
#include <physics/mam/mam_coupling.hpp>
#include "physics/mam/readfiles/tracer_reader_utils.hpp"
#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"

namespace scream
{

STRATLinoz::STRATLinoz(const ekat::Comm &comm, const ekat::ParameterList &params)
 : AtmosphereProcess(comm, params)
{
  // LINOZ namelist parameters
    m_config.o3_lbl = m_params.get<int>("mam4_o3_lbl");
    m_config.o3_tau = m_params.get<double>("mam4_o3_tau");
    m_config.o3_sfc = m_params.get<double>("mam4_o3_sfc");
    m_config.psc_T  = m_params.get<double>("mam4_psc_T");

    m_orbital_year = m_params.get<int>("orbital_year", -9999);

    // Get orbital parameters from yaml file
    m_orbital_eccen = m_params.get<double>("orbital_eccentricity", -9999);
    m_orbital_obliq = m_params.get<double>("orbital_obliquity", -9999);
    m_orbital_mvelp = m_params.get<double>("orbital_mvelp", -9999);

    use_prescribed_ozone_=true;

}
void
STRATLinoz::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  grid_                 = grids_manager->get_grid("physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();      // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // number of levels per column
  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);
  // Creating a Linoz reader and setting Linoz parameters involves reading data
  // from a file and configuring the necessary parameters for the Linoz model.
  m_var_names_linoz = {
        "o3_clim",  "o3col_clim", "t_clim",      "PmL_clim",
        "dPmL_dO3", "dPmL_dT",    "dPmL_dO3col", "cariolle_pscs"};

  constexpr auto nondim = Units::nondimensional();
  // The DataInterpolation class uses Field. We save these fields in FM.
  for(const auto &field_name : m_var_names_linoz) {
      add_field<Computed>(field_name, scalar3d_mid, nondim, grid_name);
  }
  // get column geometry and locations
  col_latitudes_ = grid_->get_geometry_data("lat").get_view<const Real *>();

  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);
  // Total pressure [Pa] at midpoints
  add_field<Required>("p_mid", scalar3d_mid, Pa, grid_name);
  // Layer thickness(pdel) [Pa] at midpoints
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);
  constexpr auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  // specific humidity [kg/kg]
  add_tracer<Required>("qv", grid_, q_unit);
  // If not using prescribed O3 (i.e., prognostic O3), we add it as a tracer
  add_tracer<Updated>("O3", grid_, q_unit);

}

// set DataInterpolation object for linoz reader.
void STRATLinoz::set_linoz_reader(){
  auto pmid = get_field_in("p_mid");
  // Beg of any year, since we use yearly periodic timeline
  util::TimeStamp ref_ts_linoz (1,1,1,0,0,0);
  const auto m_linoz_file_name = m_params.get<std::string>("mam4_linoz_file_name");
  const std::string linoz_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");
  std::vector<Field> linoz_fields;
  for(const auto &field_name : m_var_names_linoz) {
      linoz_fields.push_back(get_field_out(field_name));
  }

  data_interp_linoz_ = std::make_shared<DataInterpolation>(grid_,linoz_fields);
  data_interp_linoz_->setup_time_database ({m_linoz_file_name},util::TimeLine::YearlyPeriodic, ref_ts_linoz);
  data_interp_linoz_->create_horiz_remappers (linoz_map_file=="none" ? "" : linoz_map_file);
  data_interp_linoz_->set_logger(m_atm_logger);

  DataInterpolation::VertRemapData remap_data_linoz;
  remap_data_linoz.vr_type = DataInterpolation::Static1D;
  // lev is the name of variables for vertical interpolation
  remap_data_linoz.pname = "lev";
  remap_data_linoz.pmid = pmid;
  // Static1D can also be employed instead of mam4xx routine.
  bool mam4_use_mam4xx_linoz_vert_remap=true;
  if (mam4_use_mam4xx_linoz_vert_remap){
    // We are using a custom remapper that invokes the MAM4XX routine
    // for vertical interpolation.
    // The type used is VertRemapType::MAM4_ZONAL.
    auto grid_after_hremap_linoz = data_interp_linoz_->get_grid_after_hremap();
    auto vertical_remapper_linoz = std::make_shared<VerticalRemapperMAM4>(grid_after_hremap_linoz, grid_,
    VerticalRemapperMAM4::VertRemapType::MAM4_ZONAL);
    remap_data_linoz.custom_remapper=vertical_remapper_linoz;
  }
  data_interp_linoz_->create_vert_remapper (remap_data_linoz);
  data_interp_linoz_->init_data_interval (start_of_step_ts());
}

void STRATLinoz::initialize_impl(const RunType run_type) {

  // climatology data for linear stratospheric chemistry
  auto ts = start_of_step_ts();
  std::string linoz_chlorine_file =
        m_params.get<std::string>("mam4_linoz_chlorine_file");
  int chlorine_loading_ymd = m_params.get<int>("mam4_chlorine_loading_ymd");
  scream::mam_coupling::create_linoz_chlorine_reader(
        linoz_chlorine_file, ts, chlorine_loading_ymd, chlorine_values_,
        chlorine_time_secs_);

  set_linoz_reader();

  acos_cosine_zenith_host_ = view_1d_host("host_acos(cosine_zenith)", ncol_);
  init_temporary_views();
}

int STRATLinoz::get_len_temporary_views() {
  int work_len              = 0;
  // m_vmr, m_o3_col_dens
  work_len += 2*ncol_ * nlev_+ ncol_ * (nlev_+1) + ncol_;
  return work_len;
}

void STRATLinoz::init_temporary_views() {
  //FIXME: use buffer instead of temporary_views
  view_1d temporary_views("temporary_views",get_len_temporary_views());
  auto work_ptr             = (Real *)temporary_views.data();
  m_o3_col_dens = view_2d(work_ptr, ncol_, nlev_);
  work_ptr += ncol_ * nlev_;
  m_vmr = view_2d(work_ptr, ncol_, nlev_);
  work_ptr += ncol_ * nlev_;
  m_o3_col_deltas = view_2d(work_ptr, ncol_, nlev_+1);
  work_ptr += ncol_ * (nlev_+1);
  acos_cosine_zenith_=view_1d(work_ptr, ncol_);
  work_ptr += ncol_;
  // Error check
  // NOTE: workspace_provided can be larger than workspace_used, but let's try
  // to use the minimum amount of memory
  const int workspace_used     = work_ptr - temporary_views.data();
  const int workspace_provided = temporary_views.extent(0);
  EKAT_REQUIRE_MSG(workspace_used == workspace_provided,
                   "Error: workspace_used (" + std::to_string(workspace_used) +
                       ") and workspace_provided (" +
                       std::to_string(workspace_provided) +
                       ") should be equal. \n");
}

void STRATLinoz::run_impl(const double dt) {
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  using PF = scream::PhysicsFunctions<DefaultDevice>;

  const auto& T_mid  = get_field_in("T_mid").get_view<const Real **>();
  const auto& p_mid  = get_field_in("p_mid").get_view<const Real **>();
  const auto& p_del  = get_field_in("pseudo_density").get_view<const Real **>();

  const int ncol = ncol_;
  const int nlev = nlev_;
  const auto policy = TPF::get_default_team_policy(ncol, nlev);
  const_view_1d &col_latitudes     = col_latitudes_;
  auto& vmr =m_vmr;
  // CHECK: I assume it is mmr wet
  const auto& mmr_o3  = get_field_out("O3").get_view<Real **>();
  const auto& qv = get_field_in("qv").get_view<const Real **>();
  constexpr Real mw_o3=  47.998200; // g/mol

  // Conversion from wet mmr to  dry vmr
  Kokkos::parallel_for(
    "STRATLinoz::run_impl::convert_to_vmr_dry", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
    const int icol     = team.league_rank();   // column index
    Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        const Real mmr_o3_dry =
        PF::calculate_drymmr_from_wetmmr(qv(icol, kk), mmr_o3(icol,kk));
        vmr(icol, kk) =  mam4::conversions::vmr_from_mmr(mmr_o3_dry,mw_o3);
       });
  });

  // climatology data for linear stratospheric chemistry
  // ozone (climatology) [vmr]
  view_2d linoz_o3_clim;
  // column o3 above box (climatology) [Dobson Units (DU)]
  view_2d linoz_o3col_clim;
  // temperature (climatology) [K]
  view_2d linoz_t_clim;
  // P minus L (climatology) [vmr/s]
  view_2d linoz_PmL_clim;
  // sensitivity of P minus L to O3 [1/s]
  view_2d linoz_dPmL_dO3;
  // sensitivity of P minus L to T3 [K]
  view_2d linoz_dPmL_dT;
  // sensitivity of P minus L to overhead O3 column [vmr/DU]
  view_2d linoz_dPmL_dO3col;
  // Cariolle parameter for PSC loss of ozone [1/s]
  view_2d linoz_cariolle_pscs;
  view_2d linoz_views[8];

  data_interp_linoz_->run(end_of_step_ts());
  for (size_t i = 0; i < m_var_names_linoz.size(); ++i) {
      linoz_views[i] = get_field_out(m_var_names_linoz[i]).get_view<Real **>();
  }
  linoz_o3_clim = linoz_views[0];
  linoz_o3col_clim = linoz_views[1];
  linoz_t_clim = linoz_views[2];
  linoz_PmL_clim = linoz_views[3];
  linoz_dPmL_dO3 = linoz_views[4];
  linoz_dPmL_dT = linoz_views[5];
  linoz_dPmL_dO3col = linoz_views[6];
  linoz_cariolle_pscs = linoz_views[7];

  /* Gather time and state information for interpolation */
  const auto ts = end_of_step_ts();
  const Real chlorine_loading = scream::mam_coupling::chlorine_loading_advance(
      ts, chlorine_values_, chlorine_time_secs_);
  m_config.chlorine_loading=chlorine_loading;

  const auto& linoz_conf=m_config;
  const auto& o3_col_dens = m_o3_col_dens;

    // Compute orbital parameters; these are used both for computing
  // the solar zenith angle.
  // Note: We are following the RRTMGP EAMxx interface to compute the zenith
  // angle. This operation is performed on the host because the routine
  // shr_orb_cosz_c2f has not been ported to C++.
  auto orbital_year = m_orbital_year;
  // Note: We need double precision because
  // shr_orb_params_c2f and shr_orb_decl_c2f only support double precision.
  double obliqr, lambm0, mvelpp;
  double eccen = m_orbital_eccen;
  double obliq = m_orbital_obliq;
  double mvelp = m_orbital_mvelp;
  // Use the orbital parameters to calculate the solar declination and
  // eccentricity factor
  double delta, eccf;
  if(eccen >= 0 && obliq >= 0 && mvelp >= 0) {
    // use fixed orbital parameters; to force this, we need to set
    // orbital_year to SHR_ORB_UNDEF_INT, which is exposed through
    // our c2f bridge as shr_orb_undef_int_c2f
    orbital_year = shr_orb_undef_int_c2f;
  } else if(orbital_year < 0) {
    // compute orbital parameters based on current year
    orbital_year = start_of_step_ts().get_year();
  }
  shr_orb_params_c2f(&orbital_year,                                       // in
                     &eccen, &obliq, &mvelp, &obliqr, &lambm0, &mvelpp);  // out

  // Want day + fraction; calday 1 == Jan 1 0Z
  auto calday = start_of_step_ts().frac_of_year_in_days() + 1;
  shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0, obliqr,  // in
                   &delta, &eccf);                         // out
  {
    const auto col_latitudes_host =
        grid_->get_geometry_data("lat").get_view<const Real *, Host>();
    const auto col_longitudes_host =
        grid_->get_geometry_data("lon").get_view<const Real *, Host>();
    // get a host copy of lat/lon
    // Determine the cosine zenith angle
    // NOTE: Since we are bridging to F90 arrays this must be done on HOST and
    // then deep copied to a device view.

    // Now use solar declination to calculate zenith angle for all points
    for(int i = 0; i < ncol; i++) {
      Real lat =
          col_latitudes_host(i) * M_PI / 180.0;  // Convert lat/lon to radians
      Real lon = col_longitudes_host(i) * M_PI / 180.0;
      // what's the aerosol microphys frequency?
      Real temp = shr_orb_cosz_c2f(calday, lat, lon, delta, dt);
      acos_cosine_zenith_host_(i) = acos(temp);
    }
    Kokkos::deep_copy(acos_cosine_zenith_, acos_cosine_zenith_host_);
  }
  const auto zenith_angle = acos_cosine_zenith_;
  // FIXME: we found a bug in the following bock of code in mam4xx. It needs to be updated.
  const auto& o3_col_deltas=m_o3_col_deltas;
  Kokkos::parallel_for(
    "MAMMicrophysics::run_impl::compute_o3_column_density", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
    const int icol     = team.league_rank();   // column index
    // calculate o3 column densities (first component of col_dens in Fortran
    // code)
    auto o3_col_dens_i = ekat::subview(o3_col_dens, icol);
    auto o3_col_deltas_icol = ekat::subview(o3_col_deltas, icol);
    // NOTE: if we need o2 column densities, set_ub_col and setcol must be changed
    const int nlev = mam4::nlev;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](const int kk) {
      // compute the change in o3 density for this column above its neighbor
      constexpr Real xfactor = 2.8704e21 / (9.80616 * 1.38044); // BAD_CONSTANT!
      o3_col_deltas_icol(kk + 1) = xfactor * p_del(icol, kk) * vmr(icol, kk);
    });
    team.team_barrier();
    // sum the o3 column deltas to densities
    mam4::mo_photo::setcol(team, o3_col_deltas.data(), // in
                         o3_col_dens_i);        // out
  });

  Kokkos::parallel_for(
    "STRATLinoz::run_impl::linoz", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol     = team.league_rank();   // column index
      const auto o3_col_dens_i = ekat::subview(o3_col_dens, icol);
      // convert column latitude to radians
      const Real rlats = col_latitudes(icol) * M_PI / 180.0;

      mam4::microphysics::LinozData linoz_data;
      linoz_data.linoz_o3_clim_icol = ekat::subview(linoz_o3_clim, icol);
      linoz_data.linoz_t_clim_icol  = ekat::subview(linoz_t_clim, icol);
      linoz_data.linoz_o3col_clim_icol = ekat::subview(linoz_o3col_clim, icol);
      linoz_data.linoz_PmL_clim_icol = ekat::subview(linoz_PmL_clim, icol);
      linoz_data.linoz_dPmL_dO3_icol = ekat::subview(linoz_dPmL_dO3, icol);
      linoz_data.linoz_dPmL_dT_icol  = ekat::subview(linoz_dPmL_dT, icol);
      linoz_data.linoz_dPmL_dO3col_icol = ekat::subview(linoz_dPmL_dO3col, icol);
      linoz_data.linoz_cariolle_pscs_icol = ekat::subview(linoz_cariolle_pscs, icol);

      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
      //-----------------
      // LINOZ chemistry
      //-----------------
      const Real temp = T_mid(icol, kk);
      const Real pmid = p_mid(icol, kk);
      const Real pdel = p_del(icol, kk);

      // the following things are diagnostics, which we're not
      // including in the first rev
      Real do3_linoz = 0, do3_linoz_psc = 0, ss_o3 = 0, o3col_du_diag = 0,
           o3clim_linoz_diag = 0, zenith_angle_degrees = 0;
      // index of "O3" in solsym array (in EAM)
      mam4::lin_strat_chem::lin_strat_chem_solve_kk(
          // in
          o3_col_dens_i(kk), temp, zenith_angle(icol), pmid, dt, rlats,
          linoz_data.linoz_o3_clim_icol(kk), linoz_data.linoz_t_clim_icol(kk),
          linoz_data.linoz_o3col_clim_icol(kk),
          linoz_data.linoz_PmL_clim_icol(kk),
          linoz_data.linoz_dPmL_dO3_icol(kk), linoz_data.linoz_dPmL_dT_icol(kk),
          linoz_data.linoz_dPmL_dO3col_icol(kk),
          linoz_data.linoz_cariolle_pscs_icol(kk), linoz_conf.chlorine_loading,
          linoz_conf.psc_T,
          // in/out
          vmr(icol, kk),
          // outputs that are not used
          do3_linoz, do3_linoz_psc, ss_o3, o3col_du_diag, o3clim_linoz_diag,
          zenith_angle_degrees);

      // Update source terms above the ozone decay threshold
      if (kk >= nlev - linoz_conf.o3_lbl) {
        const Real o3l_vmr_old = vmr(icol, kk);
        Real do3mass = 0;
        const Real o3l_vmr_new =
            mam4::lin_strat_chem::lin_strat_sfcsink_kk(dt, pdel,          // in
                                                       o3l_vmr_old,       // in
                                                       linoz_conf.o3_sfc, // in
                                                       linoz_conf.o3_tau, // in
                                                       do3mass);          // out
        // Update the mixing ratio (vmr) for O3
        vmr(icol, kk) = o3l_vmr_new;
      }
        });
        });

    // Conversion from mmr to vmr
  Kokkos::parallel_for(
    "STRATLinoz::run_impl::convert_to_mmr_wet", policy,
    KOKKOS_LAMBDA(const ThreadTeam &team) {
    const int icol     = team.league_rank();   // column index
    Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, nlev),
       [&](const int kk) {
        const Real qv_dry =
        PF::calculate_drymmr_from_wetmmr(qv(icol, kk), qv(icol, kk));
        const Real mmr_dry = mam4::conversions::mmr_from_vmr(vmr(icol, kk),mw_o3);
        mmr_o3(icol,kk) = PF::calculate_wetmmr_from_drymmr(
          mmr_dry, qv_dry);

       });
  });

      }
}  // namespace scream
