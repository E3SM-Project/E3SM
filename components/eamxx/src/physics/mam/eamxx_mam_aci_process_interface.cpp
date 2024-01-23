#include <physics/mam/eamxx_mam_aci_process_interface.hpp>
#include "ekat/util/ekat_units.hpp"
#include "mam4xx/aero_config.hpp"
#include "mam4xx/ndrop.hpp"

namespace scream
{

namespace
{
void copy_scream_array_to_mam4xx(
    haero::ThreadTeamPolicy team_policy,
    MAMAci::view_2d mam4xx_view[], 
    MAMAci::const_view_3d scream_view,
    const int nlev,
    const int num_views) 
{
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, nlev), KOKKOS_LAMBDA(int kk) { 
      for (int i=0; i<num_views; ++i) {
        mam4xx_view[i](icol,kk) = scream_view(icol, kk, i);
      }
    });
  });
}
void compute_w0_and_rho(
    haero::ThreadTeamPolicy team_policy,
    MAMAci::view_2d w0, 
    MAMAci::view_2d rho, 
    MAMAci::const_view_2d p_mid,
    MAMAci::const_view_2d T_mid,
    MAMAci::const_view_2d omega,
    const int top_lev,
    const int nlev) 
{
  //Get physical constants
  using C  = physics::Constants<Real>;
  static constexpr auto gravit = C::gravit; // Gravity [m/s2]
  static constexpr auto rair   = C::Rair;   // Gas constant for dry air [J/(kg*K) or J/Kg/K]
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, top_lev), KOKKOS_LAMBDA(int kk) { 
      w0(icol,kk) = 0;
      rho(icol, kk) = -999.0;
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, top_lev, nlev), KOKKOS_LAMBDA(int kk) { 
      rho(icol,kk) = p_mid(icol,kk)/(rair*T_mid(icol,kk));
      w0(icol,kk) = -1.0*omega(icol,kk)/(rho(icol,kk)*gravit);
    });
  });
}
void compute_tke_using_w_sec(
    haero::ThreadTeamPolicy team_policy,
    MAMAci::view_2d tke, 
    MAMAci::const_view_2d w_sec,
    const int nlev) 
{
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    // FIXME Is this the correct boundary condition for tke at the surface?
    // TKE seems to be at interfaces but w_sec is at cell centers so this
    // descrepensy needs to be worked out.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, nlev+1), KOKKOS_LAMBDA(int kk) { 
      tke(icol,kk) = (3.0/2.0)*w_sec(icol,kk);
    });
  });
}
void compute_subgrid_scale_velocities(
    haero::ThreadTeamPolicy team_policy,
    MAMAci::view_2d wsub, 
    MAMAci::view_2d wsubice, 
    MAMAci::view_2d wsig, 
    MAMAci::const_view_2d tke,
    const Real wsubmin,
    const int top_lev,
    const int nlev) 
{
  // More refined computation of sub-grid vertical velocity
  // Set to be zero at the surface by initialization.
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, top_lev), KOKKOS_LAMBDA(int kk) { 
      wsub(icol,kk)  = wsubmin;
      wsubice(icol,kk) = 0.001;
      wsig(icol,kk)  = 0.001;
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, top_lev, nlev), KOKKOS_LAMBDA(int kk) { 
      wsub(icol,kk)  = haero::sqrt(0.5*(tke(icol,kk) + tke(icol,kk+1))*(2.0/3.0));
      wsig(icol,kk)  = mam4::utils::min_max_bound(0.001, 10.0, wsub(icol,kk));
      wsubice(icol,kk) = mam4::utils::min_max_bound(0.2, 10.0, wsub(icol,kk));
      wsub(icol,kk)  = haero::max(wsubmin, wsub(icol,kk));
    });
  });
}


KOKKOS_INLINE_FUNCTION
Real subgrid_mean_updraft(const Real w0, const Real wsig)
{
/* ---------------------------------------------------------------------------------
 Purpose: Calculate the mean updraft velocity inside a GCM grid assuming the
          vertical velocity distribution is Gaussian and peaks at the
          GCM resolved large-scale vertical velocity.
          When icenul_wsub_scheme = 2, the model uses the mean updraft velocity as the
          characteristic updraft velocity to calculate the ice nucleation rate.
 Author:  Kai Zhang (kai.zhang@pnnl.gov)
 Last Modified: Oct, 2015
-------------------------------------------------------------------------------- */

// interface

//   in :: wsig  standard deviation (m/s)
//   in :: w0  large scale vertical velocity (m/s)
//   out::   mean updraft velocity(m/s) -> characteristic w*

  // FIXME should nbin be a user parameter?
  const int nbin = 50;

  
  using C  = physics::Constants<Real>;
  constexpr Real pi       = C::Pi;
  Real zz[nbin], wa = 0, ww = 0;
  int kp = 0;

  const Real sigma  = haero::max(0.001, wsig);
  const Real wlarge = w0;

  const Real xx = 6.0 * sigma / nbin;

  for (int ibin = 0; ibin < nbin; ++ibin) {
    Real yy = wlarge - 3.0*sigma + 0.5*xx;
    yy += (ibin-1)*xx;
    // wbar = integrator < w * f(w) * dw >
    zz[ibin] = yy * haero::exp(-1.0*haero::square(yy-wlarge)/(2*sigma*sigma))/(sigma*haero::sqrt(2*pi))*xx;
  }
  for (int ibin = 0; ibin < nbin; ++ibin) {
    if (zz[ibin] > 0) ++kp, wa += zz[ibin];
  }
  if (kp) {
    // wbar = integrator < w * f(w) * dw >
    ww = wa;
  } else {
    ww = 0.001;
  }
  return ww;
}
void compute_subgrid_mean_updraft_velocities(
    haero::ThreadTeamPolicy team_policy,
    MAMAci::view_2d w2, 
    MAMAci::const_view_2d w0,
    MAMAci::const_view_2d wsig,
    const int nlev) 
{
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, nlev), KOKKOS_LAMBDA(int kk) { 
      w2(icol,kk) = subgrid_mean_updraft(w0(icol,kk), wsig(icol,kk));
    });
  });
}
void compute_aitken_dry_diameter(
    haero::ThreadTeamPolicy team_policy,
    MAMAci::view_2d aitken_dry_dia, 
    MAMAci::const_view_3d dgnum,
    const int top_lev) 
{
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    const int aitken_idx = static_cast<int>(mam4::ModeIndex::Aitken);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, top_lev), KOKKOS_LAMBDA(int kk) { 
      aitken_dry_dia(icol,kk)  = dgnum(icol, kk, aitken_idx);
    });
  });
}

void compute_nucleate_ice_tendencies(
    haero::ThreadTeamPolicy team_policy,
    const int nlev) 
{
}

}


//FIXME: The following variables are namelist variables
const Real wsubmin = 1;

MAMAci::MAMAci(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params){
}


//Return type of the process
AtmosphereProcessType MAMAci::type() const{
  return AtmosphereProcessType::Physics;
}

//return name of the process
std::string MAMAci::name() const{
  return "mam4_aci";
  }

//set grid for all the inputs and outputs
void MAMAci::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  m_atm_logger->log(ekat::logger::LogLevel::info,"Calling ACI set grid");

  grid_ = grids_manager->get_grid("Physics"); //Use physics grid
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs(); // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // Number of levels per column

  Kokkos::resize(rho_,     ncol_, nlev_);
  Kokkos::resize(w0_,      ncol_, nlev_);
  Kokkos::resize(tke_,     ncol_, nlev_+1);
  Kokkos::resize(wsub_,    ncol_, nlev_);
  Kokkos::resize(wsubice_, ncol_, nlev_);
  Kokkos::resize(wsig_,    ncol_, nlev_);
  Kokkos::resize(w2_,      ncol_, nlev_);
  Kokkos::resize(lcldn_,   ncol_, nlev_);
  Kokkos::resize(lcldo_,   ncol_, nlev_);
  Kokkos::resize(aitken_dry_dia_, ncol_, nlev_);
  Kokkos::resize(rpdel_,   ncol_, nlev_);

  for (int  i=0; i<15; ++i) {
     Kokkos::resize(dropmixnuc_scratch_mem_[i],ncol_, nlev_);
  }
  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variables
  FieldLayout scalar3d_layout_mid{ {COL, LEV}, {ncol_, nlev_} }; // mid points
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {ncol_, nlev_+1} }; //interfaces

  using namespace ekat::units;
  auto q_unit = kg/kg; // units of mass mixing ratios of tracers
  q_unit.set_string("kg/kg");
  auto n_unit = 1/kg; // units of number mixing ratios of tracers
  n_unit.set_string("#/kg");
  auto nondim = Units::nondimensional();
  // atmospheric quantities
  add_field<Required>("qc",             scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud liquid mass mixing ratio [kg/kg]
  add_field<Required>("qi",             scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud ice mass mixing ratio [kg/kg]
  add_field<Required>("nc",             scalar3d_layout_mid, n_unit, grid_name, "tracers");// cloud liquid number mixing ratio [1/kg]
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,      grid_name); // Temperature[K] at midpoints
  add_field<Required>("omega",          scalar3d_layout_mid, Pa/s,   grid_name); // Vertical pressure velocity [Pa/s] at midpoints
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa,     grid_name); // Total pressure [Pa] at midpoints
  add_field<Required>("p_int",          scalar3d_layout_int, Pa,     grid_name); // Total pressure [Pa] at interfaces
  add_field<Required>("qv",             scalar3d_layout_mid, q_unit, grid_name); // Water vapor mixing ratio [kg vapor / kg dry air]
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,     grid_name); // Layer thickness(pdel) [Pa] at midpoints
										 //
										 //
  add_field<Computed>("stratiform_cloud_fraction", scalar3d_layout_mid, nondim, grid_name); // Layer thickness(pdel) [Pa] at midpoints
  add_field<Computed>("activation_fraction_accum", scalar3d_layout_mid, nondim, grid_name); // Layer thickness(pdel) [Pa] at midpoints
  add_field<Computed>("activation_fraction_coarse", scalar3d_layout_mid, nondim, grid_name); // Layer thickness(pdel) [Pa] at midpoints

  
  // FIXME: These should all come from the intput file.  Might need to be added to the input file?
  add_field<Required>("qi_coarse_dst",   scalar3d_layout_mid, q_unit, grid_name); // Dust mixing ratio for coarse mode [kg/kg]
  add_field<Required>("qi_coarse_nacl",  scalar3d_layout_mid, q_unit, grid_name); // Salt mixing ratio for coarse mode [kg/kg]
  add_field<Required>("qi_coarse_so4",   scalar3d_layout_mid, q_unit, grid_name); // Sulfuric Acid mixing ratio for coarse mode [kg/kg]
  add_field<Required>("qi_coarse_mom",   scalar3d_layout_mid, q_unit, grid_name); // Marine Organic Matter mixing ratio for coarse mode [kg/kg]
  add_field<Required>("qi_coarse_bc",    scalar3d_layout_mid, q_unit, grid_name); // Black Carbon mixing ratio for coarse mode [kg/kg]
  add_field<Required>("qi_coarse_pom",   scalar3d_layout_mid, q_unit, grid_name); // Primary Organic Matter mixing ratio for coarse mode [kg/kg]
  add_field<Required>("qi_coarse_soa",   scalar3d_layout_mid, q_unit, grid_name); // Secondary Organic Aerosol mixing ratio for coarse mode [kg/kg]
  //add_field<Required>("qi_accum_bc",    scalar3d_layout_mid, q_unit, grid_name); // Black Carbon mixing ratio for accum mode [kg/kg]
  //add_field<Required>("qi_accum_pom",   scalar3d_layout_mid, q_unit, grid_name); // Primary Organic Matter mixing ratio for accum mode [kg/kg]
  //add_field<Required>("qi_accum_soa",   scalar3d_layout_mid, q_unit, grid_name); // Secondary Organic Aerosol mixing ratio for accum mode [kg/kg]
  //add_field<Required>("qi_accumn_bc",    scalar3d_layout_mid, q_unit, grid_name); // Black Carbon mixing ratio for primary carbon mode [kg/kg]
  //add_field<Required>("qi_accumn_mom",   scalar3d_layout_mid, q_unit, grid_name); // Marine Organic Matter mixing ratio for primary carbon mode [kg/kg]]
  //add_field<Required>("qi_accumn_pom",   scalar3d_layout_mid, q_unit, grid_name); // Primary Organic Matter mixing ratio for primary carbon mode [kg/kg]kg]
  //add_field<Required>("qi_pcarbon_bc",    scalar3d_layout_mid, q_unit, grid_name); // Black Carbon mixing ratio for primary carbon mode [kg/kg]
  //add_field<Required>("qi_pcarbon_mom",   scalar3d_layout_mid, q_unit, grid_name); // Marine Organic Matter mixing ratio for primary carbon mode [kg/kg]]
  //add_field<Required>("qi_pcarbon_pom",   scalar3d_layout_mid, q_unit, grid_name); // Primary Organic Matter mixing ratio for primary carbon mode [kg/kg]kg]
  //add_field<Required>("ni_accum",        scalar3d_layout_mid, 1/kg,   grid_name); // Coarse mode number mixing ratio [1/kg dry air]
  //add_field<Required>("ni_coarse",       scalar3d_layout_mid, 1/kg,   grid_name); // Coarse mode number mixing ratio [1/kg dry air]
  //add_field<Required>("ni_aitken",       scalar3d_layout_mid, 1/kg,   grid_name); // Aitken mode number mixing ratio [1/kg dry air]
  //add_field<Required>("zm",             scalar3d_layout_mid, m,      grid_name); // geopotential height of level (m)
  //add_field<Required>("state_q",FieldLayout{ {COL, LEV, CMP}, {ncol_, nlev_, mam4::ndrop::nvars} } , q_unit, grid_name); // aerosol mmrs [kg/kg]
  //add_field<Required>("ncldwtr", , n_unit, grid_name); // initial droplet number mixing ratio [#/kg]
  //add_field<Required>("cldo", , unitless, grid_name); // cloud fraction on previous time step [fraction]
  //
  add_field<Required>("w_updraft",      scalar3d_layout_mid, q_unit, grid_name); // updraft velocity [m/s]
  add_field<Required>("qqcw", FieldLayout{ {COL, LEV, CMP}, {ncol_, nlev_, mam4::ndrop::ncnst_tot} } , q_unit, grid_name); // cloud-borne aerosol mass, number mixing ratios [#/kg or kg/kg]
  add_field<Required>("cldfrac_tot",    scalar3d_layout_mid, nondim, grid_name); // cloud fraction [nondimentional]


  add_field<Computed>("nihf",scalar3d_layout_mid , 1/m/m/m,   grid_name); // number conc of ice nuclei due to heterogeneous freezing [1/m3]
  add_field<Computed>("niim",scalar3d_layout_mid , 1/m/m/m,   grid_name); // number conc of ice nuclei due to immersion freezing (hetero nuc) [1/m3]
  add_field<Computed>("nidep",scalar3d_layout_mid , 1/m/m/m,   grid_name); // number conc of ice nuclei due to deposition nucleation (hetero nuc)[1/m3]
  add_field<Computed>("nimey",scalar3d_layout_mid , 1/m/m/m,   grid_name); // number conc of ice nuclei due to meyers deposition [1/m3]
  add_field<Computed>("naai_hom",scalar3d_layout_mid , n_unit, grid_name); // number of activated aerosol for ice nucleation (homogeneous freezing only) [#/kg]
  add_field<Computed>("naai",scalar3d_layout_mid , n_unit, grid_name); // number of activated aerosol for ice nucleation[#/kg]
  add_field<Computed>("qcld",scalar3d_layout_mid , n_unit, grid_name); // cloud droplet number mixing ratio [#/kg]
  add_field<Computed>("ptend_q", FieldLayout{ {COL, LEV, CMP}, {ncol_, nlev_, mam4::ndrop::nvar_ptend_q}}, n_unit, grid_name); // tendencies for interstitial and cloud borne aerosols [#/kg]
  add_field<Computed>("tendnd",scalar3d_layout_mid , n_unit/s, grid_name); // tendency in droplet number mixing ratio [#/kg/s]
  add_field<Computed>("factnum", FieldLayout{ {COL, LEV, CMP}, {ncol_, nlev_, mam4::AeroConfig::num_modes()}}, nondim, grid_name); // activation fraction for aerosol number [fraction]
  add_field<Computed>("ndropcol",scalar3d_layout_mid , n_unit/s, grid_name); // 
  add_field<Computed>("ndropmix",scalar3d_layout_mid , n_unit/s, grid_name); // droplet number mixing ratio tendency due to mixing [#/kg/s]
  add_field<Computed>("nsource",scalar3d_layout_mid , n_unit/s, grid_name); // droplet number mixing ratio source tendency [#/kg/s]
  add_field<Computed>("wtke",scalar3d_layout_mid , n_unit/s, grid_name); // 
  add_field<Computed>("ccn", FieldLayout{ {COL, LEV, CMP}, {ncol_, nlev_, mam4::ndrop::psat}}, n_unit, grid_name); //number conc of aerosols activated at supersat [#/m^3]
                   //      note:  activation fraction fluxes are defined as
                   //     fluxn = [flux of activated aero. number into cloud [#/m^2/s]]
                   //           / [aero. number conc. in updraft, just below cloudbase [#/m^3]]
  add_field<Computed>("coltend", FieldLayout{ {COL, LEV, CMP}, {ncol_, nlev_, mam4::ndrop::ncnst_tot} } , nondim, grid_name); // column tendency for diagnostic output
  add_field<Computed>("coltend_cw", FieldLayout{ {COL, LEV, CMP}, {ncol_, nlev_, mam4::ndrop::ncnst_tot} } , nondim, grid_name); // column tendency 


  //MUST FIXME: The aerosols has a wet mixing ratio, we should convert that to dry

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    //interstitial aerosol tracers of interest: number (n) mixing ratios
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit, grid_name, "tracers");

    //cloudborne aerosol tracers of interest: number (n) mixing ratios
    const char* cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);

    //NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are NOT advected
    add_field<Updated>(cld_nmr_field_name, scalar3d_layout_mid, n_unit, grid_name); 

    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char* int_mmr_field_name = mam_coupling::int_aero_mmr_field_name(m, a);
      if (strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name, "tracers");
      }
      
      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const char* cld_mmr_field_name = mam_coupling::cld_aero_mmr_field_name(m, a);
      if (strlen(cld_mmr_field_name) > 0) {
        //NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are NOT advected
        add_field<Updated>(cld_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name);
      }
    }
  }

  //Inputs (atmospheric quantities) for aci codes that existed in PBUF in EAM
  //These outputs should come from the cloud macrophysics process (e.g., SHOC)
  auto m2 = m*m;
  m2.set_string("m^2");
  auto s2 = s*s;
  s2.set_string("s^2");

  //MUST FIXME: w_sec,  is at OLD time step; strat_cld_frac and liq_strat_cld_frac may also need OLD time
  add_field<Required>("w_sec",              scalar3d_layout_int, m2/s2,  grid_name); // Vertical velocity variance (wp2) at midpoints

  add_field<Required>("strat_cld_frac",     scalar3d_layout_mid, nondim, grid_name); // Stratiform cloud fraction at midpoints
  add_field<Required>("liq_strat_cld_frac", scalar3d_layout_mid, nondim, grid_name); // Liquid stratiform cloud fraction  at midpoints
  add_field<Required>("kvh",                scalar3d_layout_int, m2/s, grid_name); // Eddy diffusivity for heat
  
  // Layout for 4D (2d horiz X 1d vertical x number of modes) variables
  num_aero_modes_ = mam_coupling::num_aero_modes();
  FieldLayout scalar4d_layout_mid{ {COL, LEV, NUM_MODES}, {ncol_, nlev_, num_aero_modes_} }; // mid points
  add_field<Required>("dgnum", scalar4d_layout_mid, m, grid_name); // dry diameter of aerosols
  

  /*NOTE on other inputs for the aci process:
  1. reciprocal of pseudo_density (rpdel): computed from the pseudo_density
  2. geopotential height at midpoints: computed geopotential height at interfaces, which inturn is computed using
     pseudo_density, p_mid, T_mid and qv_mid (see dry_static_energy.cpp's "compute_diagnostic_impl" function).
     qv_mid can be obtained from "get_field_in" call*/

}

void MAMAci::initialize_impl(const RunType run_type) {
    m_atm_logger->log(ekat::logger::LogLevel::info,"Calling ACI init");
    /*
      NOTE: All other inputs should follow the way "pseudo_density" is initialized
    */
    // set atmosphere state data
    pdel_ = get_field_in("pseudo_density").get_view<const Real**>();
    omega_ = get_field_in("omega").get_view<const Real**>();
    p_mid_ = get_field_in("p_mid").get_view<const Real**>();
    p_int_ = get_field_in("p_int").get_view<const Real**>();
    T_mid_ = get_field_in("T_mid").get_view<const Real**>();
    w_sec_ = get_field_in("w_sec").get_view<const Real**>();
    qv_dry_ = get_field_in("qv").get_view<const Real**>();
    w_updraft_ = get_field_in("w_updraft").get_view<const Real**>(); // updraft velocity [m/s]
								      //
    stratiform_cloud_fraction_ = get_field_out("stratiform_cloud_fraction").get_view<Real**>();
    activation_fraction_accum_idx_ = get_field_out("activation_fraction_accum").get_view<Real**>();
    activation_fraction_coarse_idx_ = get_field_out("activation_fraction_coarse").get_view<Real**>();

    // FIXME: These should all come from the intput file.  Might need to be added to the input file?
    Kokkos::resize(qc_accum_dst_, ncol_, nlev_); // qc_accum_dst_ = get_field_in("qc_accum_dst").get_view<Real**>();
    Kokkos::resize(qc_accum_so4_, ncol_, nlev_); // qc_accum_so4_ = get_field_in("qc_accum_so4").get_view<Real**>();
    Kokkos::resize(qc_accum_mom_, ncol_, nlev_); // qc_accum_mom_ = get_field_in("qc_accum_mom").get_view<Real**>();
    Kokkos::resize(qc_accum_nacl_, ncol_, nlev_); // qc_accum_mom_ = get_field_in("qc_accum_mom").get_view<Real**>();
    Kokkos::resize(qc_accum_bc_, ncol_, nlev_); // qc_accum_bc_ = get_field_in("qc_accum_bc").get_view<Real**>();
    Kokkos::resize(qc_accum_pom_, ncol_, nlev_); // qc_accum_pom_ = get_field_in("qc_accum_pom").get_view<Real**>();
    Kokkos::resize(qc_accum_soa_, ncol_, nlev_); // qc_accum_soa_ = get_field_in("qc_accum_soa").get_view<Real**>();
    Kokkos::resize(qi_accum_dst_, ncol_, nlev_); // qi_accum_dst_ = get_field_in("qi_accum_dst").get_view<Real**>();
    Kokkos::resize(qi_accum_so4_, ncol_, nlev_); // qi_accum_so4_ = get_field_in("qi_accum_so4").get_view<Real**>();
    Kokkos::resize(qi_accum_mom_, ncol_, nlev_); // qi_accum_mom_ = get_field_in("qi_accum_mom").get_view<Real**>();
    Kokkos::resize(qi_accum_bc_, ncol_, nlev_); // qi_accum_bc_ = get_field_in("qi_accum_bc").get_view<Real**>();
    Kokkos::resize(qi_accum_pom_, ncol_, nlev_); // qi_accum_pom_ = get_field_in("qi_accum_pom").get_view<Real**>();
    Kokkos::resize(qi_accum_soa_, ncol_, nlev_); // qi_accum_soa_ = get_field_in("qi_accum_soa").get_view<Real**>();
    Kokkos::resize(qc_coarse_dst_, ncol_, nlev_); // qc_coarse_dst_ = get_field_in("qc_coarse_dst").get_view<Real**>();
    Kokkos::resize(qc_coarse_nacl_, ncol_, nlev_); // qc_coarse_nacl_ = get_field_in("qc_coarse_nacl").get_view<Real**>();
    Kokkos::resize(qc_coarse_mom_, ncol_, nlev_); // ci_coarse_mom_ = get_field_in("qc_coarse_mom").get_view<Real**>();
    Kokkos::resize(qc_coarse_bc_, ncol_, nlev_); // qc_coarse_bc_ = get_field_in("qc_coarse_bc").get_view<Real**>();
    Kokkos::resize(qc_coarse_pom_, ncol_, nlev_); // qc_coarse_pom_ = get_field_in("qc_coarse_pom").get_view<Real**>();
    Kokkos::resize(qc_coarse_soa_, ncol_, nlev_); // qc_coarse_soa_ = get_field_in("qc_coarse_soa").get_view<Real**>();
						  //
    qi_coarse_dst_ = get_field_in("qi_coarse_dst").get_view<Real**>();
    qi_coarse_nacl_ = get_field_in("qi_coarse_nacl").get_view<Real**>();
    qi_coarse_so4_ = get_field_in("qi_coarse_so4").get_view<Real**>();
    qi_coarse_mom_ = get_field_in("qi_coarse_mom").get_view<Real**>();
    qi_coarse_bc_ = get_field_in("qi_coarse_bc").get_view<Real**>();
    qi_coarse_pom_ = get_field_in("qi_coarse_pom").get_view<Real**>();
    qi_coarse_soa_ = get_field_in("qi_coarse_soa").get_view<Real**>();

    Kokkos::resize(qi_pcarbon_bc_, ncol_, nlev_); // qi_pcarbon_bc_ = get_field_in("qi_pcarbon_bc").get_view<Real**>();
    Kokkos::resize(qi_pcarbon_mom_, ncol_, nlev_); // qi_pcarbon_mom_ = get_field_in("qi_pcarbon_mom").get_view<Real**>();
    Kokkos::resize(qi_pcarbon_pom_, ncol_, nlev_); // qi_pcarbon_pom_ = get_field_in("qi_pcarbon_pom").get_view<Real**>();
    Kokkos::resize(nc_accum_, ncol_, nlev_); // nc_accum_ = get_field_in("nc_accum").get_view<Real**>();
    Kokkos::resize(ni_accum_, ncol_, nlev_); // ni_accum_ = get_field_in("ni_accum").get_view<Real**>();
    Kokkos::resize(ni_coarse_, ncol_, nlev_); // ni_coarse_ = get_field_in("ni_coarse").get_view<Real**>();
    Kokkos::resize(ni_aitken_, ncol_, nlev_); // ni_aitken_ = get_field_in("ni_aitken").get_view<Real**>();
    Kokkos::resize(zm_, ncol_, nlev_); // zm_ =  get_field_in("zm").get_view<const Real**>();
    Kokkos::resize(state_q_, ncol_, nlev_, mam4::ndrop::nvars); // state_q_ =  get_field_in("state_q").get_view<const Real**>();
    Kokkos::resize(ncldwtr_, ncol_, nlev_); // ncldwtr_ =  get_field_in("ncldwtr").get_view<const Real**>();
    Kokkos::resize(cldo_, ncol_, nlev_); // cldo_ =  get_field_in("cldo").get_view<const Real**>();
    qqcw_input_ =  get_field_in("qqcw").get_view<const Real***>();

    for (int i=0; i<mam4::ndrop::ncnst_tot; ++i) {
      Kokkos::resize(qqcw_[i], ncol_, nlev_);
      Kokkos::resize(coltend_[i],  ncol_, nlev_);
      Kokkos::resize(coltend_cw_[i],  ncol_, nlev_);
    }
    for (int i=0; i<mam4::ndrop::nvar_ptend_q; ++i) {
      Kokkos::resize(ptend_q_[i], ncol_, nlev_);
    }
    for (int i=0; i<mam4::ndrop::pver; ++i) {
      for (int j=0; j<2; ++j) {
        Kokkos::resize(raercol_cw_[i][j], ncol_, mam4::ndrop::ncnst_tot);
        Kokkos::resize(raercol_[i][j], ncol_, mam4::ndrop::ncnst_tot);
      }
    }
    Kokkos::resize(hetfrz_immersion_nucleation_tend_, ncol_, nlev_);
    Kokkos::resize(hetfrz_contact_nucleation_tend_, ncol_, nlev_);
    Kokkos::resize(hetfrz_depostion_nucleation_tend_, ncol_, nlev_);

    for (int i=0; i<42; ++i) 
      Kokkos::resize(diagnostic_scratch_[i],ncol_, nlev_);

    Kokkos::resize(nact_,  ncol_, nlev_, mam4::AeroConfig::num_modes());
    Kokkos::resize(mact_,  ncol_, nlev_, mam4::AeroConfig::num_modes());

    cldfrac_ = get_field_in("cldfrac_tot").get_view<const Real**>();
    dgnum_ = get_field_out("dgnum").get_view<const Real***>();
    nihf_ = get_field_out("nihf").get_view<Real**>();
    niim_ = get_field_out("niim").get_view<Real**>();
    nidep_ = get_field_out("nidep").get_view<Real**>();
    nimey_ = get_field_out("nimey").get_view<Real**>();
    naai_hom_ = get_field_out("naai_hom").get_view<Real**>();
    naai_ = get_field_out("naai").get_view<Real**>();
    qcld_ = get_field_out("qcld").get_view<Real**>();
    ptend_q_output_ = get_field_out("ptend_q").get_view<Real***>();
    tendnd_ = get_field_out("tendnd").get_view<Real**>();
    factnum_ = get_field_out("factnum").get_view<Real***>();
    ndropcol_ = get_field_out("ndropcol").get_view<Real**>();
    ndropmix_ = get_field_out("ndropmix").get_view<Real**>();
    nsource_ = get_field_out("nsource").get_view<Real**>();
    wtke_ = get_field_out("wtke").get_view<Real**>();
    ccn_ = get_field_out("ccn").get_view<Real***>();
    coltend_outp_ = get_field_out("coltend").get_view<Real***>();
    coltend_cw_outp_ = get_field_out("coltend_cw").get_view<Real***>();



    liqcldf_ = get_field_in("liq_strat_cld_frac").get_view<const Real**>();
    qc_ =  get_field_in("qc").get_view<const Real**>();
    qi_ =  get_field_in("qi").get_view<const Real**>();
    nc_ =  get_field_in("nc").get_view<const Real**>();
    kvh_ =  get_field_in("kvh").get_view<const Real**>();


   mam4::AeroConfig aero_config;
   // configure the nucleation parameterization
   mam4::NucleateIce::Config nucleate_ice_config;
   nucleate_ice_.init(aero_config, nucleate_ice_config);

   // configure the heterogeneous freezing parameterization
   mam4::Hetfrz::Config hetfrz_config;
   hetfrz_.init(aero_config, hetfrz_config);
    /*
      NOTE: All derived variables (like rpdel and geopotential height) should be computed in
      preprocess struct
    */
    //preprocess_.set_variables(pdel_);

}

void MAMAci::run_impl(const double dt) {
  m_atm_logger->log(ekat::logger::LogLevel::info,"calling ACI run");

  //----------------------------------------------------------
  // Convert from omega to w (vertical velocity)
  // Negative omega means rising motion
  //---------------------------------------------------------

  
  // Alias member variables
  auto rho = rho_;
  auto tke = tke_;
  auto T_mid = T_mid_;
  auto p_mid = p_mid_;
  auto p_int = p_int_;
  auto pdel = pdel_;
  auto wsub = wsub_;
  auto wsubice = wsubice_;
  auto qv_dry = qv_dry_;
  auto cldfrac = cldfrac_;
  auto w_updraft = w_updraft_;
  const auto nlev = nlev_;
  const auto top_lev = top_lev_;

  auto qc_accum_dst = qc_accum_dst_;
  auto qc_accum_so4 = qc_accum_so4_;
  auto qc_accum_mom = qc_accum_mom_;
  auto qc_accum_nacl = qc_accum_nacl_;
  auto qc_accum_bc = qc_accum_bc_;
  auto qc_accum_pom = qc_accum_pom_;
  auto qc_accum_soa = qc_accum_soa_;

  auto qi_accum_dst = qi_accum_dst_;
  auto qi_accum_so4 = qi_accum_so4_;
  auto qi_accum_mom = qi_accum_mom_;
  auto qi_accum_bc = qi_accum_bc_;
  auto qi_accum_pom = qi_accum_pom_;
  auto qi_accum_soa = qi_accum_soa_;

  auto qc_coarse_dst = qc_coarse_dst_;
  auto qc_coarse_nacl = qc_coarse_nacl_;
  auto qc_coarse_mom = qc_coarse_mom_;
  auto qc_coarse_bc = qc_coarse_bc_;
  auto qc_coarse_pom = qc_coarse_pom_;
  auto qc_coarse_soa = qc_coarse_soa_;
  auto qi_coarse_dst = qi_coarse_dst_;
  auto qi_coarse_nacl = qi_coarse_nacl_;
  auto qi_coarse_so4 = qi_coarse_so4_;
  auto qi_coarse_mom = qi_coarse_mom_;
  auto qi_coarse_bc = qi_coarse_bc_;
  auto qi_coarse_pom = qi_coarse_pom_;
  auto qi_coarse_soa = qi_coarse_soa_;
  auto qi_pcarbon_bc = qi_pcarbon_bc_;
  auto qi_pcarbon_mom = qi_pcarbon_mom_;
  auto qi_pcarbon_pom = qi_pcarbon_pom_;
  auto nc_accum = nc_accum_;
  auto ni_accum = ni_accum_;
  auto ni_coarse = ni_coarse_;
  auto ni_aitken = ni_aitken_;
  auto nihf = nihf_;
  auto niim = niim_;
  auto nidep = nidep_;
  auto nimey = nimey_;
  auto naai_hom = naai_hom_;
  auto naai = naai_;
  auto liqcldf = liqcldf_;
  auto qc = qc_;
  auto qi = qi_;
  auto nc = nc_;
  auto lcldn = lcldn_;
  auto lcldo = lcldo_;
  auto aitken_dry_dia = aitken_dry_dia_;
  auto rpdel = rpdel_;
  auto zm = zm_;
  auto state_q = state_q_;
  auto ncldwtr = ncldwtr_;
  auto kvh = kvh_;
  auto qcld = qcld_;
  auto tendnd = tendnd_;
  auto factnum = factnum_;
  auto ndropcol =  ndropcol_ ;
  auto ndropmix =  ndropmix_ ;
  auto nsource = nsource_ ;
  auto wtke =   wtke_ ;
  auto cldo = cldo_;
  auto ccn = ccn_;
  auto coltend_outp = coltend_outp_;
  auto coltend_cw_outp = coltend_cw_outp_;
  auto nact = nact_;
  auto mact = mact_;
  auto eddy_diff = dropmixnuc_scratch_mem_[0];
  auto zn = dropmixnuc_scratch_mem_[1];
  auto csbot = dropmixnuc_scratch_mem_[2];
  auto zs = dropmixnuc_scratch_mem_[3];
  auto overlapp = dropmixnuc_scratch_mem_[4];
  auto overlapm = dropmixnuc_scratch_mem_[5];
  auto eddy_diff_kp = dropmixnuc_scratch_mem_[6];
  auto eddy_diff_km = dropmixnuc_scratch_mem_[7];
  auto qncld = dropmixnuc_scratch_mem_[8];
  auto srcn = dropmixnuc_scratch_mem_[9];
  auto source = dropmixnuc_scratch_mem_[10];
  auto dz = dropmixnuc_scratch_mem_[11];
  auto csbot_cscen = dropmixnuc_scratch_mem_[12];
  auto raertend = dropmixnuc_scratch_mem_[13];
  auto qqcwtend = dropmixnuc_scratch_mem_[14];

  auto hetfrz_immersion_nucleation_tend =   hetfrz_immersion_nucleation_tend_; 
  auto hetfrz_contact_nucleation_tend =     hetfrz_contact_nucleation_tend_; 
  auto hetfrz_depostion_nucleation_tend =   hetfrz_depostion_nucleation_tend_; 
  auto stratiform_cloud_fraction = stratiform_cloud_fraction_;
  auto activation_fraction_accum_idx = activation_fraction_accum_idx_;
  auto activation_fraction_coarse_idx = activation_fraction_coarse_idx_;
  auto diagnostic_scratch = diagnostic_scratch_; 

  haero::ThreadTeamPolicy team_policy(ncol_, Kokkos::AUTO);

  copy_scream_array_to_mam4xx(team_policy, qqcw_, qqcw_input_, nlev_, mam4::ndrop::ncnst_tot);

  // All the inputs are available to compute w0 and rho
  compute_w0_and_rho(team_policy, w0_, rho_, p_mid_, T_mid_, omega_, top_lev_, nlev_);

  compute_tke_using_w_sec(team_policy, tke_, w_sec_, nlev_);

  compute_subgrid_scale_velocities(team_policy, wsub_, wsubice_, wsig_, tke_, wsubmin, top_lev_, nlev_);
  
  compute_subgrid_mean_updraft_velocities(team_policy, w2_, w0_, wsig_, nlev_);
  
  compute_aitken_dry_diameter(team_policy, aitken_dry_dia_, dgnum_, top_lev_);

  compute_nucleate_ice_tendencies(team_policy, nlev_);
      //T_mid_, p_mid_, qv_dry_, cldfrac_, w_updraft_
    // qi_coarse_dst_, qi_coarse_nacl_, qi_coarse_so4_, qi_coarse_mom_, qi_coarse_bc_, qi_coarse_pom_, qi_coarse_soa_, 

  //-------------------------------------------------------------
  // Get number of activated aerosol for ice nucleation (naai)
  // from ice nucleation
  //-------------------------------------------------------------
  using view_1d = typename KokkosTypes<DefaultDevice>::template view_1d<Real>;
  view_1d dummy("DummyView", nlev);
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();

    // Set up an atmosphere, surface, diagnostics, pronostics and tendencies class.
    Real pblh = 0;
    haero::Atmosphere atmos(
      nlev, 
      ekat::subview(T_mid, icol),
      ekat::subview(p_mid, icol),
      ekat::subview(qv_dry, icol),
      dummy, dummy, dummy, dummy, dummy, dummy,
      ekat::subview(cldfrac, icol),
      ekat::subview(w_updraft, icol), pblh);
    // set surface state data
    haero::Surface surf{};
    mam4::Prognostics progs(nlev);
    const int coarse_idx = static_cast<int>(mam4::ModeIndex::Coarse);
    const int aitken_idx = static_cast<int>(mam4::ModeIndex::Aitken);
    const int dst_idx = static_cast<int>(mam4::AeroId::DST);

    progs.q_aero_i[coarse_idx][dst_idx] = ekat::subview(qi_coarse_dst, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::NaCl)] = ekat::subview(qi_coarse_nacl, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::SO4)] = ekat::subview(qi_coarse_so4, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::MOM)] = ekat::subview(qi_coarse_mom, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::BC)] = ekat::subview(qi_coarse_bc, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::POM)] = ekat::subview(qi_coarse_pom, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::SOA)] = ekat::subview(qi_coarse_soa, icol);
    progs.n_mode_i[coarse_idx] = ekat::subview(ni_coarse, icol);
    progs.n_mode_i[aitken_idx] = ekat::subview(ni_aitken, icol);

    // nucleation doesn't use any diagnostics, so it's okay to leave this alone
    // for now
    mam4::Diagnostics diags(nlev);
    diags.dry_geometric_mean_diameter_i[aitken_idx] = ekat::subview(aitken_dry_dia, icol);

    // These are the fields that are updated. Taking subviews means that
    // the nihf, niim, nidep, nimey, naai_hom, and naai filds are updated 
    // in nucleate_ice_.compute_tendencies.
    diags.icenuc_num_hetfrz = ekat::subview(nihf, icol);
    diags.icenuc_num_immfrz = ekat::subview(niim, icol);
    diags.icenuc_num_depnuc = ekat::subview(nidep, icol);
    diags.icenuc_num_meydep = ekat::subview(nimey, icol);

    // naai and naai_hom are the outputs needed for nucleate_ice and these are not tendencies.
    diags.num_act_aerosol_ice_nucle_hom = ekat::subview(naai_hom, icol);
    diags.num_act_aerosol_ice_nucle = ekat::subview(naai, icol);

    // grab views from the buffer to store tendencies, not used as all values are store in diags above.
    const mam4::Tendencies tends(nlev);
    const mam4::AeroConfig aero_config;
    const Real t=0, dt=0;
    nucleate_ice_.compute_tendencies(aero_config, team, t, dt, atmos, surf, progs, diags, tends);

    /*
      NOTE:"state_q" is a combination of subset of tracers added by "int_mmr_field_name" and "int_nmr_field_name".
      Only output we care about is "naai", "naai_hom" is never used anywhere

      Fortran code:
      call nucleate_ice_cam_calc(ncol, lchnk, temperature, state_q, pmid, &      ! input
           rho, wsubice, strat_cld_frac, dgnum, &          ! input
           naai, naai_hom)                                 ! output
    */
  });
    //-------------------------------------------------------------
    // Get old and new liquid cloud fractions when amount of cloud 
    // is above qsmall threshold value
    //-------------------------------------------------------------

  const Real dtmicro = dtmicro_;
   static constexpr auto qsmall = 1e-18; //cut-off for cloud amount (ice or liquid)
   
    /*
    MUST FIXME NOTE: We need old and new liquid cloud fractions here.
    We have the new liquid cloud fraction (liq_strat_cld_frac) but we need to 
    store the old (liq_strat_cld_frac_old) before we call SHOC. For now, we will make
    a note of it and use the new cloud fraction for the 
    old cloud fraction.
    */
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, top_lev), KOKKOS_LAMBDA(int kk) { 
      const Real qcld = qc(icol,kk) + qi(icol,kk);
      if (qcld > qsmall) {
        lcldn(icol,kk)=liqcldf(icol,kk);
        lcldo(icol,kk)=liqcldf(icol,kk); // FIXME should be liqcldf_old
      }
    });
  });
    /*
    Fortran code:
    liqcldf(:ncol,:pver) = liq_strat_cld_frac(:ncol,:pver)
    lcldn = 0._r8
    lcldo = 0._r8
    do kk = top_lev, pver
       do icol = 1, ncol
          qcld = qc(icol,kk) + qi(icol,kk)
          if (qcld > qsmall) then
             lcldn(icol,kk)=liqcldf(icol,kk)
             lcldo(icol,kk)=liqcldfo(icol,kk)
          end if
       end do
    end do
    */

    //-------------------------------------------------------------
    // Save cloud borne aerosols to be used in the heterozenous 
    // freezing before they are changed by the droplet activation
    // process. This is only a select subset of cloud borne 
    // aerosols, not all the cloud borne aerosols.
    //-------------------------------------------------------------

    /*NOTE: We probably need to store indices for the select few
    cloud borne aerosols
  
    Fortran code:
    lchnk_zb = lchnk - begchunk
    ! save copy of cloud borne aerosols for use in heterogeneous freezing before
    !we change it in dropmixnuc
    do ispec = 1, ncnst
       call pbuf_get_field(pbuf, hetfrz_aer_spec_idx(ispec), ptr2d)
       aer_cb(:ncol,:,ispec,lchnk_zb) = ptr2d(:ncol,:)
       aer_cb(:ncol,:,ispec,lchnk_zb) = aer_cb(:ncol,:,ispec,lchnk_zb) * rho(:ncol,:)
    enddo
    */


    //-------------------------------------------------------------
    // Compute activated fraction of aerosols
    //-------------------------------------------------------------
    const int ntot_amode = mam4::AeroConfig::num_modes();
    const int maxd_aspectype = mam4::ndrop::maxd_aspectype;
    const int nspec_max = mam4::ndrop::nspec_max;
    int nspec_amode[ntot_amode];
    int lspectype_amode[maxd_aspectype][ntot_amode];
    int lmassptr_amode[maxd_aspectype][ntot_amode];
    Real specdens_amode[maxd_aspectype];
    Real spechygro[maxd_aspectype];
    int numptr_amode[ntot_amode];
    int mam_idx[ntot_amode][nspec_max];
    int mam_cnst_idx[ntot_amode][nspec_max];
    mam4::ndrop::get_e3sm_parameters(
        nspec_amode, lspectype_amode, lmassptr_amode, numptr_amode,
        specdens_amode, spechygro, mam_idx, mam_cnst_idx);
    Real exp45logsig[ntot_amode],
        alogsig[ntot_amode],
        num2vol_ratio_min_nmodes[ntot_amode],
        num2vol_ratio_max_nmodes[ntot_amode] = {};
    Real aten = 0;
    mam4::ndrop::ndrop_init(exp45logsig, alogsig, aten,
                      num2vol_ratio_min_nmodes,  // voltonumbhi_amode
                      num2vol_ratio_max_nmodes); // voltonumblo_amode
    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) { 
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, nlev), KOKKOS_LAMBDA(int kk) { 
        rpdel(icol,kk) = 1/pdel(icol,kk);
      });
    });



    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
      const int icol = team.league_rank();    

      mam4::ndrop::View1D raercol_cw[mam4::ndrop::pver][2];
      mam4::ndrop::View1D raercol[mam4::ndrop::pver][2];
      for (int i=0; i<mam4::ndrop::pver; ++i) {
        for (int j=0; j<2; ++j) {
          raercol_cw[i][j] = ekat::subview(raercol_cw_[i][j], icol);
          raercol[i][j] = ekat::subview(raercol_[i][j], icol);
        }
      }
      mam4::ColumnView qqcw[mam4::ndrop::ncnst_tot];
      for (int i=0; i<mam4::ndrop::ncnst_tot; ++i) {
        qqcw[i] = ekat::subview(qqcw_[i], icol);
      }
      mam4::ColumnView ptend_q[mam4::ndrop::nvar_ptend_q];
      for (int i=0; i<mam4::ndrop::nvar_ptend_q; ++i) {
        ptend_q[i] = ekat::subview(ptend_q_[i], icol);
      }
      
      mam4::ColumnView coltend[mam4::ndrop::ncnst_tot], coltend_cw[mam4::ndrop::ncnst_tot];
      for (int i=0; i<mam4::ndrop::ncnst_tot; ++i) {
        coltend[i] = ekat::subview(coltend_[i], icol);
        coltend_cw[i] = ekat::subview(coltend_cw_[i], icol);
      }


      mam4::ndrop::dropmixnuc(
       team, dtmicro,
       ekat::subview(T_mid, icol), 
       ekat::subview(p_mid, icol),
       ekat::subview(p_int, icol),
       ekat::subview(pdel, icol),
       ekat::subview(rpdel, icol),
       ekat::subview(zm, icol), //  ! in zm[kk] - zm[kk+1], for pver zm[kk-1] - zm[kk]
       ekat::subview(state_q, icol),
       ekat::subview(ncldwtr, icol),
       ekat::subview(kvh, icol), // kvh[kk+1]
       ekat::subview(cldfrac, icol),
       lspectype_amode, specdens_amode, spechygro, lmassptr_amode,
       num2vol_ratio_min_nmodes, num2vol_ratio_max_nmodes, numptr_amode,
       nspec_amode, exp45logsig, alogsig, aten, mam_idx, mam_cnst_idx,
       ekat::subview(qcld, icol),// out
       ekat::subview(wsub, icol), // in
       ekat::subview(cldo, icol), // in
       qqcw, // inout
       ptend_q, 
       ekat::subview(tendnd, icol), 
       ekat::subview(factnum,  icol), 
       ekat::subview(ndropcol, icol), 
       ekat::subview(ndropmix, icol), 
       ekat::subview(nsource, icol), 
       ekat::subview(wtke, icol), 
       ekat::subview(ccn, icol), 
       coltend, 
       coltend_cw, 
       raercol_cw, 
       raercol, 
       ekat::subview(nact, icol),
       ekat::subview(mact, icol),
       ekat::subview(eddy_diff, icol), 
       // work arrays
       ekat::subview(zn,  icol), 
       ekat::subview(csbot,  icol), 
       ekat::subview(zs,  icol), 
       ekat::subview(overlapp,  icol), 
       ekat::subview(overlapm,  icol), 
       ekat::subview(eddy_diff_kp,  icol), 
       ekat::subview(eddy_diff_km,  icol), 
       ekat::subview(qncld,  icol), 
       ekat::subview(srcn, icol), 
       ekat::subview(source,  icol), 
       ekat::subview(dz,  icol), 
       ekat::subview(csbot_cscen,  icol), 
       ekat::subview(raertend,  icol), 
       ekat::subview(qqcwtend,  icol) );
  });

  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, nlev), KOKKOS_LAMBDA(int kk) { 
      for (int i=0; i<mam4::ndrop::nvar_ptend_q; ++i) {
         ptend_q_output_(icol, kk, i) = ptend_q_[i](icol,kk);
      }
      for (int i=0; i<mam4::ndrop::ncnst_tot; ++i) {
        coltend_outp(icol, kk, i) = coltend_[i](icol,kk);
        coltend_cw_outp(icol, kk, i) = coltend_cw_[i](icol,kk);
      }
    });
  });

    std::cout<<"pdel_ in run_impl is:"<<pdel_(0,0)<<std::endl;
    /*
    NOTE: "deltain" is the model time step. "state_q" is a combination of tracers
    fields with "int_mmr_field_name" and "int_nmr_field_name". "z_mid" is computed.
    "qqcw" is the combination of cld_mmr_field_name and cld_nmr_field_name.
    The output "ptend" will have tendencies for interstitial and cloud borne aerosols.

    Fortan code:
      call dropmixnuc(lchnk, ncol, deltatin, T_mid, p_mid, p_int, p_del, rpdel, z_mid, &  ! in
         state_q, nc, kvh, wsub, lcldn, lcldo, &  ! in
         qqcw, &  ! inout
         ptend, nctend_mixnuc, factnum)  !out
    */


    //-------------------------------------------------------------
    // Heterogeneous freezing
    // frzimm, frzcnt, frzdep are the outputs of 
    // hetfrz_classnuc_cam_calc used by the microphysics (e.g. p3)
    //-------------------------------------------------------------

  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();

    // Set up an atmosphere, surface, diagnostics, pronostics and tendencies class.
    Real pblh = 0;
    haero::Atmosphere atmos(
      nlev, 
      ekat::subview(T_mid, icol),
      ekat::subview(p_mid, icol),
      dummy, 
      ekat::subview(qc, icol),
      ekat::subview(nc, icol),
      dummy, dummy, dummy, dummy, dummy, dummy,
      pblh);
    // set surface state data
    haero::Surface surf{};
    mam4::Prognostics progs(nlev);
    const int accum_idx  = static_cast<int>(mam4::ModeIndex::Accumulation);
    const int coarse_idx = static_cast<int>(mam4::ModeIndex::Coarse);
    const int pcarbon_idx = static_cast<int>(mam4::ModeIndex::PrimaryCarbon);

    progs.q_aero_c[accum_idx][static_cast<int>(mam4::AeroId::BC)] = ekat::subview(qc_accum_bc, icol);
    progs.q_aero_c[accum_idx][static_cast<int>(mam4::AeroId::DST)] = ekat::subview(qc_accum_dst, icol);
    progs.q_aero_c[accum_idx][static_cast<int>(mam4::AeroId::MOM)] = ekat::subview(qc_accum_mom, icol);
    progs.q_aero_c[accum_idx][static_cast<int>(mam4::AeroId::POM)] = ekat::subview(qc_accum_pom, icol);
    progs.q_aero_c[accum_idx][static_cast<int>(mam4::AeroId::NaCl)] = ekat::subview(qc_accum_nacl, icol);
    progs.q_aero_c[accum_idx][static_cast<int>(mam4::AeroId::SO4)] = ekat::subview(qc_accum_so4, icol);
    progs.q_aero_c[accum_idx][static_cast<int>(mam4::AeroId::SOA)] = ekat::subview(qc_accum_soa, icol);

    progs.q_aero_i[accum_idx][static_cast<int>(mam4::AeroId::BC)] = ekat::subview(qi_accum_bc, icol);
    progs.q_aero_i[accum_idx][static_cast<int>(mam4::AeroId::DST)] = ekat::subview(qi_accum_dst, icol);
    progs.q_aero_i[accum_idx][static_cast<int>(mam4::AeroId::MOM)] = ekat::subview(qi_accum_mom, icol);
    progs.q_aero_i[accum_idx][static_cast<int>(mam4::AeroId::POM)] = ekat::subview(qi_accum_pom, icol);
    progs.q_aero_i[accum_idx][static_cast<int>(mam4::AeroId::SO4)] = ekat::subview(qi_accum_so4, icol);
    progs.q_aero_i[accum_idx][static_cast<int>(mam4::AeroId::SOA)] = ekat::subview(qi_accum_soa, icol);

    progs.q_aero_c[coarse_idx][static_cast<int>(mam4::AeroId::BC)] = ekat::subview(qc_coarse_bc, icol);
    progs.q_aero_c[coarse_idx][static_cast<int>(mam4::AeroId::DST)] = ekat::subview(qc_coarse_dst, icol);
    progs.q_aero_c[coarse_idx][static_cast<int>(mam4::AeroId::MOM)] = ekat::subview(qc_coarse_mom, icol);
    progs.q_aero_c[coarse_idx][static_cast<int>(mam4::AeroId::NaCl)] = ekat::subview(qc_coarse_nacl, icol);
    progs.q_aero_c[coarse_idx][static_cast<int>(mam4::AeroId::POM)] = ekat::subview(qc_coarse_pom, icol);
    progs.q_aero_c[coarse_idx][static_cast<int>(mam4::AeroId::SOA)] = ekat::subview(qc_coarse_soa, icol);

    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::BC)] = ekat::subview(qi_coarse_bc, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::DST)] = ekat::subview(qi_coarse_dst, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::MOM)] = ekat::subview(qi_coarse_mom, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::NaCl)] = ekat::subview(qi_coarse_nacl, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::POM)] = ekat::subview(qi_coarse_pom, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::SO4)] = ekat::subview(qi_coarse_so4, icol);
    progs.q_aero_i[coarse_idx][static_cast<int>(mam4::AeroId::SOA)] = ekat::subview(qi_coarse_soa, icol);

    progs.q_aero_i[pcarbon_idx][static_cast<int>(mam4::AeroId::BC)] = ekat::subview(qi_pcarbon_bc, icol);
    progs.q_aero_i[pcarbon_idx][static_cast<int>(mam4::AeroId::MOM)] = ekat::subview(qi_pcarbon_mom, icol);
    progs.q_aero_i[pcarbon_idx][static_cast<int>(mam4::AeroId::POM)] = ekat::subview(qi_pcarbon_pom, icol);

    progs.n_mode_c[accum_idx ] = ekat::subview(nc_accum, icol);

    progs.n_mode_i[accum_idx ] = ekat::subview(ni_accum, icol);
    progs.n_mode_i[coarse_idx] = ekat::subview(ni_coarse, icol);

    mam4::Diagnostics diags(nlev);
    diags.stratiform_cloud_fraction = ekat::subview(stratiform_cloud_fraction, icol);
    diags.activation_fraction[accum_idx] = ekat::subview(activation_fraction_accum_idx, icol);
    diags.activation_fraction[coarse_idx] = ekat::subview(activation_fraction_coarse_idx, icol);

    diags.hetfrz_immersion_nucleation_tend = ekat::subview(hetfrz_immersion_nucleation_tend, icol);
    diags.hetfrz_contact_nucleation_tend   = ekat::subview(hetfrz_contact_nucleation_tend, icol);
    diags.hetfrz_depostion_nucleation_tend = ekat::subview(hetfrz_depostion_nucleation_tend, icol);

    diags.bc_num = ekat::subview(diagnostic_scratch[0], icol);
    diags.dst1_num = ekat::subview(diagnostic_scratch[1], icol);
    diags.dst3_num = ekat::subview(diagnostic_scratch[2], icol);
    diags.bcc_num = ekat::subview(diagnostic_scratch[3], icol);
    diags.dst1c_num = ekat::subview(diagnostic_scratch[4], icol);
    diags.dst3c_num = ekat::subview(diagnostic_scratch[5], icol);
    diags.bcuc_num = ekat::subview(diagnostic_scratch[6], icol);
    diags.dst1uc_num = ekat::subview(diagnostic_scratch[7], icol);
    diags.dst3uc_num = ekat::subview(diagnostic_scratch[8], icol);
    diags.bc_a1_num = ekat::subview(diagnostic_scratch[0], icol);
    diags.dst_a1_num = ekat::subview(diagnostic_scratch[10], icol);
    diags.dst_a3_num = ekat::subview(diagnostic_scratch[11], icol);
    diags.bc_c1_num = ekat::subview(diagnostic_scratch[12], icol);
    diags.dst_c1_num = ekat::subview(diagnostic_scratch[13], icol);
    diags.dst_c3_num = ekat::subview(diagnostic_scratch[14], icol);
    diags.fn_bc_c1_num = ekat::subview(diagnostic_scratch[15], icol);
    diags.fn_dst_c1_num = ekat::subview(diagnostic_scratch[16], icol);
    diags.fn_dst_c3_num = ekat::subview(diagnostic_scratch[17], icol);
    diags.na500 = ekat::subview(diagnostic_scratch[18], icol);
    diags.totna500 = ekat::subview(diagnostic_scratch[19], icol);
    diags.freqimm = ekat::subview(diagnostic_scratch[20], icol);
    diags.freqcnt = ekat::subview(diagnostic_scratch[21], icol);
    diags.freqdep = ekat::subview(diagnostic_scratch[22], icol);
    diags.freqmix = ekat::subview(diagnostic_scratch[23], icol);
    diags.dstfrezimm = ekat::subview(diagnostic_scratch[24], icol);
    diags.dstfrezcnt = ekat::subview(diagnostic_scratch[25], icol);
    diags.dstfrezdep = ekat::subview(diagnostic_scratch[26], icol);
    diags.bcfrezimm = ekat::subview(diagnostic_scratch[27], icol);
    diags.bcfrezcnt = ekat::subview(diagnostic_scratch[28], icol);
    diags.bcfrezdep = ekat::subview(diagnostic_scratch[19], icol);
    diags.nimix_imm = ekat::subview(diagnostic_scratch[30], icol);
    diags.nimix_cnt = ekat::subview(diagnostic_scratch[31], icol);
    diags.nimix_dep = ekat::subview(diagnostic_scratch[32], icol);
    diags.dstnidep = ekat::subview(diagnostic_scratch[33], icol);
    diags.dstnicnt = ekat::subview(diagnostic_scratch[34], icol);
    diags.dstniimm = ekat::subview(diagnostic_scratch[35], icol);
    diags.bcnidep = ekat::subview(diagnostic_scratch[36], icol);
    diags.bcnicnt = ekat::subview(diagnostic_scratch[37], icol);
    diags.bcniimm = ekat::subview(diagnostic_scratch[38], icol);
    diags.numice10s = ekat::subview(diagnostic_scratch[39], icol);
    diags.numimm10sdst = ekat::subview(diagnostic_scratch[40], icol);
    diags.numimm10sbc = ekat::subview(diagnostic_scratch[41], icol);

    // naai and naai_hom are the outputs needed for nucleate_ice and these are not tendencies.
    diags.num_act_aerosol_ice_nucle_hom = ekat::subview(naai_hom, icol);
    diags.num_act_aerosol_ice_nucle = ekat::subview(naai, icol);

    // grab views from the buffer to store tendencies, not used as all values are store in diags above.
    const mam4::Tendencies tends(nlev);
    const mam4::AeroConfig aero_config;
    const Real t=0, dt=0;
    hetfrz_.compute_tendencies(aero_config, team, t, dt, atmos, surf, progs, diags, tends);

  });
    
    /*
    call hetfrz_classnuc_cam_calc(ncol, lchnk, temperature, pmid, rho, ast, &   ! in
         qc, nc, state_q, aer_cb(:,:,:,lchnk_zb), deltatin, factnum, & ! in
         frzimm, frzcnt, frzdep)  
    */       

}

void MAMAci::finalize_impl(){
  m_atm_logger->log(ekat::logger::LogLevel::info,"calling ACI final");
}

} // namespace scream
