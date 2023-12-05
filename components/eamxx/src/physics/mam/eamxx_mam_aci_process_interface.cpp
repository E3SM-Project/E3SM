#include <physics/mam/eamxx_mam_aci_process_interface.hpp>
namespace scream
{

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

  //FIXME: The following variables are namelist variables
  Real wsubmin = 1;

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

  Kokkos::resize(rho_, ncol_, nlev_);
  Kokkos::resize(w0_, ncol_, nlev_);
  Kokkos::resize(tke_, ncol_, nlev_+1);
  Kokkos::resize(wsub_, ncol_, nlev_);
  Kokkos::resize(wsubice_, ncol_, nlev_);
  Kokkos::resize(wsig_, ncol_, nlev_);
  Kokkos::resize(w2_, ncol_, nlev_);

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

  // atmospheric quantities
  add_field<Required>("qc",             scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud liquid mass mixing ratio [kg/kg]
  add_field<Required>("qi",             scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud ice mass mixing ratio [kg/kg]
  add_field<Required>("ni",             scalar3d_layout_mid, n_unit, grid_name, "tracers");// cloud liquid mass mixing ratio [kg/kg]
  add_field<Required>("T_mid",          scalar3d_layout_mid, K,      grid_name); // Temperature[K] at midpoints
  add_field<Required>("omega",          scalar3d_layout_mid, Pa/s,   grid_name); // Vertical pressure velocity [Pa/s] at midpoints
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa,     grid_name); // Total pressure [Pa] at midpoints
  add_field<Required>("p_int",          scalar3d_layout_int, Pa,     grid_name); // Total pressure [Pa] at interfaces
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,     grid_name); // Layer thickness(pdel) [Pa] at midpoints

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
  add_field<Required>("w_sec",              scalar3d_layout_mid, m2/s2,  grid_name); // Vertical velocity variance (wp2) at midpoints

  auto nondim = Units::nondimensional();
  add_field<Required>("strat_cld_frac",     scalar3d_layout_mid, nondim, grid_name); // Stratiform cloud fraction at midpoints
  add_field<Required>("liq_strat_cld_frac", scalar3d_layout_mid, nondim, grid_name); // Liquid stratiform cloud fraction  at midpoints
  add_field<Required>("kvh",                scalar3d_layout_mid, m2/s, grid_name); // Eddy diffusivity for heat
  
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
    T_mid_ = get_field_in("T_mid").get_view<const Real**>();
    w_sec_ = get_field_in("w_sec").get_view<const Real**>();

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

  //Get physical constants
  using C  = physics::Constants<Real>;
  
  static constexpr auto gravit = C::gravit; // Gravity [m/s2]
  static constexpr auto rair   = C::Rair;   // Gas constant for dry air [J/(kg*K) or J/Kg/K]

  // NOTE: All the inputs are available to compute w0
  auto team_policy = haero::ThreadTeamPolicy(ncol_, Kokkos::AUTO);
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, top_lev_), KOKKOS_LAMBDA(int kk) { 
      w0_(icol,kk) = 0;
      rho_(icol, kk) = -999.0;
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, top_lev_, nlev_), KOKKOS_LAMBDA(int kk) { 
      rho_(icol,kk) = p_mid_(icol,kk)/(rair*T_mid_(icol,kk));
      w0_(icol,kk) = -1.0*omega_(icol,kk)/(rho_(icol,kk)*gravit);
    });
  });

  //---------------------------------------------------------
  //Compute TKE using "w_sec"
  //---------------------------------------------------------
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    // FIXME Is this the correct boundary condition for tke at the surface?
    tke_(icol,nlev_) = 0;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, nlev_), KOKKOS_LAMBDA(int kk) { 
      tke_(icol,kk) = (3.0/2.0)*w_sec_(icol,kk);
    });
  });
   //---------------------------------------------------------
   // Compute subgrid scale velocities(wsub, wsig and wsubice)
   //---------------------------------------------------------
  
  // More refined computation of sub-grid vertical velocity
  // Set to be zero at the surface by initialization.
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, top_lev_), KOKKOS_LAMBDA(int kk) { 
      wsub_(icol,kk)  = wsubmin;
      wsubice_(icol,kk) = 0.001;
      wsig_(icol,kk)  = 0.001;
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, top_lev_, nlev_), KOKKOS_LAMBDA(int kk) { 
      wsub_(icol,kk)  = haero::sqrt(0.5*(tke_(icol,kk) + tke_(icol,kk+1))*(2.0/3.0));
      wsig_(icol,kk)  = mam4::utils::min_max_bound(0.001, 10.0, wsub_(icol,kk));
      wsubice_(icol,kk) = mam4::utils::min_max_bound(0.2, 10.0, wsub_(icol,kk));
      wsub_(icol,kk)  = haero::max(wsubmin, wsub_(icol,kk));
    });
  });

    //---------------------------------------------------------
    // Compute subgrid mean updraft velocity (w2)
    //---------------------------------------------------------
  
    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const haero::ThreadTeam &team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, nlev_), KOKKOS_LAMBDA(int kk) { 
        w2_(icol,kk) = subgrid_mean_updraft(w0_(icol,kk), wsig_(icol,kk));
      });
  });

    //-------------------------------------------------------------
    // Get number of activated aerosol for ice nucleation (naai)
    // from ice nucleation
    //-------------------------------------------------------------

    /*
    NOTE:"state_q" is a combination of subset of tracers added by "int_mmr_field_name" and "int_nmr_field_name".
    Only output we care about is "naai", "naai_hom" is never used anywhere

    Fortran code:
    call nucleate_ice_cam_calc(ncol, lchnk, temperature, state_q, pmid, &      ! input
         rho, wsubice, strat_cld_frac, dgnum, &                     ! input
         naai, naai_hom)                                 ! output
    */

    //-------------------------------------------------------------
    // Get old and new liquid cloud fractions when amount of cloud 
    // is above qsmall threshold value
    //-------------------------------------------------------------

   static constexpr auto qsmall = 1e-18; //cut-off for cloud amount (ice or liquid)
   
    /*
    MUST FIXME NOTE: We need old and new liquid cloud fractions here.
    We have the new liquid cloud fraction (liq_strat_cld_frac) but we need to 
    store the old (liq_strat_cld_frac_old) before we call SHOC. For now, we will make
    a note of it and use the new cloud fraction for the 
    old cloud fraction.

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
