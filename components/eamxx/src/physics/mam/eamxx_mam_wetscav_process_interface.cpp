#include "physics/mam/eamxx_mam_wetscav_process_interface.hpp"

#include <ekat_team_policy_utils.hpp>

// =========================================================================================
// Helper function to initialize scratch arrays for convection processing
// =========================================================================================
namespace {

// Initialize scratch1Dviews array from work array for convective processing
// This must be called for each column to set up per-column scratch workspace
KOKKOS_INLINE_FUNCTION
void initialize_scratch1d_views(
    mam4::wetdep::ConvProc::Col1DView scratch1Dviews[mam4::wetdep::ConvProc::Col1DViewInd::NumScratch],
    Real* work_ptr,
    const int nlev)
{
  constexpr int pcnst = mam4::aero_model::pcnst;
  constexpr int pcnst_extd = mam4::wetdep::ConvProc::pcnst_extd;
  using ConvProc = mam4::wetdep::ConvProc;
  
  // Allocate each scratch array from the work array
  // Index 0: q - tracer mixing ratios
  scratch1Dviews[ConvProc::Col1DViewInd::q] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst);
  work_ptr += nlev * pcnst;
  
  // Indices 1-2: mu, md - updraft/downdraft mass fluxes
  scratch1Dviews[ConvProc::Col1DViewInd::mu] = 
      ConvProc::Col1DView(work_ptr, nlev + 1);
  work_ptr += nlev + 1;
  
  scratch1Dviews[ConvProc::Col1DViewInd::md] = 
      ConvProc::Col1DView(work_ptr, nlev + 1);
  work_ptr += nlev + 1;
  
  // Indices 3-6: eudp, dudp, eddp, dddp - entrainment/detrainment × dp
  scratch1Dviews[ConvProc::Col1DViewInd::eudp] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
  
  scratch1Dviews[ConvProc::Col1DViewInd::dudp] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
  
  scratch1Dviews[ConvProc::Col1DViewInd::eddp] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
  
  scratch1Dviews[ConvProc::Col1DViewInd::dddp] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
  
  // Index 7: rhoair - air density
  scratch1Dviews[ConvProc::Col1DViewInd::rhoair] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
  
  // Index 8: zmagl - height above ground level
  scratch1Dviews[ConvProc::Col1DViewInd::zmagl] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
  
  // Indices 9-12: gath, chat, conu, cond - gathered tracers and concentrations
  scratch1Dviews[ConvProc::Col1DViewInd::gath] = 
      ConvProc::Col1DView(work_ptr, (nlev + 1) * pcnst_extd);
  work_ptr += (nlev + 1) * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::chat] = 
      ConvProc::Col1DView(work_ptr, (nlev + 1) * pcnst_extd);
  work_ptr += (nlev + 1) * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::conu] = 
      ConvProc::Col1DView(work_ptr, (nlev + 1) * pcnst_extd);
  work_ptr += (nlev + 1) * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::cond] = 
      ConvProc::Col1DView(work_ptr, (nlev + 1) * pcnst_extd);
  work_ptr += (nlev + 1) * pcnst_extd;
  
  // Indices 13-14: dconudt_wetdep, dconudt_activa - updraft tendencies
  scratch1Dviews[ConvProc::Col1DViewInd::dconudt_wetdep] = 
      ConvProc::Col1DView(work_ptr, (nlev + 1) * pcnst_extd);
  work_ptr += (nlev + 1) * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::dconudt_activa] = 
      ConvProc::Col1DView(work_ptr, (nlev + 1) * pcnst_extd);
  work_ptr += (nlev + 1) * pcnst_extd;
  
  // Index 15: fa_u - updraft fractional area
  scratch1Dviews[ConvProc::Col1DViewInd::fa_u] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
  
  // Indices 16-20: dcondt arrays - downdraft tendencies
  scratch1Dviews[ConvProc::Col1DViewInd::dcondt] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst_extd);
  work_ptr += nlev * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::dcondt_wetdep] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst_extd);
  work_ptr += nlev * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::dcondt_prevap] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst_extd);
  work_ptr += nlev * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::dcondt_prevap_hist] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst_extd);
  work_ptr += nlev * pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::dcondt_resusp] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst_extd);
  work_ptr += nlev * pcnst_extd;
  
  // Indices 21-27: wd_flux and sum arrays
  scratch1Dviews[ConvProc::Col1DViewInd::wd_flux] = 
      ConvProc::Col1DView(work_ptr, pcnst_extd);
  work_ptr += pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::sumactiva] = 
      ConvProc::Col1DView(work_ptr, pcnst_extd);
  work_ptr += pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::sumaqchem] = 
      ConvProc::Col1DView(work_ptr, pcnst_extd);
  work_ptr += pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::sumprevap] = 
      ConvProc::Col1DView(work_ptr, pcnst_extd);
  work_ptr += pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::sumprevap_hist] = 
      ConvProc::Col1DView(work_ptr, pcnst_extd);
  work_ptr += pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::sumresusp] = 
      ConvProc::Col1DView(work_ptr, pcnst_extd);
  work_ptr += pcnst_extd;
  
  scratch1Dviews[ConvProc::Col1DViewInd::sumwetdep] = 
      ConvProc::Col1DView(work_ptr, pcnst_extd);
  work_ptr += pcnst_extd;
  
  // Indices 28-29: dqdt, qnew - tracer tendency and updated values
  scratch1Dviews[ConvProc::Col1DViewInd::dqdt] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst);
  work_ptr += nlev * pcnst;
  
  scratch1Dviews[ConvProc::Col1DViewInd::qnew] = 
      ConvProc::Col1DView(work_ptr, nlev * pcnst);
  work_ptr += nlev * pcnst;
  
  // Index 30: dlfdp - detrainment × dp
  scratch1Dviews[ConvProc::Col1DViewInd::dlfdp] = 
      ConvProc::Col1DView(work_ptr, nlev);
  work_ptr += nlev;
}

} // anonymous namespace

/*
-----------------------------------------------------------------
NOTES:
1. We should connect surface fluxes and add code to update the fluxes
2. Identify diagnostic variables and remove them from FM
3. Add assert statements to check output ranges
*/

namespace scream
{

// =========================================================================================
MAMWetscav::MAMWetscav(const ekat::Comm &comm, const ekat::ParameterList &params)
 : MAMGenericInterface(comm, params)
{
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
  check_fields_intervals_          = m_params.get<bool>("create_fields_interval_checks", false);
  scav_fraction_in_cloud_strat_    = m_params.get<Real>("scav_fraction_in_cloud_strat", 1.00);
  scav_fraction_in_cloud_conv_     = m_params.get<Real>("scav_fraction_in_cloud_conv", 0.00);
  scav_fraction_below_cloud_strat_ = m_params.get<Real>("scav_fraction_below_cloud_strat", 0.03);
  activation_fraction_in_cloud_conv_ =
      m_params.get<Real>("activation_fraction_in_cloud_conv", 0.40);
  convproc_do_aer_ = m_params.get<bool>("convproc_do_aer", false);
  convproc_do_gas_ = m_params.get<bool>("convproc_do_gas", false);
}

// ================================================================
//  SET_GRIDS
// ================================================================
void
MAMWetscav::create_requests()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  grid_                 = m_grids_manager->get_grid("physics");
  const auto &grid_name = grid_->name();

  ncol_                = grid_->get_num_local_dofs();      // Number of columns on this rank
  nlev_                = grid_->get_num_vertical_levels(); // Number of levels per column
  len_temporary_views_ = get_len_temporary_views();
  buffer_.set_len_temporary_views(len_temporary_views_);
  buffer_.set_num_scratch(num_2d_scratch_);

  const int nmodes    = mam4::AeroConfig::num_modes(); // Number of modes
  constexpr int pcnst = mam4::aero_model::pcnst;

  // layout for 3D (2d horiz X 1d vertical) variables at level
  // midpoints/interfaces
  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(LEV);
  FieldLayout scalar3d_int = grid_->get_3d_scalar_layout(ILEV);

  // layout for 2D (1d horiz X 1d vertical) variables
  FieldLayout scalar2d = grid_->get_2d_scalar_layout();

  // layout for 3D (ncol, nmodes, nlevs)
  FieldLayout scalar3d_mid_nmodes =
      grid_->get_3d_vector_layout(LEV, nmodes, mam_coupling::num_modes_tag_name());

  // layout for 2D (ncol, pcnst)
  FieldLayout scalar2d_pcnst = grid_->get_2d_vector_layout(pcnst, "num_phys_constituents");

  // --------------------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // --------------------------------------------------------------------------

  // ----------- Atmospheric quantities -------------
  // Specific humidity [kg/kg]
  add_tracers_wet_atm();
  add_fields_dry_atm();

  // cloud liquid number mixing ratio [1/kg]
  auto n_unit = 1 / kg; // units of number mixing ratios of tracers
  add_tracer<Required>("nc", grid_, n_unit);

  static constexpr auto m2 = m * m;
  static constexpr auto s2 = s * s;

  // Surface geopotential [m2/s2]
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);
  add_field<Required>("pseudo_density_dry", scalar3d_mid, Pa, grid_name);


  //----------- Variables from microphysics scheme -------------

  // Evaporation from stratiform rain [kg/kg/s]
  add_field<Required>("nevapr", scalar3d_mid, kg / kg / s, grid_name);

  // Stratiform rain production rate [kg/kg/s]
  add_field<Required>("precip_total_tend", scalar3d_mid, kg / kg / s, grid_name);

  //----------- Variables from macrophysics scheme -------------

  // Total cloud fraction [fraction]
  add_field<Required>("cldfrac_liq", scalar3d_mid, none, grid_name);

  //----------- Variables from ZM deep convection scheme -------------
  // Deep convection cloud water detrainment [kg/kg/s]
  add_field<Required>("zm_detr_qc", scalar3d_mid, kg / kg / s, grid_name);
  add_field<Required>("zm_detr_qi", scalar3d_mid, kg / kg / s, grid_name);

  // Cloud top and base indices from ZM deep convection
  add_field<Required>("zm_jt", scalar2d, none, grid_name);      // Cloud top level index
  add_field<Required>("zm_jcbot", scalar2d, none, grid_name);   // Cloud base level index

  //----------- Variables from ZM convection scheme -------------
  // Convective mass fluxes and entrainment/detrainment rates
  add_field<Required>("zm_mflx_up", scalar3d_mid, kg / m2 / s, grid_name);  // Updraft mass flux
  add_field<Required>("zm_mflx_dn", scalar3d_mid, kg / m2 / s, grid_name);  // Downdraft mass flux
  add_field<Required>("zm_entr_up", scalar3d_mid, 1 / s, grid_name);         // Updraft entrainment rate
  add_field<Required>("zm_detr_up", scalar3d_mid, 1 / s, grid_name);         // Updraft detrainment rate
  add_field<Required>("zm_entr_dn", scalar3d_mid, 1 / s, grid_name);         // Downdraft entrainment rate

  // ---------------------------------------------------------------------
  // These variables are "updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------
  // FIXME: we have not added code to update the surface fluxes.
  // -- surface fluxes (input/outpts) for the coupler's cam_out data struture
  // for the land model

  // Wet deposition of hydrophilic black carbon [kg/m2/s]
  add_field<Updated>("wetdep_hydrophilic_bc", scalar3d_mid, kg / m2 / s, grid_name);

  // Dry deposition of hydrophilic black carbon [kg/m2/s]
  add_field<Updated>("drydep_hydrophilic_bc", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of hydrophilic organic carbon [kg/m2/s]
  add_field<Updated>("wetdep_hydrophilic_oc", scalar3d_mid, kg / m2 / s, grid_name);

  // Dry deposition of hydrophilic organic carbon [kg/m2/s]
  add_field<Updated>("drydep_hydrophilic_oc", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of dust (bin1) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin1", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of dust (bin2) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin2", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of dust (bin3) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin3", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of dust (bin4) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin4", scalar3d_mid, kg / m2 / s, grid_name);

  // Interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios

  // NOTE: Interstitial aerosols are updated in the interface using the
  // "tendencies" from the wetscavenging process
  // add tracers, e.g., num_a1, soa_a1
  add_tracers_interstitial_aerosol();
  // add tracer gases, e.g., O3
  add_tracers_gases();
  // add fields e.g., num_c1, soa_c1
  add_fields_cloudborne_aerosol();

  // -------------------------------------------------------------
  // These variables are "Computed" or outputs for the process
  // -------------------------------------------------------------
  static constexpr auto m3 = m * m * m;

  // Aerosol dry particle diameter [m]
  add_field<Computed>("dgnum", scalar3d_mid_nmodes, m, grid_name);

  // Wet aerosol density [kg/m3]
  add_field<Computed>("wetdens", scalar3d_mid_nmodes, kg / m3, grid_name);

  // Aerosol water [kg/kg]
  add_field<Computed>("qaerwat", scalar3d_mid_nmodes, kg / kg, grid_name);

  // Wet aerosol diameter [m]
  add_field<Computed>("dgnumwet", scalar3d_mid_nmodes, m, grid_name);

  // Fraction of transported species that are insoluble [fraction]
  add_field<Computed>("fracis", scalar3d_mid, none, grid_name);

  // Aerosol wet deposition (interstitial) [kg/m2/s]
  add_field<Computed>("aerdepwetis", scalar2d_pcnst, kg / m2 / s, grid_name);

  // Aerosol wet deposition (cloud water) [kg/m2/s]
  add_field<Computed>("aerdepwetcw", scalar2d_pcnst, kg / m2 / s, grid_name);
}

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to store
// intermediate (dry) quantities on the given number of columns with the given
// number of vertical levels. Returns the number of bytes allocated.
void
MAMWetscav::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");
  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  std::cout << "used_mem " << used_mem << "requested_buffer_size_in_bytes() "
            << requested_buffer_size_in_bytes() << "\n";
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMWetscav.");
}

int MAMWetscav::get_len_temporary_views() {
  const int work_len = mam4::wetdep::get_aero_model_wetdep_work_len() * ncol_;
  return work_len;
}

void MAMWetscav::init_temporary_views() {
  auto work_ptr = (Real *)buffer_.temporary_views.data();
  // Allocate work array
  const int work_len = mam4::wetdep::get_aero_model_wetdep_work_len();
  work_              = view_2d(work_ptr, ncol_, work_len);
  work_ptr += ncol_ * work_len;

  /// error check
  // NOTE: workspace_provided can be larger than workspace_used, but let's try
  // to use the minimum amount of memory
  const int workspace_used     = work_ptr - buffer_.temporary_views.data();
  const int workspace_provided = buffer_.temporary_views.extent(0);
  EKAT_REQUIRE_MSG(workspace_used == workspace_provided,
                   "Error: workspace_used (" + std::to_string(workspace_used) +
                       ") and workspace_provided (" +
                       std::to_string(workspace_provided) +
                       ") should be equal. \n");
}
// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMWetscav::initialize_impl(const RunType run_type) {
  // Check the interval values for the following fields used by this interface.
  // NOTE: We do not include aerosol and gas species, e.g., soa_a1, num_a1,
  // because we automatically added these fields.
  const std::map<std::string, std::pair<Real, Real>> ranges_wetscav = {
      // wetscav
      {"drydep_hydrophilic_bc", {-1e10, 1e10}},  // FIXME
      {"drydep_hydrophilic_oc", {-1e10, 1e10}},  // FIXME
      {"wetdep_dust_bin1", {-1e10, 1e10}},       // FIXME
      {"wetdep_dust_bin2", {-1e10, 1e10}},       // FIXME
      {"wetdep_dust_bin3", {-1e10, 1e10}},       // FIXME
      {"wetdep_dust_bin4", {-1e10, 1e10}},       // FIXME
      {"wetdep_hydrophilic_bc", {-1e10, 1e10}},  // FIXME
      {"wetdep_hydrophilic_oc", {-1e10, 1e10}},  // FIXME
      {"precip_total_tend", {-1e10, 1e10}},      // FIXME
      {"aerdepwetcw", {-1e10, 1e10}},            // FIXME
      {"aerdepwetis", {-1e10, 1e10}},            // FIXME
      {"fracis", {-1e10, 1e10}},                 // FIXME
      {"qaerwat", {-1e10, 1e10}}                 // FIXME
  };

  set_ranges_process(ranges_wetscav);
  add_interval_checks();
  // ---------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ---------------------------------------------------------------

  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);
  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  // It populates wet_aero struct (wet_aero_) with:
  // interstitial aerosol, e.g., soa_a_1
  populate_interstitial_wet_aero(wet_aero_);
  // gases, e.g., O3
  populate_gases_wet_aero(wet_aero_);
  // cloudborne aerosol, e.g., soa_c_1
  populate_cloudborne_wet_aero(wet_aero_);
  // It populates dry_aero struct (dry_aero_) with:
  // interstitial aerosol, e.g., soa_a_1
  populate_interstitial_dry_aero(dry_aero_, buffer_);
  // gases, e.g., O3
  populate_gases_dry_aero(dry_aero_, buffer_);
  // cloudborne aerosol, e.g., soa_c_1
  populate_cloudborne_dry_aero(dry_aero_, buffer_);

  dry_atm_.phis = get_field_in("phis").get_view<const Real *>();
  // NOTE: In populate_dry_atm we use cldfrac_tot
  dry_atm_.cldfrac = get_field_in("cldfrac_liq").get_view<const Real **>();

  //---------------------------------------------------------------------------------
  // Allocate memory
  //---------------------------------------------------------------------------------
  // Alllocate aerosol-related gas tendencies
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    set_field_w_scratch_buffer(dry_aero_tends_.gas_mmr[g], buffer_, true);
  }

  // Allocate aerosol state tendencies (interstitial aerosols only)
  for(int imode = 0; imode < mam_coupling::num_aero_modes(); ++imode) {
    set_field_w_scratch_buffer(dry_aero_tends_.int_aero_nmr[imode], buffer_,
                               true);

    for(int ispec = 0; ispec < mam_coupling::num_aero_species(); ++ispec) {
      set_field_w_scratch_buffer(dry_aero_tends_.int_aero_mmr[imode][ispec],
                                 buffer_, true);
    }
  }
  // Allocate work array
  init_temporary_views();
  isprx_ = int_view_2d("isprx", ncol_, nlev_);
  // TODO: Following variables are from convective parameterization (not
  // implemented yet in EAMxx), so should be zero for now
  // NOTE:If we use buffer_ to set the following inputs,
  // we must set these views to zero at every time step.
  sh_frac_ = view_2d("sh_frac_", ncol_, nlev_);

  // Deep convective cloud fraction [fraction]
  dp_frac_ = view_2d("dp_frac_", ncol_, nlev_);

  // Evaporation rate of shallow convective precipitation >=0. [kg/kg/s]
  evapcsh_ = view_2d("evapcsh_", ncol_, nlev_);

  // Evaporation rate of deep convective precipitation >=0. [kg/kg/s]
  evapcdp_ = view_2d("evapcdp_", ncol_, nlev_);

  // Rain production, shallow convection [kg/kg/s]
  rprdsh_ = view_2d("rprdsh_", ncol_, nlev_);

  // Rain production, deep convection [kg/kg/s]
  rprddp_ = view_2d("rprddp_", ncol_, nlev_);

  // In cloud water mixing ratio, deep convection
  icwmrdp_ = view_2d("icwmrdp_", ncol_, nlev_);

  // In cloud water mixing ratio, shallow convection
  icwmrsh_ = view_2d("icwmrsh_", ncol_, nlev_);
  
  // Detraining cloud water from deep convection [kg/kg/s]
  // This will be computed from ZM scheme outputs: dlf = zm_detr_qc + zm_detr_qi
  dlf_ = view_2d("dlf_", ncol_, nlev_);

  // Shallow convection detrainment [kg/kg/s]
  // NOTE: EAMxx does not yet have a shallow convection scheme implemented.
  // This is initialized to zero and will be updated when shallow convection is added.
  dlfsh_ = view_2d("dlfsh_", ncol_, nlev_);
  
  // Shallow convection entrainment/(entrainment+detrainment) ratio
  // NOTE: Also not available without a shallow convection scheme.
  // Initialized to zero for now.
  sh_e_ed_ratio_ = view_2d("sh_e_ed_ratio_", ncol_, nlev_);

  // Temporary view for pressure thickness in mb (converted from Pa)
  dp_tmp_ = view_2d("dp_tmp", ncol_, nlev_);
  
  calcsize_data_.initialize();
  // wetscav uses update_mmr=true;
  calcsize_data_.set_update_mmr(true);

  view_2d_host scavimptblvol_host("scavimptblvol_host",
                                  mam4::aero_model::nimptblgrow_total,
                                  mam4::AeroConfig::num_modes());
  view_2d_host scavimptblnum_host("scavimptblnum_host",
                                  mam4::aero_model::nimptblgrow_total,
                                  mam4::AeroConfig::num_modes());

  mam4::wetdep::init_scavimptbl(aero_config_, scavimptblvol_host, scavimptblnum_host);

  scavimptblnum_ = view_2d("scavimptblnum", mam4::aero_model::nimptblgrow_total,
                           mam4::AeroConfig::num_modes());
  scavimptblvol_ = view_2d("scavimptblvol", mam4::aero_model::nimptblgrow_total,
                           mam4::AeroConfig::num_modes());
  Kokkos::deep_copy(scavimptblnum_, scavimptblnum_host);
  Kokkos::deep_copy(scavimptblvol_, scavimptblvol_host);

  // Initialize ConvProc configuration for species classification and resuspension mapping
  // The default constructor initializes species_class and mmtoo_prevap_resusp arrays
  // with reference values from mam4xx/convproc.hpp
  convproc_config_ = mam4::ConvProc::Config();

  // Initialize ptend_lq flags based on species to process
  // This determines which species undergo convective processing based on:
  //   - species_class (aerosol vs gas vs other) from convproc_config_
  //   - convproc_do_aer and convproc_do_gas runtime flags from namelist
  // NOTE: This assumes convproc_do_aer_ and convproc_do_gas_ don't change during simulation
  for (int i = 0; i < mam4::aero_model::pcnst; ++i) {
    ptend_lq_[i] = (convproc_config_.species_class[i] == mam4::ConvProc::species_class::aerosol && 
                    convproc_do_aer_) ||
                   (convproc_config_.species_class[i] == mam4::ConvProc::species_class::gas && 
                    convproc_do_gas_);
  }
}

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMWetscav::run_impl(const double dt) {
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of all variables
  // needed by this process or setting up MAM4xx classes and their objects
  pre_process(wet_aero_, dry_aero_, wet_atm_, dry_atm_);
  Kokkos::fence();

  const mam_coupling::DryAtmosphere &dry_atm = dry_atm_;
  const auto &dry_aero                       = dry_aero_;
  const auto &work                           = work_;
  const auto &isprx                          = isprx_;
  const auto &dry_aero_tends                 = dry_aero_tends_;

  // ---------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // ---------------------------------------------------------------

  //----------- Variables from convective scheme -------------

  // TODO: Following variables are from convective parameterization (not
  // implemented yet in EAMxx), so should be zero for now

  auto sh_frac = sh_frac_;

  // Deep convective cloud fraction [fraction]
  auto dp_frac = dp_frac_;

  // Evaporation rate of shallow convective precipitation >=0. [kg/kg/s]
  auto evapcsh = evapcsh_;

  // Evaporation rate of deep convective precipitation >=0. [kg/kg/s]
  auto evapcdp = evapcdp_;

  // Rain production, shallow convection [kg/kg/s]
  auto rprdsh = rprdsh_;

  // Rain production, deep convection [kg/kg/s]
  auto rprddp = rprddp_;

  // In cloud water mixing ratio, deep convection
  auto icwmrdp = icwmrdp_;

  // In cloud water mixing ratio, shallow convection
  auto icwmrsh = icwmrsh_;

  // Shallow convection detrainment (initialized to zero in initialize_impl)
  auto dlfsh = dlfsh_;

  // Shallow convection entrainment/(entrainment+detrainment) ratio (initialized to zero)
  auto sh_e_ed_ratio = sh_e_ed_ratio_;

  //----------- Variables from ZM deep convection scheme -------------
  // Deep convection cloud water detrainment [kg/kg/s]
  // zm_detr_qc and zm_detr_qi are the liquid and ice detrainment tendencies from ZM
  // dlf = total cloud water detrainment (liquid + ice)
  auto zm_detr_qc_in = get_field_in("zm_detr_qc").get_view<const Real **>();
  auto zm_detr_qi_in = get_field_in("zm_detr_qi").get_view<const Real **>();
  
  // Compute total detrainment (dlf = zm_detr_qc + zm_detr_qi)
  Kokkos::parallel_for("compute_dlf", 
    policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
        dlf_(icol, kk) = zm_detr_qc_in(icol, kk) + zm_detr_qi_in(icol, kk);
      });
    });
  Kokkos::fence();

  // Detraining cld H20 from deep convection [kg/kg/s]
  auto dlf = dlf_;

  //----------- Variables from macrophysics scheme -------------
  // Total cloud fraction
  auto cldt = get_field_in("cldfrac_liq").get_view<const Real **>();

  //----------- Variables from microphysics scheme -------------

  // Evaporation from stratiform rain [kg/kg/s]
  auto nevapr = get_field_in("nevapr").get_view<const Real **>();

  // Stratiform rain production rate [kg/kg/s]
  auto prain = get_field_in("precip_total_tend").get_view<const Real **>();
  
  // Convection mass fluxes and rates from ZM scheme
  auto zm_mflx_up = get_field_in("zm_mflx_up").get_view<const Real **>();
  auto zm_mflx_dn = get_field_in("zm_mflx_dn").get_view<const Real **>();
  auto zm_entr_up = get_field_in("zm_entr_up").get_view<const Real **>();
  auto zm_detr_up = get_field_in("zm_detr_up").get_view<const Real **>();
  auto zm_entr_dn = get_field_in("zm_entr_dn").get_view<const Real **>();
  auto pseudo_density_dry = get_field_in("pseudo_density_dry").get_view<const Real **>();
  
  // ------------------------------------------------------------------
  // These variables are "Computed" or pure outputs for the process
  // ------------------------------------------------------------------

  const auto aerdepwetis = get_field_out("aerdepwetis").get_view<Real **>();
  const auto aerdepwetcw = get_field_out("aerdepwetcw").get_view<Real **>();

  const auto wet_geometric_mean_diameter_i =
      get_field_out("dgnumwet").get_view<Real ***>();
  const auto dry_geometric_mean_diameter_i =
      get_field_out("dgnum").get_view<Real ***>();
  const auto qaerwat = get_field_out("qaerwat").get_view<Real ***>();
  const auto wetdens = get_field_out("wetdens").get_view<Real ***>();

  const auto policy = TPF::get_default_team_policy(ncol_, nlev_);

  // Making a local copy of 'nlev_' because we cannot use a member of a class
  // inside a parallel_for.
  const int nlev = nlev_;

  // Zero out tendencies otherwise, they are initialized to junk values
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    Kokkos::deep_copy(dry_aero_tends_.int_aero_nmr[m], 0);
    for(int a = 0; a < mam4::num_species_mode(m); ++a) {
      Kokkos::deep_copy(dry_aero_tends_.int_aero_mmr[m][a], 0);
    }
  }

  const auto &calsize_data  = calcsize_data_;
  const auto &scavimptblnum = scavimptblnum_;
  const auto &scavimptblvol = scavimptblvol_;
  const Real scav_fraction_in_cloud_strat  = scav_fraction_in_cloud_strat_;   
  const Real scav_fraction_in_cloud_conv = scav_fraction_in_cloud_conv_;
  const Real scav_fraction_below_cloud_strat  = scav_fraction_below_cloud_strat_;
  const Real activation_fraction_in_cloud_conv = activation_fraction_in_cloud_conv_;
  const auto &aero_config = aero_config_;
  const auto &dp_tmp = dp_tmp_;
 
  // Get cloud top and base indices from ZM convection scheme
  auto zm_jt_in    = get_field_in("zm_jt").get_view<const Real *>();
  auto zm_jcbot_in = get_field_in("zm_jcbot").get_view<const Real *>();
 
  // Get convective processing flags from namelist parameters
  const bool convproc_do_aer = convproc_do_aer_;
  const bool convproc_do_gas = convproc_do_gas_;
  
  // Species classification and resuspension mapping arrays
  // Get pointers to species classification and resuspension mapping arrays
  // These were initialized in initialize_impl from mam4xx ConvProc defaults
  const int* species_class = convproc_config_.species_class;
  const int* mmtoo_prevap_resusp = convproc_config_.mmtoo_prevap_resusp;

  // Get pointer to species processing flags
  // These were initialized in initialize_impl based on species_class and runtime flags
  const bool* ptend_lq = ptend_lq_;
  
  // Loop over atmosphere columns
  Kokkos::parallel_for("MAMWetscav::run_impl::aero_model_wetdep",
      policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol = team.league_rank();  // column index

        const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
        // set surface state data
        // fetch column-specific subviews into aerosol prognostics
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);
        // fetch column-specific subviews into aerosol tendencies
        // Note: we are only updating interstitial aerosols.
        mam4::Tendencies tends =
            mam_coupling::interstitial_aerosols_tendencies_for_column(
                dry_aero_tends, icol);
        /// shallow_convective_precipitation_production
        const auto rprdsh_icol = ekat::subview(rprdsh, icol);
        // deep_convective_precipitation_production
        const auto rprddp_icol = ekat::subview(rprddp, icol);
        // deep_convective_precipitation_evaporation
        const auto evapcdp_icol = ekat::subview(evapcdp, icol);
        // shallow_convective_precipitation_evaporation =
        const auto evapcsh_icol = ekat::subview(evapcsh, icol);
        // deep_convective_cloud_fraction
        const auto dp_frac_icol = ekat::subview(dp_frac, icol);
        // shallow_convective_cloud_fraction    =
        const auto sh_frac_icol = ekat::subview(sh_frac, icol);

        const auto icwmrdp_col  = ekat::subview(icwmrdp, icol);
        const auto icwmrsh_icol = ekat::subview(icwmrsh, icol);
        const auto nevapr_icol  = ekat::subview(nevapr, icol);
        const auto cldt_icol    = ekat::subview(cldt, icol);

        const auto dlf_icol   = ekat::subview(dlf, icol);
        auto aerdepwetis_icol = ekat::subview(aerdepwetis, icol);
        auto aerdepwetcw_icol = ekat::subview(aerdepwetcw, icol);
        auto work_icol        = ekat::subview(work, icol);
        auto wet_diameter_icol =
            ekat::subview(wet_geometric_mean_diameter_i, icol);
        auto dry_diameter_icol =
            ekat::subview(dry_geometric_mean_diameter_i, icol);
        auto qaerwat_icol     = ekat::subview(qaerwat, icol);
        auto wetdens_icol     = ekat::subview(wetdens, icol);
        const auto prain_icol = ekat::subview(prain, icol);

        auto isprx_icol = ekat::subview(isprx, icol);
        
        // Scratch arrays for convection processing (per column)
        // These are required by ma_convproc_intr when convproc_do_aer or convproc_do_gas is true
        mam4::wetdep::ConvProc::Col1DView scratch1Dviews[mam4::wetdep::ConvProc::Col1DViewInd::NumScratch];
        
        // Initialize scratch arrays from work array for convective processing
        initialize_scratch1d_views(scratch1Dviews, work_icol.data(), nlev);
        
        // Convection mass flux fields from ZM scheme
        // Extract convection fields for this column
        const auto mu_icol = ekat::subview(zm_mflx_up, icol);  // Updraft mass flux
        const auto md_icol = ekat::subview(zm_mflx_dn, icol);  // Downdraft mass flux
        const auto eu_icol = ekat::subview(zm_entr_up, icol);  // Updraft entrainment rate
        const auto du_icol = ekat::subview(zm_detr_up, icol);  // Updraft detrainment rate
        const auto ed_icol = ekat::subview(zm_entr_dn, icol);  // Downdraft entrainment rate
        
        // Get pressure thickness from atmosphere and convert from Pa to mb
        // dp (total) is from dry_atm.p_del, dpdry (dry) is from pseudo_density_dry
        // Both need to be converted from Pa to mb: 1 mb = 100 Pa, so multiply by 0.01
        constexpr Real pa_to_mb = 0.01;
        
        // Get the p_del column for this icol (total pressure thickness)
        const auto p_del_icol = ekat::subview(dry_atm.p_del, icol);
        
        // Get the pseudo_density_dry column for this icol (dry pressure thickness)
        const auto p_del_dry_icol = ekat::subview(pseudo_density_dry, icol);
        
        
        // Get temporary dp view for this column
        auto dp_tmp_icol = ekat::subview(dp_tmp, icol);
        // Convert pressure thickness from Pa to mb
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, mam4::nlev), [&](int kk) {
          dp_tmp_icol(kk) = p_del_icol(kk) * pa_to_mb;
        });
        team.team_barrier();
        
        const auto dp_icol = dp_tmp_icol;
        // Extract shallow convection fields for this column
        // NOTE: These are initialized to zero in initialize_impl since shallow convection
        // is not yet implemented in EAMxx.
        const auto dlfsh_icol = ekat::subview(dlfsh, icol);
        const auto sh_e_ed_ratio_icol = ekat::subview(sh_e_ed_ratio, icol);

        // Get cloud top and base indices from ZM convection for this column
        const int ktop = static_cast<int>(zm_jt_in(icol));
        const int kbot = static_cast<int>(zm_jcbot_in(icol));

        mam4::wetdep::aero_model_wetdep(
            team, atm, progs, tends, dt,
            scav_fraction_in_cloud_strat,
            scav_fraction_in_cloud_conv,
            scav_fraction_below_cloud_strat,
            activation_fraction_in_cloud_conv,
            // inputs
            cldt_icol, rprdsh_icol, rprddp_icol, evapcdp_icol, evapcsh_icol,
            dp_frac_icol, sh_frac_icol, icwmrdp_col, icwmrsh_icol, nevapr_icol,
            dlf_icol, prain_icol, scavimptblnum, scavimptblvol, calsize_data,
            // in/out calcsize and water_uptake
            wet_diameter_icol, dry_diameter_icol, qaerwat_icol, wetdens_icol,
            // output
            aerdepwetis_icol, aerdepwetcw_icol, work_icol, isprx_icol,
           // Convection mass flux parameters
           scratch1Dviews,
           mu_icol, md_icol, du_icol, eu_icol, ed_icol,
           dp_icol, p_del_dry_icol, dlfsh_icol, sh_e_ed_ratio_icol,
           ktop, kbot,
           convproc_do_aer, convproc_do_gas,
           species_class, mmtoo_prevap_resusp,
           ptend_lq,
           aero_config);
        team.team_barrier();
        // update interstitial aerosol state
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
          for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
            const auto n_mode_i       = progs.n_mode_i[m];
            const auto tends_n_mode_i = tends.n_mode_i[m];
            n_mode_i(kk) += tends_n_mode_i(kk) * dt;
            for(int a = 0; a < mam4::num_species_mode(m); ++a) {
              const auto q_aero_i       = progs.q_aero_i[m][a];
              const auto tends_q_aero_i = tends.q_aero_i[m][a];
              q_aero_i(kk) += tends_q_aero_i(kk) * dt;
            }
          }
        });  // parallel_for for update interstitial aerosol state
      });    // icol parallel_for loop

  // call post processing to convert dry mixing ratios to wet mixing ratios
  // and update the state
  Kokkos::fence();
  post_process(wet_aero_, dry_aero_, dry_atm_);
  Kokkos::fence();  // wait before returning to calling function
}

// =========================================================================================
}  // namespace scream
