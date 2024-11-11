#include <physics/mam/eamxx_mam_constituent_fluxes_functions.hpp>
#include <physics/mam/eamxx_mam_constituent_fluxes_interface.hpp>
namespace scream {

// ================================================================
//  Constructor
// ================================================================
MAMConstituentFluxes::MAMConstituentFluxes(const ekat::Comm &comm,
                                           const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMConstituentFluxes::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  using namespace ekat::units;
  auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  auto n_unit = 1 / kg;   // units of number mixing ratios of tracers
  auto nondim = ekat::units::Units::nondimensional();

  FieldLayout scalar2d     = grid_->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = grid_->get_3d_scalar_layout(false);

  static constexpr int pcnst = mam4::aero_model::pcnst;

  const FieldLayout scalar2d_pcnct =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constituents");

  // --------------------------------------------------------------------------
  // These variables are "Required" or pure inputs for the process
  // --------------------------------------------------------------------------
  // ----------- Atmospheric quantities -------------
  // Specific humidity [kg/kg](Require only for building DS)
  add_tracer<Required>("qv", grid_, q_unit);

  // Cloud liquid mass mixing ratio [kg/kg](Require only for building DS)
  add_tracer<Required>("qc", grid_, q_unit);

  // Cloud ice mass mixing ratio [kg/kg](Require only for building DS)
  add_tracer<Required>("qi", grid_, q_unit);

  // Cloud liquid number mixing ratio [1/kg](Require only for building DS)
  add_tracer<Required>("nc", grid_, n_unit);

  // Cloud ice number mixing ratio [1/kg](Require only for building DS)
  add_tracer<Required>("ni", grid_, n_unit);

  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);

  // Vertical pressure velocity [Pa/s] at midpoints (Require only for building
  // DS)
  add_field<Required>("omega", scalar3d_mid, Pa / s, grid_name);

  // Total pressure [Pa] at midpoints
  add_field<Required>("p_mid", scalar3d_mid, Pa, grid_name);

  // Total pressure [Pa] at interfaces
  add_field<Required>("p_int", scalar3d_int, Pa, grid_name);

  // Layer thickness(pdel) [Pa] at midpoints
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);

  // Planetary boundary layer height [m] (Require only for building DS)
  add_field<Required>("pbl_height", scalar2d, m, grid_name);

  // cloud fraction [nondimensional] computed by eamxx_cld_fraction_process
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);

  static constexpr Units m2(m * m, "m2");
  static constexpr Units s2(s * s, "s2");

  // Surface geopotential [m2/s2] (Require only for building DS)
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  // Constituent fluxes at the surface (gasses and aerosols)
  //[units: kg/m2/s (mass) or #/m2/s (number)]
  add_field<Required>("constituent_fluxes", scalar2d_pcnct, kg / m2 / s,
                      grid_name);

  // ---------------------------------------------------------------------
  // These variables are "Updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------
  // NOTE: Cloud borne aerosols are not updated in this process but are included
  // to create data structures.

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int mode = 0; mode < mam_coupling::num_aero_modes(); ++mode) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const std::string int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(mode);
    add_tracer<Updated>(int_nmr_field_name, grid_, n_unit);

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
    // NOT advected
    const std::string cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(mode);
    add_field<Updated>(cld_nmr_field_name, scalar3d_mid, n_unit, grid_name);

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(mode, a);
      if(not int_mmr_field_name.empty()) {
        add_tracer<Updated>(int_mmr_field_name, grid_, q_unit);
      }
      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
      // NOT advected
      const std::string cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(mode, a);
      if(not cld_mmr_field_name.empty()) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_mid, q_unit, grid_name);
      }
    }  // end for loop num species
  }    // end for loop for num modes

  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const std::string gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_tracer<Updated>(gas_mmr_field_name, grid_, q_unit);
  }  // end for loop num gases

}  // set_grid

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by the above
// Buffer type given the number of columns and vertical levels
size_t MAMConstituentFluxes::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializes the Buffer type with sufficient memory to store
// intermediate (dry) quantities on the given number of columns with the given
// number of vertical levels. Returns the number of bytes allocated.
void MAMConstituentFluxes::init_buffers(
    const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(
      used_mem == requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for MAMConstituentFluxes.");
}

// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMConstituentFluxes::initialize_impl(const RunType run_type) {
  // ---------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ---------------------------------------------------------------

  // Populate the wet atmosphere state with views from fields
  // FIMXE: specifically look which among these are actually used by the process
  wet_atm_.qv = get_field_in("qv").get_view<const Real **>();

  // Following wet_atm vars are required only for building DS
  wet_atm_.qc = get_field_in("qc").get_view<const Real **>();
  wet_atm_.nc = get_field_in("nc").get_view<const Real **>();
  wet_atm_.qi = get_field_in("qi").get_view<const Real **>();
  wet_atm_.ni = get_field_in("ni").get_view<const Real **>();

  // Populate the dry atmosphere state with views from fields
  dry_atm_.T_mid = get_field_in("T_mid").get_view<const Real **>();
  dry_atm_.p_mid = get_field_in("p_mid").get_view<const Real **>();
  dry_atm_.p_del = get_field_in("pseudo_density").get_view<const Real **>();
  dry_atm_.p_int = get_field_in("p_int").get_view<const Real **>();

  // Following dry_atm vars are required only for building DS
  dry_atm_.cldfrac = get_field_in("cldfrac_tot").get_view<const Real **>();
  dry_atm_.pblh    = get_field_in("pbl_height").get_view<const Real *>();
  dry_atm_.omega   = get_field_in("omega").get_view<const Real **>();
  dry_atm_.p_del   = get_field_in("pseudo_density").get_view<const Real **>();

  // store fields converted to dry mmr from wet mmr in dry_atm_
  dry_atm_.z_mid     = buffer_.z_mid;
  dry_atm_.z_iface   = buffer_.z_iface;
  dry_atm_.dz        = buffer_.dz;
  dry_atm_.qv        = buffer_.qv_dry;
  dry_atm_.qc        = buffer_.qc_dry;
  dry_atm_.nc        = buffer_.nc_dry;
  dry_atm_.qi        = buffer_.qi_dry;
  dry_atm_.ni        = buffer_.ni_dry;
  dry_atm_.w_updraft = buffer_.w_updraft;
  dry_atm_.z_surf    = 0.0;  // FIXME: for now

  // Constituent fluxes at the surface (gasses and aerosols) [kg/m2/s]
  constituent_fluxes_ =
      get_field_in("constituent_fluxes").get_view<const Real **>();

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const std::string int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(m);
    wet_aero_.int_aero_nmr[m] =
        get_field_out(int_nmr_field_name).get_view<Real **>();

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    const std::string cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(m);
    wet_aero_.cld_aero_nmr[m] =
        get_field_out(cld_nmr_field_name).get_view<Real **>();

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(not int_mmr_field_name.empty()) {
        wet_aero_.int_aero_mmr[m][a] =
            get_field_out(int_mmr_field_name).get_view<Real **>();
      }

      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const std::string cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if(not cld_mmr_field_name.empty()) {
        wet_aero_.cld_aero_mmr[m][a] =
            get_field_out(cld_mmr_field_name).get_view<Real **>();
      }
    }
  }
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const std::string gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] =
        get_field_out(gas_mmr_field_name).get_view<Real **>();
  }

}  // end initialize_impl()

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMConstituentFluxes::run_impl(const double dt) {
  // -------------------------------------------------------------------
  // (LONG) NOTE: The following code is an adaptation of cflx.F90 code in
  // E3SM. In EAMxx, all constituents are considered "wet" (or have wet
  // mixing ratios), we are *not* doing any wet to dry conversions in the
  // for this process. We are simply updating the MAM4xx tracers using the
  // "constituent fluxes".
  // We are converting wet atm to dry atm. Since we do not use or update
  // any of the water constituents (qc, qv, qi etc.), we should be okay
  // to do this conversion. We need to do this conversion as our function
  // are built following HAERO data structures.
  // -------------------------------------------------------------------

  // Compute vertical layer heights and updraft velocity. We need these to fully
  // populate dry_atm_, so that we can form a HAERO atmosphere object. HAERO
  // atmosphere object is used to for state%q like array.
  // NOTE: We cannot pass a member of the interface class (in this case, MAMConstituentFluxes)
  // inside a parallel_for. Instead, we must create a soft copy of each member.
  const auto & wet_atm = wet_atm_;
  const auto & dry_atm= dry_atm_;

  auto lambda =
      KOKKOS_LAMBDA(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) {
    const int icol = team.league_rank();       // column index
    compute_dry_mixing_ratios(team, wet_atm,  // in
                              dry_atm,        // out
                              icol);           // in
    team.team_barrier();
    // vertical heights has to be computed after computing dry mixing ratios
    // for atmosphere
    compute_vertical_layer_heights(team,        // in
                                   dry_atm,    // out
                                   icol);       // in
    compute_updraft_velocities(team, wet_atm,  // in
                               dry_atm,        // out
                               icol);           // in
  };
  // policy
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  Kokkos::parallel_for("mam_cfi_compute_updraft", scan_policy, lambda);

  update_gas_aerosols_using_constituents(ncol_, nlev_, dt, dry_atm_,
                                         constituent_fluxes_,
                                         // output
                                         wet_aero_);
}  // run_impl ends

// =============================================================================
}  // namespace scream
