#include <physics/mam/eamxx_mam_constituent_fluxes_functions.hpp>
#include <physics/mam/eamxx_mam_constituent_fluxes_interface.hpp>
#include <physics/mam/physical_limits.hpp>
#include <share/property_checks/field_within_interval_check.hpp>
namespace scream {

// ================================================================
//  Constructor
// ================================================================
MAMConstituentFluxes::MAMConstituentFluxes(const ekat::Comm &comm,
                                           const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);
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
  static constexpr int pcnst = mam4::aero_model::pcnst;

  const FieldLayout scalar2d_pcnct =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constituents");

  add_tracers_wet_atm();
  add_fields_dry_atm();
  static constexpr Units m2(m * m, "m2");
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
  // add tracers, e.g., num_a1, soa_a1
  add_tracers_interstitial_aerosol();
  // add tracer gases, e.g., O3
  add_tracers_gases();
  // add fields e.g., num_c1, soa_c1
  add_fields_cloudborne_aerosol();
}  // set_grid

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by the above
// Buffer type given the number of columns and vertical levels
size_t MAMConstituentFluxes::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_, 0, 0);
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
  // Check the interval values for the following fields used by this interface.
  // NOTE: We do not include aerosol and gas species, e.g., soa_a1, num_a1,
  // because we automatically added these fields.
  const std::map<std::string, std::pair<Real, Real>> ranges_cons_fluxes = {
      {"constituent_fluxes", {0, 1e10}}  // FIXME
  };
  set_ranges_process(ranges_cons_fluxes);
  add_interval_checks();

  // Populate the wet atmosphere state with views from fields
  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);

  // Constituent fluxes at the surface (gasses and aerosols) [kg/m2/s]
  constituent_fluxes_ =
      get_field_in("constituent_fluxes").get_view<const Real **>();

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
  // NOTE: We cannot pass a member of the interface class (in this case,
  // MAMConstituentFluxes) inside a parallel_for. Instead, we must create a soft
  // copy of each member.
  const auto &wet_atm = wet_atm_;
  const auto &dry_atm = dry_atm_;

  auto lambda =
      KOKKOS_LAMBDA(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) {
    const int icol = team.league_rank();      // column index
    compute_dry_mixing_ratios(team, wet_atm,  // in
                              dry_atm,        // out
                              icol);          // in
    team.team_barrier();
    // vertical heights has to be computed after computing dry mixing ratios
    // for atmosphere
    compute_vertical_layer_heights(team,       // in
                                   dry_atm,    // out
                                   icol);      // in
    compute_updraft_velocities(team, wet_atm,  // in
                               dry_atm,        // out
                               icol);          // in
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
