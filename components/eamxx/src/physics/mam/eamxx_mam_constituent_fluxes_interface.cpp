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

  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);

  static constexpr int pcnst = mam4::aero_model::pcnst;

  const FieldLayout scalar2d_pcnct =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constituents");

  // --------------------------------------------------------------------------
  // These variables are "Required" or pure inputs for the process
  // --------------------------------------------------------------------------
  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);

  // Layer thickness(pdel) [Pa] at midpoints
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);

  static constexpr Units m2(m * m, "m2");
  add_field<Required>("constituent_fluxes", scalar2d_pcnct, kg / m2 / s,
                      grid_name);

  // ---------------------------------------------------------------------
  // These variables are "Updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);

    add_field<Updated>(int_nmr_field_name, scalar3d_mid, n_unit, grid_name,
                       "tracers");
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_mid, q_unit, grid_name,
                           "tracers");
      }
    }
  }

  // aerosol-related gases: mass mixing ratios
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_mid, q_unit, grid_name,
                       "tracers");
  }

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
  // Inputs
  // ---------------------------------------------------------------
  constituent_fluxes_ =
      get_field_in("constituent_fluxes").get_view<const Real **>();
  dry_atm_.p_del = get_field_in("pseudo_density").get_view<const Real **>();

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    wet_aero_.int_aero_nmr[m] =
        get_field_out(int_nmr_field_name).get_view<Real **>();
    dry_aero_.int_aero_nmr[m] = buffer_.dry_int_aero_nmr[m];

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(strlen(int_mmr_field_name) > 0) {
        wet_aero_.int_aero_mmr[m][a] =
            get_field_out(int_mmr_field_name).get_view<Real **>();
        dry_aero_.int_aero_mmr[m][a] = buffer_.dry_int_aero_mmr[m][a];
      }
    }
  }
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] =
        get_field_out(gas_mmr_field_name).get_view<Real **>();
    dry_aero_.gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

  //-----------------------------------------------------------------
  // Allocate memory
  //-----------------------------------------------------------------

  rpdel_ = view_2d("rpdel_", ncol_, nlev_);
  //-----------------------------------------------------------------
  // Setup preprocessing and post processing
  //-----------------------------------------------------------------
  // preprocess_.initialize(constituent_fluxes_);

}  // end initialize_impl()

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMConstituentFluxes::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  // Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  // Kokkos::fence();
  using C                      = physics::Constants<Real>;
  static constexpr auto gravit = C::gravit;  // Gravity [m/s2]
  /*FORTRAN CODE:
    What we need:
    1. ncol
    2. all tracers
    3. timestep (verify what is rztodt??)
    4. gravit, rpdel
    5. cflx
    6. Ensure dry/et-wet/dry isconsistent....comment on that!!


    ncol = state%ncol

    !-------------------------------------------------------
    ! Assume 'wet' mixing ratios in surface diffusion code.
    ! don't convert co2 tracers to wet mixing ratios

    cnst_type_loc(:) = cnst_type(:)
    call set_dry_to_wet(state, cnst_type_loc)

    !-------------------------------------------------------
    ! Initialize ptend

    lq(:) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'clubb_srf', lq=lq)

    !-------------------------------------------------------
    ! Calculate tracer mixing ratio tendencies from cflx

    rztodt                 = 1._r8/ztodt
    ptend%q(:ncol,:pver,:) = state%q(:ncol,:pver,:)
    tmp1(:ncol)            = ztodt * gravit * state%rpdel(:ncol,pver)

    do m = 2, pcnst
      ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m) + tmp1(:ncol) *
    cam_in%cflx(:ncol,m) enddo

    ptend%q(:ncol,:pver,:) = (ptend%q(:ncol,:pver,:) - state%q(:ncol,:pver,:)) *
    rztodt

    ! Convert tendencies of dry constituents to dry basis.
    do m = 1,pcnst
       if (cnst_type(m).eq.'dry') then
          ptend%q(:ncol,:pver,m) =
    ptend%q(:ncol,:pver,m)*state%pdel(:ncol,:pver)/state%pdeldry(:ncol,:pver)
       endif
    end do

    !-------------------------------------------------------
    ! convert wet mmr back to dry before conservation check
    ! avoid converting co2 tracers again

    cnst_type_loc(:) = cnst_type(:)
    call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
    call set_wet_to_dry(state, cnst_type_loc)

  */
}

// =============================================================================
}  // namespace scream
