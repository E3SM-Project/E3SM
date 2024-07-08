#include "physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp"

/*
Future work:
Wirte comments
write in/outs for all variables clearly
*/

namespace scream {

// =========================================================================================
MAMSrfOnlineEmiss::MAMSrfOnlineEmiss(const ekat::Comm &comm,
                                     const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMSrfOnlineEmiss::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  // set grid for all the inputs and outputs
  // use physics grid
  grid_ = grids_manager->get_grid("Physics");

  // Name of the grid
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variables
  // mid points
  FieldLayout scalar3d_layout_mid{{COL, LEV}, {ncol_, nlev_}};
  // interfaces
  FieldLayout scalar3d_layout_int{{COL, ILEV}, {ncol_, nlev_ + 1}};

  // layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{{COL}, {ncol_}};

  using namespace ekat::units;
  auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  auto n_unit = 1 / kg;   // units of number mixing ratios of tracers

  auto nondim = ekat::units::Units::nondimensional();

  // -------------------------------------------------------------------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // -------------------------------------------------------------------------------------------------------------------------
  // atmospheric quantities
  // specific humidity [kg/kg]
  add_field<Required>("qv", scalar3d_layout_mid, q_unit, grid_name, "tracers");

  // cloud liquid mass mixing ratio [kg/kg]
  add_field<Required>("qc", scalar3d_layout_mid, q_unit, grid_name, "tracers");

  // cloud ice mass mixing ratio [kg/kg]
  add_field<Required>("qi", scalar3d_layout_mid, q_unit, grid_name, "tracers");

  // cloud liquid number mixing ratio [1/kg]
  add_field<Required>("nc", scalar3d_layout_mid, n_unit, grid_name, "tracers");

  // cloud ice number mixing ratio [1/kg]
  add_field<Required>("ni", scalar3d_layout_mid, n_unit, grid_name, "tracers");

  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name);

  // Vertical pressure velocity [Pa/s] at midpoints
  add_field<Required>("omega", scalar3d_layout_mid, Pa / s, grid_name);

  // Total pressure [Pa] at midpoints
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name);

  // Total pressure [Pa] at interfaces
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid_name);

  // Layer thickness(pdel) [Pa] at midpoints
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name);

  // planetary boundary layer height
  add_field<Required>("pbl_height", scalar2d_layout_col, m, grid_name);

  // ========================================================================
  // Output from this whole process
  // ========================================================================

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int mode = 0; mode < mam_coupling::num_aero_modes(); ++mode) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char *int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(mode);
    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name, "tracers");

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
    // NOT advected
    const char *cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(mode);
    add_field<Updated>(cld_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name);

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(mode, a);
      if(strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, q_unit,
                           grid_name, "tracers");
      }
      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
      // NOT advected
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(mode, a);
      if(strlen(cld_mmr_field_name) > 0) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_layout_mid, q_unit,
                           grid_name);
      }
    }  // end for loop num species
  }    // end for loop for num modes

  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_layout_mid, q_unit,
                       grid_name, "tracers");
  }  // end for loop num gases

}  // function set_grids ends

// ================================================================
//  INIT_BUFFERS
// ================================================================
size_t MAMSrfOnlineEmiss::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

void MAMSrfOnlineEmiss::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(
      used_mem == requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for MAMSrfOnlineEmiss.");
}

// ================================================================
//  INITIALIZE_IMPL
// ================================================================

void MAMSrfOnlineEmiss::initialize_impl(const RunType run_type) {
  // Gather runtime options
  //(e.g.) runtime_options.lambda_low    = m_params.get<double>("lambda_low");

  wet_atm_.qv = get_field_in("qv").get_view<const Real **>();
  wet_atm_.qc = get_field_in("qc").get_view<const Real **>();
  wet_atm_.nc = get_field_in("nc").get_view<const Real **>();
  wet_atm_.qi = get_field_in("qi").get_view<const Real **>();
  wet_atm_.ni = get_field_in("ni").get_view<const Real **>();

  dry_atm_.T_mid = get_field_in("T_mid").get_view<const Real **>();
  dry_atm_.p_mid = get_field_in("p_mid").get_view<const Real **>();
  dry_atm_.p_del = get_field_in("pseudo_density").get_view<const Real **>();
  dry_atm_.omega = get_field_in("omega").get_view<const Real **>();

  // store fields converted to dry mmr from wet mmr in dry_atm_
  dry_atm_.qv = buffer_.qv_dry;
  dry_atm_.qc = buffer_.qc_dry;
  dry_atm_.nc = buffer_.nc_dry;
  dry_atm_.qi = buffer_.qi_dry;
  dry_atm_.ni = buffer_.ni_dry;

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    wet_aero_.int_aero_nmr[m] =
        get_field_out(int_nmr_field_name).get_view<Real **>();
    dry_aero_.int_aero_nmr[m] = buffer_.dry_int_aero_nmr[m];

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    const char *cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);
    wet_aero_.cld_aero_nmr[m] =
        get_field_out(cld_nmr_field_name).get_view<Real **>();
    dry_aero_.cld_aero_nmr[m] = buffer_.dry_cld_aero_nmr[m];

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);
      if(strlen(int_mmr_field_name) > 0) {
        wet_aero_.int_aero_mmr[m][a] =
            get_field_out(int_mmr_field_name).get_view<Real **>();
        dry_aero_.int_aero_mmr[m][a] = buffer_.dry_int_aero_mmr[m][a];
      }

      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if(strlen(cld_mmr_field_name) > 0) {
        wet_aero_.cld_aero_mmr[m][a] =
            get_field_out(cld_mmr_field_name).get_view<Real **>();
        dry_aero_.cld_aero_mmr[m][a] = buffer_.dry_cld_aero_mmr[m][a];
      }
    }
  }
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] =
        get_field_out(gas_mmr_field_name).get_view<Real **>();
    dry_aero_.gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

  // set up our preprocess functor
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                         dry_aero_);
}

// =========================================================================================
void MAMSrfOnlineEmiss::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  /* Rough notes:

  Here we should implement or port the chem_emissions subroutine in
  chemistry.F90. Basically call two subroutines, aero_model_emissions and
  set_srf_emissions.

  Here is the code:

    ! initialize chemistry constituent surface fluxes to zero
    do m = 2,pcnst
       n = map2chm(m)
       if (n>0) cam_in%cflx(:,m) = 0._r8
    enddo

    ! aerosol emissions ...
    call aero_model_emissions( state, & ! in
                               cam_in ) ! out

    ! prescribed emissions from file ...

    !-----------------------------------------------------------------------
    !        ... Set surface emissions
    !-----------------------------------------------------------------------
    call set_srf_emissions( lchnk, ncol, sflx(:,:) )

    do m = 1,pcnst
       n = map2chm(m)
       if ( n /= h2o_ndx .and. n > 0 ) then
          cam_in%cflx(:ncol,m) = cam_in%cflx(:ncol,m) + sflx(:ncol,n)
          call outfld( sflxnam(m), cam_in%cflx(:ncol,m), ncol,lchnk )
       endif
    enddo


  */

  std::cout << "End of surface emissions run" << std::endl;
}

// =========================================================================================
}  // namespace scream
