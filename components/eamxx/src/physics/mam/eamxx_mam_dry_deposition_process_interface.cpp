#include "physics/mam/eamxx_mam_dry_deposition_process_interface.hpp"

// Drydep functions are stored in the following hpp file
#include <physics/mam/eamxx_mam_dry_deposition_functions.hpp>

/*
-----------------------------------------------------------------
NOTES:
1. Add a CIME test and multi-process tests
2. Ensure that the submodule for MAM4xx is the main branch
3. Read file for fractional landuse
-----------------------------------------------------------------
*/
namespace scream {

MAMDryDep::MAMDryDep(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMDryDep::set_grids(
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

  // Layout for 2D (2d horiz) variable
  const FieldLayout scalar2d{{COL}, {ncol_}};

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_mid{{COL, LEV}, {ncol_, nlev_}};
  const FieldLayout scalar3d_int{{COL, ILEV}, {ncol_, nlev_ + 1}};

  // Layout for 4D (2d horiz X 1d vertical x number of modes) variables
  // at mid points
  auto make_layout = [](const std::vector<int> &extents,
                        const std::vector<std::string> &names) {
    std::vector<FieldTag> tags(extents.size(), CMP);
    return FieldLayout(tags, extents, names);
  };
  const int num_aero_modes = mam_coupling::num_aero_modes();
  FieldLayout scalar4d_mid =
      make_layout({ncol_, num_aero_modes, nlev_}, {"COL", "num_modes", "lev"});

  using namespace ekat::units;

  auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  auto n_unit = 1 / kg;   // units of number mixing ratios of tracers

  auto nondim = ekat::units::Units::nondimensional();

  auto m3 = m * m * m;  // meter cubed

  // --------------------------------------------------------------------------
  // These variables are "Required" or pure inputs for the process
  // --------------------------------------------------------------------------

  // ----------- Atmospheric quantities -------------
  // Specific humidity [kg/kg](Require only for building DS)
  add_field<Required>("qv", scalar3d_mid, q_unit, grid_name, "tracers");

  // Cloud liquid mass mixing ratio [kg/kg](Require only for building DS)
  add_field<Required>("qc", scalar3d_mid, q_unit, grid_name, "tracers");

  // Cloud ice mass mixing ratio [kg/kg](Require only for building DS)
  add_field<Required>("qi", scalar3d_mid, q_unit, grid_name, "tracers");

  // Cloud liquid number mixing ratio [1/kg](Require only for building DS)
  add_field<Required>("nc", scalar3d_mid, n_unit, grid_name, "tracers");

  // Cloud ice number mixing ratio [1/kg](Require only for building DS)
  add_field<Required>("ni", scalar3d_mid, n_unit, grid_name, "tracers");

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

  static constexpr auto m2 = m * m;
  static constexpr auto s2 = s * s;

  // Surface geopotential [m2/s2] (Require only for building DS)
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  //----------- Variables from microphysics scheme -------------

  // Total cloud fraction [fraction] (Require only for building DS)
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);

  //----------- Variables from coupler (land component)---------
  // Obukhov length [m]
  add_field<Required>("Obukhov_length", scalar2d, m, grid_name);

  // Surface friction velocty or ustar[m/s]
  add_field<Required>("surface_friction_velocty", scalar2d, m / s, grid_name);

  // Land fraction [fraction]
  add_field<Required>("land_fraction", scalar2d, nondim, grid_name);

  // Friction velocity from land model [m/s]
  add_field<Required>("friction_velocity", scalar2d, m / s, grid_name);

  // Aerodynamical resistance from land model [s/m]
  add_field<Required>("aerodynamical_resistance", scalar2d, s / m, grid_name);

  //----------- Variables from coupler (ice component)---------

  // Ice fraction [unitless]
  add_field<Required>("ice_fraction", scalar2d, nondim, grid_name);

  //----------- Variables from coupler (ocean component)---------
  // Ocean fraction [unitless]
  add_field<Required>("ocean_fraction", scalar2d, nondim, grid_name);

  //----------- Variables from other mam4xx processes ------------
  // Geometric mean wet diameter for number distribution [m]
  add_field<Required>("dgnumwet", scalar4d_mid, m, grid_name);

  // Wet density of interstitial aerosol [kg/m3]
  add_field<Required>("wetdens", scalar4d_mid, kg / m3, grid_name);

  // ---------------------------------------------------------------------
  // These variables are "updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  for(int m = 0; m < num_aero_modes; ++m) {
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
  // (cloud) aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for(int m = 0; m < num_aero_modes; ++m) {
    const char *cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);

    add_field<Updated>(cld_nmr_field_name, scalar3d_mid, n_unit, grid_name);
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);

      if(strlen(cld_mmr_field_name) > 0) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_mid, q_unit, grid_name);
      }
    }
  }

  // aerosol-related gases: mass mixing ratios
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_mid, q_unit, grid_name,
                       "tracers");
  }

  // -------------------------------------------------------------
  // These variables are "Computed" or outputs for the process
  // -------------------------------------------------------------
  // FIXME: These are diagnostics, remove them from FM after initial evaluation
  // surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]
  add_field<Computed>("deposition_flux_of_cloud_borne_aerosols", scalar3d_mid,
                      1 / m2 / s, grid_name);
  // surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
  add_field<Computed>("deposition_flux_of_interstitial_aerosols", scalar3d_mid,
                      1 / m2 / s, grid_name);
}  // set_grids

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by
// the above. Buffer type given the number of columns and vertical
// levels
size_t MAMDryDep::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}  // requested_buffer_size_in_bytes

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to
// store intermediate (dry) quantities on the given number of
// columns with the given number of vertical levels. Returns the
// number of bytes allocated.
void MAMDryDep::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMDryDep.");
}  // init_buffers

// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMDryDep::initialize_impl(const RunType run_type) {
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

  // ---- set wet/dry aerosol-related gas state data
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

  //-----------------------------------------------------------------
  // Allocate memory
  //-----------------------------------------------------------------
  const int pcnst = mam4::aero_model::pcnst;
  // FIXME: comment what they are and units.....
  qtracers_ = view_3d("qtracers_", ncol_, nlev_, pcnst);

  rho_     = view_2d("rho", ncol_, nlev_);
  ptend_q_ = view_3d("ptend_q_", ncol_, nlev_, pcnst);

  for(int i = 0; i < mam4::AeroConfig::num_modes(); ++i) {
    for(int j = 0; j < aerosol_categories_; ++j) {
      vlc_dry_[i][j] = view_2d("vlc_dry_[i][j]", ncol_, nlev_);
      vlc_grv_[i][j] = view_2d("vlc_grv_[i][j]", ncol_, nlev_);
      vlc_trb_[i][j] = view_1d("vlc_trb_[i][j]", ncol_);
    }
  }

  for(int i = 0; i < pcnst; ++i) {
    qqcw_[i]     = view_2d("qqcw_[i]", ncol_, nlev_);
    dqdt_tmp_[i] = view_2d("dqdt_tmp_[i]", ncol_, nlev_);
  }

  static constexpr int n_land_type = mam4::DryDeposition::n_land_type;
  for(int i = 0; i < n_land_type; ++i) {
    // FIXME: This should come from a file reading
    fraction_landuse_[i] = view_1d("fraction_landuse_[i]", ncol_);
  }

  //-----------------------------------------------------------------
  // Setup preprocessing and post processing
  //-----------------------------------------------------------------
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                         dry_aero_);
  postprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                          dry_aero_);
}  // initialize_impl

// =========================================================================================
void MAMDryDep::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  // -------------------------------------------------------------
  // Inputs fields for the process
  // -------------------------------------------------------------

  // Geometric mean wet diameter for number distribution [m]
  auto dgncur_awet_ = get_field_in("dgnumwet").get_view<const Real ***>();
  // Wet density of interstitial aerosol [kg/m3]
  auto wet_dens_ = get_field_in("wetdens").get_view<const Real ***>();
  // Obukhov length [m]
  auto obukhov_length_ =
      get_field_in("Obukhov_length").get_view<const Real *>();
  // Land fraction [unitless]
  auto land_fraction_ = get_field_in("land_fraction").get_view<const Real *>();
  // Ice fraction [unitless]
  auto ice_fraction_ = get_field_in("ice_fraction").get_view<const Real *>();
  // Ocean fraction [unitless]
  auto ocean_fraction_ =
      get_field_in("ocean_fraction").get_view<const Real *>();
  // Friction velocity from land model [m/s]
  auto friction_velocity_ =
      get_field_in("friction_velocity").get_view<const Real *>();
  // Aerodynamical resistance from land model [s/m]
  auto aerodynamical_resistance_ =
      get_field_in("aerodynamical_resistance").get_view<const Real *>();
  //  Sfc friction velocity or ustar [m/s]
  auto surface_friction_velocty_ =
      get_field_in("surface_friction_velocty").get_view<const Real *>();

  // -------------------------------------------------------------
  // Output fields for the process
  // -------------------------------------------------------------
  // Surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]
  auto aerdepdrycw_ = get_field_out("deposition_flux_of_cloud_borne_aerosols")
                          .get_view<Real **>();
  // Surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
  auto aerdepdryis_ = get_field_out("deposition_flux_of_interstitial_aerosols")
                          .get_view<Real **>();

  // FIXME: remove it if it read from a file
  populated_fraction_landuse(fraction_landuse_, ncol_);

  // Call drydeposition and get tendencies
  compute_tendencies(ncol_, nlev_, dt, obukhov_length_,
                     surface_friction_velocty_, land_fraction_, ice_fraction_,
                     ocean_fraction_, friction_velocity_,
                     aerodynamical_resistance_, qtracers_, fraction_landuse_,
                     dgncur_awet_, wet_dens_, dry_atm_, dry_aero_,
                     // Inouts-outputs
                     qqcw_,
                     // Outputs
                     ptend_q_, aerdepdrycw_, aerdepdryis_,
                     // work arrays
                     rho_, vlc_dry_, vlc_trb_, vlc_grv_, dqdt_tmp_);
  Kokkos::fence();

  // Update the interstitial aerosols using ptend.
  update_interstitial_mmrs(ptend_q_, dt, ncol_, nlev_,  // inputs
                           dry_aero_);                  // output

  // Update the interstitial aerosols
  update_cloudborne_mmrs(qqcw_, dt, nlev_,  // inputs
                         dry_aero_);        // output

  // call post processing to convert dry mixing ratios to wet mixing ratios
  // and update the state
  Kokkos::parallel_for("postprocess", scan_policy, postprocess_);
  Kokkos::fence();  // wait before returning to calling function
}  // run_impl
}  // namespace scream
