#include <physics/mam/eamxx_mam_generic_process_interface.hpp>
#include <physics/mam/physical_limits.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

namespace scream {
// ================================================================
//  Constructor
// ================================================================

MAMGenericInterface::MAMGenericInterface(const ekat::Comm &comm,
                                         const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}
// ================================================================
void MAMGenericInterface::set_aerosol_and_gas_ranges() {
  // NOTE: Using only one range for all num variables.
  // std::map<std::string, std::pair<Real, Real>> limits_aerosol_gas_tracers_;

  const std::string nmr_label = "nmr";
  const std::string mmr_label = "mmr";

  for(int mode = 0; mode < mam_coupling::num_aero_modes(); ++mode) {
    const std::string int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(mode);
    limits_aerosol_gas_tracers_[int_nmr_field_name] =
        mam_coupling::physical_min_max(nmr_label);

    const std::string cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(mode);
    limits_aerosol_gas_tracers_[cld_nmr_field_name] =
        mam_coupling::physical_min_max(nmr_label);

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(mode, a);
      if(not int_mmr_field_name.empty()) {
        limits_aerosol_gas_tracers_[int_mmr_field_name] =
            mam_coupling::physical_min_max(mmr_label);
      }
      const std::string cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(mode, a);
      if(not cld_mmr_field_name.empty()) {
        limits_aerosol_gas_tracers_[cld_mmr_field_name] =
            mam_coupling::physical_min_max(mmr_label);
      }
    }  // end for loop num species
  }

  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const std::string gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    limits_aerosol_gas_tracers_[gas_mmr_field_name] =
        mam_coupling::physical_min_max(mmr_label);
  }  // end for loop num gases
}
void MAMGenericInterface::set_ranges_process(
    const std::map<std::string, std::pair<Real, Real>> &max_min_process) {
  // NOTE: We are using the same range (mmr) for all aerosols and gases.
  // And aerosol numbers (nmr).
  // set ranges for aerosol and gases
  // populates limits_aerosol_gas_tracers_
  set_aerosol_and_gas_ranges();
  // set ranges for other variables
  max_min_process_ = max_min_process;
  // Ensure that we populate the maps.
  set_ranges_ = true;
}
const std::pair<Real, Real> MAMGenericInterface::get_ranges(
    const std::string &field_name) {
  EKAT_ASSERT_MSG(
      set_ranges_,
      "Error: min_max ranges are not set. Please invoke set_ranges_process.");

  std::pair<Real, Real> min_max;
  // We obtain the minimum and maximum values for aerosol and gas species.
  auto it = limits_aerosol_gas_tracers_.find(field_name);
  if(it != limits_aerosol_gas_tracers_.end()) {
    min_max = it->second;
  } else {
    // Next, we proceed to check the other variables.
    auto it_process = max_min_process_.find(field_name);
    if(it_process != max_min_process_.end()) {
      min_max = it_process->second;
    } else {
      // If we do not find a variable name in the previous maps,
      // we return a pair (-1, -1) that will bypass the interval check.
      min_max = std::make_pair(-1, -1);
    }
  }
  return min_max;
}

// ================================================================
void MAMGenericInterface::add_fields_cloudborne_aerosol() {
  using namespace ekat::units;
  auto q_unit           = kg / kg;  // units of mass mixing ratios of tracers
  auto n_unit           = 1 / kg;   // units of number mixing ratios of tracers
  const auto &grid_name = grid_->name();

  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);

  // ---------------------------------------------------------------------
  // These variables are "Updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------
  // NOTE: Cloud borne aerosols are not updated in this process but are included
  // to create data structures.

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int mode = 0; mode < mam_coupling::num_aero_modes(); ++mode) {
    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
    // NOT advected
    const std::string cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(mode);
    add_field<Updated>(cld_nmr_field_name, scalar3d_mid, n_unit, grid_name);

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
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
}

void MAMGenericInterface::add_tracers_interstitial_aerosol() {
  using namespace ekat::units;
  auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  auto n_unit = 1 / kg;   // units of number mixing ratios of tracers

  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);

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
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(mode, a);
      if(not int_mmr_field_name.empty()) {
        add_tracer<Updated>(int_mmr_field_name, grid_, q_unit);
      }
    }  // end for loop num species
  }    // end for loop for num modes
}
// ================================================================

void MAMGenericInterface::add_tracers_gases() {
  using namespace ekat::units;
  auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const std::string gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_tracer<Updated>(gas_mmr_field_name, grid_, q_unit);
  }  // end for loop num gases
}
// ================================================================
void MAMGenericInterface::populate_cloudborne_wet_aero(
    mam_coupling::AerosolState &wet_aero) {
  // cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    const std::string cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(m);
    wet_aero.cld_aero_nmr[m] =
        get_field_out(cld_nmr_field_name).get_view<Real **>();

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const std::string cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if(not cld_mmr_field_name.empty()) {
        wet_aero.cld_aero_mmr[m][a] =
            get_field_out(cld_mmr_field_name).get_view<Real **>();
      }
    }
  }
}

// ================================================================
void MAMGenericInterface::populate_cloudborne_dry_aero(
    mam_coupling::AerosolState &dry_aero, mam_coupling::Buffer &buffer) {
  // cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    dry_aero.cld_aero_nmr[m] = buffer.dry_cld_aero_nmr[m];

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const std::string cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if(not cld_mmr_field_name.empty()) {
        dry_aero.cld_aero_mmr[m][a] = buffer.dry_cld_aero_mmr[m][a];
      }
    }
  }
}
// ================================================================
void MAMGenericInterface::set_field_w_scratch_buffer(mam_coupling::view_2d& var,
  mam_coupling::Buffer &buffer, const bool set_to_zero)
{
  var = buffer.scratch[i_scratch_vars_];
  i_scratch_vars_++;
  EKAT_REQUIRE_MSG(
      i_scratch_vars_ < buffer.max_num_2d_scratch ,
      "Error! Insufficient number of scratch size in mam buffer.\n");
  if (set_to_zero){
    Kokkos::deep_copy(var, 0.0);
  }
}
// ================================================================
void MAMGenericInterface::populate_gases_dry_aero(
    mam_coupling::AerosolState &dry_aero, mam_coupling::Buffer &buffer) {
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    dry_aero.gas_mmr[g] = buffer.dry_gas_mmr[g];
  }
}

// ================================================================
void MAMGenericInterface::populate_gases_wet_aero(
    mam_coupling::AerosolState &wet_aero) {
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const std::string gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero.gas_mmr[g] = get_field_out(gas_mmr_field_name).get_view<Real **>();
  }
}

// ================================================================
void MAMGenericInterface::populate_interstitial_dry_aero(
    mam_coupling::AerosolState &dry_aero, mam_coupling::Buffer &buffer) {
  // interstitial aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    dry_aero.int_aero_nmr[m] = buffer.dry_int_aero_nmr[m];

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(not int_mmr_field_name.empty()) {
        dry_aero.int_aero_mmr[m][a] = buffer.dry_int_aero_mmr[m][a];
      }
    }
  }
}

// ================================================================
void MAMGenericInterface::populate_interstitial_wet_aero(
    mam_coupling::AerosolState &wet_aero) {
  // interstitial aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const std::string int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(m);
    wet_aero.int_aero_nmr[m] =
        get_field_out(int_nmr_field_name).get_view<Real **>();
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(not int_mmr_field_name.empty()) {
        wet_aero.int_aero_mmr[m][a] =
            get_field_out(int_mmr_field_name).get_view<Real **>();
      }
    }
  }
}

void MAMGenericInterface::populate_wet_atm(
    mam_coupling::WetAtmosphere &wet_atm) {
  // store fields only to be converted to dry mmrs in wet_atm_
  wet_atm.qv = get_field_in("qv").get_view<const Real **>();
  wet_atm.qc = get_field_in("qc").get_view<const Real **>();
  wet_atm.nc = get_field_in("nc").get_view<const Real **>();
  wet_atm.qi = get_field_in("qi").get_view<const Real **>();
  wet_atm.ni = get_field_in("ni").get_view<const Real **>();
}
void MAMGenericInterface::populate_dry_atm(mam_coupling::DryAtmosphere &dry_atm,
                                           mam_coupling::Buffer &buffer) {
  // store rest fo the atm fields in dry_atm_in
  dry_atm.z_surf = 0;
  dry_atm.T_mid  = get_field_in("T_mid").get_view<const Real **>();
  dry_atm.p_mid  = get_field_in("p_mid").get_view<const Real **>();
  dry_atm.p_int  = get_field_in("p_int").get_view<const Real **>();
  dry_atm.p_del  = get_field_in("pseudo_density").get_view<const Real **>();
  dry_atm.omega  = get_field_in("omega").get_view<const Real **>();

  // store fields converted to dry mmr from wet mmr in dry_atm
  dry_atm.qv = buffer.qv_dry;
  dry_atm.qc = buffer.qc_dry;
  dry_atm.nc = buffer.nc_dry;
  dry_atm.qi = buffer.qi_dry;
  dry_atm.ni = buffer.ni_dry;

  // pbl_height
  dry_atm.pblh = get_field_in("pbl_height").get_view<const Real *>();

  // geometric thickness of layers (m)
  dry_atm.dz = buffer.dz;

  // geopotential height above surface at interface levels (m)
  dry_atm.z_iface = buffer.z_iface;

  // geopotential height above surface at mid levels (m)
  dry_atm.z_mid = buffer.z_mid;

  // total cloud fraction
  dry_atm.cldfrac = get_field_in("cldfrac_tot").get_view<const Real **>();

  // computed updraft velocity
  dry_atm.w_updraft = buffer.w_updraft;
}

void MAMGenericInterface::add_tracers_wet_atm() {
  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;
  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  using namespace ekat::units;
  constexpr auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  constexpr auto n_unit = 1 / kg;   // units of number mixing ratios of tracers

  // atmospheric quantities
  // specific humidity [kg/kg]
  add_tracer<Required>("qv", grid_, q_unit);

  // cloud liquid mass mixing ratio [kg/kg]
  add_tracer<Required>("qc", grid_, q_unit);

  // cloud ice mass mixing ratio [kg/kg]
  add_tracer<Required>("qi", grid_, q_unit);

  // cloud liquid number mixing ratio [1/kg]
  add_tracer<Required>("nc", grid_, n_unit);

  // cloud ice number mixing ratio [1/kg]
  add_tracer<Required>("ni", grid_, n_unit);
}

void MAMGenericInterface::add_fields_dry_atm() {
  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;
  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const auto &grid_name = grid_->name();
  const int ncol =
      grid_->get_num_local_dofs();  // Number of columns on this rank
  const FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);
  const FieldLayout scalar3d_int = grid_->get_3d_scalar_layout(false);
  // layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{{COL}, {ncol}};
  using namespace ekat::units;
  constexpr auto nondim = Units::nondimensional();

  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);

  // Vertical pressure velocity [Pa/s] at midpoints
  add_field<Required>("omega", scalar3d_mid, Pa / s, grid_name);

  // Total pressure [Pa] at midpoints
  add_field<Required>("p_mid", scalar3d_mid, Pa, grid_name);

  // Total pressure [Pa] at interfaces
  add_field<Required>("p_int", scalar3d_int, Pa, grid_name);

  // Layer thickness(pdel) [Pa] at midpoints
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);

  // planetary boundary layer height
  add_field<Required>("pbl_height", scalar2d_layout_col, m, grid_name);

  // cloud fraction [nondimensional] computed by eamxx_cld_fraction_process
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);
}

// ================================================================
void MAMGenericInterface::add_interval_checks() {
  if(check_fields_intervals_) {
    const auto &in_fields = get_fields_in();
    for(const auto &field : in_fields) {
      const auto ranges    = get_ranges(field.name());
      const auto min_value = ranges.first;
      const auto max_value = ranges.second;
      // Only check the range for variables that are specified in
      // physical_limits.hpp. Some variables from get_fields_in may not be part
      // of physical_min_max.
      if(min_value != -1 && max_value != -1)
        add_precondition_check<FieldWithinIntervalCheck>(
            field, grid_, min_value, max_value, false);
    }

    const auto &out_fields = get_fields_out();
    for(const auto &field : out_fields) {
      const auto ranges    = get_ranges(field.name());
      const auto min_value = ranges.first;
      const auto max_value = ranges.second;
      // Only check the range for variables that are specified in
      // physical_limits.hpp. Some variables from get_fields_out may not be part
      // of physical_min_max.
      if(min_value != -1 && max_value != -1)
        add_postcondition_check<FieldWithinIntervalCheck>(
            field, grid_, min_value, max_value, false);
    }
  }
}

void MAMGenericInterface::pre_process(mam_coupling::AerosolState &wet_aero,
                                      mam_coupling::AerosolState &dry_aero,
                                      mam_coupling::WetAtmosphere &wet_atm,
                                      mam_coupling::DryAtmosphere &dry_atm) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);
  Kokkos::parallel_for(
      scan_policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int i = team.league_rank();  // column index

        mam_coupling::compute_dry_mixing_ratios(team, wet_atm, dry_atm, i);
        mam_coupling::compute_dry_mixing_ratios(team, wet_atm, wet_aero,
                                                dry_aero, i);
        team.team_barrier();
        // vertical heights has to be computed after computing dry mixing ratios
        // for atmosphere
        mam_coupling::compute_vertical_layer_heights(team, dry_atm, i);
        mam_coupling::compute_updraft_velocities(team, wet_atm, dry_atm, i);
        // allows kernels below to use layer heights operator()
        team.team_barrier();
      });
}

void MAMGenericInterface::post_process(mam_coupling::AerosolState &wet_aero,
                                       mam_coupling::AerosolState &dry_aero,
                                       mam_coupling::DryAtmosphere &dry_atm) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);
  Kokkos::parallel_for(
      scan_policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int i = team.league_rank();  // column index
        compute_wet_mixing_ratios(team, dry_atm, dry_aero, wet_aero, i);
      });
}
}  // namespace scream
