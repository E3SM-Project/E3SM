#include <physics/mam/eamxx_mam_generic_process_interface.hpp>

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
void MAMGenericInterface::get_aerosol_gas_map()
{
// NOTE: Using only one range for all num variables.
  // std::map<std::string, std::pair<Real, Real>> limits_aerosol_gas_tracers_;

  const std::string nmr_label="nmr";
  const std::string mmr_label="mmr";

  for(int mode = 0; mode < mam_coupling::num_aero_modes(); ++mode) {
    const std::string int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(mode);
    limits_aerosol_gas_tracers_[int_nmr_field_name]=mam_coupling::physical_min_max(nmr_label);

    const std::string cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(mode);
    limits_aerosol_gas_tracers_[cld_nmr_field_name]=mam_coupling::physical_min_max(nmr_label);

   for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const std::string int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(mode, a);
      if(not int_mmr_field_name.empty()) {
        limits_aerosol_gas_tracers_[int_mmr_field_name]=mam_coupling::physical_min_max(mmr_label);
      }
      const std::string cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(mode, a);
      if(not cld_mmr_field_name.empty()) {
        limits_aerosol_gas_tracers_[cld_mmr_field_name]=mam_coupling::physical_min_max(mmr_label);
      }
    }  // end for loop num species
  }

  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const std::string gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    limits_aerosol_gas_tracers_[gas_mmr_field_name]=mam_coupling::physical_min_max(mmr_label);
   }  // end for loop num gases
}

const std::pair<Real, Real> MAMGenericInterface::get_range(const std::string &field_name)
{
  get_aerosol_gas_map();
  std::pair<Real, Real> min_max;

  auto it = limits_aerosol_gas_tracers_.find(field_name);
  if (it != limits_aerosol_gas_tracers_.end()) {
        min_max = it->second;
  } else {
        min_max = mam_coupling::physical_min_max(field_name);
        // std::cout << "Key not found" << std::endl;
  }
  return min_max;
}

// ================================================================
void MAMGenericInterface::add_aerosol_tracers()
{

  using namespace ekat::units;
  auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  auto n_unit = 1 / kg;   // units of number mixing ratios of tracers
  const auto &grid_name = grid_->name();
  std::cout << "grid_name" << grid_name << "\n";
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
}

// ================================================================
void MAMGenericInterface::populate_wet_and_dry_aero()
{
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


}

void MAMGenericInterface::print_fields_names()
{
  const auto& in_fields = get_fields_in();
  std::cout << "Checking interval pre for..." << "\n";
  for(const auto &item : in_fields) {
    auto& field_name = item.name();
    std::cout << field_name << "\n";
  }
  std::cout << "Checking interval post for..." << "\n";
  const auto& out_fields = get_fields_out();
  for(const auto &item : out_fields) {
    auto& field_name = item.name();
    std::cout << field_name << "\n";
  }
}
// ================================================================
void MAMGenericInterface::add_interval_checks()
{
  if (check_fields_intervals_) {
  const auto& in_fields = get_fields_in();
  // std::cout << "Checking interval pre for..." << "\n";
  for(const auto &item : in_fields) {
    auto& field_name = item.name();
    // std::cout << field_name<< "\n";
    const auto ranges = get_range(field_name);
    const auto min_value = ranges.first;
    const auto max_value = ranges.second;
    add_precondition_check<FieldWithinIntervalCheck>(
      get_field_in(field_name),grid_,min_value,max_value,false);
  }

  const auto& out_fields = get_fields_out();
  // std::cout << "Checking interval post for..." << "\n";
  for(const auto &item : out_fields) {
    auto& field_name = item.name();
    //  std::cout << field_name<< "\n";
    const auto ranges = get_range(field_name);
    const auto min_value = ranges.first;
    const auto max_value = ranges.second;
    add_postcondition_check<FieldWithinIntervalCheck>(
      get_field_out(field_name),grid_,min_value,max_value,false);
  }
  }

}

}  // namespace scream
