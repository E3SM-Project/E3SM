#ifndef EAMXX_MAM_GENERIC_PROCESS_HPP
#define EAMXX_MAM_GENERIC_PROCESS_HPP
// For declaring contituent fluxes class derived from atm process
// class
#include <share/atm_process/atmosphere_process.hpp>
#include <share/property_checks/field_within_interval_check.hpp>
// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>
#include <physics/mam/physical_limits.hpp>
#include <string>
/* We implemented the MAMGenericInterface class to eliminate duplicate code in
the MAM4xx processes. Consequently, all MAM4xx processes must derive from this
class.
*/

namespace scream {
class MAMGenericInterface : public scream::AtmosphereProcess {
 public:
  using KT = ekat::KokkosTypes<DefaultDevice>;
  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

  // Constructor
  MAMGenericInterface(const ekat::Comm &comm,
                      const ekat::ParameterList &params);
  // Add tracers needed for aerosols and gases."
  void add_tracers_interstitial_aerosol_and_gases();
  void add_tracers_cloudborne_aerosol();
  void add_tracers_aerosol_and_gases();
  // Perform interval checks for all MAM4xx fields.
  // The limits are declared in physical_limits.
  void add_interval_checks();
  // Print all fields that are added in a MAM4xx process.
  void print_fields_names();
  //
  void populate_interstitial_wet_and_dry_aero();
  void populate_cloudborne_wet_and_dry_aero();
  // Populate the wet_aero and dry_aero structs.
  void populate_wet_and_dry_aero();
  // Populate the wet_atm and dry_atm struct.
  void populate_wet_and_dry_atm();
  // Add tracers that are needed by the wet_atm and dry_atm.
  void add_tracers_wet_and_dry_atm();

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  void pre_process();
  // Atmosphere processes often have a post-processing step prepares output
  // from this process for the Field Manager. This functor implements this
  // step, which is called during run_impl.
  // Postprocessing functor
  void post_process();
  // Physics grid for column information.
  std::shared_ptr<const AbstractGrid> grid_;
  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;
  bool check_fields_intervals_{false};
    // number of horizontal columns and vertical levels
  int ncol_, nlev_;
 private:
  // The type of subcomponent
  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }
  std::map<std::string, std::pair<Real, Real>> limits_aerosol_gas_tracers_;
  void get_aerosol_gas_map();
  const std::pair<Real, Real> get_range(const std::string &field_name);

};  // MAMGenericInterface
}  // namespace scream

#endif  // ifdef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
