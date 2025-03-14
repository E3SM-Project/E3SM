#ifndef EAMXX_MAM_GENERIC_PROCESS_HPP
#define EAMXX_MAM_GENERIC_PROCESS_HPP
// For declaring contituent fluxes class derived from atm process
// class
#include <share/atm_process/atmosphere_process.hpp>
// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>
#include <string>
/* We implemented the MAMGenericInterface class to eliminate duplicate code in
the MAM4xx processes. Consequently, all MAM4xx processes must derive from this
class.
*/

namespace scream {
class MAMGenericInterface : public scream::AtmosphereProcess {
 public:
  // Constructor
  MAMGenericInterface(const ekat::Comm &comm,
                      const ekat::ParameterList &params);

  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }

 protected:
  using KT = ekat::KokkosTypes<DefaultDevice>;
  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;
  // Add tracers needed for aerosols and gases."
  void add_tracers_interstitial_aerosol();
  void add_tracers_gases();
  void add_fields_cloudborne_aerosol();
  // Perform interval checks for all MAM4xx fields.
  // The limits are declared in physical_limits.
  void add_interval_checks();
  // Populate the wet_aero and dry_aero structs.
  void populate_interstitial_wet_aero(mam_coupling::AerosolState &wet_aero);
  void populate_interstitial_dry_aero(mam_coupling::AerosolState &dry_aero,
                                      mam_coupling::Buffer &buffer);
  void populate_gases_dry_aero(mam_coupling::AerosolState &dry_aero,
                               mam_coupling::Buffer &buffer);
  void populate_gases_wet_aero(mam_coupling::AerosolState &wet_aero);
  void populate_cloudborne_wet_aero(mam_coupling::AerosolState &wet_aero);
  void populate_cloudborne_dry_aero(mam_coupling::AerosolState &dry_aero,
                                    mam_coupling::Buffer &buffer);

  // Populate the wet_atm and dry_atm struct.
  void populate_dry_atm(mam_coupling::DryAtmosphere &dry_atm,
                        mam_coupling::Buffer &buffer);
  void populate_wet_atm(mam_coupling::WetAtmosphere &wet_atm);
  void add_fields_dry_atm();
  void add_tracers_wet_atm();

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  void pre_process(mam_coupling::AerosolState &wet_aero,
                   mam_coupling::AerosolState &dry_aero,
                   mam_coupling::WetAtmosphere &wet_atm,
                   mam_coupling::DryAtmosphere &dry_atm);
  // Atmosphere processes often have a post-processing step prepares output
  // from this process for the Field Manager. This functor implements this
  // step, which is called during run_impl.
  // Postprocessing functor
  void post_process(mam_coupling::AerosolState &wet_aero,
                    mam_coupling::AerosolState &dry_aero,
                    mam_coupling::DryAtmosphere &dry_atm);
  // Physics grid for column information.
  void set_field_w_scratch_buffer(mam_coupling::view_2d& var,
  mam_coupling::Buffer &buffer, const bool set_to_zero);
  std::shared_ptr<const AbstractGrid> grid_;
  bool check_fields_intervals_{false};
  // number of horizontal columns and vertical levels
  int ncol_, nlev_;
  void set_ranges_process(
      const std::map<std::string, std::pair<Real, Real>> &max_min_process);

 private:
  // The type of subcomponent
  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------
  std::map<std::string, std::pair<Real, Real>> limits_aerosol_gas_tracers_;
  void set_aerosol_and_gas_ranges();
  const std::pair<Real, Real> get_ranges(const std::string &field_name);
  std::map<std::string, std::pair<Real, Real>> max_min_process_;
  bool set_ranges_{false};
  int i_scratch_vars_{0};

};  // MAMGenericInterface
}  // namespace scream

#endif  // ifdef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
