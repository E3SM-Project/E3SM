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

namespace scream {
class MAMGenericInterface : public scream::AtmosphereProcess {
 public:
   // Constructor
  MAMGenericInterface(const ekat::Comm &comm,
                       const ekat::ParameterList &params);

  // void add_invariant_check_for_aerosol();
  void add_aerosol_tracers();
  void add_interval_checks();
  void print_fields_names();
    // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
  std::vector<std::string> wet_atm_names_ = {"qv", "qc", "nc", "qi", "ni"};
  std::vector<std::string> dry_atm_names_ = {
        "T_mid",
        "p_mid",
        "p_int",
        "pseudo_density",
        "omega",
        "pbl_height",
        "cldfrac_tot"
    };
  bool check_fields_intervals_{false};

 private:
   // The type of subcomponent
  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }
  std::map<std::string, std::pair<Real, Real>>  limits_aerosol_gas_tracers_;
  void get_aerosol_gas_map();
  const std::pair<Real, Real> get_range(const std::string &field_name);

};  // MAMGenericInterface
}  // namespace scream

#endif  // ifdef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
