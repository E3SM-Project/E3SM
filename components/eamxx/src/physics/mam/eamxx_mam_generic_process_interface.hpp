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

  void add_invariant_check_for_aerosol();
  void add_aerosol_tracers();
    // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

 private:
   // The type of subcomponent
  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }



};  // MAMGenericInterface
}  // namespace scream

#endif  // ifdef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
