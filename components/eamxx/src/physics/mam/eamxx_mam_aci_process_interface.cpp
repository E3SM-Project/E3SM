#include <physics/mam/eamxx_mam_aci_process_interface.hpp>
#include "scream_config.h" // for SCREAM_CIME_BUILD
#include <share/property_checks/field_lower_bound_check.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

#include "scream_config.h" // for SCREAM_CIME_BUILD

#include <ekat/ekat_assert.hpp>

namespace scream
{

MAMAci::MAMAci(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params){
}

AtmosphereProcessType MAMAci::type() const {
  return AtmosphereProcessType::Physics;
}
//return name of the process
std::string MAMAci::name() const{
  return "mam4_aci";
  }


//grid
void MAMAci::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
}

void MAMAci::initialize_impl(const RunType run_type) {
}

void MAMAci::run_impl(const double dt) {
}

void MAMAci::finalize_impl(){
}

} // namespace scream
