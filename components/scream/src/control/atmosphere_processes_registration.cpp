// Inlcude all atmosphere processes headers here
#include <control/surface_coupling.hpp>
#include <control/atmosphere_process_group.hpp>
#include <dynamics/atmosphere_dynamics.hpp>

namespace scream
{

// This routine registers all the atmosphere processes in the AtmosphereProcessFactory,
// so that we can later create them by invoking the create method of the factory.
void register_all_atm_processes ()
{
  bool success = true;

  success &= AtmosphereProcessFactory::instance().register_product("Surface Coupling", &create_surface_coupling);
  success &= AtmosphereProcessFactory::instance().register_product("Process Group",    &create_process_group);
  success &= AtmosphereProcessFactory::instance().register_product("Dynamics",         &create_atmosphere_dynamics);

  error::runtime_check(success, "Error! Something went wrong while registering atmosphere processes.\n");
}

} // namespace scream
