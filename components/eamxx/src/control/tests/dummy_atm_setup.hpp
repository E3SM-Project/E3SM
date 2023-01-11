#ifndef SCREAM_DUMMY_ATM_SETUP_HPP
#define SCREAM_DUMMY_ATM_SETUP_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

#include "dummy_atm_proc.hpp"

namespace scream {

inline void dummy_atm_init () {
  using namespace scream;

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("Dummy",&create_atmosphere_process<DummyProcess>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);
}

inline void dummy_atm_cleanup () {
  AtmosphereProcessFactory::instance().clean_up();
  GridsManagerFactory::instance().clean_up();
}

} // namespace scream

#endif // SCREAM_DUMMY_ATM_SETUP_HPP
