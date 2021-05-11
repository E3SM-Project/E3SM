#ifndef SCREAM_REGISTER_DYNAMICS_PROCESS_HPP
#define SCREAM_REGISTER_DYNAMICS_PROCESS_HPP

#include "share/atm_process/atmosphere_process.hpp"

#ifdef SCREAM_HAS_HOMME
#include "homme/atmosphere_dynamics.hpp"
#include "homme/dynamics_driven_grids_manager.hpp"
#endif

namespace scream {

inline void register_dynamics () {
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory   = GridsManagerFactory::instance();

#ifdef SCREAM_HAS_HOMME
  proc_factory.register_product("Homme",&create_atmosphere_process<HommeDynamics>);

  gm_factory.register_product("Dynamics Driven",&create_dynamics_driven_grids_manager);
#endif
  (void) proc_factory;
  (void) gm_factory;
}

} // namespace scream

#endif // SCREAM_REGISTER_DYNAMICS_PROCESS_HPP
