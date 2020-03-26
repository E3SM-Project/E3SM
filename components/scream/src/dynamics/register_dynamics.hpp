#ifndef SCREAM_REGISTER_DYNAMICS_PROCESS_HPP
#define SCREAM_REGISTER_DYNAMICS_PROCESS_HPP

#include "share/atmosphere_process.hpp"

#ifdef SCREAM_HAS_HOMME
#include "homme/atmosphere_dynamics.hpp"
#endif

namespace scream {

inline void register_dynamics () {
  auto& proc_factory = AtmosphereProcessFactory::instance();
#ifdef SCREAM_HAS_HOMME
  proc_factory.register_product("Homme",&create_atmosphere_process<HommeDynamics>);
#endif
  (void) proc_factory;
}

} // namespace scream

#endif // SCREAM_REGISTER_DYNAMICS_PROCESS_HPP
