#ifndef SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
#define SCREAM_REGISTER_PHYSICS_PROCESSES_HPP

#include "share/atm_process/atmosphere_process.hpp"

// TODO: in the future, you may want to guard these headers,
//       in case we add more options for each parametrization,
//       and we want to only register the ones built,
//       without hardcoding all of them.

#include "physics/p3/atmosphere_microphysics.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"

#include "physics/common/physics_only_grids_manager.hpp"

namespace scream {

inline void register_physics () {
  auto& proc_factory = AtmosphereProcessFactory::instance();

  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
  proc_factory.register_product("SHOC",&create_atmosphere_process<SHOCMacrophysics>);

  // A physics-only grids manager, in case we run a physics-only test
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Physics Only",&create_physics_only_grids_manager);
}

} // namespace scream

#endif // SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
