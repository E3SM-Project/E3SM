#ifndef SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
#define SCREAM_REGISTER_PHYSICS_PROCESSES_HPP

#include "share/atm_process/atmosphere_process.hpp"

// TODO: in the future, you may want to guard these headers,
//       in case we add more options for each parametrization,
//       and we want to only register the ones built,
//       without hardcoding all of them.

#include "physics/p3/eamxx_p3.hpp"
#include "physics/shoc/eamxx_shoc.hpp"
#include "physics/cld_fraction/eamxx_cld_fraction.hpp"
#include "physics/rrtmgp/eamxx_rrtmgp.hpp"
#include "physics/spa/eamxx_spa.hpp"
#include "physics/nudging/eamxx_nudging.hpp"

namespace scream {

inline void register_physics () {
  auto& proc_factory = AtmosphereProcessFactory::instance();

  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
  proc_factory.register_product("SHOC",&create_atmosphere_process<SHOCMacrophysics>);
  proc_factory.register_product("CldFraction",&create_atmosphere_process<CldFraction>);
  proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);
  proc_factory.register_product("SPA",&create_atmosphere_process<SPA>);
  proc_factory.register_product("Nudging",&create_atmosphere_process<Nudging>);
}

} // namespace scream

#endif // SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
