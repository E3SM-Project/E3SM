#ifndef SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
#define SCREAM_REGISTER_PHYSICS_PROCESSES_HPP

#include "share/atm_process/atmosphere_process.hpp"

// Only include headers/register processes which
// have been built.

#ifdef EAMXX_HAS_P3
#include "physics/p3/eamxx_p3_process_interface.hpp"
#endif
#ifdef EAMXX_HAS_SHOC
#include "physics/shoc/eamxx_shoc_process_interface.hpp"
#endif
#ifdef EAMXX_HAS_CLD_FRACTION
#include "physics/cld_fraction/eamxx_cld_fraction_process_interface.hpp"
#endif
#ifdef EAMXX_HAS_RRTMGP
#include "physics/rrtmgp/eamxx_rrtmgp_process_interface.hpp"
#endif
#ifdef EAMXX_HAS_SPA
#include "physics/spa/eamxx_spa_process_interface.hpp"
#endif
#ifdef EAMXX_HAS_NUDGING
#include "physics/nudging/eamxx_nudging_process_interface.hpp"
#endif

namespace scream {

inline void register_physics () {
  auto& proc_factory = AtmosphereProcessFactory::instance();
#ifdef EAMXX_HAS_P3
  proc_factory.register_product("p3",&create_atmosphere_process<P3Microphysics>);
#endif
#ifdef EAMXX_HAS_SHOC
  proc_factory.register_product("SHOC",&create_atmosphere_process<SHOCMacrophysics>);
#endif
#ifdef EAMXX_HAS_CLD_FRACTION
  proc_factory.register_product("CldFraction",&create_atmosphere_process<CldFraction>);
#endif
#ifdef EAMXX_HAS_RRTMGP
  proc_factory.register_product("RRTMGP",&create_atmosphere_process<RRTMGPRadiation>);
#endif
#ifdef EAMXX_HAS_SPA
  proc_factory.register_product("SPA",&create_atmosphere_process<SPA>);
#endif
#ifdef EAMXX_HAS_NUDGING
  proc_factory.register_product("Nudging",&create_atmosphere_process<Nudging>);
#endif
}

} // namespace scream

#endif // SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
