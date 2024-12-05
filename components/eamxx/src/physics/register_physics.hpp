#ifndef SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
#define SCREAM_REGISTER_PHYSICS_PROCESSES_HPP

#include "share/atm_process/atmosphere_process.hpp"

// Only include headers and register processes for libs that have been linked in

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
#ifdef EAMXX_HAS_MAM
#include "physics/mam/eamxx_mam_microphysics_process_interface.hpp"
#include "physics/mam/eamxx_mam_optics_process_interface.hpp"
#include "physics/mam/eamxx_mam_dry_deposition_process_interface.hpp"
#include "physics/mam/eamxx_mam_aci_process_interface.hpp"
#include "physics/mam/eamxx_mam_wetscav_process_interface.hpp"
#include "physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp"
#include "physics/mam/eamxx_mam_constituent_fluxes_interface.hpp"
#endif
#ifdef EAMXX_HAS_COSP
#include "physics/cosp/eamxx_cosp.hpp"
#endif
#ifdef EAMXX_HAS_TMS
#include "physics/tms/eamxx_tms_process_interface.hpp"
#endif
#ifdef EAMXX_HAS_ML_CORRECTION
#include "physics/ml_correction/eamxx_ml_correction_process_interface.hpp"
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
#ifdef EAMXX_HAS_MAM
  proc_factory.register_product("mam4_aero_microphys",&create_atmosphere_process<MAMMicrophysics>);
  proc_factory.register_product("mam4_optics",&create_atmosphere_process<MAMOptics>);
  proc_factory.register_product("mam4_drydep",&create_atmosphere_process<MAMDryDep>);
  proc_factory.register_product("mam4_aci",&create_atmosphere_process<MAMAci>);
  proc_factory.register_product("mam4_wetscav",&create_atmosphere_process<MAMWetscav>);
  proc_factory.register_product("mam4_srf_online_emiss",&create_atmosphere_process<MAMSrfOnlineEmiss>);
  proc_factory.register_product("mam4_constituent_fluxes",&create_atmosphere_process<MAMConstituentFluxes>);
#endif
#ifdef EAMXX_HAS_COSP
  proc_factory.register_product("Cosp",&create_atmosphere_process<Cosp>);
#endif
#ifdef EAMXX_HAS_TMS
  proc_factory.register_product("tms",&create_atmosphere_process<TurbulentMountainStress>);
#endif
#ifdef EAMXX_HAS_ML_CORRECTION
  proc_factory.register_product("MLCorrection",&create_atmosphere_process<MLCorrection>);
#endif

  // If no physics was enabled, silence compile warning about unused var
  (void) proc_factory;
}

} // namespace scream

#endif // SCREAM_REGISTER_PHYSICS_PROCESSES_HPP
