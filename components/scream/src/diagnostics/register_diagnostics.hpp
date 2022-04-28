#ifndef SCREAM_REGISTER_DIAGNOSTICS_HPP
#define SCREAM_REGISTER_DIAGNOSTICS_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
// Include all diagnostics
#include "diagnostics/potential_temperature.hpp"
#include "diagnostics/atm_density.hpp"
#include "diagnostics/exner.hpp"

namespace scream {

inline void register_diagnostics () {
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("PotentialTemperature",&create_atmosphere_diagnostic<PotentialTemperatureDiagnostic>);
  diag_factory.register_product("AtmosphereDensity",&create_atmosphere_diagnostic<AtmDensityDiagnostic>);
  diag_factory.register_product("Exner",&create_atmosphere_diagnostic<ExnerDiagnostic>);
}

} // namespace scream
#endif // SCREAM_REGISTER_DIAGNOSTICS_HPP
