#ifndef SCREAM_REGISTER_DIAGNOSTICS_HPP
#define SCREAM_REGISTER_DIAGNOSTICS_HPP

// Include all diagnostics
#include "diagnostics/field_at_level.hpp"
#include "diagnostics/field_at_height.hpp"
#include "diagnostics/potential_temperature.hpp"
#include "diagnostics/atm_density.hpp"
#include "diagnostics/exner.hpp"
#include "diagnostics/virtual_temperature.hpp"
#include "diagnostics/vertical_layer.hpp"
#include "diagnostics/dry_static_energy.hpp"
#include "diagnostics/sea_level_pressure.hpp"
#include "diagnostics/water_path.hpp"
#include "diagnostics/shortwave_cloud_forcing.hpp"
#include "diagnostics/longwave_cloud_forcing.hpp"
#include "diagnostics/relative_humidity.hpp"
#include "diagnostics/vapor_flux.hpp"
#include "diagnostics/field_at_pressure_level.hpp"
#include "diagnostics/precip_surf_mass_flux.hpp"
#include "diagnostics/surf_upward_latent_heat_flux.hpp"
#include "diagnostics/wind_speed.hpp"
#include "diagnostics/aodvis.hpp"
#include "diagnostics/number_path.hpp"
#include "diagnostics/aerocom_cld.hpp"
#include "diagnostics/atm_backtend.hpp"
#include "diagnostics/horiz_avg.hpp"

namespace scream {

inline void register_diagnostics () {
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("PotentialTemperature",&create_atmosphere_diagnostic<PotentialTemperatureDiagnostic>);
  diag_factory.register_product("FieldAtLevel",&create_atmosphere_diagnostic<FieldAtLevel>);
  diag_factory.register_product("FieldAtHeight",&create_atmosphere_diagnostic<FieldAtHeight>);
  diag_factory.register_product("FieldAtPressureLevel",&create_atmosphere_diagnostic<FieldAtPressureLevel>);
  diag_factory.register_product("AtmosphereDensity",&create_atmosphere_diagnostic<AtmDensityDiagnostic>);
  diag_factory.register_product("Exner",&create_atmosphere_diagnostic<ExnerDiagnostic>);
  diag_factory.register_product("VirtualTemperature",&create_atmosphere_diagnostic<VirtualTemperatureDiagnostic>);
  diag_factory.register_product("DryStaticEnergy",&create_atmosphere_diagnostic<DryStaticEnergyDiagnostic>);
  diag_factory.register_product("SeaLevelPressure",&create_atmosphere_diagnostic<SeaLevelPressureDiagnostic>);
  diag_factory.register_product("WaterPath",&create_atmosphere_diagnostic<WaterPathDiagnostic>);
  diag_factory.register_product("ShortwaveCloudForcing",&create_atmosphere_diagnostic<ShortwaveCloudForcingDiagnostic>);
  diag_factory.register_product("LongwaveCloudForcing",&create_atmosphere_diagnostic<LongwaveCloudForcingDiagnostic>);
  diag_factory.register_product("RelativeHumidity",&create_atmosphere_diagnostic<RelativeHumidityDiagnostic>);
  diag_factory.register_product("VaporFlux",&create_atmosphere_diagnostic<VaporFluxDiagnostic>);
  diag_factory.register_product("VerticalLayer",&create_atmosphere_diagnostic<VerticalLayerDiagnostic>);
  diag_factory.register_product("precip_surf_mass_flux",&create_atmosphere_diagnostic<PrecipSurfMassFlux>);
  diag_factory.register_product("surface_upward_latent_heat_flux",&create_atmosphere_diagnostic<SurfaceUpwardLatentHeatFlux>);
  diag_factory.register_product("wind_speed",&create_atmosphere_diagnostic<WindSpeed>);
  diag_factory.register_product("AerosolOpticalDepth550nm",&create_atmosphere_diagnostic<AODVis>);
  diag_factory.register_product("NumberPath",&create_atmosphere_diagnostic<NumberPathDiagnostic>);
  diag_factory.register_product("AeroComCld",&create_atmosphere_diagnostic<AeroComCld>);
  diag_factory.register_product("AtmBackTendDiag",&create_atmosphere_diagnostic<AtmBackTendDiag>);
  diag_factory.register_product("HorizAvgDiag",&create_atmosphere_diagnostic<HorizAvgDiag>);
}

} // namespace scream

#endif // SCREAM_REGISTER_DIAGNOSTICS_HPP
