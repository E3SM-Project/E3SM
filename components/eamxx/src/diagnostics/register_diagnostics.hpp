#ifndef SCREAM_REGISTER_DIAGNOSTICS_HPP
#define SCREAM_REGISTER_DIAGNOSTICS_HPP

// Include all diagnostics
#include "diagnostics/field_at_level.hpp"
#include "diagnostics/potential_temperature.hpp"
#include "diagnostics/atm_density.hpp"
#include "diagnostics/exner.hpp"
#include "diagnostics/virtual_temperature.hpp"
#include "diagnostics/vertical_layer.hpp"
#include "diagnostics/dry_static_energy.hpp"
#include "diagnostics/sea_level_pressure.hpp"
#include "diagnostics/liquid_water_path.hpp"
#include "diagnostics/ice_water_path.hpp"
#include "diagnostics/rime_water_path.hpp"
#include "diagnostics/vapor_water_path.hpp"
#include "diagnostics/rain_water_path.hpp"
#include "diagnostics/shortwave_cloud_forcing.hpp"
#include "diagnostics/longwave_cloud_forcing.hpp"
#include "diagnostics/relative_humidity.hpp"
#include "diagnostics/zonal_vapor_flux.hpp"
#include "diagnostics/meridional_vapor_flux.hpp"
#include "diagnostics/field_at_pressure_level.hpp"
#include "diagnostics/precip_liq_surf_mass_flux.hpp"
#include "diagnostics/precip_ice_surf_mass_flux.hpp"
#include "diagnostics/precip_total_surf_mass_flux.hpp"

namespace scream {

inline void register_diagnostics () {
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("PotentialTemperature",&create_atmosphere_diagnostic<PotentialTemperatureDiagnostic>);
  diag_factory.register_product("FieldAtLevel",&create_atmosphere_diagnostic<FieldAtLevel>);
  diag_factory.register_product("FieldAtPressureLevel",&create_atmosphere_diagnostic<FieldAtPressureLevel>);
  diag_factory.register_product("AtmosphereDensity",&create_atmosphere_diagnostic<AtmDensityDiagnostic>);
  diag_factory.register_product("Exner",&create_atmosphere_diagnostic<ExnerDiagnostic>);
  diag_factory.register_product("VirtualTemperature",&create_atmosphere_diagnostic<VirtualTemperatureDiagnostic>);
  diag_factory.register_product("z_int",&create_atmosphere_diagnostic<VerticalLayerDiagnostic>);
  diag_factory.register_product("geopotential_int",&create_atmosphere_diagnostic<VerticalLayerDiagnostic>);
  diag_factory.register_product("z_mid",&create_atmosphere_diagnostic<VerticalLayerDiagnostic>);
  diag_factory.register_product("geopotential_mid",&create_atmosphere_diagnostic<VerticalLayerDiagnostic>);
  diag_factory.register_product("dz",&create_atmosphere_diagnostic<VerticalLayerDiagnostic>);
  diag_factory.register_product("DryStaticEnergy",&create_atmosphere_diagnostic<DryStaticEnergyDiagnostic>);
  diag_factory.register_product("SeaLevelPressure",&create_atmosphere_diagnostic<SeaLevelPressureDiagnostic>);
  diag_factory.register_product("LiqWaterPath",&create_atmosphere_diagnostic<LiqWaterPathDiagnostic>);
  diag_factory.register_product("IceWaterPath",&create_atmosphere_diagnostic<IceWaterPathDiagnostic>);
  diag_factory.register_product("VapWaterPath",&create_atmosphere_diagnostic<VapWaterPathDiagnostic>);
  diag_factory.register_product("RainWaterPath",&create_atmosphere_diagnostic<RainWaterPathDiagnostic>);
  diag_factory.register_product("RimeWaterPath",&create_atmosphere_diagnostic<RimeWaterPathDiagnostic>);
  diag_factory.register_product("ShortwaveCloudForcing",&create_atmosphere_diagnostic<ShortwaveCloudForcingDiagnostic>);
  diag_factory.register_product("LongwaveCloudForcing",&create_atmosphere_diagnostic<LongwaveCloudForcingDiagnostic>);
  diag_factory.register_product("RelativeHumidity",&create_atmosphere_diagnostic<RelativeHumidityDiagnostic>);
  diag_factory.register_product("ZonalVapFlux",&create_atmosphere_diagnostic<ZonalVapFluxDiagnostic>);
  diag_factory.register_product("MeridionalVapFlux",&create_atmosphere_diagnostic<MeridionalVapFluxDiagnostic>);
  diag_factory.register_product("PrecipLiqSurfMassFlux",&create_atmosphere_diagnostic<PrecipLiqSurfMassFluxDiagnostic>);
  diag_factory.register_product("PrecipIceSurfMassFlux",&create_atmosphere_diagnostic<PrecipIceSurfMassFluxDiagnostic>);
  diag_factory.register_product("PrecipTotalSurfMassFlux",&create_atmosphere_diagnostic<PrecipTotalSurfMassFluxDiagnostic>);
}

} // namespace scream
#endif // SCREAM_REGISTER_DIAGNOSTICS_HPP
