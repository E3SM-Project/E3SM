#ifndef SCREAM_REGISTER_DIAGNOSTICS_HPP
#define SCREAM_REGISTER_DIAGNOSTICS_HPP

// Include all diagnostics
#include "field_at_level.hpp"
#include "field_at_height.hpp"
#include "potential_temperature.hpp"
#include "atm_density.hpp"
#include "exner.hpp"
#include "virtual_temperature.hpp"
#include "vertical_layer.hpp"
#include "dry_static_energy.hpp"
#include "sea_level_pressure.hpp"
#include "water_path.hpp"
#include "shortwave_cloud_forcing.hpp"
#include "longwave_cloud_forcing.hpp"
#include "relative_humidity.hpp"
#include "vapor_flux.hpp"
#include "field_at_pressure_level.hpp"
#include "precip_surf_mass_flux.hpp"
#include "surf_upward_latent_heat_flux.hpp"
#include "wind_speed.hpp"
#include "aodvis.hpp"
#include "number_path.hpp"
#include "aerocom_cld.hpp"
#include "atm_backtend.hpp"
#include "horiz_avg.hpp"
#include "vert_contract.hpp"
#include "vert_derivative.hpp"
#include "zonal_avg.hpp"
#include "conditional_sampling.hpp"
#include "binary_ops.hpp"
#include "histogram.hpp"

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
  diag_factory.register_product("VertContractDiag",&create_atmosphere_diagnostic<VertContractDiag>);
  diag_factory.register_product("VertDerivativeDiag",&create_atmosphere_diagnostic<VertDerivativeDiag>);
  diag_factory.register_product("ZonalAvgDiag",&create_atmosphere_diagnostic<ZonalAvgDiag>);
  diag_factory.register_product("ConditionalSampling",&create_atmosphere_diagnostic<ConditionalSampling>);
  diag_factory.register_product("BinaryOpsDiag", &create_atmosphere_diagnostic<BinaryOpsDiag>);
  diag_factory.register_product("HistogramDiag",&create_atmosphere_diagnostic<HistogramDiag>);
}

} // namespace scream

#endif // SCREAM_REGISTER_DIAGNOSTICS_HPP
