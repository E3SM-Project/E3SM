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
#include "field_prev.hpp"
#include "field_over_dt.hpp"
#include "horiz_avg.hpp"
#include "vert_contract.hpp"
#include "vert_derivative.hpp"
#include "zonal_avg.hpp"
#include "conditional_sampling.hpp"
#include "binary_op.hpp"
#include "histogram.hpp"

namespace scream {

inline void register_diagnostics () {
  auto& diag_factory = DiagnosticFactory::instance();
  diag_factory.register_product("PotentialTemperature",&create_diagnostic<PotentialTemperature>);
  diag_factory.register_product("FieldAtLevel",&create_diagnostic<FieldAtLevel>);
  diag_factory.register_product("FieldAtHeight",&create_diagnostic<FieldAtHeight>);
  diag_factory.register_product("FieldAtPressureLevel",&create_diagnostic<FieldAtPressureLevel>);
  diag_factory.register_product("AtmosphereDensity",&create_diagnostic<AtmDensity>);
  diag_factory.register_product("Exner",&create_diagnostic<Exner>);
  diag_factory.register_product("VirtualTemperature",&create_diagnostic<VirtualTemperature>);
  diag_factory.register_product("DryStaticEnergy",&create_diagnostic<DryStaticEnergy>);
  diag_factory.register_product("SeaLevelPressure",&create_diagnostic<SeaLevelPressure>);
  diag_factory.register_product("WaterPath",&create_diagnostic<WaterPath>);
  diag_factory.register_product("ShortwaveCloudForcing",&create_diagnostic<ShortwaveCloudForcing>);
  diag_factory.register_product("LongwaveCloudForcing",&create_diagnostic<LongwaveCloudForcing>);
  diag_factory.register_product("RelativeHumidity",&create_diagnostic<RelativeHumidity>);
  diag_factory.register_product("VaporFlux",&create_diagnostic<VaporFlux>);
  diag_factory.register_product("VerticalLayer",&create_diagnostic<VerticalLayer>);
  diag_factory.register_product("precip_surf_mass_flux",&create_diagnostic<PrecipSurfMassFlux>);
  diag_factory.register_product("surface_upward_latent_heat_flux",&create_diagnostic<SurfaceUpwardLatentHeatFlux>);
  diag_factory.register_product("wind_speed",&create_diagnostic<WindSpeed>);
  diag_factory.register_product("AerosolOpticalDepth550nm",&create_diagnostic<AODVis>);
  diag_factory.register_product("NumberPath",&create_diagnostic<NumberPath>);
  diag_factory.register_product("AeroComCld",&create_diagnostic<AeroComCld>);
  diag_factory.register_product("FieldPrev",&create_diagnostic<FieldPrev>);
  diag_factory.register_product("FieldOverDt",&create_diagnostic<FieldOverDt>);
  diag_factory.register_product("HorizAvg",&create_diagnostic<HorizAvg>);
  diag_factory.register_product("VertContract",&create_diagnostic<VertContract>);
  diag_factory.register_product("VertDerivative",&create_diagnostic<VertDerivative>);
  diag_factory.register_product("ZonalAvg",&create_diagnostic<ZonalAvg>);
  diag_factory.register_product("ConditionalSampling",&create_diagnostic<ConditionalSampling>);
  diag_factory.register_product("BinaryOp", &create_diagnostic<BinaryOp>);
  diag_factory.register_product("Histogram",&create_diagnostic<Histogram>);
}

} // namespace scream

#endif // SCREAM_REGISTER_DIAGNOSTICS_HPP
