#ifndef SCREAM_IO_UTILS_HPP
#define SCREAM_IO_UTILS_HPP

#include "scream_io_control.hpp"
#include "share/util/scream_time_stamp.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <string>

namespace scream
{

enum class FileType {
  ModelOutput,
  ModelRestart,
  HistoryRestart,
  Unset
};

inline std::string e2str(const FileType avg) {
  using FT = FileType;
  switch (avg) {
    case FT::ModelOutput:     return "model-output";
    case FT::ModelRestart:    return "model-restart";
    case FT::HistoryRestart:  return "history-restart";
    default:                  return "UNSET";
  }
}

enum class OutputAvgType {
  Instant,
  Max,
  Min,
  Average,
  Invalid
};

inline std::string e2str(const OutputAvgType avg) {
  using OAT = OutputAvgType;
  switch (avg) {
    case OAT::Instant:  return "INSTANT";
    case OAT::Max:      return "MAX";
    case OAT::Min:      return "MIN";
    case OAT::Average:  return "AVERAGE";
    default:            return "INVALID";
  }
}

inline OutputAvgType str2avg (const std::string& s) {
  auto s_ci = ekat::upper_case(s);
  using OAT = OutputAvgType;
  for (auto e : {OAT::Instant, OAT::Max, OAT::Min, OAT::Average}) {
    if (s_ci==e2str(e)) {
      return e;
    }
  }

  return OAT::Invalid;
}

// The AD will pass a default constructed control, since it doesn't know the values
// of REST_N/REST_OPTION used in the previous run
// Output streams MUST pass a valid control structure, cause we need to differentiate
// between, e.g., streams with same filename prefix, but different output freq specs
std::string find_filename_in_rpointer (
    const std::string& filename_prefix,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0,
    const OutputAvgType avg_type = OutputAvgType::Instant,
    const IOControl& control = {});

struct DefaultMetadata {

  std::string get_longname (const std::string& name) {
    if (name_2_longname.count(name)>0) {
      return name_2_longname.at(name);
    } else {
      // TODO: Do we want to print a Warning message?  I'm not sure if its needed.
      return name;
    }
  }

  std::string get_standardname (const std::string& name) {
    if (name_2_standardname.count(name)>0) {
      return name_2_standardname.at(name);
    } else {
      // TODO: Do we want to print a Warning message?  I'm not sure if its needed.
      return name;
    }
  }

  // Create map of longnames, can be added to as developers see fit.
  std::map<std::string,std::string> name_2_longname = {
    {"lev","hybrid level at midpoints (1000*(A+B))"},
    {"ilev","hybrid level at interfaces (1000*(A+B))"},
    {"hyai","hybrid A coefficient at layer interfaces"},
    {"hybi","hybrid B coefficient at layer interfaces"},
    {"hyam","hybrid A coefficient at layer midpoints"},
    {"hybm","hybrid B coefficient at layer midpoints"}
  };

  // Create map of longnames, can be added to as developers see fit.
  std::map<std::string,std::string> name_2_standardname = {
    {"p_mid"                                                       , "air_pressure"},
    {"p_mid_at_cldtop"                                             , "air_pressure_at_cloud_top"},
    {"T_2m"                                                        , "air_temperature"},
    {"T_mid"                                                       , "air_temperature"},
    {"T_mid_at_cldtop"                                             , "air_temperature_at_cloud_top"},
    {"aero_g_sw"                                                   , "asymmetry_factor_of_ambient_aerosol_particles"},
    {"pbl_height"                                                  , "atmosphere_boundary_layer_thickness"},
    {"precip_liq_surf_mass"                                        , "atmosphere_mass_content_of_liquid_precipitation"},
    {"cldlow"                                                      , "low_type_cloud_area_fraction"},
    {"cldmed"                                                      , "medium_type_cloud_area_fraction"},
    {"cldhgh"                                                      , "high_type_cloud_area_fraction"},
    {"cldtot"                                                      , "cloud_area_fraction"},
    {"cldfrac_tot_at_cldtop"                                       , "cloud_area_fraction"},
    {"cldfrac_tot"                                                 , "cloud_area_fraction_in_atmosphere_layer"},
    {"cldfrac_tot_for_analysis"                                    , "cloud_area_fraction_in_atmosphere_layer"},
    {"cldfrac_rad"                                                 , "cloud_area_fraction_in_atmosphere_layer"},
    {"qi"                                                          , "cloud_ice_mixing_ratio"},
    {"qc"                                                          , "cloud_liquid_water_mixing_ratio"},
    {"U"                                                           , "eastward_wind"},
    {"eff_radius_qi"                                               , "effective_radius_of_cloud_ice_particles"},
    {"eff_radius_qc"                                               , "effective_radius_of_cloud_liquid_water_particles"},
    {"eff_radius_qc_at_cldtop"                                     , "effective_radius_of_cloud_liquid_water_particles_at_liquid_water_cloud_top"},
    {"eff_radius_qr"                                               , "effective_radius_of_cloud_rain_particles"},
    {"qv"                                                          , "humidity_mixing_ratio"},
    {"cldfrac_ice_at_cldtop"                                       , "ice_cloud_area_fraction"},
    {"cldfrac_ice"                                                 , "ice_cloud_area_fraction_in_atmosphere_layer"},
    {"omega"                                                       , "lagrangian_tendency_of_air_pressure"},
    {"landfrac"                                                    , "land_area_fraction"},
    {"latitude"                                                    , "latitude"},
    {"cldfrac_liq_at_cldtop"                                       , "liquid_water_cloud_area_fraction"},
    {"cldfrac_liq"                                                 , "liquid_water_cloud_area_fraction_in_atmosphere_layer"},
    {"longitude"                                                   , "longitude"},
    {"rainfrac"                                                    , "mass_fraction_of_liquid_precipitation_in_air"},
    {"V"                                                           , "northward_wind"},
    {"nc"                                                          , "number_concentration_of_cloud_liquid_water_particles_in_air"},
    {"cdnc_at_cldtop"                                              , "number_concentration_of_cloud_liquid_water_particles_in_air_at_liquid_water_cloud_top"},
    {"ni"                                                          , "number_concentration_of_ice_crystals_in_air"},
    {"aero_tau_sw"                                                 , "optical_thickness_of_atmosphere_layer_due_to_ambient_aerosol_particles"},
    {"aero_tau_lw"                                                 , "optical_thickness_of_atmosphere_layer_due_to_ambient_aerosol_particles"},
    {"aero_ssa_sw"                                                 , "single_scattering_albedo_in_air_due_to_ambient_aerosol_particles"},
    {"sunlit"                                                      , "sunlit_binary_mask"},
    {"ps"                                                          , "surface_air_pressure"},
    {"LW_flux_dn_at_model_bot"                                     , "surface_downwelling_longwave_flux_in_air"},
    {"SW_flux_dn_at_model_bot"                                     , "surface_downwelling_shortwave_flux_in_air"},
    {"SW_clrsky_flux_dn_at_model_bot"                              , "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky"},
    {"phis"                                                        , "surface_geopotential"},
    {"surf_radiative_T"                                            , "surface_temperature"},
    {"surf_sens_flux"                                              , "surface_upward_sensible_heat_flux"},
    {"SW_flux_dn_at_model_top"                                     , "toa_incoming_shortwave_flux"},
    {"LW_flux_up_at_model_top"                                     , "toa_outgoing_longwave_flux"},
    {"LW_clrsky_flux_up_at_model_top"                              , "toa_outgoing_longwave_flux_assuming_clear_sky"},
    {"surf_evap"                                                   , "water_evapotranspiration_flux"},
    {"AtmosphereDensity"                                           , "air_density"},
    {"PotentialTemperature"                                        , "air_potential_temperature"},
    {"SeaLevelPressure"                                            , "air_pressure_at_mean_sea_level"},
    {"IceWaterPath"                                                , "atmosphere_mass_content_of_cloud_ice"},
    {"LiqWaterPath"                                                , "atmosphere_mass_content_of_cloud_liquid_water"},
    {"VapWaterPath"                                                , "atmosphere_mass_content_of_water_vapor"},
    {"AerosolOpticalDepth550nm"                                    , "atmosphere_optical_thickness_due_to_ambient_aerosol_particles"},
    {"Exner"                                                       , "dimensionless_exner_function"},
    {"z_mid"                                                       , "geopotential_height"},
    {"geopotential_mid"                                            , "geopotential_height"},
    {"RelativeHumidity"                                            , "relative_humidity"},
    {"surface_upward_latent_heat_flux"                             , "surface_upward_latent_heat_flux"},
    {"LongwaveCloudForcing"                                        , "toa_longwave_cloud_radiative_effect"},
    {"ShortwaveCloudForcing"                                       , "toa_shortwave_cloud_radiative_effect"},
    {"VirtualTemperature"                                          , "virtual_temperature"},
    {"VaporFlux"                                                   , "water_evapotranspiration_flux"},
    {"wind_speed"                                                  , "wind_speed"}
  };
  
};

// Shortcut to write/read to/from YYYYMMDD/HHMMSS attributes in the NC file
void write_timestamp (const std::string& filename, const std::string& ts_name,
                      const util::TimeStamp& ts, const bool write_nsteps = false);
util::TimeStamp read_timestamp (const std::string& filename,
                                const std::string& ts_name,
                                const bool read_nsteps = false);

} // namespace scream
#endif // SCREAM_IO_UTILS_HPP
