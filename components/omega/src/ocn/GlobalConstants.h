#ifndef OMEGA_GLOBALCONSTANTS_H
#define OMEGA_GLOBALCONSTANTS_H
//===-- ocn/GlobalConstants.h - Global Constants --------------------*- C++
//-*-===//
//
/// \file
/// \brief Contains global constants for use across Omega
///
/// This header defines constants to be used across Omega (to be replaced by a
/// shared E3SM constants file when that becomes available)
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "pcd_const.h"
#include <cmath>

namespace OMEGA {

// Mathematical constants
constexpr Real Pi = pcd::pi; // Pi (from Physical Constants Dictionary)
constexpr Real TwoPi =
    2.0 * pcd::pi; // 2*Pi (from Physical Constants Dictionary)
constexpr Real Radian =
    pcd::radian; // Radians in one degree (from Physical Constants Dictionary)
constexpr Real Degree =
    pcd::degree; // Degrees in one radian (from Physical Constants Dictionary)
constexpr Real SqrtTwo =
    pcd::square_root_of_2; // Square root of 2 (from Physical Constants
                           // Dictionary)
constexpr Real SqrtThree =
    pcd::square_root_of_3; // Square root of 3 (from Physical Constants
                           // Dictionary)

// Earth constants
constexpr Real Gravity =
    pcd::standard_acceleration_of_gravity; // Acceleration due to gravity ~
                                           // m/s^2 (from Physical Constants
                                           // Dictionary)
constexpr Real CDay = 86400.0;             // Seconds in a calendar day ~ sec
constexpr Real Omega =
    pcd::angular_velocity; // Angular velocity of the Earth ~ rad/sec (from
                           // Physical Constants Dictionary)
constexpr Real SDay   = TwoPi / Omega;    // Seconds in a sidereal day ~ sec
constexpr Real REarth = pcd::mean_radius; // Mean radius of the Earth ~ m (from
                                          // Physical Constants Dictionary)

/// Physical constants
constexpr Real TkTrip =
    pcd::water_triple_point_temperature; // Triple point of fresh water ~ K
                                         // (from Physical Constants Dictionary)
constexpr Real TkFrz = pcd::pure_water_freezing_temperature_reference;
// Freezing point of fresh water ~ K (from Physical Constants Dictionary)
constexpr Real TkFrzSw = TkFrz - 1.8; // Freezing point of seawater ~ K
constexpr Real RhoAir =
    pcd::dry_air_density_at_standard_temperature_and_pressure;
// Density of air ~ kg/m^3 (from Physical Constants Dictionary)
constexpr Real RhoFw = pcd::pure_water_density_reference;
// Density of fresh water ~ kg/m^3 (from Physical Constants Dictionary)
constexpr Real RhoSw = pcd::seawater_density_reference;
// Density of seawater ~ kg/m^3 (from Physical Constants Dictionary)
constexpr Real RhoIce = pcd::sea_ice_density_reference;
// Density of ice ~ kg/m^3 (from Physical Constants Dictionary)
constexpr Real CpAir = pcd::dry_air_specific_heat_capacity_reference;
// Specific heat capacity of air ~ J/(kg*K) (from Physical Constants Dictionary)
constexpr Real CpFw = pcd::pure_water_specific_heat_capacity_reference;
// Specific heat capacity of fresh water ~ J/(kg*K) (from Physical Constants
// Dictionary)
constexpr Real CpSw = pcd::seawater_specific_heat_capacity_reference;
// Specific heat capacity of seawater ~ J/(kg*K) (from Physical Constants
// Dictionary)
constexpr Real Cp0Sw =
    3991.86795711963; // Specific heat capacity of seawater for TEOS-10
constexpr Real CpIce = pcd::sea_ice_specific_heat_capacity_reference;
// Specific heat capacity of ice ~ J/(kg*K) (from Physical Constants Dictionary)
constexpr Real LatIce = pcd::latent_heat_of_fusion_reference;
// Latent heat of fusion ~ J/kg (from Physical Constants Dictionary)
constexpr Real LatVap = pcd::latent_heat_of_vaporization_reference;
// Latent heat of vaporization ~ J/kg (from Physical Constants Dictionary)
constexpr Real LatSub = pcd::latent_heat_of_sublimation_reference;
// Latent heat of sublimation ~ J/kg (from Physical Constants Dictionary)
constexpr Real CondIce = pcd::sea_ice_thermal_conductivity_reference;
// Thermal conductivity of ice ~ W/(m*K) (from Physical Constants Dictionary)
constexpr Real OcnRefSal = pcd::ocean_reference_salinity;
// Reference ocean salinity ~ g/kg (from Physical Constants Dictionary)
constexpr Real IceRefSal = pcd::sea_ice_reference_salinity;
// Reference ice salinity ~ g/kg (from Physical Constants Dictionary)
constexpr Real Sound = pcd::speed_of_sound_in_seawater_reference;
// Speed of sound in seawater ~ m/s (from Physical Constants Dictionary)
constexpr Real VonKar = 0.4; // Von Karman constant ~ dimensionless
constexpr Real Emiss  = 1.0; // Emissivity ~ dimensionless
constexpr Real AtmRefP =
    pcd::standard_atmosphere; // Reference atmospheric pressure ~ Pa (from
                              // Physical Constants Dictionary)

// Conversion factors
constexpr Real Sec2Day = 1.0 / 86400.0; // Seconds to days
constexpr Real Day2Sec = 86400.0;       // Days to seconds
constexpr Real Rad2Deg =
    pcd::radian; // Radians to degrees (from Physical Constants Dictionary)
constexpr Real Deg2Rad =
    pcd::degree; // Degrees to radians (from Physical Constants Dictionary)
constexpr Real Salt2PPt = 1000.0;   // Salinity (kg/kg) to parts per thousand
constexpr Real SS0      = 35.16504; // Standard ocean reference salinity (g/kg)
constexpr Real Psu2Gpkg =
    SS0 / 35.0; // Unit conversion factor for salinities (psu -> g/kg)
constexpr Real Sfac       = 0.0248826675584615; // 1 / (40 * Psu2Gpkg)
constexpr Real PPt2Salt   = 1.0e-3;    // Parts per thousand to salinity (kg/kg)
constexpr Real Mass2Sv    = 1.0e-12;   // Mass flux (kg/s) to Sverdrup
constexpr Real Heat2Pw    = 4.186e-15; // Heat flux (W) to Petawatt
constexpr Real Salt2SvPpt = 1.0e-9;    // Salt flux (kg/s) to Sv*ppt
constexpr Real Salt2MmDay = 3.1536e+5; // Salt flux to water flux (mm/day)
constexpr Real Db2Pa      = 1.0e4;     // Decibar to Pascal
constexpr Real Pa2Db      = 1.0e-4;    // Pascal to Decibar
constexpr Real Cm2M       = 1.0e-2;    // Centimeters to meters
constexpr Real M2Cm       = 1.0e2;     // Meters to centimeters
constexpr Real HFluxFac =
    1.0 / (RhoSw * CpSw);         // Heat flux (W/m^2) to temp flux (C*m/s)
constexpr Real FwFluxFac = 1.e-6; // Fw flux (kg/m^2/s) to salt((msu/psu)*m/s)
constexpr Real SaltFac =
    -OcnRefSal * FwFluxFac;    // Fw flux (kg/m^2/s) to salt flux (msu*m/s)
constexpr Real SFluxFac = 1.0; // Salt flux (kg/m^2/s) to salt flux (msu*m/s)

} // namespace OMEGA
#endif
