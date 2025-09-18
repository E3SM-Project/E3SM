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

namespace OMEGA {

// Earth constants
constexpr Real Gravity = 9.80616; // Acceleration due to gravity ~ m/s^2
constexpr Real Pi      = 3.14159265358979323846; // Pi
constexpr Real TwoPi   = 2.0 * Pi;               // 2*Pi
constexpr Real CDay    = 86400.0; // Seconds in a calendar day ~ sec
constexpr Real SDay    = 86164.0; // Seconds in a sidereal day ~ sec
constexpr Real Omega =
    2.0 * Pi / SDay;               // Angular velocity of the Earth ~ rad/sec
constexpr Real REarth = 6.37122e6; // Mean radius of the Earth ~ m

/// Physical constants
constexpr Real TkTrip  = 273.16;      // Triple point of fresh water ~ K
constexpr Real TkFrz   = 273.15;      // Freezing point of fresh water ~ K
constexpr Real TkFrzSw = TkFrz - 1.8; // Freezing point of seawater ~ K
constexpr Real RhoAir  = 1.2;         // Density of air ~ kg/m^3
constexpr Real RhoFw   = 1.000e3;     // Density of fresh water ~ kg/m^3
constexpr Real RhoSw   = 1.026e3;     // Density of seawater ~ kg/m^3
constexpr Real RhoIce  = 0.917e3;     // Density of ice ~ kg/m^3
constexpr Real CpAir   = 1.005e3; // Specific heat capacity of air ~ J/(kg*K)
constexpr Real CpFw =
    4.188e3; // Specific heat capacity of fresh water ~ J/(kg*K)
constexpr Real CpSw  = 3.996e3; // Specific heat capacity of seawater ~ J/(kg*K)
constexpr Real CpIce = 2.108e3; // Specific heat capacity of ice ~ J/(kg*K)
constexpr Real LatIce    = 3.337e5; // Latent heat of fusion ~ J/kg
constexpr Real LatVap    = 2.501e6; // Latent heat of vaporization ~ J/kg
constexpr Real LatSub    = LatIce + LatVap; // Latent heat of sublimation ~ J/kg
constexpr Real CondIce   = 2.1;      // Universal gas constant ~ J/(mol*K)
constexpr Real OcnRefSal = 34.7;     // Reference ocean salinity ~ psu
constexpr Real IceRefSal = 4.0;      // Reference ice salinity ~ psu
constexpr Real Sound     = 1.5e2;    // Speed of sound ~ m/s
constexpr Real VonKar    = 0.4;      // Von Karman constant ~ dimensionless
constexpr Real Emiss     = 1.0;      // Emissivity ~ dimensionless
constexpr Real AtmRefP   = 101325.0; // Reference atmospheric pressure ~ Pa

// Conversion factors
constexpr Real Sec2Day    = 1.0 / 86400.0; // Seconds to days
constexpr Real Day2Sec    = 86400.0;       // Days to seconds
constexpr Real Salt2PPt   = 1000.0;    // Salinity (kg/kg) to parts per thousand
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
