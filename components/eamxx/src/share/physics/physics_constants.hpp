#ifndef EAMXX_PHYSICS_CONSTANTS_HPP
#define EAMXX_PHYSICS_CONSTANTS_HPP

#include "share/util/eamxx_quantity.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_string_utils.hpp>
#include <ekat_math_utils.hpp>
#include <ekat_units.hpp>

namespace scream {
namespace physics {

/*
 * Physical/math constants used in EAMxx.
 *
 * Note that a potential optimization could be to change the type of
 * Scalar constants that have integer values to int.
 */

template <typename Scalar>
struct Constants
{
  using quantity_t = Quantity<Scalar>;
  using ci_string  = ekat::CaseInsensitiveString;
  using Units = ekat::units::Units;

  // Some predefined units
  static constexpr auto m   = ekat::units::m;
  static constexpr auto s   = ekat::units::s;
  static constexpr auto mol = ekat::units::mol;
  static constexpr auto kg  = ekat::units::kg;
  static constexpr auto g   = ekat::units::g;
  static constexpr auto Pa  = ekat::units::Pa;
  static constexpr auto K   = ekat::units::K;
  static constexpr auto J   = ekat::units::J;
  static constexpr auto W   = ekat::units::W;

  // Some shorthand units
  static constexpr auto m2     = Units(m*m,"m2");
  static constexpr auto m3     = Units(m*m*m,"m3");
  static constexpr auto s2     = Units(s*s,"s2");
  static constexpr auto nondim = Units::nondimensional();
  static constexpr auto rad    = Units(nondim,"rad");

  // Commonly used numeric constants
  static constexpr Scalar ZERO    = 0.0;
  static constexpr Scalar ONE     = 1;
  static constexpr Scalar THIRD   = 1.0/3.0;
  static constexpr Scalar SXTH    = 1.0/6.0;
  static constexpr Scalar Pi      = 3.14159265358979323;
  static constexpr Scalar PIOV3   = Pi*THIRD;
  static constexpr Scalar PIOV6   = Pi*SXTH;
  static constexpr Scalar macheps = std::numeric_limits<Real>::epsilon();

  // Physical constants
  static constexpr auto Cpair           = quantity_t(1004.64, J/kg/K);
  static constexpr auto CP              = Cpair;          // heat constant of air at constant pressure, J/kg
  static constexpr auto INV_CP          = ONE/CP;
  static constexpr auto Rair            = quantity_t(287.042, J/kg/K); // gas constant for dry air
  static constexpr auto RD              = Rair;
  static constexpr auto RH2O            = quantity_t(461.505, J/kg/K); // Water vapor gas constant
  static constexpr auto RV              = RH2O;

  static constexpr auto RHO_H2O         = quantity_t(1000.0,kg/m3);
  static constexpr auto INV_RHO_H2O     = 1/RHO_H2O;
  static constexpr auto RHOW            = RHO_H2O;
  static constexpr auto INV_RHOW        = 1/RHOW;
  static constexpr auto RhoIce          = quantity_t(917.0,kg/m3);     // Ice density at 0 C from Wallace+Hobbes 1977

  static constexpr auto MWH2O           = quantity_t(18.016, g/mol);   // water molar mass
  static constexpr auto MWWV            = MWH2O;
  static constexpr auto MWdry           = quantity_t(28.966, g/mol);   // dry air molar mass
  static constexpr auto ep_2            = MWH2O / MWdry;               // dimensionless (MWH2O/MWdry)
  static constexpr auto o2mmr           = quantity_t(0.23143,nondim);  // o2 mass mixing ratio (dimensionless)

  static constexpr auto gravit          = quantity_t(9.80616, m/s2);

  static constexpr auto LatVap          = quantity_t(2501000.0, J/kg);
  static constexpr auto LatIce          = quantity_t(333700.0, J/kg);
  static constexpr auto CpLiq           = quantity_t(4188.0, J/kg/K);

  static constexpr auto Tmelt           = quantity_t(273.15, K);
  static constexpr auto T_zerodegc      = Tmelt;
  static constexpr auto T_homogfrz      = quantity_t(Tmelt.value - 40,K);
  static constexpr auto T_rainfrz       = quantity_t(Tmelt.value - 4,K);

  static constexpr Scalar RHO_RIMEMIN   =  50.0;  //Min limit for rime density [kg m-3]
  static constexpr Scalar RHO_RIMEMAX   = 900.0;  //Max limit for rime density [kg m-3]
  static constexpr Scalar INV_RHO_RIMEMAX  =  1.0/RHO_RIMEMAX; // Inverse for limits for rime density [kg m-3]
  static constexpr Scalar BIMM          = 2.0;
  static constexpr Scalar CONS1         = PIOV6*RHOW.value;
  static constexpr Scalar CONS2         = 4.*PIOV3*RHOW.value;
  static constexpr Scalar CONS3         =  1.0/(CONS2*1.562500000000000e-14); // 1./(CONS2*pow(25.e-6,3.0));
  static constexpr Scalar CONS5         = PIOV6*BIMM;
  static constexpr Scalar CONS6         = PIOV6*PIOV6*RHOW.value*BIMM;
  static constexpr Scalar CONS7         = 4.*PIOV3*RHOW.value*1.e-18;
  static constexpr Scalar QSMALL        = 1.e-14;
  static constexpr Scalar QTENDSMALL    = 1e-20;
  static constexpr Scalar BSMALL        = 1.e-15;
  static constexpr Scalar NSMALL        = 1.e-16;
  static constexpr auto P0              = quantity_t(100000.0,Pa);        // reference pressure, Pa
  static constexpr auto RHOSUR          = P0/(RD*Tmelt);
  static constexpr Scalar rhosui        = 60000/(RD.value*253.15);
  static constexpr auto RHO_1000MB      = P0/(RD*Tmelt);
  static constexpr Scalar RHO_600MB     = 60000/(RD.value*253.15);
  //  static constexpr Scalar Tol           = ekat::is_single_precision<Real>::value ? 2e-5 : 1e-14;
  static constexpr Scalar dt_left_tol   = 1.e-4;
  static constexpr Scalar bcn           = 2.;
  static constexpr Scalar dropmass      = 5.2e-7;
  static constexpr Scalar NCCNST        = 200.0e+6;
  static constexpr Scalar incloud_limit = 5.1e-3;
  static constexpr Scalar precip_limit  = 1.0e-2;
  static constexpr Scalar Karman        = 0.4;
  static constexpr auto Avogad          = quantity_t (6.02214e26,1/mol);
  static constexpr auto Boltz           = quantity_t (1.38065e-23,J/K);
  static constexpr auto Rgas            = Avogad * Boltz;
  static constexpr auto RWV             = Rgas / MWWV;
  static constexpr Scalar ZVIR          = (RWV.value / Rair.value) - 1.0;
  static constexpr Scalar f1r           = 0.78;
  static constexpr Scalar f2r           = 0.32;
  static constexpr Scalar nmltratio     = 1.0; // ratio of rain number produced to ice number loss from melting
  static constexpr Scalar basetemp      = 300.0;
  static constexpr auto r_earth         = quantity_t(6.376e6,m); // Radius of the earth in m
  static constexpr auto stebol          = quantity_t(5.670374419e-8,W/m2/pow(K,4)); // Stefan-Boltzmann's constant (W/m^2/K^4)
  static constexpr auto omega           = quantity_t(7.292e-5,rad/s); // Earth's rotation (rad/sec)

  // Table dimension constants
  static constexpr int VTABLE_DIM0    = 300;
  static constexpr int VTABLE_DIM1    = 10;
  static constexpr int MU_R_TABLE_DIM = 150;

  // Turbulent Mountain Stress constants
  static constexpr Scalar orocnst     = 1;     // Converts from standard deviation to height [ no unit ]
  static constexpr Scalar z0fac       = 0.075; // Factor determining z_0 from orographic standard deviation [ no unit ]

  // switch for warm-rain parameterization
  // = 1 Seifert and Beheng 2001
  // = 2 Beheng 1994
  // = 3 Khairoutdinov and Kogan 2000
  static constexpr int IPARAM         = 3;

  // Gases
  static Scalar get_gas_mol_weight(ci_string gas_name);

  // For use in converting area to length for a column cell
  // World Geodetic System 1984 (WGS84)
  static constexpr Scalar earth_ellipsoid1 = 111132.92; // first coefficient, meters per degree longitude at equator
  static constexpr Scalar earth_ellipsoid2 = 559.82;    // second expansion coefficient for WGS84 ellipsoid
  static constexpr Scalar earth_ellipsoid3 = 1.175;     // third expansion coefficient for WGS84 ellipsoid
};

// Gases
// Define the molecular weight for each gas, which can then be
// used to determine the volume mixing ratio for each gas.
template <typename Scalar>
Scalar Constants<Scalar>::get_gas_mol_weight(ci_string gas_name) {
  //TODO: Possible improvement would be to design a device friendly function
  if        (gas_name == "h2o") {
    return MWH2O.value;
  } else if (gas_name == "co2") {
    return 44.0095;
  } else if (gas_name == "o3" ) {
    return 47.9982;
  } else if (gas_name == "n2o") {
    return 44.0128;
  } else if (gas_name == "co" ) {
    return 28.0101;
  } else if (gas_name == "ch4") {
    return 16.04246;
  } else if (gas_name == "o2" ) {
    return 31.998;
  } else if (gas_name == "n2" ) {
    return 28.0134;
  } else if (gas_name == "cfc11" ) {
    return 136.;
  } else if (gas_name == "cfc12" ) {
    return 120.;
  }
  return ekat::invalid<Scalar>();
}

template <typename Scalar>
constexpr Scalar Constants<Scalar>::NSMALL;

template <typename Scalar>
constexpr Scalar Constants<Scalar>::QSMALL;

template <typename Scalar>
constexpr Scalar Constants<Scalar>::QTENDSMALL;

} // namespace physics
} // namespace scream

#endif // EAMXX_PHYSICS_CONSTANTS_HPP
