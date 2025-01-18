#ifndef PHYSICS_CONSTANTS_HPP
#define PHYSICS_CONSTANTS_HPP

#include "share/scream_types.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/logging/ekat_logger.hpp"

#include <vector>

namespace scream {
namespace physics {

/*
 * Mathematical constants used by atmosphere processes.
 *
 * Note that a potential optimization could be to change the type of
 * Scalar constants that have integer values to int.
 */

template <typename Scalar>
struct Constants
{
  using ci_string      = ekat::CaseInsensitiveString;

  static constexpr Scalar Cpair         = 1004.64;
  static constexpr Scalar Rair          = 287.042;
  static constexpr Scalar RH2O          = 461.505;
  static constexpr Scalar RV            = RH2O;         // Water vapor gas constant ~ J/K/kg     !461.51
  static constexpr Scalar RHO_H2O       = 1000.0;
  static constexpr Scalar INV_RHO_H2O   = 1.0/RHO_H2O;
  static constexpr Scalar RhoIce        = 917.0;        // Ice density at 0 C from Wallace+Hobbes 1977
  static constexpr Scalar MWH2O         = 18.016;
  static constexpr Scalar MWdry         = 28.966;
  static constexpr Scalar o2mmr         = 0.23143;      // o2 mass mixing ratio
  static constexpr Scalar ep_2          = MWH2O/MWdry;  // ratio of molecular mass of water to the molecular mass of dry air !0.622
  static constexpr Scalar gravit        = 9.80616;
  static constexpr Scalar LatVap        = 2501000.0;
  static constexpr Scalar LatIce        = 333700.0;
  static constexpr Scalar CpLiq         = 4188.0;
  static constexpr Scalar Tmelt         = 273.15;
  static constexpr Scalar T_zerodegc    = Tmelt;
  static constexpr Scalar T_homogfrz    = Tmelt - 40;
  static constexpr Scalar T_rainfrz     = Tmelt - 4;
  static constexpr Scalar Pi            = 3.14159265358979323;
  static constexpr Scalar RHOW          = RHO_H2O;
  static constexpr Scalar INV_RHOW      = 1.0/RHOW;
  static constexpr Scalar RHO_RIMEMIN   =  50.0;  //Min limit for rime density [kg m-3]
  static constexpr Scalar RHO_RIMEMAX   = 900.0;  //Max limit for rime density [kg m-3]
  static constexpr Scalar INV_RHO_RIMEMAX  =  1.0/RHO_RIMEMAX; // Inverse for limits for rime density [kg m-3]
  static constexpr Scalar THIRD         = 1.0/3.0;
  static constexpr Scalar SXTH          = 1.0/6.0;
  static constexpr Scalar PIOV3         = Pi*THIRD;
  static constexpr Scalar PIOV6         = Pi*SXTH;
  static constexpr Scalar BIMM          = 2.0;
  static constexpr Scalar CONS1         = PIOV6*RHOW;
  static constexpr Scalar CONS2         = 4.*PIOV3*RHOW;
  static constexpr Scalar CONS3         =  1.0/(CONS2*1.562500000000000e-14); // 1./(CONS2*pow(25.e-6,3.0));
  static constexpr Scalar CONS5         = PIOV6*BIMM;
  static constexpr Scalar CONS6         = PIOV6*PIOV6*RHOW*BIMM;
  static constexpr Scalar CONS7         = 4.*PIOV3*RHOW*1.e-18;
  static constexpr Scalar QSMALL        = 1.e-14;
  static constexpr Scalar QTENDSMALL    = 1e-20;
  static constexpr Scalar BSMALL        = 1.e-15;
  static constexpr Scalar NSMALL        = 1.e-16;
  static constexpr Scalar ZERO          = 0.0;
  static constexpr Scalar ONE           = 1.0;
  static constexpr Scalar P0            = 100000.0;        // reference pressure, Pa
  static constexpr Scalar RD            = Rair;          // gas constant for dry air, J/kg/K
  static constexpr Scalar RHOSUR        = P0/(RD*Tmelt);
  static constexpr Scalar rhosui        = 60000/(RD*253.15);
  static constexpr Scalar RHO_1000MB    = P0/(RD*Tmelt);
  static constexpr Scalar RHO_600MB     = 60000/(RD*253.15);
  static constexpr Scalar CP            = Cpair;          // heat constant of air at constant pressure, J/kg
  static constexpr Scalar INV_CP        = 1.0/CP;
  //  static constexpr Scalar Tol           = ekat::is_single_precision<Real>::value ? 2e-5 : 1e-14;
  static constexpr Scalar macheps = std::numeric_limits<Real>::epsilon();
  static constexpr Scalar dt_left_tol   = 1.e-4;
  static constexpr Scalar bcn           = 2.;
  static constexpr Scalar dropmass      = 5.2e-7;
  static constexpr Scalar NCCNST        = 200.0e+6;
  static constexpr Scalar incloud_limit = 5.1e-3;
  static constexpr Scalar precip_limit  = 1.0e-2;
  static constexpr Scalar Karman        = 0.4;
  static constexpr Scalar Avogad        = 6.02214e26;
  static constexpr Scalar Boltz         = 1.38065e-23;
  static constexpr Scalar Rgas          = Avogad * Boltz;
  static constexpr Scalar MWWV          = MWH2O;
  static constexpr Scalar RWV           = Rgas / MWWV;
  static constexpr Scalar ZVIR          = (RWV / Rair) - 1.0;
  static constexpr Scalar f1r           = 0.78;
  static constexpr Scalar f2r           = 0.32;
  static constexpr Scalar nmltratio     = 1.0; // ratio of rain number produced to ice number loss from melting
  static constexpr Scalar basetemp      = 300.0;
  static constexpr Scalar r_earth       = 6.376e6; // Radius of the earth in m
  static constexpr Scalar stebol        = 5.67e-8; // Stefan-Boltzmann's constant (W/m^2/K^4)
  static constexpr Scalar omega         = 7.292e-5; // Earth's rotation (rad/sec)

  // Table dimension constants
  static constexpr int VTABLE_DIM0    = 300;
  static constexpr int VTABLE_DIM1    = 10;
  static constexpr int MU_R_TABLE_DIM = 150;

  // Turbulent Mountain Stress constants
  static constexpr Scalar orocnst = 1;     // Converts from standard deviation to height [ no unit ]
  static constexpr Scalar z0fac   = 0.075; // Factor determining z_0 from orographic standard deviation [ no unit ]

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
    return Scalar(Constants<Scalar>::MWH2O);
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
  return ekat::ScalarTraits<Scalar>::invalid();
}

template <typename Scalar>
constexpr Scalar Constants<Scalar>::NSMALL;

template <typename Scalar>
constexpr Scalar Constants<Scalar>::QSMALL;

template <typename Scalar>
constexpr Scalar Constants<Scalar>::QTENDSMALL;

template<typename Scalar>
constexpr Scalar Constants<Scalar>::ZERO;

template <typename Scalar>
constexpr Scalar Constants<Scalar>::Tmelt;

template <typename Scalar>
constexpr Scalar Constants<Scalar>::macheps;

} // namespace physics
} // namespace scream

#endif
