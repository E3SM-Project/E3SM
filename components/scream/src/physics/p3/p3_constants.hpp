#ifndef P3_CONSTANTS_HPP
#define P3_CONSTANTS_HPP

#include "share/scream_types.hpp"

#include <vector>

namespace scream {
namespace p3 {

/*
 * Mathematical constants used by p3.
 */

template <typename Scalar>
struct Constants
{
  static constexpr Scalar Cpair       = 1004.64;
  static constexpr Scalar Rair        = 287.042;
  static constexpr Scalar RH2O        = 461.505;
  static constexpr Scalar RhoH2O      = 1000.0;
  static constexpr Scalar MWH2O       = 18.016;
  static constexpr Scalar MWdry       = 28.966;
  static constexpr Scalar gravit      = 9.80616;
  static constexpr Scalar LatVap      = 2501000.0;
  static constexpr Scalar LatIce      = 333700.0;
  static constexpr Scalar CpLiq       = 4188.0;
  static constexpr Scalar Tmelt       = 273.15;
  static constexpr Scalar Pi          = 3.14159265;
  static constexpr long long int    iulog       = 98;
  static constexpr bool   masterproc  = true;
  static constexpr Scalar RHOW        = RhoH2O;
  static constexpr Scalar INV_RHOW    = 1.0/RHOW;
  static constexpr Scalar THIRD       = 1.0/3.0;
  static constexpr Scalar SXTH        = 1.0/6.0;
  static constexpr Scalar PIOV6       = Pi*SXTH;
  static constexpr Scalar CONS1       = PIOV6*RHOW;
  static constexpr Scalar QSMALL      = 1.e-14;
  static constexpr Scalar NSMALL      = 1.e-16;
  static constexpr Scalar P0          = 100000.0;        // reference pressure, Pa
  static constexpr Scalar RD          = 287.15;          // gas constant for dry air, J/kg/K
  static constexpr Scalar RHOSUR      = P0/(RD*273.15);
  static constexpr Scalar CP          = Cpair;          // heat constant of air at constant pressure, J/kg
  static constexpr Scalar INV_CP      = 1.0/CP;

  // Constants for ice lookup tables
  static constexpr int    DENSIZE     = 5;
  static constexpr int    RIMSIZE     = 4;
  static constexpr int    ISIZE       = 50;
  static constexpr int    TABSIZE     = 12; // number of quantities used from lookup table
  static constexpr int    RCOLLSIZE   = 30;
  static constexpr int    COLTABSIZE  = 2;  // number of ice-rain collection  quantities used from lookup table

  static constexpr Scalar LOOKUP_TABLE_1A_DUM1_C = 1.0/(0.1*std::log10(261.7));

  static constexpr const char* P3_VERSION = "2.8.2";
};

template <typename Scalar>
constexpr Scalar Constants<Scalar>::NSMALL;

template <typename Scalar>
constexpr int Constants<Scalar>::ISIZE;

template <typename Scalar>
using vector_2d_t = std::vector<std::vector<Scalar> >;

template <typename Scalar>
struct Globals
{
  static constexpr int VTABLE_DIM0 = 300;
  static constexpr int VTABLE_DIM1 = 10;
  static constexpr int MU_R_TABLE_DIM = 150;

  static vector_2d_t<Scalar> VN_TABLE, VM_TABLE;
  static std::vector<Scalar> MU_R_TABLE;
};

template <typename Scalar>
vector_2d_t<Scalar> Globals<Scalar>::VN_TABLE(VTABLE_DIM0, std::vector<Scalar>(VTABLE_DIM1));

template <typename Scalar>
vector_2d_t<Scalar> Globals<Scalar>::VM_TABLE(VTABLE_DIM0, std::vector<Scalar>(VTABLE_DIM1));

template <typename Scalar>
std::vector<Scalar> Globals<Scalar>::MU_R_TABLE(MU_R_TABLE_DIM);

} // namespace p3
} // namespace scream

#endif
