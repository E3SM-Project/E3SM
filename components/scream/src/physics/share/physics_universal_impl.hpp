#ifndef PHYSICS_UNIVERSAL_IMPL_HPP
#define PHYSICS_UNIVERSAL_IMPL_HPP

#include "physics_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_constants.hpp"

namespace scream {
namespace physics {

/*
 * Implementation of universal physics functions. Clients should NOT #include
 * this file, #include physics_functions.hpp instead.
 */

//-----------------------------------------------------------------------------------------------//
// Applies Exners Function which follows:
//   Exner = (P/P0)^(Rd/Cp),
// where,
//   P  is the pressure at this location, Pa
//   P0 is a reference pressure, Pa
//   Rd is the gas constant, J/K
//   Cp is heat capacity of dry air, J/Ki
// All universal constants, P0, Rd, and Cp, are defined in physics_constants.hpp
// Note: Another experssion for Exner is,
//   Exner = T/th
// whre,
//   T  is the temperature, K
//   th is the potential temperature, K
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::get_exner(const Spack& P, const Smask& range_mask)
{
  Spack result;
  
  static constexpr Scalar p0     = C::P0;
  static constexpr Scalar rd     = C::RD;
  static constexpr Scalar inv_cp = C::INV_CP;
  const Spack exner = pow( P/p0, rd*inv_cp );
  // Check that there are no obvious errors in the result.
  EKAT_KERNEL_ASSERT_MSG(!((isnan(exner) && range_mask).any()), "Error in get_exner, Exner has NaN values.\n"); // exit with an error message
  EKAT_KERNEL_ASSERT_MSG(!(((exner <= 0) && range_mask).any()), "Error in get_exner, Exner has negative values.\n"); // exit with an error message
  // Set the values of the result
  result.set(range_mask,exner);
  return result;
}
//-----------------------------------------------------------------------------------------------//
// Converts temperature to potential temperature using Exners function:
//   th_atm = T_atm/exner,
// where
//   th_atm is the potential temperature, K
//   T_atm  is the temperature, K
//   exner  is the exners formula, see definition above, unitless
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::T_to_th(const Spack& T_atm, const Spack& exner, const Smask& range_mask)
{
  Spack result;
  const Spack th_atm = T_atm/exner;
  // Check that there are no obvious errors in the result.
  check_temperature(th_atm,"T_to_th",range_mask);
  // Set the values of the result
  result.set(range_mask,th_atm);
  return result;
}
//-----------------------------------------------------------------------------------------------//
// Converts potential temperature to temperature using Exners function:
//   T_atm = th_atm*exner,
// where
//   th_atm is the potential temperature, K
//   T_atm  is the temperature, K
//   exner  is the exners formula, see definition above, unitless
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::th_to_T(const Spack& th_atm, const Spack& exner, const Smask& range_mask)
{
  Spack result;
  const Spack T_atm = th_atm*exner;
  // Check that there are no obvious errors in the result.
  check_temperature(T_atm,"th_to_T",range_mask);
  // Set the values of the result
  result.set(range_mask,T_atm);
  return result;
}
//-----------------------------------------------------------------------------------------------//
// Determines the vertical layer thickness given the interface heights:
//   dz = zi_top-zi_bot,
// where
//   dz     is the vertical layer thickness, m
//   zi_top is the above surface height of the top of the layer, m
//   zi_bot is the above surface height of the bottom of the layer, m
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::get_dz(const Spack& zi_top, const Spack& zi_bot, const Smask& range_mask)
{
  Spack result;
  const Spack dz = zi_top-zi_bot;
  // Check that there are no obvious errors in the result.
  EKAT_KERNEL_ASSERT_MSG(!((isnan(dz) && range_mask).any()), "Error in get_dz, dz has NaN values.\n"); // exit with an error message
  EKAT_KERNEL_ASSERT_MSG(!(((dz <= 0) && range_mask).any()), "Error in get_dz, dz has negative values.\n"); // exit with an error message
  // Set the values of the result
  result.set(range_mask,dz);
  return result;
}
//-----------------------------------------------------------------------------------------------//
// Compute dry static energy (DSE).
// The result unit is in J/kg
// The inputs are
//   T_mid is the atmospheric temperature. Units in K.
//   z_mid is the geopotential height above surface at midpoints. Units in m.
//   surf_geopotential is the surface geopotential height. Units in m.
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::get_dse(const Spack& T_mid, const Spack& z_mid, const Real surf_geopotential, const Smask& range_mask)
{
  static constexpr Scalar cp  = C::CP;
  static constexpr Scalar ggr = C::gravit;

  const Spack dse(range_mask, cp*T_mid + ggr*z_mid + surf_geopotential);

  // Check that there are no obvious errors in the result.
  EKAT_KERNEL_ASSERT_MSG((isnan(T_mid) && range_mask).none(), "Error in get_dse, T_mid has NaN values.\n"); // exit with an error message
  EKAT_KERNEL_ASSERT_MSG((isnan(z_mid) && range_mask).none(), "Error in get_dse, z_mid has NaN values.\n"); // exit with an error message
  EKAT_KERNEL_ASSERT_MSG(!isnan(surf_geopotential), "Error in get_dse, surf_geopotential has NaN values.\n"); // exit with an error message
  EKAT_KERNEL_ASSERT_MSG(((dse <= 0) && range_mask).none(), "Error in get_dse, dse has negative values.\n"); // exit with an error message

  return dse;
}
//-----------------------------------------------------------------------------------------------//
  // Compute virtual temperature
  // The result unit is in K
  // The inputs are
  //   T_mid is the atmospheric temperature.  Units in K.
  //   qv    is the water vapor mass mixing ratio.  Units in kg/kg
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::get_virtual_temperature(const Spack& T_mid, const Spack& qv, const Smask& range_mask)
{
  static constexpr Scalar ep_2 = C::ep_2;

  const Spack T_virt(range_mask, T_mid*(qv+ep_2)/(ep_2*(1.0+qv)));

// Check that there are no obvious errors in the result.
  EKAT_KERNEL_ASSERT_MSG((isnan(T_mid) && range_mask).none(),  "Error in get_virtual_temperature, T_mid has NaN values.\n"); // exit with an error message
  EKAT_KERNEL_ASSERT_MSG((isnan(qv) && range_mask).none(),     "Error in get_virtual_temperature, qv has NaN values.\n"); // exit with an error message
  EKAT_KERNEL_ASSERT_MSG((isnan(T_virt) && range_mask).none(), "Error in get_virtual_temperature, T_virt has NaN values.\n"); // exit with an error message

  return T_virt;
}

//-----------------------------------------------------------------------------------------------//
} // namespace physics
} // namespace scream

#endif // PHYSICS_UNIVERSAL_IMPL_HPP
