#ifndef COMMON_PHYSICS_IMPL_HPP
#define COMMON_PHYSICS_IMPL_HPP

#include "common_physics_functions.hpp"
#include "physics_constants.hpp"

namespace scream {
namespace physics {

  
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
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::get_exner(const ScalarT& p_mid)
{
  using C = scream::physics::Constants<ScalarT>;
  static constexpr ScalarT p0(C::P0);
  static constexpr ScalarT rd(C::RD);
  static constexpr ScalarT inv_cp(C::INV_CP);

  ScalarT exner = pow( p_mid/p0, rd*inv_cp );
  // Simple check to make sure no obvious errors occurred
  EKAT_REQUIRE_MSG(!isnan(exner),"Error in get_exner, exner(p = " + std::to_string(p_mid) + ") = NaN"); 
  EKAT_REQUIRE_MSG(!(exner<0),   "Error in get_exner, exner(p = " + std::to_string(p_mid) + ") = " + std::to_string(exner) + " < 0");
  // Return exner 
  return exner;
}
// For operations on a full column at a time:
template<typename DeviceT>
template<typename ScalarT, typename InputProviderP>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_exner(const MemberType& team,
                      const InputProviderP& p_mid, 
                      const view_1d<ScalarT>& exner)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,exner.extent(0)),
                       [&] (const int k) {
    exner(k) = get_exner(p_mid(k)); 
  });
}

//-----------------------------------------------------------------------------------------------//
// Converts temperature to potential temperature using Exners function:
//   theta = T_atm/exner(p_mid),
// where
//   theta  is the potential temperature, K
//   T_mid  is the temperature, K
//   p_mid  is the pressure, Pa.  Used for exners formula, see definition above, unitless
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::get_theta_from_T(const ScalarT& T_mid, const ScalarT& p_mid)
{
  ScalarT exner = get_exner(p_mid);
  ScalarT theta = T_mid/exner;
  // Simple check to make sure no obvious errors occurred
  EKAT_REQUIRE_MSG(!isnan(theta),"Error in get_theta_from_T, theta(T = " + std::to_string(T_mid) + ", p = " + std::to_string(p_mid) + ") = NaN"); 
  EKAT_REQUIRE_MSG(!(theta<0),"Error in get_theta_from_T, theta(T = " + std::to_string(T_mid) + ", p = " + std::to_string(p_mid) + ") = " + std::to_string(theta) + " < 0"); 
  // Return theta
  return theta;
}
// For operations on a full column at a time:
template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_theta_from_T(const MemberType& team,
                                                 const InputProviderT& T_mid,
                                                 const InputProviderP& p_mid,
                                                 const view_1d<ScalarT>& theta)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,theta.extent(0)),
                       [&] (const int k) {
    theta(k) = get_theta_from_T(T_mid(k),p_mid(k)); 
  });
}

//-----------------------------------------------------------------------------------------------//
// Converts potential temperature to temperature using Exners function:
//   T_mid = theta*exner(p_mid),
// where
//   theta  is the potential temperature, K
//   T_mid  is the temperature, K
//   p_mid  is the pressure, Pa.  Used for exners formula, see definition above, unitless
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::get_T_from_theta(const ScalarT& theta, const ScalarT& p_mid)
{
  ScalarT exner = get_exner(p_mid);
  ScalarT T_mid = theta*exner;
  // Simple check to make sure no obvious errors occurred
  EKAT_REQUIRE_MSG(!isnan(T_mid),"Error in get_T_from_theta, T(theta = " + std::to_string(theta) + ", p = " + std::to_string(p_mid) + ") = NaN"); 
  EKAT_REQUIRE_MSG(!(T_mid<0),"Error in get_T_from_theta, T(theta = " + std::to_string(theta) + ", p = " + std::to_string(p_mid) + ") = " + std::to_string(T_mid) + " < 0"); 
  // Return T
  return T_mid;
}
// For operations on a full column at a time:
template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_T_from_theta(const MemberType& team,
                                                 const InputProviderT& theta,
                                                 const InputProviderP& p_mid,
                                                 const view_1d<ScalarT>& T_mid)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_mid.extent(0)),
                       [&] (const int k) {
    T_mid(k) = get_T_from_theta(theta(k),p_mid(k)); 
  });
}

//-----------------------------------------------------------------------------------------------//
  // Compute virtual temperature
  // The result unit is in K
  // The inputs are
  //   T_mid is the atmospheric temperature.  Units in K.
  //   qv    is the water vapor mass mixing ratio.  Units in kg/kg
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::get_virtual_temperature(const ScalarT& T_mid, const ScalarT& qv)
{
  using C = scream::physics::Constants<ScalarT>;
  static constexpr ScalarT ep_2(C::ep_2);
  ScalarT T_virtual = T_mid*(qv+ep_2)/(ep_2*(1.0+qv));
  // Simple check to make sure no obvious errors occurred
  EKAT_REQUIRE_MSG(!isnan(T_virtual),"Error in get_virtual_temperature, T_virtual(T = " + std::to_string(T_mid) + ", qv = " + std::to_string(qv) + ") = NaN"); 
  EKAT_REQUIRE_MSG(!(T_virtual<0),"Error in get_virtual_temperature, T_virtual(T = " + std::to_string(T_mid) + ", qc = " + std::to_string(qv) + ") = " + std::to_string(T_virtual) + " < 0");
  // Return T_virtual
  return T_virtual;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_virtual_temperature(const MemberType& team,
                                                        const InputProviderT& T_mid,
                                                        const InputProviderQ& qv,
                                                        const view_1d<ScalarT>& T_virtual)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_virtual.extent(0)),
                       [&] (const int k) {
    T_virtual(k) = get_virtual_temperature(T_mid(k),qv(k)); 
  });
}

//-----------------------------------------------------------------------------------------------//
// Compute dry static energy (DSE).
// The result unit is in J/kg
// The inputs are
//   T_mid is the atmospheric temperature. Units in K.
//   z_mid is the geopotential height above surface at midpoints. Units in m.
//   surf_geopotential is the surface geopotential height. Units in m.
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::get_dse(const ScalarT& T_mid, const ScalarT& z_mid, const Real surf_geopotential)
{
  using C = scream::physics::Constants<ScalarT>;
  static constexpr ScalarT cp (C::CP);
  static constexpr ScalarT ggr(C::gravit);

  ScalarT dse = cp*T_mid + ggr*z_mid + surf_geopotential;
  return dse;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_dse(const MemberType& team,
                                        const InputProviderT& T_mid,
                                        const InputProviderZ& z_mid,
                                        const Real surf_geopotential,
                                        const view_1d<ScalarT>& dse)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dse.extent(0)),
                       [&] (const int k) {
    dse(k) = get_dse(T_mid(k),z_mid(k),surf_geopotential); 
  });
}
//-----------------------------------------------------------------------------------------------//
} // namespace physics
} // namespace scream

#endif // COMMON_PHYSICS_IMPL_HPP
