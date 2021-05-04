#ifndef COMMON_PHYSICS_IMPL_HPP
#define COMMON_PHYSICS_IMPL_HPP

#include "share/util/scream_common_physics_functions.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {
  
//-----------------------------------------------------------------------------------------------//
// Computes exner function
// The result is exners formula, and is [dimensionless]
// The input is mid-level pressure, and has units of [Pa]
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::exner_function(const ScalarT& p_mid)
{
  using C = scream::physics::Constants<Real>;
  static const ScalarT p0(C::P0);
  static const ScalarT rd(C::RD);
  static const ScalarT inv_cp(C::INV_CP);

  ScalarT exner = pow( p_mid/p0, rd*inv_cp );
  // Return exner 
  return exner;
}
// For operations on a full column at a time:
template<typename DeviceT>
template<typename ScalarT, typename InputProviderP>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::exner_function(const MemberType& team,
                                               const InputProviderP& p_mid, 
                                               const view_1d<ScalarT>& exner)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,exner.extent(0)),
                       [&] (const int k) {
    exner(k) = exner_function(p_mid(k)); 
  });
}

//-----------------------------------------------------------------------------------------------//
// Converts temperature into potential temperature
// The result is the potential temperature, units in [K]
// The inputs are
//   T_mid is the atmospheric temperature, units in [K]
//   p_mid is the atmospheric pressure, units in [Pa].  p_mid is used in Exners function using `exner` defined above.
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_theta_from_T(const ScalarT& T_mid, const ScalarT& p_mid)
{
  ScalarT exner = exner_function(p_mid);
  ScalarT theta = T_mid/exner;
  // Return theta
  return theta;
}
// For operations on a full column at a time:
template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_theta_from_T(const MemberType& team,
                                                 const InputProviderT& T_mid,
                                                 const InputProviderP& p_mid,
                                                 const view_1d<ScalarT>& theta)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,theta.extent(0)),
                       [&] (const int k) {
    theta(k) = calculate_theta_from_T(T_mid(k),p_mid(k)); 
  });
}

//-----------------------------------------------------------------------------------------------//
// Converts potential temperature into temperature
// The result is the temperature, units in [K]
// The inputs are
//   theta is the potential temperature, units in [K]
//   p_mid is the atmospheric pressure, units in [Pa].  p_mid is used in Exners function using `exner` defined above.
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_T_from_theta(const ScalarT& theta, const ScalarT& p_mid)
{
  ScalarT exner = exner_function(p_mid);
  ScalarT T_mid = theta*exner;
  // Return T
  return T_mid;
}
// For operations on a full column at a time:
template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_T_from_theta(const MemberType& team,
                                                 const InputProviderT& theta,
                                                 const InputProviderP& p_mid,
                                                 const view_1d<ScalarT>& T_mid)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_mid.extent(0)),
                       [&] (const int k) {
    T_mid(k) = calculate_T_from_theta(theta(k),p_mid(k)); 
  });
}
//-----------------------------------------------------------------------------------------------//
// Compute temperature from virtual temperature
// The result unit is in [K]
// The inputs are
//   T_virtual is the virtual temperature.  Units in [K].
//   qv        is the water vapor mass mixing ratio.  Units in [kg/kg]
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_temperature_from_virtual_temperature(const ScalarT& T_virtual, const ScalarT& qv)
{
  using C = scream::physics::Constants<Real>;
  static const ScalarT ep_2(C::ep_2);
  ScalarT T_mid = T_virtual*(ep_2*(1.0+qv))/(qv+ep_2);
  // Return T_mid
  return T_mid;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_temperature_from_virtual_temperature(const MemberType& team,
                                                                         const InputProviderT& T_virtual,
                                                                         const InputProviderQ& qv,
                                                                         const view_1d<ScalarT>& T_mid)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_mid.extent(0)),
                       [&] (const int k) {
    T_mid(k) = calculate_temperature_from_virtual_temperature(T_virtual(k),qv(k)); 
  });
}


//-----------------------------------------------------------------------------------------------//
// Compute virtual temperature
// The result unit is in [K]
// The inputs are
//   T_mid is the atmospheric temperature.  Units in [K].
//   qv    is the water vapor mass mixing ratio.  Units in [kg/kg]
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_virtual_temperature(const ScalarT& T_mid, const ScalarT& qv)
{
  using C = scream::physics::Constants<Real>;
  static const ScalarT ep_2(C::ep_2);
  ScalarT T_virtual = T_mid*(qv+ep_2)/(ep_2*(1.0+qv));
  // Return T_virtual
  return T_virtual;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_virtual_temperature(const MemberType& team,
                                                        const InputProviderT& T_mid,
                                                        const InputProviderQ& qv,
                                                        const view_1d<ScalarT>& T_virtual)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_virtual.extent(0)),
                       [&] (const int k) {
    T_virtual(k) = calculate_virtual_temperature(T_mid(k),qv(k)); 
  });
}

//-----------------------------------------------------------------------------------------------//
// Compute dry static energy (DSE).
// The result unit is in [J/kg]
// The inputs are
//   T_mid is the atmospheric temperature. Units in [K].
//   z_mid is the geopotential height above surface at midpoints. Units in [m].
//   surf_geopotential is the surface geopotential height. Units in [m].
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_dse(const ScalarT& T_mid, const ScalarT& z_mid, const Real surf_geopotential)
{
  using C = scream::physics::Constants<Real>;
  static const ScalarT cp (C::CP);
  static const ScalarT ggr(C::gravit);

  ScalarT dse = cp*T_mid + ggr*z_mid + surf_geopotential;
  return dse;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_dse(const MemberType& team,
                                        const InputProviderT& T_mid,
                                        const InputProviderZ& z_mid,
                                        const Real surf_geopotential,
                                        const view_1d<ScalarT>& dse)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dse.extent(0)),
                       [&] (const int k) {
    dse(k) = calculate_dse(T_mid(k),z_mid(k),surf_geopotential); 
  });
}
//-----------------------------------------------------------------------------------------------//
// Determine the physical thickness of a vertical layer
// The result is dz, units in [m]
// The inputs are
//   pseudo_density is the pressure level thickness, [Pa]
//   p_mid          is the avgerage atmosphere pressure over the level, [Pa]
//   T_mid          is the atmospheric temperature, [K] - needed for T_virtual
//   qv             is the water vapor mass mixing ratio, [kg/kg] - needed for T_virtual
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_dz(const ScalarT& pseudo_density, const ScalarT& p_mid, const ScalarT& T_mid, const ScalarT& qv)
{
  using C = scream::physics::Constants<Real>;
  static const ScalarT Rd  (C::RD);
  static const ScalarT ggr (C::gravit);
  // Need to first back out virtual temperature
  ScalarT T_virtual = calculate_virtual_temperature(T_mid,qv);
  // Now can back out the vertical layer thickness
  ScalarT dz = pseudo_density*Rd*T_virtual / (p_mid*ggr);
  return dz;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderPD, typename InputProviderP, typename InputProviderT, typename InputProviderQ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_dz(const MemberType& team, 
                                       const InputProviderPD& pseudo_density,
                                       const InputProviderP& p_mid,
                                       const InputProviderT& T_mid,
                                       const InputProviderQ& qv,
                                       const view_1d<ScalarT>& dz)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dz.extent(0)),
                       [&] (const int k) {
    dz(k) = calculate_dz(pseudo_density(k),p_mid(k),T_mid(k),qv(k)); 
  });
}
//-----------------------------------------------------------------------------------------------//
// Determine the geopotential height of level interfaces
// The result is z_int, units in [m]
// The input is
//   dz the vertical level thickness, [m]
// Note: Only applicable over an entire column due to the need to integrate over dz.
template<typename DeviceT>
template<typename ScalarT, typename InputProviderZ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_z_int(const MemberType& team,
                                          const int num_levs, 
                                          const InputProviderZ& dz,
                                          const view_1d<ScalarT>& z_int)
{
  using column_ops  = ColumnOps<DeviceT,Real>;
  // Use the column ops scan function.  Note, we set FromTop to false so that we scan from the true bottom
  // of the column to the top.
  constexpr bool FromTop = false;
  Real zbot = 0.0; // TODO: question, is it true zbot is actually always 0?  What about topography?
  column_ops::column_scan<FromTop>(team,num_levs,dz,z_int,zbot);
}
//-----------------------------------------------------------------------------------------------//
} // namespace scream

#endif // COMMON_PHYSICS_IMPL_HPP
