#ifndef COMMON_PHYSICS_IMPL_HPP
#define COMMON_PHYSICS_IMPL_HPP

#include "share/util/scream_common_physics_functions.hpp"
#include "physics/share/physics_constants.hpp"

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
  // Compute temperature from virtual temperature
  // The result unit is in K
  // The inputs are
  //   T_virtual is the virtual temperature.  Units in K.
  //   qv        is the water vapor mass mixing ratio.  Units in kg/kg
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::get_temperature_from_virtual_temperature(const ScalarT& T_virtual, const ScalarT& qv)
{
  using C = scream::physics::Constants<ScalarT>;
  static constexpr ScalarT ep_2(C::ep_2);
  ScalarT T_mid = T_virtual*(ep_2*(1.0+qv))/(qv+ep_2);
  // Return T_mid
  return T_mid;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_temperature_from_virtual_temperature(const MemberType& team,
                                                                         const InputProviderT& T_virtual,
                                                                         const InputProviderQ& qv,
                                                                         const view_1d<ScalarT>& T_mid)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_mid.extent(0)),
                       [&] (const int k) {
    T_mid(k) = get_temperature_from_virtual_temperature(T_virtual(k),qv(k)); 
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
// Determines the vertical layer thickness using the equation of state:
//   dz = - (-psuedo_density)*Rd*T_virtual / (p_mid*g)
//     note the extra negative sign because the psuedo_density in the model is measured in the positive direction.
// where
//   dz             is the vertical layer thickness, m
//   psuedo_density is the pressure level thickness, Pa
//   T_virtual      is the virtual temperature - calculated using a separate function from this suite, K
//   p_mid          is the avgerage atmosphere pressure over the level, Pa
//   g              is the graviational constant, m s-2
//   Rd             is the universal gas constant for dry air, J/kg/K
//   T_mid          is the atmospheric temperature, K - needed for T_virtual
//   qv             is the water vapor mass mixing ratio, kg/kg - needed for T_virtual
template<typename DeviceT>
template<typename ScalarT>
KOKKOS_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::get_dz(const ScalarT& psuedo_density, const ScalarT& p_mid, const ScalarT& T_mid, const ScalarT& qv)
{
  using C = scream::physics::Constants<ScalarT>;
  static constexpr ScalarT Rd  (C::RD);
  static constexpr ScalarT ggr (C::gravit);
  // Need to first back out virtual temperature
  ScalarT T_virtual = get_virtual_temperature(T_mid,qv);
  // Now can back out the vertical layer thickness
  ScalarT dz = psuedo_density*Rd*T_virtual / (p_mid*ggr);
  return dz;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderPD, typename InputProviderP, typename InputProviderT, typename InputProviderQ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_dz(const MemberType& team, 
                                       const InputProviderPD& psuedo_density,
                                       const InputProviderP& p_mid,
                                       const InputProviderT& T_mid,
                                       const InputProviderQ& qv,
                                       const view_1d<ScalarT>& dz)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dz.extent(0)),
                       [&] (const int k) {
    dz(k) = get_dz(psuedo_density(k),p_mid(k),T_mid(k),qv(k)); 
  });
}
//-----------------------------------------------------------------------------------------------//
// Determines the vertical layer interface height from the vertical layer thicknesses:
//   z_int = int_0^z(dz)
// where
//   dz             is the vertical layer thickness, m
// Note: because this function does an integral it cannot be run just on a single level.  It requires
// the full column wise integration.
template<typename DeviceT>
template<typename ScalarT, typename InputProviderZ>
KOKKOS_FUNCTION
void PhysicsFunctions<DeviceT>::get_z_int(const MemberType& team, 
                                          const InputProviderZ& dz,
                                          const view_1d<ScalarT>& z_int)
{
  using column_ops  = ColumnOps<DeviceT,ScalarT>;
  using pack_type = typename ekat::Pack<ScalarT,1>;
  // Use the column ops scan function.  Note, we set FromTop to false so that we scan from the true bottom
  // of the column to the top.
  constexpr bool FromTop = false;
  int num_levs = z_int.extent(0)-1;
  ScalarT zbot = 0.0;
  column_ops::column_scan<FromTop>(team,num_levs,dz,z_int,zbot);
}
//-----------------------------------------------------------------------------------------------//
} // namespace physics
} // namespace scream

#endif // COMMON_PHYSICS_IMPL_HPP
