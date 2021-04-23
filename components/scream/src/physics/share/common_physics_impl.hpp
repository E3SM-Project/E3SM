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

  return pow( p_mid/p0, rd*inv_cp );
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

} // namespace physics
} // namespace scream

#endif // COMMON_PHYSICS_IMPL_HPP
