#ifndef SCREAM_COMMON_PHYSICS_IMPL_HPP
#define SCREAM_COMMON_PHYSICS_IMPL_HPP

#include "physics/share/physics_constants.hpp"
#include "share/util/scream_column_ops.hpp"

namespace scream {

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::exner_function(const ScalarT& pressure)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto p0 = C::P0;
  static constexpr auto rd = C::RD;
  static constexpr auto inv_cp = C::INV_CP;

  return pow( pressure/p0, rd*inv_cp );
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderP>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::exner_function(const MemberType& team,
                                               const InputProviderP& pressure,
                                               const view_1d<ScalarT>& exner)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,exner.extent(0)),
                       [&] (const int k) {
    exner(k) = exner_function(pressure(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_theta_from_T(const ScalarT& temperature, const ScalarT& pressure)
{
  return temperature/exner_function(pressure);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_theta_from_T(const MemberType& team,
                                                       const InputProviderT& temperature,
                                                       const InputProviderP& pressure,
                                                       const view_1d<ScalarT>& theta)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,theta.extent(0)),
                       [&] (const int k) {
    theta(k) = calculate_theta_from_T(temperature(k),pressure(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_T_from_theta(const ScalarT& theta, const ScalarT& pressure)
{
  return theta*exner_function(pressure);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_T_from_theta(const MemberType& team,
                                                       const InputProviderT& theta,
                                                       const InputProviderP& pressure,
                                                       const view_1d<ScalarT>& temperature)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,temperature.extent(0)),
                       [&] (const int k) {
    temperature(k) = calculate_T_from_theta(theta(k),pressure(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_temperature_from_virtual_temperature(const ScalarT& T_virtual, const ScalarT& qv)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto ep_2 = C::ep_2;
  static constexpr auto one  = C::ONE;
  static constexpr auto c1   = - one + one/ep_2;

  return T_virtual / ( one + c1*qv );
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::
calculate_temperature_from_virtual_temperature(const MemberType& team,
                                               const InputProviderT& T_virtual,
                                               const InputProviderQ& qv,
                                               const view_1d<ScalarT>& temperature)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,temperature.extent(0)),
                       [&] (const int k) {
    temperature(k) = calculate_temperature_from_virtual_temperature(T_virtual(k),qv(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_virtual_temperature(const ScalarT& temperature, const ScalarT& qv)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto ep_2 = C::ep_2;
  static constexpr auto one  = C::ONE;
  static constexpr auto c1   = - one + one/ep_2;

  return temperature * ( one + c1*qv );
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::
calculate_virtual_temperature(const MemberType& team,
                              const InputProviderT& temperature,
                              const InputProviderQ& qv,
                              const view_1d<ScalarT>& T_virtual)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_virtual.extent(0)),
                       [&] (const int k) {
    T_virtual(k) = calculate_virtual_temperature(temperature(k),qv(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_dse(const ScalarT& temperature, const ScalarT& z, const Real surf_geopotential)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto cp = C::CP;
  static constexpr auto g  = C::gravit;

  return cp*temperature + g*z + surf_geopotential;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_dse(const MemberType& team,
                                              const InputProviderT& temperature,
                                              const InputProviderZ& z,
                                              const Real surf_geopotential,
                                              const view_1d<ScalarT>& dse)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dse.extent(0)),
                       [&] (const int k) {
    dse(k) = calculate_dse(temperature(k),z(k),surf_geopotential);
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_wetmmr_from_drymmr(const ScalarT& drymmr, const ScalarT& qv)
{
  return drymmr*(1-qv);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderX, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_wetmmr_from_drymmr(const MemberType& team,
                                                             const InputProviderX& drymmr,
                                                             const InputProviderQ& qv,
                                                             const view_1d<ScalarT>& wetmmr)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,wetmmr.extent(0)),
                       [&] (const int k) {
                         wetmmr(k) = calculate_wetmmr_from_drymmr(drymmr(k),qv(k));
                       });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_drymmr_from_wetmmr(const ScalarT& wetmmr, const ScalarT& qv)
{
  return wetmmr/(1-qv);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderX, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_drymmr_from_wetmmr(const MemberType& team,
                                                             const InputProviderX& wetmmr,
                                                             const InputProviderQ& qv,
                                                             const view_1d<ScalarT>& drymmr)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,drymmr.extent(0)),
                       [&] (const int k) {
                         drymmr(k) = calculate_drymmr_from_wetmmr(wetmmr(k),qv(k));
                       });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_dz(const ScalarT& pseudo_density, const ScalarT& p_mid, const ScalarT& T_mid, const ScalarT& qv)
{
  using C = scream::physics::Constants<Real>;

  const ScalarT& T_virtual = calculate_virtual_temperature(T_mid,qv);

  static constexpr auto rd = C::RD;
  static constexpr auto g  = C::gravit;
  return (rd/g)*pseudo_density*T_virtual / p_mid;
}

template<typename DeviceT>
template<typename ScalarT,
         typename InputProviderPD, typename InputProviderP,
         typename InputProviderT,  typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
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

template<typename DeviceT>
template<typename ScalarT, typename InputProviderZ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_z_int(const MemberType& team,
                                                const int num_levs,
                                                const InputProviderZ& dz,
                                                const Real z_surf,
                                                const view_1d<ScalarT>& z_int)
{
  using column_ops  = ColumnOps<DeviceT,Real>;
  // Note, we set FromTop to false since we are prescribing the *bottom* elevation.
  constexpr bool FromTop = false;
  column_ops::template column_scan<FromTop>(team,num_levs,dz,z_int,z_surf);
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_vmr_from_mmr(const Real& gas_mol_weight, const ScalarT& qv, const ScalarT& mmr)
{
  using C = scream::physics::Constants<Real>;
  constexpr Real air_mol_weight   = C::MWdry;

  return mmr / (1.0 - qv) * air_mol_weight/gas_mol_weight;

}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderQ, typename InputProviderX>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_vmr_from_mmr(const MemberType& team,
                                                       const Real gas_mol_weight,
                                                       const InputProviderQ& qv,
                                                       const InputProviderX& mmr,
                                                       const view_1d<ScalarT>& vmr)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,vmr.extent(0)),
                       [&] (const int k) {
    vmr(k) = calculate_vmr_from_mmr(gas_mol_weight,qv(k),mmr(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_mmr_from_vmr(const Real& gas_mol_weight, const ScalarT& qv, const ScalarT& vmr)
{
  using C = scream::physics::Constants<Real>;
  constexpr Real air_mol_weight   = C::MWdry;
  const Real mol_weight_ratio = gas_mol_weight/air_mol_weight;

  return mol_weight_ratio * vmr * (1.0 - qv); 

}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderQ, typename InputProviderX>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_mmr_from_vmr(const MemberType& team,
                                                       const Real gas_mol_weight,
                                                       const InputProviderQ& qv,
                                                       const InputProviderX& vmr,
                                                       const view_1d<ScalarT>& mmr)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,mmr.extent(0)),
                       [&] (const int k) {
    mmr(k) = calculate_mmr_from_vmr(gas_mol_weight,qv(k),vmr(k));
  });
}

} // namespace scream

#endif // SCREAM_COMMON_PHYSICS_IMPL_HPP
