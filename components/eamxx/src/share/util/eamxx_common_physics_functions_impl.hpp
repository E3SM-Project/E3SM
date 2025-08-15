#ifndef SCREAM_COMMON_PHYSICS_IMPL_HPP
#define SCREAM_COMMON_PHYSICS_IMPL_HPP

#include "physics/share/physics_constants.hpp"
#include "share/util/eamxx_column_ops.hpp"

namespace scream {

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_dx_from_area(const ScalarT& area, const ScalarT& lat)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto coeff_1 = C::earth_ellipsoid1;
  static constexpr auto coeff_2 = C::earth_ellipsoid2;
  static constexpr auto coeff_3 = C::earth_ellipsoid3;
  static constexpr auto pi      = C::Pi;

  // Compute latitude in radians
  auto lat_in_rad = lat*(pi/180.0);

  // Now find meters per degree latitude
  // Below equation finds distance between two points on an ellipsoid, derived from expansion
  // taking into account ellipsoid using World Geodetic System (WGS84) reference
  auto m_per_degree_lat = coeff_1 - coeff_2 * std::cos(2.0*lat_in_rad) + coeff_3 * std::cos(4.0*lat_in_rad);
  // Note, for the formula we need to convert area from radians to degrees.
  return m_per_degree_lat * std::sqrt(area)*(180.0/pi);

}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_density(const ScalarT& pseudo_density, const ScalarT& dz)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto g = C::gravit;

  return pseudo_density/dz/g;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderP, typename InputProviderZ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_density(const MemberType& team,
                                                  const InputProviderP& pseudo_density,
                                                  const InputProviderZ& dz,
                                                  const view_1d<ScalarT>& density)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,density.extent(0)),
                       [&] (const int k) {
    density(k) = calculate_density(pseudo_density(k),dz(k));
  });
}

template <typename DeviceT>
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_vertical_velocity(const ScalarT& omega, const ScalarT& density)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto g = C::gravit;

  return -omega/(density * g);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderOmega, typename InputProviderRho>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_vertical_velocity(const MemberType& team,
                                                  const InputProviderOmega& omega,
                                                  const InputProviderRho& rho,
                                                  const view_1d<ScalarT>& w)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, w.extent(0)),
    [&] (const int k) {
      w(k) = calculate_vertical_velocity(omega(k), rho(k));
    });
}

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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,exner.extent(0)),
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
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_thetal_from_theta(const ScalarT& theta, const ScalarT& temperature, const ScalarT& qc)
{
  using C = scream::physics::Constants<Real>;

  return theta - (theta / temperature) * (C::LatVap/C::Cpair) * qc;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_theta_from_T(const MemberType& team,
                                                       const InputProviderT& temperature,
                                                       const InputProviderP& pressure,
                                                       const view_1d<ScalarT>& theta)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,theta.extent(0)),
                       [&] (const int k) {
    theta(k) = calculate_theta_from_T(temperature(k),pressure(k));
  });
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderTheta, typename InputProviderT, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_thetal_from_theta(const MemberType& team,
                                                            const InputProviderTheta& theta,
                                                            const InputProviderT& temperature,
                                                            const InputProviderQ& qc,
                                                            const view_1d<ScalarT>& thetal)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,thetal.extent(0)),
                       [&] (const int k) {
    thetal(k) = calculate_thetal_from_theta(theta(k),temperature(k),qc(k));
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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,temperature.extent(0)),
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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,temperature.extent(0)),
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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,T_virtual.extent(0)),
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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dse.extent(0)),
                       [&] (const int k) {
    dse(k) = calculate_dse(temperature(k),z(k),surf_geopotential);
  });
}

template<typename DeviceT>
template<typename ScalarT>

KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_temperature_from_dse(const ScalarT& dse, const ScalarT& z, const Real surf_geopotential)
{
  using C = scream::physics::Constants<Real>;

  static constexpr auto cp = C::CP;
  static constexpr auto g  = C::gravit;

  return (dse - g*z - surf_geopotential)/cp;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::
calculate_temperature_from_dse(const MemberType& team,
                               const InputProviderT& dse,
                               const InputProviderZ& z,
                               const Real surf_geopotential,
                               const view_1d<ScalarT>& temperature)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dse.extent(0)),
                       [&] (const int k) {
    temperature(k) = calculate_temperature_from_dse(dse(k),z(k),surf_geopotential);
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_wetmmr_from_drymmr(const ScalarT& drymmr, const ScalarT& qv_dry)
{
  return drymmr/(1+qv_dry);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderX, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_wetmmr_from_drymmr(const MemberType& team,
                                                             const InputProviderX& drymmr,
                                                             const InputProviderQ& qv_dry,
                                                             const view_1d<ScalarT>& wetmmr)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,wetmmr.extent(0)),
                       [&] (const int k) {
                         wetmmr(k) = calculate_wetmmr_from_drymmr(drymmr(k),qv_dry(k));
                       });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_wetmmr_from_drymmr_dp_based(const ScalarT& drymmr,
                                      const ScalarT& pseudo_density, const ScalarT& pseudo_density_dry)
{
  return drymmr*pseudo_density_dry/pseudo_density;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderX, typename InputProviderPD>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_wetmmr_from_drymmr_dp_based(const MemberType& team,
                                                             const InputProviderX& drymmr,
                                                             const InputProviderPD& pseudo_density,
                                                             const InputProviderPD& pseudo_density_dry,
                                                             const view_1d<ScalarT>& wetmmr)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,wetmmr.extent(0)),
                       [&] (const int k) {
                         wetmmr(k) = calculate_wetmmr_from_drymmr_dp_based(drymmr(k),pseudo_density(k),pseudo_density_dry(k));
                       });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_drymmr_from_wetmmr(const ScalarT& wetmmr, const ScalarT& qv_wet)
{
  return wetmmr/(1-qv_wet);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderX, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_drymmr_from_wetmmr(const MemberType& team,
                                                             const InputProviderX& wetmmr,
                                                             const InputProviderQ& qv_wet,
                                                             const view_1d<ScalarT>& drymmr)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,drymmr.extent(0)),
                       [&] (const int k) {
                         drymmr(k) = calculate_drymmr_from_wetmmr(wetmmr(k),qv_wet(k));
                       });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_drymmr_from_wetmmr_dp_based(const ScalarT& wetmmr,
                             const ScalarT& pseudo_density, const ScalarT& pseudo_density_dry)
{
  return wetmmr*pseudo_density/pseudo_density_dry;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderX, typename InputProviderPD>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_drymmr_from_wetmmr_dp_based(const MemberType& team,
                                                             const InputProviderX& wetmmr,
                                                             const InputProviderPD& pseudo_density,
                                                             const InputProviderPD& pseudo_density_dry,
                                                             const view_1d<ScalarT>& drymmr)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,drymmr.extent(0)),
                       [&] (const int k) {
                         drymmr(k) = calculate_drymmr_from_wetmmr_dp_based(wetmmr(k),pseudo_density(k),pseudo_density_dry(k));
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
         typename InputProviderT,  typename InputProviderQ,
         typename MT>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_dz(const MemberType& team,
                                             const InputProviderPD& pseudo_density,
                                             const InputProviderP& p_mid,
                                             const InputProviderT& T_mid,
                                             const InputProviderQ& qv,
                                             const view_1d<ScalarT, MT>& dz)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dz.extent(0)),
                       [&] (const int k) {
    dz(k) = calculate_dz(pseudo_density(k),p_mid(k),T_mid(k),qv(k));
  });
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderZ, typename MT>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_z_int(const MemberType& team,
                                                const int num_levs,
                                                const InputProviderZ& dz,
                                                const Real z_surf,
                                                const view_1d<ScalarT, MT>& z_int)
{
  using column_ops  = ColumnOps<DeviceT,Real>;
  // Note, we set FromTop to false since we are prescribing the *bottom* elevation.
  constexpr bool FromTop = false;
  column_ops::template column_scan<FromTop>(team,num_levs,dz,z_int,z_surf);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProvider, typename MT>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_z_mid(const MemberType& team,
                                                const int num_levs,
                                                const InputProvider& z_int,
                                                const view_1d<ScalarT, MT>& z_mid)
{
  using column_ops  = ColumnOps<DeviceT,Real>;
  column_ops::compute_midpoint_values(team,num_levs,z_int,z_mid);
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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,vmr.extent(0)),
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
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team,mmr.extent(0)),
                       [&] (const int k) {
    mmr(k) = calculate_mmr_from_vmr(gas_mol_weight,qv(k),vmr(k));
  });
}

template<typename DeviceT>
KOKKOS_INLINE_FUNCTION
Real PhysicsFunctions<DeviceT>::calculate_surface_air_T(const Real& T_mid_bot, const Real& z_mid_bot)
{
  /*Compute temperature at the bottom of the gridcell closest to the ground. The implementation here
    is really clunky and is meant to provide calculate_psl with the values used by CESM... Think
    twice before using it for anything else. Inputs are T at midpoint of layer closest to the surface (K) and
    the geometric height at that point (m). Note that z_mid_bot is distance from the surface rather than from
    sea level!
  */

  //Old version extrapolated off lowest 2 midpoint values (needs different function arguments).
  //Ditching this version for fear that weird lowest layer T relationships could yield strange surface values.
  //const Real T_weighting = ( p_int_i(num_levs) - p_mid_i(last_entry))/(p_mid_i(last_entry - 1) - p_mid_i(last_entry) );
  //return T_mid_i(last_entry - 1)*T_weighting + T_mid_i(last_entry)*(1-T_weighting);

  //Assume 6.5 K/km lapse rate between cell's midpoint and its bottom edge
  return T_mid_bot + sp(0.0065)*z_mid_bot;
}

template<typename DeviceT>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::lapse_T_for_psl(const Real& T_ground, const Real& phi_ground,
					        Real& lapse, Real& T_ground_tmp )
{
  /*
    Choose lapse rate and effective ground temperature to use for sea-level pressure calculation.
    This function should only be used by calculate_psl and is separated from that function solely
    to improve specificity of unit testing.
 */

  using C = scream::physics::Constants<Real>;
  constexpr Real gravit = C::gravit;

  //Get preliminary surface and sea level temperature to decide on lapse rate
  auto T_sl = T_ground + sp(0.0065)*phi_ground/gravit; //start by assuming lapse rate is 6.5 K/km
  T_ground_tmp= T_ground; //make copy b/c may need to modify later

  if (T_ground<=290.5 && T_sl>290.5){
    lapse=gravit/phi_ground*(sp(290.5)-T_ground); //choose lapse rate to make T at sea level 290.5K
  } else if (T_ground>290.5 && T_sl>290.5){
    lapse=0; //choose lapse rate = 0 to not make things any worse.
    T_ground_tmp=sp(0.5)*(sp(290.5)+T_ground); //reduce effective T in a smooth way to avoid unrealistic result
  } else if (T_ground<255){
  //Following EAM's treatment for cold cases even though it seems overly crude: what if T_sl>255?
  //Or phi_ground<0 so positive lapse makes T_sl even colder? We should eventually ditch this entire scheme.
  lapse=0.0065;
  T_ground_tmp=sp(0.5)*(255+T_ground);
  } else{
    //note lack of "elif T_ground>290.5 and T_sl<290.5" case (phi_ground<0 and hot) is missing on purpose
    //because 6.5K/km lapse rate will cool T_sl in that case and that's what we want.
    lapse = 0.0065; //assume 6.5K/km lapse rate for reasonable temperatures
  }
}

template<typename DeviceT>
KOKKOS_INLINE_FUNCTION
Real PhysicsFunctions<DeviceT>::calculate_psl(const Real& T_ground, const Real& p_ground, const Real& phi_ground)
{
  /*
     Compute sea level pressure (psl) assuming atmosphere below the land surface is dry and has a lapse
     rate of 6.5K/km unless conditions are very warm. See components/eamxx/docs/tech_doc/physics/psl/
     for a description. Note that all input/out variables are only defined at the surface rather than
     being 3d variables so no need to template on InputProvider.
 */


  using C = scream::physics::Constants<Real>;
  constexpr Real gravit = C::gravit;
  constexpr Real Rair = C::Rair;
  Real psl;

  // if phi_ground is very close to sea level already, set psl to existing p_ground
  if (std::abs(phi_ground/gravit) < 1e-4){
    psl = p_ground;
  //otherwise compute
  }else{
    Real lapse;
    Real T_ground_tmp;
    lapse_T_for_psl(T_ground,phi_ground, lapse, T_ground_tmp);
    Real alpha=lapse*Rair/gravit;
    Real beta=phi_ground/(Rair*T_ground_tmp);
    psl = p_ground*std::exp(beta*( 1 - alpha*beta/2 + std::pow(alpha*beta,2)/3 ) );
  }

  return psl;
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::apply_rayleigh_friction(const Real dt, const ScalarT& otau,
                                                        ScalarT& u_wind, ScalarT& v_wind, ScalarT& T_mid)
{
  using C = scream::physics::Constants<Real>;
  constexpr Real cp = C::CP;

  const Real dt_inv = 1.0/dt;

  const ScalarT c2 = 1.0/(1.0 + otau*dt);
  const ScalarT c1 = -1.0*otau*c2;
  const ScalarT c3 = 0.5*(1.0 - c2*c2)*dt_inv;

  const ScalarT u2 = u_wind*u_wind;
  const ScalarT v2 = v_wind*v_wind;

  u_wind += c1*u_wind;
  v_wind += c1*v_wind;
  T_mid  += c3*(u2 + v2)/cp;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderOtau, typename MT>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::apply_rayleigh_friction (const MemberType& team,
                                                         const Real dt,
                                                         const InputProviderOtau& otau,
                                                         const view_1d<ScalarT, MT>& u_wind,
                                                         const view_1d<ScalarT, MT>& v_wind,
                                                         const view_1d<ScalarT, MT>& T_mid)
{
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, T_mid.extent(0)),
                       [&] (const int k) {
    apply_rayleigh_friction(dt, otau(k), u_wind(k), v_wind(k), T_mid(k));
  });
}

} // namespace scream

#endif // SCREAM_COMMON_PHYSICS_IMPL_HPP
