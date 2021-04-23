#ifndef COMMON_PHYSICS_FUNCTIONS_HPP
#define COMMON_PHYSICS_FUNCTIONS_HPP

#include "physics_constants.hpp"

namespace scream {
namespace physics {

template <typename DeviceT>
struct PhysicsFunctions
{

  using Device = DeviceT;
  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;

  //--------------------------------------------------------------------//
  //  Functions
  //--------------------------------------------------------------------//

  // Computes exner function
  // The result is exners formula, and is dimensionless
  // The input is mid-level pressure, and has units of Pa
  template<typename ScalarT>
  KOKKOS_FUNCTION
  static ScalarT get_exner(const ScalarT& p_mid);

  template<typename ScalarT, typename InputProviderP>
  KOKKOS_FUNCTION
  static void get_exner(const MemberType& team,
                        const InputProviderP& p_mid, 
                        const view_1d<ScalarT>& exner);

  // Converts temperature into potential temperature
  // The result is the potential temperature, units in K
  // The inputs are
  //   T_mid is the atmospheric temperature, units in K
  //   p_mid is the atmospheric pressure, units in Pa.  p_mid is used in Exners function using `get_exner` defined above.
  template<typename ScalarT>
  KOKKOS_FUNCTION
  static ScalarT get_theta_from_T(const ScalarT& T_mid, const ScalarT& p_mid);

  template<typename ScalarT, typename InputProviderT, typename InputProviderP>
  KOKKOS_FUNCTION
  static void get_theta_from_T(const MemberType& team,
                               const InputProviderT& T_mid,
                               const InputProviderP& p_mid,
                               const view_1d<ScalarT>& theta);

  // Converts potential temperature into temperature
  // The result is the temperature, units in K
  // The inputs are
  //   theta is the potential temperature, units in K
  //   p_mid is the atmospheric pressure, units in Pa.  p_mid is used in Exners function using `get_exner` defined above.
  template<typename ScalarT>
  KOKKOS_FUNCTION
  static ScalarT get_T_from_theta(const ScalarT& theta, const ScalarT& p_mid);

  template<typename ScalarT, typename InputProviderT, typename InputProviderP>
  KOKKOS_FUNCTION
  static void get_T_from_theta(const MemberType& team,
                               const InputProviderT& theta,
                               const InputProviderP& p_mid,
                               const view_1d<ScalarT>& T_mid);

  // Compute virtual temperature
  // The result unit is in K
  // The inputs are
  //   T_mid is the atmospheric temperature.  Units in K.
  //   qv    is the water vapor mass mixing ratio.  Units in kg/kg
  template<typename ScalarT>
  KOKKOS_FUNCTION
  static ScalarT get_virtual_temperature(const ScalarT& T_mid, const ScalarT& qv);

  template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
  KOKKOS_FUNCTION
  static void get_virtual_temperature(const MemberType& team,
                                      const InputProviderT& T_mid,
                                      const InputProviderQ& qv,
                                      const view_1d<ScalarT>& T_virtual);

  // Compute dry static energy (DSE).
  // The result unit is in J/kg
  // The inputs are
  //   T_mid is the atmospheric temperature. Units in K.
  //   z_mid is the geopotential height above surface at midpoints. Units in m.
  //   surf_geopotential is the surface geopotential height. Units in m.
  template<typename ScalarT>
  KOKKOS_FUNCTION
  static ScalarT get_dse(const ScalarT& T_mid, const ScalarT& z_mid, const Real surf_geopotential);

  template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
  KOKKOS_FUNCTION
  static ScalarT get_dse(const MemberType& team,
                         const InputProviderT& T_mid,
                         const InputProviderZ& z_mid,
                         const Real surf_geopotential,
                         const view_1d<ScalarT>& dse);


}; // struct PhysicsFunctions

} // namespace physics
} // namespace scream

#endif // COMMON_PHYSICS_FUNCTIONS_HPP
