#ifndef COMMON_PHYSICS_FUNCTIONS_HPP
#define COMMON_PHYSICS_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"
#include "share/util/scream_column_ops.hpp"

namespace scream {

template <typename DeviceT>
struct PhysicsFunctions
{

  using Device = DeviceT;
  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;

//-----------------------------------------------------------------------------------------------//
//                                          Functions                                            //
//-----------------------------------------------------------------------------------------------//
// Applies Exners Function which follows:
//   Exner = (P/P0)^(Rd/Cp),
// where,
//   P  is the pressure at this location, [Pa]
//   P0 is a reference pressure, [Pa]
//   Rd is the gas constant, [J/K]
//   Cp is heat capacity of dry air, [J/K]
// All universal constants, P0, Rd, and Cp, are defined in physics_constants.hpp
// Note: Another experssion for Exner is,
//   Exner = T/th
// where,
//   T  is the temperature, [K]
//   th is the potential temperature, [K]
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT exner_function(const ScalarT& pressure);

  template<typename ScalarT, typename InputProviderP>
  KOKKOS_INLINE_FUNCTION
  static void exner_function(const MemberType& team,
                        const InputProviderP& pressure, 
                        const view_1d<ScalarT>& exner);

//-----------------------------------------------------------------------------------------------//
// Converts temperature to potential temperature using Exners function:
//   theta = T_atm/exner(pressure),
// where
//   theta  is the potential temperature, [K]
//   temperature  is the temperature, [K]
//   pressure  is the pressure, [Pa].  Used for exners formula, see definition above, unitless
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_theta_from_T(const ScalarT& temperature, const ScalarT& pressure);

  template<typename ScalarT, typename InputProviderT, typename InputProviderP>
  KOKKOS_INLINE_FUNCTION
  static void calculate_theta_from_T(const MemberType& team,
                               const InputProviderT& temperature,
                               const InputProviderP& pressure,
                               const view_1d<ScalarT>& theta);

//-----------------------------------------------------------------------------------------------//
// Converts potential temperature to temperature using Exners function:
//   temperature = theta*exner(pressure),
// where
//   theta  is the potential temperature, [K]
//   temperature  is the temperature, [K]
//   pressure  is the pressure, [Pa].  Used for exners formula, see definition above, unitless
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_T_from_theta(const ScalarT& theta, const ScalarT& pressure);

  template<typename ScalarT, typename InputProviderT, typename InputProviderP>
  KOKKOS_INLINE_FUNCTION
  static void calculate_T_from_theta(const MemberType& team,
                               const InputProviderT& theta,
                               const InputProviderP& pressure,
                               const view_1d<ScalarT>& temperature);

//-----------------------------------------------------------------------------------------------//
// Compute temperature from virtual temperature
//   temperature = T_virtual * (ep_2*(1+qv)/(qv+ep_2)
// where
//   ep_2        is ratio of molecular mass of water to the molecular mass of dry air 
//   T_virtual   is the virtual temperature.  Units in [K].
//   qv          is the water vapor mass mixing ratio.  Units in [kg/kg]
//   temperature is the atmospheric temperature.  Units in [K].
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_temperature_from_virtual_temperature(const ScalarT& T_virtual, const ScalarT& qv);

  template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_temperature_from_virtual_temperature(const MemberType& team,
                                                       const InputProviderT& T_virtual,
                                                       const InputProviderQ& qv,
                                                       const view_1d<ScalarT>& temperature);

//-----------------------------------------------------------------------------------------------//
// Compute virtual temperature
//   T_virtual = temperature * (qv+ep_2)/(qv+1)
// where
//   ep_2        is ratio of molecular mass of water to the molecular mass of dry air 
//   temperature is the atmospheric temperature.  Units in [K].
//   qv          is the water vapor mass mixing ratio.  Units in [kg/kg]
//   T_virtual   is the virtual temperature.  Units in [K].
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_virtual_temperature(const ScalarT& temperature, const ScalarT& qv);

  template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_virtual_temperature(const MemberType& team,
                                      const InputProviderT& temperature,
                                      const InputProviderQ& qv,
                                      const view_1d<ScalarT>& T_virtual);

//-----------------------------------------------------------------------------------------------//
// Compute dry static energy (DSE).
//   DSE = Cp*temperature + ggr*z + surf_geopotential
// where
//   Cp                is the heat constant of air at constant pressure [J/kg]
//   ggr               is the gravitational constant [m s-2] 
//   temperature       is the atmospheric temperature. Units in [K].
//   z                 is the geopotential height above surface at midpoints. Units in [m].
//   surf_geopotential is the surface geopotential height. Units in [m].
//   DSE               is the dry static energy.  Units in [J/kg]
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_dse(const ScalarT& temperature, const ScalarT& z, const Real surf_geopotential);

  template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_dse(const MemberType& team,
                      const InputProviderT& temperature,
                      const InputProviderZ& z,
                      const Real surf_geopotential,
                      const view_1d<ScalarT>& dse);

//-----------------------------------------------------------------------------------------------//
// Determines the vertical layer thickness using the equation of state:
//   dz = - (-pseudo_density)*Rd*T_virtual / (p_mid*g)
//     note the extra negative sign because the pseudo_density in the model is measured in the positive direction.
// where
//   dz             is the vertical layer thickness, [m]
//   pseudo_density is the pressure level thickness, [Pa]
//   T_virtual      is the virtual temperature - calculated using a separate function from this suite, [K]
//   p_mid          is the avgerage atmosphere pressure over the level, [Pa]
//   g              is the graviational constant, [m s-2]
//   Rd             is the universal gas constant for dry air, [J/kg/K]
//   T_mid          is the atmospheric temperature, [K] - needed for T_virtual
//   qv             is the water vapor mass mixing ratio, [kg/kg] - needed for T_virtual
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_dz(const ScalarT& pseudo_density, const ScalarT& p_mid, const ScalarT& T_mid, const ScalarT& qv);

  template<typename ScalarT, typename InputProviderPD, typename InputProviderP, typename InputProviderT, typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_dz(const MemberType& team, 
                     const InputProviderPD& pseudo_density,
                     const InputProviderP& p_mid,
                     const InputProviderT& T_mid,
                     const InputProviderQ& qv,
                     const view_1d<ScalarT>& dz);

//-----------------------------------------------------------------------------------------------//
// Determines the vertical layer interface height from the vertical layer thicknesses:
//   z_int = int_0^z(dz)
// where
//   dz             is the vertical layer thickness, [m]
// Note: because this function does an integral it cannot be run just on a single level.  It requires
// the full column wise integration.
  template<typename ScalarT, typename InputProviderZ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_z_int(const MemberType& team,
                        const int num_levs, 
                        const InputProviderZ& dz,
                        const view_1d<ScalarT>& z_int);



}; // struct PhysicsFunctions

} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "scream_common_physics_impl.hpp"
#endif

#endif // COMMON_PHYSICS_FUNCTIONS_HPP
