#ifndef SCREAM_COMMON_PHYSICS_FUNCTIONS_HPP
#define SCREAM_COMMON_PHYSICS_FUNCTIONS_HPP

#include "share/scream_types.hpp"

namespace scream {

template <typename DeviceT>
struct PhysicsFunctions
{
  // ---------------------------------------------------------------- //
  //                     Single-scalar Functions                      //
  // ---------------------------------------------------------------- //
  //                                                                  //
  // These functions act on a single scalar. They implement laws and  //
  // formulas that can be used across the whole atm component.        //
  // The functions are templated on the scalar type, so that they     //
  // can be used with builtin scalars inputs (double, float) as well  //
  // as with derived types (e.g., ekat::Pack<T,N>), provided that the //
  // required operators/functions are overloaded for such types.      //
  //                                                                  //
  // NOTE: these functions are technically independent from DeviceT,  //
  //       so they should probably be moved out of the class, and     //
  //       into a namespace.                                          //
  // ---------------------------------------------------------------- //

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
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT exner_function(const ScalarT& pressure);

  //-----------------------------------------------------------------------------------------------//
  // Converts temperature to potential temperature using Exners function:
  //   theta = T_atm/exner(pressure),
  // where
  //   theta  is the potential temperature, [K]
  //   temperature  is the temperature, [K]
  //   pressure  is the pressure, [Pa].  Used for exners formula, see definition above, unitless
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_theta_from_T(const ScalarT& temperature, const ScalarT& pressure);

  //-----------------------------------------------------------------------------------------------//
  // Converts potential temperature to temperature using Exners function:
  //   temperature = theta*exner(pressure),
  // where
  //   theta  is the potential temperature, [K]
  //   temperature  is the temperature, [K]
  //   pressure  is the pressure, [Pa].  Used for exners formula, see definition above, unitless
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_T_from_theta(const ScalarT& theta, const ScalarT& pressure);

  //-----------------------------------------------------------------------------------------------//
  // Compute temperature from virtual temperature
  //   temperature = T_virtual * (ep_2*(1+qv)/(qv+ep_2)
  // where
  //   ep_2        is ratio of molecular mass of water to the molecular mass of dry air
  //   T_virtual   is the virtual temperature.  Units in [K].
  //   qv          is the water vapor mass mixing ratio.  Units in [kg/kg]
  //   temperature is the atmospheric temperature.  Units in [K].
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_temperature_from_virtual_temperature(const ScalarT& T_virtual, const ScalarT& qv);

  //-----------------------------------------------------------------------------------------------//
  // Compute virtual temperature
  //   T_virtual = temperature * (qv+ep_2)/(qv+1)
  // where
  //   ep_2        is ratio of molecular mass of water to the molecular mass of dry air
  //   temperature is the atmospheric temperature.  Units in [K].
  //   qv          is the water vapor mass mixing ratio.  Units in [kg/kg]
  //   T_virtual   is the virtual temperature.  Units in [K].
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_virtual_temperature(const ScalarT& temperature, const ScalarT& qv);

  //-----------------------------------------------------------------------------------------------//
  // Compute dry static energy (DSE).
  //   DSE = Cp*temperature + g*z + surf_geopotential
  // where
  //   Cp                is the heat constant of air at constant pressure [J/kg]
  //   g                 is the gravitational constant [m s-2]
  //   temperature       is the atmospheric temperature. Units in [K].
  //   z                 is the geopotential height above surface at midpoints. Units in [m].
  //   surf_geopotential is the surface geopotential height. Units in [m].
  //   DSE               is the dry static energy.  Units in [J/kg]
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_dse(const ScalarT& temperature, const ScalarT& z, const Real surf_geopotential);

  //-----------------------------------------------------------------------------------------------//
  // Determines the vertical layer thickness using the equation of state:
  //   dz = - (-pseudo_density)*Rd*T_virtual / (p_mid*g)
  // where
  //   dz             is the vertical layer thickness, [m]
  //   pseudo_density is the pressure level thickness, [Pa]
  //   T_virtual      is the virtual temperature - calculated using a separate function from this suite, [K]
  //   p_mid          is the avgerage atmosphere pressure over the level, [Pa]
  //   g              is the graviational constant, [m s-2]
  //   Rd             is the universal gas constant for dry air, [J/kg/K]
  //   T_mid          is the atmospheric temperature, [K] - needed for T_virtual
  //   qv             is the water vapor mass mixing ratio, [kg/kg] - needed for T_virtual
  //
  // Note: the extra negative sign is due to the fact that the pseudo_density
  // in the model is measured in the positive direction.
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_dz (const ScalarT& pseudo_density, const ScalarT& p_mid,
                               const ScalarT& T_mid, const ScalarT& qv);

  //-----------------------------------------------------------------------------------------------//
  // Calculate the volume mixing ratio given the wet mass mixing ratio:
  //   X_vmr = X_mmr / (1 - qv) * mol_weight_air/mol_weight_X
  // where
  //   X_vmr          is the volume mixing ratio X
  //   X_mmr          is the mass mixing ratio of X
  //   qv             is the water vapor mass mixing ratio
  //   mol_weight_air is the molecular weight of dry air
  //   mol_weight_X   is the molecular weight of X
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_vmr_from_mmr(const Real& gas_mol_weight, const ScalarT& qv, const ScalarT& mmr);

  //-----------------------------------------------------------------------------------------------//
  // Calculate wet mass mixing ratio the given the volume mixing ratio:
  //   X_mmr = a*X_vmr * (1 -qv)
  // where
  //   X_vmr          is the volume mixing ratio X
  //   X_mmr          is the mass mixing ratio of X
  //   qv             is the water vapor mass mixing ratio
  //   mol_weight_air is the molecular weight of dry air
  //   mol_weight_X   is the molecular weight of X
  //   a = mol_weight_X/mol_weight_air is the ratio of the molecular weight of the gas to the molecular weight of dry air
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_mmr_from_vmr(const Real& gas_mol_weight, const ScalarT& qv, const ScalarT& vmr);

  // ---------------------------------------------------------------- //
  //                     Whole column Functions                       //
  // ---------------------------------------------------------------- //
  //                                                                  //
  // These functions act on a 1d view, representing a whole column    //
  // worth of scalars. As for the single-scalar versions, they can be //
  // used for builtin scalars as well as Pack-type scalars.           //
  // These functions are meant to be called from within a parallel    //
  // region; more precisely, they should be called from within the    //
  // outermost parallel_for, which should have been dispatched with a //
  // TeamPolicy. In other words, these routines will dispatch a       //
  // parallel_for using a TeamThreadRange policy. No third layer of   //
  // parallelism will be used (no ThreadVectorRange for loop), so the //
  // team policy on GPU should be built with vector_length=1.         //
  // The routines are templated on the type of their inputs, so that  //
  // these can be 1d views as well as other routines (e.g. lambdas).  //
  // The only requirement is that the provider supports op()(int k),  //
  // which should return scalar corresponding to level k (which can   //
  // be a physical level or a "packed level", depending on the scalar //
  // type.                                                            //
  //                                                                  //
  // Most of these routines simply call the homonymous routine that   //
  // act on a single scalar, except for calculate_z_int, which only   //
  // makes sense for a whole column (i.e., there is not single-scalar //
  // version of compute_z_int).                                       //
  // ---------------------------------------------------------------- //

  using Device = DeviceT;
  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;

  template<typename ScalarT, typename InputProviderP>
  KOKKOS_INLINE_FUNCTION
  static void exner_function (const MemberType& team,
                              const InputProviderP& pressure,
                              const view_1d<ScalarT>& exner);


  template<typename ScalarT, typename InputProviderT, typename InputProviderP>
  KOKKOS_INLINE_FUNCTION
  static void calculate_theta_from_T (const MemberType& team,
                                      const InputProviderT& temperature,
                                      const InputProviderP& pressure,
                                      const view_1d<ScalarT>& theta);

  template<typename ScalarT, typename InputProviderT, typename InputProviderP>
  KOKKOS_INLINE_FUNCTION
  static void calculate_T_from_theta (const MemberType& team,
                                      const InputProviderT& theta,
                                      const InputProviderP& pressure,
                                      const view_1d<ScalarT>& temperature);

  template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_temperature_from_virtual_temperature (const MemberType& team,
                                                              const InputProviderT& T_virtual,
                                                              const InputProviderQ& qv,
                                                              const view_1d<ScalarT>& temperature);

  template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_virtual_temperature(const MemberType& team,
                                            const InputProviderT& temperature,
                                            const InputProviderQ& qv,
                                            const view_1d<ScalarT>& T_virtual);


  template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_dse (const MemberType& team,
                             const InputProviderT& temperature,
                             const InputProviderZ& z,
                             const Real surf_geopotential,
                             const view_1d<ScalarT>& dse);

  template<typename ScalarT,
           typename InputProviderPD, typename InputProviderP,
           typename InputProviderT,  typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_dz (const MemberType& team,
                            const InputProviderPD& pseudo_density,
                            const InputProviderP& p_mid,
                            const InputProviderT& T_mid,
                            const InputProviderQ& qv,
                            const view_1d<ScalarT>& dz);

  template<typename ScalarT, typename InputProviderQ, typename InputProviderX>
  KOKKOS_INLINE_FUNCTION
  static void calculate_vmr_from_mmr(const MemberType& team,
                                     const Real gas_mol_weight,
                                     const InputProviderQ& qv,
                                     const InputProviderX& mmr,
                                     const view_1d<ScalarT>& vmr);

  template<typename ScalarT, typename InputProviderQ, typename InputProviderX>
  KOKKOS_INLINE_FUNCTION
  static void calculate_mmr_from_vmr(const MemberType& team,
                                     const Real gas_mol_weight,
                                     const InputProviderQ& qv,
                                     const InputProviderX& vmr,
                                     const view_1d<ScalarT>& mmr);

  //-----------------------------------------------------------------------------------------------//
  // Determines the vertical layer interface height from the vertical layer thicknesses:
  //   z_int = int_0^z(dz)
  // where
  //   dz             is the vertical layer thickness, [m]
  // Note: because this function does an integral it cannot be run just on a single level.  It requires
  // the full column wise integration.
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT, typename InputProviderZ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_z_int (const MemberType& team,
                               const int num_levs,
                               const InputProviderZ& dz,
                               const Real z_surf,
                               const view_1d<ScalarT>& z_int);

  

}; // struct PhysicsFunctions

} // namespace scream

#endif // SCREAM_COMMON_PHYSICS_FUNCTIONS_HPP

// We don't do ETI, since we don't know some of the concrete types.
// E.g., we don't know InputProvider, or ScalarT (although we could
// ETI the "common" cases, where the provider is a view_1d, and 
// Scalar=Real or Scalar=Pack<Real,SCREAM_PACK_SIZE>).
# include "scream_common_physics_functions_impl.hpp"
