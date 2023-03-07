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
  // Determine the length of a square column cell at a specifc latitude given the area in radians
  //   grid_dx = mpdeglat * area
  // where,
  //   mpdeglat is the distance between two points on an ellipsoid
  //   area     is the area of the column cell in radians.
  //   lat      is the latitude of the grid column in radians.
  // NOTE - Here we assume that the column area is a SQUARE so dx=dy.  We will need a different
  //        routine for a rectangular (or other shape) area.
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_dx_from_area(const ScalarT& area, const ScalarT& lat);

  //-----------------------------------------------------------------------------------------------//
  // Determines the density given the definition of pseudo_density passed by the dycore
  //   rho = pseudo_density/dz/g
  // where,
  //   rho            is the density of air, [kg/m3]
  //   pseudo_density is the pressure level thickness given a shallow atmosphere, [Pa]
  //   dz             is the geopotential thickness of the layer, [m]
  //   g              is the gravitational constant, [m/s2] - defined in physics_constants.hpp
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_density(const ScalarT& pseudo_density, const ScalarT& dz);

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
  // Compute temperature from dry static energy (DSE).
  //   temperature = (DSE - g*z - surf_geopotential)/Cp
  // where
  //   Cp                is the heat constant of air at constant pressure [J/kg]
  //   g                 is the gravitational constant [m s-2]
  //   DSE               is the dry static energy.  Units in [J/kg]
  //   z                 is the geopotential height above surface at midpoints. Units in [m].
  //   surf_geopotential is the surface geopotential height. Units in [m].
  //   temperature       is the atmospheric temperature. Units in [K].
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_temperature_from_dse(const ScalarT& dse, const ScalarT& z, const Real surf_geopotential);

  //-----------------------------------------------------------------------------------------------//
  // Computes drymmr (mass of a constituent divided by mass of dry air; commonly known as mixing ratio)
  // for any wetmmr constituent (mass of a constituent divided by mass of dry air plus water
  // vapor) using qv_wet (mass of water vapor divided by mass of dry air plus
  // water vapor; see specific humidity):
  //   drymmr = wetmmr / (1 - qv_wet)
  // where
  //   drymmr         is the dry mass mixing ratio of a species
  //   wetmmr         is the wet mass mixing ratio of a species
  //   qv_wet         is water vapor wet mass mixing ratio
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_drymmr_from_wetmmr(const ScalarT& wetmmr, const ScalarT& qv_wet);

  //-----------------------------------------------------------------------------------------------//
  // Computes wetmmr (mass of a constituent divided by mass of dry air plus water vapor)
  // for any drymmr constituent (mass of a constituent divided by mass of dry air;
  // commonly known as mixing ratio) using qv_dry (mass of water vapor divided by mass
  // of dry air):
  //   wetmmr = drymmr / (1 + qv_dry)
  // where
  //   wetmmr         is the wet mass mixing ratio of a species
  //   drymmr         is the dry mass mixing ratio of a species
  //   qv_dry         is specific humidity of water vapor
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT calculate_wetmmr_from_drymmr(const ScalarT& drymmr, const ScalarT& qv_dry);

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

  //-----------------------------------------------------------------------------------------------
  // Calculate T at the bottom of the grid cell closest to the surface for use in PSL computation.
  // This is done assuming a 6.5 K/km lapse rate, which is a horrible assumption but avoids problems that
  // computing lapse rate via extrapolation might produce strange answers. It is also what CESM has done
  // for the last 20 yrs so seems to be sufficient. Don't assume this method is appropriate for any other use.
  // INPUTS:
  // T_mid_bot
  // z_mid_bot
  // RETURNS:
  // T at the bottom of the cell nearest the surface (K)
  //-----------------------------------------------------------------------------------------------
  KOKKOS_INLINE_FUNCTION
  static Real calculate_surface_air_T(const Real& T_mid_bot, const Real& z_mid_bot);

  //-----------------------------------------------------------------------------------------------//
  // Compute the lapse rate and effective ground temperature for use in calculating psl. This function should only
  // be used by calculate_psl.
  // INPUTS:
  // T_ground is the air temperature at the bottom of the cell closest to the surface (aka T_int[nlev+1]; K)
  // phi_ground is the geopotential at surface (aka surf_geopotential; m2/s2)
  // OUTPUTS:
  // lapse (K/m) is the lapse rate
  // T_ground_tmp is the effective ground temperature (K)
  //-----------------------------------------------------------------------------------------------//
  KOKKOS_INLINE_FUNCTION
  static void lapse_T_for_psl(const Real& T_ground, const Real& phi_ground,
                              Real& lapse, Real& T_ground_tmp );

  //-----------------------------------------------------------------------------------------------//
  // Calculate sea level pressure assuming dry air between ground and sea level and using a lapse
  // rate of 6.5K/km except in very warm conditions. See docs/tech_doc/physics/psl/psl_doc.tex for details
  // INPUTS:
  // T_ground is the air temperature at the bottom of the cell closest to the surface (aka T_int[nlev+1]; K)
  // p_ground is the pressure at the bottom of the cell closest to the surface (Pa)
  // phi_ground is the geopotential at surface (aka surf_geopotential; m2/s2)
  // OUTPUTS:
  // psl is the sea level pressure (Pa)
  //-----------------------------------------------------------------------------------------------//
  KOKKOS_INLINE_FUNCTION
  static Real calculate_psl(const Real& T_ground, const Real& p_ground, const Real& phi_ground);

  //-----------------------------------------------------------------------------------------------//
  // Apply rayleigh friction. Given the decay rate profile, we compute the tendencies in u
  // and v components of the horizontal wind using an Euler backward scheme, and then apply
  // the negative of the kinetic energy tendency to the dry static energy.
  // Note: We don't actually calculate dse since this is simply a tendancy of cp*T_mid.
  // INPUTS:
  // dt is the physics timestep
  // otau is the decay rate
  // INPUT/OUTPUTS:
  // u_wind is u component of the horizontal wind (m/s)
  // v_wind is v component of the horizontal wind (m/s)
  // T_mid is the atmospheric temperature at the midpoints [K]
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static void apply_rayleigh_friction(const Real dt, const ScalarT& otau,
                                      ScalarT& u_wind, ScalarT& v_wind, ScalarT& T_mid);

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
  // parallel_for using a TeamVectorRange policy with no third layer  //
  // of parallelism (no ThreadVectorRange for loop).                  //
  // The routines are templated on the type of their inputs, so that  //
  // these can be 1d views as well as other routines (e.g. lambdas).  //
  // The only requirement is that the provider supports op()(int k),  //
  // which should return scalar corresponding to level k (which can   //
  // be a physical level or a "packed level", depending on the scalar //
  // type.                                                            //
  //                                                                  //
  // Most of these routines simply call the homonymous routine that   //
  // act on a single scalar, except for calculate_z_int/mid, which    //
  // only makes sense for a whole column (i.e., there is not          //
  // single-scalar version of compute_z_int/mid).                     //
  // ---------------------------------------------------------------- //

  using Device = DeviceT;
  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template<typename ScalarT, typename MT = Kokkos::MemoryManaged>
  using view_1d = typename KT::template view_1d<ScalarT, MT>;

  template<typename ScalarT, typename InputProviderP, typename InputProviderZ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_density (const MemberType& team,
                                 const InputProviderP& pseudo_density,
                                 const InputProviderZ& dz,
                                 const view_1d<ScalarT>& density);

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

  template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_temperature_from_dse (const MemberType& team,
                                              const InputProviderT& dse,
                                              const InputProviderZ& z,
                                              const Real surf_geopotential,
                                              const view_1d<ScalarT>& temperature);

  template<typename ScalarT, typename InputProviderX, typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_wetmmr_from_drymmr (const MemberType& team,
                             const InputProviderX& drymmr,
                             const InputProviderQ& qv,
                             const view_1d<ScalarT>& wetmmr);

  template<typename ScalarT, typename InputProviderX, typename InputProviderQ>
  KOKKOS_INLINE_FUNCTION
  static void calculate_drymmr_from_wetmmr (const MemberType& team,
                             const InputProviderX& wetmmr,
                             const InputProviderQ& qv_wet,
                             const view_1d<ScalarT>& drymmr);

  template<typename ScalarT,
           typename InputProviderPD, typename InputProviderP,
           typename InputProviderT,  typename InputProviderQ,
           typename MT = Kokkos::MemoryManaged>
  KOKKOS_INLINE_FUNCTION
  static void calculate_dz (const MemberType& team,
                            const InputProviderPD& pseudo_density,
                            const InputProviderP& p_mid,
                            const InputProviderT& T_mid,
                            const InputProviderQ& qv,
                            const view_1d<ScalarT, MT>& dz);

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

  template<typename ScalarT, typename InputProviderOtau, typename MT = Kokkos::MemoryManaged>
  KOKKOS_INLINE_FUNCTION
  static void apply_rayleigh_friction (const MemberType& team,
                                       const Real dt,
                                       const InputProviderOtau& otau,
                                       const view_1d<ScalarT, MT>& u_wind,
                                       const view_1d<ScalarT, MT>& v_wind,
                                       const view_1d<ScalarT, MT>& T_mid);

  //-----------------------------------------------------------------------------------------------//
  // Determines the vertical layer interface height from the vertical layer thicknesses:
  //   z_int = int_0^z(dz)
  // where
  //   dz             is the vertical layer thickness, [m]
  // Note: because this function does an integral it cannot be run just on a single level.  It requires
  // the full column wise integration.
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT, typename InputProviderZ, typename MT = Kokkos::MemoryManaged>
  KOKKOS_INLINE_FUNCTION
  static void calculate_z_int (const MemberType& team,
                               const int num_levs,
                               const InputProviderZ& dz,
                               const Real z_surf,
                               const view_1d<ScalarT, MT>& z_int);

  //-----------------------------------------------------------------------------------------------//
  // Determines the vertical layer height on mid points from the vertical layer interface height:
  //   z_mid(i,k) = (z_int(i,k) + z_int(i,k+1))/2.0
  // where
  //   z_int is the vertical layer interface height, [m]
  //-----------------------------------------------------------------------------------------------//
  template<typename ScalarT, typename InputProviderZ, typename MT = Kokkos::MemoryManaged>
  KOKKOS_INLINE_FUNCTION
  static void calculate_z_mid (const MemberType& team,
                               const int num_levs,
                               const InputProviderZ& z_int,
                               const view_1d<ScalarT, MT>& z_mid);

}; // struct PhysicsFunctions

} // namespace scream

#endif // SCREAM_COMMON_PHYSICS_FUNCTIONS_HPP

// We don't do ETI, since we don't know some of the concrete types.
// E.g., we don't know InputProvider, or ScalarT (although we could
// ETI the "common" cases, where the provider is a view_1d, and
// Scalar=Real or Scalar=Pack<Real,SCREAM_PACK_SIZE>).
# include "scream_common_physics_functions_impl.hpp"
