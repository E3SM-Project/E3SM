#ifndef RAYLEIGH_FRICTION_FUNCTIONS_HPP
#define RAYLEIGH_FRICTION_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_kokkos_types.hpp>

namespace scream {
namespace rayleigh_friction {

/*
 * Functions is a stateless struct used to encapsulate the functions
 * of the Rayleigh friction scheme, ported from EAM (rayleigh_friction.F90).
 *
 * Rayleigh friction is applied in the region of the model top. A decay rate
 * profile is specified that is largest at the model top and drops off
 * vertically using a hyperbolic tangent profile:
 *
 *   otau(k) = otau0 * [1 + tanh(x)] / 2,   x = (rayk0 - k) / krange
 *
 * where k is the 1-based level index (k=1 is the model top). The tendencies
 * in u and v are computed with an Euler backward scheme, and the kinetic
 * energy lost is converted to heat (i.e., applied to the dry static energy).
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using KT = ekat::KokkosTypes<Device>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;

  //
  // --------- Functions ---------
  //

  // Compute the decay rate profile otau [1/s].
  //   nlevs     : number of vertical levels
  //   rayk0     : 1-based vertical level at which the profile is centered
  //   raykrange : width of the profile in levels (if 0, set to (rayk0-1)/2)
  //   raytau0   : approximate decay time at model top [days] (if 0, otau=0)
  static void rayleigh_friction_init(
    const int&             nlevs,
    const int&             rayk0,
    const Scalar&          raykrange,
    const Scalar&          raytau0,
    const view_1d<Scalar>& otau);

  // Apply Rayleigh friction to the horizontal winds and temperature in place.
  //   dt          : physics timestep [s]
  //   otau        : decay rate profile [1/s]
  //   horiz_winds : (ncols, 2, nlevs) horizontal winds [m/s]
  //   T_mid       : (ncols, nlevs) midpoint temperature [K]
  static void rayleigh_friction_tend(
    const int&                   ncols,
    const int&                   nlevs,
    const Scalar&                dt,
    const view_1d<const Scalar>& otau,
    const view_3d<Scalar>&       horiz_winds,
    const view_2d<Scalar>&       T_mid);

}; // struct Functions

} // namespace rayleigh_friction
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)  \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)

# include "rayleigh_friction_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE

#endif // RAYLEIGH_FRICTION_FUNCTIONS_HPP
