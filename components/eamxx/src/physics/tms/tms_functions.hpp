#ifndef TMS_FUNCTIONS_HPP
#define TMS_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/eamxx_types.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace tms {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for a process. Currently only one
 * function exists: compute_tms().
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
  static void compute_tms(
    const int&                   ncols,
    const int&                   nlevs,
    const view_3d<const Scalar>& horiz_wind,
    const view_2d<const Scalar>& t_mid,
    const view_2d<const Scalar>& p_mid,
    const view_2d<const Scalar>& exner,
    const view_2d<const Scalar>& z_mid,
    const view_1d<const Scalar>& sgh,
    const view_1d<const Scalar>& landfrac,
    const view_1d<Scalar>&       ksrf,
    const view_2d<Scalar>&       tau_tms);

}; // struct tms

} // namespace tms
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)  \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)

# include "compute_tms_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE

#endif // TMS_FUNCTIONS_HPP
