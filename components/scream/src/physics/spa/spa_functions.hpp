#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace spa {

template <typename ScalarT, typename DeviceT>
struct SPAFunctions
{

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  /* ------------------------------------------------------------------------------------------- */
  // SPA structures to help manage all of the variables:
  struct SPATimeState {
    SPATimeState() = default;
    // The current month
    Int current_month;
    // Julian Date for the beginning of the month
    Real t_beg_month;
    // Number of days in the current month, cast as a Real
    Real days_this_month;
  }; // SPATimeState

  struct SPAPressureState {
    SPAPressureState() = default;
    // Number of vertical levels for the data
    Int nlevs;
    // Surface pressure for data at the beginning of the month
    view_1d<const Real> ps_this_month;
    // Surface pressure for data at the beginning of next month
    view_1d<const Real> ps_next_month;
    // Hybrid coordinate values
    view_1d<const Spack> hyam, hybm;
  }; // SPAPressureState

  struct SPAData {
    SPAData() = default;
    // CCN3
    view_2d<const Spack> CCN3;
    // AER_G_SW - 14 bands
    view_3d<const Spack> AER_G_SW;
    // AER_SSA_SW - 14 bands
    view_3d<const Spack> AER_SSA_SW;
    // AER_TAU_LW - 16 bands
    view_3d<const Spack> AER_TAU_LW;
    // AER_TAU_SW - 14 bands
    view_3d<const Spack> AER_TAU_SW;
  }; // SPAPrescribedAerosolData
  /* ------------------------------------------------------------------------------------------- */
  // SPA routines
  KOKKOS_FUNCTION
  static void reconstruct_pressure_profile(
    const Int ncols,
    const Int nlevs,
    const view_1d<const Spack>& hya,
    const view_1d<const Spack>& hyb,
    const view_1d<const Real>&  PS,
    const view_2d<Spack>&       pres);

  KOKKOS_FUNCTION
  static void aero_vertical_remap(
    const Int ncols,
    const Int nlevs_src,
    const Int nlevs_tgt,
    const view_2d<const Spack>& pres_src,
    const view_2d<const Spack>& pres_tgt,
    const view_2d<const Spack>& aero_src,
    const view_2d<Spack>& aero_tgt); 

  // TODO: This function should really be templated to work with Scalars and Packed views
  KOKKOS_FUNCTION
  static void aero_time_interp(
    const Real& t0,
    const Real& ts,
    const Real& tlen,
    const Real& y0,
    const Real& y1,
          Real& y_out);

}; // struct Functions

} // namespace spa 
} // namespace scream

#endif // SPA_FUNCTIONS_HPP
