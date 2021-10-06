#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"
#include "ekat/mpi/ekat_comm.hpp"

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

  using WorkspaceManager = typename ekat::WorkspaceManager<Spack, Device>;
  using Workspace        = typename WorkspaceManager::Workspace;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  /* ------------------------------------------------------------------------------------------- */
  // SPA structures to help manage all of the variables:
  struct SPATimeState {
    SPATimeState() = default;
    // The current month
    Int current_month;
    // Julian Date for the beginning of the month, as defined in
    //           /src/share/util/scream_time_stamp.hpp
    // See this file for definition of Julian Date.
    Real t_beg_month;
    // Current simulation Julian Date
    Real t_now;
    // Number of days in the current month, cast as a Real
    Real days_this_month;
  }; // SPATimeState

  struct SPAPressureState {
    SPAPressureState() = default;
    // Number of horizontal columns and vertical levels for the data
    Int ncols;
    Int nlevs;
    // Surface pressure for data at the beginning of the month
    view_1d<const Real> ps_this_month;
    // Surface pressure for data at the beginning of next month
    view_1d<const Real> ps_next_month;
    // Hybrid coordinate values
    view_1d<const Spack> hyam, hybm;
    // Current simulation pressure levels
    view_2d<const Spack> pmid;
  }; // SPAPressureState

  struct SPAData {
    SPAData() = default;
    // Basic spatial dimensions of the data
    Int ncols;
    Int nlevs;
    // CCN3: CCN concentration at S=0.1%, units = #/cm3, dimensions = (ncol,nlev)
    view_2d<Spack> CCN3;
    // AER_Gnit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_G_SW;
    // AER_S unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_SSA_SW;
    // AER_T unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_TAU_SW;
    // AER_T unit = #/cm3, dimensions = (ncol,nswband=16,nlev)
    view_3d<Spack> AER_TAU_LW;
  }; // SPAPrescribedAerosolData

  struct SPAOutput {
    SPAOutput() = default;
    // CCN3: CCN concentration at S=0.1%, units = #/cm3, dimensions = (ncol,nlev)
    view_2d<Spack> CCN3;
    // AER_G_SW - 14 bands
    // AER_G_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_G_SW;
    // AER_SSA_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_SSA_SW;
    // AER_TAU_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_TAU_SW;
    // AER_TAU_LW: unit = #/cm3, dimensions = (ncol,nswband=16,nlev)
    view_3d<Spack> AER_TAU_LW;
  }; // SPAPrescribedAerosolData

  struct SPAHorizInterp {
    // This structure stores the information need by SPA to conduct horizontal
    // interpolation from a set of source data to horizontal locations in the
    // simulation grid.
    // The source_grid_loc stores the column index in the source data,
    // The target_grid_loc stores the column index in the target data that will be mapped to
    // The weights stores the remapping weight to be applied to the source grid data for this location
    //   in the target data.
    SPAHorizInterp() = default;
    // Number of weights in remap data
    Int length;
    // Number of columns and levels on source grid
    // Note, the number of columns and levels on target grid should already be known.
    //       these are given based on the simulation grid.
    Int source_grid_ncols, source_grid_nlevs;
    // 1D index of weights.  Needs decoder indexing, see below
    view_1d<Real> weights;
    // 1D index of source grid column.
    view_1d<Int> source_grid_loc;
    // 1D index of target grid column.
    view_1d<Int> target_grid_loc;
  }; // SPAHorizInterp
  /* ------------------------------------------------------------------------------------------- */
  // SPA routines
  static void spa_main(
    const SPATimeState& time_state,
    const SPAPressureState& pressure_state,
    const SPAData&   data_beg,
    const SPAData&   data_end,
    const SPAOutput& data_out,
    Int ncols_scream,
    Int nlevs_scream,
    Int nswbands,
    Int nlwbands);

  static void get_remap_weights_from_file(
    const std::string& remap_file_name,
          SPAHorizInterp& spa_horiz_interp,
    const Int ncols_scream);

}; // struct Functions

} // namespace spa 
} // namespace scream

#endif // SPA_FUNCTIONS_HPP

// We don't do ETI, since we don't know some of the concrete types.
// E.g., we don't know InputProvider, or ScalarT (although we could
// ETI the "common" cases, where the provider is a view_1d, and
// Scalar=Real or Scalar=Pack<Real,SCREAM_PACK_SIZE>).
# include "spa_functions_impl.hpp"
