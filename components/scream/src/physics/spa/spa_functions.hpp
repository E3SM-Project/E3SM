#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"
#include "ekat/mpi/ekat_comm.hpp"

#include <numeric>

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

  using gid_type = AbstractGrid::gid_type;

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

  template <int N>
  using view_Nd_host = typename KT::template view_ND<Real,N>::HostMirror;
  using view_1d_host = view_Nd_host<1>;
  /* ------------------------------------------------------------------------------------------- */
  // SPA structures to help manage all of the variables:
  struct SPATimeState {
    SPATimeState() = default;
    // Whether the timestate has been initialized.
    bool inited = false;
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

  struct SPAData {
    SPAData() = default;
    SPAData(const int ncol_, const int nlev_, const int nswbands_, const int nlwbands_) :
      ncols(ncol_)
      ,nlevs(nlev_)
      ,nswbands(nswbands_)
      ,nlwbands(nlwbands_)
    {
      hyam       = view_1d<Spack>("",nlevs);
      hybm       = view_1d<Spack>("",nlevs);
      PS         = view_1d<Real>("",ncols);
      CCN3       = view_2d<Spack>("",ncols,nlevs);
      AER_G_SW   = view_3d<Spack>("",ncols,nswbands,nlevs); 
      AER_SSA_SW = view_3d<Spack>("",ncols,nswbands,nlevs); 
      AER_TAU_SW = view_3d<Spack>("",ncols,nswbands,nlevs); 
      AER_TAU_LW = view_3d<Spack>("",ncols,nlwbands,nlevs); 
    }
    // Basic spatial dimensions of the data
    Int ncols;
    Int nlevs;
    Int nswbands;
    Int nlwbands;
    // Hybrid Coordinates
    view_1d<Spack> hyam, hybm;
    // PS unit = Pa: Surface pressure
    view_1d<Real> PS;
    // CCN3: CCN concentration at S=0.1%, units = #/cm3, dimensions = (ncol,nlev)
    view_2d<Spack> CCN3;
    // AER_G_SW unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_G_SW;
    // AER_SSA_SW unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_SSA_SW;
    // AER_TAU_SW unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_TAU_SW;
    // AER_TAU_LW unit = #/cm3, dimensions = (ncol,nswband=16,nlev)
    view_3d<Spack> AER_TAU_LW;

  }; // SPAPrescribedAerosolData

  struct SPAOutput {
    SPAOutput() = default;
    SPAOutput(const int ncol_, const int nlev_, const int nswbands_, const int nlwbands_) :
      ncols(ncol_)
      ,nlevs(nlev_)
      ,nswbands(nswbands_)
      ,nlwbands(nlwbands_)
    {
      CCN3       = view_2d<Spack>("",ncols,nlevs);
      AER_G_SW   = view_3d<Spack>("",ncols,nswbands,nlevs); 
      AER_SSA_SW = view_3d<Spack>("",ncols,nswbands,nlevs); 
      AER_TAU_SW = view_3d<Spack>("",ncols,nswbands,nlevs); 
      AER_TAU_LW = view_3d<Spack>("",ncols,nlwbands,nlevs); 
    }
    // Basic spatial dimensions of the data
    Int ncols;
    Int nlevs;
    Int nswbands;
    Int nlwbands;
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
    SPAHorizInterp(const int length_)
    {
      length = length_;
      weights = view_1d<Real>("",length_);
      source_grid_loc = view_1d<Int>("",length_);
      target_grid_loc = view_1d<Int>("",length_);
    }
    // Comm group used for SPA
    ekat::Comm m_comm;
    // Number of weights in remap data
    Int length;
    // Number of columns and levels on source grid
    // Note, the number of columns on target grid should already be known.
    //       these are given based on the simulation grid.
    Int source_grid_ncols;
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
    const view_2d<const Spack>& p_tgt,
    const SPAData&   data_beg,
    const SPAData&   data_end,
    const SPAOutput& data_out,
    Int ncols_scream,
    Int nlevs_scream,
    Int nswbands,
    Int nlwbands);

  static void get_remap_weights_from_file(
    const std::string&       remap_file_name,
    const Int                ncols_scream,
    gid_type                 min_dof,
    const view_1d<gid_type>& dofs_gids,
          SPAHorizInterp&    spa_horiz_interp);

  static void set_remap_weights_one_to_one(
    const Int                ncols_scream,
    gid_type                 min_dof,
    const view_1d<gid_type>& dofs_gids,
          SPAHorizInterp&    spa_horiz_interp);

  static void update_spa_data_from_file(
    const std::string&    spa_data_file_name,
    const Int             time_index,
    const Int             nswbands,
    const Int             nlwbands,
          SPAHorizInterp& spa_horiz_interp,
          SPAData&        spa_data);

  static void update_spa_timestate(
    const std::string&     spa_data_file_name,
    const Int              nswbands,
    const Int              nlwbands,
    const util::TimeStamp& ts,
          SPAHorizInterp&  spa_horiz_interp,
          SPATimeState&    time_state, 
          SPAData&         spa_beg,
          SPAData&         spa_end);
    

}; // struct Functions

} // namespace spa 
} // namespace scream

#endif // SPA_FUNCTIONS_HPP

// We don't do ETI, since we don't know some of the concrete types.
// E.g., we don't know InputProvider, or ScalarT (although we could
// ETI the "common" cases, where the provider is a view_1d, and
// Scalar=Real or Scalar=Pack<Real,SCREAM_PACK_SIZE>).
# include "spa_functions_impl.hpp"
