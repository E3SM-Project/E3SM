#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/horizontal_remap_utility.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_workspace.hpp"
#include "ekat/mpi/ekat_comm.hpp"

namespace scream {
namespace spa {

template <typename ScalarType, typename DeviceType>
struct SPAFunctions
{

  //
  // ------- Types --------
  //

  using Scalar = ScalarType;
  using Device = DeviceType;

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

  template <typename S, int N>
  using view_Nd_host = typename KT::template view_ND<S,N>::HostMirror;

  template <typename S>
  using view_1d_host = view_Nd_host<S,1>;
  /* ------------------------------------------------------------------------------------------- */
  // SPA structures to help manage all of the variables:
  struct SPATimeState {
    SPATimeState() = default;
    // Whether the timestate has been initialized.
    bool inited = false;
    // The current month
    int current_month = -1;
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
    SPAData(const int ncol_, const int nlev_, const int nswbands_, const int nlwbands_)
    {
      init (ncol_,nlev_,nswbands_,nlwbands_,true);
    }

    void init(const int ncol_, const int nlev_, const int nswbands_, const int nlwbands_, const bool allocate)
    {
      ncols = ncol_;
      nlevs = nlev_;
      nswbands = nswbands_;
      nlwbands = nlwbands_;

      if (allocate) {
        const int npacks = ekat::PackInfo<Spack::n>::num_packs(nlevs);

        CCN3       = view_2d<Spack>("",ncols,npacks);
        AER_G_SW   = view_3d<Spack>("",ncols,nswbands,npacks);
        AER_SSA_SW = view_3d<Spack>("",ncols,nswbands,npacks);
        AER_TAU_SW = view_3d<Spack>("",ncols,nswbands,npacks);
        AER_TAU_LW = view_3d<Spack>("",ncols,nlwbands,npacks);
      }
    }

    // Basic spatial dimensions of the data
    int ncols;
    int nlevs;
    int nswbands;
    int nlwbands;

    view_2d<Spack> CCN3;        // CCN concentration at S=0.1: units = #/cm3, dimensions = (ncols,nlevs)
    view_3d<Spack> AER_G_SW;    // AER_G_SW:   units = #/cm3, dimensions = (ncols,nswbands,nlevs)
    view_3d<Spack> AER_SSA_SW;  // AER_SSA_SW: units = #/cm3, dimensions = (ncols,nswbands,nlevs)
    view_3d<Spack> AER_TAU_SW;  // AER_TAU_SW: units = #/cm3, dimensions = (ncols,nswbands,nlevs)
    view_3d<Spack> AER_TAU_LW;  // AER_TAU_LW: units = #/cm3, dimensions = (ncols,nlwbands,nlevs)
  }; // SPAData

  struct SPAInput {
    SPAInput() = default;
    SPAInput(const int ncols_, const int nlevs_, const int nswbands_, const int nlwbands_)
    {
      init(ncols_,nlevs_,nswbands_,nlwbands_);
    }

    void init(const int ncols_, const int nlevs_, const int nswbands_, const int nlwbands_)
    {
      const int npacks = ekat::PackInfo<Spack::n>::num_packs(nlevs_);

      data.init(ncols_,nlevs_,nswbands_,nlwbands_,true);
      PS  = view_1d<Real>("", ncols_);
      hyam = view_1d<Spack>("", npacks);
      hybm = view_1d<Spack>("", npacks);
    }

    view_1d<Spack>  hyam, hybm;   // Hybrid Coordinates
    view_1d<Real>   PS;           // PS unit = Pa: Surface pressure

    SPAData         data;         // All spa fields
  }; // SPAInput

  // The output is really just SPAData, but for clarity it might
  // help to see a SPAOutput along a SPAInput in functions signatures
  using SPAOutput = SPAData;

  struct SPAHorizInterp {
    // This structure stores the information need by SPA to conduct horizontal
    // interpolation from a set of source data to horizontal locations in the
    // simulation grid.
    // The source_grid_loc stores the column index in the source data,
    // The target_grid_loc stores the column index in the target data that will be mapped to
    // The weights stores the remapping weight to be applied to the source grid data for this location
    //   in the target data.
    SPAHorizInterp() = default;
    explicit SPAHorizInterp(const ekat::Comm& comm)
    {
      m_comm = comm;
    }
    // Horizontal Remap
    HorizontalMap horiz_map;
    // Comm group used for SPA
    ekat::Comm m_comm;

  }; // SPAHorizInterp
  /* ------------------------------------------------------------------------------------------- */
  // SPA routines
  static void spa_main(
    const SPATimeState& time_state,
    const view_2d<const Spack>& p_tgt,
    const view_2d<      Spack>& p_src,  // Temporary
    const SPAInput&   data_beg,
    const SPAInput&   data_end,
    const SPAInput&   data_tmp,         // Temporary
    const SPAOutput&  data_out);

  static void get_remap_weights_from_file(
    const std::string&             remap_file_name,
    const gid_type                 min_dof,
    const view_1d<const gid_type>& dofs_gids,
          SPAHorizInterp&          spa_horiz_interp);

  static void set_remap_weights_one_to_one(
    gid_type                       min_dof,
    const view_1d<const gid_type>& dofs_gids,
          SPAHorizInterp&          spa_horiz_interp);

  static void update_spa_data_from_file(
    const std::string&    spa_data_file_name,
    const int             time_index,
    const int             nswbands,
    const int             nlwbands,
          SPAHorizInterp& spa_horiz_interp,
          SPAInput&       spa_data);

  static void update_spa_timestate(
    const std::string&     spa_data_file_name,
    const int              nswbands,
    const int              nlwbands,
    const util::TimeStamp& ts,
          SPAHorizInterp&  spa_horiz_interp,
          SPATimeState&    time_state,
          SPAInput&        spa_beg,
          SPAInput&        spa_end);

  // The following three are called during spa_main
  static void perform_time_interpolation (
      const SPATimeState& time_state,
      const SPAInput&  data_beg,
      const SPAInput&  data_end,
      const SPAInput&  data_out);

  static void compute_source_pressure_levels (
      const view_1d<const Real>& ps_src,
      const view_2d<      Spack>& p_src,
      const view_1d<const Spack>& hyam,
      const view_1d<const Spack>& hybm);

  static void perform_vertical_interpolation (
      const view_2d<const Spack>& p_src,
      const view_2d<const Spack>& p_tgt,
      const SPAData&  data_in,
      const SPAData&  data_out);

  // Return the subcolumn of the proper variable, where ivar
  // is a condensed idx for var and possibly band. In particular:
  //  - ivar=0: return CCN
  //  - else, jvar = ivar-1, and
  //     - jvar=[0...nswbands): return aer_g_sw
  //     - jvar=[nswbands...2*nswbands): return aer_ssa_sw
  //     - jvar=[2*nswbands...3*nswbands): return aer_tau_sw
  //     - jvar=[3*nswbands...3*nswbands+nlbands): return aer_tau_lw
  KOKKOS_INLINE_FUNCTION
  static view_1d<Spack> get_var_column (const SPAData& data, const int icol, const int ivar);

  // Performs convex interpolation of x0 and x1 at point t
  template<typename ScalarX,typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarX linear_interp(const ScalarX& x0, const ScalarX& x1, const ScalarT& t);

}; // struct Functions

} // namespace spa
} // namespace scream

#endif // SPA_FUNCTIONS_HPP

// We don't do ETI, since we don't know some of the concrete types.
// E.g., we don't know InputProvider, or ScalarType (although we could
// ETI the "common" cases, where the provider is a view_1d, and
// Scalar=Real or Scalar=Pack<Real,SCREAM_PACK_SIZE>).
# include "spa_functions_impl.hpp"
