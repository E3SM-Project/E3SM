#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "control/intensive_observation_period.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"

#include <ekat/ekat_pack_utils.hpp>

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

  using Spack = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;

  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  using gid_type = AbstractGrid::gid_type;

  using iop_ptr_type = std::shared_ptr<control::IntensiveObservationPeriod>;

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

  struct IOPReader {
    IOPReader (iop_ptr_type& iop_,
               const std::string file_name_,
               const std::vector<Field>& io_fields_,
               const std::shared_ptr<const AbstractGrid>& io_grid_)
      : iop(iop_), file_name(file_name_)
    {
      field_mgr = std::make_shared<FieldManager>(io_grid_);
      for (auto& f : io_fields_) {
        field_mgr->add_field(f);
        field_names.push_back(f.name());
      }

      // Set IO info for this grid and file in IOP object
      iop->setup_io_info(file_name, io_grid_);
    }

    void read_variables(const int time_index, const util::TimeStamp& ts) {
      iop->read_fields_from_file_for_iop(file_name, field_names, ts, field_mgr, time_index);
    }

    iop_ptr_type iop;
    std::string file_name;
    std::vector<std::string> field_names;
    std::shared_ptr<FieldManager> field_mgr;
  };

  // The output is really just SPAData, but for clarity it might
  // help to see a SPAOutput along a SPAInput in functions signatures
  using SPAOutput = SPAData;

  /* ------------------------------------------------------------------------------------------- */
  // SPA routines

  static std::shared_ptr<AbstractRemapper>
  create_horiz_remapper (
      const std::shared_ptr<const AbstractGrid>& model_grid,
      const std::string& spa_data_file,
      const std::string& map_file,
      const bool use_iop = false);

  static std::shared_ptr<AtmosphereInput>
  create_spa_data_reader (
      const std::shared_ptr<AbstractRemapper>& horiz_remapper,
      const std::string& spa_data_file);

  static std::shared_ptr<IOPReader>
  create_spa_data_reader (
      iop_ptr_type& iop,
      const std::shared_ptr<AbstractRemapper>& horiz_remapper,
      const std::string& spa_data_file);

  static void spa_main(
    const SPATimeState& time_state,
    const view_2d<const Spack>& p_tgt,
    const view_2d<      Spack>& p_src,  // Temporary
    const SPAInput&   data_beg,
    const SPAInput&   data_end,
    const SPAInput&   data_tmp,         // Temporary
    const SPAOutput&  data_out);

  static void update_spa_data_from_file(
    std::shared_ptr<AtmosphereInput>& scorpio_reader,
    std::shared_ptr<IOPReader>&       iop_reader,
    const util::TimeStamp&            ts,
    const int                         time_index, // zero-based
    AbstractRemapper&                 spa_horiz_interp,
    SPAInput&                         spa_input);

  static void update_spa_timestate(
    std::shared_ptr<AtmosphereInput>& scorpio_reader,
    std::shared_ptr<IOPReader>&       iop_reader,
    const util::TimeStamp&            ts,
    AbstractRemapper&                 spa_horiz_interp,
    SPATimeState&                     time_state,
    SPAInput&                         spa_beg,
    SPAInput&                         spa_end);

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
