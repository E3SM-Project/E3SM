#ifndef EAMXX_MAM_OPTICS_HPP
#define EAMXX_MAM_OPTICS_HPP

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>
#include <physics/mam/mam_aerosol_optics_read_tables.hpp>
#include <physics/mam/mam_coupling.hpp>
#include <share/atm_process/ATMBufferManager.hpp>
#include <share/atm_process/atmosphere_process.hpp>
#include <share/util/eamxx_common_physics_functions.hpp>
#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected_except_cuda public
#define private_except_cuda public
#else
#define protected_except_cuda protected
#define private_except_cuda private
#endif

namespace scream {

// The process responsible for handling MAM4 aerosol optical properties. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMOptics final : public scream::AtmosphereProcess {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // a quantity stored in a single vertical column with a single index
  using ColumnView = mam4::ColumnView;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

 public:
  // Constructor
  MAMOptics(const ekat::Comm &comm, const ekat::ParameterList &params);

  protected_except_cuda :

      // --------------------------------------------------------------------------
      // AtmosphereProcess overrides (see
      // share/atm_process/atmosphere_process.hpp)
      // --------------------------------------------------------------------------

      // process metadata
      AtmosphereProcessType
      type() const override;
  std::string name() const override;

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  private_except_cuda :
      // FIXME: duplicate code from microphysics: ask it can be moved to place
      // where other process can see it. Atmosphere processes often have a
      // pre-processing step that constructs required variables from the set of
      // fields stored in the field manager. This functor implements this step,
      // which is called during run_impl.
      struct Preprocess {
    Preprocess() = default;

    // on host: initializes preprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere &wet_atm,
                    const mam_coupling::AerosolState &wet_aero,
                    const mam_coupling::DryAtmosphere &dry_atm,
                    const mam_coupling::AerosolState &dry_aero) {
      ncol_     = ncol;
      nlev_     = nlev;
      wet_atm_  = wet_atm;
      wet_aero_ = wet_aero;
      dry_atm_  = dry_atm;
      dry_aero_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) const {
      const int i = team.league_rank();  // column index
      // first, compute dry fields
      compute_dry_mixing_ratios(team, wet_atm_, dry_atm_, i);
      compute_dry_mixing_ratios(team, wet_atm_, wet_aero_, dry_aero_, i);
      team.team_barrier();
      // second, we can use dry fields to compute dz, zmin, zint
      compute_vertical_layer_heights(team, dry_atm_, i);
      compute_updraft_velocities(team, wet_atm_, dry_atm_, i);
    }  // operator()

    // number of horizontal columns and vertical levels
    int ncol_, nlev_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_;
    mam_coupling::DryAtmosphere dry_atm_;
    mam_coupling::AerosolState wet_aero_, dry_aero_;

  };  // MAMMicrophysics::Preprocess

  // Postprocessing functor
  struct Postprocess {
    Postprocess() = default;

    // on host: initializes postprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere &wet_atm,
                    const mam_coupling::AerosolState &wet_aero,
                    const mam_coupling::DryAtmosphere &dry_atm,
                    const mam_coupling::AerosolState &dry_aero) {
      ncol_     = ncol;
      nlev_     = nlev;
      wet_atm_  = wet_atm;
      wet_aero_ = wet_aero;
      dry_atm_  = dry_atm;
      dry_aero_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) const {
      const int i = team.league_rank();  // column index
      compute_wet_mixing_ratios(team, dry_atm_, dry_aero_, wet_aero_, i);
    }  // operator()

    // number of horizontal columns and vertical levels
    int ncol_, nlev_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_;
    mam_coupling::DryAtmosphere dry_atm_;
    mam_coupling::AerosolState wet_aero_, dry_aero_;
  };  // MAMMicrophysics::Postprocess

  // pre- and postprocessing scratch pads (for wet <-> dry conversions)
  Preprocess preprocess_;
  Postprocess postprocess_;

  // state variable
  // mam_coupling::view_3d state_q_,  qqcw_;// odap_aer_,

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // number of shortwave and longwave radiation bands
  int nswbands_, nlwbands_;

  // FIXME: move these values to mam_coupling
  mam_coupling::const_view_2d p_int_, p_del_;

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // atmospheric and aerosol state variables
  // atmospheric and aerosol state variables
  mam_coupling::WetAtmosphere wet_atm_;
  mam_coupling::DryAtmosphere dry_atm_;
  mam_coupling::AerosolState wet_aero_, dry_aero_;

  mam_coupling::view_3d ssa_cmip6_sw_, af_cmip6_sw_, ext_cmip6_sw_;
  // long wave extinction in the units of [1/km]
  mam_coupling::view_3d ext_cmip6_lw_;
  mam4::modal_aer_opt::AerosolOpticsDeviceData aerosol_optics_device_data_;
  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
  mam_coupling::view_2d work_;
  mam_coupling::view_3d tau_ssa_g_sw_, tau_ssa_sw_, tau_sw_, tau_f_sw_;
  // Mapping from old RRTMG sw bands to new band ordering in RRTMGP
  //  given old index swband (RRTMG) return new index swband RRTMGP
  mam_coupling::view_int_1d get_idx_rrtmgp_from_rrtmg_swbands_;

  mam_coupling::Buffer buffer_;
};  // MAMOptics

}  // namespace scream

#endif  // EAMXX_MAM_OPTICS_HPP
