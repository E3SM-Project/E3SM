#ifndef EAMXX_MAM_MICROPHYSICS_HPP
#define EAMXX_MAM_MICROPHYSICS_HPP

#include <physics/mam/mam_coupling.hpp>
#include <share/atm_process/atmosphere_process.hpp>
#include <share/util/scream_common_physics_functions.hpp>

#include "impl/mam4_amicphys.cpp" // mam4xx top-level microphysics function(s)

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected_except_cuda public
#define private_except_cuda public
#else
#define protected_except_cuda protected
#define private_except_cuda private
#endif

namespace scream
{

// The process responsible for handling MAM4 aerosol microphysics. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMMicrophysics final : public scream::AtmosphereProcess {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // views for single- and multi-column data
  using view_1d_int   = typename KT::template view_1d<int>;
  using view_1d       = typename KT::template view_1d<Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using const_view_1d = typename KT::template view_1d<const Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  // unmanaged views (for buffer and workspace manager)
  using uview_1d = Unmanaged<typename KT::template view_1d<Real>>;
  using uview_2d = Unmanaged<typename KT::template view_2d<Real>>;

  // a quantity stored in a single vertical column with a single index
  using ColumnView = mam4::ColumnView;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

public:

  // Constructor
  MAMMicrophysics(const ekat::Comm& comm, const ekat::ParameterList& params);

protected_except_cuda:

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // set aerosol microphysics configuration parameters (called by constructor)
  void configure(const ekat::ParameterList& params);

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  // performs some checks on the tracers group
  void set_computed_group_impl(const FieldGroup& group) override;

private_except_cuda:

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // configuration data (for the moment, we plan to be able to move this to
  // the device, so we can't use C++ strings)
  struct Config {
    // photolysis parameters
    struct {
      char rsf_file[MAX_FILENAME_LEN];
      char xs_long_file[MAX_FILENAME_LEN];
    } photolysis;

    // stratospheric chemistry parameters
    struct {
      int o3_lbl; // number of layers with ozone decay from the surface
      int o3_sfc; // set from namelist input linoz_sfc
      int o3_tau; // set from namelist input linoz_tau
      Real psc_T; // set from namelist input linoz_psc_T
      char chlorine_loading_file[MAX_FILENAME_LEN];
    } linoz;

    // aqueous chemistry parameters
    mam4::mo_setsox::Config setsox;

    // aero microphysics configuration (see impl/mam4_amicphys.cpp)
    impl::AmicPhysConfig amicphys;

    // dry deposition parameters
    struct {
      char srf_file[MAX_FILENAME_LEN];
    } drydep;
  };
  Config config_;

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct Preprocess {
    Preprocess() = default;

    // on host: initializes preprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere& wet_atm,
                    const mam_coupling::AerosolState& wet_aero,
                    const mam_coupling::DryAtmosphere& dry_atm,
                    const mam_coupling::AerosolState& dry_aero) {
      ncol_ = ncol;
      nlev_ = nlev;
      wet_atm_ = wet_atm;
      wet_aero_ = wet_aero;
      dry_atm_ = dry_atm;
      dry_aero_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank(); // column index

      compute_vertical_layer_heights(team, dry_atm_, i);
      team.team_barrier(); // allows kernels below to use layer heights
      compute_updraft_velocities(team, wet_atm_, dry_atm_, i);
      compute_dry_mixing_ratios(team, wet_atm_, dry_atm_, i);
      compute_dry_mixing_ratios(team, wet_atm_, wet_aero_, dry_aero_, i);
      team.team_barrier();
    } // operator()

    // number of horizontal columns and vertical levels
    int ncol_, nlev_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_;
    mam_coupling::DryAtmosphere dry_atm_;
    mam_coupling::AerosolState  wet_aero_, dry_aero_;

  }; // MAMMicrophysics::Preprocess

  // Postprocessing functor
  struct Postprocess {
    Postprocess() = default;

    // on host: initializes postprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere& wet_atm,
                    const mam_coupling::AerosolState& wet_aero,
                    const mam_coupling::DryAtmosphere& dry_atm,
                    const mam_coupling::AerosolState& dry_aero) {
      ncol_ = ncol;
      nlev_ = nlev;
      wet_atm_ = wet_atm;
      wet_aero_ = wet_aero;
      dry_atm_ = dry_atm;
      dry_aero_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank(); // column index
      compute_wet_mixing_ratios(team, dry_atm_, dry_aero_, wet_aero_, i);
      team.team_barrier();
    } // operator()

    // number of horizontal columns and vertical levels
    int ncol_, nlev_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_;
    mam_coupling::DryAtmosphere dry_atm_;
    mam_coupling::AerosolState  wet_aero_, dry_aero_;
  }; // MAMMicrophysics::Postprocess

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // pre- and postprocessing scratch pads (for wet <-> dry conversions)
  Preprocess preprocess_;
  Postprocess postprocess_;

  // atmospheric and aerosol state variables
  mam_coupling::WetAtmosphere wet_atm_;
  mam_coupling::DryAtmosphere dry_atm_;
  mam_coupling::AerosolState  wet_aero_, dry_aero_;

  // photolysis rate table (column-independent)
  mam4::mo_photo::PhotoTableData photo_table_;

  // column areas, latitudes, longitudes
  const_view_1d col_areas_, col_latitudes_, col_longitudes_;

  // time step number
  int step_;

  // workspace manager for internal local variables
  //ekat::WorkspaceManager<Real, KT::Device> workspace_mgr_;
  mam_coupling::Buffer buffer_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

  // sets defaults for "namelist parameters"
  void set_defaults_();

}; // MAMMicrophysics

} // namespace scream

#endif // EAMXX_MAM_MICROPHYSICS_HPP
