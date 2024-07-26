#ifndef EAMXX_MAM_SRF_ONLINE_EMISS_HPP
#define EAMXX_MAM_SRF_ONLINE_EMISS_HPP

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>
#include <physics/mam/mam_emissions_utils.hpp>
// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>
#include <physics/mam/srf_emission.hpp>
#include <share/atm_process/ATMBufferManager.hpp>
// For declaring surface and online emission class derived from atm process
// class
#include <share/atm_process/atmosphere_process.hpp>
// #include <share/util/scream_common_physics_functions.hpp>
#include <string>

namespace scream {

// The process responsible for handling MAM4 surface and online emissions. The
// AD stores exactly ONE instance of this class in its list of subcomponents.
class MAMSrfOnlineEmiss final : public scream::AtmosphereProcess {
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // Wet and dry states of atmosphere
  mam_coupling::WetAtmosphere wet_atm_;
  mam_coupling::DryAtmosphere dry_atm_;

  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;

  // buffer for sotring temporary variables
  mam_coupling::Buffer buffer_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

 public:
  using srfEmissFunc = mam_coupling::srfEmissFunctions<Real, DefaultDevice>;

  template <typename ScalarT>
  using uview_2d = Unmanaged<typename KT::template view_2d<ScalarT>>;

  // Constructor
  MAMSrfOnlineEmiss(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The type of subcomponent
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name() const { return "mam_srf_online_emissions"; }

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // Initialize variables
  void initialize_impl(const RunType run_type) override;

  // Run the process by one time step
  void run_impl(const double dt) override;

  // Finalize
  void finalize_impl(){/*Do nothing*/};
  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct Preprocess {
    Preprocess() = default;
    // on host: initializes preprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere &wet_atm,
                    const mam_coupling::AerosolState &wet_aero,
                    const mam_coupling::DryAtmosphere &dry_atm,
                    const mam_coupling::AerosolState &dry_aero) {
      ncol_pre_     = ncol;
      nlev_pre_     = nlev;
      wet_atm_pre_  = wet_atm;
      wet_aero_pre_ = wet_aero;
      dry_atm_pre_  = dry_atm;
      dry_aero_pre_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) const {
      const int i = team.league_rank();  // column index

      compute_dry_mixing_ratios(team, wet_atm_pre_, dry_atm_pre_, i);
      compute_dry_mixing_ratios(team, wet_atm_pre_, wet_aero_pre_,
                                dry_aero_pre_, i);
      team.team_barrier();

    }  // operator()

    // local variables for preprocess struct
    // number of horizontal columns and vertical levels
    int ncol_pre_, nlev_pre_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_pre_;
    mam_coupling::DryAtmosphere dry_atm_pre_;
    mam_coupling::AerosolState wet_aero_pre_, dry_aero_pre_;
  };  // MAMAci::Preprocess

 private:
  // preprocessing scratch pad
  Preprocess preprocess_;

  // IO structure to read in data for standard grids

  std::shared_ptr<AbstractRemapper> dmsSrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> dmsSrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState dmsSrfEmissTimeState_;
  srfEmissFunc::srfEmissInput dmsSrfEmissData_start_, dmsSrfEmissData_end_;
  srfEmissFunc::srfEmissOutput dmsSrfEmissData_out_;

  // Structures to store the data used for interpolation
  std::shared_ptr<AbstractRemapper> so2SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> so2SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState so2SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput so2SrfEmissData_start_, so2SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput so2SrfEmissData_out_;

  std::shared_ptr<AbstractRemapper> bc_a4SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> bc_a4SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState bc_a4SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput bc_a4SrfEmissData_start_, bc_a4SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput bc_a4SrfEmissData_out_;

  std::shared_ptr<AbstractRemapper> num_a1SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> num_a1SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState num_a1SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput num_a1SrfEmissData_start_,
      num_a1SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput num_a1SrfEmissData_out_;

  std::shared_ptr<AbstractRemapper> num_a2SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> num_a2SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState num_a2SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput num_a2SrfEmissData_start_,
      num_a2SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput num_a2SrfEmissData_out_;

  std::shared_ptr<AbstractRemapper> num_a4SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> num_a4SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState num_a4SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput num_a4SrfEmissData_start_,
      num_a4SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput num_a4SrfEmissData_out_;

  std::shared_ptr<AbstractRemapper> pom_a4SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> pom_a4SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState pom_a4SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput pom_a4SrfEmissData_start_,
      pom_a4SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput pom_a4SrfEmissData_out_;

  std::shared_ptr<AbstractRemapper> so4_a1SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> so4_a1SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState so4_a1SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput so4_a1SrfEmissData_start_,
      so4_a1SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput so4_a1SrfEmissData_out_;

  std::shared_ptr<AbstractRemapper> so4_a2SrfEmissHorizInterp_;
  std::shared_ptr<AtmosphereInput> so4_a2SrfEmissDataReader_;
  srfEmissFunc::srfEmissTimeState so4_a2SrfEmissTimeState_;
  srfEmissFunc::srfEmissInput so4_a2SrfEmissData_start_,
      so4_a2SrfEmissData_end_;
  srfEmissFunc::srfEmissOutput so4_a2SrfEmissData_out_;

};  // MAMSrfOnlineEmiss

}  // namespace scream

#endif  // EAMXX_MAM_SRF_ONLINE_EMISS_HPP
