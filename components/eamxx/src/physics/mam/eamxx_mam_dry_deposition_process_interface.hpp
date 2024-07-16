#ifndef EAMXX_MAM_DRYDEP_HPP
#define EAMXX_MAM_DRYDEP_HPP

// For declaring dry deposition class derived from atm process class
#include <share/atm_process/atmosphere_process.hpp>

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

// For component name
#include <string>

namespace scream {

// The process responsible for handling MAM4 dry deposition. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMDryDep final : public scream::AtmosphereProcess {
 public:
  using view_1d       = Field::view_dev_t<Real *>;
  using view_2d       = Field::view_dev_t<Real **>;
  using view_3d       = Field::view_dev_t<Real ***>;
  using const_view_1d = Field::view_dev_t<const Real *>;
  using const_view_3d = Field::view_dev_t<const Real ***>;

  static constexpr int num_aero_modes = mam_coupling::num_aero_modes();
  static constexpr int aerosol_categories_ =
      mam4::DryDeposition::aerosol_categories;
  static constexpr int n_land_type = mam4::DryDeposition::n_land_type;

 private:
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

  view_3d qtracers_;
  view_3d ptend_q_;

  // inputs
  // FIXME: collect all inputs and outputs together
  view_1d fraction_landuse_[mam4::DryDeposition::n_land_type];
  view_1d vlc_trb_[mam4::AeroConfig::num_modes()][aerosol_categories_];
  view_2d rho_;

  view_2d vlc_dry_[mam4::AeroConfig::num_modes()][aerosol_categories_];
  view_2d vlc_grv_[mam4::AeroConfig::num_modes()][aerosol_categories_];
  view_2d dqdt_tmp_[mam4::aero_model::pcnst];

  view_2d qqcw_[mam4::aero_model::pcnst];

 public:
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // Constructor
  MAMDryDep(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The type of subcomponent
  AtmosphereProcessType type() const override {
    return AtmosphereProcessType::Physics;
  }

  // The name of the subcomponent
  std::string name() const override { return "mam_dry_deposition"; }

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
  void finalize_impl() override{/*Do nothing*/};

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
      // vertical heights has to be computed after computing dry mixing ratios
      // for atmosphere
      compute_vertical_layer_heights(team, dry_atm_pre_, i);
      compute_updraft_velocities(team, wet_atm_pre_, dry_atm_pre_, i);
    }  // Preprocess operator()

    // local variables for preprocess struct
    // number of horizontal columns and vertical levels
    int ncol_pre_, nlev_pre_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_pre_;
    mam_coupling::DryAtmosphere dry_atm_pre_;
    mam_coupling::AerosolState wet_aero_pre_, dry_aero_pre_;
  };  // Preprocess

  // Atmosphere processes often have a post-processing step prepares output
  // from this process for the Field Manager. This functor implements this
  // step, which is called during run_impl.
  // Postprocessing functor
  struct Postprocess {
    Postprocess() = default;

    // on host: initializes postprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere &wet_atm,
                    const mam_coupling::AerosolState &wet_aero,
                    const mam_coupling::DryAtmosphere &dry_atm,
                    const mam_coupling::AerosolState &dry_aero) {
      ncol_post_     = ncol;
      nlev_post_     = nlev;
      wet_atm_post_  = wet_atm;
      wet_aero_post_ = wet_aero;
      dry_atm_post_  = dry_atm;
      dry_aero_post_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) const {
      const int i = team.league_rank();  // column index
      compute_wet_mixing_ratios(team, dry_atm_post_, dry_aero_post_,
                                wet_aero_post_, i);
    }  // operator() Postprocess

    // number of horizontal columns and vertical levels
    int ncol_post_, nlev_post_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_post_;
    mam_coupling::DryAtmosphere dry_atm_post_;
    mam_coupling::AerosolState wet_aero_post_, dry_aero_post_;
  };  // Postprocess

 private:
  // pre- and postprocessing scratch pads
  Preprocess preprocess_;
  Postprocess postprocess_;
};  // MAMDryDep

}  // namespace scream

#endif  // EAMXX_MAM_DRYDEP_HPP
