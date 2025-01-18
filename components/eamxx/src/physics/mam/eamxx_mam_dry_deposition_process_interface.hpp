#ifndef EAMXX_MAM_DRYDEP_HPP
#define EAMXX_MAM_DRYDEP_HPP

// For declaring dry deposition class derived from atm process class
#include <share/atm_process/atmosphere_process.hpp>

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

// For AtmosphereInput
#include "share/io/scorpio_input.hpp"

// For component name
#include <string>

namespace scream {

// The process responsible for handling MAM4 dry deposition. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMDryDep final : public scream::AtmosphereProcess {
 public:
  static constexpr int num_aero_modes = mam_coupling::num_aero_modes();
  static constexpr int aerosol_categories_ =
      mam4::DryDeposition::aerosol_categories;
  static constexpr int n_land_type = mam4::DryDeposition::n_land_type;

  using view_1d       = Field::view_dev_t<Real *>;
  using view_2d       = Field::view_dev_t<Real **>;
  using view_3d       = Field::view_dev_t<Real ***>;
  using view_4d       = Field::view_dev_t<Real ****>;
  using const_view_1d = Field::view_dev_t<const Real *>;
  using const_view_2d = Field::view_dev_t<const Real **>;
  using const_view_3d = Field::view_dev_t<const Real ***>;

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

  /* Note on mam4::DryDeposition::aerosol_categories = 4
     used in deposition velocity dimension defined below. These
     correspond to the two attachment states and two moments:
     0 - interstitial aerosol, 0th moment (i.e., number)
     1 - interstitial aerosol, 3rd moment (i.e., volume/mass)
     2 - cloud-borne aerosol,  0th moment (i.e., number)
     3 - cloud-borne aerosol,  3rd moment (i.e., volume/mass)
     see comments in the DryDeposition class in mam4xx.
  */
  // Output deposition velocity of turbulent dry deposition [m/s]
  // Dimensions
  //    [numer of modes, aerosol_categories_, num columns]
  view_3d vlc_trb_;

  // Output deposition velocity of gravitational settling [m/s]
  // Dimensions
  //   [num_modes, aerosol_categories_, num columns, num levels]
  view_4d vlc_grv_;

  // Output deposition velocity, [m/s]
  // fraction landuse weighted sum of vlc_grv and vlc_trb
  // Dimensions
  //   [num_modes, aerosol_categories_, num columns, num levels]
  view_4d vlc_dry_;

  // Output of the the mixing ratio tendencies [kg/kg/s or 1/kg/s]
  // Dimensions
  //   [num columns, num levels, mam4::aero_model::pcnst]
  // Packed the same way qtracers_ is layed out.
  view_3d ptend_q_;

  // Work array to hold the mixing ratios [kg/kg or 1/kg]
  // Dimensions
  //   [num columns, num levels, mam4::aero_model::pcnst]
  // Packs AerosolState::int_aero_nmr
  // and   AerosolState::int_aero_nmr
  // into one array, hence is mixed kg/kg and 1/kg.
  view_3d qtracers_;

  // Work array to hold the air density [kg/m3]
  // Dimensions
  //   [num columns, num levels]
  // Calculated from air pressure at layer midpoint,
  // Constants::r_gas_dry_air and air temperture.
  view_2d rho_;

  // Work array to hold tendency for 1 species [kg/kg/s] or [1/kg/s]
  // Dimensions
  //   [mam4::aero_model::pcnst, num column, num level]
  view_3d dqdt_tmp_;

  // Work array to hold cloud borne aerosols mixing ratios [kg/kg or 1/kg]
  // Dimensions
  //   [mam4::aero_model::pcnst, num column, num level]
  // Filled with Prognostics::n_mode_c and Prognostics::q_aero_c
  view_3d qqcw_;

  // For reading fractional land use file
  std::shared_ptr<AbstractRemapper> horizInterp_;
  std::shared_ptr<AtmosphereInput> dataReader_;
  const_view_2d frac_landuse_;

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
