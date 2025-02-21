#ifndef EAMXX_MAM_ACI_HPP
#define EAMXX_MAM_ACI_HPP

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

// For declaring ACI class derived from atm process class
#include <share/atm_process/atmosphere_process.hpp>

// For aerosol configuration
#include "mam4xx/aero_config.hpp"

// For calling ndrop functions
#include "mam4xx/ndrop.hpp"

namespace scream {

class MAMAci final : public scream::AtmosphereProcess {
 public:
  // declare some constant scratch space lengths
  static constexpr int hetro_scratch_   = 43;
  static constexpr int dropmix_scratch_ = 15;

  // views for multi-column data
  using view_2d       = scream::mam_coupling::view_2d;
  using const_view_2d = scream::mam_coupling::const_view_2d;
  using view_3d       = scream::mam_coupling::view_3d;
  using const_view_3d = scream::mam_coupling::const_view_3d;

 private:
  using KT = ekat::KokkosTypes<DefaultDevice>;

  mam4::NucleateIce nucleate_ice_;
  mam4::Hetfrz hetfrz_;

  // views for single-column data
  using view_1d       = scream::mam_coupling::view_1d;
  using const_view_1d = scream::mam_coupling::const_view_1d;

  //------------------------------------------------------------------------
  // ACI runtime ( or namelist) options
  //------------------------------------------------------------------------

  Real wsubmin_;                   // Minimum subgrid vertical velocity
  bool enable_aero_vertical_mix_;  // To enable vertical mixing of aerosols
  int top_lev_;                    // Top level for MAM4xx

  //------------------------------------------------------------------------
  // END: ACI runtime ( or namelist) options
  //------------------------------------------------------------------------

  // rho is air density [kg/m3]
  view_2d rho_;

  // w0_ is large scale velocity (m/s)
  view_2d w0_;

  // turbulent kinetic energy  [m^2/s^2]
  view_2d tke_;

  // Dry diameter of the aitken mode (for ice nucleation)
  view_2d aitken_dry_dia_;

  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;

  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;

  // aerosol dry diameter
  const_view_3d dgnum_;

  // ice nucleation diagnostic variables
  view_2d nihf_;
  view_2d niim_;
  view_2d nidep_;
  view_2d nimey_;
  view_2d naai_hom_;

  // ice nucleation output for FM
  view_2d naai_;

  // droplet activation inputs and outputs
  view_2d kvh_int_;  // Eddy diffusivity of heat at the interfaces
  const_view_2d liqcldf_;
  const_view_2d liqcldf_prev_;
  const_view_2d kvh_mid_;  // Eddy diffusivity of heat at the midpoints

  view_2d cloud_frac_;
  view_2d cloud_frac_prev_;
  view_2d qcld_;
  view_2d ptend_q_[mam4::aero_model::pcnst];
  view_3d factnum_;
  const_view_3d qqcw_input_;
  view_2d qqcw_[mam4::ndrop::ncnst_tot];
  view_2d ndropcol_;
  view_2d ndropmix_;
  view_2d nsource_;
  view_2d nc_inp_to_aci_;  // FIXME: TEMPORARY output
  view_2d ccn_0p02_;       // FIXME: TEMPORARY output
  view_2d ccn_0p05_;       // FIXME: TEMPORARY output
  view_2d ccn_0p1_;        // FIXME: TEMPORARY output
  view_2d ccn_0p2_;        // FIXME: TEMPORARY output
  view_2d ccn_0p5_;        // FIXME: TEMPORARY output
  view_2d ccn_1p0_;        // FIXME: TEMPORARY output
  view_2d wtke_;
  view_3d ccn_;
  view_2d coltend_[mam4::ndrop::ncnst_tot];
  view_2d coltend_cw_[mam4::ndrop::ncnst_tot];

  // raercol_cw_ and raercol_ are work arrays for dropmixnuc, allocated on the
  // stack.
  view_2d raercol_cw_[mam4::ndrop::pver][2];
  view_2d raercol_[mam4::ndrop::pver][2];

  view_3d nact_;
  view_3d mact_;
  view_2d dropmixnuc_scratch_mem_[dropmix_scratch_];

  // droplet activation output for the FM
  view_2d tendnd_;

  // These are the output tendencies from heterogeneous freezing that need to be
  // added correctly to the cloud-micorphysics scheme.
  view_2d hetfrz_immersion_nucleation_tend_;
  view_2d hetfrz_contact_nucleation_tend_;
  view_2d hetfrz_deposition_nucleation_tend_;

  view_2d diagnostic_scratch_[hetro_scratch_];

  // Subgrid scale velocities
  view_2d wsub_, wsubice_, wsig_, w2_;

  // local atmospheric state column variables
  const_view_2d pdel_;       // pressure thickess of layer [Pa]
  view_2d rpdel_;            // Inverse of pdel_
  const_view_2d w_sec_mid_;  // Vertical velocity variance at midpoints
  view_2d w_sec_int_;        // Vertical velocity variance at interfaces

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;

  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

  // A view array to carry cloud borne aerosol mmrs/nmrs
  view_2d qqcw_fld_work_[mam4::ndrop::ncnst_tot];

  // A view to carry interstitial aerosol mmrs/nmrs
  view_3d state_q_work_;

 public:
  // Constructor
  MAMAci(const ekat::Comm &comm, const ekat::ParameterList &params);

  // Process metadata: Return type of the process
  AtmosphereProcessType type() const override {
    return AtmosphereProcessType::Physics;
  }

  // Return name of the process
  std::string name() const override { return "mam4_aci"; }

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override {
    return mam_coupling::buffer_size(ncol_, nlev_);
  }

  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override{/*DO NOTHING*/};

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
      set_min_background_mmr(team, dry_aero_pre_,
                             i);  // dry_atm_pre_ is the output
    }                             // operator()

    // local variables for preprocess struct
    // number of horizontal columns and vertical levels
    int ncol_pre_, nlev_pre_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_pre_;
    mam_coupling::DryAtmosphere dry_atm_pre_;
    mam_coupling::AerosolState wet_aero_pre_, dry_aero_pre_;
  };  // MAMAci::Preprocess

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
    }  // operator()

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

};  // MAMAci

}  // namespace scream

#endif  // EAMXX_MAM_ACI_HPP
