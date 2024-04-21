#ifndef EAMXX_MAM_ACI_HPP
#define EAMXX_MAM_ACI_HPP

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

// For declaring ACI class derived from atm process class
#include <share/atm_process/atmosphere_process.hpp>

// For EKAT units package
#include "ekat/util/ekat_units.hpp"

// For aerosol configuration
#include "mam4xx/aero_config.hpp"

// For calling ndrop functions
#include "mam4xx/ndrop.hpp"

namespace scream {

class MAMAci final : public scream::AtmosphereProcess {
  using KT = ekat::KokkosTypes<DefaultDevice>;

  mam4::NucleateIce nucleate_ice_;
  mam4::Hetfrz hetfrz_;

  // views for single- and multi-column data
  using view_1d       = scream::mam_coupling::view_1d;
  using view_2d       = scream::mam_coupling::view_2d;
  using view_3d       = scream::mam_coupling::view_3d;
  using const_view_1d = scream::mam_coupling::const_view_1d;
  using const_view_2d = scream::mam_coupling::const_view_2d;
  using const_view_3d = scream::mam_coupling::const_view_3d;

  // FIXME: Should the following variables be public? They are like that in
  // micriphysics and optics codes

  // ACI runtime ( or namelist) options
  Real wsubmin_;

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

  // FIXME: Make sure these output are used somewhere and add comments about
  // these
  view_2d nihf_;
  view_2d niim_;
  view_2d nidep_;
  view_2d nimey_;
  view_2d naai_hom_;
  view_2d naai_;
  view_2d kvh_int_;  // Eddy diffusivity of heat at the interfaces
  const_view_2d liqcldf_;
  const_view_2d liqcldf_prev_;
  const_view_2d kvh_mid_;  // Eddy diffusivity of heat at the midpoints

  view_2d cloud_frac_;
  view_2d cloud_frac_prev_;
  view_2d qcld_;
  view_2d tendnd_;
  // ptend_q_ is just ptend_q_output_ reformatted.
  view_2d ptend_q_[mam4::aero_model::pcnst];
  view_3d ptend_q_output_;
  view_3d factnum_;
  const_view_3d qqcw_input_;
  view_2d qqcw_[mam4::ndrop::ncnst_tot];
  view_2d ndropcol_;
  view_2d ndropmix_;
  view_2d nsource_;
  view_2d wtke_;
  view_3d ccn_;
  view_3d coltend_outp_;
  view_2d coltend_[mam4::ndrop::ncnst_tot];
  view_3d coltend_cw_outp_;
  view_2d coltend_cw_[mam4::ndrop::ncnst_tot];

  // raercol_cw_ and raercol_ are work arrays for dropmixnuc, allocated on the
  // stack.
  view_2d raercol_cw_[mam4::ndrop::pver][2];
  view_2d raercol_[mam4::ndrop::pver][2];

  view_3d nact_;
  view_3d mact_;

  static constexpr int dropmix_scratch_ = 15;
  view_2d dropmixnuc_scratch_mem_[dropmix_scratch_];

  view_2d stratiform_cloud_fraction_;
  view_2d activation_fraction_accum_idx_;
  view_2d activation_fraction_coarse_idx_;

  // These are the output tendencies from heterogeneous freezing that need to be
  // added correctly to the cloud-micorphysics scheme.
  view_2d hetfrz_immersion_nucleation_tend_;
  view_2d hetfrz_contact_nucleation_tend_;
  view_2d hetfrz_depostion_nucleation_tend_;
  view_2d diagnostic_scratch_[42];

  // Subgrid scale velocities
  view_2d wsub_, wsubice_, wsig_, w2_;
  // Top level for troposphere cloud physics
  // FIXME: This should be read in to make user selectable.
  const int top_lev_ = 6;

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
  // FIXME: 25 should be a const int
  view_2d qqcw_fld_work_[25];

  // A view to carry interstitial aerosol mmrs/nmrs
  view_3d state_q_work_;

 public:
  // Constructor
  MAMAci(const ekat::Comm &comm, const ekat::ParameterList &params);

  // Process metadata: Return type of the process
  AtmosphereProcessType MAMAci::type() const {
    return AtmosphereProcessType::Physics;
  }

  // Return name of the process
  std::string MAMAci::name() const { return "mam4_aci"; }

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t MAMAci::requested_buffer_size_in_bytes() const {
    return mam_coupling::buffer_size(ncol_, nlev_);
  }

  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

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
  // pre- and postprocessing scratch pads
  Preprocess preprocess_;

};  // MAMAci

}  // namespace scream

#endif  // EAMXX_MAM_ACI_HPP
