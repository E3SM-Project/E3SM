#ifndef EAMXX_MAM_ACI_HPP
#define EAMXX_MAM_ACI_HPP

//For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

//For declaring ACI class derived from atm process class
#include <share/atm_process/atmosphere_process.hpp>

//For physical constants
#include "physics/share/physics_constants.hpp"

namespace scream
{

class MAMAci final : public scream::AtmosphereProcess {

public:

  using KT = ekat::KokkosTypes<DefaultDevice>;

  mam4::NucleateIce nucleate_ice_;
  mam4::Hetfrz hetfrz_;

  // views for single- and multi-column data
  using view_1d       = typename KT::template view_1d<Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using view_3d       = typename KT::template view_3d<Real>;
  using const_view_3d = typename KT::template view_3d<const Real>;

  template <typename Scalar, typename MemoryTraits = Kokkos::MemoryManaged>
  using view_4d = KT::view<Scalar****,MemoryTraits>;

  // FIXME the time step for microphysics [s] need to get from the input
  const Real dtmicro_ = .0001;

private:
  // rho is air density [kg/m3]
  view_2d rho_;

  // w0_ is large scale velocity (m/s)
  view_2d w0_;

  // turbulent kinetic energy  [m^2/s^2]
  view_2d tke_;

  const_view_2d qv_dry_;
  const_view_2d cldfrac_;
  const_view_2d w_updraft_;

  view_2d aitken_dry_dia_;
  view_2d qc_coarse_bc_;
  view_2d qc_coarse_dst_;
  view_2d qc_coarse_mom_;
  view_2d qc_coarse_nacl_;
  view_2d qc_coarse_pom_;
  view_2d qc_coarse_soa_;

  view_2d qc_accum_bc_;
  view_2d qc_accum_dst_;
  view_2d qc_accum_mom_;
  view_2d qc_accum_nacl_;
  view_2d qc_accum_pom_;
  view_2d qc_accum_so4_;
  view_2d qc_accum_soa_;

  view_2d qi_accum_dst_;
  view_2d qi_accum_so4_;
  view_2d qi_accum_mom_;
  view_2d qi_accum_bc_;
  view_2d qi_accum_pom_;
  view_2d qi_accum_soa_;

  view_2d qi_coarse_dst_;
  view_2d qi_coarse_nacl_;
  view_2d qi_coarse_so4_;
  view_2d qi_coarse_mom_;
  view_2d qi_coarse_bc_;
  view_2d qi_coarse_pom_;
  view_2d qi_coarse_soa_;

  view_2d qi_pcarbon_bc_;
  view_2d qi_pcarbon_mom_;
  view_2d qi_pcarbon_pom_;
  view_2d nc_accum_;
  view_2d ni_accum_;
  view_2d ni_coarse_;
  view_2d ni_aitken_;
  const_view_3d dgnum_;
  view_2d nihf_;
  view_2d niim_;
  view_2d nidep_;
  view_2d nimey_;
  view_2d naai_hom_;
  view_2d naai_;
  const_view_2d liqcldf_;
  const_view_2d qc_;
  const_view_2d qi_;  
  const_view_2d nc_;  
  const_view_2d kvh_;  
  view_2d cldo_;

  view_2d zm_;
  view_3d state_q_;
  view_2d ncldwtr_;

  view_2d lcldn_;
  view_2d lcldo_;
  view_2d qcld_;
  view_2d tendnd_;
  view_2d ptend_q_[mam4::ndrop::nvar_ptend_q];
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
  view_2d raercol_cw_[mam4::ndrop::pver][2];
  view_2d raercol_[mam4::ndrop::pver][2];
  view_3d nact_;
  view_3d mact_;
  view_2d dropmixnuc_scratch_mem_[15];

  view_2d stratiform_cloud_fraction_;
  view_2d activation_fraction_accum_idx_;
  view_2d activation_fraction_coarse_idx_;

  view_2d hetfrz_immersion_nucleation_tend_;
  view_2d hetfrz_contact_nucleation_tend_;
  view_2d hetfrz_depostion_nucleation_tend_;
  view_2d diagnostic_scratch_[42];

  // Subgrid scale velocities
  view_2d wsub_, wsubice_, wsig_, w2_;
  // Top level for troposphere cloud physics
  // FIXME: This should be read in to make user selectable.
  const int top_lev_ = 6;

public:
  // Constructor
  MAMAci(const ekat::Comm& comm, const ekat::ParameterList& params);
  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;


  //Local variables
  
  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // number of aerosol modes
  int num_aero_modes_;
  

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct Preprocess {
    Preprocess() = default;


    //const_view_2d pdel_;    // hydrostatic "pressure thickness" at grid
                            // interfaces [Pa]
    
    // assigns local variables
    void set_variables(const const_view_2d&     pdel) {
      //p1del_ = pdel;
    } // set_variables
  }; // MAMAci::Preprocess


  // pre- and postprocessing scratch pads
  Preprocess preprocess_;

  // local atmospheric state column variables
  const_view_2d omega_; // Vertical pressure velocity [Pa/s] at midpoints
  const_view_2d p_mid_; // Total pressure [Pa] at midpoints
  const_view_2d p_int_; // Total pressure [Pa] at interfaces
  const_view_2d T_mid_; // Temperature[K] at midpoints
  const_view_2d pdel_;  // pressure thickess of layer [Pa]
  view_2d rpdel_; // Inverse of pdel_
  const_view_2d w_sec_; // Vertical velocity variance

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMAci

} // namespace scream


#endif // EAMXX_MAM_ACI_HPP
