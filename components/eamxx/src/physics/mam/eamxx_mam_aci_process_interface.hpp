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


  using KT = ekat::KokkosTypes<DefaultDevice>;

  mam4::NucleateIce nucleate_ice_;

  // views for single- and multi-column data
  using view_1d       = typename KT::template view_1d<Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using view_3d       = typename KT::template view_3d<Real>;
  using const_view_3d = typename KT::template view_3d<const Real>;

  // rho is air density [kg/m3]
  view_2d rho_;

  // w0_ is large scale velocity (m/s)
  view_2d w0_;

  // turbulent kinetic energy  [m^2/s^2]
  view_2d tke_;

  const_view_2d qv_dry_;
  const_view_2d cldfrac_;
  //const_view_2d w_updraft_;
  view_2d w_updraft_;

  view_2d aitken_dry_dia_;
  view_2d q_coarse_dst_;
  view_2d q_coarse_nacl_;
  view_2d q_coarse_so4_;
  view_2d q_coarse_mom_;
  view_2d q_coarse_bc_;
  view_2d q_coarse_pom_;
  view_2d q_coarse_soa_;
  view_2d n_coarse_;
  view_2d n_aitken_;
  view_3d dgnum_;
  view_2d nihf_;
  view_2d niim_;
  view_2d nidep_;
  view_2d nimey_;
  view_2d naai_hom_;
  view_2d naai_;
  const_view_2d liqcldf_;
  const_view_2d qc_;
  const_view_2d qi_;
  view_2d lcldn_;
  view_2d lcldo_;

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
  const_view_2d pdel_;    // hydrostatic "pressure thickness" at grid
                          // interfaces [Pa]
  const_view_2d omega_; // Vertical pressure velocity [Pa/s] at midpoints
  const_view_2d p_mid_; // Total pressure [Pa] at midpoints
  const_view_2d T_mid_; // Temperature[K] at midpoints
  const_view_2d w_sec_; // Vertical velocity variance

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMAci

} // namespace scream


#endif // EAMXX_MAM_ACI_HPP
