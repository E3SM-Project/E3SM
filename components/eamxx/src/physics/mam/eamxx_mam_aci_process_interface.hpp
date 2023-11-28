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

  // views for single- and multi-column data
  using const_view_2d = typename KT::template view_2d<const Real>;
  

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


  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMAci

} // namespace scream


#endif // EAMXX_MAM_ACI_HPP
