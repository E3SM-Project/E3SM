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


  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
}; // MAMAci

} // namespace scream


#endif // EAMXX_MAM_ACI_HPP
