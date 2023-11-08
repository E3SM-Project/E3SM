#ifndef EAMXX_MAM_ACI_HPP
#define EAMXX_MAM_ACI_HPP

#include <share/atm_process/atmosphere_process.hpp>

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
}; // MAMAci

} // namespace scream


#endif // EAMXX_MAM_ACI_HPP
