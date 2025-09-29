#ifndef EAMXX_BELOW_OR_ABOVE_INTERFACE
#define EAMXX_BELOW_OR_ABOVE_INTERFACE

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will compute the visible aerosol optical depth.
 */

class BelowOrAboveInterface : public AtmosphereDiagnostic {
 public:
  // Constructors
  BelowOrAboveInterface(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "BelowOrAboveInterface"; }

  // Set the grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  void initialize_impl(const RunType /*run_type*/) override;

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl() override;

  int m_num_cols;
  int m_num_levs;

  std::string m_name;
  std::string m_type;
};

}  // namespace scream

#endif  // EAMXX_BELOW_OR_ABOVE_INTERFACE
