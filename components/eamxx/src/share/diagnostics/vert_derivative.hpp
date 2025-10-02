
#ifndef EAMXX_VERT_DERIVATIVE_HPP
#define EAMXX_VERT_DERIVATIVE_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will calculate the derivative of a field X in
 * the vertical direction such that dX/dp or dX/dz
 */

class VertDerivativeDiag : public AtmosphereDiagnostic {
public:
  // Constructors
  VertDerivativeDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "VertDerivativeDiag"; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl();
  void initialize_impl(const RunType /*run_type*/);

  // Name of each field (because the diagnostic impl is generic)
  std::string m_diag_name;
  // Name of derivative method (differential dp or dz)
  std::string m_derivative_method;
};

} // namespace scream

#endif // EAMXX_VERT_DERIVATIVE_HPP
