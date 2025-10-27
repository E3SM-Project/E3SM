#ifndef EAMXX_UNARY_OPS_HPP
#define EAMXX_UNARY_OPS_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will compute the visible aerosol optical depth.
 */

class UnaryOpsDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  UnaryOpsDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "UnaryOpsDiag"; }

  // Set the grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  void initialize_impl(const RunType /*run_type*/) override;

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl() override;

  std::string m_fn;
  std::string m_op;
};

}  // namespace scream

#endif  // EAMXX_UNARY_OPS_HPP
