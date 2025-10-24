#ifndef EAMXX_BINARY_OPS_DIAG_HPP
#define EAMXX_BINARY_OPS_DIAG_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will perform binary ops
 * like +, -, *, ÷ on two input fields.
 */

class BinaryOpsDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  BinaryOpsDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "BinaryOpsDiag"; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

  void initialize_impl(const RunType /*run_type*/) override;

  std::string m_field_1;
  std::string m_field_2;
  std::string m_binary_op;

};  // class BinaryOpsDiag

}  // namespace scream

#endif  // EAMXX_BINARY_OPS_DIAG_HPP
