#ifndef EAMXX_CONDITIONAL_SAMPLING_HPP
#define EAMXX_CONDITIONAL_SAMPLING_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic implements conditional sampling by applying a condition
 * to a field and outputting the input field values where the condition is met,
 * and fill values where the condition is not met.
 */

class ConditionalSampling : public AtmosphereDiagnostic {
public:
  // Constructors
  ConditionalSampling(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "ConditionalSampling"; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl();

protected:
  void initialize_impl(const RunType /*run_type*/);

  Real m_mask_val;

  std::string m_diag_name;    // X_where_Y_comp_VAL
  std::string m_input_f;      // X
  std::string m_condition_f;  // Y
  Real m_condition_v;         // VAL
  std::string m_condition_op; // comp
};

} // namespace scream

#endif // EAMXX_CONDITIONAL_SAMPLING_HPP
