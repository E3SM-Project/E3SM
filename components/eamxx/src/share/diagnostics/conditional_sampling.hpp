#ifndef EAMXX_CONDITIONAL_SAMPLING_HPP
#define EAMXX_CONDITIONAL_SAMPLING_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic implements conditional sampling by applying a condition
 * to a field and outputting the input field values where the condition is met,
 * and fill values where the condition is not met.
 *
 * The general diag name is X_where_Y_cmp_Z, where:
 *  - cmp is a string among "eq", "ne", "lt", "le", "gt", "ge"
 *  - X, Y, and Z are one of the following combinations:
 *    - X,Y,Z are fields with the same layout (e.g., T_mid_where_qc_gt_qr)
 *    - X,Y are fields with same layout and Z is a number (e.g., T_mid_where_qc_gt_1e-6)
 *    - X is a field, and Y is the string "lev" and Z is an integer (e.g., T_mid_where_lev_gt_10)
 *    - any of the above, with X="mask": compute the mask field of where the condition holds
 *      E.g., mask_where_T_mid_gt_273 gives diag=1 where T_mid>273 and 0 elsewhere
 */

class ConditionalSampling : public AbstractDiagnostic {
public:
  // Constructors
  ConditionalSampling(const ekat::Comm &comm, const ekat::ParameterList &params,
                      const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "ConditionalSampling"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl();

  void initialize_impl();

protected:

  std::string m_diag_name;      // X_where_Y_comp_Z
  std::string m_input_f;        // X
  std::string m_condition_lhs;  // Y
  std::string m_condition_rhs;  // Z

  bool            m_lhs_is_lev;     // true if we sample w.r.t. level
  bool            m_rhs_is_field;   // true if Z is a field (not a value)
  bool            m_diag_is_mask;   // true if X="mask"

  Comparison      m_condition_cmp;
  ScalarWrapper   m_rhs_value;    // Only if m_rhs_is_field=false

  Field m_lev_mask; // Only used if m_lhs_is_lev=true, both for mask diagnostics
                    // and for broadcasting in non-mask diagnostics
};

} // namespace scream

#endif // EAMXX_CONDITIONAL_SAMPLING_HPP
