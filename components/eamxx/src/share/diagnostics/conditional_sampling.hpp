#ifndef EAMXX_CONDITIONAL_SAMPLING_HPP
#define EAMXX_CONDITIONAL_SAMPLING_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic implements conditional sampling by applying a condition
 * to a field and outputting the input field values where the condition is met,
 * and fill values where the condition is not met.
 *
 * It is requested as X.where(Y cmp Z), where:
 *  - cmp is one of ==, !=, <, <=, >, >=
 *  - X, Y, and Z are one of the following combinations:
 *    - X,Y,Z are fields with the same layout (e.g., T_mid.where(qc>qr))
 *    - X,Y are fields with same layout and Z is a number (e.g., T_mid.where(qc>1e-6))
 *    - X is a field, and Y is the placeholder "lev" and Z is an integer (e.g., T_mid.where(lev>10))
 *    - any of the above, with X="mask": compute the mask field of where the condition holds
 *      E.g., mask.where(T_mid>273) gives diag=1 where T_mid>273 and 0 elsewhere
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
  Field m_lev_mask_bcast; // m_lev_mask broadcast to the output layout (non-mask diagnostics only)
};

} // namespace scream

#endif // EAMXX_CONDITIONAL_SAMPLING_HPP
