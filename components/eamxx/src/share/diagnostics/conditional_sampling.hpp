#ifndef EAMXX_CONDITIONAL_SAMPLING_HPP
#define EAMXX_CONDITIONAL_SAMPLING_HPP

#include "share/diagnostics/atmosphere_diagnostic.hpp"

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
  void create_requests ();

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl();

  void initialize_impl();

protected:

  // General syntax is X_where_Y_comp_Z
  // where either:
  //  - X,Y,Z are fields with the same layout
  //  - X,Y are fields with same layout and Z is a number
  //  - X is a field, and Y is the string "lev" and Z is an integer
  //  - all of the above, but with X="mask", which computes the mask field of where the condition holds
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
