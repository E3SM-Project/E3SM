
#ifndef EAMXX_VERT_CONTRACT_HPP
#define EAMXX_VERT_CONTRACT_HPP

#include "share/diagnostics/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will calculate the dp- or dz-weighted sum or average of a
 * field across the LEV tag dimension, producing an N-1 dimensional field.
 *
 * For "sum":   out = sum_lev(weight * f * mask)
 * For "avg":   out = sum_lev(weight * f * mask) / sum_lev(weight * mask)
 *
 * If the input field has a valid_mask, entries with mask==0 are excluded from
 * both numerator and denominator.  Output entries where the denominator is zero
 * are filled with fill_value and the output field's valid_mask is set to 0.
 */

class VertContractDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  VertContractDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "VertContractDiag"; }

  // Set the grid
  void create_requests ();

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();
  void initialize_impl();

  // Name of each field (because the diagnostic impl is generic)
  std::string m_diag_name;
  // Name of contraction method (avg, sum)
  std::string m_contract_method;
  // Name of weighting method (dp, dz, none)
  std::string m_weighting_method;

  // Weight field (dp/g, dz, or all-ones depending on weighting_method)
  Field m_weight;

  // For avg mode: denominator scratch field (same layout as m_diagnostic_output).
  Field m_weight_sum;

  // Ones field with the same layout as the weight field (for "avg" only)
  Field m_ones;
};

}  // namespace scream

#endif  // EAMXX_VERT_CONTRACT_HPP
