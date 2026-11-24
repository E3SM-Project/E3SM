#ifndef EAMXX_VERT_CONTRACT_HPP
#define EAMXX_VERT_CONTRACT_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

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
 * are invalid: the output field's valid_mask is set to 0 at those locations.
 *
 * The weight is dp/g ("dp"), dz ("dz"), 1 ("none"), or, with 'weight_name',
 * any input field, as in X.sum('lev',weights=W). With 'weight_name' the
 * semantics are those of xarray's weighted reductions:
 *  - sum: invalid entries of f are skipped, so a column with no valid level
 *         gives 0 (not invalid), and the output has no mask;
 *  - avg: the denominator sums W over the valid entries of f only, and the
 *         output is invalid wherever that is zero, masked input or not.
 * W itself must not be masked (xarray rejects missing weights), and must be
 * (lev) or f's layout without CMP.
 */

class VertContract : public AbstractDiagnostic {
 public:
  // Constructors
  VertContract(const ekat::Comm &comm, const ekat::ParameterList &params,
               const std::shared_ptr<const AbstractGrid> &grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "VertContract"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_impl();
  void initialize_impl();

  std::string m_field_name; // Input field name

  // Name of contraction method (avg, sum)
  std::string m_contract_method;
  // Name of weighting method (dp, dz, none)
  std::string m_weighting_method;

  // Name of the input field used as weight, if any (see class description)
  std::string m_weight_name;

  // Weight field (dp/g, dz, or all-ones depending on weighting_method)
  Field m_weight;

  // For avg mode: denominator scratch field (same layout as get()).
  Field m_weight_sum;

  // Ones field with the same layout as the weight field (for "avg" only)
  Field m_ones;
};

}  // namespace scream

#endif  // EAMXX_VERT_CONTRACT_HPP
