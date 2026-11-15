#ifndef EAMXX_HORIZ_AVERAGE_HPP
#define EAMXX_HORIZ_AVERAGE_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will calculate the weighted average of a field across the
 * COL tag dimension, producing an N-1 dimensional field:
 *   sum(w*f*mask) / sum(w*mask)
 * where mask is the input's valid_mask, if any. Output entries where the
 * denominator is 0 are set to fill_value, and the output field's valid_mask
 * is set to 0 at those locations.
 *
 * The weight w is
 *  - the grid's cell area (weighting_method='area', the default);
 *  - 1 (weighting_method='none');
 *  - or, with 'weight_name', the field of that name, as in xarray's
 *    X.weighted(W).mean('col'). The name 'area' is the grid's cell area, which
 *    is grid geometry, not a field. A weight field must have layout (COL) and
 *    no mask; the output then has a mask, 0 where the weights sum to 0.
 */

class HorizAvg : public AbstractDiagnostic {
 public:
  // Constructors
  HorizAvg(const ekat::Comm &comm, const ekat::ParameterList &params,
               const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "HorizAvg"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_impl();

 protected:
  void initialize_impl();

  // The weight: the grid's area, ones, or an input field (see above)
  Field m_area;

  // The name of the input field used as weight, if any
  std::string m_weight_name;
  bool m_weight_is_field = false;

  // Utility fields to compute the (scalar) denominator
  Field m_denom;
  Field m_ones;
};

}  // namespace scream

#endif  // EAMXX_HORIZ_AVERAGE_HPP
