#ifndef EAMXX_HORIZ_AVERAGE_HPP
#define EAMXX_HORIZ_AVERAGE_HPP

#include "share/diagnostics/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will calculate the area-weighted average of a field
 * across the COL tag dimension, producing an N-1 dimensional field
 * that is area-weighted average of the input field.
 *
 * If the input field has a valid_mask, the average is computed only over
 * valid (mask != 0) entries:  sum(area*f*mask) / sum(area*mask).
 * Output entries where sum(area*mask)==0 are set to fill_value,
 * and the output field's valid_mask is set to 0 at those locations.
 */

class HorizAvgDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  HorizAvgDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "HorizAvg"; }

  // Set the grid
  void create_requests ();

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

 protected:
  void initialize_impl();

  // Area field
  Field m_area;

  // Utility fields to compute the (scalar) denominator
  Field m_denom;
  Field m_ones;
};

}  // namespace scream

#endif  // EAMXX_HORIZ_AVERAGE_HPP
