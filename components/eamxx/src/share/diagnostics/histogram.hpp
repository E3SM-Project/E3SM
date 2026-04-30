#ifndef EAMXX_HISTOGRAM_HPP
#define EAMXX_HISTOGRAM_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {
/*
 * This diagnostic will calculate a histogram of a field across all dimensions
 * producing a one dimensional field, with CMP tag dimension named "bin", that
 * indicates how many times a field value in the specified range occurred.
 */

class Histogram : public AbstractDiagnostic {

public:
  // Constructors
  Histogram(const ekat::Comm &comm, const ekat::ParameterList &params,
            const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic
  std::string name() const { return "Histogram"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void initialize_impl();
  void compute_impl();

protected:
  std::string m_field_name;

  std::vector<Real> m_bin_reals;
  Field m_bin_values;
};

} // namespace scream

#endif // EAMXX_HISTOGRAM_HPP
