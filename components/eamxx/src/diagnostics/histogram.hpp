#ifndef EAMXX_HISTOTGRAM_HPP
#define EAMXX_HISTOGRAM_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {
/*
 * This diagnostic will calculate a histogram of a field across all dimensions
 * producing a one dimensional field, with CMP tag dimension named "bin", that
 * indicates how many times a field value in the specified range occurred.
 */

class HistogramDiag : public AtmosphereDiagnostic {

public:
  // Constructors
  HistogramDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic
  std::string name() const { return m_diag_name; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void initialize_impl(const RunType /*run_type*/);
  void compute_diagnostic_impl();

protected:
  std::string m_diag_name;
  std::vector<Real> m_bin_reals;
  Field m_bin_values;
};

} // namespace scream

#endif // EAMXX_HISTOGRAM_HPP
