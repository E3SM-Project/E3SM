#ifndef EAMXX_HISTOGRAM_HPP
#define EAMXX_HISTOGRAM_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {
/*
 * This diagnostic will calculate a histogram of a field across all dimensions
 * producing a one dimensional field, with CMP tag dimension named "bin", that
 * indicates how many times a field value in the specified range occurred.
 *
 * Examples:
 *  - If the input field is T_mid(ncol,nlev) and the bin configuration string is 100_200_300_400_500,
 *    the diag will have layout (cmp), with extent 5 (5 bins). diag(1) will count
 *    the number of (icol,ilev) entries where 200<T_mid<300.
 * Notes:
 *  - we do add a bin (-inf,100) and (500,inf)  before/after the provided ones, to catch the tails.
 *  - the bins endpoints MUST be listed in strictly increasing order (we error out if they are not).
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
