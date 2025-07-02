#ifndef EAMXX_ZONAL_AVERAGE_HPP
#define EAMXX_ZONAL_AVERAGE_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {
/*
 * This diagnostic will calculate area-weighted zonal averages of a field across
 * the COL tag dimension producing an N dimensional field, where the COL tag
 * dimension is replaced by a CMP tag dimension named "bin" that indicates which
 * zonal band the average value corresponds to.
 */

class ZonalAvgDiag : public AtmosphereDiagnostic {

public:
  // Constructors
  ZonalAvgDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

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
  int m_num_zonal_bins;

  Field m_lat;
  Field m_scaled_area;
  Field m_bin_to_cols;

};

} // namespace scream

#endif // EAMXX_ZONAL_AVERAGE_HPP
