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

  // TODO: make it a local function in the cpp file
  // Utility to compute the contraction of a field along its column dimension.
  // This is equivalent to f_out = einsum('i,i...k->...k', weight, f_in).
  // The implementation is such that:
  // - all Field objects must be allocated
  // - the first dimension for field, weight, and lat is for the columns (COL)
  // - the first dimension for result is for the zonal bins (CMP,"bin")
  // - field and result must be the same dimension, up to 3
  // TODO: make it a local function in the cpp file
  static void compute_zonal_sum(const Field &result, const Field &field, const Field &weight,
                                const Field &lat, const ekat::Comm *comm = nullptr);

protected:
  std::string m_diag_name;
  int m_num_zonal_bins;

  Field m_lat;
  Field m_scaled_area;
};

} // namespace scream

#endif // EAMXX_ZONAL_AVERAGE_HPP
