#ifndef EAMXX_ZONAL_AVERAGE_HPP
#define EAMXX_ZONAL_AVERAGE_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {
// TODO: Update this comment
/*
 * This diagnostic will calculate the area-weighted average of a field
 * across the COL tag dimension, producing an N-1 dimensional field
 * that is area-weighted average of the input field.
 */

class ZonalAvgDiag : public AtmosphereDiagnostic {

  // TODO: comment this, noting it's a utility function that could exist elsewhere
  //       also note the assumptions that field is arranged column first and
  //       result is arranged "bin" first
  static void compute_zonal_sum(const Field &result, const Field &field,
    const Field &weight, const Field &lat, const ekat::Comm *comm = nullptr);

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
  void compute_diagnostic_impl();

 protected:
  void initialize_impl(const RunType /*run_type*/);

  std::string m_diag_name;
  std::string m_bin_dim_name;
  int m_num_zonal_bins;

  Field m_lat;
  Field m_scaled_area;

};

}  // namespace scream

#endif  // EAMXX_ZONAL_AVERAGE_HPP
