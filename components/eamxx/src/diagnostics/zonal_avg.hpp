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

  template <bool USE_WEIGHT>
  KOKKOS_INLINE_FUNCTION Real evaluate_weight(
    const Field::get_view_type<const scream::Real *, scream::Device> &weight,
    const int &i);

  template <bool USE_WEIGHT=true>
  void compute_zonal_sum(const Field &field, const Field &weight,
    const Field &lat, const Field &result, const ekat::Comm &comm);

  void compute_zonal_area(const Field &field, const Field &lat,
    const Field &result, const ekat::Comm &comm);

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

  Field m_area;
  Field m_lat;
  int m_lat_num;

  Field m_zonal_area;

};

}  // namespace scream

#endif  // EAMXX_ZONAL_AVERAGE_HPP
