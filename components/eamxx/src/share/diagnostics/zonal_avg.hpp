#ifndef EAMXX_ZONAL_AVERAGE_HPP
#define EAMXX_ZONAL_AVERAGE_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {
/*
 * This diagnostic will calculate area-weighted zonal averages of a field across
 * the COL tag dimension producing an N dimensional field, where the COL tag
 * dimension is replaced by a CMP tag dimension named "bin" that indicates which
 * zonal band the average value corresponds to. All bins are "closed" at the
 * lower value and "open" at the upper value (lat_lower <= lat < lat_upper),
 * with the exception of the "last" bin that is closed at both ends to capture
 * any column that is centered at the northern pole (lat_lower <= lat <= 90).
 */

class ZonalAvg : public AbstractDiagnostic {

public:
  // Constructors
  ZonalAvg(const ekat::Comm &comm, const ekat::ParameterList &params,
           const std::shared_ptr<const AbstractGrid> &grid);

  // The name of the diagnostic
  std::string name() const { return "ZonalAvg"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void initialize_impl();
  void compute_diagnostic_impl();

protected:

  std::string m_field_name;
  int m_num_zonal_bins;

  Field m_lat;
  Field m_bin_to_cols;

  // Weight field if input is NOT masked
  Field m_scaled_area;

  // Weight field if input IS masked
  Field m_area;

  // Fields used to compute masked area zonal sum at runtime if field is masked
  Field m_zonal_area;
  Field m_ones;
};

} // namespace scream

#endif // EAMXX_ZONAL_AVERAGE_HPP
