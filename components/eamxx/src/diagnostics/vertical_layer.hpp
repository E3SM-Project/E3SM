#ifndef EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP
#define EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce data related to the vertical layer based on
 * the parameter "diag_name" (required). The following can be produced:
 *   - diag_name = "dz":                       Vertical layer thickness for each column and level
 *   - diag_name = "z_int":                    Vertical layer height at each column and interface level.
 *                                             Values are computed from sea level (i.e., surf_geopotential=0).
 *   - diag_name = "geopotential_int":         Same as z_int, but computed from topography data
 *                                             (i.e., surf_geopotential=phis).
 *   - diag_name = "z_mid"/"geopotential_mid": Same as z_int/geopotential_int but at midpoint levels.
 */

class VerticalLayerDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  VerticalLayerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "VerticalLayer"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
  void compute_diagnostic_impl ();
  void initialize_impl (const RunType /* run_type */);
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int PackSize>
  void do_compute_diagnostic_impl ();
protected:

  // Keep track of field dimensions
  Int m_num_cols;
  Int m_num_levs;

  // Temporaries to use for calculation of dz, z_int, and z_mid
  Field m_tmp_interface;
  Field m_tmp_midpoint;

  // The diagnostic name.
  std::string m_diag_name;

  // Store if the diagnostic output field exists on interface values
  bool m_is_interface_layout;

  // True z_mid/int and geopotential_mid/int, false for height_mid/int. Unused for dz
  bool m_from_sea_level;

  // If true, output is a geopotential (units m2/s2), otherwise an elevation
  bool m_geopotential;
};

} //namespace scream

#endif // EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP
