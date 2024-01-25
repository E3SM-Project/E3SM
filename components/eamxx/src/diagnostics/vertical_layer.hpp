#ifndef EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP
#define EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

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
  using Pack          = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using PF            = scream::PhysicsFunctions<DefaultDevice>;
  using KT            = KokkosTypes<DefaultDevice>;
  using MemberType    = typename KT::MemberType;
  using view_1d_const = typename KT::template view_1d<const Real>;
  using view_2d       = typename KT::template view_2d<Pack>;

  // Constructors
  VerticalLayerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic.
  std::string name () const { return m_diag_name; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  // Keep track of field dimensions
  Int m_num_cols;
  Int m_num_levs;

  // Temporary view to set dz, z_mid, and z_int
  view_2d m_tmp_interface_view;
  view_2d m_tmp_midpoint_view;

  // The diagnostic name. This will dictate which
  // field in the computation is output (dz, z_int, or z_mid).
  std::string m_diag_name;

  // Store if we only need to compute dz to save computation/memory requirements.
  bool m_only_compute_dz;

  // Store if the diagnostic output field exists on interface values
  bool m_is_interface_layout;

  // If z_int or z_mid is computed, determine whether the BC
  // is from sea level or not (from topography data).
  bool m_from_sea_level;

}; // class VerticalLayerDiagnostic

} //namespace scream

#endif // EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP
