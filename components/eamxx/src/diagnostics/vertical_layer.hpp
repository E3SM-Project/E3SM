#ifndef EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP
#define EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream
{

/*
 * This diagnostic will produce the vertical layer height at midpoint.
 */

class VerticalLayerDiagnostic : public AtmosphereDiagnostic
{
public:
  using Pack          = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using PF            = scream::PhysicsFunctions<DefaultDevice>;
  using KT            = KokkosTypes<DefaultDevice>;
  using MemberType    = typename KT::MemberType;
  using view_2d       = typename KT::template view_2d<Pack>;

  // Constructors
  VerticalLayerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

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

  // Temporary view to set z_int and z_mid in compute diagnostic
  view_2d m_z_int;
  view_2d m_z_mid;

  // Store the diagnostic name. This will dictate which
  // field in the computation is output (dz, z_int, or z_mid).
  std::string m_diag_name;

  // If z_int or z_mid is computed, determine whether the BC
  // is from sea level or not (from topography data).
  bool m_from_sea_level;

  // Store if the diagnostic output field exists on interface values
  bool m_is_interface_layout;

}; // class VerticalLayerDiagnostic

} //namespace scream

#endif // EAMXX_VERTICAL_LAY_MID_DIAGNOSTIC_HPP
