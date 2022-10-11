#ifndef EAMXX_DRY_STATIC_ENERGY_DIAGNOSTIC_HPP
#define EAMXX_DRY_STATIC_ENERGY_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class DryStaticEnergyDiagnostic : public AtmosphereDiagnostic
{
public:
  using Pack          = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using PF            = scream::PhysicsFunctions<DefaultDevice>;

  using KT            = KokkosTypes<DefaultDevice>;
  using MemberType    = typename KT::MemberType;
  using view_1d       = typename KT::template view_1d<Pack>;

  // Constructors
  DryStaticEnergyDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // The name of the diagnostic
  std::string name () const { return "DryStaticEnergy"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef HOMMEXX_ENABLE_GPU
public:
#endif
  void compute_diagnostic_impl ();
protected:

  // Keep track of field dimensions
  Int m_num_cols;
  Int m_num_levs;

}; // class DryStaticEnergyDiagnostic

} //namespace scream

#endif // EAMXX_DRY_STATIC_ENERGY_DIAGNOSTIC_HPP
