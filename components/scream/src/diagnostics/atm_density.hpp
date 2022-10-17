#ifndef EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP
#define EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class AtmDensityDiagnostic : public AtmosphereDiagnostic
{
public:
  using Pack          = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using PF            = scream::PhysicsFunctions<DefaultDevice>;

  // Constructors
  AtmDensityDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // The name of the diagnostic
  std::string name () const { return "AtmosphereDensity"; }

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

}; // class AtmDensityDiagnostic

} //namespace scream

#endif // EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP
