#ifndef EAMXX_PRECIP_TOTAL_SURF_MASS_DIAGNOSTIC_HPP
#define EAMXX_PRECIP_TOTAL_SURF_MASS_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream
{

/*
 * This diagnostic will produce the total precipitation surface mass flux.
 */

class PrecipTotalSurfMassFluxDiagnostic : public AtmosphereDiagnostic
{
public:
  using PC = scream::physics::Constants<Real>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // Constructors
  PrecipTotalSurfMassFluxDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // The name of the diagnostic
  std::string name () const { return "PrecipTotalSurfMassFlux"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  Int m_num_cols;
}; // class PrecipTotalSurfMassFluxDiagnostic

} //namespace scream

#endif // EAMXX_PRECIP_TOTAL_SURF_MASS_DIAGNOSTIC_HPP
