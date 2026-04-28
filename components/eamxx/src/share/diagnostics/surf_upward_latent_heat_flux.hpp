#ifndef EAMXX_SURF_LATENT_HEAT_FLUX_HPP
#define EAMXX_SURF_LATENT_HEAT_FLUX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the surface latent heat flux.
 */

class SurfaceUpwardLatentHeatFlux : public AbstractDiagnostic
{
public:
  // Constructors
  SurfaceUpwardLatentHeatFlux (const ekat::Comm& comm, const ekat::ParameterList& params,
                                const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "SurfaceUpwardLatentHeatFlux"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
};

} // namespace scream
#endif
