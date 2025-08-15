#ifndef EAMXX_SURF_LATENT_HEAT_FLUX_HPP
#define EAMXX_SURF_LATENT_HEAT_FLUX_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the surface latent heat flux.
 */

class SurfaceUpwardLatentHeatFlux : public AtmosphereDiagnostic
{
public:
  // Constructors
  SurfaceUpwardLatentHeatFlux (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "SurfaceUpwardLatentHeatFlux"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  Int m_num_cols;

  int m_type;
  std::string m_name;
  std::string cf_long_name;
};

} // namespace scream
#endif
