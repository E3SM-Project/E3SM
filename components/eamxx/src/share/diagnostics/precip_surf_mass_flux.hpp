#ifndef EAMXX_PRECIP_SURF_MASS_FLUX_HPP
#define EAMXX_PRECIP_SURF_MASS_FLUX_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the precipitation (ice) surface mass flux.
 */

class PrecipSurfMassFlux : public AtmosphereDiagnostic
{
public:
  // Constructors
  PrecipSurfMassFlux (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "PrecipSurfMassFlux"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  Int m_num_cols;

  static constexpr int s_ice = 1;
  static constexpr int s_liq = 2;

  int m_type;
  std::string m_name;
};

}

#endif // EAMXX_PRECIP_SURF_MASS_FLUX_HPP
