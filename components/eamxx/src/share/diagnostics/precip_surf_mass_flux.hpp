#ifndef EAMXX_PRECIP_SURF_MASS_FLUX_HPP
#define EAMXX_PRECIP_SURF_MASS_FLUX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the precipitation (ice) surface mass flux.
 */

class PrecipSurfMassFlux : public AbstractDiagnostic
{
public:
  // Constructors
  PrecipSurfMassFlux (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "PrecipSurfMassFlux"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  static constexpr int s_ice = 1;
  static constexpr int s_liq = 2;

  int m_type;
  std::string m_name;
};

}

#endif // EAMXX_PRECIP_SURF_MASS_FLUX_HPP
