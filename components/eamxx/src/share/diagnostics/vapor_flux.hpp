#ifndef EAMXX_VAPOR_FLUX_HPP
#define EAMXX_VAPOR_FLUX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the zonal or meridional water vapor flux.
 */

class VaporFlux : public AbstractDiagnostic
{
public:
  // Constructors
  VaporFlux (const ekat::Comm& comm, const ekat::ParameterList& params,
             const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const override { return "VaporFlux"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl () override;
protected:

  int m_component;

  std::string m_name;
};

} //namespace scream

#endif // EAMXX_VAPOR_FLUX_HPP
