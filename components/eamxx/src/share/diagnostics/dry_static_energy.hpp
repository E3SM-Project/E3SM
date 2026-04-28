#ifndef EAMXX_DRY_STATIC_ENERGY_HPP
#define EAMXX_DRY_STATIC_ENERGY_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

class DryStaticEnergy : public AbstractDiagnostic
{
public:
  // Constructors
  DryStaticEnergy (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "DryStaticEnergy"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
};

} //namespace scream

#endif // EAMXX_DRY_STATIC_ENERGY_HPP
