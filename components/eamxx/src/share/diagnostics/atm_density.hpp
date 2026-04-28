#ifndef EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP
#define EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

class AtmDensity : public AbstractDiagnostic
{
public:
  // Constructors
  AtmDensity (const ekat::Comm& comm, const ekat::ParameterList& params,
              const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "AtmosphereDensity"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
  void initialize_impl ();
};

} //namespace scream

#endif // EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP
