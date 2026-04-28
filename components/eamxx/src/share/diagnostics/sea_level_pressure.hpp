#ifndef EAMXX_SEA_LEVEL_PRESSURE_HPP
#define EAMXX_SEA_LEVEL_PRESSURE_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class SeaLevelPressure : public AbstractDiagnostic
{
public:
  // Constructors
  SeaLevelPressure (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "SeaLevelPressure"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
};

} //namespace scream

#endif // EAMXX_SEA_LEVEL_PRESSURE_HPP
