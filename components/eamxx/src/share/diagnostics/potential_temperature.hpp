#ifndef EAMXX_POTENTIAL_TEMP_HPP
#define EAMXX_POTENTIAL_TEMP_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class PotentialTemperature : public AbstractDiagnostic
{
public:
  // Constructors
  PotentialTemperature (const ekat::Comm& comm, const ekat::ParameterList& params,
                        const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const override { return "PotentialTemperature"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl () override;
protected:

  // What type of potential temperature to compute
  std::string m_ptype;
};

} //namespace scream

#endif // EAMXX_POTENTIAL_TEMP_HPP
