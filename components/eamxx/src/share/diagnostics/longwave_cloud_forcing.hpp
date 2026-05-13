#ifndef EAMXX_LONGWAVE_CLOUD_FORCING_HPP
#define EAMXX_LONGWAVE_CLOUD_FORCING_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce longwave cloud forcing.
 */

class LongwaveCloudForcing : public AbstractDiagnostic
{
public:
  // Constructors
  LongwaveCloudForcing (const ekat::Comm& comm, const ekat::ParameterList& params,
                        const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "LongwaveCloudForcing"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl ();
};

} //namespace scream

#endif // EAMXX_LONGWAVE_CLOUD_FORCING_HPP
