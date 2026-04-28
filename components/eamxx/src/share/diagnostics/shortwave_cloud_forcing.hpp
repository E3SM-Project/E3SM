#ifndef EAMXX_SHORTWAVE_CLOUD_FORCING_HPP
#define EAMXX_SHORTWAVE_CLOUD_FORCING_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the shortwave cloud forcing.
 */

class ShortwaveCloudForcing : public AbstractDiagnostic
{
public:
  // Constructors
  ShortwaveCloudForcing (const ekat::Comm& comm, const ekat::ParameterList& params,
                         const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "ShortwaveCloudForcing"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
};

} //namespace scream

#endif // EAMXX_SHORTWAVE_CLOUD_FORCING_HPP
