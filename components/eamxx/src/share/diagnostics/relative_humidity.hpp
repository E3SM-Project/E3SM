#ifndef EAMXX_RELATIVE_HUMIDITY_HPP
#define EAMXX_RELATIVE_HUMIDITY_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

class RelativeHumidity : public AbstractDiagnostic
{
public:
  // Constructors
  RelativeHumidity (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "RelativeHumidity"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
};

} //namespace scream

#endif // EAMXX_RELATIVE_HUMIDITY_HPP
