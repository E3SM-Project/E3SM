#ifndef EAMXX_EXNER_DIAGNOSTIC_HPP
#define EAMXX_EXNER_DIAGNOSTIC_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the exner function.
 */

class Exner : public AbstractDiagnostic
{
public:
  // Constructors
  Exner (const ekat::Comm& comm, const ekat::ParameterList& params,
         const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "Exner"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void initialize_impl ();
  void compute_diagnostic_impl ();
};

} //namespace scream

#endif // EAMXX_EXNER_DIAGNOSTIC_HPP
