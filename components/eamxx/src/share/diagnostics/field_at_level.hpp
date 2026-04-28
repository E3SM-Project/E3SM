#ifndef EAMXX_FIELD_AT_LEVEL_HPP
#define EAMXX_FIELD_AT_LEVEL_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce a slice of a field at a given vertical level index
 */

class FieldAtLevel : public AbstractDiagnostic
{
public:
  // Constructors
  FieldAtLevel (const ekat::Comm& comm, const ekat::ParameterList& params,
                const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "FieldAtLevel"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void initialize_impl ();

  std::string   m_field_name;
  int           m_field_level;
}; // class FieldAtLevel

} //namespace scream

#endif // EAMXX_FIELD_AT_LEVEL_HPP
