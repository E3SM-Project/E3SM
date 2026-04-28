#ifndef EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP
#define EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce a slice of a field at a given pressure level
 */

class FieldAtPressureLevel : public AbstractDiagnostic
{
public:

  // Constructors
  FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params,
                        const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "FieldAtPressureLevel"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void initialize_impl ();

  std::string         m_pressure_name;
  std::string         m_field_name;
  std::string         m_diag_name;

  Real                m_pressure_level;
};

} //namespace scream

#endif // EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP
