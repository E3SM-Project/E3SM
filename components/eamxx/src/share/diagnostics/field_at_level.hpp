#ifndef EAMXX_FIELD_AT_LEVEL_HPP
#define EAMXX_FIELD_AT_LEVEL_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce a slice of a field at a given vertical level index
 */

class FieldAtLevel : public AtmosphereDiagnostic
{
public:
  // Constructors
  FieldAtLevel (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "FieldAtLevel"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void initialize_impl (const RunType /*run_type*/);

  std::string   m_diag_name;
  int           m_field_level;
}; // class FieldAtLevel

} //namespace scream

#endif // EAMXX_FIELD_AT_LEVEL_HPP
