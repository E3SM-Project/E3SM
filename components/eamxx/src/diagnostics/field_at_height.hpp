#ifndef EAMXX_FIELD_AT_HEIGHT_HPP
#define EAMXX_FIELD_AT_HEIGHT_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce a slice of a field at a given height (above surface)
 */

class FieldAtHeight : public AtmosphereDiagnostic
{
public:

  // Constructors
  FieldAtHeight (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "FieldAtHeight"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void initialize_impl (const RunType /*run_type*/);

  std::string         m_diag_name;
  std::string         m_z_name;
  std::string         m_z_suffix;
  std::string         m_field_name;

  Real                m_z;
};

} //namespace scream

#endif // EAMXX_FIELD_AT_HEIGHT_HPP
