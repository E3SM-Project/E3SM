#ifndef EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP
#define EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

#include <ekat/ekat_pack.hpp>

namespace scream
{

/*
 * This diagnostic will produce a slice of a field at a given pressure level
 */

class FieldAtPressureLevel : public AtmosphereDiagnostic
{
public:

  // Constructors
  FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "FieldAtPressureLevel"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:
  void initialize_impl (const RunType /*run_type*/);

  std::string         m_pressure_name;
  std::string         m_field_name;
  std::string         m_diag_name;

  Real                m_pressure_level;
  int                 m_num_levs;
  Real                m_mask_val;

}; // class FieldAtPressureLevel

} //namespace scream

#endif // EAMXX_FIELD_AT_PRESSURE_LEVEL_HPP
