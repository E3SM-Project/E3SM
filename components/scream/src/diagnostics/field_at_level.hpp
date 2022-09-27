#ifndef EAMXX_FIELD_AT_LEVEL_HPP
#define EAMXX_FIELD_AT_LEVEL_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class FieldAtLevel : public AtmosphereDiagnostic
{
public:
  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  // Constructors
  FieldAtLevel (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const { return m_field_name + "@lev_" + std::to_string(m_field_level); }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  void set_required_field_impl (const Field& f);

  // Keep track of field dimensions
  std::string         m_field_name;
  FieldLayout         m_field_layout;
  ekat::units::Units  m_field_units;

  int                 m_field_level;

}; // class FieldAtLevel

} //namespace scream

#endif // EAMXX_FIELD_AT_LEVEL_HPP
