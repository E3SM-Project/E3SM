#ifndef EAMXX_SEA_LEVEL_PRESSURE_DIAGNOSTIC_HPP
#define EAMXX_SEA_LEVEL_PRESSURE_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class SeaLevelPressureDiagnostic : public AtmosphereDiagnostic
{
public:

  // Constructors
  SeaLevelPressureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic
  std::string name () const { return "SeaLevelPressure"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  // Keep track of field dimensions
  Int m_num_cols;
  Int m_num_levs;

}; // class SeaLevelPressureDiagnostic

} //namespace scream

#endif // EAMXX_SEA_LEVEL_PRESSURE_DIAGNOSTIC_HPP
