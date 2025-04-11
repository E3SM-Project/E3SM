#ifndef EAMXX_VIRTUAL_TEMP_DIAGNOSTIC_HPP
#define EAMXX_VIRTUAL_TEMP_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the virtual temperature.
 */

class VirtualTemperatureDiagnostic : public AtmosphereDiagnostic
{
public:

  // Constructors
  VirtualTemperatureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "VirtualTemperature"; }

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

}; // class VirtualTemperatureDiagnostic

} //namespace scream

#endif // EAMXX_VIRTUAL_TEMP_DIAGNOSTIC_HPP
