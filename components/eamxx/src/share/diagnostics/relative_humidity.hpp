#ifndef EAMXX_RELATIVE_HUMIDITY_DIAGNOSTIC_HPP
#define EAMXX_RELATIVE_HUMIDITY_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

class RelativeHumidityDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  RelativeHumidityDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "RelativeHumidity"; }

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

}; // class RelativeHumidityDiagnostic

} //namespace scream

#endif // EAMXX_RELATIVE_HUMIDITY_DIAGNOSTIC_HPP
