#ifndef EAMXX_WATER_PATH_DIAGNOSTIC_HPP
#define EAMXX_WATER_PATH_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class WaterPathDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  WaterPathDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const override { return "WaterPath"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl ();
protected:

  // Keep track of field dimensions
  int m_num_cols;
  int m_num_levs;

  std::string m_qname;
  std::string m_kind;
}; // class WaterPathDiagnostic

} //namespace scream

#endif // EAMXX_WATER_PATH_DIAGNOSTIC_HPP
