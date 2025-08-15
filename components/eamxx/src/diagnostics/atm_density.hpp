#ifndef EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP
#define EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

class AtmDensityDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  AtmDensityDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "AtmosphereDensity"; }

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

}; // class AtmDensityDiagnostic

} //namespace scream

#endif // EAMXX_ATM_DENSITY_DIAGNOSTIC_HPP
