#ifndef EAMXX_LONGWAVE_CLOUD_FORCING_DIAGNOSTIC_HPP
#define EAMXX_LONGWAVE_CLOUD_FORCING_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce longwave cloud forcing.
 */

class LongwaveCloudForcingDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  LongwaveCloudForcingDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "LongwaveCloudForcing"; }

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

}; // class LongwaveCloudForcingDiagnostic

} //namespace scream

#endif // EAMXX_LONGWAVE_CLOUD_FORCING_DIAGNOSTIC_HPP
