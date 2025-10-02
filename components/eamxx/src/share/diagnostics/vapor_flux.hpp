#ifndef EAMXX_VAPOR_FLUX_DIAGNOSTIC_HPP
#define EAMXX_VAPOR_FLUX_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the zonal or meridional water vapor flux.
 */

class VaporFluxDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  VaporFluxDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const override { return "VaporFlux"; }

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

  int m_component;

  std::string m_name;
};

} //namespace scream

#endif // EAMXX_VAPOR_FLUX_HPP
