#ifndef EAMXX_EXNER_DIAGNOSTIC_HPP
#define EAMXX_EXNER_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will produce the exner function.
 */

class ExnerDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  ExnerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "Exner"; }

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

}; // class ExnerDiagnostic

} //namespace scream

#endif // EAMXX_EXNER_DIAGNOSTIC_HPP
