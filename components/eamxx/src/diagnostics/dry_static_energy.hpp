#ifndef EAMXX_DRY_STATIC_ENERGY_DIAGNOSTIC_HPP
#define EAMXX_DRY_STATIC_ENERGY_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

#include "share/eamxx_types.hpp"

namespace scream
{

class DryStaticEnergyDiagnostic : public AtmosphereDiagnostic
{
public:
  using KT            = KokkosTypes<DefaultDevice>;
  using view_2d       = typename KT::template view_2d<Real>;

  // Constructors
  DryStaticEnergyDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "DryStaticEnergy"; }

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

  // Temporary view to set dz in compute diagnostic
  view_2d m_tmp_mid;
  view_2d m_tmp_int;

}; // class DryStaticEnergyDiagnostic

} //namespace scream

#endif // EAMXX_DRY_STATIC_ENERGY_DIAGNOSTIC_HPP
