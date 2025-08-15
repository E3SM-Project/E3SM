#ifndef EAMXX_WIND_SPEED_HPP
#define EAMXX_WIND_SPEED_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will compute the magnitute of the horiz_winds vector
 */

class WindSpeed : public AtmosphereDiagnostic
{
public:
  // Constructors
  WindSpeed (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const override { return "wind_speed"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) override;

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_diagnostic_impl () override;

  int m_ncols;
  int m_nlevs;
};

} //namespace scream

#endif // EAMXX_WIND_SPEED_HPP
