#ifndef EAMXX_WIND_SPEED_HPP
#define EAMXX_WIND_SPEED_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic will compute the magnitute of the horiz_winds vector
 */

class WindSpeed : public AbstractDiagnostic
{
public:
  // Constructors
  WindSpeed (const ekat::Comm& comm, const ekat::ParameterList& params,
             const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const override { return "wind_speed"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl () override;
};

} //namespace scream

#endif // EAMXX_WIND_SPEED_HPP
