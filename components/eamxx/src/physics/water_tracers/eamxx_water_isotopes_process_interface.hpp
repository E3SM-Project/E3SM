#ifndef SCREAM_WATER_ISOTOPES_HPP
#define SCREAM_WATER_ISOTOPES_HPP

#include "eamxx_water_tracers_process_interface.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
{
/*
 * The class responsible to handle water isotope transport and fractionation
 *
 * This process extends WaterTracers to add equilibrium and kinetic fractionation
 * during phase changes for water isotope species (e.g., HDO, H2-18O, HTO).
 *
 * By inheriting from WaterTracers, this class reuses all tracer field handling
 * and only needs to override specific fractionation hooks.
 *
 * Note: This is a stub implementation that registers the process. Fractionation
 * physics will be added in later specs of the water isotope campaign.
*/

class WaterIsotopes : public WaterTracers
{
public:
  // Constructors
  WaterIsotopes (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Override the name to distinguish from base WaterTracers
  std::string name () const override { return "water_isotopes"; }

protected:

  // Override run_impl to add fractionation physics
  // For now, just calls base class implementation
  void run_impl (const double dt) override;

  // TODO (later campaign specs): Add virtual hooks for fractionation processes
  // e.g., apply_equilibrium_fractionation(), apply_kinetic_fractionation()

}; // class WaterIsotopes

} // namespace scream

#endif // SCREAM_WATER_ISOTOPES_HPP
