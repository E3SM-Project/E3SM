#ifndef EAMXX_GW_EMULATOR_HPP
#define EAMXX_GW_EMULATOR_HPP

#include "share/atm_process/atmosphere_process.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{
/*
 * A single-column emulator of the Gravity-Wave parametrization
 *
*/

class GWEmulator : public AtmosphereProcess
{
 public:
  // Constructors
  GWEmulator (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const override { return "gw-emulator"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager) override;

protected:

  void initialize_impl (const RunType run_type) override;
  void run_impl        (const double dt);
  void finalize_impl   ();

  std::shared_ptr<const AbstractGrid> m_grid;
};

} // namespace scream

#endif // EAMXX_GW_EMULATOR_HPP
