#ifndef SCREAM_WATER_TRACERS_HPP
#define SCREAM_WATER_TRACERS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
{
/*
 * The class responsible to handle water tracer transport through the atmosphere
 *
 * This process manages additional water species that track through the model
 * without undergoing fractionation. Water isotopes (a special case with
 * fractionation) are handled by the WaterIsotopes subclass.
 *
 * Note: This is a stub implementation that registers the process and sets up
 * the basic infrastructure. Field definitions and physics implementation are
 * deferred to later specs in the water isotope campaign.
*/

class WaterTracers : public AtmosphereProcess
{
public:
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // Constructors
  WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  virtual std::string name () const { return "water_tracers"; }

  // Create grid-dependent field requests
  void create_requests ();

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   ();

  // Keep track of field dimensions
  int m_num_cols;
  int m_num_levs;

  // Number of tracers to track (from parameter list)
  int m_tracer_count;

  std::shared_ptr<const AbstractGrid> m_grid;

}; // class WaterTracers

} // namespace scream

#endif // SCREAM_WATER_TRACERS_HPP
