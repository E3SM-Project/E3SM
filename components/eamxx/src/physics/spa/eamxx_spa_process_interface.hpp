#ifndef SCREAM_PRESCRIBED_AEROSOL_HPP
#define SCREAM_PRESCRIBED_AEROSOL_HPP

#include "share/atm_process/atmosphere_process.hpp"

namespace scream
{

class DataInterpolation;

/*
 * The class responsible to handle the calculation of the subgrid cloud fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class SPA : public AtmosphereProcess
{
public:
  // Constructors
  SPA (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "spa"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   () { /* Nothing to do */ }

  std::shared_ptr<const AbstractGrid>   m_model_grid;

  std::shared_ptr<DataInterpolation>    m_data_interpolation;
}; // class SPA

} // namespace scream

#endif // SCREAM_SPA_HPP
