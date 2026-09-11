#ifndef SCREAM_PRESCRIBED_AEROSOL_HPP
#define SCREAM_PRESCRIBED_AEROSOL_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

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
  using KT           = ekat::KokkosTypes<DefaultDevice>;
  using view_2d      = typename KT::template view_2d<Real>;
public:
  // Constructors
  SPA (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "spa"; }

  // Set the grid
  void create_requests ();

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   () { /* Nothing to do */ }

  protected:

  std::shared_ptr<const AbstractGrid>   m_model_grid;

  std::shared_ptr<DataInterpolation>    m_data_interpolation;

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // number of shortwave and longwave bands
  int nswbands_, nlwbands_;

  // layer thickness
  view_2d dz_;

}; // class SPA

} // namespace scream

#endif // SCREAM_SPA_HPP
