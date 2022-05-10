#ifndef EAMXX_POTENTIAL_TEMP_DIAGNOSTIC_HPP
#define EAMXX_POTENTIAL_TEMP_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class PotentialTemperatureDiagnostic : public AtmosphereDiagnostic
{
public:
  // Constructors
  PotentialTemperatureDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // The name of the diagnostic
  std::string name () const { return "Potential Temperature"; } 

  // Get the required grid for the diagnostic
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // run_impl needs to be public so we can use a Kokkos LAMBDA in it.
  void run_impl        (const int dt);
protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void finalize_impl   ();

  // Keep track of field dimensions
  Int m_num_cols; 
  Int m_num_levs;

}; // class PotentialTemperatureDiagnostic

} //namespace scream

#endif // EAMXX_POTENTIAL_TEMP_DIAGNOSTIC_HPP
