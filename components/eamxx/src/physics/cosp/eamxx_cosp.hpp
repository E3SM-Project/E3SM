#ifndef SCREAM_COSP_HPP
#define SCREAM_COSP_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the calculation of COSP diagnostics
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class Cosp : public AtmosphereProcess
{
public:

  // Constructors
  Cosp (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "cosp"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   ();

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_subcols = 50;
  Int m_num_levs;
  Int m_num_tau;
  Int m_num_ctp;
  Int m_num_cth;

  std::shared_ptr<const AbstractGrid> m_grid;
}; // class Cosp

} // namespace scream

#endif // SCREAM_COSP_HPP
