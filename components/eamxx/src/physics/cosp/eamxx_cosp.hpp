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

  inline bool cosp_do(const int icosp, const int nstep) {
      // If icosp == 0, then never do cosp;
      // Otherwise, we always call cosp at the first step,
      // and afterwards we do cosp if the timestep is divisible
      // by icosp
      if (icosp == 0) {
          return false;
      } else {
          return ( (nstep == 0) || (nstep % icosp == 0) );
      }
  }


protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   ();

  // cosp frequency; positive is interpreted as number of steps, negative as number of hours
  int m_cosp_frequency;
  ekat::CaseInsensitiveString m_cosp_frequency_units;

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_subcols;
  Int m_num_levs;
  Int m_num_isccptau = 7;
  Int m_num_isccpctp = 7;
  Int m_num_cth = 16;

  std::shared_ptr<const AbstractGrid> m_grid;

}; // class Cosp

} // namespace scream

#endif // SCREAM_COSP_HPP
