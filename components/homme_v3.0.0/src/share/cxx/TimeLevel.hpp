/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_TIME_LEVEL_HPP
#define HOMMEXX_TIME_LEVEL_HPP

#include "ErrorDefs.hpp"
#include "HommexxEnums.hpp"

namespace Homme
{

struct TimeLevel
{
  // Dynamics relative time levels
  int nm1;
  int n0;
  int np1;

  // Absolute time level since simulation start
  int nstep;

  // Absolute time level of first complete leapfrog timestep
  int nstep0;

  // Tracers relative time levels
  int n0_qdp;
  int np1_qdp;

  // Time passed since the start of the simulation
  // TODO: I think this is used only with when CAM is defined
  double tevolve;

  void update_dynamics_levels (UpdateType type) {
    int tmp;
    switch(type) {
      case UpdateType::LEAPFROG:
        tmp = np1;
        np1 = nm1;
        nm1 = n0;
        n0  = tmp;
        break;
      case UpdateType::FORWARD:
        tmp = np1;
        np1 = n0;
        n0  = tmp;
        break;
      default:
        Errors::runtime_abort("Unknown time level update type",
                              Errors::err_unknown_option);
    }
    ++nstep;
  }

  void update_tracers_levels (const int qsplit) {
    int i_temp = nstep/qsplit;
    if (i_temp%2 == 0) {
      n0_qdp  = 0;
      np1_qdp = 1;
    } else {
      n0_qdp  = 1;
      np1_qdp = 0;
    }
  }
};

} // namespace Homme

#endif // HOMMEXX_TIME_LEVEL_HPP
