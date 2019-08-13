/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_SIMULATION_PARAMS_HPP
#define HOMMEXX_SIMULATION_PARAMS_HPP

#include "HommexxEnums.hpp"

namespace Homme
{

/*
 * A struct to hold simulation parameters.
 *
 */
struct SimulationParams
{
  SimulationParams() : params_set(false), ftype(ForcingAlg::FORCING_OFF) {}

  int       time_step_type; // TODO: convert to enum
  int       rsplit;
  int       qsplit;
  int       qsize;

  MoistDry  moisture;
  RemapAlg  remap_alg;
  TestCase  test_case;
  ForcingAlg ftype;

  int       limiter_option; // TODO: convert to enum

  bool      prescribed_wind;

  int       state_frequency;
  bool      disable_diagnostics;
  bool      use_semi_lagrangian_transport;
  bool      use_cpstar;

  double    nu;
  double    nu_p;
  double    nu_q;
  double    nu_s;
  double    nu_top;
  double    nu_div;
  int       hypervis_order;
  int       hypervis_subcycle;
  double    hypervis_scaling;
  double    nu_ratio1, nu_ratio2; //control balance between div and vort components in vector laplace

  // Use this member to check whether the struct has been initialized
  bool      params_set;
};

} // namespace Homme

#endif // HOMMEXX_SIMULATION_PARAMS_HPP
