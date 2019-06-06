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
  SimulationParams() : ftype(ForcingAlg::FORCING_OFF), params_set(false) {}

  TimeStepType  time_step_type;
  MoistDry      moisture;
  RemapAlg      remap_alg;
  TestCase      test_case;
  ForcingAlg    ftype;
  AdvectionForm theta_adv_form; // Only for theta model

  int           rsplit;
  int           qsplit;
  int           qsize;


  int       limiter_option; // TODO: convert to enum

  bool      prescribed_wind;

  int       state_frequency;
  bool      disable_diagnostics;
  bool      use_semi_lagrangian_transport;
  bool      use_cpstar;
  bool      theta_hydrostatic_mode;   // Only for theta model

  double    dcmip16_mu;               // Only for theta model
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
