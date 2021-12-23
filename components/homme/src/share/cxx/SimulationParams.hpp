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

  void print();

  TimeStepType  time_step_type;
  MoistDry      moisture;
  RemapAlg      remap_alg;
  TestCase      test_case;
  ForcingAlg    ftype;
  AdvectionForm theta_adv_form; // Only for theta model

  int           rsplit, dt_remap_factor;
  int           qsplit, dt_tracer_factor;
  int           qsize;


  int       limiter_option; // TODO: convert to enum

  bool      prescribed_wind;

  int       state_frequency;
  bool      disable_diagnostics;
  int       transport_alg;
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
  int       hypervis_subcycle_tom;
  double    hypervis_scaling;
  double    nu_ratio1, nu_ratio2; //control balance between div and vort components in vector laplace
  double    rearth; //propagated then to Geometry and SphereOps

  // Use this member to check whether the struct has been initialized
  bool      params_set;
};

inline void SimulationParams::print () {

  printf ("\n************** CXX SimulationParams **********************\n\n");
  printf ("   time_step_type: %d\n", etoi(time_step_type));
  printf ("   moisture: %s\n", moisture==MoistDry::DRY ? "dry" : "moist");
  printf ("   remap_alg: %d\n", etoi(remap_alg));
  printf ("   test case: %d\n", etoi(test_case));
  printf ("   ftype: %d\n", etoi(ftype));
  printf ("   theta_adv_form: %d\n", etoi(theta_adv_form));
  printf ("   rsplit: %d\n", rsplit);
  printf ("   qsplit: %d\n", qsplit);
  printf ("   qsize: %d\n", qsize);
  printf ("   limiter_option: %d\n", limiter_option);
  printf ("   state_frequency: %d\n", state_frequency);
  printf ("   dcmip16_mu: %f\n", dcmip16_mu);
  printf ("   nu: %f\n", nu);
  printf ("   nu_p: %f\n", nu_p);
  printf ("   nu_q: %f\n", nu_q);
  printf ("   nu_s: %f\n", nu_s);
  printf ("   nu_top: %f\n", nu_top);
  printf ("   nu_div: %f\n", nu_div);
  printf ("   hypervis_order: %d\n", hypervis_order);
  printf ("   hypervis_subcycle: %d\n", hypervis_subcycle);
  printf ("   hypervis_subcycle_tom: %d\n", hypervis_subcycle_tom);
  printf ("   hypervis_scaling: %f\n", hypervis_scaling);
  printf ("   nu_ratio1: %f\n", nu_ratio1);
  printf ("   nu_ratio2: %f\n", nu_ratio2);
  printf ("   use_cpstar: %s\n", (use_cpstar ? "yes" : "no"));
  printf ("   transport_alg: %d\n", transport_alg);
  printf ("   disable_diagnostics: %s\n", (disable_diagnostics ? "yes" : "no"));
  printf ("   theta_hydrostatic_mode: %s\n", (theta_hydrostatic_mode ? "yes" : "no"));
  printf ("   prescribed_wind: %s\n", (prescribed_wind ? "yes" : "no"));
  printf ("   rearth: %f\n", rearth);
  printf ("\n**********************************************************\n");
}

} // namespace Homme

#endif // HOMMEXX_SIMULATION_PARAMS_HPP
