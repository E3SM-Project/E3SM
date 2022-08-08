/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_SIMULATION_PARAMS_HPP
#define HOMMEXX_SIMULATION_PARAMS_HPP

#include "HommexxEnums.hpp"

#include <iostream>

namespace Homme
{

/*
 * A struct to hold simulation parameters.
 *
 */
struct SimulationParams
{
  void print(std::ostream& out = std::cout);

  TimeStepType  time_step_type;
  MoistDry      moisture;
  RemapAlg      remap_alg;
  TestCase      test_case;
  ForcingAlg    ftype = ForcingAlg::FORCING_OFF;
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
  int       nsplit = 0;
  int       nsplit_iteration;
  double    rearth; //propagated then to Geometry and SphereOps

  // Use this member to check whether the struct has been initialized
  bool      params_set = false;
};

inline void SimulationParams::print (std::ostream& out) {

  out << "\n************** CXX SimulationParams **********************\n\n";
  out << "   time_step_type: %d\n", etoi(time_step_type);
  out << "   moisture: %s\n", moisture==MoistDry::DRY ? "dry" : "moist";
  out << "   remap_alg: %d\n", etoi(remap_alg);
  out << "   test case: %d\n", etoi(test_case);
  out << "   ftype: %d\n", etoi(ftype);
  out << "   theta_adv_form: %d\n", etoi(theta_adv_form);
  out << "   rsplit: %d\n", rsplit;
  out << "   qsplit: %d\n", qsplit;
  out << "   qsize: %d\n", qsize;
  out << "   limiter_option: %d\n", limiter_option;
  out << "   state_frequency: %d\n", state_frequency;
  out << "   dcmip16_mu: %f\n", dcmip16_mu;
  out << "   nu: %f\n", nu;
  out << "   nu_p: %f\n", nu_p;
  out << "   nu_q: %f\n", nu_q;
  out << "   nu_s: %f\n", nu_s;
  out << "   nu_top: %f\n", nu_top;
  out << "   nu_div: %f\n", nu_div;
  out << "   hypervis_order: %d\n", hypervis_order;
  out << "   hypervis_subcycle: %d\n", hypervis_subcycle;
  out << "   hypervis_subcycle_tom: %d\n", hypervis_subcycle_tom;
  out << "   hypervis_scaling: %f\n", hypervis_scaling;
  out << "   nu_ratio1: %f\n", nu_ratio1;
  out << "   nu_ratio2: %f\n", nu_ratio2;
  out << "   use_cpstar: %s\n", (use_cpstar ? "yes" : "no");
  out << "   transport_alg: %d\n", transport_alg;
  out << "   disable_diagnostics: %s\n", (disable_diagnostics ? "yes" : "no");
  out << "   theta_hydrostatic_mode: %s\n", (theta_hydrostatic_mode ? "yes" : "no");
  out << "   prescribed_wind: %s\n", (prescribed_wind ? "yes" : "no");
  out << "   nsplit: %d\n", nsplit;
  out << "   rearth: %f\n", rearth;
  out << "\n**********************************************************\n";
}

} // namespace Homme

#endif // HOMMEXX_SIMULATION_PARAMS_HPP
