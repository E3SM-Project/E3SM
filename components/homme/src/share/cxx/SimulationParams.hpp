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
  double    nu_ratio1, nu_ratio2; // control balance between div and vort components in vector laplace
  int       nsplit = 0;
  int       nsplit_iteration;
  double    scale_factor; // radius of Earth in sphere case; propagated then to Geometry and SphereOps
  double    laplacian_rigid_factor; // propagated to SphereOps
  bool      pgrad_correction;

  double    dp3d_thresh;
  double    vtheta_thresh;

  // Use this member to check whether the struct has been initialized
  bool      params_set = false;
};

inline void SimulationParams::print (std::ostream& out) {

  out << "\n************** CXX SimulationParams **********************\n\n";
  out << "   time_step_type: " << etoi(time_step_type) << "\n";
  out << "   moisture: " << (moisture==MoistDry::DRY ? "dry" : "moist") << "\n";
  out << "   remap_alg: " << etoi(remap_alg) << "\n";
  out << "   test case: " << etoi(test_case) << "\n";
  out << "   ftype: " << etoi(ftype) << "\n";
  out << "   theta_adv_form: " << etoi(theta_adv_form) << "\n";
  out << "   rsplit: " << rsplit << "\n";
  out << "   qsplit: " << qsplit << "\n";
  out << "   dt_remap_factor: " << dt_remap_factor << "\n";
  out << "   dt_tracer_factor: " << dt_tracer_factor << "\n";
  out << "   qsize: " << qsize << "\n";
  out << "   limiter_option: " << limiter_option << "\n";
  out << "   state_frequency: " << state_frequency << "\n";
  out << "   dcmip16_mu: " << dcmip16_mu << "\n";
  out << "   nu: " << nu << "\n";
  out << "   nu_p: " << nu_p << "\n";
  out << "   nu_q: " << nu_q << "\n";
  out << "   nu_s: " << nu_s << "\n";
  out << "   nu_top: " << nu_top << "\n";
  out << "   nu_div: " << nu_div << "\n";
  out << "   hypervis_order: " << hypervis_order << "\n";
  out << "   hypervis_subcycle: " << hypervis_subcycle << "\n";
  out << "   hypervis_subcycle_tom: " << hypervis_subcycle_tom << "\n";
  out << "   hypervis_scaling: " << hypervis_scaling << "\n";
  out << "   nu_ratio1: " << nu_ratio1 << "\n";
  out << "   nu_ratio2: " << nu_ratio2 << "\n";
  out << "   use_cpstar: " << (use_cpstar ? "yes" : "no") << "\n";
  out << "   transport_alg: " << transport_alg << "\n";
  out << "   disable_diagnostics: " << (disable_diagnostics ? "yes" : "no") << "\n";
  out << "   theta_hydrostatic_mode: " << (theta_hydrostatic_mode ? "yes" : "no") << "\n";
  out << "   prescribed_wind: " << (prescribed_wind ? "yes" : "no") << "\n";
  out << "   nsplit: " << nsplit << "\n";
  out << "   scale_factor: " << scale_factor << "\n";
  out << "   laplacian_rigid_factor: " << laplacian_rigid_factor << "\n";
  out << "   dp3d_thresh: " << dp3d_thresh << "\n";
  out << "   vtheta_thresh: " << vtheta_thresh << "\n";
  out << "\n**********************************************************\n";
}

} // namespace Homme

#endif // HOMMEXX_SIMULATION_PARAMS_HPP
