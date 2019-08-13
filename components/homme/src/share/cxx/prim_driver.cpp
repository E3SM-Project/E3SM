/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/


#include "Context.hpp"
#include "Diagnostics.hpp"
#include "Elements.hpp"
#include "Tracers.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "ErrorDefs.hpp"
#include "CamForcing.hpp"
#include "profiling.hpp"

#include <iostream>

namespace Homme
{

void prim_step (const Real, const bool);
void vertical_remap (const Real);
void apply_test_forcing ();
void update_q (const int np1_qdp, const int np1);

extern "C" {

void prim_run_subcycle_c (const Real& dt, int& nstep, int& nm1, int& n0, int& np1, const int& next_output_step)
{
  GPTLstart("tl-sc prim_run_subcycle_c");

  // Get simulation params
  SimulationParams& params = Context::singleton().get_simulation_params();
  assert(params.params_set);

  // Get time info and compute dt for tracers and remap
  TimeLevel& tl = Context::singleton().get_time_level();
  const Real dt_q = dt*params.qsplit;
  Real dt_remap = dt_q;
  int nstep_end = tl.nstep + params.qsplit;
  if (params.rsplit>0) {
    dt_remap  = dt_q*params.rsplit;
    nstep_end = tl.nstep + params.qsplit*params.rsplit;
  }

  // Check if needed to compute diagnostics or energy
  bool compute_diagnostics = false;
  if (nstep_end%params.state_frequency==0 || nstep_end==tl.nstep0 ||
      nstep_end>=next_output_step) {
    compute_diagnostics = true;
  }

  if (params.disable_diagnostics) {
    compute_diagnostics = false;
  }

  if (compute_diagnostics) {
    Diagnostics& diags = Context::singleton().get_diagnostics();
    diags.prim_diag_scalars(true,2);
    diags.prim_energy_halftimes(true,2);
  }

  tl.update_tracers_levels(params.qsplit);

#ifndef CAM
  apply_test_forcing ();
#endif


  // Apply forcing.
  // In standalone mode, params.ftype == ForcingAlg::FORCING_DEBUG
  // Corresponds to ftype == 0 in Fortran
  if(params.ftype == ForcingAlg::FORCING_DEBUG) {
    apply_cam_forcing(dt_remap);
  }
  // Corresponds to ftype == 2 in Fortran
  else if(params.ftype == ForcingAlg::FORCING_2) {
    apply_cam_forcing_dynamics(dt_remap);
  }

  if (compute_diagnostics) {
    Diagnostics& diags = Context::singleton().get_diagnostics();
    diags.prim_energy_halftimes(true,0);
    diags.prim_diag_scalars(true,0);
  }

  // Initialize dp3d from ps
  GPTLstart("tl-sc dp3d-from-ps");
  Elements& elements = Context::singleton().get_elements();
  HybridVCoord& hvcoord = Context::singleton().get_hvcoord();
  const auto hybrid_ai_delta = hvcoord.hybrid_ai_delta;
  const auto hybrid_bi_delta = hvcoord.hybrid_bi_delta;
  const auto ps0 = hvcoord.ps0;
  const auto ps_v = elements.m_ps_v;
  {
    const auto dp3d = elements.m_dp3d;
    const auto n0 = tl.n0;
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace> (0,elements.num_elems()*NP*NP*NUM_LEV),
                         KOKKOS_LAMBDA(const int idx) {
      const int ie   = ((idx / NUM_LEV) / NP) / NP;
      const int igp  = ((idx / NUM_LEV) / NP) % NP;
      const int jgp  =  (idx / NUM_LEV) % NP;
      const int ilev =   idx % NUM_LEV;

      dp3d(ie,n0,igp,jgp,ilev) = hybrid_ai_delta[ilev]*ps0
                               + hybrid_bi_delta[ilev]*ps_v(ie,n0,igp,jgp);
    });
  }
  ExecSpace::fence();
  GPTLstop("tl-sc dp3d-from-ps");

  // Loop over rsplit vertically lagrangian timesteps
  GPTLstart("tl-sc prim_step-loop");
  prim_step(dt,compute_diagnostics);
  for (int r=1; r<params.rsplit; ++r) {
    tl.update_dynamics_levels(UpdateType::LEAPFROG);
    prim_step(dt,false);
  }
  GPTLstop("tl-sc prim_step-loop");

  ////////////////////////////////////////////////////////////////////////
  // apply vertical remap
  // always for tracers
  // if rsplit>0:  also remap dynamics and compute reference level ps_v
  ////////////////////////////////////////////////////////////////////////
  tl.update_tracers_levels(params.qsplit);
  GPTLstart("tl-sc vertical_remap");
  vertical_remap(dt_remap);
  GPTLstop("tl-sc vertical_remap");

  ////////////////////////////////////////////////////////////////////////
  // time step is complete.  update some diagnostic variables:
  // Q    (mixing ratio)
  ////////////////////////////////////////////////////////////////////////
  update_q(tl.np1_qdp,tl.np1);

  if (compute_diagnostics) {
    Diagnostics& diags = Context::singleton().get_diagnostics();
    diags.prim_diag_scalars(false,1);
    diags.prim_energy_halftimes(false,1);
  }

  // Update dynamics time levels
  tl.update_dynamics_levels(UpdateType::LEAPFROG);

  // Update the timelevel info to pass back to fortran
  nstep = tl.nstep;
  nm1   = tl.nm1;
  n0    = tl.n0;
  np1   = tl.np1;

  GPTLstop("tl-sc prim_run_subcycle_c");
}

} // extern "C"

void apply_test_forcing () {
  // Get simulation params
  SimulationParams& params = Context::singleton().get_simulation_params();

  if (params.test_case==TestCase::DCMIP2012_TEST2_1 ||
      params.test_case==TestCase::DCMIP2012_TEST2_2) {
    Errors::runtime_abort("Test case not yet available in C++ build.\n",
                          Errors::err_not_implemented);
  }
}

void update_q (const int np1_qdp, const int np1)
{
  // Get simulation params
  SimulationParams& params = Context::singleton().get_simulation_params();
  assert(params.params_set);

  // Get hybrid vertical coordinate
  HybridVCoord& hvcoord = Context::singleton().get_hvcoord();
  auto hyai_delta = hvcoord.hybrid_ai_delta;
  auto hybi_delta = hvcoord.hybrid_bi_delta;
  const Real ps0 = hvcoord.ps0;

  // Get ps_v from Elements
  Elements& elements = Context::singleton().get_elements();
  auto ps_v = elements.m_ps_v;

  // Get the tracers concentration and mass from Tracers
  Tracers& tracers = Context::singleton().get_tracers();
  auto qdp = tracers.qdp;
  auto Q = tracers.Q;

  // Update the device copy of Q, stored in Tracers
  const int num_elems = elements.num_elems();
  const int qsize = params.qsize;
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,num_elems*qsize*NP*NP*NUM_LEV),
                       KOKKOS_LAMBDA(const int idx) {
    const int ie    =  idx / (qsize*NP*NP*NUM_LEV);
    const int iq    = (idx / (NP*NP*NUM_LEV)) % qsize;
    const int igp   = (idx / (NP*NUM_LEV)) % NP;
    const int jgp   = (idx / NUM_LEV) % NP;
    const int level =  idx % NUM_LEV;

    const Scalar dp = hyai_delta(level)*ps0 + hybi_delta(level)*ps_v(ie,np1,igp,jgp);
    Q(ie,iq,igp,jgp,level) = qdp(ie,np1_qdp,iq,igp,jgp,level)/dp;
  });
}

} // namespace Homme
