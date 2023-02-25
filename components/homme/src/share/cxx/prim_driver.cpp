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

namespace Homme
{

void prim_step (const Real, const bool);
void prim_step_flexible (const Real, const bool);
void vertical_remap (const Real);
void update_q (const int np1_qdp, const int np1);

extern "C" {

void initialize_dp3d_from_ps_c () {
  // Initialize dp3d from ps
  GPTLstart("tl-sc dp3d-from-ps");

  auto& context = Context::singleton();
  auto& tl = context.get<TimeLevel>();

  Elements& elements = context.get<Elements>();
  HybridVCoord& hvcoord = context.get<HybridVCoord>();
  const auto hybrid_ai_delta = hvcoord.hybrid_ai_delta;
  const auto hybrid_bi_delta = hvcoord.hybrid_bi_delta;
  const auto ps0 = hvcoord.ps0;
  const auto ps_v = elements.m_state.m_ps_v;
  {
    const auto dp3d = elements.m_state.m_dp3d;
    const auto tln0 = tl.n0;
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace> (0,elements.num_elems()*NP*NP*NUM_LEV),
                         KOKKOS_LAMBDA(const int idx) {
      const int ie   = ((idx / NUM_LEV) / NP) / NP;
      const int igp  = ((idx / NUM_LEV) / NP) % NP;
      const int jgp  =  (idx / NUM_LEV) % NP;
      const int ilev =   idx % NUM_LEV;

      dp3d(ie,tln0,igp,jgp,ilev) = hybrid_ai_delta[ilev]*ps0
                                 + hybrid_bi_delta[ilev]*ps_v(ie,tln0,igp,jgp);
    });
  }
  Kokkos::fence();
  GPTLstop("tl-sc dp3d-from-ps");
}

void prim_run_subcycle_c (const Real& dt, int& nstep, int& nm1, int& n0, int& np1, 
                          const int& next_output_step, const int& nsplit_iteration)
{
  GPTLstart("tl-sc prim_run_subcycle_c");

  auto& context = Context::singleton();

  // Get simulation params
  SimulationParams& params = context.get<SimulationParams>();
  params.nsplit_iteration = nsplit_iteration;
  assert(params.params_set);

  const bool independent_time_steps = (params.transport_alg > 0 &&
                                       params.dt_remap_factor < params.dt_tracer_factor);

  // Get time info and compute dt for tracers and remap
  TimeLevel& tl = context.get<TimeLevel>();
  const Real dt_q = dt*params.dt_tracer_factor;
  Real dt_remap;
  int nstep_end; // nstep at end of this routine
  if (params.dt_remap_factor == 0) {
    // dt_remap_factor = 0 means use eulerian code, not vert. lagrange
    dt_remap = dt_q;
    nstep_end = tl.nstep + params.dt_tracer_factor;
  } else {
    dt_remap = dt*params.dt_remap_factor;
    nstep_end = tl.nstep + (std::max(params.dt_remap_factor, params.dt_tracer_factor));
  }

  // Check if we need to compute diagnostics or energy.
  const bool compute_diagnostics =
    ( ! params.disable_diagnostics &&
      ( // periodic display to stdout
        nstep_end % params.state_frequency == 0 ||
        // first two time steps
        tl.nstep <= tl.nstep0 + (nstep_end - tl.nstep) ));

  if (compute_diagnostics) {
    Diagnostics& diags = context.get<Diagnostics>();
    diags.run_diagnostics(true,2);
  }

  if ( ! independent_time_steps) {
    tl.update_tracers_levels(params.dt_tracer_factor);

    // Apply forcing.
#if defined(CAM) || defined(SCREAM)
    // CAM and SCREAM support only ftype 0 and 2.
    if (params.ftype == ForcingAlg::FORCING_0){
      apply_cam_forcing_tracers(dt_remap);
    }
    if (params.ftype == ForcingAlg::FORCING_2 && params.nsplit_iteration == 1 ){
      apply_cam_forcing_tracers(dt_remap*params.nsplit);
    }

    apply_cam_forcing_dynamics(dt_remap);
    
#else
    // standalone homme, support ftype0 and ftype2
    // ftype0  = ftype2 if dt_remap>=dt_tracer, but
    // ftype0 != ftype2 for dt_remap<dt_tracer
    if (params.ftype == ForcingAlg::FORCING_0 || 
        params.ftype == ForcingAlg::FORCING_2) {
      apply_cam_forcing(dt_remap);
    }
#endif

    if (compute_diagnostics) {
      Diagnostics& diags = context.get<Diagnostics>();
      diags.run_diagnostics(true,0);
    }

    // Loop over rsplit vertically lagrangian timesteps
    GPTLstart("tl-sc prim_step-loop");
    prim_step(dt,compute_diagnostics);
    for (int r=1; r<params.rsplit; ++r) {
      tl.update_dynamics_levels(UpdateType::LEAPFROG);
      prim_step(dt,false);
    }
    GPTLstop("tl-sc prim_step-loop");

    tl.update_tracers_levels(params.dt_tracer_factor);

    if (compute_diagnostics) {
      Diagnostics& diags = context.get<Diagnostics>();
      diags.run_diagnostics(false,3);
    }

    ////////////////////////////////////////////////////////////////////////
    // apply vertical remap
    // always for tracers
    // if rsplit>0:  also remap dynamics and compute reference level ps_v
    ////////////////////////////////////////////////////////////////////////
    GPTLstart("tl-sc vertical_remap");
    vertical_remap(dt_remap);
    GPTLstop("tl-sc vertical_remap");

    ////////////////////////////////////////////////////////////////////////
    // time step is complete.  update some diagnostic variables:
    // Q    (mixing ratio)
    ////////////////////////////////////////////////////////////////////////
    update_q(tl.np1_qdp,tl.np1);
  } else { // independent_time_steps
    prim_step_flexible(dt, compute_diagnostics);
  }

  if (compute_diagnostics) {
    Diagnostics& diags = context.get<Diagnostics>();
    diags.run_diagnostics(false,1);
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

void update_q (const int np1_qdp, const int np1)
{
  auto& context = Context::singleton();

  // Get simulation params
  SimulationParams& params = context.get<SimulationParams>();
  assert(params.params_set);

  // Get hybrid vertical coordinate
  HybridVCoord& hvcoord = context.get<HybridVCoord>();
  auto hyai_delta = hvcoord.hybrid_ai_delta;
  auto hybi_delta = hvcoord.hybrid_bi_delta;
  const Real ps0 = hvcoord.ps0;

  // Get ps_v from Elements
  Elements& elements = context.get<Elements>();
  auto ps_v = elements.m_state.m_ps_v;

  // Get the tracers concentration and mass from Tracers
  Tracers& tracers = context.get<Tracers>();
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
