/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"
#include "Elements.hpp"
#include "TimeLevel.hpp"
#include "SimulationParams.hpp"
#include "CamForcing.hpp"
#include "Diagnostics.hpp"
#include "ComposeTransport.hpp"
#include "profiling.hpp"

namespace Homme
{

void prim_advance_exp (TimeLevel& tl, const Real dt, const bool compute_diagnostics);
void prim_advec_tracers_remap (const Real);
void vertical_remap (const Real);

static void set_tracer_transport_derived_values (
  const SimulationParams& params, const Elements& elements, const TimeLevel& tl)
{
  // ===============
  // initialize mean flux accumulation variables and save some variables at n0
  // for use by advection
  // ===============
  GPTLstart("tl-s deep_copy+derived_dp");
  {
    const auto eta_dot_dpdn = elements.m_derived.m_eta_dot_dpdn;
    const auto derived_vn0 = elements.m_derived.m_vn0;
    const auto omega_p = elements.m_derived.m_omega_p;
    const auto derived_dpdiss_ave = elements.m_derived.m_dpdiss_ave;
    const auto derived_dpdiss_biharmonic = elements.m_derived.m_dpdiss_biharmonic;
    const auto derived_dp = elements.m_derived.m_dp;
    const auto dp3d = elements.m_state.m_dp3d;
    const auto vstar = elements.m_derived.m_vstar;
    const auto v = elements.m_state.m_v;
    const auto n0 = tl.n0;
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace> (0,elements.num_elems()*NP*NP*NUM_LEV),
                         KOKKOS_LAMBDA(const int idx) {
      const int ie   = ((idx / NUM_LEV) / NP) / NP;
      const int igp  = ((idx / NUM_LEV) / NP) % NP;
      const int jgp  =  (idx / NUM_LEV) % NP;
      const int ilev =   idx % NUM_LEV;
      eta_dot_dpdn(ie,igp,jgp,ilev) = 0;
      derived_vn0(ie,0,igp,jgp,ilev) = 0;
      derived_vn0(ie,1,igp,jgp,ilev) = 0;
      omega_p(ie,igp,jgp,ilev) = 0;
      if (params.nu_p>0) {
        derived_dpdiss_ave(ie,igp,jgp,ilev) = 0;
        derived_dpdiss_biharmonic(ie,igp,jgp,ilev) = 0;
      }
      derived_dp(ie,igp,jgp,ilev) = dp3d(ie,n0,igp,jgp,ilev);
      if (params.transport_alg > 0) {
        vstar(ie,0,igp,jgp,ilev) = v(ie,n0,0,igp,jgp,ilev);
        vstar(ie,1,igp,jgp,ilev) = v(ie,n0,1,igp,jgp,ilev);
      }
    });
  }
  Kokkos::fence();
  GPTLstop("tl-s deep_copy+derived_dp");  
}

void prim_step (const Real dt, const bool compute_diagnostics)
{
  GPTLstart("tl-s prim_step");
  // Get control and simulation params
  SimulationParams& params = Context::singleton().get<SimulationParams>();
  assert(params.params_set);

  // Get the elements structure
  Elements& elements = Context::singleton().get<Elements>();

  // Get the time level info
  TimeLevel& tl = Context::singleton().get<TimeLevel>();

  set_tracer_transport_derived_values(params, elements, tl);

  // ===============
  // Dynamical Step
  // ===============
  GPTLstart("tl-s prim_advance_exp-loop");
  prim_advance_exp(tl,dt,compute_diagnostics);
  tl.tevolve += dt;
  for (int n=1; n<params.dt_tracer_factor; ++n) {
    tl.update_dynamics_levels(UpdateType::LEAPFROG);
    prim_advance_exp(tl,dt,false);
    tl.tevolve += dt;
  }
  GPTLstop("tl-s prim_advance_exp-loop");

  // ===============
  // Tracer Advection.
  // in addition, this routine will apply the DSS to:
  //        derived%eta_dot_dpdn    =  mean vertical velocity (used for remap below)
  //        derived%omega           =
  // Tracers are always vertically lagrangian.
  // For rsplit=0:
  //   if tracer scheme needs v on lagrangian levels it has to vertically interpolate
  //   if tracer scheme needs dp3d, it needs to derive it from ps_v
  // ===============
  // Advect tracers if their count is > 0.
  // not be advected.  This will be cleaned up when the physgrid is merged into CAM trunk
  // Currently advecting all species
  GPTLstart("tl-s prim_advec_tracers_remap");
  if (params.qsize>0) {
    prim_advec_tracers_remap(dt*params.dt_tracer_factor);
  }
  GPTLstop("tl-s prim_advec_tracers_remap");
  GPTLstop("tl-s prim_step");
}

void prim_step_flexible (const Real dt, const bool compute_diagnostics) {
#ifdef MODEL_THETA_L
  GPTLstart("tl-s prim_step_flexible");
  const auto& context = Context::singleton();
  const SimulationParams& params = context.get<SimulationParams>();
  assert(params.params_set);
  Elements& elements = context.get<Elements>();
  TimeLevel& tl = context.get<TimeLevel>();

  const auto dt_remap = params.dt_remap_factor == 0 ? dt : dt*params.dt_remap_factor;

  tl.update_tracers_levels(params.dt_tracer_factor);
  
  const bool forcing_0or2 = (params.ftype == ForcingAlg::FORCING_0 ||
                             params.ftype == ForcingAlg::FORCING_2);
  bool apply_forcing;

  const auto dt_q =        dt * params.dt_tracer_factor;
  // In standalone HOMME, nsplit=1 always.
  const auto dt_q_nsplit = dt * params.dt_tracer_factor * params.nsplit;
  
  // Decide on tracer forcing.
#if defined(CAM) || defined(SCREAM)
  // CAM + xx supports only ftype 0 and 2.
  apply_forcing = (params.ftype == ForcingAlg::FORCING_0) ||
                  (params.ftype == ForcingAlg::FORCING_2 && params.nsplit_iteration == 1);
#else
  apply_forcing = forcing_0or2;
#endif

  if (apply_forcing) {
    if (params.ftype == ForcingAlg::FORCING_0) apply_cam_forcing_tracers(dt_q);
    if (params.ftype == ForcingAlg::FORCING_2) apply_cam_forcing_tracers(dt_q_nsplit);
  }

  set_tracer_transport_derived_values(params, elements, tl);

  for (int n = 0; n < params.dt_tracer_factor; ++n) {
    const bool compute_diagnostics_it = compute_diagnostics && n == 0;
    if (n > 0) tl.update_dynamics_levels(UpdateType::LEAPFROG);

    if (forcing_0or2) {
      // Apply dynamics forcing over the dynamics (vertically Eulerian) or
      // vertical remap time step if we're at reference levels.
      apply_forcing = (params.dt_remap_factor == 0 ||
                       n % params.dt_remap_factor == 0);
      if (apply_forcing) {
        apply_cam_forcing_dynamics(dt_remap);
        if (compute_diagnostics_it)
          context.get<Diagnostics>().run_diagnostics(true, 0);
      }
    }

    prim_advance_exp(tl, dt, compute_diagnostics);
    tl.tevolve += dt;

    if (params.dt_remap_factor == 0) {
      // Since dt_remap == 0, the only part of vertical_remap that is active is
      // the updates to ps_v(:,:,np1) and dp3d(:,:,:,np1).
      vertical_remap(dt_remap);
    } else if ((n+1) % params.dt_remap_factor == 0) {
      if (compute_diagnostics)
        context.get<Diagnostics>().run_diagnostics(false, 3);
      if (params.prescribed_wind) {
        // Prescribed winds are evaluated on reference levels, not floating
        // levels, so don't remap, just update dp3d.
        Errors::runtime_abort("prim_step_flexible: need to impl prescribed_wind\n",
                              Errors::err_not_implemented);
      } else {
        // Remap dynamics variables but not tracers.
        GPTLstart("tl-sc vertical_remap");
        vertical_remap(dt_remap);
        GPTLstop("tl-sc vertical_remap");
      }
    }
  }

  if (params.qsize > 0)
    prim_advec_tracers_remap(dt_q);

  if (params.dt_remap_factor == 0 && compute_diagnostics)
    context.get<Diagnostics>().run_diagnostics(false, 3);

  // Remap tracers.
#ifdef HOMME_ENABLE_COMPOSE	  
  if (params.qsize > 0)
    Context::singleton().get<ComposeTransport>().remap_q(tl);
#endif

  GPTLstop("tl-s prim_step_flexible");
#else
  Errors::runtime_abort("prim_step_flexible not supported in non-theta-l builds.");
#endif
}

} // namespace Homme
