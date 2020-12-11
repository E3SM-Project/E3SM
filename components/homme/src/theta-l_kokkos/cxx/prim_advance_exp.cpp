/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "CaarFunctor.hpp"
#include "DirkFunctor.hpp"
#include "Context.hpp"
#include "Diagnostics.hpp"
#include "Elements.hpp"
#include "HyperviscosityFunctor.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"

#include "profiling.hpp"

namespace Homme
{

// Declare all the timestepping schemes routines
void RK2_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w);
void imex_KG254_explicit_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w);
void u3_5stage_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w);
void imex_KG243_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w);
void imex_KG254_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w);
void imex_KG255_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w);

// -------------- IMPLEMENTATIONS -------------- //

void prim_advance_exp (TimeLevel& tl, const Real dt, const bool compute_diagnostics)
{
  GPTLstart("tl-ae prim_advance_exp");

#ifdef ARKODE
  Errors::runtime_abort("'ARKODE' support not yet available in C++ build.\n",
                         Errors::err_not_implemented);
#endif

  auto& context = Context::singleton();

  // Get simulation params
  SimulationParams& params = context.get<SimulationParams>();

  // Determine the tracers time level
  tl.update_tracers_levels(params.qsplit);

  // Set eta_ave_w
  Real eta_ave_w = 1.0/params.qsplit;

  // From f90 code: "this should not be needed, but in case physics update u without updating w b.c."
  if (!params.theta_hydrostatic_mode) {
    auto e = context.get<Elements>();
    auto w_i = e.m_state.m_w_i;
    auto v = e.m_state.m_v;
    auto gradphis = e.m_geometry.m_gradphis;
    auto n0 = tl.n0;
    constexpr auto LAST_LEV_P = ColInfo<NUM_INTERFACE_LEV>::LastPack;
    constexpr auto LAST_LEV   = ColInfo<NUM_PHYSICAL_LEV>::LastPack;
    constexpr auto LAST_INTERFACE_VEC_IDX = ColInfo<NUM_INTERFACE_LEV>::LastPackEnd;
    constexpr auto LAST_MIDPOINT_VEC_IDX  = ColInfo<NUM_PHYSICAL_LEV>::LastPackEnd;
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,NP*NP*e.m_geometry.num_elems()),
                         KOKKOS_LAMBDA(const int idx) {
      const int ie  = idx / (NP*NP);
      const int igp = (idx / NP) % NP;
      const int jgp = idx % NP;
      w_i(ie,n0,igp,jgp,LAST_LEV_P)[LAST_INTERFACE_VEC_IDX] = 
                      (v(ie,n0,0,igp,jgp,LAST_LEV)[LAST_MIDPOINT_VEC_IDX]*gradphis(ie,0,igp,jgp) +
                       v(ie,n0,1,igp,jgp,LAST_LEV)[LAST_MIDPOINT_VEC_IDX]*gradphis(ie,1,igp,jgp))/PhysicalConstants::g;
    });
  }

#ifndef CAM
  // if "prescribed wind" set dynamics explicitly and skip time-integration
  if (params.prescribed_wind) {
    Errors::runtime_abort("'prescribed wind' functionality not yet available in C++ build.\n",
                           Errors::err_not_implemented);
  }
#endif

  switch (params.time_step_type) {
    case TimeStepType::RK2:
      RK2_timestep (tl, dt, eta_ave_w);
      break;
    case TimeStepType::IMEX_KG254_EX:
      imex_KG254_explicit_timestep (tl, dt, eta_ave_w);
      break;
    case TimeStepType::ULLRICH_RK35:
      u3_5stage_timestep (tl, dt, eta_ave_w);
      break;
    case TimeStepType::IMEX_KG243:
      imex_KG243_timestep (tl, dt, eta_ave_w);
      break;
    case TimeStepType::IMEX_KG254:
      imex_KG254_timestep (tl, dt, eta_ave_w);
      break;
    case TimeStepType::IMEX_KG255:
      imex_KG255_timestep (tl, dt, eta_ave_w);
      break;
    default:
      {
        std::string msg = "[prim_advance_exp]:";
        msg += "missing some code for time step method ";
        msg += etoi(params.time_step_type);
        msg += ".\n";
        Errors::runtime_abort(msg,Errors::err_not_implemented);
      }
  }

  if (compute_diagnostics) {
    auto& diags = context.get<Diagnostics>();
    diags.run_diagnostics(false,4);
  }

  if (params.hypervis_order==2 && params.nu>0) {
    HyperviscosityFunctor& functor = context.get<HyperviscosityFunctor>();
    GPTLstart("tl-ae advance_hypervis_dp");
    functor.run(tl.np1,dt,eta_ave_w);
    GPTLstop("tl-ae advance_hypervis_dp");
  }

  if (params.dcmip16_mu>0) {
    Errors::runtime_abort("'dcmip16_mu>0' functionality not yet available in C++ build.\n",
                           Errors::err_not_implemented);
  }

  if (compute_diagnostics) {
    auto& diags = context.get<Diagnostics>();
    diags.run_diagnostics(false,5);
  }

  GPTLstop("tl-ae prim_advance_exp");
}

// Implementations of timestep schemes, in terms of CaarFunctor runs

void RK2_timestep(const TimeLevel& /* tl */,
                  const Real /* dt */,
                  const Real /* eta_ave_w */)
{
  // TODO
  assert(false);
}

void imex_KG254_explicit_timestep(const TimeLevel& /* tl */,
                                  const Real /* dt */,
                                  const Real /* eta_ave_w */)
{
  // TODO
  assert(false);
}

void u3_5stage_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w)
{
  GPTLstart("tl-ae U3-5stage_timestep");
  // Get elements structure
  Elements& elements = Context::singleton().get<Elements>();
  SimulationParams& params = Context::singleton().get<SimulationParams>();

  // Create the functor
  CaarFunctor& functor = Context::singleton().get<CaarFunctor>();

  const int nm1 = tl.nm1;
  const int n0  = tl.n0;
  const int np1 = tl.np1;
  const int qn0 = tl.n0_qdp;

  // ===================== RK STAGES ===================== //

  // Stage 1: u1 = u0 + dt/5 RHS(u0),          t_rhs = t
  functor.run(RKStageData(n0, n0, nm1, qn0, dt/5.0, eta_ave_w/4.0));

  // Stage 2: u2 = u0 + dt/5 RHS(u1),          t_rhs = t + dt/5
  functor.run(RKStageData(n0, nm1, np1, qn0, dt/5.0, 0.0));

  // Stage 3: u3 = u0 + dt/3 RHS(u2),          t_rhs = t + dt/5 + dt/5
  functor.run(RKStageData(n0, np1, np1, qn0, dt/3.0, 0.0));

  // Stage 4: u4 = u0 + 2dt/3 RHS(u3),         t_rhs = t + dt/5 + dt/5 + dt/3
  functor.run(RKStageData(n0, np1, np1, qn0, 2.0*dt/3.0, 0.0));

  // Compute (5u1-u0)/4 and store it in timelevel nm1
  {
    const auto v         = elements.m_state.m_v;
    const auto w         = elements.m_state.m_w_i;
    const auto vtheta_dp = elements.m_state.m_vtheta_dp;
    const auto phinh     = elements.m_state.m_phinh_i;
    const auto dp3d      = elements.m_state.m_dp3d;
    const auto hydrostatic_mode = params.theta_hydrostatic_mode;

    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, elements.num_elems()*NP*NP*NUM_LEV),
      KOKKOS_LAMBDA(const int it) {
        const int ie = it / (NP*NP*NUM_LEV);
        const int igp = (it / (NP*NUM_LEV)) % NP;
        const int jgp = (it / NUM_LEV) % NP;
        const int ilev = it % NUM_LEV;
        v(ie,nm1,0,igp,jgp,ilev) = (5.0*v(ie,nm1,0,igp,jgp,ilev)-v(ie,n0,0,igp,jgp,ilev))/4.0;
        v(ie,nm1,1,igp,jgp,ilev) = (5.0*v(ie,nm1,1,igp,jgp,ilev)-v(ie,n0,1,igp,jgp,ilev))/4.0;
        vtheta_dp(ie,nm1,igp,jgp,ilev) = (5.0*vtheta_dp(ie,nm1,igp,jgp,ilev)-vtheta_dp(ie,n0,igp,jgp,ilev))/4.0;
        dp3d(ie,nm1,igp,jgp,ilev) = (5.0*dp3d(ie,nm1,igp,jgp,ilev)-dp3d(ie,n0,igp,jgp,ilev))/4.0;
        if (!hydrostatic_mode) {
          w(ie,nm1,igp,jgp,ilev) = (5.0*w(ie,nm1,igp,jgp,ilev)-w(ie,n0,igp,jgp,ilev))/4.0;
          phinh(ie,nm1,igp,jgp,ilev) = (5.0*phinh(ie,nm1,igp,jgp,ilev)-phinh(ie,n0,igp,jgp,ilev))/4.0;
        }
    });
    // If NUM_LEV==NUM_LEV_P, the code above will take care also of the last interface
    if (NUM_LEV_P>NUM_LEV && !hydrostatic_mode) {
      const int LAST_INT = NUM_LEV_P-1;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, elements.num_elems()*NP*NP),
        KOKKOS_LAMBDA(const int it) {
           const int ie  =  it / (NP*NP);
           const int igp = (it / NP) % NP;
           const int jgp =  it % NP;
           w(ie,nm1,igp,jgp,LAST_INT) = (5.0*w(ie,nm1,igp,jgp,LAST_INT)-w(ie,n0,igp,jgp,LAST_INT))/4.0;
      });
    }
  }
  ExecSpace::impl_static_fence();

  // Stage 5: u5 = (5u1-u0)/4 + 3dt/4 RHS(u4), t_rhs = t + dt/5 + dt/5 + dt/3 + 2dt/3
  functor.run(RKStageData(nm1, np1, np1, qn0, 3.0*dt/4.0, 3.0*eta_ave_w/4.0));
  GPTLstop("tl-ae U3-5stage_timestep");
}

void imex_KG243_timestep(const TimeLevel& /* tl */,
                         const Real /* dt */,
                         const Real /* eta_ave_w */)
{
  // TODO
  assert(false);
}

void imex_KG254_timestep(const TimeLevel& /* tl */,
                         const Real /* dt */,
                         const Real /* eta_ave_w */)
{
  // TODO
  assert(false);
}


void imex_KG255_timestep(const TimeLevel& tl,
                         const Real dt_dyn,
                         const Real eta_ave_w)
{
  GPTLstart("IMEX_KG255");

  // The context
  const auto& c = Context::singleton();

  // Get elements, hvcoord, and functors
  auto& elements = c.get<Elements>();
  auto& hvcoord  = c.get<HybridVCoord>();
  auto& dirk     = c.get<DirkFunctor>();
  auto& caar     = c.get<CaarFunctor>();

  const int nm1 = tl.nm1;
  const int n0  = tl.n0;
  const int np1 = tl.np1;
  const int qn0 = tl.n0_qdp;

  // ===================== IMEX STAGES ===================== //

/////////////////////
//  Time level indices
//    caar: nm1, n0, np1
//    dirk: nm1, n0, np1
/////////////////////

  // Stage 1
  Real dt = dt_dyn/4.0;

  caar.run(RKStageData(n0, n0, nm1, qn0, dt, 0.0, 1.0, 0.0, 1.0));
  dirk.run(nm1, 0.0, n0, 0.0, nm1, dt, elements, hvcoord);

  // Stage 2
  dt = dt_dyn/6.0;

  caar.run(RKStageData(n0, nm1, np1, qn0, dt, 0.0, 1.0, 0.0, 1.0));
  dirk.run(nm1, 0.0, n0, 0.0, np1, dt, elements, hvcoord);

  // Stage 3
  dt = 3.0*dt_dyn/8.0;

  caar.run(RKStageData(n0, np1, np1, qn0, dt, 0.0, 1.0, 0.0, 1.0));
  dirk.run(nm1, 0.0, n0, 0.0, np1, dt, elements, hvcoord);

  // Stage 4
  dt = dt_dyn/2.0;

  caar.run(RKStageData(n0, np1, np1, qn0, dt, 0.0, 1.0, 0.0, 1.0));
  dirk.run(nm1, 0.0, n0, 0.0, np1, dt, elements, hvcoord);

  // Stage 5
  Real a1 = 0.24362;
  Real a2 = 0.34184;
  Real a3 = 1-(a1+a2);
  dt = dt_dyn;

  caar.run(RKStageData(n0, np1, np1, qn0, dt, eta_ave_w, 1.0, 0.0, 1.0));
  dirk.run(nm1, a2*dt, n0, a1*dt, np1, a3*dt, elements, hvcoord);

  GPTLstop("IMEX_KG255");
}

} // namespace Homme
