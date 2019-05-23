/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "CaarFunctor.hpp"
#include "Context.hpp"
#include "Elements.hpp"
#include "HyperviscosityFunctor.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"

#include "profiling.hpp"

namespace Homme
{

void RK2_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics);
void imex_KG254_explicit_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics);
void u3_5stage_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics);
void imex_KGNO243_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics);
void imex_KGO254_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics);

// -------------- IMPLEMENTATIONS -------------- //

void prim_advance_exp (TimeLevel& tl, const Real dt, const bool compute_diagnostics)
{
  GPTLstart("tl-ae prim_advance_exp");

#ifdef ARKODE
  Errors::runtime_abort("'ARKODE' support not yet available in C++ build.\n",
                         Errors::err_not_implemented);
#endif

  // Get simulation params
  SimulationParams& params = Context::singleton().get<SimulationParams>();

  // Note: In the following, all the checks are superfluous, since we already check that
  //       the options are supported when we init the simulation params. However, this way
  //       we remind ourselves that in these cases there is some missing code to convert from Fortran

  // Determine the tracers time level
  tl.update_tracers_levels(params.qsplit);

  // Set eta_ave_w
  Real eta_ave_w = 1.0/params.qsplit;

  // From f90 code: "this should not be needed, but in case physics update u without updating w b.c."
  {
    auto e = Context::singleton().get<Elements>();
    auto w_i = e.m_state.m_w_i;
    auto v = e.m_state.m_v;
    auto gradphis = e.m_geometry.m_gradphis;
    auto n0 = tl.n0;
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,NP*NP*e.m_geometry.num_elems()),
                         KOKKOS_LAMBDA(const int idx) {
      const int ie  = idx / (NP*NP);
      const int igp = idx / NP;
      const int jgp = idx % NP;
      w_i(ie,n0,igp,jgp,NUM_PHYSICAL_LEV) = (v(ie,n0,0,igp,jgp,NUM_PHYSICAL_LEV-1)*gradphis(ie,0,igp,jgp) +
                                             v(ie,n0,1,igp,jgp,NUM_PHYSICAL_LEV-1)*gradphis(ie,1,igp,jgp))/PhysicalConstants::g;
    // do ie=nets,nete
    //    elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
    //                                        elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
    // enddo
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
    case 1:
      RK2_timestep (tl, dt, eta_ave_w, compute_diagnostics);
      break;
    case 4:
      imex_KG254_explicit_timestep (tl, dt, eta_ave_w, compute_diagnostics);
      break;
    case 5:
      u3_5stage_timestep (tl, dt, eta_ave_w, compute_diagnostics);
      break;
    case 6:
      imex_KGNO243_timestep (tl, dt, eta_ave_w, compute_diagnostics);
      break;
    case 7:
      imex_KGO254_timestep (tl, dt, eta_ave_w, compute_diagnostics);
      break;
    default:
      {
        std::string msg = "[prim_advance_exp]:";
        msg += "missing some code for time step method ";
        msg += params.time_step_type;
        msg += ".\n";
        Errors::runtime_abort(msg,Errors::err_not_implemented);
      }
  }

  if (params.hypervis_order==2 && params.nu>0) {
    HyperviscosityFunctor& functor = Context::singleton().get<HyperviscosityFunctor>();
    GPTLstart("tl-ae advance_hypervis_dp");
    functor.run(tl.np1,dt,eta_ave_w);
    GPTLstop("tl-ae advance_hypervis_dp");
  }

  if (params.dcmip16_mu>0) {
    Errors::runtime_abort("'dcmip16_mu>0' functionality not yet available in C++ build.\n",
                           Errors::err_not_implemented);
  }

  GPTLstop("tl-ae prim_advance_exp");
}

// Implementations of timestep schemes, in terms of CaarFunctor runs

void RK2_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics)
{
  // TODO
  assert(false);
}

void imex_KG254_explicit_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics)
{
  // TODO
  assert(false);
}

void u3_5stage_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics)
{
  GPTLstart("tl-ae U3-5stage_timestep");
  // Get elements structure
  Elements& elements = Context::singleton().get<Elements>();

  // Create the functor
  CaarFunctor& functor = Context::singleton().get<CaarFunctor>();

  // ===================== RK STAGES ===================== //

  // Stage 1: u1 = u0 + dt/5 RHS(u0),          t_rhs = t
  functor.run(RKStageData(tl.n0,tl.n0,tl.nm1,tl.n0_qdp,dt/5.0,eta_ave_w/4.0,compute_diagnostics));

  // Stage 2: u2 = u0 + dt/5 RHS(u1),          t_rhs = t + dt/5
  functor.run(RKStageData(tl.n0,tl.nm1,tl.np1,tl.n0_qdp,dt/5.0,0.0,false));

  // Stage 3: u3 = u0 + dt/3 RHS(u2),          t_rhs = t + dt/5 + dt/5
  functor.run(RKStageData(tl.n0,tl.np1,tl.np1,tl.n0_qdp,dt/3.0,0.0,false));

  // Stage 4: u4 = u0 + 2dt/3 RHS(u3),         t_rhs = t + dt/5 + dt/5 + dt/3
  functor.run(RKStageData(tl.n0,tl.np1,tl.np1,tl.n0_qdp,2.0*dt/3.0,0.0,false));

  // Compute (5u1-u0)/4 and store it in timelevel nm1
  {
    // const auto t    = elements.m_state.m_t;
    const auto v    = elements.m_state.m_v;
    const auto dp3d = elements.m_state.m_dp3d;
    const auto nm1  = tl.nm1;
    const auto n0   = tl.n0;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, elements.num_elems()*NP*NP*NUM_LEV),
      KOKKOS_LAMBDA(const int it) {
         const int ie = it / (NP*NP*NUM_LEV);
         const int igp = (it / (NP*NUM_LEV)) % NP;
         const int jgp = (it / NUM_LEV) % NP;
         const int ilev = it % NUM_LEV;
         // t(ie,nm1,igp,jgp,ilev) = (5.0*t(ie,nm1,igp,jgp,ilev)-t(ie,n0,igp,jgp,ilev))/4.0;
         v(ie,nm1,0,igp,jgp,ilev) = (5.0*v(ie,nm1,0,igp,jgp,ilev)-v(ie,n0,0,igp,jgp,ilev))/4.0;
         v(ie,nm1,1,igp,jgp,ilev) = (5.0*v(ie,nm1,1,igp,jgp,ilev)-v(ie,n0,1,igp,jgp,ilev))/4.0;
         dp3d(ie,nm1,igp,jgp,ilev) = (5.0*dp3d(ie,nm1,igp,jgp,ilev)-dp3d(ie,n0,igp,jgp,ilev))/4.0;
    });
  }
  ExecSpace::fence();

  // Stage 5: u5 = (5u1-u0)/4 + 3dt/4 RHS(u4), t_rhs = t + dt/5 + dt/5 + dt/3 + 2dt/3
  functor.run(RKStageData(tl.nm1,tl.np1,tl.np1,tl.n0_qdp,3.0*dt/4.0,3.0*eta_ave_w/4.0,false));
  GPTLstop("tl-ae U3-5stage_timestep");
}

void imex_KGNO243_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics)
{
  // TODO
  assert(false);
}

void imex_KGO254_timestep(const TimeLevel& tl, const Real dt, const Real eta_ave_w, const bool compute_diagnostics)
{
  // TODO
  assert(false);
}

} // namespace Homme
