/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "CaarFunctor.hpp"
#include "CaarFunctorImpl.hpp"
#include "Context.hpp"
#include "SimulationParams.hpp"

#include "profiling.hpp"
#include "ErrorDefs.hpp"

#include <assert.h>
#include <type_traits>


namespace Homme {

CaarFunctor::CaarFunctor()
 : m_policy (Homme::get_default_team_policy<ExecSpace>(Context::singleton().get_elements().num_elems()))
{
  Elements&        elements   = Context::singleton().get_elements();
  Tracers&         tracers    = Context::singleton().get_tracers();
  Derivative&      derivative = Context::singleton().get_derivative();
  HybridVCoord&    hvcoord    = Context::singleton().get_hvcoord();
  SphereOperators& sphere_ops = Context::singleton().get_sphere_operators();
  const int        rsplit     = Context::singleton().get_simulation_params().rsplit;

  // Build functor impl
  m_caar_impl.reset(new CaarFunctorImpl(elements,tracers,derivative,hvcoord,sphere_ops,rsplit));
  m_caar_impl->m_sphere_ops.allocate_buffers(m_policy);
}

CaarFunctor::CaarFunctor(const Elements &elements, const Tracers &tracers,
                         const Derivative &derivative,
                         const HybridVCoord &hvcoord,
                         const SphereOperators &sphere_ops, 
                         const int rsplit)
    : m_policy(
          Homme::get_default_team_policy<ExecSpace>(elements.num_elems())) {
  // Build functor impl
  m_caar_impl.reset(new CaarFunctorImpl(elements, tracers, derivative, hvcoord,
                                        sphere_ops, rsplit));
  m_caar_impl->m_sphere_ops.allocate_buffers(m_policy);
}

CaarFunctor::~CaarFunctor ()
{
  // This empty destructor (where CaarFunctorImpl type is completely known)
  // is necessary for pimpl idiom to work with unique_ptr. The issue is the
  // deleter, which needs to know the size of the stored type, and which
  // would be called from the implicitly declared default destructor, which
  // would be in the header file, where CaarFunctorImpl type is incomplete.
}

void CaarFunctor::init_boundary_exchanges (const std::shared_ptr<BuffersManager>& bm_exchange) {
  assert (m_caar_impl);
  m_caar_impl->init_boundary_exchanges(bm_exchange);
}

void CaarFunctor::set_n0_qdp (const int n0_qdp)
{
  // Sanity check (should NEVER happen)
  assert (m_caar_impl);

  // Forward input to impl
  m_caar_impl->set_n0_qdp(n0_qdp);
}

void CaarFunctor::set_rk_stage_data (const int nm1, const int n0,   const int np1,
                                     const Real dt, const Real eta_ave_w,
                                     const bool compute_diagnostics)
{
  // Sanity check (should NEVER happen)
  assert (m_caar_impl);

  // Forward inputs to impl
  m_caar_impl->set_rk_stage_data(nm1,n0,np1,dt,eta_ave_w,compute_diagnostics);
}

void CaarFunctor::run ()
{
  // Sanity check (should NEVER happen)
  assert (m_caar_impl);

  // Run functor
  profiling_resume();
  GPTLstart("caar compute");
  Kokkos::parallel_for("caar loop pre-boundary exchange", m_policy, *m_caar_impl);
  ExecSpace::fence();
  GPTLstop("caar compute");
  profiling_pause();
}

void CaarFunctor::run (const int nm1, const int n0, const int np1,
                       const Real dt, const Real eta_ave_w,
                       const bool compute_diagnostics)
{
  // Sanity check (should NEVER happen)
  assert (m_caar_impl);

  // Forward inputs to impl
  m_caar_impl->set_rk_stage_data(nm1,n0,np1,dt,eta_ave_w,compute_diagnostics);

  // Run functor
  profiling_resume();
  GPTLstart("caar compute");
  Kokkos::parallel_for("caar loop pre-boundary exchange", m_policy, *m_caar_impl);
  ExecSpace::fence();
  GPTLstop("caar compute");
  start_timer("caar_bexchV");
  m_caar_impl->m_bes[np1]->exchange(m_caar_impl->m_elements.m_rspheremp);
  stop_timer("caar_bexchV");
  profiling_pause();
}

} // Namespace Homme
