/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_EULER_STEP_FUNCTOR_IMPL_HPP
#define HOMMEXX_EULER_STEP_FUNCTOR_IMPL_HPP

#include "Context.hpp"
#include "Derivative.hpp"
#include "Elements.hpp"
#include "ErrorDefs.hpp"
#include "EulerStepFunctor.hpp"
#include "HommexxEnums.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "SphereOperators.hpp"
#include "Tracers.hpp"
#include "profiling.hpp"
#include "mpi/BoundaryExchange.hpp"
#include "mpi/MpiContext.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/VectorUtils.hpp"
#include "vector/vector_pragmas.hpp"


namespace Kokkos {
struct Real2 {
  Homme::Real v[2];
  KOKKOS_FORCEINLINE_FUNCTION Real2 () { v[0] = v[1] = 0; }
  KOKKOS_FORCEINLINE_FUNCTION Real2 operator+= (const Real2& o) const {
    Real2 r;
    r.v[0] = v[0] + o.v[0];
    r.v[1] = v[1] + o.v[1];
    return r;
  }
};

template<> struct reduction_identity<Real2> {
  KOKKOS_FORCEINLINE_FUNCTION static Real2 sum() { return Real2(); }
};
}

namespace Homme {

// On older machines, low memory b/w is the most important performance-influence
// characteristic. On these machines, a tracer will be processed in serial. Take
// advantage of this to optimize memory access.
template <typename ExecSpace>
struct SerialLimiter {
  template <int limiter_option, typename ArrayGll, typename ArrayGllLvl, typename Array2Lvl,
            typename Array2GllLvl>
  KOKKOS_INLINE_FUNCTION static void
  run(const ArrayGll& sphweights, const ArrayGllLvl& idpmass,
      const Array2Lvl& iqlim, const ArrayGllLvl& iptens,
      const Array2GllLvl& irwrk);
};
// GPU doesn't have a serial impl.
#if defined KOKKOS_ENABLE_CUDA
template <>
struct SerialLimiter<Kokkos::Cuda> {
  template <int limiter_option, typename ArrayGll, typename ArrayGllLvl, typename Array2Lvl,
            typename Array2GllLvl>
  KOKKOS_INLINE_FUNCTION static void
  run (const ArrayGll& sphweights, const ArrayGllLvl& idpmass,
       const Array2Lvl& iqlim, const ArrayGllLvl& iptens,
       const Array2GllLvl& irwrk) {
    Kokkos::abort("SerialLimiter::run: Should not be called on GPU.");
  }
};
#endif

class EulerStepFunctorImpl {
  struct EulerStepData {
    EulerStepData ()
      : qsize(-1), limiter_option(0), nu_p(0), nu_q(0), consthv(1)
    {}

    int   qsize;
    int   limiter_option;
    Real  rhs_viss;
    Real  rhs_multiplier;

    Real  nu_p;
    Real  nu_q;

    Real  dt;
    int   np1_qdp;
    int   n0_qdp;

    DSSOption   DSSopt;

    bool consthv;
  };

  const Elements      m_elements;
  const Tracers       m_tracers;
  const Derivative    m_deriv;
  const HybridVCoord  m_hvcoord;
  EulerStepData       m_data;
  SphereOperators     m_sphere_ops;

  bool                m_kernel_will_run_limiters;

  std::shared_ptr<BoundaryExchange> m_mm_be, m_mmqb_be;
  Kokkos::Array<std::shared_ptr<BoundaryExchange>, 3*Q_NUM_TIME_LEVELS> m_bes;

  enum { m_mem_per_team = 2 * NP * NP * sizeof(Real) };

public:

  EulerStepFunctorImpl ()
   : m_elements   (Context::singleton().get_elements())
   , m_tracers    (Context::singleton().get_tracers())
   , m_deriv      (Context::singleton().get_derivative())
   , m_hvcoord    (Context::singleton().get_hvcoord())
   , m_sphere_ops (Context::singleton().get_sphere_operators())
  {
  }

  void reset (const SimulationParams& params) {
    m_data.rhs_viss = 0.0;
    m_data.qsize = params.qsize;
    m_data.limiter_option = params.limiter_option;
    m_data.nu_p = params.nu_p;
    m_data.nu_q = params.nu_q;
    m_data.consthv = (params.hypervis_scaling == 0);

    if (m_data.limiter_option == 4) {
      std::string msg = "[EulerStepFunctorImpl::reset]:";
      msg += "limiter_option=4 is not yet supported in C++. ";
      msg += "The program should have errored out earlier though. Plese, investigate.";
      Errors::runtime_abort(msg,Errors::err_not_implemented);
    }

    // This will fit the needs of all calls to sphere operators.
    m_sphere_ops.allocate_buffers(Homme::get_default_team_policy<ExecSpace>(m_elements.num_elems()*m_data.qsize));
  }

  void init_boundary_exchanges () {
    assert(m_data.qsize >= 0); // after reset() called

    auto bm_exchange = MpiContext::singleton().get_buffers_manager(MPI_EXCHANGE);
    for (int np1_qdp = 0, k = 0; np1_qdp < Q_NUM_TIME_LEVELS; ++np1_qdp) {
      for (int dssi = 0; dssi < 3; ++dssi, ++k) {
        m_bes[k] = std::make_shared<BoundaryExchange>();
        BoundaryExchange& be = *m_bes[k];
        be.set_buffers_manager(bm_exchange);
        be.set_num_fields(0, 0, m_data.qsize+1);
        be.register_field(m_tracers.qdp, np1_qdp, m_data.qsize, 0);
        be.register_field(dssi == 0 ? m_elements.m_eta_dot_dpdn :
                          dssi == 1 ? m_elements.m_omega_p :
                          m_elements.m_derived_divdp_proj);
        be.registration_completed();
      }
    }

    {
      m_mmqb_be = std::make_shared<BoundaryExchange>();
      m_mmqb_be->set_buffers_manager(bm_exchange);
      m_mmqb_be->set_num_fields(0, 0, m_data.qsize);
      m_mmqb_be->register_field(m_tracers.qtens_biharmonic, m_data.qsize, 0);
      m_mmqb_be->registration_completed();
    }

    {
      auto bm_exchange_minmax = MpiContext::singleton().get_buffers_manager(MPI_EXCHANGE_MIN_MAX);
      m_mm_be = std::make_shared<BoundaryExchange>();
      BoundaryExchange& be = *m_mm_be;
      be.set_buffers_manager(bm_exchange_minmax);
      be.set_num_fields(m_data.qsize, 0, 0);
      be.register_min_max_fields(m_tracers.qlim, m_data.qsize, 0);
      be.registration_completed();
    }
  }

  static size_t limiter_team_shmem_size (const int team_size) {
    return Memory<ExecSpace>::on_gpu ?
      (team_size * m_mem_per_team) :
      0;
  }

  size_t team_shmem_size (const int team_size) const {
    return m_kernel_will_run_limiters ? limiter_team_shmem_size(team_size) : 0;
  }

  struct BIHPreNup {};
  struct BIHPreNoNup {};
  struct BIHPostConstHV {};
  struct BIHPostTensorHV {};

  /*
    ! get new min/max values, and also compute biharmonic mixing term

    ! two scalings depending on nu_p:
    ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscosity only to q)
    ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
   */
  void compute_biharmonic_pre() {
    profiling_resume();
    assert(m_data.rhs_multiplier == 2.0);
    m_data.rhs_viss = 3.0;

    if(m_data.nu_p > 0){
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, BIHPreNup>(
                             m_elements.num_elems() * m_data.qsize),
                         *this);
    }else{
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, BIHPreNoNup>(
                             m_elements.num_elems() * m_data.qsize),
                         *this);

    }

    ExecSpace::fence();
    profiling_pause();
  }

  void compute_biharmonic_post() {
    profiling_resume();
    assert(m_data.rhs_multiplier == 2.0);

    if(m_data.consthv){
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, BIHPostConstHV>(
                             m_elements.num_elems() * m_data.qsize),
                         *this);
    }else{
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, BIHPostTensorHV>(
                             m_elements.num_elems() * m_data.qsize),
                         *this);
    }
    ExecSpace::fence();
    profiling_pause();
  }

//case when nu_p > 0
  KOKKOS_INLINE_FUNCTION
  void operator() (const BIHPreNup&, const TeamMember& team) const {
    KernelVariables kv(team,m_data.qsize);
    const auto qtens_biharmonic = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    dpdiss_adjustment(kv, team);
    m_sphere_ops.laplace_simple(kv, qtens_biharmonic, qtens_biharmonic);
  }

//case when nu_p == 0
  KOKKOS_INLINE_FUNCTION
  void operator() (const BIHPreNoNup&, const TeamMember& team) const {
    KernelVariables kv(team,m_data.qsize);
    const auto qtens_biharmonic = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    m_sphere_ops.laplace_simple(kv, qtens_biharmonic, qtens_biharmonic);
  }

  KOKKOS_INLINE_FUNCTION
  void dpdiss_adjustment (KernelVariables & kv, const TeamMember& team) const {

    const auto qtens_biharmonic = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
      const auto dpdiss_ave = Homme::subview(m_elements.m_derived_dpdiss_ave, kv.ie);
      Kokkos::parallel_for (
        Kokkos::TeamThreadRange(team, NP*NP),
        [&] (const int loop_idx) {
          const int i = loop_idx / NP;
          const int j = loop_idx % NP;
          Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(team, NUM_LEV),
            [&] (const int& k) {
              qtens_biharmonic(i,j,k) = qtens_biharmonic(i,j,k) * dpdiss_ave(i,j,k) / m_hvcoord.dp0(k);
            });
        });
      team.team_barrier();
  }



  KOKKOS_INLINE_FUNCTION
  void operator() (const BIHPostConstHV&, const TeamMember& team) const {
    KernelVariables kv(team,m_data.qsize);
    const auto qtens_biharmonic = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    team.team_barrier();
    m_sphere_ops.laplace_simple(kv, qtens_biharmonic, qtens_biharmonic);
    // laplace_simple provides the barrier.
   rhsviss_adjustment(kv, team);
  }//end of BIHPostConstHV ()

  KOKKOS_INLINE_FUNCTION
  void operator() (const BIHPostTensorHV&, const TeamMember& team) const {
    KernelVariables kv(team,m_data.qsize);
    const auto& e = m_elements;
    const auto qtens_biharmonic = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    const auto tensor = Homme::subview(e.m_tensorvisc, kv.ie);
    team.team_barrier();
    m_sphere_ops.laplace_tensor(kv, tensor, qtens_biharmonic, qtens_biharmonic);
    // divergence_sphere_wk provides the barrier.
    rhsviss_adjustment(kv, team);
    }//end of BIHPostTensorHV ()

  KOKKOS_INLINE_FUNCTION
  void rhsviss_adjustment (KernelVariables & kv, const TeamMember & team) const {
    const auto qtens_biharmonic = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    const auto f = -m_data.rhs_viss * m_data.dt * m_data.nu_q;
    const auto spheremp = Homme::subview(m_elements.m_spheremp, kv.ie);
      Kokkos::parallel_for (
        Kokkos::TeamThreadRange(team, NP*NP),
        [&] (const int loop_idx) {
          const int i = loop_idx / NP;
          const int j = loop_idx % NP;
          Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(team, NUM_LEV),
            [&] (const int& k) {
              qtens_biharmonic(i,j,k) = (f * m_hvcoord.dp0(k) * qtens_biharmonic(i,j,k) /
                                         spheremp(i,j));
            });
        });
    }



  struct AALSetupPhase {};
  struct AALTracerPhase {};

  void advect_and_limit() {
    profiling_resume();
    Kokkos::parallel_for(
      Homme::get_default_team_policy<ExecSpace, AALSetupPhase>(
        m_elements.num_elems()),
      *this);
    ExecSpace::fence();
    m_kernel_will_run_limiters = true;
    Kokkos::parallel_for(
      Homme::get_default_team_policy<ExecSpace, AALTracerPhase>(
        m_elements.num_elems() * m_data.qsize),
      *this);
    ExecSpace::fence();
    m_kernel_will_run_limiters = false;
    profiling_pause();
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const AALSetupPhase&, const TeamMember& team) const {
    KernelVariables kv(team);
    run_setup_phase(kv);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const AALTracerPhase&, const TeamMember& team) const {
    KernelVariables kv(team, m_data.qsize);
    run_tracer_phase(kv);
  }

  struct PrecomputeDivDp {};

  void precompute_divdp() {
    assert(m_data.qsize >= 0); // reset() already called
    profiling_resume();

    Kokkos::parallel_for(
        Homme::get_default_team_policy<ExecSpace, PrecomputeDivDp>(
            m_elements.num_elems()),
        *this);

    ExecSpace::fence();
    profiling_pause();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const PrecomputeDivDp &, const TeamMember &team) const {
    KernelVariables kv(team);
    m_sphere_ops.divergence_sphere(kv,
                      Homme::subview(m_elements.m_derived_vn0, kv.ie),
                      Homme::subview(m_elements.m_derived_divdp, kv.ie));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, NUM_LEV),
                           [&](const int ilev) {
        m_elements.m_derived_divdp_proj(kv.ie, igp, jgp, ilev) =
            m_elements.m_derived_divdp(kv.ie, igp, jgp, ilev);
      });
    });
  }

  void qdp_time_avg (const int n0_qdp, const int np1_qdp) {
    const int qsize = m_data.qsize;
    const auto qdp = m_tracers.qdp;
    const Real rkstage = 3.0;
    Kokkos::parallel_for(
      Homme::get_default_team_policy<ExecSpace>(m_elements.num_elems()*m_data.qsize),
      KOKKOS_LAMBDA(const TeamMember& team) {
        KernelVariables kv(team, qsize);
        const auto qdp_n0 = Homme::subview(qdp, kv.ie, n0_qdp, kv.iq);
        const auto qdp_np1 = Homme::subview(qdp, kv.ie, np1_qdp, kv.iq);
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(kv.team, NP*NP),
          [&] (const int& idx) {
            const int i = idx / NP;
            const int j = idx % NP;
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
              [&] (const int& ilev) {
               qdp_np1(i,j,ilev) =
                 (qdp_n0(i,j,ilev) + (rkstage-1)*qdp_np1(i,j,ilev)) /
                 rkstage;
            });
          });
      });
  }

  // TODO make GPUable.
  void compute_dp () {
    const auto& e = m_elements;
    const auto& c = m_data;
    const auto dp = e.m_derived_dp;
    const auto divdp_proj = e.m_derived_divdp_proj;
    const auto rhsmdt = c.rhs_multiplier * c.dt;
    const auto buf = e.buffers.pressure;
    Kokkos::parallel_for(
      Homme::get_default_team_policy<ExecSpace>(m_elements.num_elems()),
      KOKKOS_LAMBDA (const TeamMember& team) {
        KernelVariables kv(team);
        Kokkos::parallel_for (
          Kokkos::TeamThreadRange(kv.team, NP*NP),
          [&] (const int loop_idx) {
            const int i = loop_idx / NP;
            const int j = loop_idx % NP;
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
              [&] (const int& k) {
                //! derived variable divdp_proj() (DSS'd version of divdp) will only
                //! be correct on 2nd and 3rd stage but that's ok because
                //! rhs_multiplier=0 on the first stage:
                // Store this in unused buffer.
                buf(kv.ie,i,j,k) =
                  dp(kv.ie,i,j,k) - rhsmdt * divdp_proj(kv.ie,i,j,k);
              });
          });
      });
  }

  void compute_qmin_qmax() {
    // Temporaries, due to issues capturing *this on device
    const int qsize = m_data.qsize;
    const Real rhs_multiplier = m_data.rhs_multiplier;
    const int n0_qdp = m_data.n0_qdp;
    const auto qdp = m_tracers.qdp;
    const auto dp = m_elements.buffers.pressure;
    const auto qtens_biharmonic = m_tracers.qtens_biharmonic;
    const auto qlim = m_tracers.qlim;
    const auto num_parallel_iterations = m_elements.num_elems() * m_data.qsize;
    ThreadPreferences tp;
    tp.max_threads_usable = NUM_LEV;
    tp.max_vectors_usable = 1;
    const auto tv =
      DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
        num_parallel_iterations, tp);
    Kokkos::TeamPolicy<ExecSpace> policy(num_parallel_iterations,
                                         tv.first, tv.second);
    Kokkos::parallel_for(
      policy,
      KOKKOS_LAMBDA (const TeamMember& team) {
        KernelVariables kv(team, qsize);
        const auto dp_t = Homme::subview(dp, kv.ie);
        const auto qdp_t = Homme::subview(qdp, kv.ie, n0_qdp, kv.iq);
        const auto qtens_biharmonic_t = Homme::subview(qtens_biharmonic, kv.ie, kv.iq);
        const auto qlim_t = Homme::subview(qlim, kv.ie, kv.iq);
        if (rhs_multiplier != 1.0) {
          Kokkos::parallel_for(
            Kokkos::TeamThreadRange(kv.team, NUM_LEV),
            [&] (const int& k) {
              const auto v = qdp_t(0,0,k) / dp_t(0,0,k);
              qtens_biharmonic_t(0,0,k) = v;
              qlim_t(0,k) = v;
              qlim_t(1,k) = v;
            });
        }
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j) {
            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(kv.team, NUM_LEV),
              [&] (const int& k) {
                const auto v = qdp_t(i,j,k) / dp_t(i,j,k);
                qtens_biharmonic_t(i,j,k) = v;
                qlim_t(0,k) = min(qlim_t(0,k), v);
                qlim_t(1,k) = max(qlim_t(1,k), v);
              });
          }
      });
    ExecSpace::fence();
  }

  void neighbor_minmax_start() {
    assert(m_mm_be->is_registration_completed());
    m_mm_be->pack_and_send_min_max();
  }

  void neighbor_minmax_finish() {
    m_mm_be->recv_and_unpack_min_max();
  }

  void minmax_and_biharmonic() {
    neighbor_minmax_start();
    compute_biharmonic_pre();
    m_mmqb_be->exchange(m_elements.m_rspheremp);
    compute_biharmonic_post();
    neighbor_minmax_finish();
  }

  void neighbor_minmax() {
    assert(m_mm_be->is_registration_completed());
    m_mm_be->exchange_min_max();
  }

  void exchange_qdp_dss_var () {
    const int idx = 3*m_data.np1_qdp + static_cast<int>(m_data.DSSopt);
    m_bes[idx]->exchange(m_elements.m_rspheremp);
  }

  void euler_step(const int np1_qdp, const int n0_qdp, const Real dt,
                  const Real rhs_multiplier, const DSSOption DSSopt) {

    m_data.n0_qdp         = n0_qdp;
    m_data.np1_qdp        = np1_qdp;
    m_data.dt             = dt;
    m_data.rhs_multiplier = rhs_multiplier;
    m_data.DSSopt         = DSSopt;

    if (EulerStepFunctor::is_quasi_monotone(m_data.limiter_option)) {
      // when running lim8, we also need to limit the biharmonic, so that term
      // needs to be included in each euler step.  three possible algorithms
      // here:
      // most expensive:
      //   compute biharmonic (which also computes qmin/qmax) during all 3
      //   stages be sure to set rhs_viss=1 cost:  3 biharmonic steps with 3 DSS

      // cheapest:
      //   compute biharmonic (which also computes qmin/qmax) only on first
      //   stage be sure to set rhs_viss=3 reuse qmin/qmax for all following
      //   stages (but update based on local qmin/qmax) cost:  1 biharmonic
      //   steps with 1 DSS main concern: viscosity

      // compromise:
      //   compute biharmonic (which also computes qmin/qmax) only on last stage
      //   be sure to set rhs_viss=3
      //   compute qmin/qmax directly on first stage
      //   reuse qmin/qmax for 2nd stage stage (but update based on local
      //   qmin/qmax) cost:  1 biharmonic steps, 2 DSS

      //  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation),
      //        we should apply dissipation to Q (not Qdp) to preserve Q=1
      //        i.e.  laplace(Qdp) ~  dp0 laplace(Q)
      //        for nu_p=nu_q>0, we need to apply dissipation to Q *diffusion_dp

      // initialize dp, and compute Q from Qdp(and store Q in Qtens_biharmonic)

      compute_dp();
      compute_qmin_qmax();
      if (m_data.rhs_multiplier == 0.0) {
        neighbor_minmax();
      } else if (m_data.rhs_multiplier == 2.0) {
        minmax_and_biharmonic();
      }
    }
    advect_and_limit();
    exchange_qdp_dss_var();
  }

private:

  KOKKOS_INLINE_FUNCTION
  void run_setup_phase (const KernelVariables& kv) const {
    compute_2d_advection_step(kv);
  }

  KOKKOS_INLINE_FUNCTION
  void run_tracer_phase (const KernelVariables& kv) const {
    compute_qtens(kv);
    kv.team_barrier();
    if (m_data.limiter_option == 8) {
      limiter_optim_iter_full(kv);
      kv.team_barrier();
    } else if (m_data.limiter_option == 9) {
      limiter_clip_and_sum(kv);
      kv.team_barrier();
    }
    apply_spheremp(kv);
  }

  KOKKOS_INLINE_FUNCTION
  void compute_2d_advection_step (const KernelVariables& kv) const {
    const auto& c = m_data;
    const auto& e = m_elements;
    const bool lim_quasi_monotone
      = EulerStepFunctor::is_quasi_monotone(c.limiter_option);
    const bool add_ps_diss = c.nu_p > 0 && c.rhs_viss != 0.0;
    const Real diss_fac = add_ps_diss ? -c.rhs_viss * c.dt * c.nu_q : 0;
    const auto& f_dss = (c.DSSopt == DSSOption::ETA ?
                         e.m_eta_dot_dpdn :
                         c.DSSopt == DSSOption::OMEGA ?
                         e.m_omega_p :
                         e.m_derived_divdp_proj);
    Kokkos::parallel_for (
      Kokkos::TeamThreadRange(kv.team, NP*NP),
      [&] (const int loop_idx) {
        const int i = loop_idx / NP;
        const int j = loop_idx % NP;
        Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
          [&] (const int& k) {
            const auto dp = e.buffers.pressure(kv.ie,i,j,k);
            e.buffers.vstar(kv.ie,0,i,j,k) = e.m_derived_vn0(kv.ie,0,i,j,k) / dp;
            e.buffers.vstar(kv.ie,1,i,j,k) = e.m_derived_vn0(kv.ie,1,i,j,k) / dp;
            if (lim_quasi_monotone) {
              //! Note that the term dpdissk is independent of Q
              //! UN-DSS'ed dp at timelevel n0+1:
              e.buffers.dpdissk(kv.ie,i,j,k) = dp - c.dt * e.m_derived_divdp(kv.ie,i,j,k);
              if (add_ps_diss) {
                //! add contribution from UN-DSS'ed PS dissipation
                //!          dpdiss(:,:) = ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) *
                //!          elem(ie)%derived%psdiss_biharmonic(:,:)
                e.buffers.dpdissk(kv.ie,i,j,k) += diss_fac *
                  e.m_derived_dpdiss_biharmonic(kv.ie,i,j,k) / e.m_spheremp(kv.ie,i,j);
              }
            }
            //! also DSS extra field
            //! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
            //! all zero so we only have to DSS 1:nlev
            f_dss(kv.ie,i,j,k) *= e.m_spheremp(kv.ie,i,j);
          });
      });
  }

  KOKKOS_INLINE_FUNCTION
  void compute_qtens (const KernelVariables& kv) const {
    m_sphere_ops.divergence_sphere_update(
      kv, -m_data.dt, m_data.rhs_viss != 0.0,
      Homme::subview(m_elements.buffers.vstar, kv.ie),
      Homme::subview(m_tracers.qdp, kv.ie, m_data.n0_qdp, kv.iq),
      // On input, qtens_biharmonic if add_hyperviscosity, undefined
      // if not; on output, qtens.
      Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq));
  }

  KOKKOS_INLINE_FUNCTION
  void limiter_optim_iter_full (const KernelVariables& kv) const {
    const auto sphweights = Homme::subview(m_elements.m_spheremp, kv.ie);
    const auto dpmass = Homme::subview(m_elements.buffers.dpdissk, kv.ie);
    const auto ptens = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    const auto qlim = Homme::subview(m_tracers.qlim, kv.ie, kv.iq);
    if ( ! OnGpu<ExecSpace>::value && kv.team.team_size() == 1)
      SerialLimiter<ExecSpace>::run<8>(
        sphweights, dpmass, qlim, ptens,
        Homme::subview(m_sphere_ops.vector_buf_ml, kv.team_idx, 0));
    else
      limiter_optim_iter_full(kv.team, sphweights, dpmass, qlim, ptens);
  }

  KOKKOS_INLINE_FUNCTION
  void limiter_clip_and_sum (const KernelVariables& kv) const {
    const auto sphweights = Homme::subview(m_elements.m_spheremp, kv.ie);
    const auto dpmass = Homme::subview(m_elements.buffers.dpdissk, kv.ie);
    const auto ptens = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    const auto qlim = Homme::subview(m_tracers.qlim, kv.ie, kv.iq);

    if ( ! OnGpu<ExecSpace>::value && kv.team.team_size() == 1)
      SerialLimiter<ExecSpace>::run<9>(
        sphweights, dpmass, qlim, ptens,
        Homme::subview(m_sphere_ops.vector_buf_ml, kv.team_idx, 0));
    else
      limiter_clip_and_sum(kv.team, sphweights, dpmass, qlim, ptens);
  }

  //! apply mass matrix, overwrite np1 with solution:
  //! dont do this earlier, since we allow np1_qdp == n0_qdp
  //! and we dont want to overwrite n0_qdp until we are done using it
  KOKKOS_INLINE_FUNCTION
  void apply_spheremp (const KernelVariables& kv) const {
    const auto qdp = Homme::subview(m_tracers.qdp, kv.ie, m_data.np1_qdp, kv.iq);
    const auto qtens = Homme::subview(m_tracers.qtens_biharmonic, kv.ie, kv.iq);
    const auto spheremp = Homme::subview(m_elements.m_spheremp, kv.ie);
    Kokkos::parallel_for (
      Kokkos::TeamThreadRange(kv.team, NP * NP),
      [&] (const int loop_idx) {
        const int igp = loop_idx / NP;
        const int jgp = loop_idx % NP;
        Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
          [&] (const int& ilev) {
            qdp(igp, jgp, ilev) = spheremp(igp, jgp) * qtens(igp, jgp, ilev);
          });
      });
  }

  // Do all the setup and teardown associated with a limiter. Call a limiter
  // functor to do the actual math given the problem data (mass, minp, maxp, c,
  // x), where the limiter possibly alters x to place it in the constraint set
  //    {x: (i) minp <= x_k <= maxp and (ii) c'x = mass }.
  template <typename Limit, typename ArrayGll, typename ArrayGllLvl, typename Array2Lvl>
  KOKKOS_INLINE_FUNCTION static void
  with_limiter_shell (const TeamMember& team, const Limit& limit,
                      const ArrayGll& sphweights, const ArrayGllLvl& dpmass,
                      const Array2Lvl& qlim, const ArrayGllLvl& ptens) {
    const int NP2 = NP * NP;

    // Size doesn't matter; just need to get a pointer to the start of the
    // shared memory.
    Real* const team_data = Memory<ExecSpace>::get_shmem<Real>(team);

    const auto f = [&] (const int ilev) {
      const int vpi = ilev / VECTOR_SIZE, vsi = ilev % VECTOR_SIZE;

      Real* const data = team_data ?
      team_data + 2 * NP2 * team.team_rank() :
      nullptr;
      Memory<ExecSpace>::AutoArray<Real, NP2> x(data), c(data + NP2);

      Dispatch<>::parallel_for_NP2(team, [&] (const int& k) {
          const int i = k / NP, j = k % NP;
          const auto& dpm = dpmass(i,j,vpi)[vsi];
          c[k] = sphweights(i,j)*dpm;
          x[k] = ptens(i,j,vpi)[vsi]/dpm;
        });

      Kokkos::Real2 sums;
      Dispatch<>::parallel_reduce_NP2(team, [&] (const int& k, Kokkos::Real2& sums) {
          sums.v[0] += x[k]*c[k];
          sums.v[1] += c[k];
        }, sums);
      if (sums.v[1] <= 0) return; //! this should never happen, but if it does, dont limit

      Real minp = qlim(0,vpi)[vsi], maxp = qlim(1,vpi)[vsi];
      // This is a slightly different spot than where this comment came from,
      // but it's logically equivalent to do it here.
      //! IMPOSE ZERO THRESHOLD.  do this here so it can be turned off for
      //! testing
      if (minp < 0)
        minp = qlim(0,vpi)[vsi] = 0;

      //! relax constraints to ensure limiter has a solution:
      //! This is only needed if running with the SSP CFL>1 or
      //! due to roundoff errors
      // This is technically a write race condition, but the same value is
      // being written, so it doesn't matter.
      if (sums.v[0] < minp*sums.v[1])
        minp = qlim(0,vpi)[vsi] = sums.v[0]/sums.v[1];
      if (sums.v[0] > maxp*sums.v[1])
        maxp = qlim(1,vpi)[vsi] = sums.v[0]/sums.v[1];

      const bool modified =
      limit(team, sums.v[0], minp, maxp, x.data(), c.data());

      if (modified)
        Dispatch<>::parallel_for_NP2(team, [&] (const int& k) {
            const int i = k / NP, j = k % NP;
            ptens(i,j,vpi)[vsi] = x[k]*dpmass(i,j,vpi)[vsi];
          });
    };

    if (OnGpu<ExecSpace>::value || team.team_size() > 1) {
      Kokkos::parallel_for (
        Kokkos::TeamThreadRange(team, NUM_PHYSICAL_LEV),
        f);
    } else {
VECTOR_SIMD_LOOP
      for (int ilev = 0; ilev < NUM_PHYSICAL_LEV; ++ilev)
        f(ilev);
    }
  }

public: // Expose for unit testing.

  // limiter_option = 8.
  template <typename ArrayGll, typename ArrayGllLvl, typename Array2Lvl>
  KOKKOS_INLINE_FUNCTION static void
  limiter_optim_iter_full (const TeamMember& team,
                           const ArrayGll& sphweights, const ArrayGllLvl& dpmass,
                           const Array2Lvl& qlim, const ArrayGllLvl& ptens) {
    struct Limit {
      KOKKOS_INLINE_FUNCTION bool
      operator() (const TeamMember& team, const Real& mass,
                  const Real& minp, const Real& maxp,
                  Real* KOKKOS_RESTRICT const x,
                  Real const* KOKKOS_RESTRICT const c) const {
        const int maxiter = NP*NP - 1;
        const Real tol_limiter = 5e-14;

        for (int iter = 0; iter < maxiter; ++iter) {
          Real addmass = 0;
          Dispatch<>::parallel_reduce_NP2(team, [&] (const int& k, Real& addmass) {
              Real delta = 0;
              if (x[k] > maxp) {
                delta = x[k] - maxp;
                x[k] = maxp;
              } else if (x[k] < minp) {
                delta = x[k] - minp;
                x[k] = minp;
              }
              addmass += delta*c[k];
            }, addmass);

          if (std::abs(addmass) <= tol_limiter*std::abs(mass))
            break;

          if (addmass > 0) {
            Real weightssum = 0;
            Dispatch<>::parallel_reduce_NP2(team, [&] (const int& k, Real& weightssum) {
                if (x[k] < maxp)
                  weightssum += c[k];
              }, weightssum);
            const auto adw = addmass/weightssum;
            Dispatch<>::parallel_for_NP2(team, [&] (const int& k) {
                x[k] += (x[k] < maxp) ? adw : 0;
              });
          } else {
            Real weightssum = 0;
            Dispatch<>::parallel_reduce_NP2(team, [&] (const int& k, Real& weightssum) {
                if (x[k] > minp)
                  weightssum += c[k];
              }, weightssum);
            const auto adw = addmass/weightssum;
            Dispatch<>::parallel_for_NP2(team, [&] (const int& k) {
                x[k] += (x[k] > minp) ? adw : 0;
              });
          }
        }
        return true;
      }
    };

    with_limiter_shell(team, Limit(), sphweights, dpmass, qlim, ptens);
  }

  // This is limiter_option = 9 in ACME master.
  template <typename ArrayGll, typename ArrayGllLvl, typename Array2Lvl>
  KOKKOS_INLINE_FUNCTION static void
  limiter_clip_and_sum (const TeamMember& team,
                        const ArrayGll& sphweights, const ArrayGllLvl& dpmass,
                        const Array2Lvl& qlim, const ArrayGllLvl& ptens) {
    struct Limit {
      KOKKOS_INLINE_FUNCTION bool
      operator() (const TeamMember& team, const Real& mass,
                  const Real& minp, const Real& maxp,
                  Real* KOKKOS_RESTRICT const x,
                  Real const* KOKKOS_RESTRICT const c) const {
        // Clip.
        Kokkos::Real2 reds;
        Dispatch<>::parallel_reduce_NP2(team, [&] (const int& k, Kokkos::Real2& reds) {
            Real delta = 0;
            if (x[k] > maxp) {
              delta = x[k] - maxp;
              x[k] = maxp;
              reds.v[1] += 1;
            } else if (x[k] < minp) {
              delta = x[k] - minp;
              x[k] = minp;
              reds.v[1] += 1;
            }
            reds.v[0] += delta*c[k];
          }, reds);
        if (reds.v[0] == 0) return false;
        const Real addmass = reds.v[0];

        if (addmass > 0) {
          Real fac = 0;
          // Get sum of weights. Don't store them; we don't want another array.
          Dispatch<>::parallel_reduce_NP2(team, [&] (const int& k, Real& fac) {
              fac += c[k]*(maxp - x[k]);
            }, fac);
          if (fac > 0) {
            // Update.
            fac = addmass/fac;
            Dispatch<>::parallel_for_NP2(team, [&] (const int& k) {
                x[k] += fac*(maxp - x[k]);
              });
          }
        } else {
          Real fac = 0;
          Dispatch<>::parallel_reduce_NP2(team, [&] (const int& k, Real& fac) {
              fac += c[k]*(x[k] - minp);
            }, fac);
          if (fac > 0) {
            fac = addmass/fac;
            Dispatch<>::parallel_for_NP2(team, [&] (const int& k) {
                x[k] += fac*(x[k] - minp);
              });
          }
        }
        return true;
      }
    };

    with_limiter_shell(team, Limit(), sphweights, dpmass, qlim, ptens);
  }
};

// Code repetition results from needing BFB and slight differences between lim 8
// and 9 Fortran impls.
template <typename ExecSpace>
template <int limiter_option, typename ArrayGll, typename ArrayGllLvl, typename Array2Lvl,
          typename Array2GllLvl>
KOKKOS_INLINE_FUNCTION void SerialLimiter<ExecSpace>
::run (const ArrayGll& sphweights, const ArrayGllLvl& idpmass,
       const Array2Lvl& iqlim, const ArrayGllLvl& iptens,
       const Array2GllLvl& irwrk) {

#define forij for (int i = 0; i < NP; ++i) for (int j = 0; j < NP; ++j)
#define forlev for (int lev = 0; lev < NUM_PHYSICAL_LEV; ++lev)

  ViewUnmanaged<const Real[NP][NP][NUM_LEV*VECTOR_SIZE]>
    dpmass(&idpmass(0,0,0)[0]);
  ViewUnmanaged<Real[NP][NP][NUM_LEV*VECTOR_SIZE]>
    c(&irwrk(0,0,0,0)[0]);
  ViewUnmanaged<Real[2][NUM_LEV*VECTOR_SIZE]>
    qlim(&iqlim(0,0)[0]);

  if (limiter_option == 8) {
    ViewUnmanaged<Real[NP][NP][NUM_LEV*VECTOR_SIZE]>
      x(&iptens(0,0,0)[0]);
    Real mass[NUM_PHYSICAL_LEV] = {0}, sumc[NUM_PHYSICAL_LEV] = {0};

    forij {
      const auto& sphij = sphweights(i,j);
      VECTOR_SIMD_LOOP forlev {
        const auto& dpm = dpmass(i,j,lev);
        c(i,j,lev) = sphij*dpm;
        x(i,j,lev) /= dpm;
        mass[lev] += c(i,j,lev)*x(i,j,lev);
        sumc[lev] += c(i,j,lev);
      }
    }

    VECTOR_SIMD_LOOP forlev {
      if (qlim(0,lev) < 0)
        qlim(0,lev) = 0;
      if (mass[lev] < qlim(0,lev)*sumc[lev])
        qlim(0,lev) = mass[lev]/sumc[lev];
      if (mass[lev] > qlim(1,lev)*sumc[lev])
        qlim(1,lev) = mass[lev]/sumc[lev];
    }

    static const int maxiter = NP*NP - 1;
    static const Real tol_limiter = 5e-14;
    int donecnt = 0;
    char done[NUM_PHYSICAL_LEV] = {0};
    for (int iter = 0; iter < maxiter; ++iter) {
      Real addmass[NUM_PHYSICAL_LEV] = {0};

      forij {
        VECTOR_SIMD_LOOP forlev {
          auto& xij = x(i,j,lev);
          Real delta = 0;
          if (xij < qlim(0,lev)) {
            delta = xij - qlim(0,lev);
            xij = qlim(0,lev);
          } else if (xij > qlim(1,lev)) {
            delta = xij - qlim(1,lev);
            xij = qlim(1,lev);
          }
          addmass[lev] += delta*c(i,j,lev);
        }
      }

      forlev {
        if (std::abs(addmass[lev]) <= tol_limiter*std::abs(mass[lev]) &&
            ! done[lev]) {
          done[lev] = 1;
          ++donecnt;
        }
      }
      if (donecnt == NUM_PHYSICAL_LEV) break;

      Real f[NUM_PHYSICAL_LEV] = {0};
      forij {
        VECTOR_SIMD_LOOP forlev {
          if (done[lev]) continue;
          if (addmass[lev] <= 0) {
            if (x(i,j,lev) > qlim(0,lev))
              f[lev] += c(i,j,lev);
          } else {
            if (x(i,j,lev) < qlim(1,lev))
              f[lev] += c(i,j,lev);
          }
        }
      }

      VECTOR_SIMD_LOOP forlev {
        if (f[lev] != 0)
          f[lev] = addmass[lev] / f[lev];
      }

      forij {
        VECTOR_SIMD_LOOP forlev {
          if (done[lev]) continue;
          if (addmass[lev] <= 0) {
            if (x(i,j,lev) > qlim(0,lev))
              x(i,j,lev) += f[lev];
          } else {
            if (x(i,j,lev) < qlim(1,lev))
              x(i,j,lev) += f[lev];
          }
        }
      }
    }

    forij {
      VECTOR_SIMD_LOOP forlev {
        x(i,j,lev) *= dpmass(i,j,lev);
      }
    }
  } else if (limiter_option == 9) {
    ViewUnmanaged<Real[NP][NP][NUM_LEV*VECTOR_SIZE]>
      ptens(&iptens(0,0,0)[0]);
    ViewUnmanaged<Real[NP][NP][NUM_LEV*VECTOR_SIZE]>
      x(&irwrk(1,0,0,0)[0]);
    Real mass[NUM_PHYSICAL_LEV] = {0}, sumc[NUM_PHYSICAL_LEV] = {0};

    forij {
      const auto& sphij = sphweights(i,j);
      VECTOR_SIMD_LOOP forlev {
        const auto& dpm = dpmass(i,j,lev);
        c(i,j,lev) = sphij*dpm;
        x(i,j,lev) = ptens(i,j,lev) / dpm;
        mass[lev] += c(i,j,lev)*x(i,j,lev);
        sumc[lev] += c(i,j,lev);
      }
    }

    VECTOR_SIMD_LOOP forlev {
      if (qlim(0,lev) < 0)
        qlim(0,lev) = 0;
      if (mass[lev] < qlim(0,lev)*sumc[lev])
        qlim(0,lev) = mass[lev]/sumc[lev];
      if (mass[lev] > qlim(1,lev)*sumc[lev])
        qlim(1,lev) = mass[lev]/sumc[lev];
    }

    Real addmass[NUM_PHYSICAL_LEV] = {0};
    // This is here to be BFB with the Fortran impl. When BFB is no longer
    // required, it is almost certainly better to remove this variable and the
    // conditionals on it.
    int modified[NUM_PHYSICAL_LEV] = {0};

    forij {
      VECTOR_SIMD_LOOP forlev {
        auto& xij = x(i,j,lev);
        Real delta = 0;
        if (xij < qlim(0,lev)) {
          delta = xij - qlim(0,lev);
          xij = qlim(0,lev);
          modified[lev] = 1;
        } else if (xij > qlim(1,lev)) {
          delta = xij - qlim(1,lev);
          xij = qlim(1,lev);
          modified[lev] = 1;
        }
        if (modified[lev])
          addmass[lev] += delta*c(i,j,lev);
      }
    }

    Real f[NUM_PHYSICAL_LEV] = {0};
    forij {
      VECTOR_SIMD_LOOP forlev {
        auto& xij = x(i,j,lev);
        if (addmass[lev] <= 0) {
          if (xij > qlim(0,lev))
            f[lev] += c(i,j,lev)*(xij - qlim(0,lev));
        } else {
          if (xij < qlim(1,lev))
            f[lev] += c(i,j,lev)*(qlim(1,lev) - xij);
        }
      }
    }
    
    VECTOR_SIMD_LOOP forlev {
      if (f[lev] != 0)
        f[lev] = addmass[lev] / f[lev];
    }

    forij {
      VECTOR_SIMD_LOOP forlev {
        auto& xij = x(i,j,lev);
        if (addmass[lev] <= 0) {
          if (xij > qlim(0,lev))
            xij += f[lev]*(xij - qlim(0,lev));
        } else {
          if (xij < qlim(1,lev))
            xij += f[lev]*(qlim(1,lev) - xij);
        }
      }
    }

    forij {
      VECTOR_SIMD_LOOP forlev {
        if (modified[lev])
          ptens(i,j,lev) = x(i,j,lev) * dpmass(i,j,lev);
      }
    }
  } else {
    Kokkos::abort("Only limiter_option 8 and 9 is impl'ed.");
  }

# undef forlev
# undef forij
} // end SerialLimiter

}

#endif
