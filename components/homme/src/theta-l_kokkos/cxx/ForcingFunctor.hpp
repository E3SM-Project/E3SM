/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"
#include "CamForcing.hpp"
#include "ColumnOps.hpp"
#include "ElementOps.hpp"
#include "ElementsForcing.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "EquationOfState.hpp"
#include "FunctorsBuffersManager.hpp"
#include "Tracers.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "KernelVariables.hpp"
#include "profiling.hpp"

#include "utilities/BfbUtils.hpp"

namespace Homme {

class ForcingFunctor
{
public:

  struct TagStates {};
  struct TagTracersPre {};
  struct TagTracers {};
  struct TagTracersPost {};

  ForcingFunctor ()
   : m_policy_states(1,1)       // Need to init in constructor. Fix sizes in the body.
   , m_policy_tracers_pre(1,1)
   , m_policy_tracers(1,1)
   , m_policy_tracers_post(1,1)
   , m_tu_states(m_policy_states)
   , m_tu_tracers_pre(m_policy_tracers_pre)
   , m_tu_tracers(m_policy_tracers)
   , m_tu_tracers_post(m_policy_tracers_post)
  {
    const auto& c = Context::singleton();
    const auto& p = c.get<SimulationParams>();

    m_hydrostatic = p.theta_hydrostatic_mode;
    m_qsize = p.qsize;

    m_state    = c.get<ElementsState>();
    m_forcing  = c.get<ElementsForcing>();
    m_geometry = c.get<ElementsGeometry>();
    m_tracers  = c.get<Tracers>();
    m_hvcoord  = c.get<HybridVCoord>();

    // TODO: this may change, depending on the simulation params
    m_adjust_ps = true;

    m_eos.init(m_hydrostatic,m_hvcoord);
    m_elem_ops.init(m_hvcoord);

    // Check everything is init-ed
    assert (m_state.num_elems()>0);
    assert (m_forcing.num_elems()>0);
    assert (m_tracers.num_elems()>0);
    assert (m_geometry.num_elems()>0);
    assert (m_hvcoord.m_inited);

    // Create the team policy
    m_policy_states = Homme::get_default_team_policy<ExecSpace,TagStates>(m_state.num_elems());
    m_policy_tracers_pre = Homme::get_default_team_policy<ExecSpace,TagTracersPre>(m_tracers.num_elems());
    m_policy_tracers = Homme::get_default_team_policy<ExecSpace,TagTracers>(m_tracers.num_elems()*m_tracers.num_tracers());
    m_policy_tracers_post = Homme::get_default_team_policy<ExecSpace,TagTracersPost>(m_tracers.num_elems());

    m_tu_states       = TeamUtils<ExecSpace>(m_policy_states);
    m_tu_tracers_pre  = TeamUtils<ExecSpace>(m_policy_tracers_pre);
    m_tu_tracers      = TeamUtils<ExecSpace>(m_policy_tracers);
    m_tu_tracers_post = TeamUtils<ExecSpace>(m_policy_tracers_post);
  }

  int requested_buffer_size () const {
    const int nteams = m_tu_tracers.get_num_concurrent_teams();
    const int nelems = m_state.num_elems();
    constexpr int mid_size = NP*NP*NUM_LEV*VECTOR_SIZE;
    constexpr int int_size = NP*NP*NUM_LEV_P*VECTOR_SIZE;

    // 3 persistent midlayers, 2 non-persistent midlayer, and 1 non-persistent interface
    return mid_size*(nelems*4+nteams) + (m_hydrostatic ? int_size*nteams : 0);
  }

  void init_buffers (const FunctorsBuffersManager& fbm) {
    const int num_teams = m_tu_tracers.get_num_concurrent_teams();
    const int num_elems = m_state.num_elems();

    constexpr int mid_size = NP*NP*NUM_LEV;

    Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());

    m_tn1 = decltype(m_tn1)(mem,num_elems);
    mem += mid_size*num_elems;

    m_dp_adj = decltype(m_dp_adj)(mem,num_elems);
    mem += mid_size*num_elems;

    m_pnh = decltype(m_pnh)(mem,num_elems);
    mem += mid_size*num_elems;

    m_exner = decltype(m_exner)(mem,num_elems);
    mem += mid_size*num_elems;

    m_Rstar = decltype(m_Rstar)(mem,num_teams);
    mem += mid_size*num_teams;

    m_pi_i = decltype(m_pi_i)(mem,num_teams);
  }

  void states_forcing (const Real dt, const int np1) {
    m_dt = dt;
    m_np1 = np1;

    Kokkos::parallel_for("compute states forcing",m_policy_states,*this);
    Kokkos::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagStates&, const TeamMember& team) const {
    constexpr int LAST_MID_PACK     = ColInfo<NUM_PHYSICAL_LEV>::LastPack;
    constexpr int LAST_MID_PACK_END = ColInfo<NUM_PHYSICAL_LEV>::LastPackEnd;
    constexpr int LAST_INT_PACK     = ColInfo<NUM_INTERFACE_LEV>::LastPack;
    constexpr int LAST_INT_PACK_END = ColInfo<NUM_INTERFACE_LEV>::LastPackEnd;

    KernelVariables kv(team, m_tu_states);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {

      const int igp  = idx / NP;
      const int jgp  = idx % NP;

      const auto fm_x = Homme::subviewConst(m_forcing.m_fm,kv.ie,0,igp,jgp);
      const auto fm_y = Homme::subviewConst(m_forcing.m_fm,kv.ie,1,igp,jgp);
      const auto fm_z = Homme::subviewConst(m_forcing.m_fm,kv.ie,2,igp,jgp);
      const auto fvtheta = Homme::subviewConst(m_forcing.m_fvtheta,kv.ie,igp,jgp);
      const auto fphi    = Homme::subviewConst(m_forcing.m_fphi,kv.ie,igp,jgp);

      auto u      = Homme::subview(m_state.m_v,kv.ie,m_np1,0,igp,jgp);
      auto v      = Homme::subview(m_state.m_v,kv.ie,m_np1,1,igp,jgp);
      auto w      = Homme::subview(m_state.m_w_i,kv.ie,m_np1,igp,jgp);
      auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_np1,igp,jgp);
      auto phi    = Homme::subview(m_state.m_phinh_i,kv.ie,m_np1,igp,jgp);

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        vtheta(ilev) += m_dt*fvtheta(ilev);
        phi(ilev) += m_dt*fphi(ilev);

        u(ilev) += m_dt*fm_x(ilev);
        v(ilev) += m_dt*fm_y(ilev);
        w(ilev) += m_dt*fm_z(ilev);
      });

      // Fix w at the surface
      w(LAST_INT_PACK)[LAST_INT_PACK_END] =
        (u(LAST_MID_PACK)[LAST_MID_PACK_END]*m_geometry.m_gradphis(kv.ie,0,igp,jgp) +
         v(LAST_MID_PACK)[LAST_MID_PACK_END]*m_geometry.m_gradphis(kv.ie,1,igp,jgp)) / PhysicalConstants::g;
    });
  }

  void tracers_forcing (const Real dt, const int np1, const int np1_qdp, const bool adjustment, const MoistDry moisture) {
    m_dt = dt;
    m_np1 = np1;
    m_np1_qdp = np1_qdp;
    m_adjustment = adjustment;
    m_moist = (moisture==MoistDry::MOIST);

    Kokkos::parallel_for("temperature, NH perturb press, FQps",m_policy_tracers_pre,*this);
    Kokkos::fence();

    Kokkos::parallel_for("apply tracers forcing", m_policy_tracers,*this);
    Kokkos::fence();

    Kokkos::parallel_for("update temperature, pressure and phi", m_policy_tracers_post,*this);
    Kokkos::fence();
  }

  KOKKOS_INLINE_FUNCTION
  Real compute_fqdt (const int& k,
                     const ExecViewUnmanaged<Scalar[NUM_LEV]>& fq,
                     const ExecViewUnmanaged<Scalar[NUM_LEV]>& qdp) const {
    const int ilev = k / VECTOR_SIZE;
    const int ivec = k % VECTOR_SIZE;
    Real fqdt = m_dt * fq(ilev)[ivec];
    const Real& qdp_s = qdp(ilev)[ivec];
    if (qdp_s + fqdt < 0.0 && fqdt < 0.0) {
      if (qdp_s < 0.0) {
        fqdt = 0.0;
      } else {
        fqdt = -qdp_s;
      }
    }
    return fqdt;
  }

  KOKKOS_INLINE_FUNCTION
  Scalar compute_fqdt_pack (const int& ilev,
                            const ExecViewUnmanaged<Scalar[NUM_LEV]>& fq,
                            const ExecViewUnmanaged<Scalar[NUM_LEV]>& qdp) const {
    Scalar fqdt = m_dt * fq(ilev);
    const Scalar& qdp_s = qdp(ilev);
    // NOTE: here is where masks for simd operations would be handy
    VECTOR_SIMD_LOOP
    for (int iv=0; iv<VECTOR_SIZE; ++iv) {
      if (qdp_s[iv] + fqdt[iv] < 0.0 && fqdt[iv] < 0.0) {
        if (qdp_s[iv] < 0.0) {
          fqdt[iv] = 0.0;
        } else {
          fqdt[iv] = -qdp_s[iv];
        }
      }
    }
    return fqdt;
  }

  KOKKOS_INLINE_FUNCTION
  void operator () (const TagTracersPre&, const TeamMember& team) const {
    KernelVariables kv(team, m_tu_tracers_pre);
    constexpr Real Rgas = PhysicalConstants::Rgas;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto dp = Homme::subview(m_state.m_dp3d,kv.ie,m_np1,igp,jgp);
      Real& ps = m_state.m_ps_v(kv.ie,m_np1,igp,jgp);

      // The hydrostatic pressure in compute_pnh_and_exner in EOS is only used
      // for the theta_hydrostatic_mode case. So only compute it then
      auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_np1,igp,jgp);
      auto phinh  = Homme::subview(m_state.m_phinh_i,kv.ie,m_np1,igp,jgp);
      auto pnh    = Homme::subview(m_pnh,kv.ie,igp,jgp);
      auto exner  = Homme::subview(m_exner,kv.ie,igp,jgp);
      if (m_hydrostatic) {
      auto p_i = Homme::subview(m_pi_i,kv.team_idx,igp,jgp);
        m_elem_ops.compute_hydrostatic_p(kv,dp,p_i,pnh);
        m_eos.compute_exner(kv,pnh,exner);
      } else {
        m_eos.compute_pnh_and_exner(kv,vtheta,phinh,pnh,exner);
      }

      if (m_moist) {
        // Compute pprime=pnh-pi, store in pnh (this is 0 for hydrostatic mode)
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                             [&](const int ilev) {
          pnh(ilev) -= m_hvcoord.ps0*m_hvcoord.hybrid_am(ilev) +
                                 ps *m_hvcoord.hybrid_bm(ilev);
        });
      }

      // Compute Rstar
      auto Rstar = Homme::subview(m_Rstar,kv.team_idx,igp,jgp);
      m_elem_ops.get_R_star (kv, m_moist,
                             Homme::subview(m_tracers.Q,kv.ie,0,igp,jgp),
                             Rstar);

      // Compute temperature
      auto tn1 = Homme::subview(m_tn1,kv.ie,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        tn1(ilev) = exner(ilev) * vtheta(ilev)*(Rgas/Rstar(ilev)) / dp(ilev);
      });

      if (m_moist) {
        auto fq = Homme::subview(m_tracers.fq,kv.ie,0,igp,jgp);
        auto q = Homme::subview(m_tracers.Q,kv.ie,0,igp,jgp);
        auto qdp = Homme::subview(m_tracers.qdp,kv.ie,m_np1_qdp,0,igp,jgp);
        auto dp_adj = Homme::subview(m_dp_adj,kv.ie,igp,jgp);
        if (m_adjustment) {
          Dispatch<ExecSpace>::parallel_reduce(
            kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
            [&](const int &k, Real &accumulator) {
              const int ilev = k / VECTOR_SIZE;
              const int ivec = k % VECTOR_SIZE;
              accumulator += dp(ilev)[ivec]*(fq(ilev)[ivec]-q(ilev)[ivec]);
            },ps);
          if (!m_adjust_ps) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                                 [&](const int ilev) {
              dp_adj(ilev) = dp(ilev) + dp(ilev)*(fq(ilev)-q(ilev));
            });
          }
        } else {
          Real ps_forcing = 0.0;
          Dispatch<ExecSpace>::parallel_reduce(
            kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
            [&](const int &k, Real &accumulator) {
              accumulator += compute_fqdt(k,fq,qdp)/m_dt;
            },ps_forcing);
          ps += ps_forcing*m_dt;
          if (!m_adjust_ps) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                                 [&](const int& ilev) {
              dp_adj(ilev) = dp(ilev) + compute_fqdt_pack(ilev,fq,qdp);
            });
          }
        }

        if (m_adjust_ps) {
          m_hvcoord.compute_dp_ref(kv,ps,dp_adj);
        }
      }
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagTracers&, const TeamMember& team) const {
    KernelVariables kv(team, m_qsize, m_tu_tracers);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto qdp    = Homme::subview(m_tracers.qdp,kv.ie,m_np1_qdp,kv.iq,igp,jgp);
      auto fq     = Homme::subview(m_tracers.fq, kv.ie,kv.iq,igp,jgp);
      auto Q      = Homme::subview(m_tracers.Q, kv.ie,kv.iq,igp,jgp);
      auto dp     = Homme::subview(m_state.m_dp3d, kv.ie,m_np1,igp,jgp);
      auto dp_adj = Homme::subview(m_dp_adj, kv.ie,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        Scalar& qdp_ilev = qdp(ilev);
        if (m_adjustment) {
          qdp_ilev = dp(ilev)*fq(ilev);
        } else {
          qdp_ilev += compute_fqdt_pack(ilev,fq,qdp);
        }

        if (m_moist) {
          dp(ilev) = dp_adj(ilev);
        }

        // Update tracers concentration
        Q(ilev) = qdp_ilev/dp(ilev);
      });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagTracersPost&, const TeamMember& team) const {
    constexpr int LAST_INT_PACK     = ColInfo<NUM_INTERFACE_LEV>::LastPack;
    constexpr int LAST_INT_PACK_END = ColInfo<NUM_INTERFACE_LEV>::LastPackEnd;

    KernelVariables kv(team, m_tu_tracers_post);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // Compute Rstar
      auto Rstar = Homme::subview(m_Rstar,kv.team_idx,igp,jgp);
      m_elem_ops.get_R_star (kv, m_moist,
                             Homme::subview(m_tracers.Q,kv.ie,0,igp,jgp),
                             Rstar);

      auto tn1    = Homme::subview(m_tn1,kv.ie,igp,jgp);
      auto pnh    = Homme::subview(m_pnh,kv.ie,igp,jgp);
      auto exner  = Homme::subview(m_exner,kv.ie,igp,jgp);
      auto dp     = Homme::subview(m_state.m_dp3d,kv.ie,m_np1,igp,jgp);
      auto dp_adj = Homme::subview(m_dp_adj,kv.ie,igp,jgp);
      auto ft     = Homme::subview(m_forcing.m_ft,kv.ie,igp,jgp);
      auto fphi   = Homme::subview(m_forcing.m_fphi,kv.ie,igp,jgp);
      const Real ps = m_state.m_ps_v(kv.ie,m_np1,igp,jgp);

      if (m_moist) {
        // Need to update pnh and exner. Currently, pnh is storing pnh-pi

        if (m_adjust_ps) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                               [&](const int ilev) {
            pnh(ilev) += m_hvcoord.ps0*m_hvcoord.hybrid_am(ilev) +
                                   ps *m_hvcoord.hybrid_bm(ilev);
          });
        } else {
          // Compute hydrostatic p from dp. Store in exner, then add to pnh
          auto p_i = Homme::subview(m_pi_i,kv.team_idx,igp,jgp);
          m_elem_ops.compute_hydrostatic_p(kv,dp,p_i,exner);
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                               [&](const int ilev) {
            pnh(ilev) += exner(ilev);
          });
        }
        
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                             [&](const int ilev) {
          using namespace PhysicalConstants;
#ifdef HOMMEXX_BFB_TESTING
          exner(ilev) = bfb_pow(pnh(ilev)/p0,kappa);
#else
          exner(ilev) = pow(pnh(ilev)/p0,kappa);
#endif
        });
      }

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {

        // Update temperature
        tn1(ilev) += m_dt*ft(ilev);

        // Compute theta, store in tn1.
        tn1(ilev) = (Rstar(ilev)/PhysicalConstants::Rgas)*tn1(ilev)*dp(ilev)/exner(ilev);
      });

      // Compute phi as fcn of new theta, store in fphi
      auto phi_i  = Homme::subview(m_state.m_phinh_i,kv.ie,m_np1,igp,jgp);
      fphi(LAST_INT_PACK)[LAST_INT_PACK_END] = phi_i(LAST_INT_PACK)[LAST_INT_PACK_END];
      auto integrand_provider = [&](const int ilev)->Scalar {
        return PhysicalConstants::Rgas*tn1(ilev)*exner(ilev)/pnh(ilev);
      };
      ColumnOps::column_scan_mid_to_int<false>(kv,integrand_provider,fphi);

      // Finally, update forcing for theta and phi_i
      auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_np1,igp,jgp);
      auto ftheta = Homme::subview(m_forcing.m_fvtheta,kv.ie,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        ftheta(ilev) = (tn1(ilev) - vtheta(ilev)) / m_dt;
        fphi(ilev) = (fphi(ilev) - phi_i(ilev)) / m_dt;
      });

      // If NUM_LEV=NUM_LEV_P, the code above already took care of the last level
      if (NUM_LEV!=NUM_LEV_P) {
        // Last level
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          fphi(LAST_INT_PACK)[LAST_INT_PACK_END] = (fphi(LAST_INT_PACK)[LAST_INT_PACK_END] - phi_i(LAST_INT_PACK)[LAST_INT_PACK_END]) / m_dt;
        });
      }
    });
  }

private:

  ElementsState     m_state;
  ElementsForcing   m_forcing;
  ElementsGeometry  m_geometry;
  Tracers           m_tracers;

  HybridVCoord      m_hvcoord;
  ElementOps        m_elem_ops;
  EquationOfState   m_eos;

  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_Rstar;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_dp_adj;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_pnh;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV_P]> m_pi_i;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_exner;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_tn1;

  int m_qsize;
  int m_np1;
  int m_np1_qdp;
  bool m_moist;
  bool m_hydrostatic;
  bool m_adjustment;
  bool m_adjust_ps;
  Real m_dt;

  Kokkos::TeamPolicy<ExecSpace,TagStates>       m_policy_states;
  Kokkos::TeamPolicy<ExecSpace,TagTracersPre>   m_policy_tracers_pre;
  Kokkos::TeamPolicy<ExecSpace,TagTracers>      m_policy_tracers;
  Kokkos::TeamPolicy<ExecSpace,TagTracersPost>  m_policy_tracers_post;
  TeamUtils<ExecSpace> m_tu_states, m_tu_tracers_pre, m_tu_tracers, m_tu_tracers_post;
};

} // namespace Homme
