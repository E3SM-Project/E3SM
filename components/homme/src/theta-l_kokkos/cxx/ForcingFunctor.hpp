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

    m_eos.init(m_hydrostatic,m_hvcoord);

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
  }

  int requested_buffer_size () const {
    const int nteams = get_num_concurrent_teams(m_policy_tracers);
    const int nelems = m_state.num_elems();
    constexpr int mid_size = NP*NP*NUM_LEV*VECTOR_SIZE;
    constexpr int int_size = NP*NP*NUM_LEV_P*VECTOR_SIZE;

    // 3 persistent midlayers, 2 non-persistent midlayer, and 1 non-persistent interface
    return mid_size*(nelems*3+nteams*2) + (m_hydrostatic ? int_size*nteams : 0);
  }

  void init_buffers (const FunctorsBuffersManager& fbm) {
    const int num_teams = get_num_concurrent_teams(m_policy_tracers);
    const int num_elems = m_state.num_elems();

    constexpr int mid_size = NP*NP*NUM_LEV;

    Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());

    m_tn1 = decltype(m_Rstar)(mem,num_elems);
    mem += mid_size*num_elems;

    m_dp = decltype(m_dp)(mem,num_elems);
    mem += mid_size*num_elems;

    m_pnh = decltype(m_pnh)(mem,num_elems);
    mem += mid_size*num_elems;

    m_exner = decltype(m_exner)(mem,num_teams);
    mem += mid_size*num_teams;

    m_Rstar = decltype(m_exner)(mem,num_teams);
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
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {

      const int igp  = idx / NP;
      const int jgp  = idx % NP;

      auto fm_x = Homme::subviewConst(m_forcing.m_fm,kv.ie,0,igp,jgp);
      auto fm_y = Homme::subviewConst(m_forcing.m_fm,kv.ie,1,igp,jgp);
      auto fm_z = Homme::subviewConst(m_forcing.m_fm,kv.ie,2,igp,jgp);
      auto fvtheta = Homme::subviewConst(m_forcing.m_fvtheta,kv.ie,igp,jgp);
      auto fphi    = Homme::subviewConst(m_forcing.m_fphi,kv.ie,igp,jgp);

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
      constexpr int last_int_pack = ColInfo<NUM_INTERFACE_LEV>::LastPack;
      constexpr int last_int_vec_end = ColInfo<NUM_INTERFACE_LEV>::LastVecEnd;
      constexpr int last_mid_pack = ColInfo<NUM_PHYSICAL_LEV>::LastPack;
      constexpr int last_mid_vec_end = ColInfo<NUM_PHYSICAL_LEV>::LastVecEnd;

      w(last_int_pack)[last_int_vec_end] =
        u(last_mid_pack)[last_mid_vec_end]*m_geometry.m_gradphis(kv.ie,0,igp,jgp) +
        v(last_mid_pack)[last_mid_vec_end]*m_geometry.m_gradphis(kv.ie,1,igp,jgp);
    });
  }

  void tracers_forcing (const Real dt, const int np1, const int np1_qdp, const MoistDry moisture) {
    m_dt = dt;
    m_np1 = np1;
    m_np1_qdp = np1_qdp;
    m_moist = moisture==MoistDry::MOIST;

    Kokkos::parallel_for("temperature, NH perturb press, FQps",m_policy_tracers_pre,*this);
    Kokkos::fence();

    Kokkos::parallel_for("apply tracers forcing", m_policy_tracers,*this);
    Kokkos::fence();

    Kokkos::parallel_for("update temperature, pressure and phi", m_policy_tracers_post,*this);
    Kokkos::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator () (const TagTracersPre&, const TeamMember& team) const {
    KernelVariables kv(team);
    constexpr Real Rgas = PhysicalConstants::Rgas;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto dp = Homme::subview(m_dp,kv.ie,igp,jgp);
      Real& ps_v = m_state.m_ps_v(kv.ie,m_np1,igp,jgp);

      // Compute dp from current ps_v
      m_hvcoord.compute_dp_ref(kv,ps_v, dp);

      // The hydrostatic pressure in compute_pnh_and_exner in EOS is only used
      // for the theta_hydrostatic_mode case. So only compute it then
      auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_np1,igp,jgp);
      auto phinh  = Homme::subview(m_state.m_phinh_i,kv.ie,m_np1,igp,jgp);
      auto pnh    = Homme::subview(m_pnh,kv.ie,igp,jgp);
      auto exner  = Homme::subview(m_exner,kv.team_idx,igp,jgp);
      if (m_hydrostatic) {
        auto p_i = Homme::subview(m_pi_i,kv.team_idx,igp,jgp);
        m_eos.compute_hydrostatic_p(kv,dp,p_i,pnh);
        m_eos.compute_exner(kv,pnh,exner);
      } else {
        m_eos.compute_pnh_and_exner(kv,vtheta,phinh,pnh,exner);
      }
      // Compute pprime, store in pnh
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        pnh(ilev) -= (m_hvcoord.ps0*m_hvcoord.hybrid_am(ilev) +
                      ps_v*m_hvcoord.hybrid_bm(ilev));
      });

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
        // This conserves the dry mass in the process of forcing tracer 0
        Real ps_v_forcing = 0.0;

        auto fq = Homme::subview(m_tracers.fq,kv.ie,0,igp,jgp);
        auto qdp = Homme::subview(m_tracers.qdp,kv.ie,m_np1_qdp,0,igp,jgp);
        Dispatch<ExecSpace>::parallel_reduce(
            kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
            [&](const int &k, Real &accumulator) {
              const int ilev = k / VECTOR_SIZE;
              const int vlev = k % VECTOR_SIZE;
              Real v1 = m_dt * fq(ilev)[vlev];
              const Real& qdp_s = qdp(ilev)[vlev];
              if (qdp_s + v1 < 0.0 && v1 < 0.0) {
                if (qdp_s < 0.0) {
                  v1 = 0.0;
                } else {
                  v1 = -qdp_s;
                }
              }
              accumulator += v1;
            },
            ps_v_forcing);
        ps_v += ps_v_forcing;

        // Now update dp
        m_hvcoord.compute_dp_ref(kv,ps_v,dp);
      }
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagTracers&, const TeamMember& team) const {
    KernelVariables kv(team,m_qsize);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto qdp = Homme::subview(m_tracers.qdp,kv.ie,m_np1_qdp,kv.iq,igp,jgp);
      auto fq  = Homme::subview(m_tracers.fq, kv.ie,kv.iq,igp,jgp);
      auto Q   = Homme::subview(m_tracers.Q, kv.ie,kv.iq,igp,jgp);
      auto dp  = Homme::subview(m_dp, kv.ie,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        Scalar v1 = m_dt*fq(ilev);
        Scalar& qdp_s = qdp(ilev);
        // TODO: here is where masked loops as done in scream::Pack would be useful
        VECTOR_SIMD_LOOP
        for (int iv=0; iv<VECTOR_SIZE; ++iv) {
          if (qdp_s[iv] + v1[iv] < 0.0 && v1[iv] < 0.0) {
            if (qdp_s[iv] < 0.0) {
              v1[iv] = 0.0;
            } else {
              v1[iv] = -qdp_s[iv];
            }
          }
        }
        // Not sure if the above loop is guaranteed to vectorize,
        // so I pull out
        qdp_s += v1;

        // Update tracers concentration
        Q(ilev) = qdp_s/dp(ilev);
      });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagTracersPost&, const TeamMember& team) const {
    KernelVariables kv(team);
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
      auto exner  = Homme::subview(m_exner,kv.team_idx,igp,jgp);
      auto dp     = Homme::subview(m_dp,kv.ie,igp,jgp);
      auto ft     = Homme::subview(m_forcing.m_ft,kv.ie,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {

        // Update temperature
        tn1(ilev) += m_dt*ft(ilev);
      
        // Update pnh and recompute exner
        pnh(ilev) += m_hvcoord.ps0*m_hvcoord.hybrid_am(ilev) + 
                     m_state.m_ps_v(kv.ie,m_np1,igp,jgp)*m_hvcoord.hybrid_bm(ilev);
        exner(ilev) = pow(pnh(ilev)/PhysicalConstants::p0,PhysicalConstants::kappa);

        // Compute theta, store in tn1.
        tn1(ilev) = (Rstar(ilev)/PhysicalConstants::Rgas)*tn1(ilev)*dp(ilev)/exner(ilev);
      });

      // Compute phi as fcn of new theta, store in fphi
      auto phi_i  = Homme::subview(m_state.m_phinh_i,kv.ie,m_np1,igp,jgp);
      auto fphi   = Homme::subview(m_forcing.m_fphi,kv.ie,igp,jgp);
      using Info = ColInfo<NUM_INTERFACE_LEV>;
      fphi(Info::LastPack)[Info::LastVecEnd] = phi_i(Info::LastPack)[Info::LastVecEnd];
      auto integrand_provider = [&](const int ilev)->Scalar {
        return PhysicalConstants::Rgas*tn1(ilev)*exner(ilev)/pnh(ilev);
      };
      m_col_ops.column_scan_mid_to_int<false>(kv,integrand_provider,fphi);

      // Finally, update forcing for theta and phi_i
      auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_np1,igp,jgp);
      auto ftheta = Homme::subview(m_forcing.m_fvtheta,kv.ie,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        ftheta(ilev) = (tn1(ilev) - vtheta(ilev)) / m_dt;
        fphi(ilev) = (fphi(ilev) - phi_i(ilev)) / m_dt;
      });

      // Last level
      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        fphi(Info::LastPack) = (fphi(Info::LastPack) - phi_i(Info::LastPack)) / m_dt;
      });
    });
  }

private:

  ElementsState     m_state;
  ElementsForcing   m_forcing;
  ElementsGeometry  m_geometry;
  Tracers           m_tracers;

  HybridVCoord      m_hvcoord;
  ColumnOps         m_col_ops;
  ElementOps        m_elem_ops;
  EquationOfState   m_eos;

  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_Rstar;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_dp;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_pnh;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV_P]> m_pi_i;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_exner;
  ExecViewUnmanaged<Scalar*[NP][NP][NUM_LEV]>   m_tn1;

  int m_qsize;
  int m_np1;
  int m_np1_qdp;
  bool m_moist;
  bool m_hydrostatic;
  Real m_dt;

  Kokkos::TeamPolicy<ExecSpace,TagStates>       m_policy_states;
  Kokkos::TeamPolicy<ExecSpace,TagTracersPre>   m_policy_tracers_pre;
  Kokkos::TeamPolicy<ExecSpace,TagTracers>      m_policy_tracers;
  Kokkos::TeamPolicy<ExecSpace,TagTracersPost>  m_policy_tracers_post;
};

} // namespace Homme
