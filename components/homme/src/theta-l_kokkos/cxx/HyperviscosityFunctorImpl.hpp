/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP
#define HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP

#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ElementsDerivedState.hpp"
#include "ColumnOps.hpp"
#include "EquationOfState.hpp"
#include "ElementOps.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "SimulationParams.hpp"
#include "SphereOperators.hpp"

#include "utilities/VectorUtils.hpp"

#include <memory>

namespace Homme
{

class BoundaryExchange;
struct FunctorsBuffersManager;

class HyperviscosityFunctorImpl
{
  // TODO: don't pass nu_ratio1/2. Instead, do like in F90: compute them from
  //       nu, nu_div, and hv_scaling
  struct HyperviscosityData {
    HyperviscosityData(const int hypervis_subcycle_in, const Real nu_ratio1_in, const Real nu_ratio2_in, const Real nu_top_in,
                       const Real nu_in, const Real nu_p_in, const Real nu_s_in,
                       const Real hypervis_scaling_in)
                      : hypervis_subcycle(hypervis_subcycle_in), nu_ratio1(nu_ratio1_in), nu_ratio2(nu_ratio2_in)
                      , nu_top(nu_top_in), nu(nu_in), nu_p(nu_p_in), nu_s(nu_s_in)
                      , consthv(hypervis_scaling_in == 0){}

    const int   hypervis_subcycle;

    Real  nu_ratio1;
    Real  nu_ratio2;

    const Real  nu_top;
    const Real  nu;
    const Real  nu_p;
    const Real  nu_s;

    int         np1; // The time-level on which to apply hv
    Real        dt;

    Real        eta_ave_w;

    bool consthv;
  };

  static constexpr int NUM_BIHARMONIC_PHYSICAL_LEVELS = 3;
  static constexpr int NUM_BIHARMONIC_LEV = ColInfo<NUM_BIHARMONIC_PHYSICAL_LEVELS>::NumPacks;

  struct Buffers {
    ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>    dptens;
    ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>    ttens;
    ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>    wtens;
    ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>    phitens;
    ExecViewManaged<Scalar * [2][NP][NP][NUM_LEV]> vtens;

    ExecViewManaged<Scalar * [NP][NP][NUM_BIHARMONIC_LEV]>    lapl_dp;
    ExecViewManaged<Scalar * [NP][NP][NUM_BIHARMONIC_LEV]>    lapl_theta;
    ExecViewManaged<Scalar * [NP][NP][NUM_BIHARMONIC_LEV]>    lapl_w;
    ExecViewManaged<Scalar * [NP][NP][NUM_BIHARMONIC_LEV]>    lapl_phi;
    ExecViewManaged<Scalar * [2][NP][NP][NUM_BIHARMONIC_LEV]> lapl_v;
  };

public:

  struct TagFirstLaplaceHV {};
  struct TagSecondLaplaceConstHV {};
  struct TagSecondLaplaceTensorHV {};
  struct TagUpdateStates {};
  struct TagApplyInvMass {};
  struct TagHyperPreExchange {};

  HyperviscosityFunctorImpl (const SimulationParams&     params,
                             const ElementsGeometry&     geometry,
                             const ElementsState&        state,
                             const ElementsDerivedState& derived);

  int requested_buffer_size () const;
  void init_buffers (const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges();

  void run (const int np1, const Real dt, const Real eta_ave_w);

  void biharmonic_wk_theta () const;

  // first iter of laplace, const hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagFirstLaplaceHV&, const TeamMember& team) const {
    using IntColumn = decltype(Homme::subview(m_state.m_w_i,0,0,0,0));

    KernelVariables kv(team, m_tu);
    // Subtract the reference states from the states
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1,igp,jgp);
      auto dp    = Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp);

      auto theta_ref = Homme::subview(m_state.m_ref_states.theta_ref,kv.ie,igp,jgp);
      auto dp_ref    = Homme::subview(m_state.m_ref_states.dp_ref,kv.ie,igp,jgp);

      IntColumn phi_i, phi_i_ref;

      if (m_process_nh_vars) {
        phi_i = Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1,igp,jgp);
        phi_i_ref = Homme::subview(m_state.m_ref_states.phi_i_ref,kv.ie,igp,jgp);
      }
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        vtheta(ilev) -= theta_ref(ilev);
        dp(ilev)     -= dp_ref(ilev);
        if (m_process_nh_vars) {
          phi_i(ilev)  -= phi_i_ref(ilev);
        }
      });

#ifndef XX_NONBFB_COMING
      // It would be fine to not even bother with the surface level, since
      // phitens is only NUM_LEV long, so all the hv stuff does not even happen
      // at NUM_LEV_P (unless NUM_LEV_P==NUM_LEV). However, removing the subtraction
      // and addition of phi_i_ref at NUM_LEV_P introduces NON BFB diffs.
      if (m_process_nh_vars && NUM_LEV!=NUM_LEV_P) {
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          phi_i(NUM_LEV_P-1) -= phi_i_ref(NUM_LEV_P-1);
        });
      }
#endif
    });

    // Laplacian of layer thickness
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1),
                   Homme::subview(m_buffers.dptens,kv.ie));
    // Laplacian of theta
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1),
                   Homme::subview(m_buffers.ttens,kv.ie));

    if (m_process_nh_vars) {
      // Laplacian of vertical velocity (do not compute last interface)
      m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                     Homme::subview(m_state.m_w_i,kv.ie,m_data.np1),
                     Homme::subview(m_buffers.wtens,kv.ie));
      // Laplacian of geopotential (do not compute last interface)
      m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                     Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1),
                     Homme::subview(m_buffers.phitens,kv.ie));
    }

    // Laplacian of velocity
    m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio1,
                              Homme::subview(m_state.m_v,kv.ie,m_data.np1),
                              Homme::subview(m_buffers.vtens,kv.ie));
  }


  //second iter of laplace, const hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagSecondLaplaceConstHV&, const TeamMember& team) const {
    KernelVariables kv(team, m_tu);
    // Laplacian of layers thickness
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_buffers.dptens,kv.ie),
                   Homme::subview(m_buffers.dptens,kv.ie));
    // Laplacian of theta
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_buffers.ttens,kv.ie),
                   Homme::subview(m_buffers.ttens,kv.ie));

    if (m_process_nh_vars) {
      // Laplacian of vertical velocity
      m_sphere_ops.laplace_simple(kv,
                     Homme::subview(m_buffers.wtens,kv.ie),
                     Homme::subview(m_buffers.wtens,kv.ie));
      // Laplacian of vertical geopotential
      m_sphere_ops.laplace_simple(kv,
                     Homme::subview(m_buffers.phitens,kv.ie),
                     Homme::subview(m_buffers.phitens,kv.ie));
    }
    // Laplacian of velocity
    m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio2,
                              Homme::subview(m_buffers.vtens,kv.ie),
                              Homme::subview(m_buffers.vtens,kv.ie));
  }

  //second iter of laplace, tensor hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagSecondLaplaceTensorHV&, const TeamMember& team) const {
    KernelVariables kv(team, m_tu);
    // Laplacian of layers thickness
    m_sphere_ops.laplace_tensor(kv,
                   Homme::subview(m_geometry.m_tensorvisc,kv.ie),
                   Homme::subview(m_buffers.dptens,kv.ie),
                   Homme::subview(m_buffers.dptens,kv.ie));
    // Laplacian of theta
    m_sphere_ops.laplace_tensor(kv,
                   Homme::subview(m_geometry.m_tensorvisc,kv.ie),
                   Homme::subview(m_buffers.ttens,kv.ie),
                   Homme::subview(m_buffers.ttens,kv.ie));

    if (m_process_nh_vars) {
      // Laplacian of vertical velocity
      m_sphere_ops.laplace_tensor(kv,
                     Homme::subview(m_geometry.m_tensorvisc,kv.ie),
                     Homme::subview(m_buffers.wtens,kv.ie),
                     Homme::subview(m_buffers.wtens,kv.ie));
      // Laplacian of geopotential
      m_sphere_ops.laplace_tensor(kv,
                     Homme::subview(m_geometry.m_tensorvisc,kv.ie),
                     Homme::subview(m_buffers.phitens,kv.ie),
                     Homme::subview(m_buffers.phitens,kv.ie));
    }

    // Laplacian of velocity
    m_sphere_ops.vlaplace_sphere_wk_cartesian(kv, 
                   Homme::subview(m_geometry.m_tensorvisc,kv.ie),
                   Homme::subview(m_geometry.m_vec_sph2cart,kv.ie),
                   Homme::subview(m_buffers.vtens,kv.ie),
                   Homme::subview(m_buffers.vtens,kv.ie));
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagUpdateStates&, const TeamMember& team) const {
    KernelVariables kv(team, m_tu);

    using MidColumn = decltype(Homme::subview(m_buffers.wtens,0,0,0));
    using IntColumn = decltype(Homme::subview(m_state.m_w_i,0,0,0,0));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // Add Xtens quantities back to the states, except for vtheta
      auto u = Homme::subview(m_state.m_v,kv.ie,m_data.np1,0,igp,jgp);
      auto v = Homme::subview(m_state.m_v,kv.ie,m_data.np1,1,igp,jgp);
      auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1,igp,jgp);
      auto dp     = Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp);

      auto utens   = Homme::subview(m_buffers.vtens,kv.ie,0,igp,jgp);
      auto vtens   = Homme::subview(m_buffers.vtens,kv.ie,1,igp,jgp);
      auto ttens   = Homme::subview(m_buffers.ttens,kv.ie,igp,jgp);
      auto dptens  = Homme::subview(m_buffers.dptens,kv.ie,igp,jgp);
      const auto& rspheremp = m_geometry.m_rspheremp(kv.ie,igp,jgp);

      MidColumn wtens, phitens;
      IntColumn w, phi_i;

      if (m_process_nh_vars) {
        wtens   = Homme::subview(m_buffers.wtens,kv.ie,igp,jgp);
        phitens = Homme::subview(m_buffers.phitens,kv.ie,igp,jgp);
        w       = Homme::subview(m_state.m_w_i,kv.ie,m_data.np1,igp,jgp);
        phi_i   = Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1,igp,jgp);
      }

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        utens(ilev)   *= m_data.dt*rspheremp;
        vtens(ilev)   *= m_data.dt*rspheremp;
        ttens(ilev)   *= m_data.dt*rspheremp;
        dptens(ilev)  *= m_data.dt*rspheremp;

        u(ilev)      += utens(ilev);
        v(ilev)      += vtens(ilev);
        vtheta(ilev) += ttens(ilev);
        dp(ilev)     += dptens(ilev);
        if (m_process_nh_vars) {
          wtens(ilev)   *= m_data.dt*rspheremp;
          phitens(ilev) *= m_data.dt*rspheremp;

          w(ilev)      += wtens(ilev);
          phi_i(ilev)  += phitens(ilev);
        }
      });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagHyperPreExchange, const TeamMember &team) const {
    using IntColumn = decltype(Homme::subview(m_state.m_w_i,0,0,0,0));

    KernelVariables kv(team, m_tu);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;

      const auto dpdiss_ave = Homme::subview(m_derived.m_dpdiss_ave,kv.ie, igp, jgp);
      const auto dpdiss_bih = Homme::subview(m_derived.m_dpdiss_biharmonic,kv.ie, igp, jgp);
      const auto dp3d = Homme::subview(m_state.m_dp3d,kv.ie, m_data.np1, igp, jgp);
      const auto dp_ref = Homme::subview(m_state.m_ref_states.dp_ref,kv.ie, igp, jgp);
      const auto theta = Homme::subview(m_state.m_vtheta_dp,kv.ie, m_data.np1, igp, jgp);
      const auto theta_ref = Homme::subview(m_state.m_ref_states.theta_ref,kv.ie, igp, jgp);
      const auto dptens = Homme::subview(m_buffers.dptens,kv.ie, igp, jgp);

      IntColumn phi, phi_ref;

      if (m_process_nh_vars) {
        phi = Homme::subview(m_state.m_phinh_i,kv.ie, m_data.np1, igp, jgp);
        phi_ref = Homme::subview(m_state.m_ref_states.phi_i_ref,kv.ie, igp, jgp);
      }

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int &ilev) {

        dp3d(ilev) += dp_ref(ilev);
        theta(ilev) += theta_ref(ilev);
        if (m_process_nh_vars) {
          phi(ilev) += phi_ref(ilev);
        }
        if (m_data.nu_p>0) {
          dpdiss_ave(ilev) += m_data.eta_ave_w*dp3d(ilev) / m_data.hypervis_subcycle;
          dpdiss_bih(ilev) += m_data.eta_ave_w*dptens(ilev) / m_data.hypervis_subcycle;
        }
      });
#ifndef XX_NONBFB_COMING
      // It would be fine to not even bother with the surface level, since
      // phitens is only NUM_LEV long, so all the hv stuff does not even happen
      // at NUM_LEV_P (unless NUM_LEV_P==NUM_LEV). However, removing the subtraction
      // and addition of phi_i_ref at NUM_LEV_P introduces NON BFB diffs.
      if (m_process_nh_vars && NUM_LEV!=NUM_LEV_P) {
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          phi(NUM_LEV_P-1) += phi_ref(NUM_LEV_P-1);
        });
      }
#endif
    });
    kv.team_barrier();

    // laplace subfunctors cannot be called from a TeamThreadRange or
    // ThreadVectorRange
    if (m_data.nu_top > 0) {

      //for top 3 levels and laplace, there is trivial nu_ratio only

      m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV,NUM_LEV>(
            kv, Homme::subview(m_state.m_dp3d, kv.ie, m_data.np1),
                Homme::subview(m_buffers.lapl_dp, kv.team_idx));

      m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV,NUM_LEV>(
            kv, Homme::subview(m_state.m_vtheta_dp, kv.ie, m_data.np1),
                Homme::subview(m_buffers.lapl_theta, kv.team_idx));

      if (m_process_nh_vars) {
        m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV,NUM_LEV_P>(
              kv, Homme::subview(m_state.m_w_i, kv.ie, m_data.np1),
                  Homme::subview(m_buffers.lapl_w, kv.team_idx));

        m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV,NUM_LEV_P>(
              kv, Homme::subview(m_state.m_phinh_i, kv.ie, m_data.np1),
                  Homme::subview(m_buffers.lapl_phi, kv.team_idx));
      }

      m_sphere_ops.vlaplace_sphere_wk_contra<NUM_BIHARMONIC_LEV,NUM_LEV>(
            kv, 1.0, Homme::subview(m_state.m_v, kv.ie, m_data.np1),
                     Homme::subview(m_buffers.lapl_v, kv.team_idx));
    }//if nu_top>0
    kv.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int &lev) {
        m_buffers.vtens(kv.ie, 0, igp, jgp, lev) *= -m_data.nu;
        m_buffers.vtens(kv.ie, 1, igp, jgp, lev) *= -m_data.nu;
        m_buffers.ttens(kv.ie, igp, jgp, lev) *= -m_data.nu;
        m_buffers.dptens(kv.ie, igp, jgp, lev) *= -m_data.nu_p;
        if (m_process_nh_vars) {
          m_buffers.wtens(kv.ie, igp, jgp, lev) *= -m_data.nu;
          m_buffers.phitens(kv.ie, igp, jgp, lev) *= -m_data.nu_s;
        }
      });

      if (m_data.nu_top > 0) {
        // Need to loop on the actual number of levels, not packs, to avoid adding more sponge layers
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, int(NUM_BIHARMONIC_PHYSICAL_LEVELS)),
                           [&](const int k) {
          const int ilev = k / VECTOR_SIZE;
          const int ivec = k % VECTOR_SIZE;

          m_buffers.vtens(kv.ie, 0, igp, jgp, ilev)[ivec] +=
              m_nu_scale_top[ilev][ivec] *
              m_buffers.lapl_v(kv.team_idx, 0, igp, jgp, ilev)[ivec];

          m_buffers.vtens(kv.ie, 1, igp, jgp, ilev)[ivec] +=
              m_nu_scale_top[ilev][ivec] *
              m_buffers.lapl_v(kv.team_idx, 1, igp, jgp, ilev)[ivec];

          m_buffers.ttens(kv.ie, igp, jgp, ilev)[ivec] +=
              m_nu_scale_top[ilev][ivec] *
              m_buffers.lapl_theta(kv.team_idx, igp, jgp, ilev)[ivec];

          m_buffers.dptens(kv.ie, igp, jgp, ilev)[ivec] +=
              m_nu_scale_top[ilev][ivec] *
              m_buffers.lapl_dp(kv.team_idx, igp, jgp, ilev)[ivec];

          if (m_process_nh_vars) {
            m_buffers.wtens(kv.ie, igp, jgp, ilev)[ivec] +=
                m_nu_scale_top[ilev][ivec] *
                m_buffers.lapl_w(kv.team_idx, igp, jgp, ilev)[ivec];

            m_buffers.phitens(kv.ie, igp, jgp, ilev)[ivec] +=
                m_nu_scale_top[ilev][ivec] *
                m_buffers.lapl_phi(kv.team_idx, igp, jgp, ilev)[ivec];
          }
        });
      }
    });
  }

protected:

  HyperviscosityData    m_data;
  ElementsState         m_state;
  ElementsDerivedState  m_derived;
  ElementsGeometry      m_geometry;
  SphereOperators       m_sphere_ops;
  ElementOps            m_elem_ops;
  EquationOfState       m_eos;
  Buffers               m_buffers;
  HybridVCoord          m_hvcoord;

  bool m_process_nh_vars;

  // Policies
  Kokkos::TeamPolicy<ExecSpace,TagUpdateStates>     m_policy_update_states;
  Kokkos::TeamPolicy<ExecSpace,TagFirstLaplaceHV>   m_policy_first_laplace;
  Kokkos::TeamPolicy<ExecSpace,TagHyperPreExchange> m_policy_pre_exchange;

  TeamUtils<ExecSpace> m_tu; // If the policies only differ by tag, just need one tu

  std::shared_ptr<BoundaryExchange> m_be;

  ExecViewManaged<Scalar[NUM_LEV]> m_nu_scale_top;
};

} // namespace Homme

#endif // HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP
