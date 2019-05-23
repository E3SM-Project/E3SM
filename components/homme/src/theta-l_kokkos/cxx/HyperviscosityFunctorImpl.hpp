/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP
#define HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP

#include "Elements.hpp"
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
  struct HyperviscosityData {
    HyperviscosityData(const int hypervis_subcycle_in, const Real nu_ratio1_in, const Real nu_ratio2_in, const Real nu_top_in,
                       const Real nu_in, const Real nu_p_in, const Real nu_s_in,
                       const Real hypervis_scaling_in)
                      : hypervis_subcycle(hypervis_subcycle_in), nu_ratio1(nu_ratio1_in), nu_ratio2(nu_ratio2_in)
                      , nu_top(nu_top_in), nu(nu_in), nu_p(nu_p_in), nu_s(nu_s_in)
                      , hypervis_scaling(hypervis_scaling_in)
                      , consthv(hypervis_scaling_in == 0){}


    const int   hypervis_subcycle;

    const Real  nu_ratio1;
    const Real  nu_ratio2;

    const Real  nu_top;
    const Real  nu;
    const Real  nu_p;
    const Real  nu_s;

    int         np1;
    Real        dt;

    Real        eta_ave_w;

    Real hypervis_scaling;
    bool consthv;
  };

  struct Buffers {
    static constexpr int num_2d_scalar_buf = 0;
    static constexpr int num_3d_scalar_mid_buf = 0;
    static constexpr int num_3d_vector_mid_buf = 0;
    static constexpr int num_3d_scalar_int_buf = 0;

    ExecViewUnmanaged<Real   * [NP][NP]> ps_ref;

    ExecViewUnmanaged<Scalar * [NP][NP][NUM_LEV]> dp_ref;
    ExecViewUnmanaged<Scalar * [NP][NP][NUM_LEV]> theta_ref;
    ExecViewUnmanaged<Scalar * [NP][NP][NUM_LEV]> thetadp_ref;

    ExecViewUnmanaged<Scalar * [NP][NP][NUM_LEV_P]> p_i;
    ExecViewUnmanaged<Scalar * [NP][NP][NUM_LEV_P]> phi_i_ref;
  };

  static constexpr Real T0 = 0.0065*288.0*PhysicalConstants::cp/PhysicalConstants::g;
  static constexpr Real T1 = 288.0-T0;

public:

  struct TagReferenceStates {};

  struct TagFirstLaplaceHV {};
  struct TagSecondLaplaceConstHV {};
  struct TagSecondLaplaceTensorHV {};
  struct TagUpdateStates {};
  struct TagApplyInvMass {};
  struct TagHyperPreExchange {};

  HyperviscosityFunctorImpl (const SimulationParams&       params,
                             const Elements&               elements);

  void request_buffers (FunctorsBuffersManager& fbm) const;
  void init_buffers (const FunctorsBuffersManager& fbm);
  void init_boundary_exchanges();

  void run (const int np1, const Real dt, const Real eta_ave_w);

  void biharmonic_wk_dp3d () const;

  // first iter of laplace, const hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagReferenceStates&, const TeamMember& team) const {
    KernelVariables kv(team);

    // Compute ps_ref and then dp_ref
    m_hvcoord.compute_ps_ref(kv,Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1),
                                Homme::subview(m_buffers.ps_ref,kv.team_idx));
    m_hvcoord.compute_dp_ref(kv,Homme::subview(m_buffers.ps_ref,kv.team_idx),
                                Homme::subview(m_buffers.dp_ref,kv.ie));

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // Compute theta_ref
      // Step 1: compute p at interfaces
      const auto& dp_real = viewAsReal(Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp));
      const auto& p_i = viewAsReal(Homme::subview(m_buffers.p_i,kv.team_idx,igp,jgp));
      p_i(0) = m_hvcoord.hybrid_ai0*m_hvcoord.ps0;
      Dispatch<>::parallel_scan(team,NUM_INTERFACE_LEV-1,
                                [&](const int ilev, Real& accumulator, const bool last) {
        accumulator += dp_real(ilev);
        if (last) {
          p_i(ilev+1) = accumulator;
        }
      });

      // Step 2: compute p at midpoints
      // Note: insteda of p_i above, re-subview the one in buffers, so that we can still operate on packs.
      const auto& theta_ref = Homme::subview(m_buffers.theta_ref,kv.ie,igp,jgp);
      const auto& thetadp_ref = Homme::subview(m_buffers.thetadp_ref,kv.ie,igp,jgp);
      const auto& dp_ref = Homme::subview(m_buffers.dp_ref,kv.ie,igp,jgp);
      m_elem_ops.compute_midpoint_values(kv,Homme::subview(m_buffers.p_i,kv.team_idx,igp,jgp),
                                            theta_ref);

      // Step 3-4: compute exner = (p/p0)^k, then theta_ref = T0/exner + T1 and thetadp_ref = theta_ref*dp_ref
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        theta_ref(ilev) /= m_hvcoord.ps0;
        pow_update(theta_ref(ilev),PhysicalConstants::kappa);

        theta_ref(ilev) = T0/theta_ref(ilev) + T1;
        thetadp_ref(ilev) = theta_ref(ilev)*dp_ref(ilev);
      });

      // Compute phi_i_ref = phis + int_top^k Rgas*vtheta_dp*(p/p0)^(k-1) / p0
      // TODO: you could exploit exner computed 5 lines up, to replace (p/p0)^(k-1) with exner*p0/p.
    // phi_i(:,:,nlevp) = phis(:,:)
    // do k=nlev,1,-1
    //    phi_i(:,:,k) = phi_i(:,:,k+1)+(Rgas*vtheta_dp(:,:,k)*(p(:,:,k)/p0)**(kappa-1))/p0
    // enddo
    });

    // // Laplacian of temperature
    // m_sphere_ops.laplace_simple(kv,
    //                Homme::subview(m_elements.m_state.m_t,kv.ie,m_data.np1),
    //                Homme::subview(m_elements.m_buffers.ttens,kv.ie));
    // // Laplacian of pressure
    // m_sphere_ops.laplace_simple(kv,
    //                Homme::subview(m_elements.m_state.m_dp3d,kv.ie,m_data.np1),
    //                Homme::subview(m_elements.m_buffers.dptens,kv.ie));

    // // Laplacian of velocity
    // m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio1,
    //                           Homme::subview(m_elements.m_state.m_v,kv.ie,m_data.np1),
    //                           Homme::subview(m_elements.m_buffers.vtens,kv.ie));
  }


  // first iter of laplace, const hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagFirstLaplaceHV&, const TeamMember& team) const {
    KernelVariables kv(team);
    // // Laplacian of temperature
    // m_sphere_ops.laplace_simple(kv,
    //                Homme::subview(m_elements.m_state.m_t,kv.ie,m_data.np1),
    //                Homme::subview(m_elements.m_buffers.ttens,kv.ie));
    // // Laplacian of pressure
    // m_sphere_ops.laplace_simple(kv,
    //                Homme::subview(m_elements.m_state.m_dp3d,kv.ie,m_data.np1),
    //                Homme::subview(m_elements.m_buffers.dptens,kv.ie));

    // // Laplacian of velocity
    // m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio1,
    //                           Homme::subview(m_elements.m_state.m_v,kv.ie,m_data.np1),
    //                           Homme::subview(m_elements.m_buffers.vtens,kv.ie));
  }


  //second iter of laplace, const hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagSecondLaplaceConstHV&, const TeamMember& team) const {
    KernelVariables kv(team);
    // // Laplacian of temperature
    // m_sphere_ops.laplace_simple(kv,
    //                Homme::subview(m_elements.m_buffers.ttens,kv.ie),
    //                Homme::subview(m_elements.m_buffers.ttens,kv.ie));
    // // Laplacian of pressure
    // m_sphere_ops.laplace_simple(kv,
    //                Homme::subview(m_elements.m_buffers.dptens,kv.ie),
    //                Homme::subview(m_elements.m_buffers.dptens,kv.ie));

    // // Laplacian of velocity
    // m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio2,
    //                           Homme::subview(m_elements.m_buffers.vtens,kv.ie),
    //                           Homme::subview(m_elements.m_buffers.vtens,kv.ie));
  }

  //second iter of laplace, tensor hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagSecondLaplaceTensorHV&, const TeamMember& team) const {
    KernelVariables kv(team);
    // // Laplacian of temperature
    // m_sphere_ops.laplace_tensor(kv,
    //                Homme::subview(m_elements.m_geometry.m_tensorvisc,kv.ie),
    //                Homme::subview(m_elements.m_buffers.ttens,kv.ie),
    //                Homme::subview(m_elements.m_buffers.ttens,kv.ie));
    // // Laplacian of pressure
    // m_sphere_ops.laplace_tensor(kv,
    //                Homme::subview(m_elements.m_geometry.m_tensorvisc,kv.ie),
    //                Homme::subview(m_elements.m_buffers.dptens,kv.ie),
    //                Homme::subview(m_elements.m_buffers.dptens,kv.ie));
    // // Laplacian of velocity
    // m_sphere_ops.vlaplace_sphere_wk_cartesian(kv, 
    //                           Homme::subview(m_elements.m_geometry.m_tensorvisc,kv.ie),
    //                           Homme::subview(m_elements.m_geometry.m_vec_sph2cart,kv.ie),
    //                           Homme::subview(m_elements.m_buffers.vtens,kv.ie),
    //                           Homme::subview(m_elements.m_buffers.vtens,kv.ie));
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const TagUpdateStates&, const int idx) const {
    // const int ie   =  idx / (NP*NP*NUM_LEV);
    // const int igp  = (idx / (NP*NUM_LEV)) % NP;
    // const int jgp  = (idx / NUM_LEV) % NP;
    // const int ilev =  idx % NUM_LEV;

    // // Apply inverse mass matrix
    // m_elements.m_buffers.vtens(ie,0,igp,jgp,ilev) = (m_data.dt * m_elements.m_buffers.vtens(ie,0,igp,jgp,ilev) *
    //                                                m_elements.m_geometry.m_rspheremp(ie,igp,jgp));
    // m_elements.m_buffers.vtens(ie,1,igp,jgp,ilev) = (m_data.dt * m_elements.m_buffers.vtens(ie,1,igp,jgp,ilev) *
    //                                                m_elements.m_geometry.m_rspheremp(ie,igp,jgp));
    // m_elements.m_state.m_v(ie,m_data.np1,0,igp,jgp,ilev) += m_elements.m_buffers.vtens(ie,0,igp,jgp,ilev);
    // m_elements.m_state.m_v(ie,m_data.np1,1,igp,jgp,ilev) += m_elements.m_buffers.vtens(ie,1,igp,jgp,ilev);

    // m_elements.m_buffers.ttens(ie,igp,jgp,ilev) = (m_data.dt*m_elements.m_buffers.ttens(ie,igp,jgp,ilev) *
    //                                              m_elements.m_geometry.m_rspheremp(ie,igp,jgp));
    // const Scalar heating = m_elements.m_buffers.vtens(ie,0,igp,jgp,ilev)*m_elements.m_state.m_v(ie,m_data.np1,0,igp,jgp,ilev)
    //                      + m_elements.m_buffers.vtens(ie,1,igp,jgp,ilev)*m_elements.m_state.m_v(ie,m_data.np1,1,igp,jgp,ilev);
    // m_elements.m_state.m_t(ie,m_data.np1,igp,jgp,ilev) =
    //   m_elements.m_state.m_t(ie,m_data.np1,igp,jgp,ilev) + m_elements.m_buffers.ttens(ie,igp,jgp,ilev) -
    //   heating/PhysicalConstants::cp;

    // m_elements.m_state.m_dp3d(ie,m_data.np1,igp,jgp,ilev) = (m_elements.m_buffers.dptens(ie,igp,jgp,ilev) *
    //                                                  m_elements.m_geometry.m_rspheremp(ie,igp,jgp));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagHyperPreExchange, const TeamMember &team) const {
    // KernelVariables kv(team);
    // Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
    //                      [&](const int &point_idx) {
    //   const int igp = point_idx / NP;
    //   const int jgp = point_idx % NP;
    //   Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
    //                        [&](const int &lev) {

    //     m_elements.m_derived.m_dpdiss_ave(kv.ie, igp, jgp, lev) +=
    //         m_data.eta_ave_w *
    //         m_elements.m_state.m_dp3d(kv.ie, m_data.np1, igp, jgp, lev) /
    //         m_data.hypervis_subcycle;
    //     m_elements.m_derived.m_dpdiss_biharmonic(kv.ie, igp, jgp, lev) +=
    //         m_data.eta_ave_w * m_elements.m_buffers.dptens(kv.ie, igp, jgp, lev) /
    //         m_data.hypervis_subcycle;
    //   });
    // });
    // kv.team_barrier();

    // // Alias these for more descriptive names
    // auto &laplace_v = m_elements.m_buffers.div_buf;
    // auto &laplace_t = m_elements.m_buffers.lapl_buf_1;
    // auto &laplace_dp3d = m_elements.m_buffers.lapl_buf_2;
    // // laplace subfunctors cannot be called from a TeamThreadRange or
    // // ThreadVectorRange
    // constexpr int NUM_BIHARMONIC_PHYSICAL_LEVELS = 3;
    // constexpr int NUM_BIHARMONIC_LEV = (NUM_BIHARMONIC_PHYSICAL_LEVELS + VECTOR_SIZE - 1) / VECTOR_SIZE;

    // if (m_data.nu_top > 0) {

    //   //for top 3 levels and laplace, there is trivial nu_ratio only
    //   m_sphere_ops.vlaplace_sphere_wk_contra<NUM_BIHARMONIC_LEV>(
    //         kv, 1.0,
    //         // input
    //         Homme::subview(m_elements.m_state.m_v, kv.ie, m_data.np1),
    //         // output
    //         Homme::subview(laplace_v, kv.ie));

    //   m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV>(
    //         kv,
    //         // input
    //         Homme::subview(m_elements.m_state.m_t, kv.ie, m_data.np1),
    //         // output
    //         Homme::subview(laplace_t, kv.ie));

    //   m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV>(
    //         kv,
    //         // input
    //         Homme::subview(m_elements.m_state.m_dp3d, kv.ie, m_data.np1),
    //         // output
    //         Homme::subview(laplace_dp3d, kv.ie));
    // }//if nu_top>0
    // kv.team_barrier();

    // Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
    //                      [&](const int &point_idx) {
    //   const int igp = point_idx / NP;
    //   const int jgp = point_idx % NP;
    //   Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
    //                        [&](const int &lev) {
    //     m_elements.m_buffers.vtens(kv.ie, 0, igp, jgp, lev) *= -m_data.nu;
    //     m_elements.m_buffers.vtens(kv.ie, 1, igp, jgp, lev) *= -m_data.nu;
    //     m_elements.m_buffers.ttens(kv.ie, igp, jgp, lev) *= -m_data.nu_s;
    //     m_elements.m_buffers.dptens(kv.ie, igp, jgp, lev) *= -m_data.nu_p;
    //   });

    //   if (m_data.nu_top > 0) {
    //     Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, int(NUM_BIHARMONIC_LEV)),
    //                        [&](const int ilev) {
    //       m_elements.m_buffers.vtens(kv.ie, 0, igp, jgp, ilev) +=
    //           m_nu_scale_top[ilev] *
    //           laplace_v(kv.ie, 0, igp, jgp, ilev);
    //       m_elements.m_buffers.vtens(kv.ie, 1, igp, jgp, ilev) +=
    //           m_nu_scale_top[ilev] *
    //           laplace_v(kv.ie, 1, igp, jgp, ilev);

    //       m_elements.m_buffers.ttens(kv.ie, igp, jgp, ilev) +=
    //           m_nu_scale_top[ilev] *
    //           laplace_t(kv.ie, igp, jgp, ilev);

    //       m_elements.m_buffers.dptens(kv.ie, igp, jgp, ilev) +=
    //           m_nu_scale_top[ilev] *
    //           laplace_dp3d(kv.ie, igp, jgp, ilev);
    //     });
    //   }
    //   // While for T and v we exchange the tendencies, for dp3d we exchange the updated state.
    //   // However, since the BE structure already has registerd the *tens quantities, we store
    //   // the updated state in dptens.
    //   Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
    //                        [&](const int &lev) {
    //       m_elements.m_buffers.dptens(kv.ie, igp, jgp, lev) *= m_data.dt;
    //       m_elements.m_buffers.dptens(kv.ie, igp, jgp, lev) += m_elements.m_state.m_dp3d(kv.ie,m_data.np1,igp,jgp,lev)
    //                                                          * m_elements.m_geometry.m_spheremp(kv.ie,igp,jgp);
    //   });
    // });
  }

private:

  HyperviscosityData  m_data;
  ElementsState       m_state;
  ElementsGeometry    m_geometry;
  SphereOperators     m_sphere_ops;
  ElementOps          m_elem_ops;
  Buffers             m_buffers;
  HybridVCoord        m_hvcoord;

  // Policies
  Kokkos::RangePolicy<ExecSpace,TagUpdateStates>    m_policy_update_states;
  Kokkos::TeamPolicy<ExecSpace,TagFirstLaplaceHV>   m_policy_first_laplace;
  Kokkos::TeamPolicy<ExecSpace,TagHyperPreExchange> m_policy_pre_exchange;

  std::shared_ptr<BoundaryExchange> m_be;

  ExecViewManaged<Scalar[NUM_LEV]> m_nu_scale_top;
};

} // namespace Homme

#endif // HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP
