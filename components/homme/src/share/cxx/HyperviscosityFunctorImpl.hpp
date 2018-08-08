/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP
#define HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP

#include "Elements.hpp"
#include "Derivative.hpp"
#include "SimulationParams.hpp"
#include "KernelVariables.hpp"
#include "SphereOperators.hpp"

#include <memory>

namespace Homme
{

class BoundaryExchange;

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

public:

  struct TagFirstLaplaceHV {};
  struct TagSecondLaplaceConstHV {};
  struct TagSecondLaplaceTensorHV {};
  struct TagUpdateStates {};
  struct TagApplyInvMass {};
  struct TagHyperPreExchange {};

  HyperviscosityFunctorImpl (const SimulationParams& params, const Elements& elements, const Derivative& deriv);

  void init_boundary_exchanges();

  void run (const int np1, const Real dt, const Real eta_ave_w);

  void biharmonic_wk_dp3d () const;

// first iter of laplace, const hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagFirstLaplaceHV&, const TeamMember& team) const {
    KernelVariables kv(team);
    // Laplacian of temperature
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_elements.m_t,kv.ie,m_data.np1),
                   Homme::subview(m_elements.buffers.ttens,kv.ie));
    // Laplacian of pressure
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_elements.m_dp3d,kv.ie,m_data.np1),
                   Homme::subview(m_elements.buffers.dptens,kv.ie));

    // Laplacian of velocity
    m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio1,
                              Homme::subview(m_elements.m_v,kv.ie,m_data.np1),
                              Homme::subview(m_elements.buffers.vtens,kv.ie));
  }


//second iter of laplace, const hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagSecondLaplaceConstHV&, const TeamMember& team) const {
    KernelVariables kv(team);
    // Laplacian of temperature
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_elements.buffers.ttens,kv.ie),
                   Homme::subview(m_elements.buffers.ttens,kv.ie));
    // Laplacian of pressure
    m_sphere_ops.laplace_simple(kv,
                   Homme::subview(m_elements.buffers.dptens,kv.ie),
                   Homme::subview(m_elements.buffers.dptens,kv.ie));

    // Laplacian of velocity
    m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio2,
                              Homme::subview(m_elements.buffers.vtens,kv.ie),
                              Homme::subview(m_elements.buffers.vtens,kv.ie));
  }

//second iter of laplace, tensor hv
  KOKKOS_INLINE_FUNCTION
  void operator() (const TagSecondLaplaceTensorHV&, const TeamMember& team) const {
    KernelVariables kv(team);
    // Laplacian of temperature
    m_sphere_ops.laplace_tensor(kv,
                   Homme::subview(m_elements.m_tensorvisc,kv.ie),
                   Homme::subview(m_elements.buffers.ttens,kv.ie),
                   Homme::subview(m_elements.buffers.ttens,kv.ie));
    // Laplacian of pressure
    m_sphere_ops.laplace_tensor(kv,
                   Homme::subview(m_elements.m_tensorvisc,kv.ie),
                   Homme::subview(m_elements.buffers.dptens,kv.ie),
                   Homme::subview(m_elements.buffers.dptens,kv.ie));
    // Laplacian of velocity
    m_sphere_ops.vlaplace_sphere_wk_cartesian(kv, 
                              Homme::subview(m_elements.m_tensorvisc,kv.ie),
                              Homme::subview(m_elements.m_vec_sph2cart,kv.ie),
                              Homme::subview(m_elements.buffers.vtens,kv.ie),
                              Homme::subview(m_elements.buffers.vtens,kv.ie));

  }


  KOKKOS_INLINE_FUNCTION
  void operator() (const TagUpdateStates&, const int idx) const {
    const int ie   =  idx / (NP*NP*NUM_LEV);
    const int igp  = (idx / (NP*NUM_LEV)) % NP;
    const int jgp  = (idx / NUM_LEV) % NP;
    const int ilev =  idx % NUM_LEV;

    // Apply inverse mass matrix
    m_elements.buffers.vtens(ie,0,igp,jgp,ilev) = (m_data.dt * m_elements.buffers.vtens(ie,0,igp,jgp,ilev) *
                                                   m_elements.m_rspheremp(ie,igp,jgp));
    m_elements.buffers.vtens(ie,1,igp,jgp,ilev) = (m_data.dt * m_elements.buffers.vtens(ie,1,igp,jgp,ilev) *
                                                   m_elements.m_rspheremp(ie,igp,jgp));
    m_elements.m_v(ie,m_data.np1,0,igp,jgp,ilev) += m_elements.buffers.vtens(ie,0,igp,jgp,ilev);
    m_elements.m_v(ie,m_data.np1,1,igp,jgp,ilev) += m_elements.buffers.vtens(ie,1,igp,jgp,ilev);

    m_elements.buffers.ttens(ie,igp,jgp,ilev) = (m_data.dt*m_elements.buffers.ttens(ie,igp,jgp,ilev) *
                                                 m_elements.m_rspheremp(ie,igp,jgp));
    const Scalar heating = m_elements.buffers.vtens(ie,0,igp,jgp,ilev)*m_elements.m_v(ie,m_data.np1,0,igp,jgp,ilev)
                         + m_elements.buffers.vtens(ie,1,igp,jgp,ilev)*m_elements.m_v(ie,m_data.np1,1,igp,jgp,ilev);
    m_elements.m_t(ie,m_data.np1,igp,jgp,ilev) =
      m_elements.m_t(ie,m_data.np1,igp,jgp,ilev) + m_elements.buffers.ttens(ie,igp,jgp,ilev) -
      heating/PhysicalConstants::cp;

    m_elements.m_dp3d(ie,m_data.np1,igp,jgp,ilev) = (m_elements.buffers.dptens(ie,igp,jgp,ilev) *
                                                     m_elements.m_rspheremp(ie,igp,jgp));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagHyperPreExchange, const TeamMember &team) const {
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int &lev) {

        m_elements.m_derived_dpdiss_ave(kv.ie, igp, jgp, lev) +=
            m_data.eta_ave_w *
            m_elements.m_dp3d(kv.ie, m_data.np1, igp, jgp, lev) /
            m_data.hypervis_subcycle;
        m_elements.m_derived_dpdiss_biharmonic(kv.ie, igp, jgp, lev) +=
            m_data.eta_ave_w * m_elements.buffers.dptens(kv.ie, igp, jgp, lev) /
            m_data.hypervis_subcycle;
      });
    });
    kv.team_barrier();

    // Alias these for more descriptive names
    auto &laplace_v = m_elements.buffers.div_buf;
    auto &laplace_t = m_elements.buffers.lapl_buf_1;
    auto &laplace_dp3d = m_elements.buffers.lapl_buf_2;
    // laplace subfunctors cannot be called from a TeamThreadRange or
    // ThreadVectorRange
    constexpr int NUM_BIHARMONIC_PHYSICAL_LEVELS = 3;
    constexpr int NUM_BIHARMONIC_LEV = (NUM_BIHARMONIC_PHYSICAL_LEVELS + VECTOR_SIZE - 1) / VECTOR_SIZE;

    if (m_data.nu_top > 0) {

      //for top 3 levels and laplace, there is trivial nu_ratio only
      m_sphere_ops.vlaplace_sphere_wk_contra<NUM_BIHARMONIC_LEV>(
            kv, 1.0,
            // input
            Homme::subview(m_elements.m_v, kv.ie, m_data.np1),
            // output
            Homme::subview(laplace_v, kv.ie));

      m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV>(
            kv,
            // input
            Homme::subview(m_elements.m_t, kv.ie, m_data.np1),
            // output
            Homme::subview(laplace_t, kv.ie));

      m_sphere_ops.laplace_simple<NUM_BIHARMONIC_LEV>(
            kv,
            // input
            Homme::subview(m_elements.m_dp3d, kv.ie, m_data.np1),
            // output
            Homme::subview(laplace_dp3d, kv.ie));
    }//if nu_top>0
    kv.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int &lev) {
        m_elements.buffers.vtens(kv.ie, 0, igp, jgp, lev) *= -m_data.nu;
        m_elements.buffers.vtens(kv.ie, 1, igp, jgp, lev) *= -m_data.nu;
        m_elements.buffers.ttens(kv.ie, igp, jgp, lev) *= -m_data.nu_s;
        m_elements.buffers.dptens(kv.ie, igp, jgp, lev) *= -m_data.nu_p;
      });

      if (m_data.nu_top > 0) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, int(NUM_BIHARMONIC_LEV)),
                           [&](const int ilev) {
          m_elements.buffers.vtens(kv.ie, 0, igp, jgp, ilev) +=
              m_nu_scale_top[ilev] *
              laplace_v(kv.ie, 0, igp, jgp, ilev);
          m_elements.buffers.vtens(kv.ie, 1, igp, jgp, ilev) +=
              m_nu_scale_top[ilev] *
              laplace_v(kv.ie, 1, igp, jgp, ilev);

          m_elements.buffers.ttens(kv.ie, igp, jgp, ilev) +=
              m_nu_scale_top[ilev] *
              laplace_t(kv.ie, igp, jgp, ilev);

          m_elements.buffers.dptens(kv.ie, igp, jgp, ilev) +=
              m_nu_scale_top[ilev] *
              laplace_dp3d(kv.ie, igp, jgp, ilev);
        });
      }
      // While for T and v we exchange the tendencies, for dp3d we exchange the updated state.
      // However, since the BE structure already has registerd the *tens quantities, we store
      // the updated state in dptens.
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int &lev) {
          m_elements.buffers.dptens(kv.ie, igp, jgp, lev) *= m_data.dt;
          m_elements.buffers.dptens(kv.ie, igp, jgp, lev) += m_elements.m_dp3d(kv.ie,m_data.np1,igp,jgp,lev)
                                                           * m_elements.m_spheremp(kv.ie,igp,jgp);
      });
    });
  }

private:

  Elements            m_elements;
  Derivative          m_deriv;
  HyperviscosityData  m_data;
  SphereOperators     m_sphere_ops;

  std::shared_ptr<BoundaryExchange> m_be;

  ExecViewManaged<Scalar[NUM_LEV]> m_nu_scale_top;
};

} // namespace Homme

#endif // HOMMEXX_HYPERVISCOSITY_FUNCTOR_IMPL_HPP
