/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CAAR_FUNCTOR_IMPL_HPP
#define HOMMEXX_CAAR_FUNCTOR_IMPL_HPP

#include "Types.hpp"
#include "Elements.hpp"
#include "Tracers.hpp"
#include "HybridVCoord.hpp"
#include "Derivative.hpp"
#include "KernelVariables.hpp"
#include "SphereOperators.hpp"

#include "mpi/BoundaryExchange.hpp"
#include "utilities/SubviewUtils.hpp"

#include "profiling.hpp"
#include "ErrorDefs.hpp"

#include <assert.h>
#include <type_traits>


namespace Homme {

struct CaarFunctorImpl {

  struct CaarData {
    CaarData (const int rsplit_in) : rsplit(rsplit_in) {}
    int       nm1;
    int       n0;
    int       np1;
    int       n0_qdp;

    Real      dt;
    Real      eta_ave_w;

    const int rsplit;

    bool      compute_diagnostics;
  };

  CaarData              m_data;
  const HybridVCoord    m_hvcoord;
  const Elements        m_elements;
  const Tracers         m_tracers;
  const Derivative      m_deriv;
  SphereOperators       m_sphere_ops;

  Kokkos::Array<std::shared_ptr<BoundaryExchange>, NUM_TIME_LEVELS> m_bes;

  CaarFunctorImpl(const Elements &elements, const Tracers &tracers,
                  const Derivative &derivative, const HybridVCoord &hvcoord,
                  const SphereOperators &sphere_ops, 
                  const int rsplit)
      : m_data(rsplit)
      , m_hvcoord(hvcoord)
      , m_elements(elements)
      , m_tracers(tracers)
      , m_deriv(derivative)
      , m_sphere_ops(sphere_ops) {
    // Nothing to be done here
  }

  void init_boundary_exchanges (const std::shared_ptr<BuffersManager>& bm_exchange) {
    for (int tl=0; tl<NUM_TIME_LEVELS; ++tl) {
      m_bes[tl] = std::make_shared<BoundaryExchange>();
      auto& be = *m_bes[tl];
      be.set_buffers_manager(bm_exchange);
      be.set_num_fields(0,0,4);
      be.register_field(m_elements.m_v,tl,2,0);
      be.register_field(m_elements.m_t,1,tl);
      be.register_field(m_elements.m_dp3d,1,tl);
      be.registration_completed();
    }
  }

  void set_n0_qdp (const int n0_qdp) { m_data.n0_qdp = n0_qdp; }

  void set_rk_stage_data (const int nm1, const int n0,   const int np1,
                          const Real dt, const Real eta_ave_w,
                          const bool compute_diagnostics) {
    m_data.nm1 = nm1;
    m_data.n0  = n0;
    m_data.np1 = np1;
    m_data.dt  = dt;

    m_data.eta_ave_w = eta_ave_w;
    m_data.compute_diagnostics = compute_diagnostics;
  }

  // Depends on PHI (after preq_hydrostatic), PECND
  // Modifies Ephi_grad
  // Computes \nabla (E + phi) + \nabla (P) * Rgas * T_v / P
  KOKKOS_INLINE_FUNCTION void compute_energy_grad(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), [&] (const int& ilev) {
        // pre-fill energy_grad with the pressure(_grad)-temperature part
        m_elements.buffers.energy_grad(kv.team_idx, 0, igp, jgp, ilev) =
            PhysicalConstants::Rgas *
            (m_elements.buffers.temperature_virt(kv.team_idx, igp, jgp, ilev) /
             m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev)) *
            m_elements.buffers.pressure_grad(kv.team_idx, 0, igp, jgp, ilev);

        m_elements.buffers.energy_grad(kv.team_idx, 1, igp, jgp, ilev) =
            PhysicalConstants::Rgas *
            (m_elements.buffers.temperature_virt(kv.team_idx, igp, jgp, ilev) /
             m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev)) *
            m_elements.buffers.pressure_grad(kv.team_idx, 1, igp, jgp, ilev);

        // Kinetic energy + PHI (geopotential energy) +
        Scalar k_energy =
            0.5 * (m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev) *
                       m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev) +
                   m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev) *
                       m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev));
        m_elements.buffers.ephi(kv.team_idx, igp, jgp, ilev) =
            k_energy + m_elements.m_phi(kv.ie, igp, jgp, ilev);
      });
    });
    kv.team_barrier();

    m_sphere_ops.gradient_sphere_update(kv,
        Homme::subview(m_elements.buffers.ephi, kv.team_idx),
        Homme::subview(m_elements.buffers.energy_grad, kv.team_idx));
  } // TESTED 1

  // Depends on pressure, PHI, U_current, V_current, METDET,
  // D, DINV, U, V, FCOR, SPHEREMP, T_v, ETA_DPDN
  KOKKOS_INLINE_FUNCTION void compute_phase_3(KernelVariables &kv) const {
    if (m_data.rsplit == 0) {
      // vertical Eulerian
      assign_zero_to_sdot_sum(kv);
      compute_eta_dot_dpdn_vertadv_euler(kv);
      preq_vertadv(kv);
      accumulate_eta_dot_dpdn(kv);
    }
    compute_omega_p(kv);
    compute_temperature_np1(kv);
    compute_velocity_np1(kv);
    compute_dp3d_np1(kv);
  } // TRIVIAL
  //is it?

  KOKKOS_INLINE_FUNCTION
  void print_debug(KernelVariables &kv, const int ie, const int which) const {
    if( kv.ie == ie ){
      for(int k = 0; k < NUM_PHYSICAL_LEV; ++k){
        const int ilev = k / VECTOR_SIZE;
        const int ivec = k % VECTOR_SIZE;
        int igp = 0, jgp = 0;
        Real val;
        if( which == 0)
          val = m_elements.m_t(ie, m_data.np1, igp, jgp, ilev)[ivec];
        if( which == 1)
          val = m_elements.m_v(ie, m_data.np1, 0, igp, jgp, ilev)[ivec];
        if( which == 2)
          val = m_elements.m_v(ie, m_data.np1, 1, igp, jgp, ilev)[ivec];
        if( which == 3)
          val = m_elements.m_dp3d(ie, m_data.np1, igp, jgp, ilev)[ivec];
        Kokkos::single(Kokkos::PerTeam(kv.team), [&] () {
            if( which == 0)
              std::printf("m_t %d (%d %d): % .17e \n", k, ilev, ivec, val);
            if( which == 1)
              std::printf("m_v(0) %d (%d %d): % .17e \n", k, ilev, ivec, val);
            if( which == 2)
              std::printf("m_v(1) %d (%d %d): % .17e \n", k, ilev, ivec, val);
            if( which == 3)
              std::printf("m_dp3d %d (%d %d): % .17e \n", k, ilev, ivec, val);
          });
      }
    }
  }

  // Depends on pressure, PHI, U_current, V_current, METDET,
  // D, DINV, U, V, FCOR, SPHEREMP, T_v
  KOKKOS_INLINE_FUNCTION
  void compute_velocity_np1(KernelVariables &kv) const {
    compute_energy_grad(kv);

    m_sphere_ops.vorticity_sphere(kv,
        Homme::subview(m_elements.m_v, kv.ie, m_data.n0),
        Homme::subview(m_elements.buffers.vorticity, kv.team_idx));

    const bool rsplit_gt0 = m_data.rsplit > 0;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), [&] (const int& ilev) {
        // Recycle vort to contain (fcor+vort)
        m_elements.buffers.vorticity(kv.team_idx, igp, jgp, ilev) +=
            m_elements.m_fcor(kv.ie, igp, jgp);

        m_elements.buffers.energy_grad(kv.team_idx, 0, igp, jgp, ilev) *= -1;

        m_elements.buffers.energy_grad(kv.team_idx, 0, igp, jgp, ilev) +=
            (rsplit_gt0 ? 0 : - m_elements.buffers.v_vadv_buf(kv.team_idx, 0, igp, jgp, ilev)) +
            m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev) *
            m_elements.buffers.vorticity(kv.team_idx, igp, jgp, ilev);

        m_elements.buffers.energy_grad(kv.team_idx, 1, igp, jgp, ilev) *= -1;

        m_elements.buffers.energy_grad(kv.team_idx, 1, igp, jgp, ilev) +=
            (rsplit_gt0 ? 0 : - m_elements.buffers.v_vadv_buf(kv.team_idx, 1, igp, jgp, ilev)) -
            m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev) *
            m_elements.buffers.vorticity(kv.team_idx, igp, jgp, ilev);

        m_elements.buffers.energy_grad(kv.team_idx, 0, igp, jgp, ilev) *= m_data.dt;
        m_elements.buffers.energy_grad(kv.team_idx, 0, igp, jgp, ilev) +=
            m_elements.m_v(kv.ie, m_data.nm1, 0, igp, jgp, ilev);
        m_elements.buffers.energy_grad(kv.team_idx, 1, igp, jgp, ilev) *= m_data.dt;
        m_elements.buffers.energy_grad(kv.team_idx, 1, igp, jgp, ilev) +=
            m_elements.m_v(kv.ie, m_data.nm1, 1, igp, jgp, ilev);

        // Velocity at np1 = spheremp * buffer
        m_elements.m_v(kv.ie, m_data.np1, 0, igp, jgp, ilev) =
            m_elements.m_spheremp(kv.ie, igp, jgp) *
            m_elements.buffers.energy_grad(kv.team_idx, 0, igp, jgp, ilev);
        m_elements.m_v(kv.ie, m_data.np1, 1, igp, jgp, ilev) =
            m_elements.m_spheremp(kv.ie, igp, jgp) *
            m_elements.buffers.energy_grad(kv.team_idx, 1, igp, jgp, ilev);
      });
    });
    kv.team_barrier();
  } // UNTESTED 2
  //og: i'd better make a test for this

  //m_eta is zeroed outside of local kernels, in prim_step
  KOKKOS_INLINE_FUNCTION
  void accumulate_eta_dot_dpdn(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         KOKKOS_LAMBDA(const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      for(int k = 0; k < NUM_LEV; ++k){
        m_elements.m_eta_dot_dpdn(kv.ie, igp, jgp, k) +=
           m_data.eta_ave_w * m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, k);
      }//k loop
    });
    kv.team_barrier();
  } //tested against caar_adjust_eta_dot_dpdn_c_int


  KOKKOS_INLINE_FUNCTION
  void assign_zero_to_sdot_sum(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         KOKKOS_LAMBDA(const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      m_elements.buffers.sdot_sum(kv.team_idx, igp, jgp) = 0;
    });
    kv.team_barrier();
  } // TRIVIAL

  KOKKOS_INLINE_FUNCTION
  void compute_eta_dot_dpdn_vertadv_euler(KernelVariables &kv) const {

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         KOKKOS_LAMBDA(const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      for(int k = 0; k < NUM_PHYSICAL_LEV-1; ++k){
        const int ilev = k / VECTOR_SIZE;
        const int ivec = k % VECTOR_SIZE;
        const int kp1 = k+1;
        const int ilevp1 = kp1 / VECTOR_SIZE;
        const int ivecp1 = kp1 % VECTOR_SIZE;
        m_elements.buffers.sdot_sum(kv.team_idx, igp, jgp) +=
           m_elements.buffers.div_vdp(kv.team_idx, igp, jgp, ilev)[ivec];
        m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilevp1)[ivecp1] =
           m_elements.buffers.sdot_sum(kv.team_idx, igp, jgp);
      }//k loop
      //one more entry for sdot, separately, cause eta_dot_ is not of size LEV+1
      {
        const int ilev = (NUM_PHYSICAL_LEV - 1) / VECTOR_SIZE;
        const int ivec = (NUM_PHYSICAL_LEV - 1) % VECTOR_SIZE;
        m_elements.buffers.sdot_sum(kv.team_idx, igp, jgp) +=
           m_elements.buffers.div_vdp(kv.team_idx, igp, jgp, ilev)[ivec];
      }

      //note that index starts from 1
      for(int k = 1; k < NUM_PHYSICAL_LEV; ++k){
        const int ilev = k / VECTOR_SIZE;
        const int ivec = k % VECTOR_SIZE;
        m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilev)[ivec] =
           m_hvcoord.hybrid_bi(k)*m_elements.buffers.sdot_sum(kv.team_idx, igp, jgp) -
           m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilev)[ivec];
      }//k loop
      m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, 0)[0] = 0.0;
    });//NP*NP loop
    kv.team_barrier();
  }//TESTED against compute_eta_dot_dpdn_vertadv_euler_c_int

  // Depends on PHIS, DP3D, PHI, pressure, T_v
  // Modifies PHI
  KOKKOS_INLINE_FUNCTION
  void preq_hydrostatic(KernelVariables &kv) const {
    preq_hydrostatic_impl<ExecSpace>(kv);
  } // TESTED 3

  // Depends on pressure, U_current, V_current, div_vdp,
  // omega_p
  KOKKOS_INLINE_FUNCTION
  void preq_omega_ps(KernelVariables &kv) const {
    preq_omega_ps_impl<ExecSpace>(kv);
  } // TESTED 4

  // Depends on DP3D
  KOKKOS_INLINE_FUNCTION
  void compute_pressure(KernelVariables &kv) const {
    compute_pressure_impl<ExecSpace>(kv);
  } // TESTED 5

  // Depends on DP3D, PHIS, DP3D, PHI, T_v
  // Modifies pressure, PHI
  KOKKOS_INLINE_FUNCTION
  void compute_scan_properties(KernelVariables &kv) const {
    compute_pressure(kv);
    preq_hydrostatic(kv);
    preq_omega_ps(kv);
  } // TRIVIAL

//should be renamed, instead of no tracer should be dry
  KOKKOS_INLINE_FUNCTION
  void compute_temperature_no_tracers_helper(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), [&] (const int& ilev) {
        m_elements.buffers.temperature_virt(kv.team_idx, igp, jgp, ilev) =
            m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilev);
      });
    });
    kv.team_barrier();
  } // TESTED 6

//should be renamed
  KOKKOS_INLINE_FUNCTION
  void compute_temperature_tracers_helper(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), [&] (const int& ilev) {
        Scalar Qt = m_tracers.qdp(kv.ie, m_data.n0_qdp, 0, igp, jgp, ilev) /
                    m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev);
        Qt *= (PhysicalConstants::Rwater_vapor / PhysicalConstants::Rgas - 1.0);
        Qt += 1.0;
        m_elements.buffers.temperature_virt(kv.team_idx, igp, jgp, ilev) =
            m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilev) * Qt;
      });
    });
    kv.team_barrier();
  } // TESTED 7

  // Depends on DERIVED_UN0, DERIVED_VN0, METDET, DINV
  // Initializes div_vdp, which is used 2 times afterwards
  // Modifies DERIVED_UN0, DERIVED_VN0
  // Requires NUM_LEV * 5 * NP * NP
  KOKKOS_INLINE_FUNCTION
  void compute_div_vdp(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), [&] (const int& ilev) {
        m_elements.buffers.vdp(kv.team_idx, 0, igp, jgp, ilev) =
            m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev) *
            m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev);

        m_elements.buffers.vdp(kv.team_idx, 1, igp, jgp, ilev) =
            m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev) *
            m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev);

        m_elements.m_derived_vn0(kv.ie, 0, igp, jgp, ilev) +=
            m_data.eta_ave_w * m_elements.buffers.vdp(kv.team_idx, 0, igp, jgp, ilev);

        m_elements.m_derived_vn0(kv.ie, 1, igp, jgp, ilev) +=
            m_data.eta_ave_w * m_elements.buffers.vdp(kv.team_idx, 1, igp, jgp, ilev);
      });
    });
    kv.team_barrier();

    m_sphere_ops.divergence_sphere(kv,
        Homme::subview(m_elements.buffers.vdp, kv.team_idx),
        Homme::subview(m_elements.buffers.div_vdp, kv.team_idx));
  } // TESTED 8

  // Depends on T_current, DERIVE_UN0, DERIVED_VN0, METDET,
  // DINV
  // Might depend on QDP, DP3D_current
  KOKKOS_INLINE_FUNCTION
  void compute_temperature_div_vdp(KernelVariables &kv) const {
    if (m_data.n0_qdp < 0) {
      compute_temperature_no_tracers_helper(kv);
    } else {
      compute_temperature_tracers_helper(kv);
    }
    compute_div_vdp(kv);
  } // TESTED 9

  KOKKOS_INLINE_FUNCTION
  void compute_omega_p(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), [&] (const int& ilev) {
        m_elements.m_omega_p(kv.ie, igp, jgp, ilev) +=
            m_data.eta_ave_w * m_elements.buffers.omega_p(kv.team_idx, igp, jgp, ilev);
      });
    });
    kv.team_barrier();
  } // TESTED 10

  // Depends on T (global), OMEGA_P (global), U (global), V
  // (global),
  // SPHEREMP (global), T_v, and omega_p
  // block_3d_scalars
  KOKKOS_INLINE_FUNCTION
  void compute_temperature_np1(KernelVariables &kv) const {

    m_sphere_ops.gradient_sphere(kv,
        Homme::subview(m_elements.m_t, kv.ie, m_data.n0),
        Homme::subview(m_elements.buffers.temperature_grad, kv.team_idx));

    const bool rsplit_gt0 = m_data.rsplit > 0;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV), [&] (const int& ilev) {
        const Scalar vgrad_t =
            m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev) *
                m_elements.buffers.temperature_grad(kv.team_idx, 0, igp, jgp, ilev) +
            m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev) *
                m_elements.buffers.temperature_grad(kv.team_idx, 1, igp, jgp, ilev);

        const Scalar ttens =
              (rsplit_gt0 ? 0 : -m_elements.buffers.t_vadv_buf(kv.team_idx, igp, jgp,
                                                             ilev)) -
            vgrad_t +
            PhysicalConstants::kappa *
                m_elements.buffers.temperature_virt(kv.team_idx, igp, jgp, ilev) *
                m_elements.buffers.omega_p(kv.team_idx, igp, jgp, ilev);
        Scalar temp_np1 = ttens * m_data.dt +
                          m_elements.m_t(kv.ie, m_data.nm1, igp, jgp, ilev);
        temp_np1 *= m_elements.m_spheremp(kv.ie, igp, jgp);
        m_elements.m_t(kv.ie, m_data.np1, igp, jgp, ilev) = temp_np1;
      });
    });
    kv.team_barrier();
  } // TESTED 11

  // Depends on DERIVED_UN0, DERIVED_VN0, U, V,
  // Modifies DERIVED_UN0, DERIVED_VN0, OMEGA_P, T, and DP3D
  KOKKOS_INLINE_FUNCTION
  void compute_dp3d_np1(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int &ilev) {
        Scalar tmp = m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilev);
        tmp.shift_left(1);
        tmp[VECTOR_SIZE - 1] = (ilev + 1 < NUM_LEV)
                                   ? m_elements.buffers.eta_dot_dpdn_buf(
                                         kv.team_idx, igp, jgp, ilev + 1)[0]
                                   : 0;
        // Add div_vdp before subtracting the previous value to eta_dot_dpdn
        // This will hopefully reduce numeric error
        tmp += m_elements.buffers.div_vdp(kv.team_idx, igp, jgp, ilev);
        tmp -= m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilev);
        tmp = m_elements.m_dp3d(kv.ie, m_data.nm1, igp, jgp, ilev) -
              tmp * m_data.dt;

        m_elements.m_dp3d(kv.ie, m_data.np1, igp, jgp, ilev) =
            m_elements.m_spheremp(kv.ie, igp, jgp) * tmp;
      });
    });
    kv.team_barrier();
  } // TESTED 12

  // depends on eta_dot_dpdn, dp3d, T, v, modifies v_vadv, t_vadv
  KOKKOS_INLINE_FUNCTION
  void preq_vertadv(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         KOKKOS_LAMBDA(const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // first level
      int k = 0;
      int ilev = k / VECTOR_SIZE;
      int ivec = k % VECTOR_SIZE;
      const int kp1 = k + 1;
      const int ilevp1 = kp1 / VECTOR_SIZE;
      const int ivecp1 = kp1 % VECTOR_SIZE;

      // lets do this 1/dp thing to make it bfb with F and follow F for extra
      // (), not clear why
      Real facp =
          (0.5 * 1 /
           m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev)[ivec]) *
          m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilevp1)[ivecp1];
      Real facm;
      m_elements.buffers.t_vadv_buf(kv.team_idx, igp, jgp, ilev)[ivec] =
          facp * (m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilevp1)[ivecp1] -
                  m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilev)[ivec]);
      m_elements.buffers.v_vadv_buf(kv.team_idx, 0, igp, jgp, ilev)[ivec] =
          facp *
          (m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilevp1)[ivecp1] -
           m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev)[ivec]);
      m_elements.buffers.v_vadv_buf(kv.team_idx, 1, igp, jgp, ilev)[ivec] =
          facp *
          (m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilevp1)[ivecp1] -
           m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev)[ivec]);

      for (int k = 1; k < NUM_PHYSICAL_LEV - 1; ++k) {
        const int ilev = k / VECTOR_SIZE;
        const int ivec = k % VECTOR_SIZE;
        const int km1 = k - 1;
        const int ilevm1 = km1 / VECTOR_SIZE;
        const int ivecm1 = km1 % VECTOR_SIZE;
        const int kp1 = k + 1;
        const int ilevp1 = kp1 / VECTOR_SIZE;
        const int ivecp1 = kp1 % VECTOR_SIZE;

        facp = 0.5 *
               (1 / m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev)[ivec]) *
               m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp,
                                                   ilevp1)[ivecp1];
        facm = 0.5 *
               (1 / m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev)[ivec]) *
               m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilev)[ivec];

        m_elements.buffers.t_vadv_buf(kv.team_idx, igp, jgp, ilev)[ivec] =
            facp * (m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilevp1)[ivecp1] -
                    m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilev)[ivec]) +
            facm * (m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilev)[ivec] -
                    m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilevm1)[ivecm1]);

        m_elements.buffers.v_vadv_buf(kv.team_idx, 0, igp, jgp, ilev)[ivec] =
            facp *
                (m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilevp1)[ivecp1] -
                 m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev)[ivec]) +
            facm *
                (m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev)[ivec] -
                 m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilevm1)[ivecm1]);

        m_elements.buffers.v_vadv_buf(kv.team_idx, 1, igp, jgp, ilev)[ivec] =
            facp *
                (m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilevp1)[ivecp1] -
                 m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev)[ivec]) +
            facm *
                (m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev)[ivec] -
                 m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilevm1)[ivecm1]);
      } // k loop

      k = NUM_PHYSICAL_LEV - 1;
      ilev = k / VECTOR_SIZE;
      ivec = k % VECTOR_SIZE;
      const int km1 = k - 1;
      const int ilevm1 = km1 / VECTOR_SIZE;
      const int ivecm1 = km1 % VECTOR_SIZE;
      // note the (), just to comply with F
      facm = (0.5 *
              (1 / m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev)[ivec])) *
             m_elements.buffers.eta_dot_dpdn_buf(kv.team_idx, igp, jgp, ilev)[ivec];

      m_elements.buffers.t_vadv_buf(kv.team_idx, igp, jgp, ilev)[ivec] =
          facm * (m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilev)[ivec] -
                  m_elements.m_t(kv.ie, m_data.n0, igp, jgp, ilevm1)[ivecm1]);

      m_elements.buffers.v_vadv_buf(kv.team_idx, 0, igp, jgp, ilev)[ivec] =
          facm *
          (m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev)[ivec] -
           m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilevm1)[ivecm1]);

      m_elements.buffers.v_vadv_buf(kv.team_idx, 1, igp, jgp, ilev)[ivec] =
          facm *
          (m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev)[ivec] -
           m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilevm1)[ivecm1]);
    }); // NP*NP
    kv.team_barrier();
  } // TESTED against preq_vertadv

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamMember &team) const {
    KernelVariables kv(team);

    compute_temperature_div_vdp(kv);
    kv.team.team_barrier();

    compute_scan_properties(kv);
    kv.team.team_barrier();

    compute_phase_3(kv);
  }

  KOKKOS_INLINE_FUNCTION
  size_t shmem_size(const int team_size) const {
    return KernelVariables::shmem_size(team_size);
  }

private:
  template <typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<
      !std::is_same<ExecSpaceType, Hommexx_Cuda>::value, void>::type
  compute_pressure_impl(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;

      Real dp_prev = 0;
      Real p_prev = m_hvcoord.hybrid_ai0 * m_hvcoord.ps0;

      for (int ilev = 0; ilev < NUM_LEV; ++ilev) {
        const int vector_end =
            (ilev == NUM_LEV - 1
                 ? ((NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE)
                 : VECTOR_SIZE - 1);

        auto p = m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev);
        const auto &dp = m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev);

        for (int iv = 0; iv <= vector_end; ++iv) {
          // p[k] = p[k-1] + 0.5*(dp[k-1] + dp[k])
          p[iv] = p_prev + 0.5 * (dp_prev + dp[iv]);
          // Update p[k-1] and dp[k-1]
          p_prev = p[iv];
          dp_prev = dp[iv];
        }
        m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev) = p;
      };
    });
    kv.team_barrier();
  }

  template <typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<
      std::is_same<ExecSpaceType, Hommexx_Cuda>::value, void>::type
  compute_pressure_impl(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;

      ExecViewUnmanaged<Real[NUM_PHYSICAL_LEV]>        p (&m_elements.buffers.pressure(kv.team_idx,igp,jgp,0)[0]);
      ExecViewUnmanaged<const Real[NUM_PHYSICAL_LEV]> dp (&m_elements.m_dp3d(kv.team_idx,m_data.n0,igp,jgp,0)[0]);

      const Real p0 = m_hvcoord.hybrid_ai0 * m_hvcoord.ps0 + 0.5*dp(0);
      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        p(0) = p0;
      });
      // Compute cumsum of (dp(k-1) + dp(k)) in [1,NUM_PHYSICAL_LEV] and update p(k)
      // For level=1 (which is k=0), first add p0=hyai0*ps0 + 0.5*dp(0).
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, NUM_PHYSICAL_LEV-1,
                                            [&](const int k, Real& accumulator, const bool last) {

        // Note: the actual level must range in [1,NUM_PHYSICAL_LEV], while k ranges in [0, NUM_PHYSICAL_LEV-1].
        // If k=0, add the p0 part: hyai0*ps0 + 0.5*dp(0)
        accumulator += (k==0 ? p0 : 0);
        accumulator += (dp(k) + dp(k+1))/2;

        if (last) {
          p(k+1) = accumulator;
        }
      });
    });
    kv.team_barrier();
  }

  template <typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<
      !std::is_same<ExecSpaceType, Hommexx_Cuda>::value, void>::type
  preq_hydrostatic_impl(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;

      // Note: we add VECTOR_SIZE-1 rather than subtracting 1 since (0-1)%N=-1
      // while (0+N-1)%N=N-1.
      constexpr int last_lvl_last_vector_idx =
          (NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE;

      Real integration = 0;
      for (int ilev = NUM_LEV - 1; ilev >= 0; --ilev) {
        const int vec_start = (ilev == (NUM_LEV - 1) ? last_lvl_last_vector_idx
                                                     : VECTOR_SIZE - 1);

        const Real phis = m_elements.m_phis(kv.ie, igp, jgp);
        auto &phi = m_elements.m_phi(kv.ie, igp, jgp, ilev);
        const auto &t_v =
            m_elements.buffers.temperature_virt(kv.team_idx, igp, jgp, ilev);
        const auto &dp3d = m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev);
        const auto &p = m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev);

        // Precompute this product as a SIMD operation
        const auto rgas_tv_dp_over_p =
            PhysicalConstants::Rgas * t_v * (dp3d * 0.5 / p);

        // Integrate
        Scalar integration_ij;
        integration_ij[vec_start] = integration;
        for (int iv = vec_start - 1; iv >= 0; --iv)
          integration_ij[iv] =
              integration_ij[iv + 1] + rgas_tv_dp_over_p[iv + 1];

        // Add integral and constant terms to phi
        phi = phis + 2.0 * integration_ij + rgas_tv_dp_over_p;
        integration = integration_ij[0] + rgas_tv_dp_over_p[0];
      }
    });
    kv.team_barrier();
  }

  static KOKKOS_INLINE_FUNCTION void assert_vector_size_1() {
#ifndef NDEBUG
    if (VECTOR_SIZE != 1)
      Kokkos::abort("This impl is for GPU, for which VECTOR_SIZE is 1. It will "
                    "not work if VECTOR_SIZE > 1. Eventually, we may get "
                    "VECTOR_SIZE > 1 on GPU, at which point the alternative to "
                    "this impl will be the one to use, anyway.");
#endif
  }

  // CUDA version
  template <typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<
      std::is_same<ExecSpaceType, Hommexx_Cuda>::value, void>::type
  preq_hydrostatic_impl(KernelVariables &kv) const {
    assert_vector_size_1();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;

      // Use currently unused buffers to store one column of data.
      const auto rgas_tv_dp_over_p = Homme::subview(m_elements.buffers.vstar, kv.team_idx, 0, igp, jgp);
      const auto integration = Homme::subview(m_elements.buffers.divergence_temp,kv.team_idx,igp,jgp);

      // Precompute this product
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int &ilev) {
        const auto &t_v =
            m_elements.buffers.temperature_virt(kv.team_idx, igp, jgp, ilev);
        const auto &dp3d = m_elements.m_dp3d(kv.ie, m_data.n0, igp, jgp, ilev);
        const auto &p = m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev);

        rgas_tv_dp_over_p(ilev) =
            PhysicalConstants::Rgas * t_v * (dp3d * 0.5 / p);

        // Init integration to 0, so it's ready for later
        integration(ilev) = 0.0;
      });

      // accumulate rgas_tv_dp_over_p in [NUM_PHYSICAL_LEV,1].
      // Store sum(NUM_PHYSICAL_LEV,k) in integration(k-1)
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, NUM_PHYSICAL_LEV-1,
                                            [&](const int k, Real& accumulator, const bool last) {
        // level must range in [NUM_PHYSICAL_LEV-1,1], while k ranges in [0, NUM_PHYSICAL_LEV-2].
        const int level = NUM_PHYSICAL_LEV-1-k;

        accumulator += rgas_tv_dp_over_p(level)[0];

        if (last) {
          integration(level-1) = accumulator;
        }
      });

      // Update phi with a parallel_for: phi(k)= phis + 2.0 * integration(k+1) + rgas_tv_dp_over_p(k);
      const Real phis = m_elements.m_phis(kv.ie, igp, jgp);
      const auto phi = Homme::subview(m_elements.m_phi,kv.ie, igp, jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int level) {
        phi(level) = phis + 2.0 * integration(level) + rgas_tv_dp_over_p(level);
      });
    });
    kv.team_barrier();
  }

  // CUDA version
  template <typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<
      std::is_same<ExecSpaceType, Hommexx_Cuda>::value, void>::type
  preq_omega_ps_impl(KernelVariables &kv) const {
    assert_vector_size_1();
#ifdef DEBUG_TRACE
    Kokkos::single(Kokkos::PerTeam(kv.team), [&]() {
      m_elements.buffers.kernel_start_times(kv.ie) = clock();
    });
#endif
    m_sphere_ops.gradient_sphere(
        kv, Homme::subview(m_elements.buffers.pressure, kv.team_idx),
        Homme::subview(m_elements.buffers.pressure_grad, kv.team_idx));

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;

      ExecViewUnmanaged<Real[NUM_PHYSICAL_LEV]> omega_p(&m_elements.buffers.omega_p(kv.team_idx, igp, jgp, 0)[0]);
      ExecViewUnmanaged<Real[NUM_PHYSICAL_LEV]> div_vdp(&m_elements.buffers.div_vdp(kv.team_idx, igp, jgp, 0)[0]);
      Dispatch<ExecSpaceType>::parallel_scan(kv.team, NUM_PHYSICAL_LEV,
                                            [&](const int level, Real& accumulator, const bool last) {
        if (last) {
          omega_p(level) = accumulator;
        }

        accumulator += div_vdp(level);
      });

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                           [&](const int ilev) {
        const Scalar vgrad_p =
            m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev) *
                m_elements.buffers.pressure_grad(kv.team_idx, 0, igp, jgp, ilev) +
            m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev) *
                m_elements.buffers.pressure_grad(kv.team_idx, 1, igp, jgp, ilev);

        const auto &p = m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev);
        m_elements.buffers.omega_p(kv.team_idx, igp, jgp, ilev) =
            (vgrad_p -
             (m_elements.buffers.omega_p(kv.team_idx, igp, jgp, ilev) +
              0.5 * m_elements.buffers.div_vdp(kv.team_idx, igp, jgp, ilev))) /
            p;
      });
    });
#ifdef DEBUG_TRACE
    Kokkos::single(Kokkos::PerTeam(kv.team), [&]() {
      m_elements.buffers.kernel_end_times(kv.ie) = clock();
    });
#endif
  }

  // Non-CUDA version
  template <typename ExecSpaceType>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<
      !std::is_same<ExecSpaceType, Hommexx_Cuda>::value, void>::type
  preq_omega_ps_impl(KernelVariables &kv) const {
    m_sphere_ops.gradient_sphere(
        kv, Homme::subview(m_elements.buffers.pressure, kv.team_idx),
        Homme::subview(m_elements.buffers.pressure_grad, kv.team_idx));

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;

      Real integration = 0;
      for (int ilev = 0; ilev < NUM_LEV; ++ilev) {
        const int vector_end =
            (ilev == NUM_LEV - 1
                 ? ((NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE)
                 : VECTOR_SIZE - 1);

        const Scalar vgrad_p =
            m_elements.m_v(kv.ie, m_data.n0, 0, igp, jgp, ilev) *
                m_elements.buffers.pressure_grad(kv.team_idx, 0, igp, jgp, ilev) +
            m_elements.m_v(kv.ie, m_data.n0, 1, igp, jgp, ilev) *
                m_elements.buffers.pressure_grad(kv.team_idx, 1, igp, jgp, ilev);
        auto &omega_p = m_elements.buffers.omega_p(kv.team_idx, igp, jgp, ilev);
        const auto &p = m_elements.buffers.pressure(kv.team_idx, igp, jgp, ilev);
        const auto &div_vdp =
            m_elements.buffers.div_vdp(kv.team_idx, igp, jgp, ilev);

        Scalar integration_ij;
        integration_ij[0] = integration;
        for (int iv = 0; iv < vector_end; ++iv)
          integration_ij[iv + 1] = integration_ij[iv] + div_vdp[iv];
        omega_p = (vgrad_p - (integration_ij + 0.5 * div_vdp)) / p;
        integration = integration_ij[vector_end] + div_vdp[vector_end];
      }
    });
  }
};

} // Namespace Homme

#endif // HOMMEXX_CAAR_FUNCTOR_IMPL_HPP
