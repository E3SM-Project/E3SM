#include "CaarFunctorImpl.hpp"

namespace Homme {

template <bool HYDROSTATIC, bool CONSERVATIVE>
void CaarFunctorImpl::epoch1_blockOps()
{
  auto buffers_dp_tens = viewAsReal(m_buffers.dp_tens);
  auto buffers_exner = viewAsReal(m_buffers.exner);
  auto buffers_phi = viewAsReal(m_buffers.phi);
  auto buffers_pnh = viewAsReal(m_buffers.pnh);
  auto buffers_theta_tens = viewAsReal(m_buffers.theta_tens);
  auto buffers_vdp = viewAsReal(m_buffers.vdp);

  const Real data_eta_ave_w = m_data.eta_ave_w;
  const int data_n0 = m_data.n0;

  auto derived_vn0 = viewAsReal(m_derived.m_vn0);

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = viewAsReal(m_state.m_dp3d);
  auto state_phinh_i = viewAsReal(m_state.m_phinh_i);
  auto state_v = viewAsReal(m_state.m_v);
  auto state_vtheta_dp= viewAsReal(m_state.m_vtheta_dp);

  Kokkos::parallel_for(
    __PRETTY_FUNCTION__,
    SphereBlockOps::policy(m_num_elems, 4),
    KOKKOS_LAMBDA(const Team &team) {

      SphereBlockOps b(sg, team);
      if (b.skip()) return;

      if (!HYDROSTATIC) {

        const Real dphi = state_phinh_i(b.e,data_n0,b.x,b.y,b.z+1) - state_phinh_i(b.e,data_n0,b.x,b.y,b.z);
        const Real vtheta_dp = state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        if ((vtheta_dp < 0) || (dphi > 0)) abort();
        EquationOfState::compute_pnh_and_exner(vtheta_dp, dphi, buffers_pnh(b.e,b.x,b.y,b.z), buffers_exner(b.e,b.x,b.y,b.z));

        buffers_phi(b.e,b.x,b.y,b.z) = 0.5 * (state_phinh_i(b.e,data_n0,b.x,b.y,b.z) + state_phinh_i(b.e,data_n0,b.x,b.y,b.z+1));
      }

      const Real v0 = state_v(b.e,data_n0,0,b.x,b.y,b.z) * state_dp3d(b.e,data_n0,b.x,b.y,b.z);
      const Real v1 = state_v(b.e,data_n0,1,b.x,b.y,b.z) * state_dp3d(b.e,data_n0,b.x,b.y,b.z);
      buffers_vdp(b.e,0,b.x,b.y,b.z) = v0;
      buffers_vdp(b.e,1,b.x,b.y,b.z) = v1;
      derived_vn0(b.e,0,b.x,b.y,b.z) += data_eta_ave_w * v0;
      derived_vn0(b.e,1,b.x,b.y,b.z) += data_eta_ave_w * v1;

      SphereBlockScratch ttmp0(b);
      SphereBlockScratch ttmp1(b);
      b.divInit(ttmp0, ttmp1, v0, v1);

      SphereBlockScratch ttmp2(b);
      SphereBlockScratch ttmp3(b);
      Real vtheta = 0;

      if (CONSERVATIVE) {
        const Real sv0 = state_v(b.e,data_n0,0,b.x,b.y,b.z) * state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        const Real sv1 = state_v(b.e,data_n0,1,b.x,b.y,b.z) * state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        b.divInit(ttmp2, ttmp3, sv0, sv1);
      } else {
        vtheta = state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z) / state_dp3d(b.e,data_n0,b.x,b.y,b.z);
        b.gradInit(ttmp2, vtheta);
      }

      b.barrier();

      const Real dvdp = b.div(ttmp0, ttmp1);
      buffers_dp_tens(b.e,b.x,b.y,b.z) = dvdp;

      if (CONSERVATIVE) {
        buffers_theta_tens(b.e,b.x,b.y,b.z) = b.div(ttmp2, ttmp3);
      } else {
        Real grad0, grad1;
        b.grad(grad0, grad1, ttmp2);
        Real theta_tens = dvdp * vtheta;
        theta_tens += grad0 * v0;
        theta_tens += grad1 * v1;
        buffers_theta_tens(b.e,b.x,b.y,b.z) = theta_tens;
      }
    });
}

template <bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch2_scanOps()
{
  auto buffers_dp_tens = viewAsReal(m_buffers.dp_tens);
  auto buffers_dp_i = viewAsReal(m_buffers.dp_i);
  auto buffers_eta_dot_dpdn = viewAsReal(m_buffers.eta_dot_dpdn);
  auto buffers_w_tens = viewAsReal(m_buffers.w_tens);

  const int data_n0 = m_data.n0;
  auto &hvcoord_hybrid_bi = m_hvcoord.hybrid_bi;
  const Real pi_i00 = m_hvcoord.ps0 * m_hvcoord.hybrid_ai0;
  auto state_dp3d = viewAsReal(m_state.m_dp3d);

  Kokkos::parallel_for(
    __PRETTY_FUNCTION__,
    SphereScanOps::policy(m_num_elems),
    KOKKOS_LAMBDA(const Team &team) {

      const SphereScanOps s(team);

      s.scan(buffers_dp_i, state_dp3d, data_n0, pi_i00);
      s.scan(buffers_w_tens, buffers_dp_tens, 0);

      if (RSPLIT_ZERO) {

        s.scan(buffers_eta_dot_dpdn, buffers_dp_tens, 0);

        const Real last = buffers_eta_dot_dpdn(s.e,s.x,s.y,NUM_PHYSICAL_LEV);

        Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(s.t, 1, NUM_PHYSICAL_LEV),
          [&](const int z) {
            Real eta_dot_dpdn = -buffers_eta_dot_dpdn(s.e,s.x,s.y,z);
            eta_dot_dpdn += hvcoord_hybrid_bi(z) * last;
            buffers_eta_dot_dpdn(s.e,s.x,s.y,z) = eta_dot_dpdn;
          });

        Kokkos::single(
          Kokkos::PerThread(team),
          [&]() {
            buffers_eta_dot_dpdn(s.e,s.x,s.y,0) = buffers_eta_dot_dpdn(s.e,s.x,s.y,NUM_PHYSICAL_LEV) = 0;
          });
      }
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch3_blockOps()
{
  auto buffers_dp_i = viewAsReal(m_buffers.dp_i);
  auto buffers_dp_tens = viewAsReal(m_buffers.dp_tens);
  auto buffers_eta_dot_dpdn = viewAsReal(m_buffers.eta_dot_dpdn);
  auto buffers_exner = viewAsReal(m_buffers.exner);
  auto buffers_omega_p = viewAsReal(m_buffers.omega_p);
  auto buffers_pnh = viewAsReal(m_buffers.pnh);
  auto buffers_temp = viewAsReal(m_buffers.temp);
  auto buffers_v_tens = viewAsReal(m_buffers.v_tens);
  auto buffers_w_tens = viewAsReal(m_buffers.w_tens);

  const Real data_eta_ave_w = m_data.eta_ave_w;
  const int data_n0 = m_data.n0;

  auto derived_eta_dot_dpdn = viewAsReal(m_derived.m_eta_dot_dpdn);
  auto derived_omega_p = viewAsReal(m_derived.m_omega_p);

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = viewAsReal(m_state.m_dp3d);
  auto state_v = viewAsReal(m_state.m_v);
  auto state_vtheta_dp= viewAsReal(m_state.m_vtheta_dp);
  auto state_w_i = viewAsReal(m_state.m_w_i);

  Kokkos::parallel_for(
    __PRETTY_FUNCTION__,
    SphereBlockOps::policy(m_num_elems, 1),
    KOKKOS_LAMBDA(const Team &team) {

      SphereBlockOps b(sg, team);
      if (b.skip()) return;

      const Real pi = 0.5 * (buffers_dp_i(b.e,b.x,b.y,b.z) + buffers_dp_i(b.e,b.x,b.y,b.z+1));
      SphereBlockScratch tmp0(b);
      b.gradInit(tmp0, pi);

      if (HYDROSTATIC) {
        Real exner = pi;
        EquationOfState::pressure_to_exner(exner);
        buffers_exner(b.e,b.x,b.y,b.z) = exner;
        buffers_pnh(b.e,b.x,b.y,b.z) = EquationOfState::compute_dphi(state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z), exner, pi);
      }

      b.barrier();

      Real grad0, grad1;
      b.grad(grad0, grad1, tmp0);

      Real omega = -0.5 * (buffers_w_tens(b.e,b.x,b.y,b.z) + buffers_w_tens(b.e,b.x,b.y,b.z+1));
      const Real uz = state_v(b.e,data_n0,0,b.x,b.y,b.z);
      const Real vz = state_v(b.e,data_n0,1,b.x,b.y,b.z);
      omega += uz * grad0 + vz * grad1;
      buffers_omega_p(b.e,b.x,b.y,b.z) = omega;

      derived_omega_p(b.e,b.x,b.y,b.z) += data_eta_ave_w * omega;

      if (RSPLIT_ZERO) {

        const Real dp = state_dp3d(b.e,data_n0,b.x,b.y,b.z);
        const Real etap = buffers_eta_dot_dpdn(b.e,b.x,b.y,b.z+1);
        const Real etaz = buffers_eta_dot_dpdn(b.e,b.x,b.y,b.z);

        buffers_dp_tens(b.e,b.x,b.y,b.z) += etap - etaz;

        derived_eta_dot_dpdn(b.e,b.x,b.y,b.z) += data_eta_ave_w * etaz;

        if (!HYDROSTATIC) {
          const Real dw = state_w_i(b.e,data_n0,b.x,b.y,b.z+1) - state_w_i(b.e,data_n0,b.x,b.y,b.z);
          const Real eta = 0.5 * (etaz + etap);
          buffers_temp(b.e,b.x,b.y,b.z) = dw * eta;
        }

        Real u = 0;
        Real v = 0;
        if (b.z < NUM_PHYSICAL_LEV-1) {
          const Real facp = 0.5 * etap / dp;
          u = facp * (state_v(b.e,data_n0,0,b.x,b.y,b.z+1) - uz);
          v = facp * (state_v(b.e,data_n0,1,b.x,b.y,b.z+1) - vz);
        }
        if (b.z > 0) {
          const Real facm = 0.5 * etaz / dp;
          u += facm * (uz - state_v(b.e,data_n0,0,b.x,b.y,b.z-1));
          v += facm * (vz - state_v(b.e,data_n0,1,b.x,b.y,b.z-1));
        }
        buffers_v_tens(b.e,0,b.x,b.y,b.z) = u;
        buffers_v_tens(b.e,1,b.x,b.y,b.z) = v;
      }
    });
}

void CaarFunctorImpl::epoch4_scanOps()
{
  auto buffers_phi = viewAsReal(m_buffers.phi);
  auto buffers_pnh = viewAsReal(m_buffers.pnh);

  const int data_n0 = m_data.n0;

  auto &geometry_phis = m_geometry.m_phis;

  auto state_phinh_i = viewAsReal(m_state.m_phinh_i);

  Kokkos::parallel_for(
    __PRETTY_FUNCTION__,
    SphereScanOps::policy(m_num_elems),
    KOKKOS_LAMBDA(const Team &team) {
      const SphereScanOps s(team);
      s.nacs(state_phinh_i, data_n0, buffers_pnh, geometry_phis);
      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(s.t, NUM_PHYSICAL_LEV),
        [&](const int z) {
          buffers_phi(s.e,s.x,s.y,z) = 0.5 * (state_phinh_i(s.e,data_n0,s.x,s.y,z) + state_phinh_i(s.e,data_n0,s.x,s.y,z+1));
        });
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch5_colOps()
{
  auto buffers_dp_i = viewAsReal(m_buffers.dp_i);
  auto buffers_dpnh_dp_i = viewAsReal(m_buffers.dpnh_dp_i);
  auto buffers_eta_dot_dpdn = viewAsReal(m_buffers.eta_dot_dpdn);
  auto buffers_exner = viewAsReal(m_buffers.exner);
  auto buffers_grad_phinh_i = viewAsReal(m_buffers.grad_phinh_i);
  auto buffers_grad_w_i = viewAsReal(m_buffers.grad_w_i);
  auto buffers_phi = viewAsReal(m_buffers.phi);
  auto buffers_phi_tens = viewAsReal(m_buffers.phi_tens);
  auto buffers_pnh = viewAsReal(m_buffers.pnh);
  auto buffers_temp = viewAsReal(m_buffers.temp);
  auto buffers_v_i = viewAsReal(m_buffers.v_i);
  auto buffers_vtheta_i = viewAsReal(m_buffers.vtheta_i);
  auto buffers_w_tens = viewAsReal(m_buffers.w_tens);

  const int data_n0 = m_data.n0;
  const Real dscale = m_data.scale1 - m_data.scale2;

  auto &geometry_gradphis = m_geometry.m_gradphis;

  const Real gscale1 = m_data.scale1 * PhysicalConstants::g;
  const Real gscale2 = m_data.scale2 * PhysicalConstants::g;

  auto hvcoord_hybrid_bi_packed = viewAsReal(m_hvcoord.hybrid_bi_packed);

  const Real ndata_scale1 = -m_data.scale1;
  const Real pi_i00 = m_hvcoord.ps0 * m_hvcoord.hybrid_ai0;

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = viewAsReal(m_state.m_dp3d);
  auto state_phinh_i = viewAsReal(m_state.m_phinh_i);
  auto state_v = viewAsReal(m_state.m_v);
  auto state_w_i = viewAsReal(m_state.m_w_i);

  Kokkos::parallel_for(
    __PRETTY_FUNCTION__,
    SphereColOps::policy(m_num_elems, NUM_INTERFACE_LEV),
    KOKKOS_LAMBDA(const Team &team) {

      const SphereColOps c(sg, team);

      c.grad(buffers_grad_phinh_i, state_phinh_i, data_n0);
      if (!HYDROSTATIC) c.grad(buffers_grad_w_i, state_w_i, data_n0); 

      const Real dm = (c.z == 0) ? 0 : state_dp3d(c.e,data_n0,c.x,c.y,c.z-1);
      const Real dz = (c.z == NUM_PHYSICAL_LEV) ? 0 : state_dp3d(c.e,data_n0,c.x,c.y,c.z);
      const Real dp_i = (c.z == 0) ? dz : (c.z == NUM_PHYSICAL_LEV) ? dm : 0.5 * (dz + dm);
      buffers_dp_i(c.e,c.x,c.y,c.z) = dp_i;

      if (!HYDROSTATIC) {

        const Real v0m = (c.z == 0) ? 0 : state_v(c.e,data_n0,0,c.x,c.y,c.z-1);
        const Real v0z = (c.z == NUM_PHYSICAL_LEV) ? 0 : state_v(c.e,data_n0,0,c.x,c.y,c.z);
        const Real v_i0 = (c.z == 0) ? v0z : (c.z == NUM_PHYSICAL_LEV) ? v0m : (dz * v0z + dm * v0m) / (dm + dz);
        buffers_v_i(c.e,0,c.x,c.y,c.z) = v_i0;

        const Real v1m = (c.z == 0) ? 0 : state_v(c.e,data_n0,1,c.x,c.y,c.z-1);
        const Real v1z = (c.z == NUM_PHYSICAL_LEV) ? 0 : state_v(c.e,data_n0,1,c.x,c.y,c.z);
        const Real v_i1 = (c.z == 0) ? v1z : (c.z == NUM_PHYSICAL_LEV) ? v1m :(dz * v1z + dm * v1m) / (dm + dz);
        buffers_v_i(c.e,1,c.x,c.y,c.z) = v_i1;

        const Real pm = (c.z == 0) ? pi_i00 : buffers_pnh(c.e,c.x,c.y,c.z-1);
        const Real pz = (c.z == NUM_PHYSICAL_LEV) ? pm + 0.5 * dm : buffers_pnh(c.e,c.x,c.y,c.z);
        buffers_dpnh_dp_i(c.e,c.x,c.y,c.z) = 2.0 * (pz - pm) / (dm + dz);
      }

      if (RSPLIT_ZERO) {

        const Real phim = (c.z == 0) ? 0 : buffers_phi(c.e,c.x,c.y,c.z-1);
        const Real phiz = (c.z == NUM_PHYSICAL_LEV) ? 0 : buffers_phi(c.e,c.x,c.y,c.z);

        if (!HYDROSTATIC) {
          const Real phi_vadv = ((c.z == 0) || (c.z == NUM_PHYSICAL_LEV)) ? 0 : (phiz - phim) * buffers_eta_dot_dpdn(c.e,c.x,c.y,c.z) / dp_i;
          buffers_phi_tens(c.e,c.x,c.y,c.z) = phi_vadv;
        }

        Real vtheta_i = phiz - phim;
        if (!(c.z == 0) && !(c.z == NUM_PHYSICAL_LEV)) {
          const Real dexner = buffers_exner(c.e,c.x,c.y,c.z) - buffers_exner(c.e,c.x,c.y,c.z-1);
          vtheta_i /= dexner;
        }
        vtheta_i /= -PhysicalConstants::cp;
        if (!HYDROSTATIC) vtheta_i *= buffers_dpnh_dp_i(c.e,c.x,c.y,c.z);
        buffers_vtheta_i(c.e,c.x,c.y,c.z) = vtheta_i;
      }

      if (!HYDROSTATIC) {

        Real w_tens = 0;
        if (RSPLIT_ZERO) {
          const Real tempm = (c.z == 0) ? 0 : buffers_temp(c.e,c.x,c.y,c.z-1);
          const Real tempz = (c.z == NUM_PHYSICAL_LEV) ? 0 : buffers_temp(c.e,c.x,c.y,c.z);
          const Real dw = (c.z == 0) ? tempz : (c.z == NUM_PHYSICAL_LEV) ? tempm : 0.5 * (tempz + tempm);
          w_tens = dw / dp_i;
        }
        w_tens += buffers_v_i(c.e,0,c.x,c.y,c.z) * buffers_grad_w_i(c.e,0,c.x,c.y,c.z) + buffers_v_i(c.e,1,c.x,c.y,c.z) * buffers_grad_w_i(c.e,1,c.x,c.y,c.z);
        w_tens *= ndata_scale1;
        const Real gscale = (c.z == NUM_PHYSICAL_LEV) ? gscale1 : gscale2;
        w_tens += (buffers_dpnh_dp_i(c.e,c.x,c.y,c.z)-Real(1)) * gscale;
        buffers_w_tens(c.e,c.x,c.y,c.z) = w_tens;

        Real phi_tens = (RSPLIT_ZERO) ? buffers_phi_tens(c.e,c.x,c.y,c.z) : 0;
        phi_tens += buffers_v_i(c.e,0,c.x,c.y,c.z) * buffers_grad_phinh_i(c.e,0,c.x,c.y,c.z) + buffers_v_i(c.e,1,c.x,c.y,c.z) * buffers_grad_phinh_i(c.e,1,c.x,c.y,c.z);
        phi_tens *= ndata_scale1;
        phi_tens += state_w_i(c.e,data_n0,c.x,c.y,c.z) * gscale;

        if (dscale) phi_tens += dscale * (buffers_v_i(c.e,0,c.x,c.y,c.z) * geometry_gradphis(c.e,0,c.x,c.y) + buffers_v_i(c.e,1,c.x,c.y,c.z) * geometry_gradphis(c.e,1,c.x,c.y)) * hvcoord_hybrid_bi_packed(c.z);

        buffers_phi_tens(c.e,c.x,c.y,c.z) = phi_tens;
      }
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO, bool PGRAD_CORRECTION>
void CaarFunctorImpl::epoch6_blockOps()
{
  auto buffers_dpnh_dp_i = viewAsReal(m_buffers.dpnh_dp_i);
  auto buffers_exner = viewAsReal(m_buffers.exner);
  auto buffers_grad_phinh_i = viewAsReal(m_buffers.grad_phinh_i);
  auto buffers_grad_w_i = viewAsReal(m_buffers.grad_w_i);
  auto buffers_v_tens = viewAsReal(m_buffers.v_tens);

  const int data_n0 = m_data.n0;

  auto &geometry_fcor = m_geometry.m_fcor;

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = viewAsReal(m_state.m_dp3d);
  auto state_v = viewAsReal(m_state.m_v);
  auto state_vtheta_dp= viewAsReal(m_state.m_vtheta_dp);
  auto state_w_i = viewAsReal(m_state.m_w_i);

  Kokkos::parallel_for(
    "caar compute_v_tens",
    SphereBlockOps::policy(m_num_elems, 6),
    KOKKOS_LAMBDA(const Team &team) {

      SphereBlockOps b(sg, team);
      if (b.skip()) return;

      SphereBlockScratch ttmp0(b);
      const Real w2 = (HYDROSTATIC) ? 0 : 0.25 * (state_w_i(b.e,data_n0,b.x,b.y,b.z) * state_w_i(b.e,data_n0,b.x,b.y,b.z) + state_w_i(b.e,data_n0,b.x,b.y,b.z+1) * state_w_i(b.e,data_n0,b.x,b.y,b.z+1));
      b.gradInit(ttmp0, w2);

      SphereBlockScratch ttmp1(b);
      const Real exneriz = buffers_exner(b.e,b.x,b.y,b.z);
      b.gradInit(ttmp1, exneriz);

      SphereBlockScratch ttmp2(b);
      const Real log_exneriz = (PGRAD_CORRECTION) ? log(exneriz) : 0;
      b.gradInit(ttmp2, log_exneriz);

      const Real v0 = state_v(b.e,data_n0,0,b.x,b.y,b.z);
      const Real v1 = state_v(b.e,data_n0,1,b.x,b.y,b.z);

      SphereBlockScratch ttmp3(b);
      SphereBlockScratch ttmp4(b);
      b.vortInit(ttmp3, ttmp4, v0, v1);

      SphereBlockScratch ttmp5(b);
      b.gradInit(ttmp5, 0.5 * (v0 * v0 + v1 * v1));

      b.barrier();

      Real grad_v0, grad_v1;
      b.grad(grad_v0, grad_v1, ttmp5);

      Real u_tens = (RSPLIT_ZERO) ? buffers_v_tens(b.e,0,b.x,b.y,b.z) : 0;
      Real v_tens = (RSPLIT_ZERO) ? buffers_v_tens(b.e,1,b.x,b.y,b.z) : 0;
      u_tens += grad_v0;
      v_tens += grad_v1;

      const Real cp_vtheta = PhysicalConstants::cp * (state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z) / state_dp3d(b.e,data_n0,b.x,b.y,b.z));

      Real grad_exner0, grad_exner1;
      b.grad(grad_exner0, grad_exner1, ttmp1);

      u_tens += cp_vtheta * grad_exner0;
      v_tens += cp_vtheta * grad_exner1;

      Real mgrad_x, mgrad_y;
      if (HYDROSTATIC) {

        mgrad_x = 0.5 * (buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z) + buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z+1));
        mgrad_y = 0.5 * (buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z) + buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z+1));

      } else {

        mgrad_x = 0.5 * (buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z) * buffers_dpnh_dp_i(b.e,b.x,b.y,b.z) + buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z+1) * buffers_dpnh_dp_i(b.e,b.x,b.y,b.z+1));
        mgrad_y = 0.5 * (buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z) * buffers_dpnh_dp_i(b.e,b.x,b.y,b.z) + buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z+1) * buffers_dpnh_dp_i(b.e,b.x,b.y,b.z+1));

      }

      if (PGRAD_CORRECTION) {

        Real grad_lexner0, grad_lexner1;
        b.grad(grad_lexner0, grad_lexner1, ttmp2);

        namespace PC = PhysicalConstants;
        constexpr Real cpt0 = PC::cp * (PC::Tref - PC::Tref_lapse_rate * PC::Tref * PC::cp / PC::g);
        mgrad_x += cpt0 * (grad_lexner0 - grad_exner0 / exneriz);
        mgrad_y += cpt0 * (grad_lexner1 - grad_exner1 / exneriz);
      }

      Real wvor_x = 0;
      Real wvor_y = 0;
      if (!HYDROSTATIC) {
        b.grad(wvor_x, wvor_y, ttmp0);
        wvor_x -= 0.5 * (buffers_grad_w_i(b.e,0,b.x,b.y,b.z) * state_w_i(b.e,data_n0,b.x,b.y,b.z) + buffers_grad_w_i(b.e,0,b.x,b.y,b.z+1) * state_w_i(b.e,data_n0,b.x,b.y,b.z+1));
        wvor_y -= 0.5 * (buffers_grad_w_i(b.e,1,b.x,b.y,b.z) * state_w_i(b.e,data_n0,b.x,b.y,b.z) + buffers_grad_w_i(b.e,1,b.x,b.y,b.z+1) * state_w_i(b.e,data_n0,b.x,b.y,b.z+1));
      }

      u_tens += mgrad_x + wvor_x;
      v_tens += mgrad_y + wvor_y;

      const Real vort = b.vort(ttmp3, ttmp4) + geometry_fcor(b.e,b.x,b.y);
      u_tens -= v1 * vort;
      v_tens += v0 * vort;

      buffers_v_tens(b.e,0,b.x,b.y,b.z) = u_tens;
      buffers_v_tens(b.e,1,b.x,b.y,b.z) = v_tens;
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch7_col()
{
  auto buffers_dp_tens = viewAsReal(m_buffers.dp_tens);
  auto buffers_eta_dot_dpdn = viewAsReal(m_buffers.eta_dot_dpdn);
  auto buffers_phi_tens = viewAsReal(m_buffers.phi_tens);
  auto buffers_theta_tens = viewAsReal(m_buffers.theta_tens);
  auto buffers_v_tens = viewAsReal(m_buffers.v_tens);
  auto buffers_vtheta_i = viewAsReal(m_buffers.vtheta_i);
  auto buffers_w_tens = viewAsReal(m_buffers.w_tens);

  const Real data_dt = m_data.dt;
  const int data_nm1 = m_data.nm1;
  const int data_np1 = m_data.np1;
  const Real data_scale3 = m_data.scale3;

  auto &geometry_spheremp = m_geometry.m_spheremp;

  const Real scale1_dt = m_data.scale1 * m_data.dt;

  auto state_dp3d = viewAsReal(m_state.m_dp3d);
  auto state_phinh_i = viewAsReal(m_state.m_phinh_i);
  auto state_v = viewAsReal(m_state.m_v);
  auto state_vtheta_dp = viewAsReal(m_state.m_vtheta_dp);
  auto state_w_i = viewAsReal(m_state.m_w_i);

  Kokkos::parallel_for(
    __PRETTY_FUNCTION__,
    SphereCol::policy(m_num_elems, NUM_PHYSICAL_LEV),
    KOKKOS_LAMBDA(const Team &team) {

      const SphereCol c(team);

      const Real spheremp = geometry_spheremp(c.e,c.x,c.y);
      const Real scale1_dt_spheremp = scale1_dt * spheremp;
      const Real scale3_spheremp = data_scale3 * spheremp;

      Real dp_tens = buffers_dp_tens(c.e,c.x,c.y,c.z);
      dp_tens *= scale1_dt_spheremp;
      Real dp_np1 = scale3_spheremp * state_dp3d(c.e,data_nm1,c.x,c.y,c.z);
      dp_np1 -= dp_tens;
      state_dp3d(c.e,data_np1,c.x,c.y,c.z) = dp_np1;

      Real theta_tens = buffers_theta_tens(c.e,c.x,c.y,c.z);
      if (RSPLIT_ZERO) {
        const Real etap = buffers_eta_dot_dpdn(c.e,c.x,c.y,c.z+1);
        const Real etaz = buffers_eta_dot_dpdn(c.e,c.x,c.y,c.z);
        const Real thetap = etap * buffers_vtheta_i(c.e,c.x,c.y,c.z+1);
        const Real thetaz = etaz * buffers_vtheta_i(c.e,c.x,c.y,c.z);
        theta_tens += thetap - thetaz;
      }
      theta_tens *= -scale1_dt_spheremp;
      Real vtheta_np1 = state_vtheta_dp(c.e,data_nm1,c.x,c.y,c.z);
      vtheta_np1 *= scale3_spheremp;
      vtheta_np1 += theta_tens;
      state_vtheta_dp(c.e,data_np1,c.x,c.y,c.z) = vtheta_np1;

      Real u_tens = buffers_v_tens(c.e,0,c.x,c.y,c.z);
      u_tens *= -scale1_dt_spheremp;
      Real u_np1 = state_v(c.e,data_nm1,0,c.x,c.y,c.z);
      u_np1 *= scale3_spheremp;
      u_np1 += u_tens;
      state_v(c.e,data_np1,0,c.x,c.y,c.z) = u_np1;

      Real v_tens = buffers_v_tens(c.e,1,c.x,c.y,c.z);
      v_tens *= -scale1_dt_spheremp;
      Real v_np1 = state_v(c.e,data_nm1,1,c.x,c.y,c.z);
      v_np1 *= scale3_spheremp;
      v_np1 += v_tens;
      state_v(c.e,data_np1,1,c.x,c.y,c.z) = v_np1;

      if (!HYDROSTATIC) {

        const Real dt_spheremp = data_dt * spheremp;

        Real phi_tens = buffers_phi_tens(c.e,c.x,c.y,c.z);
        phi_tens *= dt_spheremp;
        Real phi_np1 = state_phinh_i(c.e,data_nm1,c.x,c.y,c.z);
        phi_np1 *= scale3_spheremp;
        phi_np1 += phi_tens;
        state_phinh_i(c.e,data_np1,c.x,c.y,c.z) = phi_np1;

        Real w_tens = buffers_w_tens(c.e,c.x,c.y,c.z);
        w_tens *= dt_spheremp;
        Real w_np1 = state_w_i(c.e,data_nm1,c.x,c.y,c.z);
        w_np1 *= scale3_spheremp;
        w_np1 += w_tens;
        state_w_i(c.e,data_np1,c.x,c.y,c.z) = w_np1;

        if (c.z == NUM_PHYSICAL_LEV-1) {
          Real w_tens = buffers_w_tens(c.e,c.x,c.y,NUM_PHYSICAL_LEV);
          w_tens *= dt_spheremp;
          buffers_w_tens(c.e,c.x,c.y,NUM_PHYSICAL_LEV) = w_tens;
          Real w_np1 = state_w_i(c.e,data_nm1,c.x,c.y,NUM_PHYSICAL_LEV);
          w_np1 *= scale3_spheremp;
          w_np1 += w_tens;
          state_w_i(c.e,data_np1,c.x,c.y,NUM_PHYSICAL_LEV) = w_np1;
        }
      }
    });
}

void CaarFunctorImpl::caar_compute() 
{
  if (m_theta_hydrostatic_mode) {
    if (m_theta_advection_form == AdvectionForm::Conservative) epoch1_blockOps<true,true>();
    else epoch1_blockOps<true,false>();
  } else {
    if (m_theta_advection_form == AdvectionForm::Conservative) epoch1_blockOps<false,true>();
    else epoch1_blockOps<false,false>();
  }

  if (m_rsplit == 0) epoch2_scanOps<true>();
  else epoch2_scanOps<false>();

  if (m_theta_hydrostatic_mode) {
    if (m_rsplit == 0) epoch3_blockOps<true,true>();
    else epoch3_blockOps<true,false>();
  } else {
    if (m_rsplit == 0) epoch3_blockOps<false,true>();
    else epoch3_blockOps<false,false>();
  }

  if (m_theta_hydrostatic_mode) epoch4_scanOps();

  if (m_theta_hydrostatic_mode) {
    if (m_rsplit == 0) epoch5_colOps<true,true>();
    else epoch5_colOps<true,false>();
  } else {
    if (m_rsplit == 0) epoch5_colOps<false,true>();
    else epoch5_colOps<false,false>();
  }

  if (m_theta_hydrostatic_mode) {
    if (m_rsplit == 0) {
      if (m_pgrad_correction) epoch6_blockOps<true,true,true>();
      else epoch6_blockOps<true,true,false>();
    } else {
      if (m_pgrad_correction) epoch6_blockOps<true,false,true>();
      else epoch6_blockOps<true,false,false>();
    }
  } else {
    if (m_rsplit == 0) {
      if (m_pgrad_correction) epoch6_blockOps<false,true,true>();
      else epoch6_blockOps<false,true,false>();
    } else {
      if (m_pgrad_correction) epoch6_blockOps<false,false,true>();
      else epoch6_blockOps<false,false,false>();
    }
  }

  if (m_theta_hydrostatic_mode) {
    if (m_rsplit == 0) epoch7_col<true,true>();
    else epoch7_col<true,false>();
  } else {
    if (m_rsplit == 0) epoch7_col<false,true>();
    else epoch7_col<false,false>();
  }
}

}
