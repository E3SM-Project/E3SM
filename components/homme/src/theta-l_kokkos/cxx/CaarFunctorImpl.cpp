#include "Config.hpp"

#ifdef HOMMEXX_ENABLE_CAAR_OPT

#include "CaarFunctorImpl-caar-opt.hpp"

namespace Homme {

template <bool HYDROSTATIC, AdvectionForm ADV_FORM>
void CaarFunctorImpl::epoch1_blockOps(const SphereOuter::Team &team)
{
  auto buffers_dp_tens = m_buffers.dp_tens;
  auto buffers_exner = m_buffers.exner;
  auto buffers_phi = m_buffers.phi;
  auto buffers_pnh = m_buffers.pnh;
  auto buffers_theta_tens = m_buffers.theta_tens;
  auto buffers_vdp = m_buffers.vdp;

  const Real data_eta_ave_w = m_data.eta_ave_w;
  const int data_n0 = m_data.n0;

  auto derived_vn0 = m_derived.m_vn0;

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = m_state.m_dp3d;
  auto state_phinh_i = m_state.m_phinh_i;
  auto state_v = m_state.m_v;
  auto state_vtheta_dp= m_state.m_vtheta_dp;

  SphereBlockOps::parallel_for(
    team, 5,
    KOKKOS_LAMBDA(const SphereBlockOps::Team &team) {
      SphereBlockOps b(sg, team);
      SphereBlockScratch ttmp0(b), ttmp1(b), ttmp2(b), ttmp3(b), ttmp4(b);

      SPHERE_BLOCK_START3(b, v0, v1, vtheta);

      if (!HYDROSTATIC) {

        const Scalar *const phinh_above = (b.z < NUM_LEV_P-1) ? &state_phinh_i(b.e,data_n0,b.x,b.y,b.z+1) : nullptr;
        Scalar dphi = b.zabove(&state_phinh_i(b.e,data_n0,b.x,b.y,b.z), phinh_above) - state_phinh_i(b.e,data_n0,b.x,b.y,b.z);
        const Scalar vtheta_dp = state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        for (int i = 0; i < VECTOR_SIZE; i++) {
          if ((vtheta_dp[i] < 0) || (dphi[i] > 0)) {
            if (b.z * VECTOR_SIZE + i < NUM_PHYSICAL_LEV) abort();
          }
        }

        // Sanitize dphi in padding slots before compute_pnh_and_exner so we never divide 0/0.
        if (NUM_PHYSICAL_LEV % VECTOR_SIZE != 0 && b.z == NUM_LEV-1) {
          for (int i = NUM_PHYSICAL_LEV % VECTOR_SIZE; i < VECTOR_SIZE; i++) dphi[i] = -1.0;
        }

        EquationOfState::compute_pnh_and_exner(vtheta_dp, dphi, buffers_pnh(b.e,b.x,b.y,b.z), buffers_exner(b.e,b.x,b.y,b.z));

        // Zero out padding slots in the last pack to prevent NaN propagation.
        // When NUM_PHYSICAL_LEV % VECTOR_SIZE != 0, the last midpoint pack (z=NUM_LEV-1)
        // has slots NUM_PHYSICAL_LEV%VECTOR_SIZE..VECTOR_SIZE-1 that are padding (vtheta_dp=0,
        // dphi=0 or -phis), causing 0/0=NaN in compute_pnh_and_exner.
        if (NUM_PHYSICAL_LEV % VECTOR_SIZE != 0 && b.z == NUM_LEV-1) {
          for (int i = NUM_PHYSICAL_LEV % VECTOR_SIZE; i < VECTOR_SIZE; i++) {
            buffers_pnh(b.e,b.x,b.y,b.z)[i] = 1.0;
            buffers_exner(b.e,b.x,b.y,b.z)[i] = 1.0;
          }
        }

        buffers_phi(b.e,b.x,b.y,b.z) = 0.5 * (state_phinh_i(b.e,data_n0,b.x,b.y,b.z) + b.zabove(&state_phinh_i(b.e,data_n0,b.x,b.y,b.z), phinh_above));
      }

      Scalar dp = state_dp3d(b.e,data_n0,b.x,b.y,b.z);
      v0 = state_v(b.e,data_n0,0,b.x,b.y,b.z) * dp;
      v1 = state_v(b.e,data_n0,1,b.x,b.y,b.z) * dp;
      if (NUM_PHYSICAL_LEV % VECTOR_SIZE != 0 && b.z == NUM_LEV-1) {
        for (int i = NUM_PHYSICAL_LEV % VECTOR_SIZE; i < VECTOR_SIZE; i++) dp[i] = 1.0;
      }
      b.divInit(ttmp0, ttmp1, v0, v1);
      buffers_vdp(b.e,0,b.x,b.y,b.z) = v0;
      buffers_vdp(b.e,1,b.x,b.y,b.z) = v1;
      derived_vn0(b.e,0,b.x,b.y,b.z) += data_eta_ave_w * v0;
      derived_vn0(b.e,1,b.x,b.y,b.z) += data_eta_ave_w * v1;

      vtheta = 0;

      if (ADV_FORM == AdvectionForm::Conservative) {
        const Scalar sv0 = state_v(b.e,data_n0,0,b.x,b.y,b.z) * state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        const Scalar sv1 = state_v(b.e,data_n0,1,b.x,b.y,b.z) * state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        b.divInit(ttmp2, ttmp3, sv0, sv1);
      } else if (ADV_FORM == AdvectionForm::NonConservative) {
        vtheta = state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z) / dp;
        b.gradInit(ttmp2, vtheta);
      } else { // AdvectionForm::Split
        const Scalar sv0 = state_v(b.e,data_n0,0,b.x,b.y,b.z) * state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        const Scalar sv1 = state_v(b.e,data_n0,1,b.x,b.y,b.z) * state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z);
        b.divInit(ttmp2, ttmp3, sv0, sv1);
        vtheta = state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z) / dp;
        b.gradInit(ttmp4, vtheta);
      }

      SPHERE_BLOCK_MIDDLE3(b, v0, v1, vtheta);

      const Scalar dvdp = b.div(ttmp0, ttmp1);
      buffers_dp_tens(b.e,b.x,b.y,b.z) = dvdp;

      if (ADV_FORM == AdvectionForm::Conservative) {
        buffers_theta_tens(b.e,b.x,b.y,b.z) = b.div(ttmp2, ttmp3);
      } else if (ADV_FORM == AdvectionForm::NonConservative) {
        Scalar grad0, grad1;
        b.grad(grad0, grad1, ttmp2);
        Scalar theta_tens = dvdp * vtheta;
        theta_tens += grad0 * v0;
        theta_tens += grad1 * v1;
        buffers_theta_tens(b.e,b.x,b.y,b.z) = theta_tens;
      } else { // AdvectionForm::Split
        const Scalar cons_tens = b.div(ttmp2, ttmp3);
        Scalar grad0, grad1;
        b.grad(grad0, grad1, ttmp4);
        const Scalar noncons_tens = dvdp * vtheta + grad0 * v0 + grad1 * v1;
        buffers_theta_tens(b.e,b.x,b.y,b.z) = 0.5 * (cons_tens + noncons_tens);
      }

      SPHERE_BLOCK_END();
    });
}

template <bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch2_scanOps(const SphereOuter::Team &team)
{
  auto buffers_dp_tens = viewAsReal(m_buffers.dp_tens);
  auto buffers_dp_i = viewAsReal(m_buffers.dp_i);
  auto buffers_eta_dot_dpdn = viewAsReal(m_buffers.eta_dot_dpdn);
  auto buffers_w_tens = viewAsReal(m_buffers.w_tens);

  const int data_n0 = m_data.n0;
  auto &hvcoord_hybrid_bi = m_hvcoord.hybrid_bi;
  const Real pi_i00 = m_hvcoord.ps0 * m_hvcoord.hybrid_ai0;
  auto state_dp3d = viewAsReal(m_state.m_dp3d);

  SphereScanOps::parallel_for(
    team,
    KOKKOS_LAMBDA(const SphereScanOps::Team &team) {
      SphereScanOps s(team);

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
          Kokkos::PerThread(s.t),
          [&]() {
            buffers_eta_dot_dpdn(s.e,s.x,s.y,0) = buffers_eta_dot_dpdn(s.e,s.x,s.y,NUM_PHYSICAL_LEV) = 0;
          });
      }
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch3_blockOps(const SphereOuter::Team &team)
{
  auto buffers_dp_i = m_buffers.dp_i;
  auto buffers_dp_tens = m_buffers.dp_tens;
  auto buffers_eta_dot_dpdn = m_buffers.eta_dot_dpdn;
  auto buffers_exner = m_buffers.exner;
  auto buffers_omega_p = m_buffers.omega_p;
  auto buffers_pnh = m_buffers.pnh;
  auto buffers_temp = m_buffers.temp;
  auto buffers_v_tens = m_buffers.v_tens;
  auto buffers_w_tens = m_buffers.w_tens;

  const Real data_eta_ave_w = m_data.eta_ave_w;
  const int data_n0 = m_data.n0;

  auto derived_eta_dot_dpdn = m_derived.m_eta_dot_dpdn;
  auto derived_omega_p = m_derived.m_omega_p;

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = m_state.m_dp3d;
  auto state_v = m_state.m_v;
  auto state_vtheta_dp= m_state.m_vtheta_dp;
  auto state_w_i = m_state.m_w_i;

  SphereBlockOps::parallel_for(
    team, 1,
    KOKKOS_LAMBDA(const SphereBlockOps::Team &team) {
      SphereBlockOps b(sg, team);
      SphereBlockScratch ttmp0(b);

      SPHERE_BLOCK_START0(b);

      const Scalar *const dp_i_above = (b.z < NUM_LEV_P-1) ? &buffers_dp_i(b.e,b.x,b.y,b.z+1) : nullptr;
      const Scalar pi = 0.5 * (buffers_dp_i(b.e,b.x,b.y,b.z) + b.zabove(&buffers_dp_i(b.e,b.x,b.y,b.z), dp_i_above));
      b.gradInit(ttmp0, pi);

      if (HYDROSTATIC) {
        Scalar p = pi;
        if (NUM_PHYSICAL_LEV % VECTOR_SIZE != 0 && b.z == NUM_LEV-1) {
          for (int i = NUM_PHYSICAL_LEV % VECTOR_SIZE; i < VECTOR_SIZE; i++) p[i] = PhysicalConstants::p0;
        }
        Scalar exner = p;
        EquationOfState::pressure_to_exner(exner);
        buffers_exner(b.e,b.x,b.y,b.z) = exner;
        buffers_pnh(b.e,b.x,b.y,b.z) = EquationOfState::compute_dphi(state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z), exner, p);
      }

      SPHERE_BLOCK_MIDDLE0(b);

      Scalar grad0, grad1;
      b.grad(grad0, grad1, ttmp0);

      const Scalar *const w_tens_above = (b.z < NUM_LEV_P-1) ? &buffers_w_tens(b.e,b.x,b.y,b.z+1) : nullptr;
      Scalar omega = -0.5 * (buffers_w_tens(b.e,b.x,b.y,b.z) + b.zabove(&buffers_w_tens(b.e,b.x,b.y,b.z), w_tens_above));
      const Scalar uz = state_v(b.e,data_n0,0,b.x,b.y,b.z);
      const Scalar vz = state_v(b.e,data_n0,1,b.x,b.y,b.z);
      omega += uz * grad0 + vz * grad1;
      buffers_omega_p(b.e,b.x,b.y,b.z) = omega;

      derived_omega_p(b.e,b.x,b.y,b.z) += data_eta_ave_w * omega;

      if (RSPLIT_ZERO) {

        Scalar dp = state_dp3d(b.e,data_n0,b.x,b.y,b.z);
        if (NUM_PHYSICAL_LEV % VECTOR_SIZE != 0 && b.z == NUM_LEV-1) {
          for (int i = NUM_PHYSICAL_LEV % VECTOR_SIZE; i < VECTOR_SIZE; i++) dp[i] = 1.0;
        }
        const Scalar *const eta_above = (b.z < NUM_LEV_P-1) ? &buffers_eta_dot_dpdn(b.e,b.x,b.y,b.z+1) : nullptr;
        const Scalar etap = b.zabove(&buffers_eta_dot_dpdn(b.e,b.x,b.y,b.z), eta_above);
        const Scalar etaz = buffers_eta_dot_dpdn(b.e,b.x,b.y,b.z);

        buffers_dp_tens(b.e,b.x,b.y,b.z) += etap - etaz;

        derived_eta_dot_dpdn(b.e,b.x,b.y,b.z) += data_eta_ave_w * etaz;

        if (!HYDROSTATIC) {
          const Scalar *const w_i_above = (b.z < NUM_LEV_P-1) ? &state_w_i(b.e,data_n0,b.x,b.y,b.z+1) : nullptr;
          const Scalar dw = b.zabove(&state_w_i(b.e,data_n0,b.x,b.y,b.z), w_i_above) - state_w_i(b.e,data_n0,b.x,b.y,b.z);
          const Scalar eta = 0.5 * (etaz + etap);
          buffers_temp(b.e,b.x,b.y,b.z) = dw * eta;
        }

        const Scalar *const v0_above = (b.z < NUM_LEV-1) ? &state_v(b.e,data_n0,0,b.x,b.y,b.z+1) : nullptr;
        const Scalar *const v1_above = (b.z < NUM_LEV-1) ? &state_v(b.e,data_n0,1,b.x,b.y,b.z+1) : nullptr;
        const Scalar *const v0_below = (b.z > 0) ? &state_v(b.e,data_n0,0,b.x,b.y,b.z-1) : nullptr;
        const Scalar *const v1_below = (b.z > 0) ? &state_v(b.e,data_n0,1,b.x,b.y,b.z-1) : nullptr;

        const Scalar facp = 0.5 * etap / dp;
        Scalar u = facp * (b.zabove(uz, &state_v(b.e,data_n0,0,b.x,b.y,b.z), v0_above) - uz);
        Scalar v = facp * (b.zabove(vz, &state_v(b.e,data_n0,1,b.x,b.y,b.z), v1_above) - vz);

        const Scalar facm = 0.5 * etaz / dp;
        u += facm * (uz - b.zbelow(uz, v0_below, &state_v(b.e,data_n0,0,b.x,b.y,b.z)));
        v += facm * (vz - b.zbelow(vz, v1_below, &state_v(b.e,data_n0,1,b.x,b.y,b.z)));
        buffers_v_tens(b.e,0,b.x,b.y,b.z) = u;
        buffers_v_tens(b.e,1,b.x,b.y,b.z) = v;
      }
      SPHERE_BLOCK_END();
    });
}

void CaarFunctorImpl::epoch4_scanOps(const SphereOuter::Team &team)
{
  auto buffers_phi = viewAsReal(m_buffers.phi);
  auto buffers_pnh = viewAsReal(m_buffers.pnh);

  const int data_n0 = m_data.n0;

  auto &geometry_phis = m_geometry.m_phis;

  auto state_phinh_i = viewAsReal(m_state.m_phinh_i);

  SphereScanOps::parallel_for(
    team,
    KOKKOS_LAMBDA(const SphereScanOps::Team &team) {
      SphereScanOps s(team);
      s.bwd_scan(state_phinh_i, data_n0, buffers_pnh, geometry_phis);
      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(s.t, NUM_PHYSICAL_LEV),
        [&](const int z) {
          buffers_phi(s.e,s.x,s.y,z) = 0.5 * (state_phinh_i(s.e,data_n0,s.x,s.y,z) + state_phinh_i(s.e,data_n0,s.x,s.y,z+1));
        });
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch5_colOps(const SphereOuter::Team &team)
{
  auto buffers_dp_i = m_buffers.dp_i;
  auto buffers_dpnh_dp_i = m_buffers.dpnh_dp_i;
  auto buffers_eta_dot_dpdn = m_buffers.eta_dot_dpdn;
  auto buffers_exner = m_buffers.exner;
  auto buffers_grad_phinh_i = m_buffers.grad_phinh_i;
  auto buffers_grad_w_i = m_buffers.grad_w_i;
  auto buffers_phi = m_buffers.phi;
  auto buffers_phi_tens = m_buffers.phi_tens;
  auto buffers_pnh = m_buffers.pnh;
  auto buffers_temp = m_buffers.temp;
  auto buffers_v_i = m_buffers.v_i;
  auto buffers_vtheta_i = m_buffers.vtheta_i;
  auto buffers_w_tens = m_buffers.w_tens;

  const int data_n0 = m_data.n0;
  const Real dscale = m_data.scale1 - m_data.scale2;

  auto &geometry_gradphis = m_geometry.m_gradphis;

  const Real gscale1 = m_data.scale1 * PhysicalConstants::g;
  const Real gscale2 = m_data.scale2 * PhysicalConstants::g;

  auto hvcoord_hybrid_bi_packed = m_hvcoord.hybrid_bi_packed;

  const Real ndata_scale1 = -m_data.scale1;
  const Real pi_i00 = m_hvcoord.ps0 * m_hvcoord.hybrid_ai0;

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = m_state.m_dp3d;
  auto state_phinh_i = m_state.m_phinh_i;
  auto state_v = m_state.m_v;
  auto state_w_i = m_state.m_w_i;

  SphereColOps::parallel_for(
    team, NUM_LEV_P,
    KOKKOS_LAMBDA(const SphereColOps::Team &team) {
      SphereColOps c(sg, team);

      c.grad(buffers_grad_phinh_i, state_phinh_i, data_n0);
      if (!HYDROSTATIC) c.grad(buffers_grad_w_i, state_w_i, data_n0);

      const Scalar *const dp3d_below = (c.z == 0) ? nullptr : &state_dp3d(c.e,data_n0,c.x,c.y,c.z-1);
      const Scalar *const dp3d_here  = (c.z == NUM_LEV) ? nullptr : &state_dp3d(c.e,data_n0,c.x,c.y,c.z);
      Scalar dm = c.zbelow(0, dp3d_below, dp3d_here);
      Scalar dz = c.zabove(0, dp3d_here);
      if (c.z == NUM_LEV_P-1) {
        for (int i = (NUM_PHYSICAL_LEV % VECTOR_SIZE) + 1; i < VECTOR_SIZE; i++) {
          dm[i] = 1.0;
          dz[i] = 1.0;
        }
      }
      const Scalar dp_i = c.ztrim(dz, dm, 0.5 * (dz + dm));
      buffers_dp_i(c.e,c.x,c.y,c.z) = dp_i;

      if (!HYDROSTATIC) {

        const Scalar *const v0_below = (c.z == 0) ? nullptr : &state_v(c.e,data_n0,0,c.x,c.y,c.z-1);
        const Scalar *const v0_here  = (c.z == NUM_LEV) ? nullptr : &state_v(c.e,data_n0,0,c.x,c.y,c.z);
        const Scalar v0m = c.zbelow(0, v0_below, v0_here);
        const Scalar v0z = c.zabove(0, v0_here);
        const Scalar v_i0 = c.ztrim(v0z, v0m, (dz * v0z + dm * v0m) / (dm + dz));
        buffers_v_i(c.e,0,c.x,c.y,c.z) = v_i0;

        const Scalar *const v1_below = (c.z == 0) ? nullptr : &state_v(c.e,data_n0,1,c.x,c.y,c.z-1);
        const Scalar *const v1_here  = (c.z == NUM_LEV) ? nullptr : &state_v(c.e,data_n0,1,c.x,c.y,c.z);
        const Scalar v1m = c.zbelow(0, v1_below, v1_here);
        const Scalar v1z = c.zabove(0, v1_here);
        const Scalar v_i1 = c.ztrim(v1z, v1m, (dz * v1z + dm * v1m) / (dm + dz));
        buffers_v_i(c.e,1,c.x,c.y,c.z) = v_i1;

        const Scalar *const pnh_below = (c.z == 0) ? nullptr : &buffers_pnh(c.e,c.x,c.y,c.z-1);
        const Scalar *const pnh_here  = (c.z == NUM_LEV) ? nullptr : &buffers_pnh(c.e,c.x,c.y,c.z);
        const Scalar pm = c.zbelow(pi_i00, pnh_below, pnh_here);
        const Scalar pz = c.zabove(pm + 0.5 * dm, pnh_here);
        buffers_dpnh_dp_i(c.e,c.x,c.y,c.z) = 2.0 * (pz - pm) / (dm + dz);
      }

      if (RSPLIT_ZERO) {

        const Scalar *const phi_below = (c.z == 0) ? nullptr : &buffers_phi(c.e,c.x,c.y,c.z-1);
        const Scalar *const phi_here  = (c.z == NUM_LEV) ? nullptr : &buffers_phi(c.e,c.x,c.y,c.z);
        const Scalar phim = c.zbelow(0, phi_below, phi_here);
        const Scalar phiz = c.zabove(0, phi_here);

        if (!HYDROSTATIC) {
          const Scalar phi_vadv = c.ztrim(0, 0, (phiz - phim) * buffers_eta_dot_dpdn(c.e,c.x,c.y,c.z) / dp_i);
          buffers_phi_tens(c.e,c.x,c.y,c.z) = phi_vadv;
        }

        Scalar vtheta_i = phiz - phim;

        const Scalar *const exner_below = (c.z == 0) ? nullptr : &buffers_exner(c.e,c.x,c.y,c.z-1);
        const Scalar *const exner_here  = (c.z == NUM_LEV) ? nullptr : &buffers_exner(c.e,c.x,c.y,c.z);
        const Scalar exnerm = c.zbelow(0, exner_below, exner_here);
        const Scalar exnerz = c.zabove(0, exner_here);
        Scalar dexner = c.ztrim(1, 1, exnerz - exnerm);
        if (c.z == NUM_LEV_P-1) {
          for (int i = (NUM_PHYSICAL_LEV % VECTOR_SIZE) + 1; i < VECTOR_SIZE; i++) dexner[i] = 1.0;
        }
        vtheta_i /= dexner;

        vtheta_i /= -PhysicalConstants::cp; // separate for BFB equality with original Fortran

        if (!HYDROSTATIC) vtheta_i *= buffers_dpnh_dp_i(c.e,c.x,c.y,c.z);
        buffers_vtheta_i(c.e,c.x,c.y,c.z) = vtheta_i;
      }

      if (!HYDROSTATIC) {

        Scalar w_tens = 0;
        if (RSPLIT_ZERO) {
          const Scalar *const temp_below = (c.z == 0) ? nullptr : &buffers_temp(c.e,c.x,c.y,c.z-1);
          const Scalar *const temp_here  = (c.z == NUM_LEV) ? nullptr : &buffers_temp(c.e,c.x,c.y,c.z);
          const Scalar tempm = c.zbelow(0, temp_below, temp_here);
          const Scalar tempz = c.zabove(0, temp_here);
          const Scalar dw = c.ztrim(tempz, tempm, 0.5 * (tempz + tempm));
          w_tens = dw / dp_i;
        }
        w_tens += buffers_v_i(c.e,0,c.x,c.y,c.z) * buffers_grad_w_i(c.e,0,c.x,c.y,c.z) + buffers_v_i(c.e,1,c.x,c.y,c.z) * buffers_grad_w_i(c.e,1,c.x,c.y,c.z);
        w_tens *= ndata_scale1;
        const Scalar gscale = c.ztop(gscale1, gscale2);
        w_tens += (buffers_dpnh_dp_i(c.e,c.x,c.y,c.z)-Scalar(1)) * gscale;
        buffers_w_tens(c.e,c.x,c.y,c.z) = w_tens;

        Scalar phi_tens = (RSPLIT_ZERO) ? buffers_phi_tens(c.e,c.x,c.y,c.z) : 0;
        phi_tens += buffers_v_i(c.e,0,c.x,c.y,c.z) * buffers_grad_phinh_i(c.e,0,c.x,c.y,c.z) + buffers_v_i(c.e,1,c.x,c.y,c.z) * buffers_grad_phinh_i(c.e,1,c.x,c.y,c.z);
        phi_tens *= ndata_scale1;
        phi_tens += state_w_i(c.e,data_n0,c.x,c.y,c.z) * gscale;

        if (dscale) phi_tens += dscale * (buffers_v_i(c.e,0,c.x,c.y,c.z) * geometry_gradphis(c.e,0,c.x,c.y) + buffers_v_i(c.e,1,c.x,c.y,c.z) * geometry_gradphis(c.e,1,c.x,c.y)) * hvcoord_hybrid_bi_packed(c.z);

        buffers_phi_tens(c.e,c.x,c.y,c.z) = phi_tens;
      }
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO, bool PGRAD_CORRECTION, bool SPLIT_FORM>
void CaarFunctorImpl::epoch6_blockOps(const SphereOuter::Team &team)
{
  auto buffers_dpnh_dp_i = m_buffers.dpnh_dp_i;
  auto buffers_exner = m_buffers.exner;
  auto buffers_grad_phinh_i = m_buffers.grad_phinh_i;
  auto buffers_grad_w_i = m_buffers.grad_w_i;
  auto buffers_v_tens = m_buffers.v_tens;

  const int data_n0 = m_data.n0;

  auto &geometry_fcor = m_geometry.m_fcor;

  const SphereGlobal sg(m_sphere_ops);

  auto state_dp3d = m_state.m_dp3d;
  auto state_v = m_state.m_v;
  auto state_vtheta_dp= m_state.m_vtheta_dp;
  auto state_w_i = m_state.m_w_i;

  SphereBlockOps::parallel_for(
    team, 8,
    KOKKOS_LAMBDA(const SphereBlockOps::Team &team) {
      SphereBlockOps b(sg, team);
      SphereBlockScratch ttmp0(b), ttmp1(b), ttmp2(b), ttmp3(b), ttmp4(b), ttmp5(b), ttmp6(b), ttmp7(b);

      SPHERE_BLOCK_START3(b, exneriz, v0, v1);

      if (!HYDROSTATIC) {
        const Scalar *const w_i_above = (b.z < NUM_LEV_P-1) ? &state_w_i(b.e,data_n0,b.x,b.y,b.z+1) : nullptr;
        const Scalar wz = b.zabove(&state_w_i(b.e,data_n0,b.x,b.y,b.z), w_i_above);
        const Scalar w2 = 0.25 * (state_w_i(b.e,data_n0,b.x,b.y,b.z) * state_w_i(b.e,data_n0,b.x,b.y,b.z) + wz * wz);
        b.gradInit(ttmp0, w2);
      }

      exneriz = buffers_exner(b.e,b.x,b.y,b.z);
      b.gradInit(ttmp1, exneriz);

      const Scalar log_exneriz = (PGRAD_CORRECTION) ? log(exneriz) : 0;
      b.gradInit(ttmp2, log_exneriz);

      v0 = state_v(b.e,data_n0,0,b.x,b.y,b.z);
      v1 = state_v(b.e,data_n0,1,b.x,b.y,b.z);

      b.vortInit(ttmp3, ttmp4, v0, v1);
      b.gradInit(ttmp5, 0.5 * (v0 * v0 + v1 * v1));

      if (SPLIT_FORM) {
        Scalar dp = state_dp3d(b.e,data_n0,b.x,b.y,b.z);
        if (NUM_PHYSICAL_LEV % VECTOR_SIZE != 0 && b.z == NUM_LEV-1) {
          for (int i = NUM_PHYSICAL_LEV % VECTOR_SIZE; i < VECTOR_SIZE; i++) dp[i] = 1.0;
        }
        const Scalar vtheta = state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z) / dp;
        b.gradInit(ttmp6, vtheta);
        b.gradInit(ttmp7, vtheta * exneriz);
      }

      SPHERE_BLOCK_MIDDLE3(b, exneriz, v0, v1);

      Scalar grad_v0, grad_v1;
      b.grad(grad_v0, grad_v1, ttmp5);

      Scalar u_tens = (RSPLIT_ZERO) ? buffers_v_tens(b.e,0,b.x,b.y,b.z) : 0;
      Scalar v_tens = (RSPLIT_ZERO) ? buffers_v_tens(b.e,1,b.x,b.y,b.z) : 0;
      u_tens += grad_v0;
      v_tens += grad_v1;

      Scalar dp = state_dp3d(b.e,data_n0,b.x,b.y,b.z);
      if (NUM_PHYSICAL_LEV % VECTOR_SIZE != 0 && b.z == NUM_LEV-1) {
        for (int i = NUM_PHYSICAL_LEV % VECTOR_SIZE; i < VECTOR_SIZE; i++) dp[i] = 1.0;
      }
      const Scalar cp_vtheta = PhysicalConstants::cp * (state_vtheta_dp(b.e,data_n0,b.x,b.y,b.z) / dp);

      Scalar grad_exner0, grad_exner1;
      b.grad(grad_exner0, grad_exner1, ttmp1);

      if (SPLIT_FORM) {
        Scalar grad_vtheta0, grad_vtheta1;
        b.grad(grad_vtheta0, grad_vtheta1, ttmp6);

        Scalar grad_vtheta_exner0, grad_vtheta_exner1;
        b.grad(grad_vtheta_exner0, grad_vtheta_exner1, ttmp7);

        u_tens += 0.5 * (cp_vtheta * grad_exner0 + PhysicalConstants::cp * (grad_vtheta_exner0 - exneriz * grad_vtheta0));
        v_tens += 0.5 * (cp_vtheta * grad_exner1 + PhysicalConstants::cp * (grad_vtheta_exner1 - exneriz * grad_vtheta1));
      } else {
        u_tens += cp_vtheta * grad_exner0;
        v_tens += cp_vtheta * grad_exner1;
      }

      Scalar mgrad_x, mgrad_y;
      const Scalar *const grad_phinh_i0_above = (b.z < NUM_LEV_P-1) ? &buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z+1) : nullptr;
      const Scalar *const grad_phinh_i1_above = (b.z < NUM_LEV_P-1) ? &buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z+1) : nullptr;
      if (HYDROSTATIC) {

        mgrad_x = 0.5 * (buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z) + b.zabove(&buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z), grad_phinh_i0_above));
        mgrad_y = 0.5 * (buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z) + b.zabove(&buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z), grad_phinh_i1_above));

      } else {

        const Scalar *const dpnh_dp_i_above = (b.z < NUM_LEV_P-1) ? &buffers_dpnh_dp_i(b.e,b.x,b.y,b.z+1) : nullptr;
        mgrad_x = 0.5 * (buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z) * buffers_dpnh_dp_i(b.e,b.x,b.y,b.z) + b.zabove(&buffers_grad_phinh_i(b.e,0,b.x,b.y,b.z), grad_phinh_i0_above) * b.zabove(&buffers_dpnh_dp_i(b.e,b.x,b.y,b.z), dpnh_dp_i_above));
        mgrad_y = 0.5 * (buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z) * buffers_dpnh_dp_i(b.e,b.x,b.y,b.z) + b.zabove(&buffers_grad_phinh_i(b.e,1,b.x,b.y,b.z), grad_phinh_i1_above) * b.zabove(&buffers_dpnh_dp_i(b.e,b.x,b.y,b.z), dpnh_dp_i_above));

      }

      if (PGRAD_CORRECTION) {

        Scalar grad_lexner0, grad_lexner1;
        b.grad(grad_lexner0, grad_lexner1, ttmp2);

        namespace PC = PhysicalConstants;
        constexpr Real cpt0 = PC::cp * (PC::Tref - PC::Tref_lapse_rate * PC::Tref * PC::cp / PC::g);
        mgrad_x += cpt0 * (grad_lexner0 - grad_exner0 / exneriz);
        mgrad_y += cpt0 * (grad_lexner1 - grad_exner1 / exneriz);
      }

      Scalar wvor_x = 0;
      Scalar wvor_y = 0;
      if (!HYDROSTATIC) {
        b.grad(wvor_x, wvor_y, ttmp0);
        const Scalar *const grad_w_i0_above = (b.z < NUM_LEV_P-1) ? &buffers_grad_w_i(b.e,0,b.x,b.y,b.z+1) : nullptr;
        const Scalar *const grad_w_i1_above = (b.z < NUM_LEV_P-1) ? &buffers_grad_w_i(b.e,1,b.x,b.y,b.z+1) : nullptr;
        const Scalar *const w_i_above = (b.z < NUM_LEV_P-1) ? &state_w_i(b.e,data_n0,b.x,b.y,b.z+1) : nullptr;
        wvor_x -= 0.5 * (buffers_grad_w_i(b.e,0,b.x,b.y,b.z) * state_w_i(b.e,data_n0,b.x,b.y,b.z) + b.zabove(&buffers_grad_w_i(b.e,0,b.x,b.y,b.z), grad_w_i0_above) * b.zabove(&state_w_i(b.e,data_n0,b.x,b.y,b.z), w_i_above));
        wvor_y -= 0.5 * (buffers_grad_w_i(b.e,1,b.x,b.y,b.z) * state_w_i(b.e,data_n0,b.x,b.y,b.z) + b.zabove(&buffers_grad_w_i(b.e,1,b.x,b.y,b.z), grad_w_i1_above) * b.zabove(&state_w_i(b.e,data_n0,b.x,b.y,b.z), w_i_above));
      }

      u_tens += mgrad_x + wvor_x;
      v_tens += mgrad_y + wvor_y;

      const Scalar vort = b.vort(ttmp3, ttmp4) + geometry_fcor(b.e,b.x,b.y);
      u_tens -= v1 * vort;
      v_tens += v0 * vort;

      buffers_v_tens(b.e,0,b.x,b.y,b.z) = u_tens;
      buffers_v_tens(b.e,1,b.x,b.y,b.z) = v_tens;

      SPHERE_BLOCK_END();
    });
}

template <bool HYDROSTATIC, bool RSPLIT_ZERO>
void CaarFunctorImpl::epoch7_col(const SphereOuter::Team &team)
{
  auto buffers_dp_tens = m_buffers.dp_tens;
  auto buffers_eta_dot_dpdn = m_buffers.eta_dot_dpdn;
  auto buffers_phi_tens = m_buffers.phi_tens;
  auto buffers_theta_tens = m_buffers.theta_tens;
  auto buffers_v_tens = m_buffers.v_tens;
  auto buffers_vtheta_i = m_buffers.vtheta_i;
  auto buffers_w_tens = m_buffers.w_tens;

  const Real data_dt = m_data.dt;
  const int data_nm1 = m_data.nm1;
  const int data_np1 = m_data.np1;
  const Real data_scale3 = m_data.scale3;

  auto &geometry_phis = m_geometry.m_phis;
  auto &geometry_spheremp = m_geometry.m_spheremp;

  const Real scale1_dt = m_data.scale1 * m_data.dt;

  auto state_dp3d = m_state.m_dp3d;
  auto state_phinh_i = m_state.m_phinh_i;
  auto state_v = m_state.m_v;
  auto state_vtheta_dp = m_state.m_vtheta_dp;
  auto state_w_i = m_state.m_w_i;

  SphereCol::parallel_for(
    team, NUM_LEV,
    KOKKOS_LAMBDA(const SphereCol::Team &team) {
      const SphereCol c(team);

      const Real spheremp = geometry_spheremp(c.e,c.x,c.y);
      const Real scale1_dt_spheremp = scale1_dt * spheremp;
      const Real scale3_spheremp = data_scale3 * spheremp;

      Scalar dp_tens = buffers_dp_tens(c.e,c.x,c.y,c.z);
      dp_tens *= scale1_dt_spheremp;
      Scalar dp_np1 = scale3_spheremp * state_dp3d(c.e,data_nm1,c.x,c.y,c.z);
      dp_np1 -= dp_tens;
      state_dp3d(c.e,data_np1,c.x,c.y,c.z) = dp_np1;

      Scalar theta_tens = buffers_theta_tens(c.e,c.x,c.y,c.z);
      if (RSPLIT_ZERO) {
        const Scalar *const eta_above = (c.z < NUM_LEV_P-1) ? &buffers_eta_dot_dpdn(c.e,c.x,c.y,c.z+1) : nullptr;
        const Scalar *const vtheta_above = (c.z < NUM_LEV_P-1) ? &buffers_vtheta_i(c.e,c.x,c.y,c.z+1) : nullptr;
        const Scalar etap = c.zplus(&buffers_eta_dot_dpdn(c.e,c.x,c.y,c.z), eta_above);
        const Scalar etaz = buffers_eta_dot_dpdn(c.e,c.x,c.y,c.z);
        const Scalar thetap = etap * c.zplus(&buffers_vtheta_i(c.e,c.x,c.y,c.z), vtheta_above);
        const Scalar thetaz = etaz * buffers_vtheta_i(c.e,c.x,c.y,c.z);
        theta_tens += thetap - thetaz;
      }
      theta_tens *= -scale1_dt_spheremp;
      Scalar vtheta_np1 = state_vtheta_dp(c.e,data_nm1,c.x,c.y,c.z);
      vtheta_np1 *= scale3_spheremp;
      vtheta_np1 += theta_tens;
      state_vtheta_dp(c.e,data_np1,c.x,c.y,c.z) = vtheta_np1;

      Scalar u_tens = buffers_v_tens(c.e,0,c.x,c.y,c.z);
      u_tens *= -scale1_dt_spheremp;
      Scalar u_np1 = state_v(c.e,data_nm1,0,c.x,c.y,c.z);
      u_np1 *= scale3_spheremp;
      u_np1 += u_tens;
      state_v(c.e,data_np1,0,c.x,c.y,c.z) = u_np1;

      Scalar v_tens = buffers_v_tens(c.e,1,c.x,c.y,c.z);
      v_tens *= -scale1_dt_spheremp;
      Scalar v_np1 = state_v(c.e,data_nm1,1,c.x,c.y,c.z);
      v_np1 *= scale3_spheremp;
      v_np1 += v_tens;
      state_v(c.e,data_np1,1,c.x,c.y,c.z) = v_np1;

      if (!HYDROSTATIC) {

        const Real dt_spheremp = data_dt * spheremp;

        Scalar phi_tens = buffers_phi_tens(c.e,c.x,c.y,c.z);
        phi_tens *= dt_spheremp;
        Scalar phi_np1 = state_phinh_i(c.e,data_nm1,c.x,c.y,c.z);
        phi_np1 *= scale3_spheremp;
        phi_np1 += phi_tens;
        state_phinh_i(c.e,data_np1,c.x,c.y,c.z) = phi_np1;

        Scalar w_tens = buffers_w_tens(c.e,c.x,c.y,c.z);
        w_tens *= dt_spheremp;
        Scalar w_np1 = state_w_i(c.e,data_nm1,c.x,c.y,c.z);
        w_np1 *= scale3_spheremp;
        w_np1 += w_tens;
        state_w_i(c.e,data_np1,c.x,c.y,c.z) = w_np1;

        if (c.z == NUM_LEV-1) {
          constexpr int z = NUM_LEV_P-1;
          constexpr int dz = NUM_PHYSICAL_LEV%VECTOR_SIZE;
          if (NUM_LEV == NUM_LEV_P) {
            state_phinh_i(c.e,data_np1,c.x,c.y,z)[dz] = geometry_phis(c.e,c.x,c.y);
          } else {
            Real w_tens = buffers_w_tens(c.e,c.x,c.y,z)[dz];
            w_tens *= dt_spheremp;
            buffers_w_tens(c.e,c.x,c.y,z)[dz] = w_tens;
            Real w_np1 = state_w_i(c.e,data_nm1,c.x,c.y,z)[dz];
            w_np1 *= scale3_spheremp;
            w_np1 += w_tens;
            state_w_i(c.e,data_np1,c.x,c.y,z)[dz] = w_np1;
          }
        }
      }
    });
}

void CaarFunctorImpl::caar_compute()
{
  SphereOuter::parallel_for(
    m_num_elems,
    KOKKOS_LAMBDA(const SphereOuter::Team &team) {

      if (m_theta_hydrostatic_mode) {
        if (m_theta_advection_form == AdvectionForm::Conservative) epoch1_blockOps<true, AdvectionForm::Conservative>(team);
        else if (m_theta_advection_form == AdvectionForm::NonConservative) epoch1_blockOps<true, AdvectionForm::NonConservative>(team);
        else epoch1_blockOps<true, AdvectionForm::Split>(team);
      } else {
        if (m_theta_advection_form == AdvectionForm::Conservative) epoch1_blockOps<false, AdvectionForm::Conservative>(team);
        else if (m_theta_advection_form == AdvectionForm::NonConservative) epoch1_blockOps<false, AdvectionForm::NonConservative>(team);
        else epoch1_blockOps<false, AdvectionForm::Split>(team);
      }

      if (m_rsplit == 0) epoch2_scanOps<true>(team);
      else epoch2_scanOps<false>(team);

      if (m_theta_hydrostatic_mode) {
        if (m_rsplit == 0) epoch3_blockOps<true,true>(team);
        else epoch3_blockOps<true,false>(team);
      } else {
        if (m_rsplit == 0) epoch3_blockOps<false,true>(team);
        else epoch3_blockOps<false,false>(team);
      }

      if (m_theta_hydrostatic_mode) epoch4_scanOps(team);

      if (m_theta_hydrostatic_mode) {
        if (m_rsplit == 0) epoch5_colOps<true,true>(team);
        else epoch5_colOps<true,false>(team);
      } else {
        if (m_rsplit == 0) epoch5_colOps<false,true>(team);
        else epoch5_colOps<false,false>(team);
      }

      const bool is_split = (m_theta_advection_form == AdvectionForm::Split);
      if (m_theta_hydrostatic_mode) {
        if (m_rsplit == 0) {
          if (m_pgrad_correction) {
            if (is_split) epoch6_blockOps<true,true,true,true>(team);
            else epoch6_blockOps<true,true,true,false>(team);
          } else {
            if (is_split) epoch6_blockOps<true,true,false,true>(team);
            else epoch6_blockOps<true,true,false,false>(team);
          }
        } else {
          if (m_pgrad_correction) {
            if (is_split) epoch6_blockOps<true,false,true,true>(team);
            else epoch6_blockOps<true,false,true,false>(team);
          } else {
            if (is_split) epoch6_blockOps<true,false,false,true>(team);
            else epoch6_blockOps<true,false,false,false>(team);
          }
        }
      } else {
        if (m_rsplit == 0) {
          if (m_pgrad_correction) {
            if (is_split) epoch6_blockOps<false,true,true,true>(team);
            else epoch6_blockOps<false,true,true,false>(team);
          } else {
            if (is_split) epoch6_blockOps<false,true,false,true>(team);
            else epoch6_blockOps<false,true,false,false>(team);
          }
        } else {
          if (m_pgrad_correction) {
            if (is_split) epoch6_blockOps<false,false,true,true>(team);
            else epoch6_blockOps<false,false,true,false>(team);
          } else {
            if (is_split) epoch6_blockOps<false,false,false,true>(team);
            else epoch6_blockOps<false,false,false,false>(team);
          }
        }
      }

      if (m_theta_hydrostatic_mode) {
        if (m_rsplit == 0) epoch7_col<true,true>(team);
        else epoch7_col<true,false>(team);
      } else {
        if (m_rsplit == 0) epoch7_col<false,true>(team);
        else epoch7_col<false,false>(team);
      }
    });
}

} // namespace Homme

#endif // HOMMEXX_ENABLE_CAAR_OPT

