/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImpl.hpp"

namespace Homme {

namespace {

constexpr bool print_tracer_sgs_diffusivity_clipping = true;
constexpr Real tracer_sgs_cfl_target = 1.00;

KOKKOS_INLINE_FUNCTION
constexpr Real get_lambda_vis_ct ()
{
  switch (NP) {
  case 2: return 12.0;
  case 3: return 30.0;
  case 4: return 91.6742;
  case 5: return 190.1176;
  case 6: return 374.7788;
  case 7: return 652.3015;
  default: return 0.0;
  }
}

KOKKOS_INLINE_FUNCTION
Real get_local_laplace_metric_ct (const Real a, const Real b, const Real c, const Real d,
                                  const Real lambda_vis, const Real scale_factor_inv)
{
  const Real s11 = a*a + c*c;
  const Real s22 = b*b + d*d;
  const Real s12 = a*b + c*d;
  const Real disc = (s11 - s22)*(s11 - s22) + 4.0*s12*s12;
  const Real max_eig = 0.5 * (s11 + s22 + std::sqrt(disc));
  const Real norm_dinv = std::sqrt(max_eig);
  return lambda_vis * (scale_factor_inv * norm_dinv) * (scale_factor_inv * norm_dinv);
}

KOKKOS_INLINE_FUNCTION
bool is_tom_sponge_level_ct (const Real eta_top, const Real eta_mid)
{
  const Real ptop_over_press = eta_top / eta_mid;
  const Real ptop_over_press_sq = ptop_over_press * ptop_over_press;
  const Real nu_scale_top = 16.0 * ptop_over_press_sq / (ptop_over_press_sq + 1.0);
  return nu_scale_top >= 0.15;
}

} // namespace

void ComposeTransportImpl::advance_horizontal_turbulent_diffusion_scalar (const Real dt_q) {
  const auto dt = dt_q / m_data.hv_subcycle_q_sgs;
  const auto hv_q = m_data.hv_q;
  const auto Qtens = m_tracers.qtens_biharmonic;
  const auto Q = m_tracers.Q;
  const auto Kh = m_derived.m_turb_diff_heat;
  const auto dinv = m_geometry.m_dinv;
  const auto spheremp = m_geometry.m_spheremp;
  const Real scale_factor_inv = 1.0 / m_geometry.m_scale_factor;
  const Real lambda_vis = get_lambda_vis_ct();
  const auto tu_ne_hv_q = m_tu_ne_hv_q;
  const auto sphere_ops = m_sphere_ops;
  const bool tom_sponge_active = Context::singleton().get<SimulationParams>().nu_top > 0;
  int num_tom_sponge_levels = 0;
  if (tom_sponge_active) {
    const auto etam_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_hvcoord.etam);
    const auto etai_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), m_hvcoord.etai);
    const Real eta_top = etai_h(0) == 0.0 ? etam_h(0)[0] : etai_h(0);
    for (int phys_lev = 0; phys_lev < NUM_PHYSICAL_LEV; ++phys_lev) {
      const int lev = phys_lev / VECTOR_SIZE;
      const int s = phys_lev % VECTOR_SIZE;
      if (is_tom_sponge_level_ct(eta_top, etam_h(lev)[s])) {
        num_tom_sponge_levels = phys_lev + 1;
      }
    }
  }

  if (print_tracer_sgs_diffusivity_clipping && lambda_vis > 0) {
    const auto kh_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kh);
    const auto dinv_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dinv);
    const int rank = Context::singleton().get<Comm>().rank();

    for (int ie = 0; ie < m_data.nelemd; ++ie) {
      for (int i = 0; i < NP; ++i) {
        for (int j = 0; j < NP; ++j) {
          const Real a = dinv_h(ie,0,0,i,j);
          const Real b = dinv_h(ie,0,1,i,j);
          const Real c = dinv_h(ie,1,0,i,j);
          const Real d = dinv_h(ie,1,1,i,j);
          const Real laplace_metric = get_local_laplace_metric_ct(a, b, c, d,
                                                                  lambda_vis, scale_factor_inv);
          if (laplace_metric <= 0) continue;

          const Real max_diffusivity = 2.0 * tracer_sgs_cfl_target / (dt * laplace_metric);
          for (int lev = 0; lev < NUM_LEV; ++lev) {
            const auto kh = kh_h(ie,i,j,lev);
            for (int s = 0; s < VECTOR_SIZE; ++s) {
              const int phys_lev = lev * VECTOR_SIZE + s;
              if (phys_lev >= NUM_PHYSICAL_LEV) continue;
              if (phys_lev < num_tom_sponge_levels) continue;
              if (kh[s] > max_diffusivity) {
                printf("Warning: rank %d clipped tracer SGS Kh at ie=%d igp=%d jgp=%d lev=%d from %.16e to %.16e.\n",
                       rank, ie, i, j, phys_lev, kh[s], max_diffusivity);
              }
            }
          }
        }
      }
    }
  }

  for (int it = 0; it < m_data.hv_subcycle_q_sgs; ++it) {
    { // Qtens = Q
      const auto f = KOKKOS_LAMBDA (const int idx) {
        int ie, q, i, j, lev;
        idx_ie_q_ij_nlev<num_lev_pack>(hv_q, idx, ie, q, i, j, lev);
        Qtens(ie,q,i,j,lev) = Q(ie,q,i,j,lev);
      };
      launch_ie_q_ij_nlev<num_lev_pack>(hv_q, f);
    }
    // biharmonic_wk_scalar
    const auto laplace_simple_Qtens = [&] () {
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        KernelVariables kv(team, hv_q, tu_ne_hv_q);
        const auto Qtens_ie = Homme::subview(Qtens, kv.ie, kv.iq);
        sphere_ops.laplace_simple(kv, Qtens_ie, Qtens_ie);
      };
      Kokkos::parallel_for(m_tp_ne_hv_q, f);
    };
    laplace_simple_Qtens();
    m_hv_dss_be[0]->exchange(m_geometry.m_rspheremp);

    Kokkos::fence();

    { // Compute Q = Q spheremp - dt Kh Qtens. N.B. spheremp is already in
      // Qtens from divergence_sphere_wk.
      const auto f = KOKKOS_LAMBDA (const int idx) {
        int ie, q, i, j, lev;
        idx_ie_q_ij_nlev<num_lev_pack>(hv_q, idx, ie, q, i, j, lev);
        auto kh_eff = Kh(ie,i,j,lev);

        if (lambda_vis > 0) {
          const Real a = dinv(ie,0,0,i,j);
          const Real b = dinv(ie,0,1,i,j);
          const Real c = dinv(ie,1,0,i,j);
          const Real d = dinv(ie,1,1,i,j);
          const Real laplace_metric = get_local_laplace_metric_ct(a, b, c, d,
                                                                  lambda_vis, scale_factor_inv);
          if (laplace_metric > 0) {
            const Real max_diffusivity = 2.0 * tracer_sgs_cfl_target / (dt * laplace_metric);
            for (int s = 0; s < VECTOR_SIZE; ++s) {
              const int phys_lev = lev * VECTOR_SIZE + s;
              if (phys_lev >= num_tom_sponge_levels &&
                  phys_lev < NUM_PHYSICAL_LEV && kh_eff[s] > max_diffusivity) {
                kh_eff[s] = max_diffusivity;
              }
            }
          }
        }
        auto q_new = Q(ie,q,i,j,lev) * spheremp(ie,i,j);
        for (int s = 0; s < VECTOR_SIZE; ++s) {
          const int phys_lev = lev * VECTOR_SIZE + s;
          if (phys_lev < NUM_PHYSICAL_LEV &&
              phys_lev >= num_tom_sponge_levels) {
            q_new[s] = (Q(ie,q,i,j,lev)[s] * spheremp(ie,i,j)
                        - dt * kh_eff[s] * Qtens(ie,q,i,j,lev)[s]);
          }
        }
        Q(ie,q,i,j,lev) = q_new;
      };
      launch_ie_q_ij_nlev<num_lev_pack>(hv_q, f);
    }
    // Halo exchange Q and apply rspheremp.
    Kokkos::fence();
    m_hv_dss_be[1]->exchange(m_geometry.m_rspheremp);
  }
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
