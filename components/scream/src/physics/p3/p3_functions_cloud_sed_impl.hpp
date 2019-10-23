#ifndef P3_FUNCTIONS_CLOUD_SED_IMPL_HPP
#define P3_FUNCTIONS_CLOUD_SED_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 cloud sedimentation function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_sedimentation(
    const uview_1d<const Spack>& qc_incld,
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& lcldm,
    const uview_1d<const Spack>& acn,
    const uview_1d<const Spack>& inv_dzq,
    const view_dnu_table& dnu,
    const MemberType& team,
    const Workspace& workspace,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& odt, const bool& log_predictNc,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& mu_c,
    const uview_1d<Spack>& lamc,
    const uview_1d<Spack>& qc_tend,
    const uview_1d<Spack>& nc_tend,
    Scalar& prt_liq)
{
  // Get temporary workspaces needed for the cloud-sed calculation
  uview_1d<Spack> V_qc, V_nc, flux_qx, flux_nx;
  workspace.template take_many<4>(
    {"V_qc", "V_nc", "flux_qx", "flux_nx"},
    {&V_qc, &V_nc, &flux_qx, &flux_nx});

  view_1d_ptr_array<Spack, 2>
    fluxes_ptr = {&flux_qx, &flux_nx},
    vs_ptr     = {&V_qc, &V_nc},
    qnr_ptr    = {&qc, &nc};

  // find top, determine qxpresent
  const auto sqc = scalarize(qc);
  bool log_qxpresent;
  const Int k_qxtop = find_top(team, sqc, C::QSMALL, kbot, ktop, kdir, log_qxpresent);

  if (log_qxpresent) {
    Scalar dt_left = dt;    // time remaining for sedi over full model (mp) time step
    Scalar prt_accum = 0.0; // precip rate for individual category

    // find bottom
    Int k_qxbot = find_bottom(team, sqc, C::QSMALL, kbot, k_qxtop, kdir, log_qxpresent);

    while (dt_left > 1.e-4) {
      Scalar Co_max = 0.0;
      Int kmin, kmax;

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, V_qc.extent(0)), [&] (Int k) {
          V_qc(k) = 0;
          if (log_predictNc) {
            V_nc(k) = 0;
          }
      });
      team.team_barrier();

      // Convert top/bot to pack indices
      util::set_min_max(k_qxbot, k_qxtop, kmin, kmax, Spack::n);

      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int pk_, Scalar& lmax) {
          const int pk = kmin + pk_;
          auto qc_gt_small = (qc_incld(pk) > C::QSMALL);
          if (qc_gt_small.any()) {
            // compute Vq, Vn
            Spack nu, cdist, cdist1, dum;
            get_cloud_dsd2(qc_gt_small, qc_incld(pk), nc_incld(pk), mu_c(pk), rho(pk), nu, dnu, lamc(pk), cdist, cdist1, lcldm(pk));
            nc(pk).set(qc_gt_small, nc_incld(pk)*lcldm(pk));
            dum = 1 / (lamc(pk)*lamc(pk)); // bcn = 2
            V_qc(pk).set(qc_gt_small, acn(pk)*pack::tgamma(4 + 2 + mu_c(pk)) * dum / (pack::tgamma(mu_c(pk)+4)));
            if (log_predictNc) {
              V_nc(pk).set(qc_gt_small, acn(pk)*pack::tgamma(1 + 2 + mu_c(pk)) * dum / (pack::tgamma(mu_c(pk)+1)));
            }

            const auto Co_max_local = max(qc_gt_small, -1,
                                          V_qc(pk) * dt_left * inv_dzq(pk));
            if (Co_max_local > lmax)
              lmax = Co_max_local;
          }
      }, Kokkos::Max<Scalar>(Co_max));
      team.team_barrier();

      generalized_sedimentation<2>(rho, inv_rho, inv_dzq, team, nk, k_qxtop, k_qxbot, kbot, kdir, Co_max, dt_left, prt_accum, fluxes_ptr, vs_ptr, qnr_ptr);
    }

    Kokkos::single(
      Kokkos::PerTeam(team), [&] () {
        prt_liq = prt_accum * C::INV_RHOW * odt;
      });
  }

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, qc_tend.extent(0)), [&] (int pk) {
      qc_tend(pk) = (qc(pk) - qc_tend(pk)) * odt; // Liq. sedimentation tendency, measure
      nc_tend(pk) = (nc(pk) - nc_tend(pk)) * odt; // Liq. # sedimentation tendency, measure
  });

  workspace.release(V_qc);
  workspace.release(V_nc);
  workspace.release(flux_qx);
  workspace.release(flux_nx);
}

} // namespace p3
} // namespace scream

#endif
