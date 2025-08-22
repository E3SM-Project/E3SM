#ifndef GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP
#define GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "util/eamxx_utils.hpp"
#include <ekat_math_utils.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_compute_stress_profiles_and_diffusivities. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_compute_stress_profiles_and_diffusivities(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const GwCommonInit& init,
  const Int& pver,
  const Int& pgwv,
  const Int& src_level,
  const uview_1d<const Real>& ubi,
  const uview_1d<const Real>& c,
  const uview_1d<const Real>& rhoi,
  const uview_1d<const Real>& ni,
  const uview_1d<const Real>& kvtt,
  const uview_1d<const Real>& t,
  const uview_1d<const Real>& ti,
  const uview_1d<const Real>& piln,
  // Inputs/Outputs
  const uview_2d<Real>& tau)
{
  static const auto ubmc2mn = GWC::ubmc2mn;

  const int num_pgwv = 2*pgwv + 1;

  // Get temporary workspaces and change them to desired dimensions
  uview_1d<Real> tausat_1d, dsat_1d, wrk1_1d, wrk2_1d;
  workspace.template take_many_contiguous_unsafe<4>(
    {"tausat_1d", "dsat_1d", "wrk1_1d", "wrk2_1d"},
    {&tausat_1d, &dsat_1d, &wrk1_1d, &wrk2_1d});
  uview_2d<Real>
    tausat(tausat_1d.data(), pver+1, num_pgwv),
    dsat(dsat_1d.data(), pver+1, num_pgwv),
    wrk1(wrk1_1d.data(), pver+1, num_pgwv),
    wrk2(wrk2_1d.data(), pver+1, num_pgwv);

  // Loop from bottom to top to get stress profiles. Instead of having 2 levels
  // of parallelism, we collapse all the parallelism into the top level by multiplying
  // the level by num_pgwv.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, init.ktop*num_pgwv, (src_level+1)*num_pgwv), [&] (const int k_pgwv) {

    const int k = k_pgwv / num_pgwv;
    const int l = k_pgwv % num_pgwv;

    // Determine the absolute value of the saturation stress.
    // Define critical levels where the sign of (u-c) changes between interfaces.
    const Real ubmc = ubi(k) - c(l);

    // Test to see if u-c has the same sign here as the level below.
    if (ubmc * (ubi(k + 1) - c(l)) > 0.0) {
      tausat(k, l) = std::abs(init.effkwv * rhoi(k) * bfb_cube(ubmc) /
                                   (2.0 * ni(k)));
      if (tausat(k, l) <= GWC::taumin) tausat(k, l) = 0.0;
    }
    else {
      tausat(k, l) = 0.0;
    }

    if (!init.do_molec_diff) {
      dsat(k, l) = bfb_square(ubmc / ni(k)) *
        (init.effkwv * bfb_square(ubmc) /
         (GWC::rog * ti(k) * ni(k)) - init.alpha(k));
    }

    if (k <= init.nbot_molec || !init.do_molec_diff) {
      const Real ubmc2 = ekat::impl::max(bfb_square(ubmc), ubmc2mn);
      const Real at = ni(k) / (2.0 * init.kwv * ubmc2);
      const Real bt = init.alpha(k);
      const Real ct = bfb_square(ni(k)) / ubmc2;
      const Real et = -2.0 * GWC::rog * t(k) * (piln(k + 1) - piln(k));
      wrk1(k, l) = at*bt*et;
      wrk2(k, l) = at*ct*et;
    }
  });

  team.team_barrier();

  // The outer loop is serial because tau(k) depends on tau(k+1), which eliminates
  // parallelism in the vertical levels. We can still parallelize over pgwvs though.
  for (Int k = src_level; k >= init.ktop; --k) {
    // Determine the diffusivity for each column.
    Real d = GWC::dback;
    if (init.do_molec_diff) {
      d += kvtt(k);
    }
    else {
      Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, -pgwv, pgwv+1), [&] (const int l, Real& lmax) {
        const int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays
        const Real dscal = ekat::impl::min((Real)1.0, tau(pl_idx, k+1) / (tausat(k, pl_idx) + GWC::taumin));
        lmax = ekat::impl::max(lmax, dscal * dsat(k, pl_idx));
      }, Kokkos::Max<Real>(d));
    }

    team.team_barrier();

    // Compute stress for each wave. The stress at this level is the min of
    // the saturation stress and the stress at the level below reduced by
    // damping. The sign of the stress must be the same as at the level below.
    //
    // If molecular diffusion is on, only do this in levels with molecular
    // diffusion. Otherwise, do it everywhere.
    if (k <= init.nbot_molec || !init.do_molec_diff) {
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, -pgwv, pgwv+1), [&] (const int l) {
        const int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays

        const Real wrk = wrk1(k, pl_idx) + wrk2(k, pl_idx) * d;

        Real taudmp;
        if (wrk >= -150.0 || !init.do_molec_diff) {
          taudmp = tau(pl_idx, k+1) * std::exp(wrk);
        } else {
          taudmp = 0.0;
        }
        if (taudmp <= GWC::taumin) taudmp = 0.0;
        tau(pl_idx, k) = ekat::impl::min(taudmp, tausat(k, pl_idx));
      });
    }
    else {
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, -pgwv, pgwv+1), [&] (const int l) {
        int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays
        tau(pl_idx, k) = ekat::impl::min(tau(pl_idx, k+1), tausat(k, pl_idx));
      });
    }
    team.team_barrier();
  }

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<4>(
    {&tausat_1d, &dsat_1d, &wrk1_1d, &wrk2_1d});
}

} // namespace gw
} // namespace scream

#endif
