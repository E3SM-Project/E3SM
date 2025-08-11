#ifndef GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP
#define GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

#define bfb_square(val) ((val)*(val))
#define bfb_cube(val)   ((val)*(val)*(val))

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
  const uview_2d<Real>& tausat,
  const uview_2d<Real>& dsat,
  const uview_2d<Real>& wrk1,
  const uview_2d<Real>& wrk2,
  const uview_2d<Real>& tau)
{
  static const auto ubmc2mn = GWC::ubmc2mn;

  // Loop from bottom to top to get stress profiles.
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, init.ktop, src_level+1), [&] (const int k) {

    // Determine the absolute value of the saturation stress.
    // Define critical levels where the sign of (u-c) changes between interfaces.
    Kokkos::parallel_for(
      Kokkos::ThreadVectorRange(team, -pgwv, pgwv+1), [&] (const int l) {

      const int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays
      const Real ubmc = ubi(k) - c(pl_idx);

      // Test to see if u-c has the same sign here as the level below.
      if (ubmc * (ubi(k + 1) - c(pl_idx)) > 0.0) {
        tausat(k, pl_idx) = std::abs(init.effkwv * rhoi(k) * bfb_cube(ubmc) /
                                     (2.0 * ni(k)));
        if (tausat(k, pl_idx) <= GWC::taumin) tausat(k, pl_idx) = 0.0;
      }
      else {
        tausat(k, pl_idx) = 0.0;
      }

      if (!init.do_molec_diff) {
        dsat(k, pl_idx) = bfb_square(ubmc / ni(k)) *
          (init.effkwv * bfb_square(ubmc) /
           (GWC::rog * ti(k) * ni(k)) - init.alpha(k));
      }

      if (k <= init.nbot_molec || !init.do_molec_diff) {
        const Real ubmc2 = ekat::impl::max(bfb_square(ubmc), ubmc2mn);
        const Real at = ni(k) / (2.0 * init.kwv * ubmc2);
        const Real bt = init.alpha(k);
        const Real ct = bfb_square(ni(k)) / ubmc2;
        const Real et = -2.0 * GWC::rog * t(k) * (piln(k + 1) - piln(k));
        wrk1(k, pl_idx) = at*bt*et;
        wrk2(k, pl_idx) = at*ct*et;
      }
    });
  });

  team.team_barrier();

  // This loop is serial, so it may as well only be performed by one thread.
  // tau(k) depends on tau(k+1), which eliminates parallelism in the vertical
  // levels.
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
        const Real dscal = ekat::impl::min(1.0, tau(pl_idx, k+1) / (tausat(k, pl_idx) + GWC::taumin));
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
}

} // namespace gw
} // namespace scream

#endif
