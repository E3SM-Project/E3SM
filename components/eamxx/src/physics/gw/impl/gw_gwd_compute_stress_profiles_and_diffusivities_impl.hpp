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
  const Int& max_level,
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
  // Local storage
  // (ub-c)
  uview_1d<Real> ubmc, tausat;
  workspace.template take_many_contiguous_unsafe<2>(
    {"ubmc", "tausat"},
    {&ubmc, &tausat});

  // Loop from bottom to top to get stress profiles.
  for (Int k = max_level; k >= init.ktop; --k) {
    if (src_level >= k) {

      // Determine the absolute value of the saturation stress.
      // Define critical levels where the sign of (u-c) changes between interfaces.
      for (Int l = -pgwv; l <= pgwv; ++l) {
        int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays
        ubmc(pl_idx) = ubi(k) - c(pl_idx);

        // Test to see if u-c has the same sign here as the level below.
        if (ubmc(pl_idx) * (ubi(k + 1) - c(pl_idx)) > 0.0) {
          tausat(pl_idx) = std::abs(init.effkwv * rhoi(k) * bfb_cube(ubmc(pl_idx)) /
                                    (2.0 * ni(k)));
          if (tausat(pl_idx) <= GWC::taumin) tausat(pl_idx) = 0.0;
        }
        else {
          tausat(pl_idx) = 0.0;
        }
      }

      // Determine the diffusivity for each column.
      Real d = GWC::dback;
      if (init.do_molec_diff) {
        d += kvtt(k);
      }
      else {
        for (Int l = -pgwv; l <= pgwv; ++l) {
          int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays
          const Real dsat = bfb_square(ubmc(pl_idx) / ni(k)) *
            (init.effkwv * bfb_square(ubmc(pl_idx)) /
             (GWC::rog * ti(k) * ni(k)) - init.alpha(k));
          const Real dscal = std::min(1.0, tau(pl_idx, k+1) / (tausat(pl_idx) + GWC::taumin));
          d = std::max(d, dscal * dsat);
        }
      }

      // Compute stress for each wave. The stress at this level is the min of
      // the saturation stress and the stress at the level below reduced by
      // damping. The sign of the stress must be the same as at the level below.
      //
      // If molecular diffusion is on, only do this in levels with molecular
      // diffusion. Otherwise, do it everywhere.
      if (k <= init.nbot_molec || !init.do_molec_diff) {
        for (Int l = -pgwv; l <= pgwv; ++l) {
          int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays
          const Real ubmc2 = std::max(bfb_square(ubmc(pl_idx)), GWC::ubmc2mn);
          const Real mi = ni(k) / (2.0 * init.kwv * ubmc2) * (init.alpha(k) + bfb_square(ni(k)) / ubmc2 * d);
          const Real wrk = -2.0 * mi * GWC::rog * t(k) * (piln(k + 1) - piln(k));
          Real taudmp;
          if (wrk >= -150.0 || !init.do_molec_diff) {
            taudmp = tau(pl_idx, k+1) * std::exp(wrk);
          } else {
            taudmp = 0.0;
          }
          if (taudmp <= GWC::taumin) taudmp = 0.0;
          tau(pl_idx, k) = std::min(taudmp, tausat(pl_idx));
        }
      }
      else {
        for (Int l = -pgwv; l <= pgwv; ++l) {
          int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays
          tau(pl_idx, k) = std::min(tau(pl_idx, k+1), tausat(pl_idx));
        }
      }
    }
  }

  workspace.template release_many_contiguous<2>(
    {&ubmc, &tausat});
}

} // namespace gw
} // namespace scream

#endif
