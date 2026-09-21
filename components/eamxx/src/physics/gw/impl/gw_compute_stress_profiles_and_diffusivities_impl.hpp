#ifndef GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP
#define GW_GWD_COMPUTE_STRESS_PROFILES_AND_DIFFUSIVITIES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "share/util/eamxx_utils.hpp"

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
    Kokkos::TeamVectorRange(team, (init.ktop+1)*num_pgwv, (src_level+1)*num_pgwv), [&] (const int k_pgwv) {

    const int k = k_pgwv / num_pgwv;
    const int l = k_pgwv % num_pgwv;

    // Determine the absolute value of the saturation stress.
    // Define critical levels where the sign of (u-c) changes between interfaces.
    const Real ubmc = ubi(k) - c(l);

    // Test to see if u-c has the same sign here as the level below.
    if (ubmc * (ubi(k + 1) - c(l)) > 0) {
      tausat(k, l) = std::abs(init.effkwv * rhoi(k) * bfb_cube(ubmc) /
                                   (2 * ni(k)));
      if (tausat(k, l) <= GWC::taumin) tausat(k, l) = 0;
    }
    else {
      tausat(k, l) = 0;
    }

    if (!init.do_molec_diff) {
      dsat(k, l) = bfb_square(ubmc / ni(k)) *
        (init.effkwv * bfb_square(ubmc) /
         (GWC::rog * ti(k) * ni(k)) - init.alpha(k));
    }

    if (k <= init.nbot_molec || !init.do_molec_diff) {
      const Real ubmc2 = ekat::impl::max(bfb_square(ubmc), ubmc2mn);
      const Real at = ni(k) / (2 * init.kwv * ubmc2);
      const Real bt = init.alpha(k);
      const Real ct = bfb_square(ni(k)) / ubmc2;
      const Real et = -2 * GWC::rog * t(k) * (piln(k + 1) - piln(k));
      wrk1(k, l) = at*bt*et;
      wrk2(k, l) = at*ct*et;
    }
  });

  team.team_barrier();

  // The outer loop is serial because tau(k) depends on tau(k+1), which eliminates
  // parallelism in the vertical levels. We can still parallelize over pgwvs though.
  for (Int k = src_level; k > init.ktop; --k) {
    // Determine the diffusivity for each column.
    Real d = GWC::dback;
    if (init.do_molec_diff) {
      d += kvtt(k);
    }
    else {
      Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, num_pgwv), [&] (const int pl_idx, Real& lmax) {
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
        Kokkos::TeamVectorRange(team, num_pgwv), [&] (const int pl_idx) {

        const Real wrk = wrk1(k, pl_idx) + wrk2(k, pl_idx) * d;

        Real taudmp;
        if (wrk >= -150 || !init.do_molec_diff) {
          taudmp = tau(pl_idx, k+1) * std::exp(wrk);
        } else {
          taudmp = 0;
        }
        if (taudmp <= GWC::taumin) taudmp = 0;
        tau(pl_idx, k) = ekat::impl::min(taudmp, tausat(k, pl_idx));
      });
    }
    else {
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, num_pgwv), [&] (const int pl_idx) {
        tau(pl_idx, k) = ekat::impl::min(tau(pl_idx, k+1), tausat(k, pl_idx));
      });
    }
    team.team_barrier();
  }

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<4>(
    {&tausat_1d, &dsat_1d, &wrk1_1d, &wrk2_1d});
}

// Serial version: follows the Fortran gwd_compute_stress_profiles_and_diffusivities
// loop structure exactly. The outer k-loop is strictly serial (tau(k) depends on
// tau(k+1)) and, together with the inner loops over the wave spectrum, runs on a
// single team thread (see the WORKAROUND note below).
//
// Unlike the parallel version above, there is no two-pass precomputation of tausat/
// dsat/wrk1/wrk2 for all levels. Instead, each quantity is computed in-place within
// the single downward-propagating k-loop, matching the Fortran order of operations,
// and nothing is held in team-shared scratch between stages, so this variant needs
// no workspace at all.
template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_compute_stress_profiles_and_diffusivities_serial(
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
  // This variant keeps nothing in team-shared scratch, so the workspace
  // argument is unused; it stays in the signature to match the parallel
  // variant.
  (void)workspace;

  // ---------------------------------------------------------------------------
  // WORKAROUND: this routine is currently called from inside a team-policy
  // parallel_for in run_impl. On the Kokkos/CUDA build used here, a long serial
  // outer loop containing any inner team-collective sync (team.team_barrier(),
  // Kokkos::parallel_reduce on TeamVectorRange, etc.) progressively loses team
  // threads over iterations: thread 0 races ahead while threads 1..N never
  // reach the loop exit. The hang then surfaces at the next team-collective
  // op (typically the workspace release). We verified that all 128 team_ranks
  // pass the barriers at the first iteration but only thread 0 reaches the
  // last iteration, so the loss is cumulative across iterations rather than
  // immediate. Root cause appears to be below the application level (Kokkos
  // / CUDA / driver interaction with team_size=128 and many consecutive
  // __syncthreads() in a serial outer loop) and is not debuggable without
  // CUDA-level tooling.
  //
  // The workaround is to run the entire k-loop on a single team thread inside
  // Kokkos::single(PerTeam), with no team-collective syncs inside the loop, and
  // a single team_barrier after it to publish tau to the rest of the team. This
  // is the same pattern vd_lu_decomp/vd_lu_solve use for their serial
  // recurrences. It preserves column-level parallelism (still one team per
  // column) but gives up intra-column parallelism over the wave spectrum.
  //
  // NOTE: an earlier version of this workaround ran the loop redundantly on
  // every team thread. That is a data race even though each thread computes
  // the same values: threads are not in lockstep, so one thread reads
  // tau(pl_idx, k+1) while another is still writing it, and unsynchronized
  // concurrent read/write of the same location is undefined behavior
  // regardless of the values involved. Every tau element must have exactly one
  // writer.
  //
  // If the underlying stack issue is fixed (Kokkos / CUDA / driver / EKAT
  // team policy), the parallel variant
  // gwd_compute_stress_profiles_and_diffusivities() above can be revisited
  // and dropped back in by changing the call site in gw_drag_prof.
  // ---------------------------------------------------------------------------

  static const auto ubmc2mn = GWC::ubmc2mn;

  const int num_pgwv = 2*pgwv + 1;

  // Saturation stress at interface k for wave pl_idx.
  //
  // Recomputed on demand rather than cached in a workspace array; this keeps
  // the routine free of team-shared scratch. The arithmetic is identical to
  // caching it, so answers are unchanged.
  auto tausat_at = [&] (const Int kk, const int pl_idx) -> Real {
    const Real ubmc = ubi(kk) - c(pl_idx);
    if (ubmc * (ubi(kk + 1) - c(pl_idx)) > 0) {
      const Real ts = Kokkos::abs(init.effkwv * rhoi(kk) * bfb_cube(ubmc) /
                                  (2 * ni(kk)));
      return (ts <= GWC::taumin) ? Real(0) : ts;
    }
    return 0;
  };

  // Single writer for every tau element: see the WORKAROUND note above.
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    // Serial outer loop from the source level upward to the model top.
    // Matches Fortran: do k = maxval(src_level)-1, ktop, -1
    for (Int k = src_level; k > init.ktop; --k) {

      // -------------------------------------------------------------------------
      // Stage 2: Diffusivity d for this level
      // -------------------------------------------------------------------------
      Real d = GWC::dback;
      if (init.do_molec_diff) {
        d += kvtt(k);
      } else {
        for (int pl_idx = 0; pl_idx < num_pgwv; ++pl_idx) {
          const Real ubmc  = ubi(k) - c(pl_idx);
          const Real dsat  = bfb_square(ubmc / ni(k)) *
            (init.effkwv * bfb_square(ubmc) /
             (GWC::rog * ti(k) * ni(k)) - init.alpha(k));
          const Real dscal = ekat::impl::min((Real)1.0,
            tau(pl_idx, k+1) / (tausat_at(k, pl_idx) + GWC::taumin));
          d = ekat::impl::max(d, dscal * dsat);
        }
      }

      // -------------------------------------------------------------------------
      // Stage 3: Stress at interface k
      // -------------------------------------------------------------------------
      if (k <= init.nbot_molec || !init.do_molec_diff) {
        for (int pl_idx = 0; pl_idx < num_pgwv; ++pl_idx) {
          const Real ubmc  = ubi(k) - c(pl_idx);
          const Real ubmc2 = ekat::impl::max(bfb_square(ubmc), ubmc2mn);
          const Real mi    = ni(k) / (2 * init.kwv * ubmc2) *
                             (init.alpha(k) + bfb_square(ni(k)) / ubmc2 * d);
          const Real wrk   = -2 * mi * GWC::rog * t(k) * (piln(k + 1) - piln(k));

          Real taudmp;
          if (wrk >= -150 || !init.do_molec_diff) {
            taudmp = tau(pl_idx, k+1) * Kokkos::exp(wrk);
          } else {
            taudmp = 0;
          }

          if (taudmp <= GWC::taumin) taudmp = 0;
          tau(pl_idx, k) = ekat::impl::min(taudmp, tausat_at(k, pl_idx));
        }
      } else {
        for (int pl_idx = 0; pl_idx < num_pgwv; ++pl_idx) {
          tau(pl_idx, k) = ekat::impl::min(tau(pl_idx, k+1), tausat_at(k, pl_idx));
        }
      }
    }
  });

  team.team_barrier();
}

} // namespace gw
} // namespace scream

#endif
