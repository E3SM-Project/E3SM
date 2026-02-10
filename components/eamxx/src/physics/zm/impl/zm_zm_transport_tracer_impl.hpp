#ifndef ZM_ZM_TRANSPORT_TRACER_IMPL_HPP
#define ZM_ZM_TRANSPORT_TRACER_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_transport_tracer. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_transport_tracer(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver,                        // number of mid-point levels
  const uview_1d<const bool>& doconvtran, // flag for doing convective transport
  const uview_2d<const Real>& q,          // tracer array (including water vapor)
  const Int& ncnst,                       // number of tracers to transport
  const uview_1d<const Real>& mu,         // mass flux up
  const uview_1d<const Real>& md,         // mass flux down
  const uview_1d<const Real>& du,         // mass detraining from updraft
  const uview_1d<const Real>& eu,         // mass entraining from updraft
  const uview_1d<const Real>& ed,         // mass entraining from downdraft
  const uview_1d<const Real>& dp,         // delta pressure between interfaces
  const Int& jt,                          // index of cloud top for this column
  const Int& mx,                          // index of cloud bottom for this column
  const Int& ktm,                         // Highest top level for any column
  const Int& kbm,                         // Highest bottom level for any column
  const uview_2d<const Real>& fracis,     // fraction of tracer that is insoluble
  const uview_1d<const Real>& dpdry,      // delta pressure between interfaces
  const Real& dt,                         // model time increment)
  // Outputs
  const uview_2d<Real>& dqdt)             // output tendency array
{
  // Allocate temporary arrays (no pcols dimension)
  uview_1d<Real> chat, cond, const_arr, fisg, conu, dcondt, dutmp, eutmp, edtmp, dptmp;
  workspace.template take_many_contiguous_unsafe<10>(
    {"chat", "cond", "const_arr", "fisg", "conu", "dcondt", "dutmp", "eutmp", "edtmp", "dptmp"},
    {&chat, &cond, &const_arr, &fisg, &conu, &dcondt, &dutmp, &eutmp, &edtmp, &dptmp});

  // Loop over each constituent (skip water vapor at m=1, Fortran indexing)
  for (Int m = 1; m < ncnst; ++m) {
    if (!doconvtran(m)) continue;

    // Initialize temporary arrays (always use moist formulation)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
      // if (cnst_get_type_byind(m).eq.'dry') then
      if (false) { // dry, assume false for now
        const Real ratio = dp(k)/dpdry(k);
        dptmp(k) = dpdry(k);
        dutmp(k) = du(k) * ratio;
        eutmp(k) = eu(k) * ratio;
        edtmp(k) = ed(k) * ratio;
      }
      else {
        dptmp(k) = dp(k);
        dutmp(k) = du(k);
        eutmp(k) = eu(k);
        edtmp(k) = ed(k);
      }
      // Gather up the constituent and set tend to zero
      const_arr(k) = q(k, m);
      fisg(k) = fracis(k, m);
    });
    team.team_barrier();

    // From now on work only with gathered data

    // Interpolate environment tracer values to interfaces
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
      const Int km1 = ekat::impl::max(0, k-1);
      const Real minc = ekat::impl::min(const_arr(km1), const_arr(k));
      const Real maxc = ekat::impl::max(const_arr(km1), const_arr(k));

      Real cdifr;
      if (minc < 0) {
        cdifr = 0.0;
      } else {
        cdifr = std::abs(const_arr(k) - const_arr(km1)) / ekat::impl::max(maxc, ZMC::small);
      }

      // If the two layers differ significantly use a geometric averaging
      if (cdifr > ZMC::cdifr_min) {
        const Real cabv = ekat::impl::max(const_arr(km1), maxc*ZMC::maxc_factor);
        const Real cbel = ekat::impl::max(const_arr(k), maxc*ZMC::maxc_factor);
        chat(k) = std::log(cabv/cbel) / (cabv-cbel) * cabv * cbel;
      } else { // Small diff, so just arithmetic mean
        chat(k) = ZMC::half * (const_arr(k) + const_arr(km1));
      }

      // Provisional updraft and downdraft values
      conu(k) = chat(k);
      cond(k) = chat(k);
      // provisional tendencies
      dcondt(k) = 0.0;
    });
    team.team_barrier();

    // Do levels adjacent to top and bottom
    const Int kk = pver - 1;
    Kokkos::single(Kokkos::PerTeam(team), [&] {
      const Real mupdudp = mu(kk) + dutmp(kk) * dptmp(kk);
      if (mupdudp > ZMC::mbsth) {
        conu(kk) = (eutmp(kk) * fisg(kk) * const_arr(kk) * dptmp(kk)) / mupdudp;
      }
      if (md(1) < -ZMC::mbsth) {
        cond(1) = (-edtmp(0) * fisg(0) * const_arr(0) * dptmp(0)) / md(1);
      }
    });
    team.team_barrier();

    // Updraft from bottom to top
    Kokkos::single(Kokkos::PerTeam(team), [&] {
      for (Int kk = pver-2; kk >= 0; --kk) {
        const Int kkp1 = ekat::impl::min(pver-1, kk+1);
        const Real mupdudp = mu(kk) + dutmp(kk) * dptmp(kk);
        if (mupdudp > ZMC::mbsth) {
          conu(kk) = (mu(kkp1) * conu(kkp1) + eutmp(kk) * fisg(kk) * const_arr(kk) * dptmp(kk)) / mupdudp;
        }
      }
    });

    // Downdraft from top to bottom
    Kokkos::single(Kokkos::PerTeam(team), [&] {
      for (Int k = 2; k < pver; ++k) {
        const Int km1 = ekat::impl::max(0, k-1);
        if (md(k) < -ZMC::mbsth) {
          cond(k) = (md(km1) * cond(km1) - edtmp(km1) * fisg(km1) * const_arr(km1) * dptmp(km1)) / md(k);
        }
      }
    });
    team.team_barrier();

    // Compute tendencies from cloud top to bottom
    Kokkos::single(Kokkos::PerTeam(team), [&] {
      for (Int k = ktm; k < pver; ++k) {
        const Int km1 = ekat::impl::max(0, k-1);
        const Int kp1 = ekat::impl::min(pver-1, k+1);

        const Real fluxin = mu(kp1) * conu(kp1) + mu(k) * ekat::impl::min(chat(k), const_arr(km1))
                          - (md(k) * cond(k) + md(kp1) * ekat::impl::min(chat(kp1), const_arr(kp1)));
        const Real fluxout = mu(k) * conu(k) + mu(kp1) * ekat::impl::min(chat(kp1), const_arr(k))
                           - (md(kp1) * cond(kp1) + md(k) * ekat::impl::min(chat(k), const_arr(k)));

        Real netflux = fluxin - fluxout;
        if (std::abs(netflux) < ekat::impl::max(fluxin, fluxout) * ZMC::flux_factor) {
          netflux = 0.0;
        }
        dcondt(k) = netflux / dptmp(k);
      }
    });
    team.team_barrier();

    // Handle cloud base levels
    Kokkos::single(Kokkos::PerTeam(team), [&] {
      for (Int k = kbm; k < pver; ++k) {
        const Int km1 = ekat::impl::max(0, k-1);
        if (k == mx) {
          // limit fluxes outside convection to mass in appropriate layer
          // these limiters are probably only safe for positive definite quantitities
          // it assumes that mu and md already satify a courant number limit of 1
          const Real fluxin = mu(k) * ekat::impl::min(chat(k), const_arr(km1)) - md(k) * cond(k);
          const Real fluxout = mu(k) * conu(k) - md(k) * ekat::impl::min(chat(k), const_arr(k));
          Real netflux = fluxin - fluxout;
          if (std::abs(netflux) < ekat::impl::max(fluxin, fluxout) * ZMC::flux_factor) {
            netflux = 0.0;
          }
          dcondt(k) = netflux / dptmp(k);
        } else if (k > mx) {
          dcondt(k) = 0.0;
        }
      }
    });
    team.team_barrier();

    // Conservation check for ZM microphysics
    if (runtime_opt.zm_microp) {
      Kokkos::single(Kokkos::PerTeam(team), [&] {
        for (Int k = jt; k <= mx; ++k) {
          if (dcondt(k) * dt + const_arr(k) < 0.0) {
            Real negadt = dcondt(k) + const_arr(k) / dt;
            dcondt(k) = -const_arr(k) / dt;

            // Try to redistribute to levels above (k+1 to mx)
            for (Int kk = k+1; kk <= mx; ++kk) {
              if (negadt < 0.0 && dcondt(kk) * dt + const_arr(kk) > 0.0) {
                Real qtmp = dcondt(kk) + negadt * dptmp(k) / dptmp(kk);
                if (qtmp * dt + const_arr(kk) > 0.0) {
                  dcondt(kk) = qtmp;
                  negadt = 0.0;
                } else {
                  negadt = negadt + (const_arr(kk) / dt + dcondt(kk)) * dptmp(kk) / dptmp(k);
                  dcondt(kk) = -const_arr(kk) / dt;
                }
              }
            }

            // Try to redistribute to levels below (k-1 down to jt)
            for (Int kk = k-1; kk >= jt; --kk) {
              if (negadt < 0.0 && dcondt(kk) * dt + const_arr(kk) > 0.0) {
                Real qtmp = dcondt(kk) + negadt * dptmp(k) / dptmp(kk);
                if (qtmp * dt + const_arr(kk) > 0.0) {
                  dcondt(kk) = qtmp;
                  negadt = 0.0;
                } else {
                  negadt = negadt + (const_arr(kk) / dt + dcondt(kk)) * dptmp(kk) / dptmp(k);
                  dcondt(kk) = -const_arr(kk) / dt;
                }
              }
            }

            // If still negative, add back to current level
            if (negadt < 0.0) {
              dcondt(k) = dcondt(k) - negadt;
            }
          }
        }
      });
      team.team_barrier();
    }

    // Scatter tendency back to output array
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
      dqdt(k, m) = dcondt(k);
    });
    team.team_barrier();
  }

  workspace.template release_many_contiguous<10>(
    {&chat, &cond, &const_arr, &fisg, &conu, &dcondt, &dutmp, &eutmp, &edtmp, &dptmp});
}

} // namespace zm
} // namespace scream

#endif
