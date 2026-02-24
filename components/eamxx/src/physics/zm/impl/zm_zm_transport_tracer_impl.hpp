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
  constexpr Real small = ZMC::small;
  // Allocate temporary arrays (2D: ncnst x pver)
  uview_1d<Real> chat1d, cond1d, const_arr1d, fisg1d, conu1d, dcondt1d, dutmp1d, eutmp1d, edtmp1d, dptmp1d;
  workspace.template take_many_contiguous_unsafe<10>(
    {"chat", "cond", "const_arr", "fisg", "conu", "dcondt", "dutmp", "eutmp", "edtmp", "dptmp"},
    {&chat1d, &cond1d, &const_arr1d, &fisg1d, &conu1d, &dcondt1d, &dutmp1d, &eutmp1d, &edtmp1d, &dptmp1d});

  uview_2d<Real>
    chat(chat1d.data(), ncnst, pver),
    cond(cond1d.data(), ncnst, pver),
    const_arr(const_arr1d.data(), ncnst, pver),
    fisg(fisg1d.data(), ncnst, pver),
    conu(conu1d.data(), ncnst, pver),
    dcondt(dcondt1d.data(), ncnst, pver),
    dutmp(dutmp1d.data(), ncnst, pver),
    eutmp(eutmp1d.data(), ncnst, pver),
    edtmp(edtmp1d.data(), ncnst, pver),
    dptmp(dptmp1d.data(), ncnst, pver);

  // Parallel loop over each constituent (skip water vapor at m=0)
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 1, ncnst), [&] (const Int& m) {
    if (!doconvtran(m)) return;

    // Initialize temporary arrays (always use moist formulation)
    for (Int k = 0; k < pver; ++k) {
      // if (cnst_get_type_byind(m).eq.'dry') then
      if (false) { // dry, assume false for now
        const Real ratio = dp(k)/dpdry(k);
        dptmp(m, k) = dpdry(k);
        dutmp(m, k) = du(k) * ratio;
        eutmp(m, k) = eu(k) * ratio;
        edtmp(m, k) = ed(k) * ratio;
      }
      else {
        dptmp(m, k) = dp(k);
        dutmp(m, k) = du(k);
        eutmp(m, k) = eu(k);
        edtmp(m, k) = ed(k);
      }
      // Gather up the constituent and set tend to zero
      const_arr(m, k) = q(k, m);
      fisg(m, k) = fracis(k, m);
    }

    // From now on work only with gathered data

    // Interpolate environment tracer values to interfaces
    for (Int k = 0; k < pver; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      const Real minc = ekat::impl::min(const_arr(m, km1), const_arr(m, k));
      const Real maxc = ekat::impl::max(const_arr(m, km1), const_arr(m, k));

      Real cdifr;
      if (minc < 0) {
        cdifr = 0.0;
      } else {
        cdifr = std::abs(const_arr(m, k) - const_arr(m, km1)) / ekat::impl::max(maxc, small);
      }

      // If the two layers differ significantly use a geometric averaging
      if (cdifr > ZMC::cdifr_min) {
        const Real cabv = ekat::impl::max(const_arr(m, km1), maxc*ZMC::maxc_factor);
        const Real cbel = ekat::impl::max(const_arr(m, k), maxc*ZMC::maxc_factor);
        chat(m, k) = std::log(cabv/cbel) / (cabv-cbel) * cabv * cbel;
      } else { // Small diff, so just arithmetic mean
        chat(m, k) = ZMC::half * (const_arr(m, k) + const_arr(m, km1));
      }

      // Provisional updraft and downdraft values
      conu(m, k) = chat(m, k);
      cond(m, k) = chat(m, k);
      // provisional tendencies
      dcondt(m, k) = 0.0;
    }

    // Do levels adjacent to top and bottom
    const Int kk = pver - 1;
    const Real mupdudp = mu(kk) + dutmp(m, kk) * dptmp(m, kk);
    if (mupdudp > ZMC::mbsth) {
      conu(m, kk) = (eutmp(m, kk) * fisg(m, kk) * const_arr(m, kk) * dptmp(m, kk)) / mupdudp;
    }
    if (md(1) < -ZMC::mbsth) {
      cond(m, 1) = (-edtmp(m, 0) * fisg(m, 0) * const_arr(m, 0) * dptmp(m, 0)) / md(1);
    }

    // Updraft from bottom to top
    for (Int kk = pver-2; kk >= 0; --kk) {
      const Int kkp1 = ekat::impl::min(pver-1, kk+1);
      const Real mupdudp_kk = mu(kk) + dutmp(m, kk) * dptmp(m, kk);
      if (mupdudp_kk > ZMC::mbsth) {
        conu(m, kk) = (mu(kkp1) * conu(m, kkp1) + eutmp(m, kk) * fisg(m, kk) * const_arr(m, kk) * dptmp(m, kk)) / mupdudp_kk;
      }
    }

    // Downdraft from top to bottom
    for (Int k = 2; k < pver; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      if (md(k) < -ZMC::mbsth) {
        cond(m, k) = (md(km1) * cond(m, km1) - edtmp(m, km1) * fisg(m, km1) * const_arr(m, km1) * dptmp(m, km1)) / md(k);
      }
    }

    // Compute tendencies from cloud top to bottom
    for (Int k = ktm; k < pver; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      const Int kp1 = ekat::impl::min(pver-1, k+1);

      const Real fluxin = mu(kp1) * conu(m, kp1) + mu(k) * ekat::impl::min(chat(m, k), const_arr(m, km1))
                        - (md(k) * cond(m, k) + md(kp1) * ekat::impl::min(chat(m, kp1), const_arr(m, kp1)));
      const Real fluxout = mu(k) * conu(m, k) + mu(kp1) * ekat::impl::min(chat(m, kp1), const_arr(m, k))
                         - (md(kp1) * cond(m, kp1) + md(k) * ekat::impl::min(chat(m, k), const_arr(m, k)));

      Real netflux = fluxin - fluxout;
      if (std::abs(netflux) < ekat::impl::max(fluxin, fluxout) * ZMC::flux_factor) {
        netflux = 0.0;
      }
      dcondt(m, k) = netflux / dptmp(m, k);
    }

    // Handle cloud base levels
    for (Int k = kbm; k < pver; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      if (k == mx) {
        // limit fluxes outside convection to mass in appropriate layer
        // these limiters are probably only safe for positive definite quantitities
        // it assumes that mu and md already satify a courant number limit of 1
        const Real fluxin = mu(k) * ekat::impl::min(chat(m, k), const_arr(m, km1)) - md(k) * cond(m, k);
        const Real fluxout = mu(k) * conu(m, k) - md(k) * ekat::impl::min(chat(m, k), const_arr(m, k));
        Real netflux = fluxin - fluxout;
        if (std::abs(netflux) < ekat::impl::max(fluxin, fluxout) * ZMC::flux_factor) {
          netflux = 0.0;
        }
        dcondt(m, k) = netflux / dptmp(m, k);
      } else if (k > mx) {
        dcondt(m, k) = 0.0;
      }
    }

    // Conservation check for ZM microphysics
    if (runtime_opt.zm_microp) {
      for (Int k = jt; k <= mx; ++k) {
        if (dcondt(m, k) * dt + const_arr(m, k) < 0.0) {
          Real negadt = dcondt(m, k) + const_arr(m, k) / dt;
          dcondt(m, k) = -const_arr(m, k) / dt;

          // Try to redistribute to levels above (k+1 to mx)
          for (Int kk = k+1; kk <= mx; ++kk) {
            if (negadt < 0.0 && dcondt(m, kk) * dt + const_arr(m, kk) > 0.0) {
              const Real qtmp = dcondt(m, kk) + negadt * dptmp(m, k) / dptmp(m, kk);
              if (qtmp * dt + const_arr(m, kk) > 0.0) {
                dcondt(m, kk) = qtmp;
                negadt = 0.0;
              } else {
                negadt = negadt + (const_arr(m, kk) / dt + dcondt(m, kk)) * dptmp(m, kk) / dptmp(m, k);
                dcondt(m, kk) = -const_arr(m, kk) / dt;
              }
            }
          }

          // Try to redistribute to levels below (k-1 down to jt)
          for (Int kk = k-1; kk >= jt; --kk) {
            if (negadt < 0.0 && (dcondt(m, kk) * dt + const_arr(m, kk) > 0.0)) {
              const Real qtmp = dcondt(m, kk) + negadt * dptmp(m, k) / dptmp(m, kk);
              if (qtmp * dt + const_arr(m, kk) > 0.0) {
                dcondt(m, kk) = qtmp;
                negadt = 0.0;
              } else {
                negadt = negadt + (const_arr(m, kk) / dt + dcondt(m, kk)) * dptmp(m, kk) / dptmp(m, k);
                dcondt(m, kk) = -const_arr(m, kk) / dt;
              }
            }
          }

          // If still negative, add back to current level
          if (negadt < 0.0) {
            dcondt(m, k) = dcondt(m, k) - negadt;
          }
        }
      }
    }

    // Scatter tendency back to output array
    for (Int k = 0; k < pver; ++k) {
      dqdt(k, m) = dcondt(m, k);
    }
  });
  team.team_barrier();

  workspace.template release_many_contiguous<10>(
    {&chat1d, &cond1d, &const_arr1d, &fisg1d, &conu1d, &dcondt1d, &dutmp1d, &eutmp1d, &edtmp1d, &dptmp1d});
}

} // namespace zm
} // namespace scream

#endif
