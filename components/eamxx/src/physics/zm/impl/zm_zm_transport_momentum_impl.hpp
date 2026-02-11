#ifndef ZM_ZM_TRANSPORT_MOMENTUM_IMPL_HPP
#define ZM_ZM_TRANSPORT_MOMENTUM_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_transport_momentum. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_transport_momentum(
  // Inputs
  const MemberType& team,
  const Int& ncol, // actual number of columns
  const Int& pver, // number of mid-point levels
  const Int& pverp, // number of interface levels
  const uview_2d<const Real>& wind_in, // input Momentum array
  const Int& nwind, // number of tracers to transport
  const uview_1d<const Real>& mu, // mass flux up
  const uview_1d<const Real>& md, // mass flux down
  const uview_1d<const Real>& du, // mass detraining from updraft
  const uview_1d<const Real>& eu, // mass entraining from updraft
  const uview_1d<const Real>& ed, // mass entraining from downdraft
  const uview_1d<const Real>& dp, // gathered pressure delta between interfaces
  const Int& jt, // index of cloud top for each column
  const Int& mx, // index of cloud top for each column
  const Int& ideep, // gathering array
  const Int& il1g, // gathered min ncol index
  const Int& il2g, // gathered max ncol index
  const Real& dt, // time step in seconds : 2*delta_t
  // Outputs
  const uview_2d<Real>& wind_tend, // output momentum tendency
  const uview_2d<Real>& pguall, // apparent force from  updraft PG
  const uview_2d<Real>& pgdall, // apparent force from  downdraft PG
  const uview_2d<Real>& icwu, // in-cloud winds in updraft
  const uview_2d<Real>& icwd, // in-cloud winds in downdraft
  const uview_1d<Real>& seten) // dry energy tendency)
{
  constexpr Real momcu = 0.4; // pressure gradient term constant for updrafts
  constexpr Real momcd = 0.4; // pressure gradient term constant for downdrafts

  // Allocate temporary arrays (2D: nwind x pver)
  uview_2d<Real> wind0, windf, wind_mid, wind_int, wind_int_d, wind_int_u, wind_tend_tmp;
  uview_2d<Real> mududp, mddudp, pgu, pgd, gseten;
  uview_3d<Real> mflux;

  workspace.template take_many_contiguous_unsafe<11>(
    {"wind0", "windf", "wind_mid", "wind_int", "wind_int_d", "wind_int_u", "wind_tend_tmp",
     "mududp", "mddudp", "pgu", "pgd"},
    {&wind0, &windf, &wind_mid, &wind_int, &wind_int_d, &wind_int_u, &wind_tend_tmp,
     &mududp, &mddudp, &pgu, &pgd});
  workspace.template take_many_contiguous_unsafe<2>(
    {"gseten", "mflux"},
    {&gseten, &mflux});

  // Initialize gseten and seten to zero
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    gseten(0, k) = 0.0;
    seten(k) = 0.0;
  });
  team.team_barrier();

  // Parallel loop over each wind component
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nwind), [&] (const Int& m) {

    // Initialize in-cloud winds to environmental wind and gather
    for (Int k = 0; k < pver; ++k) {
      wind_mid(m, k) = wind_in(k, m);
      wind0(m, k) = wind_mid(m, k);
      icwu(k, m) = wind_in(k, m);
      icwd(k, m) = wind_in(k, m);
      pguall(k, m) = 0.0;
      pgdall(k, m) = 0.0;
    }

    for (Int k = 0; k <= pver; ++k) {
      mflux(m, k, 0) = 0.0;
    }

    // Interpolate winds to interfaces
    for (Int k = 0; k < pver; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      // use arithmetic mean
      wind_int(m, k) = ZMC::half * (wind_mid(m, k) + wind_mid(m, km1));
      // Provisional up and down draft values
      wind_int_u(m, k) = wind_int(m, k);
      wind_int_d(m, k) = wind_int(m, k);
      // provisional tendency
      wind_tend_tmp(m, k) = 0.0;
    }

    // Calculate pressure perturbation terms

    // upper boundary - assume mu is zero
    pgu(m, 0) = 0.0;
    pgd(m, 0) = 0.0;

    // interior points
    for (Int k = 1; k < pver-1; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      const Int kp1 = ekat::impl::min(pver-1, k+1);

      mududp(m, k) = (mu(k)   * (wind_mid(m, k)   - wind_mid(m, km1)) / dp(km1)
                    + mu(kp1) * (wind_mid(m, kp1) - wind_mid(m, k))   / dp(k));
      mddudp(m, k) = (md(k)   * (wind_mid(m, k)   - wind_mid(m, km1)) / dp(km1)
                    + md(kp1) * (wind_mid(m, kp1) - wind_mid(m, k))   / dp(k));
      pgu(m, k) = -momcu * ZMC::half * mududp(m, k);
      pgd(m, k) = -momcd * ZMC::half * mddudp(m, k);
    }

    // bottom boundary
    {
      const Int k = pver - 1;
      const Int km1 = ekat::impl::max(0, k-1);
      mududp(m, k) = mu(k) * (wind_mid(m, k) - wind_mid(m, km1)) / dp(km1);
      mddudp(m, k) = md(k) * (wind_mid(m, k) - wind_mid(m, km1)) / dp(km1);
      pgu(m, k) = -momcu * mududp(m, k);
      pgd(m, k) = -momcd * mddudp(m, k);
    }

    // Calculate in-cloud velocity

    // levels adjacent to top and bottom
    {
      const Int kk = pver - 1;
      const Int kkm1 = ekat::impl::max(0, kk-1);
      Real mupdudp = mu(kk) + du(kk) * dp(kk);
      if (mupdudp > ZMC::mbsth) {
        wind_int_u(m, kk) = (eu(kk) * wind_mid(m, kk) * dp(kk) + pgu(m, kk) * dp(kk)) / mupdudp;
      }
      if (md(1) < -ZMC::mbsth) {
        wind_int_d(m, 1) = (-ed(0) * wind_mid(m, 0) * dp(0) - pgd(m, 0) * dp(0)) / md(1);
      }
    }

    // Updraft from bottom to top
    for (Int kk = pver-2; kk >= 0; --kk) {
      const Int kkp1 = ekat::impl::min(pver-1, kk+1);
      Real mupdudp = mu(kk) + du(kk) * dp(kk);
      if (mupdudp > ZMC::mbsth) {
        wind_int_u(m, kk) = (mu(kkp1) * wind_int_u(m, kkp1) + eu(kk) * wind_mid(m, kk) * dp(kk)
                            + pgu(m, kk) * dp(kk)) / mupdudp;
      }
    }

    // Downdraft from top to bottom
    for (Int k = 2; k < pver; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      if (md(k) < -ZMC::mbsth) {
        wind_int_d(m, k) = (md(km1) * wind_int_d(m, km1) - ed(km1) * wind_mid(m, km1) * dp(km1)
                           - pgd(m, km1) * dp(km1)) / md(k);
      }
    }

    // Calculate momentum tendency
    for (Int k = ktm; k < pver; ++k) {
      const Int kp1 = ekat::impl::min(pver-1, k+1);
      wind_tend_tmp(m, k) = (mu(kp1) * (wind_int_u(m, kp1) - wind_int(m, kp1))
                           - mu(k)   * (wind_int_u(m, k)   - wind_int(m, k))
                           + md(kp1) * (wind_int_d(m, kp1) - wind_int(m, kp1))
                           - md(k)   * (wind_int_d(m, k)   - wind_int(m, k))) / dp(k);
    }

    // tendency for bottom layer
    for (Int k = kbm; k < pver; ++k) {
      if (k == mx) {
        wind_tend_tmp(m, k) = (-mu(k) * (wind_int_u(m, k) - wind_int(m, k))
                              - md(k) * (wind_int_d(m, k) - wind_int(m, k))) / dp(k);
      }
    }

    // Scatter tendency and related outputs back to full array
    for (Int k = 0; k < pver; ++k) {
      wind_tend(k, m) = wind_tend_tmp(m, k);
      pguall(k, m) = -pgu(m, k);
      pgdall(k, m) = -pgd(m, k);
      icwu(k, m) = wind_int_u(m, k);
      icwd(k, m) = wind_int_d(m, k);
    }

    // Calculate momentum flux in units of mb*m/s2
    for (Int k = ktm; k < pver; ++k) {
      mflux(m, k, 0) = -mu(k) * (wind_int_u(m, k) - wind_int(m, k))
                       -md(k) * (wind_int_d(m, k) - wind_int(m, k));
    }

    // Calculate winds at the end of the time step
    for (Int k = ktm; k < pver; ++k) {
      const Int kp1 = k + 1;
      windf(m, k) = wind_mid(m, k) - (mflux(m, kp1, 0) - mflux(m, k, 0)) * dt / dp(k);
    }
  });
  team.team_barrier();

  // Need to add an energy fix to account for the dissipation of kinetic energy
  // Formulation follows from Boville and Bretherton (2003) - modified by Phil Rasch
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    for (Int k = ktm; k < pver; ++k) {
      const Int km1 = ekat::impl::max(0, k-1);
      const Int kp1 = ekat::impl::min(pver-1, k+1);

      // calculate the KE fluxes at top and bot of layer
      // based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
      const Real utop = (wind0(0, k) + wind0(0, km1)) * ZMC::half;
      const Real vtop = (wind0(1, k) + wind0(1, km1)) * ZMC::half;
      const Real ubot = (wind0(0, kp1) + wind0(0, k)) * ZMC::half;
      const Real vbot = (wind0(1, kp1) + wind0(1, k)) * ZMC::half;
      const Real fket = utop * mflux(0, k, 0) + vtop * mflux(1, k, 0);   // top of layer
      const Real fkeb = ubot * mflux(0, k+1, 0) + vbot * mflux(1, k+1, 0); // bot of layer

      // divergence of these fluxes should give a conservative redistribution of KE
      const Real ketend_cons = (fket - fkeb) / dp(k);

      // tendency in kinetic energy resulting from the momentum transport
      const Real ketend = ((windf(0, k) * windf(0, k) + windf(1, k) * windf(1, k))
                         - (wind0(0, k) * wind0(0, k) + wind0(1, k) * wind0(1, k))) * ZMC::half / dt;

      // the difference should be the dissipation
      gseten(0, k) = ketend_cons - ketend;
      seten(k) = gseten(0, k);
    }
  });
  team.team_barrier();

  workspace.template release_many_contiguous<11>(
    {&wind0, &windf, &wind_mid, &wind_int, &wind_int_d, &wind_int_u, &wind_tend_tmp,
     &mududp, &mddudp, &pgu, &pgd});
  workspace.template release_many_contiguous<2>(
    {&gseten, &mflux});
}

} // namespace zm
} // namespace scream

#endif
