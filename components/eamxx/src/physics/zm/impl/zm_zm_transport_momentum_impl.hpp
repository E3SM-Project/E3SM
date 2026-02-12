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
  const Workspace& workspace,
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
  const Int& ktm, // Highest top level for any column
  const Int& kbm, // Highest bottom level for any column
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

  // Allocate temporary arrays (2D: pver x nwind or pverp x nwind)
  uview_1d<Real> wind0_1d, windf_1d, wind_mid_1d, wind_int_1d, wind_int_d_1d, wind_int_u_1d;
  uview_1d<Real> wind_tend_tmp_1d, mududp_1d, mddudp_1d, pgu_1d, pgd_1d, mflux_1d;
  workspace.template take_many_contiguous_unsafe<12>(
    {"wind0", "windf", "wind_mid", "wind_int", "wind_int_d", "wind_int_u",
     "wind_tend_tmp", "mududp", "mddudp", "pgu", "pgd", "mflux"},
    {&wind0_1d, &windf_1d, &wind_mid_1d, &wind_int_1d, &wind_int_d_1d, &wind_int_u_1d,
     &wind_tend_tmp_1d, &mududp_1d, &mddudp_1d, &pgu_1d, &pgd_1d, &mflux_1d});

  uview_2d<Real>
    wind0(wind0_1d.data(), nwind, pver),
    windf(windf_1d.data(), nwind, pver),
    wind_mid(wind_mid_1d.data(), nwind, pver),
    wind_int(wind_int_1d.data(), nwind, pver),
    wind_int_d(wind_int_d_1d.data(), nwind, pver),
    wind_int_u(wind_int_u_1d.data(), nwind, pver),
    wind_tend_tmp(wind_tend_tmp_1d.data(), nwind, pver),
    mududp(mududp_1d.data(), nwind, pver),
    mddudp(mddudp_1d.data(), nwind, pver),
    pgu(pgu_1d.data(), nwind, pver),
    pgd(pgd_1d.data(), nwind, pver),
    mflux(mflux_1d.data(), nwind, pverp);

  // Initialize outputs
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver*nwind), [&] (const Int& idx) {
    const Int k = idx / nwind;
    const Int m = idx % nwind;
    wind_tend(k,m) = 0.0;
    pguall(k,m) = 0.0;
    pgdall(k,m) = 0.0;
    icwu(k,m) = wind_in(k,m);
    icwd(k,m) = wind_in(k,m);
    wind0(m,k) = 0.0;
    windf(m,k) = 0.0;
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    seten(k) = 0.0;
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pverp*nwind), [&] (const Int& idx) {
    mflux_1d(idx) = 0.0;
  });

  team.team_barrier();

  // Loop over each wind component using team parallelism
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nwind), [&] (const Int& m) {

    // Gather up the winds
    for (Int k = 0; k < pver; ++k) {
      wind_mid(m,k) = wind_in(k,m);
      wind0(m,k) = wind_mid(m,k);
    }

    // Interpolate winds to interfaces
    for (Int k = 0; k < pver; ++k) {
      const Int km1 = (k == 0) ? 0 : k-1;
      wind_int(m,k) = 0.5 * (wind_mid(m,k) + wind_mid(m,km1));
      wind_int_u(m,k) = wind_int(m,k);
      wind_int_d(m,k) = wind_int(m,k);
      wind_tend_tmp(m,k) = 0.0;
    }

    // Calculate pressure perturbation terms
    // upper boundary
    pgu(m,0) = 0.0;
    pgd(m,0) = 0.0;

    // interior points
    for (Int k = 1; k < pver-1; ++k) {
      const Int km1 = k-1;
      const Int kp1 = k+1;
      mududp(m,k) = (mu(k)  *(wind_mid(m,k)  -wind_mid(m,km1))/dp(km1) +
                     mu(kp1)*(wind_mid(m,kp1)-wind_mid(m,k)  )/dp(k));
      mddudp(m,k) = (md(k)  *(wind_mid(m,k)  -wind_mid(m,km1))/dp(km1) +
                     md(kp1)*(wind_mid(m,kp1)-wind_mid(m,k)  )/dp(k));
      pgu(m,k) = -momcu * 0.5 * mududp(m,k);
      pgd(m,k) = -momcd * 0.5 * mddudp(m,k);
    }

    // bottom boundary
    {
      const Int k = pver-1;
      const Int km1 = k-1;
      mududp(m,k) = mu(k) * (wind_mid(m,k)-wind_mid(m,km1))/dp(km1);
      mddudp(m,k) = md(k) * (wind_mid(m,k)-wind_mid(m,km1))/dp(km1);
      pgu(m,k) = -momcu * mududp(m,k);
      pgd(m,k) = -momcd * mddudp(m,k);
    }

    // Calculate in-cloud velocity
    // levels adjacent to top and bottom
    {
      const Int k = 1;
      const Int km1 = 0;
      const Int kk = pver-1;
      const Real mupdudp = mu(kk) + du(kk)*dp(kk);
      if (mupdudp > ZMC::mbsth) {
        wind_int_u(m,kk) = (eu(kk)*wind_mid(m,kk)*dp(kk) + pgu(m,kk)*dp(kk)) / mupdudp;
      }
      if (md(k) < -ZMC::mbsth) {
        wind_int_d(m,k) = (-ed(km1)*wind_mid(m,km1)*dp(km1) - pgd(m,km1)*dp(km1)) / md(k);
      }
    }

    // Updraft from bottom to top
    for (Int kk = pver-2; kk >= 0; --kk) {
      const Int kkp1 = kk+1;
      const Real mupdudp = mu(kk) + du(kk)*dp(kk);
      if (mupdudp > ZMC::mbsth) {
        wind_int_u(m,kk) = (mu(kkp1)*wind_int_u(m,kkp1) + eu(kk)*wind_mid(m,kk)*dp(kk) + pgu(m,kk)*dp(kk)) / mupdudp;
      }
    }

    // Downdraft from top to bottom
    for (Int k = 2; k < pver; ++k) {
      const Int km1 = k-1;
      if (md(k) < -ZMC::mbsth) {
        wind_int_d(m,k) = (md(km1)*wind_int_d(m,km1) - ed(km1)*wind_mid(m,km1)*dp(km1) - pgd(m,km1)*dp(km1)) / md(k);
      }
    }

    // Calculate momentum tendency
    for (Int k = ktm; k < pver; ++k) {
      const Int kp1 = k+1;
      wind_tend_tmp(m,k) = (mu(kp1)*(wind_int_u(m,kp1)-wind_int(m,kp1)) -
                            mu(k)  *(wind_int_u(m,k)  -wind_int(m,k)  ) +
                            md(kp1)*(wind_int_d(m,kp1)-wind_int(m,kp1)) -
                            md(k)  *(wind_int_d(m,k)  -wind_int(m,k)  )) / dp(k);
    }

    // dcont for bottom layer
    for (Int k = kbm; k < pver; ++k) {
      if (k == mx) {
        wind_tend_tmp(m,k) = (-mu(k)*(wind_int_u(m,k)-wind_int(m,k)) -
                              md(k)*(wind_int_d(m,k)-wind_int(m,k))) / dp(k);
      }
    }

    // Set outputs
    for (Int k = 0; k < pver; ++k) {
      wind_tend(k,m) = wind_tend_tmp(m,k);
      pguall(k,m) = -pgu(m,k);
      pgdall(k,m) = -pgd(m,k);
      icwu(k,m) = wind_int_u(m,k);
      icwd(k,m) = wind_int_d(m,k);
    }

    // Calculate momentum flux
    for (Int k = ktm; k < pver; ++k) {
      mflux(m,k) = -mu(k)*(wind_int_u(m,k) - wind_int(m,k)) -
                   md(k)*(wind_int_d(m,k) - wind_int(m,k));
    }

    // Calculate winds at the end of the time step
    for (Int k = ktm; k < pver; ++k) {
      const Int kp1 = k+1;
      windf(m,k) = wind_mid(m,k) - (mflux(m,kp1) - mflux(m,k)) * dt/dp(k);
    }

  }); // end parallel_for over nwind

  team.team_barrier();

  // Energy fix to account for dissipation of kinetic energy
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    if (k >= ktm) {
      const Int km1 = (k == 0) ? 0 : k-1;
      const Int kp1 = (k == pver-1) ? pver-1 : k+1;

      // Calculate KE fluxes at top and bot of layer
      const Real utop = (wind0(0,k) + wind0(0,km1)) / 2.0;
      const Real vtop = (wind0(1,k) + wind0(1,km1)) / 2.0;
      const Real ubot = (wind0(0,kp1) + wind0(0,k)) / 2.0;
      const Real vbot = (wind0(1,kp1) + wind0(1,k)) / 2.0;
      const Real fket = utop*mflux(0,k) + vtop*mflux(1,k);
      const Real fkeb = ubot*mflux(0,k+1) + vbot*mflux(1,k+1);

      // Divergence of fluxes gives conservative redistribution of KE
      const Real ketend_cons = (fket-fkeb)/dp(k);

      // Tendency in kinetic energy from momentum transport
      const Real ketend = ((windf(0,k)*windf(0,k) + windf(1,k)*windf(1,k)) -
                          (wind0(0,k)*wind0(0,k) + wind0(1,k)*wind0(1,k))) * 0.5 / dt;

      // The difference is the dissipation
      seten(k) = ketend_cons - ketend;
    }
  });

  team.team_barrier();

  workspace.template release_many_contiguous<12>(
    {&wind0_1d, &windf_1d, &wind_mid_1d, &wind_int_1d, &wind_int_d_1d, &wind_int_u_1d,
     &wind_tend_tmp_1d, &mududp_1d, &mddudp_1d, &pgu_1d, &pgd_1d, &mflux_1d});
}

} // namespace zm
} // namespace scream

#endif
