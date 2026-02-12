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

  // Get workspace manager views (1d views)
  const auto ws_size = pver * nwind;
  uview_1d<Real> wind0_ws        = workspace.take_macro_block("wind0",        ws_size);
  uview_1d<Real> windf_ws        = workspace.take_macro_block("windf",        ws_size);
  uview_1d<Real> wind_mid_ws     = workspace.take_macro_block("wind_mid",     pver*nwind);
  uview_1d<Real> wind_int_ws     = workspace.take_macro_block("wind_int",     pver*nwind);
  uview_1d<Real> wind_int_d_ws   = workspace.take_macro_block("wind_int_d",   pver*nwind);
  uview_1d<Real> wind_int_u_ws   = workspace.take_macro_block("wind_int_u",   pver*nwind);
  uview_1d<Real> wind_tend_tmp_ws= workspace.take_macro_block("wind_tend_tmp",pver*nwind);
  uview_1d<Real> mududp_ws       = workspace.take_macro_block("mududp",       pver*nwind);
  uview_1d<Real> mddudp_ws       = workspace.take_macro_block("mddudp",       pver*nwind);
  uview_1d<Real> pgu_ws          = workspace.take_macro_block("pgu",          pver*nwind);
  uview_1d<Real> pgd_ws          = workspace.take_macro_block("pgd",          pver*nwind);
  uview_1d<Real> mflux_ws        = workspace.take_macro_block("mflux",        pverp*nwind);

  // Convert 1d workspace views to 2d views (pver, nwind)
  auto wind0 = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(wind0_ws.data(), pver, nwind);
  auto windf = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(windf_ws.data(), pver, nwind);
  auto wind_mid = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(wind_mid_ws.data(), pver, nwind);
  auto wind_int = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(wind_int_ws.data(), pver, nwind);
  auto wind_int_d = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(wind_int_d_ws.data(), pver, nwind);
  auto wind_int_u = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(wind_int_u_ws.data(), pver, nwind);
  auto wind_tend_tmp = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(wind_tend_tmp_ws.data(), pver, nwind);
  auto mududp = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(mududp_ws.data(), pver, nwind);
  auto mddudp = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(mddudp_ws.data(), pver, nwind);
  auto pgu = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(pgu_ws.data(), pver, nwind);
  auto pgd = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(pgd_ws.data(), pver, nwind);
  auto mflux = ekat::screamify<Real[::][::], ekat::MemLayout::Right>(mflux_ws.data(), pverp, nwind);

  // Initialize outputs
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver*nwind), [&] (const Int& idx) {
    const Int k = idx / nwind;
    const Int m = idx % nwind;
    wind_tend(k,m) = 0.0;
    pguall(k,m) = 0.0;
    pgdall(k,m) = 0.0;
    icwu(k,m) = wind_in(k,m);
    icwd(k,m) = wind_in(k,m);
    wind0(k,m) = 0.0;
    windf(k,m) = 0.0;
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    seten(k) = 0.0;
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pverp*nwind), [&] (const Int& idx) {
    mflux_ws(idx) = 0.0;
  });

  team.team_barrier();

  // Loop over each wind component using team parallelism
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nwind), [&] (const Int& m) {

    // Gather up the winds
    for (Int k = 0; k < pver; ++k) {
      wind_mid(k,m) = wind_in(k,m);
      wind0(k,m) = wind_mid(k,m);
    }

    // Interpolate winds to interfaces
    for (Int k = 0; k < pver; ++k) {
      const Int km1 = (k == 0) ? 0 : k-1;
      wind_int(k,m) = 0.5 * (wind_mid(k,m) + wind_mid(km1,m));
      wind_int_u(k,m) = wind_int(k,m);
      wind_int_d(k,m) = wind_int(k,m);
      wind_tend_tmp(k,m) = 0.0;
    }

    // Calculate pressure perturbation terms
    // upper boundary
    pgu(0,m) = 0.0;
    pgd(0,m) = 0.0;

    // interior points
    for (Int k = 1; k < pver-1; ++k) {
      const Int km1 = k-1;
      const Int kp1 = k+1;
      mududp(k,m) = (mu(k)  *(wind_mid(k,m)  -wind_mid(km1,m))/dp(km1) +
                     mu(kp1)*(wind_mid(kp1,m)-wind_mid(k,m)  )/dp(k));
      mddudp(k,m) = (md(k)  *(wind_mid(k,m)  -wind_mid(km1,m))/dp(km1) +
                     md(kp1)*(wind_mid(kp1,m)-wind_mid(k,m)  )/dp(k));
      pgu(k,m) = -momcu * 0.5 * mududp(k,m);
      pgd(k,m) = -momcd * 0.5 * mddudp(k,m);
    }

    // bottom boundary
    {
      const Int k = pver-1;
      const Int km1 = k-1;
      mududp(k,m) = mu(k) * (wind_mid(k,m)-wind_mid(km1,m))/dp(km1);
      mddudp(k,m) = md(k) * (wind_mid(k,m)-wind_mid(km1,m))/dp(km1);
      pgu(k,m) = -momcu * mududp(k,m);
      pgd(k,m) = -momcd * mddudp(k,m);
    }

    // Calculate in-cloud velocity
    // levels adjacent to top and bottom
    {
      const Int k = 1;
      const Int km1 = 0;
      const Int kk = pver-1;
      const Real mupdudp = mu(kk) + du(kk)*dp(kk);
      if (mupdudp > mbsth) {
        wind_int_u(kk,m) = (eu(kk)*wind_mid(kk,m)*dp(kk) + pgu(kk,m)*dp(kk)) / mupdudp;
      }
      if (md(k) < -mbsth) {
        wind_int_d(k,m) = (-ed(km1)*wind_mid(km1,m)*dp(km1) - pgd(km1,m)*dp(km1)) / md(k);
      }
    }

    // Updraft from bottom to top
    for (Int kk = pver-2; kk >= 0; --kk) {
      const Int kkp1 = kk+1;
      const Real mupdudp = mu(kk) + du(kk)*dp(kk);
      if (mupdudp > mbsth) {
        wind_int_u(kk,m) = (mu(kkp1)*wind_int_u(kkp1,m) + eu(kk)*wind_mid(kk,m)*dp(kk) + pgu(kk,m)*dp(kk)) / mupdudp;
      }
    }

    // Downdraft from top to bottom
    for (Int k = 2; k < pver; ++k) {
      const Int km1 = k-1;
      if (md(k) < -mbsth) {
        wind_int_d(k,m) = (md(km1)*wind_int_d(km1,m) - ed(km1)*wind_mid(km1,m)*dp(km1) - pgd(km1,m)*dp(km1)) / md(k);
      }
    }

    // Calculate momentum tendency
    for (Int k = jt; k < pver; ++k) {
      const Int kp1 = k+1;
      wind_tend_tmp(k,m) = (mu(kp1)*(wind_int_u(kp1,m)-wind_int(kp1,m)) -
                            mu(k)  *(wind_int_u(k,m)  -wind_int(k,m)  ) +
                            md(kp1)*(wind_int_d(kp1,m)-wind_int(kp1,m)) -
                            md(k)  *(wind_int_d(k,m)  -wind_int(k,m)  )) / dp(k);
    }

    // dcont for bottom layer
    if (mx >= mx && mx < pver) {
      const Int k = mx;
      wind_tend_tmp(k,m) = (-mu(k)*(wind_int_u(k,m)-wind_int(k,m)) -
                            md(k)*(wind_int_d(k,m)-wind_int(k,m))) / dp(k);
    }

    // Set outputs
    for (Int k = 0; k < pver; ++k) {
      wind_tend(k,m) = wind_tend_tmp(k,m);
      pguall(k,m) = -pgu(k,m);
      pgdall(k,m) = -pgd(k,m);
      icwu(k,m) = wind_int_u(k,m);
      icwd(k,m) = wind_int_d(k,m);
    }

    // Calculate momentum flux
    for (Int k = jt; k < pver; ++k) {
      mflux(k,m) = -mu(k)*(wind_int_u(k,m) - wind_int(k,m)) -
                   md(k)*(wind_int_d(k,m) - wind_int(k,m));
    }

    // Calculate winds at the end of the time step
    for (Int k = jt; k < pver; ++k) {
      const Int kp1 = k+1;
      windf(k,m) = wind_mid(k,m) - (mflux(kp1,m) - mflux(k,m)) * dt/dp(k);
    }

  }); // end parallel_for over nwind

  team.team_barrier();

  // Energy fix to account for dissipation of kinetic energy
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    if (k >= jt) {
      const Int km1 = (k == 0) ? 0 : k-1;
      const Int kp1 = (k == pver-1) ? pver-1 : k+1;
      
      // Calculate KE fluxes at top and bot of layer
      const Real utop = (wind0(k,0) + wind0(km1,0)) / 2.0;
      const Real vtop = (wind0(k,1) + wind0(km1,1)) / 2.0;
      const Real ubot = (wind0(kp1,0) + wind0(k,0)) / 2.0;
      const Real vbot = (wind0(kp1,1) + wind0(k,1)) / 2.0;
      const Real fket = utop*mflux(k,0) + vtop*mflux(k,1);
      const Real fkeb = ubot*mflux(k+1,0) + vbot*mflux(k+1,1);
      
      // Divergence of fluxes gives conservative redistribution of KE
      const Real ketend_cons = (fket-fkeb)/dp(k);
      
      // Tendency in kinetic energy from momentum transport
      const Real ketend = ((windf(k,0)*windf(k,0) + windf(k,1)*windf(k,1)) -
                          (wind0(k,0)*wind0(k,0) + wind0(k,1)*wind0(k,1))) * 0.5 / dt;
      
      // The difference is the dissipation
      seten(k) = ketend_cons - ketend;
    }
  });

  team.team_barrier();

  // Release workspace
  workspace.release_macro_block(wind0_ws,        ws_size);
  workspace.release_macro_block(windf_ws,        ws_size);
  workspace.release_macro_block(wind_mid_ws,     pver*nwind);
  workspace.release_macro_block(wind_int_ws,     pver*nwind);
  workspace.release_macro_block(wind_int_d_ws,   pver*nwind);
  workspace.release_macro_block(wind_int_u_ws,   pver*nwind);
  workspace.release_macro_block(wind_tend_tmp_ws,pver*nwind);
  workspace.release_macro_block(mududp_ws,       pver*nwind);
  workspace.release_macro_block(mddudp_ws,       pver*nwind);
  workspace.release_macro_block(pgu_ws,          pver*nwind);
  workspace.release_macro_block(pgd_ws,          pver*nwind);
  workspace.release_macro_block(mflux_ws,        pverp*nwind);
}

} // namespace zm
} // namespace scream

#endif
