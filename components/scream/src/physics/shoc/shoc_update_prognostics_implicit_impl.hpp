#ifndef SHOC_UPDATE_PROGNOSTICS_IMPLICIT_IMPL_HPP
#define SHOC_UPDATE_PROGNOSTICS_IMPLICIT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc update_prognostics_implicit. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::update_prognostics_implicit(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Int&                   num_tracer,
  const Scalar&                dtime,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& rho_zt,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<const Spack>& tk,
  const uview_1d<const Spack>& tkh,
  const Scalar&                uw_sfc,
  const Scalar&                vw_sfc,
  const Scalar&                wthl_sfc,
  const Scalar&                wqw_sfc,
  const uview_1d<const Spack>& wtracer_sfc,
  const uview_1d<Spack>&       rdp_zt,
  const uview_1d<Spack>&       tmpi,
  const uview_1d<Spack>&       tkh_zi,
  const uview_1d<Spack>&       tk_zi,
  const uview_1d<Spack>&       rho_zi,
  const uview_1d<Scalar>&      du,
  const uview_1d<Scalar>&      dl,
  const uview_1d<Scalar>&      d,
  const uview_2d<Spack>&       X1,
  const uview_1d<Spack>&       thetal,
  const uview_1d<Spack>&       qw,
  const uview_2d<Spack>&       tracer,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       u_wind,
  const uview_1d<Spack>&       v_wind)
{
  const auto last_nlev_pack = (nlev-1)/Spack::n;
  const auto last_nlev_indx = (nlev-1)%Spack::n;
  const auto last_nlevi_pack = (nlevi-1)/Spack::n;
  const auto last_nlevi_indx = (nlevi-1)%Spack::n;

  // linearly interpolate tkh, tk, and air density onto the interface grids
  linear_interp(team,zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,0);

  // Define the tmpi variable, which is really dt*(g*rho)**2/dp
  // at interfaces. Substitue dp = g*rho*dz in the above equation
  compute_tmpi(team, nlevi, dtime, rho_zi, dz_zi, tmpi);

  // compute 1/dp term, needed in diffusion solver
  dp_inverse(team, nlev, rho_zt, dz_zt, rdp_zt);
  team.team_barrier();

  // compute terms needed for the implicit surface stress (ksrf)
  // and tke flux calc (wtke_sfc)
  Scalar ksrf, wtke_sfc;
  {
    const Scalar wsmin = 1;
    const Scalar ksrfmin = 1e-4;
    const Scalar ustarmin = 0.01;

    const Scalar rho = rho_zi(last_nlevi_pack)[last_nlevi_indx];
    const Scalar uw = uw_sfc;
    const Scalar vw = vw_sfc;

    const Scalar taux = rho*uw;
    const Scalar tauy = rho*vw;

    const Scalar u_wind_sfc = u_wind(last_nlev_pack)[last_nlev_indx];
    const Scalar v_wind_sfc = v_wind(last_nlev_pack)[last_nlev_indx];

    const Scalar ws = ekat::impl::max(std::sqrt((u_wind_sfc*u_wind_sfc) + v_wind_sfc*v_wind_sfc), wsmin);
    const Scalar tau = std::sqrt(taux*taux + tauy*tauy);
    ksrf = ekat::impl::max(tau/ws, ksrfmin);

    const Scalar ustar = ekat::impl::max(std::sqrt(std::sqrt(uw*uw + vw*vw)), ustarmin);
    wtke_sfc = ustar*ustar*ustar;
  }

  // compute surface fluxes for liq. potential temp, water and tke
  {
    const auto cmnfac = dtime*(C::gravit*rho_zi(last_nlevi_pack)[last_nlevi_indx]*rdp_zt(last_nlev_pack)[last_nlev_indx]);
    Kokkos::single(Kokkos::PerTeam(team), [&] () {
      thetal(last_nlev_pack)[last_nlev_indx] += cmnfac*wthl_sfc;
      qw(last_nlev_pack)[last_nlev_indx] += cmnfac*wqw_sfc;
      tke(last_nlev_pack)[last_nlev_indx] += cmnfac*wtke_sfc;
    });

    const auto tracer_sfc = Kokkos::subview(tracer, nlev-1, Kokkos::ALL());
    const Int num_tracer_pack = ekat::npack<Spack>(num_tracer);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_tracer_pack), [&] (const Int& p) {
      tracer_sfc(p) += cmnfac*wtracer_sfc(p);
    });
  }

  // Store RHS values in X1 and tracer for 1st and 2nd solve respectively
  team.team_barrier();
  const auto s_u_wind = ekat::scalarize(u_wind);
  const auto s_v_wind = ekat::scalarize(v_wind);
  const auto s_thetal = ekat::scalarize(thetal);
  const auto s_qw = ekat::scalarize(qw);
  const auto s_tke = ekat::scalarize(tke);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
    X1(k,0)[0] = s_u_wind(k);
    X1(k,1/Spack::n)[1%Spack::n] = s_v_wind(k);

    tracer(k,(num_tracer)/Spack::n)[(num_tracer)%Spack::n] = s_thetal(k);
    tracer(k,(num_tracer+1)/Spack::n)[(num_tracer+1)%Spack::n] = s_qw(k);
    tracer(k,(num_tracer+2)/Spack::n)[(num_tracer+2)%Spack::n] = s_tke(k);
  });

  // march u_wind and v_wind one step forward using implicit solver
  {
    // Call decomp for momentum variables
    vd_shoc_decomp(team, nlev, tk_zi, tmpi, rdp_zt, dtime, ksrf, du, dl, d);

    // Solve
    team.team_barrier();
    vd_shoc_solve(team, du, dl, d, X1);
  }

  // march temperature, total water, tke,and tracers one step forward using implicit solver
  {
    // Call decomp for thermo variables. Fluxes applied explicitly, so zero
    // fluxes out for implicit solver decomposition.
    team.team_barrier();
    vd_shoc_decomp(team, nlev, tkh_zi, tmpi, rdp_zt, dtime, 0, du, dl, d);

    // Solve
    team.team_barrier();
    vd_shoc_solve(team, du, dl, d, tracer);
  }

  // Copy RHS values from X1 and X2
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
    s_u_wind(k) = X1(k,0)[0];
    s_v_wind(k) = X1(k,1/Spack::n)[1%Spack::n];

    s_thetal(k) = tracer(k,(num_tracer)/Spack::n)[(num_tracer)%Spack::n];
    s_qw(k) = tracer(k,(num_tracer+1)/Spack::n)[(num_tracer+1)%Spack::n];
    s_tke(k) = tracer(k,(num_tracer+2)/Spack::n)[(num_tracer+2)%Spack::n];
  });
}

} // namespace shoc
} // namespace scream

#endif
