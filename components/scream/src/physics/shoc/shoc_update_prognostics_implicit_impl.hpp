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
  const Int&                   num_qtracers,
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
  const Workspace&             workspace,
  const uview_1d<Spack>&       thetal,
  const uview_1d<Spack>&       qw,
  const uview_2d<Spack>&       qtracers,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       u_wind,
  const uview_1d<Spack>&       v_wind)
{
  const int num_wind_transpose_packs = ekat::npack<Spack>(2);
  const int num_qtracers_transpose_packs = ekat::npack<Spack>(num_qtracers+3);

  const auto nlev_v  = (nlev-1)/Spack::n;
  const auto nlev_p  = (nlev-1)%Spack::n;
  const auto nlevi_v = (nlevi-1)/Spack::n;
  const auto nlevi_p = (nlevi-1)%Spack::n;

  // Define temporary variables via the WorkspaceManager

  // 1d allocations
  uview_1d<Spack> tmpi, tkh_zi,
                  tk_zi, rho_zi,
                  rdp_zt;
  uview_1d<Scalar> du_workspace, dl_workspace, d_workspace;

  workspace.template take_many_contiguous_unsafe<5>(
    {"tmpi", "tkh_zi", "tk_zi", "rho_zi", "rdp_zt"},
    {&tmpi, &tkh_zi, &tk_zi, &rho_zi, &rdp_zt});

  workspace.template take_many_contiguous_unsafe<3, Scalar>(
    {"du_workspace", "dl_workspace", "d_workspace"},
    {&du_workspace, &dl_workspace, &d_workspace});
  auto du = Kokkos::subview(du_workspace, Kokkos::make_pair(0,nlev));
  auto dl = Kokkos::subview(dl_workspace, Kokkos::make_pair(0,nlev));
  auto d  = Kokkos::subview(d_workspace,  Kokkos::make_pair(0,nlev));

  // 2d allocations for solver RHS
  const int n_wind_slots = num_wind_transpose_packs*Spack::n;
  const int n_trac_slots = num_qtracers_transpose_packs*Spack::n;
  const auto wind_slot    = workspace.template take_macro_block<Scalar>("wind_slot",n_wind_slots);
  const auto tracers_slot = workspace.template take_macro_block<Scalar>("tracers_slot",n_trac_slots);

  // Reshape 2d views
  const auto wind_transpose     = uview_2d<Spack>(reinterpret_cast<Spack*>(wind_slot.data()),
                                                  nlev,ekat::npack<Spack>(2));
  const auto qtracers_transpose = uview_2d<Spack>(reinterpret_cast<Spack*>(tracers_slot.data()),
                                                  nlev,ekat::npack<Spack>(num_qtracers+3));

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

    const Scalar rho = rho_zi(nlevi_v)[nlevi_p];
    const Scalar uw = uw_sfc;
    const Scalar vw = vw_sfc;

    const Scalar taux = rho*uw;
    const Scalar tauy = rho*vw;

    const Scalar u_wind_sfc = u_wind(nlev_v)[nlev_p];
    const Scalar v_wind_sfc = v_wind(nlev_v)[nlev_p];

    const Scalar ws = ekat::impl::max(std::sqrt((u_wind_sfc*u_wind_sfc) + v_wind_sfc*v_wind_sfc), wsmin);
    const Scalar tau = std::sqrt(taux*taux + tauy*tauy);
    ksrf = ekat::impl::max(tau/ws, ksrfmin);

    const Scalar ustar = ekat::impl::max(std::sqrt(std::sqrt(uw*uw + vw*vw)), ustarmin);
    wtke_sfc = ustar*ustar*ustar;
  }

  // compute surface fluxes for liq. potential temp, water and tke
  {
    const auto cmnfac = dtime*(C::gravit*rho_zi(nlevi_v)[nlevi_p]*rdp_zt(nlev_v)[nlev_p]);
    Kokkos::single(Kokkos::PerTeam(team), [&] () {
      thetal(nlev_v)[nlev_p] += cmnfac*wthl_sfc;
      qw(nlev_v)[nlev_p] += cmnfac*wqw_sfc;
      tke(nlev_v)[nlev_p] += cmnfac*wtke_sfc;
    });

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_qtracers), [&] (const Int& q) {
      qtracers(q,nlev_v)[nlev_p] += cmnfac*wtracer_sfc(q/Spack::n)[q%Spack::n];
    });
  }

  // Store RHS values in X1 and tracer for 1st and 2nd solve respectively
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
    const int k_v = k/Spack::n;
    const int k_p = k%Spack::n;

    wind_transpose(k,0/Spack::n)[0%Spack::n] = u_wind(k_v)[k_p];
    wind_transpose(k,1/Spack::n)[1%Spack::n] = v_wind(k_v)[k_p];

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_qtracers), [&] (const Int& q) {
      qtracers_transpose(k, q/Spack::n)[q%Spack::n] = qtracers(q,k_v)[k_p];
    });
    qtracers_transpose(k, num_qtracers/Spack::n)[num_qtracers%Spack::n]         = thetal(k_v)[k_p];
    qtracers_transpose(k, (num_qtracers+1)/Spack::n)[(num_qtracers+1)%Spack::n] = qw(k_v)[k_p];
    qtracers_transpose(k, (num_qtracers+2)/Spack::n)[(num_qtracers+2)%Spack::n] = tke(k_v)[k_p];
  });

  // march u_wind and v_wind one step forward using implicit solver
  {
    // Call decomp for momentum variables
    vd_shoc_decomp(team, nlev, tk_zi, tmpi, rdp_zt, dtime, ksrf, du, dl, d);

    // Solve
    team.team_barrier();
    vd_shoc_solve(team, du, dl, d, wind_transpose);
  }

  // march temperature, total water, tke,and tracers one step forward using implicit solver
  {
    // Call decomp for thermo variables. Fluxes applied explicitly, so zero
    // fluxes out for implicit solver decomposition.
    team.team_barrier();
    vd_shoc_decomp(team, nlev, tkh_zi, tmpi, rdp_zt, dtime, 0, du, dl, d);

    // Solve
    team.team_barrier();
    vd_shoc_solve(team, du, dl, d, qtracers_transpose);
  }

  // Copy RHS values from X1 and X2
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
    const int k_v = k/Spack::n;
    const int k_p = k%Spack::n;

    u_wind(k_v)[k_p] = wind_transpose(k,0/Spack::n)[0%Spack::n];
    v_wind(k_v)[k_p] = wind_transpose(k,1/Spack::n)[1%Spack::n];

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_qtracers), [&] (const Int& q) {
      qtracers(q,k_v)[k_p] = qtracers_transpose(k, q/Spack::n)[q%Spack::n];
    });
    thetal(k_v)[k_p] = qtracers_transpose(k, num_qtracers/Spack::n)[num_qtracers%Spack::n];
    qw(k_v)[k_p]     = qtracers_transpose(k, (num_qtracers+1)/Spack::n)[(num_qtracers+1)%Spack::n];
    tke(k_v)[k_p]    = qtracers_transpose(k, (num_qtracers+2)/Spack::n)[(num_qtracers+2)%Spack::n];
  });


  // Release temporary variables from the workspace
  team.team_barrier();
  workspace.template release_macro_block<Scalar>(tracers_slot,n_trac_slots);
  workspace.template release_macro_block<Scalar>(wind_slot,n_wind_slots);
  workspace.template release_many_contiguous<3,Scalar>(
    {&du_workspace, &dl_workspace, &d_workspace});
  workspace.template release_many_contiguous<5>(
    {&tmpi, &tkh_zi, &tk_zi, &rho_zi, &rdp_zt});
}

} // namespace shoc
} // namespace scream

#endif
