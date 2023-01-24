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
  const int num_wind_transpose_packs = ekat::npack<Spack>(2);
  const int num_qtracers_transpose_packs = ekat::npack<Spack>(num_qtracers+3);

  const int n_wind_slots = num_wind_transpose_packs*Spack::n;
  const int n_trac_slots = num_qtracers_transpose_packs*Spack::n;

  const auto wind_slot    = workspace.template take_macro_block<Scalar>("wind_slot",n_wind_slots);
  const auto tracers_slot = workspace.template take_macro_block<Scalar>("tracers_slot",n_trac_slots);

  // Reshape 2d views
  const auto wind_rhs     = uview_2d<Spack>(reinterpret_cast<Spack*>(wind_slot.data()),
                                            nlev, num_wind_transpose_packs);
  const auto qtracers_rhs  = uview_2d<Spack>(reinterpret_cast<Spack*>(tracers_slot.data()),
                                            nlev, num_qtracers_transpose_packs);

  // scalarized versions of some views will be needed
  const auto rdp_zt_s       = ekat::scalarize(rdp_zt);
  const auto rho_zi_s       = ekat::scalarize(rho_zi);
  const auto u_wind_s       = ekat::scalarize(u_wind);
  const auto v_wind_s       = ekat::scalarize(v_wind);
  const auto wind_rhs_s     = ekat::scalarize(wind_rhs);
  const auto thetal_s       = ekat::scalarize(thetal);
  const auto qw_s           = ekat::scalarize(qw);
  const auto tke_s          = ekat::scalarize(tke);
  const auto qtracers_s     = ekat::scalarize(qtracers);
  const auto qtracers_rhs_s = ekat::scalarize(qtracers_rhs);
  const auto wtracer_sfc_s  = ekat::scalarize(wtracer_sfc);

  // linearly interpolate tkh, tk, and air density onto the interface grids
  linear_interp(team,zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,0);
  // A note on team_barrier after linear_interp: One does not need a barrier
  // between two parallel_for loops with the same index set: e.g., two
  // [0,nlev_pack) parallel_for ranges. But one does if e.g. a linear_interp is
  // over nlevi_pack and the next parallel_for is over nlev_pack.
  team.team_barrier();

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

    const Scalar rho = rho_zi_s(nlevi-1);
    const Scalar uw = uw_sfc;
    const Scalar vw = vw_sfc;

    const Scalar taux = rho*uw;
    const Scalar tauy = rho*vw;

    const Scalar u_wind_sfc = u_wind_s(nlev-1);
    const Scalar v_wind_sfc = v_wind_s(nlev-1);

    const Scalar ws = ekat::impl::max(std::sqrt((u_wind_sfc*u_wind_sfc) + v_wind_sfc*v_wind_sfc), wsmin);
    const Scalar tau = std::sqrt(taux*taux + tauy*tauy);
    ksrf = ekat::impl::max(tau/ws, ksrfmin);

    const Scalar ustar = ekat::impl::max(std::sqrt(std::sqrt(uw*uw + vw*vw)), ustarmin);
    wtke_sfc = ustar*ustar*ustar;
  }

  // compute surface fluxes for liq. potential temp, water and tke
  {
    const auto cmnfac = dtime*(C::gravit*rho_zi_s(nlevi-1)*rdp_zt_s(nlev-1));
    Kokkos::single(Kokkos::PerTeam(team), [&] () {
      thetal_s(nlev-1) += cmnfac*wthl_sfc;
      qw_s(nlev-1)     += cmnfac*wqw_sfc;
      tke_s(nlev-1)    += cmnfac*wtke_sfc;
    });

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, num_qtracers), [&] (const Int& q) {
      qtracers_s(q, nlev-1) += cmnfac*wtracer_sfc_s(q);
    });
  }

  // Store RHS values in wind_rhs and qtracers_rhs for 1st and 2nd solve respectively
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
    wind_rhs_s(k,0) = u_wind_s(k);
    wind_rhs_s(k,1) = v_wind_s(k);

    // The rhs version of the tracers is the transpose of the input/output layout
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_qtracers), [&] (const Int& q) {
      qtracers_rhs_s(k, q) = qtracers_s(q, k);
    });
    qtracers_rhs_s(k, num_qtracers)   = thetal_s(k);
    qtracers_rhs_s(k, num_qtracers+1) = qw_s(k);
    qtracers_rhs_s(k, num_qtracers+2) = tke_s(k);
  });

  // march u_wind and v_wind one step forward using implicit solver
  {
    // Call decomp for momentum variables
    vd_shoc_decomp(team, nlev, tk_zi, tmpi, rdp_zt, dtime, ksrf, du, dl, d);

    // Solve
    team.team_barrier();
    vd_shoc_solve(team, du, dl, d, wind_rhs);
  }

  // march temperature, total water, tke,and tracers one step forward using implicit solver
  {
    // Call decomp for thermo variables. Fluxes applied explicitly, so zero
    // fluxes out for implicit solver decomposition.
    team.team_barrier();
    vd_shoc_decomp(team, nlev, tkh_zi, tmpi, rdp_zt, dtime, 0, du, dl, d);

    // Solve
    team.team_barrier();
    vd_shoc_solve(team, du, dl, d, qtracers_rhs);
  }

  // Copy RHS values back into output variables
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev), [&] (const Int& k) {
    u_wind_s(k) = wind_rhs_s(k, 0);
    v_wind_s(k) = wind_rhs_s(k, 1);

    // Transpose tracers back to  input/output layout
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, num_qtracers), [&] (const Int& q) {
      qtracers_s(q, k) = qtracers_rhs_s(k, q);
    });
    thetal_s(k) = qtracers_rhs_s(k, num_qtracers);
    qw_s(k)     = qtracers_rhs_s(k, num_qtracers+1);
    tke_s(k)    = qtracers_rhs_s(k, num_qtracers+2);
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
