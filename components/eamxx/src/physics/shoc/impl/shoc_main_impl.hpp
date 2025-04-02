#ifndef SHOC_MAIN_IMPL_HPP
#define SHOC_MAIN_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include "ekat/kokkos/ekat_subview_utils.hpp"

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_main. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
Int Functions<S,D>::shoc_init(
  const Int&                  nbot_shoc,
  const Int&                  ntop_shoc,
  const view_1d<const Spack>& pref_mid)
{
  // This function calculates the maximum number of levels
  // in pbl from surface

  using ExeSpace = typename KT::ExeSpace;
  view_1d<Int> npbl_d("npbl",1);

  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Scalar pblmaxp = SC::pblmaxp;

    Int npbl_val = 1;

    const int begin_pack_indx = ntop_shoc/Spack::n;
    const int end_pack_indx   = nbot_shoc/Spack::n+1;
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, begin_pack_indx, end_pack_indx),
                                                    [&] (const Int& k, Int& local_max) {
      auto range = ekat::range<IntSmallPack>(k*Spack::n);
      auto condition = (range >= ntop_shoc && range < nbot_shoc);
      if (condition.any()) {
        condition = condition && pref_mid(k) >= pblmaxp;
      }

      auto levels_from_surface = nbot_shoc - range;
      levels_from_surface.set(!condition, 1);

      if (local_max < ekat::max(levels_from_surface))
        local_max = ekat::max(levels_from_surface);

    }, Kokkos::Max<Int>(npbl_val));

    npbl_d(0) = npbl_val;
  });

  const auto host_view = Kokkos::create_mirror_view(npbl_d);
  Kokkos::deep_copy(host_view, npbl_d);

  return host_view(0);
}

#ifndef SCREAM_SHOC_SMALL_KERNELS
template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_main_internal(
  const MemberType&            team,
  const Int&                   nlev,         // Number of levels
  const Int&                   nlevi,        // Number of levels on interface grid
  const Int&                   npbl,         // Maximum number of levels in pbl from surface
  const Int&                   nadv,         // Number of times to loop SHOC
  const Int&                   num_qtracers, // Number of tracers
  const Scalar&                dtime,        // SHOC timestep [s]
  // Runtime Parameters
  const Scalar&                lambda_low,
  const Scalar&                lambda_high,
  const Scalar&                lambda_slope,
  const Scalar&                lambda_thresh,
  const Scalar&                thl2tune,
  const Scalar&                qw2tune,
  const Scalar&                qwthl2tune,
  const Scalar&                w2tune,
  const Scalar&                length_fac,
  const Scalar&                c_diag_3rd_mom,
  const Scalar&                Ckh,
  const Scalar&                Ckm,
  // Input Variables
  const Scalar&                dx,
  const Scalar&                dy,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& presi,
  const uview_1d<const Spack>& pdel,
  const uview_1d<const Spack>& thv,
  const uview_1d<const Spack>& w_field,
  const Scalar&                wthl_sfc,
  const Scalar&                wqw_sfc,
  const Scalar&                uw_sfc,
  const Scalar&                vw_sfc,
  const uview_1d<const Spack>& wtracer_sfc,
  const uview_1d<const Spack>& inv_exner,
  const Scalar&                phis,
  // Workspace/Local Variables
  const Workspace&             workspace,
  // Input/Output Variables
  const uview_1d<Spack>&       host_dse,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       thetal,
  const uview_1d<Spack>&       qw,
  const uview_1d<Spack>&       u_wind,
  const uview_1d<Spack>&       v_wind,
  const uview_1d<Spack>&       wthv_sec,
  const uview_2d_strided<Spack>& qtracers,
  const uview_1d<Spack>&       tk,
  const uview_1d<Spack>&       shoc_cldfrac,
  const uview_1d<Spack>&       shoc_ql,
  // Output Variables
  Scalar&                      pblh,
  Scalar&                      ustar,
  Scalar&                      obklen,
  const uview_1d<Spack>&       shoc_ql2,
  const uview_1d<Spack>&       tkh,
  // Diagnostic Output Variables
  const uview_1d<Spack>&       shoc_mix,
  const uview_1d<Spack>&       shoc_cond,
  const uview_1d<Spack>&       shoc_evap,
  const uview_1d<Spack>&       w_sec,
  const uview_1d<Spack>&       thl_sec,
  const uview_1d<Spack>&       qw_sec,
  const uview_1d<Spack>&       qwthl_sec,
  const uview_1d<Spack>&       wthl_sec,
  const uview_1d<Spack>&       wqw_sec,
  const uview_1d<Spack>&       wtke_sec,
  const uview_1d<Spack>&       uw_sec,
  const uview_1d<Spack>&       vw_sec,
  const uview_1d<Spack>&       w3,
  const uview_1d<Spack>&       wqls_sec,
  const uview_1d<Spack>&       brunt,
  const uview_1d<Spack>&       isotropy)
{

  // Define temporary variables
  uview_1d<Spack> rho_zt, shoc_qv, shoc_tabs, dz_zt, dz_zi;
  workspace.template take_many_and_reset<5>(
    {"rho_zt", "shoc_qv", "shoc_tabs", "dz_zt", "dz_zi"},
    {&rho_zt, &shoc_qv, &shoc_tabs, &dz_zt, &dz_zi});

  // Local scalars
  Scalar se_b{0},   ke_b{0},   wv_b{0}, wl_b{0},
         se_a{0},   ke_a{0},   wv_a{0}, wl_a{0},
         kbfs{0},   ustar2{0}, wstar{0};

  // Scalarize some views for single entry access
  const auto s_thetal  = ekat::scalarize(thetal);
  const auto s_shoc_ql = ekat::scalarize(shoc_ql);
  const auto s_shoc_qv = ekat::scalarize(shoc_qv);

  // Compute integrals of static energy, kinetic energy, water vapor, and liquid water
  // for the computation of total energy before SHOC is called.  This is for an
  // effort to conserve energy since liquid water potential temperature (which SHOC
  // conserves) and static energy (which E3SM conserves) are not exactly equal.
  shoc_energy_integrals(team,nlev,host_dse,pdel,qw,shoc_ql,u_wind,v_wind, // Input
                        se_b,ke_b,wv_b,wl_b);                             // Output

  for (Int t=0; t<nadv; ++t) {
    // Check TKE to make sure values lie within acceptable
    // bounds after host model performs horizontal advection
    check_tke(team,nlev, // Input
              tke);      // Input/Output

    // Define vertical grid arrays needed for
    // vertical derivatives in SHOC, also
    // define air density (rho_zt)
    shoc_grid(team,nlev,nlevi,      // Input
              zt_grid,zi_grid,pdel, // Input
              dz_zt,dz_zi,rho_zt);  // Output

    // Compute the planetary boundary layer height, which is an
    // input needed for the length scale calculation.

    // Update SHOC water vapor,
    // to be used by the next two routines
    compute_shoc_vapor(team,nlev,qw,shoc_ql, // Input
                       shoc_qv);             // Output

    // Update SHOC temperature
    compute_shoc_temperature(team,nlev,thetal,  // Input
                             shoc_ql,inv_exner, // Input
                             shoc_tabs);        // Output

    team.team_barrier();
    shoc_diag_obklen(uw_sfc,vw_sfc,     // Input
                     wthl_sfc, wqw_sfc, // Input
                     s_thetal(nlev-1),  // Input
                     s_shoc_ql(nlev-1), // Input
                     s_shoc_qv(nlev-1), // Input
                     ustar,kbfs,obklen); // Output

    pblintd(team,nlev,nlevi,npbl,     // Input
            zt_grid,zi_grid,thetal,   // Input
            shoc_ql,shoc_qv,u_wind,   // Input
            v_wind,ustar,obklen,kbfs, // Input
            shoc_cldfrac,             // Input
            workspace,                // Workspace
            pblh);                    // Output

    // Update the turbulent length scale
    shoc_length(team,nlev,nlevi,       // Input
                length_fac,            // Runtime Options
                dx,dy,                 // Input
                zt_grid,zi_grid,dz_zt, // Input
                tke,thv,               // Input
                workspace,             // Workspace
                brunt,shoc_mix);       // Output

    // Advance the SGS TKE equation
    shoc_tke(team,nlev,nlevi,dtime,              // Input
	     lambda_low,lambda_high,lambda_slope, // Runtime options
	     lambda_thresh,Ckh,Ckm,              // Runtime options
	     wthv_sec,                           // Input
             shoc_mix,dz_zi,dz_zt,pres,shoc_tabs,// Input
             u_wind,v_wind,brunt,zt_grid,        // Input
             zi_grid,pblh,                       // Input
             workspace,                          // Workspace
             tke,tk,tkh,                         // Input/Output
             isotropy);                          // Output

    // Update SHOC prognostic variables here
    // via implicit diffusion solver
    team.team_barrier();
    update_prognostics_implicit(team,nlev,nlevi,num_qtracers,dtime,dz_zt,   // Input
                                dz_zi,rho_zt,zt_grid,zi_grid,tk,tkh,uw_sfc, // Input
                                vw_sfc,wthl_sfc,wqw_sfc,wtracer_sfc,        // Input
                                workspace,                                  // Workspace
                                thetal,qw,qtracers,tke,u_wind,v_wind);   // Input/Output

    // Diagnose the second order moments
    diag_second_shoc_moments(team,nlev,nlevi,
                             thl2tune, qw2tune, qwthl2tune, w2tune,     // Runtime options
                             thetal,qw,u_wind,v_wind,                   // Input
                             tke,isotropy,tkh,tk,dz_zi,zt_grid,zi_grid, // Input
                             shoc_mix,wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,   // Input
                             ustar2,wstar,                              // Input/Output
                             workspace,                                 // Workspace
                             thl_sec,qw_sec,wthl_sec,wqw_sec,qwthl_sec, // Output
                             uw_sec,vw_sec,wtke_sec,w_sec);             // Output

    // Diagnose the third moment of vertical velocity,
    //  needed for the PDF closure
    diag_third_shoc_moments(team,nlev,nlevi,
                            c_diag_3rd_mom,                         // Runtime options
                            w_sec,thl_sec,wthl_sec,                 // Input
                            isotropy,brunt,thetal,tke,dz_zt,dz_zi,  // Input
                            zt_grid,zi_grid,                        // Input
                            workspace,                              // Workspace
                            w3);                                    // Output

    // Call the PDF to close on SGS cloud and turbulence
    team.team_barrier();
    shoc_assumed_pdf(team,nlev,nlevi,dtime,thetal,qw,w_field,thl_sec,qw_sec, // Input
                     wthl_sec,w_sec,wqw_sec,qwthl_sec,w3,pres,         // Input
                     zt_grid, zi_grid,                                 // Input
                     workspace,                                        // Workspace
                     shoc_cldfrac,shoc_ql,wqls_sec,wthv_sec,shoc_ql2,shoc_cond,shoc_evap); // Ouptut

    // Check TKE to make sure values lie within acceptable
    // bounds after vertical advection, etc.
    check_tke(team,nlev,tke);
  }

  // End SHOC parameterization

  // Use SHOC outputs to update the host model
  // temperature
  update_host_dse(team,nlev,thetal,shoc_ql, // Input
                  inv_exner,zt_grid,phis,   // Input
                  host_dse);                // Output

  team.team_barrier();
  shoc_energy_integrals(team,nlev,host_dse,pdel,  // Input
                        qw,shoc_ql,u_wind,v_wind, // Input
                        se_a,ke_a,wv_a,wl_a);     // Output

  shoc_energy_fixer(team,nlev,nlevi,dtime,nadv,zt_grid,zi_grid, // Input
                    se_b,ke_b,wv_b,wl_b,se_a,ke_a,wv_a,wl_a,    // Input
                    wthl_sfc,wqw_sfc,rho_zt,tke,presi,          // Input
                    workspace,                                  // Workspace
                    host_dse);                                  // Output

  // Remaining code is to diagnose certain quantities
  // related to PBL.  No answer changing subroutines
  // should be placed at this point onward.

  // Update PBLH, as other routines outside of SHOC
  // may require this variable.

  // Update SHOC water vapor, to be used by the next two routines
  compute_shoc_vapor(team,nlev,qw,shoc_ql, // Input
                     shoc_qv);             // Output

  team.team_barrier();
  shoc_diag_obklen(uw_sfc,vw_sfc,      // Input
                   wthl_sfc,wqw_sfc,   // Input
                   s_thetal(nlev-1),   // Input
                   s_shoc_ql(nlev-1),  // Input
                   s_shoc_qv(nlev-1),  // Input
                   ustar,kbfs,obklen); // Output

  pblintd(team,nlev,nlevi,npbl,zt_grid,   // Input
          zi_grid,thetal,shoc_ql,shoc_qv, // Input
          u_wind,v_wind,ustar,obklen,     // Input
          kbfs,shoc_cldfrac,              // Input
          workspace,                      // Workspace
          pblh);                          // Output

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<5>(
    {&rho_zt, &shoc_qv, &shoc_tabs, &dz_zt, &dz_zi});
}
#else
template<typename S, typename D>
void Functions<S,D>::shoc_main_internal(
  const Int&                   shcol,        // Number of columns
  const Int&                   nlev,         // Number of levels
  const Int&                   nlevi,        // Number of levels on interface grid
  const Int&                   npbl,         // Maximum number of levels in pbl from surface
  const Int&                   nadv,         // Number of times to loop SHOC
  const Int&                   num_qtracers, // Number of tracers
  const Scalar&                dtime,        // SHOC timestep [s]
  // Runtime Parameters
  const Scalar&                lambda_low,
  const Scalar&                lambda_high,
  const Scalar&                lambda_slope,
  const Scalar&                lambda_thresh,
  const Scalar&                thl2tune,
  const Scalar&                qw2tune,
  const Scalar&                qwthl2tune,
  const Scalar&                w2tune,
  const Scalar&                length_fac,
  const Scalar&                c_diag_3rd_mom,
  const Scalar&                Ckh,
  const Scalar&                Ckm,
  // Input Variables
  const view_1d<const Scalar>& dx,
  const view_1d<const Scalar>& dy,
  const view_2d<const Spack>& zt_grid,
  const view_2d<const Spack>& zi_grid,
  const view_2d<const Spack>& pres,
  const view_2d<const Spack>& presi,
  const view_2d<const Spack>& pdel,
  const view_2d<const Spack>& thv,
  const view_2d<const Spack>& w_field,
  const view_1d<const Scalar>& wthl_sfc,
  const view_1d<const Scalar>& wqw_sfc,
  const view_1d<const Scalar>& uw_sfc,
  const view_1d<const Scalar>& vw_sfc,
  const view_2d<const Spack>& wtracer_sfc,
  const view_2d<const Spack>& inv_exner,
  const view_1d<const Scalar>& phis,
  // Workspace Manager
  WorkspaceMgr&      workspace_mgr,
  // Input/Output Variables
  const view_2d<Spack>&       host_dse,
  const view_2d<Spack>&       tke,
  const view_2d<Spack>&       thetal,
  const view_2d<Spack>&       qw,
  const uview_2d<Spack>&      u_wind,
  const uview_2d<Spack>&      v_wind,
  const view_2d<Spack>&       wthv_sec,
  const view_3d_strided<Spack>& qtracers,
  const view_2d<Spack>&       tk,
  const view_2d<Spack>&       shoc_cldfrac,
  const view_2d<Spack>&       shoc_ql,
  // Output Variables
  const view_1d<Scalar>&      pblh,
  const view_1d<Scalar>&      ustar,
  const view_1d<Scalar>&      obklen,
  const view_2d<Spack>&       shoc_ql2,
  const view_2d<Spack>&       tkh,
  // Diagnostic Output Variables
  const view_2d<Spack>&       shoc_mix,
  const view_2d<Spack>&       shoc_cond,
  const view_2d<Spack>&       shoc_evap,
  const view_2d<Spack>&       w_sec,
  const view_2d<Spack>&       thl_sec,
  const view_2d<Spack>&       qw_sec,
  const view_2d<Spack>&       qwthl_sec,
  const view_2d<Spack>&       wthl_sec,
  const view_2d<Spack>&       wqw_sec,
  const view_2d<Spack>&       wtke_sec,
  const view_2d<Spack>&       uw_sec,
  const view_2d<Spack>&       vw_sec,
  const view_2d<Spack>&       w3,
  const view_2d<Spack>&       wqls_sec,
  const view_2d<Spack>&       brunt,
  const view_2d<Spack>&       isotropy,
  // Temporaries
  const view_1d<Scalar>& se_b,
  const view_1d<Scalar>& ke_b,
  const view_1d<Scalar>& wv_b,
  const view_1d<Scalar>& wl_b,
  const view_1d<Scalar>& se_a,
  const view_1d<Scalar>& ke_a,
  const view_1d<Scalar>& wv_a,
  const view_1d<Scalar>& wl_a,
  const view_1d<Scalar>& kbfs,
  const view_1d<Scalar>& ustar2,
  const view_1d<Scalar>& wstar,
  const view_2d<Spack>& rho_zt,
  const view_2d<Spack>& shoc_qv,
  const view_2d<Spack>& shoc_tabs,
  const view_2d<Spack>& dz_zt,
  const view_2d<Spack>& dz_zi)
{
  // Scalarize some views for single entry access
  const auto s_thetal  = ekat::scalarize(thetal);
  const auto s_shoc_ql = ekat::scalarize(shoc_ql);
  const auto s_shoc_qv = ekat::scalarize(shoc_qv);

  // Compute integrals of static energy, kinetic energy, water vapor, and liquid water
  // for the computation of total energy before SHOC is called.  This is for an
  // effort to conserve energy since liquid water potential temperature (which SHOC
  // conserves) and static energy (which E3SM conserves) are not exactly equal.
  shoc_energy_integrals_disp(shcol,nlev,host_dse,pdel,qw,shoc_ql,u_wind,v_wind,
                             se_b, ke_b, wv_b, wl_b); // Input

  for (Int t=0; t<nadv; ++t) {
    // Check TKE to make sure values lie within acceptable
    // bounds after host model performs horizontal advection
    check_tke_disp(shcol,nlev, // Input
                   tke);      // Input/Output

    // Define vertical grid arrays needed for
    // vertical derivatives in SHOC, also
    // define air density (rho_zt)
    shoc_grid_disp(shcol,nlev,nlevi,      // Input
                   zt_grid,zi_grid,pdel, // Input
                   dz_zt,dz_zi,rho_zt);  // Output

    // Compute the planetary boundary layer height, which is an
    // input needed for the length scale calculation.

    // Update SHOC water vapor,
    // to be used by the next two routines
    compute_shoc_vapor_disp(shcol,nlev,qw,shoc_ql, // Input
                            shoc_qv);             // Output

    // Update SHOC temperature
    compute_shoc_temperature_disp(shcol,nlev,thetal,  // Input
                                  shoc_ql,inv_exner, // Input
                                  shoc_tabs);        // Output

    shoc_diag_obklen_disp(shcol, nlev,
                          uw_sfc,vw_sfc,     // Input
                          wthl_sfc, wqw_sfc, // Input
                          s_thetal,  // Input
                          s_shoc_ql, // Input
                          s_shoc_qv, // Input
                          ustar,kbfs,obklen); // Output

    pblintd_disp(shcol,nlev,nlevi,npbl,    // Input
                 zt_grid,zi_grid,thetal,   // Input
                 shoc_ql,shoc_qv,u_wind,   // Input
                 v_wind,ustar,obklen,kbfs, // Input
                 shoc_cldfrac,             // Input
                 workspace_mgr,            // Workspace mgr
                 pblh);                    // Output

    // Update the turbulent length scale
    shoc_length_disp(shcol,nlev,nlevi,      // Input
                     length_fac,            // Runtime Options
                     dx,dy,                 // Input
                     zt_grid,zi_grid,dz_zt, // Input
                     tke,thv,               // Input
                     workspace_mgr,         // Workspace mgr
                     brunt,shoc_mix);       // Output

    // Advance the SGS TKE equation
    shoc_tke_disp(shcol,nlev,nlevi,dtime,             // Input
	          lambda_low,lambda_high,lambda_slope, // Runtime options
		  lambda_thresh,Ckh,Ckm,              // Runtime options
                  wthv_sec,                           // Input
                  shoc_mix,dz_zi,dz_zt,pres,shoc_tabs,// Input
                  u_wind,v_wind,brunt,zt_grid,        // Input
                  zi_grid,pblh,                       // Input
                  workspace_mgr,                      // Workspace mgr
                  tke,tk,tkh,                         // Input/Output
                  isotropy);                          // Output

    // Update SHOC prognostic variables here
    // via implicit diffusion solver
    update_prognostics_implicit_disp(shcol,nlev,nlevi,num_qtracers,dtime,dz_zt,  // Input
                                     dz_zi,rho_zt,zt_grid,zi_grid,tk,tkh,uw_sfc, // Input
                                     vw_sfc,wthl_sfc,wqw_sfc,wtracer_sfc,        // Input
                                     workspace_mgr,                              // Workspace mgr
                                     thetal,qw,qtracers,tke,u_wind,v_wind);      // Input/Output

    // Diagnose the second order moments
    diag_second_shoc_moments_disp(shcol,nlev,nlevi,
                                  thl2tune, qw2tune, qwthl2tune, w2tune,     // Runtime options
                                  thetal,qw,u_wind,v_wind,                   // Input
                                  tke,isotropy,tkh,tk,dz_zi,zt_grid,zi_grid, // Input
                                  shoc_mix,wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,   // Input
                                  ustar2,wstar,                              // Input/Output
                                  workspace_mgr,                             // Workspace
                                  thl_sec,qw_sec,wthl_sec,wqw_sec,qwthl_sec, // Output
                                  uw_sec,vw_sec,wtke_sec,w_sec);             // Output

    // Diagnose the third moment of vertical velocity,
    //  needed for the PDF closure
    diag_third_shoc_moments_disp(shcol,nlev,nlevi,
                                 c_diag_3rd_mom,                         // Runtime options
                                 w_sec,thl_sec,wthl_sec,                 // Input
                                 isotropy,brunt,thetal,tke,dz_zt,dz_zi,  // Input
                                 zt_grid,zi_grid,                        // Input
                                 workspace_mgr,                          // Workspace mgr
                                 w3);                                    // Output

    // Call the PDF to close on SGS cloud and turbulence
    shoc_assumed_pdf_disp(shcol,nlev,nlevi,thetal,qw,w_field,thl_sec,qw_sec, // Input
                          wthl_sec,w_sec,wqw_sec,qwthl_sec,w3,pres,         // Input
                          zt_grid, zi_grid,                                 // Input
                          workspace_mgr,                                    // Workspace mgr
                          shoc_cldfrac,shoc_ql,wqls_sec,wthv_sec,shoc_ql2); // Ouptut

    // Check TKE to make sure values lie within acceptable
    // bounds after vertical advection, etc.
    check_tke_disp(shcol,nlev,tke);
  }

  // End SHOC parameterization

  // Use SHOC outputs to update the host model
  // temperature
  update_host_dse_disp(shcol,nlev,thetal,shoc_ql, // Input
                       inv_exner,zt_grid,phis,   // Input
                       host_dse);                // Output

  shoc_energy_integrals_disp(shcol,nlev,host_dse,pdel,  // Input
                        qw,shoc_ql,u_wind,v_wind, // Input
                        se_a,ke_a,wv_a,wl_a);     // Output

  shoc_energy_fixer_disp(shcol,nlev,nlevi,dtime,nadv,zt_grid,zi_grid, // Input
                         se_b,ke_b,wv_b,wl_b,se_a,ke_a,wv_a,wl_a,    // Input
                         wthl_sfc,wqw_sfc,rho_zt,tke,presi,          // Input
                         workspace_mgr,                              // Workspace
                         host_dse);                                  // Output

  // Remaining code is to diagnose certain quantities
  // related to PBL.  No answer changing subroutines
  // should be placed at this point onward.

  // Update PBLH, as other routines outside of SHOC
  // may require this variable.

  // Update SHOC water vapor, to be used by the next two routines
  compute_shoc_vapor_disp(shcol,nlev,qw,shoc_ql, // Input
                          shoc_qv);             // Output

  shoc_diag_obklen_disp(shcol, nlev, uw_sfc,vw_sfc,      // Input
                        wthl_sfc,wqw_sfc,   // Input
                        s_thetal,   // Input
                        s_shoc_ql,  // Input
                        s_shoc_qv,  // Input
                        ustar,kbfs,obklen); // Output

  pblintd_disp(shcol,nlev,nlevi,npbl,zt_grid,   // Input
               zi_grid,thetal,shoc_ql,shoc_qv, // Input
               u_wind,v_wind,ustar,obklen,     // Input
               kbfs,shoc_cldfrac,              // Input
               workspace_mgr,                  // Workspace mgr
               pblh);                          // Output
}
#endif

template<typename S, typename D>
Int Functions<S,D>::shoc_main(
  const Int&               shcol,               // Number of SHOC columns in the array
  const Int&               nlev,                // Number of levels
  const Int&               nlevi,               // Number of levels on interface grid
  const Int&               npbl,                // Maximum number of levels in pbl from surface
  const Int&               nadv,                // Number of times to loop SHOC
  const Int&               num_qtracers,        // Number of tracers
  const Scalar&            dtime,               // SHOC timestep [s]
  WorkspaceMgr&            workspace_mgr,       // WorkspaceManager for local variables
  const SHOCRuntime&       shoc_runtime,        // Runtime Options
  const SHOCInput&         shoc_input,          // Input
  const SHOCInputOutput&   shoc_input_output,   // Input/Output
  const SHOCOutput&        shoc_output,         // Output
  const SHOCHistoryOutput& shoc_history_output  // Output (diagnostic)
#ifdef SCREAM_SHOC_SMALL_KERNELS
  , const SHOCTemporaries& shoc_temporaries     // Temporaries for small kernels
#endif
                              )
{
  // Start timer
  auto start = std::chrono::steady_clock::now();

  // Runtime options
  const Scalar lambda_low    = shoc_runtime.lambda_low;
  const Scalar lambda_high   = shoc_runtime.lambda_high;
  const Scalar lambda_slope  = shoc_runtime.lambda_slope;
  const Scalar lambda_thresh = shoc_runtime.lambda_thresh;
  const Scalar thl2tune      = shoc_runtime.thl2tune;
  const Scalar qw2tune       = shoc_runtime.qw2tune;
  const Scalar qwthl2tune    = shoc_runtime.qwthl2tune;
  const Scalar w2tune        = shoc_runtime.w2tune;
  const Scalar length_fac    = shoc_runtime.length_fac;
  const Scalar c_diag_3rd_mom = shoc_runtime.c_diag_3rd_mom;
  const Scalar Ckh           = shoc_runtime.Ckh;
  const Scalar Ckm           = shoc_runtime.Ckm;

#ifndef SCREAM_SHOC_SMALL_KERNELS
  using ExeSpace = typename KT::ExeSpace;

  // SHOC main loop
  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const Scalar dx_s{shoc_input.dx(i)};
    const Scalar dy_s{shoc_input.dy(i)};
    const Scalar wthl_sfc_s{shoc_input.wthl_sfc(i)};
    const Scalar wqw_sfc_s{shoc_input.wqw_sfc(i)};
    const Scalar uw_sfc_s{shoc_input.uw_sfc(i)};
    const Scalar vw_sfc_s{shoc_input.vw_sfc(i)};
    const Scalar phis_s{shoc_input.phis(i)};
    Scalar pblh_s{0};
    Scalar ustar_s{0};
    Scalar obklen_s{0};

    const auto zt_grid_s      = ekat::subview(shoc_input.zt_grid, i);
    const auto zi_grid_s      = ekat::subview(shoc_input.zi_grid, i);
    const auto pres_s         = ekat::subview(shoc_input.pres, i);
    const auto presi_s        = ekat::subview(shoc_input.presi, i);
    const auto pdel_s         = ekat::subview(shoc_input.pdel, i);
    const auto thv_s          = ekat::subview(shoc_input.thv, i);
    const auto w_field_s      = ekat::subview(shoc_input.w_field, i);
    const auto wtracer_sfc_s  = ekat::subview(shoc_input.wtracer_sfc, i);
    const auto inv_exner_s    = ekat::subview(shoc_input.inv_exner, i);
    const auto host_dse_s     = ekat::subview(shoc_input_output.host_dse, i);
    const auto tke_s          = ekat::subview(shoc_input_output.tke, i);
    const auto thetal_s       = ekat::subview(shoc_input_output.thetal, i);
    const auto qw_s           = ekat::subview(shoc_input_output.qw, i);
    const auto wthv_sec_s     = ekat::subview(shoc_input_output.wthv_sec, i);
    const auto tk_s           = ekat::subview(shoc_input_output.tk, i);
    const auto shoc_cldfrac_s = ekat::subview(shoc_input_output.shoc_cldfrac, i);
    const auto shoc_ql_s      = ekat::subview(shoc_input_output.shoc_ql, i);
    const auto shoc_ql2_s     = ekat::subview(shoc_output.shoc_ql2, i);
    const auto tkh_s          = ekat::subview(shoc_output.tkh, i);
    const auto shoc_mix_s     = ekat::subview(shoc_history_output.shoc_mix, i);
    const auto shoc_cond_s     = ekat::subview(shoc_history_output.shoc_cond, i);
    const auto shoc_evap_s     = ekat::subview(shoc_history_output.shoc_evap, i);
    const auto w_sec_s        = ekat::subview(shoc_history_output.w_sec, i);
    const auto thl_sec_s      = ekat::subview(shoc_history_output.thl_sec, i);
    const auto qw_sec_s       = ekat::subview(shoc_history_output.qw_sec, i);
    const auto qwthl_sec_s    = ekat::subview(shoc_history_output.qwthl_sec, i);
    const auto wthl_sec_s     = ekat::subview(shoc_history_output.wthl_sec, i);
    const auto wqw_sec_s      = ekat::subview(shoc_history_output.wqw_sec, i);
    const auto wtke_sec_s     = ekat::subview(shoc_history_output.wtke_sec, i);
    const auto uw_sec_s       = ekat::subview(shoc_history_output.uw_sec, i);
    const auto vw_sec_s       = ekat::subview(shoc_history_output.vw_sec, i);
    const auto w3_s           = ekat::subview(shoc_history_output.w3, i);
    const auto wqls_sec_s     = ekat::subview(shoc_history_output.wqls_sec, i);
    const auto brunt_s        = ekat::subview(shoc_history_output.brunt, i);
    const auto isotropy_s     = ekat::subview(shoc_history_output.isotropy, i);

    const auto u_wind_s   = Kokkos::subview(shoc_input_output.horiz_wind, i, 0, Kokkos::ALL());
    const auto v_wind_s   = Kokkos::subview(shoc_input_output.horiz_wind, i, 1, Kokkos::ALL());
    const auto qtracers_s = Kokkos::subview(shoc_input_output.qtracers, i, Kokkos::ALL(), Kokkos::ALL());

    shoc_main_internal(team, nlev, nlevi, npbl, nadv, num_qtracers, dtime,
	               lambda_low, lambda_high, lambda_slope, lambda_thresh,  // Runtime options
                       thl2tune, qw2tune, qwthl2tune, w2tune, length_fac,     // Runtime options
                       c_diag_3rd_mom, Ckh, Ckm,                              // Runtime options
                       dx_s, dy_s, zt_grid_s, zi_grid_s,                      // Input
                       pres_s, presi_s, pdel_s, thv_s, w_field_s,             // Input
                       wthl_sfc_s, wqw_sfc_s, uw_sfc_s, vw_sfc_s,             // Input
                       wtracer_sfc_s, inv_exner_s, phis_s,                    // Input
                       workspace,                                             // Workspace
                       host_dse_s, tke_s, thetal_s, qw_s, u_wind_s, v_wind_s, // Input/Output
                       wthv_sec_s, qtracers_s, tk_s, shoc_cldfrac_s,          // Input/Output
                       shoc_ql_s,                                             // Input/Output
                       pblh_s, ustar_s, obklen_s, shoc_ql2_s, tkh_s,          // Output
                       shoc_mix_s,shoc_cond_s, shoc_evap_s, w_sec_s, thl_sec_s, qw_sec_s, qwthl_sec_s, // Diagnostic Output Variables
                       wthl_sec_s, wqw_sec_s, wtke_sec_s, uw_sec_s, vw_sec_s, // Diagnostic Output Variables
                       w3_s, wqls_sec_s, brunt_s, isotropy_s);                // Diagnostic Output Variables

    shoc_output.pblh(i) = pblh_s;
    shoc_output.ustar(i) = ustar_s;
    shoc_output.obklen(i) = obklen_s;
  });
  Kokkos::fence();
#else
  const auto u_wind_s   = Kokkos::subview(shoc_input_output.horiz_wind, Kokkos::ALL(), 0, Kokkos::ALL());
  const auto v_wind_s   = Kokkos::subview(shoc_input_output.horiz_wind, Kokkos::ALL(), 1, Kokkos::ALL());

  shoc_main_internal(shcol, nlev, nlevi, npbl, nadv, num_qtracers, dtime,
    lambda_low, lambda_high, lambda_slope, lambda_thresh,  // Runtime options
    thl2tune, qw2tune, qwthl2tune, w2tune, length_fac,     // Runtime options
    c_diag_3rd_mom, Ckh, Ckm,                              // Runtime options
    shoc_input.dx, shoc_input.dy, shoc_input.zt_grid, shoc_input.zi_grid, // Input
    shoc_input.pres, shoc_input.presi, shoc_input.pdel, shoc_input.thv, shoc_input.w_field, // Input
    shoc_input.wthl_sfc, shoc_input.wqw_sfc, shoc_input.uw_sfc, shoc_input.vw_sfc, // Input
    shoc_input.wtracer_sfc, shoc_input.inv_exner, shoc_input.phis, // Input
    workspace_mgr, // Workspace Manager
    shoc_input_output.host_dse, shoc_input_output.tke, shoc_input_output.thetal, shoc_input_output.qw, u_wind_s, v_wind_s, // Input/Output
    shoc_input_output.wthv_sec, shoc_input_output.qtracers, shoc_input_output.tk, shoc_input_output.shoc_cldfrac, // Input/Output
    shoc_input_output.shoc_ql, // Input/Output
    shoc_output.pblh, shoc_output.ustar, shoc_output.obklen, shoc_output.shoc_ql2, shoc_output.tkh, // Output
    shoc_history_output.shoc_mix, shoc_history_output.shoc_cond, shoc_history_output.shoc_evap, shoc_history_output.w_sec, shoc_history_output.thl_sec, shoc_history_output.qw_sec, shoc_history_output.qwthl_sec, // Diagnostic Output Variables
    shoc_history_output.wthl_sec, shoc_history_output.wqw_sec, shoc_history_output.wtke_sec, shoc_history_output.uw_sec, shoc_history_output.vw_sec, // Diagnostic Output Variables
    shoc_history_output.w3, shoc_history_output.wqls_sec, shoc_history_output.brunt, shoc_history_output.isotropy, // Diagnostic Output Variables
    // Temporaries
    shoc_temporaries.se_b, shoc_temporaries.ke_b, shoc_temporaries.wv_b, shoc_temporaries.wl_b,
    shoc_temporaries.se_a, shoc_temporaries.ke_a, shoc_temporaries.wv_a, shoc_temporaries.wl_a,
    shoc_temporaries.kbfs, shoc_temporaries.ustar2,
    shoc_temporaries.wstar, shoc_temporaries.rho_zt, shoc_temporaries.shoc_qv,
    shoc_temporaries.tabs, shoc_temporaries.dz_zt, shoc_temporaries.dz_zi);
#endif

  auto finish = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
  return duration.count();
}

} // namespace shoc
} // namespace scream

#endif
