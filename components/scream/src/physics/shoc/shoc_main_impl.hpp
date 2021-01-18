#ifndef SHOC_MAIN_IMPL_HPP
#define SHOC_MAIN_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_main. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

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
  // Input Variables
  const Scalar&                host_dx,
  const Scalar&                host_dy,
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
  const uview_1d<const Spack>& exner,
  const Scalar&                phis,
  // Local Variables
  const uview_1d<Spack>&       rho_zt,
  const uview_1d<Spack>&       shoc_qv,
  const uview_1d<Spack>&       dz_zt,
  const uview_1d<Spack>&       dz_zi,
  const uview_1d<Spack>&       thv_zi,
  const uview_1d<Spack>&       sterm,
  const uview_1d<Spack>&       sterm_zt,
  const uview_1d<Spack>&       a_diss,
  const uview_1d<Spack>&       rdp_zt,
  const uview_1d<Spack>&       tmpi,
  const uview_1d<Spack>&       tkh_zi,
  const uview_1d<Spack>&       tk_zi,
  const uview_1d<Spack>&       rho_zi,
  const uview_1d<Scalar>&      du,
  const uview_1d<Scalar>&      dl,
  const uview_1d<Scalar>&      d,
  const uview_2d<Spack>&       X1,
  const uview_1d<Spack>&       isotropy_zi,
  const uview_1d<Spack>&       w_sec_zi,
  const uview_1d<Spack>&       brunt_zi,
  const uview_1d<Spack>&       thetal_zi,
  const uview_1d<Spack>&       wthl_sec_zt,
  const uview_1d<Spack>&       wqw_sec_zt,
  const uview_1d<Spack>&       w3_zt,
  const uview_1d<Spack>&       thl_sec_zt,
  const uview_1d<Spack>&       qwthl_sec_zt,
  const uview_1d<Spack>&       qw_sec_zt,
  const uview_1d<Spack>&       pblintd_thv,
  const uview_1d<Spack>&       rino,
  // Input/Output Variables
  const uview_1d<Spack>&       host_dse,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       thetal,
  const uview_1d<Spack>&       qw,
  const uview_1d<Spack>&       u_wind,
  const uview_1d<Spack>&       v_wind,
  const uview_1d<Spack>&       wthv_sec,
  const uview_2d<Spack>&       qtracers,
  const uview_1d<Spack>&       tk,
  const uview_1d<Spack>&       tkh,
  const uview_1d<Spack>&       shoc_cldfrac,
  const uview_1d<Spack>&       shoc_ql,
  // Output Variables
  Scalar&                      pblh,
  const uview_1d<Spack>&       shoc_ql2,
  // Diagnostic Output Variables
  const uview_1d<Spack>&       shoc_mix,
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
  // View/pack indices for nlev, nlevi
  const Int nlev_v = (nlev-1)/Spack::n;
  const Int nlev_p = (nlev-1)%Spack::n;

  // Local variables
  Scalar se_b{0},   ke_b{0}, wv_b{0},   wl_b{0},
         se_a{0},   ke_a{0}, wv_a{0},   wl_a{0},
         ustar{0},  kbfs{0}, obklen{0}, ustar2{0}, wstar{0};

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
    // define air density
    shoc_grid(team,nlev,nlevi,      // Input
              zt_grid,zi_grid,pdel, // Input
              dz_zt,dz_zi,rho_zt);  // Output

    // Compute the planetary boundary layer height, which is an
    // input needed for the length scale calculation.

    // Update SHOC water vapor, to be used by the next two routines
    compute_shoc_vapor(team,nlev,qw,shoc_ql, // Input
                       shoc_qv);             // Output

    team.team_barrier();
    shoc_diag_obklen(uw_sfc,vw_sfc,          // Input
                     wthl_sfc, wqw_sfc,      // Input
                    thetal(nlev_v)[nlev_p],  // Input
                    shoc_ql(nlev_v)[nlev_p], // Input
                    shoc_qv(nlev_v)[nlev_p], // Input
                    ustar,kbfs,obklen);      // Output

    pblintd(team,nlev,nlevi,npbl,          // Input
            zt_grid,zi_grid,thetal,        // Input
            shoc_ql,shoc_qv,u_wind,        // Input
            v_wind,ustar,obklen,kbfs,      // Input
            shoc_cldfrac,rino,pblintd_thv, // Input
            pblh);                         // Output

    // Update the turbulent length scale
    shoc_length(team,nlev,nlevi,host_dx,host_dy, // Input
                pblh,tke,zt_grid,zi_grid,dz_zt,  // Input
                wthv_sec,thv,thv_zi,             // Input
                brunt,shoc_mix);                 // Output

    // Advance the SGS TKE equation
    shoc_tke(team,nlev,nlevi,dtime,wthv_sec,      // Input
             shoc_mix,dz_zi,dz_zt,pres,u_wind,    // Input
             v_wind,brunt,obklen,zt_grid,zi_grid, // Input
             pblh,sterm,sterm_zt,a_diss,          // Input
             tke,tk,tkh,                          // Input/Output
             isotropy);                           // Output

    // Update SHOC prognostic variables here
    // via implicit diffusion solver
    team.team_barrier();
    update_prognostics_implicit(team,nlev,nlevi,num_qtracers,dtime,dz_zt,   // Input
                                dz_zi,rho_zt,zt_grid,zi_grid,tk,tkh,uw_sfc, // Input
                                vw_sfc,wthl_sfc,wqw_sfc,wtracer_sfc,rdp_zt, // Input
                                tmpi,tkh_zi,tk_zi,rho_zi,du,dl,d,X1,        // Input
                               thetal, qw,qtracers,tke,u_wind,v_wind);      // Input/Output

    // Diagnose the second order moments
    diag_second_shoc_moments(team,nlev,nlevi,thetal,qw,u_wind,v_wind,  // Input
                            tke,isotropy,tkh,tk,dz_zi,zt_grid,zi_grid, // Input
                            shoc_mix,wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,   // Input
                            ustar2,wstar,isotropy_zi,tkh_zi,tk_zi,     // Input
                            thl_sec,qw_sec,wthl_sec,wqw_sec,qwthl_sec, // Output
                            uw_sec,vw_sec,wtke_sec,w_sec);             // Output

    // Diagnose the third moment of vertical velocity,
    //  needed for the PDF closure
    diag_third_shoc_moments(team,nlev,nlevi,w_sec,thl_sec,wthl_sec, // Input
                            isotropy,brunt,thetal,tke,dz_zt,dz_zi,  // Input
                            zt_grid,zi_grid,w_sec_zi,isotropy_zi,   // Input
                            brunt_zi,thetal_zi,                     // Input
                            w3);                                    // Output

    // Call the PDF to close on SGS cloud and turbulence
    team.team_barrier();
    shoc_assumed_pdf(team,nlev,nlevi,thetal,qw,w_field,thl_sec,qw_sec, // Input
                     wthl_sec,w_sec,wqw_sec,qwthl_sec,w3,pres,zt_grid, // Input
                     zi_grid,wthl_sec_zt,wqw_sec_zt,w3_zt,thl_sec_zt,  // Input
                     qwthl_sec_zt,qw_sec_zt,                           // Input
                    shoc_cldfrac,shoc_ql,wqls_sec,wthv_sec,shoc_ql2);  // Ouptut

    // Check TKE to make sure values lie within acceptable
    // bounds after vertical advection, etc.
    check_tke(team,nlev,tke);
  }

  // End SHOC parameterization

  // Use SHOC outputs to update the host model
  // temperature
  update_host_dse(team,nlev,thetal,shoc_ql, // Input
                  exner,zt_grid,phis,       // Input
                  host_dse);                // Output

  team.team_barrier();
  shoc_energy_integrals(team,nlev,host_dse,pdel,  // Input
                        qw,shoc_ql,u_wind,v_wind, // Input
                        se_a,ke_a,wv_a,wl_a);     // Output

  shoc_energy_fixer(team,nlev,nlevi,dtime,nadv,zt_grid,zi_grid, // Input
                    se_b,ke_b,wv_b,wl_b,se_a,ke_a,wv_a,wl_a,    // Input
                    wthl_sfc,wqw_sfc,rho_zt,tke,presi,rho_zi,   // Input
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
  shoc_diag_obklen(uw_sfc,vw_sfc,           // Input
                   wthl_sfc,wqw_sfc,        // Input
                   thetal(nlev_v)[nlev_p],  // Input
                   shoc_ql(nlev_v)[nlev_p], // Input
                   shoc_qv(nlev_v)[nlev_p], // Input
                   ustar,kbfs,obklen);      // Output

  pblintd(team,nlev,nlevi,npbl,zt_grid,       // Input
          zi_grid,thetal,shoc_ql,shoc_qv,     // Input
          u_wind,v_wind,ustar,obklen,         // Input
          kbfs,shoc_cldfrac,rino,pblintd_thv, // Input
          pblh);                              // Output
}


template<typename S, typename D>
Int Functions<S,D>::shoc_main(
  const Int&               shcol,        // Number of SHOC columns in the array
  const Int&               nlev,         // Number of levels
  const Int&               nlevi,        // Number of levels on interface grid
  const Int&               npbl,         // Maximum number of levels in pbl from surface
  const Int&               nadv,         // Number of times to loop SHOC
  const Int&               num_qtracers, // Number of tracers
  const Scalar&            dtime,        // SHOC timestep [s]
  const SHOCInput&         shoc_input,
  const SHOCInputOutput&   shoc_input_output,
  const SHOCOutput&        shoc_output,
  const SHOCHistoryOutput& shoc_history_output)
{
  using ExeSpace = typename KT::ExeSpace;

  // Number of packs for nlev, nlevi
  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto nlevi_packs = ekat::npack<Spack>(nlevi);

  // Local variables
  view_2d<Spack>
    rho_zt_d("rho_zt",shcol,nlev_packs),
    shoc_qv_d("rho_zt",shcol,nlev_packs),
    dz_zt_d("rho_zt",shcol,nlev_packs),
    dz_zi_d("rho_zt",shcol,nlevi_packs),
    thv_zi_d("thv_zi",shcol,nlevi_packs),
    sterm_d("sterm",shcol,nlevi_packs),
    sterm_zt_d("sterm_zt",shcol,nlev_packs),
    a_diss_d("a_diss",shcol,nlev_packs),
    rdp_zt_d("rdp_zt", shcol, nlev_packs),
    tmpi_d("tmpi", shcol, nlevi_packs),
    tkh_zi_d("tkh_zi", shcol, nlevi_packs),
    tk_zi_d("tk_zi", shcol, nlevi_packs),
    rho_zi_d("rho_zi", shcol, nlevi_packs),
    isotropy_zi_d("isotropy_zi", shcol, nlevi_packs),
    w_sec_zi_d("w_sec_zi", shcol, nlevi_packs),
    brunt_zi_d("brunt_zi", shcol, nlevi_packs),
    thetal_zi_d("thetal_zi", shcol, nlevi_packs),
    wthl_sec_zt_d("wthl_sec_zt", shcol, nlev_packs),
    wqw_sec_zt_d("wqw_sec_zt", shcol, nlev_packs),
    w3_zt_d("w3_zt", shcol, nlev_packs),
    thl_sec_zt_d("thl_sec_zt", shcol, nlev_packs),
    qwthl_sec_zt_d("qwthl_sec_zt", shcol, nlev_packs),
    qw_sec_zt_d("qw_sec_zt", shcol, nlev_packs),
    pblintd_thv_d("pblintd_thv", shcol, nlev_packs),
    pblintd_rino_d("rino", shcol, nlev_packs);

  view_2d<Scalar>
    du_d("du", shcol, nlev),
    dl_d("dl", shcol, nlev),
    d_d("d", shcol, nlev);

  view_3d<Spack>
    X1_d("X1",shcol,nlev,ekat::npack<Spack>(2));

  // Start timer
  auto start = std::chrono::steady_clock::now();

  // SHOC main loop
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    const Scalar host_dx_s{shoc_input.host_dx(i)[0]};
    const Scalar host_dy_s{shoc_input.host_dy(i)[0]};
    const Scalar wthl_sfc_s{shoc_input.wthl_sfc(i)[0]};
    const Scalar wqw_sfc_s{shoc_input.wqw_sfc(i)[0]};
    const Scalar uw_sfc_s{shoc_input.uw_sfc(i)[0]};
    const Scalar vw_sfc_s{shoc_input.vw_sfc(i)[0]};
    const Scalar phis_s{shoc_input.phis(i)[0]};
    Scalar pblh_s{0};

    const auto zt_grid_s      = ekat::subview(shoc_input.zt_grid, i);
    const auto zi_grid_s      = ekat::subview(shoc_input.zi_grid, i);
    const auto pres_s         = ekat::subview(shoc_input.pres, i);
    const auto presi_s        = ekat::subview(shoc_input.presi, i);
    const auto pdel_s         = ekat::subview(shoc_input.pdel, i);
    const auto thv_s          = ekat::subview(shoc_input.thv, i);
    const auto w_field_s      = ekat::subview(shoc_input.w_field, i);
    const auto wtracer_sfc_s  = ekat::subview(shoc_input.wtracer_sfc, i);
    const auto exner_s        = ekat::subview(shoc_input.exner, i);

    const auto rho_zt_s       = ekat::subview(rho_zt_d, i);
    const auto shoc_qv_s      = ekat::subview(shoc_qv_d, i);
    const auto dz_zt_s        = ekat::subview(dz_zt_d, i);
    const auto dz_zi_s        = ekat::subview(dz_zi_d, i);
    const auto thv_zi_s       = ekat::subview(thv_zi_d, i);
    const auto sterm_s        = ekat::subview(sterm_d, i);
    const auto sterm_zt_s     = ekat::subview(sterm_zt_d, i);
    const auto a_diss_s       = ekat::subview(a_diss_d, i);
    const auto rdp_zt_s       = ekat::subview(rdp_zt_d, i);
    const auto tmpi_s         = ekat::subview(tmpi_d, i);
    const auto tkh_zi_s       = ekat::subview(tkh_zi_d, i);
    const auto tk_zi_s        = ekat::subview(tk_zi_d, i);
    const auto rho_zi_s       = ekat::subview(rho_zi_d, i);
    const auto du_s           = ekat::subview(du_d, i);
    const auto dl_s           = ekat::subview(dl_d, i);
    const auto d_s            = ekat::subview(d_d, i);
    const auto isotropy_zi_s  = ekat::subview(isotropy_zi_d, i);
    const auto w_sec_zi_s     = ekat::subview(w_sec_zi_d, i);
    const auto brunt_zi_s     = ekat::subview(brunt_zi_d, i);
    const auto thetal_zi_s    = ekat::subview(thetal_zi_d, i);
    const auto wthl_sec_zt_s  = ekat::subview(wthl_sec_zt_d, i);
    const auto wqw_sec_zt_s   = ekat::subview(wqw_sec_zt_d, i);
    const auto w3_zt_s        = ekat::subview(w3_zt_d, i);
    const auto thl_sec_zt_s   = ekat::subview(thl_sec_zt_d, i);
    const auto qwthl_sec_zt_s = ekat::subview(qwthl_sec_zt_d, i);
    const auto qw_sec_zt_s    = ekat::subview(qw_sec_zt_d, i);
    const auto pblintd_thv_s  = ekat::subview(pblintd_thv_d, i);
    const auto pblintd_rino_s = ekat::subview(pblintd_rino_d, i);

    const auto host_dse_s     = ekat::subview(shoc_input_output.host_dse, i);
    const auto tke_s          = ekat::subview(shoc_input_output.tke, i);
    const auto thetal_s       = ekat::subview(shoc_input_output.thetal, i);
    const auto qw_s           = ekat::subview(shoc_input_output.qw, i);
    const auto u_wind_s       = ekat::subview(shoc_input_output.u_wind, i);
    const auto v_wind_s       = ekat::subview(shoc_input_output.v_wind, i);
    const auto wthv_sec_s     = ekat::subview(shoc_input_output.wthv_sec, i);
    const auto tk_s           = ekat::subview(shoc_input_output.tk, i);
    const auto tkh_s          = ekat::subview(shoc_input_output.tkh, i);
    const auto shoc_cldfrac_s = ekat::subview(shoc_input_output.shoc_cldfrac, i);
    const auto shoc_ql_s      = ekat::subview(shoc_input_output.shoc_ql, i);
    const auto shoc_ql2_s     = ekat::subview(shoc_output.shoc_ql2, i);
    const auto shoc_mix_s     = ekat::subview(shoc_history_output.shoc_mix, i);
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

    const auto X1_s       = Kokkos::subview(X1_d, i, Kokkos::ALL(), Kokkos::ALL());
    const auto qtracers_s = Kokkos::subview(shoc_input_output.qtracers, i, Kokkos::ALL(), Kokkos::ALL());

    shoc_main_internal(team, nlev, nlevi, npbl, nadv, num_qtracers, dtime,
                       host_dx_s, host_dy_s, zt_grid_s, zi_grid_s,
                       pres_s, presi_s, pdel_s, thv_s, w_field_s,
                       wthl_sfc_s, wqw_sfc_s, uw_sfc_s, vw_sfc_s,
                       wtracer_sfc_s, exner_s, phis_s, rho_zt_s,
                       shoc_qv_s, dz_zt_s, dz_zi_s, thv_zi_s, sterm_s,
                       sterm_zt_s, a_diss_s, rdp_zt_s, tmpi_s, tkh_zi_s,
                       tk_zi_s, rho_zi_s, du_s, dl_s, d_s, X1_s, isotropy_zi_s,
                       w_sec_zi_s, brunt_zi_s, thetal_zi_s, wthl_sec_zt_s,
                       wqw_sec_zt_s, w3_zt_s, thl_sec_zt_s, qwthl_sec_zt_s,
                       qw_sec_zt_s, pblintd_thv_s, pblintd_rino_s, host_dse_s, tke_s,
                       thetal_s, qw_s, u_wind_s, v_wind_s, wthv_sec_s, qtracers_s,
                       tk_s, tkh_s, shoc_cldfrac_s, shoc_ql_s, pblh_s, shoc_ql2_s,
                       shoc_mix_s, w_sec_s, thl_sec_s, qw_sec_s, qwthl_sec_s,
                       wthl_sec_s, wqw_sec_s, wtke_sec_s, uw_sec_s, vw_sec_s, w3_s,
                       wqls_sec_s, brunt_s, isotropy_s);

    shoc_output.pblh(i)[0] = pblh_s;
  });

  auto finish = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
  return duration.count();
}

} // namespace shoc
} // namespace scream

#endif
