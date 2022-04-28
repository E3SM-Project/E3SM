
#include "shoc.h"

// temporarily disable all this code for now
#ifdef use_shoc

using namespace scream;
using namespace scream::shoc;

// initialize shoc
void shoc_initialize() {
  YAKL_SCOPE( sgs_field      , ::sgs_field);
  YAKL_SCOPE( sgs_field_diag , ::sgs_field_diag);
  int  constexpr perturb_num_layers = 5;

  parallel_for( SimpleBounds<5>(nsgs_fields,nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int l, int k, int j, int i, int icrm) {
    if (k<=perturb_num_layers) {
      real perturb_k_scaling = ((real)perturb_num_layers-k) / (real)perturb_num_layers;
      sgs_field(l,k,j,i,icrm) = 0.01 * perturb_k_scaling + shoc_min_tke;
    } else {
      sgs_field(l,k,j,i,icrm) = shoc_min_tke;
    }
  });

  parallel_for( SimpleBounds<5>(nsgs_fields_diag,nzm,dimy_d,dimx_d,ncrms) , YAKL_LAMBDA (int l, int k, int j, int i, int icrm) {
    sgs_field_diag(l,k,j,i,icrm) = 0.0;
  });
}

void shoc_proc() {
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Spack      = typename SHOC::Spack;
  using view_1d    = typename SHOC::view_1d<Scalar>;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using view_3d    = typename SHOC::view_3d<Spack>;
  using ExeSpace   = typename SHOC::KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;
  YAKL_SCOPE( dtime          , :: dt);
  YAKL_SCOPE( ncrms          , :: ncrms);
  YAKL_SCOPE( dx             , :: dx);
  YAKL_SCOPE( dy             , :: dy);
  YAKL_SCOPE( dz             , :: dz);
  YAKL_SCOPE( u              , :: u);
  YAKL_SCOPE( v              , :: v);
  YAKL_SCOPE( w              , :: w);
  YAKL_SCOPE( fluxbu         , :: fluxbu);
  YAKL_SCOPE( fluxbv         , :: fluxbv);
  YAKL_SCOPE( gamaz          , :: gamaz);
  YAKL_SCOPE( t              , :: t);
  YAKL_SCOPE( tabs           , :: tabs);
  YAKL_SCOPE( qc             , :: qc);
  YAKL_SCOPE( qv             , :: qv);
  YAKL_SCOPE( pmid_in        , :: pres);
  YAKL_SCOPE( pint_in        , :: presi);
  YAKL_SCOPE( pdel_in        , :: pdel);
  YAKL_SCOPE( phis           , :: phis);
  YAKL_SCOPE( CF3D           , :: CF3D);
  YAKL_SCOPE( sgs_field      , :: sgs_field);
  YAKL_SCOPE( sgs_field_diag , :: sgs_field_diag);
  YAKL_SCOPE( z              , :: z);
  YAKL_SCOPE( zi             , :: zi);
  YAKL_SCOPE( qcl            , :: qcl);
  YAKL_SCOPE( qci            , :: qci);
  YAKL_SCOPE( qpl            , :: qpl);
  YAKL_SCOPE( qpi            , :: qpi); 
  YAKL_SCOPE( micro_field    , :: micro_field);

  const int nlev   = nzm;
  const int nlevi  = nzm+1;
  const int ncol   = ncrms*nx*ny;
  const int npack  = ekat::npack<Spack>(nlev);
  const int nipack = ekat::npack<Spack>(nlevi);

  // keep these constants here for references
  // const Real theta0        = 300.;   // Reference temperature                     [K]
  // const Real ts_nudge      = 86400.; // Time scale for u/v nudging (not used)     [s]
  // const Real p0_shoc       = 100000.;
  // const Real shoc_tk1      = 268.15;
  // const Real shoc_tk2      = 238.15;
  // const Real shoc_liq_deep = 8.e-6;
  // const Real shoc_liq_sh   = 10.e-6;
  // const Real shoc_ice_deep = 25.e-6;
  // const Real shoc_ice_sh   = 50.e-6;

  real2d inv_exner("inv_exner", ncol, nlev);
  real2d shoc_qw("shoc_qw", ncol, nlev); // total water (vapor + cloud liquid)
  real2d shoc_ql("shoc_ql", ncol, nlev); // cloud liquid water
  real2d thlm("thlm", ncol, nlev);
  real2d thv("thv", ncol, nlev);
  real2d dz_g("dz_g", ncol, nlev);
  real2d zt_g("zt_g", ncol, nlev);
  real2d zi_g("zi_g", ncol, nlevi);
  real2d rrho("rrho", ncol, nlev);
  real2d wm_zt("wm_zt", ncol, nlev);
  real2d tke("tke", ncol, nlev);
  real2d tk("tk", ncol, nlev);
  real2d tkh("tkh", ncol, nlev);
  real2d wthv("wthv", ncol, nlev);
  real2d pmid("pmid", ncol, nlev);
  real2d pdel("pdel", ncol, nlev);
  real2d pint("pint", ncol, nlevi);
  real2d shoc_dse("shoc_dse", ncol, nlev);
  real2d shoc_cldfrac("shoc_cldfrac", ncol, nlev);
  real3d shoc_hwind("shoc_hwind", ncol, 2, nlev);
  real3d qtracers("qtracers", ncol, num_shoc_tracers, nlev);

  // ------------------------------------------------- 
  // Set input state for SHOC
  // -------------------------------------------------
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = nzm-(k+1);

    pmid(icol,ilev)  = 100.*pmid_in(k,icrm);
    pdel(icol,ilev)  = 100.*pdel_in(k,icrm);

    inv_exner(icol,ilev) = 1./std::pow((pmid_in(k,icrm)*1.0e-3), (rgas/cp));

    dz_g(icol,ilev)   = dz(icrm);
    zt_g(icol,ilev)   = z(k,icrm);

    rrho(icol,ilev)   = (1./ggr)*(pdel(icol,ilev)/dz_g(icol,ilev));
    wm_zt(icol,ilev)  = w(k,j+offy_w,i+offx_w,icrm);

    shoc_qw(icol,ilev) = qv(k,j,i,icrm) + qcl(k,j,i,icrm);
    shoc_ql(icol,ilev) = qcl(k,j,i,icrm);

    shoc_hwind(icol,0,ilev) = u(k,j+offy_u,i+offx_u,icrm);
    shoc_hwind(icol,1,ilev) = v(k,j+offy_v,i+offx_v,icrm);

    tabs(k,j,i,icrm) = t(k,j+offy_s,i+offx_s,icrm) - gamaz(k,icrm)
                      + fac_cond *( qcl(k,j,i,icrm) + qpl(k,j,i,icrm) ) 
                      + fac_sub  *( qci(k,j,i,icrm) + qpi(k,j,i,icrm) );

    real theta_zt = tabs(k,j,i,icrm) * inv_exner(icol, ilev);
    thlm(icol,ilev) = theta_zt - (theta_zt/tabs(k,j,i,icrm)) * (latvap/cp)*qcl(k,j,i,icrm);
    thv(icol,ilev)  = theta_zt*(1 + zvir*qv(k,j,i,icrm) - qcl(k,j,i,icrm));
    shoc_dse(icol,ilev) = cp*tabs(k,j,i,icrm) + ggr*z(k,icrm) + phis(icrm);

    tke (icol,ilev) = max(shoc_min_tke, sgs_field(0,k,j+offy_s,i+offx_s,icrm) ); // enforce min TKE value - same as SCREAM
    tk  (icol,ilev) = sgs_field_diag(0,k,j+offy_d,i+offx_d,icrm);
    tkh (icol,ilev) = sgs_field_diag(1,k,j+offy_d,i+offx_d,icrm);
    wthv(icol,ilev) = sgs_field_diag(2,k,j+offy_d,i+offx_d,icrm);

    // Cloud fraction needs to be initialized for first PBL height calculation
    shoc_cldfrac(icol,ilev) = CF3D(k,j,i,icrm);

    qtracers(icol,shoc_idx_qv,ilev) = micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_nc,ilev) = micro_field(idx_nc,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_qr,ilev) = micro_field(idx_qr,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_nr,ilev) = micro_field(idx_nr,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_qi,ilev) = micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_qm,ilev) = micro_field(idx_qm,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_ni,ilev) = micro_field(idx_ni,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_bm,ilev) = micro_field(idx_bm,k,j+offy_s,i+offx_s,icrm);
    qtracers(icol,shoc_idx_qc,ilev) = micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm);
    // qtracers(icol,shoc_idx_tke,ilev)= sgs_field(0,k,offy_s+j,offx_s+i,icrm);
    qtracers(icol,shoc_idx_tke,ilev)= tke (icol,ilev);

    // if (use_ESMT) {
    //   qtracers(icol,shoc_idx_esmt_u,ilev) = u_esmt(k,j,i,icrm)
    //   qtracers(icol,shoc_idx_esmt_v,ilev) = v_esmt(k,j,i,icrm)
    // }
  });

  parallel_for( SimpleBounds<4>(nz, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = nz-(k+1);
    zi_g(icol,ilev) = zi(k,icrm);
    pint(icol,ilev) = 100.*pint_in(k,icrm);
  });

  // -------------------------------------------------
  // Set surface geopotential and fluxes
  // -------------------------------------------------
  view_1d host_dx_1d("host_dx", ncol),
          host_dy_1d("host_dy", ncol),
          wthl_sfc_1d("wthl_sfc", ncol),  // Surface sensible heat flux [K m/s]
          wqw_sfc_1d("wqw_sfc", ncol),    // Surface latent heat flux [kg/kg m/s]
          uw_sfc_1d("uw_sfc", ncol),      // Surface momentum flux (u-direction) [m2/s2]
          vw_sfc_1d("vw_sfc", ncol),      // Surface momentum flux (v-direction) [m2/s2]
          phis_1d("phis_1d", ncol);       // Host model surface geopotential height


  Kokkos::parallel_for("host grid", ncol, KOKKOS_LAMBDA (const int& i) {
    host_dx_1d(i) = dx;
    host_dy_1d(i) = dy;
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ny, nx, ncrms}), KOKKOS_LAMBDA(int j, int i, int icrm) {
    // Surface fluxes provided by host model
    int icol = i+nx*(j+ny*icrm);
    // wthl_sfc_1d(icol) = shf.myData[icrm]/(cp*rrho.myData[icol*rrho.dimension[1]+nlev-1]);
    // wqw_sfc_1d(icol)  = cflx.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev-1];
    // uw_sfc_1d(icol)   = wsx.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev-1];
    // vw_sfc_1d(icol)   = wsy.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev-1];

    // wthl_sfc_1d(icol) = fluxbt(j,i,icrm);
    // wqw_sfc_1d(icol)  = fluxbq(j,i,icrm);
    // uw_sfc_1d(icol)   = fluxbu(j,i,icrm);
    // vw_sfc_1d(icol)   = fluxbv(j,i,icrm);

    wthl_sfc_1d(icol) = 0.;
    wqw_sfc_1d(icol)  = 0.;
    uw_sfc_1d(icol)   = 0.;
    vw_sfc_1d(icol)   = 0.;

    phis_1d(icol)     = phis.myData[icrm];//*100.0;
  });

  view_2d zt_grid_2d("zt_grid", ncol, npack),                 // heights, for thermo grid [m]
          zi_grid_2d("zi_grid", ncol, nipack),                // heights, for interface grid [m]
          pmid_2d("pmid", ncol, npack),                       // pressure levels on thermo grid [Pa]
          pint_2d("pint", ncol, nipack),                      // pressure levels on interface grid [Pa]
          pdel_2d("pdel", ncol, npack),                       // Differences in pressure levels [Pa]
          thv_2d("thv", ncol, npack),                         // virtual potential temperature [K]
          w_field_2d("w_field", ncol, npack),                 // large scale vertical velocity [m/s]
          wtracer_sfc_2d("wtracer", ncol, num_shoc_tracers),  // Surface flux for tracers [varies]
          inv_exner_2d("inv_exner", ncol, npack);

  array_to_view(zt_g.myData,      ncol, nlev,  zt_grid_2d);
  array_to_view(zi_g.myData,      ncol, nlevi, zi_grid_2d);
  array_to_view(pmid.myData,      ncol, nlev,  pmid_2d);
  array_to_view(pint.myData,      ncol, nlevi, pint_2d);
  array_to_view(pdel.myData,      ncol, nlev,  pdel_2d);
  array_to_view(inv_exner.myData, ncol, nlev,  inv_exner_2d);
  array_to_view(thv.myData,       ncol, nlev,  thv_2d);
  array_to_view(wm_zt.myData,     ncol, nlev,  w_field_2d);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, num_shoc_tracers}), KOKKOS_LAMBDA(int icol, int itrc) {
    wtracer_sfc_2d(icol,itrc) = 0.;
  });

  SHOC::SHOCInput shoc_input{host_dx_1d, host_dy_1d, zt_grid_2d, zi_grid_2d,
                            pmid_2d, pint_2d, pdel_2d, thv_2d,
                            w_field_2d, wthl_sfc_1d, wqw_sfc_1d, uw_sfc_1d,
                            vw_sfc_1d, wtracer_sfc_2d, inv_exner_2d, phis_1d};

  view_2d host_dse_2d("host_dse", ncol, npack);                 // dry static energy [J/kg] : dse = Cp*T + g*z + phis
  view_2d tke_2d("tke", ncol, npack);                           // turbulent kinetic energy [m2/s2]
  view_2d thetal_2d("thetal", ncol, npack);                     // liquid water potential temperature [K]
  view_2d shoc_qw_2d("shoc_qw", ncol, npack);                   // total water mixing ratio [kg/kg]
  view_2d shoc_ql_2d("shoc_ql", ncol, npack);                   // Vector-valued wind (u,v) [m/s]
  view_2d wthv_sec_2d("wthv_sec", ncol, npack);                 // buoyancy flux [K m/s]
  view_2d tk_2d("tk", ncol, npack);                             // tracers [varies]
  view_2d tkh_2d("tkh", ncol, npack);                           // eddy coefficient for momentum [m2/s]
  view_2d shoc_cldfrac_2d("shoc_cldfrac", ncol, npack);         // eddy heat conductivity
  view_3d shoc_hwind_3d("shoc_hwind",ncol,2,npack);             // Cloud fraction [-]
  view_3d qtracers_3d("qtracers",ncol,num_shoc_tracers,npack);  // cloud liquid mixing ratio [kg/kg]

  array_to_view(shoc_dse.myData,     ncol, nlev, host_dse_2d);
  array_to_view(tke.myData,          ncol, nlev, tke_2d);
  array_to_view(thlm.myData,         ncol, nlev, thetal_2d);
  array_to_view(shoc_qw.myData,      ncol, nlev, shoc_qw_2d);
  array_to_view(shoc_ql.myData,      ncol, nlev, shoc_ql_2d);
  array_to_view(tk.myData,           ncol, nlev, tk_2d);
  array_to_view(tkh.myData,          ncol, nlev, tkh_2d);
  array_to_view(wthv.myData,         ncol, nlev, wthv_sec_2d);
  array_to_view(shoc_cldfrac.myData, ncol, nlev, shoc_cldfrac_2d);
  array_to_view(shoc_hwind.myData, ncol, 2, nlev, shoc_hwind_3d);
  array_to_view(qtracers.myData, ncol, num_shoc_tracers, nlev, qtracers_3d);

  SHOC::SHOCInputOutput shoc_input_output{host_dse_2d, tke_2d, thetal_2d, shoc_qw_2d,
                                         shoc_hwind_3d, wthv_sec_2d, qtracers_3d,
                                         tk_2d, tkh_2d, shoc_cldfrac_2d, shoc_ql_2d};

  view_1d pblh_1d("pblh",ncol);
  view_2d shoc_ql2_2d("shoc_ql2",ncol,npack);
  SHOC::SHOCOutput shoc_output{pblh_1d, shoc_ql2_2d};

  view_2d shoc_mix_2d("shoc_mix", ncol, npack),     // Turbulent length scale [m]
          w_sec_2d("w_sec", ncol, npack),           // vertical velocity variance [m2/s2]
          thl_sec_2d("thl_sec", ncol, nipack),      // temperature variance [K^2]
          qw_sec_2d("qw_sec", ncol, nipack),        // moisture variance [kg2/kg2]
          qwthl_sec_2d("qwthl_sec", ncol, nipack),  // temp moisture covariance [K kg/kg]
          wthl_sec_2d("wthl_sec", ncol, nipack),    // vertical heat flux [K m/s]
          wqw_sec_2d("wqw_sec", ncol, nipack),      // vertical moisture flux [K m/s]
          wtke_sec_2d("wtke_sec", ncol, nipack),    // vertical tke flux [m3/s3]
          uw_sec_2d("uw_sec", ncol, nipack),        // vertical zonal momentum flux [m2/s2]
          vw_sec_2d("vw_sec", ncol, nipack),        // vertical meridional momentum flux [m2/s2]
          w3_2d("w3", ncol, nipack),                // third moment vertical velocity [m3/s3]
          wqls_sec_2d("wqls_sec", ncol, npack),     // liquid water flux [kg/kg m/s]
          brunt_2d("brunt", ncol, npack),           // brunt vaisala frequency [s-1]
          isotropy_2d("isotropy", ncol, npack);     // return to isotropic timescale [s]

  SHOC::SHOCHistoryOutput shoc_history_output{shoc_mix_2d, w_sec_2d, thl_sec_2d, qw_sec_2d,
                                              qwthl_sec_2d, wthl_sec_2d, wqw_sec_2d, wtke_sec_2d,
                                              uw_sec_2d, vw_sec_2d, w3_2d, wqls_sec_2d, brunt_2d, isotropy_2d};

  const int nwind = ekat::npack<Spack>(2)*Spack::n;
  const int ntrac = ekat::npack<Spack>(num_shoc_tracers+3)*Spack::n;
  const auto policy = ekat::ExeSpaceUtils<SHOC::KT::ExeSpace>::get_default_team_policy(ncol, npack);
  ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 128+(nwind+ntrac), policy);

  const auto elapsed_microsec = SHOC::shoc_main(ncol, nlev, nlevi, nlev, 1, num_shoc_tracers, dtime, workspace_mgr,
                                                shoc_input, shoc_input_output, shoc_output, shoc_history_output);

  // get SHOC output back to CRM 
  view_to_array(shoc_input_output.tk,   ncol, nlev, tk);
  view_to_array(shoc_input_output.tkh,  ncol, nlev, tkh);
  view_to_array(shoc_input_output.wthv_sec, ncol, nlev, wthv);
  // view_to_array(shoc_input_output.tke, ncol, nlev, tke);
  view_to_array(shoc_input_output.shoc_cldfrac, ncol, nlev, shoc_cldfrac);
  view_to_array(shoc_input_output.horiz_wind, ncol, 2, nlev, shoc_hwind);
  view_to_array(shoc_input_output.host_dse, ncol, nlev, shoc_dse);
  view_to_array(shoc_input_output.qtracers, ncol, num_shoc_tracers, nlev, qtracers);


  // update tracers and TKE
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = nzm-(k+1);
    micro_field(idx_qv,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_qv,ilev);
    micro_field(idx_nc,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_nc,ilev);
    micro_field(idx_qr,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_qr,ilev);
    micro_field(idx_nr,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_nr,ilev);
    micro_field(idx_qi,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_qi,ilev);
    micro_field(idx_qm,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_qm,ilev);
    micro_field(idx_ni,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_ni,ilev);
    micro_field(idx_bm,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_bm,ilev);
    micro_field(idx_qc,k,j+offy_s,i+offx_s,icrm) = qtracers(icol,shoc_idx_qc,ilev);
    sgs_field(0,k,offy_s+j,offx_s+i,icrm)        = qtracers(icol,shoc_idx_tke,ilev);
  });

  // update diagnostic micro fields based on micro scheme
  if (is_same_str(microphysics_scheme, "sam1mom") == 0) { micro_diagnose(); }
  if (is_same_str(microphysics_scheme, "p3")      == 0) { micro_p3_diagnose(); }

  // update other variables
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = nzm-(k+1);
    tabs(k,j,i,icrm) = ( shoc_dse(icol,ilev) - ggr*z(k,icrm) - phis(icrm) )/cp;
    t(k,j+offy_s,i+offx_s,icrm) = tabs(k,j,i,icrm) + gamaz(k,icrm)
                  - fac_cond *( qcl(k,j,i,icrm) - qpl(k,j,i,icrm) )
                  - fac_sub  *( qci(k,j,i,icrm) - qpi(k,j,i,icrm) );
    u(k,j+offy_u,i+offx_u,icrm) = shoc_hwind(icol,0,ilev);
    v(k,j+offy_v,i+offx_v,icrm) = shoc_hwind(icol,1,ilev);
    // sgs_field(0,k,offy_s+j,offx_s+i,icrm)      = tke(icol,nlev-(ilev+1));
    sgs_field_diag(0,k,offy_d+j,offx_d+i,icrm) = tk(icol,ilev);
    sgs_field_diag(1,k,offy_d+j,offx_d+i,icrm) = tkh(icol,ilev);
    sgs_field_diag(2,k,offy_d+j,offx_d+i,icrm) = wthv(icol,ilev);
    CF3D(k,j,i,icrm) = shoc_cldfrac(icol,ilev);

    // if (use_ESMT) {
    //   u_esmt(k,j,i,icrm) = qtracers(icol,shoc_idx_esmt_u,ilev)
    //   v_esmt(k,j,i,icrm) = qtracers(icol,shoc_idx_esmt_v,ilev)
    // }

  });

}

#endif
