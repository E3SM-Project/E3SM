
#include "shoc.h"

using namespace scream;
using namespace scream::shoc;

void shoc_proc() {
  using SHOC       = Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Spack      = typename SHOC::Spack;
  using view_1d    = typename SHOC::view_1d<Scalar>;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using view_3d    = typename SHOC::view_3d<Spack>;
  using ExeSpace   = typename SHOC::KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;
  auto &u                  = :: u;
  auto &v                  = :: v;
  auto &w                  = :: w;
  auto &rho                = :: rho;
  auto &dtime              = :: dt;
  auto &tabs               = :: tabs;
  auto &ncrms              = :: ncrms;
  auto &dx                 = :: dx;
  auto &dy                 = :: dy;
  auto &dz                 = :: dz;
  auto &z0                 = :: z0;
  auto &zi                 = :: zi;
  auto &zm                 = :: zm;
  auto &z                  = :: z;
  auto &sl                 = :: sl;
  auto &omega              = :: omega;
  auto &crm_state_qv       = :: crm_state_qv; // ixq
  auto &crm_state_qc       = :: crm_state_qc; // ixcldliq
  auto &crm_state_qi       = :: crm_state_qi; // ixcldice
  auto &crm_state_nc       = :: crm_state_ni; // ixnumliq
  auto &crm_state_ni       = :: crm_state_ni; // ixnumice
  auto &phis               = :: phis;
  auto &shf                = :: crm_input_shf;
  auto &cflx               = :: crm_input_cflx;
  auto &wsx                = :: crm_input_wsx;
  auto &wsy                = :: crm_input_wsy;
  auto &crm_input_tke      = :: crm_input_tke_zt;
  auto &crm_input_wthv     = :: crm_input_wthv;
  auto &crm_input_tkh      = :: crm_input_tkh;
  auto &crm_input_tk       = :: crm_input_tk;
  auto &crm_input_alst     = :: crm_input_alst;
  auto &crm_input_qtracers = :: crm_input_qtracers;
  auto &crm_input_tl       = :: crm_input_tl;
  auto &crm_input_zint     = :: crm_input_zint;
  auto &crm_input_zmid     = :: crm_input_zmid;
  auto &crm_input_pmid     = :: crm_input_pmid;
  auto &crm_input_pint     = :: crm_input_pint;
  auto &crm_input_pdel     = :: crm_input_pdel;
  auto &crm_input_phis     = :: crm_input_phis;
  auto &crm_input_sl       = :: crm_input_sl;
  auto &crm_input_omega    = :: crm_input_omega;
  auto &crm_input_qccl     = :: crm_input_qccl;
  auto &crm_input_ql       = :: crm_input_ql;
  auto &crm_input_ul       = :: crm_input_ul;
  auto &crm_input_vl       = :: crm_input_vl;
  auto &sgs_field          = :: sgs_field;
  auto &sgs_field_diag     = :: sgs_field_diag;
 
  const int nlev   = nzm;
  const int nlevi  = nzm+1;
  const int ncol   = ncrms*nx*ny;
  const int npack  = ekat::npack<Spack>(nlev);
  const int nipack = ekat::npack<Spack>(nlevi);

  const Real theta0        = 300.;   // Reference temperature                     [K]
  const Real ts_nudge      = 86400.; // Time scale for u/v nudging (not used)     [s]
  const Real p0_shoc       = 100000.;
  const Real shoc_tk1      = 268.15;
  const Real shoc_tk2      = 238.15;
  const Real shoc_liq_deep = 8.e-6;
  const Real shoc_liq_sh   = 10.e-6;
  const Real shoc_ice_deep = 25.e-6;
  const Real shoc_ice_sh   = 50.e-6;

  real2d inv_exner("inv_exner", ncol, nlev);
  real2d rvm("rvm", ncol, nlev);
  real2d rcm("rcm", ncol, nlev);
  real2d rtm("rtm", ncol, nlev);
  real2d um("um", ncol, nlev);
  real2d vm("vm", ncol, nlev);
  real2d thlm("thlm", ncol, nlev);
  real2d thv("thv", ncol, nlev);
  real2d cloud_frac("cloud_frac", ncol, nlev);
  real2d dz_g("dz_g", ncol, nlev);
  real2d zt_g("zt_g", ncol, nlev);
  real2d zi_g("zi_g", ncol, nlevi);
  real2d rrho("rrho", ncol, nlev);
  real2d wm_zt("wm_zt", ncol, nlev);
  real2d shoc_s("shoc_s", ncol, nlev);
  real2d tke("tke", ncol, nlev);
  real2d tk("tk", ncol, nlev);
  real2d tkh("tkh", ncol, nlev);
  real2d wthv("wthv", ncol, nlev);
  real2d alst("alst", ncol, nlev);
  real2d qv("qv", ncol, nlev);
  real2d qc("qc", ncol, nlev);
  real2d qi("qi", ncol, nlev);
  real2d nc("nc", ncol, nlev);
  real2d ni("ni", ncol, nlev);
  real2d pres("pres", ncol, nlev);
  real2d pdels("pdel", ncol, nlev);
  real2d presi("presi", ncol, nlevi);
  real3d qtracers("qtracers", ncol, 10, nlev);

  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    inv_exner(icol, ilev) = 1./std::pow((crm_input_pmid(plev-nzm+k,icrm)/p0_shoc), (rgas/cp));
  });

  // At each SHOC call, initialize mean momentum  and thermo SHOC state 
  //  from the E3SM state
  //do k=1,pver   ! loop over levels
  //  do i=1,ncol ! loop over columns
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;

    rvm(icol,ilev) = crm_input_ql(plev-nzm+k,icrm);
    rcm(icol,ilev) = crm_input_qccl(plev-nzm+k,icrm);
    rtm(icol,ilev) = rvm(icol,ilev) + rcm(icol,ilev);
    um(icol,ilev)  = crm_input_ul(plev-nzm+k,icrm);
    vm(icol,ilev)  = crm_input_vl(plev-nzm+k,icrm);

    thlm(icol,ilev)  = crm_input_tl(plev-nzm+k,icrm)*inv_exner(icol,ilev)-(latvap/cp)*crm_state_qc(k,j,i,icrm);
    thv(icol,ilev)   = crm_input_tl(plev-nzm+k,icrm)*inv_exner(icol,ilev)
                       *(1.0+zvir*crm_state_qv(k,j,i,icrm)-crm_state_qc(k,j,i,icrm));

    tke(icol,ilev)   = crm_input_tke(plev-nzm+k,icrm);
    tkh(icol,ilev)   = crm_input_tkh(plev-nzm+k,icrm);
    tk (icol,ilev)   = crm_input_tk(plev-nzm+k,icrm);
    wthv(icol,ilev)  = crm_input_wthv(plev-nzm+k,icrm);
    alst(icol,ilev)  = crm_input_alst(plev-nzm+k,icrm);
    pres(icol,ilev)  = crm_input_pmid(plev-nzm+k,icrm);
    pdels(icol,ilev) = crm_input_pdel(plev-nzm+k,icrm);

    // Cloud fraction needs to be initialized for first 
    // PBL height calculation call
    cloud_frac(icol,ilev) = crm_input_alst(plev-nzm+k,icrm);

    for (auto q=0; q<10; ++q)  {
       qtracers(icol,q,ilev) = crm_input_qtracers(plev-nzm+k+1,icrm,q);
    }
  });

  parallel_for( SimpleBounds<4>(nz, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    presi(icol,ilev) = crm_input_pint(plev+1-nz+k,icrm);
  });

  // ------------------------------------------------- 
  // Prepare inputs for SHOC call                      
  //------------------------------------------------- 
  //do k=1,pver
  //  do i=1,ncol
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
//    dz_g(icol,ilev) = zi(k+1,icrm)-zi(k,icrm);   //compute thickness
      dz_g(icol,ilev) = crm_input_zint(plev-nzm+k,icrm)-crm_input_zint(plev-nzm+k+1,icrm);
  });


  //  Define the SHOC thermodynamic grid (in units of m)
//   wm_zt(:,pver) = 0._r8
//   do k=1,pver
//     do i=1,ncol
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
//    zt_g(icol,ilev)   = z(k,icrm);
    zt_g(icol,ilev)   = crm_input_zmid(plev-nzm+k,icrm)-crm_input_zint(plev+1,icrm);
    rrho(icol,ilev)   = (1./ggr)*(pdels(icol,ilev)/dz_g(icol,ilev));
    wm_zt(icol,ilev)  = -1.*crm_input_omega(plev-nzm+k,i)/(rrho(icol,ilev)*ggr);
    shoc_s(icol,ilev) = crm_input_sl(plev-nzm+k,icrm);
  });

//   do k=1,pverp
//     do i=1,ncol
  parallel_for( SimpleBounds<4>(nz, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
//    zi_g(icol,ilev) = zi(k,icrm);
    zi_g(icol,ilev) = crm_input_zint(plev+1-nz+k,icrm)-crm_input_zint(plev+1,icrm);
  });

  view_1d host_dx_1d("host_dx", ncol),
          host_dy_1d("host_dy", ncol),
          wthl_sfc_1d("wthl_sfc", ncol), // wpthlp_sfc(1:ncol)
          wqw_sfc_1d("wqw_sfc", ncol),   // wprtp_sfc(1:ncol)
          uw_sfc_1d("uw_sfc", ncol),     // upwp_sfc(1:ncol)
          vw_sfc_1d("vw_sfc", ncol),     // vpwp_sfc(1:ncol)
          phis_1d("phis_1d", ncol);

  Kokkos::parallel_for("host grid", ncol, KOKKOS_LAMBDA (const int& i) {
    host_dx_1d(i) = dx;
    host_dy_1d(i) = dy;
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ny, nx, ncrms}), KOKKOS_LAMBDA(int j, int i, int icrm) {
    // Surface fluxes provided by host model
    int icol = i+nx*(j+ny*icrm);
    wthl_sfc_1d(icol) = shf.myData[icrm]/(cp*rrho.myData[icol*rrho.dimension[1]+nlev]);
    wqw_sfc_1d(icol)  = cflx.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev];
    uw_sfc_1d(icol)   = wsx.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev];
    vw_sfc_1d(icol)   = wsy.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev];
    phis_1d(icol)     = phis.myData[icrm];
  });

  view_2d zt_grid_2d("zt_grid", ncol, npack),
          zi_grid_2d("zi_grid", ncol, nipack),
          pres_2d("pres", ncol, npack),
          presi_2d("presi", ncol, nipack),
          pdel_2d("pdel", ncol, npack),
          thv_2d("thv", ncol, npack),
          w_field_2d("w_field", ncol, npack),  // wm_zt
          wtracer_sfc_2d("wtracer", ncol, 10),
          inv_exner_2d("inv_exner", ncol, npack);

  array_to_view(zt_g.myData, ncol, nlev, zt_grid_2d);
  array_to_view(zi_g.myData, ncol, nlevi, zi_grid_2d);
  array_to_view(pres.myData, ncol, nlev, pres_2d);
  array_to_view(presi.myData, ncol, nlevi, presi_2d);
  array_to_view(pdels.myData, ncol, nlev, pdel_2d);
  array_to_view(inv_exner.myData, ncol, nlev, inv_exner_2d);
  array_to_view(thv.myData, ncol, nlev, thv_2d);
  array_to_view(wm_zt.myData, ncol, nlev, w_field_2d);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, 10}), KOKKOS_LAMBDA(int k, int j) {
    wtracer_sfc_2d(j,k) = 0.;
  });

  SHOC::SHOCInput shoc_input{host_dx_1d, host_dy_1d, zt_grid_2d, zi_grid_2d,
                             pres_2d, presi_2d, pdel_2d, thv_2d,
                             w_field_2d, wthl_sfc_1d, wqw_sfc_1d, uw_sfc_1d,
                             vw_sfc_1d, wtracer_sfc_2d, inv_exner_2d, phis_1d};

  view_2d host_dse_2d("host_dse", ncol, npack), // shoc_s
          tke_2d("tke", ncol, npack),
          thetal_2d("thetal", ncol, npack),  // thlm
          qw_2d("qw", ncol, npack),  // rtm
          wthv_sec_2d("wthv_sec", ncol, npack),
          tk_2d("tk", ncol, npack), 
          shoc_cldfrac_2d("shoc_cldfrac", ncol, npack), //cloud_frac
          shoc_ql_2d("shoc_ql", ncol, npack), //rcm
          u_wind_2d("u_wind", ncol, npack),
          v_wind_2d("v_wind", ncol, npack);

  array_to_view(shoc_s.myData, ncol, nlev, host_dse_2d);
  array_to_view(tke.myData, ncol, nlev, tke_2d);
  array_to_view(thlm.myData, ncol, nlev, thetal_2d);
  array_to_view(rtm.myData, ncol, nlev, qw_2d);
  array_to_view(tk.myData, ncol, nlev, tk_2d);
  array_to_view(wthv.myData, ncol, nlev, wthv_sec_2d);
  array_to_view(cloud_frac.myData, ncol, nlev, shoc_cldfrac_2d);
  array_to_view(rcm.myData, ncol, nlev, shoc_ql_2d);
  array_to_view(um.myData, ncol, nlev, u_wind_2d);
  array_to_view(vm.myData, ncol, nlev, v_wind_2d);

  int num_qtracers {10};
  view_3d horiz_wind_3d("horiz_wind",ncol,2,npack);
  view_3d qtracers_3d("qtracers",ncol,num_qtracers,npack);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, npack}), KOKKOS_LAMBDA(int k, int j) {
     horiz_wind_3d(k,0,j) = u_wind_2d(k,j);
     horiz_wind_3d(k,1,j) = v_wind_2d(k,j);
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, num_qtracers, npack}), KOKKOS_LAMBDA(int k, int j, int i) {
    qtracers_3d(k,j,i)[0] = qtracers.myData[(k*qtracers.dimension[1]+j)*qtracers.dimension[2]+i];   
  });

  SHOC::SHOCInputOutput shoc_input_output{host_dse_2d, tke_2d, thetal_2d, qw_2d,
                                         horiz_wind_3d, wthv_sec_2d, qtracers_3d,
                                         tk_2d, shoc_cldfrac_2d, shoc_ql_2d};

  view_1d pblh_1d("pblh",ncol);
  view_2d shoc_ql2_2d("shoc_ql2",ncol,npack);
  SHOC::SHOCOutput shoc_output{pblh_1d, shoc_ql2_2d};

  view_2d shoc_mix_2d("shoc_mix", ncol, npack),
          w_sec_2d("w_sec", ncol, npack),
          thl_sec_2d("thl_sec", ncol, nipack),
          qw_sec_2d("qw_sec", ncol, nipack),
          qwthl_sec_2d("qwthl_sec", ncol, nipack),
          wthl_sec_2d("wthl_sec", ncol, nipack),
          wqw_sec_2d("wqw_sec", ncol, nipack),
          wtke_sec_2d("wtke_sec", ncol, nipack),
          uw_sec_2d("uw_sec", ncol, nipack),
          vw_sec_2d("vw_sec", ncol, nipack),
          w3_2d("w3", ncol, nipack),
          wqls_sec_2d("wqls_sec", ncol, npack),
          brunt_2d("brunt", ncol, npack),
          isotropy_2d("isotropy", ncol, npack);

  SHOC::SHOCHistoryOutput shoc_history_output{shoc_mix_2d, w_sec_2d, thl_sec_2d, qw_sec_2d,
                                              qwthl_sec_2d, wthl_sec_2d, wqw_sec_2d, wtke_sec_2d,
                                              uw_sec_2d, vw_sec_2d, w3_2d, wqls_sec_2d, brunt_2d, isotropy_2d};

  const int nwind = ekat::npack<Spack>(2)*Spack::n;
  const int ntrac = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
  const auto policy = ekat::ExeSpaceUtils<SHOC::KT::ExeSpace>::get_default_team_policy(ncol, npack);
  ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 23+(nwind+ntrac), policy);

  const auto elapsed_microsec = SHOC::shoc_main(ncol, nlev, nlevi, nlev, 1, num_qtracers, dtime, workspace_mgr,
                                               shoc_input, shoc_input_output, shoc_output, shoc_history_output);

  // get SHOC output back to CRM 
  view_to_array(shoc_input_output.tke, ncol, nlev, tke);
  view_to_array(shoc_input_output.tk, ncol, nlev, tk);

  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    sgs_field(1,k,offy_s+j,offx_s+i,icrm)      = tke(icol,nlev-(ilev+1));
    sgs_field(2,k,offy_s+j,offx_s+i,icrm)      = tk(icol,nlev-(ilev+1));
    sgs_field_diag(1,k,offy_d+j,offx_d+i,icrm) = tke(icol,nlev-(ilev+1));
    sgs_field_diag(2,k,offy_d+j,offx_d+i,icrm) = tk(icol,nlev-(ilev+1));
  });
}
