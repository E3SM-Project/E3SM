
#include "shoc.h"

using namespace scream;
using namespace scream::shoc;

// initialize shoc
void shoc_initialize() {
  auto &sgs_field      = ::sgs_field;
  auto &sgs_field_diag = ::sgs_field_diag;

  parallel_for( SimpleBounds<5>(nsgs_fields,nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int l, int k, int j, int i, int icrm) {
    sgs_field(l,k,j,i,icrm) = 0.0;
  });

  parallel_for( SimpleBounds<5>(nsgs_fields_diag,nzm,dimy_d,dimx_d,ncrms) , YAKL_LAMBDA (int l, int k, int j, int i, int icrm) {
    sgs_field_diag(l,k,j,i,icrm) = 0.0;
  });
}

// update precipitation
template <typename ArrayT>
void shoc_update_precipitation(const ArrayT& qtracers) {
    auto &qc   = :: qc;
    auto &qv   = :: qv;
    auto &qcl  = :: qcl;
    auto &qpl  = :: qpl;
    auto &qci  = :: qci;
    auto &qpi  = :: qpi;

    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int icol = i+nx*(j+ny*icrm);
      int ilev = k;
      int m = nzm-(k+1);

      micro_field(idx_qt,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_qt,ilev) + qtracers(icol,idx_qc,ilev);
      micro_field(idx_nc,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_nc,ilev);
      micro_field(idx_qr,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_qr,ilev);
      micro_field(idx_nr,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_nr,ilev);
      micro_field(idx_qi,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_qi,ilev);
      micro_field(idx_qm,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_qm,ilev);
      micro_field(idx_ni,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_ni,ilev);
      micro_field(idx_bm,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_bm,ilev);
      micro_field(idx_qc,m,j+offy_s,i+offx_s,icrm) = qtracers(icol,idx_qc,ilev);
      sgs_field(0,m,offy_s+j,offx_s+i,icrm)        = qtracers(icol,idx_tke,ilev);

      qc(m,j,i,icrm)  = qtracers(icol,idx_qc,ilev);
      qv(m,j,i,icrm)  = qtracers(icol,idx_qt,ilev);
      qcl(m,j,i,icrm) = qc(m,j,i,icrm);
      qpl(m,j,i,icrm) = qtracers(icol,idx_qr,ilev);
      qci(m,j,i,icrm) = qtracers(icol,idx_qi,ilev);
      qpi(m,j,i,icrm) = 0.;
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
  auto &dtime      = :: dt;
  auto &ncrms      = :: ncrms;
  auto &dx         = :: dx;
  auto &dy         = :: dy;
  auto &dz         = :: dz;
  auto &u          = :: u;
  auto &v          = :: v;
  auto &w          = :: w;
  auto &fluxbu     = :: fluxbu;
  auto &fluxbv     = :: fluxbv;
  auto &gamaz      = :: gamaz;
  auto &t          = :: t;
  auto &tabs       = :: tabs;
  auto &qc         = :: qc;
  auto &qv         = :: qv;
  auto &presin     = :: pres;
  auto &pdelin     = :: pdel;
  auto &pint       = :: presi;
  auto &prespot    = :: prespot;
  auto &phis       = :: phis;
  auto &CF3D       = :: CF3D;
  auto &sgs_field  = :: sgs_field;
  auto &sgs_field_diag = :: sgs_field_diag;
 
  const int nlev   = nzm;
  const int nlevi  = nzm+1;
  const int ncol   = ncrms*nx*ny;
  const int npack  = ekat::npack<Spack>(nlev);
  const int nipack = ekat::npack<Spack>(nlevi);

  // keep these constants here for references
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
  real2d pres("pres", ncol, nlev);
  real2d pdels("pdel", ncol, nlev);
  real2d presi("presi", ncol, nlevi);
  real3d qtracers("qtracers", ncol, num_tracers, nlev);

  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    inv_exner(icol, ilev) = 1./std::pow((pdel(nzm-(k+1),icrm)/1000.), (rgas/cp));
  });

  // At each SHOC call, initialize mean momentum  and thermo SHOC state 
  //  from the E3SM state
  //do k=1,pver   ! loop over levels
  //  do i=1,ncol ! loop over columns
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;

    rvm(icol,ilev) = qv(nzm-(k+1),j,i,icrm);
    rcm(icol,ilev) = qc(nzm-(k+1),j,i,icrm);
    rtm(icol,ilev) = rvm(icol,ilev) + rcm(icol,ilev);

    um(icol,ilev)  = u(nzm-(k+1),j+offy_u,i+offx_u,icrm);
    vm(icol,ilev)  = v(nzm-(k+1),j+offy_v,i+offx_v,icrm);

    thlm(icol,ilev)  = tabs(nzm-(k+1),j,i,icrm)*inv_exner(icol,ilev)-(latvap/cp)*qc(nzm-(k+1),j,i,icrm);
    thv(icol,ilev)   = tabs(nzm-(k+1),j,i,icrm)*inv_exner(icol,ilev)
                       *(1.0+zvir*rtm(icol,ilev)-rcm(icol,ilev));

    tke(icol,ilev)   = sgs_field(0,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    tkh(icol,ilev)   = sgs_field_diag(1,nzm-(k+1),j+offy_s,i+offx_s,icrm); 
    tk (icol,ilev)   = sgs_field_diag(0,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    wthv(icol,ilev)  = 0.; //crm_input_wthv(plev-nzm+k,icrm);
    pres(icol,ilev)  = 100.*presin(nzm-(k+1),icrm); 
    pdels(icol,ilev) = 100.*pdelin(nzm-(k+1),icrm);

    // Cloud fraction needs to be initialized for first 
    // PBL height calculation call
    cloud_frac(icol,ilev) = CF3D(nzm-(k+1),j,i,icrm);

    qtracers(icol,idx_qc,ilev) = micro_field(idx_qc,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_qt,ilev) = micro_field(idx_qt,nzm-(k+1),j+offy_s,i+offx_s,icrm) - qtracers(icol,idx_qc,ilev);
    qtracers(icol,idx_nc,ilev) = micro_field(idx_nc,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_qr,ilev) = micro_field(idx_qr,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_nr,ilev) = micro_field(idx_nr,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_qi,ilev) = micro_field(idx_qi,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_qm,ilev) = micro_field(idx_qm,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_ni,ilev) = micro_field(idx_ni,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_bm,ilev) = micro_field(idx_bm,nzm-(k+1),j+offy_s,i+offx_s,icrm);
    qtracers(icol,idx_tke,ilev)= sgs_field(0,nzm-(k+1),offy_s+j,offx_s+i,icrm);
  });

  parallel_for( SimpleBounds<4>(nz, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    presi(icol,ilev) = 100.*pint(nz-(k+1),icrm);
  });

  // ------------------------------------------------- 
  // Prepare inputs for SHOC call                      
  //------------------------------------------------- 
  //do k=1,pver
  //  do i=1,ncol
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
   //   dz_g(icol,ilev) = zi(k+1,icrm)-zi(k,icrm);   //compute thickness
      dz_g(icol,ilev) = dz(icrm);
  });


  //  Define the SHOC thermodynamic grid (in units of m)
//   wm_zt(:,pver) = 0._r8
//   do k=1,pver
//     do i=1,ncol
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    zt_g(icol,ilev)   = z(nzm-(k+1),icrm);
    rrho(icol,ilev)   = (1./ggr)*(pdels(icol,ilev)/dz_g(icol,ilev));
    wm_zt(icol,ilev)  = w(nzm-(k+1),j+offy_w,i+offx_w,icrm);
    shoc_s(icol,ilev) = cp*thlm(icol,ilev)+ggr*zt_g(icol,ilev)+phis(icrm);
  });

//   do k=1,pverp
//     do i=1,ncol
  parallel_for( SimpleBounds<4>(nz, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    zi_g(icol,ilev) = zi(nz-(k+1),icrm);
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
    wthl_sfc_1d(icol) = 0.; //shf.myData[icrm]/(cp*rrho.myData[icol*rrho.dimension[1]+nlev-1]);
    wqw_sfc_1d(icol)  = 0.; //cflx.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev-1];
    uw_sfc_1d(icol)   = 0.; //fluxbu(j,i,icrm); //crm_input_ul.myData[icrm]; //wsx.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev-1];
    vw_sfc_1d(icol)   = 0.; //fluxbu(j,i,icrm); //crm_input_vl.myData[icrm]; //wsy.myData[icrm]/rrho.myData[icol*rrho.dimension[1]+nlev-1];
    phis_1d(icol)     = 100.*phis.myData[icrm];
  });

  view_2d zt_grid_2d("zt_grid", ncol, npack),
          zi_grid_2d("zi_grid", ncol, nipack),
          pres_2d("pres", ncol, npack),
          presi_2d("presi", ncol, nipack),
          pdel_2d("pdel", ncol, npack),
          thv_2d("thv", ncol, npack),
          w_field_2d("w_field", ncol, npack),  // wm_zt
          wtracer_sfc_2d("wtracer", ncol, num_tracers),
          inv_exner_2d("inv_exner", ncol, npack);

  array_to_view(zt_g.myData, ncol, nlev, zt_grid_2d);
  array_to_view(zi_g.myData, ncol, nlevi, zi_grid_2d);
  array_to_view(pres.myData, ncol, nlev, pres_2d);
  array_to_view(presi.myData, ncol, nlevi, presi_2d);
  array_to_view(pdels.myData, ncol, nlev, pdel_2d);
  array_to_view(inv_exner.myData, ncol, nlev, inv_exner_2d);
  array_to_view(thv.myData, ncol, nlev, thv_2d);
  array_to_view(wm_zt.myData, ncol, nlev, w_field_2d);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, num_tracers}), KOKKOS_LAMBDA(int k, int j) {
    wtracer_sfc_2d(k,j) = 0.;
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
          tkh_2d("tkh", ncol, npack), 
          shoc_cldfrac_2d("shoc_cldfrac", ncol, npack), //cloud_frac
          shoc_ql_2d("shoc_ql", ncol, npack), //rcm
          u_wind_2d("u_wind", ncol, npack),
          v_wind_2d("v_wind", ncol, npack);

  array_to_view(shoc_s.myData, ncol, nlev, host_dse_2d);
  array_to_view(tke.myData, ncol, nlev, tke_2d);
  array_to_view(thlm.myData, ncol, nlev, thetal_2d);
  array_to_view(rtm.myData, ncol, nlev, qw_2d);
  array_to_view(tk.myData, ncol, nlev, tk_2d);
  array_to_view(tkh.myData, ncol, nlev, tkh_2d);
  array_to_view(wthv.myData, ncol, nlev, wthv_sec_2d);
  array_to_view(cloud_frac.myData, ncol, nlev, shoc_cldfrac_2d);
  array_to_view(rcm.myData, ncol, nlev, shoc_ql_2d);
  array_to_view(um.myData, ncol, nlev, u_wind_2d);
  array_to_view(vm.myData, ncol, nlev, v_wind_2d);

  view_3d horiz_wind_3d("horiz_wind",ncol,2,npack);
  view_3d qtracers_3d("qtracers",ncol,num_tracers,npack);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, npack}), KOKKOS_LAMBDA(int icol, int ilev) {
     horiz_wind_3d(icol,0,ilev) = u_wind_2d(icol,ilev);
     horiz_wind_3d(icol,1,ilev) = v_wind_2d(icol,ilev);
  });

  array_to_view(qtracers.myData, ncol, num_tracers, npack, qtracers_3d);

  SHOC::SHOCInputOutput shoc_input_output{host_dse_2d, tke_2d, thetal_2d, qw_2d,
                                         horiz_wind_3d, wthv_sec_2d, qtracers_3d,
                                         tk_2d, tkh_2d, shoc_cldfrac_2d, shoc_ql_2d};

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
  const int ntrac = ekat::npack<Spack>(num_tracers+3)*Spack::n;
  const auto policy = ekat::ExeSpaceUtils<SHOC::KT::ExeSpace>::get_default_team_policy(ncol, npack);
  ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 128+(nwind+ntrac), policy);

  const auto elapsed_microsec = SHOC::shoc_main(ncol, nlev, nlevi, nlev, 1, num_tracers, dtime, workspace_mgr,
                                                shoc_input, shoc_input_output, shoc_output, shoc_history_output);

  // get SHOC output back to CRM 
  view_to_array(shoc_input_output.tke, ncol, nlev, tke);
  view_to_array(shoc_input_output.tk, ncol, nlev, tk);
  view_to_array(shoc_input_output.tkh, ncol, nlev, tkh);
  view_to_array(shoc_input_output.qtracers, ncol, num_tracers, nlev, qtracers);

  // get precipitation variable updated
  shoc_update_precipitation(qtracers);

{auto fp=fopen("tke.txt","w");
  // update sgs_field tke
  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    sgs_field(0,k,offy_s+j,offx_s+i,icrm)      = tke(icol,nlev-(ilev+1));
    sgs_field_diag(0,k,offy_d+j,offx_d+i,icrm) = tk(icol,nlev-(ilev+1));
    sgs_field_diag(1,k,offy_d+j,offx_d+i,icrm) = tkh(icol,nlev-(ilev+1));
 
      fprintf(fp,"%d, %d, %d, %d %13.6e, %13.6e, %13.6e\n",k,j,i+nx*icrm,icrm,
      tke(icol,nlev-(ilev+1)), tk(icol,nlev-(ilev+1)), tkh(icol,nlev-(ilev+1)));
  });
fclose(fp);}

{auto fp=fopen("shoc_temp.txt","w");
  // output data 
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
    int i    = icol%nx;
    int j    = (icol/nx)%ny;
    int icrm = (icol/nx)/ny;
    int k    = ilev*Spack::n + s;
    if (k < nlev) {
        CF3D(k,j,i,icrm) = shoc_input_output.shoc_cldfrac(icol,nlev-(ilev+1))[s];
        u(k,j+offy_u,i+offx_u,icrm) = shoc_input_output.horiz_wind(icol,0,nlev-(ilev+1))[s];
        v(k,j+offy_v,i+offx_v,icrm) = shoc_input_output.horiz_wind(icol,1,nlev-(ilev+1))[s];
        tabs(k,j,i,icrm)            = (shoc_input_output.thetal(icol,nlev-(ilev+1))[s] + (latvap/cp)*qc(k,j,i,icrm))
                                      /inv_exner(icol,nlev-(ilev+1));
        t(k,j+offy_s,i+offx_s,icrm) = tabs(k,j,i,icrm) + gamaz(k,icrm)
                                    - fac_cond *( qcl(k,j,i,icrm) - qpl(k,j,i,icrm) )
                                    - fac_sub  *( qci(k,j,i,icrm) - qpi(k,j,i,icrm) );

      fprintf(fp,"%d, %d, %d, %d %13.6e, %13.6e, %13.6e, %13.6e\n",k,j,i+nx*icrm,icrm,
      t(k,j+offy_s,i+offx_s,icrm),tabs(k,j,i,icrm),inv_exner(icol,ilev),shoc_input_output.horiz_wind(icol,0,ilev)[s]);

    }
  });
fclose(fp);}
}
