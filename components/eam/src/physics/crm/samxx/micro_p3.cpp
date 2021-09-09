
#include "micro_p3.h"

using namespace scream;
using namespace scream::p3;

void micro_p3_diagnose() {
  auto &qv          = :: qv;
  auto &qcl         = :: qcl;
  auto &qci         = :: qci;
  auto &qpl         = :: qpl;
  auto &qpi         = :: qpi;
  auto &micro_field = :: micro_field;

  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    qv(k,j,i,icrm)  = micro_field(ixqv,    k,j+offy_s,i+offx_s,icrm);
    qcl(k,j,i,icrm) = micro_field(ixcldliq,k,j+offy_s,i+offx_s,icrm);
    qci(k,j,i,icrm) = micro_field(ixcldice,k,j+offy_s,i+offx_s,icrm);
    qpl(k,j,i,icrm) = micro_field(ixrain,  k,j+offy_s,i+offx_s,icrm);
    qpi(k,j,i,icrm) = micro_field(ixcldrim,k,j+offy_s,i+offx_s,icrm);
  });
}


void get_cloud_fraction(int its, int ite, int kts, int kte, 
                       real2d& ast, real2d& qc, real2d& qr, real2d& qi, 
                       std::string& method, 
                       real2d& cld_frac_i, real2d& cld_frac_l, real2d& cld_frac_r, 
                       real2d& cldm) 
{
  int nk = kte-kts+1;
  int ni = ite-its+1;

  parallel_for( SimpleBounds<2>(ni,nk) , YAKL_LAMBDA (int i, int k) {
      int i0 = i+its;
      int k0 = kte-kts-k;

      cldm(i0,k0)  = max(ast(i0,k0), mincld);
      cld_frac_i(i0,k0) = max(ast(i0,k0), mincld);
      cld_frac_l(i0,k0) = max(ast(i0,k0), mincld);
      cld_frac_r(i0,k0) = cldm(i0,k0);
  });

  if (method == "in_cloud") {

     parallel_for( SimpleBounds<2>(ni,nk) , YAKL_LAMBDA (int i, int k) {
       int i0 = i+its;
       int k0 = kte-kts-k;

       // in_cloud means that precip_frac (cld_frac_r) = cloud (cldm) frac when cloud mass
       // is present. Below cloud, precip frac is equal to the cloud
       // fraction from the last layer that had cloud. Since presence or
       // absence of cloud is defined as mass > qsmall, sub-cloud precip
       // frac for the in_cloud method tends to be very small and is
       // very sensitive to tiny changes in condensate near cloud base.
       if (qc(i0,k0) < qsmall && qi(i0,k0) < qsmall) {
         // max(cld_frac_r above and cld_frac_r for this layer) is taken here
         // because code is incapable of handling cld_frac_r<cldm for a
         // given grid cell
         cld_frac_r(i0,k0) = max(cld_frac_r(i0,k0-1),cld_frac_r(i0,k0));
       }
     });

  } else if (method == "max_overlap") {
      // max overlap is the max cloud fraction in all layers above which are
      // connected to this one by a continuous band of precip mass. If
      // there's no precip mass falling into a cell, it's precip frac is equal
      // to the cloud frac, which is probably ~zero.

      // IF rain or ice mix ratios are smaller than threshold,
      // then leave cld_frac_r as cloud fraction at current level

      parallel_for( SimpleBounds<2>(ni,nk) , YAKL_LAMBDA (int i, int k) {
        int i0 = i+its;
        int k0 = kte-kts-k;

        if (qr(i0,k0-1) >= qsmall || qi(i0,k0-1) >= qsmall) {
           cld_frac_r(i0,k0) = max(cld_frac_r(i0,k0-1),cld_frac_r(i0,k0));
        }
      });
 }
}

//
// main micro_p3 microproc
//
void micro_p3_proc() {
 using P3F        = Functions<Real, DefaultDevice>;
 using KT         = typename P3F::KT;
 using ExeSpace   = typename KT::ExeSpace;
 using MemberType = typename P3F::MemberType;
 using Spack      = typename P3F::Spack;
 using view_1d    = typename P3F::view_1d<Spack>;
 using uview_1d   = typename P3F::uview_1d<Spack>;
 using view_2d    = typename P3F::view_2d<Spack>;
 using sview_1d   = typename P3F::view_1d<Real>;
 using sview_2d   = typename P3F::view_2d<Real>;

 auto &micro_field        = :: micro_field;
 auto &rho                = :: rho;
 auto &dt                 = :: dt;
 auto &pres               = :: pres;
 auto &pdel               = :: crm_input_pdel;
 auto &tabs               = :: tabs;
 auto &ncrms              = :: ncrms;
 auto &dz                 = :: dz;
 auto &longitude0         = :: longitude0;
 auto &latitude0          = :: latitude0;
 auto &z0                 = :: z0;
 auto &nc_nuceat_tend     = :: nc_nuceat_tend;
 auto &nccn_prescribed    = :: nccn_prescribed;
 auto &ni_activated       = :: ni_activated;
 auto &relvar             = :: relvar;
 auto &diag_eff_radius_qc = :: diag_eff_radius_qc;
 auto &diag_eff_radius_qi = :: diag_eff_radius_qi;
 auto &precip_total_tend  = :: precip_total_tend;
 auto &nevapr             = :: nevapr;
 auto &qr_evap_tend       = :: qr_evap_tend;
 auto &mu                 = :: mu;
 auto &lambdac            = :: lambdac;
 auto &t_prev             = :: t_prev;
 auto &qv_prev            = :: qv_prev;
 auto &docloud            = :: docloud;
 auto &ast                = :: ast;

 // output 
 auto &qv2qi_depos_tend   = :: qv2qi_depos_tend;
 auto &precip_liq_surf    = :: precip_liq_surf;
 auto &precip_ice_surf    = :: precip_ice_surf;
 auto &rho_qi             = :: rho_qi;
 auto &precip_liq_flux    = :: precip_liq_flux;
 auto &precip_ice_flux    = :: precip_ice_flux;

 auto &liq_ice_exchange  =  :: liq_ice_exchange;
 auto &vap_liq_exchange  =  :: vap_liq_exchange;
 auto &vap_ice_exchange  =  :: vap_ice_exchange;

 // output
 auto &precsfc           = :: precsfc;
 auto &prec_xy           = :: prec_xy;

 const int nlev  = nzm;
 const int ncol  = ncrms*nx*ny;
 const int npack = ekat::npack<Spack>(nlev);

 real2d qc_in("qc",ncol, nlev);
 real2d nc_in("nc",ncol, nlev);
 real2d qr_in("qr",ncol, nlev);
 real2d nr_in("nr",ncol, nlev);
 real2d qi_in("qi",ncol, nlev);
 real2d qm_in("qm",ncol, nlev);
 real2d ni_in("ni",ncol, nlev);
 real2d bm_in("bm",ncol, nlev);
 real2d qv_in("qv",ncol, nlev);
 real2d th_in("th",ncol, nlev);
 
 real2d nc_nuceat_tend_in("nuceat",ncol, nlev);
 real2d nccn_prescribed_in("nccn_prescribed",ncol, nlev);
 real2d ni_activated_in("ni_act",ncol, nlev);
 real2d inv_qc_relvar_in("inv_qc",ncol, nlev);
 real2d cld_frac_i_in("cld_frac_i",ncol, nlev);
 real2d cld_frac_l_in("cld_frac_l",ncol, nlev);
 real2d cld_frac_r_in("cld_frac_r",ncol, nlev);
 real2d pres_in("pres", ncol, nlev);
 real2d dz_in("dz", ncol, nlev);
 real2d dpres_in("dpres",ncol, nlev);
 real2d exner_in("exner",ncol, nlev);
 real2d qv_prev_in("qv_prev",ncol, nlev);
 real2d t_prev_in("t_prev",ncol, nlev);
 
 real2d ast_in("ast_in",ncol, nlev);
 real2d cldm_in("cldm_in", ncol, nlev);
 
 real2d qv2qi_depos_tend_in("qv2qi",ncol, nlev);
 real2d precip_liq_surf_in("precip_liq",ncol, nlev);
 real2d precip_ice_surf_in("precip_ice",ncol, nlev);
 real2d diag_eff_radius_qc_in("diag_eff_qc",ncol, nlev);
 real2d diag_eff_radius_qi_in("diag_eff_qi",ncol, nlev);
 real2d rho_qi_in("rho_qi",ncol, nlev);
 real2d precip_liq_flux_in("precip_liq_flux",ncol, nlev);
 real2d precip_ice_flux_in("precip_ice_flux",ncol, nlev);

 // p3 output variables
 real2d qv2qi_depos_tend_out("qv2qi_depos_tend", ncol, nlev);
 real2d diag_eff_radius_qc_out("diag_eff_radius_qc", ncol, nlev);
 real2d diag_eff_radius_qi_out("diag_eff_radius_qi", ncol, nlev);
 real2d rho_qi_out("rho_qi", nzm, nx);
 real2d precip_liq_flux_out("precip_liq_flux", ncol, nlev);
 real2d precip_ice_flux_out("precip_ice_flux", ncol, nlev);

 real1d precip_liq_surf_out("precip_liq_surf_d", nlev);
 real1d precip_ice_surf_out("precip_ice_surf_d", nlev);

printf("micro_p3_00: nx=%d, ny=%d, nzm=%d, ncrms=%d, ncol=%d, nlev=%d\n",nx,ny,nzm,ncrms,ncol,nlev);

 parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
   int icol = i+nx*(j+ny*icrm);
   int ilev = k;
   exner_in(icol, ilev) = 1./std::pow((pres(k,icrm)*1.0e-5), (rgas/cp));
   th_in(icol, ilev) = tabs(k,j,i,icrm)*exner_in(icol, ilev);
 });

 parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+ny*icrm);
    int ilev = k;
    qc_in(icol,ilev) = micro_field(ixcldliq, k,j+offy_s,i+offx_s,icrm); 
    nc_in(icol,ilev) = micro_field(ixnumliq, k,j+offy_s,i+offx_s,icrm);
    qr_in(icol,ilev) = micro_field(ixrain,   k,j+offy_s,i+offx_s,icrm);
    nr_in(icol,ilev) = micro_field(ixnumrain,k,j+offy_s,i+offx_s,icrm);
    qi_in(icol,ilev) = micro_field(ixcldice, k,j+offy_s,i+offx_s,icrm);
    qm_in(icol,ilev) = micro_field(ixcldrim, k,j+offy_s,i+offx_s,icrm);
    ni_in(icol,ilev) = micro_field(ixnumice, k,j+offy_s,i+offx_s,icrm);
    bm_in(icol,ilev) = micro_field(ixrimvol, k,j+offy_s,i+offx_s,icrm);
    qv_in(icol,ilev) = micro_field(ixqv,     k,j+offy_s,i+offx_s,icrm);
  });

  view_2d qc_d("qc", ncol, npack),
          nc_d("nc", ncol, npack),
          qr_d("qr", ncol, npack),
          nr_d("nr", ncol, npack),
          qi_d("qi", ncol, npack),
          qm_d("qm", ncol, npack),
          ni_d("ni", ncol, npack),
          bm_d("bm", ncol, npack),
          qv_d("qv", ncol, npack),
          th_atm_d("th", ncol, npack);

  array_to_view(qc_in.myData, ncol, nlev, qc_d);
  array_to_view(nc_in.myData, ncol, nlev, nc_d);
  array_to_view(qr_in.myData, ncol, nlev, qr_d);
  array_to_view(nr_in.myData, ncol, nlev, nr_d);
  array_to_view(qi_in.myData, ncol, nlev, qi_d);
  array_to_view(qm_in.myData, ncol, nlev, qm_d);
  array_to_view(ni_in.myData, ncol, nlev, ni_d);
  array_to_view(bm_in.myData, ncol, nlev, bm_d);
  array_to_view(qv_in.myData, ncol, nlev, qv_d);
  array_to_view(th_in.myData, ncol, nlev, th_atm_d);
   
  P3F::P3PrognosticState prog_state{qc_d, nc_d, qr_d, nr_d, qi_d, qm_d,
                                    ni_d, bm_d, qv_d, th_atm_d};

  parallel_for( SimpleBounds<4>(nzm, ny, nx, ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    int icol = i+nx*(j+icrm*ny);
    int ilev = k;
    nccn_prescribed_in(icol,ilev) = nccn_prescribed(k,icrm)*0.0; // TODO: we zero out the nccn_prescribed because of the missing of model
    nc_nuceat_tend_in(icol,ilev)  = nc_nuceat_tend(k,icrm)*0.0;  // TODO: we zero out the nc_nuceat_tend because of the missing of model
    ni_activated_in(icol,ilev)    = ni_activated(k,icrm)*0.0;    // TODO: we zero out the ni_activated because of the missing of model
    inv_qc_relvar_in(icol,ilev)   = relvar(k,icrm);
    pres_in(icol,ilev)            = pres(k,icrm);
    dz_in(icol,ilev)              = dz(icrm);
    dpres_in(icol,ilev)           = pdel(k,icrm);
    ast_in(icol,ilev)             = ast(k,icrm);
    qv_prev_in(icol,ilev)         = qv_prev(k,j,i,icrm);
    t_prev_in(icol,ilev)          = t_prev(k,j,i,icrm);
  });

  std::string method("in_cloud");
  get_cloud_fraction(0, ncol, 0, nlev,
                     ast_in, qc_in, qr_in, qi_in, method,
                     cld_frac_i_in, cld_frac_l_in, cld_frac_r_in, cldm_in);

  view_2d nc_nuceat_tend_d("nc_nuceat_tend", ncol, npack),
          nccn_prescribed_d("nccn_prescribed", ncol, npack),
          ni_activated_d("ni_activated", ncol, npack),
          inv_qc_relvar_d("inv_qc_relvar", ncol, npack),
          pres_d("pres", ncol, npack),
          dz_d("dz", ncol, npack),
          dpres_d("dpres", ncol, npack),
          exner_d("exner", ncol, npack),
          t_prev_d("t_prev", ncol, npack),
          qv_prev_d("qv_prev", ncol, npack),
          cld_frac_i_d("cld_frac_i", ncol, npack),
          cld_frac_l_d("cld_frac_l", ncol, npack),
          cld_frac_r_d("cld_frac_r", ncol, npack);

  array_to_view(nc_nuceat_tend_in.myData, ncol, nlev, nc_nuceat_tend_d);
  array_to_view(nccn_prescribed_in.myData, ncol, nlev, nccn_prescribed_d);
  array_to_view(ni_activated_in.myData, ncol, nlev, ni_activated_d);
  array_to_view(inv_qc_relvar_in.myData, ncol, nlev, inv_qc_relvar_d);
  array_to_view(pres_in.myData, ncol, nlev, pres_d);
  array_to_view(dz_in.myData, ncol, nlev, dz_d);
  array_to_view(dpres_in.myData, ncol, nlev, dpres_d);
  array_to_view(exner_in.myData, ncol, nlev, exner_d);
  array_to_view(t_prev_in.myData, ncol, nlev, t_prev_d);
  array_to_view(qv_prev_in.myData, ncol, nlev, qv_prev_d);
  array_to_view(cld_frac_i_in.myData, ncol, nlev, cld_frac_i_d);
  array_to_view(cld_frac_l_in.myData, ncol, nlev, cld_frac_l_d);
  array_to_view(cld_frac_r_in.myData, ncol, nlev, cld_frac_r_d);

  P3F::P3DiagnosticInputs diag_inputs{nc_nuceat_tend_d, nccn_prescribed_d, ni_activated_d, inv_qc_relvar_d, cld_frac_i_d,
                                      cld_frac_l_d, cld_frac_r_d, pres_d, dz_d, dpres_d,
                                      exner_d, qv_prev_d, t_prev_d};

  view_2d qv2qi_depos_tend_d("qv2qi_depos_tend", ncol, npack),
          diag_eff_radius_qc_d("diag_eff_radius_qc", ncol, npack),
          diag_eff_radius_qi_d("diag_eff_radius_qi", ncol, npack),
          rho_qi_d("rho_qi", ncol, npack),
          precip_liq_flux_d("precip_liq_flux", ncol, npack),
          precip_ice_flux_d("precip_ice_flux", ncol, npack);

  sview_1d precip_liq_surf_d("precip_liq_surf_d", npack), 
           precip_ice_surf_d("precip_ice_surf_d", npack);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, npack}), KOKKOS_LAMBDA(int k, int j) {
     qv2qi_depos_tend_d(k,j)[0] = 0.;
     diag_eff_radius_qc_d(k,j)[0] = 0.;
     diag_eff_radius_qi_d(k,j)[0] = 0.;
     rho_qi_d(k,j)[0] = 0.;
     precip_liq_flux_d(k,j)[0] = 0.;
     precip_ice_flux_d(k,j)[0] = 0.;
  });

  Kokkos::parallel_for("precip", nlev, KOKKOS_LAMBDA (const int& i) {
    precip_liq_surf_d(i) = 0.;
    precip_ice_surf_d(i) = 0.;
  });

  P3F::P3DiagnosticOutputs diag_outputs {qv2qi_depos_tend_d, precip_liq_surf_d,
                                         precip_ice_surf_d, diag_eff_radius_qc_d, diag_eff_radius_qi_d,
                                         rho_qi_d,precip_liq_flux_d, precip_ice_flux_d};

  int it;
  int its{0};
  int ite{ncrms*crm_nx*crm_ny};
  int kts{0}; 
  int kte{nzm}; 
  bool do_predict_nc, do_prescribed_CCN;
  sview_2d col_location_d("col_location_d", ncol, 3);

  do_predict_nc = false;
  do_prescribed_CCN = false;

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0}, {nzm, ny, nx, ncrms}), KOKKOS_LAMBDA(int k, int j, int i, int icrm) {
// Kokkos::parallel_for("col_location", ncol, KOKKOS_LAMBDA (const int& icol) {
     int icol = i+nx*(j+icrm*ny);
     col_location_d(icol, 1) = z0.myData[icrm];
     col_location_d(icol, 2) = longitude0.myData[icrm];
     col_location_d(icol, 3) = latitude0.myData[icrm];
  });

  P3F::P3Infrastructure infrastructure{dt, it, its, ite, kts, kte,
                                       do_predict_nc, do_prescribed_CCN, col_location_d};

  view_2d liq_ice_exchange_d("liq_ice_exchange_d", ncol, npack),
          vap_liq_exchange_d("vap_liq_exchange_d", ncol, npack),
          vap_ice_exchange_d("vap_ice_exchange_d", ncol, npack);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, npack}), KOKKOS_LAMBDA(int k, int j) {
     liq_ice_exchange_d(k,j) = 0.;
     vap_liq_exchange_d(k,j) = 0.;
     vap_ice_exchange_d(k,j) = 0.;
  });

  P3F::P3HistoryOnly history_only {liq_ice_exchange_d, vap_liq_exchange_d,
                                   vap_ice_exchange_d};

  const int nlev_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol, nlev_pack);
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlev_pack, 52, policy);


  auto elapsed_time = P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
                                   history_only, workspace_mgr, ncol, nlev);

     
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, nlev, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
     int i    = icol%nx;
     int j    = (icol/nx)%ny;
     int icrm = (icol/nx)/ny;
     int k    = ilev*Spack::n+s;
     if (k <= nlev) {
         precsfc(j,i,icrm) = precsfc(j,i,icrm)+(diag_outputs.precip_liq_surf(ilev)+diag_outputs.precip_ice_surf(ilev))*1000.0*dt/dz.myData[icrm];
         prec_xy(j,i,icrm) = prec_xy(j,i,icrm)+(diag_outputs.precip_liq_surf(ilev)+diag_outputs.precip_ice_surf(ilev))*1000.0*dt/dz.myData[icrm];
       }
  });

  // update microfield
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
     int i    = icol%nx;
     int j    = (icol/nx)%ny;
     int icrm = (icol/nx)/ny;
     int k    = ilev*Spack::n + s;
     if (k <= nlev) {
        int d1 = micro_field.dimension[1];
        int d2 = micro_field.dimension[2];
        int d3 = micro_field.dimension[3];
        int d4 = micro_field.dimension[4];
        int d5 = micro_field.dimension[5];
        micro_field.myData[((((ixcldliq*d1  +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.qc(icol,ilev)[s];
        micro_field.myData[((((ixnumliq*d1  +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.nc(icol,ilev)[s];
        micro_field.myData[((((ixrain*d1    +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.qr(icol,ilev)[s];
        micro_field.myData[((((ixnumrain*d1 +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.nr(icol,ilev)[s];
        micro_field.myData[((((ixcldice*d1  +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.qi(icol,ilev)[s];
        micro_field.myData[((((ixcldrim*d1  +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.qm(icol,ilev)[s];
        micro_field.myData[((((ixnumice*d1  +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.ni(icol,ilev)[s];
        micro_field.myData[((((ixrimvol*d1  +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.bm(icol,ilev)[s];
        micro_field.myData[((((ixqv*d1      +k)*d2+j+offy_s)*d3+i+offx_s)*d4+k)*d5+icrm] = prog_state.qv(icol,ilev)[s];        
     }
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
     int i    = icol%nx;
     int j    = (icol/nx)%ny;
     int icrm = (icol/nx)/ny;
     int k    = j*Spack::n + s;
     if (k <= nlev) {
        qv2qi_depos_tend(k,icrm)   = diag_outputs.qv2qi_depos_tend(icol,ilev)[s];
        diag_eff_radius_qc(k,icrm) = diag_outputs.diag_eff_radius_qc(icol,ilev)[s];
        diag_eff_radius_qi(k,icrm) = diag_outputs.diag_eff_radius_qi(icol,ilev)[s];
        rho_qi(k,icrm)             = diag_outputs.rho_qi(icol,ilev)[s];
        precip_liq_flux(k,icrm)    = diag_outputs.precip_liq_flux(icol,ilev)[s];
        precip_ice_flux(k,icrm)    = diag_outputs.precip_ice_flux(icol,ilev)[s];
     }
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
     int i    = icol%nx;
     int j    = (icol/nx)%ny;
     int icrm = (icol/nx)/ny;
     int k    = ilev*Spack::n + s;
     if (k <= nlev) {
         liq_ice_exchange(k,icrm) = history_only.liq_ice_exchange(icol,ilev)[s];
         vap_liq_exchange(k,icrm) = history_only.vap_liq_exchange(icol,ilev)[s];
         vap_ice_exchange(k,icrm) = history_only.vap_ice_exchange(icol,ilev)[s];
     }
  });
     
  if (docloud)  micro_p3_diagnose();   // leave this line here
}

