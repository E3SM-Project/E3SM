
#include "ice_fall.h"

void ice_fall() {
  auto &qcl           = :: qcl;
  auto &qci           = :: qci;
  auto &tabs          = :: tabs;
  auto &qifall        = :: qifall;
  auto &tlatqi        = :: tlatqi;
  auto &dtn           = :: dtn;
  auto &adz           = :: adz;
  auto &dz            = :: dz;
  auto &rho           = :: rho;
  auto &micro_field   = :: micro_field;
  auto &t             = :: t;
  auto &ncrms         = :: ncrms;
  auto &precsfc       = :: precsfc;
  auto &precssfc      = :: precssfc;

  int1d  kmax("kmax",ncrms);
  int1d  kmin("kmin",ncrms);
  real4d fz  ("fz"  ,nz,ny,nx,ncrms);

  // for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( ncrms , YAKL_LAMBDA (int icrm) {
    kmax(icrm) = -1;
    kmin(icrm) = nzm;
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    for(int k=0; k < nzm; k++) {
      if(qcl(k,j,i,icrm)+qci(k,j,i,icrm) > 0.0 && tabs(k,j,i,icrm) < 273.15) {
        yakl::atomicMin(kmin(icrm),k);
        yakl::atomicMax(kmax(icrm),k);
      }
    }
  });

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    qifall(k,icrm) = 0.0;
    tlatqi(k,icrm) = 0.0;
  });

  if(index_cloud_ice == -1) { return;}

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nz,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    fz(k,j,i,icrm) = 0.0;
  });


  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nz,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (k >= max(0,kmin(icrm)-1) && k <= kmax(icrm) ) {
      // Set up indices for x-y planes above and below current plane.
      int kc = min(k+1,nzm-1);
      int kb = max(k-1,0    );

      // CFL number based on grid spacing interpolated to interface i,j,k-1/2
      real coef = dtn/(0.5*(adz(kb,icrm)+adz(k,icrm))*dz(icrm));

      // Compute cloud ice density in this cell and the ones above/below.
      // Since cloud ice is falling, the above cell is u(icrm,upwind),
      // this cell is c (center) and the one below is d (downwind).
      real qiu = rho(kc,icrm)*qci(kc,j,i,icrm);
      real qic = rho(k ,icrm)*qci(k ,j,i,icrm);
      real qid = rho(kb,icrm)*qci(kb,j,i,icrm);

      // Ice sedimentation velocity depends on ice content. The fiting is
      // based on the data by Heymsfield (JAS,2003). -Marat
      real vt_ice = min( 0.4 , 8.66 * pow( (max(0.,qic)+1.e-10) , 0.24) );   // Heymsfield (JAS, 2003, p.2607)

      // Use MC flux limiter in computation of flux correction.
      // (MC = monotonized centered difference).
      //         if (qic.eq.qid) then
      real tmp_phi;
      if ( abs(qic-qid) < 1.0e-25 ) {  // when qic, and qid is very small, qic_qid can still be zero
        // even if qic is not equal to qid. so add a fix here +++mhwang
        tmp_phi = 0.;
      } else {
        real tmp_theta = (qiu-qic) / (qic-qid);
        tmp_phi = max( 0. , min( 0.5*(1.+tmp_theta) , min( 2. , 2.*tmp_theta ) ) );
      }

      // Compute limited flux.
      // Since falling cloud ice is a 1D advection problem, this
      // flux-limited advection scheme is monotonic.
      fz(k,j,i,icrm) = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid));
    }
  });

  // for (int j=0; j<ny; j++) {
  //  for (int i=0; i<nx; i++) {
  //    for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    fz(nz-1,j,i,icrm) = 0.0;
  });

  int constexpr ici = index_cloud_ice;

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nz,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if ( k >= max(0,kmin(icrm)-2) && k <= kmax(icrm) ) {
      real coef = dtn/(dz(icrm)*adz(k,icrm)*rho(k,icrm));
      // The cloud ice increment is the difference of the fluxes.
      real dqi  = coef*(fz(k,j,i,icrm)-fz(k+1,j,i,icrm));
      // Add this increment to both non-precipitating and total water.
      micro_field(ici,k,offy_s+j,offx_s+i,icrm) += dqi;
      // Include this effect in the total moisture budget.
      //$acc atomic update
      yakl::atomicAdd(qifall(k,icrm),dqi);

      // The latent heat flux induced by the falling cloud ice enters
      // the liquid-ice static energy budget in the same way as the
      // precipitation.  Note: use latent heat of sublimation.
      real lat_heat = (fac_cond+fac_fus)*dqi;
      // Add divergence of latent heat flux to liquid-ice static energy.
      t(k,offy_s+j,offx_s+i,icrm) -= lat_heat;
      // Add divergence to liquid-ice static energy budget.
      yakl::atomicAdd(tlatqi(k,icrm),-lat_heat);
    }
  });

  // for (int j=0; j<ny; j++) {
  //    for (int i=0; i<nx; i++) {
  //      for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
    real coef = dtn/dz(icrm);
    real dqi = -coef*fz(0,j,i,icrm);
    precsfc (j,i,icrm) = precsfc (j,i,icrm)+dqi;
    precssfc(j,i,icrm) = precssfc(j,i,icrm)+dqi;
  });
}


