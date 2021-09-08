
#include "sgs.h"

void kurant_sgs(real &cfl) {
  auto &sgs_field_diag = :: sgs_field_diag;
  auto &dz             = :: dz;
  auto &dy             = :: dy;
  auto &dx             = :: dx;
  auto &dt             = :: dt;
  auto &adzw           = :: adzw;
  auto &grdf_x         = :: grdf_x;
  auto &grdf_y         = :: grdf_y;
  auto &grdf_z         = :: grdf_z;
  auto &ncrms          = :: ncrms;

  real2d tkhmax("tkhmax",nzm,ncrms);

  // for (int k=0; k<nzm; k++) {
  //   for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    tkhmax(k,icrm) = 0.;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    yakl::atomicMax( tkhmax(k,icrm) , sgs_field_diag(1,k,offy_d+j,offx_d+i,icrm) );
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int icrm=0; icrm < ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    real dztmp = dz(icrm)*adzw(k,icrm);
    real xdir = 0.5*tkhmax(k,icrm)*grdf_x(k,icrm)*dt/(dx*dx);
    real ydir = 0.5*tkhmax(k,icrm)*grdf_y(k,icrm)*dt/(dy*dy)*YES3D;
    real zdir = 0.5*tkhmax(k,icrm)*grdf_z(k,icrm)*dt/(dztmp*dztmp);
    tkhmax(k,icrm) = max( max( xdir , ydir ) , zdir );
  });

  // Perform a max reduction over tkhmax
  yakl::ParallelMax<real,yakl::memDevice> pmax( nzm*ncrms );
  real cfl_loc = pmax( tkhmax.data() );
  cfl = max(cfl , cfl_loc);
}


void sgs_proc() {
  auto &sgs_field      = :: sgs_field;
  auto &sgs_field_diag = :: sgs_field_diag;
  auto &tke2           = :: tke2;
  auto &tk2            = :: tk2;
  auto &ncrms          = :: ncrms;

  if (dosgs) {
    tke_full(sgs_field, 0, sgs_field_diag, 0, sgs_field_diag, 1);
  }

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<dimy_s; j++) {
  //     for (int i=0; i<dimx_s; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    tke2(k,j,i,icrm) = sgs_field(0,k,j,i,icrm);
  });
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<dimy_d; j++) {
  //     for (int i=0; i<dimx_d; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,dimy_d,dimx_d,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    tk2(k,j,i,icrm) = sgs_field_diag(0,k,j,i,icrm);
  });
}


void sgs_mom() {
  diffuse_mom();
}


void sgs_scalars() {
  real2d dummy("dummy", nz, ncrms);

  diffuse_scalar(sgs_field_diag,1,t,fluxbt,fluxtt,tdiff,twsb);

  if (advect_sgs) {
    diffuse_scalar(sgs_field_diag,1,sgs_field,0,fzero,fzero,dummy,dummy);
  }

  micro_flux();

  for (int k=0; k<nmicro_fields; k++) {
    if (k==index_water_vapor || (docloud && flag_precip(k)!=1) || (doprecip && flag_precip(k)==1)) {
      diffuse_scalar(sgs_field_diag,1,micro_field,k,fluxbmk,k,fluxtmk,k,mkdiff,k,mkwsb,k);
    }
  }

#if defined(MMF_ESMT)
  diffuse_scalar(sgs_field_diag,1,u_esmt,fluxb_u_esmt,fluxt_u_esmt,u_esmt_diff,u_esmt_sgs);
  diffuse_scalar(sgs_field_diag,1,v_esmt,fluxb_v_esmt,fluxt_v_esmt,v_esmt_diff,v_esmt_sgs);
#endif
}


void sgs_init() {
  auto &sgs_field      = ::sgs_field;
  auto &sgs_field_diag = ::sgs_field_diag;
  auto &grdf_x         = ::grdf_x;
  auto &grdf_y         = ::grdf_y;
  auto &grdf_z         = ::grdf_z;
  auto &adz            = ::adz;
  auto &dx             = ::dx;
  auto &dy             = ::dy;
  auto &dz             = ::dz;
  auto &ncrms          = ::ncrms;

  if (nrestart == 0) {
    // for (int l=0; l<nsgs_fields; l++) {
    //  for (int k=0; k<nzm; k++) {
    //    for (int j=0; j<dimy_s; j++) {
    //      for (int i=0; i<dimx_s; i++) {
    //        for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<5>(nsgs_fields,nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int l, int k, int j, int i, int icrm) {
      sgs_field(l,k,j,i,icrm) = 0.0;
    });

    // for (int k=0; k<nsgs_fields_diag; k++) {
    //  for (int k=0; k<nzm; k++) {
    //    for (int j=0; j<dimy_d; j++) {
    //      for (int i=0; i<dimx_d; i++) {
    //        for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<5>(nsgs_fields_diag,nzm,dimy_d,dimx_d,ncrms) , 
                  YAKL_LAMBDA (int l, int k, int j, int i, int icrm) {
      sgs_field_diag(l,k,j,i,icrm) = 0.0;
    });
  }

  if (les) {
    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      real tmp1 = (adz(k,icrm)*dz(icrm));
      real tmp2 = (adz(k,icrm)*dz(icrm));
      grdf_x(k,icrm) = dx*dx/(tmp1*tmp1);
      grdf_y(k,icrm) = dy*dy/(tmp2*tmp2);
      grdf_z(k,icrm) = 1.0;
    });
  } else {
    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      real tmp1 = (adz(k,icrm)*dz(icrm));
      real tmp2 = (adz(k,icrm)*dz(icrm));
      grdf_x(k,icrm) = min( 16.0, dx*dx/(tmp1*tmp1));
      grdf_y(k,icrm) = min( 16.0, dy*dy/(tmp2*tmp2));
      grdf_z(k,icrm) = 1.0;
    });
  }
}
