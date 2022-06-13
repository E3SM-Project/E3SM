
#include "periodic.h"

void periodic(int flag) {
  YAKL_SCOPE( w     , ::w );
  YAKL_SCOPE( sstxy , ::sstxy );
  YAKL_SCOPE( ncrms , ::ncrms );

  if (flag == 0) {
    bound_exchange(u,nzm,1,1,1,1, 1);
    bound_exchange(v,nzm,1,1,1,1, 2);
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(ny,nx,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
      w(nz-1,j+offy_w,i+offx_w,icrm) = sstxy(j+offy_sstxy,i+offx_sstxy,icrm);
    });

    bound_exchange(w,nz,1,1,1,1, 3);
    //   for (int j=0; j<ny+2*YES3D; j++) {
    //     for (int i=0; i<nxp3; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(ny+2*YES3D,nxp3,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
      int jStart = 1-YES3D-1;
      int jInd = 1-YES3D-1 +j;
      int jEnd = ny+YES3D-1;
      int iStart = 0-1;
      int iInd = 0-1+i;
      int iEnd = nx+1 -1;
      if ( iInd>= iStart && i<= (nx-1) && jInd >= jStart && jInd <= (ny-1)) {
        sstxy(jInd+offy_sstxy,iInd+offx_sstxy,icrm) = w(nz-1,jInd+offy_w,iInd+offx_w,icrm);
      }
    });

    //   for (int j=0; j<ny+2*YES3D; j++) {
    //     for (int i=0; i<nxp3; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(ny+2*YES3D,nxp3,ncrms) , YAKL_LAMBDA (int j, int i, int icrm) {
      int jStart = 1-YES3D-1;
      int jInd = 1-YES3D-1+j;
      int jEnd = ny+YES3D-1;
      int iStart = 0-1;
      int iInd = 0-1+i;
      int iEnd = nx+1-1;
      if( iInd>= iStart && iInd<= (nx-1) && jInd >= jStart && jInd <= ny+YES3D-1) {
        w(nz-1,jInd+offy_w,iInd+offx_w,icrm) = 0.0;
      }
    });
  }

  if (flag == 2) {
    bound_exchange(u,nzm,2,3,2,2, 1);
    bound_exchange(v,nzm,2,2,2,3, 2);
    bound_exchange(w,nz,2,2,2,2, 3);
    bound_exchange(t,nzm,3,3,3,3, 4);

    for (int i=0; i<nsgs_fields; i++) {
      if (dosgs && advect_sgs) {
        bound_exchange(sgs_field,i,nzm,3,3,3,3, 4);
      }
    }

    for (int i=0; i<nmicro_fields; i++) {
      if (i == index_water_vapor || (docloud && flag_precip(i)!=1) || (doprecip && flag_precip(i)==1)) {
        bound_exchange(micro_field,i,nzm,3,3,3,3, 4);
      }
    }
#ifdef MMF_ESMT
    bound_exchange(u_esmt, nzm, 3, 3, 3, 3, 4);
    bound_exchange(v_esmt, nzm, 3, 3, 3, 3, 4);
#endif
  }

  if (flag == 3) {
    bound_exchange(t,nzm,1,1,1,1, 4);
    for (int i=0; i<nsgs_fields; i++) {
      if (dosgs && advect_sgs) {
        bound_exchange(sgs_field,i,nzm,1,1,1,1, 4);
      }
    }

    for (int i=0; i<nmicro_fields; i++) {
      if ( i==index_water_vapor || (docloud && flag_precip(i)!=1) || (doprecip && flag_precip(i)==1)) {
          bound_exchange(micro_field,i,nzm,1,1,1,1, 4);
      }
    }
#ifdef MMF_ESMT
    bound_exchange(u_esmt, nzm, 1, 1, 1, 1, 4);
    bound_exchange(v_esmt, nzm, 1, 1, 1, 1, 4);
#endif
  }

  if (flag == 4) {
    for (int i=0; i<nsgs_fields_diag; i++) {
      if (dosgs && do_sgsdiag_bound) {
        bound_exchange(sgs_field_diag,i,nzm,1,1,1,1, 5);
      }
    }
  }

}
