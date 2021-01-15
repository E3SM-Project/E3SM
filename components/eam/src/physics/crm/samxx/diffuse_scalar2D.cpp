
#include "diffuse_scalar2D.h"

void diffuse_scalar2D(real4d &field, real3d &fluxb, real3d &fluxt, real5d &tkh,
                      int ind_tkh, real2d &flux) {
  auto &dx     = ::dx;
  auto &rhow   = ::rhow;
  auto &adzw   = ::adzw;
  auto &adz    = ::adz; 
  auto &dz     = ::dz; 
  auto &dtn    = ::dtn;
  auto &rho    = ::rho;
  auto &grdf_x = ::grdf_x;
  auto &grdf_z = ::grdf_z;
  auto &ncrms  = ::ncrms;

  if (dosgs || docolumn) {
    real rdx2=1.0/(dx*dx);
    int constexpr j=0;
    int constexpr offx_flx = 1;
    int constexpr offz_flx = 1;

    real4d flx("flx", nzm+1, 1, nx+1, ncrms);
    real4d dfdt("dfdt", nzm, ny, nx, ncrms);

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      dfdt(k,j,i,icrm)=0.0;
    });

    if (!docolumn) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<nx+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        real rdx5=0.5*rdx2 * grdf_x(k,icrm);
        int ic=i+1;
        real tkx=rdx5*(tkh(ind_tkh,k,j,i+offx_d-1,icrm)+tkh(ind_tkh,k,j,ic+offx_d-1,icrm));
        flx(k+offz_flx,j,i,icrm)=-tkx*(field(k,j,ic+offx_s-1,icrm)-field(k,j,i+offx_s-1,icrm));
      });

      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<nx; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        int ib=i-1;
        dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx(k+offz_flx,j,i+offx_flx,icrm)-flx(k+offz_flx,j,ib+offx_flx,icrm));
      });
    }

    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(k,icrm) = 0.0;
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      if (k <= nzm-2) {
        int kc=k+1;
        real rhoi = rhow(kc,icrm)/adzw(kc,icrm);
        real rdz2=1.0/(dz(icrm)*dz(icrm));
        real rdz5=0.5*rdz2 * grdf_z(k,icrm);
        real tkz=rdz5*(tkh(ind_tkh,k,j,i+offx_d,icrm)+tkh(ind_tkh,kc,j,i+offx_d,icrm));
        flx(k+offz_flx,j,i+offx_flx,icrm)=-tkz*(field(kc,j,i+offx_s,icrm)-field(k,j,i+offx_s,icrm))*rhoi;
        yakl::atomicAdd(flux(kc,icrm), flx(k+offz_flx,j,i+offx_flx,icrm));
      } else if (k == nzm-1) {
        real tmp=1.0/adzw(nz-1,icrm);
        real rdz=1.0/dz(icrm);
        flx(0,j,i+offx_flx,icrm)=fluxb(j,i,icrm)*rdz*rhow(0,icrm);
        flx(nzm-1+offz_flx,j,i+offx_flx,icrm)=fluxt(j,i,icrm)*rdz*tmp*rhow(nz-1,icrm);
        yakl::atomicAdd(flux(0,icrm),flx(0,j,i+offx_flx,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kb=k-1;
      real rhoi = 1.0/(adz(k,icrm)*rho(k,icrm));
      dfdt(k,j,i,icrm)=dtn*(dfdt(k,j,i,icrm)-(flx(k+offz_flx,j,i+offx_flx,icrm)-flx(kb+offz_flx,j,i+offx_flx,icrm))*rhoi);
      field(k,j,i+offx_s,icrm)=field(k,j,i+offx_s,icrm) + dfdt(k,j,i,icrm);
    });
  }
}


void diffuse_scalar2D(real5d &field, int ind_field, real3d &fluxb, real3d &fluxt,
                      real5d &tkh, int ind_tkh, real2d &flux) {
  auto &dx     = ::dx;
  auto &rhow   = ::rhow;
  auto &adzw   = ::adzw;
  auto &adz    = ::adz; 
  auto &dz     = ::dz; 
  auto &dtn    = ::dtn;
  auto &rho    = ::rho;
  auto &grdf_x = ::grdf_x;
  auto &grdf_z = ::grdf_z;
  auto &ncrms  = ::ncrms;

  if (dosgs || docolumn) {
    real rdx2=1.0/(dx*dx);
    int constexpr j=0;
    int constexpr offx_flx = 1;
    int constexpr offz_flx = 1;

    real4d flx("flx", nzm+1, 1, nx+1, ncrms);
    real4d dfdt("dfdt", nzm, ny, nx, ncrms);

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      dfdt(k,j,i,icrm)=0.0;
    });

    if (!docolumn) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<nx+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        real rdx5=0.5*rdx2 * grdf_x(k,icrm);
        int ic=i+1;
        real tkx=rdx5*(tkh(ind_tkh,k,j,i+offx_d-1,icrm)+tkh(ind_tkh,k,j,ic+offx_d-1,icrm));
        flx(k+offz_flx,j,i,icrm)=-tkx*(field(ind_field,k,j,ic+offx_s-1,icrm)-field(ind_field,k,j,i+offx_s-1,icrm));
      });

      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<nx; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        int ib=i-1;
        dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx(k+offz_flx,j,i+offx_flx,icrm)-flx(k+offz_flx,j,ib+offx_flx,icrm));
      });
    }

    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(k,icrm) = 0.0;
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      if (k <= nzm-2) {
        int kc=k+1;
        real rhoi = rhow(kc,icrm)/adzw(kc,icrm);
        real rdz2=1.0/(dz(icrm)*dz(icrm));
        real rdz5=0.5*rdz2 * grdf_z(k,icrm);
        real tkz=rdz5*(tkh(ind_tkh,k,j,i+offx_d,icrm)+tkh(ind_tkh,kc,j,i+offx_d,icrm));
        flx(k+offz_flx,j,i+offx_flx,icrm)=-tkz*(field(ind_field,kc,j,i+offx_s,icrm)-
                                                field(ind_field,k,j,i+offx_s,icrm))*rhoi;
        yakl::atomicAdd(flux(kc,icrm), flx(k+offz_flx,j,i+offx_flx,icrm));
      } else if (k == nzm-1) {
        real tmp=1.0/adzw(nz-1,icrm);
        real rdz=1.0/dz(icrm);
        flx(0,j,i+offx_flx,icrm)=fluxb(j,i,icrm)*rdz*rhow(0,icrm);
        flx(nzm-1+offz_flx,j,i+offx_flx,icrm)=fluxt(j,i,icrm)*rdz*tmp*rhow(nz-1,icrm);
        yakl::atomicAdd(flux(0,icrm),flx(0,j,i+offx_flx,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kb=k-1;
      real rhoi = 1.0/(adz(k,icrm)*rho(k,icrm));
      dfdt(k,j,i,icrm)=dtn*(dfdt(k,j,i,icrm)-(flx(k+offz_flx,j,i+offx_flx,icrm)-flx(kb+offz_flx,j,i+offx_flx,icrm))*rhoi);
      field(ind_field,k,j,i+offx_s,icrm)=field(ind_field,k,j,i+offx_s,icrm) + dfdt(k,j,i,icrm);
    });
  }
}

void diffuse_scalar2D(real5d &field, int ind_field, real4d &fluxb, int ind_fluxb, real4d &fluxt,
                      int ind_fluxt, real5d &tkh, int ind_tkh, real3d &flux, int ind_flux) {
  auto &dx            = :: dx;
  auto &rhow          = :: rhow;
  auto &adzw          = :: adzw;
  auto &adz           = :: adz; 
  auto &dz            = :: dz; 
  auto &dtn           = :: dtn;
  auto &rho           = :: rho;
  auto &grdf_x        = :: grdf_x;
  auto &grdf_z        = :: grdf_z;
  auto &ncrms         = :: ncrms;

  if (dosgs || docolumn) {
    real rdx2=1.0/(dx*dx);
    int constexpr j=0;
    int constexpr offx_flx = 1;
    int constexpr offz_flx = 1;

    real4d flx("flx", nzm+1, 1, nx+1, ncrms);
    real4d dfdt("dfdt", nzm, ny, nx, ncrms);

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      dfdt(k,j,i,icrm)=0.0;
    });

    if (!docolumn) {
      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<nx+1; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        real rdx5=0.5*rdx2 * grdf_x(k,icrm);
        int ic=i+1;
        real tkx=rdx5*(tkh(ind_tkh,k,j,i+offx_d-1,icrm)+tkh(ind_tkh,k,j,ic+offx_d-1,icrm));
        flx(k+offz_flx,j,i,icrm)=-tkx*(field(ind_field,k,j,ic+offx_s-1,icrm)-field(ind_field,k,j,i+offx_s-1,icrm));
      });

      // for (int k=0; k<nzm; k++) {
      //  for (int i=0; i<nx; i++) {
      //    for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
        int ib=i-1;
        dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx(k+offz_flx,j,i+offx_flx,icrm)-flx(k+offz_flx,j,ib+offx_flx,icrm));
      });
    }

    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(ind_flux,k,icrm) = 0.0;
    });

    // for (int k=0; k<nzm; k++) {
    //  for (int i=0; i<nx; i++) {
    //    for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      if (k <= nzm-2) {
        int kc=k+1;
        real rhoi = rhow(kc,icrm)/adzw(kc,icrm);
        real rdz2=1.0/(dz(icrm)*dz(icrm));
        real rdz5=0.5*rdz2 * grdf_z(k,icrm);
        real tkz=rdz5*(tkh(ind_tkh,k,j,i+offx_d,icrm)+tkh(ind_tkh,kc,j,i+offx_d,icrm));
        flx(k+offz_flx,j,i+offx_flx,icrm)=-tkz*(field(ind_field,kc,j,i+offx_s,icrm)-
                                                field(ind_field,k,j,i+offx_s,icrm))*rhoi;
        yakl::atomicAdd(flux(ind_flux,kc,icrm), flx(k+offz_flx,j,i+offx_flx,icrm));
      }
      else if (k == nzm-1) {
        real tmp=1.0/adzw(nz-1,icrm);
        real rdz=1.0/dz(icrm);
        flx(0,j,i+offx_flx,icrm)=fluxb(ind_fluxb,j,i,icrm)*rdz*rhow(0,icrm);
        flx(nzm-1+offz_flx,j,i+offx_flx,icrm)=fluxt(ind_fluxt,j,i,icrm)*rdz*tmp*rhow(nz-1,icrm);
        yakl::atomicAdd(flux(ind_flux,0,icrm),flx(0,j,i+offx_flx,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int kb=k-1;
      real rhoi = 1.0/(adz(k,icrm)*rho(k,icrm));
      dfdt(k,j,i,icrm)=dtn*(dfdt(k,j,i,icrm)-(flx(k+offz_flx,j,i+offx_flx,icrm)-flx(kb+offz_flx,j,i+offx_flx,icrm))*rhoi);
      field(ind_field,k,j,i+offx_s,icrm)=field(ind_field,k,j,i+offx_s,icrm) + dfdt(k,j,i,icrm);
    });
  }
}


