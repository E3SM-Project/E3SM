#include "diffuse_scalar3D.h"

void diffuse_scalar3D(real4d &field, real3d &fluxb, real3d &fluxt, real5d &tkh,
                      int ind_tkh, real2d &flux) {
  auto &dx     = ::dx;
  auto &dy     = ::dy;
  auto &rhow   = ::rhow;
  auto &adzw   = ::adzw;
  auto &adz    = ::adz;
  auto &dz     = ::dz; 
  auto &dtn    = ::dtn;
  auto &rho    = ::rho;
  auto &grdf_x = ::grdf_x;
  auto &grdf_y = ::grdf_y;
  auto &grdf_z = ::grdf_z;
  auto &ncrms  = ::ncrms;

  if (dosgs) {
    real4d flx_x("flx_x", nzm+1, ny+1, nx+1, ncrms);
    real4d flx_y("flx_y", nzm+1, ny+1, nx+1, ncrms);
    real4d flx_z("flx_z", nzm+1, ny+1, nx+1, ncrms);
    real4d dfdt("dfdt", nz, ny, nx, ncrms);

    int constexpr offx_flx = 1;
    int constexpr offy_flx = 1;
    int constexpr offz_flx = 1;

    real rdx2=1.0/(dx*dx);
    real rdy2=1.0/(dy*dy);
    real dxy=dx/dy;
    real dyx=dy/dx;

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      dfdt(k,j,i,icrm)=0.0;
    });

    //  Horizontal diffusion:
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+1; j++) {
    //     for (int i=0; i<nx+1; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+1,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (j >= 1) {
        int ic=i+1;
        real rdx5=0.5*rdx2 * grdf_x(k,icrm);
        real tkx=rdx5*(tkh(ind_tkh,k,j+offy_d-1,i+offx_d-1,icrm)+tkh(ind_tkh,k,j+offy_d-1,ic+offx_d-1,icrm));
        flx_x(k+offz_flx,j,i,icrm)=-tkx*(field(k,j+offy_s-1,ic+offx_s-1,icrm)-field(k,j+offy_s-1,i+offx_s-1,icrm));
      }
      if (i >= 1) {
        int jc=j+1;
        real rdy5=0.5*rdy2 * grdf_y(k,icrm);
        real tky=rdy5*(tkh(ind_tkh,k,j+offy_d-1,i+offx_d-1,icrm)+tkh(ind_tkh,k,jc+offy_d-1,i+offx_d-1,icrm));
        flx_y(k+offz_flx,j,i,icrm)=-tky*(field(k,jc+offy_s-1,i+offx_s-1,icrm)-field(k,j+offy_s-1,i+offx_s-1,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int ib=i-1;
      dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx_x(k+offz_flx,j+offy_flx,i +offx_flx,icrm)-
                                         flx_x(k+offz_flx,j+offy_flx,ib+offx_flx,icrm));
      int jb=j-1;
      dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx_y(k+offz_flx,j +offy_flx,i+offx_flx,icrm)-
                                         flx_y(k+offz_flx,jb+offy_flx,i+offx_flx,icrm));
    });

    //  Vertical diffusion:
    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(k,icrm) = 0.0;
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (k <= nzm-2) {
        int kc=k+1;
        real rhoi = rhow(kc,icrm)/adzw(kc,icrm);
        real rdz2 = 1.0/(dz(icrm)*dz(icrm));
        real rdz5 = 0.5*rdz2 * grdf_z(k,icrm);
        real tkz = rdz5*(tkh(ind_tkh,k,j+offy_d,i+offx_d,icrm)+tkh(ind_tkh,kc,j+offy_d,i+offx_d,icrm));
        flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm)=-tkz*(field(kc,j+offy_s,i+offx_s,icrm)-
                                                           field(k ,j+offy_s,i+offx_s,icrm))*rhoi;
        yakl::atomicAdd(flux(kc,icrm), flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm));
      }
      else if (k == nzm-1) {
        real tmp=1.0/adzw(nz-1,icrm);
        real rdz=1.0/dz(icrm);
        flx_z(0,j+offy_flx,i+offx_flx,icrm)=fluxb(j,i,icrm)*rdz*rhow(0,icrm);
        flx_z(nzm-1+offz_flx,j+offy_flx,i+offx_flx,icrm)=fluxt(j,i,icrm)*rdz*tmp*rhow(nz-1,icrm);
        yakl::atomicAdd(flux(0,icrm),flx_z(0,j+offy_flx,i+offx_flx,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kb=k-1;
      real rhoi = 1.0/(adz(k,icrm)*rho(k,icrm));
      dfdt(k,j,i,icrm)=dtn*(dfdt(k,j,i,icrm)-(flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm)-
                                              flx_z(kb+offz_flx,j+offy_flx,i+offx_flx,icrm))*rhoi);
      field(k,j+offy_s,i+offx_s,icrm)=field(k,j+offy_s,i+offx_s,icrm)+dfdt(k,j,i,icrm);
    });
  }
}


void diffuse_scalar3D(real5d &field, int ind_field, real3d &fluxb, real3d &fluxt, real5d &tkh,
                      int ind_tkh, real2d &flux) {
  auto &dx     = ::dx;
  auto &dy     = ::dy;
  auto &rhow   = ::rhow;
  auto &adzw   = ::adzw;
  auto &adz    = ::adz;
  auto &dz     = ::dz; 
  auto &dtn    = ::dtn;
  auto &rho    = ::rho;
  auto &grdf_x = ::grdf_x;
  auto &grdf_y = ::grdf_y;
  auto &grdf_z = ::grdf_z;
  auto &ncrms  = ::ncrms;
  
  if (dosgs) {
    real4d flx_x("flx_x", nzm+1, ny+1, nx+1, ncrms);
    real4d flx_y("flx_y", nzm+1, ny+1, nx+1, ncrms);
    real4d flx_z("flx_z", nzm+1, ny+1, nx+1, ncrms);
    real4d dfdt("dfdt", nz, ny, nx, ncrms);
    int constexpr offx_flx = 1;
    int constexpr offy_flx = 1;
    int constexpr offz_flx = 1;

    real rdx2=1.0/(dx*dx);
    real rdy2=1.0/(dy*dy);
    real dxy=dx/dy;
    real dyx=dy/dx;

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      dfdt(k,j,i,icrm)=0.0;
    });

    //  Horizontal diffusion:
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+1; j++) {
    //     for (int i=0; i<nx+1; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+1,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (j >= 1) {
        int ic=i+1;
        real rdx5=0.5*rdx2 * grdf_x(k,icrm);
        real tkx=rdx5*(tkh(ind_tkh,k,j+offy_d-1,i+offx_d-1,icrm)+tkh(ind_tkh,k,j+offy_d-1,ic+offx_d-1,icrm));
        flx_x(k+offz_flx,j,i,icrm)=-tkx*(field(ind_field,k,j+offy_s-1,ic+offx_s-1,icrm)-
                                         field(ind_field,k,j+offy_s-1,i+offx_s-1,icrm));
      }
      if (i >= 1) {
        int jc=j+1;
        real rdy5=0.5*rdy2 * grdf_y(k,icrm);
        real tky=rdy5*(tkh(ind_tkh,k,j+offy_d-1,i+offx_d-1,icrm)+tkh(ind_tkh,k,jc+offy_d-1,i+offx_d-1,icrm));
        flx_y(k+offz_flx,j,i,icrm)=-tky*(field(ind_field,k,jc+offy_s-1,i+offx_s-1,icrm)-
                                         field(ind_field,k,j+offy_s-1,i+offx_s-1,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int ib=i-1;
      dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx_x(k+offz_flx,j+offy_flx,i+offx_flx,icrm)-
                                         flx_x(k+offz_flx,j+offy_flx,ib+offx_flx,icrm));
      int jb=j-1;
      dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx_y(k+offz_flx,j+offy_flx,i+offx_flx,icrm)-
                                         flx_y(k+offz_flx,jb+offy_flx,i+offx_flx,icrm));
    });

    //  Vertical diffusion:
    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(k,icrm) = 0.0;
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (k <= nzm-2) {
        int kc=k+1;
        real rhoi = rhow(kc,icrm)/adzw(kc,icrm);
        real rdz2 = 1.0/(dz(icrm)*dz(icrm));
        real rdz5 = 0.5*rdz2 * grdf_z(k,icrm);
        real tkz = rdz5*(tkh(ind_tkh,k,j+offy_d,i+offx_d,icrm)+tkh(ind_tkh,kc,j+offy_d,i+offx_d,icrm));
        flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm)=-tkz*(field(ind_field,kc,j+offy_s,i+offx_s,icrm)-
                                                           field(ind_field,k,j+offy_s,i+offx_s,icrm))*rhoi;
        yakl::atomicAdd(flux(kc,icrm), flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm));
      } else if (k == nzm-1) {
        real tmp=1.0/adzw(nz-1,icrm);
        real rdz=1.0/dz(icrm);
        flx_z(0,j+offy_flx,i+offx_flx,icrm)=fluxb(j,i,icrm)*rdz*rhow(0,icrm);
        flx_z(nzm-1+offz_flx,j+offy_flx,i+offx_flx,icrm)=fluxt(j,i,icrm)*rdz*tmp*rhow(nz-1,icrm);
        yakl::atomicAdd(flux(0,icrm),flx_z(0,j+offy_flx,i+offx_flx,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kb=k-1;
      real rhoi = 1.0/(adz(k,icrm)*rho(k,icrm));
      dfdt(k,j,i,icrm)=dtn*(dfdt(k,j,i,icrm)-(flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm)-
                                              flx_z(kb+offz_flx,j+offy_flx,i+offx_flx,icrm))*rhoi);
      field(ind_field,k,j+offy_s,i+offx_s,icrm)=field(ind_field,k,j+offy_s,i+offx_s,icrm)+dfdt(k,j,i,icrm);
    });

  }
}

void diffuse_scalar3D(real5d &field, int ind_field, real4d &fluxb, int ind_fluxb, real4d &fluxt,
                      int ind_fluxt, real5d &tkh, int ind_tkh, real3d &flux, int ind_flux) {
  auto &dx     = ::dx;
  auto &dy     = ::dy;
  auto &rhow   = ::rhow;
  auto &adzw   = ::adzw;
  auto &adz    = ::adz;
  auto &dz     = ::dz; 
  auto &dtn    = ::dtn;
  auto &rho    = ::rho;
  auto &grdf_x = ::grdf_x;
  auto &grdf_y = ::grdf_y;
  auto &grdf_z = ::grdf_z;
  auto &ncrms  = ::ncrms;
  
  if (dosgs) {
    real4d flx_x("flx_x", nzm+1, ny+1, nx+1, ncrms);
    real4d flx_y("flx_y", nzm+1, ny+1, nx+1, ncrms);
    real4d flx_z("flx_z", nzm+1, ny+1, nx+1, ncrms);
    real4d dfdt("dfdt", nz, ny, nx, ncrms);

    int constexpr offx_flx = 1;
    int constexpr offy_flx = 1;
    int constexpr offz_flx = 1;

    real rdx2=1.0/(dx*dx);
    real rdy2=1.0/(dy*dy);
    real dxy=dx/dy;
    real dyx=dy/dx;

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      dfdt(k,j,i,icrm)=0.0;
    });

    //  Horizontal diffusion:
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny+1; j++) {
    //     for (int i=0; i<nx+1; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny+1,nx+1,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (j >= 1) {
        int ic=i+1;
        real rdx5=0.5*rdx2 * grdf_x(k,icrm);
        real tkx=rdx5*(tkh(ind_tkh,k,j+offy_d-1,i+offx_d-1,icrm)+tkh(ind_tkh,k,j+offy_d-1,ic+offx_d-1,icrm));
        flx_x(k+offz_flx,j,i,icrm)=-tkx*(field(ind_field,k,j+offy_s-1,ic+offx_s-1,icrm)-
                                         field(ind_field,k,j+offy_s-1,i+offx_s-1,icrm));
      }
      if (i >= 1) {
        int jc=j+1;
        real rdy5=0.5*rdy2 * grdf_y(k,icrm);
        real tky=rdy5*(tkh(ind_tkh,k,j+offy_d-1,i+offx_d-1,icrm)+tkh(ind_tkh,k,jc+offy_d-1,i+offx_d-1,icrm));
        flx_y(k+offz_flx,j,i,icrm)=-tky*(field(ind_field,k,jc+offy_s-1,i+offx_s-1,icrm)-
                                         field(ind_field,k,j+offy_s-1,i+offx_s-1,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int ib=i-1;
      dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx_x(k+offz_flx,j+offy_flx,i+offx_flx,icrm)-
                                         flx_x(k+offz_flx,j+offy_flx,ib+offx_flx,icrm));
      int jb=j-1;
      dfdt(k,j,i,icrm)=dfdt(k,j,i,icrm)-(flx_y(k+offz_flx,j+offy_flx,i+offx_flx,icrm)-
                                         flx_y(k+offz_flx,jb+offy_flx,i+offx_flx,icrm));
    });

    //  Vertical diffusion:
    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(ind_flux,k,icrm) = 0.0;
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (k <= nzm-2) {
        int kc=k+1;
        real rhoi = rhow(kc,icrm)/adzw(kc,icrm);
        real rdz2 = 1.0/(dz(icrm)*dz(icrm));
        real rdz5 = 0.5*rdz2 * grdf_z(k,icrm);
        real tkz = rdz5*(tkh(ind_tkh,k,j+offy_d,i+offx_d,icrm)+tkh(ind_tkh,kc,j+offy_d,i+offx_d,icrm));
        flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm)=-tkz*(field(ind_field,kc,j+offy_s,i+offx_s,icrm)-
                                                           field(ind_field,k,j+offy_s,i+offx_s,icrm))*rhoi;
        yakl::atomicAdd(flux(ind_flux,kc,icrm), flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm));
      } else if (k == nzm-1) {
        real tmp=1.0/adzw(nz-1,icrm);
        real rdz=1.0/dz(icrm);
        flx_z(0,j+offy_flx,i+offx_flx,icrm)=fluxb(ind_fluxb,j,i,icrm)*rdz*rhow(0,icrm);
        flx_z(nzm-1+offz_flx,j+offy_flx,i+offx_flx,icrm)=fluxt(ind_fluxt,j,i,icrm)*rdz*tmp*rhow(nz-1,icrm);
        yakl::atomicAdd(flux(ind_flux,0,icrm),flx_z(0,j+offy_flx,i+offx_flx,icrm));
      }
    });

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int kb=k-1;
      real rhoi = 1.0/(adz(k,icrm)*rho(k,icrm));
      dfdt(k,j,i,icrm)=dtn*(dfdt(k,j,i,icrm)-(flx_z(k+offz_flx,j+offy_flx,i+offx_flx,icrm)-
                                              flx_z(kb+offz_flx,j+offy_flx,i+offx_flx,icrm))*rhoi);
      field(ind_field,k,j+offy_s,i+offx_s,icrm)=field(ind_field,k,j+offy_s,i+offx_s,icrm)+dfdt(k,j,i,icrm);
    });

  }
}


