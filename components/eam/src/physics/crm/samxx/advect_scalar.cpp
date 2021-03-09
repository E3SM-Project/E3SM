#include "advect_scalar.h"

void advect_scalar(real4d &f, real2d &fadv, real2d &flux) {
  auto &ncrms = ::ncrms;

  real4d f0("f0", nzm, dimy_s, dimx_s, ncrms);

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  if (docolumn) {

    parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(k,icrm) = 0.0;
    });

  } else {

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<dimy_s; j++) {
    //     for (int i=0; i<dimx_s; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      f0(k,j,i,icrm) = f(k,j,i,icrm);
    });

    if (RUN3D) {
      advect_scalar3D(f,flux);
    } else {
      advect_scalar2D(f,flux);
    }

    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      fadv(k,icrm)=0.0;
    });
    
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      real tmp = f(k,j+offy_s,i+offx_s,icrm)-f0(k,j+offy_s,i+offx_s,icrm);
      yakl::atomicAdd(fadv(k,icrm),tmp);
    });

  }  

}

void advect_scalar(real5d &f, int ind_f, real2d &fadv, real2d &flux) {
  auto &ncrms         = :: ncrms;

  real4d f0("f0", nzm, dimy_s, dimx_s, ncrms);

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  if(docolumn) {

    parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(k,icrm) = 0.0;
    });

  } else {

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<dimy_s; j++) {
    //     for (int i=0; i<dimx_s; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      f0(k,j,i,icrm) = f(ind_f,k,j,i,icrm);
    });

    if(RUN3D) {
      advect_scalar3D(f,ind_f,flux);
    } else {
      advect_scalar2D(f,ind_f,flux);
    }

    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      fadv(k,icrm)=0.0;
    });
    
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      real tmp = f(ind_f,k,j+offy_s,i+offx_s,icrm)-f0(k,j+offy_s,i+offx_s,icrm);
      yakl::atomicAdd(fadv(k,icrm),tmp);
    });

  }  

}

void advect_scalar(real5d &f, int ind_f, real3d &fadv, int ind_fadv, real3d &flux, int ind_flux) {
  auto &ncrms         = :: ncrms;

  real4d f0("f0", nzm, dimy_s, dimx_s, ncrms);

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  if (docolumn) {
    parallel_for( SimpleBounds<2>(nz,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      flux(ind_flux,k,icrm) = 0.0;
    });

  } else {

    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<dimy_s; j++) {
    //     for (int i=0; i<dimx_s; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      f0(k,j,i,icrm) = f(ind_f,k,j,i,icrm);
    });

    if(RUN3D) {
      advect_scalar3D(f,ind_f,flux,ind_flux);
    } else {
      advect_scalar2D(f,ind_f,flux,ind_flux);
    }

    // for (int k=0; k<nzm; k++) {
    //  for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      fadv(ind_fadv,k,icrm)=0.0;
    });
    
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      real tmp = f(ind_f,k,j+offy_s,i+offx_s,icrm)-f0(k,j+offy_s,i+offx_s,icrm);
      yakl::atomicAdd(fadv(ind_fadv,k,icrm),tmp);
    });

  }  

}
