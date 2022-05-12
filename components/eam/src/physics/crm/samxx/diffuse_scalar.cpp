#include "diffuse_scalar.h"


void diffuse_scalar(real5d &tkh, int ind_tkh, real4d &f, real3d &fluxb, real3d &fluxt, real2d &fdiff, real2d &flux) {
  YAKL_SCOPE( ncrms , ::ncrms );
  real4d df("df", nzm, dimy_s, dimx_s, ncrms);
  
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<dimy_s; j++) {
  //     for (int i=0; i<dimx_s; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    df(k,j,i,icrm) = f(k,j,i,icrm);
  });

  if (RUN3D) {
    diffuse_scalar3D(f,fluxb,fluxt,tkh,ind_tkh,flux);
  } else {
    diffuse_scalar2D(f,fluxb,fluxt,tkh,ind_tkh,flux);
  }

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    fdiff(k,icrm) = 0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real tmp = f(k,j+offy_s,i+offx_s,icrm)-df(k,j+offy_s,i+offx_s,icrm);
    yakl::atomicAdd(fdiff(k,icrm),tmp);
  });
}

void diffuse_scalar(real5d &tkh, int ind_tkh, real5d &f, int ind_f, real3d &fluxb,
                    real3d &fluxt, real2d &fdiff, real2d &flux) {
  YAKL_SCOPE( ncrms , ::ncrms );
  real4d df("df", nzm, dimy_s, dimx_s, ncrms);
  
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<dimy_s; j++) {
  //     for (int i=0; i<dimx_s; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    df(k,j,i,icrm) = f(ind_f,k,j,i,icrm);
  });

  if (RUN3D) {
    diffuse_scalar3D(f,ind_f,fluxb,fluxt,tkh,ind_tkh,flux);
  } else {
    diffuse_scalar2D(f,ind_f,fluxb,fluxt,tkh,ind_tkh,flux);
  }

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    fdiff(k,icrm) = 0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real tmp = f(ind_f,k,j+offy_s,i+offx_s,icrm)-df(k,j+offy_s,i+offx_s,icrm);
    yakl::atomicAdd(fdiff(k,icrm),tmp);
  });
}

void diffuse_scalar(real5d &tkh, int ind_tkh, real5d &f, int ind_f, real4d &fluxb, int ind_fluxb,
                    real4d &fluxt, int ind_fluxt, real3d &fdiff, int ind_fdiff, real3d &flux, int ind_flux) {
  YAKL_SCOPE( ncrms , ::ncrms );
  real4d df("df", nzm, dimy_s, dimx_s, ncrms);
  
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<dimy_s; j++) {
  //     for (int i=0; i<dimx_s; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,dimy_s,dimx_s,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    df(k,j,i,icrm) = f(ind_f,k,j,i,icrm);
  });

  if (RUN3D) {
    diffuse_scalar3D(f,ind_f,fluxb,ind_fluxb,fluxt,ind_fluxt,tkh,ind_tkh,flux,ind_flux);
  } else {
    diffuse_scalar2D(f,ind_f,fluxb,ind_fluxb,fluxt,ind_fluxt,tkh,ind_tkh,flux,ind_flux);
  }

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    fdiff(ind_fdiff,k,icrm) = 0.0;
  });

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real tmp = f(ind_f,k,j+offy_s,i+offx_s,icrm)-df(k,j+offy_s,i+offx_s,icrm);
    yakl::atomicAdd(fdiff(ind_fdiff,k,icrm),tmp);
  });
}

