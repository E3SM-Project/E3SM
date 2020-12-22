#include "test_fft_yakl.h"

//==============================================================================
//==============================================================================

void VT_filter(real4d &f_in, real4d &f_out) {
  // local variables
  int nx2 = nx+2;
  int ny2 = ny+2*YES3D;
  real4d fft_out ("fft_out" , nzm, ny2, nx2, ncrms);

  int constexpr fftySize = ny > 4 ? ny : 4;

  int nwx = nx2-(filter_wn_max+1)*2;
  int nwy = ny2-(filter_wn_max+1)*2;
  
  yakl::FFT<nx> fftx;
  yakl::FFT<fftySize> ffty;

  //----------------------------------------------------------------------------
  // Forward Fourier transform
  
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
    real ftmp[nx+2]; real tmp [nx];
    for (int i=0; i<nx ; i++) { ftmp[i] = f_in(k,j,i,icrm); }
    fftx.forward(ftmp, tmp);
    for (int i=0; i<nx2; i++) { fft_out(k,j,i,icrm) = ftmp[i]; }
  });

  if (RUN3D) {
    // for (int k=0; k<nzm; k++) {
    //   for (int i=0; j<nx+1; i++) {
    //     for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      real ftmp[ny+2]; real tmp [ny];
      for (int j=0; j<ny ; j++) { ftmp[j] = fft_out(k,j,i,icrm); }
      ffty.forward(ftmp, tmp);
      for (int j=0; j<ny2; j++) { fft_out(k,j,i,icrm) = ftmp[j]; }
    });
  }

  //----------------------------------------------------------------------------
  // Zero out the higher modes

  if (not RUN3D) {
    parallel_for( SimpleBounds<4>(nzm,ny2,nx2/2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      printf("  i:%2d  j:%2d  k:%2d    fft_out: %+12.6f  %+12.6f \n",i*2, j, k, fft_out(k,j,i*2,icrm), fft_out(k,j,i*2+1,icrm) );
    });
  }

  if (RUN3D) {
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<nwy; j++) {
    //     for (int i=0; i<nwx+1; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<4>(nzm,nwy,nwx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      int ii = i + 2*(filter_wn_max+1) ;
      int jj = j + 2*(filter_wn_max+1) ;
      fft_out(k,jj,ii,icrm) = 0.0;
    });
  } else {
    // for (int k=0; k<nzm; k++) {
    //   for (int i=0; i<nwx+1; i++) {
    //     for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nwx,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      int ii = i + 2*(filter_wn_max+1) ;
      fft_out(k,0,ii,icrm) = 0.0;
    });
  }

  //----------------------------------------------------------------------------
  // Backward Fourier transform

  if (RUN3D) {
    // for (int k=0; k<nzm; k++) {
    //   for (int i=0; i<nx+1; i++) {
    //     for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      real ftmp[ny+2]; real tmp [ny];
      for(int j=0; j<ny+2; j++) { ftmp[j] = fft_out(k,j,i,icrm); }
      ffty.inverse(ftmp,tmp);
      for(int j=0; j<ny  ; j++) { fft_out(k,j,i,icrm) = ftmp[j]; } 
    });
  }

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; i<ny; i++) {
  //     for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
    real ftmp[nx+2]; real tmp [nx];
    for(int i=0; i<nx+2; i++) { ftmp[i] = fft_out(k,j,i,icrm); }
    fftx.inverse(ftmp,tmp);
    for(int i=0; i<nx  ; i++) { f_out(k,j,i,icrm) = ftmp[i]; }
  });


}

//==============================================================================
//==============================================================================

int main() {
  
  yakl::init();

  real4d data1("data" , nzm, ny, nx, ncrms);
  real4d data2("data" , nzm, ny, nx, ncrms);

  real constexpr pi = 3.14159265359;

  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    data2(k,j,i,icrm) = 0;
    // data1(k,j,i,icrm) = k*100.0 + i;
    real angle = (2*pi)*i/nx; real phase = (2*pi)*0.4;
    data1(k,j,i,icrm) = 0;

    // data1(k,j,i,icrm) = data1(k,j,i,icrm) + 8*cos(1*angle) ;
    // data1(k,j,i,icrm) = data1(k,j,i,icrm) + 2*sin(1*angle) ;
    // data1(k,j,i,icrm) = data1(k,j,i,icrm) + 3*sin(2*angle-phase) ; 
    // data1(k,j,i,icrm) = data1(k,j,i,icrm) + 10*cos(6*angle) ;
    // if (ny>1) {
    //   angle = (2*pi)*(j-1)/nx;
    //   data1(k,j,i,icrm) = data1(k,j,i,icrm) + 5*cos(6*angle) ;
    //   // data1(k,j,i,icrm) = data1(k,j,i,icrm) + 10*cos(6*angle) ;
    // }

    for(int w=1; w<nx/2+1; w++) { 
      data1(k,j,i,icrm) = data1(k,j,i,icrm) + 10*cos(w*angle) ;
      data1(k,j,i,icrm) = data1(k,j,i,icrm) + 10*sin(w*angle) ;
    }

  });

  VT_filter(data1,data2);

  real max_diff = 0.0;
  
  printf("fft input\n");
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    if (i==0 && k>0) { printf("\n"); }
    real diff = data1(k,j,i,icrm) - data2(k,j,i,icrm) ; 
    printf("  i:%2d  j:%2d  k:%2d    data1: %+12.4e   data2: %+12.4e    diff: %+12.4e \n",
            i, j, k, data1(k,j,i,icrm), data2(k,j,i,icrm), diff );
    max_diff = max(max_diff,diff);
  });
  printf("\n");


  printf("  max diff: %+12.4e \n\n",max_diff);

  // print coeffecients again
  if (not RUN3D) { VT_filter(data2,data2); }

  return 0;
}