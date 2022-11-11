#include "crm_variance_transport.h"

//==============================================================================
//==============================================================================

void VT_filter(int filter_wn_max, real4d &f_in, real4d &f_out) {
  // local variables
  int nx2 = nx+2;
  int ny2 = ny+2*YES3D;
  real4d fft_out ("fft_out" , nzm, ny2, nx2, ncrms);

  int constexpr fftySize = ny > 4 ? ny : 4;
  
  int nwx = nx2-(filter_wn_max+1)*2;
  int nwy = ny2-(filter_wn_max+1)*2;
  
  yakl::RealFFT1D<nx> fftx;
  yakl::RealFFT1D<fftySize> ffty;
  fftx.init();
  ffty.init();

  //----------------------------------------------------------------------------
  // Forward Fourier transform
  
  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
    SArray<real,1,nx+2> ftmp;
    for (int i=0; i<nx ; i++) { ftmp(i) = f_in(k,j,i,icrm); }
    fftx.forward(ftmp, fftx.trig, yakl::FFT_SCALE_ECMWF);
    for (int i=0; i<nx2; i++) { fft_out(k,j,i,icrm) = ftmp(i); }
  });

  if (RUN3D) {
    // for (int k=0; k<nzm; k++) {
    //   for (int i=0; j<nx+1; i++) {
    //     for (int icrm=0; icrm<ncrms; icrm++) {
    parallel_for( SimpleBounds<3>(nzm,nx+1,ncrms) , YAKL_LAMBDA (int k, int i, int icrm) {
      SArray<real,1,ny+2> ftmp;
      for (int j=0; j<ny ; j++) { ftmp(j) = fft_out(k,j,i,icrm); }
      ffty.forward(ftmp, ffty.trig, yakl::FFT_SCALE_ECMWF);
      for (int j=0; j<ny2; j++) { fft_out(k,j,i,icrm) = ftmp(j); }
    });
  }

  //----------------------------------------------------------------------------
  // Zero out the higher modes

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
      SArray<real,1,ny+2> ftmp;
      for(int j=0; j<ny+2; j++) { ftmp(j) = fft_out(k,j,i,icrm); }
      ffty.inverse(ftmp, ffty.trig, yakl::FFT_SCALE_ECMWF);
      for(int j=0; j<ny  ; j++) { fft_out(k,j,i,icrm) = ftmp(j); }
    });
  }

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; i<ny; i++) {
  //     for (int icrm=0; icrm<ncrms; icrm++) {
  parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
    SArray<real,1,nx+2> ftmp;
    for(int i=0; i<nx+2; i++) { ftmp(i) = fft_out(k,j,i,icrm); }
    fftx.inverse(ftmp, fftx.trig, yakl::FFT_SCALE_ECMWF);
    for(int i=0; i<nx  ; i++) { f_out(k,j,i,icrm) = ftmp(i); }
  });

}

//==============================================================================
//==============================================================================

void VT_diagnose() {
  YAKL_SCOPE( t             , :: t);
  YAKL_SCOPE( micro_field   , :: micro_field);
  YAKL_SCOPE( factor_xy     , :: factor_xy);
  YAKL_SCOPE( t_vt_pert     , :: t_vt_pert);
  YAKL_SCOPE( q_vt_pert     , :: q_vt_pert);
  YAKL_SCOPE( t_vt          , :: t_vt);
  YAKL_SCOPE( q_vt          , :: q_vt);
  YAKL_SCOPE( ncrms         , :: ncrms);
  YAKL_SCOPE( u             , :: u);
  YAKL_SCOPE( u_vt_pert     , :: u_vt_pert);
  YAKL_SCOPE( u_vt          , :: u_vt);

  // local variables
  real2d t_mean("t_mean", nzm, ncrms);
  real2d q_mean("q_mean", nzm, ncrms);
  real2d u_mean("u_mean", nzm, ncrms);

  int idx_qt = index_water_vapor;

  //----------------------------------------------------------------------------
  // calculate horizontal mean
  //----------------------------------------------------------------------------
  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      t_mean(k,icrm) = 0.0;
      q_mean(k,icrm) = 0.0;
      u_mean(k,icrm) = 0.0;
      t_vt(k,icrm) = 0.0;
      q_vt(k,icrm) = 0.0;
      u_vt(k,icrm) = 0.0;
  });

  // do k = 1,nzm
  //  do j = 1,ny
  //    do i = 1,nx
  //      do icrm = 1,ncrms
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    yakl::atomicAdd( t_mean(k,icrm) , t(k,j+offy_s,i+offx_s,icrm) );
    yakl::atomicAdd( q_mean(k,icrm) , micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) );
    yakl::atomicAdd( u_mean(k,icrm) , u(k,j+offy_u,i+offx_u,icrm) );
  });

  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      t_mean(k,icrm) = t_mean(k,icrm) * factor_xy ;
      q_mean(k,icrm) = q_mean(k,icrm) * factor_xy ;
      u_mean(k,icrm) = u_mean(k,icrm) * factor_xy ;
  });

  //----------------------------------------------------------------------------
  // calculate fluctuations - either from horz mean or with a low-pass filter
  //----------------------------------------------------------------------------
  if (VT_wn_max>0) { // use filtered state for fluctuations
  

    real4d tmp_t("tmp_t", nzm, ny, nx, ncrms);
    real4d tmp_q("tmp_q", nzm, ny, nx, ncrms);
    real4d tmp_u("tmp_u", nzm, ny, nx, ncrms);

    // do k = 1,nzm
    //   do j = 1,ny
    //     do i = 1,nx
    //       do icrm = 1,ncrms
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      tmp_t(k,j,i,icrm) = t(k,j+offy_s,i+offx_s,icrm);
      tmp_q(k,j,i,icrm) = micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm);
      tmp_t(k,j,i,icrm) = tmp_t(k,j,i,icrm) - t_mean(k,icrm);
      tmp_q(k,j,i,icrm) = tmp_q(k,j,i,icrm) - q_mean(k,icrm);
      tmp_u(k,j,i,icrm) = u(k,j+offy_u,i+offx_u,icrm);
      tmp_u(k,j,i,icrm) = tmp_u(k,j,i,icrm) - u_mean(k,icrm);
    });

    VT_filter( VT_wn_max, tmp_t, t_vt_pert );
    VT_filter( VT_wn_max, tmp_q, q_vt_pert );
    VT_filter( VT_wn_max, tmp_u, u_vt_pert );

  } else { // use total variance

    // do k = 1,nzm
    //   do j = 1,ny
    //     do i = 1,nx
    //       do icrm = 1,ncrms
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      t_vt_pert(k,j,i,icrm) = t(k,j+offy_s,i+offx_s,icrm) - t_mean(k,icrm);
      q_vt_pert(k,j,i,icrm) = micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) - q_mean(k,icrm);
      u_vt_pert(k,j,i,icrm) = u(k,j+offy_u,i+offx_u,icrm) - u_mean(k,icrm);
    });
    
  }

  //----------------------------------------------------------------------------
  // calculate variance
  //----------------------------------------------------------------------------

  // do k = 1,nzm
  //   do j = 1,ny
  //     do i = 1,nx
  //       do icrm = 1,ncrms
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    yakl::atomicAdd( t_vt(k,icrm) , t_vt_pert(k,j,i,icrm) * t_vt_pert(k,j,i,icrm) );
    yakl::atomicAdd( q_vt(k,icrm) , q_vt_pert(k,j,i,icrm) * q_vt_pert(k,j,i,icrm) );
    yakl::atomicAdd( u_vt(k,icrm) , u_vt_pert(k,j,i,icrm) * u_vt_pert(k,j,i,icrm) );
  });


  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    t_vt(k,icrm) = t_vt(k,icrm) * factor_xy ;
    q_vt(k,icrm) = q_vt(k,icrm) * factor_xy ;
    u_vt(k,icrm) = u_vt(k,icrm) * factor_xy ;
  });

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
}

//==============================================================================
//==============================================================================

void VT_forcing() {
  YAKL_SCOPE( t            , :: t);
  YAKL_SCOPE( micro_field  , :: micro_field);
  YAKL_SCOPE( t_vt_tend    , :: t_vt_tend);
  YAKL_SCOPE( q_vt_tend    , :: q_vt_tend);
  YAKL_SCOPE( t_vt_pert    , :: t_vt_pert);
  YAKL_SCOPE( q_vt_pert    , :: q_vt_pert);
  YAKL_SCOPE( t_vt         , :: t_vt);
  YAKL_SCOPE( q_vt         , :: q_vt);
  YAKL_SCOPE( ncrms        , :: ncrms);
  YAKL_SCOPE( dtn          , :: dtn);
  YAKL_SCOPE( u            , :: u);
  YAKL_SCOPE( u_vt_pert    , :: u_vt_pert);
  YAKL_SCOPE( u_vt         , :: u_vt);
  YAKL_SCOPE( u_vt_tend    , :: u_vt_tend);

  // local variables
  real2d t_pert_scale("t_pert_scale", nzm, ncrms);
  real2d q_pert_scale("q_pert_scale", nzm, ncrms);
  real2d u_pert_scale("u_pert_scale", nzm, ncrms);

  int idx_qt = index_water_vapor;

  // min and max perturbation scaling values are used to limit the 
  // large-scale forcing from variance transport. This is meant to 
  // protect against creating unstable situations, although 
  // problematic scenarios were extremely rare in testing.
  // A scaling limit of +/- 10% was found to be adequate.
  real constexpr pert_scale_min = 1.0 - 0.1;
  real constexpr pert_scale_max = 1.0 + 0.1;

  //----------------------------------------------------------------------------
  // calculate scaling factor for local perturbations
  //----------------------------------------------------------------------------
  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    // initialize scaling factors to 1.0
    t_pert_scale(k,icrm) = 1.0;
    q_pert_scale(k,icrm) = 1.0;
    u_pert_scale(k,icrm) = 1.0;
    real tmp_t_scale = -1.0;
    real tmp_q_scale = -1.0;
    real tmp_u_scale = -1.0;
    // set scaling factors as long as there are perturbations to scale
    if (t_vt(k,icrm)>0.0) { tmp_t_scale = 1.0 + dtn * t_vt_tend(k,icrm) / t_vt(k,icrm); }
    if (q_vt(k,icrm)>0.0) { tmp_q_scale = 1.0 + dtn * q_vt_tend(k,icrm) / q_vt(k,icrm); }
    if (u_vt(k,icrm)>0.0) { tmp_u_scale = 1.0 + dtn * u_vt_tend(k,icrm) / u_vt(k,icrm); }
    if (tmp_t_scale>0.0) { t_pert_scale(k,icrm) = sqrt( tmp_t_scale ); }
    if (tmp_q_scale>0.0) { q_pert_scale(k,icrm) = sqrt( tmp_q_scale ); }
    if (tmp_u_scale>0.0) { u_pert_scale(k,icrm) = sqrt( tmp_u_scale ); }
    // enforce minimum scaling
    t_pert_scale(k,icrm) = max( t_pert_scale(k,icrm), pert_scale_min );
    q_pert_scale(k,icrm) = max( q_pert_scale(k,icrm), pert_scale_min );
    u_pert_scale(k,icrm) = max( u_pert_scale(k,icrm), pert_scale_min );
    // enforce maximum scaling
    t_pert_scale(k,icrm) = min( t_pert_scale(k,icrm), pert_scale_max );
    q_pert_scale(k,icrm) = min( q_pert_scale(k,icrm), pert_scale_max );
    u_pert_scale(k,icrm) = min( u_pert_scale(k,icrm), pert_scale_max );
  });

  //----------------------------------------------------------------------------
  // apply variance forcing tendency
  //----------------------------------------------------------------------------
  // do k = 1,nzm
  //   do j = 1,ny
  //     do i = 1,nx
  //       do icrm = 1,ncrms
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real ttend_loc = ( t_pert_scale(k,icrm) * t_vt_pert(k,j,i,icrm) - t_vt_pert(k,j,i,icrm) ) / dtn;
    real qtend_loc = ( q_pert_scale(k,icrm) * q_vt_pert(k,j,i,icrm) - q_vt_pert(k,j,i,icrm) ) / dtn;
    t(k,j+offy_s,i+offx_s,icrm)                  = t(k,j+offy_s,i+offx_s,icrm)                  + ttend_loc * dtn;
    micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) = micro_field(idx_qt,k,j+offy_s,i+offx_s,icrm) + qtend_loc * dtn;
    real utend_loc = ( u_pert_scale(k,icrm) * u_vt_pert(k,j,i,icrm) - u_vt_pert(k,j,i,icrm) ) / dtn;
    u(k,j+offy_u,i+offx_u,icrm) = u(k,j+offy_u,i+offx_u,icrm) + utend_loc * dtn;
  });

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
}

//==============================================================================
//==============================================================================
