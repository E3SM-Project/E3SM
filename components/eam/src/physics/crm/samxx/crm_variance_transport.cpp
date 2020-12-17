#include "crm_variance_transport.h"

//==============================================================================
//==============================================================================

// void VT_filter() {
//   // local variables
//   integer, parameter :: lensav = nx+15 ! must be at least N + INT(LOG(REAL(N))) + 4.
//   real(crm_rknd), dimension(nx)    :: fft_out   ! for FFT input and output
//   real(crm_rknd), dimension(nx)    :: work      ! work array
//   real(crm_rknd), dimension(lensav):: wsave     ! prime factors of N and certain trig values used in rfft1f
//   ! real(crm_rknd), dimension(nx)    :: wave_num    ! only for debugging
//   integer :: i, j, k, icrm   ! loop iterators
//   integer :: ier             ! FFT error return code
//   //----------------------------------------------------------------------------
//   // initialization for FFT
//   // call rfft1i(nx,wsave,lensav,ier)
//   // if(ier /= 0) write(0,*) 'ERROR: rfftmi(): VT_filter - FFT initialization error ',ier

//   // do k = 1,nzm
//   //   do j = 1,ny
//   //     do icrm = 1,ncrms
//   parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
          
//           // initialize FFT input
//           parallel_for( SimpleBounds<1>(nx) , YAKL_LAMBDA (int k, int j, int icrm) {
//              fft_out(i) = f_in(icrm,i,j,k)
//           end do

//           // do the forward transform
//           call rfft1f( nx, 1, fft_out(:), nx, wsave, lensav, work(:), nx, ier )
//           if (ier/=0) write(0,*) 'ERROR: rfftmf(): VT_filter - forward FFT error ',ier

//           // filter out high frequencies
//           fft_out(2*(filter_wn_max+1):) = 0

//           // transform back
//           call rfft1b( nx, 1, fft_out(:), nx, wsave, lensav, work(:), nx, ier )
//           if(ier /= 0) write(0,*) 'ERROR: rfftmb(): VT_filter - backward FFT error ',ier

//           // copy to output
//           do i = 1,nx
//              f_out(icrm,i,j,k) = fft_out(i)
//           end do

//   });

// }

//==============================================================================
//==============================================================================

void VT_diagnose() {
  auto &t            = :: t;
  auto &qv           = :: qv;
  auto &qcl          = :: qcl;
  auto &qci          = :: qci;
  auto &factor_xy    = :: factor_xy;
  auto &t_vt_pert   = :: t_vt_pert;
  auto &q_vt_pert   = :: q_vt_pert;
  auto &t_vt        = :: t_vt;
  auto &q_vt        = :: q_vt;
  auto &ncrms        = :: ncrms;

  // local variables
  real2d t_mean("t_mean", nzm, ncrms);
  real2d q_mean("q_mean", nzm, ncrms);
  // real4d tmp_t("tmp_t", nzm, ny, nx, ncrms);
  // real4d tmp_q("tmp_q", nzm, ny, nx, ncrms);

  //----------------------------------------------------------------------------
  // calculate horizontal mean
  //----------------------------------------------------------------------------
  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      t_mean(k,icrm) = 0.0;
      q_mean(k,icrm) = 0.0;
      t_vt(k,icrm)  = 0.0;
      q_vt(k,icrm)  = 0.0;
  });

  // do k = 1,nzm
  //  do j = 1,ny
  //    do i = 1,nx
  //      do icrm = 1,ncrms
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real t_tmp = t(k,j+offy_s,i+offx_s,icrm);
    real q_tmp = qv(k,j,i,icrm) + qcl(k,j,i,icrm) + qci(k,j,i,icrm) ;
    yakl::atomicAdd( t_mean(k,icrm) , t_tmp);
    yakl::atomicAdd( q_mean(k,icrm) , q_tmp);
  });

  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      t_mean(k,icrm) = t_mean(k,icrm) * factor_xy ;
      q_mean(k,icrm) = q_mean(k,icrm) * factor_xy ;
  });

  //----------------------------------------------------------------------------
  // calculate anomalies
  //----------------------------------------------------------------------------
  // if (filter_wn_max>0) {
  // // use filtered state for anomalies

  //   // do k = 1,nzm
  //   //   do j = 1,ny
  //   //     do i = 1,nx
  //   //       do icrm = 1,ncrms
  //   parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
  //     tmp_t(k,j,i,icrm) = t(k,j+offy_s,i+offx_s,icrm)
  //     tmp_q(k,j,i,icrm) = qv(k,j,i,icrm) + qcl(k,j,i,icrm) + qci(k,j,i,icrm)
  //     tmp_t(k,j,i,icrm) = tmp_t(k,j,i,icrm) - t_mean(k,icrm)
  //     tmp_q(k,j,i,icrm) = tmp_q(k,j,i,icrm) - q_mean(k,icrm)
  //   });

  //   VT_filter( tmp_t, t_vt_pert )
  //   VT_filter( tmp_q, q_vt_pert )

  // else { 
  // // use total variance

    // do k = 1,nzm
    //   do j = 1,ny
    //     do i = 1,nx
    //       do icrm = 1,ncrms
    parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      t_vt_pert(k,j,i,icrm) = t(k,j+offy_s,i+offx_s,icrm) - t_mean(k,icrm);
      q_vt_pert(k,j,i,icrm) = qv(k,j,i,icrm) + qcl(k,j,i,icrm) + qci(k,j,i,icrm) - q_mean(k,icrm);
    });
    
  // }

  //----------------------------------------------------------------------------
  // calculate variance
  //----------------------------------------------------------------------------

  // do k = 1,nzm
  //   do j = 1,ny
  //     do i = 1,nx
  //       do icrm = 1,ncrms
  parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    real t_tmp = t_vt_pert(k,j,i,icrm) * t_vt_pert(k,j,i,icrm) ;
    real q_tmp = q_vt_pert(k,j,i,icrm) * q_vt_pert(k,j,i,icrm) ;
    yakl::atomicAdd( t_vt(k,icrm) , t_tmp);
    yakl::atomicAdd( q_vt(k,icrm) , q_tmp);
  });


  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    t_vt(k,icrm) = t_vt(k,icrm) * factor_xy ;
    q_vt(k,icrm) = q_vt(k,icrm) * factor_xy ;
  });

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
}

//==============================================================================
//==============================================================================

void VT_forcing() {
  auto &t            = :: t;
  auto &micro_field  = :: micro_field;
  auto &t_vt_tend   = :: t_vt_tend;
  auto &q_vt_tend   = :: q_vt_tend;
  auto &t_vt_pert   = :: t_vt_pert;
  auto &q_vt_pert   = :: q_vt_pert;
  auto &t_vt        = :: t_vt;
  auto &q_vt        = :: q_vt;
  auto &ncrms        = :: ncrms;
  auto &dtn          = :: dtn;

  // local variables
  real2d t_pert_scale("t_pert_scale", nzm, ncrms);
  real2d q_pert_scale("q_pert_scale", nzm, ncrms);
  real4d ttend_loc("ttend_loc", nzm, ny, nx, ncrms);
  real4d qtend_loc("qtend_loc", nzm, ny, nx, ncrms);

  int idx_qt = index_water_vapor;

  real pert_scale_min = 1.0 - 0.1;
  real pert_scale_max = 1.0 + 0.1;

  //----------------------------------------------------------------------------
  // calculate scaling factor for local perturbations
  //----------------------------------------------------------------------------
  // do k = 1,nzm
  //   do icrm = 1,ncrms
  parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
    // initialize scaling factors to 1.0
    t_pert_scale(k,icrm) = 1.0;
    q_pert_scale(k,icrm) = 1.0;
    real tmp_t_scale = -1;
    real tmp_q_scale = -1;
    // set scaling factors as long as there are perturbations to scale
    if (t_vt(k,icrm)>0.0) { tmp_t_scale = 1.0 + dtn * t_vt_tend(k,icrm) / t_vt(k,icrm); }
    if (q_vt(k,icrm)>0.0) { tmp_q_scale = 1.0 + dtn * q_vt_tend(k,icrm) / q_vt(k,icrm); }
    if (tmp_t_scale>0.0) { t_pert_scale(k,icrm) = sqrt( tmp_t_scale ); }
    if (tmp_q_scale>0.0) { q_pert_scale(k,icrm) = sqrt( tmp_q_scale ); }
    // enforce minimum scaling
    t_pert_scale(k,icrm) = max( t_pert_scale(k,icrm), pert_scale_min );
    q_pert_scale(k,icrm) = max( q_pert_scale(k,icrm), pert_scale_min );
    // enforce maximum scaling
    t_pert_scale(k,icrm) = min( t_pert_scale(k,icrm), pert_scale_max );
    q_pert_scale(k,icrm) = min( q_pert_scale(k,icrm), pert_scale_max );
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
  });

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
}

//==============================================================================
//==============================================================================