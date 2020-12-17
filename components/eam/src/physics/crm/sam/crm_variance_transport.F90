module variance_transport_mod
!-------------------------------------------------------------------------------
! This module supports the capability of transporting the CRM's internal 
! variance on the large-scale host model grid. This is done by calculating a 
! total or filtered variance and using this to couple to a mass-less tracer 
! on the GCM grid.
!
! The variance transport only applies to the prognostic scalars of the CRM. 
! Extending this to momentum may be possible, but it does not seem necessary.
! 
! Author: Walter Hannah - Lawrence Livermore National Lab
!-------------------------------------------------------------------------------
   use params_kind,  only: crm_rknd
   use grid,         only: nx,ny,nzm,dtn
   use params,       only: pi
   use openacc_utils

   implicit none

   public allocate_VT
   public deallocate_VT
   public VT_diagnose
   public VT_forcing

   real(crm_rknd), allocatable :: t_vt_tend(:,:)     ! forcing tendency of LSE amp per wavenumber
   real(crm_rknd), allocatable :: q_vt_tend(:,:)     ! forcing tendency of QT  amp per wavenumber

   real(crm_rknd), allocatable :: t_vt(:,:)          ! LSE variance tracer
   real(crm_rknd), allocatable :: q_vt(:,:)          ! QT  variance tracer

   real(crm_rknd), allocatable :: t_vt_pert(:,:,:,:) ! LSE perturbation from horizontal mean
   real(crm_rknd), allocatable :: q_vt_pert(:,:,:,:) ! QT  perturbation from horizontal mean

#if defined(MMF_VT_KMAX)
   integer, parameter :: filter_wn_max = MMF_VT_KMAX
#else
   integer, parameter :: filter_wn_max = 0
#endif

contains

!===============================================================================
!===============================================================================
subroutine allocate_VT(ncrms)
   !----------------------------------------------------------------------------
   ! Purpose: Allocate and initialize VT variables
   !----------------------------------------------------------------------------
   implicit none
   integer, intent(in) :: ncrms

   allocate( t_vt_tend( ncrms, nzm ) )
   allocate( q_vt_tend( ncrms, nzm ) )

   allocate( t_vt( ncrms, nzm ) )
   allocate( q_vt( ncrms, nzm ) )   

   t_vt_tend(:,:) = 0.0_crm_rknd
   q_vt_tend(:,:) = 0.0_crm_rknd

   t_vt(:,:)      = 0.0_crm_rknd
   q_vt(:,:)      = 0.0_crm_rknd

   allocate( t_vt_pert( ncrms, nx, ny, nzm ) )
   allocate( q_vt_pert( ncrms, nx, ny, nzm ) )
   t_vt_pert(:,:,:,:)  = 0.0_crm_rknd
   q_vt_pert(:,:,:,:)  = 0.0_crm_rknd

end subroutine allocate_VT

!===============================================================================
!===============================================================================
subroutine deallocate_VT()
   !----------------------------------------------------------------------------
   ! Purpose: Deallocate VT variables
   !----------------------------------------------------------------------------
   deallocate( t_vt_tend )
   deallocate( q_vt_tend )
   deallocate( t_vt )
   deallocate( q_vt )   
   deallocate( t_vt_pert )
   deallocate( q_vt_pert )

end subroutine deallocate_VT

!===============================================================================
!===============================================================================
subroutine VT_filter(ncrms,f_in,f_out)
   !----------------------------------------------------------------------------
   ! Purpose: use FFT to filter out high frequency modes
   !----------------------------------------------------------------------------
   use crmdims,   only: crm_dx
   use fftpack51D
   implicit none

   ! interface arguments
   integer, intent(in) :: ncrms
   real(crm_rknd), dimension(ncrms,nx,ny,nzm), intent(in ) :: f_in
   real(crm_rknd), dimension(ncrms,nx,ny,nzm), intent(out) :: f_out

   ! local variables
   integer, parameter :: lensav = nx+15 ! must be at least N + INT(LOG(REAL(N))) + 4.
   real(crm_rknd), dimension(nx)    :: fft_out   ! for FFT input and output
   real(crm_rknd), dimension(nx)    :: work      ! work array
   real(crm_rknd), dimension(lensav):: wsave     ! prime factors of N and certain trig values used in rfft1f
   ! real(crm_rknd), dimension(nx)    :: wave_num    ! only for debugging
   integer :: i, j, k, icrm   ! loop iterators
   integer :: ier             ! FFT error return code
   !----------------------------------------------------------------------------
   ! initialization for FFT
   call rfft1i(nx,wsave,lensav,ier)
   if(ier /= 0) write(0,*) 'ERROR: rfftmi(): VT_filter - FFT initialization error ',ier
   
   !$acc parallel loop collapse(3) async(asyncid)
   do k = 1,nzm
      do j = 1,ny
         do icrm = 1,ncrms
            
            ! initialize FFT input
            do i = 1,nx
               fft_out(i) = f_in(icrm,i,j,k)
            end do

            ! do the forward transform
            call rfft1f( nx, 1, fft_out(:), nx, wsave, lensav, work(:), nx, ier )
            if (ier/=0) write(0,*) 'ERROR: rfftmf(): VT_filter - forward FFT error ',ier

            ! filter out high frequencies
            fft_out(2*(filter_wn_max+1):) = 0

            ! transform back
            call rfft1b( nx, 1, fft_out(:), nx, wsave, lensav, work(:), nx, ier )
            if(ier /= 0) write(0,*) 'ERROR: rfftmb(): VT_filter - backward FFT error ',ier

            ! copy to output
            do i = 1,nx
               f_out(icrm,i,j,k) = fft_out(i)
            end do

         end do
      end do
   end do

end subroutine VT_filter

!===============================================================================
!===============================================================================
subroutine VT_diagnose(ncrms)
   !----------------------------------------------------------------------------
   ! Purpose: Diagnose amplitude for each wavenumber for variance transport
   !----------------------------------------------------------------------------
   use crmdims,   only: crm_dx
   use vars,      only: t,qv,qcl,qci
   implicit none

   ! interface arguments
   integer, intent(in) :: ncrms

   ! local variables
   real(crm_rknd), allocatable :: t_mean(:,:)
   real(crm_rknd), allocatable :: q_mean(:,:)
   real(crm_rknd), allocatable :: tmp_t(:,:,:,:)
   real(crm_rknd), allocatable :: tmp_q(:,:,:,:)
   real(crm_rknd) :: tmp
   real(crm_rknd) :: factor_xy
   integer :: i, j, k, icrm   ! loop iterators
   !----------------------------------------------------------------------------

   allocate( t_mean( ncrms, nzm ) )
   allocate( q_mean( ncrms, nzm ) )

   allocate( tmp_t( ncrms, nx, ny, nzm ) )
   allocate( tmp_q( ncrms, nx, ny, nzm ) )

   factor_xy = 1._crm_rknd/dble(nx*ny)

   !----------------------------------------------------------------------------
   ! calculate horizontal mean
   !----------------------------------------------------------------------------
   !$acc parallel loop collapse(2) async(asyncid)
   do k = 1,nzm
      do icrm = 1,ncrms
         t_mean(icrm,k) = 0.
         q_mean(icrm,k) = 0.
         t_vt(icrm,k) = 0.
         q_vt(icrm,k) = 0.
      end do
   end do

   !$acc parallel loop collapse(4) async(asyncid)
   do k = 1,nzm
      do j = 1,ny
         do i = 1,nx
            do icrm = 1,ncrms
               !$acc atomic update
               t_mean(icrm,k) = t_mean(icrm,k) + t(icrm,i,j,k)
               tmp = qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)
               !$acc atomic update
               q_mean(icrm,k) = q_mean(icrm,k) + tmp
            end do
         end do
      end do
   end do

   !$acc parallel loop collapse(2) async(asyncid)
   do k = 1,nzm
      do icrm = 1,ncrms
         t_mean(icrm,k) = t_mean(icrm,k) * factor_xy
         q_mean(icrm,k) = q_mean(icrm,k) * factor_xy
      end do
   end do

   !----------------------------------------------------------------------------
   ! calculate anomalies
   !----------------------------------------------------------------------------
   if (filter_wn_max>0) then ! use filtered state for anomalies

      !$acc parallel loop collapse(4) async(asyncid)
      do k = 1,nzm
         do j = 1,ny
            do i = 1,nx
               do icrm = 1,ncrms
                  tmp_t(icrm,i,j,k) = t(icrm,i,j,k) 
                  tmp_q(icrm,i,j,k) = qv(icrm,i,j,k) + qcl(icrm,i,j,k) + qci(icrm,i,j,k)
                  tmp_t(icrm,i,j,k) = tmp_t(icrm,i,j,k) - t_mean(icrm,k)
                  tmp_q(icrm,i,j,k) = tmp_q(icrm,i,j,k) - q_mean(icrm,k)
               end do
            end do
         end do
      end do

      call VT_filter( ncrms, tmp_t, t_vt_pert )
      call VT_filter( ncrms, tmp_q, q_vt_pert )

   else ! use total variance

      !$acc parallel loop collapse(4) async(asyncid)
      do k = 1,nzm
         do j = 1,ny
            do i = 1,nx
               do icrm = 1,ncrms
                  t_vt_pert(icrm,i,j,k) = t(icrm,i,j,k) - t_mean(icrm,k)
                  tmp = qv(icrm,i,j,k) + qcl(icrm,i,j,k) + qci(icrm,i,j,k)
                  q_vt_pert(icrm,i,j,k) = tmp           - q_mean(icrm,k)
               end do
            end do
         end do
      end do

   end if

   !----------------------------------------------------------------------------
   ! calculate variance
   !----------------------------------------------------------------------------
   !$acc parallel loop collapse(4) async(asyncid)
   do k = 1,nzm
      do j = 1,ny
         do i = 1,nx
            do icrm = 1,ncrms
               t_vt(icrm,k) = t_vt(icrm,k) + t_vt_pert(icrm,i,j,k) * t_vt_pert(icrm,i,j,k)
               q_vt(icrm,k) = q_vt(icrm,k) + q_vt_pert(icrm,i,j,k) * q_vt_pert(icrm,i,j,k)
            end do
         end do
      end do
   end do

   !$acc parallel loop collapse(2) async(asyncid)
   do k = 1,nzm
      do icrm = 1,ncrms
         t_vt(icrm,k) = t_vt(icrm,k) * factor_xy
         q_vt(icrm,k) = q_vt(icrm,k) * factor_xy
      end do
   end do

   !----------------------------------------------------------------------------
   ! clean up
   !----------------------------------------------------------------------------
   deallocate( t_mean )
   deallocate( q_mean )

   deallocate( tmp_t )
   deallocate( tmp_q )

end subroutine VT_diagnose

!===============================================================================
!===============================================================================
subroutine VT_forcing(ncrms)
   !----------------------------------------------------------------------------
   ! Purpose: Calculate forcing for variance injection and apply limiters
   !----------------------------------------------------------------------------
   use crmdims,      only: crm_dx
   use vars,         only: t
   use microphysics, only: micro_field, index_water_vapor
   implicit none

   ! interface arguments
   integer, intent(in) :: ncrms

   ! local variables
   real(crm_rknd), parameter   :: pert_scale_min = 1.0 - 0.1
   real(crm_rknd), parameter   :: pert_scale_max = 1.0 + 0.1
   real(crm_rknd), allocatable :: t_pert_scale(:,:)
   real(crm_rknd), allocatable :: q_pert_scale(:,:)
   real(crm_rknd) :: ttend_loc, qtend_loc
   real(crm_rknd) :: tmp_t_scale, tmp_q_scale
   integer :: i, j, k, icrm   ! loop iterators
   integer :: idx_qt
   !----------------------------------------------------------------------------

   allocate( t_pert_scale( ncrms, nzm ) )
   allocate( q_pert_scale( ncrms, nzm ) )

   idx_qt = index_water_vapor

   !----------------------------------------------------------------------------
   ! calculate local tendencies scaled by local perturbations
   !----------------------------------------------------------------------------
   !$acc parallel loop collapse(2) async(asyncid)
   do k=1,nzm
      do icrm = 1 , ncrms
         ! initialize scaling factors to 1.0
         t_pert_scale(icrm,k) = 1.0_crm_rknd
         q_pert_scale(icrm,k) = 1.0_crm_rknd
         tmp_t_scale = -1
         tmp_q_scale = -1
         ! set scaling factors as long as there are perturbations to scale
         if (t_vt(icrm,k)>0) tmp_t_scale = 1.0_crm_rknd + dtn * t_vt_tend(icrm,k) / t_vt(icrm,k)
         if (q_vt(icrm,k)>0) tmp_q_scale = 1.0_crm_rknd + dtn * q_vt_tend(icrm,k) / q_vt(icrm,k)
         if (tmp_t_scale>0) t_pert_scale(icrm,k) = sqrt( tmp_t_scale )
         if (tmp_q_scale>0) q_pert_scale(icrm,k) = sqrt( tmp_q_scale )
         ! enforce minimum scaling
         t_pert_scale(icrm,k) = max( t_pert_scale(icrm,k), pert_scale_min )
         q_pert_scale(icrm,k) = max( q_pert_scale(icrm,k), pert_scale_min )
         ! enforce maximum scaling
         t_pert_scale(icrm,k) = min( t_pert_scale(icrm,k), pert_scale_max )
         q_pert_scale(icrm,k) = min( q_pert_scale(icrm,k), pert_scale_max )
      end do
   end do

   !----------------------------------------------------------------------------
   ! apply tendencies
   !----------------------------------------------------------------------------
   !$acc parallel loop collapse(4) async(asyncid)
   do k = 1,nzm
      do j = 1,ny
         do i = 1,nx
            do icrm = 1,ncrms
               ttend_loc = ( t_pert_scale(icrm,k) * t_vt_pert(icrm,i,j,k) - t_vt_pert(icrm,i,j,k) ) / dtn
               qtend_loc = ( q_pert_scale(icrm,k) * q_vt_pert(icrm,i,j,k) - q_vt_pert(icrm,i,j,k) ) / dtn

               t(icrm,i,j,k)                  = t(icrm,i,j,k)                  + ttend_loc * dtn
               micro_field(icrm,i,j,k,idx_qt) = micro_field(icrm,i,j,k,idx_qt) + qtend_loc * dtn
            end do
         end do
      end do
   end do

   !----------------------------------------------------------------------------
   ! clean up
   !----------------------------------------------------------------------------
   deallocate( t_pert_scale )
   deallocate( q_pert_scale )

end subroutine VT_forcing

!===============================================================================
!===============================================================================
end module variance_transport_mod
