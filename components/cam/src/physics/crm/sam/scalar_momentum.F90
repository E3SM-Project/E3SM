module scalar_momentum_mod
!-----------------------------------------------------------------------
! Purpose:
!
! Explicit Scalar Momentum Transport (ESMT)
! Transport large-scale momentum as non-conserved scalars,
! including estimate of 3D pressure gradient force. This was
! implemented to help avoid using any physics modules outside
! of CRM and radiation routines. These non-SP routines can be
! problematic when changing the vertical grid.
! See Tulich (2015) for further details on ESMT
!
! Revision history:
! Nov, 2017 - Walter Hannah - Lawrence Livermore National Lab
!              initial version based on crmtracers.F90
!              Possoin solver and fft routines provided by Stefan Tulich
!
!---------------------------------------------------------------------------
   use params
   use grid, only: nx,ny,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s
   use openacc_utils

   implicit none

   public allocate_scalar_momentum
   public deallocate_scalar_momentum
   public scalar_momentum_tend 

   real(crm_rknd), allocatable :: u_esmt(:,:,:,:)       ! scalar zonal velocity
   real(crm_rknd), allocatable :: v_esmt(:,:,:,:)       ! scalar meridonal velocity

   real(crm_rknd), allocatable :: u_esmt_sgs (:,:)      ! SGS vertical flux
   real(crm_rknd), allocatable :: v_esmt_sgs (:,:)      !
   real(crm_rknd), allocatable :: u_esmt_diff(:,:)      ! large-scale tendency due to vertical diffusion
   real(crm_rknd), allocatable :: v_esmt_diff(:,:)      !

   real(crm_rknd), allocatable :: fluxb_u_esmt(:,:,:)   ! flux of u_esmt at surface    (normally set to zero)
   real(crm_rknd), allocatable :: fluxb_v_esmt(:,:,:)   ! flux of v_esmt at surface    (normally set to zero)
   real(crm_rknd), allocatable :: fluxt_u_esmt(:,:,:)   ! flux of u_esmt at model top  (normally set to zero)
   real(crm_rknd), allocatable :: fluxt_v_esmt(:,:,:)   ! flux of v_esmt at model top  (normally set to zero)

   character*30 u_esmt_name
   character*30 v_esmt_name
   character*10 esmt_units

contains

!========================================================================================
!========================================================================================
subroutine allocate_scalar_momentum(ncrms)
   !------------------------------------------------------------------
   ! Purpose: Allocate and initialize variables for ESMT
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   !------------------------------------------------------------------
   implicit none
   integer, intent(in) :: ncrms

   allocate( u_esmt(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )
   allocate( v_esmt(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )
   allocate( fluxb_u_esmt (nx, ny,ncrms) )
   allocate( fluxb_v_esmt (nx, ny,ncrms) )
   allocate( fluxt_u_esmt (nx, ny,ncrms) )
   allocate( fluxt_v_esmt (nx, ny,ncrms) )
   allocate( u_esmt_sgs   (nz,ncrms)  )
   allocate( v_esmt_sgs   (nz,ncrms)  )
   allocate( u_esmt_diff  (nz,ncrms)  )
   allocate( v_esmt_diff  (nz,ncrms)  )

   call prefetch( u_esmt )
   call prefetch( v_esmt )
   call prefetch( fluxb_u_esmt )
   call prefetch( fluxb_v_esmt )
   call prefetch( fluxt_u_esmt )
   call prefetch( fluxt_v_esmt )
   call prefetch( u_esmt_sgs )
   call prefetch( v_esmt_sgs )
   call prefetch( u_esmt_diff )
   call prefetch( v_esmt_diff )

   u_esmt       = 0.0_crm_rknd
   v_esmt       = 0.0_crm_rknd
   fluxb_u_esmt = 0.0_crm_rknd
   fluxb_v_esmt = 0.0_crm_rknd
   fluxt_u_esmt = 0.0_crm_rknd
   fluxt_v_esmt = 0.0_crm_rknd
   u_esmt_sgs   = 0.0_crm_rknd
   u_esmt_diff  = 0.0_crm_rknd
   v_esmt_sgs   = 0.0_crm_rknd
   v_esmt_diff  = 0.0_crm_rknd

   u_esmt_name = 'Zonal Velocity'
   v_esmt_name = 'Meridonal Velocity'
   esmt_units  = 'm/s'

end subroutine allocate_scalar_momentum

!========================================================================================
!========================================================================================
subroutine deallocate_scalar_momentum()
   !------------------------------------------------------------------
   ! Purpose: Deallocate ESMT variables
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   !------------------------------------------------------------------
   implicit none
   deallocate( u_esmt       )
   deallocate( v_esmt       )
   deallocate( fluxb_u_esmt )
   deallocate( fluxb_v_esmt )
   deallocate( fluxt_u_esmt )
   deallocate( fluxt_v_esmt )
   deallocate( u_esmt_sgs   )
   deallocate( v_esmt_sgs   )
   deallocate( u_esmt_diff  )
   deallocate( v_esmt_diff  )
end subroutine deallocate_scalar_momentum

!========================================================================================
!========================================================================================
subroutine scalar_momentum_tend(ncrms)
   !------------------------------------------------------------------
   ! Purpose: Calculate pressure gradient effects on scalar momentum
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   !------------------------------------------------------------------
   use grid
   implicit none
   integer, intent(in) :: ncrms

   real(crm_rknd), dimension(nx,ny,nzm) :: u_esmt_pgf_3D
   real(crm_rknd), dimension(nx,ny,nzm) :: v_esmt_pgf_3D
   real(crm_rknd) factor_xy
   integer :: i,j,k,icrm

   do icrm = 1 , ncrms
      factor_xy = 1._crm_rknd/real(nx*ny,crm_rknd)

      call scalar_momentum_pgf(ncrms,icrm,u_esmt(icrm,:,:,:),u_esmt_pgf_3D)
      call scalar_momentum_pgf(ncrms,icrm,v_esmt(icrm,:,:,:),v_esmt_pgf_3D)

      ! Add PGF tendency
      do k=1,nzm
         do j=1,ny
            do i=1,nx
               u_esmt(icrm,i,j,k) = u_esmt(icrm,i,j,k) + u_esmt_pgf_3D(i,j,k)*dtn
               v_esmt(icrm,i,j,k) = v_esmt(icrm,i,j,k) + v_esmt_pgf_3D(i,j,k)*dtn
            end do
         end do
      end do

   end do

end subroutine scalar_momentum_tend

!========================================================================================
!========================================================================================

subroutine scalar_momentum_pgf( ncrms, icrm, u_s, tend )
   !------------------------------------------------------------------
   ! Purpose: calculate pgf for scalar momentum transport
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   ! adapted from SP-WRF code by Stefan Tulich
   !------------------------------------------------------------------
   use grid,    only: nx,ny,nz,nzm,z,pres,zi
   use crmdims, only: crm_dx
   use vars,    only: w, rho

   implicit none
   integer, intent(in) :: ncrms,icrm
   !------------------------------------------------------------------
   ! interface variables
   !------------------------------------------------------------------
   real(crm_rknd), dimension(nx,ny,nzm), intent(in ) :: u_s      ! scalar momentum
   real(crm_rknd), dimension(nx,ny,nzm), intent(out) :: tend     ! output tendency
   !------------------------------------------------------------------
   ! local variables
   !------------------------------------------------------------------
   integer :: i,j,k
   real(crm_rknd) :: dampwt
   real(crm_rknd), dimension(nx)       :: k_arr      ! "zonal" wavelength from forward FFT
   real(crm_rknd), dimension(nzm)      :: u_s_avg    ! horizonal average of u_s
   real(crm_rknd), dimension(nzm)      :: a,b,c      ! Poisson solver boundary conditions
   real(crm_rknd), dimension(nzm)      :: rhs        ! right-hand side of pressure equation
   real(crm_rknd), dimension(nzm)      :: shr        ! vertical shear of scalar momentum (du/dz)
   real(crm_rknd), dimension(nx,nzm)   :: w_i        ! w interpolated to scalar levels
   real(crm_rknd), dimension(nx,nzm)   :: w_hat      ! w after Fourier Transform
   real(crm_rknd), dimension(nx,nzm)   :: pgf_hat    ! pressure grad force (Fourier transform space)
   real(crm_rknd), dimension(nx,nzm)   :: pgf        ! pressure grad force for final tendency

   real(crm_rknd), dimension(nzm+1)    :: dz         ! layer thickness
   !------------------------------------------------------------------------
   !------------------------------------------------------------------------

   ! The loop over "y" points is mostly unessary, since ESMT
   ! is for 2D CRMs, but it is useful for directly comparing
   ! ESMT tendencies to fully resolved 3D momentum transport

   ! Calculate layer thickness
   do k = 1,nzm
      dz(k) = zi(icrm,k+1)-zi(icrm,k)
   enddo
   dz(nzm+1) = dz(nzm)

   do j=1,ny

      !-----------------------------------------
      ! Initialize stuff for averaging
      !-----------------------------------------
      u_s_avg(:) = 0.0
      shr(:) = 0.0_crm_rknd

      !-----------------------------------------
      ! Calculate shear of domain average profile
      ! defined on scalar levels
      !-----------------------------------------
      do k = 1,nzm
         do i = 1,nx
            u_s_avg(k) = u_s_avg(k) + u_s(i,j,k)
            ! note that w is on interface levels
            w_i(i,k) = ( w(icrm,i,j,k) + w(icrm,i,j,k+1) )/2.0_crm_rknd
         end do
         u_s_avg(k) = u_s_avg(k) / real(nx,crm_rknd)
      end do

      shr(1) = ( u_s_avg(2) - u_s_avg(1) ) / ( z(icrm,2) - z(icrm,1) )
      do k = 2,nzm-1
         shr(k) = ( u_s_avg(k+1) - u_s_avg(k-1) )/(z(icrm,k+1)-z(icrm,k-1))
      end do

      !------------------------------------------------------------------------
      !------------------------------------------------------------------------
      ! Use Poisson solver to calculate pressure gradient force (PGF)
      ! pgf is diagnosed from w * du_si/dz using the poisson equation
      ! (see Wu and Yanai 1994)
      !------------------------------------------------------------------------
      !------------------------------------------------------------------------

      !-----------------------------------------
      ! compute forward fft of w
      !-----------------------------------------
      call esmt_fft_forward(nx,nzm,crm_dx,w_i,k_arr,w_hat)

      !-----------------------------------------
      ! solve vertical structure equation
      ! for each zonal wavelength (k_arr)
      ! solution method involves constructing
      ! a tridiagonal matrix
      !-----------------------------------------

      pgf_hat(:,:) = 0.

      ! Loop through wavelengths
      do i = 2,nx

         do k = 1,nzm
            a(k) = dz(k+1) / ( dz(k+1) + dz(k) )
            ! the factor of 1.25 crudely accounts for difference between 2D and 3D updraft geometry
            b(k) = -0.5_crm_rknd * 1.25_crm_rknd * k_arr(i)**2.0_crm_rknd * dz(k) * dz(k+1) - 1.0_crm_rknd 
            c(k) = dz(k) / ( dz(k+1) + dz(k) )
            rhs(k) = k_arr(i)**2 * w_hat(i,k) * shr(k) * dz(k) * dz(k+1)
         end do ! k

         !lower boundary condition (symmetric)
         b(1) = b(1) + a(1)
         a(1) = 0.0_crm_rknd

         !upper boundary condition (symmetric)
         b(nzm) = b(nzm) + c(nzm)
         c(nzm) = 0.0_crm_rknd

         ! gaussian elimination with no pivoting
         do k = 1,nzm-1
            b(k+1) = b(k+1) - a(k+1) / b(k) * c(k)
            rhs(k+1) = rhs(k+1) - a(k+1) / b(k) * rhs(k)
         end do ! k

         ! backward substitution
         rhs(nzm) = rhs(nzm) / b(nzm)
         do k=nzm-1,1,-1
            rhs(k) = ( rhs(k) - c(k) * rhs(k+1) ) / b(k)
         end do

         do k = 1,nzm
            pgf_hat(i,k) = rhs(k)
         end do ! k
      end do ! i - zonal wavelength

      ! Note sure what this part does... 
      ! something about the Nyquist freq and whether nx is odd or even
      if (mod(nx,2) == 0) then
         do k=1,nzm
            pgf_hat(nx,k) = pgf_hat(nx,k) / 2.0_crm_rknd
         end do ! k
      end if

      !-----------------------------------------
      ! invert fft of pgf_hat to get pgf
      !-----------------------------------------
      call esmt_fft_backward(nx,nzm,pgf_hat,pgf)

      !-----------------------------------------
      ! Compute final tendency
      !-----------------------------------------
      do k = 1,nzm
         do i = 1,nx
            if (k.eq.1) then
               tend(i,j,k) = 0.0_crm_rknd
            else
               tend(i,j,k) = -1.0_crm_rknd * pgf(i,k) * rho(icrm,k)
            end if
         enddo ! i
      end do ! k

      !------------------------------------------------------------------------
      !------------------------------------------------------------------------
      ! End of Possion solver
      !------------------------------------------------------------------------
      !------------------------------------------------------------------------

   end do ! j

end subroutine  scalar_momentum_pgf

!========================================================================================
!========================================================================================

subroutine esmt_fft_forward(nx,nzm,dx,arr_in,k_out,arr_out)
   !------------------------------------------------------------------
   ! Purpose: calculate forward FFT transform
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   ! adapted from SP-WRF code provided by Stefan Tulich
   !------------------------------------------------------------------
   use fftpack51D
   implicit none
   integer                          , intent(in ) :: nx
   integer                          , intent(in ) :: nzm
   real(crm_rknd)                   , intent(in ) :: dx
   real(crm_rknd), dimension(nx,nzm), intent(in ) :: arr_in
   real(crm_rknd), dimension(nx)    , intent(out) :: k_out
   real(crm_rknd), dimension(nx,nzm), intent(out) :: arr_out
   ! local variables
   integer :: lensave, ier, nh, n1, i, j, k
   integer :: lot, jump, n, inc, lenr, lensav, lenwrk
   real(crm_rknd) :: pi
   real(crm_rknd), dimension(nx+15)  :: wsave
   real(crm_rknd), dimension(nx,nzm) :: work

   ! naming convention follows fftpack5 routines
   pi = 2.*asin(1.0)
   n = nx
   lot = 1
   lensav = n+15
   inc = 1
   lenr = nx
   jump = nx
   lenwrk = lenr
   ! initialization for FFT
   call rfft1i(n,wsave,lensav,ier)
   if(ier /= 0) write(0,*) 'ERROR: rfftmi(): ESMT - FFT initialization error ',ier
   !  do the forward transform
   do k = 1,nzm
      do i = 1,nx
         arr_out(i,k) = arr_in(i,k)
      enddo
      call rfft1f( n, inc, arr_out(:,k), lenr, wsave, lensav, work(:,k), lenwrk, ier )
      if(ier /= 0) write(0,*) 'ERROR: rfftmf(): ESMT - Forward FFT error ',ier
   enddo

   if(mod(n,2) == 0) then
      nh = n/2 - 1
   else
      nh = (n-1)/2
   endif

   k_out(1)=0.0
   do j = 1,nh
      k_out(2*j)   = 2.*pi*real(j,crm_rknd)/(real(n,crm_rknd)*dx)   !cos
      k_out(2*j+1) = 2.*pi*real(j,crm_rknd)/(real(n,crm_rknd)*dx)   !sin
   enddo
   if (mod(n,2) == 0) then
      k_out(n) =  2.*pi/(2.*dx)                  !nyquist wavelength for even n
   end if
end subroutine esmt_fft_forward

!========================================================================================
!========================================================================================

subroutine esmt_fft_backward(nx,nzm,arr_in,arr_out)
   !------------------------------------------------------------------
   ! Purpose: calculate backward FFT transform
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   ! adapted from SP-WRF code provided by Stefan Tulich
   !------------------------------------------------------------------
   use fftpack51D
   implicit none
   integer                         , intent(in ) :: nx
   integer                         , intent(in ) :: nzm
   real(crm_rknd), dimension(nx,nzm), intent(in ) :: arr_in
   real(crm_rknd), dimension(nx,nzm), intent(out) :: arr_out
   ! local variables
   integer :: lensave, ier, nh, n1, i, k
   integer :: lot, jump, n, inc, lenr, lensav, lenwrk
   real(crm_rknd), dimension(nx+15) :: wsave
   real(crm_rknd), dimension(nx,nzm) :: work

   ! naming convention follows fftpack5 routines
   n = nx
   lot = 1
   lensav = n+15
   inc = 1
   lenr = nx
   jump = nx
   lenwrk = lenr
   ! initialization for FFT
   call rfft1i(n,wsave,lensav,ier)
   if(ier /= 0) write(0,*) 'ERROR: rfftmi(): ESMT - FFT initialization error ',ier
   !  do the backward transform
   do k = 1,nzm
      do i = 1,nx
         arr_out(i,k) = arr_in(i,k)
      enddo
      call rfft1b( n, inc, arr_out(:,k), lenr, wsave, lensav, work(:,k), lenwrk, ier )
      if(ier /= 0) write(0,*) 'ERROR: rfftmb(): ESMT - backward FFT error ',ier
   enddo
end subroutine esmt_fft_backward

end module scalar_momentum_mod
