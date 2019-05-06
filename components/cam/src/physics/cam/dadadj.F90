module dadadj
!----------------------------------------------------------------------- 
! 
! Purpose: 
! GFDL style dry adiabatic adjustment
! 
! Method: 
! if stratification is unstable, adjustment to the dry adiabatic lapse
! rate is forced subject to the condition that enthalpy is conserved.
! 
! Author: J.Hack
! 
!-----------------------------------------------------------------------

use shr_kind_mod,    only: r8 => shr_kind_r8

implicit none
private
save

public :: &
   dadadj_initial, &
   dadadj_calc

integer  :: nlvdry  ! number of layers from top of model to apply the adjustment
integer  :: niter   ! number of iterations for convergence

!===============================================================================
contains
!===============================================================================

subroutine dadadj_initial(nlvdry_in, niter_in)

   integer,  intent(in) :: nlvdry_in
   integer,  intent(in) :: niter_in

   nlvdry = nlvdry_in
   niter = niter_in

end subroutine dadadj_initial

!===============================================================================

subroutine dadadj_calc( &
   ncol, pmid, pint, pdel, cappav, t, &
   q, dadpdf, icol_err)

   ! Arguments

   integer, intent(in) :: ncol                ! number of atmospheric columns

   real(r8), intent(in) :: pmid(:,:)   ! pressure at model levels
   real(r8), intent(in) :: pint(:,:)   ! pressure at model interfaces
   real(r8), intent(in) :: pdel(:,:)   ! vertical delta-p
   real(r8), intent(in) :: cappav(:,:) ! variable Kappa

   real(r8), intent(inout) :: t(:,:)   ! temperature (K)
   real(r8), intent(inout) :: q(:,:)   ! specific humidity
   
   real(r8), intent(out) :: dadpdf(:,:)  ! PDF of where adjustments happened

   integer,  intent(out) :: icol_err ! index of column in which error occurred

   !---------------------------Local workspace-----------------------------

   integer :: i,k             ! longitude, level indices
   integer :: jiter           ! iteration index

   real(r8), allocatable :: c1dad(:) ! intermediate constant
   real(r8), allocatable :: c2dad(:) ! intermediate constant
   real(r8), allocatable :: c3dad(:) ! intermediate constant
   real(r8), allocatable :: c4dad(:) ! intermediate constant
   real(r8) :: gammad    ! dry adiabatic lapse rate (deg/Pa)
   real(r8) :: zeps      ! convergence criterion (deg/Pa)
   real(r8) :: rdenom    ! reciprocal of denominator of expression
   real(r8) :: dtdp      ! delta-t/delta-p
   real(r8) :: zepsdp    ! zeps*delta-p
   real(r8) :: zgamma    ! intermediate constant
   real(r8) :: qave      ! mean q between levels
   real(r8) :: cappa     ! Kappa at level intefaces

   logical :: ilconv          ! .TRUE. ==> convergence was attained
   logical :: dodad(ncol)     ! .TRUE. ==> do dry adjustment

   !-----------------------------------------------------------------------

   icol_err = 0
   zeps = 2.0e-5_r8           ! set convergence criteria

   allocate(c1dad(nlvdry), c2dad(nlvdry), c3dad(nlvdry), c4dad(nlvdry))

   ! Find gridpoints with unstable stratification

   do i = 1, ncol
      cappa = 0.5_r8*(cappav(i,2) + cappav(i,1))
      gammad = cappa*0.5_r8*(t(i,2) + t(i,1))/pint(i,2)
      dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
      dodad(i) = (dtdp + zeps) .gt. gammad
   end do
   
   dadpdf(:ncol,:) = 0._r8
   do k= 2, nlvdry
      do i = 1, ncol
         cappa = 0.5_r8*(cappav(i,k+1) + cappav(i,k))
         gammad = cappa*0.5_r8*(t(i,k+1) + t(i,k))/pint(i,k+1)
         dtdp = (t(i,k+1) - t(i,k))/(pmid(i,k+1) - pmid(i,k))
         dodad(i) = dodad(i) .or. (dtdp + zeps).gt.gammad
         if ((dtdp + zeps).gt.gammad) then
           dadpdf(i,k) = 1._r8
         end if
      end do
   end do

   ! Make a dry adiabatic adjustment
   ! Note: nlvdry ****MUST**** be < pver

   COL: do i = 1, ncol

      if (dodad(i)) then

         zeps = 2.0e-5_r8

         do k = 1, nlvdry
            c1dad(k) = cappa*0.5_r8*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
            c2dad(k) = (1._r8 - c1dad(k))/(1._r8 + c1dad(k))
            rdenom = 1._r8/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
            c3dad(k) = rdenom*pdel(i,k)
            c4dad(k) = rdenom*pdel(i,k+1)
         end do

50       continue

         do jiter = 1, niter
            ilconv = .true.

            do k = 1, nlvdry
               zepsdp = zeps*(pmid(i,k+1) - pmid(i,k))
               zgamma = c1dad(k)*(t(i,k) + t(i,k+1))

               if ((t(i,k+1)-t(i,k)) >= (zgamma+zepsdp)) then
                  ilconv = .false.
                  t(i,k+1) = t(i,k)*c3dad(k) + t(i,k+1)*c4dad(k)
                  t(i,k) = c2dad(k)*t(i,k+1)
                  qave = (pdel(i,k+1)*q(i,k+1) + pdel(i,k)*q(i,k))/(pdel(i,k+1)+ pdel(i,k))
                  q(i,k+1) = qave
                  q(i,k) = qave
               end if

            end do

            if (ilconv) cycle COL ! convergence => next longitude
         end do

         ! Double convergence criterion if no convergence in niter iterations

         zeps = zeps + zeps
         if (zeps > 1.e-4_r8) then
            icol_err = i
            return                ! error return
         else
            go to 50
         end if

      end if

   end do COL

   deallocate(c1dad, c2dad, c3dad, c4dad)

end subroutine dadadj_calc

!===============================================================================

end module dadadj
