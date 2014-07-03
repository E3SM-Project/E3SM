
subroutine limdx(pidim   ,ibeg    ,len     ,dx      ,f       ,&
                 fxl     ,fxr     )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Limit the derivative estimates for data on an equally spaced grid
! so they satisfy the SCM0 condition, that is, the spline will be 
! monotonic, but only C0 continuous on the domain
! 
! Method: 
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use scanslt,      only: plond
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog
!-----------------------------------------------------------------------
   implicit none
!---------------------------Local parameters----------------------------
!
   integer pbpts             ! (length of latitude slice)*fields
   parameter(pbpts = plond) 
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: pidim             ! vector dimension
   integer, intent(in) :: ibeg              ! index of vector to begin computation
   integer, intent(in) :: len               ! length of vector to compute
!
   real(r8), intent(in) :: dx                   ! length of grid inteval
   real(r8), intent(in) :: f(pidim)             ! field
!
! Input/output arguments
!
   real(r8), intent(inout) :: fxl(pidim)           ! x-derivs at left  edge of interval
   real(r8), intent(inout) :: fxr(pidim)           ! x-derivs at right edge of interval
!
!-----------------------------------------------------------------------
!
!  pidim   Length of f, fxl, and fxr.
!  ibeg    First interval of grid for which derivatives are computed.
!  len     Number of grid intervals for which derivatives are computed.
!          (There are pidim - 1 intervals between the pidim gridpoints
!          represented in f, fxl, and fxr.)
!  dx      Value of grid spacing.
!  f       Values on equally spaced grid from which derivatives fxl and
!          fxr were computed.
!  fxl     fxl(i) is the limited derivative at the left  edge of
!          interval
!  fxr     fxr(i) is the limited derivative at the right edge of
!          interval
!
!---------------------------Local variables-----------------------------
!
   integer i                 ! index
   integer iend              ! index to end work on vector
!
   real(r8) rdx                  ! 1./dx
   real(r8) deli(pbpts)          ! simple linear derivative
!
!-----------------------------------------------------------------------
!
   if(pidim .gt. pbpts) then
      write(iulog,9000) pidim
      call endrun
   end if
!
   iend = ibeg + len - 1
   rdx  = 1._r8/dx
!
   do i = ibeg,iend
      deli(i) = ( f(i+1) - f(i) )*rdx
   end do
!
! Limiter
!
   call scm0(len     ,deli(ibeg),fxl(ibeg),fxr(ibeg))
!
   return
9000 format('LIMDX: Local work array DELI not dimensioned large enough' &
             ,/'       Increase local parameter pbpts to ',i5)
end subroutine limdx

