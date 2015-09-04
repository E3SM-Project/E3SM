
subroutine scm0(n       ,deli    ,df1     ,df2     )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Apply SCM0 limiter to derivative estimates.
! See Rasch and Williamson (1990)
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)    :: n         ! length of vectors
  real(r8), intent(in)    :: deli(n)   ! discrete derivative
  real(r8), intent(inout) :: df1(n)    ! limited left -edge derivative
  real(r8), intent(inout) :: df2(n)    ! limited right-edge derivative
!
! n      Dimension of input arrays.
! deli   deli(i) is the discrete derivative on interval i, i.e.,
!        deli(i) = ( f(i+1) - f(i) )/( x(i+1) - x(i) ).
! df1    df1(i) is the limited derivative at the left  edge of interval
! df2    df2(i) is the limited derivative at the right edge of interval
!-----------------------------------------------------------------------


!---------------------------Local variables-----------------------------
  integer i                 ! index
  real(r8) fac              ! factor applied in limiter
  real(r8) tmp1             ! derivative factor
  real(r8) tmp2             ! abs(tmp1)
!-----------------------------------------------------------------------
!
  fac = 3._r8*(1._r8 - 10._r8*epsilon(fac))
  do i = 1,n
     tmp1 = fac*deli(i)
     tmp2 = abs( tmp1 )
     if( deli(i)*df1(i)   <= 0.0_r8  ) df1(i) = 0._r8
     if( deli(i)*df2(i)   <= 0.0_r8  ) df2(i) = 0._r8
     if( abs(    df1(i) ) >  tmp2 ) df1(i) = tmp1
     if( abs(    df2(i) ) >  tmp2 ) df2(i) = tmp1
  end do

  return
end subroutine scm0

