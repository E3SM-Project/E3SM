
subroutine reordp(irow    ,iy      ,zalp    ,zdalp   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Renormalize associated Legendre polynomials and their derivatives.
! 
! Method: 
! Reorder associated Legendre polynomials and their derivatives from
! column rectangular storage to diagonal pentagonal storage. The
! reordered polynomials and derivatives are returned via common/comspe/
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pspect
  use comspe
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: irow            ! latitude pair index
  integer , intent(in)  :: iy              ! dimension of input polynomials
  real(r8), intent(in)  :: zalp(iy)        ! Legendre polynomial
  real(r8), intent(in)  :: zdalp(iy)       ! Legendre polynomial derivative
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer mr              ! spectral index
  integer m               ! index along diagonal and row
  integer n               ! index of diagonal
  real(r8) sqrt2              ! sqrt(2)
!-----------------------------------------------------------------------
!
! Multiply ALP and DALP by SQRT(2.) in order to get proper
! normalization. DALP is multiplied by -1 to correct for - sign
! in Copenhagen definition.
!
  sqrt2 = sqrt(2._r8)
  do m=1,pmmax
     mr = nstart(m)
     do n=1,nlen(m)
        alp(mr+n,irow) = zalp((m-1)*pmax + n)*sqrt2
        dalp(mr+n,irow) = -zdalp((m-1)*pmax + n)*sqrt2
     end do
  end do

  return
end subroutine reordp

