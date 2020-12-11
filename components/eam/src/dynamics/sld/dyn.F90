subroutine dyn(irow    ,grlps1  ,grt1    ,grq1    ,grz1    , &
               grd1    ,grfu1   ,grfv1   ,grlps2  ,grt2    , &
               grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Combine undifferentiated and longitudinally differentiated Fourier
! coefficient terms for later use in the Gaussian quadrature
!
! Computational note: Index "2*m-1" refers to the real part of the
! complex coefficient, and "2*m" to the imaginary.
!
! The naming convention is as follows:
!  - t, q, d, z refer to temperature, specific humidity, divergence
!     and vorticity
!  - "1" suffix to an array => symmetric component of current latitude
!     pair
!  - "2" suffix to an array => antisymmetric component
!
! Author:   J. Rosinski
! Modified: P. Worley, October 2002
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use rgrid
  use comspe
  use commap
  use physconst, only: rearth
  use spmd_utils, only: iam
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: irow                 ! latitude pair index
  real(r8), intent(in)   :: grlps1(2*maxm)      ! sym. undifferentiated term in Ps eqn.
  real(r8), intent(in)   :: grt1  (2*maxm,plev) ! sym. undifferentiated term in t eqn.
  real(r8), intent(in)   :: grq1  (2*maxm,plev) ! sym. undifferentiated term in q
  real(r8), intent(out)  :: grz1  (2*maxm,plev) ! sym. undifferentiated term in z eqn.
  real(r8), intent(out)  :: grd1  (2*maxm,plev) ! sym. undifferentiated term in d eqn.
  real(r8), intent(in)   :: grfu1 (2*maxm,plev) ! sym. nonlinear terms in u eqn.
  real(r8), intent(in)   :: grfv1 (2*maxm,plev) ! sym. nonlinear terms in v eqn.
  real(r8), intent(in)   :: grlps2(2*maxm)      ! antisym. undifferentiated term in Ps eq
  real(r8), intent(in)   :: grt2  (2*maxm,plev) ! antisym. undifferentiated term in t eq
  real(r8), intent(in)   :: grq2  (2*maxm,plev) ! antisym. undifferentiated term in q
  real(r8), intent(out)  :: grz2  (2*maxm,plev) ! antisym. undifferentiated term in z eq
  real(r8), intent(out)  :: grd2  (2*maxm,plev) ! antisym. undifferentiated term in d eq
  real(r8), intent(in)   :: grfu2 (2*maxm,plev) ! antisym. nonlinear terms in u eqn.
  real(r8), intent(in)   :: grfv2 (2*maxm,plev) ! antisym. nonlinear terms in v eqn.
!
!---------------------------Local workspace-----------------------------
!
  real(r8) tmp1
  real(r8) tmp2
  real(r8) zxm(pmmax)         ! m*2dt/(a*cos(lat)**2)
  real(r8) zrcsj              ! 1./(a*cos(lat)**2)
  integer lm, mlength         ! local Fourier wavenumber index
                              !  and number of local indices
  integer m                   ! Fourier index
  integer k                   ! level index
!
!-----------------------------------------------------------------------
!
! Set constants
!
  mlength = numm(iam)
  zrcsj = 1._r8/(cs(irow)*rearth)
!
! Combine constants with Fourier wavenumber m
!
  do lm=1,mlength
     zxm(lm) = zrcsj*xm(locm(lm,iam))
  end do
!
! Combine undifferentiated and longitudinal derivative terms for
! later use in Gaussian quadrature
!
  do k=1,plev
     do lm=1,mlength
        grd1(2*lm-1,k) = - zxm(lm)*grfu1(2*lm,k)
        grd1(2*lm,k)   =   zxm(lm)*grfu1(2*lm-1,k)
        grz1(2*lm-1,k) = - zxm(lm)*grfv1(2*lm,k)
        grz1(2*lm,k)   =   zxm(lm)*grfv1(2*lm-1,k)
!
        grd2(2*lm-1,k) = - zxm(lm)*grfu2(2*lm,k)
        grd2(2*lm,k)   =   zxm(lm)*grfu2(2*lm-1,k)
        grz2(2*lm-1,k) = - zxm(lm)*grfv2(2*lm,k)
        grz2(2*lm,k)   =   zxm(lm)*grfv2(2*lm-1,k)
     end do
  end do
  return
end subroutine dyn
