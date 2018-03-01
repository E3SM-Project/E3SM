   subroutine dyn(irow    ,grlps1  ,grt1    ,grz1    ,grd1    ,  &
                  grfu1   ,grfv1   ,grut1   ,grvt1   ,grrh1   ,  &
                  grlps2  ,grt2    ,grz2    ,grd2    ,grfu2   ,  &
                  grfv2   ,grut2   ,grvt2   ,grrh2, ztodt   )
!-----------------------------------------------------------------------
!
! Combine undifferentiated and longitudinally differentiated Fourier
! coefficient terms for later use in the Gaussian quadrature
!
! Computational note: Index "2*m-1" refers to the real part of the
! complex coefficient, and "2*m" to the imaginary.
!
! The naming convention is as follows:
!  - t, q, d, z refer to temperature, specific humidity, divergence
!     and vorticity
!  - "1" suffix to an array => symmetric component of current latitude pair
!  - "2" suffix to an array => antisymmetric component
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, B. Boville, J. Hack, August 1992
! Reviewed:          D. Williamson, March 1996
! Modified:          P. Worley, September 2002
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
      use time_manager, only: get_step_size, is_first_step
      use spmd_utils, only: iam
      implicit none

!
! Input arguments
!
      integer irow                ! latitude pair index
!
! Input/output arguments
!
      real(r8) grlps1(2*maxm)        ! sym. surface pressure equation term
      real(r8) grt1(2*maxm,plev)     ! sym. undifferentiated term in t eqn.
      real(r8) grz1(2*maxm,plev)     ! sym. undifferentiated term in z eqn.
      real(r8) grd1(2*maxm,plev)     ! sym. undifferentiated term in d eqn.
      real(r8) grfu1(2*maxm,plev)    ! sym. nonlinear terms in u eqn.
      real(r8) grfv1(2*maxm,plev)    ! sym. nonlinear terms in v eqn.
      real(r8) grut1(2*maxm,plev)    ! sym. lambda derivative term in t eqn.
      real(r8) grvt1(2*maxm,plev)    ! sym. mu derivative term in t eqn.
      real(r8) grrh1(2*maxm,plev)    ! sym. RHS of divergence eqn (del^2 term)
      real(r8) grlps2(2*maxm)        ! antisym. surface pressure equation term
      real(r8) grt2(2*maxm,plev)     ! antisym. undifferentiated term in t eqn.
      real(r8) grz2(2*maxm,plev)     ! antisym. undifferentiated term in z eqn.
      real(r8) grd2(2*maxm,plev)     ! antisym. undifferentiated term in d eqn.
      real(r8) grfu2(2*maxm,plev)    ! antisym. nonlinear terms in u eqn.
      real(r8) grfv2(2*maxm,plev)    ! antisym. nonlinear terms in v eqn.
      real(r8) grut2(2*maxm,plev)    ! antisym. lambda derivative term in t eqn.
      real(r8) grvt2(2*maxm,plev)    ! antisym. mu derivative term in t eqn.
      real(r8) grrh2(2*maxm,plev)    ! antisym. RHS of divergence eqn (del^2 term)
      real(r8) ztodt
!
!---------------------------Local workspace-----------------------------
!
      real(r8) tmp1,tmp2              ! temporaries
      real(r8) zxm(pmmax)             ! m*2dt/(a*cos(lat)**2)
      real(r8) zrcsj                  ! 1./(a*cos(lat)**2)
!      real(r8) dtime                  ! timestep size [seconds]
      real(r8) ztdtrc                 ! 2dt/(a*cos(lat)**2)  1dt/..... at nstep=0
      integer lm, mlength             ! local Fourier wavenumber index
                                      !  and number of local indices
      integer k                       ! level index
!DIR$ NOSTREAM
!
! Set constants
!
      mlength = numm(iam)
!      dtime = get_step_size()

      zrcsj = 1._r8/(cs(irow)*rearth)
      ztdtrc = ztodt*zrcsj

!      if (is_first_step()) then
!         ztdtrc = dtime*zrcsj
!      else
!         ztdtrc = 2.0_r8*dtime*zrcsj
!      end if
!
! Combine constants with Fourier wavenumber m
!
!cdir altcode=(loopcnt)
      do lm=1,mlength
         zxm(lm) = ztdtrc*xm(locm(lm,iam))
      end do
!
! Combine undifferentiated and longitudinal derivative terms for
! later use in Gaussian quadrature
!
      do k=1,plev
!cdir loopchg
         do lm=1,mlength
            grt1(2*lm-1,k) = grt1(2*lm-1,k) + zxm(lm)*grut1(2*lm,k)
            grt1(2*lm,k)   = grt1(2*lm,k)   - zxm(lm)*grut1(2*lm-1,k)
            grd1(2*lm-1,k) = grd1(2*lm-1,k) - zxm(lm)*grfu1(2*lm,k)
            grd1(2*lm,k)   = grd1(2*lm,k)   + zxm(lm)*grfu1(2*lm-1,k)
            grz1(2*lm-1,k) = grz1(2*lm-1,k) - zxm(lm)*grfv1(2*lm,k)
            grz1(2*lm,k)   = grz1(2*lm,k)   + zxm(lm)*grfv1(2*lm-1,k)
!
            grt2(2*lm-1,k) = grt2(2*lm-1,k) + zxm(lm)*grut2(2*lm,k)
            grt2(2*lm,k)   = grt2(2*lm,k)   - zxm(lm)*grut2(2*lm-1,k)
            grd2(2*lm-1,k) = grd2(2*lm-1,k) - zxm(lm)*grfu2(2*lm,k)
            grd2(2*lm,k)   = grd2(2*lm,k)   + zxm(lm)*grfu2(2*lm-1,k)
            grz2(2*lm-1,k) = grz2(2*lm-1,k) - zxm(lm)*grfv2(2*lm,k)
            grz2(2*lm,k)   = grz2(2*lm,k)   + zxm(lm)*grfv2(2*lm-1,k)
         end do
      end do

      return
   end subroutine dyn

