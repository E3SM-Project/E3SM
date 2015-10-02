
subroutine settau(zdt)
!-----------------------------------------------------------------------
!
! Purpose:
! Set time invariant hydrostatic matrices, which depend on the reference
! temperature and pressure in the semi-implicit time step. Note that
! this subroutine is actually called twice, because the effective time
! step changes between step 0 and step 1.
!
! zdt = delta t for next semi-implicit time step.
!
! Original version:  CCM1
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
  use comspe
  use commap
  use scanslt,      only: epssld
  use physconst,    only: cappa, rair, gravit, omega
  use spmd_utils,   only: masterproc
  use hycoef,       only: hypd, hypi, hybi
  use cam_logfile,  only: iulog
  implicit none

!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: zdt  ! time step (or dt/2 at time 0)
!
!---------------------------Local workspace-----------------------------
!
  real(r8) zci(plev)             ! dummy, used to print phase speeds
  real(r8) zdt2                  ! zdt**2
  real(r8) factor                ! intermediate workspace
  real(r8) zdt0u                 ! vertical diff. of ref. temp (above)
  real(r8) zshu                  ! interface "sigma" (above)
  real(r8) zr2ds                 ! 1./(2.*hypd(k))
  real(r8) zdt0d                 ! vertical diff. of ref. temp (below)
  real(r8) zshd                  ! interface "sigma" (below)
  real(r8) ztd                   ! temporary accumulator
  real(r8) zb(plev,plev)         ! semi-implicit matrix in d equation
  real(r8) zcr1(plev)            ! real(r8) eigenvalues of semi-impl. matrix
  real(r8) onepeps               ! (1 + epssld)
  real(r8) onepepss              ! (1 + epssld)**2
  real(r8) dnml                  ! tmp variable
  real(r8) dpnml                 ! tmp variable
  real(r8) fm                    ! tmp index for computation purposes
  real(r8) fn                    ! tmp index for computation purposes
  real(r8) alphal                ! (1 + eps)*omega*2*deltat
  integer k,kk,kkk               ! level indices
  integer m                      ! fourier wavenumber index
  integer n                      ! n-wavenumber index
  integer fnp1                   ! fn+1
  integer mr,mc                  ! real and imaginary spectral indices
  integer ir,ii                  ! real and imaginary spectral indices
!
!-----------------------------------------------------------------------
!
  onepeps  = 1._r8 + epssld
  onepepss = onepeps**2
  zdt2 = zdt*zdt
!
! Set mean temperature
! NOTE: Making t0 an actual function of height ***DOES NOT WORK***
!
  do k=1,plev
     t0(k) = 350._r8
  end do
!
! Calculate thermodynamic matrix tau.  1st index = column; 2nd index =
! row of matrix.
!
  zdt0u = 0._r8
  zshu = 0._r8
  do k=1,plev
     zr2ds = 1._r8/(2._r8*hypd(k))
     if (k.lt.plev) then
        zdt0d = t0(k+1) - t0(k)
        zshd = hybi(k+1)
     else
        zdt0d = 0._r8
        zshd = 0._r8
     end if
!
     factor = ((zdt0u*zshu + zdt0d*zshd) - (zdt0d + zdt0u))*zr2ds
     do kk=1,k-1
        tau(kk,k) = factor*hypd(kk) + cappa*t0(k)*ecref(kk,k)
     end do
!
     factor = (zdt0u*zshu + zdt0d*zshd - zdt0d)*zr2ds
     tau(k,k) = factor*hypd(k) + cappa*t0(k)*ecref(k,k)
!
     factor = (zdt0u*zshu + zdt0d*zshd)*zr2ds
     do kk=k+1,plev
        tau(kk,k) = factor*hypd(kk)
     end do
     zdt0u = zdt0d
     zshu = zshd
  end do
!
! Vector for linear surface pressure term in divergence
! Pressure gradient and diagonal term of hydrostatic components
!
  do k=1,plev
     bps(k) = t0(k)
     bps(k) = bps(k)*rair
  end do
  do k=1,plev
     do kk=1,plev
        ztd = bps(k) * hypd(kk)/hypi(plevp)
        do kkk=1,plev
           ztd = ztd + href(kkk,k)*tau(kk,kkk)
        end do
        zb(kk,k) = ztd
     end do
  end do
!
! Determine eigenvalues/vectors of hydrostatic matrix (reference atm)
! for use in computing vertical normal modes.
!
  call nmmatrix(zb      ,zcr1    ,bm1     ,bmi     )
!
! Compute and print gravity wave equivalent depths and phase speeds
!
  do k=1,plev
     zci(k) = sqrt(zcr1(k))
  end do

  if (masterproc) then
     write(iulog,910) (t0(k),k=1,plev)
     write(iulog,920) (zci(k),k=1,plev)
  end if

  do k=1,plev
     zci(k) = zcr1(k) / gravit
  end do

  if (masterproc) then
     write(iulog,930) (zci(k),k=1,plev)
  end if
!
! Compute zcr(n)=(1+e)**2*sq*g*D(l)*delt**2
!
  do k=1,plev
     do n=2,pnmax
        zcr(n,k) = onepepss*zcr1(k)*zdt2*sq(n)
        zcr(n,k) = onepepss*zcr1(k)*zdt2*sq(n)
     end do
  end do
!
! Overwrite zcr for n=1
!
  do k=1,plev
     zcr(1,k) = 0._r8
  end do
!
! Compute coefficients to be used in normal mode space.
! NOTE:  Storage will be sequential along columns ("N")
!
  alphal = onepeps*omega*2._r8*zdt
  do m = 1,pmmax
     fm = m - 1
     mr = nstart(m)
     mc = 2*mr
     do n = 1,nlen(m)
        ir = mc + 2*n - 1
        ii = ir + 1
        fn    = fm + n - 1
        fnp1  = fn + 1
        dnml  = sqrt( (fn  *fn   - fm*fm)/(4._r8*fn  *fn   - 1._r8) )
        dpnml = sqrt( (fnp1*fnp1 - fm*fm)/(4._r8*fnp1*fnp1 - 1._r8) )
!
        a0nm(ir) = 1._r8
        bpnm(ir) = fn*alphal/(fnp1) * dpnml
        if(fn .eq. 0._r8) then
           a0nm(ii) = 0._r8
           bmnm(ir) = 0._r8
        else
           a0nm(ii) = -fm*alphal/(fn*fnp1)
           bmnm(ir) =  fnp1*alphal/(fn) * dnml
        endif
     end do
!
! Compute coefficients to be used in normal mode space to solve tri-
! diagonal matrix.
! NOTE:  Storage will be sequential along columns ("N")
!
     call tricoef(nlen(m) ,a0nm(mc+1),bpnm(mc+1),bmnm(mc+1),atri(mc+1), &
                  btri(mc+1),ctri(mc+1) )
  end do
!
  return
!
! Formats
! 
910 format(' REFERENCE TEMPERATURES FOR SEMI-IMPLICIT SCHEME = '  /(1x,12f9.3))
920 format(' GRAVITY WAVE PHASE SPEEDS (M/S) FOR MEAN STATE = '   /(1x,12f9.3))
930 format(' GRAVITY WAVE EQUIVALENT DEPTHS (M) FOR MEAN STATE = '/(1x,12f9.3))
end subroutine settau

