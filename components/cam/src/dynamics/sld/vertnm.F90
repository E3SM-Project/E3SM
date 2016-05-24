subroutine vertnm(lm)
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the system of semi-implicit divergence/vorticity
! equations.  Equations have been de-coupled in the vertical by
! transforming the spectral terms into vertical normal modes.
!
! Author:  J. Olson
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
  use spmd_utils, only : iam
  
  implicit none

!------------------------------Arguments--------------------------------
!
  integer, intent(in) :: lm ! local Fourier wavenumber index 
!
!---------------------------Local workspace-----------------------------
!
  integer k                 ! vertical index
  integer m                 ! diagonal element of cmplx array
  integer n                 ! wave number
  integer nr                ! n   (real)
  integer ni                ! n   (imag)
  integer np1r              ! n+1 (real)
  integer np1i              ! n+1 (imag)
  integer nm1r              ! n-1 (real)
  integer nm1i              ! n-1 (imag)
  integer mr                ! real spectral index
  integer mc                ! imaginary spectral index
  real(r8) tmp  (2)         ! real/imaginary temp spaces
  real(r8) tmp1 (2)         ! real/imaginary temp spaces
  real(r8) tmp2 (2)         ! real/imaginary temp spaces
  real(r8) btrie(2*pmax)    ! "btri" + eigenvalue of ref atmosphere
  real(r8) dtri (2*pmax)    ! RHS in the solution of the tri-diagonal matrix
  real(r8) denom            ! denominator
!
!-----------------------------------------------------------------------
!
! Complete computation of btri and compute dtri (complex arithmetic)
!
  m  = locm(lm,iam)
  mr = nstart(m)
  mc = 2*mr
  do k = 1,plev
     do n = 1,nlen(m)
        nr   = mc   + 2*n     - 1
        ni   = nr   + 1
        np1r = mc   + 2*(n+1) - 1
        np1i = np1r + 1
        nm1r = mc   + 2*(n-1) - 1
        nm1i = nm1r + 1
!
        btrie(2*n-1) = btri(nr) + zcr(m+n-1,k)
        btrie(2*n  ) = btri(ni)
!
        tmp1(1) = 0._r8
        tmp1(2) = 0._r8 
        tmp2(1) = 0._r8
        tmp2(2) = 0._r8 
        if(n .ne. nlen(m)) then
           tmp(1)  =  bpnm(nr  )*vznm(np1r,k)
           tmp(2)  =  bpnm(nr  )*vznm(np1i,k)
           denom   =  a0nm(np1r)*a0nm(np1r) + a0nm(np1i)*a0nm(np1i)
           tmp1(1) = (a0nm(np1r)*tmp(1) + a0nm(np1i)*tmp(2))/denom
           tmp1(2) = (a0nm(np1r)*tmp(2) - a0nm(np1i)*tmp(1))/denom
        endif
        if(n .ne. 1    ) then
           tmp(1)  =  bmnm(nr  )*vznm(nm1r,k)
           tmp(2)  =  bmnm(nr  )*vznm(nm1i,k)
           denom   =  a0nm(nm1r)*a0nm(nm1r) + a0nm(nm1i)*a0nm(nm1i)
           tmp2(1) = (a0nm(nm1r)*tmp(1) + a0nm(nm1i)*tmp(2))/denom
           tmp2(2) = (a0nm(nm1r)*tmp(2) - a0nm(nm1i)*tmp(1))/denom
        endif
        dtri(2*n-1) = dsnm(nr,k) + hsnm(nr,k) + tmp1(1) + tmp2(1)
        dtri(2*n  ) = dsnm(ni,k) + hsnm(ni,k) + tmp1(2) + tmp2(2)
     end do
!
! Solve tridiagonal matrix:  call once for ODDs and once for EVENs
!
     if(mod(nlen(m),2) .eq. 0) n = nlen(m) - 1
     if(mod(nlen(m),2) .ne. 0) n = nlen(m)
     call trisolve(n       ,atri(mc+1),btrie(1),ctri(mc+1),dtri(1),dnm(mc+1,k)    )
     if(mod(nlen(m),2) .eq. 0) n = nlen(m) - 1
     if(mod(nlen(m),2) .ne. 0) n = nlen(m) - 2
     if(n .gt. 0) then
        call trisolve(n       ,atri(mc+3),btrie(3),ctri(mc+3),dtri(3),dnm(mc+3,k)    )
     endif
!
! Solve for vorticity
!
     do n = 1,nlen(m)
        nr   = mc   + 2*n     - 1
        ni   = nr   + 1
        np1r = mc   + 2*(n+1) - 1
        np1i = np1r + 1
        nm1r = mc   + 2*(n-1) - 1
        nm1i = nm1r + 1
!
        tmp1(1) = 0._r8
        tmp1(2) = 0._r8 
        tmp2(1) = 0._r8
        tmp2(2) = 0._r8 
        if(n .ne. nlen(m)) then
           tmp1(1) = bpnm(nr)*dnm(np1r,k)
           tmp1(2) = bpnm(nr)*dnm(np1i,k)
        endif
        if(n .ne. 1      ) then
           tmp2(1) = bmnm(nr)*dnm(nm1r,k)
           tmp2(2) = bmnm(nr)*dnm(nm1i,k)
        endif
        vznm(nr,k) = vznm(nr,k) - tmp1(1) - tmp2(1)
        vznm(ni,k) = vznm(ni,k) - tmp1(2) - tmp2(2)
        denom  =  a0nm(nr)*a0nm(nr)   + a0nm(ni)*a0nm(ni)
        tmp(1) = (a0nm(nr)*vznm(nr,k) + a0nm(ni)*vznm(ni,k))/denom
        tmp(2) = (a0nm(nr)*vznm(ni,k) - a0nm(ni)*vznm(nr,k))/denom
        vznm(nr,k) = tmp(1)
        vznm(ni,k) = tmp(2)
     end do
  end do
!
  return
end subroutine vertnm
