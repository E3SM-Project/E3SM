subroutine hdinti(rearth, deltat)
!-----------------------------------------------------------------------
!
! Purpose:
! Time independent initialization for the horizontal diffusion.
!
! Author:  D. Williamson
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use pmgrid
  use pspect
  use sld_control_mod
  implicit none
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: rearth ! radius of the earth
  real(r8), intent(in)   :: deltat ! time step
!
!---------------------------Local workspace-----------------------------
!
  integer     k                    ! level index
  integer(i8) n                    ! n-wavenumber index
!
!-----------------------------------------------------------------------
!
! Top level for del**4 diffusion, set for 18-level model
!
  kmnhd4 = 5
!
! Bottom level for increased del**2 diffusion (kmxhd2 < kmnhd4)
!
  kmxhd2 = 3
!
! Initialize physical constants for courant number based spect truncation
!
  nmaxhd = ptrk
  cnlim  = 0.999_r8          ! maximum allowable Courant number
  cnfac  = deltat*real(nmaxhd,r8)/rearth
!
! Initialize arrays used for courant number based spectral truncation
!
  do k=1,plev
     nindex(k) = 2*nmaxhd
  end do
!
! Set the Del^2 and Del^4 diffusion coefficients for each wavenumber
!
  hdfst2(1) = 0._r8
  hdfsd2(1) = 0._r8
!
  hdfst4(1) = 0._r8
  hdfsd4(1) = 0._r8
  do n=2,pnmax
     hdfst2(n) = dif2 * (n*(n-1)  ) / rearth**2
     hdfsd2(n) = dif2 * (n*(n-1)-2) / rearth**2

     hdfst4(n) = dif4 * (n*(n-1)*n*(n-1)  ) / rearth**4
     hdfsd4(n) = dif4 * (n*(n-1)*n*(n-1)-4) / rearth**4
  end do
!
  return
end subroutine hdinti

