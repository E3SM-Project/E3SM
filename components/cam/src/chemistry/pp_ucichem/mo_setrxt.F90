
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,27) = 1.800E-12_r8
      rate(:,:,40) = 0.05E+00_r8
      rate(:,:,41) = 0.03E+00_r8
      rate(:,:,42) = 0.10E+00_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,17) = 1.630E-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,18) = 1.200E-10_r8 * exp_fac(:,:)
      rate(:,:,19) = 1.750E-10_r8 * exp_fac(:,:)
      rate(:,:,21) = 1.800E-12_r8 * exp_fac(:,:)
      rate(:,:,27) = 1.800E-12_r8 * exp_fac(:,:)
      rate(:,:,31) = 2.20e-11_r8 * exp_fac(:,:)
      rate(:,:,40) = 0.05E+00_r8 * exp_fac(:,:)
      rate(:,:,41) = 0.03E+00_r8 * exp_fac(:,:)
      rate(:,:,42) = 0.10E+00_r8 * exp_fac(:,:)
      rate(:,:,49) = 1.50e-13_r8 * exp_fac(:,:)
      rate(:,:,62) = 6.800E-14_r8 * exp_fac(:,:)
      rate(:,:,63) = 2.500E-14_r8 * exp_fac(:,:)
      rate(:,:,66) = 8.010E-12_r8 * exp_fac(:,:)
      rate(:,:,72) = 1.33e-13_r8 * exp_fac(:,:)
      rate(:,:,76) = 1.00E-13_r8 * exp_fac(:,:)
      rate(:,:,79) = 1.00E-13_r8 * exp_fac(:,:)
      rate(:,:,20) = 1.700E-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,22) = 1.000E-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,23) = 4.800E-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:,24) = 2.100E-33_r8 * exp( -920._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 460._r8 * itemp(:,:) )
      rate(:,:,25) = 3.000E-13_r8 * exp_fac(:,:)
      rate(:,:,45) = 2.400E-14_r8 * exp_fac(:,:)
      rate(:,:,28) = 3.000E-12_r8 * exp( -1500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,29) = 3.300E-12_r8 * exp_fac(:,:)
      rate(:,:,69) = 8.10E-12_r8 * exp_fac(:,:)
      rate(:,:,30) = 1.200E-13_r8 * exp( -2450._r8 * itemp(:,:) )
      rate(:,:,32) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,37) = 2.10E-27_r8 * exp( 10900._r8 * itemp(:,:) )
      rate(:,:,38) = 5.80E-27_r8 * exp( 10840._r8 * itemp(:,:) )
      rate(:,:,39) = 9.00E-29_r8 * exp( 14000._r8 * itemp(:,:) )
      rate(:,:,43) = 1.300E-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,44) = 2.700E-17_r8 * exp( 2199._r8 * itemp(:,:) )
      rate(:,:,46) = 2.800E-11_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,47) = 2.450E-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,50) = 7.660E-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,51) = 8.700E-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,53) = 1.200E-14_r8 * exp( -2630._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,54) = 1.10e-14_r8 * exp_fac(:,:)
      rate(:,:,73) = 3.80e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,55) = 3.000E-11_r8 * exp_fac(:,:)
      rate(:,:,74) = 2.70E-12_r8 * exp_fac(:,:)
      rate(:,:,77) = 2.70E-12_r8 * exp_fac(:,:)
      rate(:,:,82) = 2.70E-12_r8 * exp_fac(:,:)
      rate(:,:,56) = 5.500E-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,57) = 4.100E-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,58) = 2.800E-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,59) = 9.500E-14_r8 * exp( 390._r8 * itemp(:,:) )
      rate(:,:,60) = 3.800E-12_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:,61) = 2.600E-12_r8 * exp( 365._r8 * itemp(:,:) )
      rate(:,:,64) = 7.500E-13_r8 * exp( 700._r8 * itemp(:,:) )
      rate(:,:,65) = 1.900E-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,67) = 4.63E-12_r8 * exp( 350._r8 * itemp(:,:) )
      rate(:,:,68) = 1.40E-12_r8 * exp( -1900._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,70) = 2.00E-12_r8 * exp_fac(:,:)
      rate(:,:,71) = 2.90E-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1300._r8 * itemp(:,:) )
      rate(:,:,75) = 1.50E-13_r8 * exp_fac(:,:)
      rate(:,:,78) = 2.05E-13_r8 * exp_fac(:,:)
      rate(:,:,83) = 1.82E-13_r8 * exp_fac(:,:)
      rate(:,:,80) = 8.50E-16_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,81) = 2.60E-12_r8 * exp( 610._r8 * itemp(:,:) )
      rate(:,:,84) = 9.600E-12_r8 * exp( -234._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,26), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.80e-11_r8
      call jpl( rate(1,1,33), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.90e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 4.00e-12_r8 * itemp(:,:)**0.3_r8
      call jpl( rate(1,1,34), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.40E-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 1.60E-12_r8 * itemp(:,:)**(-0.1_r8)
      call jpl( rate(1,1,35), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.70E-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30E-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,36), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90E-33_r8 * itemp(:,:)**(-1.0_r8)
      kinf(:,:) = 1.10E-12_r8 * itemp(:,:)**1.30_r8
      call jpl( rate(1,1,48), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.10E-28_r8 * itemp(:,:)**(-3.50_r8)
      kinf(:,:) = 8.40E-12_r8 * itemp(:,:)**(-1.75_r8)
      call jpl( rate(1,1,52), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
