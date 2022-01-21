
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

      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,20) = 1.630E-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,21) = 1.200E-10_r8 * exp_fac(:,:)
      rate(:,:,22) = 1.750E-10_r8 * exp_fac(:,:)
      rate(:,:,23) = 1.500E-13_r8 * exp_fac(:,:)
      rate(:,:,35) = 1.800E-12_r8 * exp_fac(:,:)
      rate(:,:,41) = 1.800E-12_r8 * exp_fac(:,:)
      rate(:,:,45) = 2.200e-11_r8 * exp_fac(:,:)
      rate(:,:,62) = 6.800E-14_r8 * exp_fac(:,:)
      rate(:,:,63) = 2.500E-14_r8 * exp_fac(:,:)
      rate(:,:,66) = 8.010E-12_r8 * exp_fac(:,:)
      rate(:,:,72) = 1.330e-13_r8 * exp_fac(:,:)
      rate(:,:,76) = 1.000E-13_r8 * exp_fac(:,:)
      rate(:,:,79) = 1.000E-13_r8 * exp_fac(:,:)
      rate(:,:,84) = 1.286E-7_r8 * exp_fac(:,:)
      rate(:,:,25) = 2.800E-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,26) = 2.450E-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,27) = 7.660E-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,28) = 8.700E-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,30) = 1.200E-14_r8 * exp( -2630._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,31) = 1.100e-14_r8 * exp_fac(:,:)
      rate(:,:,73) = 3.820e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,32) = 3.000E-11_r8 * exp_fac(:,:)
      rate(:,:,74) = 2.700E-12_r8 * exp_fac(:,:)
      rate(:,:,77) = 2.700E-12_r8 * exp_fac(:,:)
      rate(:,:,82) = 2.700E-12_r8 * exp_fac(:,:)
      rate(:,:,33) = 5.500E-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,34) = 1.700E-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,36) = 1.000E-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,37) = 4.800E-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:,38) = 2.100E-33_r8 * exp( 920._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 460._r8 * itemp(:,:) )
      rate(:,:,39) = 3.000E-13_r8 * exp_fac(:,:)
      rate(:,:,48) = 2.400E-14_r8 * exp_fac(:,:)
      rate(:,:,42) = 3.000E-12_r8 * exp( -1500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,43) = 3.300E-12_r8 * exp_fac(:,:)
      rate(:,:,69) = 8.100E-12_r8 * exp_fac(:,:)
      rate(:,:,44) = 1.200E-13_r8 * exp( -2450._r8 * itemp(:,:) )
      rate(:,:,46) = 1.500e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,47) = 1.300E-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,54) = 2.10E-27_r8 * exp( 10900._r8 * itemp(:,:) )
      rate(:,:,55) = 5.80E-27_r8 * exp( 10840._r8 * itemp(:,:) )
      rate(:,:,56) = 9.00E-29_r8 * exp( 14000._r8 * itemp(:,:) )
      rate(:,:,57) = 4.100E-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,58) = 2.800E-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,59) = 9.500E-14_r8 * exp( 390._r8 * itemp(:,:) )
      rate(:,:,60) = 3.800E-12_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:,61) = 2.600E-12_r8 * exp( 365._r8 * itemp(:,:) )
      rate(:,:,64) = 7.500E-13_r8 * exp( 700._r8 * itemp(:,:) )
      rate(:,:,65) = 1.900E-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,67) = 4.630E-12_r8 * exp( 350._r8 * itemp(:,:) )
      rate(:,:,68) = 1.400E-12_r8 * exp( -1900._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,70) = 2.000E-12_r8 * exp_fac(:,:)
      rate(:,:,71) = 2.900E-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1300._r8 * itemp(:,:) )
      rate(:,:,75) = 1.500E-13_r8 * exp_fac(:,:)
      rate(:,:,78) = 2.050E-13_r8 * exp_fac(:,:)
      rate(:,:,83) = 1.820E-13_r8 * exp_fac(:,:)
      rate(:,:,80) = 8.500E-16_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,81) = 2.600E-12_r8 * exp( 610._r8 * itemp(:,:) )
      rate(:,:,85) = 9.600E-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,88) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 1.10E-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,24), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.10E-28_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 8.40E-12_r8 * itemp(:,:)**1.75_r8
      call jpl( rate(1,1,29), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,40), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.80e-11_r8
      call jpl( rate(1,1,50), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.90e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 4.00e-12_r8 * itemp(:,:)**0.3_r8
      call jpl( rate(1,1,51), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.40E-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 1.60E-12_r8 * itemp(:,:)**(-0.1_r8)
      call jpl( rate(1,1,52), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.70E-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30E-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,53), m, 0.6_r8, ko, kinf, n )

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
