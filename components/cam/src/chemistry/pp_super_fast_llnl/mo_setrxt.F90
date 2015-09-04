
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

      rate(:,:,11) = 1.800E-12_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,7) = 1.700E-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,8) = 1.000E-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,9) = 4.800E-11_r8 * exp_fac(:,:)
      rate(:,:,13) = 3.500E-12_r8 * exp_fac(:,:)
      rate(:,:,12) = 3.000E-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:,15) = 2.450E-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,17) = 5.500E-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,18) = 4.100E-13_r8 * exp( 750._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,19) = 2.700E-12_r8 * exp_fac(:,:)
      rate(:,:,20) = 1.100E-12_r8 * exp_fac(:,:)
      rate(:,:,21) = 2.800E-12_r8 * exp( 300._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 390._r8 * itemp(:,:) )
      rate(:,:,22) = 9.500E-14_r8 * exp_fac(:,:)
      rate(:,:,29) = 2.700E-11_r8 * exp_fac(:,:)
      rate(:,:,24) = 1.100E-11_r8 * exp( -240._r8 * itemp(:,:) )
      rate(:,:,30) = 5.590E-15_r8 * exp( -1814._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 1.800E-30_r8 * itemp(:,:)**3.00_r8
      kinf(:,:) = 2.800E-11_r8
      call jpl( rate(1,1,14), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 3.300E-31_r8 * itemp(:,:)**4.30_r8
      kinf(:,:) = 1.600E-12_r8
      call jpl( rate(1,1,26), m, 0.6_r8, ko, kinf, n )

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
