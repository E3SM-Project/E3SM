      module mo_indprd
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: indprd
      contains
      subroutine indprd( class, prod, nprod, y, extfrc, rxt, ncol )
      use chem_mods, only : gas_pcnst, extcnt, rxntot
      use ppgrid, only : pver
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(ncol,pver,gas_pcnst)
      real(r8), intent(in) :: rxt(ncol,pver,rxntot)
      real(r8), intent(in) :: extfrc(ncol,pver,extcnt)
      real(r8), intent(inout) :: prod(ncol,pver,nprod)
!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,:,1) = (rxt(:,:,4) +rxt(:,:,5) +rxt(:,:,17)*y(:,:,2))*y(:,:,9) &
                  +.050_r8*rxt(:,:,30)*y(:,:,15)*y(:,:,1) + extfrc(:,:,3)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,10) = 0._r8
         prod(:,:,13) = 0._r8
         prod(:,:,11) = 0._r8
         prod(:,:,4) = 0._r8
         prod(:,:,12) = + extfrc(:,:,1)
         prod(:,:,9) = + extfrc(:,:,2)
         prod(:,:,1) = 0._r8
         prod(:,:,7) = 0._r8
         prod(:,:,14) = 0._r8
         prod(:,:,8) = 0._r8
         prod(:,:,3) = 0._r8
         prod(:,:,5) = + extfrc(:,:,4)
         prod(:,:,2) = + extfrc(:,:,5)
         prod(:,:,6) = 0._r8
      end if
      end subroutine indprd
      end module mo_indprd
