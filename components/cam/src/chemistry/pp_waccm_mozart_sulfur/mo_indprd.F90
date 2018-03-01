




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
         prod(:,:,1) = 0._r8

         prod(:,:,2) =rxt(:,:,118)*y(:,:,10)*y(:,:,8)

         prod(:,:,3) = (rxt(:,:,17) +rxt(:,:,18) +rxt(:,:,140)*y(:,:,11) + &
                 rxt(:,:,141)*y(:,:,22) +rxt(:,:,142)*y(:,:,2) + &
                 rxt(:,:,166)*y(:,:,27) +rxt(:,:,189)*y(:,:,36))*y(:,:,18) &
                  + (rxt(:,:,52) +rxt(:,:,208)*y(:,:,2))*y(:,:,59) + extfrc(:,:,3)

         prod(:,:,4) =rxt(:,:,18)*y(:,:,18) +rxt(:,:,148)*y(:,:,23)*y(:,:,21) &
                  +rxt(:,:,20)*y(:,:,66)

         prod(:,:,5) = 0._r8

         prod(:,:,6) = 0._r8

         prod(:,:,7) = 0._r8

         prod(:,:,8) = 0._r8

         prod(:,:,9) = 0._r8

         prod(:,:,10) = 0._r8

         prod(:,:,11) = 0._r8

         prod(:,:,12) = 0._r8

         prod(:,:,13) = 0._r8

         prod(:,:,14) = 0._r8

         prod(:,:,15) = 0._r8

         prod(:,:,16) = 0._r8

         prod(:,:,17) = 0._r8

!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,42) = 0._r8

         prod(:,:,39) = (rxt(:,:,46) +rxt(:,:,79))*y(:,:,51) +.180_r8*rxt(:,:,48) &
                 *y(:,:,15)

         prod(:,:,45) =rxt(:,:,5)*y(:,:,7)

         prod(:,:,49) = 0._r8

         prod(:,:,3) = 0._r8

         prod(:,:,2) = 0._r8

         prod(:,:,27) = (rxt(:,:,63) +.800_r8*rxt(:,:,66) +rxt(:,:,75) + &
                 .800_r8*rxt(:,:,78)) + extfrc(:,:,10)

         prod(:,:,38) = + extfrc(:,:,1)

         prod(:,:,40) = + extfrc(:,:,2)

         prod(:,:,43) =.660_r8*rxt(:,:,48)*y(:,:,15) + extfrc(:,:,12)

         prod(:,:,33) = 0._r8

         prod(:,:,22) = 0._r8

         prod(:,:,13) = 0._r8

         prod(:,:,8) = 0._r8

         prod(:,:,28) =rxt(:,:,47)*y(:,:,15) +rxt(:,:,36)*y(:,:,41) +rxt(:,:,43) &
                 *y(:,:,42)

         prod(:,:,11) = 0._r8

         prod(:,:,46) =.180_r8*rxt(:,:,48)*y(:,:,15)

         prod(:,:,32) =rxt(:,:,47)*y(:,:,15)

         prod(:,:,47) = 0._r8

         prod(:,:,18) = 0._r8

         prod(:,:,44) =.050_r8*rxt(:,:,48)*y(:,:,15)

         prod(:,:,48) =rxt(:,:,36)*y(:,:,41) +3.000_r8*rxt(:,:,39)*y(:,:,43) &
                  +2.000_r8*rxt(:,:,40)*y(:,:,44) +3.000_r8*rxt(:,:,41)*y(:,:,45) &
                  +rxt(:,:,42)*y(:,:,46) +4.000_r8*rxt(:,:,37)*y(:,:,47) &
                  +3.000_r8*rxt(:,:,38)*y(:,:,48) +rxt(:,:,45)*y(:,:,50)

         prod(:,:,5) = 0._r8

         prod(:,:,37) = 0._r8

         prod(:,:,9) = 0._r8

         prod(:,:,1) = 0._r8

         prod(:,:,36) = 0._r8

         prod(:,:,29) = 0._r8

         prod(:,:,30) = 0._r8

         prod(:,:,6) = 0._r8

         prod(:,:,35) =rxt(:,:,43)*y(:,:,42) +rxt(:,:,44)*y(:,:,49) +rxt(:,:,45) &
                 *y(:,:,50)

         prod(:,:,41) = 0._r8

         prod(:,:,20) = 0._r8

         prod(:,:,23) = 0._r8

         prod(:,:,17) = 0._r8

         prod(:,:,19) = (rxt(:,:,59) +rxt(:,:,71)) + extfrc(:,:,8)

         prod(:,:,24) = + extfrc(:,:,6)

         prod(:,:,14) = (rxt(:,:,63) +rxt(:,:,64) +rxt(:,:,75) +rxt(:,:,76)) &
                  + extfrc(:,:,7)

         prod(:,:,15) = + extfrc(:,:,5)

         prod(:,:,25) = 0._r8

         prod(:,:,16) = (rxt(:,:,64) +1.200_r8*rxt(:,:,66) +rxt(:,:,76) + &
                 1.200_r8*rxt(:,:,78)) + extfrc(:,:,9)

         prod(:,:,26) = (rxt(:,:,59) +rxt(:,:,63) +rxt(:,:,64) +rxt(:,:,71) + &
                 rxt(:,:,75) +rxt(:,:,76)) + extfrc(:,:,11)

         prod(:,:,10) = 0._r8

         prod(:,:,21) = 0._r8

         prod(:,:,34) = 0._r8

         prod(:,:,31) = + extfrc(:,:,4)

         prod(:,:,7) = 0._r8

         prod(:,:,12) = 0._r8

         prod(:,:,4) = 0._r8

      end if

      end subroutine indprd

      end module mo_indprd
