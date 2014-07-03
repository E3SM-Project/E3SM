




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
         prod(:,:,1) =.080_r8*rxt(:,:,116)*y(:,:,33)*y(:,:,1)

         prod(:,:,2) = 0._r8

         prod(:,:,3) = (.500_r8*rxt(:,:,91)*y(:,:,22) + &
                 .560_r8*rxt(:,:,116)*y(:,:,33) +.300_r8*rxt(:,:,147)*y(:,:,50) + &
                 .050_r8*rxt(:,:,159)*y(:,:,53) +.200_r8*rxt(:,:,161)*y(:,:,54)) &
                 *y(:,:,1) + (rxt(:,:,79)*y(:,:,19) +rxt(:,:,133)*y(:,:,60) + &
                 .220_r8*rxt(:,:,164)*y(:,:,55) +.500_r8*rxt(:,:,181)*y(:,:,61)) &
                 *y(:,:,7) + (rxt(:,:,80)*y(:,:,19) +.350_r8*rxt(:,:,89)*y(:,:,101) + &
                 rxt(:,:,112)*y(:,:,72) +rxt(:,:,132)*y(:,:,60))*y(:,:,12) &
                  + (.220_r8*rxt(:,:,162)*y(:,:,5) +.110_r8*rxt(:,:,166)*y(:,:,17) + &
                 .220_r8*rxt(:,:,167)*y(:,:,31))*y(:,:,55) &
                  + (.500_r8*rxt(:,:,180)*y(:,:,5) +.200_r8*rxt(:,:,183)*y(:,:,17) + &
                 .500_r8*rxt(:,:,184)*y(:,:,31))*y(:,:,61) + (rxt(:,:,13) + &
                 rxt(:,:,14))*y(:,:,19) + (.670_r8*rxt(:,:,21) +.670_r8*rxt(:,:,22)) &
                 *y(:,:,54) +rxt(:,:,33)*y(:,:,26) +rxt(:,:,16)*y(:,:,30) &
                  +.700_r8*rxt(:,:,23)*y(:,:,53) +rxt(:,:,28)*y(:,:,60) &
                  +.450_r8*rxt(:,:,35)*y(:,:,71) +2.000_r8*rxt(:,:,36)*y(:,:,72) &
                  +rxt(:,:,30)*y(:,:,77) + extfrc(:,:,2)

         prod(:,:,4) = 0._r8

         prod(:,:,5) = 0._r8

         prod(:,:,6) =rxt(:,:,14)*y(:,:,19)

         prod(:,:,7) = 0._r8

         prod(:,:,8) = 0._r8

!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,91) = 0._r8

         prod(:,:,30) =rxt(:,:,4)*y(:,:,4)

         prod(:,:,59) =2.000_r8*rxt(:,:,1)

         prod(:,:,88) = + extfrc(:,:,1)

         prod(:,:,93) = 0._r8

         prod(:,:,94) = 0._r8

         prod(:,:,58) = 0._r8

         prod(:,:,41) = 0._r8

         prod(:,:,31) = 0._r8

         prod(:,:,89) = 0._r8

         prod(:,:,90) = 0._r8

         prod(:,:,26) = 0._r8

         prod(:,:,92) = 0._r8

         prod(:,:,43) = 0._r8

         prod(:,:,82) = 0._r8

         prod(:,:,62) = 0._r8

         prod(:,:,33) = 0._r8

         prod(:,:,38) = 0._r8

         prod(:,:,27) = 0._r8

         prod(:,:,45) = 0._r8

         prod(:,:,50) = 0._r8

         prod(:,:,63) = 0._r8

         prod(:,:,18) = 0._r8

         prod(:,:,71) = 0._r8

         prod(:,:,42) = 0._r8

         prod(:,:,76) = 0._r8

         prod(:,:,95) = 0._r8

         prod(:,:,49) = 0._r8

         prod(:,:,75) = 0._r8

         prod(:,:,19) = 0._r8

         prod(:,:,72) = 0._r8

         prod(:,:,40) = 0._r8

         prod(:,:,67) = 0._r8

         prod(:,:,52) = 0._r8

         prod(:,:,68) = 0._r8

         prod(:,:,78) = 0._r8

         prod(:,:,34) = 0._r8

         prod(:,:,20) = 0._r8

         prod(:,:,35) = 0._r8

         prod(:,:,21) = 0._r8

         prod(:,:,70) = 0._r8

         prod(:,:,65) = 0._r8

         prod(:,:,51) = 0._r8

         prod(:,:,66) = 0._r8

         prod(:,:,36) = 0._r8

         prod(:,:,74) = 0._r8

         prod(:,:,86) = 0._r8

         prod(:,:,61) = 0._r8

         prod(:,:,84) = 0._r8

         prod(:,:,80) = 0._r8

         prod(:,:,85) = 0._r8

         prod(:,:,37) = 0._r8

         prod(:,:,87) = 0._r8

         prod(:,:,44) = 0._r8

         prod(:,:,79) = 0._r8

         prod(:,:,81) = 0._r8

         prod(:,:,83) = 0._r8

         prod(:,:,28) = 0._r8

         prod(:,:,64) = 0._r8

         prod(:,:,77) = 0._r8

         prod(:,:,53) = 0._r8

         prod(:,:,23) = 0._r8

         prod(:,:,24) = 0._r8

         prod(:,:,55) = 0._r8

         prod(:,:,46) = 0._r8

         prod(:,:,29) = 0._r8

         prod(:,:,56) = 0._r8

         prod(:,:,54) = 0._r8

         prod(:,:,57) = 0._r8

         prod(:,:,47) = 0._r8

         prod(:,:,60) = 0._r8

         prod(:,:,69) = 0._r8

         prod(:,:,73) = 0._r8

         prod(:,:,22) = + extfrc(:,:,3)

         prod(:,:,32) = 0._r8

         prod(:,:,1) = 0._r8

         prod(:,:,17) = 0._r8

         prod(:,:,2) = 0._r8

         prod(:,:,3) = 0._r8

         prod(:,:,4) = 0._r8

         prod(:,:,5) = 0._r8

         prod(:,:,6) = 0._r8

         prod(:,:,7) = 0._r8

         prod(:,:,8) = 0._r8

         prod(:,:,25) = 0._r8

         prod(:,:,39) = 0._r8

         prod(:,:,48) = 0._r8

         prod(:,:,9) = 0._r8

         prod(:,:,10) = 0._r8

         prod(:,:,11) = 0._r8

         prod(:,:,12) = 0._r8

         prod(:,:,13) = 0._r8

         prod(:,:,14) = 0._r8

         prod(:,:,15) = 0._r8

         prod(:,:,16) = 0._r8

      end if

      end subroutine indprd

      end module mo_indprd
