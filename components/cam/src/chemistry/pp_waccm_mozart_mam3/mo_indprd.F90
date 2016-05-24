
      module mo_indprd

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: indprd

      contains

      subroutine indprd( class, prod, nprod, y, extfrc, rxt, ncol )

      use chem_mods, only : gas_pcnst, extcnt, rxntot
      use ppgrid,    only : pver

      implicit none

!--------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      real(r8), intent(in)    :: y(ncol,pver,gas_pcnst)
      real(r8), intent(in)    :: rxt(ncol,pver,rxntot)
      real(r8), intent(in)    :: extfrc(ncol,pver,extcnt)
      real(r8), intent(inout) :: prod(ncol,pver,nprod)

!--------------------------------------------------------------------
!       ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,:,1) = 0._r8
                                                                                          
         prod(:,:,2) =rxt(:,:,113)*y(:,:,10)*y(:,:,8)
                                                                                          
         prod(:,:,3) = (rxt(:,:,17) +rxt(:,:,18) +rxt(:,:,135)*y(:,:,11) + &
                 rxt(:,:,136)*y(:,:,22) +rxt(:,:,137)*y(:,:,2) + &
                 rxt(:,:,161)*y(:,:,27) +rxt(:,:,184)*y(:,:,36))*y(:,:,18) &
                  + extfrc(:,:,3)
                                                                                          
         prod(:,:,4) =rxt(:,:,18)*y(:,:,18) +rxt(:,:,143)*y(:,:,23)*y(:,:,21) &
                  +rxt(:,:,20)*y(:,:,78)
                                                                                          
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
!       ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,59) = 0._r8
                                                                                          
         prod(:,:,61) = (rxt(:,:,46) +rxt(:,:,74))*y(:,:,51) +.180_r8*rxt(:,:,48) &
                 *y(:,:,15)
                                                                                          
         prod(:,:,58) =rxt(:,:,5)*y(:,:,7)
                                                                                          
         prod(:,:,48) = 0._r8
                                                                                          
         prod(:,:,22) = 0._r8
                                                                                          
         prod(:,:,21) = 0._r8
                                                                                          
         prod(:,:,40) = (rxt(:,:,58) +.800_r8*rxt(:,:,61) +rxt(:,:,70) + &
                 .800_r8*rxt(:,:,73)) + extfrc(:,:,16)
                                                                                          
         prod(:,:,52) = + extfrc(:,:,1)
                                                                                          
         prod(:,:,53) = + extfrc(:,:,2)
                                                                                          
         prod(:,:,54) =.660_r8*rxt(:,:,48)*y(:,:,15) + extfrc(:,:,18)
                                                                                          
         prod(:,:,55) = 0._r8
                                                                                          
         prod(:,:,36) = 0._r8
                                                                                          
         prod(:,:,29) = 0._r8
                                                                                          
         prod(:,:,26) = 0._r8
                                                                                          
         prod(:,:,42) =rxt(:,:,47)*y(:,:,15) +rxt(:,:,36)*y(:,:,41) +rxt(:,:,43) &
                 *y(:,:,42)
                                                                                          
         prod(:,:,27) = 0._r8
                                                                                          
         prod(:,:,46) =.180_r8*rxt(:,:,48)*y(:,:,15)
                                                                                          
         prod(:,:,47) =rxt(:,:,47)*y(:,:,15)
                                                                                          
         prod(:,:,51) = 0._r8
                                                                                          
         prod(:,:,34) = 0._r8
                                                                                          
         prod(:,:,45) =.050_r8*rxt(:,:,48)*y(:,:,15)
                                                                                          
         prod(:,:,56) =rxt(:,:,36)*y(:,:,41) +3.000_r8*rxt(:,:,39)*y(:,:,43) &
                  +2.000_r8*rxt(:,:,40)*y(:,:,44) +3.000_r8*rxt(:,:,41)*y(:,:,45) &
                  +rxt(:,:,42)*y(:,:,46) +4.000_r8*rxt(:,:,37)*y(:,:,47) &
                  +3.000_r8*rxt(:,:,38)*y(:,:,48) +rxt(:,:,45)*y(:,:,50)
                                                                                          
         prod(:,:,23) = 0._r8
                                                                                          
         prod(:,:,50) = 0._r8
                                                                                          
         prod(:,:,20) = 0._r8
                                                                                          
         prod(:,:,18) = 0._r8
                                                                                          
         prod(:,:,60) = 0._r8
                                                                                          
         prod(:,:,43) = 0._r8
                                                                                          
         prod(:,:,44) = 0._r8
                                                                                          
         prod(:,:,24) = 0._r8
                                                                                          
         prod(:,:,57) =rxt(:,:,43)*y(:,:,42) +rxt(:,:,44)*y(:,:,49) +rxt(:,:,45) &
                 *y(:,:,50)
                                                                                          
         prod(:,:,49) = 0._r8
                                                                                          
         prod(:,:,35) = 0._r8
                                                                                          
         prod(:,:,41) = 0._r8
                                                                                          
         prod(:,:,32) = 0._r8
                                                                                          
         prod(:,:,33) = (rxt(:,:,54) +rxt(:,:,66)) + extfrc(:,:,14)
                                                                                          
         prod(:,:,37) = + extfrc(:,:,12)
                                                                                          
         prod(:,:,28) = (rxt(:,:,58) +rxt(:,:,59) +rxt(:,:,70) +rxt(:,:,71)) &
                  + extfrc(:,:,13)
                                                                                          
         prod(:,:,30) = + extfrc(:,:,11)
                                                                                          
         prod(:,:,38) = 0._r8
                                                                                          
         prod(:,:,31) = (rxt(:,:,59) +1.200_r8*rxt(:,:,61) +rxt(:,:,71) + &
                 1.200_r8*rxt(:,:,73)) + extfrc(:,:,15)
                                                                                          
         prod(:,:,39) = (rxt(:,:,54) +rxt(:,:,58) +rxt(:,:,59) +rxt(:,:,66) + &
                 rxt(:,:,70) +rxt(:,:,71)) + extfrc(:,:,17)
                                                                                          
         prod(:,:,19) = + extfrc(:,:,4)
                                                                                          
         prod(:,:,25) = 0._r8
                                                                                          
         prod(:,:,1) = 0._r8
                                                                                          
         prod(:,:,2) = 0._r8
                                                                                          
         prod(:,:,3) = + extfrc(:,:,5)
                                                                                          
         prod(:,:,4) = + extfrc(:,:,7)
                                                                                          
         prod(:,:,5) = 0._r8
                                                                                          
         prod(:,:,6) = + extfrc(:,:,8)
                                                                                          
         prod(:,:,7) = 0._r8
                                                                                          
         prod(:,:,8) = 0._r8
                                                                                          
         prod(:,:,9) = + extfrc(:,:,9)
                                                                                          
         prod(:,:,10) = + extfrc(:,:,6)
                                                                                          
         prod(:,:,11) = 0._r8
                                                                                          
         prod(:,:,12) = 0._r8
                                                                                          
         prod(:,:,13) = + extfrc(:,:,10)
                                                                                          
         prod(:,:,14) = 0._r8
                                                                                          
         prod(:,:,15) = 0._r8
                                                                                          
         prod(:,:,16) = 0._r8
                                                                                          
         prod(:,:,17) = 0._r8
                                                                                          
      end if                                                                              
                                                                                          
      end subroutine indprd                                                               
                                                                                          
      end module mo_indprd                                                                
