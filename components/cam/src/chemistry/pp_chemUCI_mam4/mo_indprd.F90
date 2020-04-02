













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
         prod(:,:,1) = (rxt(:,:,4) +rxt(:,:,5) +rxt(:,:,30)*y(:,:,2))*y(:,:,5) &
                  + (2.000_r8*rxt(:,:,74)*y(:,:,27) +rxt(:,:,79)*y(:,:,29))*y(:,:,3) &
                  + (rxt(:,:,13) +.500_r8*rxt(:,:,63)*y(:,:,2))*y(:,:,24) &
                  +rxt(:,:,27)*y(:,:,19)*y(:,:,1) +1.500_r8*rxt(:,:,16)*y(:,:,28)
                                                                                          
         prod(:,:,2) = 0._r8
                                                                                          
         prod(:,:,3) = 0._r8
                                                                                          
         prod(:,:,4) = 0._r8
                                                                                          
         prod(:,:,5) = 0._r8
                                                                                          
         prod(:,:,6) = 0._r8
                                                                                          
         prod(:,:,7) = 0._r8
                                                                                          
!--------------------------------------------------------------------
!       ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,53) =2.000_r8*rxt(:,:,2)
                                                                                          
         prod(:,:,52) = 0._r8
                                                                                          
         prod(:,:,51) = 0._r8
                                                                                          
         prod(:,:,30) = 0._r8
                                                                                          
         prod(:,:,39) = 0._r8
                                                                                          
         prod(:,:,48) =1.330_r8*rxt(:,:,15)*y(:,:,21)
                                                                                          
         prod(:,:,37) = 0._r8
                                                                                          
         prod(:,:,47) = + extfrc(:,:,1)
                                                                                          
         prod(:,:,49) = + extfrc(:,:,2)
                                                                                          
         prod(:,:,50) = 0._r8
                                                                                          
         prod(:,:,32) = 0._r8
                                                                                          
         prod(:,:,35) = 0._r8
                                                                                          
         prod(:,:,36) = 0._r8
                                                                                          
         prod(:,:,31) = 0._r8
                                                                                          
         prod(:,:,42) = 0._r8
                                                                                          
         prod(:,:,38) = 0._r8
                                                                                          
         prod(:,:,45) = 0._r8
                                                                                          
         prod(:,:,46) =.670_r8*rxt(:,:,15)*y(:,:,21)
                                                                                          
         prod(:,:,33) = 0._r8
                                                                                          
         prod(:,:,41) = 0._r8
                                                                                          
         prod(:,:,34) = 0._r8
                                                                                          
         prod(:,:,43) = 0._r8
                                                                                          
         prod(:,:,44) = 0._r8
                                                                                          
         prod(:,:,40) = 0._r8
                                                                                          
         prod(:,:,29) = 0._r8
                                                                                          
         prod(:,:,28) = + extfrc(:,:,3)
                                                                                          
         prod(:,:,1) = 0._r8
                                                                                          
         prod(:,:,2) = 0._r8
                                                                                          
         prod(:,:,3) = + extfrc(:,:,4)
                                                                                          
         prod(:,:,4) = + extfrc(:,:,5)
                                                                                          
         prod(:,:,5) = 0._r8
                                                                                          
         prod(:,:,6) = + extfrc(:,:,6)
                                                                                          
         prod(:,:,7) = 0._r8
                                                                                          
         prod(:,:,8) = 0._r8
                                                                                          
         prod(:,:,9) = 0._r8
                                                                                          
         prod(:,:,10) = 0._r8
                                                                                          
         prod(:,:,11) = 0._r8
                                                                                          
         prod(:,:,12) = + extfrc(:,:,7)
                                                                                          
         prod(:,:,13) = 0._r8
                                                                                          
         prod(:,:,14) = 0._r8
                                                                                          
         prod(:,:,15) = 0._r8
                                                                                          
         prod(:,:,16) = 0._r8
                                                                                          
         prod(:,:,17) = 0._r8
                                                                                          
         prod(:,:,18) = 0._r8
                                                                                          
         prod(:,:,19) = 0._r8
                                                                                          
         prod(:,:,20) = 0._r8
                                                                                          
         prod(:,:,21) = 0._r8
                                                                                          
         prod(:,:,22) = 0._r8
                                                                                          
         prod(:,:,23) = 0._r8
                                                                                          
         prod(:,:,24) = + extfrc(:,:,8)
                                                                                          
         prod(:,:,25) = + extfrc(:,:,9)
                                                                                          
         prod(:,:,26) = 0._r8
                                                                                          
         prod(:,:,27) = 0._r8
                                                                                          
      end if                                                                              
                                                                                          
      end subroutine indprd                                                               
                                                                                          
      end module mo_indprd                                                                
