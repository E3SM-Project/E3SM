













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
         prod(:,:,1) = (rxt(:,:,4) +rxt(:,:,5) +rxt(:,:,17)*y(:,:,2))*y(:,:,9) &
                  +.050_r8*rxt(:,:,28)*y(:,:,14)*y(:,:,1)
                                                                                          
!--------------------------------------------------------------------
!       ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,40) = 0._r8
                                                                                          
         prod(:,:,38) = 0._r8
                                                                                          
         prod(:,:,39) = 0._r8
                                                                                          
         prod(:,:,31) = 0._r8
                                                                                          
         prod(:,:,36) = + extfrc(:,:,1)
                                                                                          
         prod(:,:,35) = + extfrc(:,:,2)
                                                                                          
         prod(:,:,1) = 0._r8
                                                                                          
         prod(:,:,33) = 0._r8
                                                                                          
         prod(:,:,37) = 0._r8
                                                                                          
         prod(:,:,34) = 0._r8
                                                                                          
         prod(:,:,30) = 0._r8
                                                                                          
         prod(:,:,29) = + extfrc(:,:,3)
                                                                                          
         prod(:,:,32) = 0._r8
                                                                                          
         prod(:,:,2) = 0._r8
                                                                                          
         prod(:,:,3) = + extfrc(:,:,11)
                                                                                          
         prod(:,:,4) = + extfrc(:,:,4)
                                                                                          
         prod(:,:,5) = 0._r8
                                                                                          
         prod(:,:,6) = 0._r8
                                                                                          
         prod(:,:,7) = 0._r8
                                                                                          
         prod(:,:,8) = 0._r8
                                                                                          
         prod(:,:,9) = 0._r8
                                                                                          
         prod(:,:,10) = 0._r8
                                                                                          
         prod(:,:,11) = + extfrc(:,:,8)
                                                                                          
         prod(:,:,12) = + extfrc(:,:,5)
                                                                                          
         prod(:,:,13) = 0._r8
                                                                                          
         prod(:,:,14) = 0._r8
                                                                                          
         prod(:,:,15) = 0._r8
                                                                                          
         prod(:,:,16) = + extfrc(:,:,9)
                                                                                          
         prod(:,:,17) = 0._r8
                                                                                          
         prod(:,:,18) = 0._r8
                                                                                          
         prod(:,:,19) = 0._r8
                                                                                          
         prod(:,:,20) = 0._r8
                                                                                          
         prod(:,:,21) = 0._r8
                                                                                          
         prod(:,:,22) = 0._r8
                                                                                          
         prod(:,:,23) = 0._r8
                                                                                          
         prod(:,:,24) = 0._r8
                                                                                          
         prod(:,:,25) = + extfrc(:,:,6)
                                                                                          
         prod(:,:,26) = + extfrc(:,:,7)
                                                                                          
         prod(:,:,27) = 0._r8
                                                                                          
         prod(:,:,28) = + extfrc(:,:,10)
                                                                                          
      end if                                                                              
                                                                                          
      end subroutine indprd                                                               
                                                                                          
      end module mo_indprd                                                                
