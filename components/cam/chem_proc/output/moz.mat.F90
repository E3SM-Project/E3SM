












      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid,       only : pver

      implicit none

!--------------------------------------------------------------------
!     ... dummy args                                                                      
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in)    ::  y(:,:,:)
      real(r8), intent(in)    ::  rxt(:,:,:)
      real(r8), intent(in)    ::  het_rates(:,:,:)



!--------------------------------------------------------------------
!       ... loss and production for Explicit method
!--------------------------------------------------------------------


         loss(:,:,1) = ( + het_rates(:,:,1))* y(:,:,1)
         prod(:,:,1) = 0._r8

      end subroutine exp_prod_loss

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid,       only : pver

      implicit none

!--------------------------------------------------------------------
!     ... dummy args                                                                      
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in)    ::  y(:)
      real(r8), intent(in)    ::  rxt(:)
      real(r8), intent(in)    ::  het_rates(:)



!--------------------------------------------------------------------
!       ... loss and production for Implicit method
!--------------------------------------------------------------------


         loss(1) = ( + rxt(1) + rxt(3) + het_rates(2))* y(2)
         prod(1) = 0._r8
         loss(2) = ( + het_rates(3))* y(3)
         prod(2) =rxt(4)*y(4)
         loss(3) = ( + rxt(4) + het_rates(4))* y(4)
         prod(3) = (rxt(5) +.500_r8*rxt(6) +rxt(7))*y(5)
         loss(4) = ( + rxt(5) + rxt(6) + rxt(7) + het_rates(5))* y(5)
         prod(4) = 0._r8
         loss(5) = ( + het_rates(6))* y(6)
         prod(5) = 0._r8
         loss(6) = ( + het_rates(7))* y(7)
         prod(6) = 0._r8
         loss(7) = ( + het_rates(8))* y(8)
         prod(7) = 0._r8
         loss(8) = ( + het_rates(9))* y(9)
         prod(8) = 0._r8
         loss(9) = ( + het_rates(10))* y(10)
         prod(9) = 0._r8
         loss(10) = ( + het_rates(12))* y(12)
         prod(10) = 0._r8
         loss(11) = ( + het_rates(13))* y(13)
         prod(11) = 0._r8
         loss(12) = ( + het_rates(11))* y(11)
         prod(12) = 0._r8
         loss(13) = ( + het_rates(14))* y(14)
         prod(13) = 0._r8
         loss(14) = ( + het_rates(15))* y(15)
         prod(14) = 0._r8
         loss(15) = ( + het_rates(16))* y(16)
         prod(15) = 0._r8
         loss(16) = ( + het_rates(17))* y(17)
         prod(16) = 0._r8
         loss(17) = ( + het_rates(18))* y(18)
         prod(17) = 0._r8
         loss(18) = ( + het_rates(20))* y(20)
         prod(18) = 0._r8
         loss(19) = ( + het_rates(21))* y(21)
         prod(19) = 0._r8
         loss(20) = ( + het_rates(19))* y(19)
         prod(20) = 0._r8
         loss(21) = ( + het_rates(22))* y(22)
         prod(21) = 0._r8
         loss(22) = ( + het_rates(23))* y(23)
         prod(22) = 0._r8
         loss(23) = ( + het_rates(24))* y(24)
         prod(23) = 0._r8
         loss(24) = ( + het_rates(25))* y(25)
         prod(24) = 0._r8
         loss(25) = ( + het_rates(26))* y(26)
         prod(25) = 0._r8
         loss(26) = ( + het_rates(27))* y(27)
         prod(26) = 0._r8
         loss(27) = ( + het_rates(29))* y(29)
         prod(27) = 0._r8
         loss(28) = ( + het_rates(30))* y(30)
         prod(28) = 0._r8
         loss(29) = ( + het_rates(28))* y(28)
         prod(29) = 0._r8
         loss(30) = ( + het_rates(31))* y(31)
         prod(30) = 0._r8
         loss(31) = ( + het_rates(32))* y(32)
         prod(31) = 0._r8
         loss(32) = ( + het_rates(33))* y(33)
         prod(32) = 0._r8
         loss(33) = ( + het_rates(34))* y(34)
         prod(33) = 0._r8

      end subroutine imp_prod_loss

      end module mo_prod_loss

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
!--------------------------------------------------------------------
!       ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,:,1) =rxt(:,:,2)
                                                                                          
         prod(:,:,2) = 0._r8
                                                                                          
         prod(:,:,3) = + extfrc(:,:,1)
                                                                                          
         prod(:,:,4) = 0._r8
                                                                                          
         prod(:,:,5) = + extfrc(:,:,9)
                                                                                          
         prod(:,:,6) = + extfrc(:,:,2)
                                                                                          
         prod(:,:,7) = 0._r8
                                                                                          
         prod(:,:,8) = 0._r8
                                                                                          
         prod(:,:,9) = 0._r8
                                                                                          
         prod(:,:,10) = 0._r8
                                                                                          
         prod(:,:,11) = 0._r8
                                                                                          
         prod(:,:,12) = 0._r8
                                                                                          
         prod(:,:,13) = 0._r8
                                                                                          
         prod(:,:,14) = + extfrc(:,:,6)
                                                                                          
         prod(:,:,15) = + extfrc(:,:,3)
                                                                                          
         prod(:,:,16) = 0._r8
                                                                                          
         prod(:,:,17) = 0._r8
                                                                                          
         prod(:,:,18) = 0._r8
                                                                                          
         prod(:,:,19) = + extfrc(:,:,7)
                                                                                          
         prod(:,:,20) = 0._r8
                                                                                          
         prod(:,:,21) = 0._r8
                                                                                          
         prod(:,:,22) = 0._r8
                                                                                          
         prod(:,:,23) = 0._r8
                                                                                          
         prod(:,:,24) = 0._r8
                                                                                          
         prod(:,:,25) = 0._r8
                                                                                          
         prod(:,:,26) = 0._r8
                                                                                          
         prod(:,:,27) = 0._r8
                                                                                          
         prod(:,:,28) = 0._r8
                                                                                          
         prod(:,:,29) = 0._r8
                                                                                          
         prod(:,:,30) = + extfrc(:,:,4)
                                                                                          
         prod(:,:,31) = + extfrc(:,:,5)
                                                                                          
         prod(:,:,32) = 0._r8
                                                                                          
         prod(:,:,33) = + extfrc(:,:,8)
                                                                                          
      end if                                                                              
                                                                                          
      end subroutine indprd                                                               
                                                                                          
      end module mo_indprd                                                                

      module mo_lin_matrix

      private
      public :: linmat

      contains

      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(in)    ::  het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) ::  mat(nzcnt)

         mat(1) = -( rxt(1) + rxt(3) + het_rates(2) )

         mat(2) = -( het_rates(3) )
         mat(3) = rxt(4)

         mat(4) = -( rxt(4) + het_rates(4) )
         mat(5) = rxt(5) + .500_r8*rxt(6) + rxt(7)

         mat(6) = -( rxt(5) + rxt(6) + rxt(7) + het_rates(5) )

         mat(7) = -( het_rates(6) )

         mat(8) = -( het_rates(7) )

         mat(9) = -( het_rates(8) )

         mat(10) = -( het_rates(9) )

         mat(11) = -( het_rates(10) )

         mat(12) = -( het_rates(12) )

         mat(13) = -( het_rates(13) )

         mat(14) = -( het_rates(11) )

         mat(15) = -( het_rates(14) )

         mat(16) = -( het_rates(15) )

         mat(17) = -( het_rates(16) )

         mat(18) = -( het_rates(17) )

         mat(19) = -( het_rates(18) )

         mat(20) = -( het_rates(20) )

         mat(21) = -( het_rates(21) )

         mat(22) = -( het_rates(19) )

         mat(23) = -( het_rates(22) )

         mat(24) = -( het_rates(23) )

         mat(25) = -( het_rates(24) )

         mat(26) = -( het_rates(25) )

         mat(27) = -( het_rates(26) )

         mat(28) = -( het_rates(27) )

         mat(29) = -( het_rates(29) )

         mat(30) = -( het_rates(30) )

         mat(31) = -( het_rates(28) )

         mat(32) = -( het_rates(31) )

         mat(33) = -( het_rates(32) )

         mat(34) = -( het_rates(33) )

         mat(35) = -( het_rates(34) )


      end subroutine linmat01

      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(in)    ::  het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) ::  mat(nzcnt)

      call linmat01( mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix

      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: nlnmat

      contains

      subroutine     nlnmat( mat, y, rxt, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  dti
      real(r8), intent(in)    ::  lmat(nzcnt)
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(inout) ::  mat(nzcnt)


































      call     nlnmat_finit( mat, lmat, dti )

      end subroutine nlnmat

      subroutine     nlnmat_finit( mat, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  dti
      real(r8), intent(in)    ::  lmat(nzcnt)
      real(r8), intent(inout) ::  mat(nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------


         mat(   1) = lmat(   1)
         mat(   2) = lmat(   2)
         mat(   3) = lmat(   3)
         mat(   4) = lmat(   4)
         mat(   5) = lmat(   5)
         mat(   6) = lmat(   6)
         mat(   7) = lmat(   7)
         mat(   8) = lmat(   8)
         mat(   9) = lmat(   9)
         mat(  10) = lmat(  10)
         mat(  11) = lmat(  11)
         mat(  12) = lmat(  12)
         mat(  13) = lmat(  13)
         mat(  14) = lmat(  14)
         mat(  15) = lmat(  15)
         mat(  16) = lmat(  16)
         mat(  17) = lmat(  17)
         mat(  18) = lmat(  18)
         mat(  19) = lmat(  19)
         mat(  20) = lmat(  20)
         mat(  21) = lmat(  21)
         mat(  22) = lmat(  22)
         mat(  23) = lmat(  23)
         mat(  24) = lmat(  24)
         mat(  25) = lmat(  25)
         mat(  26) = lmat(  26)
         mat(  27) = lmat(  27)
         mat(  28) = lmat(  28)
         mat(  29) = lmat(  29)
         mat(  30) = lmat(  30)
         mat(  31) = lmat(  31)
         mat(  32) = lmat(  32)
         mat(  33) = lmat(  33)
         mat(  34) = lmat(  34)
         mat(  35) = lmat(  35)
         mat(   1) = mat(   1) - dti
         mat(   2) = mat(   2) - dti
         mat(   4) = mat(   4) - dti
         mat(   6) = mat(   6) - dti
         mat(   7) = mat(   7) - dti
         mat(   8) = mat(   8) - dti
         mat(   9) = mat(   9) - dti
         mat(  10) = mat(  10) - dti
         mat(  11) = mat(  11) - dti
         mat(  12) = mat(  12) - dti
         mat(  13) = mat(  13) - dti
         mat(  14) = mat(  14) - dti
         mat(  15) = mat(  15) - dti
         mat(  16) = mat(  16) - dti
         mat(  17) = mat(  17) - dti
         mat(  18) = mat(  18) - dti
         mat(  19) = mat(  19) - dti
         mat(  20) = mat(  20) - dti
         mat(  21) = mat(  21) - dti
         mat(  22) = mat(  22) - dti
         mat(  23) = mat(  23) - dti
         mat(  24) = mat(  24) - dti
         mat(  25) = mat(  25) - dti
         mat(  26) = mat(  26) - dti
         mat(  27) = mat(  27) - dti
         mat(  28) = mat(  28) - dti
         mat(  29) = mat(  29) - dti
         mat(  30) = mat(  30) - dti
         mat(  31) = mat(  31) - dti
         mat(  32) = mat(  32) - dti
         mat(  33) = mat(  33) - dti
         mat(  34) = mat(  34) - dti
         mat(  35) = mat(  35) - dti

      end subroutine nlnmat_finit

      end module mo_nln_matrix

      module mo_lu_factor

      private
      public :: lu_fac

      contains
                                                                        
      subroutine lu_fac01( lu )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) ::   lu(:)
                                                                        
         lu(1) = 1._r8 / lu(1)
                                                                        
         lu(2) = 1._r8 / lu(2)
                                                                        
         lu(4) = 1._r8 / lu(4)
                                                                        
         lu(6) = 1._r8 / lu(6)
                                                                        
         lu(7) = 1._r8 / lu(7)
                                                                        
         lu(8) = 1._r8 / lu(8)
                                                                        
         lu(9) = 1._r8 / lu(9)
                                                                        
         lu(10) = 1._r8 / lu(10)
                                                                        
         lu(11) = 1._r8 / lu(11)
                                                                        
         lu(12) = 1._r8 / lu(12)
                                                                        
         lu(13) = 1._r8 / lu(13)
                                                                        
         lu(14) = 1._r8 / lu(14)
                                                                        
         lu(15) = 1._r8 / lu(15)
                                                                        
         lu(16) = 1._r8 / lu(16)
                                                                        
         lu(17) = 1._r8 / lu(17)
                                                                        
         lu(18) = 1._r8 / lu(18)
                                                                        
         lu(19) = 1._r8 / lu(19)
                                                                        
         lu(20) = 1._r8 / lu(20)
                                                                        
         lu(21) = 1._r8 / lu(21)
                                                                        
         lu(22) = 1._r8 / lu(22)
                                                                        
         lu(23) = 1._r8 / lu(23)
                                                                        
         lu(24) = 1._r8 / lu(24)
                                                                        
         lu(25) = 1._r8 / lu(25)
                                                                        
         lu(26) = 1._r8 / lu(26)
                                                                        
         lu(27) = 1._r8 / lu(27)
                                                                        
         lu(28) = 1._r8 / lu(28)
                                                                        
         lu(29) = 1._r8 / lu(29)
                                                                        
         lu(30) = 1._r8 / lu(30)
                                                                        
         lu(31) = 1._r8 / lu(31)
                                                                        
         lu(32) = 1._r8 / lu(32)
                                                                        
         lu(33) = 1._r8 / lu(33)
                                                                        
         lu(34) = 1._r8 / lu(34)
                                                                        
         lu(35) = 1._r8 / lu(35)
                                                                        
                                                                        
      end subroutine lu_fac01
                                                                        
      subroutine lu_fac( lu )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) ::   lu(:)
                                                                        
      call lu_fac01( lu )
                                                                        
      end subroutine lu_fac
                                                                        
      end module mo_lu_factor

      module mo_lu_solve

      private
      public :: lu_slv

      contains
                                                                        
      subroutine lu_slv01( lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real(r8), intent(in) ::   lu(:)
      real(r8), intent(inout) ::   b(:)
                                                                        
!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
                                                                        
!-----------------------------------------------------------------------
!       ... solve L * y = b
!-----------------------------------------------------------------------
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!-----------------------------------------------------------------------
!       ... Solve U * x = y
!-----------------------------------------------------------------------
         b(33) = b(33) * lu(35)
                                                                        
         b(32) = b(32) * lu(34)
                                                                        
         b(31) = b(31) * lu(33)
                                                                        
         b(30) = b(30) * lu(32)
                                                                        
         b(29) = b(29) * lu(31)
                                                                        
         b(28) = b(28) * lu(30)
                                                                        
         b(27) = b(27) * lu(29)
                                                                        
         b(26) = b(26) * lu(28)
                                                                        
         b(25) = b(25) * lu(27)
                                                                        
         b(24) = b(24) * lu(26)
                                                                        
         b(23) = b(23) * lu(25)
                                                                        
         b(22) = b(22) * lu(24)
                                                                        
         b(21) = b(21) * lu(23)
                                                                        
         b(20) = b(20) * lu(22)
                                                                        
         b(19) = b(19) * lu(21)
                                                                        
         b(18) = b(18) * lu(20)
                                                                        
         b(17) = b(17) * lu(19)
                                                                        
         b(16) = b(16) * lu(18)
                                                                        
         b(15) = b(15) * lu(17)
                                                                        
         b(14) = b(14) * lu(16)
                                                                        
         b(13) = b(13) * lu(15)
                                                                        
         b(12) = b(12) * lu(14)
                                                                        
         b(11) = b(11) * lu(13)
                                                                        
         b(10) = b(10) * lu(12)
                                                                        
         b(9) = b(9) * lu(11)
                                                                        
         b(8) = b(8) * lu(10)
                                                                        
         b(7) = b(7) * lu(9)
                                                                        
         b(6) = b(6) * lu(8)
                                                                        
         b(5) = b(5) * lu(7)
                                                                        
         b(4) = b(4) * lu(6)
         b(3) = b(3) - lu(5) * b(4)
                                                                        
         b(3) = b(3) * lu(4)
         b(2) = b(2) - lu(3) * b(3)
                                                                        
         b(2) = b(2) * lu(2)
                                                                        
         b(1) = b(1) * lu(1)
                                                                        
                                                                        
      end subroutine lu_slv01
                                                                        
      subroutine lu_slv( lu, b )
                                                                        

      use shr_kind_mod, only : r8 => shr_kind_r8
                                                                        
      implicit none
                                                                        
!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      real(r8), intent(in) ::   lu(:)
      real(r8), intent(inout) ::   b(:)
                                                                        
      call lu_slv01( lu, b )
                                                                        
      end subroutine lu_slv
                                                                        
      end module mo_lu_solve
