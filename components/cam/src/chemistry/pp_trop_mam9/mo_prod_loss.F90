





      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:,:,:)
      real(r8), intent(in) :: rxt(:,:,:)
      real(r8), intent(in) :: het_rates(:,:,:)


      end subroutine exp_prod_loss

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:)
      real(r8), intent(in) :: rxt(:)
      real(r8), intent(in) :: het_rates(:)



!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------


         loss(1) = ( + rxt(1) + rxt(3) + het_rates(1))* y(1)
         prod(1) = 0._r8
         loss(2) = ( + het_rates(2))* y(2)
         prod(2) =rxt(4)*y(3)
         loss(3) = ( + rxt(4) + het_rates(3))* y(3)
         prod(3) = (rxt(5) +.500_r8*rxt(6) +rxt(7))*y(4)
         loss(4) = ( + rxt(5) + rxt(6) + rxt(7) + het_rates(4))* y(4)
         prod(4) = 0._r8
         loss(5) = ( + rxt(8) + het_rates(5))* y(5)
         prod(5) = 0._r8
         loss(6) = ( + het_rates(6))* y(6)
         prod(6) = 0._r8
         loss(7) = ( + het_rates(7))* y(7)
         prod(7) = 0._r8
         loss(8) = ( + het_rates(8))* y(8)
         prod(8) = 0._r8
         loss(9) = ( + het_rates(9))* y(9)
         prod(9) = 0._r8
         loss(10) = ( + het_rates(10))* y(10)
         prod(10) = 0._r8
         loss(11) = ( + het_rates(11))* y(11)
         prod(11) = 0._r8
         loss(12) = ( + het_rates(12))* y(12)
         prod(12) = 0._r8
         loss(13) = ( + het_rates(13))* y(13)
         prod(13) = 0._r8
         loss(14) = ( + het_rates(14))* y(14)
         prod(14) = 0._r8
         loss(15) = ( + het_rates(15))* y(15)
         prod(15) = 0._r8
         loss(16) = ( + het_rates(16))* y(16)
         prod(16) = 0._r8
         loss(17) = ( + het_rates(17))* y(17)
         prod(17) = 0._r8
         loss(18) = ( + het_rates(18))* y(18)
         prod(18) = 0._r8
         loss(19) = ( + het_rates(19))* y(19)
         prod(19) = 0._r8
         loss(20) = ( + het_rates(20))* y(20)
         prod(20) = 0._r8
         loss(21) = ( + het_rates(21))* y(21)
         prod(21) = 0._r8
         loss(22) = ( + het_rates(22))* y(22)
         prod(22) = 0._r8
         loss(23) = ( + het_rates(23))* y(23)
         prod(23) = 0._r8
         loss(24) = ( + het_rates(24))* y(24)
         prod(24) = 0._r8
         loss(25) = ( + het_rates(25))* y(25)
         prod(25) = 0._r8
         loss(26) = ( + het_rates(26))* y(26)
         prod(26) = 0._r8
         loss(27) = ( + het_rates(27))* y(27)
         prod(27) = 0._r8
         loss(28) = ( + het_rates(28))* y(28)
         prod(28) = 0._r8
         loss(29) = ( + het_rates(29))* y(29)
         prod(29) = 0._r8
         loss(30) = ( + het_rates(30))* y(30)
         prod(30) = 0._r8
         loss(31) = ( + het_rates(31))* y(31)
         prod(31) = 0._r8
         loss(32) = ( + het_rates(32))* y(32)
         prod(32) = 0._r8
         loss(33) = ( + het_rates(33))* y(33)
         prod(33) = 0._r8
         loss(34) = ( + het_rates(34))* y(34)
         prod(34) = 0._r8
         loss(35) = ( + het_rates(35))* y(35)
         prod(35) = 0._r8
         loss(36) = ( + het_rates(36))* y(36)
         prod(36) = 0._r8
         loss(37) = ( + het_rates(37))* y(37)
         prod(37) = 0._r8
         loss(38) = ( + het_rates(38))* y(38)
         prod(38) = 0._r8
         loss(39) = ( + het_rates(39))* y(39)
         prod(39) = 0._r8
         loss(40) = ( + het_rates(40))* y(40)
         prod(40) = 0._r8
         loss(41) = ( + het_rates(41))* y(41)
         prod(41) = 0._r8
         loss(42) = ( + het_rates(42))* y(42)
         prod(42) = 0._r8
         loss(43) = ( + het_rates(43))* y(43)
         prod(43) = 0._r8
         loss(44) = ( + het_rates(44))* y(44)
         prod(44) = 0._r8
         loss(45) = ( + het_rates(45))* y(45)
         prod(45) = 0._r8
         loss(46) = ( + het_rates(46))* y(46)
         prod(46) = 0._r8
         loss(47) = ( + het_rates(47))* y(47)
         prod(47) = 0._r8
         loss(48) = ( + het_rates(48))* y(48)
         prod(48) = 0._r8
         loss(49) = ( + het_rates(49))* y(49)
         prod(49) = 0._r8
         loss(50) = ( + het_rates(50))* y(50)
         prod(50) = 0._r8
         loss(51) = ( + het_rates(51))* y(51)
         prod(51) = 0._r8
         loss(52) = ( + het_rates(52))* y(52)
         prod(52) = 0._r8
         loss(53) = ( + het_rates(53))* y(53)
         prod(53) = 0._r8
         loss(54) = ( + het_rates(54))* y(54)
         prod(54) = 0._r8

      end subroutine imp_prod_loss

      end module mo_prod_loss
