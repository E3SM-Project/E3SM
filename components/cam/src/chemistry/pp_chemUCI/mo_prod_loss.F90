












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


         loss(:,:,1) = ((rxt(:,:,20) +rxt(:,:,21))* y(:,:,2) + het_rates(:,:,15)) &
                 * y(:,:,15)
         prod(:,:,1) =.330_r8*rxt(:,:,15)*y(:,:,21)
         loss(:,:,2) = (rxt(:,:,19)* y(:,:,1) +rxt(:,:,23)* y(:,:,2) &
                  + het_rates(:,:,16))* y(:,:,16)
         prod(:,:,2) = 0._r8
         loss(:,:,3) = (rxt(:,:,24)* y(:,:,2) + het_rates(:,:,17))* y(:,:,17)
         prod(:,:,3) = 0._r8
         loss(:,:,4) = (rxt(:,:,25)* y(:,:,2) + het_rates(:,:,18))* y(:,:,18)
         prod(:,:,4) = 0._r8
         loss(:,:,5) = ((rxt(:,:,68) +rxt(:,:,69))* y(:,:,2) + rxt(:,:,15) &
                  + het_rates(:,:,21))* y(:,:,21)
         prod(:,:,5) =.800_r8*rxt(:,:,25)*y(:,:,18)*y(:,:,2)
         loss(:,:,6) = ( + rxt(:,:,83) + het_rates(:,:,30))* y(:,:,30)
         prod(:,:,6) = 0._r8
         loss(:,:,7) = ( + het_rates(:,:,31))* y(:,:,31)
         prod(:,:,7) = 0._r8

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


         loss(19) = (rxt(31)* y(2) +rxt(33)* y(3) +rxt(39)* y(8) +rxt(41)* y(9) &
                  +rxt(19)* y(16) +rxt(27)* y(19) +rxt(28)* y(26) +rxt(76)* y(28) &
                  + rxt(17) + rxt(18) + het_rates(1))* y(1)
         prod(19) =rxt(32)*y(2)*y(2) +rxt(8)*y(9) +.886_r8*rxt(9)*y(10)
         loss(22) = (rxt(31)* y(1) + 2._r8*(rxt(32) +rxt(37))* y(2) +rxt(34)* y(3) &
                  +rxt(38)* y(4) +rxt(30)* y(5) +rxt(56)* y(7) +rxt(46)* y(9) +rxt(42) &
                 * y(10) +rxt(45)* y(12) +rxt(44)* y(13) + (rxt(20) +rxt(21))* y(15) &
                  +rxt(23)* y(16) +rxt(24)* y(17) +rxt(25)* y(18) +rxt(26)* y(19) &
                  + (rxt(68) +rxt(69))* y(21) +rxt(61)* y(23) +rxt(63)* y(24) +rxt(29) &
                 * y(26) +rxt(77)* y(28) + rxt(22) + het_rates(2))* y(2)
         prod(22) = (2.000_r8*rxt(17) +rxt(18) +rxt(19)*y(16) +rxt(28)*y(26) + &
                 rxt(33)*y(3))*y(1) + (rxt(40)*y(8) +rxt(71)*y(20))*y(3) + (rxt(6) + &
                 .300_r8*rxt(56)*y(2))*y(7) +2.000_r8*rxt(3)*y(4) +rxt(11)*y(12) &
                  +.330_r8*rxt(12)*y(13) +rxt(7)*y(23)
         loss(24) = (rxt(33)* y(1) +rxt(34)* y(2) + 2._r8*(rxt(35) +rxt(36))* y(3) &
                  +rxt(53)* y(6) +rxt(40)* y(8) +rxt(47)* y(9) +rxt(71)* y(20) &
                  +rxt(60)* y(22) +rxt(74)* y(27) +rxt(79)* y(29) + rxt(82) &
                  + het_rates(3))* y(3)
         prod(24) = (rxt(22) +rxt(21)*y(15) +rxt(25)*y(18) +rxt(30)*y(5) + &
                 rxt(31)*y(1) +rxt(38)*y(4) +rxt(42)*y(10))*y(2) + (rxt(54)*y(8) + &
                 2.000_r8*rxt(55)*y(6) +2.000_r8*rxt(59)*y(22) +rxt(66)*y(25) + &
                 rxt(72)*y(20) +2.000_r8*rxt(75)*y(27))*y(6) + (.670_r8*rxt(12) + &
                 rxt(50))*y(13) + (rxt(57)*y(8) +2.000_r8*rxt(58)*y(22))*y(22) &
                  +rxt(18)*y(1) +2.000_r8*rxt(4)*y(5) +rxt(6)*y(7) +rxt(7)*y(23) &
                  +rxt(13)*y(24) +.500_r8*rxt(16)*y(28)
         loss(1) = (rxt(38)* y(2) + rxt(3) + het_rates(4))* y(4)
         prod(1) = (rxt(35)*y(3) +rxt(36)*y(3))*y(3) +rxt(37)*y(2)*y(2)
         loss(10) = (rxt(30)* y(2) + rxt(4) + rxt(5) + het_rates(5))* y(5)
         prod(10) = (rxt(54)*y(8) +2.000_r8*rxt(55)*y(6) +rxt(59)*y(22) + &
                 rxt(66)*y(25) +2.000_r8*rxt(72)*y(20) +2.000_r8*rxt(75)*y(27))*y(6) &
                  + (rxt(27)*y(19) +rxt(28)*y(26) +rxt(76)*y(28))*y(1) &
                  + (rxt(68)*y(21) +rxt(69)*y(21) +.300_r8*rxt(56)*y(7))*y(2) &
                  + (rxt(70)*y(20) +rxt(73)*y(27) +rxt(78)*y(29))*y(8) +rxt(71)*y(20) &
                 *y(3) +rxt(6)*y(7) +2.000_r8*rxt(16)*y(28)
         loss(20) = (rxt(53)* y(3) + 2._r8*rxt(55)* y(6) +rxt(54)* y(8) +rxt(72) &
                 * y(20) +rxt(59)* y(22) +rxt(75)* y(27) + het_rates(6))* y(6)
         prod(20) = (rxt(23)*y(16) +.700_r8*rxt(56)*y(7) +.500_r8*rxt(63)*y(24))*y(2) &
                  + (rxt(65)*y(8) +2.000_r8*rxt(67)*y(25))*y(25) + (.500_r8*rxt(16) + &
                 rxt(76)*y(1))*y(28) +rxt(13)*y(24)
         loss(6) = (rxt(56)* y(2) + rxt(6) + het_rates(7))* y(7)
         prod(6) =rxt(53)*y(6)*y(3)
         loss(23) = (rxt(39)* y(1) +rxt(40)* y(3) +rxt(54)* y(6) +rxt(43)* y(10) &
                  +rxt(70)* y(20) +rxt(57)* y(22) +rxt(65)* y(25) +rxt(73)* y(27) &
                  +rxt(78)* y(29) + het_rates(8))* y(8)
         prod(23) =rxt(8)*y(9) +.114_r8*rxt(9)*y(10)
         loss(18) = (rxt(41)* y(1) +rxt(46)* y(2) +rxt(47)* y(3) +rxt(48)* y(10) &
                  +rxt(49)* y(25) + rxt(8) + het_rates(9))* y(9)
         prod(18) = (rxt(39)*y(1) +rxt(40)*y(3) +2.000_r8*rxt(43)*y(10) + &
                 rxt(54)*y(6) +rxt(57)*y(22) +rxt(65)*y(25) +rxt(70)*y(20) + &
                 rxt(73)*y(27) +rxt(78)*y(29))*y(8) + (.670_r8*rxt(12) +rxt(50) + &
                 rxt(44)*y(2))*y(13) + (.886_r8*rxt(9) +rxt(42)*y(2))*y(10) &
                  + (rxt(10) +rxt(51))*y(11) + (rxt(14) +rxt(52))*y(14) +rxt(11)*y(12)
         loss(21) = (rxt(42)* y(2) +rxt(43)* y(8) +rxt(48)* y(9) +rxt(64)* y(24) &
                  + rxt(9) + rxt(81) + het_rates(10))* y(10)
         prod(21) = (rxt(10) +rxt(51))*y(11) +rxt(41)*y(9)*y(1) +rxt(45)*y(12)*y(2) &
                  +.330_r8*rxt(12)*y(13)
         loss(3) = ( + rxt(10) + rxt(51) + rxt(80) + het_rates(11))* y(11)
         prod(3) =rxt(48)*y(10)*y(9)
         loss(7) = (rxt(45)* y(2) + rxt(11) + het_rates(12))* y(12)
         prod(7) = (rxt(81) +rxt(64)*y(24))*y(10) +rxt(46)*y(9)*y(2) +2.000_r8*rxt(80) &
                 *y(11)
         loss(8) = (rxt(44)* y(2) + rxt(12) + rxt(50) + het_rates(13))* y(13)
         prod(8) =rxt(47)*y(9)*y(3)
         loss(2) = ( + rxt(14) + rxt(52) + het_rates(14))* y(14)
         prod(2) =rxt(49)*y(25)*y(9)
         loss(13) = (rxt(60)* y(3) +rxt(59)* y(6) +rxt(57)* y(8) + 2._r8*rxt(58) &
                 * y(22) + het_rates(22))* y(22)
         prod(13) = (rxt(24)*y(17) +rxt(61)*y(23))*y(2)
         loss(9) = ((rxt(61) +rxt(62))* y(2) + rxt(7) + het_rates(23))* y(23)
         prod(9) = (rxt(60)*y(22) +rxt(74)*y(27) +rxt(79)*y(29))*y(3)
         loss(16) = (rxt(63)* y(2) +rxt(64)* y(10) + rxt(13) + het_rates(24))* y(24)
         prod(16) = (rxt(57)*y(8) +2.000_r8*rxt(58)*y(22) +rxt(59)*y(6))*y(22) &
                  + (.200_r8*rxt(25)*y(18) +rxt(62)*y(23))*y(2) + (rxt(70)*y(8) + &
                 rxt(71)*y(3))*y(20) +.500_r8*rxt(27)*y(19)*y(1) +rxt(7)*y(23)
         loss(17) = (rxt(66)* y(6) +rxt(65)* y(8) +rxt(49)* y(9) + 2._r8*rxt(67) &
                 * y(25) + het_rates(25))* y(25)
         prod(17) = (rxt(68)*y(21) +rxt(69)*y(21) +.500_r8*rxt(63)*y(24))*y(2) &
                  + (rxt(14) +rxt(52))*y(14) +.500_r8*rxt(76)*y(28)*y(1) &
                  +.500_r8*rxt(78)*y(29)*y(8) +rxt(64)*y(24)*y(10)
         loss(4) = (rxt(27)* y(1) +rxt(26)* y(2) + het_rates(19))* y(19)
         prod(4) = 0._r8
         loss(12) = (rxt(71)* y(3) +rxt(72)* y(6) +rxt(70)* y(8) + het_rates(20)) &
                 * y(20)
         prod(12) =rxt(26)*y(19)*y(2)
         loss(5) = (rxt(28)* y(1) +rxt(29)* y(2) + het_rates(26))* y(26)
         prod(5) = 0._r8
         loss(14) = (rxt(74)* y(3) +rxt(75)* y(6) +rxt(73)* y(8) + het_rates(27)) &
                 * y(27)
         prod(14) =rxt(29)*y(26)*y(2)
         loss(15) = (rxt(76)* y(1) +rxt(77)* y(2) + rxt(16) + het_rates(28))* y(28)
         prod(15) = (rxt(73)*y(8) +rxt(75)*y(6))*y(27) +rxt(28)*y(26)*y(1)
         loss(11) = (rxt(79)* y(3) +rxt(78)* y(8) + het_rates(29))* y(29)
         prod(11) =rxt(77)*y(28)*y(2)

      end subroutine imp_prod_loss

      end module mo_prod_loss
