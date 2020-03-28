












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


         loss(:,:,1) = ((rxt(:,:,19) +rxt(:,:,20))* y(:,:,2) + het_rates(:,:,15)) &
                 * y(:,:,15)
         prod(:,:,1) =.330_r8*rxt(:,:,17)*y(:,:,21)
         loss(:,:,2) = (rxt(:,:,22)* y(:,:,2) + het_rates(:,:,16))* y(:,:,16)
         prod(:,:,2) = 0._r8
         loss(:,:,3) = (rxt(:,:,23)* y(:,:,2) + het_rates(:,:,17))* y(:,:,17)
         prod(:,:,3) = 0._r8
         loss(:,:,4) = (rxt(:,:,24)* y(:,:,2) + het_rates(:,:,18))* y(:,:,18)
         prod(:,:,4) = 0._r8
         loss(:,:,5) = ((rxt(:,:,67) +rxt(:,:,68))* y(:,:,2) + rxt(:,:,17) &
                  + het_rates(:,:,21))* y(:,:,21)
         prod(:,:,5) =.800_r8*rxt(:,:,24)*y(:,:,18)*y(:,:,2)
         loss(:,:,6) = ( + rxt(:,:,86) + het_rates(:,:,30))* y(:,:,30)
         prod(:,:,6) = 0._r8
         loss(:,:,7) = ( + het_rates(:,:,31))* y(:,:,31)
         prod(:,:,7) = 0._r8
         loss(:,:,8) = (rxt(:,:,81)* y(:,:,3) +rxt(:,:,80)* y(:,:,10) +rxt(:,:,79) &
                 * y(:,:,11) + het_rates(:,:,32))* y(:,:,32)
         prod(:,:,8) = 0._r8

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


         loss(53) = (rxt(30)* y(2) +rxt(32)* y(3) +rxt(38)* y(8) +rxt(40)* y(9) &
                  +rxt(26)* y(19) +rxt(27)* y(26) +rxt(75)* y(28) + rxt(1) + rxt(2) &
                  + rxt(3) + het_rates(1))* y(1)
         prod(53) =rxt(31)*y(2)*y(2) +rxt(10)*y(9) +.886_r8*rxt(11)*y(10)
         loss(52) = (rxt(30)* y(1) + 2._r8*(rxt(31) +rxt(36))* y(2) +rxt(33)* y(3) &
                  +rxt(37)* y(4) +rxt(29)* y(5) +rxt(55)* y(7) +rxt(45)* y(9) +rxt(41) &
                 * y(10) +rxt(44)* y(12) + (rxt(43) +rxt(82))* y(13) + (rxt(19) + &
                 rxt(20))* y(15) +rxt(22)* y(16) +rxt(23)* y(17) +rxt(24)* y(18) &
                  +rxt(25)* y(19) + (rxt(67) +rxt(68))* y(21) +rxt(60)* y(23) +rxt(62) &
                 * y(24) +rxt(28)* y(26) +rxt(76)* y(28) + (rxt(83) +rxt(84))* y(33) &
                  +rxt(85)* y(34) + rxt(21) + het_rates(2))* y(2)
         prod(52) = (2.000_r8*rxt(1) +rxt(2) +rxt(3) +rxt(27)*y(26) +rxt(32)*y(3)) &
                 *y(1) + (rxt(39)*y(8) +rxt(70)*y(20))*y(3) + (rxt(8) + &
                 .300_r8*rxt(55)*y(2))*y(7) +2.000_r8*rxt(5)*y(4) +rxt(13)*y(12) &
                  +.330_r8*rxt(14)*y(13) +rxt(9)*y(23)
         loss(51) = (rxt(32)* y(1) +rxt(33)* y(2) + 2._r8*(rxt(34) +rxt(35))* y(3) &
                  +rxt(52)* y(6) +rxt(39)* y(8) +rxt(46)* y(9) +rxt(70)* y(20) &
                  +rxt(59)* y(22) +rxt(73)* y(27) +rxt(78)* y(29) +rxt(81)* y(32) &
                  + het_rates(3))* y(3)
         prod(51) = (rxt(21) +rxt(20)*y(15) +rxt(24)*y(18) +rxt(29)*y(5) + &
                 rxt(30)*y(1) +rxt(37)*y(4) +rxt(41)*y(10) +.500_r8*rxt(84)*y(33)) &
                 *y(2) + (rxt(53)*y(8) +2.000_r8*rxt(54)*y(6) + &
                 2.000_r8*rxt(58)*y(22) +rxt(65)*y(25) +rxt(71)*y(20) + &
                 2.000_r8*rxt(74)*y(27))*y(6) + (.670_r8*rxt(14) +rxt(49))*y(13) &
                  + (rxt(56)*y(8) +2.000_r8*rxt(57)*y(22))*y(22) +rxt(2)*y(1) &
                  +2.000_r8*rxt(6)*y(5) +rxt(8)*y(7) +rxt(9)*y(23) +rxt(15)*y(24) &
                  +.500_r8*rxt(18)*y(28)
         loss(30) = (rxt(37)* y(2) + rxt(5) + het_rates(4))* y(4)
         prod(30) = (rxt(34)*y(3) +rxt(35)*y(3))*y(3) +rxt(36)*y(2)*y(2)
         loss(39) = (rxt(29)* y(2) + rxt(6) + rxt(7) + het_rates(5))* y(5)
         prod(39) = (rxt(53)*y(8) +2.000_r8*rxt(54)*y(6) +rxt(58)*y(22) + &
                 rxt(65)*y(25) +2.000_r8*rxt(71)*y(20) +2.000_r8*rxt(74)*y(27))*y(6) &
                  + (rxt(26)*y(19) +rxt(27)*y(26) +rxt(75)*y(28))*y(1) &
                  + (rxt(67)*y(21) +rxt(68)*y(21) +.300_r8*rxt(55)*y(7))*y(2) &
                  + (rxt(69)*y(20) +rxt(72)*y(27) +rxt(77)*y(29))*y(8) +rxt(70)*y(20) &
                 *y(3) +rxt(8)*y(7) +2.000_r8*rxt(18)*y(28)
         loss(48) = (rxt(52)* y(3) + 2._r8*rxt(54)* y(6) +rxt(53)* y(8) +rxt(71) &
                 * y(20) +rxt(58)* y(22) +rxt(74)* y(27) + het_rates(6))* y(6)
         prod(48) = (rxt(22)*y(16) +.700_r8*rxt(55)*y(7) +.500_r8*rxt(62)*y(24))*y(2) &
                  + (rxt(64)*y(8) +2.000_r8*rxt(66)*y(25))*y(25) + (.500_r8*rxt(18) + &
                 rxt(75)*y(1))*y(28) +rxt(15)*y(24)
         loss(37) = (rxt(55)* y(2) + rxt(8) + het_rates(7))* y(7)
         prod(37) =rxt(52)*y(6)*y(3)
         loss(47) = (rxt(38)* y(1) +rxt(39)* y(3) +rxt(53)* y(6) +rxt(42)* y(10) &
                  +rxt(69)* y(20) +rxt(56)* y(22) +rxt(64)* y(25) +rxt(72)* y(27) &
                  +rxt(77)* y(29) + het_rates(8))* y(8)
         prod(47) =rxt(10)*y(9) +.114_r8*rxt(11)*y(10)
         loss(49) = (rxt(40)* y(1) +rxt(45)* y(2) +rxt(46)* y(3) +rxt(47)* y(10) &
                  +rxt(48)* y(25) + rxt(10) + het_rates(9))* y(9)
         prod(49) = (rxt(38)*y(1) +rxt(39)*y(3) +2.000_r8*rxt(42)*y(10) + &
                 rxt(53)*y(6) +rxt(56)*y(22) +rxt(64)*y(25) +rxt(69)*y(20) + &
                 rxt(72)*y(27) +rxt(77)*y(29))*y(8) + (.670_r8*rxt(14) +rxt(49) + &
                 rxt(43)*y(2) +rxt(82)*y(2))*y(13) + (.886_r8*rxt(11) +rxt(41)*y(2)) &
                 *y(10) + (rxt(12) +rxt(50))*y(11) + (rxt(16) +rxt(51))*y(14) +rxt(13) &
                 *y(12)
         loss(50) = (rxt(41)* y(2) +rxt(42)* y(8) +rxt(47)* y(9) +rxt(63)* y(24) &
                  +rxt(80)* y(32) + rxt(11) + het_rates(10))* y(10)
         prod(50) = (rxt(12) +rxt(50))*y(11) +rxt(40)*y(9)*y(1) +rxt(44)*y(12)*y(2) &
                  +.330_r8*rxt(14)*y(13)
         loss(32) = (rxt(79)* y(32) + rxt(12) + rxt(50) + het_rates(11))* y(11)
         prod(32) =rxt(47)*y(10)*y(9)
         loss(35) = (rxt(44)* y(2) + rxt(13) + het_rates(12))* y(12)
         prod(35) = (rxt(80)*y(32) +rxt(63)*y(24))*y(10) +rxt(45)*y(9)*y(2) &
                  +2.000_r8*rxt(79)*y(32)*y(11)
         loss(36) = ((rxt(43) +rxt(82))* y(2) + rxt(14) + rxt(49) + het_rates(13)) &
                 * y(13)
         prod(36) =rxt(46)*y(9)*y(3)
         loss(31) = ( + rxt(16) + rxt(51) + het_rates(14))* y(14)
         prod(31) =rxt(48)*y(25)*y(9)
         loss(42) = (rxt(59)* y(3) +rxt(58)* y(6) +rxt(56)* y(8) + 2._r8*rxt(57) &
                 * y(22) + het_rates(22))* y(22)
         prod(42) = (rxt(23)*y(17) +rxt(60)*y(23))*y(2)
         loss(38) = ((rxt(60) +rxt(61))* y(2) + rxt(9) + het_rates(23))* y(23)
         prod(38) = (rxt(59)*y(22) +rxt(73)*y(27) +rxt(78)*y(29))*y(3)
         loss(45) = (rxt(62)* y(2) +rxt(63)* y(10) + rxt(15) + het_rates(24))* y(24)
         prod(45) = (rxt(56)*y(8) +2.000_r8*rxt(57)*y(22) +rxt(58)*y(6))*y(22) &
                  + (.200_r8*rxt(24)*y(18) +rxt(61)*y(23))*y(2) + (rxt(69)*y(8) + &
                 rxt(70)*y(3))*y(20) +.500_r8*rxt(26)*y(19)*y(1) +rxt(9)*y(23)
         loss(46) = (rxt(65)* y(6) +rxt(64)* y(8) +rxt(48)* y(9) + 2._r8*rxt(66) &
                 * y(25) + het_rates(25))* y(25)
         prod(46) = (rxt(67)*y(21) +rxt(68)*y(21) +.500_r8*rxt(62)*y(24))*y(2) &
                  + (rxt(16) +rxt(51))*y(14) +.500_r8*rxt(75)*y(28)*y(1) &
                  +.500_r8*rxt(77)*y(29)*y(8) +rxt(63)*y(24)*y(10)
         loss(33) = (rxt(26)* y(1) +rxt(25)* y(2) + het_rates(19))* y(19)
         prod(33) = 0._r8
         loss(41) = (rxt(70)* y(3) +rxt(71)* y(6) +rxt(69)* y(8) + het_rates(20)) &
                 * y(20)
         prod(41) =rxt(25)*y(19)*y(2)
         loss(34) = (rxt(27)* y(1) +rxt(28)* y(2) + het_rates(26))* y(26)
         prod(34) = 0._r8
         loss(43) = (rxt(73)* y(3) +rxt(74)* y(6) +rxt(72)* y(8) + het_rates(27)) &
                 * y(27)
         prod(43) =rxt(28)*y(26)*y(2)
         loss(44) = (rxt(75)* y(1) +rxt(76)* y(2) + rxt(18) + het_rates(28))* y(28)
         prod(44) = (rxt(72)*y(8) +rxt(74)*y(6))*y(27) +rxt(27)*y(26)*y(1)
         loss(40) = (rxt(78)* y(3) +rxt(77)* y(8) + het_rates(29))* y(29)
         prod(40) =rxt(76)*y(28)*y(2)
         loss(29) = ((rxt(83) +rxt(84))* y(2) + het_rates(33))* y(33)
         prod(29) = 0._r8
         loss(28) = (rxt(85)* y(2) + het_rates(34))* y(34)
         prod(28) = (rxt(83)*y(33) +.500_r8*rxt(84)*y(33))*y(2)
         loss(1) = ( + het_rates(35))* y(35)
         prod(1) =rxt(85)*y(34)*y(2)
         loss(2) = ( + het_rates(36))* y(36)
         prod(2) = 0._r8
         loss(3) = ( + het_rates(37))* y(37)
         prod(3) = 0._r8
         loss(4) = ( + het_rates(38))* y(38)
         prod(4) = 0._r8
         loss(5) = ( + het_rates(39))* y(39)
         prod(5) = 0._r8
         loss(6) = ( + het_rates(40))* y(40)
         prod(6) = 0._r8
         loss(7) = ( + het_rates(41))* y(41)
         prod(7) = 0._r8
         loss(8) = ( + het_rates(42))* y(42)
         prod(8) = 0._r8
         loss(9) = ( + het_rates(43))* y(43)
         prod(9) = 0._r8
         loss(10) = ( + het_rates(44))* y(44)
         prod(10) = 0._r8
         loss(11) = ( + het_rates(45))* y(45)
         prod(11) = 0._r8
         loss(12) = ( + het_rates(46))* y(46)
         prod(12) = 0._r8
         loss(13) = ( + het_rates(47))* y(47)
         prod(13) = 0._r8
         loss(14) = ( + het_rates(48))* y(48)
         prod(14) = 0._r8
         loss(15) = ( + het_rates(49))* y(49)
         prod(15) = 0._r8
         loss(16) = ( + het_rates(50))* y(50)
         prod(16) = 0._r8
         loss(17) = ( + het_rates(51))* y(51)
         prod(17) = 0._r8
         loss(18) = ( + het_rates(52))* y(52)
         prod(18) = 0._r8
         loss(19) = ( + het_rates(53))* y(53)
         prod(19) = 0._r8
         loss(20) = ( + het_rates(54))* y(54)
         prod(20) = 0._r8
         loss(21) = ( + het_rates(55))* y(55)
         prod(21) = 0._r8
         loss(22) = ( + het_rates(56))* y(56)
         prod(22) = 0._r8
         loss(23) = ( + het_rates(57))* y(57)
         prod(23) = 0._r8
         loss(24) = ( + het_rates(58))* y(58)
         prod(24) = 0._r8
         loss(25) = ( + het_rates(59))* y(59)
         prod(25) = 0._r8
         loss(26) = ( + het_rates(60))* y(60)
         prod(26) = 0._r8
         loss(27) = ( + het_rates(61))* y(61)
         prod(27) = 0._r8

      end subroutine imp_prod_loss

      end module mo_prod_loss
