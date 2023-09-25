










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


         loss(:,:,1) = ((rxt(:,:,26) +rxt(:,:,27))* y(:,:,2) + het_rates(:,:,15)) &
                 * y(:,:,15)
         prod(:,:,1) =.330_r8*rxt(:,:,18)*y(:,:,20)
         loss(:,:,2) = (rxt(:,:,30)* y(:,:,2) + het_rates(:,:,16))* y(:,:,16)
         prod(:,:,2) = 0._r8
         loss(:,:,3) = (rxt(:,:,31)* y(:,:,2) + het_rates(:,:,17))* y(:,:,17)
         prod(:,:,3) = 0._r8
         loss(:,:,4) = ((rxt(:,:,75) +rxt(:,:,76))* y(:,:,2) + rxt(:,:,18) &
                  + het_rates(:,:,20))* y(:,:,20)
         prod(:,:,4) =.800_r8*rxt(:,:,31)*y(:,:,17)*y(:,:,2)
         loss(:,:,5) = ( + rxt(:,:,87) + het_rates(:,:,31))* y(:,:,31)
         prod(:,:,5) = 0._r8
         loss(:,:,6) = ( + het_rates(:,:,32))* y(:,:,32)
         prod(:,:,6) = 0._r8
         loss(:,:,7) = ( + het_rates(:,:,33))* y(:,:,33)
         prod(:,:,7) = 0._r8
         loss(:,:,8) = (rxt(:,:,25)* y(:,:,1) +rxt(:,:,29)* y(:,:,2) &
                  + het_rates(:,:,34))* y(:,:,34)
         prod(:,:,8) = 0._r8
         loss(:,:,9) = ( + het_rates(:,:,35))* y(:,:,35)
         prod(:,:,9) = 0._r8
         loss(:,:,10) = ((rxt(:,:,88) +rxt(:,:,89))* y(:,:,2) +rxt(:,:,91)* y(:,:,10) &
                  + het_rates(:,:,36))* y(:,:,36)
         prod(:,:,10) = 0._r8
         loss(:,:,11) = (rxt(:,:,90)* y(:,:,2) + het_rates(:,:,37))* y(:,:,37)
         prod(:,:,11) = (rxt(:,:,88)*y(:,:,2) +.500_r8*rxt(:,:,89)*y(:,:,2) + &
                 rxt(:,:,91)*y(:,:,10))*y(:,:,36)
         loss(:,:,12) = ( + het_rates(:,:,38))* y(:,:,38)
         prod(:,:,12) =rxt(:,:,90)*y(:,:,37)*y(:,:,2)
         loss(:,:,13) = ( + het_rates(:,:,47))* y(:,:,47)
         prod(:,:,13) = 0._r8
         loss(:,:,14) = ( + het_rates(:,:,48))* y(:,:,48)
         prod(:,:,14) = 0._r8
         loss(:,:,15) = ( + het_rates(:,:,49))* y(:,:,49)
         prod(:,:,15) = 0._r8
         loss(:,:,16) = ( + het_rates(:,:,50))* y(:,:,50)
         prod(:,:,16) = 0._r8
         loss(:,:,17) = ( + het_rates(:,:,51))* y(:,:,51)
         prod(:,:,17) = 0._r8
         loss(:,:,18) = ( + het_rates(:,:,52))* y(:,:,52)
         prod(:,:,18) = 0._r8
         loss(:,:,19) = ( + het_rates(:,:,53))* y(:,:,53)
         prod(:,:,19) = 0._r8
         loss(:,:,20) = ( + het_rates(:,:,57))* y(:,:,57)
         prod(:,:,20) = 0._r8
         loss(:,:,21) = ( + het_rates(:,:,58))* y(:,:,58)
         prod(:,:,21) = 0._r8
         loss(:,:,22) = ( + het_rates(:,:,59))* y(:,:,59)
         prod(:,:,22) = 0._r8
         loss(:,:,23) = ( + het_rates(:,:,60))* y(:,:,60)
         prod(:,:,23) = 0._r8
         loss(:,:,24) = ( + het_rates(:,:,61))* y(:,:,61)
         prod(:,:,24) = 0._r8
         loss(:,:,25) = ( + het_rates(:,:,62))* y(:,:,62)
         prod(:,:,25) = 0._r8
         loss(:,:,26) = ( + het_rates(:,:,63))* y(:,:,63)
         prod(:,:,26) = 0._r8
         loss(:,:,27) = ( + het_rates(:,:,64))* y(:,:,64)
         prod(:,:,27) = 0._r8
         loss(:,:,28) = ( + het_rates(:,:,65))* y(:,:,65)
         prod(:,:,28) = 0._r8
         loss(:,:,29) = ( + het_rates(:,:,66))* y(:,:,66)
         prod(:,:,29) = 0._r8
         loss(:,:,30) = ( + het_rates(:,:,67))* y(:,:,67)
         prod(:,:,30) = 0._r8
         loss(:,:,31) = ( + het_rates(:,:,68))* y(:,:,68)
         prod(:,:,31) = 0._r8
         loss(:,:,32) = ( + het_rates(:,:,69))* y(:,:,69)
         prod(:,:,32) = 0._r8
         loss(:,:,33) = ( + het_rates(:,:,70))* y(:,:,70)
         prod(:,:,33) = 0._r8
         loss(:,:,34) = ( + het_rates(:,:,71))* y(:,:,71)
         prod(:,:,34) = 0._r8
         loss(:,:,35) = ( + het_rates(:,:,72))* y(:,:,72)
         prod(:,:,35) = 0._r8
         loss(:,:,36) = ( + het_rates(:,:,73))* y(:,:,73)
         prod(:,:,36) = 0._r8

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


         loss(35) = (rxt(37)* y(2) +rxt(39)* y(3) +rxt(45)* y(8) +rxt(47)* y(9) &
                  +rxt(33)* y(18) +rxt(34)* y(25) +rxt(83)* y(29) +rxt(25)* y(34) &
                  + rxt(23) + rxt(24) + het_rates(1))* y(1)
         prod(35) =rxt(38)*y(2)*y(2) +rxt(8)*y(9) +rxt(9)*y(10) +rxt(12)*y(11)
         loss(36) = (rxt(37)* y(1) + 2._r8*(rxt(38) +rxt(43))* y(2) +rxt(40)* y(3) &
                  +rxt(44)* y(4) +rxt(36)* y(5) +rxt(63)* y(7) +rxt(53)* y(9) +rxt(48) &
                 * y(10) + (rxt(51) +rxt(52))* y(12) +rxt(50)* y(13) + (rxt(26) + &
                 rxt(27))* y(15) +rxt(30)* y(16) +rxt(31)* y(17) +rxt(32)* y(18) &
                  + (rxt(75) +rxt(76))* y(20) +rxt(68)* y(22) +rxt(70)* y(23) +rxt(35) &
                 * y(25) +rxt(84)* y(29) +rxt(29)* y(34) + (rxt(88) +rxt(89))* y(36) &
                  +rxt(90)* y(37) + rxt(28) + het_rates(2))* y(2)
         prod(36) = (2.000_r8*rxt(23) +rxt(24) +rxt(25)*y(34) +rxt(34)*y(25) + &
                 rxt(39)*y(3))*y(1) + (rxt(46)*y(8) +rxt(78)*y(19))*y(3) + (rxt(6) + &
                 .300_r8*rxt(63)*y(2))*y(7) +2.000_r8*rxt(3)*y(4) +rxt(13)*y(12) &
                  +rxt(14)*y(13) +rxt(7)*y(22)
         loss(31) = (rxt(39)* y(1) +rxt(40)* y(2) + 2._r8*(rxt(41) +rxt(42))* y(3) &
                  +rxt(60)* y(6) +rxt(46)* y(8) +rxt(54)* y(9) +rxt(78)* y(19) &
                  +rxt(67)* y(21) +rxt(81)* y(27) +rxt(86)* y(30) + het_rates(3)) &
                 * y(3)
         prod(31) = (rxt(28) +rxt(27)*y(15) +rxt(31)*y(17) +.500_r8*rxt(89)*y(36) + &
                 rxt(36)*y(5) +rxt(37)*y(1) +rxt(44)*y(4) +rxt(48)*y(10))*y(2) &
                  + (rxt(61)*y(8) +2.000_r8*rxt(62)*y(6) +2.000_r8*rxt(66)*y(21) + &
                 rxt(73)*y(24) +rxt(79)*y(19) +2.000_r8*rxt(82)*y(27))*y(6) &
                  + (rxt(15) +rxt(57))*y(13) + (rxt(64)*y(8) +2.000_r8*rxt(65)*y(21)) &
                 *y(21) +rxt(24)*y(1) +2.000_r8*rxt(4)*y(5) +rxt(6)*y(7) +rxt(7)*y(22) &
                  +rxt(16)*y(23) +.500_r8*rxt(19)*y(29)
         loss(15) = (rxt(44)* y(2) + rxt(3) + het_rates(4))* y(4)
         prod(15) = (rxt(41)*y(3) +rxt(42)*y(3))*y(3) +rxt(43)*y(2)*y(2)
         loss(23) = (rxt(36)* y(2) + rxt(4) + rxt(5) + het_rates(5))* y(5)
         prod(23) = (rxt(61)*y(8) +2.000_r8*rxt(62)*y(6) +rxt(66)*y(21) + &
                 rxt(73)*y(24) +2.000_r8*rxt(79)*y(19) +2.000_r8*rxt(82)*y(27))*y(6) &
                  + (rxt(33)*y(18) +rxt(34)*y(25) +rxt(83)*y(29))*y(1) &
                  + (rxt(75)*y(20) +rxt(76)*y(20) +.300_r8*rxt(63)*y(7))*y(2) &
                  + (rxt(77)*y(19) +rxt(80)*y(27) +rxt(85)*y(30))*y(8) +rxt(78)*y(19) &
                 *y(3) +rxt(6)*y(7) +2.000_r8*rxt(19)*y(29)
         loss(37) = (rxt(60)* y(3) + 2._r8*rxt(62)* y(6) +rxt(61)* y(8) +rxt(79) &
                 * y(19) +rxt(66)* y(21) +rxt(82)* y(27) + het_rates(6))* y(6)
         prod(37) = (rxt(29)*y(34) +.700_r8*rxt(63)*y(7) +.500_r8*rxt(70)*y(23))*y(2) &
                  + (rxt(72)*y(8) +2.000_r8*rxt(74)*y(24))*y(24) + (.500_r8*rxt(19) + &
                 rxt(83)*y(1))*y(29) +rxt(16)*y(23)
         loss(19) = (rxt(63)* y(2) + rxt(6) + het_rates(7))* y(7)
         prod(19) =rxt(60)*y(6)*y(3)
         loss(32) = (rxt(45)* y(1) +rxt(46)* y(3) +rxt(61)* y(6) +rxt(49)* y(10) &
                  +rxt(77)* y(19) +rxt(64)* y(21) +rxt(72)* y(24) +rxt(80)* y(27) &
                  +rxt(85)* y(30) + het_rates(8))* y(8)
         prod(32) =rxt(8)*y(9) +rxt(10)*y(10) +rxt(12)*y(11)
         loss(33) = (rxt(47)* y(1) +rxt(53)* y(2) +rxt(54)* y(3) +rxt(55)* y(10) &
                  +rxt(56)* y(24) + rxt(8) + het_rates(9))* y(9)
         prod(33) = (rxt(45)*y(1) +rxt(46)*y(3) +2.000_r8*rxt(49)*y(10) + &
                 rxt(61)*y(6) +rxt(64)*y(21) +rxt(72)*y(24) +rxt(77)*y(19) + &
                 rxt(80)*y(27) +rxt(85)*y(30))*y(8) + (rxt(15) +rxt(57) +rxt(50)*y(2)) &
                 *y(13) + (rxt(9) +rxt(48)*y(2))*y(10) + (rxt(11) +rxt(58))*y(11) &
                  + (rxt(17) +rxt(59))*y(14) +rxt(13)*y(12)
         loss(34) = (rxt(48)* y(2) +rxt(49)* y(8) +rxt(55)* y(9) +rxt(71)* y(23) &
                  +rxt(91)* y(36) + rxt(9) + rxt(10) + het_rates(10))* y(10)
         prod(34) = (rxt(11) +rxt(12) +rxt(58))*y(11) + (rxt(51)*y(12) +rxt(52)*y(12)) &
                 *y(2) +rxt(47)*y(9)*y(1) +rxt(14)*y(13)
         loss(16) = ( + rxt(11) + rxt(12) + rxt(58) + het_rates(11))* y(11)
         prod(16) =rxt(55)*y(10)*y(9)
         loss(20) = ((rxt(51) +rxt(52))* y(2) + rxt(13) + het_rates(12))* y(12)
         prod(20) = (rxt(91)*y(36) +rxt(71)*y(23))*y(10) +rxt(53)*y(9)*y(2)
         loss(21) = (rxt(50)* y(2) + rxt(14) + rxt(15) + rxt(57) + het_rates(13)) &
                 * y(13)
         prod(21) =rxt(54)*y(9)*y(3)
         loss(14) = ( + rxt(17) + rxt(59) + het_rates(14))* y(14)
         prod(14) =rxt(56)*y(24)*y(9)
         loss(26) = (rxt(67)* y(3) +rxt(66)* y(6) +rxt(64)* y(8) + 2._r8*rxt(65) &
                 * y(21) + het_rates(21))* y(21)
         prod(26) = (rxt(30)*y(16) +rxt(68)*y(22))*y(2)
         loss(22) = ((rxt(68) +rxt(69))* y(2) + rxt(7) + het_rates(22))* y(22)
         prod(22) = (rxt(67)*y(21) +rxt(81)*y(27) +rxt(86)*y(30))*y(3)
         loss(29) = (rxt(70)* y(2) +rxt(71)* y(10) + rxt(16) + het_rates(23))* y(23)
         prod(29) = (rxt(64)*y(8) +2.000_r8*rxt(65)*y(21) +rxt(66)*y(6))*y(21) &
                  + (.200_r8*rxt(31)*y(17) +rxt(69)*y(22))*y(2) + (rxt(77)*y(8) + &
                 rxt(78)*y(3))*y(19) +.500_r8*rxt(33)*y(18)*y(1) +rxt(7)*y(22)
         loss(30) = (rxt(73)* y(6) +rxt(72)* y(8) +rxt(56)* y(9) + 2._r8*rxt(74) &
                 * y(24) + het_rates(24))* y(24)
         prod(30) = (rxt(75)*y(20) +rxt(76)*y(20) +.500_r8*rxt(70)*y(23))*y(2) &
                  + (rxt(17) +rxt(59))*y(14) +.500_r8*rxt(83)*y(29)*y(1) &
                  +.500_r8*rxt(85)*y(30)*y(8) +rxt(71)*y(23)*y(10)
         loss(17) = (rxt(33)* y(1) +rxt(32)* y(2) + het_rates(18))* y(18)
         prod(17) = 0._r8
         loss(25) = (rxt(78)* y(3) +rxt(79)* y(6) +rxt(77)* y(8) + het_rates(19)) &
                 * y(19)
         prod(25) =rxt(32)*y(18)*y(2)
         loss(18) = (rxt(34)* y(1) +rxt(35)* y(2) + het_rates(25))* y(25)
         prod(18) = 0._r8
         loss(1) = ( + rxt(95) + rxt(97) + rxt(99) + het_rates(26))* y(26)
         prod(1) = 0._r8
         loss(27) = (rxt(81)* y(3) +rxt(82)* y(6) +rxt(80)* y(8) + het_rates(27)) &
                 * y(27)
         prod(27) =rxt(35)*y(25)*y(2)
         loss(2) = ( + rxt(96) + rxt(98) + rxt(100) + het_rates(28))* y(28)
         prod(2) = 0._r8
         loss(28) = (rxt(83)* y(1) +rxt(84)* y(2) + rxt(19) + het_rates(29))* y(29)
         prod(28) = (rxt(80)*y(8) +rxt(82)*y(6))*y(27) +rxt(34)*y(25)*y(1)
         loss(24) = (rxt(86)* y(3) +rxt(85)* y(8) + het_rates(30))* y(30)
         prod(24) =rxt(84)*y(29)*y(2)
         loss(3) = ( + rxt(20) + het_rates(54))* y(54)
         prod(3) = 0._r8
         loss(4) = ( + rxt(21) + het_rates(55))* y(55)
         prod(4) = 0._r8
         loss(5) = ( + rxt(22) + het_rates(56))* y(56)
         prod(5) = 0._r8
         loss(6) = ( + rxt(92) + het_rates(39))* y(39)
         prod(6) = 0._r8
         loss(7) = ( + rxt(93) + het_rates(40))* y(40)
         prod(7) =1.150_r8*rxt(92)*y(39)
         loss(8) = ( + rxt(94) + het_rates(41))* y(41)
         prod(8) =1.150_r8*rxt(93)*y(40)
         loss(9) = ( + het_rates(42))* y(42)
         prod(9) =.0436005_r8*rxt(98)*y(28) +.460_r8*rxt(104)*y(43)
         loss(10) = ( + rxt(104) + het_rates(43))* y(43)
         prod(10) = (.0027205_r8*rxt(95) +.0002705_r8*rxt(99))*y(26) &
                  + (.1498605_r8*rxt(96) +.0103505_r8*rxt(98))*y(28) +.460_r8*rxt(103) &
                 *y(44)
         loss(11) = ( + rxt(103) + het_rates(44))* y(44)
         prod(11) = (.0098105_r8*rxt(95) +.0032705_r8*rxt(97) +.0406005_r8*rxt(99)) &
                 *y(26) + (.0545005_r8*rxt(96) +.0980905_r8*rxt(98) + &
                 .1749305_r8*rxt(100))*y(28) +.460_r8*rxt(94)*y(41) +.460_r8*rxt(102) &
                 *y(45)
         loss(12) = ( + rxt(102) + het_rates(45))* y(45)
         prod(12) = (.0926405_r8*rxt(96) +.0163505_r8*rxt(98) +.5901905_r8*rxt(100)) &
                 *y(28) + (.0021805_r8*rxt(95) +.0645805_r8*rxt(99))*y(26) &
                  +.460_r8*rxt(101)*y(46)
         loss(13) = ( + rxt(101) + het_rates(46))* y(46)
         prod(13) =.500_r8*rxt(94)*y(41) +.500_r8*rxt(104)*y(43) +.500_r8*rxt(103) &
                 *y(44) +.500_r8*rxt(102)*y(45) +.500_r8*rxt(101)*y(46)

      end subroutine imp_prod_loss

      end module mo_prod_loss
