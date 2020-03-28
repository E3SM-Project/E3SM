













      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: nlnmat

      contains

      subroutine     nlnmat01( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(inout) ::  mat(nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------


         mat(256) = -(rxt(26)*y(19) + rxt(27)*y(26) + rxt(30)*y(2) + rxt(32)*y(3) &
                      + rxt(38)*y(8) + rxt(40)*y(9) + rxt(75)*y(28))
         mat(50) = -rxt(26)*y(1)
         mat(56) = -rxt(27)*y(1)
         mat(241) = -rxt(30)*y(1)
         mat(216) = -rxt(32)*y(1)
         mat(159) = -rxt(38)*y(1)
         mat(186) = -rxt(40)*y(1)
         mat(128) = -rxt(75)*y(1)

         mat(241) = mat(241) + 2.000_r8*rxt(31)*y(2)

         mat(240) = -(rxt(25)*y(19) + rxt(28)*y(26) + rxt(29)*y(5) + rxt(30)*y(1) &
                      + (4._r8*rxt(31) + 4._r8*rxt(36)) * y(2) + rxt(33)*y(3) + rxt(37) &
                      *y(4) + rxt(41)*y(10) + (rxt(43) + rxt(82)) * y(13) + rxt(44) &
                      *y(12) + rxt(45)*y(9) + rxt(55)*y(7) + rxt(60)*y(23) + rxt(62) &
                      *y(24) + rxt(76)*y(28) + (rxt(83) + rxt(84)) * y(33) + rxt(85) &
                      *y(34))
         mat(49) = -rxt(25)*y(2)
         mat(55) = -rxt(28)*y(2)
         mat(78) = -rxt(29)*y(2)
         mat(255) = -rxt(30)*y(2)
         mat(215) = -rxt(33)*y(2)
         mat(37) = -rxt(37)*y(2)
         mat(196) = -rxt(41)*y(2)
         mat(65) = -(rxt(43) + rxt(82)) * y(2)
         mat(60) = -rxt(44)*y(2)
         mat(185) = -rxt(45)*y(2)
         mat(70) = -rxt(55)*y(2)
         mat(75) = -rxt(60)*y(2)
         mat(136) = -rxt(62)*y(2)
         mat(127) = -rxt(76)*y(2)
         mat(34) = -(rxt(83) + rxt(84)) * y(2)
         mat(30) = -rxt(85)*y(2)

         mat(255) = mat(255) + rxt(32)*y(3) + rxt(27)*y(26)
         mat(240) = mat(240) + .300_r8*rxt(55)*y(7)
         mat(215) = mat(215) + rxt(32)*y(1) + rxt(39)*y(8) + rxt(70)*y(20)
         mat(70) = mat(70) + .300_r8*rxt(55)*y(2)
         mat(158) = rxt(39)*y(3)
         mat(96) = rxt(70)*y(3)
         mat(55) = mat(55) + rxt(27)*y(1)

         mat(214) = -(rxt(32)*y(1) + rxt(33)*y(2) + (4._r8*rxt(34) + 4._r8*rxt(35) &
                      ) * y(3) + rxt(39)*y(8) + rxt(46)*y(9) + rxt(52)*y(6) + rxt(59) &
                      *y(22) + rxt(70)*y(20) + rxt(73)*y(27) + rxt(78)*y(29))
         mat(254) = -rxt(32)*y(3)
         mat(239) = -rxt(33)*y(3)
         mat(157) = -rxt(39)*y(3)
         mat(184) = -rxt(46)*y(3)
         mat(172) = -rxt(52)*y(3)
         mat(104) = -rxt(59)*y(3)
         mat(95) = -rxt(70)*y(3)
         mat(115) = -rxt(73)*y(3)
         mat(87) = -rxt(78)*y(3)

         mat(254) = mat(254) + rxt(30)*y(2)
         mat(239) = mat(239) + rxt(30)*y(1) + rxt(37)*y(4) + rxt(29)*y(5) + rxt(41) &
                      *y(10) + .500_r8*rxt(84)*y(33)
         mat(36) = rxt(37)*y(2)
         mat(77) = rxt(29)*y(2)
         mat(172) = mat(172) + 4.000_r8*rxt(54)*y(6) + rxt(53)*y(8) + 2.000_r8*rxt(58) &
                      *y(22) + rxt(65)*y(25) + rxt(71)*y(20) + 2.000_r8*rxt(74)*y(27)
         mat(157) = mat(157) + rxt(53)*y(6) + rxt(56)*y(22)
         mat(195) = rxt(41)*y(2)
         mat(104) = mat(104) + 2.000_r8*rxt(58)*y(6) + rxt(56)*y(8) + 4.000_r8*rxt(57) &
                      *y(22)
         mat(143) = rxt(65)*y(6)
         mat(95) = mat(95) + rxt(71)*y(6)
         mat(115) = mat(115) + 2.000_r8*rxt(74)*y(6)
         mat(33) = .500_r8*rxt(84)*y(2)

         mat(35) = -(rxt(37)*y(2))
         mat(220) = -rxt(37)*y(4)

         mat(220) = mat(220) + 2.000_r8*rxt(36)*y(2)
         mat(198) = (2.000_r8*rxt(34)+2.000_r8*rxt(35))*y(3)

         mat(76) = -(rxt(29)*y(2))
         mat(227) = -rxt(29)*y(5)

         mat(244) = rxt(26)*y(19) + rxt(27)*y(26) + rxt(75)*y(28)
         mat(227) = mat(227) + .300_r8*rxt(55)*y(7)
         mat(202) = rxt(70)*y(20)
         mat(161) = 4.000_r8*rxt(54)*y(6) + rxt(53)*y(8) + rxt(58)*y(22) + rxt(65) &
                      *y(25) + 2.000_r8*rxt(71)*y(20) + 2.000_r8*rxt(74)*y(27)
         mat(67) = .300_r8*rxt(55)*y(2)
         mat(145) = rxt(53)*y(6) + rxt(69)*y(20) + rxt(72)*y(27) + rxt(77)*y(29)
         mat(98) = rxt(58)*y(6)
         mat(138) = rxt(65)*y(6)
         mat(46) = rxt(26)*y(1)
         mat(89) = rxt(70)*y(3) + 2.000_r8*rxt(71)*y(6) + rxt(69)*y(8)
         mat(52) = rxt(27)*y(1)
         mat(107) = 2.000_r8*rxt(74)*y(6) + rxt(72)*y(8)
         mat(117) = rxt(75)*y(1)
         mat(80) = rxt(77)*y(8)

         mat(169) = -(rxt(52)*y(3) + rxt(53)*y(8) + 4._r8*rxt(54)*y(6) + rxt(58)*y(22) &
                      + rxt(71)*y(20) + rxt(74)*y(27))
         mat(211) = -rxt(52)*y(6)
         mat(154) = -rxt(53)*y(6)
         mat(102) = -rxt(58)*y(6)
         mat(93) = -rxt(71)*y(6)
         mat(113) = -rxt(74)*y(6)

         mat(251) = rxt(75)*y(28)
         mat(236) = .700_r8*rxt(55)*y(7) + .500_r8*rxt(62)*y(24)
         mat(68) = .700_r8*rxt(55)*y(2)
         mat(154) = mat(154) + rxt(64)*y(25)
         mat(132) = .500_r8*rxt(62)*y(2)
         mat(141) = rxt(64)*y(8) + 4.000_r8*rxt(66)*y(25)
         mat(124) = rxt(75)*y(1)

         mat(66) = -(rxt(55)*y(2))
         mat(225) = -rxt(55)*y(7)

         mat(200) = rxt(52)*y(6)
         mat(160) = rxt(52)*y(3)

         mat(153) = -(rxt(38)*y(1) + rxt(39)*y(3) + rxt(42)*y(10) + rxt(53)*y(6) &
                      + rxt(56)*y(22) + rxt(64)*y(25) + rxt(69)*y(20) + rxt(72)*y(27) &
                      + rxt(77)*y(29))
         mat(250) = -rxt(38)*y(8)
         mat(210) = -rxt(39)*y(8)
         mat(191) = -rxt(42)*y(8)
         mat(168) = -rxt(53)*y(8)
         mat(101) = -rxt(56)*y(8)
         mat(140) = -rxt(64)*y(8)
         mat(92) = -rxt(69)*y(8)
         mat(112) = -rxt(72)*y(8)
         mat(85) = -rxt(77)*y(8)

         mat(182) = -(rxt(40)*y(1) + rxt(45)*y(2) + rxt(46)*y(3) + rxt(47)*y(10) &
                      + rxt(48)*y(25))
         mat(252) = -rxt(40)*y(9)
         mat(237) = -rxt(45)*y(9)
         mat(212) = -rxt(46)*y(9)
         mat(193) = -rxt(47)*y(9)
         mat(142) = -rxt(48)*y(9)

         mat(252) = mat(252) + rxt(38)*y(8)
         mat(237) = mat(237) + rxt(41)*y(10) + (rxt(43)+rxt(82))*y(13)
         mat(212) = mat(212) + rxt(39)*y(8)
         mat(170) = rxt(53)*y(8)
         mat(155) = rxt(38)*y(1) + rxt(39)*y(3) + rxt(53)*y(6) + 2.000_r8*rxt(42) &
                      *y(10) + rxt(56)*y(22) + rxt(64)*y(25) + rxt(69)*y(20) + rxt(72) &
                      *y(27) + rxt(77)*y(29)
         mat(193) = mat(193) + rxt(41)*y(2) + 2.000_r8*rxt(42)*y(8)
         mat(62) = (rxt(43)+rxt(82))*y(2)
         mat(103) = rxt(56)*y(8)
         mat(142) = mat(142) + rxt(64)*y(8)
         mat(94) = rxt(69)*y(8)
         mat(114) = rxt(72)*y(8)
         mat(86) = rxt(77)*y(8)

         mat(194) = -(rxt(41)*y(2) + rxt(42)*y(8) + rxt(47)*y(9) + rxt(63)*y(24))
         mat(238) = -rxt(41)*y(10)
         mat(156) = -rxt(42)*y(10)
         mat(183) = -rxt(47)*y(10)
         mat(134) = -rxt(63)*y(10)

         mat(253) = rxt(40)*y(9)
         mat(238) = mat(238) + rxt(44)*y(12)
         mat(183) = mat(183) + rxt(40)*y(1)
         mat(59) = rxt(44)*y(2)


         mat(176) = rxt(47)*y(10)
         mat(187) = rxt(47)*y(9)

         mat(57) = -(rxt(44)*y(2))
         mat(223) = -rxt(44)*y(12)

         mat(223) = mat(223) + rxt(45)*y(9)
         mat(177) = rxt(45)*y(2)
         mat(188) = rxt(63)*y(24)
         mat(129) = rxt(63)*y(10)

         mat(61) = -((rxt(43) + rxt(82)) * y(2))
         mat(224) = -(rxt(43) + rxt(82)) * y(13)

         mat(199) = rxt(46)*y(9)
         mat(178) = rxt(46)*y(3)


         mat(175) = rxt(48)*y(25)
         mat(137) = rxt(48)*y(9)

         mat(99) = -(rxt(56)*y(8) + 4._r8*rxt(57)*y(22) + rxt(58)*y(6) + rxt(59)*y(3))
         mat(148) = -rxt(56)*y(22)
         mat(163) = -rxt(58)*y(22)
         mat(205) = -rxt(59)*y(22)

         mat(230) = rxt(60)*y(23)
         mat(72) = rxt(60)*y(2)

         mat(71) = -((rxt(60) + rxt(61)) * y(2))
         mat(226) = -(rxt(60) + rxt(61)) * y(23)

         mat(201) = rxt(59)*y(22) + rxt(73)*y(27) + rxt(78)*y(29)
         mat(97) = rxt(59)*y(3)
         mat(106) = rxt(73)*y(3)
         mat(79) = rxt(78)*y(3)

         mat(130) = -(rxt(62)*y(2) + rxt(63)*y(10))
         mat(233) = -rxt(62)*y(24)
         mat(189) = -rxt(63)*y(24)

         mat(248) = .500_r8*rxt(26)*y(19)
         mat(233) = mat(233) + rxt(61)*y(23)
         mat(208) = rxt(70)*y(20)
         mat(166) = rxt(58)*y(22)
         mat(151) = rxt(56)*y(22) + rxt(69)*y(20)
         mat(100) = rxt(58)*y(6) + rxt(56)*y(8) + 4.000_r8*rxt(57)*y(22)
         mat(73) = rxt(61)*y(2)
         mat(48) = .500_r8*rxt(26)*y(1)
         mat(91) = rxt(70)*y(3) + rxt(69)*y(8)

         mat(139) = -(rxt(48)*y(9) + rxt(64)*y(8) + rxt(65)*y(6) + 4._r8*rxt(66)*y(25))
         mat(179) = -rxt(48)*y(25)
         mat(152) = -rxt(64)*y(25)
         mat(167) = -rxt(65)*y(25)

         mat(249) = .500_r8*rxt(75)*y(28)
         mat(234) = .500_r8*rxt(62)*y(24)
         mat(152) = mat(152) + .500_r8*rxt(77)*y(29)
         mat(190) = rxt(63)*y(24)
         mat(131) = .500_r8*rxt(62)*y(2) + rxt(63)*y(10)
         mat(122) = .500_r8*rxt(75)*y(1)
         mat(84) = .500_r8*rxt(77)*y(8)

         mat(45) = -(rxt(25)*y(2) + rxt(26)*y(1))
         mat(221) = -rxt(25)*y(19)
         mat(242) = -rxt(26)*y(19)

         mat(90) = -(rxt(69)*y(8) + rxt(70)*y(3) + rxt(71)*y(6))
         mat(147) = -rxt(69)*y(20)
         mat(204) = -rxt(70)*y(20)
         mat(162) = -rxt(71)*y(20)

         mat(229) = rxt(25)*y(19)
         mat(47) = rxt(25)*y(2)


      end subroutine     nlnmat01

      subroutine     nlnmat02( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(inout) ::  mat(nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------


         mat(51) = -(rxt(27)*y(1) + rxt(28)*y(2))
         mat(243) = -rxt(27)*y(26)
         mat(222) = -rxt(28)*y(26)

         mat(109) = -(rxt(72)*y(8) + rxt(73)*y(3) + rxt(74)*y(6))
         mat(149) = -rxt(72)*y(27)
         mat(206) = -rxt(73)*y(27)
         mat(164) = -rxt(74)*y(27)

         mat(231) = rxt(28)*y(26)
         mat(53) = rxt(28)*y(2)

         mat(120) = -(rxt(75)*y(1) + rxt(76)*y(2))
         mat(247) = -rxt(75)*y(28)
         mat(232) = -rxt(76)*y(28)

         mat(247) = mat(247) + rxt(27)*y(26)
         mat(165) = rxt(74)*y(27)
         mat(150) = rxt(72)*y(27)
         mat(54) = rxt(27)*y(1)
         mat(110) = rxt(74)*y(6) + rxt(72)*y(8)

         mat(81) = -(rxt(77)*y(8) + rxt(78)*y(3))
         mat(146) = -rxt(77)*y(29)
         mat(203) = -rxt(78)*y(29)

         mat(228) = rxt(76)*y(28)
         mat(118) = rxt(76)*y(2)

         mat(32) = -((rxt(83) + rxt(84)) * y(2))
         mat(219) = -(rxt(83) + rxt(84)) * y(33)

         mat(29) = -(rxt(85)*y(2))
         mat(218) = -rxt(85)*y(34)

         mat(218) = mat(218) + (rxt(83)+.500_r8*rxt(84))*y(33)
         mat(31) = (rxt(83)+.500_r8*rxt(84))*y(2)


         mat(217) = rxt(85)*y(34)
         mat(28) = rxt(85)*y(2)




























      end subroutine     nlnmat02

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
         mat(  29) = mat(  29) + lmat(  29)
         mat(  32) = mat(  32) + lmat(  32)
         mat(  35) = mat(  35) + lmat(  35)
         mat(  37) = mat(  37) + lmat(  37)
         mat(  38) = lmat(  38)
         mat(  39) = lmat(  39)
         mat(  40) = lmat(  40)
         mat(  41) = lmat(  41)
         mat(  42) = lmat(  42)
         mat(  43) = lmat(  43)
         mat(  44) = lmat(  44)
         mat(  45) = mat(  45) + lmat(  45)
         mat(  51) = mat(  51) + lmat(  51)
         mat(  57) = mat(  57) + lmat(  57)
         mat(  58) = lmat(  58)
         mat(  60) = mat(  60) + lmat(  60)
         mat(  61) = mat(  61) + lmat(  61)
         mat(  62) = mat(  62) + lmat(  62)
         mat(  63) = lmat(  63)
         mat(  64) = lmat(  64)
         mat(  65) = mat(  65) + lmat(  65)
         mat(  66) = mat(  66) + lmat(  66)
         mat(  67) = mat(  67) + lmat(  67)
         mat(  69) = lmat(  69)
         mat(  70) = mat(  70) + lmat(  70)
         mat(  71) = mat(  71) + lmat(  71)
         mat(  73) = mat(  73) + lmat(  73)
         mat(  74) = lmat(  74)
         mat(  75) = mat(  75) + lmat(  75)
         mat(  76) = mat(  76) + lmat(  76)
         mat(  77) = mat(  77) + lmat(  77)
         mat(  81) = mat(  81) + lmat(  81)
         mat(  90) = mat(  90) + lmat(  90)
         mat(  99) = mat(  99) + lmat(  99)
         mat( 109) = mat( 109) + lmat( 109)
         mat( 117) = mat( 117) + lmat( 117)
         mat( 120) = mat( 120) + lmat( 120)
         mat( 124) = mat( 124) + lmat( 124)
         mat( 126) = lmat( 126)
         mat( 130) = mat( 130) + lmat( 130)
         mat( 132) = mat( 132) + lmat( 132)
         mat( 135) = lmat( 135)
         mat( 139) = mat( 139) + lmat( 139)
         mat( 153) = mat( 153) + lmat( 153)
         mat( 169) = mat( 169) + lmat( 169)
         mat( 180) = lmat( 180)
         mat( 182) = mat( 182) + lmat( 182)
         mat( 186) = mat( 186) + lmat( 186)
         mat( 188) = mat( 188) + lmat( 188)
         mat( 191) = mat( 191) + lmat( 191)
         mat( 193) = mat( 193) + lmat( 193)
         mat( 194) = mat( 194) + lmat( 194)
         mat( 197) = lmat( 197)
         mat( 214) = mat( 214) + lmat( 214)
         mat( 227) = mat( 227) + lmat( 227)
         mat( 230) = mat( 230) + lmat( 230)
         mat( 233) = mat( 233) + lmat( 233)
         mat( 234) = mat( 234) + lmat( 234)
         mat( 236) = mat( 236) + lmat( 236)
         mat( 239) = mat( 239) + lmat( 239)
         mat( 240) = mat( 240) + lmat( 240)
         mat( 254) = mat( 254) + lmat( 254)
         mat( 255) = mat( 255) + lmat( 255)
         mat( 256) = mat( 256) + lmat( 256)
         mat(  82) = 0._r8
         mat(  83) = 0._r8
         mat(  88) = 0._r8
         mat( 105) = 0._r8
         mat( 108) = 0._r8
         mat( 111) = 0._r8
         mat( 116) = 0._r8
         mat( 119) = 0._r8
         mat( 121) = 0._r8
         mat( 123) = 0._r8
         mat( 125) = 0._r8
         mat( 133) = 0._r8
         mat( 144) = 0._r8
         mat( 171) = 0._r8
         mat( 173) = 0._r8
         mat( 174) = 0._r8
         mat( 181) = 0._r8
         mat( 192) = 0._r8
         mat( 207) = 0._r8
         mat( 209) = 0._r8
         mat( 213) = 0._r8
         mat( 235) = 0._r8
         mat( 245) = 0._r8
         mat( 246) = 0._r8
         mat(   1) = mat(   1) - dti
         mat(   2) = mat(   2) - dti
         mat(   3) = mat(   3) - dti
         mat(   4) = mat(   4) - dti
         mat(   5) = mat(   5) - dti
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
         mat(  29) = mat(  29) - dti
         mat(  32) = mat(  32) - dti
         mat(  35) = mat(  35) - dti
         mat(  38) = mat(  38) - dti
         mat(  41) = mat(  41) - dti
         mat(  45) = mat(  45) - dti
         mat(  51) = mat(  51) - dti
         mat(  57) = mat(  57) - dti
         mat(  61) = mat(  61) - dti
         mat(  66) = mat(  66) - dti
         mat(  71) = mat(  71) - dti
         mat(  76) = mat(  76) - dti
         mat(  81) = mat(  81) - dti
         mat(  90) = mat(  90) - dti
         mat(  99) = mat(  99) - dti
         mat( 109) = mat( 109) - dti
         mat( 120) = mat( 120) - dti
         mat( 130) = mat( 130) - dti
         mat( 139) = mat( 139) - dti
         mat( 153) = mat( 153) - dti
         mat( 169) = mat( 169) - dti
         mat( 182) = mat( 182) - dti
         mat( 194) = mat( 194) - dti
         mat( 214) = mat( 214) - dti
         mat( 240) = mat( 240) - dti
         mat( 256) = mat( 256) - dti

      end subroutine nlnmat_finit

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

      call     nlnmat01( mat, y, rxt )
      call     nlnmat02( mat, y, rxt )
      call     nlnmat_finit( mat, lmat, dti )

      end subroutine nlnmat

      end module mo_nln_matrix
