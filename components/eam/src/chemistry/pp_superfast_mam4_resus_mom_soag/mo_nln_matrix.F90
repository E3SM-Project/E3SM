













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


         mat(261) = -(rxt(30)*y(18) + rxt(31)*y(25) + rxt(34)*y(2) + rxt(36)*y(3) &
                      + rxt(42)*y(8) + rxt(44)*y(9) + rxt(80)*y(27))
         mat(53) = -rxt(30)*y(1)
         mat(59) = -rxt(31)*y(1)
         mat(246) = -rxt(34)*y(1)
         mat(206) = -rxt(36)*y(1)
         mat(162) = -rxt(42)*y(1)
         mat(187) = -rxt(44)*y(1)
         mat(131) = -rxt(80)*y(1)

         mat(246) = mat(246) + 2.000_r8*rxt(35)*y(2)

         mat(245) = -(rxt(29)*y(18) + rxt(32)*y(25) + rxt(33)*y(5) + rxt(34)*y(1) &
                      + (4._r8*rxt(35) + 4._r8*rxt(40)) * y(2) + rxt(37)*y(3) + rxt(41) &
                      *y(4) + rxt(45)*y(10) + rxt(47)*y(13) + (rxt(48) + rxt(49) &
                      ) * y(12) + rxt(50)*y(9) + rxt(60)*y(7) + rxt(65)*y(22) + rxt(67) &
                      *y(23) + rxt(81)*y(27) + (rxt(85) + rxt(86)) * y(31) + rxt(87) &
                      *y(32))
         mat(52) = -rxt(29)*y(2)
         mat(58) = -rxt(32)*y(2)
         mat(81) = -rxt(33)*y(2)
         mat(260) = -rxt(34)*y(2)
         mat(205) = -rxt(37)*y(2)
         mat(33) = -rxt(41)*y(2)
         mat(174) = -rxt(45)*y(2)
         mat(68) = -rxt(47)*y(2)
         mat(63) = -(rxt(48) + rxt(49)) * y(2)
         mat(186) = -rxt(50)*y(2)
         mat(73) = -rxt(60)*y(2)
         mat(78) = -rxt(65)*y(2)
         mat(139) = -rxt(67)*y(2)
         mat(130) = -rxt(81)*y(2)
         mat(47) = -(rxt(85) + rxt(86)) * y(2)
         mat(30) = -rxt(87)*y(2)

         mat(260) = mat(260) + rxt(36)*y(3) + rxt(31)*y(25)
         mat(245) = mat(245) + .300_r8*rxt(60)*y(7)
         mat(205) = mat(205) + rxt(36)*y(1) + rxt(43)*y(8) + rxt(75)*y(19)
         mat(73) = mat(73) + .300_r8*rxt(60)*y(2)
         mat(161) = rxt(43)*y(3)
         mat(99) = rxt(75)*y(3)
         mat(58) = mat(58) + rxt(31)*y(1)

         mat(203) = -(rxt(36)*y(1) + rxt(37)*y(2) + (4._r8*rxt(38) + 4._r8*rxt(39) &
                      ) * y(3) + rxt(43)*y(8) + rxt(51)*y(9) + rxt(57)*y(6) + rxt(64) &
                      *y(21) + rxt(75)*y(19) + rxt(78)*y(26) + rxt(83)*y(28))
         mat(258) = -rxt(36)*y(3)
         mat(243) = -rxt(37)*y(3)
         mat(159) = -rxt(43)*y(3)
         mat(184) = -rxt(51)*y(3)
         mat(218) = -rxt(57)*y(3)
         mat(106) = -rxt(64)*y(3)
         mat(97) = -rxt(75)*y(3)
         mat(117) = -rxt(78)*y(3)
         mat(90) = -rxt(83)*y(3)

         mat(258) = mat(258) + rxt(34)*y(2)
         mat(243) = mat(243) + rxt(34)*y(1) + rxt(41)*y(4) + rxt(33)*y(5) + rxt(45) &
                      *y(10) + .500_r8*rxt(86)*y(31)
         mat(32) = rxt(41)*y(2)
         mat(80) = rxt(33)*y(2)
         mat(218) = mat(218) + 4.000_r8*rxt(59)*y(6) + rxt(58)*y(8) + 2.000_r8*rxt(63) &
                      *y(21) + rxt(70)*y(24) + rxt(76)*y(19) + 2.000_r8*rxt(79)*y(26)
         mat(159) = mat(159) + rxt(58)*y(6) + rxt(61)*y(21)
         mat(172) = rxt(45)*y(2)
         mat(106) = mat(106) + 2.000_r8*rxt(63)*y(6) + rxt(61)*y(8) + 4.000_r8*rxt(62) &
                      *y(21)
         mat(145) = rxt(70)*y(6)
         mat(97) = mat(97) + rxt(76)*y(6)
         mat(117) = mat(117) + 2.000_r8*rxt(79)*y(6)
         mat(46) = .500_r8*rxt(86)*y(2)

         mat(31) = -(rxt(41)*y(2))
         mat(224) = -rxt(41)*y(4)

         mat(224) = mat(224) + 2.000_r8*rxt(40)*y(2)
         mat(188) = (2.000_r8*rxt(38)+2.000_r8*rxt(39))*y(3)

         mat(79) = -(rxt(33)*y(2))
         mat(232) = -rxt(33)*y(5)

         mat(249) = rxt(30)*y(18) + rxt(31)*y(25) + rxt(80)*y(27)
         mat(232) = mat(232) + .300_r8*rxt(60)*y(7)
         mat(192) = rxt(75)*y(19)
         mat(208) = 4.000_r8*rxt(59)*y(6) + rxt(58)*y(8) + rxt(63)*y(21) + rxt(70) &
                      *y(24) + 2.000_r8*rxt(76)*y(19) + 2.000_r8*rxt(79)*y(26)
         mat(70) = .300_r8*rxt(60)*y(2)
         mat(148) = rxt(58)*y(6) + rxt(74)*y(19) + rxt(77)*y(26) + rxt(82)*y(28)
         mat(101) = rxt(63)*y(6)
         mat(141) = rxt(70)*y(6)
         mat(49) = rxt(30)*y(1)
         mat(92) = rxt(75)*y(3) + 2.000_r8*rxt(76)*y(6) + rxt(74)*y(8)
         mat(55) = rxt(31)*y(1)
         mat(110) = 2.000_r8*rxt(79)*y(6) + rxt(77)*y(8)
         mat(120) = rxt(80)*y(1)
         mat(83) = rxt(82)*y(8)

         mat(219) = -(rxt(57)*y(3) + rxt(58)*y(8) + 4._r8*rxt(59)*y(6) + rxt(63)*y(21) &
                      + rxt(76)*y(19) + rxt(79)*y(26))
         mat(204) = -rxt(57)*y(6)
         mat(160) = -rxt(58)*y(6)
         mat(107) = -rxt(63)*y(6)
         mat(98) = -rxt(76)*y(6)
         mat(118) = -rxt(79)*y(6)

         mat(259) = rxt(80)*y(27)
         mat(244) = .700_r8*rxt(60)*y(7) + .500_r8*rxt(67)*y(23)
         mat(72) = .700_r8*rxt(60)*y(2)
         mat(160) = mat(160) + rxt(69)*y(24)
         mat(138) = .500_r8*rxt(67)*y(2)
         mat(146) = rxt(69)*y(8) + 4.000_r8*rxt(71)*y(24)
         mat(129) = rxt(80)*y(1)

         mat(69) = -(rxt(60)*y(2))
         mat(230) = -rxt(60)*y(7)

         mat(190) = rxt(57)*y(6)
         mat(207) = rxt(57)*y(3)

         mat(156) = -(rxt(42)*y(1) + rxt(43)*y(3) + rxt(46)*y(10) + rxt(58)*y(6) &
                      + rxt(61)*y(21) + rxt(69)*y(24) + rxt(74)*y(19) + rxt(77)*y(26) &
                      + rxt(82)*y(28))
         mat(255) = -rxt(42)*y(8)
         mat(200) = -rxt(43)*y(8)
         mat(169) = -rxt(46)*y(8)
         mat(215) = -rxt(58)*y(8)
         mat(104) = -rxt(61)*y(8)
         mat(143) = -rxt(69)*y(8)
         mat(95) = -rxt(74)*y(8)
         mat(115) = -rxt(77)*y(8)
         mat(88) = -rxt(82)*y(8)

         mat(183) = -(rxt(44)*y(1) + rxt(50)*y(2) + rxt(51)*y(3) + rxt(52)*y(10) &
                      + rxt(53)*y(24))
         mat(257) = -rxt(44)*y(9)
         mat(242) = -rxt(50)*y(9)
         mat(202) = -rxt(51)*y(9)
         mat(171) = -rxt(52)*y(9)
         mat(144) = -rxt(53)*y(9)

         mat(257) = mat(257) + rxt(42)*y(8)
         mat(242) = mat(242) + rxt(45)*y(10) + rxt(47)*y(13)
         mat(202) = mat(202) + rxt(43)*y(8)
         mat(217) = rxt(58)*y(8)
         mat(158) = rxt(42)*y(1) + rxt(43)*y(3) + rxt(58)*y(6) + 2.000_r8*rxt(46) &
                      *y(10) + rxt(61)*y(21) + rxt(69)*y(24) + rxt(74)*y(19) + rxt(77) &
                      *y(26) + rxt(82)*y(28)
         mat(171) = mat(171) + rxt(45)*y(2) + 2.000_r8*rxt(46)*y(8)
         mat(66) = rxt(47)*y(2)
         mat(105) = rxt(61)*y(8)
         mat(144) = mat(144) + rxt(69)*y(8)
         mat(96) = rxt(74)*y(8)
         mat(116) = rxt(77)*y(8)
         mat(89) = rxt(82)*y(8)

         mat(170) = -(rxt(45)*y(2) + rxt(46)*y(8) + rxt(52)*y(9) + rxt(68)*y(23) &
                      + rxt(88)*y(31))
         mat(241) = -rxt(45)*y(10)
         mat(157) = -rxt(46)*y(10)
         mat(182) = -rxt(52)*y(10)
         mat(135) = -rxt(68)*y(10)
         mat(45) = -rxt(88)*y(10)

         mat(256) = rxt(44)*y(9)
         mat(241) = mat(241) + (rxt(48)+rxt(49))*y(12)
         mat(182) = mat(182) + rxt(44)*y(1)
         mat(61) = (rxt(48)+rxt(49))*y(2)


         mat(177) = rxt(52)*y(10)
         mat(164) = rxt(52)*y(9)

         mat(60) = -((rxt(48) + rxt(49)) * y(2))
         mat(228) = -(rxt(48) + rxt(49)) * y(12)

         mat(228) = mat(228) + rxt(50)*y(9)
         mat(178) = rxt(50)*y(2)
         mat(166) = rxt(68)*y(23) + rxt(88)*y(31)
         mat(132) = rxt(68)*y(10)
         mat(44) = rxt(88)*y(10)

         mat(64) = -(rxt(47)*y(2))
         mat(229) = -rxt(47)*y(13)

         mat(189) = rxt(51)*y(9)
         mat(179) = rxt(51)*y(3)


         mat(176) = rxt(53)*y(24)
         mat(140) = rxt(53)*y(9)

         mat(102) = -(rxt(61)*y(8) + 4._r8*rxt(62)*y(21) + rxt(63)*y(6) + rxt(64)*y(3))
         mat(151) = -rxt(61)*y(21)
         mat(210) = -rxt(63)*y(21)
         mat(195) = -rxt(64)*y(21)

         mat(235) = rxt(65)*y(22)
         mat(75) = rxt(65)*y(2)

         mat(74) = -((rxt(65) + rxt(66)) * y(2))
         mat(231) = -(rxt(65) + rxt(66)) * y(22)

         mat(191) = rxt(64)*y(21) + rxt(78)*y(26) + rxt(83)*y(28)
         mat(100) = rxt(64)*y(3)
         mat(109) = rxt(78)*y(3)
         mat(82) = rxt(83)*y(3)

         mat(133) = -(rxt(67)*y(2) + rxt(68)*y(10))
         mat(238) = -rxt(67)*y(23)
         mat(167) = -rxt(68)*y(23)

         mat(253) = .500_r8*rxt(30)*y(18)
         mat(238) = mat(238) + rxt(66)*y(22)
         mat(198) = rxt(75)*y(19)
         mat(213) = rxt(63)*y(21)
         mat(154) = rxt(61)*y(21) + rxt(74)*y(19)
         mat(103) = rxt(63)*y(6) + rxt(61)*y(8) + 4.000_r8*rxt(62)*y(21)
         mat(76) = rxt(66)*y(2)
         mat(51) = .500_r8*rxt(30)*y(1)
         mat(94) = rxt(75)*y(3) + rxt(74)*y(8)

         mat(142) = -(rxt(53)*y(9) + rxt(69)*y(8) + rxt(70)*y(6) + 4._r8*rxt(71)*y(24))
         mat(180) = -rxt(53)*y(24)
         mat(155) = -rxt(69)*y(24)
         mat(214) = -rxt(70)*y(24)

         mat(254) = .500_r8*rxt(80)*y(27)
         mat(239) = .500_r8*rxt(67)*y(23)
         mat(155) = mat(155) + .500_r8*rxt(82)*y(28)
         mat(168) = rxt(68)*y(23)
         mat(134) = .500_r8*rxt(67)*y(2) + rxt(68)*y(10)
         mat(125) = .500_r8*rxt(80)*y(1)
         mat(87) = .500_r8*rxt(82)*y(8)

         mat(48) = -(rxt(29)*y(2) + rxt(30)*y(1))
         mat(226) = -rxt(29)*y(18)
         mat(247) = -rxt(30)*y(18)

         mat(93) = -(rxt(74)*y(8) + rxt(75)*y(3) + rxt(76)*y(6))
         mat(150) = -rxt(74)*y(19)
         mat(194) = -rxt(75)*y(19)
         mat(209) = -rxt(76)*y(19)

         mat(234) = rxt(29)*y(18)
         mat(50) = rxt(29)*y(2)


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


         mat(54) = -(rxt(31)*y(1) + rxt(32)*y(2))
         mat(248) = -rxt(31)*y(25)
         mat(227) = -rxt(32)*y(25)

         mat(112) = -(rxt(77)*y(8) + rxt(78)*y(3) + rxt(79)*y(6))
         mat(152) = -rxt(77)*y(26)
         mat(196) = -rxt(78)*y(26)
         mat(211) = -rxt(79)*y(26)

         mat(236) = rxt(32)*y(25)
         mat(56) = rxt(32)*y(2)

         mat(123) = -(rxt(80)*y(1) + rxt(81)*y(2))
         mat(252) = -rxt(80)*y(27)
         mat(237) = -rxt(81)*y(27)

         mat(252) = mat(252) + rxt(31)*y(25)
         mat(212) = rxt(79)*y(26)
         mat(153) = rxt(77)*y(26)
         mat(57) = rxt(31)*y(1)
         mat(113) = rxt(79)*y(6) + rxt(77)*y(8)

         mat(84) = -(rxt(82)*y(8) + rxt(83)*y(3))
         mat(149) = -rxt(82)*y(28)
         mat(193) = -rxt(83)*y(28)

         mat(233) = rxt(81)*y(27)
         mat(121) = rxt(81)*y(2)

         mat(43) = -((rxt(85) + rxt(86)) * y(2) + rxt(88)*y(10))
         mat(225) = -(rxt(85) + rxt(86)) * y(31)
         mat(165) = -rxt(88)*y(31)

         mat(29) = -(rxt(87)*y(2))
         mat(223) = -rxt(87)*y(32)

         mat(223) = mat(223) + (rxt(85)+.500_r8*rxt(86))*y(31)
         mat(163) = rxt(88)*y(31)
         mat(42) = (rxt(85)+.500_r8*rxt(86))*y(2) + rxt(88)*y(10)


         mat(222) = rxt(87)*y(32)
         mat(28) = rxt(87)*y(2)




























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
         mat(  31) = mat(  31) + lmat(  31)
         mat(  33) = mat(  33) + lmat(  33)
         mat(  34) = lmat(  34)
         mat(  35) = lmat(  35)
         mat(  36) = lmat(  36)
         mat(  37) = lmat(  37)
         mat(  38) = lmat(  38)
         mat(  39) = lmat(  39)
         mat(  40) = lmat(  40)
         mat(  41) = lmat(  41)
         mat(  43) = mat(  43) + lmat(  43)
         mat(  48) = mat(  48) + lmat(  48)
         mat(  54) = mat(  54) + lmat(  54)
         mat(  60) = mat(  60) + lmat(  60)
         mat(  62) = lmat(  62)
         mat(  63) = mat(  63) + lmat(  63)
         mat(  64) = mat(  64) + lmat(  64)
         mat(  65) = lmat(  65)
         mat(  66) = mat(  66) + lmat(  66)
         mat(  67) = lmat(  67)
         mat(  68) = mat(  68) + lmat(  68)
         mat(  69) = mat(  69) + lmat(  69)
         mat(  70) = mat(  70) + lmat(  70)
         mat(  71) = lmat(  71)
         mat(  73) = mat(  73) + lmat(  73)
         mat(  74) = mat(  74) + lmat(  74)
         mat(  76) = mat(  76) + lmat(  76)
         mat(  77) = lmat(  77)
         mat(  78) = mat(  78) + lmat(  78)
         mat(  79) = mat(  79) + lmat(  79)
         mat(  80) = mat(  80) + lmat(  80)
         mat(  84) = mat(  84) + lmat(  84)
         mat(  93) = mat(  93) + lmat(  93)
         mat( 102) = mat( 102) + lmat( 102)
         mat( 112) = mat( 112) + lmat( 112)
         mat( 120) = mat( 120) + lmat( 120)
         mat( 123) = mat( 123) + lmat( 123)
         mat( 128) = lmat( 128)
         mat( 129) = mat( 129) + lmat( 129)
         mat( 133) = mat( 133) + lmat( 133)
         mat( 137) = lmat( 137)
         mat( 138) = mat( 138) + lmat( 138)
         mat( 142) = mat( 142) + lmat( 142)
         mat( 156) = mat( 156) + lmat( 156)
         mat( 169) = mat( 169) + lmat( 169)
         mat( 170) = mat( 170) + lmat( 170)
         mat( 171) = mat( 171) + lmat( 171)
         mat( 175) = lmat( 175)
         mat( 181) = lmat( 181)
         mat( 183) = mat( 183) + lmat( 183)
         mat( 187) = mat( 187) + lmat( 187)
         mat( 203) = mat( 203) + lmat( 203)
         mat( 219) = mat( 219) + lmat( 219)
         mat( 232) = mat( 232) + lmat( 232)
         mat( 235) = mat( 235) + lmat( 235)
         mat( 238) = mat( 238) + lmat( 238)
         mat( 239) = mat( 239) + lmat( 239)
         mat( 243) = mat( 243) + lmat( 243)
         mat( 244) = mat( 244) + lmat( 244)
         mat( 245) = mat( 245) + lmat( 245)
         mat( 258) = mat( 258) + lmat( 258)
         mat( 260) = mat( 260) + lmat( 260)
         mat( 261) = mat( 261) + lmat( 261)
         mat(  85) = 0._r8
         mat(  86) = 0._r8
         mat(  91) = 0._r8
         mat( 108) = 0._r8
         mat( 111) = 0._r8
         mat( 114) = 0._r8
         mat( 119) = 0._r8
         mat( 122) = 0._r8
         mat( 124) = 0._r8
         mat( 126) = 0._r8
         mat( 127) = 0._r8
         mat( 136) = 0._r8
         mat( 147) = 0._r8
         mat( 173) = 0._r8
         mat( 185) = 0._r8
         mat( 197) = 0._r8
         mat( 199) = 0._r8
         mat( 201) = 0._r8
         mat( 216) = 0._r8
         mat( 220) = 0._r8
         mat( 221) = 0._r8
         mat( 240) = 0._r8
         mat( 250) = 0._r8
         mat( 251) = 0._r8
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
         mat(  31) = mat(  31) - dti
         mat(  34) = mat(  34) - dti
         mat(  37) = mat(  37) - dti
         mat(  43) = mat(  43) - dti
         mat(  48) = mat(  48) - dti
         mat(  54) = mat(  54) - dti
         mat(  60) = mat(  60) - dti
         mat(  64) = mat(  64) - dti
         mat(  69) = mat(  69) - dti
         mat(  74) = mat(  74) - dti
         mat(  79) = mat(  79) - dti
         mat(  84) = mat(  84) - dti
         mat(  93) = mat(  93) - dti
         mat( 102) = mat( 102) - dti
         mat( 112) = mat( 112) - dti
         mat( 123) = mat( 123) - dti
         mat( 133) = mat( 133) - dti
         mat( 142) = mat( 142) - dti
         mat( 156) = mat( 156) - dti
         mat( 170) = mat( 170) - dti
         mat( 183) = mat( 183) - dti
         mat( 203) = mat( 203) - dti
         mat( 219) = mat( 219) - dti
         mat( 245) = mat( 245) - dti
         mat( 261) = mat( 261) - dti

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
