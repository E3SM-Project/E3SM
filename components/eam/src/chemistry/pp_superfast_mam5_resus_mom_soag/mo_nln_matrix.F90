











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


         mat(212) = -(rxt(33)*y(18) + rxt(34)*y(25) + rxt(37)*y(2) + rxt(39)*y(3) &
                      + rxt(45)*y(8) + rxt(47)*y(9) + rxt(83)*y(29))
         mat(47) = -rxt(33)*y(1)
         mat(53) = -rxt(34)*y(1)
         mat(234) = -rxt(37)*y(1)
         mat(159) = -rxt(39)*y(1)
         mat(174) = -rxt(45)*y(1)
         mat(186) = -rxt(47)*y(1)
         mat(124) = -rxt(83)*y(1)

         mat(234) = mat(234) + 2.000_r8*rxt(38)*y(2)

         mat(235) = -(rxt(32)*y(18) + rxt(35)*y(25) + rxt(36)*y(5) + rxt(37)*y(1) &
                      + (4._r8*rxt(38) + 4._r8*rxt(43)) * y(2) + rxt(40)*y(3) + rxt(44) &
                      *y(4) + rxt(48)*y(10) + rxt(50)*y(13) + (rxt(51) + rxt(52) &
                      ) * y(12) + rxt(53)*y(9) + rxt(63)*y(7) + rxt(68)*y(22) + rxt(70) &
                      *y(23) + rxt(84)*y(29))
         mat(48) = -rxt(32)*y(2)
         mat(54) = -rxt(35)*y(2)
         mat(76) = -rxt(36)*y(2)
         mat(213) = -rxt(37)*y(2)
         mat(160) = -rxt(40)*y(2)
         mat(37) = -rxt(44)*y(2)
         mat(198) = -rxt(48)*y(2)
         mat(68) = -rxt(50)*y(2)
         mat(63) = -(rxt(51) + rxt(52)) * y(2)
         mat(187) = -rxt(53)*y(2)
         mat(58) = -rxt(63)*y(2)
         mat(73) = -rxt(68)*y(2)
         mat(133) = -rxt(70)*y(2)
         mat(125) = -rxt(84)*y(2)

         mat(213) = mat(213) + rxt(39)*y(3) + rxt(34)*y(25)
         mat(235) = mat(235) + .300_r8*rxt(63)*y(7)
         mat(160) = mat(160) + rxt(39)*y(1) + rxt(46)*y(8) + rxt(78)*y(19)
         mat(58) = mat(58) + .300_r8*rxt(63)*y(2)
         mat(175) = rxt(46)*y(3)
         mat(93) = rxt(78)*y(3)
         mat(54) = mat(54) + rxt(34)*y(1)

         mat(155) = -(rxt(39)*y(1) + rxt(40)*y(2) + (4._r8*rxt(41) + 4._r8*rxt(42) &
                      ) * y(3) + rxt(46)*y(8) + rxt(54)*y(9) + rxt(60)*y(6) + rxt(67) &
                      *y(21) + rxt(78)*y(19) + rxt(81)*y(27) + rxt(86)*y(30))
         mat(208) = -rxt(39)*y(3)
         mat(230) = -rxt(40)*y(3)
         mat(170) = -rxt(46)*y(3)
         mat(182) = -rxt(54)*y(3)
         mat(245) = -rxt(60)*y(3)
         mat(99) = -rxt(67)*y(3)
         mat(90) = -rxt(78)*y(3)
         mat(110) = -rxt(81)*y(3)
         mat(83) = -rxt(86)*y(3)

         mat(208) = mat(208) + rxt(37)*y(2)
         mat(230) = mat(230) + rxt(37)*y(1) + rxt(44)*y(4) + rxt(36)*y(5) + rxt(48) &
                      *y(10)
         mat(36) = rxt(44)*y(2)
         mat(75) = rxt(36)*y(2)
         mat(245) = mat(245) + 4.000_r8*rxt(62)*y(6) + rxt(61)*y(8) + 2.000_r8*rxt(66) &
                      *y(21) + rxt(73)*y(24) + rxt(79)*y(19) + 2.000_r8*rxt(82)*y(27)
         mat(170) = mat(170) + rxt(61)*y(6) + rxt(64)*y(21)
         mat(193) = rxt(48)*y(2)
         mat(99) = mat(99) + 2.000_r8*rxt(66)*y(6) + rxt(64)*y(8) + 4.000_r8*rxt(65) &
                      *y(21)
         mat(138) = rxt(73)*y(6)
         mat(90) = mat(90) + rxt(79)*y(6)
         mat(110) = mat(110) + 2.000_r8*rxt(82)*y(6)

         mat(35) = -(rxt(44)*y(2))
         mat(215) = -rxt(44)*y(4)

         mat(215) = mat(215) + 2.000_r8*rxt(43)*y(2)
         mat(143) = (2.000_r8*rxt(41)+2.000_r8*rxt(42))*y(3)

         mat(74) = -(rxt(36)*y(2))
         mat(222) = -rxt(36)*y(5)

         mat(202) = rxt(33)*y(18) + rxt(34)*y(25) + rxt(83)*y(29)
         mat(222) = mat(222) + .300_r8*rxt(63)*y(7)
         mat(147) = rxt(78)*y(19)
         mat(238) = 4.000_r8*rxt(62)*y(6) + rxt(61)*y(8) + rxt(66)*y(21) + rxt(73) &
                      *y(24) + 2.000_r8*rxt(79)*y(19) + 2.000_r8*rxt(82)*y(27)
         mat(56) = .300_r8*rxt(63)*y(2)
         mat(162) = rxt(61)*y(6) + rxt(77)*y(19) + rxt(80)*y(27) + rxt(85)*y(30)
         mat(96) = rxt(66)*y(6)
         mat(136) = rxt(73)*y(6)
         mat(44) = rxt(33)*y(1)
         mat(87) = rxt(78)*y(3) + 2.000_r8*rxt(79)*y(6) + rxt(77)*y(8)
         mat(50) = rxt(34)*y(1)
         mat(105) = 2.000_r8*rxt(82)*y(6) + rxt(80)*y(8)
         mat(115) = rxt(83)*y(1)
         mat(78) = rxt(85)*y(8)

         mat(251) = -(rxt(60)*y(3) + rxt(61)*y(8) + 4._r8*rxt(62)*y(6) + rxt(66)*y(21) &
                      + rxt(79)*y(19) + rxt(82)*y(27))
         mat(161) = -rxt(60)*y(6)
         mat(176) = -rxt(61)*y(6)
         mat(103) = -rxt(66)*y(6)
         mat(94) = -rxt(79)*y(6)
         mat(114) = -rxt(82)*y(6)

         mat(214) = rxt(83)*y(29)
         mat(236) = .700_r8*rxt(63)*y(7) + .500_r8*rxt(70)*y(23)
         mat(59) = .700_r8*rxt(63)*y(2)
         mat(176) = mat(176) + rxt(72)*y(24)
         mat(134) = .500_r8*rxt(70)*y(2)
         mat(142) = rxt(72)*y(8) + 4.000_r8*rxt(74)*y(24)
         mat(126) = rxt(83)*y(1)

         mat(55) = -(rxt(63)*y(2))
         mat(218) = -rxt(63)*y(7)

         mat(144) = rxt(60)*y(6)
         mat(237) = rxt(60)*y(3)

         mat(171) = -(rxt(45)*y(1) + rxt(46)*y(3) + rxt(49)*y(10) + rxt(61)*y(6) &
                      + rxt(64)*y(21) + rxt(72)*y(24) + rxt(77)*y(19) + rxt(80)*y(27) &
                      + rxt(85)*y(30))
         mat(209) = -rxt(45)*y(8)
         mat(156) = -rxt(46)*y(8)
         mat(194) = -rxt(49)*y(8)
         mat(246) = -rxt(61)*y(8)
         mat(100) = -rxt(64)*y(8)
         mat(139) = -rxt(72)*y(8)
         mat(91) = -rxt(77)*y(8)
         mat(111) = -rxt(80)*y(8)
         mat(84) = -rxt(85)*y(8)

         mat(184) = -(rxt(47)*y(1) + rxt(53)*y(2) + rxt(54)*y(3) + rxt(55)*y(10) &
                      + rxt(56)*y(24))
         mat(210) = -rxt(47)*y(9)
         mat(232) = -rxt(53)*y(9)
         mat(157) = -rxt(54)*y(9)
         mat(195) = -rxt(55)*y(9)
         mat(140) = -rxt(56)*y(9)

         mat(210) = mat(210) + rxt(45)*y(8)
         mat(232) = mat(232) + rxt(48)*y(10) + rxt(50)*y(13)
         mat(157) = mat(157) + rxt(46)*y(8)
         mat(247) = rxt(61)*y(8)
         mat(172) = rxt(45)*y(1) + rxt(46)*y(3) + rxt(61)*y(6) + 2.000_r8*rxt(49) &
                      *y(10) + rxt(64)*y(21) + rxt(72)*y(24) + rxt(77)*y(19) + rxt(80) &
                      *y(27) + rxt(85)*y(30)
         mat(195) = mat(195) + rxt(48)*y(2) + 2.000_r8*rxt(49)*y(8)
         mat(66) = rxt(50)*y(2)
         mat(101) = rxt(64)*y(8)
         mat(140) = mat(140) + rxt(72)*y(8)
         mat(92) = rxt(77)*y(8)
         mat(112) = rxt(80)*y(8)
         mat(85) = rxt(85)*y(8)

         mat(196) = -(rxt(48)*y(2) + rxt(49)*y(8) + rxt(55)*y(9) + rxt(71)*y(23))
         mat(233) = -rxt(48)*y(10)
         mat(173) = -rxt(49)*y(10)
         mat(185) = -rxt(55)*y(10)
         mat(132) = -rxt(71)*y(10)

         mat(211) = rxt(47)*y(9)
         mat(233) = mat(233) + (rxt(51)+rxt(52))*y(12)
         mat(185) = mat(185) + rxt(47)*y(1)
         mat(62) = (rxt(51)+rxt(52))*y(2)


         mat(178) = rxt(55)*y(10)
         mat(189) = rxt(55)*y(9)

         mat(60) = -((rxt(51) + rxt(52)) * y(2))
         mat(219) = -(rxt(51) + rxt(52)) * y(12)

         mat(219) = mat(219) + rxt(53)*y(9)
         mat(179) = rxt(53)*y(2)
         mat(190) = rxt(71)*y(23)
         mat(127) = rxt(71)*y(10)

         mat(64) = -(rxt(50)*y(2))
         mat(220) = -rxt(50)*y(13)

         mat(145) = rxt(54)*y(9)
         mat(180) = rxt(54)*y(3)


         mat(177) = rxt(56)*y(24)
         mat(135) = rxt(56)*y(9)

         mat(97) = -(rxt(64)*y(8) + 4._r8*rxt(65)*y(21) + rxt(66)*y(6) + rxt(67)*y(3))
         mat(165) = -rxt(64)*y(21)
         mat(240) = -rxt(66)*y(21)
         mat(150) = -rxt(67)*y(21)

         mat(225) = rxt(68)*y(22)
         mat(70) = rxt(68)*y(2)

         mat(69) = -((rxt(68) + rxt(69)) * y(2))
         mat(221) = -(rxt(68) + rxt(69)) * y(22)

         mat(146) = rxt(67)*y(21) + rxt(81)*y(27) + rxt(86)*y(30)
         mat(95) = rxt(67)*y(3)
         mat(104) = rxt(81)*y(3)
         mat(77) = rxt(86)*y(3)

         mat(128) = -(rxt(70)*y(2) + rxt(71)*y(10))
         mat(228) = -rxt(70)*y(23)
         mat(191) = -rxt(71)*y(23)

         mat(206) = .500_r8*rxt(33)*y(18)
         mat(228) = mat(228) + rxt(69)*y(22)
         mat(153) = rxt(78)*y(19)
         mat(243) = rxt(66)*y(21)
         mat(168) = rxt(64)*y(21) + rxt(77)*y(19)
         mat(98) = rxt(66)*y(6) + rxt(64)*y(8) + 4.000_r8*rxt(65)*y(21)
         mat(71) = rxt(69)*y(2)
         mat(46) = .500_r8*rxt(33)*y(1)
         mat(89) = rxt(78)*y(3) + rxt(77)*y(8)

         mat(137) = -(rxt(56)*y(9) + rxt(72)*y(8) + rxt(73)*y(6) + 4._r8*rxt(74)*y(24))
         mat(181) = -rxt(56)*y(24)
         mat(169) = -rxt(72)*y(24)
         mat(244) = -rxt(73)*y(24)

         mat(207) = .500_r8*rxt(83)*y(29)
         mat(229) = .500_r8*rxt(70)*y(23)
         mat(169) = mat(169) + .500_r8*rxt(85)*y(30)
         mat(192) = rxt(71)*y(23)
         mat(129) = .500_r8*rxt(70)*y(2) + rxt(71)*y(10)
         mat(120) = .500_r8*rxt(83)*y(1)
         mat(82) = .500_r8*rxt(85)*y(8)

         mat(43) = -(rxt(32)*y(2) + rxt(33)*y(1))
         mat(216) = -rxt(32)*y(18)
         mat(200) = -rxt(33)*y(18)

         mat(88) = -(rxt(77)*y(8) + rxt(78)*y(3) + rxt(79)*y(6))
         mat(164) = -rxt(77)*y(19)
         mat(149) = -rxt(78)*y(19)
         mat(239) = -rxt(79)*y(19)

         mat(224) = rxt(32)*y(18)
         mat(45) = rxt(32)*y(2)


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


         mat(49) = -(rxt(34)*y(1) + rxt(35)*y(2))
         mat(201) = -rxt(34)*y(25)
         mat(217) = -rxt(35)*y(25)


         mat(107) = -(rxt(80)*y(8) + rxt(81)*y(3) + rxt(82)*y(6))
         mat(166) = -rxt(80)*y(27)
         mat(151) = -rxt(81)*y(27)
         mat(241) = -rxt(82)*y(27)

         mat(226) = rxt(35)*y(25)
         mat(51) = rxt(35)*y(2)


         mat(118) = -(rxt(83)*y(1) + rxt(84)*y(2))
         mat(205) = -rxt(83)*y(29)
         mat(227) = -rxt(84)*y(29)

         mat(205) = mat(205) + rxt(34)*y(25)
         mat(242) = rxt(82)*y(27)
         mat(167) = rxt(80)*y(27)
         mat(52) = rxt(34)*y(1)
         mat(108) = rxt(82)*y(6) + rxt(80)*y(8)

         mat(79) = -(rxt(85)*y(8) + rxt(86)*y(3))
         mat(163) = -rxt(85)*y(30)
         mat(148) = -rxt(86)*y(30)

         mat(223) = rxt(84)*y(29)
         mat(116) = rxt(84)*y(2)













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
         mat(  28) = lmat(  28)
         mat(  29) = lmat(  29)
         mat(  30) = lmat(  30)
         mat(  31) = lmat(  31)
         mat(  32) = lmat(  32)
         mat(  33) = lmat(  33)
         mat(  34) = lmat(  34)
         mat(  35) = mat(  35) + lmat(  35)
         mat(  37) = mat(  37) + lmat(  37)
         mat(  38) = lmat(  38)
         mat(  39) = lmat(  39)
         mat(  40) = lmat(  40)
         mat(  41) = lmat(  41)
         mat(  42) = lmat(  42)
         mat(  43) = mat(  43) + lmat(  43)
         mat(  49) = mat(  49) + lmat(  49)
         mat(  55) = mat(  55) + lmat(  55)
         mat(  56) = mat(  56) + lmat(  56)
         mat(  57) = lmat(  57)
         mat(  58) = mat(  58) + lmat(  58)
         mat(  60) = mat(  60) + lmat(  60)
         mat(  61) = lmat(  61)
         mat(  63) = mat(  63) + lmat(  63)
         mat(  64) = mat(  64) + lmat(  64)
         mat(  65) = lmat(  65)
         mat(  66) = mat(  66) + lmat(  66)
         mat(  67) = lmat(  67)
         mat(  68) = mat(  68) + lmat(  68)
         mat(  69) = mat(  69) + lmat(  69)
         mat(  71) = mat(  71) + lmat(  71)
         mat(  72) = lmat(  72)
         mat(  73) = mat(  73) + lmat(  73)
         mat(  74) = mat(  74) + lmat(  74)
         mat(  75) = mat(  75) + lmat(  75)
         mat(  79) = mat(  79) + lmat(  79)
         mat(  88) = mat(  88) + lmat(  88)
         mat(  97) = mat(  97) + lmat(  97)
         mat( 107) = mat( 107) + lmat( 107)
         mat( 115) = mat( 115) + lmat( 115)
         mat( 118) = mat( 118) + lmat( 118)
         mat( 121) = lmat( 121)
         mat( 126) = mat( 126) + lmat( 126)
         mat( 128) = mat( 128) + lmat( 128)
         mat( 130) = lmat( 130)
         mat( 134) = mat( 134) + lmat( 134)
         mat( 137) = mat( 137) + lmat( 137)
         mat( 155) = mat( 155) + lmat( 155)
         mat( 171) = mat( 171) + lmat( 171)
         mat( 183) = lmat( 183)
         mat( 184) = mat( 184) + lmat( 184)
         mat( 186) = mat( 186) + lmat( 186)
         mat( 190) = mat( 190) + lmat( 190)
         mat( 194) = mat( 194) + lmat( 194)
         mat( 195) = mat( 195) + lmat( 195)
         mat( 196) = mat( 196) + lmat( 196)
         mat( 197) = lmat( 197)
         mat( 208) = mat( 208) + lmat( 208)
         mat( 212) = mat( 212) + lmat( 212)
         mat( 213) = mat( 213) + lmat( 213)
         mat( 222) = mat( 222) + lmat( 222)
         mat( 225) = mat( 225) + lmat( 225)
         mat( 228) = mat( 228) + lmat( 228)
         mat( 229) = mat( 229) + lmat( 229)
         mat( 230) = mat( 230) + lmat( 230)
         mat( 235) = mat( 235) + lmat( 235)
         mat( 236) = mat( 236) + lmat( 236)
         mat( 251) = mat( 251) + lmat( 251)
         mat(  80) = 0._r8
         mat(  81) = 0._r8
         mat(  86) = 0._r8
         mat( 102) = 0._r8
         mat( 106) = 0._r8
         mat( 109) = 0._r8
         mat( 113) = 0._r8
         mat( 117) = 0._r8
         mat( 119) = 0._r8
         mat( 122) = 0._r8
         mat( 123) = 0._r8
         mat( 131) = 0._r8
         mat( 141) = 0._r8
         mat( 152) = 0._r8
         mat( 154) = 0._r8
         mat( 158) = 0._r8
         mat( 188) = 0._r8
         mat( 199) = 0._r8
         mat( 203) = 0._r8
         mat( 204) = 0._r8
         mat( 231) = 0._r8
         mat( 248) = 0._r8
         mat( 249) = 0._r8
         mat( 250) = 0._r8
         mat(   1) = mat(   1) - dti
         mat(   5) = mat(   5) - dti
         mat(  10) = mat(  10) - dti
         mat(  11) = mat(  11) - dti
         mat(  12) = mat(  12) - dti
         mat(  13) = mat(  13) - dti
         mat(  15) = mat(  15) - dti
         mat(  17) = mat(  17) - dti
         mat(  20) = mat(  20) - dti
         mat(  22) = mat(  22) - dti
         mat(  25) = mat(  25) - dti
         mat(  28) = mat(  28) - dti
         mat(  31) = mat(  31) - dti
         mat(  32) = mat(  32) - dti
         mat(  35) = mat(  35) - dti
         mat(  38) = mat(  38) - dti
         mat(  43) = mat(  43) - dti
         mat(  49) = mat(  49) - dti
         mat(  55) = mat(  55) - dti
         mat(  60) = mat(  60) - dti
         mat(  64) = mat(  64) - dti
         mat(  69) = mat(  69) - dti
         mat(  74) = mat(  74) - dti
         mat(  79) = mat(  79) - dti
         mat(  88) = mat(  88) - dti
         mat(  97) = mat(  97) - dti
         mat( 107) = mat( 107) - dti
         mat( 118) = mat( 118) - dti
         mat( 128) = mat( 128) - dti
         mat( 137) = mat( 137) - dti
         mat( 155) = mat( 155) - dti
         mat( 171) = mat( 171) - dti
         mat( 184) = mat( 184) - dti
         mat( 196) = mat( 196) - dti
         mat( 212) = mat( 212) - dti
         mat( 235) = mat( 235) - dti
         mat( 251) = mat( 251) - dti

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
