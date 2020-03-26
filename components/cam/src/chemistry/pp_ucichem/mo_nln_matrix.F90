













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


         mat(171) = -(rxt(20)*y(2) + rxt(22)*y(3) + rxt(28)*y(8) + rxt(30)*y(9) + rxt(53) &
                      *y(19) + rxt(54)*y(25) + rxt(80)*y(27))
         mat(230) = -rxt(20)*y(1)
         mat(190) = -rxt(22)*y(1)
         mat(145) = -rxt(28)*y(1)
         mat(242) = -rxt(30)*y(1)
         mat(39) = -rxt(53)*y(1)
         mat(45) = -rxt(54)*y(1)
         mat(114) = -rxt(80)*y(1)

         mat(230) = mat(230) + 2.000_r8*rxt(21)*y(2)

         mat(233) = -(rxt(20)*y(1) + (4._r8*rxt(21) + 4._r8*rxt(26)) * y(2) + rxt(23) &
                      *y(3) + rxt(27)*y(4) + rxt(31)*y(10) + rxt(33)*y(9) + rxt(43) &
                      *y(13) + (rxt(44) + rxt(45)) * y(12) + rxt(52)*y(19) + rxt(55) &
                      *y(25) + rxt(56)*y(5) + rxt(60)*y(7) + rxt(65)*y(22) + rxt(67) &
                      *y(23) + rxt(81)*y(27) + (rxt(84) + rxt(85)) * y(36) + rxt(86) &
                      *y(35))
         mat(174) = -rxt(20)*y(2)
         mat(193) = -rxt(23)*y(2)
         mat(27) = -rxt(27)*y(2)
         mat(159) = -rxt(31)*y(2)
         mat(245) = -rxt(33)*y(2)
         mat(59) = -rxt(43)*y(2)
         mat(54) = -(rxt(44) + rxt(45)) * y(2)
         mat(40) = -rxt(52)*y(2)
         mat(46) = -rxt(55)*y(2)
         mat(68) = -rxt(56)*y(2)
         mat(51) = -rxt(60)*y(2)
         mat(65) = -rxt(65)*y(2)
         mat(125) = -rxt(67)*y(2)
         mat(117) = -rxt(81)*y(2)
         mat(24) = -(rxt(84) + rxt(85)) * y(2)
         mat(20) = -rxt(86)*y(2)

         mat(174) = mat(174) + rxt(22)*y(3) + rxt(54)*y(25)
         mat(233) = mat(233) + .300_r8*rxt(60)*y(7)
         mat(193) = mat(193) + rxt(22)*y(1) + rxt(29)*y(8) + rxt(75)*y(29)
         mat(51) = mat(51) + .300_r8*rxt(60)*y(2)
         mat(148) = rxt(29)*y(3)
         mat(85) = rxt(75)*y(3)
         mat(46) = mat(46) + rxt(54)*y(1)

         mat(191) = -(rxt(22)*y(1) + rxt(23)*y(2) + (4._r8*rxt(24) + 4._r8*rxt(25) &
                      ) * y(3) + rxt(29)*y(8) + rxt(34)*y(9) + rxt(57)*y(6) + rxt(64) &
                      *y(21) + rxt(75)*y(29) + rxt(78)*y(26) + rxt(83)*y(28))
         mat(172) = -rxt(22)*y(3)
         mat(231) = -rxt(23)*y(3)
         mat(146) = -rxt(29)*y(3)
         mat(243) = -rxt(34)*y(3)
         mat(206) = -rxt(57)*y(3)
         mat(92) = -rxt(64)*y(3)
         mat(83) = -rxt(75)*y(3)
         mat(103) = -rxt(78)*y(3)
         mat(76) = -rxt(83)*y(3)

         mat(172) = mat(172) + rxt(20)*y(2)
         mat(231) = mat(231) + rxt(20)*y(1) + rxt(27)*y(4) + rxt(56)*y(5) + rxt(31) &
                      *y(10) + .500_r8*rxt(85)*y(36)
         mat(26) = rxt(27)*y(2)
         mat(67) = rxt(56)*y(2)
         mat(206) = mat(206) + 4.000_r8*rxt(59)*y(6) + rxt(58)*y(8) + 2.000_r8*rxt(63) &
                      *y(21) + rxt(70)*y(24) + rxt(76)*y(29) + 2.000_r8*rxt(79)*y(26)
         mat(146) = mat(146) + rxt(58)*y(6) + rxt(61)*y(21)
         mat(157) = rxt(31)*y(2)
         mat(92) = mat(92) + 2.000_r8*rxt(63)*y(6) + rxt(61)*y(8) + 4.000_r8*rxt(62) &
                      *y(21)
         mat(131) = rxt(70)*y(6)
         mat(83) = mat(83) + rxt(76)*y(6)
         mat(103) = mat(103) + 2.000_r8*rxt(79)*y(6)
         mat(23) = .500_r8*rxt(85)*y(2)

         mat(25) = -(rxt(27)*y(2))
         mat(213) = -rxt(27)*y(4)

         mat(213) = mat(213) + 2.000_r8*rxt(26)*y(2)
         mat(176) = (2.000_r8*rxt(24)+2.000_r8*rxt(25))*y(3)

         mat(66) = -(rxt(56)*y(2))
         mat(220) = -rxt(56)*y(5)

         mat(163) = rxt(53)*y(19) + rxt(54)*y(25) + rxt(80)*y(27)
         mat(220) = mat(220) + .300_r8*rxt(60)*y(7)
         mat(180) = rxt(75)*y(29)
         mat(196) = 4.000_r8*rxt(59)*y(6) + rxt(58)*y(8) + rxt(63)*y(21) + rxt(70) &
                      *y(24) + 2.000_r8*rxt(76)*y(29) + 2.000_r8*rxt(79)*y(26)
         mat(48) = .300_r8*rxt(60)*y(2)
         mat(135) = rxt(58)*y(6) + rxt(74)*y(29) + rxt(77)*y(26) + rxt(82)*y(28)
         mat(88) = rxt(63)*y(6)
         mat(128) = rxt(70)*y(6)
         mat(36) = rxt(53)*y(1)
         mat(79) = rxt(75)*y(3) + 2.000_r8*rxt(76)*y(6) + rxt(74)*y(8)
         mat(42) = rxt(54)*y(1)
         mat(97) = 2.000_r8*rxt(79)*y(6) + rxt(77)*y(8)
         mat(107) = rxt(80)*y(1)
         mat(70) = rxt(82)*y(8)

         mat(207) = -(rxt(57)*y(3) + rxt(58)*y(8) + 4._r8*rxt(59)*y(6) + rxt(63)*y(21) &
                      + rxt(76)*y(29) + rxt(79)*y(26))
         mat(192) = -rxt(57)*y(6)
         mat(147) = -rxt(58)*y(6)
         mat(93) = -rxt(63)*y(6)
         mat(84) = -rxt(76)*y(6)
         mat(104) = -rxt(79)*y(6)

         mat(173) = rxt(80)*y(27)
         mat(232) = .700_r8*rxt(60)*y(7) + .500_r8*rxt(67)*y(23)
         mat(50) = .700_r8*rxt(60)*y(2)
         mat(147) = mat(147) + rxt(69)*y(24)
         mat(124) = .500_r8*rxt(67)*y(2)
         mat(132) = rxt(69)*y(8) + 4.000_r8*rxt(71)*y(24)
         mat(116) = rxt(80)*y(1)

         mat(47) = -(rxt(60)*y(2))
         mat(216) = -rxt(60)*y(7)

         mat(177) = rxt(57)*y(6)
         mat(195) = rxt(57)*y(3)

         mat(143) = -(rxt(28)*y(1) + rxt(29)*y(3) + rxt(32)*y(10) + rxt(58)*y(6) &
                      + rxt(61)*y(21) + rxt(69)*y(24) + rxt(74)*y(29) + rxt(77)*y(26) &
                      + rxt(82)*y(28))
         mat(169) = -rxt(28)*y(8)
         mat(188) = -rxt(29)*y(8)
         mat(154) = -rxt(32)*y(8)
         mat(203) = -rxt(58)*y(8)
         mat(91) = -rxt(61)*y(8)
         mat(130) = -rxt(69)*y(8)
         mat(82) = -rxt(74)*y(8)
         mat(102) = -rxt(77)*y(8)
         mat(75) = -rxt(82)*y(8)

         mat(246) = -(rxt(30)*y(1) + rxt(33)*y(2) + rxt(34)*y(3) + rxt(35)*y(10) &
                      + rxt(36)*y(24))
         mat(175) = -rxt(30)*y(9)
         mat(234) = -rxt(33)*y(9)
         mat(194) = -rxt(34)*y(9)
         mat(160) = -rxt(35)*y(9)
         mat(134) = -rxt(36)*y(9)

         mat(175) = mat(175) + rxt(28)*y(8)
         mat(234) = mat(234) + rxt(31)*y(10) + rxt(43)*y(13)
         mat(194) = mat(194) + rxt(29)*y(8)
         mat(209) = rxt(58)*y(8)
         mat(149) = rxt(28)*y(1) + rxt(29)*y(3) + rxt(58)*y(6) + 2.000_r8*rxt(32) &
                      *y(10) + rxt(61)*y(21) + rxt(69)*y(24) + rxt(74)*y(29) + rxt(77) &
                      *y(26) + rxt(82)*y(28)
         mat(160) = mat(160) + rxt(31)*y(2) + 2.000_r8*rxt(32)*y(8)
         mat(60) = rxt(43)*y(2)
         mat(95) = rxt(61)*y(8)
         mat(134) = mat(134) + rxt(69)*y(8)
         mat(86) = rxt(74)*y(8)
         mat(106) = rxt(77)*y(8)
         mat(78) = rxt(82)*y(8)

         mat(155) = -(rxt(31)*y(2) + rxt(32)*y(8) + rxt(35)*y(9) + rxt(68)*y(23))
         mat(229) = -rxt(31)*y(10)
         mat(144) = -rxt(32)*y(10)
         mat(241) = -rxt(35)*y(10)
         mat(122) = -rxt(68)*y(10)

         mat(170) = rxt(30)*y(9)
         mat(229) = mat(229) + (rxt(44)+rxt(45))*y(12)
         mat(241) = mat(241) + rxt(30)*y(1)
         mat(53) = (rxt(44)+rxt(45))*y(2)


         mat(236) = rxt(35)*y(10)
         mat(150) = rxt(35)*y(9)

         mat(52) = -((rxt(44) + rxt(45)) * y(2))
         mat(217) = -(rxt(44) + rxt(45)) * y(12)

         mat(217) = mat(217) + rxt(33)*y(9)
         mat(237) = rxt(33)*y(2)
         mat(151) = rxt(68)*y(23)
         mat(119) = rxt(68)*y(10)

         mat(56) = -(rxt(43)*y(2))
         mat(218) = -rxt(43)*y(13)

         mat(178) = rxt(34)*y(9)
         mat(238) = rxt(34)*y(3)


         mat(235) = rxt(36)*y(24)
         mat(127) = rxt(36)*y(9)

         mat(89) = -(rxt(61)*y(8) + 4._r8*rxt(62)*y(21) + rxt(63)*y(6) + rxt(64)*y(3))
         mat(138) = -rxt(61)*y(21)
         mat(198) = -rxt(63)*y(21)
         mat(183) = -rxt(64)*y(21)

         mat(223) = rxt(65)*y(22)
         mat(62) = rxt(65)*y(2)

         mat(61) = -((rxt(65) + rxt(66)) * y(2))
         mat(219) = -(rxt(65) + rxt(66)) * y(22)

         mat(179) = rxt(64)*y(21) + rxt(78)*y(26) + rxt(83)*y(28)
         mat(87) = rxt(64)*y(3)
         mat(96) = rxt(78)*y(3)
         mat(69) = rxt(83)*y(3)

         mat(120) = -(rxt(67)*y(2) + rxt(68)*y(10))
         mat(226) = -rxt(67)*y(23)
         mat(152) = -rxt(68)*y(23)

         mat(167) = .500_r8*rxt(53)*y(19)
         mat(226) = mat(226) + rxt(66)*y(22)
         mat(186) = rxt(75)*y(29)
         mat(201) = rxt(63)*y(21)
         mat(141) = rxt(61)*y(21) + rxt(74)*y(29)
         mat(90) = rxt(63)*y(6) + rxt(61)*y(8) + 4.000_r8*rxt(62)*y(21)
         mat(63) = rxt(66)*y(2)
         mat(38) = .500_r8*rxt(53)*y(1)
         mat(81) = rxt(75)*y(3) + rxt(74)*y(8)

         mat(129) = -(rxt(36)*y(9) + rxt(69)*y(8) + rxt(70)*y(6) + 4._r8*rxt(71)*y(24))
         mat(239) = -rxt(36)*y(24)
         mat(142) = -rxt(69)*y(24)
         mat(202) = -rxt(70)*y(24)

         mat(168) = .500_r8*rxt(80)*y(27)
         mat(227) = .500_r8*rxt(67)*y(23)
         mat(142) = mat(142) + .500_r8*rxt(82)*y(28)
         mat(153) = rxt(68)*y(23)
         mat(121) = .500_r8*rxt(67)*y(2) + rxt(68)*y(10)
         mat(112) = .500_r8*rxt(80)*y(1)
         mat(74) = .500_r8*rxt(82)*y(8)

         mat(35) = -(rxt(52)*y(2) + rxt(53)*y(1))
         mat(214) = -rxt(52)*y(19)
         mat(161) = -rxt(53)*y(19)

         mat(80) = -(rxt(74)*y(8) + rxt(75)*y(3) + rxt(76)*y(6))
         mat(137) = -rxt(74)*y(29)
         mat(182) = -rxt(75)*y(29)
         mat(197) = -rxt(76)*y(29)

         mat(222) = rxt(52)*y(19)
         mat(37) = rxt(52)*y(2)


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


         mat(41) = -(rxt(54)*y(1) + rxt(55)*y(2))
         mat(162) = -rxt(54)*y(25)
         mat(215) = -rxt(55)*y(25)

         mat(99) = -(rxt(77)*y(8) + rxt(78)*y(3) + rxt(79)*y(6))
         mat(139) = -rxt(77)*y(26)
         mat(184) = -rxt(78)*y(26)
         mat(199) = -rxt(79)*y(26)

         mat(224) = rxt(55)*y(25)
         mat(43) = rxt(55)*y(2)

         mat(110) = -(rxt(80)*y(1) + rxt(81)*y(2))
         mat(166) = -rxt(80)*y(27)
         mat(225) = -rxt(81)*y(27)

         mat(166) = mat(166) + rxt(54)*y(25)
         mat(200) = rxt(79)*y(26)
         mat(140) = rxt(77)*y(26)
         mat(44) = rxt(54)*y(1)
         mat(100) = rxt(79)*y(6) + rxt(77)*y(8)

         mat(71) = -(rxt(82)*y(8) + rxt(83)*y(3))
         mat(136) = -rxt(82)*y(28)
         mat(181) = -rxt(83)*y(28)

         mat(221) = rxt(81)*y(27)
         mat(108) = rxt(81)*y(2)

         mat(22) = -((rxt(84) + rxt(85)) * y(2))
         mat(212) = -(rxt(84) + rxt(85)) * y(36)

         mat(19) = -(rxt(86)*y(2))
         mat(211) = -rxt(86)*y(35)

         mat(211) = mat(211) + (rxt(84)+.500_r8*rxt(85))*y(36)
         mat(21) = (rxt(84)+.500_r8*rxt(85))*y(2)


         mat(210) = rxt(86)*y(35)
         mat(18) = rxt(86)*y(2)


















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
         mat(  19) = mat(  19) + lmat(  19)
         mat(  22) = mat(  22) + lmat(  22)
         mat(  25) = mat(  25) + lmat(  25)
         mat(  27) = mat(  27) + lmat(  27)
         mat(  28) = lmat(  28)
         mat(  29) = lmat(  29)
         mat(  30) = lmat(  30)
         mat(  31) = lmat(  31)
         mat(  32) = lmat(  32)
         mat(  33) = lmat(  33)
         mat(  34) = lmat(  34)
         mat(  35) = mat(  35) + lmat(  35)
         mat(  41) = mat(  41) + lmat(  41)
         mat(  47) = mat(  47) + lmat(  47)
         mat(  48) = mat(  48) + lmat(  48)
         mat(  49) = lmat(  49)
         mat(  51) = mat(  51) + lmat(  51)
         mat(  52) = mat(  52) + lmat(  52)
         mat(  54) = mat(  54) + lmat(  54)
         mat(  55) = lmat(  55)
         mat(  56) = mat(  56) + lmat(  56)
         mat(  57) = lmat(  57)
         mat(  58) = lmat(  58)
         mat(  59) = mat(  59) + lmat(  59)
         mat(  60) = mat(  60) + lmat(  60)
         mat(  61) = mat(  61) + lmat(  61)
         mat(  63) = mat(  63) + lmat(  63)
         mat(  64) = lmat(  64)
         mat(  65) = mat(  65) + lmat(  65)
         mat(  66) = mat(  66) + lmat(  66)
         mat(  67) = mat(  67) + lmat(  67)
         mat(  71) = mat(  71) + lmat(  71)
         mat(  80) = mat(  80) + lmat(  80)
         mat(  89) = mat(  89) + lmat(  89)
         mat(  99) = mat(  99) + lmat(  99)
         mat( 107) = mat( 107) + lmat( 107)
         mat( 110) = mat( 110) + lmat( 110)
         mat( 115) = lmat( 115)
         mat( 116) = mat( 116) + lmat( 116)
         mat( 120) = mat( 120) + lmat( 120)
         mat( 123) = lmat( 123)
         mat( 124) = mat( 124) + lmat( 124)
         mat( 129) = mat( 129) + lmat( 129)
         mat( 143) = mat( 143) + lmat( 143)
         mat( 151) = mat( 151) + lmat( 151)
         mat( 154) = mat( 154) + lmat( 154)
         mat( 155) = mat( 155) + lmat( 155)
         mat( 156) = lmat( 156)
         mat( 160) = mat( 160) + lmat( 160)
         mat( 171) = mat( 171) + lmat( 171)
         mat( 172) = mat( 172) + lmat( 172)
         mat( 174) = mat( 174) + lmat( 174)
         mat( 191) = mat( 191) + lmat( 191)
         mat( 207) = mat( 207) + lmat( 207)
         mat( 220) = mat( 220) + lmat( 220)
         mat( 223) = mat( 223) + lmat( 223)
         mat( 226) = mat( 226) + lmat( 226)
         mat( 227) = mat( 227) + lmat( 227)
         mat( 231) = mat( 231) + lmat( 231)
         mat( 232) = mat( 232) + lmat( 232)
         mat( 233) = mat( 233) + lmat( 233)
         mat( 240) = lmat( 240)
         mat( 242) = mat( 242) + lmat( 242)
         mat( 246) = mat( 246) + lmat( 246)
         mat(  72) = 0._r8
         mat(  73) = 0._r8
         mat(  77) = 0._r8
         mat(  94) = 0._r8
         mat(  98) = 0._r8
         mat( 101) = 0._r8
         mat( 105) = 0._r8
         mat( 109) = 0._r8
         mat( 111) = 0._r8
         mat( 113) = 0._r8
         mat( 118) = 0._r8
         mat( 126) = 0._r8
         mat( 133) = 0._r8
         mat( 158) = 0._r8
         mat( 164) = 0._r8
         mat( 165) = 0._r8
         mat( 185) = 0._r8
         mat( 187) = 0._r8
         mat( 189) = 0._r8
         mat( 204) = 0._r8
         mat( 205) = 0._r8
         mat( 208) = 0._r8
         mat( 228) = 0._r8
         mat( 244) = 0._r8
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
         mat(  19) = mat(  19) - dti
         mat(  22) = mat(  22) - dti
         mat(  25) = mat(  25) - dti
         mat(  28) = mat(  28) - dti
         mat(  31) = mat(  31) - dti
         mat(  35) = mat(  35) - dti
         mat(  41) = mat(  41) - dti
         mat(  47) = mat(  47) - dti
         mat(  52) = mat(  52) - dti
         mat(  56) = mat(  56) - dti
         mat(  61) = mat(  61) - dti
         mat(  66) = mat(  66) - dti
         mat(  71) = mat(  71) - dti
         mat(  80) = mat(  80) - dti
         mat(  89) = mat(  89) - dti
         mat(  99) = mat(  99) - dti
         mat( 110) = mat( 110) - dti
         mat( 120) = mat( 120) - dti
         mat( 129) = mat( 129) - dti
         mat( 143) = mat( 143) - dti
         mat( 155) = mat( 155) - dti
         mat( 171) = mat( 171) - dti
         mat( 191) = mat( 191) - dti
         mat( 207) = mat( 207) - dti
         mat( 233) = mat( 233) - dti
         mat( 246) = mat( 246) - dti

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
