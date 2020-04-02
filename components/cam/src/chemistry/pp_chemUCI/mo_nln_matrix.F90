













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


         mat(132) = -(rxt(27)*y(19) + rxt(28)*y(26) + rxt(31)*y(2) + rxt(33)*y(3) &
                      + rxt(39)*y(8) + rxt(41)*y(9) + rxt(76)*y(28))
         mat(15) = -rxt(27)*y(1)
         mat(21) = -rxt(28)*y(1)
         mat(180) = -rxt(31)*y(1)
         mat(214) = -rxt(33)*y(1)
         mat(195) = -rxt(39)*y(1)
         mat(117) = -rxt(41)*y(1)
         mat(90) = -rxt(76)*y(1)

         mat(180) = mat(180) + 2.000_r8*rxt(32)*y(2)

         mat(183) = -(rxt(26)*y(19) + rxt(29)*y(26) + rxt(30)*y(5) + rxt(31)*y(1) &
                      + (4._r8*rxt(32) + 4._r8*rxt(37)) * y(2) + rxt(34)*y(3) + rxt(38) &
                      *y(4) + rxt(42)*y(10) + rxt(44)*y(13) + rxt(45)*y(12) + rxt(46) &
                      *y(9) + rxt(56)*y(7) + rxt(61)*y(23) + rxt(63)*y(24) + rxt(77) &
                      *y(28))
         mat(16) = -rxt(26)*y(2)
         mat(22) = -rxt(29)*y(2)
         mat(43) = -rxt(30)*y(2)
         mat(135) = -rxt(31)*y(2)
         mat(217) = -rxt(34)*y(2)
         mat(2) = -rxt(38)*y(2)
         mat(161) = -rxt(42)*y(2)
         mat(35) = -rxt(44)*y(2)
         mat(31) = -rxt(45)*y(2)
         mat(120) = -rxt(46)*y(2)
         mat(26) = -rxt(56)*y(2)
         mat(40) = -rxt(61)*y(2)
         mat(101) = -rxt(63)*y(2)
         mat(92) = -rxt(77)*y(2)

         mat(135) = mat(135) + rxt(33)*y(3) + rxt(28)*y(26)
         mat(183) = mat(183) + .300_r8*rxt(56)*y(7)
         mat(217) = mat(217) + rxt(33)*y(1) + rxt(40)*y(8) + rxt(71)*y(20)
         mat(26) = mat(26) + .300_r8*rxt(56)*y(2)
         mat(198) = rxt(40)*y(3)
         mat(60) = rxt(71)*y(3)
         mat(22) = mat(22) + rxt(28)*y(1)

         mat(219) = -(rxt(33)*y(1) + rxt(34)*y(2) + (4._r8*rxt(35) + 4._r8*rxt(36) &
                      ) * y(3) + rxt(40)*y(8) + rxt(47)*y(9) + rxt(53)*y(6) + rxt(60) &
                      *y(22) + rxt(71)*y(20) + rxt(74)*y(27) + rxt(79)*y(29))
         mat(137) = -rxt(33)*y(3)
         mat(185) = -rxt(34)*y(3)
         mat(200) = -rxt(40)*y(3)
         mat(122) = -rxt(47)*y(3)
         mat(152) = -rxt(53)*y(3)
         mat(71) = -rxt(60)*y(3)
         mat(62) = -rxt(71)*y(3)
         mat(82) = -rxt(74)*y(3)
         mat(54) = -rxt(79)*y(3)

         mat(137) = mat(137) + rxt(31)*y(2)
         mat(185) = mat(185) + rxt(31)*y(1) + rxt(38)*y(4) + rxt(30)*y(5) + rxt(42) &
                      *y(10)
         mat(3) = rxt(38)*y(2)
         mat(44) = rxt(30)*y(2)
         mat(152) = mat(152) + 4.000_r8*rxt(55)*y(6) + rxt(54)*y(8) + 2.000_r8*rxt(59) &
                      *y(22) + rxt(66)*y(25) + rxt(72)*y(20) + 2.000_r8*rxt(75)*y(27)
         mat(200) = mat(200) + rxt(54)*y(6) + rxt(57)*y(22)
         mat(163) = rxt(42)*y(2)
         mat(71) = mat(71) + 2.000_r8*rxt(59)*y(6) + rxt(57)*y(8) + 4.000_r8*rxt(58) &
                      *y(22)
         mat(110) = rxt(66)*y(6)
         mat(62) = mat(62) + rxt(72)*y(6)
         mat(82) = mat(82) + 2.000_r8*rxt(75)*y(6)

         mat(1) = -(rxt(38)*y(2))
         mat(164) = -rxt(38)*y(4)

         mat(164) = mat(164) + 2.000_r8*rxt(37)*y(2)
         mat(201) = (2.000_r8*rxt(35)+2.000_r8*rxt(36))*y(3)

         mat(42) = -(rxt(30)*y(2))
         mat(171) = -rxt(30)*y(5)

         mat(125) = rxt(27)*y(19) + rxt(28)*y(26) + rxt(76)*y(28)
         mat(171) = mat(171) + .300_r8*rxt(56)*y(7)
         mat(205) = rxt(71)*y(20)
         mat(139) = 4.000_r8*rxt(55)*y(6) + rxt(54)*y(8) + rxt(59)*y(22) + rxt(66) &
                      *y(25) + 2.000_r8*rxt(72)*y(20) + 2.000_r8*rxt(75)*y(27)
         mat(24) = .300_r8*rxt(56)*y(2)
         mat(186) = rxt(54)*y(6) + rxt(70)*y(20) + rxt(73)*y(27) + rxt(78)*y(29)
         mat(64) = rxt(59)*y(6)
         mat(104) = rxt(66)*y(6)
         mat(12) = rxt(27)*y(1)
         mat(55) = rxt(71)*y(3) + 2.000_r8*rxt(72)*y(6) + rxt(70)*y(8)
         mat(18) = rxt(28)*y(1)
         mat(73) = 2.000_r8*rxt(75)*y(6) + rxt(73)*y(8)
         mat(83) = rxt(76)*y(1)
         mat(46) = rxt(78)*y(8)

         mat(148) = -(rxt(53)*y(3) + rxt(54)*y(8) + 4._r8*rxt(55)*y(6) + rxt(59)*y(22) &
                      + rxt(72)*y(20) + rxt(75)*y(27))
         mat(215) = -rxt(53)*y(6)
         mat(196) = -rxt(54)*y(6)
         mat(68) = -rxt(59)*y(6)
         mat(59) = -rxt(72)*y(6)
         mat(79) = -rxt(75)*y(6)

         mat(133) = rxt(76)*y(28)
         mat(181) = .700_r8*rxt(56)*y(7) + .500_r8*rxt(63)*y(24)
         mat(25) = .700_r8*rxt(56)*y(2)
         mat(196) = mat(196) + rxt(65)*y(25)
         mat(99) = .500_r8*rxt(63)*y(2)
         mat(107) = rxt(65)*y(8) + 4.000_r8*rxt(67)*y(25)
         mat(91) = rxt(76)*y(1)

         mat(23) = -(rxt(56)*y(2))
         mat(167) = -rxt(56)*y(7)

         mat(202) = rxt(53)*y(6)
         mat(138) = rxt(53)*y(3)

         mat(199) = -(rxt(39)*y(1) + rxt(40)*y(3) + rxt(43)*y(10) + rxt(54)*y(6) &
                      + rxt(57)*y(22) + rxt(65)*y(25) + rxt(70)*y(20) + rxt(73)*y(27) &
                      + rxt(78)*y(29))
         mat(136) = -rxt(39)*y(8)
         mat(218) = -rxt(40)*y(8)
         mat(162) = -rxt(43)*y(8)
         mat(151) = -rxt(54)*y(8)
         mat(70) = -rxt(57)*y(8)
         mat(109) = -rxt(65)*y(8)
         mat(61) = -rxt(70)*y(8)
         mat(81) = -rxt(73)*y(8)
         mat(53) = -rxt(78)*y(8)

         mat(116) = -(rxt(41)*y(1) + rxt(46)*y(2) + rxt(47)*y(3) + rxt(48)*y(10) &
                      + rxt(49)*y(25))
         mat(131) = -rxt(41)*y(9)
         mat(179) = -rxt(46)*y(9)
         mat(213) = -rxt(47)*y(9)
         mat(157) = -rxt(48)*y(9)
         mat(106) = -rxt(49)*y(9)

         mat(131) = mat(131) + rxt(39)*y(8)
         mat(179) = mat(179) + rxt(42)*y(10) + rxt(44)*y(13)
         mat(213) = mat(213) + rxt(40)*y(8)
         mat(146) = rxt(54)*y(8)
         mat(194) = rxt(39)*y(1) + rxt(40)*y(3) + rxt(54)*y(6) + 2.000_r8*rxt(43) &
                      *y(10) + rxt(57)*y(22) + rxt(65)*y(25) + rxt(70)*y(20) + rxt(73) &
                      *y(27) + rxt(78)*y(29)
         mat(157) = mat(157) + rxt(42)*y(2) + 2.000_r8*rxt(43)*y(8)
         mat(33) = rxt(44)*y(2)
         mat(67) = rxt(57)*y(8)
         mat(106) = mat(106) + rxt(65)*y(8)
         mat(58) = rxt(70)*y(8)
         mat(78) = rxt(73)*y(8)
         mat(51) = rxt(78)*y(8)

         mat(160) = -(rxt(42)*y(2) + rxt(43)*y(8) + rxt(48)*y(9) + rxt(64)*y(24))
         mat(182) = -rxt(42)*y(10)
         mat(197) = -rxt(43)*y(10)
         mat(119) = -rxt(48)*y(10)
         mat(100) = -rxt(64)*y(10)

         mat(134) = rxt(41)*y(9)
         mat(182) = mat(182) + rxt(45)*y(12)
         mat(119) = mat(119) + rxt(41)*y(1)
         mat(30) = rxt(45)*y(2)


         mat(112) = rxt(48)*y(10)
         mat(153) = rxt(48)*y(9)

         mat(28) = -(rxt(45)*y(2))
         mat(168) = -rxt(45)*y(12)

         mat(168) = mat(168) + rxt(46)*y(9)
         mat(113) = rxt(46)*y(2)
         mat(154) = rxt(64)*y(24)
         mat(95) = rxt(64)*y(10)

         mat(32) = -(rxt(44)*y(2))
         mat(169) = -rxt(44)*y(13)

         mat(203) = rxt(47)*y(9)
         mat(114) = rxt(47)*y(3)


         mat(111) = rxt(49)*y(25)
         mat(103) = rxt(49)*y(9)

         mat(65) = -(rxt(57)*y(8) + 4._r8*rxt(58)*y(22) + rxt(59)*y(6) + rxt(60)*y(3))
         mat(189) = -rxt(57)*y(22)
         mat(141) = -rxt(59)*y(22)
         mat(208) = -rxt(60)*y(22)

         mat(174) = rxt(61)*y(23)
         mat(38) = rxt(61)*y(2)

         mat(37) = -((rxt(61) + rxt(62)) * y(2))
         mat(170) = -(rxt(61) + rxt(62)) * y(23)

         mat(204) = rxt(60)*y(22) + rxt(74)*y(27) + rxt(79)*y(29)
         mat(63) = rxt(60)*y(3)
         mat(72) = rxt(74)*y(3)
         mat(45) = rxt(79)*y(3)

         mat(96) = -(rxt(63)*y(2) + rxt(64)*y(10))
         mat(177) = -rxt(63)*y(24)
         mat(155) = -rxt(64)*y(24)

         mat(129) = .500_r8*rxt(27)*y(19)
         mat(177) = mat(177) + rxt(62)*y(23)
         mat(211) = rxt(71)*y(20)
         mat(144) = rxt(59)*y(22)
         mat(192) = rxt(57)*y(22) + rxt(70)*y(20)
         mat(66) = rxt(59)*y(6) + rxt(57)*y(8) + 4.000_r8*rxt(58)*y(22)
         mat(39) = rxt(62)*y(2)
         mat(14) = .500_r8*rxt(27)*y(1)
         mat(57) = rxt(71)*y(3) + rxt(70)*y(8)

         mat(105) = -(rxt(49)*y(9) + rxt(65)*y(8) + rxt(66)*y(6) + 4._r8*rxt(67)*y(25))
         mat(115) = -rxt(49)*y(25)
         mat(193) = -rxt(65)*y(25)
         mat(145) = -rxt(66)*y(25)

         mat(130) = .500_r8*rxt(76)*y(28)
         mat(178) = .500_r8*rxt(63)*y(24)
         mat(193) = mat(193) + .500_r8*rxt(78)*y(29)
         mat(156) = rxt(64)*y(24)
         mat(97) = .500_r8*rxt(63)*y(2) + rxt(64)*y(10)
         mat(88) = .500_r8*rxt(76)*y(1)
         mat(50) = .500_r8*rxt(78)*y(8)

         mat(11) = -(rxt(26)*y(2) + rxt(27)*y(1))
         mat(165) = -rxt(26)*y(19)
         mat(123) = -rxt(27)*y(19)

         mat(56) = -(rxt(70)*y(8) + rxt(71)*y(3) + rxt(72)*y(6))
         mat(188) = -rxt(70)*y(20)
         mat(207) = -rxt(71)*y(20)
         mat(140) = -rxt(72)*y(20)

         mat(173) = rxt(26)*y(19)
         mat(13) = rxt(26)*y(2)


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


         mat(17) = -(rxt(28)*y(1) + rxt(29)*y(2))
         mat(124) = -rxt(28)*y(26)
         mat(166) = -rxt(29)*y(26)

         mat(75) = -(rxt(73)*y(8) + rxt(74)*y(3) + rxt(75)*y(6))
         mat(190) = -rxt(73)*y(27)
         mat(209) = -rxt(74)*y(27)
         mat(142) = -rxt(75)*y(27)

         mat(175) = rxt(29)*y(26)
         mat(19) = rxt(29)*y(2)

         mat(86) = -(rxt(76)*y(1) + rxt(77)*y(2))
         mat(128) = -rxt(76)*y(28)
         mat(176) = -rxt(77)*y(28)

         mat(128) = mat(128) + rxt(28)*y(26)
         mat(143) = rxt(75)*y(27)
         mat(191) = rxt(73)*y(27)
         mat(20) = rxt(28)*y(1)
         mat(76) = rxt(75)*y(6) + rxt(73)*y(8)

         mat(47) = -(rxt(78)*y(8) + rxt(79)*y(3))
         mat(187) = -rxt(78)*y(29)
         mat(206) = -rxt(79)*y(29)

         mat(172) = rxt(77)*y(28)
         mat(84) = rxt(77)*y(2)


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


         mat(   1) = mat(   1) + lmat(   1)
         mat(   2) = mat(   2) + lmat(   2)
         mat(   4) = lmat(   4)
         mat(   5) = lmat(   5)
         mat(   6) = lmat(   6)
         mat(   7) = lmat(   7)
         mat(   8) = lmat(   8)
         mat(   9) = lmat(   9)
         mat(  10) = lmat(  10)
         mat(  11) = mat(  11) + lmat(  11)
         mat(  17) = mat(  17) + lmat(  17)
         mat(  23) = mat(  23) + lmat(  23)
         mat(  24) = mat(  24) + lmat(  24)
         mat(  26) = mat(  26) + lmat(  26)
         mat(  27) = lmat(  27)
         mat(  28) = mat(  28) + lmat(  28)
         mat(  29) = lmat(  29)
         mat(  31) = mat(  31) + lmat(  31)
         mat(  32) = mat(  32) + lmat(  32)
         mat(  33) = mat(  33) + lmat(  33)
         mat(  34) = lmat(  34)
         mat(  35) = mat(  35) + lmat(  35)
         mat(  36) = lmat(  36)
         mat(  37) = mat(  37) + lmat(  37)
         mat(  39) = mat(  39) + lmat(  39)
         mat(  40) = mat(  40) + lmat(  40)
         mat(  41) = lmat(  41)
         mat(  42) = mat(  42) + lmat(  42)
         mat(  44) = mat(  44) + lmat(  44)
         mat(  47) = mat(  47) + lmat(  47)
         mat(  56) = mat(  56) + lmat(  56)
         mat(  65) = mat(  65) + lmat(  65)
         mat(  75) = mat(  75) + lmat(  75)
         mat(  83) = mat(  83) + lmat(  83)
         mat(  86) = mat(  86) + lmat(  86)
         mat(  91) = mat(  91) + lmat(  91)
         mat(  94) = lmat(  94)
         mat(  96) = mat(  96) + lmat(  96)
         mat(  99) = mat(  99) + lmat(  99)
         mat( 102) = lmat( 102)
         mat( 105) = mat( 105) + lmat( 105)
         mat( 116) = mat( 116) + lmat( 116)
         mat( 117) = mat( 117) + lmat( 117)
         mat( 121) = lmat( 121)
         mat( 132) = mat( 132) + lmat( 132)
         mat( 135) = mat( 135) + lmat( 135)
         mat( 137) = mat( 137) + lmat( 137)
         mat( 148) = mat( 148) + lmat( 148)
         mat( 154) = mat( 154) + lmat( 154)
         mat( 157) = mat( 157) + lmat( 157)
         mat( 158) = lmat( 158)
         mat( 160) = mat( 160) + lmat( 160)
         mat( 162) = mat( 162) + lmat( 162)
         mat( 171) = mat( 171) + lmat( 171)
         mat( 174) = mat( 174) + lmat( 174)
         mat( 177) = mat( 177) + lmat( 177)
         mat( 178) = mat( 178) + lmat( 178)
         mat( 181) = mat( 181) + lmat( 181)
         mat( 183) = mat( 183) + lmat( 183)
         mat( 185) = mat( 185) + lmat( 185)
         mat( 199) = mat( 199) + lmat( 199)
         mat( 219) = mat( 219) + lmat( 219)
         mat(  48) = 0._r8
         mat(  49) = 0._r8
         mat(  52) = 0._r8
         mat(  69) = 0._r8
         mat(  74) = 0._r8
         mat(  77) = 0._r8
         mat(  80) = 0._r8
         mat(  85) = 0._r8
         mat(  87) = 0._r8
         mat(  89) = 0._r8
         mat(  93) = 0._r8
         mat(  98) = 0._r8
         mat( 108) = 0._r8
         mat( 118) = 0._r8
         mat( 126) = 0._r8
         mat( 127) = 0._r8
         mat( 147) = 0._r8
         mat( 149) = 0._r8
         mat( 150) = 0._r8
         mat( 159) = 0._r8
         mat( 184) = 0._r8
         mat( 210) = 0._r8
         mat( 212) = 0._r8
         mat( 216) = 0._r8
         mat(   1) = mat(   1) - dti
         mat(   4) = mat(   4) - dti
         mat(   7) = mat(   7) - dti
         mat(  11) = mat(  11) - dti
         mat(  17) = mat(  17) - dti
         mat(  23) = mat(  23) - dti
         mat(  28) = mat(  28) - dti
         mat(  32) = mat(  32) - dti
         mat(  37) = mat(  37) - dti
         mat(  42) = mat(  42) - dti
         mat(  47) = mat(  47) - dti
         mat(  56) = mat(  56) - dti
         mat(  65) = mat(  65) - dti
         mat(  75) = mat(  75) - dti
         mat(  86) = mat(  86) - dti
         mat(  96) = mat(  96) - dti
         mat( 105) = mat( 105) - dti
         mat( 116) = mat( 116) - dti
         mat( 132) = mat( 132) - dti
         mat( 148) = mat( 148) - dti
         mat( 160) = mat( 160) - dti
         mat( 183) = mat( 183) - dti
         mat( 199) = mat( 199) - dti
         mat( 219) = mat( 219) - dti

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
