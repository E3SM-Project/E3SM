




      module mo_lin_matrix

      private
      public :: linmat

      contains

      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

         mat(820) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(530) = rxt(43)

         mat(526) = -( rxt(43) + het_rates(2) )
         mat(811) = rxt(3)
         mat(836) = rxt(5)
         mat(104) = rxt(7)
         mat(956) = rxt(9)
         mat(447) = rxt(45) + rxt(46)

         mat(446) = -( rxt(45) + rxt(46) + rxt(48) + rxt(50)*y(4) + rxt(51)*y(4) &
                      + rxt(52)*y(16) + rxt(53)*y(16) + rxt(54)*y(16) + het_rates(3) )
         mat(804) = rxt(2)

         mat(890) = -( het_rates(17) )
         mat(729) = rxt(14) + rxt(15)
         mat(432) = rxt(17)
         mat(497) = 1.340_r8*rxt(23)
         mat(549) = .700_r8*rxt(24)
         mat(504) = rxt(30)
         mat(441) = rxt(32)
         mat(362) = rxt(35)
         mat(256) = .450_r8*rxt(37)
         mat(295) = 2.000_r8*rxt(38)

         mat(278) = -( het_rates(11) )
         mat(719) = rxt(15)
         mat(444) = rxt(54)*y(16)

         mat(77) = -( het_rates(116) )

         mat(50) = -( het_rates(117) )

         mat(391) = -( rxt(56) + het_rates(15) )
         mat(146) = rxt(13)
         mat(445) = rxt(53)*y(16)

         mat(935) = -( het_rates(5) )
         mat(846) = rxt(5) + .500_r8*rxt(229)
         mat(106) = rxt(7)
         mat(969) = rxt(10)
         mat(453) = 2.000_r8*rxt(50)*y(4)

         mat(843) = -( rxt(5) + rxt(229) + het_rates(6) )
         mat(105) = rxt(6) + rxt(82)
         mat(225) = rxt(8)
         mat(966) = rxt(9)
         mat(136) = rxt(12) + rxt(91)
         mat(220) = .600_r8*rxt(20) + rxt(134)
         mat(240) = rxt(21) + rxt(180)
         mat(440) = rxt(32)

         mat(970) = -( rxt(9) + rxt(10) + rxt(228) + het_rates(7) )
         mat(107) = rxt(6) + rxt(7) + rxt(82)
         mat(137) = rxt(11)
         mat(222) = .400_r8*rxt(20)

         mat(223) = -( rxt(8) + het_rates(8) )
         mat(103) = 2.000_r8*rxt(227)
         mat(942) = rxt(228)
         mat(830) = .500_r8*rxt(229)

         mat(133) = -( rxt(11) + rxt(12) + rxt(91) + het_rates(9) )

         mat(102) = -( rxt(6) + rxt(7) + rxt(82) + rxt(227) + het_rates(10) )

         mat(708) = -( rxt(92)*y(16) + het_rates(12) )
         mat(224) = rxt(8)
         mat(134) = rxt(11)
         mat(147) = rxt(13)
         mat(92) = 2.000_r8*rxt(16)
         mat(231) = rxt(18)
         mat(195) = rxt(19)
         mat(131) = rxt(25)
         mat(66) = rxt(26)
         mat(126) = rxt(27)
         mat(121) = rxt(28)
         mat(83) = rxt(31)
         mat(304) = rxt(39)
         mat(112) = rxt(40)
         mat(169) = rxt(41)
         mat(206) = rxt(42)
         mat(448) = 2.000_r8*rxt(48) + rxt(52)*y(16)
         mat(839) = .500_r8*rxt(229)

         mat(787) = -( rxt(237) + het_rates(13) )
         mat(135) = rxt(12) + rxt(91)
         mat(429) = rxt(17)
         mat(233) = rxt(18)
         mat(495) = 1.340_r8*rxt(22) + .660_r8*rxt(23)
         mat(132) = rxt(25)
         mat(127) = rxt(27)
         mat(502) = rxt(30)
         mat(439) = rxt(32)
         mat(264) = rxt(33)
         mat(487) = rxt(34)
         mat(361) = 2.000_r8*rxt(35)
         mat(255) = .560_r8*rxt(37)
         mat(294) = 2.000_r8*rxt(38)
         mat(306) = .900_r8*rxt(39)
         mat(207) = rxt(42)
         mat(395) = rxt(56)
         mat(185) = rxt(106)
         mat(101) = rxt(114) + rxt(115)
         mat(450) = rxt(53)*y(16)

         mat(90) = -( rxt(16) + het_rates(14) )
         mat(736) = .500_r8*rxt(237)

         mat(877) = -( het_rates(18) )
         mat(431) = rxt(17)
         mat(197) = rxt(19)
         mat(221) = .400_r8*rxt(20)
         mat(548) = .300_r8*rxt(24)
         mat(319) = rxt(29)
         mat(452) = rxt(52)*y(16)
         mat(713) = rxt(92)*y(16)

         mat(145) = -( rxt(13) + het_rates(19) )

         mat(724) = -( rxt(14) + rxt(15) + het_rates(20) )
         mat(148) = rxt(13)
         mat(232) = rxt(18)
         mat(494) = 1.340_r8*rxt(22)
         mat(122) = rxt(28)
         mat(438) = rxt(32)
         mat(263) = .690_r8*rxt(33)
         mat(486) = rxt(34)
         mat(360) = rxt(35)
         mat(305) = .100_r8*rxt(39)
         mat(184) = rxt(106)
         mat(100) = 2.000_r8*rxt(115)
         mat(449) = rxt(53)*y(16) + rxt(54)*y(16)

         mat(265) = -( het_rates(21) )

         mat(94) = -( het_rates(22) )

         mat(64) = -( rxt(26) + het_rates(28) )

         mat(154) = -( het_rates(23) )

         mat(98) = -( rxt(114) + rxt(115) + het_rates(24) )
         mat(65) = rxt(26)

         mat(285) = -( het_rates(25) )

         mat(198) = -( het_rates(26) )

         mat(358) = -( rxt(35) + het_rates(27) )
         mat(99) = rxt(114)

         mat(32) = -( het_rates(29) )

         mat(349) = -( het_rates(30) )
         mat(190) = rxt(36)

         mat(128) = -( rxt(25) + het_rates(31) )

         mat(426) = -( rxt(17) + het_rates(32) )
         mat(229) = rxt(18)
         mat(130) = rxt(25)
         mat(302) = .400_r8*rxt(39)
         mat(110) = rxt(40)

         mat(621) = -( het_rates(33) )
         mat(217) = .600_r8*rxt(20) + rxt(134)
         mat(492) = 1.340_r8*rxt(22)
         mat(542) = .300_r8*rxt(24)
         mat(120) = rxt(28)
         mat(317) = rxt(29)
         mat(500) = rxt(30)
         mat(484) = rxt(34)
         mat(191) = rxt(36)
         mat(254) = .130_r8*rxt(37)
         mat(111) = rxt(40)

         mat(193) = -( rxt(19) + het_rates(34) )

         mat(412) = -( het_rates(35) )
         mat(536) = .700_r8*rxt(24)

         mat(35) = -( het_rates(36) )

         mat(365) = -( het_rates(37) )

         mat(123) = -( rxt(27) + het_rates(38) )

         mat(321) = -( het_rates(39) )

         mat(227) = -( rxt(18) + het_rates(40) )

         mat(315) = -( rxt(29) + het_rates(41) )
         mat(124) = .820_r8*rxt(27)
         mat(299) = .250_r8*rxt(39)
         mat(202) = .100_r8*rxt(42)

         mat(399) = -( het_rates(42) )

         mat(118) = -( rxt(28) + het_rates(43) )

         mat(38) = -( het_rates(44) )

         mat(138) = -( het_rates(45) )

         mat(41) = -( het_rates(49) )

         mat(334) = -( het_rates(50) )

         mat(297) = -( rxt(39) + het_rates(51) )

         mat(188) = -( rxt(36) + het_rates(46) )
         mat(296) = .800_r8*rxt(39)

         mat(308) = -( het_rates(47) )

         mat(108) = -( rxt(40) + het_rates(48) )

         mat(376) = -( het_rates(52) )

         mat(579) = -( het_rates(53) )

         mat(257) = -( rxt(33) + het_rates(54) )

         mat(540) = -( rxt(24) + het_rates(55) )
         mat(260) = .402_r8*rxt(33)
         mat(205) = rxt(42)

         mat(488) = -( rxt(22) + rxt(23) + het_rates(56) )
         mat(258) = .288_r8*rxt(33)
         mat(204) = rxt(42)

         mat(559) = -( het_rates(57) )

         mat(113) = -( het_rates(58) )

         mat(596) = -( het_rates(59) )
         mat(236) = rxt(21) + rxt(180)
         mat(491) = .660_r8*rxt(22)

         mat(161) = -( het_rates(60) )

         mat(482) = -( rxt(34) + het_rates(61) )

         mat(499) = -( rxt(30) + het_rates(62) )
         mat(253) = .180_r8*rxt(37)
         mat(168) = .450_r8*rxt(41)

         mat(512) = -( het_rates(63) )

         mat(81) = -( rxt(31) + het_rates(64) )

         mat(269) = -( het_rates(65) )

         mat(470) = -( het_rates(66) )

         mat(201) = -( rxt(42) + het_rates(67) )

         mat(56) = -( het_rates(68) )

         mat(61) = -( het_rates(69) )

         mat(243) = -( het_rates(70) )

         mat(164) = -( rxt(41) + het_rates(71) )

         mat(73) = -( het_rates(78) )

         mat(251) = -( rxt(37) + het_rates(79) )
         mat(166) = .900_r8*rxt(41)

         mat(292) = -( rxt(38) + het_rates(80) )
         mat(252) = .130_r8*rxt(37)
         mat(167) = .450_r8*rxt(41)


      end subroutine linmat01

      subroutine linmat02( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

         mat(44) = -( het_rates(72) )

         mat(171) = -( het_rates(73) )

         mat(1) = -( het_rates(74) )

         mat(47) = -( het_rates(75) )

         mat(209) = -( het_rates(76) )

         mat(2) = -( het_rates(77) )

         mat(216) = -( rxt(20) + rxt(134) + het_rates(81) )

         mat(177) = -( het_rates(82) )

         mat(234) = -( rxt(21) + rxt(180) + het_rates(83) )

         mat(456) = -( het_rates(84) )

         mat(435) = -( rxt(32) + het_rates(85) )

         mat(54) = -( het_rates(100) )

         mat(85) = -( het_rates(101) )

         mat(3) = -( het_rates(102) )

         mat(30) = -( het_rates(103) )

         mat(4) = -( het_rates(104) )

         mat(5) = -( het_rates(105) )

         mat(6) = -( het_rates(90) )

         mat(7) = -( het_rates(91) )

         mat(8) = -( het_rates(92) )

         mat(9) = -( het_rates(93) )

         mat(10) = -( het_rates(94) )

         mat(11) = -( het_rates(95) )

         mat(12) = -( het_rates(96) )

         mat(13) = -( het_rates(97) )

         mat(14) = -( het_rates(98) )

         mat(15) = -( het_rates(99) )

         mat(16) = -( rxt(230) + het_rates(86) )

         mat(18) = -( het_rates(87) )
         mat(17) = rxt(230)

         mat(19) = -( rxt(236) + het_rates(88) )

         mat(21) = -( het_rates(89) )
         mat(20) = rxt(236)

         mat(67) = -( het_rates(118) )

         mat(150) = -( het_rates(119) )

         mat(182) = -( rxt(106) + het_rates(120) )

         mat(22) = -( het_rates(106) )

         mat(23) = -( het_rates(107) )

         mat(24) = -( het_rates(108) )

         mat(25) = -( het_rates(109) )

         mat(26) = -( het_rates(110) )

         mat(27) = -( het_rates(111) )

         mat(28) = -( het_rates(112) )

         mat(29) = -( het_rates(113) )


      end subroutine linmat02

      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

      call linmat01( mat, y, rxt, het_rates )
      call linmat02( mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix
