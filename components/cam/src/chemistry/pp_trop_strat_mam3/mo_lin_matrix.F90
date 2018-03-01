




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

         mat(1007) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(1032) = rxt(82)

         mat(1033) = -( rxt(82) + het_rates(2) )
         mat(1008) = rxt(3)
         mat(1192) = rxt(5)
         mat(798) = rxt(6)
         mat(103) = rxt(8)
         mat(930) = rxt(10)
         mat(682) = rxt(19)
         mat(1123) = rxt(22)
         mat(54) = rxt(23)
         mat(708) = rxt(30)
         mat(949) = rxt(85) + rxt(86)
         mat(70) = rxt(133)

         mat(946) = -( rxt(85) + rxt(86) + rxt(88)*y(4) + rxt(89)*y(4) + rxt(91)*y(101) &
                      + rxt(92)*y(102) + rxt(93)*y(103) + rxt(94)*y(111) + rxt(95)*y(112) &
                      + rxt(96)*y(104) + rxt(97)*y(109) + rxt(98)*y(110) + rxt(99)*y(105) &
                      + rxt(100)*y(100) + rxt(101)*y(108) + rxt(102)*y(107) &
                      + rxt(103)*y(113) + rxt(104)*y(114) + rxt(105)*y(115) &
                      + rxt(106)*y(116) + rxt(107)*y(12) + rxt(108)*y(12) + rxt(109)*y(12) &
                 + het_rates(3) )
         mat(1005) = rxt(2)
         mat(680) = rxt(18)

         mat(487) = -( het_rates(18) )
         mat(754) = rxt(16)
         mat(676) = rxt(18)
         mat(938) = rxt(109)*y(12)

         mat(495) = -( het_rates(17) )
         mat(755) = rxt(15) + rxt(16)
         mat(522) = rxt(56)
         mat(532) = 1.340_r8*rxt(62)
         mat(573) = .700_r8*rxt(63)
         mat(545) = rxt(69)
         mat(436) = rxt(71)
         mat(369) = rxt(74)
         mat(224) = .450_r8*rxt(76)
         mat(294) = 2.000_r8*rxt(77)
         mat(1205) = rxt(198)*y(99)

         mat(68) = -( rxt(133) + het_rates(5) )
         mat(1147) = rxt(5)

         mat(1196) = -( rxt(5) + het_rates(6) )
         mat(802) = rxt(6) + .500_r8*rxt(342)
         mat(104) = rxt(8)
         mat(934) = rxt(11)
         mat(71) = rxt(133)
         mat(953) = 2.000_r8*rxt(88)*y(4)

         mat(792) = -( rxt(6) + rxt(342) + het_rates(7) )
         mat(101) = rxt(7) + rxt(145)
         mat(395) = rxt(9)
         mat(924) = rxt(10)
         mat(170) = rxt(13) + rxt(154)
         mat(427) = rxt(28)
         mat(260) = rxt(33)
         mat(196) = .600_r8*rxt(59) + rxt(253)
         mat(232) = rxt(60) + rxt(299)
         mat(439) = rxt(71)

         mat(884) = -( rxt(199)*y(99) + rxt(200)*y(106) + rxt(201)*y(104) &
                      + rxt(202)*y(100) + rxt(204)*y(109) + rxt(205)*y(110) &
                      + rxt(206)*y(116) + rxt(207)*y(115) + rxt(210)*y(12) + het_rates(20) &
       )
         mat(396) = rxt(9)
         mat(171) = rxt(12)
         mat(179) = rxt(14)
         mat(679) = rxt(17)
         mat(266) = 2.000_r8*rxt(20)
         mat(417) = rxt(25)
         mat(362) = rxt(31)
         mat(246) = rxt(57)
         mat(212) = rxt(58)
         mat(125) = rxt(64)
         mat(52) = rxt(65)
         mat(152) = rxt(66)
         mat(159) = rxt(67)
         mat(130) = rxt(70)
         mat(306) = rxt(78)
         mat(116) = rxt(79)
         mat(147) = rxt(80)
         mat(190) = rxt(81)
         mat(793) = .500_r8*rxt(342)
         mat(944) = rxt(107)*y(12)

         mat(926) = -( rxt(10) + rxt(11) + rxt(341) + het_rates(8) )
         mat(102) = rxt(7) + rxt(8) + rxt(145)
         mat(172) = rxt(12)
         mat(429) = rxt(27)
         mat(261) = rxt(32)
         mat(198) = .400_r8*rxt(59)

         mat(393) = -( rxt(9) + het_rates(9) )
         mat(100) = 2.000_r8*rxt(340) + 2.000_r8*rxt(349) + 2.000_r8*rxt(355) &
                      + 2.000_r8*rxt(360)
         mat(903) = rxt(341)
         mat(780) = .500_r8*rxt(342)
         mat(423) = rxt(350) + rxt(356) + rxt(361)
         mat(258) = rxt(351) + rxt(359) + rxt(362)

         mat(168) = -( rxt(12) + rxt(13) + rxt(154) + het_rates(10) )

         mat(99) = -( rxt(7) + rxt(8) + rxt(145) + rxt(340) + rxt(349) + rxt(355) &
                      + rxt(360) + het_rates(11) )

         mat(738) = -( het_rates(13) )
         mat(527) = rxt(56)
         mat(210) = rxt(58)
         mat(194) = .400_r8*rxt(59)
         mat(582) = .300_r8*rxt(63)
         mat(335) = rxt(68)
         mat(941) = rxt(107)*y(12)
         mat(1210) = rxt(161)*y(12)
         mat(881) = rxt(210)*y(12)

         mat(174) = -( rxt(14) + het_rates(14) )

         mat(72) = -( het_rates(39) )

         mat(32) = -( het_rates(40) )

         mat(758) = -( rxt(15) + rxt(16) + het_rates(16) )
         mat(178) = rxt(14)
         mat(245) = rxt(57)
         mat(540) = 1.340_r8*rxt(61)
         mat(158) = rxt(67)
         mat(438) = rxt(71)
         mat(253) = .690_r8*rxt(72)
         mat(518) = rxt(73)
         mat(370) = rxt(74)
         mat(305) = .100_r8*rxt(78)
         mat(202) = rxt(224)
         mat(91) = 2.000_r8*rxt(234)
         mat(942) = rxt(108)*y(12) + rxt(109)*y(12)

         mat(686) = -( rxt(114) + het_rates(19) )
         mat(176) = rxt(14)
         mat(757) = 2.000_r8*rxt(15)
         mat(678) = rxt(17) + 2.000_r8*rxt(19)
         mat(962) = rxt(26)
         mat(375) = rxt(34)
         mat(940) = rxt(108)*y(12)

         mat(1101) = -( rxt(348) + het_rates(21) )
         mat(173) = rxt(13) + rxt(154)
         mat(531) = rxt(56)
         mat(247) = rxt(57)
         mat(543) = 1.340_r8*rxt(61) + .660_r8*rxt(62)
         mat(126) = rxt(64)
         mat(153) = rxt(66)
         mat(553) = rxt(69)
         mat(442) = rxt(71)
         mat(255) = rxt(72)
         mat(520) = rxt(73)
         mat(372) = 2.000_r8*rxt(74)
         mat(227) = .560_r8*rxt(76)
         mat(296) = 2.000_r8*rxt(77)
         mat(307) = .900_r8*rxt(78)
         mat(191) = rxt(81)
         mat(692) = rxt(114)
         mat(205) = rxt(224)
         mat(92) = rxt(233) + rxt(234)
         mat(950) = rxt(108)*y(12)
         mat(1219) = rxt(198)*y(99) + rxt(203)*y(100)
         mat(890) = rxt(199)*y(99) + rxt(202)*y(100)

         mat(264) = -( rxt(20) + het_rates(22) )
         mat(1058) = .500_r8*rxt(348)

         mat(677) = -( rxt(17) + rxt(18) + rxt(19) + het_rates(145) )
         mat(878) = rxt(199)*y(99) + rxt(200)*y(106) + rxt(201)*y(104) + rxt(202)*y(100) &
                      + rxt(206)*y(116) + rxt(210)*y(12)

         mat(1223) = -( rxt(161)*y(12) + rxt(198)*y(99) + rxt(203)*y(100) &
                      + rxt(208)*y(116) + rxt(209)*y(115) + het_rates(25) )
         mat(63) = 2.000_r8*rxt(21)
         mat(1128) = rxt(22)
         mat(22) = 2.000_r8*rxt(24)
         mat(421) = rxt(25)
         mat(975) = rxt(26)
         mat(433) = rxt(27)
         mat(78) = rxt(29)
         mat(954) = 3.000_r8*rxt(91)*y(101) + 2.000_r8*rxt(92)*y(102) &
                      + 3.000_r8*rxt(93)*y(103) + 2.000_r8*rxt(94)*y(111) + rxt(95)*y(112) &
                      + rxt(96)*y(104) + 2.000_r8*rxt(97)*y(109) + rxt(98)*y(110) &
                      + 4.000_r8*rxt(99)*y(105) + rxt(101)*y(108)
         mat(894) = rxt(199)*y(99) + 3.000_r8*rxt(200)*y(106) + rxt(201)*y(104) &
                      + 2.000_r8*rxt(204)*y(109) + rxt(205)*y(110)

         mat(62) = -( rxt(21) + het_rates(26) )

         mat(1125) = -( rxt(22) + het_rates(27) )
         mat(55) = rxt(23)
         mat(432) = rxt(28)
         mat(21) = 2.000_r8*rxt(173)

         mat(53) = -( rxt(23) + het_rates(28) )

         mat(20) = -( rxt(24) + rxt(173) + het_rates(29) )

         mat(968) = -( rxt(26) + het_rates(30) )
         mat(1216) = rxt(161)*y(12) + 2.000_r8*rxt(198)*y(99) + rxt(203)*y(100) &
                      + rxt(208)*y(116) + rxt(209)*y(115)

         mat(415) = -( rxt(25) + het_rates(31) )
         mat(424) = rxt(350) + rxt(356) + rxt(361)

         mat(425) = -( rxt(27) + rxt(28) + rxt(350) + rxt(356) + rxt(361) + het_rates(32) &
       )

         mat(76) = -( rxt(29) + het_rates(33) )

         mat(1144) = -( het_rates(34) )
         mat(77) = rxt(29)
         mat(711) = rxt(30)
         mat(365) = rxt(31)
         mat(263) = rxt(32)
         mat(379) = rxt(34)
         mat(952) = rxt(100)*y(100) + rxt(101)*y(108) + rxt(102)*y(107) &
                      + 2.000_r8*rxt(103)*y(113) + 2.000_r8*rxt(104)*y(114) &
                      + 3.000_r8*rxt(105)*y(115) + 2.000_r8*rxt(106)*y(116)
         mat(892) = rxt(202)*y(100) + 2.000_r8*rxt(206)*y(116) + 3.000_r8*rxt(207)*y(115)
         mat(1221) = rxt(203)*y(100) + 2.000_r8*rxt(208)*y(116) + 3.000_r8*rxt(209)*y(115)

         mat(701) = -( rxt(30) + het_rates(35) )
         mat(259) = rxt(33)

         mat(373) = -( rxt(34) + het_rates(36) )

         mat(359) = -( rxt(31) + het_rates(37) )
         mat(257) = rxt(351) + rxt(359) + rxt(362)


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

         mat(256) = -( rxt(32) + rxt(33) + rxt(351) + rxt(359) + rxt(362) + het_rates(38) &
       )

         mat(459) = -( het_rates(56) )
         mat(572) = .700_r8*rxt(63)

         mat(399) = -( het_rates(80) )

         mat(338) = -( het_rates(61) )

         mat(523) = -( rxt(56) + het_rates(47) )
         mat(243) = rxt(57)
         mat(124) = rxt(64)
         mat(303) = .400_r8*rxt(78)
         mat(114) = rxt(79)

         mat(289) = -( het_rates(46) )

         mat(240) = -( rxt(57) + het_rates(62) )

         mat(665) = -( het_rates(45) )
         mat(193) = .600_r8*rxt(59) + rxt(253)
         mat(537) = 1.340_r8*rxt(61)
         mat(579) = .300_r8*rxt(63)
         mat(156) = rxt(67)
         mat(333) = rxt(68)
         mat(547) = rxt(69)
         mat(517) = rxt(73)
         mat(183) = rxt(75)
         mat(226) = .130_r8*rxt(76)
         mat(115) = rxt(79)

         mat(207) = -( rxt(58) + het_rates(51) )

         mat(192) = -( rxt(59) + rxt(253) + het_rates(55) )

         mat(160) = -( het_rates(79) )

         mat(93) = -( het_rates(42) )

         mat(131) = -( het_rates(41) )

         mat(23) = -( het_rates(68) )

         mat(228) = -( rxt(60) + rxt(299) + het_rates(78) )

         mat(26) = -( het_rates(67) )

         mat(105) = -( het_rates(70) )

         mat(320) = -( het_rates(81) )

         mat(298) = -( rxt(78) + het_rates(82) )

         mat(180) = -( rxt(75) + het_rates(69) )
         mat(297) = .800_r8*rxt(78)

         mat(309) = -( het_rates(71) )

         mat(112) = -( rxt(79) + het_rates(72) )

         mat(42) = -( het_rates(91) )

         mat(47) = -( het_rates(92) )

         mat(214) = -( het_rates(93) )

         mat(142) = -( rxt(80) + het_rates(94) )

         mat(64) = -( het_rates(95) )

         mat(502) = -( het_rates(97) )

         mat(185) = -( rxt(81) + het_rates(98) )

         mat(222) = -( rxt(76) + het_rates(83) )
         mat(144) = .900_r8*rxt(80)

         mat(293) = -( rxt(77) + het_rates(50) )
         mat(223) = .130_r8*rxt(76)
         mat(145) = .450_r8*rxt(80)

         mat(620) = -( het_rates(85) )

         mat(577) = -( rxt(63) + het_rates(74) )
         mat(251) = .402_r8*rxt(72)
         mat(189) = rxt(81)

         mat(533) = -( rxt(61) + rxt(62) + het_rates(75) )
         mat(249) = .288_r8*rxt(72)
         mat(188) = rxt(81)

         mat(598) = -( het_rates(76) )

         mat(117) = -( het_rates(77) )

         mat(638) = -( het_rates(73) )
         mat(230) = rxt(60) + rxt(299)
         mat(536) = .660_r8*rxt(61)

         mat(350) = -( het_rates(43) )
         mat(182) = rxt(75)

         mat(122) = -( rxt(64) + het_rates(44) )

         mat(271) = -( het_rates(96) )

         mat(35) = -( het_rates(57) )

         mat(382) = -( het_rates(58) )

         mat(148) = -( rxt(66) + het_rates(59) )

         mat(331) = -( rxt(68) + het_rates(60) )
         mat(149) = .820_r8*rxt(66)
         mat(301) = .250_r8*rxt(78)
         mat(186) = .100_r8*rxt(81)

         mat(154) = -( rxt(67) + het_rates(66) )

         mat(236) = -( het_rates(15) )

         mat(85) = -( het_rates(48) )

         mat(368) = -( rxt(74) + het_rates(49) )
         mat(90) = rxt(233)

         mat(515) = -( rxt(73) + het_rates(63) )

         mat(282) = -( het_rates(52) )

         mat(89) = -( rxt(233) + rxt(234) + het_rates(53) )
         mat(51) = rxt(65)

         mat(50) = -( rxt(65) + het_rates(54) )

         mat(139) = -( het_rates(84) )

         mat(445) = -( het_rates(64) )

         mat(546) = -( rxt(69) + het_rates(65) )
         mat(225) = .180_r8*rxt(76)
         mat(146) = .450_r8*rxt(80)

         mat(475) = -( het_rates(86) )

         mat(435) = -( rxt(71) + het_rates(87) )

         mat(561) = -( het_rates(88) )

         mat(127) = -( rxt(70) + het_rates(89) )

         mat(248) = -( rxt(72) + het_rates(90) )

         mat(56) = -( het_rates(118) )

         mat(164) = -( het_rates(119) )

         mat(200) = -( rxt(224) + het_rates(120) )

         mat(40) = -( het_rates(121) )

         mat(80) = -( het_rates(122) )

         mat(29) = -( het_rates(123) )

         mat(1) = -( het_rates(124) )

         mat(2) = -( het_rates(129) )

         mat(3) = -( het_rates(127) )

         mat(4) = -( het_rates(128) )

         mat(5) = -( het_rates(130) )

         mat(6) = -( het_rates(131) )

         mat(7) = -( het_rates(132) )

         mat(8) = -( het_rates(133) )

         mat(9) = -( het_rates(134) )

         mat(10) = -( het_rates(135) )

         mat(11) = -( het_rates(136) )

         mat(12) = -( het_rates(137) )

         mat(13) = -( het_rates(138) )

         mat(14) = -( het_rates(139) )

         mat(15) = -( het_rates(140) )

         mat(16) = -( het_rates(141) )

         mat(17) = -( het_rates(142) )

         mat(18) = -( het_rates(143) )

         mat(19) = -( het_rates(144) )


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
