




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

         mat(1197) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(1132) = rxt(82)

         mat(1130) = -( rxt(82) + het_rates(2) )
         mat(1195) = rxt(3)
         mat(879) = rxt(5)
         mat(826) = rxt(6)
         mat(120) = rxt(8)
         mat(921) = rxt(10)
         mat(708) = rxt(19)
         mat(1218) = rxt(22)
         mat(70) = rxt(23)
         mat(773) = rxt(30)
         mat(1237) = rxt(85) + rxt(86)
         mat(87) = rxt(133)

         mat(1241) = -( rxt(85) + rxt(86) + rxt(88)*y(4) + rxt(89)*y(4) + rxt(91)*y(107) &
                      + rxt(92)*y(108) + rxt(93)*y(109) + rxt(94)*y(117) + rxt(95)*y(118) &
                      + rxt(96)*y(110) + rxt(97)*y(115) + rxt(98)*y(116) + rxt(99)*y(111) &
                      + rxt(100)*y(106) + rxt(101)*y(114) + rxt(102)*y(113) &
                      + rxt(103)*y(119) + rxt(104)*y(120) + rxt(105)*y(121) &
                      + rxt(106)*y(122) + rxt(107)*y(12) + rxt(108)*y(12) + rxt(109)*y(12) &
                 + het_rates(3) )
         mat(1199) = rxt(2)
         mat(710) = rxt(18)

         mat(518) = -( het_rates(18) )
         mat(782) = rxt(16)
         mat(704) = rxt(18)
         mat(1226) = rxt(109)*y(12)

         mat(526) = -( het_rates(17) )
         mat(783) = rxt(15) + rxt(16)
         mat(550) = rxt(56)
         mat(560) = 1.340_r8*rxt(62)
         mat(601) = .700_r8*rxt(63)
         mat(573) = rxt(69)
         mat(467) = rxt(71)
         mat(400) = rxt(74)
         mat(283) = .450_r8*rxt(76)
         mat(325) = 2.000_r8*rxt(77)
         mat(1143) = rxt(198)*y(105)

         mat(84) = -( rxt(133) + het_rates(5) )
         mat(832) = rxt(5)

         mat(874) = -( rxt(5) + het_rates(6) )
         mat(821) = rxt(6) + .500_r8*rxt(348)
         mat(118) = rxt(8)
         mat(916) = rxt(11)
         mat(86) = rxt(133)
         mat(1232) = 2.000_r8*rxt(88)*y(4)

         mat(820) = -( rxt(6) + rxt(348) + het_rates(7) )
         mat(117) = rxt(7) + rxt(145)
         mat(413) = rxt(9)
         mat(915) = rxt(10)
         mat(193) = rxt(13) + rxt(154)
         mat(458) = rxt(28)
         mat(260) = rxt(34)
         mat(235) = .600_r8*rxt(59) + rxt(253)
         mat(268) = rxt(60) + rxt(299)
         mat(470) = rxt(71)

         mat(1013) = -( rxt(199)*y(105) + rxt(200)*y(112) + rxt(201)*y(110) &
                      + rxt(202)*y(106) + rxt(204)*y(115) + rxt(205)*y(116) &
                      + rxt(206)*y(122) + rxt(207)*y(121) + rxt(210)*y(12) + het_rates(20) &
       )
         mat(415) = rxt(9)
         mat(195) = rxt(12)
         mat(202) = rxt(14)
         mat(707) = rxt(17)
         mat(319) = 2.000_r8*rxt(20)
         mat(448) = rxt(25)
         mat(382) = rxt(31)
         mat(293) = rxt(57)
         mat(243) = rxt(58)
         mat(141) = rxt(64)
         mat(68) = rxt(65)
         mat(179) = rxt(66)
         mat(186) = rxt(67)
         mat(146) = rxt(70)
         mat(337) = rxt(78)
         mat(132) = rxt(79)
         mat(163) = rxt(80)
         mat(221) = rxt(81)
         mat(823) = .500_r8*rxt(348)
         mat(1234) = rxt(107)*y(12)

         mat(917) = -( rxt(10) + rxt(11) + rxt(347) + het_rates(8) )
         mat(119) = rxt(7) + rxt(8) + rxt(145)
         mat(194) = rxt(12)
         mat(459) = rxt(27)
         mat(261) = rxt(33)
         mat(236) = .400_r8*rxt(59)

         mat(411) = -( rxt(9) + het_rates(9) )
         mat(116) = 2.000_r8*rxt(346) + 2.000_r8*rxt(357) + 2.000_r8*rxt(363) &
                      + 2.000_r8*rxt(368)
         mat(893) = rxt(347)
         mat(808) = .500_r8*rxt(348)
         mat(454) = rxt(358) + rxt(364) + rxt(369)
         mat(258) = rxt(359) + rxt(367) + rxt(370)

         mat(191) = -( rxt(12) + rxt(13) + rxt(154) + het_rates(10) )

         mat(115) = -( rxt(7) + rxt(8) + rxt(145) + rxt(346) + rxt(357) + rxt(363) &
                      + rxt(368) + het_rates(11) )

         mat(746) = -( het_rates(13) )
         mat(555) = rxt(56)
         mat(241) = rxt(58)
         mat(233) = .400_r8*rxt(59)
         mat(610) = .300_r8*rxt(63)
         mat(366) = rxt(68)
         mat(1229) = rxt(107)*y(12)
         mat(1148) = rxt(161)*y(12)
         mat(1007) = rxt(210)*y(12)

         mat(197) = -( rxt(14) + het_rates(14) )

         mat(88) = -( het_rates(39) )

         mat(45) = -( het_rates(40) )

         mat(786) = -( rxt(15) + rxt(16) + het_rates(16) )
         mat(201) = rxt(14)
         mat(292) = rxt(57)
         mat(568) = 1.340_r8*rxt(61)
         mat(185) = rxt(67)
         mat(469) = rxt(71)
         mat(253) = .690_r8*rxt(72)
         mat(546) = rxt(73)
         mat(401) = rxt(74)
         mat(336) = .100_r8*rxt(78)
         mat(226) = rxt(224)
         mat(107) = 2.000_r8*rxt(234)
         mat(1230) = rxt(108)*y(12) + rxt(109)*y(12)

         mat(714) = -( rxt(114) + het_rates(19) )
         mat(199) = rxt(14)
         mat(785) = 2.000_r8*rxt(15)
         mat(706) = rxt(17) + 2.000_r8*rxt(19)
         mat(1250) = rxt(26)
         mat(406) = rxt(32)
         mat(1228) = rxt(108)*y(12)

         mat(1085) = -( rxt(356) + het_rates(21) )
         mat(196) = rxt(13) + rxt(154)
         mat(559) = rxt(56)
         mat(294) = rxt(57)
         mat(570) = 1.340_r8*rxt(61) + .660_r8*rxt(62)
         mat(142) = rxt(64)
         mat(180) = rxt(66)
         mat(581) = rxt(69)
         mat(473) = rxt(71)
         mat(255) = rxt(72)
         mat(548) = rxt(73)
         mat(403) = 2.000_r8*rxt(74)
         mat(286) = .560_r8*rxt(76)
         mat(327) = 2.000_r8*rxt(77)
         mat(338) = .900_r8*rxt(78)
         mat(222) = rxt(81)
         mat(716) = rxt(114)
         mat(230) = rxt(224)
         mat(108) = rxt(233) + rxt(234)
         mat(1235) = rxt(108)*y(12)
         mat(1154) = rxt(198)*y(105) + rxt(203)*y(106)
         mat(1014) = rxt(199)*y(105) + rxt(202)*y(106)

         mat(317) = -( rxt(20) + het_rates(22) )
         mat(1047) = .500_r8*rxt(356)

         mat(705) = -( rxt(17) + rxt(18) + rxt(19) + het_rates(124) )
         mat(1005) = rxt(199)*y(105) + rxt(200)*y(112) + rxt(201)*y(110) + rxt(202)*y(106) &
                      + rxt(206)*y(122) + rxt(210)*y(12)

         mat(1157) = -( rxt(161)*y(12) + rxt(198)*y(105) + rxt(203)*y(106) &
                      + rxt(208)*y(122) + rxt(209)*y(121) + het_rates(25) )
         mat(79) = 2.000_r8*rxt(21)
         mat(1219) = rxt(22)
         mat(49) = 2.000_r8*rxt(24)
         mat(450) = rxt(25)
         mat(1259) = rxt(26)
         mat(462) = rxt(27)
         mat(94) = rxt(29)
         mat(1238) = 3.000_r8*rxt(91)*y(107) + 2.000_r8*rxt(92)*y(108) &
                      + 3.000_r8*rxt(93)*y(109) + 2.000_r8*rxt(94)*y(117) + rxt(95)*y(118) &
                      + rxt(96)*y(110) + 2.000_r8*rxt(97)*y(115) + rxt(98)*y(116) &
                      + 4.000_r8*rxt(99)*y(111) + rxt(101)*y(114)
         mat(1017) = rxt(199)*y(105) + 3.000_r8*rxt(200)*y(112) + rxt(201)*y(110) &
                      + 2.000_r8*rxt(204)*y(115) + rxt(205)*y(116)

         mat(78) = -( rxt(21) + het_rates(26) )

         mat(1221) = -( rxt(22) + het_rates(27) )
         mat(71) = rxt(23)
         mat(463) = rxt(28)
         mat(50) = 2.000_r8*rxt(173)

         mat(69) = -( rxt(23) + het_rates(28) )

         mat(48) = -( rxt(24) + rxt(173) + het_rates(29) )

         mat(1263) = -( rxt(26) + het_rates(30) )
         mat(1161) = rxt(161)*y(12) + 2.000_r8*rxt(198)*y(105) + rxt(203)*y(106) &
                      + rxt(208)*y(122) + rxt(209)*y(121)

         mat(446) = -( rxt(25) + het_rates(31) )
         mat(455) = rxt(358) + rxt(364) + rxt(369)

         mat(456) = -( rxt(27) + rxt(28) + rxt(358) + rxt(364) + rxt(369) + het_rates(32) &
       )

         mat(92) = -( rxt(29) + het_rates(33) )

         mat(1104) = -( het_rates(34) )
         mat(93) = rxt(29)
         mat(772) = rxt(30)
         mat(383) = rxt(31)
         mat(408) = rxt(32)
         mat(262) = rxt(33)
         mat(1236) = rxt(100)*y(106) + rxt(101)*y(114) + rxt(102)*y(113) &
                      + 2.000_r8*rxt(103)*y(119) + 2.000_r8*rxt(104)*y(120) &
                      + 3.000_r8*rxt(105)*y(121) + 2.000_r8*rxt(106)*y(122)
         mat(1015) = rxt(202)*y(106) + 2.000_r8*rxt(206)*y(122) + 3.000_r8*rxt(207)*y(121)
         mat(1155) = rxt(203)*y(106) + 2.000_r8*rxt(208)*y(122) + 3.000_r8*rxt(209)*y(121)

         mat(766) = -( rxt(30) + het_rates(35) )
         mat(259) = rxt(34)

         mat(404) = -( rxt(32) + het_rates(36) )

         mat(379) = -( rxt(31) + het_rates(37) )
         mat(257) = rxt(359) + rxt(367) + rxt(370)


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

         mat(256) = -( rxt(33) + rxt(34) + rxt(359) + rxt(367) + rxt(370) + het_rates(38) &
       )

         mat(490) = -( het_rates(56) )
         mat(600) = .700_r8*rxt(63)

         mat(430) = -( het_rates(80) )

         mat(369) = -( het_rates(61) )

         mat(551) = -( rxt(56) + het_rates(47) )
         mat(290) = rxt(57)
         mat(140) = rxt(64)
         mat(334) = .400_r8*rxt(78)
         mat(130) = rxt(79)

         mat(313) = -( het_rates(46) )

         mat(287) = -( rxt(57) + het_rates(62) )

         mat(693) = -( het_rates(45) )
         mat(232) = .600_r8*rxt(59) + rxt(253)
         mat(565) = 1.340_r8*rxt(61)
         mat(607) = .300_r8*rxt(63)
         mat(183) = rxt(67)
         mat(364) = rxt(68)
         mat(575) = rxt(69)
         mat(545) = rxt(73)
         mat(206) = rxt(75)
         mat(285) = .130_r8*rxt(76)
         mat(131) = rxt(79)

         mat(238) = -( rxt(58) + het_rates(51) )

         mat(231) = -( rxt(59) + rxt(253) + het_rates(55) )

         mat(171) = -( het_rates(79) )

         mat(109) = -( het_rates(42) )

         mat(147) = -( het_rates(41) )

         mat(30) = -( het_rates(68) )

         mat(264) = -( rxt(60) + rxt(299) + het_rates(78) )

         mat(33) = -( het_rates(67) )

         mat(121) = -( het_rates(70) )

         mat(351) = -( het_rates(81) )

         mat(329) = -( rxt(78) + het_rates(82) )

         mat(203) = -( rxt(75) + het_rates(69) )
         mat(328) = .800_r8*rxt(78)

         mat(340) = -( het_rates(71) )

         mat(128) = -( rxt(79) + het_rates(72) )

         mat(58) = -( het_rates(91) )

         mat(63) = -( het_rates(92) )

         mat(273) = -( het_rates(93) )

         mat(158) = -( rxt(80) + het_rates(94) )

         mat(80) = -( het_rates(95) )

         mat(476) = -( het_rates(103) )

         mat(216) = -( rxt(81) + het_rates(104) )

         mat(281) = -( rxt(76) + het_rates(83) )
         mat(160) = .900_r8*rxt(80)

         mat(324) = -( rxt(77) + het_rates(50) )
         mat(282) = .130_r8*rxt(76)
         mat(161) = .450_r8*rxt(80)

         mat(36) = -( het_rates(96) )

         mat(165) = -( het_rates(97) )

         mat(1) = -( het_rates(98) )

         mat(39) = -( het_rates(99) )

         mat(209) = -( het_rates(100) )

         mat(2) = -( het_rates(101) )

         mat(648) = -( het_rates(85) )

         mat(605) = -( rxt(63) + het_rates(74) )
         mat(251) = .402_r8*rxt(72)
         mat(220) = rxt(81)

         mat(561) = -( rxt(61) + rxt(62) + het_rates(75) )
         mat(249) = .288_r8*rxt(72)
         mat(219) = rxt(81)

         mat(626) = -( het_rates(76) )

         mat(133) = -( het_rates(77) )

         mat(666) = -( het_rates(73) )
         mat(266) = rxt(60) + rxt(299)
         mat(564) = .660_r8*rxt(61)

         mat(390) = -( het_rates(43) )
         mat(205) = rxt(75)

         mat(138) = -( rxt(64) + het_rates(44) )

         mat(295) = -( het_rates(102) )

         mat(51) = -( het_rates(57) )

         mat(418) = -( het_rates(58) )

         mat(175) = -( rxt(66) + het_rates(59) )

         mat(362) = -( rxt(68) + het_rates(60) )
         mat(176) = .820_r8*rxt(66)
         mat(332) = .250_r8*rxt(78)
         mat(217) = .100_r8*rxt(81)

         mat(181) = -( rxt(67) + het_rates(66) )

         mat(244) = -( het_rates(15) )

         mat(101) = -( het_rates(48) )

         mat(399) = -( rxt(74) + het_rates(49) )
         mat(106) = rxt(233)

         mat(543) = -( rxt(73) + het_rates(63) )

         mat(306) = -( het_rates(52) )

         mat(105) = -( rxt(233) + rxt(234) + het_rates(53) )
         mat(67) = rxt(65)

         mat(66) = -( rxt(65) + het_rates(54) )

         mat(155) = -( het_rates(84) )

         mat(532) = -( het_rates(64) )

         mat(574) = -( rxt(69) + het_rates(65) )
         mat(284) = .180_r8*rxt(76)
         mat(162) = .450_r8*rxt(80)

         mat(506) = -( het_rates(86) )

         mat(466) = -( rxt(71) + het_rates(87) )

         mat(589) = -( het_rates(88) )

         mat(143) = -( rxt(70) + het_rates(89) )

         mat(248) = -( rxt(72) + het_rates(90) )

         mat(72) = -( het_rates(125) )

         mat(187) = -( het_rates(126) )

         mat(224) = -( rxt(224) + het_rates(127) )

         mat(56) = -( het_rates(142) )

         mat(96) = -( het_rates(143) )

         mat(3) = -( het_rates(144) )

         mat(42) = -( het_rates(145) )

         mat(4) = -( het_rates(146) )

         mat(5) = -( het_rates(147) )

         mat(6) = -( het_rates(132) )

         mat(7) = -( het_rates(133) )

         mat(8) = -( het_rates(134) )

         mat(9) = -( het_rates(135) )

         mat(10) = -( het_rates(136) )

         mat(11) = -( het_rates(137) )

         mat(12) = -( het_rates(138) )

         mat(13) = -( het_rates(139) )

         mat(14) = -( het_rates(140) )

         mat(15) = -( het_rates(141) )

         mat(16) = -( rxt(349) + het_rates(128) )

         mat(18) = -( het_rates(129) )
         mat(17) = rxt(349)

         mat(19) = -( rxt(355) + het_rates(130) )

         mat(21) = -( het_rates(131) )
         mat(20) = rxt(355)

         mat(22) = -( het_rates(148) )

         mat(23) = -( het_rates(149) )

         mat(24) = -( het_rates(150) )

         mat(25) = -( het_rates(151) )

         mat(26) = -( het_rates(152) )

         mat(27) = -( het_rates(153) )

         mat(28) = -( het_rates(154) )

         mat(29) = -( het_rates(155) )


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
