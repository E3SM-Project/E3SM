




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

         mat(951) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(1176) = rxt(82)

         mat(1181) = -( rxt(82) + het_rates(2) )
         mat(956) = rxt(3)
         mat(852) = rxt(5)
         mat(1236) = rxt(6)
         mat(112) = rxt(8)
         mat(917) = rxt(10)
         mat(695) = rxt(19)
         mat(875) = rxt(22)
         mat(70) = rxt(23)
         mat(1202) = rxt(30)
         mat(1089) = rxt(85) + rxt(86)
         mat(85) = rxt(133)

         mat(1087) = -( rxt(85) + rxt(86) + rxt(88)*y(4) + rxt(89)*y(4) + rxt(91)*y(101) &
                      + rxt(92)*y(102) + rxt(93)*y(103) + rxt(94)*y(111) + rxt(95)*y(112) &
                      + rxt(96)*y(104) + rxt(97)*y(109) + rxt(98)*y(110) + rxt(99)*y(105) &
                      + rxt(100)*y(100) + rxt(101)*y(108) + rxt(102)*y(107) &
                      + rxt(103)*y(113) + rxt(104)*y(114) + rxt(105)*y(115) &
                      + rxt(106)*y(116) + rxt(107)*y(12) + rxt(108)*y(12) + rxt(109)*y(12) &
                 + het_rates(3) )
         mat(954) = rxt(2)
         mat(694) = rxt(18)

         mat(502) = -( het_rates(18) )
         mat(760) = rxt(16)
         mat(688) = rxt(18)
         mat(1074) = rxt(109)*y(12)

         mat(510) = -( het_rates(17) )
         mat(761) = rxt(15) + rxt(16)
         mat(534) = rxt(56)
         mat(544) = 1.340_r8*rxt(62)
         mat(585) = .700_r8*rxt(63)
         mat(557) = rxt(69)
         mat(451) = rxt(71)
         mat(384) = rxt(74)
         mat(243) = .450_r8*rxt(76)
         mat(309) = 2.000_r8*rxt(77)
         mat(784) = rxt(198)*y(99)

         mat(83) = -( rxt(133) + het_rates(5) )
         mat(804) = rxt(5)

         mat(844) = -( rxt(5) + het_rates(6) )
         mat(1228) = rxt(6) + .500_r8*rxt(342)
         mat(110) = rxt(8)
         mat(909) = rxt(11)
         mat(84) = rxt(133)
         mat(1081) = 2.000_r8*rxt(88)*y(4)

         mat(1238) = -( rxt(6) + rxt(342) + het_rates(7) )
         mat(113) = rxt(7) + rxt(145)
         mat(412) = rxt(9)
         mat(919) = rxt(10)
         mat(182) = rxt(13) + rxt(154)
         mat(448) = rxt(28)
         mat(262) = rxt(33)
         mat(221) = .600_r8*rxt(59) + rxt(253)
         mat(278) = rxt(60) + rxt(299)
         mat(457) = rxt(71)

         mat(1066) = -( rxt(199)*y(99) + rxt(200)*y(106) + rxt(201)*y(104) &
                      + rxt(202)*y(100) + rxt(204)*y(109) + rxt(205)*y(110) &
                      + rxt(206)*y(116) + rxt(207)*y(115) + rxt(210)*y(12) + het_rates(20) &
       )
         mat(411) = rxt(9)
         mat(180) = rxt(12)
         mat(188) = rxt(14)
         mat(693) = rxt(17)
         mat(305) = 2.000_r8*rxt(20)
         mat(435) = rxt(25)
         mat(379) = rxt(31)
         mat(269) = rxt(57)
         mat(227) = rxt(58)
         mat(133) = rxt(64)
         mat(67) = rxt(65)
         mat(161) = rxt(66)
         mat(168) = rxt(67)
         mat(138) = rxt(70)
         mat(321) = rxt(78)
         mat(124) = rxt(79)
         mat(194) = rxt(80)
         mat(205) = rxt(81)
         mat(1233) = .500_r8*rxt(342)
         mat(1086) = rxt(107)*y(12)

         mat(911) = -( rxt(10) + rxt(11) + rxt(341) + het_rates(8) )
         mat(111) = rxt(7) + rxt(8) + rxt(145)
         mat(179) = rxt(12)
         mat(444) = rxt(27)
         mat(259) = rxt(32)
         mat(219) = .400_r8*rxt(59)

         mat(408) = -( rxt(9) + het_rates(9) )
         mat(109) = 2.000_r8*rxt(340) + 2.000_r8*rxt(349) + 2.000_r8*rxt(355) &
                      + 2.000_r8*rxt(360)
         mat(886) = rxt(341)
         mat(1215) = .500_r8*rxt(342)
         mat(438) = rxt(350) + rxt(356) + rxt(361)
         mat(257) = rxt(351) + rxt(359) + rxt(362)

         mat(177) = -( rxt(12) + rxt(13) + rxt(154) + het_rates(10) )

         mat(108) = -( rxt(7) + rxt(8) + rxt(145) + rxt(340) + rxt(349) + rxt(355) &
                      + rxt(360) + het_rates(11) )

         mat(730) = -( het_rates(13) )
         mat(539) = rxt(56)
         mat(225) = rxt(58)
         mat(217) = .400_r8*rxt(59)
         mat(594) = .300_r8*rxt(63)
         mat(350) = rxt(68)
         mat(1077) = rxt(107)*y(12)
         mat(789) = rxt(161)*y(12)
         mat(1057) = rxt(210)*y(12)

         mat(183) = -( rxt(14) + het_rates(14) )

         mat(87) = -( het_rates(39) )

         mat(44) = -( het_rates(40) )

         mat(765) = -( rxt(15) + rxt(16) + het_rates(16) )
         mat(187) = rxt(14)
         mat(268) = rxt(57)
         mat(552) = 1.340_r8*rxt(61)
         mat(167) = rxt(67)
         mat(453) = rxt(71)
         mat(252) = .690_r8*rxt(72)
         mat(530) = rxt(73)
         mat(385) = rxt(74)
         mat(320) = .100_r8*rxt(78)
         mat(210) = rxt(224)
         mat(106) = 2.000_r8*rxt(234)
         mat(1079) = rxt(108)*y(12) + rxt(109)*y(12)

         mat(698) = -( rxt(114) + het_rates(19) )
         mat(185) = rxt(14)
         mat(763) = 2.000_r8*rxt(15)
         mat(690) = rxt(17) + 2.000_r8*rxt(19)
         mat(966) = rxt(26)
         mat(390) = rxt(34)
         mat(1076) = rxt(108)*y(12)

         mat(1155) = -( rxt(348) + het_rates(21) )
         mat(181) = rxt(13) + rxt(154)
         mat(542) = rxt(56)
         mat(270) = rxt(57)
         mat(555) = 1.340_r8*rxt(61) + .660_r8*rxt(62)
         mat(134) = rxt(64)
         mat(162) = rxt(66)
         mat(564) = rxt(69)
         mat(456) = rxt(71)
         mat(254) = rxt(72)
         mat(532) = rxt(73)
         mat(387) = 2.000_r8*rxt(74)
         mat(246) = .560_r8*rxt(76)
         mat(311) = 2.000_r8*rxt(77)
         mat(322) = .900_r8*rxt(78)
         mat(206) = rxt(81)
         mat(704) = rxt(114)
         mat(213) = rxt(224)
         mat(107) = rxt(233) + rxt(234)
         mat(1088) = rxt(108)*y(12)
         mat(800) = rxt(198)*y(99) + rxt(203)*y(100)
         mat(1068) = rxt(199)*y(99) + rxt(202)*y(100)

         mat(301) = -( rxt(20) + het_rates(22) )
         mat(1113) = .500_r8*rxt(348)

         mat(689) = -( rxt(17) + rxt(18) + rxt(19) + het_rates(160) )
         mat(1055) = rxt(199)*y(99) + rxt(200)*y(106) + rxt(201)*y(104) + rxt(202)*y(100) &
                      + rxt(206)*y(116) + rxt(210)*y(12)

         mat(792) = -( rxt(161)*y(12) + rxt(198)*y(99) + rxt(203)*y(100) + rxt(208)*y(116) &
                      + rxt(209)*y(115) + het_rates(25) )
         mat(78) = 2.000_r8*rxt(21)
         mat(866) = rxt(22)
         mat(48) = 2.000_r8*rxt(24)
         mat(432) = rxt(25)
         mat(969) = rxt(26)
         mat(442) = rxt(27)
         mat(99) = rxt(29)
         mat(1080) = 3.000_r8*rxt(91)*y(101) + 2.000_r8*rxt(92)*y(102) &
                      + 3.000_r8*rxt(93)*y(103) + 2.000_r8*rxt(94)*y(111) + rxt(95)*y(112) &
                      + rxt(96)*y(104) + 2.000_r8*rxt(97)*y(109) + rxt(98)*y(110) &
                      + 4.000_r8*rxt(99)*y(105) + rxt(101)*y(108)
         mat(1060) = rxt(199)*y(99) + 3.000_r8*rxt(200)*y(106) + rxt(201)*y(104) &
                      + 2.000_r8*rxt(204)*y(109) + rxt(205)*y(110)

         mat(77) = -( rxt(21) + het_rates(26) )

         mat(868) = -( rxt(22) + het_rates(27) )
         mat(69) = rxt(23)
         mat(443) = rxt(28)
         mat(49) = 2.000_r8*rxt(173)

         mat(68) = -( rxt(23) + het_rates(28) )

         mat(47) = -( rxt(24) + rxt(173) + het_rates(29) )

         mat(974) = -( rxt(26) + het_rates(30) )
         mat(797) = rxt(161)*y(12) + 2.000_r8*rxt(198)*y(99) + rxt(203)*y(100) &
                      + rxt(208)*y(116) + rxt(209)*y(115)

         mat(430) = -( rxt(25) + het_rates(31) )
         mat(439) = rxt(350) + rxt(356) + rxt(361)

         mat(440) = -( rxt(27) + rxt(28) + rxt(350) + rxt(356) + rxt(361) + het_rates(32) &
       )

         mat(97) = -( rxt(29) + het_rates(33) )

         mat(747) = -( het_rates(34) )
         mat(98) = rxt(29)
         mat(1191) = rxt(30)
         mat(376) = rxt(31)
         mat(258) = rxt(32)
         mat(391) = rxt(34)
         mat(1078) = rxt(100)*y(100) + rxt(101)*y(108) + rxt(102)*y(107) &
                      + 2.000_r8*rxt(103)*y(113) + 2.000_r8*rxt(104)*y(114) &
                      + 3.000_r8*rxt(105)*y(115) + 2.000_r8*rxt(106)*y(116)
         mat(1058) = rxt(202)*y(100) + 2.000_r8*rxt(206)*y(116) + 3.000_r8*rxt(207)*y(115)
         mat(790) = rxt(203)*y(100) + 2.000_r8*rxt(208)*y(116) + 3.000_r8*rxt(209)*y(115)

         mat(1203) = -( rxt(30) + het_rates(35) )
         mat(261) = rxt(33)

         mat(388) = -( rxt(34) + het_rates(36) )

         mat(374) = -( rxt(31) + het_rates(37) )
         mat(256) = rxt(351) + rxt(359) + rxt(362)


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

         mat(255) = -( rxt(32) + rxt(33) + rxt(351) + rxt(359) + rxt(362) + het_rates(38) &
       )

         mat(474) = -( het_rates(56) )
         mat(584) = .700_r8*rxt(63)

         mat(414) = -( het_rates(80) )

         mat(353) = -( het_rates(61) )

         mat(535) = -( rxt(56) + het_rates(47) )
         mat(266) = rxt(57)
         mat(132) = rxt(64)
         mat(318) = .400_r8*rxt(78)
         mat(122) = rxt(79)

         mat(297) = -( het_rates(46) )

         mat(263) = -( rxt(57) + het_rates(62) )

         mat(677) = -( het_rates(45) )
         mat(216) = .600_r8*rxt(59) + rxt(253)
         mat(549) = 1.340_r8*rxt(61)
         mat(591) = .300_r8*rxt(63)
         mat(165) = rxt(67)
         mat(348) = rxt(68)
         mat(559) = rxt(69)
         mat(529) = rxt(73)
         mat(198) = rxt(75)
         mat(245) = .130_r8*rxt(76)
         mat(123) = rxt(79)

         mat(222) = -( rxt(58) + het_rates(51) )

         mat(215) = -( rxt(59) + rxt(253) + het_rates(55) )

         mat(169) = -( het_rates(79) )

         mat(114) = -( het_rates(42) )

         mat(146) = -( het_rates(41) )

         mat(35) = -( het_rates(68) )

         mat(271) = -( rxt(60) + rxt(299) + het_rates(78) )

         mat(38) = -( het_rates(67) )

         mat(139) = -( het_rates(70) )

         mat(335) = -( het_rates(81) )

         mat(313) = -( rxt(78) + het_rates(82) )

         mat(195) = -( rxt(75) + het_rates(69) )
         mat(312) = .800_r8*rxt(78)

         mat(324) = -( het_rates(71) )

         mat(120) = -( rxt(79) + het_rates(72) )

         mat(57) = -( het_rates(91) )

         mat(62) = -( het_rates(92) )

         mat(233) = -( het_rates(93) )

         mat(189) = -( rxt(80) + het_rates(94) )

         mat(79) = -( het_rates(95) )

         mat(460) = -( het_rates(97) )

         mat(200) = -( rxt(81) + het_rates(98) )

         mat(241) = -( rxt(76) + het_rates(83) )
         mat(191) = .900_r8*rxt(80)

         mat(308) = -( rxt(77) + het_rates(50) )
         mat(242) = .130_r8*rxt(76)
         mat(192) = .450_r8*rxt(80)

         mat(632) = -( het_rates(85) )

         mat(589) = -( rxt(63) + het_rates(74) )
         mat(250) = .402_r8*rxt(72)
         mat(204) = rxt(81)

         mat(545) = -( rxt(61) + rxt(62) + het_rates(75) )
         mat(248) = .288_r8*rxt(72)
         mat(203) = rxt(81)

         mat(610) = -( het_rates(76) )

         mat(125) = -( het_rates(77) )

         mat(650) = -( het_rates(73) )
         mat(273) = rxt(60) + rxt(299)
         mat(548) = .660_r8*rxt(61)

         mat(365) = -( het_rates(43) )
         mat(197) = rxt(75)

         mat(130) = -( rxt(64) + het_rates(44) )

         mat(288) = -( het_rates(96) )

         mat(50) = -( het_rates(57) )

         mat(397) = -( het_rates(58) )

         mat(157) = -( rxt(66) + het_rates(59) )

         mat(346) = -( rxt(68) + het_rates(60) )
         mat(158) = .820_r8*rxt(66)
         mat(316) = .250_r8*rxt(78)
         mat(201) = .100_r8*rxt(81)

         mat(163) = -( rxt(67) + het_rates(66) )

         mat(228) = -( het_rates(15) )

         mat(100) = -( het_rates(48) )

         mat(383) = -( rxt(74) + het_rates(49) )
         mat(105) = rxt(233)

         mat(527) = -( rxt(73) + het_rates(63) )

         mat(281) = -( het_rates(52) )

         mat(104) = -( rxt(233) + rxt(234) + het_rates(53) )
         mat(66) = rxt(65)

         mat(65) = -( rxt(65) + het_rates(54) )

         mat(154) = -( het_rates(84) )

         mat(516) = -( het_rates(64) )

         mat(558) = -( rxt(69) + het_rates(65) )
         mat(244) = .180_r8*rxt(76)
         mat(193) = .450_r8*rxt(80)

         mat(490) = -( het_rates(86) )

         mat(450) = -( rxt(71) + het_rates(87) )

         mat(573) = -( het_rates(88) )

         mat(135) = -( rxt(70) + het_rates(89) )

         mat(247) = -( rxt(72) + het_rates(90) )

         mat(71) = -( het_rates(118) )

         mat(173) = -( het_rates(119) )

         mat(208) = -( rxt(224) + het_rates(120) )

         mat(55) = -( het_rates(121) )

         mat(92) = -( het_rates(122) )

         mat(41) = -( het_rates(123) )

         mat(1) = -( het_rates(128) )

         mat(2) = -( het_rates(126) )

         mat(3) = -( het_rates(127) )

         mat(4) = -( het_rates(129) )

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

         mat(20) = -( het_rates(145) )

         mat(21) = -( het_rates(146) )

         mat(22) = -( het_rates(147) )

         mat(23) = -( het_rates(148) )

         mat(24) = -( het_rates(149) )

         mat(25) = -( het_rates(150) )

         mat(26) = -( het_rates(151) )

         mat(27) = -( het_rates(152) )

         mat(28) = -( het_rates(153) )

         mat(29) = -( het_rates(154) )

         mat(30) = -( het_rates(155) )

         mat(31) = -( het_rates(156) )

         mat(32) = -( het_rates(157) )

         mat(33) = -( het_rates(158) )

         mat(34) = -( het_rates(159) )


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
