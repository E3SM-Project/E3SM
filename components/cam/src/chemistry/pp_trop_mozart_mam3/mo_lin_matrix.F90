




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

         mat(657) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(494) = rxt(43)

         mat(492) = -( rxt(43) + het_rates(2) )
         mat(650) = rxt(3)
         mat(717) = rxt(5)
         mat(88) = rxt(7)
         mat(617) = rxt(9)
         mat(398) = rxt(45) + rxt(46)

         mat(397) = -( rxt(45) + rxt(46) + rxt(48) + rxt(50)*y(4) + rxt(51)*y(4) &
                      + rxt(52)*y(16) + rxt(53)*y(16) + rxt(54)*y(16) + het_rates(3) )
         mat(642) = rxt(2)

         mat(928) = -( het_rates(17) )
         mat(883) = rxt(14) + rxt(15)
         mat(413) = rxt(17)
         mat(463) = 1.340_r8*rxt(23)
         mat(517) = .700_r8*rxt(24)
         mat(471) = rxt(30)
         mat(422) = rxt(32)
         mat(331) = rxt(35)
         mat(221) = .450_r8*rxt(37)
         mat(264) = 2.000_r8*rxt(38)

         mat(256) = -( het_rates(11) )
         mat(871) = rxt(15)
         mat(395) = rxt(54)*y(16)

         mat(61) = -( het_rates(86) )

         mat(34) = -( het_rates(87) )

         mat(360) = -( rxt(56) + het_rates(15) )
         mat(141) = rxt(13)
         mat(396) = rxt(53)*y(16)

         mat(700) = -( het_rates(5) )
         mat(722) = rxt(5) + .500_r8*rxt(223)
         mat(90) = rxt(7)
         mat(625) = rxt(10)
         mat(400) = 2.000_r8*rxt(50)*y(4)

         mat(723) = -( rxt(5) + rxt(223) + het_rates(6) )
         mat(91) = rxt(6) + rxt(82)
         mat(194) = rxt(8)
         mat(626) = rxt(9)
         mat(132) = rxt(12) + rxt(91)
         mat(188) = .600_r8*rxt(20) + rxt(134)
         mat(226) = rxt(21) + rxt(180)
         mat(418) = rxt(32)

         mat(623) = -( rxt(9) + rxt(10) + rxt(222) + het_rates(7) )
         mat(89) = rxt(6) + rxt(7) + rxt(82)
         mat(131) = rxt(11)
         mat(187) = .400_r8*rxt(20)

         mat(192) = -( rxt(8) + het_rates(8) )
         mat(87) = 2.000_r8*rxt(221)
         mat(602) = rxt(222)
         mat(711) = .500_r8*rxt(223)

         mat(130) = -( rxt(11) + rxt(12) + rxt(91) + het_rates(9) )

         mat(86) = -( rxt(6) + rxt(7) + rxt(82) + rxt(221) + het_rates(10) )

         mat(806) = -( rxt(92)*y(16) + het_rates(12) )
         mat(195) = rxt(8)
         mat(133) = rxt(11)
         mat(142) = rxt(13)
         mat(84) = 2.000_r8*rxt(16)
         mat(200) = rxt(18)
         mat(172) = rxt(19)
         mat(117) = rxt(25)
         mat(42) = rxt(26)
         mat(138) = rxt(27)
         mat(95) = rxt(28)
         mat(67) = rxt(31)
         mat(273) = rxt(39)
         mat(108) = rxt(40)
         mat(153) = rxt(41)
         mat(183) = rxt(42)
         mat(401) = 2.000_r8*rxt(48) + rxt(52)*y(16)
         mat(724) = .500_r8*rxt(223)

         mat(865) = -( rxt(229) + het_rates(13) )
         mat(134) = rxt(12) + rxt(91)
         mat(411) = rxt(17)
         mat(201) = rxt(18)
         mat(461) = 1.340_r8*rxt(22) + .660_r8*rxt(23)
         mat(118) = rxt(25)
         mat(139) = rxt(27)
         mat(470) = rxt(30)
         mat(420) = rxt(32)
         mat(236) = rxt(33)
         mat(452) = rxt(34)
         mat(329) = 2.000_r8*rxt(35)
         mat(220) = .560_r8*rxt(37)
         mat(263) = 2.000_r8*rxt(38)
         mat(274) = .900_r8*rxt(39)
         mat(184) = rxt(42)
         mat(365) = rxt(56)
         mat(163) = rxt(106)
         mat(80) = rxt(114) + rxt(115)
         mat(402) = rxt(53)*y(16)

         mat(82) = -( rxt(16) + het_rates(14) )
         mat(814) = .500_r8*rxt(229)

         mat(915) = -( het_rates(18) )
         mat(412) = rxt(17)
         mat(174) = rxt(19)
         mat(191) = .400_r8*rxt(20)
         mat(516) = .300_r8*rxt(24)
         mat(297) = rxt(29)
         mat(404) = rxt(52)*y(16)
         mat(809) = rxt(92)*y(16)

         mat(140) = -( rxt(13) + het_rates(19) )

         mat(881) = -( rxt(14) + rxt(15) + het_rates(20) )
         mat(143) = rxt(13)
         mat(202) = rxt(18)
         mat(462) = 1.340_r8*rxt(22)
         mat(96) = rxt(28)
         mat(421) = rxt(32)
         mat(237) = .690_r8*rxt(33)
         mat(453) = rxt(34)
         mat(330) = rxt(35)
         mat(275) = .100_r8*rxt(39)
         mat(164) = rxt(106)
         mat(81) = 2.000_r8*rxt(115)
         mat(403) = rxt(53)*y(16) + rxt(54)*y(16)

         mat(203) = -( het_rates(21) )

         mat(74) = -( het_rates(22) )

         mat(40) = -( rxt(26) + het_rates(28) )

         mat(123) = -( het_rates(23) )

         mat(78) = -( rxt(114) + rxt(115) + het_rates(24) )
         mat(41) = rxt(26)

         mat(249) = -( het_rates(25) )

         mat(175) = -( het_rates(26) )

         mat(327) = -( rxt(35) + het_rates(27) )
         mat(79) = rxt(114)

         mat(22) = -( het_rates(29) )

         mat(318) = -( het_rates(30) )
         mat(167) = rxt(36)

         mat(114) = -( rxt(25) + het_rates(31) )

         mat(406) = -( rxt(17) + het_rates(32) )
         mat(198) = rxt(18)
         mat(116) = rxt(25)
         mat(271) = .400_r8*rxt(39)
         mat(106) = rxt(40)

         mat(587) = -( het_rates(33) )
         mat(186) = .600_r8*rxt(20) + rxt(134)
         mat(458) = 1.340_r8*rxt(22)
         mat(508) = .300_r8*rxt(24)
         mat(94) = rxt(28)
         mat(295) = rxt(29)
         mat(466) = rxt(30)
         mat(450) = rxt(34)
         mat(168) = rxt(36)
         mat(219) = .130_r8*rxt(37)
         mat(107) = rxt(40)

         mat(170) = -( rxt(19) + het_rates(34) )

         mat(381) = -( het_rates(35) )
         mat(502) = .700_r8*rxt(24)

         mat(25) = -( het_rates(36) )

         mat(334) = -( het_rates(37) )

         mat(135) = -( rxt(27) + het_rates(38) )

         mat(285) = -( het_rates(39) )

         mat(196) = -( rxt(18) + het_rates(40) )

         mat(293) = -( rxt(29) + het_rates(41) )
         mat(136) = .820_r8*rxt(27)
         mat(268) = .250_r8*rxt(39)
         mat(179) = .100_r8*rxt(42)

         mat(438) = -( het_rates(42) )

         mat(92) = -( rxt(28) + het_rates(43) )

         mat(28) = -( het_rates(44) )

         mat(97) = -( het_rates(45) )

         mat(31) = -( het_rates(49) )

         mat(303) = -( het_rates(50) )

         mat(266) = -( rxt(39) + het_rates(51) )

         mat(165) = -( rxt(36) + het_rates(46) )
         mat(265) = .800_r8*rxt(39)

         mat(277) = -( het_rates(47) )

         mat(104) = -( rxt(40) + het_rates(48) )

         mat(345) = -( het_rates(52) )

         mat(545) = -( het_rates(53) )

         mat(230) = -( rxt(33) + het_rates(54) )

         mat(506) = -( rxt(24) + het_rates(55) )
         mat(233) = .402_r8*rxt(33)
         mat(182) = rxt(42)

         mat(454) = -( rxt(22) + rxt(23) + het_rates(56) )
         mat(231) = .288_r8*rxt(33)
         mat(181) = rxt(42)

         mat(525) = -( het_rates(57) )

         mat(109) = -( het_rates(58) )

         mat(562) = -( het_rates(59) )
         mat(224) = rxt(21) + rxt(180)
         mat(457) = .660_r8*rxt(22)

         mat(145) = -( het_rates(60) )

         mat(448) = -( rxt(34) + het_rates(61) )

         mat(465) = -( rxt(30) + het_rates(62) )
         mat(218) = .180_r8*rxt(37)
         mat(152) = .450_r8*rxt(41)

         mat(478) = -( het_rates(63) )

         mat(65) = -( rxt(31) + het_rates(64) )

         mat(238) = -( het_rates(65) )

         mat(368) = -( het_rates(66) )

         mat(178) = -( rxt(42) + het_rates(67) )

         mat(43) = -( het_rates(68) )

         mat(48) = -( het_rates(69) )

         mat(208) = -( het_rates(70) )

         mat(148) = -( rxt(41) + het_rates(71) )

         mat(57) = -( het_rates(72) )

         mat(216) = -( rxt(37) + het_rates(73) )
         mat(150) = .900_r8*rxt(41)

         mat(261) = -( rxt(38) + het_rates(74) )
         mat(217) = .130_r8*rxt(37)
         mat(151) = .450_r8*rxt(41)


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

         mat(185) = -( rxt(20) + rxt(134) + het_rates(75) )

         mat(154) = -( het_rates(76) )

         mat(222) = -( rxt(21) + rxt(180) + het_rates(77) )

         mat(425) = -( het_rates(78) )

         mat(415) = -( rxt(32) + het_rates(79) )

         mat(38) = -( het_rates(80) )

         mat(69) = -( het_rates(81) )

         mat(20) = -( het_rates(82) )

         mat(1) = -( het_rates(83) )

         mat(2) = -( het_rates(93) )

         mat(51) = -( het_rates(88) )

         mat(119) = -( het_rates(89) )

         mat(159) = -( rxt(106) + het_rates(90) )

         mat(3) = -( het_rates(91) )

         mat(4) = -( het_rates(92) )

         mat(5) = -( het_rates(94) )

         mat(6) = -( het_rates(95) )

         mat(7) = -( het_rates(96) )

         mat(8) = -( het_rates(97) )

         mat(9) = -( het_rates(98) )

         mat(10) = -( het_rates(99) )

         mat(11) = -( het_rates(100) )

         mat(12) = -( het_rates(101) )

         mat(13) = -( het_rates(102) )

         mat(14) = -( het_rates(103) )

         mat(15) = -( het_rates(104) )

         mat(16) = -( het_rates(105) )

         mat(17) = -( het_rates(106) )

         mat(18) = -( het_rates(107) )

         mat(19) = -( het_rates(108) )


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
