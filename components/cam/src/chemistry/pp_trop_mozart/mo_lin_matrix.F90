




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

         mat(716) = -( rxt(2) + rxt(3) + het_rates(1) )
         mat(798) = .890_r8*rxt(8) + .890_r8*rxt(9)
         mat(223) = rxt(41)

         mat(63) = -( rxt(43) + rxt(44) + rxt(45) + rxt(46)*y(11) + rxt(57)*y(4) &
                      + rxt(58)*y(4) + rxt(73)*y(15) + het_rates(3) )
         mat(691) = rxt(2)

         mat(219) = -( rxt(41) + het_rates(2) )
         mat(696) = rxt(3)
         mat(756) = rxt(5)
         mat(64) = rxt(43) + rxt(44)

         mat(558) = -( het_rates(5) )
         mat(762) = rxt(5) + .500_r8*rxt(202)
         mat(795) = .110_r8*rxt(8) + .110_r8*rxt(9)
         mat(66) = 2.000_r8*rxt(58)*y(4)

         mat(767) = -( rxt(5) + rxt(202) + het_rates(6) )
         mat(72) = rxt(6) + rxt(65)
         mat(217) = rxt(7)
         mat(800) = .890_r8*rxt(8) + .890_r8*rxt(9)
         mat(124) = .660_r8*rxt(10) + .660_r8*rxt(11) + rxt(71)
         mat(212) = .600_r8*rxt(19) + rxt(109)
         mat(231) = rxt(20) + rxt(176)
         mat(347) = rxt(30)

         mat(801) = -( rxt(8) + rxt(9) + rxt(201) + het_rates(7) )
         mat(73) = rxt(6) + rxt(65)
         mat(125) = .330_r8*rxt(10) + .330_r8*rxt(11)
         mat(213) = .400_r8*rxt(19)

         mat(215) = -( rxt(7) + het_rates(8) )
         mat(71) = 2.000_r8*rxt(200)
         mat(776) = rxt(201)
         mat(755) = .500_r8*rxt(202)

         mat(121) = -( rxt(10) + rxt(11) + rxt(71) + het_rates(9) )

         mat(70) = -( rxt(6) + rxt(65) + rxt(200) + het_rates(10) )

         mat(634) = -( rxt(47)*y(11) + rxt(72)*y(15) + rxt(81)*y(16) + rxt(82)*y(16) &
                      + rxt(211)*y(99) + rxt(212)*y(100) + het_rates(12) )
         mat(216) = rxt(7)
         mat(122) = .330_r8*rxt(10) + .330_r8*rxt(11)
         mat(133) = rxt(12)
         mat(50) = 2.000_r8*rxt(15)
         mat(181) = rxt(17)
         mat(165) = rxt(18)
         mat(425) = .330_r8*rxt(21) + .330_r8*rxt(22)
         mat(129) = rxt(24)
         mat(119) = rxt(25)
         mat(87) = rxt(26)
         mat(58) = rxt(29)
         mat(267) = rxt(37)
         mat(99) = rxt(38)
         mat(151) = rxt(39)
         mat(188) = rxt(40)
         mat(67) = 2.000_r8*rxt(45) + rxt(46)*y(11) + .750_r8*rxt(73)*y(15)
         mat(763) = .500_r8*rxt(202)

         mat(685) = -( rxt(210) + het_rates(13) )
         mat(123) = .660_r8*rxt(10) + .660_r8*rxt(11) + rxt(71)
         mat(134) = rxt(12)
         mat(441) = 2.000_r8*rxt(13)
         mat(385) = rxt(16)
         mat(182) = rxt(17)
         mat(426) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(130) = rxt(24)
         mat(120) = rxt(25)
         mat(432) = rxt(28)
         mat(346) = rxt(30)
         mat(240) = rxt(31)
         mat(418) = rxt(32)
         mat(249) = 2.000_r8*rxt(33)
         mat(206) = .560_r8*rxt(35)
         mat(192) = 2.000_r8*rxt(36)
         mat(268) = .900_r8*rxt(37)
         mat(189) = rxt(40)
         mat(161) = rxt(86)
         mat(55) = rxt(93) + rxt(94)
         mat(68) = rxt(46)*y(11) + .400_r8*rxt(73)*y(15)
         mat(635) = rxt(47)*y(11) + rxt(81)*y(16) + rxt(82)*y(16) + rxt(211)*y(99) &
                      + rxt(212)*y(100)

         mat(49) = -( rxt(15) + het_rates(14) )
         mat(641) = .500_r8*rxt(210)

         mat(746) = -( het_rates(17) )
         mat(386) = rxt(16)
         mat(166) = rxt(18)
         mat(211) = .400_r8*rxt(19)
         mat(470) = .300_r8*rxt(23)
         mat(290) = rxt(27)
         mat(637) = rxt(72)*y(15)
         mat(69) = .750_r8*rxt(73)*y(15)

         mat(131) = -( rxt(12) + het_rates(18) )

         mat(438) = -( rxt(13) + rxt(14) + het_rates(19) )
         mat(132) = rxt(12)
         mat(180) = rxt(17)
         mat(422) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(86) = rxt(26)
         mat(343) = rxt(30)
         mat(235) = .690_r8*rxt(31)
         mat(416) = rxt(32)
         mat(247) = rxt(33)
         mat(266) = .100_r8*rxt(37)
         mat(158) = rxt(86)
         mat(54) = 2.000_r8*rxt(94)
         mat(65) = .250_r8*rxt(73)*y(15)

         mat(241) = -( het_rates(20) )

         mat(80) = -( het_rates(21) )

         mat(106) = -( het_rates(22) )

         mat(52) = -( rxt(93) + rxt(94) + het_rates(23) )

         mat(140) = -( het_rates(24) )

         mat(168) = -( het_rates(25) )

         mat(246) = -( rxt(33) + het_rates(26) )
         mat(53) = rxt(93)

         mat(21) = -( het_rates(27) )

         mat(321) = -( het_rates(28) )
         mat(173) = rxt(34)

         mat(126) = -( rxt(24) + het_rates(29) )

         mat(383) = -( rxt(16) + het_rates(30) )
         mat(178) = rxt(17)
         mat(128) = rxt(24)
         mat(265) = .400_r8*rxt(37)
         mat(98) = rxt(38)

         mat(824) = -( het_rates(31) )
         mat(214) = .600_r8*rxt(19) + rxt(109)
         mat(428) = .670_r8*rxt(21) + .670_r8*rxt(22)
         mat(473) = .300_r8*rxt(23)
         mat(88) = rxt(26)
         mat(291) = rxt(27)
         mat(435) = rxt(28)
         mat(419) = rxt(32)
         mat(175) = rxt(34)
         mat(207) = .130_r8*rxt(35)
         mat(100) = rxt(38)

         mat(163) = -( rxt(18) + het_rates(32) )

         mat(370) = -( het_rates(33) )
         mat(459) = .700_r8*rxt(23)

         mat(24) = -( het_rates(34) )

         mat(331) = -( het_rates(35) )

         mat(116) = -( rxt(25) + het_rates(36) )

         mat(279) = -( het_rates(37) )

         mat(176) = -( rxt(17) + het_rates(38) )

         mat(287) = -( rxt(27) + het_rates(39) )
         mat(117) = .820_r8*rxt(25)
         mat(262) = .250_r8*rxt(37)
         mat(184) = .100_r8*rxt(40)

         mat(404) = -( het_rates(40) )

         mat(84) = -( rxt(26) + het_rates(41) )

         mat(27) = -( het_rates(42) )

         mat(89) = -( het_rates(43) )

         mat(30) = -( het_rates(47) )

         mat(306) = -( het_rates(48) )

         mat(260) = -( rxt(37) + het_rates(49) )

         mat(171) = -( rxt(34) + het_rates(44) )
         mat(259) = .800_r8*rxt(37)

         mat(271) = -( het_rates(45) )

         mat(96) = -( rxt(38) + het_rates(46) )

         mat(352) = -( het_rates(50) )

         mat(501) = -( het_rates(51) )

         mat(233) = -( rxt(31) + het_rates(52) )

         mat(464) = -( rxt(23) + het_rates(53) )
         mat(237) = .402_r8*rxt(31)
         mat(187) = rxt(40)

         mat(420) = -( rxt(21) + rxt(22) + het_rates(54) )
         mat(234) = .288_r8*rxt(31)
         mat(186) = rxt(40)

         mat(482) = -( het_rates(55) )

         mat(101) = -( het_rates(56) )

         mat(517) = -( het_rates(57) )
         mat(228) = rxt(20) + rxt(176)
         mat(424) = .330_r8*rxt(21) + .330_r8*rxt(22)

         mat(136) = -( het_rates(58) )

         mat(414) = -( rxt(32) + het_rates(59) )

         mat(430) = -( rxt(28) + het_rates(60) )
         mat(204) = .180_r8*rxt(35)
         mat(150) = .450_r8*rxt(39)

         mat(451) = -( het_rates(61) )

         mat(56) = -( rxt(29) + het_rates(62) )

         mat(250) = -( het_rates(63) )

         mat(392) = -( het_rates(64) )

         mat(183) = -( rxt(40) + het_rates(65) )

         mat(36) = -( het_rates(66) )

         mat(41) = -( het_rates(67) )

         mat(195) = -( het_rates(68) )

         mat(146) = -( rxt(39) + het_rates(69) )

         mat(59) = -( het_rates(70) )

         mat(203) = -( rxt(35) + het_rates(71) )
         mat(149) = .900_r8*rxt(39)

         mat(190) = -( rxt(36) + het_rates(72) )
         mat(202) = .130_r8*rxt(35)
         mat(147) = .450_r8*rxt(39)

         mat(208) = -( rxt(19) + rxt(109) + het_rates(73) )

         mat(152) = -( het_rates(74) )

         mat(225) = -( rxt(20) + rxt(176) + het_rates(75) )

         mat(292) = -( het_rates(76) )

         mat(342) = -( rxt(30) + het_rates(77) )

         mat(34) = -( het_rates(83) )

         mat(75) = -( het_rates(84) )

         mat(1) = -( het_rates(85) )

         mat(19) = -( het_rates(86) )

         mat(2) = -( het_rates(87) )

         mat(3) = -( het_rates(88) )

         mat(4) = -( het_rates(82) )

         mat(5) = -( rxt(203) + het_rates(78) )

         mat(7) = -( het_rates(79) )
         mat(6) = rxt(203)

         mat(8) = -( rxt(209) + het_rates(80) )

         mat(10) = -( het_rates(81) )
         mat(9) = rxt(209)


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

         mat(44) = -( het_rates(101) )

         mat(113) = -( het_rates(102) )

         mat(157) = -( rxt(86) + het_rates(103) )

         mat(11) = -( het_rates(89) )

         mat(12) = -( het_rates(90) )

         mat(13) = -( het_rates(91) )

         mat(14) = -( het_rates(92) )

         mat(15) = -( het_rates(93) )

         mat(16) = -( het_rates(94) )

         mat(17) = -( het_rates(95) )

         mat(18) = -( het_rates(96) )


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
