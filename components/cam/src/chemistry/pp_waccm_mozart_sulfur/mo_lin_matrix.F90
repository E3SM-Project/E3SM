




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

         mat(419) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(347) = -( rxt(54) + rxt(55) + rxt(56) + rxt(67) + rxt(68) + rxt(69) &
                 + het_rates(2) )
         mat(578) = rxt(1) + 2.000_r8*rxt(2) + rxt(60) + rxt(61) + rxt(62) &
                      + 2.000_r8*rxt(65) + rxt(72) + rxt(73) + rxt(74) + 2.000_r8*rxt(77)
         mat(416) = rxt(4)
         mat(309) = rxt(6)
         mat(372) = rxt(8)
         mat(26) = rxt(10)
         mat(217) = rxt(12)
         mat(465) = rxt(21)
         mat(288) = rxt(24)
         mat(32) = rxt(25)
         mat(395) = rxt(32)
         mat(201) = rxt(50)
         mat(20) = rxt(51)
         mat(232) = rxt(53)
         mat(484) = rxt(93)

         mat(490) = -( rxt(93) + rxt(97)*y(7) + rxt(98)*y(7) + rxt(100)*y(43) &
                      + rxt(101)*y(44) + rxt(102)*y(45) + rxt(103)*y(46) + rxt(104)*y(47) &
                      + rxt(105)*y(42) + rxt(106)*y(50) + rxt(107)*y(49) + rxt(108)*y(15) &
                      + rxt(109)*y(15) + rxt(110)*y(15) + rxt(111)*y(20) + het_rates(3) )
         mat(584) = rxt(1)
         mat(422) = rxt(3)
         mat(471) = rxt(20)

         mat(588) = -( rxt(1) + rxt(2) + rxt(58) + rxt(60) + rxt(61) + rxt(62) + rxt(65) &
                      + rxt(70) + rxt(72) + rxt(73) + rxt(74) + rxt(77) + het_rates(4) )
         mat(426) = rxt(4)
         mat(223) = rxt(13)
         mat(8) = rxt(88)
         mat(5) = rxt(91) + rxt(92)
         mat(494) = rxt(98)*y(7)

         mat(7) = -( rxt(85) + rxt(88) + rxt(87)*y(51) + het_rates(5) )

         mat(4) = -( rxt(91) + rxt(92) + het_rates(6) )
         mat(406) = rxt(3)
         mat(6) = rxt(85) + rxt(87)*y(51)

         mat(162) = -( rxt(57) + het_rates(8) )
         mat(302) = rxt(6)
         mat(70) = rxt(241)

         mat(308) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(371) = rxt(8)
         mat(25) = rxt(10)
         mat(216) = rxt(13)
         mat(135) = rxt(250)
         mat(483) = 2.000_r8*rxt(97)*y(7)

         mat(373) = -( rxt(8) + het_rates(10) )
         mat(27) = rxt(9) + rxt(126)
         mat(119) = rxt(11)
         mat(218) = rxt(12)
         mat(55) = rxt(15) + rxt(135)
         mat(194) = rxt(30)
         mat(84) = rxt(35)

         mat(450) = -( rxt(136)*y(15) + rxt(143)*y(19) + rxt(144)*y(19) + rxt(155)*y(20) &
                      + rxt(204)*y(41) + rxt(205)*y(48) + rxt(206)*y(46) + rxt(207)*y(42) &
                 + het_rates(22) )
         mat(120) = rxt(11)
         mat(56) = rxt(14)
         mat(43) = rxt(16)
         mat(469) = rxt(19)
         mat(89) = 2.000_r8*rxt(22)
         mat(183) = rxt(27)
         mat(128) = rxt(33)
         mat(488) = rxt(108)*y(15) + rxt(111)*y(20)

         mat(215) = -( rxt(12) + rxt(13) + het_rates(11) )
         mat(24) = rxt(9) + rxt(10) + rxt(126)
         mat(54) = rxt(14)
         mat(190) = rxt(29)
         mat(81) = rxt(34)

         mat(117) = -( rxt(11) + het_rates(12) )
         mat(23) = 2.000_r8*rxt(223) + 2.000_r8*rxt(229) + 2.000_r8*rxt(234)
         mat(187) = rxt(224) + rxt(230) + rxt(235)
         mat(79) = rxt(225) + rxt(233) + rxt(236)

         mat(53) = -( rxt(14) + rxt(15) + rxt(135) + het_rates(13) )

         mat(22) = -( rxt(9) + rxt(10) + rxt(126) + rxt(223) + rxt(229) + rxt(234) &
                 + het_rates(14) )

         mat(169) = -( het_rates(16) )
         mat(478) = rxt(108)*y(15)
         mat(435) = rxt(136)*y(15)
         mat(540) = rxt(167)*y(15)

         mat(40) = -( rxt(16) + het_rates(17) )

         mat(510) = -( rxt(17) + rxt(18) + het_rates(18) )
         mat(45) = rxt(16)
         mat(491) = rxt(109)*y(15) + rxt(110)*y(15)

         mat(206) = -( het_rates(21) )
         mat(42) = rxt(16)
         mat(497) = 2.000_r8*rxt(17)
         mat(460) = rxt(19) + 2.000_r8*rxt(21)
         mat(258) = rxt(28)
         mat(479) = rxt(109)*y(15) + rxt(111)*y(20)
         mat(439) = rxt(144)*y(19) + rxt(155)*y(20)
         mat(543) = rxt(162)*y(20)

         mat(535) = -( het_rates(23) )
         mat(58) = rxt(15) + rxt(135)
         mat(492) = rxt(109)*y(15)
         mat(454) = rxt(143)*y(19) + rxt(204)*y(41) + rxt(207)*y(42)
         mat(556) = rxt(203)*y(41)

         mat(86) = -( rxt(22) + het_rates(24) )

         mat(470) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(66) )
         mat(11) = rxt(49)
         mat(451) = rxt(136)*y(15) + rxt(155)*y(20) + rxt(204)*y(41) + rxt(205)*y(48) &
                      + rxt(206)*y(46) + rxt(207)*y(42)

         mat(557) = -( rxt(162)*y(20) + rxt(167)*y(15) + rxt(203)*y(41) + het_rates(27) )
         mat(13) = 2.000_r8*rxt(23)
         mat(297) = rxt(24)
         mat(3) = 2.000_r8*rxt(26)
         mat(185) = rxt(27)
         mat(273) = rxt(28)
         mat(197) = rxt(29)
         mat(16) = rxt(31)
         mat(493) = 3.000_r8*rxt(100)*y(43) + 2.000_r8*rxt(101)*y(44) &
                      + 3.000_r8*rxt(102)*y(45) + rxt(103)*y(46) + 4.000_r8*rxt(104)*y(47)
         mat(455) = rxt(204)*y(41) + 3.000_r8*rxt(205)*y(48) + rxt(206)*y(46)

         mat(12) = -( rxt(23) + het_rates(28) )

         mat(286) = -( rxt(24) + het_rates(29) )
         mat(31) = rxt(25)
         mat(192) = rxt(30)
         mat(2) = 2.000_r8*rxt(178)

         mat(28) = -( rxt(25) + het_rates(30) )

         mat(1) = -( rxt(26) + rxt(178) + het_rates(31) )

         mat(261) = -( rxt(28) + het_rates(32) )
         mat(545) = rxt(162)*y(20) + rxt(167)*y(15) + 2.000_r8*rxt(203)*y(41)

         mat(179) = -( rxt(27) + het_rates(33) )
         mat(188) = rxt(224) + rxt(230) + rxt(235)

         mat(189) = -( rxt(29) + rxt(30) + rxt(224) + rxt(230) + rxt(235) + het_rates(34) &
       )

         mat(14) = -( rxt(31) + het_rates(35) )

         mat(242) = -( het_rates(36) )
         mat(15) = rxt(31)
         mat(391) = rxt(32)
         mat(124) = rxt(33)
         mat(82) = rxt(34)
         mat(480) = rxt(105)*y(42) + rxt(106)*y(50) + rxt(107)*y(49)
         mat(442) = rxt(207)*y(42)

         mat(397) = -( rxt(32) + het_rates(37) )
         mat(85) = rxt(35)

         mat(104) = -( het_rates(38) )

         mat(123) = -( rxt(33) + het_rates(39) )
         mat(80) = rxt(225) + rxt(233) + rxt(236)

         mat(78) = -( rxt(34) + rxt(35) + rxt(225) + rxt(233) + rxt(236) + het_rates(40) &
       )

         mat(95) = -( het_rates(52) )

         mat(131) = -( rxt(250) + het_rates(53) )
         mat(568) = rxt(58) + rxt(70)
         mat(68) = rxt(243)*y(51)

         mat(60) = -( het_rates(54) )
         mat(157) = rxt(57)

         mat(67) = -( rxt(241) + rxt(243)*y(51) + het_rates(55) )
         mat(324) = rxt(54) + rxt(55) + rxt(56) + rxt(67) + rxt(68) + rxt(69)
         mat(564) = rxt(60) + rxt(61) + rxt(62) + rxt(72) + rxt(73) + rxt(74)

         mat(140) = -( het_rates(56) )
         mat(300) = rxt(7)
         mat(69) = rxt(241)
         mat(132) = rxt(250)

         mat(73) = -( het_rates(58) )

         mat(151) = -( het_rates(57) )
         mat(301) = rxt(7)
         mat(335) = rxt(54) + rxt(55) + rxt(56) + rxt(67) + rxt(68) + rxt(69)
         mat(161) = rxt(57)
         mat(570) = rxt(58) + rxt(60) + rxt(61) + rxt(62) + rxt(70) + rxt(72) + rxt(73) &
                      + rxt(74)

         mat(33) = -( rxt(52) + het_rates(59) )

         mat(110) = -( het_rates(60) )
         mat(34) = rxt(52)
         mat(225) = rxt(53)

         mat(228) = -( rxt(53) + het_rates(61) )
         mat(200) = rxt(50)

         mat(199) = -( rxt(50) + het_rates(62) )
         mat(19) = rxt(51)

         mat(18) = -( rxt(51) + het_rates(63) )
         mat(10) = rxt(49)

         mat(47) = -( het_rates(64) )

         mat(9) = -( rxt(49) + het_rates(65) )


      end subroutine linmat01

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

      end subroutine linmat

      end module mo_lin_matrix
