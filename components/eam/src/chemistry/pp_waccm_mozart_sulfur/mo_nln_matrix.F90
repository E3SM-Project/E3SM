




      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: nlnmat

      contains

      subroutine nlnmat01( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)


!----------------------------------------------
! ... local variables
!----------------------------------------------

!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------


         mat(419) = -(rxt(81)*y(2) + rxt(99)*y(3) + rxt(121)*y(9) + rxt(124)*y(10) &
                      + rxt(146)*y(21) + rxt(151)*y(22) + rxt(158)*y(23) + rxt(161) &
                      *y(27) + rxt(187)*y(36) + rxt(212)*y(60) + rxt(215)*y(61))
         mat(350) = -rxt(81)*y(1)
         mat(487) = -rxt(99)*y(1)
         mat(312) = -rxt(121)*y(1)
         mat(375) = -rxt(124)*y(1)
         mat(208) = -rxt(146)*y(1)
         mat(449) = -rxt(151)*y(1)
         mat(530) = -rxt(158)*y(1)
         mat(551) = -rxt(161)*y(1)
         mat(245) = -rxt(187)*y(1)
         mat(114) = -rxt(212)*y(1)
         mat(235) = -rxt(215)*y(1)

         mat(350) = mat(350) + rxt(80)*y(4)
         mat(581) = rxt(80)*y(2)

         mat(347) = -(rxt(80)*y(4) + rxt(81)*y(1) + 4._r8*rxt(82)*y(2) + rxt(119)*y(9) &
                      + (rxt(122) + rxt(123)) * y(10) + rxt(130)*y(11) + rxt(142) &
                      *y(18) + rxt(150)*y(22) + rxt(157)*y(23) + rxt(160)*y(24) &
                      + rxt(168)*y(29) + rxt(180)*y(32) + rxt(181)*y(33) + rxt(184) &
                      *y(34) + rxt(190)*y(37) + rxt(200)*y(38) + rxt(201)*y(39) &
                      + rxt(202)*y(40) + rxt(208)*y(59) + (rxt(242) + rxt(251) &
                      ) * y(52) + rxt(248)*y(54))
         mat(578) = -rxt(80)*y(2)
         mat(416) = -rxt(81)*y(2)
         mat(309) = -rxt(119)*y(2)
         mat(372) = -(rxt(122) + rxt(123)) * y(2)
         mat(217) = -rxt(130)*y(2)
         mat(503) = -rxt(142)*y(2)
         mat(446) = -rxt(150)*y(2)
         mat(527) = -rxt(157)*y(2)
         mat(88) = -rxt(160)*y(2)
         mat(288) = -rxt(168)*y(2)
         mat(264) = -rxt(180)*y(2)
         mat(182) = -rxt(181)*y(2)
         mat(193) = -rxt(184)*y(2)
         mat(395) = -rxt(190)*y(2)
         mat(106) = -rxt(200)*y(2)
         mat(126) = -rxt(201)*y(2)
         mat(83) = -rxt(202)*y(2)
         mat(38) = -rxt(208)*y(2)
         mat(101) = -(rxt(242) + rxt(251)) * y(2)
         mat(65) = -rxt(248)*y(2)

         mat(484) = (rxt(94)+rxt(95))*y(4)
         mat(578) = mat(578) + (rxt(94)+rxt(95))*y(3) + rxt(116)*y(8) + rxt(247)*y(54) &
                      + rxt(240)*y(55) + rxt(211)*y(60) + rxt(214)*y(61)
         mat(164) = rxt(116)*y(4) + rxt(117)*y(9) + rxt(118)*y(10) + rxt(244)*y(53)
         mat(309) = mat(309) + rxt(117)*y(8)
         mat(372) = mat(372) + rxt(118)*y(8)
         mat(446) = mat(446) + 2.000_r8*rxt(153)*y(22)
         mat(207) = rxt(149)*y(23)
         mat(527) = mat(527) + rxt(149)*y(21)
         mat(136) = rxt(244)*y(8) + 1.150_r8*rxt(253)*y(57)
         mat(65) = mat(65) + rxt(247)*y(4)
         mat(71) = rxt(240)*y(4)
         mat(144) = rxt(252)*y(57)
         mat(154) = 1.150_r8*rxt(253)*y(53) + rxt(252)*y(56)
         mat(113) = rxt(211)*y(4)
         mat(232) = rxt(214)*y(4)

         mat(490) = -((rxt(94) + rxt(95)) * y(4) + rxt(96)*y(66) + rxt(99)*y(1) &
                      + rxt(112)*y(32) + rxt(113)*y(38))
         mat(584) = -(rxt(94) + rxt(95)) * y(3)
         mat(471) = -rxt(96)*y(3)
         mat(422) = -rxt(99)*y(3)
         mat(270) = -rxt(112)*y(3)
         mat(109) = -rxt(113)*y(3)

         mat(584) = mat(584) + rxt(114)*y(58)
         mat(137) = .850_r8*rxt(253)*y(57)
         mat(76) = rxt(114)*y(4)
         mat(155) = .850_r8*rxt(253)*y(53)

         mat(588) = -(rxt(80)*y(2) + rxt(90)*y(6) + rxt(94)*y(3) + rxt(114)*y(58) &
                      + rxt(116)*y(8) + rxt(145)*y(21) + rxt(211)*y(60) + rxt(214) &
                      *y(61) + rxt(221)*y(64) + rxt(240)*y(55) + (rxt(246) + rxt(247) &
                      ) * y(54) + rxt(249)*y(52))
         mat(357) = -rxt(80)*y(4)
         mat(5) = -rxt(90)*y(4)
         mat(494) = -rxt(94)*y(4)
         mat(77) = -rxt(114)*y(4)
         mat(167) = -rxt(116)*y(4)
         mat(212) = -rxt(145)*y(4)
         mat(116) = -rxt(211)*y(4)
         mat(240) = -rxt(214)*y(4)
         mat(52) = -rxt(221)*y(4)
         mat(72) = -rxt(240)*y(4)
         mat(66) = -(rxt(246) + rxt(247)) * y(4)
         mat(103) = -rxt(249)*y(4)

         mat(426) = 2.000_r8*rxt(81)*y(2) + 2.000_r8*rxt(99)*y(3) + rxt(121)*y(9) &
                      + rxt(124)*y(10) + rxt(151)*y(22) + rxt(146)*y(21) &
                      + 2.000_r8*rxt(158)*y(23) + rxt(161)*y(27) + rxt(187)*y(36) &
                      + rxt(212)*y(60) + rxt(215)*y(61)
         mat(357) = mat(357) + 2.000_r8*rxt(81)*y(1) + 2.000_r8*rxt(82)*y(2) + rxt(89) &
                      *y(6) + rxt(122)*y(10) + rxt(150)*y(22) + rxt(130)*y(11) &
                      + rxt(157)*y(23) + rxt(168)*y(29) + rxt(190)*y(37)
         mat(494) = mat(494) + 2.000_r8*rxt(99)*y(1)
         mat(588) = mat(588) + 2.000_r8*rxt(90)*y(6)
         mat(5) = mat(5) + rxt(89)*y(2) + 2.000_r8*rxt(90)*y(4)
         mat(319) = rxt(121)*y(1) + rxt(245)*y(53)
         mat(382) = rxt(124)*y(1) + rxt(122)*y(2)
         mat(456) = rxt(151)*y(1) + rxt(150)*y(2) + rxt(134)*y(13) + rxt(152)*y(23) &
                      + rxt(170)*y(29)
         mat(223) = rxt(130)*y(2) + rxt(132)*y(23)
         mat(59) = rxt(134)*y(22)
         mat(177) = rxt(138)*y(23)
         mat(212) = mat(212) + rxt(146)*y(1) + rxt(148)*y(23)
         mat(537) = 2.000_r8*rxt(158)*y(1) + rxt(157)*y(2) + rxt(152)*y(22) + rxt(132) &
                      *y(11) + rxt(138)*y(16) + rxt(148)*y(21) + 2.000_r8*rxt(159) &
                      *y(23) + rxt(164)*y(27) + rxt(171)*y(29) + rxt(188)*y(36) &
                      + rxt(192)*y(37)
         mat(558) = rxt(161)*y(1) + rxt(164)*y(23)
         mat(298) = rxt(168)*y(2) + rxt(170)*y(22) + rxt(171)*y(23) + ( &
                      + 2.000_r8*rxt(174)+2.000_r8*rxt(175))*y(29) + (rxt(196) &
                       +rxt(197))*y(37)
         mat(251) = rxt(187)*y(1) + rxt(188)*y(23)
         mat(405) = rxt(190)*y(2) + rxt(192)*y(23) + (rxt(196)+rxt(197))*y(29) &
                      + 2.000_r8*rxt(198)*y(37)
         mat(138) = rxt(245)*y(9)
         mat(116) = mat(116) + rxt(212)*y(1)
         mat(240) = mat(240) + rxt(215)*y(1)

         mat(7) = -(rxt(83)*y(2) + rxt(84)*y(4) + rxt(86)*y(1))
         mat(321) = -rxt(83)*y(5)
         mat(560) = -rxt(84)*y(5)
         mat(407) = -rxt(86)*y(5)

         mat(476) = rxt(94)*y(4)
         mat(560) = mat(560) + rxt(94)*y(3)

         mat(4) = -(rxt(89)*y(2) + rxt(90)*y(4))
         mat(320) = -rxt(89)*y(6)
         mat(559) = -rxt(90)*y(6)

         mat(406) = rxt(86)*y(5)
         mat(320) = mat(320) + rxt(83)*y(5)
         mat(559) = mat(559) + rxt(84)*y(5)
         mat(6) = rxt(86)*y(1) + rxt(83)*y(2) + rxt(84)*y(4)

         mat(162) = -(rxt(116)*y(4) + rxt(117)*y(9) + rxt(118)*y(10) + rxt(244)*y(53))
         mat(571) = -rxt(116)*y(8)
         mat(302) = -rxt(117)*y(8)
         mat(363) = -rxt(118)*y(8)
         mat(134) = -rxt(244)*y(8)

         mat(336) = rxt(248)*y(54) + rxt(115)*y(58)
         mat(571) = mat(571) + rxt(246)*y(54)
         mat(99) = 1.100_r8*rxt(254)*y(57)
         mat(64) = rxt(248)*y(2) + rxt(246)*y(4)
         mat(142) = .200_r8*rxt(252)*y(57)
         mat(74) = rxt(115)*y(2)
         mat(152) = 1.100_r8*rxt(254)*y(52) + .200_r8*rxt(252)*y(56)

         mat(308) = -(rxt(117)*y(8) + rxt(119)*y(2) + rxt(120)*y(23) + rxt(121)*y(1) &
                      + rxt(129)*y(11) + rxt(137)*y(16) + rxt(172)*y(29) + rxt(193) &
                      *y(37) + rxt(245)*y(53))
         mat(163) = -rxt(117)*y(9)
         mat(346) = -rxt(119)*y(9)
         mat(526) = -rxt(120)*y(9)
         mat(415) = -rxt(121)*y(9)
         mat(216) = -rxt(129)*y(9)
         mat(171) = -rxt(137)*y(9)
         mat(287) = -rxt(172)*y(9)
         mat(394) = -rxt(193)*y(9)
         mat(135) = -rxt(245)*y(9)

         mat(346) = mat(346) + rxt(122)*y(10)
         mat(577) = rxt(116)*y(8) + rxt(114)*y(58)
         mat(163) = mat(163) + rxt(116)*y(4)
         mat(371) = rxt(122)*y(2) + rxt(216)*y(61)
         mat(75) = rxt(114)*y(4)
         mat(231) = rxt(216)*y(10)

         mat(373) = -(rxt(118)*y(8) + (rxt(122) + rxt(123)) * y(2) + rxt(124)*y(1) &
                      + rxt(125)*y(11) + rxt(127)*y(22) + rxt(133)*y(23) + rxt(173) &
                      *y(29) + rxt(194)*y(37) + rxt(216)*y(61))
         mat(165) = -rxt(118)*y(10)
         mat(348) = -(rxt(122) + rxt(123)) * y(10)
         mat(417) = -rxt(124)*y(10)
         mat(218) = -rxt(125)*y(10)
         mat(447) = -rxt(127)*y(10)
         mat(528) = -rxt(133)*y(10)
         mat(289) = -rxt(173)*y(10)
         mat(396) = -rxt(194)*y(10)
         mat(233) = -rxt(216)*y(10)

         mat(417) = mat(417) + rxt(121)*y(9)
         mat(348) = mat(348) + rxt(119)*y(9) + rxt(130)*y(11)
         mat(310) = rxt(121)*y(1) + rxt(119)*y(2) + 2.000_r8*rxt(129)*y(11) + rxt(137) &
                      *y(16) + rxt(120)*y(23) + rxt(172)*y(29) + rxt(193)*y(37)
         mat(447) = mat(447) + rxt(131)*y(11) + rxt(134)*y(13)
         mat(218) = mat(218) + rxt(130)*y(2) + 2.000_r8*rxt(129)*y(9) + rxt(131)*y(22) &
                      + rxt(132)*y(23)
         mat(55) = rxt(134)*y(22)
         mat(172) = rxt(137)*y(9)
         mat(528) = mat(528) + rxt(120)*y(9) + rxt(132)*y(11)
         mat(289) = mat(289) + rxt(172)*y(9)
         mat(396) = mat(396) + rxt(193)*y(9)

         mat(450) = -(rxt(127)*y(10) + rxt(128)*y(12) + rxt(131)*y(11) + rxt(134) &
                      *y(13) + rxt(139)*y(17) + rxt(141)*y(18) + rxt(150)*y(2) + rxt(151) &
                      *y(1) + rxt(152)*y(23) + (4._r8*rxt(153) + 4._r8*rxt(154) &
                      ) * y(22) + rxt(156)*y(24) + (rxt(169) + rxt(170)) * y(29) &
                      + rxt(179)*y(32) + rxt(183)*y(33) + rxt(185)*y(34) + rxt(191) &
                      *y(37) + rxt(199)*y(38) + rxt(209)*y(59) + rxt(210)*y(60) &
                      + rxt(213)*y(61) + rxt(220)*y(62))
         mat(376) = -rxt(127)*y(22)
         mat(120) = -rxt(128)*y(22)
         mat(219) = -rxt(131)*y(22)
         mat(56) = -rxt(134)*y(22)
         mat(43) = -rxt(139)*y(22)
         mat(507) = -rxt(141)*y(22)
         mat(351) = -rxt(150)*y(22)
         mat(420) = -rxt(151)*y(22)
         mat(531) = -rxt(152)*y(22)
         mat(89) = -rxt(156)*y(22)
         mat(292) = -(rxt(169) + rxt(170)) * y(22)
         mat(268) = -rxt(179)*y(22)
         mat(183) = -rxt(183)*y(22)
         mat(195) = -rxt(185)*y(22)
         mat(399) = -rxt(191)*y(22)
         mat(107) = -rxt(199)*y(22)
         mat(39) = -rxt(209)*y(22)
         mat(115) = -rxt(210)*y(22)
         mat(236) = -rxt(213)*y(22)
         mat(202) = -rxt(220)*y(22)

         mat(420) = mat(420) + rxt(146)*y(21) + rxt(158)*y(23)
         mat(351) = mat(351) + rxt(142)*y(18) + rxt(157)*y(23) + rxt(160)*y(24) &
                      + rxt(180)*y(32) + rxt(181)*y(33) + rxt(200)*y(38) + rxt(201) &
                      *y(39)
         mat(488) = 2.000_r8*rxt(96)*y(66) + rxt(112)*y(32) + rxt(113)*y(38)
         mat(313) = rxt(120)*y(23)
         mat(219) = mat(219) + rxt(132)*y(23)
         mat(507) = mat(507) + rxt(142)*y(2)
         mat(209) = rxt(146)*y(1) + 2.000_r8*rxt(147)*y(23)
         mat(531) = mat(531) + rxt(158)*y(1) + rxt(157)*y(2) + rxt(120)*y(9) &
                      + rxt(132)*y(11) + 2.000_r8*rxt(147)*y(21) + rxt(165)*y(27)
         mat(89) = mat(89) + rxt(160)*y(2)
         mat(469) = 2.000_r8*rxt(96)*y(3)
         mat(552) = rxt(165)*y(23)
         mat(268) = mat(268) + rxt(180)*y(2) + rxt(112)*y(3)
         mat(183) = mat(183) + rxt(181)*y(2)
         mat(107) = mat(107) + rxt(200)*y(2) + rxt(113)*y(3)
         mat(128) = rxt(201)*y(2)


      end subroutine nlnmat01

      subroutine nlnmat02( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)


!----------------------------------------------
! ... local variables
!----------------------------------------------

!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------


         mat(215) = -(rxt(125)*y(10) + rxt(129)*y(9) + rxt(130)*y(2) + rxt(131)*y(22) &
                      + rxt(132)*y(23) + rxt(140)*y(18))
         mat(366) = -rxt(125)*y(11)
         mat(305) = -rxt(129)*y(11)
         mat(341) = -rxt(130)*y(11)
         mat(440) = -rxt(131)*y(11)
         mat(522) = -rxt(132)*y(11)
         mat(498) = -rxt(140)*y(11)

         mat(411) = rxt(124)*y(10)
         mat(341) = mat(341) + rxt(123)*y(10) + rxt(184)*y(34) + rxt(202)*y(40)
         mat(366) = mat(366) + rxt(124)*y(1) + rxt(123)*y(2)
         mat(440) = mat(440) + rxt(128)*y(12) + rxt(185)*y(34)
         mat(118) = rxt(128)*y(22)
         mat(544) = rxt(186)*y(34)
         mat(190) = rxt(184)*y(2) + rxt(185)*y(22) + rxt(186)*y(27)
         mat(81) = rxt(202)*y(2)

         mat(117) = -(rxt(128)*y(22))
         mat(434) = -rxt(128)*y(12)

         mat(361) = rxt(127)*y(22)
         mat(434) = mat(434) + rxt(127)*y(10)
         mat(214) = rxt(140)*y(18)
         mat(496) = rxt(140)*y(11)
         mat(254) = (rxt(226)+rxt(231)+rxt(237))*y(34)
         mat(187) = (rxt(226)+rxt(231)+rxt(237))*y(32)

         mat(53) = -(rxt(134)*y(22))
         mat(430) = -rxt(134)*y(13)

         mat(359) = rxt(133)*y(23)
         mat(515) = rxt(133)*y(10)


         mat(358) = rxt(125)*y(11)
         mat(213) = rxt(125)*y(10)

         mat(169) = -(rxt(137)*y(9) + rxt(138)*y(23))
         mat(303) = -rxt(137)*y(16)
         mat(519) = -rxt(138)*y(16)

         mat(435) = rxt(139)*y(17)
         mat(41) = rxt(139)*y(22)

         mat(40) = -(rxt(139)*y(22))
         mat(428) = -rxt(139)*y(17)

         mat(168) = rxt(138)*y(23)
         mat(514) = rxt(138)*y(16)

         mat(510) = -(rxt(140)*y(11) + rxt(141)*y(22) + rxt(142)*y(2) + rxt(166)*y(27) &
                      + rxt(189)*y(36))
         mat(221) = -rxt(140)*y(18)
         mat(453) = -rxt(141)*y(18)
         mat(354) = -rxt(142)*y(18)
         mat(555) = -rxt(166)*y(18)
         mat(249) = -rxt(189)*y(18)

         mat(316) = rxt(137)*y(16)
         mat(175) = rxt(137)*y(9)

         mat(206) = -(rxt(145)*y(4) + rxt(146)*y(1) + (rxt(147) + rxt(148) + rxt(149) &
                      ) * y(23))
         mat(573) = -rxt(145)*y(21)
         mat(410) = -rxt(146)*y(21)
         mat(521) = -(rxt(147) + rxt(148) + rxt(149)) * y(21)

         mat(340) = rxt(150)*y(22)
         mat(439) = rxt(150)*y(2) + rxt(141)*y(18) + rxt(209)*y(59) + rxt(210)*y(60) &
                      + rxt(213)*y(61)
         mat(497) = rxt(141)*y(22)
         mat(36) = rxt(209)*y(22)
         mat(111) = rxt(210)*y(22)
         mat(227) = rxt(213)*y(22)

         mat(535) = -(rxt(120)*y(9) + rxt(132)*y(11) + rxt(133)*y(10) + rxt(138)*y(16) &
                      + (rxt(147) + rxt(148) + rxt(149)) * y(21) + rxt(152)*y(22) &
                      + rxt(157)*y(2) + rxt(158)*y(1) + 4._r8*rxt(159)*y(23) + (rxt(164) &
                      + rxt(165)) * y(27) + rxt(171)*y(29) + rxt(188)*y(36) + rxt(192) &
                      *y(37))
         mat(317) = -rxt(120)*y(23)
         mat(222) = -rxt(132)*y(23)
         mat(380) = -rxt(133)*y(23)
         mat(176) = -rxt(138)*y(23)
         mat(211) = -(rxt(147) + rxt(148) + rxt(149)) * y(23)
         mat(454) = -rxt(152)*y(23)
         mat(355) = -rxt(157)*y(23)
         mat(424) = -rxt(158)*y(23)
         mat(556) = -(rxt(164) + rxt(165)) * y(23)
         mat(296) = -rxt(171)*y(23)
         mat(250) = -rxt(188)*y(23)
         mat(403) = -rxt(192)*y(23)

         mat(424) = mat(424) + rxt(151)*y(22)
         mat(355) = mat(355) + rxt(142)*y(18) + rxt(160)*y(24)
         mat(586) = rxt(145)*y(21) + rxt(221)*y(64)
         mat(317) = mat(317) + rxt(137)*y(16)
         mat(454) = mat(454) + rxt(151)*y(1) + rxt(131)*y(11) + rxt(156)*y(24) &
                      + rxt(169)*y(29) + rxt(191)*y(37)
         mat(222) = mat(222) + rxt(131)*y(22) + rxt(140)*y(18)
         mat(176) = mat(176) + rxt(137)*y(9)
         mat(511) = rxt(142)*y(2) + rxt(140)*y(11) + rxt(166)*y(27) + rxt(189)*y(36)
         mat(211) = mat(211) + rxt(145)*y(4)
         mat(91) = rxt(160)*y(2) + rxt(156)*y(22) + rxt(163)*y(27)
         mat(556) = mat(556) + rxt(166)*y(18) + rxt(163)*y(24)
         mat(296) = mat(296) + rxt(169)*y(22)
         mat(250) = mat(250) + rxt(189)*y(18)
         mat(403) = mat(403) + rxt(191)*y(22)
         mat(51) = rxt(221)*y(4)

         mat(86) = -(rxt(156)*y(22) + rxt(160)*y(2) + rxt(163)*y(27))
         mat(431) = -rxt(156)*y(24)
         mat(327) = -rxt(160)*y(24)
         mat(539) = -rxt(163)*y(24)

         mat(431) = mat(431) + 2.000_r8*rxt(154)*y(22)
         mat(516) = 2.000_r8*rxt(159)*y(23)

         mat(470) = -(rxt(96)*y(3) + rxt(222)*y(63))
         mat(489) = -rxt(96)*y(66)
         mat(21) = -rxt(222)*y(66)

         mat(451) = 2.000_r8*rxt(153)*y(22) + rxt(128)*y(12) + rxt(134)*y(13) &
                      + rxt(139)*y(17) + rxt(141)*y(18) + rxt(152)*y(23) + rxt(156) &
                      *y(24) + rxt(179)*y(32) + rxt(183)*y(33) + rxt(199)*y(38)
         mat(121) = rxt(128)*y(22)
         mat(57) = rxt(134)*y(22)
         mat(44) = rxt(139)*y(22)
         mat(508) = rxt(141)*y(22)
         mat(210) = rxt(149)*y(23)
         mat(532) = rxt(152)*y(22) + rxt(149)*y(21)
         mat(90) = rxt(156)*y(22)
         mat(269) = rxt(179)*y(22) + (rxt(227)+rxt(232)+rxt(238))*y(33) + (rxt(228) &
                       +rxt(239))*y(39)
         mat(184) = rxt(183)*y(22) + (rxt(227)+rxt(232)+rxt(238))*y(32)
         mat(108) = rxt(199)*y(22)
         mat(129) = (rxt(228)+rxt(239))*y(32)

         mat(557) = -(rxt(161)*y(1) + rxt(163)*y(24) + (rxt(164) + rxt(165)) * y(23) &
                      + rxt(166)*y(18) + rxt(182)*y(33) + rxt(186)*y(34))
         mat(425) = -rxt(161)*y(27)
         mat(92) = -rxt(163)*y(27)
         mat(536) = -(rxt(164) + rxt(165)) * y(27)
         mat(512) = -rxt(166)*y(27)
         mat(185) = -rxt(182)*y(27)
         mat(197) = -rxt(186)*y(27)

         mat(356) = rxt(168)*y(29) + rxt(180)*y(32)
         mat(493) = rxt(112)*y(32)
         mat(318) = rxt(172)*y(29)
         mat(455) = rxt(169)*y(29) + rxt(179)*y(32)
         mat(297) = rxt(168)*y(2) + rxt(172)*y(9) + rxt(169)*y(22) + ( &
                      + 4.000_r8*rxt(174)+2.000_r8*rxt(176))*y(29) + rxt(196)*y(37) &
                      + rxt(217)*y(61)
         mat(273) = rxt(180)*y(2) + rxt(112)*y(3) + rxt(179)*y(22)
         mat(404) = rxt(196)*y(29)
         mat(239) = rxt(217)*y(29)


         mat(538) = rxt(186)*y(34)
         mat(276) = 2.000_r8*rxt(175)*y(29)
         mat(252) = (rxt(227)+rxt(232)+rxt(238))*y(33) + (rxt(226)+rxt(231)+rxt(237)) &
                      *y(34)
         mat(178) = (rxt(227)+rxt(232)+rxt(238))*y(32)
         mat(186) = rxt(186)*y(27) + (rxt(226)+rxt(231)+rxt(237))*y(32)

         mat(286) = -(rxt(168)*y(2) + (rxt(169) + rxt(170)) * y(22) + rxt(171)*y(23) &
                      + rxt(172)*y(9) + rxt(173)*y(10) + (4._r8*rxt(174) + 4._r8*rxt(175) &
                      + 4._r8*rxt(176) + 4._r8*rxt(177)) * y(29) + (rxt(195) + rxt(196) &
                      + rxt(197)) * y(37) + rxt(217)*y(61))
         mat(345) = -rxt(168)*y(29)
         mat(444) = -(rxt(169) + rxt(170)) * y(29)
         mat(525) = -rxt(171)*y(29)
         mat(307) = -rxt(172)*y(29)
         mat(370) = -rxt(173)*y(29)
         mat(393) = -(rxt(195) + rxt(196) + rxt(197)) * y(29)
         mat(230) = -rxt(217)*y(29)

         mat(414) = rxt(161)*y(27)
         mat(345) = mat(345) + rxt(181)*y(33) + rxt(184)*y(34)
         mat(444) = mat(444) + rxt(183)*y(33)
         mat(525) = mat(525) + rxt(165)*y(27)
         mat(546) = rxt(161)*y(1) + rxt(165)*y(23) + rxt(182)*y(33)
         mat(31) = rxt(219)*y(61)
         mat(181) = rxt(181)*y(2) + rxt(183)*y(22) + rxt(182)*y(27)
         mat(192) = rxt(184)*y(2)
         mat(230) = mat(230) + rxt(219)*y(30)

         mat(28) = -(rxt(219)*y(61))
         mat(224) = -rxt(219)*y(30)

         mat(278) = 2.000_r8*rxt(176)*y(29) + rxt(195)*y(37)
         mat(384) = rxt(195)*y(29)


         mat(275) = 2.000_r8*rxt(177)*y(29)

         mat(261) = -(rxt(112)*y(3) + rxt(179)*y(22) + rxt(180)*y(2) + (rxt(226) &
                      + rxt(231) + rxt(237)) * y(34) + (rxt(227) + rxt(232) + rxt(238) &
                      ) * y(33) + (rxt(228) + rxt(239)) * y(39))
         mat(481) = -rxt(112)*y(32)
         mat(443) = -rxt(179)*y(32)
         mat(344) = -rxt(180)*y(32)
         mat(191) = -(rxt(226) + rxt(231) + rxt(237)) * y(32)
         mat(180) = -(rxt(227) + rxt(232) + rxt(238)) * y(32)
         mat(125) = -(rxt(228) + rxt(239)) * y(32)

         mat(443) = mat(443) + rxt(170)*y(29)
         mat(500) = rxt(166)*y(27)
         mat(524) = rxt(164)*y(27)
         mat(87) = rxt(163)*y(27)
         mat(545) = rxt(166)*y(18) + rxt(164)*y(23) + rxt(163)*y(24) + rxt(182)*y(33)
         mat(285) = rxt(170)*y(22)
         mat(180) = mat(180) + rxt(182)*y(27)

         mat(179) = -(rxt(181)*y(2) + rxt(182)*y(27) + rxt(183)*y(22) + (rxt(227) &
                      + rxt(232) + rxt(238)) * y(32))
         mat(337) = -rxt(181)*y(33)
         mat(541) = -rxt(182)*y(33)
         mat(436) = -rxt(183)*y(33)
         mat(256) = -(rxt(227) + rxt(232) + rxt(238)) * y(33)

         mat(436) = mat(436) + rxt(185)*y(34)
         mat(520) = rxt(171)*y(29)
         mat(279) = rxt(171)*y(23)
         mat(188) = rxt(185)*y(22)

         mat(189) = -(rxt(184)*y(2) + rxt(185)*y(22) + rxt(186)*y(27) + (rxt(226) &
                      + rxt(231) + rxt(237)) * y(32))
         mat(338) = -rxt(184)*y(34)
         mat(437) = -rxt(185)*y(34)
         mat(542) = -rxt(186)*y(34)
         mat(257) = -(rxt(226) + rxt(231) + rxt(237)) * y(34)

         mat(364) = rxt(173)*y(29)
         mat(280) = rxt(173)*y(10)


         mat(277) = rxt(197)*y(37)
         mat(253) = (rxt(228)+rxt(239))*y(39)
         mat(383) = rxt(197)*y(29)
         mat(122) = (rxt(228)+rxt(239))*y(32)


      end subroutine nlnmat02

      subroutine nlnmat03( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)


!----------------------------------------------
! ... local variables
!----------------------------------------------

!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------


         mat(242) = -(rxt(187)*y(1) + rxt(188)*y(23) + rxt(189)*y(18))
         mat(413) = -rxt(187)*y(36)
         mat(523) = -rxt(188)*y(36)
         mat(499) = -rxt(189)*y(36)

         mat(343) = rxt(190)*y(37) + rxt(200)*y(38)
         mat(480) = rxt(113)*y(38)
         mat(306) = rxt(193)*y(37)
         mat(442) = rxt(191)*y(37) + rxt(199)*y(38)
         mat(284) = (rxt(195)+rxt(196))*y(37)
         mat(391) = rxt(190)*y(2) + rxt(193)*y(9) + rxt(191)*y(22) + (rxt(195) &
                       +rxt(196))*y(29) + 4.000_r8*rxt(198)*y(37) + rxt(218)*y(61)
         mat(105) = rxt(200)*y(2) + rxt(113)*y(3) + rxt(199)*y(22)
         mat(229) = rxt(218)*y(37)

         mat(397) = -(rxt(190)*y(2) + rxt(191)*y(22) + rxt(192)*y(23) + rxt(193)*y(9) &
                      + rxt(194)*y(10) + (rxt(195) + rxt(196) + rxt(197)) * y(29) &
                      + 4._r8*rxt(198)*y(37) + rxt(218)*y(61))
         mat(349) = -rxt(190)*y(37)
         mat(448) = -rxt(191)*y(37)
         mat(529) = -rxt(192)*y(37)
         mat(311) = -rxt(193)*y(37)
         mat(374) = -rxt(194)*y(37)
         mat(290) = -(rxt(195) + rxt(196) + rxt(197)) * y(37)
         mat(234) = -rxt(218)*y(37)

         mat(418) = rxt(187)*y(36)
         mat(349) = mat(349) + rxt(201)*y(39) + rxt(202)*y(40)
         mat(244) = rxt(187)*y(1)
         mat(127) = rxt(201)*y(2)
         mat(85) = rxt(202)*y(2)

         mat(104) = -(rxt(113)*y(3) + rxt(199)*y(22) + rxt(200)*y(2))
         mat(477) = -rxt(113)*y(38)
         mat(432) = -rxt(199)*y(38)
         mat(329) = -rxt(200)*y(38)

         mat(495) = rxt(189)*y(36)
         mat(517) = rxt(188)*y(36)
         mat(241) = rxt(189)*y(18) + rxt(188)*y(23)

         mat(123) = -(rxt(201)*y(2) + (rxt(228) + rxt(239)) * y(32))
         mat(332) = -rxt(201)*y(39)
         mat(255) = -(rxt(228) + rxt(239)) * y(39)

         mat(518) = rxt(192)*y(37)
         mat(387) = rxt(192)*y(23)

         mat(78) = -(rxt(202)*y(2))
         mat(326) = -rxt(202)*y(40)

         mat(360) = rxt(194)*y(37)
         mat(385) = rxt(194)*y(10)

         mat(95) = -((rxt(242) + rxt(251)) * y(2) + rxt(249)*y(4) + rxt(254)*y(57))
         mat(328) = -(rxt(242) + rxt(251)) * y(52)
         mat(566) = -rxt(249)*y(52)
         mat(148) = -rxt(254)*y(52)

         mat(131) = -(rxt(244)*y(8) + rxt(245)*y(9) + rxt(253)*y(57))
         mat(159) = -rxt(244)*y(53)
         mat(299) = -rxt(245)*y(53)
         mat(149) = -rxt(253)*y(53)

         mat(568) = rxt(249)*y(52) + rxt(246)*y(54) + rxt(240)*y(55)
         mat(96) = rxt(249)*y(4)
         mat(62) = rxt(246)*y(4)
         mat(68) = rxt(240)*y(4)

         mat(60) = -((rxt(246) + rxt(247)) * y(4) + rxt(248)*y(2))
         mat(563) = -(rxt(246) + rxt(247)) * y(54)
         mat(323) = -rxt(248)*y(54)

         mat(67) = -(rxt(240)*y(4))
         mat(564) = -rxt(240)*y(55)

         mat(324) = rxt(251)*y(52) + rxt(248)*y(54)
         mat(93) = rxt(251)*y(2)
         mat(61) = rxt(248)*y(2)

         mat(140) = -(rxt(252)*y(57))
         mat(150) = -rxt(252)*y(56)

         mat(334) = rxt(242)*y(52)
         mat(569) = rxt(247)*y(54)
         mat(160) = rxt(244)*y(53)
         mat(300) = rxt(245)*y(53)
         mat(97) = rxt(242)*y(2)
         mat(132) = rxt(244)*y(8) + rxt(245)*y(9)
         mat(63) = rxt(247)*y(4)

         mat(73) = -(rxt(114)*y(4) + rxt(115)*y(2))
         mat(565) = -rxt(114)*y(58)
         mat(325) = -rxt(115)*y(58)

         mat(325) = mat(325) + rxt(242)*y(52)
         mat(94) = rxt(242)*y(2) + .900_r8*rxt(254)*y(57)
         mat(139) = .800_r8*rxt(252)*y(57)
         mat(147) = .900_r8*rxt(254)*y(52) + .800_r8*rxt(252)*y(56)

         mat(151) = -(rxt(252)*y(56) + rxt(253)*y(53) + rxt(254)*y(52))
         mat(141) = -rxt(252)*y(57)
         mat(133) = -rxt(253)*y(57)
         mat(98) = -rxt(254)*y(57)

         mat(33) = -(rxt(208)*y(2) + rxt(209)*y(22))
         mat(322) = -rxt(208)*y(59)
         mat(427) = -rxt(209)*y(59)

         mat(110) = -(rxt(210)*y(22) + rxt(211)*y(4) + rxt(212)*y(1))
         mat(433) = -rxt(210)*y(60)
         mat(567) = -rxt(211)*y(60)
         mat(408) = -rxt(212)*y(60)

         mat(228) = -(rxt(213)*y(22) + rxt(214)*y(4) + rxt(215)*y(1) + rxt(216)*y(10) &
                      + rxt(217)*y(29) + rxt(218)*y(37) + rxt(219)*y(30))
         mat(441) = -rxt(213)*y(61)
         mat(574) = -rxt(214)*y(61)
         mat(412) = -rxt(215)*y(61)
         mat(367) = -rxt(216)*y(61)
         mat(283) = -rxt(217)*y(61)
         mat(390) = -rxt(218)*y(61)
         mat(30) = -rxt(219)*y(61)

         mat(412) = mat(412) + rxt(212)*y(60)
         mat(342) = rxt(208)*y(59)
         mat(574) = mat(574) + rxt(211)*y(60)
         mat(441) = mat(441) + rxt(210)*y(60)
         mat(37) = rxt(208)*y(2)
         mat(112) = rxt(212)*y(1) + rxt(211)*y(4) + rxt(210)*y(22)

         mat(199) = -(rxt(220)*y(22))
         mat(438) = -rxt(220)*y(62)

         mat(409) = rxt(215)*y(61)
         mat(572) = rxt(214)*y(61)
         mat(365) = rxt(216)*y(61)
         mat(438) = mat(438) + rxt(209)*y(59) + rxt(213)*y(61)
         mat(281) = rxt(217)*y(61)
         mat(29) = rxt(219)*y(61)
         mat(388) = rxt(218)*y(61)
         mat(35) = rxt(209)*y(22)
         mat(226) = rxt(215)*y(1) + rxt(214)*y(4) + rxt(216)*y(10) + rxt(213)*y(22) &
                      + rxt(217)*y(29) + rxt(219)*y(30) + rxt(218)*y(37)

         mat(18) = -(rxt(222)*y(66))
         mat(458) = -rxt(222)*y(63)

         mat(561) = rxt(221)*y(64)
         mat(46) = rxt(221)*y(4)

         mat(47) = -(rxt(221)*y(4))
         mat(562) = -rxt(221)*y(64)

         mat(429) = rxt(220)*y(62)
         mat(198) = rxt(220)*y(22)


         mat(457) = rxt(222)*y(63)
         mat(17) = rxt(222)*y(66)


      end subroutine nlnmat03

      subroutine nlnmat_finit( mat, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(inout) :: mat(nzcnt)


!----------------------------------------------
! ... local variables
!----------------------------------------------

!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------


         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = lmat( 3)
         mat( 4) = mat( 4) + lmat( 4)
         mat( 5) = mat( 5) + lmat( 5)
         mat( 6) = mat( 6) + lmat( 6)
         mat( 7) = mat( 7) + lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = lmat( 12)
         mat( 13) = lmat( 13)
         mat( 14) = lmat( 14)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 18) = mat( 18) + lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 22) = lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = lmat( 24)
         mat( 25) = lmat( 25)
         mat( 26) = lmat( 26)
         mat( 27) = lmat( 27)
         mat( 28) = mat( 28) + lmat( 28)
         mat( 31) = mat( 31) + lmat( 31)
         mat( 32) = lmat( 32)
         mat( 33) = mat( 33) + lmat( 33)
         mat( 34) = lmat( 34)
         mat( 40) = mat( 40) + lmat( 40)
         mat( 42) = lmat( 42)
         mat( 43) = mat( 43) + lmat( 43)
         mat( 45) = lmat( 45)
         mat( 47) = mat( 47) + lmat( 47)
         mat( 53) = mat( 53) + lmat( 53)
         mat( 54) = lmat( 54)
         mat( 55) = mat( 55) + lmat( 55)
         mat( 56) = mat( 56) + lmat( 56)
         mat( 58) = lmat( 58)
         mat( 60) = mat( 60) + lmat( 60)
         mat( 67) = mat( 67) + lmat( 67)
         mat( 68) = mat( 68) + lmat( 68)
         mat( 69) = lmat( 69)
         mat( 70) = lmat( 70)
         mat( 73) = mat( 73) + lmat( 73)
         mat( 78) = mat( 78) + lmat( 78)
         mat( 79) = lmat( 79)
         mat( 80) = lmat( 80)
         mat( 81) = mat( 81) + lmat( 81)
         mat( 82) = lmat( 82)
         mat( 84) = lmat( 84)
         mat( 85) = mat( 85) + lmat( 85)
         mat( 86) = mat( 86) + lmat( 86)
         mat( 89) = mat( 89) + lmat( 89)
         mat( 95) = mat( 95) + lmat( 95)
         mat( 104) = mat( 104) + lmat( 104)
         mat( 110) = mat( 110) + lmat( 110)
         mat( 117) = mat( 117) + lmat( 117)
         mat( 119) = lmat( 119)
         mat( 120) = mat( 120) + lmat( 120)
         mat( 123) = mat( 123) + lmat( 123)
         mat( 124) = lmat( 124)
         mat( 128) = mat( 128) + lmat( 128)
         mat( 131) = mat( 131) + lmat( 131)
         mat( 132) = mat( 132) + lmat( 132)
         mat( 135) = mat( 135) + lmat( 135)
         mat( 140) = mat( 140) + lmat( 140)
         mat( 151) = mat( 151) + lmat( 151)
         mat( 157) = lmat( 157)
         mat( 161) = lmat( 161)
         mat( 162) = mat( 162) + lmat( 162)
         mat( 169) = mat( 169) + lmat( 169)
         mat( 179) = mat( 179) + lmat( 179)
         mat( 183) = mat( 183) + lmat( 183)
         mat( 185) = mat( 185) + lmat( 185)
         mat( 187) = mat( 187) + lmat( 187)
         mat( 188) = mat( 188) + lmat( 188)
         mat( 189) = mat( 189) + lmat( 189)
         mat( 190) = mat( 190) + lmat( 190)
         mat( 192) = mat( 192) + lmat( 192)
         mat( 194) = lmat( 194)
         mat( 197) = mat( 197) + lmat( 197)
         mat( 199) = mat( 199) + lmat( 199)
         mat( 200) = lmat( 200)
         mat( 201) = lmat( 201)
         mat( 206) = mat( 206) + lmat( 206)
         mat( 215) = mat( 215) + lmat( 215)
         mat( 216) = mat( 216) + lmat( 216)
         mat( 217) = mat( 217) + lmat( 217)
         mat( 218) = mat( 218) + lmat( 218)
         mat( 223) = mat( 223) + lmat( 223)
         mat( 225) = lmat( 225)
         mat( 228) = mat( 228) + lmat( 228)
         mat( 232) = mat( 232) + lmat( 232)
         mat( 242) = mat( 242) + lmat( 242)
         mat( 258) = lmat( 258)
         mat( 261) = mat( 261) + lmat( 261)
         mat( 273) = mat( 273) + lmat( 273)
         mat( 286) = mat( 286) + lmat( 286)
         mat( 288) = mat( 288) + lmat( 288)
         mat( 297) = mat( 297) + lmat( 297)
         mat( 300) = mat( 300) + lmat( 300)
         mat( 301) = lmat( 301)
         mat( 302) = mat( 302) + lmat( 302)
         mat( 308) = mat( 308) + lmat( 308)
         mat( 309) = mat( 309) + lmat( 309)
         mat( 324) = mat( 324) + lmat( 324)
         mat( 335) = lmat( 335)
         mat( 347) = mat( 347) + lmat( 347)
         mat( 371) = mat( 371) + lmat( 371)
         mat( 372) = mat( 372) + lmat( 372)
         mat( 373) = mat( 373) + lmat( 373)
         mat( 391) = mat( 391) + lmat( 391)
         mat( 395) = mat( 395) + lmat( 395)
         mat( 397) = mat( 397) + lmat( 397)
         mat( 406) = mat( 406) + lmat( 406)
         mat( 416) = mat( 416) + lmat( 416)
         mat( 419) = mat( 419) + lmat( 419)
         mat( 422) = mat( 422) + lmat( 422)
         mat( 426) = mat( 426) + lmat( 426)
         mat( 435) = mat( 435) + lmat( 435)
         mat( 439) = mat( 439) + lmat( 439)
         mat( 442) = mat( 442) + lmat( 442)
         mat( 450) = mat( 450) + lmat( 450)
         mat( 451) = mat( 451) + lmat( 451)
         mat( 454) = mat( 454) + lmat( 454)
         mat( 455) = mat( 455) + lmat( 455)
         mat( 460) = lmat( 460)
         mat( 465) = lmat( 465)
         mat( 469) = mat( 469) + lmat( 469)
         mat( 470) = mat( 470) + lmat( 470)
         mat( 471) = mat( 471) + lmat( 471)
         mat( 478) = lmat( 478)
         mat( 479) = lmat( 479)
         mat( 480) = mat( 480) + lmat( 480)
         mat( 483) = lmat( 483)
         mat( 484) = mat( 484) + lmat( 484)
         mat( 488) = mat( 488) + lmat( 488)
         mat( 490) = mat( 490) + lmat( 490)
         mat( 491) = lmat( 491)
         mat( 492) = lmat( 492)
         mat( 493) = mat( 493) + lmat( 493)
         mat( 494) = mat( 494) + lmat( 494)
         mat( 497) = mat( 497) + lmat( 497)
         mat( 510) = mat( 510) + lmat( 510)
         mat( 535) = mat( 535) + lmat( 535)
         mat( 540) = lmat( 540)
         mat( 543) = lmat( 543)
         mat( 545) = mat( 545) + lmat( 545)
         mat( 556) = mat( 556) + lmat( 556)
         mat( 557) = mat( 557) + lmat( 557)
         mat( 564) = mat( 564) + lmat( 564)
         mat( 568) = mat( 568) + lmat( 568)
         mat( 570) = lmat( 570)
         mat( 578) = mat( 578) + lmat( 578)
         mat( 584) = mat( 584) + lmat( 584)
         mat( 588) = mat( 588) + lmat( 588)
         mat( 48) = 0._r8
         mat( 49) = 0._r8
         mat( 50) = 0._r8
         mat( 100) = 0._r8
         mat( 102) = 0._r8
         mat( 130) = 0._r8
         mat( 143) = 0._r8
         mat( 145) = 0._r8
         mat( 146) = 0._r8
         mat( 153) = 0._r8
         mat( 156) = 0._r8
         mat( 158) = 0._r8
         mat( 166) = 0._r8
         mat( 170) = 0._r8
         mat( 173) = 0._r8
         mat( 174) = 0._r8
         mat( 196) = 0._r8
         mat( 203) = 0._r8
         mat( 204) = 0._r8
         mat( 205) = 0._r8
         mat( 220) = 0._r8
         mat( 237) = 0._r8
         mat( 238) = 0._r8
         mat( 243) = 0._r8
         mat( 246) = 0._r8
         mat( 247) = 0._r8
         mat( 248) = 0._r8
         mat( 259) = 0._r8
         mat( 260) = 0._r8
         mat( 262) = 0._r8
         mat( 263) = 0._r8
         mat( 265) = 0._r8
         mat( 266) = 0._r8
         mat( 267) = 0._r8
         mat( 271) = 0._r8
         mat( 272) = 0._r8
         mat( 274) = 0._r8
         mat( 282) = 0._r8
         mat( 291) = 0._r8
         mat( 293) = 0._r8
         mat( 294) = 0._r8
         mat( 295) = 0._r8
         mat( 304) = 0._r8
         mat( 314) = 0._r8
         mat( 315) = 0._r8
         mat( 330) = 0._r8
         mat( 331) = 0._r8
         mat( 333) = 0._r8
         mat( 339) = 0._r8
         mat( 352) = 0._r8
         mat( 353) = 0._r8
         mat( 362) = 0._r8
         mat( 368) = 0._r8
         mat( 369) = 0._r8
         mat( 377) = 0._r8
         mat( 378) = 0._r8
         mat( 379) = 0._r8
         mat( 381) = 0._r8
         mat( 386) = 0._r8
         mat( 389) = 0._r8
         mat( 392) = 0._r8
         mat( 398) = 0._r8
         mat( 400) = 0._r8
         mat( 401) = 0._r8
         mat( 402) = 0._r8
         mat( 421) = 0._r8
         mat( 423) = 0._r8
         mat( 445) = 0._r8
         mat( 452) = 0._r8
         mat( 459) = 0._r8
         mat( 461) = 0._r8
         mat( 462) = 0._r8
         mat( 463) = 0._r8
         mat( 464) = 0._r8
         mat( 466) = 0._r8
         mat( 467) = 0._r8
         mat( 468) = 0._r8
         mat( 472) = 0._r8
         mat( 473) = 0._r8
         mat( 474) = 0._r8
         mat( 475) = 0._r8
         mat( 482) = 0._r8
         mat( 485) = 0._r8
         mat( 486) = 0._r8
         mat( 501) = 0._r8
         mat( 502) = 0._r8
         mat( 504) = 0._r8
         mat( 505) = 0._r8
         mat( 506) = 0._r8
         mat( 509) = 0._r8
         mat( 513) = 0._r8
         mat( 533) = 0._r8
         mat( 534) = 0._r8
         mat( 547) = 0._r8
         mat( 548) = 0._r8
         mat( 549) = 0._r8
         mat( 550) = 0._r8
         mat( 553) = 0._r8
         mat( 554) = 0._r8
         mat( 575) = 0._r8
         mat( 576) = 0._r8
         mat( 579) = 0._r8
         mat( 580) = 0._r8
         mat( 582) = 0._r8
         mat( 583) = 0._r8
         mat( 585) = 0._r8
         mat( 587) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 4) = mat( 4) - dti
         mat( 7) = mat( 7) - dti
         mat( 9) = mat( 9) - dti
         mat( 12) = mat( 12) - dti
         mat( 14) = mat( 14) - dti
         mat( 18) = mat( 18) - dti
         mat( 22) = mat( 22) - dti
         mat( 28) = mat( 28) - dti
         mat( 33) = mat( 33) - dti
         mat( 40) = mat( 40) - dti
         mat( 47) = mat( 47) - dti
         mat( 53) = mat( 53) - dti
         mat( 60) = mat( 60) - dti
         mat( 67) = mat( 67) - dti
         mat( 73) = mat( 73) - dti
         mat( 78) = mat( 78) - dti
         mat( 86) = mat( 86) - dti
         mat( 95) = mat( 95) - dti
         mat( 104) = mat( 104) - dti
         mat( 110) = mat( 110) - dti
         mat( 117) = mat( 117) - dti
         mat( 123) = mat( 123) - dti
         mat( 131) = mat( 131) - dti
         mat( 140) = mat( 140) - dti
         mat( 151) = mat( 151) - dti
         mat( 162) = mat( 162) - dti
         mat( 169) = mat( 169) - dti
         mat( 179) = mat( 179) - dti
         mat( 189) = mat( 189) - dti
         mat( 199) = mat( 199) - dti
         mat( 206) = mat( 206) - dti
         mat( 215) = mat( 215) - dti
         mat( 228) = mat( 228) - dti
         mat( 242) = mat( 242) - dti
         mat( 261) = mat( 261) - dti
         mat( 286) = mat( 286) - dti
         mat( 308) = mat( 308) - dti
         mat( 347) = mat( 347) - dti
         mat( 373) = mat( 373) - dti
         mat( 397) = mat( 397) - dti
         mat( 419) = mat( 419) - dti
         mat( 450) = mat( 450) - dti
         mat( 470) = mat( 470) - dti
         mat( 490) = mat( 490) - dti
         mat( 510) = mat( 510) - dti
         mat( 535) = mat( 535) - dti
         mat( 557) = mat( 557) - dti
         mat( 588) = mat( 588) - dti

      end subroutine nlnmat_finit

      subroutine nlnmat( mat, y, rxt, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)

      call nlnmat01( mat, y, rxt )
      call nlnmat02( mat, y, rxt )
      call nlnmat03( mat, y, rxt )
      call nlnmat_finit( mat, lmat, dti )

      end subroutine nlnmat

      end module mo_nln_matrix
