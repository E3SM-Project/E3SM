
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


         mat(469) = -(rxt(76)*y(2) + rxt(94)*y(3) + rxt(116)*y(9) + rxt(119)*y(10) &
                      + rxt(141)*y(21) + rxt(146)*y(22) + rxt(153)*y(23) + rxt(156) &
                      *y(27) + rxt(182)*y(36))
         mat(525) = -rxt(76)*y(1)
         mat(452) = -rxt(94)*y(1)
         mat(325) = -rxt(116)*y(1)
         mat(348) = -rxt(119)*y(1)
         mat(218) = -rxt(141)*y(1)
         mat(376) = -rxt(146)*y(1)
         mat(303) = -rxt(153)*y(1)
         mat(416) = -rxt(156)*y(1)
         mat(434) = -rxt(182)*y(1)

         mat(525) = mat(525) + rxt(75)*y(4)
         mat(237) = rxt(75)*y(2)

         mat(527) = -(rxt(75)*y(4) + rxt(76)*y(1) + 4._r8*rxt(77)*y(2) + rxt(114)*y(9) &
                      + (rxt(117) + rxt(118)) * y(10) + rxt(125)*y(11) + rxt(137) &
                      *y(18) + rxt(145)*y(22) + rxt(152)*y(23) + rxt(155)*y(24) &
                      + rxt(163)*y(29) + rxt(175)*y(32) + rxt(176)*y(33) + rxt(179) &
                      *y(34) + rxt(185)*y(37) + rxt(195)*y(38) + rxt(196)*y(39) &
                      + rxt(197)*y(40) + (rxt(230) + rxt(239)) * y(52) + rxt(236) &
                      *y(54))
         mat(238) = -rxt(75)*y(2)
         mat(471) = -rxt(76)*y(2)
         mat(327) = -rxt(114)*y(2)
         mat(350) = -(rxt(117) + rxt(118)) * y(2)
         mat(397) = -rxt(125)*y(2)
         mat(211) = -rxt(137)*y(2)
         mat(378) = -rxt(145)*y(2)
         mat(305) = -rxt(152)*y(2)
         mat(105) = -rxt(155)*y(2)
         mat(281) = -rxt(163)*y(2)
         mat(493) = -rxt(175)*y(2)
         mat(180) = -rxt(176)*y(2)
         mat(192) = -rxt(179)*y(2)
         mat(259) = -rxt(185)*y(2)
         mat(111) = -rxt(195)*y(2)
         mat(162) = -rxt(196)*y(2)
         mat(87) = -rxt(197)*y(2)
         mat(98) = -(rxt(230) + rxt(239)) * y(2)
         mat(61) = -rxt(236)*y(2)

         mat(454) = (rxt(89)+rxt(90))*y(4)
         mat(238) = mat(238) + (rxt(89)+rxt(90))*y(3) + rxt(111)*y(8) + rxt(235)*y(54)  &
                      + rxt(228)*y(55)
         mat(153) = rxt(111)*y(4) + rxt(112)*y(9) + rxt(113)*y(10) + rxt(232)*y(53)
         mat(327) = mat(327) + rxt(112)*y(8)
         mat(350) = mat(350) + rxt(113)*y(8)
         mat(378) = mat(378) + 2.000_r8*rxt(148)*y(22)
         mat(219) = rxt(144)*y(23)
         mat(305) = mat(305) + rxt(144)*y(21)
         mat(124) = rxt(232)*y(8) + 1.150_r8*rxt(241)*y(57)
         mat(61) = mat(61) + rxt(235)*y(4)
         mat(74) = rxt(228)*y(4)
         mat(132) = rxt(240)*y(57)
         mat(142) = 1.150_r8*rxt(241)*y(53) + rxt(240)*y(56)

         mat(451) = -((rxt(89) + rxt(90)) * y(4) + rxt(91)*y(78) + rxt(94)*y(1) &
                      + rxt(107)*y(32) + rxt(108)*y(38))
         mat(236) = -(rxt(89) + rxt(90)) * y(3)
         mat(196) = -rxt(91)*y(3)
         mat(468) = -rxt(94)*y(3)
         mat(490) = -rxt(107)*y(3)
         mat(110) = -rxt(108)*y(3)

         mat(236) = mat(236) + rxt(109)*y(58)
         mat(123) = .850_r8*rxt(241)*y(57)
         mat(79) = rxt(109)*y(4)
         mat(141) = .850_r8*rxt(241)*y(53)

         mat(231) = -(rxt(75)*y(2) + rxt(85)*y(6) + rxt(89)*y(3) + rxt(109)*y(58) &
                      + rxt(111)*y(8) + rxt(140)*y(21) + rxt(228)*y(55) + (rxt(234) &
                      + rxt(235)) * y(54) + rxt(237)*y(52))
         mat(514) = -rxt(75)*y(4)
         mat(28) = -rxt(85)*y(4)
         mat(443) = -rxt(89)*y(4)
         mat(77) = -rxt(109)*y(4)
         mat(149) = -rxt(111)*y(4)
         mat(214) = -rxt(140)*y(4)
         mat(73) = -rxt(228)*y(4)
         mat(60) = -(rxt(234) + rxt(235)) * y(4)
         mat(95) = -rxt(237)*y(4)

         mat(458) = 2.000_r8*rxt(76)*y(2) + 2.000_r8*rxt(94)*y(3) + rxt(116)*y(9)  &
                      + rxt(119)*y(10) + rxt(146)*y(22) + rxt(141)*y(21)  &
                      + 2.000_r8*rxt(153)*y(23) + rxt(156)*y(27) + rxt(182)*y(36)
         mat(514) = mat(514) + 2.000_r8*rxt(76)*y(1) + 2.000_r8*rxt(77)*y(2) + rxt(84) &
                      *y(6) + rxt(117)*y(10) + rxt(145)*y(22) + rxt(125)*y(11)  &
                      + rxt(152)*y(23) + rxt(163)*y(29) + rxt(185)*y(37)
         mat(443) = mat(443) + 2.000_r8*rxt(94)*y(1)
         mat(231) = mat(231) + 2.000_r8*rxt(85)*y(6)
         mat(28) = mat(28) + rxt(84)*y(2) + 2.000_r8*rxt(85)*y(4)
         mat(314) = rxt(116)*y(1) + rxt(233)*y(53)
         mat(337) = rxt(119)*y(1) + rxt(117)*y(2)
         mat(365) = rxt(146)*y(1) + rxt(145)*y(2) + rxt(129)*y(13) + rxt(147)*y(23)  &
                      + rxt(165)*y(29)
         mat(386) = rxt(125)*y(2) + rxt(127)*y(23)
         mat(64) = rxt(129)*y(22)
         mat(168) = rxt(133)*y(23)
         mat(214) = mat(214) + rxt(141)*y(1) + rxt(143)*y(23)
         mat(292) = 2.000_r8*rxt(153)*y(1) + rxt(152)*y(2) + rxt(147)*y(22) + rxt(127) &
                      *y(11) + rxt(133)*y(16) + rxt(143)*y(21) + 2.000_r8*rxt(154) &
                      *y(23) + rxt(159)*y(27) + rxt(166)*y(29) + rxt(183)*y(36)  &
                      + rxt(187)*y(37)
         mat(406) = rxt(156)*y(1) + rxt(159)*y(23)
         mat(268) = rxt(163)*y(2) + rxt(165)*y(22) + rxt(166)*y(23) + ( &
                      + 2.000_r8*rxt(169)+2.000_r8*rxt(170))*y(29) + (rxt(191) &
                       +rxt(192))*y(37)
         mat(423) = rxt(182)*y(1) + rxt(183)*y(23)
         mat(246) = rxt(185)*y(2) + rxt(187)*y(23) + (rxt(191)+rxt(192))*y(29)  &
                      + 2.000_r8*rxt(193)*y(37)
         mat(121) = rxt(233)*y(9)

         mat(30) = -(rxt(78)*y(2) + rxt(79)*y(4) + rxt(81)*y(1))
         mat(495) = -rxt(78)*y(5)
         mat(221) = -rxt(79)*y(5)
         mat(456) = -rxt(81)*y(5)

         mat(437) = rxt(89)*y(4)
         mat(221) = mat(221) + rxt(89)*y(3)

         mat(27) = -(rxt(84)*y(2) + rxt(85)*y(4))
         mat(494) = -rxt(84)*y(6)
         mat(220) = -rxt(85)*y(6)

         mat(455) = rxt(81)*y(5)
         mat(494) = mat(494) + rxt(78)*y(5)
         mat(220) = mat(220) + rxt(79)*y(5)
         mat(29) = rxt(81)*y(1) + rxt(78)*y(2) + rxt(79)*y(4)

         mat(148) = -(rxt(111)*y(4) + rxt(112)*y(9) + rxt(113)*y(10) + rxt(232)*y(53))
         mat(229) = -rxt(111)*y(8)
         mat(309) = -rxt(112)*y(8)
         mat(332) = -rxt(113)*y(8)
         mat(120) = -rxt(232)*y(8)

         mat(507) = rxt(236)*y(54) + rxt(110)*y(58)
         mat(229) = mat(229) + rxt(234)*y(54)
         mat(94) = 1.100_r8*rxt(242)*y(57)
         mat(59) = rxt(236)*y(2) + rxt(234)*y(4)
         mat(128) = .200_r8*rxt(240)*y(57)
         mat(76) = rxt(110)*y(2)
         mat(138) = 1.100_r8*rxt(242)*y(52) + .200_r8*rxt(240)*y(56)

         mat(318) = -(rxt(112)*y(8) + rxt(114)*y(2) + rxt(115)*y(23) + rxt(116)*y(1) &
                      + rxt(124)*y(11) + rxt(132)*y(16) + rxt(167)*y(29) + rxt(188) &
                      *y(37) + rxt(233)*y(53))
         mat(150) = -rxt(112)*y(9)
         mat(518) = -rxt(114)*y(9)
         mat(296) = -rxt(115)*y(9)
         mat(462) = -rxt(116)*y(9)
         mat(388) = -rxt(124)*y(9)
         mat(170) = -rxt(132)*y(9)
         mat(272) = -rxt(167)*y(9)
         mat(250) = -rxt(188)*y(9)
         mat(122) = -rxt(233)*y(9)

         mat(518) = mat(518) + rxt(117)*y(10)
         mat(233) = rxt(111)*y(8) + rxt(109)*y(58)
         mat(150) = mat(150) + rxt(111)*y(4)
         mat(341) = rxt(117)*y(2)
         mat(78) = rxt(109)*y(4)

         mat(342) = -(rxt(113)*y(8) + (rxt(117) + rxt(118)) * y(2) + rxt(119)*y(1) &
                      + rxt(120)*y(11) + rxt(122)*y(22) + rxt(128)*y(23) + rxt(168) &
                      *y(29) + rxt(189)*y(37))
         mat(151) = -rxt(113)*y(10)
         mat(519) = -(rxt(117) + rxt(118)) * y(10)
         mat(463) = -rxt(119)*y(10)
         mat(389) = -rxt(120)*y(10)
         mat(370) = -rxt(122)*y(10)
         mat(297) = -rxt(128)*y(10)
         mat(273) = -rxt(168)*y(10)
         mat(251) = -rxt(189)*y(10)

         mat(463) = mat(463) + rxt(116)*y(9)
         mat(519) = mat(519) + rxt(114)*y(9) + rxt(125)*y(11)
         mat(319) = rxt(116)*y(1) + rxt(114)*y(2) + 2.000_r8*rxt(124)*y(11) + rxt(132) &
                      *y(16) + rxt(115)*y(23) + rxt(167)*y(29) + rxt(188)*y(37)
         mat(370) = mat(370) + rxt(126)*y(11) + rxt(129)*y(13)
         mat(389) = mat(389) + rxt(125)*y(2) + 2.000_r8*rxt(124)*y(9) + rxt(126)*y(22)  &
                      + rxt(127)*y(23)
         mat(66) = rxt(129)*y(22)
         mat(171) = rxt(132)*y(9)
         mat(297) = mat(297) + rxt(115)*y(9) + rxt(127)*y(11)
         mat(273) = mat(273) + rxt(167)*y(9)
         mat(251) = mat(251) + rxt(188)*y(9)

         mat(371) = -(rxt(122)*y(10) + rxt(123)*y(12) + rxt(126)*y(11) + rxt(129) &
                      *y(13) + rxt(134)*y(17) + rxt(136)*y(18) + rxt(145)*y(2) + rxt(146) &
                      *y(1) + rxt(147)*y(23) + (4._r8*rxt(148) + 4._r8*rxt(149) &
                      ) * y(22) + rxt(151)*y(24) + (rxt(164) + rxt(165)) * y(29) &
                      + rxt(174)*y(32) + rxt(178)*y(33) + rxt(180)*y(34) + rxt(186) &
                      *y(37) + rxt(194)*y(38) + rxt(201)*y(59) + (rxt(202) + rxt(203) &
                      ) * y(60))
         mat(343) = -rxt(122)*y(22)
         mat(115) = -rxt(123)*y(22)
         mat(390) = -rxt(126)*y(22)
         mat(67) = -rxt(129)*y(22)
         mat(54) = -rxt(134)*y(22)
         mat(205) = -rxt(136)*y(22)
         mat(520) = -rxt(145)*y(22)
         mat(464) = -rxt(146)*y(22)
         mat(298) = -rxt(147)*y(22)
         mat(102) = -rxt(151)*y(22)
         mat(274) = -(rxt(164) + rxt(165)) * y(22)
         mat(486) = -rxt(174)*y(22)
         mat(177) = -rxt(178)*y(22)
         mat(188) = -rxt(180)*y(22)
         mat(252) = -rxt(186)*y(22)
         mat(108) = -rxt(194)*y(22)
         mat(23) = -rxt(201)*y(22)
         mat(41) = -(rxt(202) + rxt(203)) * y(22)

         mat(464) = mat(464) + rxt(141)*y(21) + rxt(153)*y(23)
         mat(520) = mat(520) + rxt(137)*y(18) + rxt(152)*y(23) + rxt(155)*y(24)  &
                      + rxt(175)*y(32) + rxt(176)*y(33) + rxt(195)*y(38) + rxt(196) &
                      *y(39)
         mat(447) = 2.000_r8*rxt(91)*y(78) + rxt(107)*y(32) + rxt(108)*y(38)
         mat(320) = rxt(115)*y(23)
         mat(390) = mat(390) + rxt(127)*y(23)
         mat(205) = mat(205) + rxt(137)*y(2)
         mat(216) = rxt(141)*y(1) + 2.000_r8*rxt(142)*y(23)
         mat(298) = mat(298) + rxt(153)*y(1) + rxt(152)*y(2) + rxt(115)*y(9)  &
                      + rxt(127)*y(11) + 2.000_r8*rxt(142)*y(21) + rxt(160)*y(27)
         mat(102) = mat(102) + rxt(155)*y(2)
         mat(195) = 2.000_r8*rxt(91)*y(3)
         mat(411) = rxt(160)*y(23)
         mat(486) = mat(486) + rxt(175)*y(2) + rxt(107)*y(3)
         mat(177) = mat(177) + rxt(176)*y(2)
         mat(108) = mat(108) + rxt(195)*y(2) + rxt(108)*y(3)
         mat(158) = rxt(196)*y(2)


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


         mat(391) = -(rxt(120)*y(10) + rxt(124)*y(9) + rxt(125)*y(2) + rxt(126)*y(22) &
                      + rxt(127)*y(23) + rxt(135)*y(18) + rxt(204)*y(60))
         mat(344) = -rxt(120)*y(11)
         mat(321) = -rxt(124)*y(11)
         mat(521) = -rxt(125)*y(11)
         mat(372) = -rxt(126)*y(11)
         mat(299) = -rxt(127)*y(11)
         mat(206) = -rxt(135)*y(11)
         mat(42) = -rxt(204)*y(11)

         mat(465) = rxt(119)*y(10)
         mat(521) = mat(521) + rxt(118)*y(10) + rxt(179)*y(34) + rxt(197)*y(40)
         mat(344) = mat(344) + rxt(119)*y(1) + rxt(118)*y(2)
         mat(372) = mat(372) + rxt(123)*y(12) + rxt(180)*y(34)
         mat(116) = rxt(123)*y(22)
         mat(412) = rxt(181)*y(34)
         mat(189) = rxt(179)*y(2) + rxt(180)*y(22) + rxt(181)*y(27)
         mat(85) = rxt(197)*y(2)

         mat(112) = -(rxt(123)*y(22))
         mat(358) = -rxt(123)*y(12)

         mat(331) = rxt(122)*y(22)
         mat(358) = mat(358) + rxt(122)*y(10)
         mat(382) = rxt(135)*y(18) + rxt(204)*y(60)
         mat(199) = rxt(135)*y(11)
         mat(474) = (rxt(214)+rxt(219)+rxt(225))*y(34)
         mat(182) = (rxt(214)+rxt(219)+rxt(225))*y(32)
         mat(39) = rxt(204)*y(11)

         mat(62) = -(rxt(129)*y(22))
         mat(355) = -rxt(129)*y(13)

         mat(329) = rxt(128)*y(23)
         mat(283) = rxt(128)*y(10)


         mat(328) = rxt(120)*y(11)
         mat(381) = rxt(120)*y(10)

         mat(164) = -(rxt(132)*y(9) + rxt(133)*y(23))
         mat(310) = -rxt(132)*y(16)
         mat(287) = -rxt(133)*y(16)

         mat(359) = rxt(134)*y(17)
         mat(50) = rxt(134)*y(22)

         mat(49) = -(rxt(134)*y(22))
         mat(354) = -rxt(134)*y(17)

         mat(163) = rxt(133)*y(23)
         mat(282) = rxt(133)*y(16)

         mat(201) = -(rxt(135)*y(11) + rxt(136)*y(22) + rxt(137)*y(2) + rxt(161)*y(27) &
                      + rxt(184)*y(36))
         mat(384) = -rxt(135)*y(18)
         mat(363) = -rxt(136)*y(18)
         mat(512) = -rxt(137)*y(18)
         mat(404) = -rxt(161)*y(18)
         mat(421) = -rxt(184)*y(18)

         mat(312) = rxt(132)*y(16)
         mat(166) = rxt(132)*y(9)

         mat(213) = -(rxt(140)*y(4) + rxt(141)*y(1) + (rxt(142) + rxt(143) + rxt(144) &
                      ) * y(23))
         mat(230) = -rxt(140)*y(21)
         mat(457) = -rxt(141)*y(21)
         mat(291) = -(rxt(142) + rxt(143) + rxt(144)) * y(21)

         mat(513) = rxt(145)*y(22)
         mat(364) = rxt(145)*y(2) + rxt(136)*y(18)
         mat(202) = rxt(136)*y(22)

         mat(295) = -(rxt(115)*y(9) + rxt(127)*y(11) + rxt(128)*y(10) + rxt(133)*y(16) &
                      + (rxt(142) + rxt(143) + rxt(144)) * y(21) + rxt(147)*y(22) &
                      + rxt(152)*y(2) + rxt(153)*y(1) + 4._r8*rxt(154)*y(23) + (rxt(159) &
                      + rxt(160)) * y(27) + rxt(166)*y(29) + rxt(183)*y(36) + rxt(187) &
                      *y(37))
         mat(317) = -rxt(115)*y(23)
         mat(387) = -rxt(127)*y(23)
         mat(340) = -rxt(128)*y(23)
         mat(169) = -rxt(133)*y(23)
         mat(215) = -(rxt(142) + rxt(143) + rxt(144)) * y(23)
         mat(368) = -rxt(147)*y(23)
         mat(517) = -rxt(152)*y(23)
         mat(461) = -rxt(153)*y(23)
         mat(408) = -(rxt(159) + rxt(160)) * y(23)
         mat(271) = -rxt(166)*y(23)
         mat(426) = -rxt(183)*y(23)
         mat(249) = -rxt(187)*y(23)

         mat(461) = mat(461) + rxt(146)*y(22)
         mat(517) = mat(517) + rxt(137)*y(18) + rxt(155)*y(24)
         mat(232) = rxt(140)*y(21)
         mat(317) = mat(317) + rxt(132)*y(16)
         mat(368) = mat(368) + rxt(146)*y(1) + rxt(126)*y(11) + rxt(151)*y(24)  &
                      + rxt(164)*y(29) + rxt(186)*y(37) + .500_r8*rxt(203)*y(60)
         mat(387) = mat(387) + rxt(126)*y(22) + rxt(135)*y(18)
         mat(169) = mat(169) + rxt(132)*y(9)
         mat(203) = rxt(137)*y(2) + rxt(135)*y(11) + rxt(161)*y(27) + rxt(184)*y(36)
         mat(215) = mat(215) + rxt(140)*y(4)
         mat(101) = rxt(155)*y(2) + rxt(151)*y(22) + rxt(158)*y(27)
         mat(408) = mat(408) + rxt(161)*y(18) + rxt(158)*y(24)
         mat(271) = mat(271) + rxt(164)*y(22)
         mat(426) = mat(426) + rxt(184)*y(18)
         mat(249) = mat(249) + rxt(186)*y(22)
         mat(40) = .500_r8*rxt(203)*y(22)

         mat(99) = -(rxt(151)*y(22) + rxt(155)*y(2) + rxt(158)*y(27))
         mat(356) = -rxt(151)*y(24)
         mat(501) = -rxt(155)*y(24)
         mat(399) = -rxt(158)*y(24)

         mat(356) = mat(356) + 2.000_r8*rxt(149)*y(22)
         mat(284) = 2.000_r8*rxt(154)*y(23)

         mat(193) = -(rxt(91)*y(3))
         mat(440) = -rxt(91)*y(78)

         mat(362) = 2.000_r8*rxt(148)*y(22) + rxt(123)*y(12) + rxt(129)*y(13)  &
                      + rxt(134)*y(17) + rxt(136)*y(18) + rxt(147)*y(23) + rxt(151) &
                      *y(24) + rxt(174)*y(32) + rxt(178)*y(33) + rxt(194)*y(38)
         mat(113) = rxt(123)*y(22)
         mat(63) = rxt(129)*y(22)
         mat(51) = rxt(134)*y(22)
         mat(200) = rxt(136)*y(22)
         mat(212) = rxt(144)*y(23)
         mat(289) = rxt(147)*y(22) + rxt(144)*y(21)
         mat(100) = rxt(151)*y(22)
         mat(478) = rxt(174)*y(22) + (rxt(215)+rxt(220)+rxt(226))*y(33) + (rxt(216) &
                       +rxt(227))*y(39)
         mat(175) = rxt(178)*y(22) + (rxt(215)+rxt(220)+rxt(226))*y(32)
         mat(107) = rxt(194)*y(22)
         mat(156) = (rxt(216)+rxt(227))*y(32)

         mat(413) = -(rxt(156)*y(1) + rxt(158)*y(24) + (rxt(159) + rxt(160)) * y(23) &
                      + rxt(161)*y(18) + rxt(177)*y(33) + rxt(181)*y(34))
         mat(466) = -rxt(156)*y(27)
         mat(103) = -rxt(158)*y(27)
         mat(300) = -(rxt(159) + rxt(160)) * y(27)
         mat(207) = -rxt(161)*y(27)
         mat(178) = -rxt(177)*y(27)
         mat(190) = -rxt(181)*y(27)

         mat(522) = rxt(163)*y(29) + rxt(175)*y(32)
         mat(449) = rxt(107)*y(32)
         mat(322) = rxt(167)*y(29)
         mat(373) = rxt(164)*y(29) + rxt(174)*y(32)
         mat(276) = rxt(163)*y(2) + rxt(167)*y(9) + rxt(164)*y(22) + ( &
                      + 4.000_r8*rxt(169)+2.000_r8*rxt(171))*y(29) + rxt(191)*y(37)
         mat(488) = rxt(175)*y(2) + rxt(107)*y(3) + rxt(174)*y(22)
         mat(254) = rxt(191)*y(29)


         mat(398) = rxt(181)*y(34)
         mat(262) = 2.000_r8*rxt(170)*y(29)
         mat(472) = (rxt(215)+rxt(220)+rxt(226))*y(33) + (rxt(214)+rxt(219)+rxt(225)) &
                      *y(34)
         mat(173) = (rxt(215)+rxt(220)+rxt(226))*y(32)
         mat(181) = rxt(181)*y(27) + (rxt(214)+rxt(219)+rxt(225))*y(32)

         mat(270) = -(rxt(163)*y(2) + (rxt(164) + rxt(165)) * y(22) + rxt(166)*y(23) &
                      + rxt(167)*y(9) + rxt(168)*y(10) + (4._r8*rxt(169) + 4._r8*rxt(170) &
                      + 4._r8*rxt(171) + 4._r8*rxt(172)) * y(29) + (rxt(190) + rxt(191) &
                      + rxt(192)) * y(37))
         mat(516) = -rxt(163)*y(29)
         mat(367) = -(rxt(164) + rxt(165)) * y(29)
         mat(294) = -rxt(166)*y(29)
         mat(316) = -rxt(167)*y(29)
         mat(339) = -rxt(168)*y(29)
         mat(248) = -(rxt(190) + rxt(191) + rxt(192)) * y(29)

         mat(460) = rxt(156)*y(27)
         mat(516) = mat(516) + rxt(176)*y(33) + rxt(179)*y(34)
         mat(367) = mat(367) + rxt(178)*y(33)
         mat(294) = mat(294) + rxt(160)*y(27)
         mat(407) = rxt(156)*y(1) + rxt(160)*y(23) + rxt(177)*y(33)
         mat(176) = rxt(176)*y(2) + rxt(178)*y(22) + rxt(177)*y(27)
         mat(186) = rxt(179)*y(2)


         mat(261) = 2.000_r8*rxt(171)*y(29) + rxt(190)*y(37)
         mat(239) = rxt(190)*y(29)


         mat(260) = 2.000_r8*rxt(172)*y(29)

         mat(492) = -(rxt(107)*y(3) + rxt(174)*y(22) + rxt(175)*y(2) + (rxt(214) &
                      + rxt(219) + rxt(225)) * y(34) + (rxt(215) + rxt(220) + rxt(226) &
                      ) * y(33) + (rxt(216) + rxt(227)) * y(39))
         mat(453) = -rxt(107)*y(32)
         mat(377) = -rxt(174)*y(32)
         mat(526) = -rxt(175)*y(32)
         mat(191) = -(rxt(214) + rxt(219) + rxt(225)) * y(32)
         mat(179) = -(rxt(215) + rxt(220) + rxt(226)) * y(32)
         mat(161) = -(rxt(216) + rxt(227)) * y(32)

         mat(377) = mat(377) + rxt(165)*y(29)
         mat(210) = rxt(161)*y(27)
         mat(304) = rxt(159)*y(27)
         mat(104) = rxt(158)*y(27)
         mat(417) = rxt(161)*y(18) + rxt(159)*y(23) + rxt(158)*y(24) + rxt(177)*y(33)
         mat(280) = rxt(165)*y(22)
         mat(179) = mat(179) + rxt(177)*y(27)

         mat(174) = -(rxt(176)*y(2) + rxt(177)*y(27) + rxt(178)*y(22) + (rxt(215) &
                      + rxt(220) + rxt(226)) * y(32))
         mat(509) = -rxt(176)*y(33)
         mat(401) = -rxt(177)*y(33)
         mat(360) = -rxt(178)*y(33)
         mat(476) = -(rxt(215) + rxt(220) + rxt(226)) * y(33)

         mat(360) = mat(360) + rxt(180)*y(34)
         mat(288) = rxt(166)*y(29)
         mat(264) = rxt(166)*y(23)
         mat(183) = rxt(180)*y(22)

         mat(184) = -(rxt(179)*y(2) + rxt(180)*y(22) + rxt(181)*y(27) + (rxt(214) &
                      + rxt(219) + rxt(225)) * y(32))
         mat(510) = -rxt(179)*y(34)
         mat(361) = -rxt(180)*y(34)
         mat(402) = -rxt(181)*y(34)
         mat(477) = -(rxt(214) + rxt(219) + rxt(225)) * y(34)

         mat(334) = rxt(168)*y(29)
         mat(265) = rxt(168)*y(10)


         mat(263) = rxt(192)*y(37)
         mat(473) = (rxt(216)+rxt(227))*y(39)
         mat(240) = rxt(192)*y(29)
         mat(154) = (rxt(216)+rxt(227))*y(32)

         mat(432) = -(rxt(182)*y(1) + rxt(183)*y(23) + rxt(184)*y(18))
         mat(467) = -rxt(182)*y(36)
         mat(301) = -rxt(183)*y(36)
         mat(208) = -rxt(184)*y(36)

         mat(523) = rxt(185)*y(37) + rxt(195)*y(38)
         mat(450) = rxt(108)*y(38)
         mat(323) = rxt(188)*y(37)
         mat(374) = rxt(186)*y(37) + rxt(194)*y(38)
         mat(277) = (rxt(190)+rxt(191))*y(37)
         mat(255) = rxt(185)*y(2) + rxt(188)*y(9) + rxt(186)*y(22) + (rxt(190) &
                       +rxt(191))*y(29) + 4.000_r8*rxt(193)*y(37)
         mat(109) = rxt(195)*y(2) + rxt(108)*y(3) + rxt(194)*y(22)


      end subroutine     nlnmat02

      subroutine     nlnmat03( mat, y, rxt )

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


         mat(247) = -(rxt(185)*y(2) + rxt(186)*y(22) + rxt(187)*y(23) + rxt(188)*y(9) &
                      + rxt(189)*y(10) + (rxt(190) + rxt(191) + rxt(192)) * y(29) &
                      + 4._r8*rxt(193)*y(37))
         mat(515) = -rxt(185)*y(37)
         mat(366) = -rxt(186)*y(37)
         mat(293) = -rxt(187)*y(37)
         mat(315) = -rxt(188)*y(37)
         mat(338) = -rxt(189)*y(37)
         mat(269) = -(rxt(190) + rxt(191) + rxt(192)) * y(37)

         mat(459) = rxt(182)*y(36)
         mat(515) = mat(515) + rxt(196)*y(39) + rxt(197)*y(40)
         mat(424) = rxt(182)*y(1)
         mat(157) = rxt(196)*y(2)
         mat(83) = rxt(197)*y(2)

         mat(106) = -(rxt(108)*y(3) + rxt(194)*y(22) + rxt(195)*y(2))
         mat(438) = -rxt(108)*y(38)
         mat(357) = -rxt(194)*y(38)
         mat(502) = -rxt(195)*y(38)

         mat(198) = rxt(184)*y(36)
         mat(285) = rxt(183)*y(36)
         mat(419) = rxt(184)*y(18) + rxt(183)*y(23)

         mat(155) = -(rxt(196)*y(2) + (rxt(216) + rxt(227)) * y(32))
         mat(508) = -rxt(196)*y(39)
         mat(475) = -(rxt(216) + rxt(227)) * y(39)

         mat(286) = rxt(187)*y(37)
         mat(243) = rxt(187)*y(23)

         mat(80) = -(rxt(197)*y(2))
         mat(499) = -rxt(197)*y(40)

         mat(330) = rxt(189)*y(37)
         mat(241) = rxt(189)*y(10)

         mat(90) = -((rxt(230) + rxt(239)) * y(2) + rxt(237)*y(4) + rxt(242)*y(57))
         mat(500) = -(rxt(230) + rxt(239)) * y(52)
         mat(225) = -rxt(237)*y(52)
         mat(134) = -rxt(242)*y(52)

         mat(117) = -(rxt(232)*y(8) + rxt(233)*y(9) + rxt(241)*y(57))
         mat(145) = -rxt(232)*y(53)
         mat(306) = -rxt(233)*y(53)
         mat(135) = -rxt(241)*y(53)

         mat(226) = rxt(237)*y(52) + rxt(234)*y(54) + rxt(228)*y(55)
         mat(91) = rxt(237)*y(4)
         mat(57) = rxt(234)*y(4)
         mat(70) = rxt(228)*y(4)

         mat(55) = -((rxt(234) + rxt(235)) * y(4) + rxt(236)*y(2))
         mat(222) = -(rxt(234) + rxt(235)) * y(54)
         mat(496) = -rxt(236)*y(54)

         mat(69) = -(rxt(228)*y(4))
         mat(223) = -rxt(228)*y(55)

         mat(497) = rxt(239)*y(52) + rxt(236)*y(54)
         mat(88) = rxt(239)*y(2)
         mat(56) = rxt(236)*y(2)

         mat(126) = -(rxt(240)*y(57))
         mat(136) = -rxt(240)*y(56)

         mat(505) = rxt(230)*y(52)
         mat(227) = rxt(235)*y(54)
         mat(146) = rxt(232)*y(53)
         mat(307) = rxt(233)*y(53)
         mat(92) = rxt(230)*y(2)
         mat(118) = rxt(232)*y(8) + rxt(233)*y(9)
         mat(58) = rxt(235)*y(4)

         mat(75) = -(rxt(109)*y(4) + rxt(110)*y(2))
         mat(224) = -rxt(109)*y(58)
         mat(498) = -rxt(110)*y(58)

         mat(498) = mat(498) + rxt(230)*y(52)
         mat(89) = rxt(230)*y(2) + .900_r8*rxt(242)*y(57)
         mat(125) = .800_r8*rxt(240)*y(57)
         mat(133) = .900_r8*rxt(242)*y(52) + .800_r8*rxt(240)*y(56)

         mat(137) = -(rxt(240)*y(56) + rxt(241)*y(53) + rxt(242)*y(52))
         mat(127) = -rxt(240)*y(57)
         mat(119) = -rxt(241)*y(57)
         mat(93) = -rxt(242)*y(57)

         mat(22) = -(rxt(201)*y(22))
         mat(352) = -rxt(201)*y(59)

         mat(352) = mat(352) + (rxt(202)+.500_r8*rxt(203))*y(60)
         mat(379) = rxt(204)*y(60)
         mat(37) = (rxt(202)+.500_r8*rxt(203))*y(22) + rxt(204)*y(11)

         mat(38) = -((rxt(202) + rxt(203)) * y(22) + rxt(204)*y(11))
         mat(353) = -(rxt(202) + rxt(203)) * y(60)
         mat(380) = -rxt(204)*y(60)


         mat(351) = rxt(201)*y(59)
         mat(21) = rxt(201)*y(22)


















      end subroutine     nlnmat03

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
         mat(  18) = lmat(  18)
         mat(  19) = lmat(  19)
         mat(  20) = lmat(  20)
         mat(  22) = mat(  22) + lmat(  22)
         mat(  24) = lmat(  24)
         mat(  25) = lmat(  25)
         mat(  26) = lmat(  26)
         mat(  27) = mat(  27) + lmat(  27)
         mat(  28) = mat(  28) + lmat(  28)
         mat(  29) = mat(  29) + lmat(  29)
         mat(  30) = mat(  30) + lmat(  30)
         mat(  31) = lmat(  31)
         mat(  32) = lmat(  32)
         mat(  33) = lmat(  33)
         mat(  34) = lmat(  34)
         mat(  35) = lmat(  35)
         mat(  36) = lmat(  36)
         mat(  38) = mat(  38) + lmat(  38)
         mat(  43) = lmat(  43)
         mat(  44) = lmat(  44)
         mat(  45) = lmat(  45)
         mat(  46) = lmat(  46)
         mat(  47) = lmat(  47)
         mat(  48) = lmat(  48)
         mat(  49) = mat(  49) + lmat(  49)
         mat(  52) = lmat(  52)
         mat(  53) = lmat(  53)
         mat(  54) = mat(  54) + lmat(  54)
         mat(  55) = mat(  55) + lmat(  55)
         mat(  62) = mat(  62) + lmat(  62)
         mat(  65) = lmat(  65)
         mat(  66) = mat(  66) + lmat(  66)
         mat(  67) = mat(  67) + lmat(  67)
         mat(  68) = lmat(  68)
         mat(  69) = mat(  69) + lmat(  69)
         mat(  70) = mat(  70) + lmat(  70)
         mat(  71) = lmat(  71)
         mat(  72) = lmat(  72)
         mat(  75) = mat(  75) + lmat(  75)
         mat(  80) = mat(  80) + lmat(  80)
         mat(  81) = lmat(  81)
         mat(  82) = lmat(  82)
         mat(  83) = mat(  83) + lmat(  83)
         mat(  84) = lmat(  84)
         mat(  85) = mat(  85) + lmat(  85)
         mat(  86) = lmat(  86)
         mat(  90) = mat(  90) + lmat(  90)
         mat(  99) = mat(  99) + lmat(  99)
         mat( 102) = mat( 102) + lmat( 102)
         mat( 106) = mat( 106) + lmat( 106)
         mat( 112) = mat( 112) + lmat( 112)
         mat( 114) = lmat( 114)
         mat( 115) = mat( 115) + lmat( 115)
         mat( 117) = mat( 117) + lmat( 117)
         mat( 118) = mat( 118) + lmat( 118)
         mat( 122) = mat( 122) + lmat( 122)
         mat( 126) = mat( 126) + lmat( 126)
         mat( 137) = mat( 137) + lmat( 137)
         mat( 143) = lmat( 143)
         mat( 147) = lmat( 147)
         mat( 148) = mat( 148) + lmat( 148)
         mat( 155) = mat( 155) + lmat( 155)
         mat( 158) = mat( 158) + lmat( 158)
         mat( 160) = lmat( 160)
         mat( 164) = mat( 164) + lmat( 164)
         mat( 174) = mat( 174) + lmat( 174)
         mat( 177) = mat( 177) + lmat( 177)
         mat( 178) = mat( 178) + lmat( 178)
         mat( 182) = mat( 182) + lmat( 182)
         mat( 183) = mat( 183) + lmat( 183)
         mat( 184) = mat( 184) + lmat( 184)
         mat( 186) = mat( 186) + lmat( 186)
         mat( 187) = lmat( 187)
         mat( 189) = mat( 189) + lmat( 189)
         mat( 190) = mat( 190) + lmat( 190)
         mat( 193) = mat( 193) + lmat( 193)
         mat( 194) = lmat( 194)
         mat( 195) = mat( 195) + lmat( 195)
         mat( 196) = mat( 196) + lmat( 196)
         mat( 197) = lmat( 197)
         mat( 201) = mat( 201) + lmat( 201)
         mat( 202) = mat( 202) + lmat( 202)
         mat( 213) = mat( 213) + lmat( 213)
         mat( 223) = mat( 223) + lmat( 223)
         mat( 226) = mat( 226) + lmat( 226)
         mat( 228) = lmat( 228)
         mat( 231) = mat( 231) + lmat( 231)
         mat( 236) = mat( 236) + lmat( 236)
         mat( 238) = mat( 238) + lmat( 238)
         mat( 247) = mat( 247) + lmat( 247)
         mat( 255) = mat( 255) + lmat( 255)
         mat( 259) = mat( 259) + lmat( 259)
         mat( 270) = mat( 270) + lmat( 270)
         mat( 276) = mat( 276) + lmat( 276)
         mat( 281) = mat( 281) + lmat( 281)
         mat( 284) = mat( 284) + lmat( 284)
         mat( 295) = mat( 295) + lmat( 295)
         mat( 307) = mat( 307) + lmat( 307)
         mat( 308) = lmat( 308)
         mat( 309) = mat( 309) + lmat( 309)
         mat( 318) = mat( 318) + lmat( 318)
         mat( 327) = mat( 327) + lmat( 327)
         mat( 331) = mat( 331) + lmat( 331)
         mat( 341) = mat( 341) + lmat( 341)
         mat( 342) = mat( 342) + lmat( 342)
         mat( 343) = mat( 343) + lmat( 343)
         mat( 350) = mat( 350) + lmat( 350)
         mat( 359) = mat( 359) + lmat( 359)
         mat( 362) = mat( 362) + lmat( 362)
         mat( 364) = mat( 364) + lmat( 364)
         mat( 368) = mat( 368) + lmat( 368)
         mat( 371) = mat( 371) + lmat( 371)
         mat( 373) = mat( 373) + lmat( 373)
         mat( 374) = mat( 374) + lmat( 374)
         mat( 382) = mat( 382) + lmat( 382)
         mat( 386) = mat( 386) + lmat( 386)
         mat( 388) = mat( 388) + lmat( 388)
         mat( 389) = mat( 389) + lmat( 389)
         mat( 391) = mat( 391) + lmat( 391)
         mat( 397) = mat( 397) + lmat( 397)
         mat( 400) = lmat( 400)
         mat( 405) = lmat( 405)
         mat( 408) = mat( 408) + lmat( 408)
         mat( 413) = mat( 413) + lmat( 413)
         mat( 417) = mat( 417) + lmat( 417)
         mat( 432) = mat( 432) + lmat( 432)
         mat( 439) = lmat( 439)
         mat( 441) = lmat( 441)
         mat( 442) = lmat( 442)
         mat( 443) = mat( 443) + lmat( 443)
         mat( 444) = lmat( 444)
         mat( 445) = lmat( 445)
         mat( 447) = mat( 447) + lmat( 447)
         mat( 449) = mat( 449) + lmat( 449)
         mat( 450) = mat( 450) + lmat( 450)
         mat( 451) = mat( 451) + lmat( 451)
         mat( 454) = mat( 454) + lmat( 454)
         mat( 455) = mat( 455) + lmat( 455)
         mat( 458) = mat( 458) + lmat( 458)
         mat( 468) = mat( 468) + lmat( 468)
         mat( 469) = mat( 469) + lmat( 469)
         mat( 471) = mat( 471) + lmat( 471)
         mat( 479) = lmat( 479)
         mat( 488) = mat( 488) + lmat( 488)
         mat( 492) = mat( 492) + lmat( 492)
         mat( 497) = mat( 497) + lmat( 497)
         mat( 506) = lmat( 506)
         mat( 527) = mat( 527) + lmat( 527)
         mat(  96) = 0._r8
         mat(  97) = 0._r8
         mat( 129) = 0._r8
         mat( 130) = 0._r8
         mat( 131) = 0._r8
         mat( 139) = 0._r8
         mat( 140) = 0._r8
         mat( 144) = 0._r8
         mat( 152) = 0._r8
         mat( 159) = 0._r8
         mat( 165) = 0._r8
         mat( 167) = 0._r8
         mat( 172) = 0._r8
         mat( 185) = 0._r8
         mat( 204) = 0._r8
         mat( 209) = 0._r8
         mat( 217) = 0._r8
         mat( 234) = 0._r8
         mat( 235) = 0._r8
         mat( 242) = 0._r8
         mat( 244) = 0._r8
         mat( 245) = 0._r8
         mat( 253) = 0._r8
         mat( 256) = 0._r8
         mat( 257) = 0._r8
         mat( 258) = 0._r8
         mat( 266) = 0._r8
         mat( 267) = 0._r8
         mat( 275) = 0._r8
         mat( 278) = 0._r8
         mat( 279) = 0._r8
         mat( 290) = 0._r8
         mat( 302) = 0._r8
         mat( 311) = 0._r8
         mat( 313) = 0._r8
         mat( 324) = 0._r8
         mat( 326) = 0._r8
         mat( 333) = 0._r8
         mat( 335) = 0._r8
         mat( 336) = 0._r8
         mat( 345) = 0._r8
         mat( 346) = 0._r8
         mat( 347) = 0._r8
         mat( 349) = 0._r8
         mat( 369) = 0._r8
         mat( 375) = 0._r8
         mat( 383) = 0._r8
         mat( 385) = 0._r8
         mat( 392) = 0._r8
         mat( 393) = 0._r8
         mat( 394) = 0._r8
         mat( 395) = 0._r8
         mat( 396) = 0._r8
         mat( 403) = 0._r8
         mat( 409) = 0._r8
         mat( 410) = 0._r8
         mat( 414) = 0._r8
         mat( 415) = 0._r8
         mat( 418) = 0._r8
         mat( 420) = 0._r8
         mat( 422) = 0._r8
         mat( 425) = 0._r8
         mat( 427) = 0._r8
         mat( 428) = 0._r8
         mat( 429) = 0._r8
         mat( 430) = 0._r8
         mat( 431) = 0._r8
         mat( 433) = 0._r8
         mat( 435) = 0._r8
         mat( 436) = 0._r8
         mat( 446) = 0._r8
         mat( 448) = 0._r8
         mat( 470) = 0._r8
         mat( 480) = 0._r8
         mat( 481) = 0._r8
         mat( 482) = 0._r8
         mat( 483) = 0._r8
         mat( 484) = 0._r8
         mat( 485) = 0._r8
         mat( 487) = 0._r8
         mat( 489) = 0._r8
         mat( 491) = 0._r8
         mat( 503) = 0._r8
         mat( 504) = 0._r8
         mat( 511) = 0._r8
         mat( 524) = 0._r8
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
         mat(  18) = mat(  18) - dti
         mat(  22) = mat(  22) - dti
         mat(  24) = mat(  24) - dti
         mat(  27) = mat(  27) - dti
         mat(  30) = mat(  30) - dti
         mat(  32) = mat(  32) - dti
         mat(  34) = mat(  34) - dti
         mat(  38) = mat(  38) - dti
         mat(  43) = mat(  43) - dti
         mat(  49) = mat(  49) - dti
         mat(  55) = mat(  55) - dti
         mat(  62) = mat(  62) - dti
         mat(  69) = mat(  69) - dti
         mat(  75) = mat(  75) - dti
         mat(  80) = mat(  80) - dti
         mat(  90) = mat(  90) - dti
         mat(  99) = mat(  99) - dti
         mat( 106) = mat( 106) - dti
         mat( 112) = mat( 112) - dti
         mat( 117) = mat( 117) - dti
         mat( 126) = mat( 126) - dti
         mat( 137) = mat( 137) - dti
         mat( 148) = mat( 148) - dti
         mat( 155) = mat( 155) - dti
         mat( 164) = mat( 164) - dti
         mat( 174) = mat( 174) - dti
         mat( 184) = mat( 184) - dti
         mat( 193) = mat( 193) - dti
         mat( 201) = mat( 201) - dti
         mat( 213) = mat( 213) - dti
         mat( 231) = mat( 231) - dti
         mat( 247) = mat( 247) - dti
         mat( 270) = mat( 270) - dti
         mat( 295) = mat( 295) - dti
         mat( 318) = mat( 318) - dti
         mat( 342) = mat( 342) - dti
         mat( 371) = mat( 371) - dti
         mat( 391) = mat( 391) - dti
         mat( 413) = mat( 413) - dti
         mat( 432) = mat( 432) - dti
         mat( 451) = mat( 451) - dti
         mat( 469) = mat( 469) - dti
         mat( 492) = mat( 492) - dti
         mat( 527) = mat( 527) - dti

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
      call     nlnmat03( mat, y, rxt )
      call     nlnmat_finit( mat, lmat, dti )

      end subroutine nlnmat

      end module mo_nln_matrix
