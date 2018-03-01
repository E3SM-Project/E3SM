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
         mat(496) = -(rxt(76)*y(2) + rxt(94)*y(3) + rxt(116)*y(9) + rxt(119)*y(10) &
                      + rxt(141)*y(21) + rxt(146)*y(22) + rxt(153)*y(23) + rxt(156) &
                      *y(27) + rxt(182)*y(36))
         mat(289) = -rxt(76)*y(1)
         mat(370) = -rxt(94)*y(1)
         mat(311) = -rxt(116)*y(1)
         mat(352) = -rxt(119)*y(1)
         mat(193) = -rxt(141)*y(1)
         mat(419) = -rxt(146)*y(1)
         mat(394) = -rxt(153)*y(1)
         mat(457) = -rxt(156)*y(1)
         mat(329) = -rxt(182)*y(1)
         mat(289) = mat(289) + rxt(75)*y(4)
         mat(212) = rxt(75)*y(2)
         mat(279) = -(rxt(75)*y(4) + rxt(76)*y(1) + 4._r8*rxt(77)*y(2) + rxt(114)*y(9) &
                      + (rxt(117) + rxt(118)) * y(10) + rxt(125)*y(11) + rxt(137) &
                      *y(18) + rxt(145)*y(22) + rxt(152)*y(23) + rxt(155)*y(24) &
                      + rxt(163)*y(29) + rxt(175)*y(32) + rxt(176)*y(33) + rxt(179) &
                      *y(34) + rxt(185)*y(37) + rxt(195)*y(38) + rxt(196)*y(39) &
                      + rxt(197)*y(40) + (rxt(222) + rxt(231)) * y(52) + rxt(228) &
                      *y(54))
         mat(206) = -rxt(75)*y(2)
         mat(486) = -rxt(76)*y(2)
         mat(301) = -rxt(114)*y(2)
         mat(342) = -(rxt(117) + rxt(118)) * y(2)
         mat(426) = -rxt(125)*y(2)
         mat(177) = -rxt(137)*y(2)
         mat(409) = -rxt(145)*y(2)
         mat(384) = -rxt(152)*y(2)
         mat(64) = -rxt(155)*y(2)
         mat(245) = -rxt(163)*y(2)
         mat(469) = -rxt(175)*y(2)
         mat(151) = -rxt(176)*y(2)
         mat(161) = -rxt(179)*y(2)
         mat(223) = -rxt(185)*y(2)
         mat(82) = -rxt(195)*y(2)
         mat(95) = -rxt(196)*y(2)
         mat(58) = -rxt(197)*y(2)
         mat(77) = -(rxt(222) + rxt(231)) * y(2)
         mat(42) = -rxt(228)*y(2)
         mat(360) = (rxt(89)+rxt(90))*y(4)
         mat(206) = mat(206) + (rxt(89)+rxt(90))*y(3) + rxt(111)*y(8) + rxt(227)*y(54) &
                      + rxt(220)*y(55)
         mat(133) = rxt(111)*y(4) + rxt(112)*y(9) + rxt(113)*y(10) + rxt(224)*y(53)
         mat(301) = mat(301) + rxt(112)*y(8)
         mat(342) = mat(342) + rxt(113)*y(8)
         mat(409) = mat(409) + 2.000_r8*rxt(148)*y(22)
         mat(189) = rxt(144)*y(23)
         mat(384) = mat(384) + rxt(144)*y(21)
         mat(105) = rxt(224)*y(8) + 1.150_r8*rxt(233)*y(57)
         mat(42) = mat(42) + rxt(227)*y(4)
         mat(48) = rxt(220)*y(4)
         mat(113) = rxt(232)*y(57)
         mat(123) = 1.150_r8*rxt(233)*y(53) + rxt(232)*y(56)
         mat(364) = -((rxt(89) + rxt(90)) * y(4) + rxt(91)*y(59) + rxt(94)*y(1) &
                      + rxt(107)*y(32) + rxt(108)*y(38))
         mat(209) = -(rxt(89) + rxt(90)) * y(3)
         mat(170) = -rxt(91)*y(3)
         mat(490) = -rxt(94)*y(3)
         mat(473) = -rxt(107)*y(3)
         mat(84) = -rxt(108)*y(3)
         mat(209) = mat(209) + rxt(109)*y(58)
         mat(107) = .850_r8*rxt(233)*y(57)
         mat(53) = rxt(109)*y(4)
         mat(125) = .850_r8*rxt(233)*y(53)
         mat(205) = -(rxt(75)*y(2) + rxt(85)*y(6) + rxt(89)*y(3) + rxt(109)*y(58) &
                      + rxt(111)*y(8) + rxt(140)*y(21) + rxt(220)*y(55) + (rxt(226) &
                      + rxt(227)) * y(54) + rxt(229)*y(52))
         mat(276) = -rxt(75)*y(4)
         mat(5) = -rxt(85)*y(4)
         mat(359) = -rxt(89)*y(4)
         mat(51) = -rxt(109)*y(4)
         mat(132) = -rxt(111)*y(4)
         mat(188) = -rxt(140)*y(4)
         mat(47) = -rxt(220)*y(4)
         mat(41) = -(rxt(226) + rxt(227)) * y(4)
         mat(76) = -rxt(229)*y(4)
         mat(483) = 2.000_r8*rxt(76)*y(2) + 2.000_r8*rxt(94)*y(3) + rxt(116)*y(9) &
                      + rxt(119)*y(10) + rxt(146)*y(22) + rxt(141)*y(21) &
                      + 2.000_r8*rxt(153)*y(23) + rxt(156)*y(27) + rxt(182)*y(36)
         mat(276) = mat(276) + 2.000_r8*rxt(76)*y(1) + 2.000_r8*rxt(77)*y(2) + rxt(84) &
                      *y(6) + rxt(117)*y(10) + rxt(145)*y(22) + rxt(125)*y(11) &
                      + rxt(152)*y(23) + rxt(163)*y(29) + rxt(185)*y(37)
         mat(359) = mat(359) + 2.000_r8*rxt(94)*y(1)
         mat(205) = mat(205) + 2.000_r8*rxt(85)*y(6)
         mat(5) = mat(5) + rxt(84)*y(2) + 2.000_r8*rxt(85)*y(4)
         mat(298) = rxt(116)*y(1) + rxt(225)*y(53)
         mat(339) = rxt(119)*y(1) + rxt(117)*y(2)
         mat(406) = rxt(146)*y(1) + rxt(145)*y(2) + rxt(129)*y(13) + rxt(147)*y(23) &
                      + rxt(165)*y(29)
         mat(425) = rxt(125)*y(2) + rxt(127)*y(23)
         mat(31) = rxt(129)*y(22)
         mat(142) = rxt(133)*y(23)
         mat(188) = mat(188) + rxt(141)*y(1) + rxt(143)*y(23)
         mat(381) = 2.000_r8*rxt(153)*y(1) + rxt(152)*y(2) + rxt(147)*y(22) + rxt(127) &
                      *y(11) + rxt(133)*y(16) + rxt(143)*y(21) + 2.000_r8*rxt(154) &
                      *y(23) + rxt(159)*y(27) + rxt(166)*y(29) + rxt(183)*y(36) &
                      + rxt(187)*y(37)
         mat(445) = rxt(156)*y(1) + rxt(159)*y(23)
         mat(242) = rxt(163)*y(2) + rxt(165)*y(22) + rxt(166)*y(23) + ( &
                      + 2.000_r8*rxt(169)+2.000_r8*rxt(170))*y(29) + (rxt(191) &
                       +rxt(192))*y(37)
         mat(316) = rxt(182)*y(1) + rxt(183)*y(23)
         mat(220) = rxt(185)*y(2) + rxt(187)*y(23) + (rxt(191)+rxt(192))*y(29) &
                      + 2.000_r8*rxt(193)*y(37)
         mat(104) = rxt(225)*y(9)
         mat(7) = -(rxt(78)*y(2) + rxt(79)*y(4) + rxt(81)*y(1))
         mat(257) = -rxt(78)*y(5)
         mat(195) = -rxt(79)*y(5)
         mat(481) = -rxt(81)*y(5)
         mat(353) = rxt(89)*y(4)
         mat(195) = mat(195) + rxt(89)*y(3)
         mat(4) = -(rxt(84)*y(2) + rxt(85)*y(4))
         mat(256) = -rxt(84)*y(6)
         mat(194) = -rxt(85)*y(6)
         mat(480) = rxt(81)*y(5)
         mat(256) = mat(256) + rxt(78)*y(5)
         mat(194) = mat(194) + rxt(79)*y(5)
         mat(6) = rxt(81)*y(1) + rxt(78)*y(2) + rxt(79)*y(4)
         mat(131) = -(rxt(111)*y(4) + rxt(112)*y(9) + rxt(113)*y(10) + rxt(224)*y(53))
         mat(203) = -rxt(111)*y(8)
         mat(293) = -rxt(112)*y(8)
         mat(335) = -rxt(113)*y(8)
         mat(103) = -rxt(224)*y(8)
         mat(270) = rxt(228)*y(54) + rxt(110)*y(58)
         mat(203) = mat(203) + rxt(226)*y(54)
         mat(75) = 1.100_r8*rxt(234)*y(57)
         mat(40) = rxt(228)*y(2) + rxt(226)*y(4)
         mat(111) = .200_r8*rxt(232)*y(57)
         mat(50) = rxt(110)*y(2)
         mat(121) = 1.100_r8*rxt(234)*y(52) + .200_r8*rxt(232)*y(56)
         mat(302) = -(rxt(112)*y(8) + rxt(114)*y(2) + rxt(115)*y(23) + rxt(116)*y(1) &
                      + rxt(124)*y(11) + rxt(132)*y(16) + rxt(167)*y(29) + rxt(188) &
                      *y(37) + rxt(225)*y(53))
         mat(134) = -rxt(112)*y(9)
         mat(280) = -rxt(114)*y(9)
         mat(385) = -rxt(115)*y(9)
         mat(487) = -rxt(116)*y(9)
         mat(427) = -rxt(124)*y(9)
         mat(143) = -rxt(132)*y(9)
         mat(246) = -rxt(167)*y(9)
         mat(224) = -rxt(188)*y(9)
         mat(106) = -rxt(225)*y(9)
         mat(280) = mat(280) + rxt(117)*y(10)
         mat(207) = rxt(111)*y(8) + rxt(109)*y(58)
         mat(134) = mat(134) + rxt(111)*y(4)
         mat(343) = rxt(117)*y(2)
         mat(52) = rxt(109)*y(4)
         mat(345) = -(rxt(113)*y(8) + (rxt(117) + rxt(118)) * y(2) + rxt(119)*y(1) &
                      + rxt(120)*y(11) + rxt(122)*y(22) + rxt(128)*y(23) + rxt(168) &
                      *y(29) + rxt(189)*y(37))
         mat(135) = -rxt(113)*y(10)
         mat(282) = -(rxt(117) + rxt(118)) * y(10)
         mat(489) = -rxt(119)*y(10)
         mat(429) = -rxt(120)*y(10)
         mat(412) = -rxt(122)*y(10)
         mat(387) = -rxt(128)*y(10)
         mat(248) = -rxt(168)*y(10)
         mat(226) = -rxt(189)*y(10)
         mat(489) = mat(489) + rxt(116)*y(9)
         mat(282) = mat(282) + rxt(114)*y(9) + rxt(125)*y(11)
         mat(304) = rxt(116)*y(1) + rxt(114)*y(2) + 2.000_r8*rxt(124)*y(11) + rxt(132) &
                      *y(16) + rxt(115)*y(23) + rxt(167)*y(29) + rxt(188)*y(37)
         mat(412) = mat(412) + rxt(126)*y(11) + rxt(129)*y(13)
         mat(429) = mat(429) + rxt(125)*y(2) + 2.000_r8*rxt(124)*y(9) + rxt(126)*y(22) &
                      + rxt(127)*y(23)
         mat(32) = rxt(129)*y(22)
         mat(144) = rxt(132)*y(9)
         mat(387) = mat(387) + rxt(115)*y(9) + rxt(127)*y(11)
         mat(248) = mat(248) + rxt(167)*y(9)
         mat(226) = mat(226) + rxt(188)*y(9)
         mat(415) = -(rxt(122)*y(10) + rxt(123)*y(12) + rxt(126)*y(11) + rxt(129) &
                      *y(13) + rxt(134)*y(17) + rxt(136)*y(18) + rxt(145)*y(2) + rxt(146) &
                      *y(1) + rxt(147)*y(23) + (4._r8*rxt(148) + 4._r8*rxt(149) &
                      ) * y(22) + rxt(151)*y(24) + (rxt(164) + rxt(165)) * y(29) &
                      + rxt(174)*y(32) + rxt(178)*y(33) + rxt(180)*y(34) + rxt(186) &
                      *y(37) + rxt(194)*y(38))
         mat(348) = -rxt(122)*y(22)
         mat(89) = -rxt(123)*y(22)
         mat(432) = -rxt(126)*y(22)
         mat(34) = -rxt(129)*y(22)
         mat(28) = -rxt(134)*y(22)
         mat(182) = -rxt(136)*y(22)
         mat(285) = -rxt(145)*y(22)
         mat(492) = -rxt(146)*y(22)
         mat(390) = -rxt(147)*y(22)
         mat(66) = -rxt(151)*y(22)
         mat(251) = -(rxt(164) + rxt(165)) * y(22)
         mat(475) = -rxt(174)*y(22)
         mat(152) = -rxt(178)*y(22)
         mat(163) = -rxt(180)*y(22)
         mat(229) = -rxt(186)*y(22)
         mat(85) = -rxt(194)*y(22)
         mat(492) = mat(492) + rxt(141)*y(21) + rxt(153)*y(23)
         mat(285) = mat(285) + rxt(137)*y(18) + rxt(152)*y(23) + rxt(155)*y(24) &
                      + rxt(175)*y(32) + rxt(176)*y(33) + rxt(195)*y(38) + rxt(196) &
                      *y(39)
         mat(366) = 2.000_r8*rxt(91)*y(59) + rxt(107)*y(32) + rxt(108)*y(38)
         mat(307) = rxt(115)*y(23)
         mat(432) = mat(432) + rxt(127)*y(23)
         mat(182) = mat(182) + rxt(137)*y(2)
         mat(192) = rxt(141)*y(1) + 2.000_r8*rxt(142)*y(23)
         mat(390) = mat(390) + rxt(153)*y(1) + rxt(152)*y(2) + rxt(115)*y(9) &
                      + rxt(127)*y(11) + 2.000_r8*rxt(142)*y(21) + rxt(160)*y(27)
         mat(66) = mat(66) + rxt(155)*y(2)
         mat(171) = 2.000_r8*rxt(91)*y(3)
         mat(453) = rxt(160)*y(23)
         mat(475) = mat(475) + rxt(175)*y(2) + rxt(107)*y(3)
         mat(152) = mat(152) + rxt(176)*y(2)
         mat(85) = mat(85) + rxt(195)*y(2) + rxt(108)*y(3)
         mat(97) = rxt(196)*y(2)
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
         mat(433) = -(rxt(120)*y(10) + rxt(124)*y(9) + rxt(125)*y(2) + rxt(126)*y(22) &
                      + rxt(127)*y(23) + rxt(135)*y(18))
         mat(349) = -rxt(120)*y(11)
         mat(308) = -rxt(124)*y(11)
         mat(286) = -rxt(125)*y(11)
         mat(416) = -rxt(126)*y(11)
         mat(391) = -rxt(127)*y(11)
         mat(183) = -rxt(135)*y(11)
         mat(493) = rxt(119)*y(10)
         mat(286) = mat(286) + rxt(118)*y(10) + rxt(179)*y(34) + rxt(197)*y(40)
         mat(349) = mat(349) + rxt(119)*y(1) + rxt(118)*y(2)
         mat(416) = mat(416) + rxt(123)*y(12) + rxt(180)*y(34)
         mat(90) = rxt(123)*y(22)
         mat(454) = rxt(181)*y(34)
         mat(164) = rxt(179)*y(2) + rxt(180)*y(22) + rxt(181)*y(27)
         mat(61) = rxt(197)*y(2)
         mat(86) = -(rxt(123)*y(22))
         mat(399) = -rxt(123)*y(12)
         mat(333) = rxt(122)*y(22)
         mat(399) = mat(399) + rxt(122)*y(10)
         mat(421) = rxt(135)*y(18)
         mat(173) = rxt(135)*y(11)
         mat(460) = (rxt(206)+rxt(211)+rxt(217))*y(34)
         mat(156) = (rxt(206)+rxt(211)+rxt(217))*y(32)
         mat(29) = -(rxt(129)*y(22))
         mat(396) = -rxt(129)*y(13)
         mat(331) = rxt(128)*y(23)
         mat(372) = rxt(128)*y(10)
         mat(330) = rxt(120)*y(11)
         mat(420) = rxt(120)*y(10)
         mat(138) = -(rxt(132)*y(9) + rxt(133)*y(23))
         mat(294) = -rxt(132)*y(16)
         mat(376) = -rxt(133)*y(16)
         mat(400) = rxt(134)*y(17)
         mat(24) = rxt(134)*y(22)
         mat(23) = -(rxt(134)*y(22))
         mat(395) = -rxt(134)*y(17)
         mat(137) = rxt(133)*y(23)
         mat(371) = rxt(133)*y(16)
         mat(175) = -(rxt(135)*y(11) + rxt(136)*y(22) + rxt(137)*y(2) + rxt(161)*y(27) &
                      + rxt(184)*y(36))
         mat(423) = -rxt(135)*y(18)
         mat(404) = -rxt(136)*y(18)
         mat(274) = -rxt(137)*y(18)
         mat(443) = -rxt(161)*y(18)
         mat(314) = -rxt(184)*y(18)
         mat(296) = rxt(132)*y(16)
         mat(140) = rxt(132)*y(9)
         mat(187) = -(rxt(140)*y(4) + rxt(141)*y(1) + (rxt(142) + rxt(143) + rxt(144) &
                      ) * y(23))
         mat(204) = -rxt(140)*y(21)
         mat(482) = -rxt(141)*y(21)
         mat(380) = -(rxt(142) + rxt(143) + rxt(144)) * y(21)
         mat(275) = rxt(145)*y(22)
         mat(405) = rxt(145)*y(2) + rxt(136)*y(18)
         mat(176) = rxt(136)*y(22)
         mat(389) = -(rxt(115)*y(9) + rxt(127)*y(11) + rxt(128)*y(10) + rxt(133)*y(16) &
                      + (rxt(142) + rxt(143) + rxt(144)) * y(21) + rxt(147)*y(22) &
                      + rxt(152)*y(2) + rxt(153)*y(1) + 4._r8*rxt(154)*y(23) + (rxt(159) &
                      + rxt(160)) * y(27) + rxt(166)*y(29) + rxt(183)*y(36) + rxt(187) &
                      *y(37))
         mat(306) = -rxt(115)*y(23)
         mat(431) = -rxt(127)*y(23)
         mat(347) = -rxt(128)*y(23)
         mat(145) = -rxt(133)*y(23)
         mat(191) = -(rxt(142) + rxt(143) + rxt(144)) * y(23)
         mat(414) = -rxt(147)*y(23)
         mat(284) = -rxt(152)*y(23)
         mat(491) = -rxt(153)*y(23)
         mat(452) = -(rxt(159) + rxt(160)) * y(23)
         mat(250) = -rxt(166)*y(23)
         mat(324) = -rxt(183)*y(23)
         mat(228) = -rxt(187)*y(23)
         mat(491) = mat(491) + rxt(146)*y(22)
         mat(284) = mat(284) + rxt(137)*y(18) + rxt(155)*y(24)
         mat(210) = rxt(140)*y(21)
         mat(306) = mat(306) + rxt(132)*y(16)
         mat(414) = mat(414) + rxt(146)*y(1) + rxt(126)*y(11) + rxt(151)*y(24) &
                      + rxt(164)*y(29) + rxt(186)*y(37)
         mat(431) = mat(431) + rxt(126)*y(22) + rxt(135)*y(18)
         mat(145) = mat(145) + rxt(132)*y(9)
         mat(181) = rxt(137)*y(2) + rxt(135)*y(11) + rxt(161)*y(27) + rxt(184)*y(36)
         mat(191) = mat(191) + rxt(140)*y(4)
         mat(65) = rxt(155)*y(2) + rxt(151)*y(22) + rxt(158)*y(27)
         mat(452) = mat(452) + rxt(161)*y(18) + rxt(158)*y(24)
         mat(250) = mat(250) + rxt(164)*y(22)
         mat(324) = mat(324) + rxt(184)*y(18)
         mat(228) = mat(228) + rxt(186)*y(22)
         mat(62) = -(rxt(151)*y(22) + rxt(155)*y(2) + rxt(158)*y(27))
         mat(397) = -rxt(151)*y(24)
         mat(262) = -rxt(155)*y(24)
         mat(438) = -rxt(158)*y(24)
         mat(397) = mat(397) + 2.000_r8*rxt(149)*y(22)
         mat(373) = 2.000_r8*rxt(154)*y(23)
         mat(167) = -(rxt(91)*y(3))
         mat(356) = -rxt(91)*y(59)
         mat(403) = 2.000_r8*rxt(148)*y(22) + rxt(123)*y(12) + rxt(129)*y(13) &
                      + rxt(134)*y(17) + rxt(136)*y(18) + rxt(147)*y(23) + rxt(151) &
                      *y(24) + rxt(174)*y(32) + rxt(178)*y(33) + rxt(194)*y(38)
         mat(87) = rxt(123)*y(22)
         mat(30) = rxt(129)*y(22)
         mat(25) = rxt(134)*y(22)
         mat(174) = rxt(136)*y(22)
         mat(186) = rxt(144)*y(23)
         mat(378) = rxt(147)*y(22) + rxt(144)*y(21)
         mat(63) = rxt(151)*y(22)
         mat(464) = rxt(174)*y(22) + (rxt(207)+rxt(212)+rxt(218))*y(33) + (rxt(208) &
                       +rxt(219))*y(39)
         mat(149) = rxt(178)*y(22) + (rxt(207)+rxt(212)+rxt(218))*y(32)
         mat(81) = rxt(194)*y(22)
         mat(93) = (rxt(208)+rxt(219))*y(32)
         mat(455) = -(rxt(156)*y(1) + rxt(158)*y(24) + (rxt(159) + rxt(160)) * y(23) &
                      + rxt(161)*y(18) + rxt(177)*y(33) + rxt(181)*y(34))
         mat(494) = -rxt(156)*y(27)
         mat(67) = -rxt(158)*y(27)
         mat(392) = -(rxt(159) + rxt(160)) * y(27)
         mat(184) = -rxt(161)*y(27)
         mat(153) = -rxt(177)*y(27)
         mat(165) = -rxt(181)*y(27)
         mat(287) = rxt(163)*y(29) + rxt(175)*y(32)
         mat(368) = rxt(107)*y(32)
         mat(309) = rxt(167)*y(29)
         mat(417) = rxt(164)*y(29) + rxt(174)*y(32)
         mat(253) = rxt(163)*y(2) + rxt(167)*y(9) + rxt(164)*y(22) + ( &
                      + 4.000_r8*rxt(169)+2.000_r8*rxt(171))*y(29) + rxt(191)*y(37)
         mat(477) = rxt(175)*y(2) + rxt(107)*y(3) + rxt(174)*y(22)
         mat(231) = rxt(191)*y(29)
         mat(437) = rxt(181)*y(34)
         mat(236) = 2.000_r8*rxt(170)*y(29)
         mat(458) = (rxt(207)+rxt(212)+rxt(218))*y(33) + (rxt(206)+rxt(211)+rxt(217)) &
                      *y(34)
         mat(147) = (rxt(207)+rxt(212)+rxt(218))*y(32)
         mat(155) = rxt(181)*y(27) + (rxt(206)+rxt(211)+rxt(217))*y(32)
         mat(244) = -(rxt(163)*y(2) + (rxt(164) + rxt(165)) * y(22) + rxt(166)*y(23) &
                      + rxt(167)*y(9) + rxt(168)*y(10) + (4._r8*rxt(169) + 4._r8*rxt(170) &
                      + 4._r8*rxt(171) + 4._r8*rxt(172)) * y(29) + (rxt(190) + rxt(191) &
                      + rxt(192)) * y(37))
         mat(278) = -rxt(163)*y(29)
         mat(408) = -(rxt(164) + rxt(165)) * y(29)
         mat(383) = -rxt(166)*y(29)
         mat(300) = -rxt(167)*y(29)
         mat(341) = -rxt(168)*y(29)
         mat(222) = -(rxt(190) + rxt(191) + rxt(192)) * y(29)
         mat(485) = rxt(156)*y(27)
         mat(278) = mat(278) + rxt(176)*y(33) + rxt(179)*y(34)
         mat(408) = mat(408) + rxt(178)*y(33)
         mat(383) = mat(383) + rxt(160)*y(27)
         mat(446) = rxt(156)*y(1) + rxt(160)*y(23) + rxt(177)*y(33)
         mat(150) = rxt(176)*y(2) + rxt(178)*y(22) + rxt(177)*y(27)
         mat(160) = rxt(179)*y(2)
         mat(235) = 2.000_r8*rxt(171)*y(29) + rxt(190)*y(37)
         mat(213) = rxt(190)*y(29)
         mat(234) = 2.000_r8*rxt(172)*y(29)
         mat(478) = -(rxt(107)*y(3) + rxt(174)*y(22) + rxt(175)*y(2) + (rxt(206) &
                      + rxt(211) + rxt(217)) * y(34) + (rxt(207) + rxt(212) + rxt(218) &
                      ) * y(33) + (rxt(208) + rxt(219)) * y(39))
         mat(369) = -rxt(107)*y(32)
         mat(418) = -rxt(174)*y(32)
         mat(288) = -rxt(175)*y(32)
         mat(166) = -(rxt(206) + rxt(211) + rxt(217)) * y(32)
         mat(154) = -(rxt(207) + rxt(212) + rxt(218)) * y(32)
         mat(99) = -(rxt(208) + rxt(219)) * y(32)
         mat(418) = mat(418) + rxt(165)*y(29)
         mat(185) = rxt(161)*y(27)
         mat(393) = rxt(159)*y(27)
         mat(68) = rxt(158)*y(27)
         mat(456) = rxt(161)*y(18) + rxt(159)*y(23) + rxt(158)*y(24) + rxt(177)*y(33)
         mat(254) = rxt(165)*y(22)
         mat(154) = mat(154) + rxt(177)*y(27)
         mat(148) = -(rxt(176)*y(2) + rxt(177)*y(27) + rxt(178)*y(22) + (rxt(207) &
                      + rxt(212) + rxt(218)) * y(32))
         mat(271) = -rxt(176)*y(33)
         mat(440) = -rxt(177)*y(33)
         mat(401) = -rxt(178)*y(33)
         mat(462) = -(rxt(207) + rxt(212) + rxt(218)) * y(33)
         mat(401) = mat(401) + rxt(180)*y(34)
         mat(377) = rxt(166)*y(29)
         mat(238) = rxt(166)*y(23)
         mat(157) = rxt(180)*y(22)
         mat(158) = -(rxt(179)*y(2) + rxt(180)*y(22) + rxt(181)*y(27) + (rxt(206) &
                      + rxt(211) + rxt(217)) * y(32))
         mat(272) = -rxt(179)*y(34)
         mat(402) = -rxt(180)*y(34)
         mat(441) = -rxt(181)*y(34)
         mat(463) = -(rxt(206) + rxt(211) + rxt(217)) * y(34)
         mat(336) = rxt(168)*y(29)
         mat(239) = rxt(168)*y(10)
         mat(237) = rxt(192)*y(37)
         mat(459) = (rxt(208)+rxt(219))*y(39)
         mat(214) = rxt(192)*y(29)
         mat(91) = (rxt(208)+rxt(219))*y(32)
         mat(321) = -(rxt(182)*y(1) + rxt(183)*y(23) + rxt(184)*y(18))
         mat(488) = -rxt(182)*y(36)
         mat(386) = -rxt(183)*y(36)
         mat(178) = -rxt(184)*y(36)
         mat(281) = rxt(185)*y(37) + rxt(195)*y(38)
         mat(362) = rxt(108)*y(38)
         mat(303) = rxt(188)*y(37)
         mat(411) = rxt(186)*y(37) + rxt(194)*y(38)
         mat(247) = (rxt(190)+rxt(191))*y(37)
         mat(225) = rxt(185)*y(2) + rxt(188)*y(9) + rxt(186)*y(22) + (rxt(190) &
                       +rxt(191))*y(29) + 4.000_r8*rxt(193)*y(37)
         mat(83) = rxt(195)*y(2) + rxt(108)*y(3) + rxt(194)*y(22)
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
         mat(221) = -(rxt(185)*y(2) + rxt(186)*y(22) + rxt(187)*y(23) + rxt(188)*y(9) &
                      + rxt(189)*y(10) + (rxt(190) + rxt(191) + rxt(192)) * y(29) &
                      + 4._r8*rxt(193)*y(37))
         mat(277) = -rxt(185)*y(37)
         mat(407) = -rxt(186)*y(37)
         mat(382) = -rxt(187)*y(37)
         mat(299) = -rxt(188)*y(37)
         mat(340) = -rxt(189)*y(37)
         mat(243) = -(rxt(190) + rxt(191) + rxt(192)) * y(37)
         mat(484) = rxt(182)*y(36)
         mat(277) = mat(277) + rxt(196)*y(39) + rxt(197)*y(40)
         mat(317) = rxt(182)*y(1)
         mat(94) = rxt(196)*y(2)
         mat(57) = rxt(197)*y(2)
         mat(80) = -(rxt(108)*y(3) + rxt(194)*y(22) + rxt(195)*y(2))
         mat(354) = -rxt(108)*y(38)
         mat(398) = -rxt(194)*y(38)
         mat(264) = -rxt(195)*y(38)
         mat(172) = rxt(184)*y(36)
         mat(374) = rxt(183)*y(36)
         mat(312) = rxt(184)*y(18) + rxt(183)*y(23)
         mat(92) = -(rxt(196)*y(2) + (rxt(208) + rxt(219)) * y(32))
         mat(266) = -rxt(196)*y(39)
         mat(461) = -(rxt(208) + rxt(219)) * y(39)
         mat(375) = rxt(187)*y(37)
         mat(217) = rxt(187)*y(23)
         mat(54) = -(rxt(197)*y(2))
         mat(261) = -rxt(197)*y(40)
         mat(332) = rxt(189)*y(37)
         mat(215) = rxt(189)*y(10)
         mat(71) = -((rxt(222) + rxt(231)) * y(2) + rxt(229)*y(4) + rxt(234)*y(57))
         mat(263) = -(rxt(222) + rxt(231)) * y(52)
         mat(199) = -rxt(229)*y(52)
         mat(117) = -rxt(234)*y(52)
         mat(100) = -(rxt(224)*y(8) + rxt(225)*y(9) + rxt(233)*y(57))
         mat(128) = -rxt(224)*y(53)
         mat(290) = -rxt(225)*y(53)
         mat(118) = -rxt(233)*y(53)
         mat(200) = rxt(229)*y(52) + rxt(226)*y(54) + rxt(220)*y(55)
         mat(72) = rxt(229)*y(4)
         mat(38) = rxt(226)*y(4)
         mat(44) = rxt(220)*y(4)
         mat(36) = -((rxt(226) + rxt(227)) * y(4) + rxt(228)*y(2))
         mat(196) = -(rxt(226) + rxt(227)) * y(54)
         mat(258) = -rxt(228)*y(54)
         mat(43) = -(rxt(220)*y(4))
         mat(197) = -rxt(220)*y(55)
         mat(259) = rxt(231)*y(52) + rxt(228)*y(54)
         mat(69) = rxt(231)*y(2)
         mat(37) = rxt(228)*y(2)
         mat(109) = -(rxt(232)*y(57))
         mat(119) = -rxt(232)*y(56)
         mat(268) = rxt(222)*y(52)
         mat(201) = rxt(227)*y(54)
         mat(129) = rxt(224)*y(53)
         mat(291) = rxt(225)*y(53)
         mat(73) = rxt(222)*y(2)
         mat(101) = rxt(224)*y(8) + rxt(225)*y(9)
         mat(39) = rxt(227)*y(4)
         mat(49) = -(rxt(109)*y(4) + rxt(110)*y(2))
         mat(198) = -rxt(109)*y(58)
         mat(260) = -rxt(110)*y(58)
         mat(260) = mat(260) + rxt(222)*y(52)
         mat(70) = rxt(222)*y(2) + .900_r8*rxt(234)*y(57)
         mat(108) = .800_r8*rxt(232)*y(57)
         mat(116) = .900_r8*rxt(234)*y(52) + .800_r8*rxt(232)*y(56)
         mat(120) = -(rxt(232)*y(56) + rxt(233)*y(53) + rxt(234)*y(52))
         mat(110) = -rxt(232)*y(57)
         mat(102) = -rxt(233)*y(57)
         mat(74) = -rxt(234)*y(57)
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
         mat( 17) = lmat( 17)
         mat( 18) = lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = mat( 23) + lmat( 23)
         mat( 26) = lmat( 26)
         mat( 27) = lmat( 27)
         mat( 28) = mat( 28) + lmat( 28)
         mat( 29) = mat( 29) + lmat( 29)
         mat( 32) = mat( 32) + lmat( 32)
         mat( 33) = lmat( 33)
         mat( 34) = mat( 34) + lmat( 34)
         mat( 35) = lmat( 35)
         mat( 36) = mat( 36) + lmat( 36)
         mat( 43) = mat( 43) + lmat( 43)
         mat( 44) = mat( 44) + lmat( 44)
         mat( 45) = lmat( 45)
         mat( 46) = lmat( 46)
         mat( 49) = mat( 49) + lmat( 49)
         mat( 54) = mat( 54) + lmat( 54)
         mat( 55) = lmat( 55)
         mat( 56) = lmat( 56)
         mat( 57) = mat( 57) + lmat( 57)
         mat( 59) = lmat( 59)
         mat( 60) = lmat( 60)
         mat( 61) = mat( 61) + lmat( 61)
         mat( 62) = mat( 62) + lmat( 62)
         mat( 66) = mat( 66) + lmat( 66)
         mat( 71) = mat( 71) + lmat( 71)
         mat( 80) = mat( 80) + lmat( 80)
         mat( 86) = mat( 86) + lmat( 86)
         mat( 88) = lmat( 88)
         mat( 89) = mat( 89) + lmat( 89)
         mat( 92) = mat( 92) + lmat( 92)
         mat( 96) = lmat( 96)
         mat( 97) = mat( 97) + lmat( 97)
         mat( 100) = mat( 100) + lmat( 100)
         mat( 101) = mat( 101) + lmat( 101)
         mat( 106) = mat( 106) + lmat( 106)
         mat( 109) = mat( 109) + lmat( 109)
         mat( 120) = mat( 120) + lmat( 120)
         mat( 126) = lmat( 126)
         mat( 130) = lmat( 130)
         mat( 131) = mat( 131) + lmat( 131)
         mat( 138) = mat( 138) + lmat( 138)
         mat( 148) = mat( 148) + lmat( 148)
         mat( 152) = mat( 152) + lmat( 152)
         mat( 153) = mat( 153) + lmat( 153)
         mat( 156) = mat( 156) + lmat( 156)
         mat( 157) = mat( 157) + lmat( 157)
         mat( 158) = mat( 158) + lmat( 158)
         mat( 160) = mat( 160) + lmat( 160)
         mat( 162) = lmat( 162)
         mat( 164) = mat( 164) + lmat( 164)
         mat( 165) = mat( 165) + lmat( 165)
         mat( 167) = mat( 167) + lmat( 167)
         mat( 168) = lmat( 168)
         mat( 169) = lmat( 169)
         mat( 170) = mat( 170) + lmat( 170)
         mat( 171) = mat( 171) + lmat( 171)
         mat( 175) = mat( 175) + lmat( 175)
         mat( 176) = mat( 176) + lmat( 176)
         mat( 187) = mat( 187) + lmat( 187)
         mat( 197) = mat( 197) + lmat( 197)
         mat( 200) = mat( 200) + lmat( 200)
         mat( 202) = lmat( 202)
         mat( 205) = mat( 205) + lmat( 205)
         mat( 206) = mat( 206) + lmat( 206)
         mat( 209) = mat( 209) + lmat( 209)
         mat( 221) = mat( 221) + lmat( 221)
         mat( 223) = mat( 223) + lmat( 223)
         mat( 225) = mat( 225) + lmat( 225)
         mat( 244) = mat( 244) + lmat( 244)
         mat( 245) = mat( 245) + lmat( 245)
         mat( 253) = mat( 253) + lmat( 253)
         mat( 259) = mat( 259) + lmat( 259)
         mat( 269) = lmat( 269)
         mat( 279) = mat( 279) + lmat( 279)
         mat( 291) = mat( 291) + lmat( 291)
         mat( 292) = lmat( 292)
         mat( 293) = mat( 293) + lmat( 293)
         mat( 301) = mat( 301) + lmat( 301)
         mat( 302) = mat( 302) + lmat( 302)
         mat( 321) = mat( 321) + lmat( 321)
         mat( 342) = mat( 342) + lmat( 342)
         mat( 343) = mat( 343) + lmat( 343)
         mat( 345) = mat( 345) + lmat( 345)
         mat( 355) = lmat( 355)
         mat( 357) = lmat( 357)
         mat( 358) = lmat( 358)
         mat( 359) = mat( 359) + lmat( 359)
         mat( 360) = mat( 360) + lmat( 360)
         mat( 361) = lmat( 361)
         mat( 362) = mat( 362) + lmat( 362)
         mat( 364) = mat( 364) + lmat( 364)
         mat( 365) = lmat( 365)
         mat( 366) = mat( 366) + lmat( 366)
         mat( 368) = mat( 368) + lmat( 368)
         mat( 389) = mat( 389) + lmat( 389)
         mat( 400) = mat( 400) + lmat( 400)
         mat( 403) = mat( 403) + lmat( 403)
         mat( 405) = mat( 405) + lmat( 405)
         mat( 411) = mat( 411) + lmat( 411)
         mat( 414) = mat( 414) + lmat( 414)
         mat( 415) = mat( 415) + lmat( 415)
         mat( 417) = mat( 417) + lmat( 417)
         mat( 425) = mat( 425) + lmat( 425)
         mat( 426) = mat( 426) + lmat( 426)
         mat( 427) = mat( 427) + lmat( 427)
         mat( 429) = mat( 429) + lmat( 429)
         mat( 433) = mat( 433) + lmat( 433)
         mat( 439) = lmat( 439)
         mat( 444) = lmat( 444)
         mat( 452) = mat( 452) + lmat( 452)
         mat( 455) = mat( 455) + lmat( 455)
         mat( 456) = mat( 456) + lmat( 456)
         mat( 465) = lmat( 465)
         mat( 477) = mat( 477) + lmat( 477)
         mat( 478) = mat( 478) + lmat( 478)
         mat( 480) = mat( 480) + lmat( 480)
         mat( 483) = mat( 483) + lmat( 483)
         mat( 486) = mat( 486) + lmat( 486)
         mat( 490) = mat( 490) + lmat( 490)
         mat( 496) = mat( 496) + lmat( 496)
         mat( 78) = 0._r8
         mat( 79) = 0._r8
         mat( 98) = 0._r8
         mat( 112) = 0._r8
         mat( 114) = 0._r8
         mat( 115) = 0._r8
         mat( 122) = 0._r8
         mat( 124) = 0._r8
         mat( 127) = 0._r8
         mat( 136) = 0._r8
         mat( 139) = 0._r8
         mat( 141) = 0._r8
         mat( 146) = 0._r8
         mat( 159) = 0._r8
         mat( 179) = 0._r8
         mat( 180) = 0._r8
         mat( 190) = 0._r8
         mat( 208) = 0._r8
         mat( 211) = 0._r8
         mat( 216) = 0._r8
         mat( 218) = 0._r8
         mat( 219) = 0._r8
         mat( 227) = 0._r8
         mat( 230) = 0._r8
         mat( 232) = 0._r8
         mat( 233) = 0._r8
         mat( 240) = 0._r8
         mat( 241) = 0._r8
         mat( 249) = 0._r8
         mat( 252) = 0._r8
         mat( 255) = 0._r8
         mat( 265) = 0._r8
         mat( 267) = 0._r8
         mat( 273) = 0._r8
         mat( 283) = 0._r8
         mat( 295) = 0._r8
         mat( 297) = 0._r8
         mat( 305) = 0._r8
         mat( 310) = 0._r8
         mat( 313) = 0._r8
         mat( 315) = 0._r8
         mat( 318) = 0._r8
         mat( 319) = 0._r8
         mat( 320) = 0._r8
         mat( 322) = 0._r8
         mat( 323) = 0._r8
         mat( 325) = 0._r8
         mat( 326) = 0._r8
         mat( 327) = 0._r8
         mat( 328) = 0._r8
         mat( 334) = 0._r8
         mat( 337) = 0._r8
         mat( 338) = 0._r8
         mat( 344) = 0._r8
         mat( 346) = 0._r8
         mat( 350) = 0._r8
         mat( 351) = 0._r8
         mat( 363) = 0._r8
         mat( 367) = 0._r8
         mat( 379) = 0._r8
         mat( 388) = 0._r8
         mat( 410) = 0._r8
         mat( 413) = 0._r8
         mat( 422) = 0._r8
         mat( 424) = 0._r8
         mat( 428) = 0._r8
         mat( 430) = 0._r8
         mat( 434) = 0._r8
         mat( 435) = 0._r8
         mat( 436) = 0._r8
         mat( 442) = 0._r8
         mat( 447) = 0._r8
         mat( 448) = 0._r8
         mat( 449) = 0._r8
         mat( 450) = 0._r8
         mat( 451) = 0._r8
         mat( 466) = 0._r8
         mat( 467) = 0._r8
         mat( 468) = 0._r8
         mat( 470) = 0._r8
         mat( 471) = 0._r8
         mat( 472) = 0._r8
         mat( 474) = 0._r8
         mat( 476) = 0._r8
         mat( 479) = 0._r8
         mat( 495) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 4) = mat( 4) - dti
         mat( 7) = mat( 7) - dti
         mat( 9) = mat( 9) - dti
         mat( 12) = mat( 12) - dti
         mat( 14) = mat( 14) - dti
         mat( 17) = mat( 17) - dti
         mat( 23) = mat( 23) - dti
         mat( 29) = mat( 29) - dti
         mat( 36) = mat( 36) - dti
         mat( 43) = mat( 43) - dti
         mat( 49) = mat( 49) - dti
         mat( 54) = mat( 54) - dti
         mat( 62) = mat( 62) - dti
         mat( 71) = mat( 71) - dti
         mat( 80) = mat( 80) - dti
         mat( 86) = mat( 86) - dti
         mat( 92) = mat( 92) - dti
         mat( 100) = mat( 100) - dti
         mat( 109) = mat( 109) - dti
         mat( 120) = mat( 120) - dti
         mat( 131) = mat( 131) - dti
         mat( 138) = mat( 138) - dti
         mat( 148) = mat( 148) - dti
         mat( 158) = mat( 158) - dti
         mat( 167) = mat( 167) - dti
         mat( 175) = mat( 175) - dti
         mat( 187) = mat( 187) - dti
         mat( 205) = mat( 205) - dti
         mat( 221) = mat( 221) - dti
         mat( 244) = mat( 244) - dti
         mat( 279) = mat( 279) - dti
         mat( 302) = mat( 302) - dti
         mat( 321) = mat( 321) - dti
         mat( 345) = mat( 345) - dti
         mat( 364) = mat( 364) - dti
         mat( 389) = mat( 389) - dti
         mat( 415) = mat( 415) - dti
         mat( 433) = mat( 433) - dti
         mat( 455) = mat( 455) - dti
         mat( 478) = mat( 478) - dti
         mat( 496) = mat( 496) - dti
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
