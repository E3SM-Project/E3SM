
      module mo_lin_matrix

      private
      public :: linmat

      contains

      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(in)    ::  het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) ::  mat(nzcnt)

         mat(469) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(527) = -( rxt(49) + rxt(50) + rxt(51) + rxt(62) + rxt(63) + rxt(64) &
                 + het_rates(2) )
         mat(238) = rxt(1) + 2.000_r8*rxt(2) + rxt(55) + rxt(56) + rxt(57) &
                      + 2.000_r8*rxt(60) + rxt(67) + rxt(68) + rxt(69) + 2.000_r8*rxt(72)
         mat(471) = rxt(4)
         mat(327) = rxt(6)
         mat(350) = rxt(8)
         mat(48) = rxt(10)
         mat(397) = rxt(12)
         mat(197) = rxt(21)
         mat(281) = rxt(24)
         mat(26) = rxt(25)
         mat(259) = rxt(32)
         mat(454) = rxt(88)

         mat(451) = -( rxt(88) + rxt(92)*y(7) + rxt(93)*y(7) + rxt(95)*y(43) &
                      + rxt(96)*y(44) + rxt(97)*y(45) + rxt(98)*y(46) + rxt(99)*y(47) &
                      + rxt(100)*y(42) + rxt(101)*y(50) + rxt(102)*y(49) + rxt(103)*y(15) &
                      + rxt(104)*y(15) + rxt(105)*y(15) + rxt(106)*y(20) + het_rates(3) )
         mat(236) = rxt(1)
         mat(468) = rxt(3)
         mat(196) = rxt(20)

         mat(231) = -( rxt(1) + rxt(2) + rxt(53) + rxt(55) + rxt(56) + rxt(57) + rxt(60) &
                      + rxt(65) + rxt(67) + rxt(68) + rxt(69) + rxt(72) + het_rates(4) )
         mat(458) = rxt(4)
         mat(386) = rxt(13)
         mat(31) = rxt(83)
         mat(28) = rxt(86) + rxt(87)
         mat(443) = rxt(93)*y(7)

         mat(30) = -( rxt(80) + rxt(83) + rxt(82)*y(51) + het_rates(5) )

         mat(27) = -( rxt(86) + rxt(87) + het_rates(6) )
         mat(455) = rxt(3)
         mat(29) = rxt(80) + rxt(82)*y(51)

         mat(148) = -( rxt(52) + het_rates(8) )
         mat(309) = rxt(6)
         mat(72) = rxt(229)

         mat(318) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(341) = rxt(8) + .500_r8*rxt(200)
         mat(45) = rxt(10)
         mat(388) = rxt(13)
         mat(122) = rxt(238)
         mat(445) = 2.000_r8*rxt(92)*y(7)

         mat(342) = -( rxt(8) + rxt(200) + het_rates(10) )
         mat(46) = rxt(9) + rxt(121)
         mat(114) = rxt(11)
         mat(389) = rxt(12)
         mat(66) = rxt(15) + rxt(130)
         mat(187) = rxt(30)
         mat(84) = rxt(35)

         mat(371) = -( rxt(131)*y(15) + rxt(138)*y(19) + rxt(139)*y(19) + rxt(150)*y(20) &
                      + rxt(207)*y(41) + rxt(208)*y(48) + rxt(209)*y(46) + rxt(210)*y(42) &
                 + het_rates(22) )
         mat(115) = rxt(11)
         mat(67) = rxt(14)
         mat(54) = rxt(16)
         mat(195) = rxt(19)
         mat(102) = 2.000_r8*rxt(22)
         mat(177) = rxt(27)
         mat(158) = rxt(33)
         mat(343) = .500_r8*rxt(200)
         mat(447) = rxt(103)*y(15) + rxt(106)*y(20)

         mat(391) = -( rxt(12) + rxt(13) + rxt(199) + het_rates(11) )
         mat(47) = rxt(9) + rxt(10) + rxt(121)
         mat(68) = rxt(14)
         mat(189) = rxt(29)
         mat(85) = rxt(34)

         mat(112) = -( rxt(11) + het_rates(12) )
         mat(44) = 2.000_r8*rxt(198) + 2.000_r8*rxt(211) + 2.000_r8*rxt(217) &
                      + 2.000_r8*rxt(222)
         mat(382) = rxt(199)
         mat(331) = .500_r8*rxt(200)
         mat(182) = rxt(212) + rxt(218) + rxt(223)
         mat(81) = rxt(213) + rxt(221) + rxt(224)

         mat(62) = -( rxt(14) + rxt(15) + rxt(130) + het_rates(13) )

         mat(43) = -( rxt(9) + rxt(10) + rxt(121) + rxt(198) + rxt(211) + rxt(217) &
                      + rxt(222) + het_rates(14) )

         mat(164) = -( het_rates(16) )
         mat(439) = rxt(103)*y(15)
         mat(359) = rxt(131)*y(15)
         mat(400) = rxt(162)*y(15)

         mat(49) = -( rxt(16) + het_rates(17) )

         mat(201) = -( rxt(17) + rxt(18) + het_rates(18) )
         mat(52) = rxt(16)
         mat(441) = rxt(104)*y(15) + rxt(105)*y(15)

         mat(213) = -( het_rates(21) )
         mat(53) = rxt(16)
         mat(202) = 2.000_r8*rxt(17)
         mat(194) = rxt(19) + 2.000_r8*rxt(21)
         mat(479) = rxt(28)
         mat(442) = rxt(104)*y(15) + rxt(106)*y(20)
         mat(364) = rxt(139)*y(19) + rxt(150)*y(20)
         mat(405) = rxt(157)*y(20)

         mat(295) = -( rxt(205) + het_rates(23) )
         mat(65) = rxt(15) + rxt(130)
         mat(444) = rxt(104)*y(15)
         mat(368) = rxt(138)*y(19) + rxt(207)*y(41) + rxt(210)*y(42)
         mat(408) = rxt(206)*y(41)

         mat(99) = -( rxt(22) + het_rates(24) )
         mat(284) = .500_r8*rxt(205)

         mat(193) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(78) )
         mat(362) = rxt(131)*y(15) + rxt(150)*y(20) + rxt(207)*y(41) + rxt(208)*y(48) &
                      + rxt(209)*y(46) + rxt(210)*y(42)

         mat(413) = -( rxt(157)*y(20) + rxt(162)*y(15) + rxt(206)*y(41) + het_rates(27) )
         mat(33) = 2.000_r8*rxt(23)
         mat(276) = rxt(24)
         mat(20) = 2.000_r8*rxt(26)
         mat(178) = rxt(27)
         mat(488) = rxt(28)
         mat(190) = rxt(29)
         mat(35) = rxt(31)
         mat(449) = 3.000_r8*rxt(95)*y(43) + 2.000_r8*rxt(96)*y(44) &
                      + 3.000_r8*rxt(97)*y(45) + rxt(98)*y(46) + 4.000_r8*rxt(99)*y(47)
         mat(373) = rxt(207)*y(41) + 3.000_r8*rxt(208)*y(48) + rxt(209)*y(46)

         mat(32) = -( rxt(23) + het_rates(28) )

         mat(270) = -( rxt(24) + het_rates(29) )
         mat(25) = rxt(25)
         mat(186) = rxt(30)
         mat(19) = 2.000_r8*rxt(173)

         mat(24) = -( rxt(25) + het_rates(30) )

         mat(18) = -( rxt(26) + rxt(173) + het_rates(31) )

         mat(492) = -( rxt(28) + het_rates(32) )
         mat(417) = rxt(157)*y(20) + rxt(162)*y(15) + 2.000_r8*rxt(206)*y(41)

         mat(174) = -( rxt(27) + het_rates(33) )
         mat(183) = rxt(212) + rxt(218) + rxt(223)

         mat(184) = -( rxt(29) + rxt(30) + rxt(212) + rxt(218) + rxt(223) + het_rates(34) &
       )

         mat(34) = -( rxt(31) + het_rates(35) )

         mat(432) = -( het_rates(36) )
         mat(36) = rxt(31)
         mat(255) = rxt(32)
         mat(160) = rxt(33)
         mat(86) = rxt(34)
         mat(450) = rxt(100)*y(42) + rxt(101)*y(50) + rxt(102)*y(49)
         mat(374) = rxt(210)*y(42)

         mat(247) = -( rxt(32) + het_rates(37) )
         mat(83) = rxt(35)

         mat(106) = -( het_rates(38) )

         mat(155) = -( rxt(33) + het_rates(39) )
         mat(82) = rxt(213) + rxt(221) + rxt(224)

         mat(80) = -( rxt(34) + rxt(35) + rxt(213) + rxt(221) + rxt(224) + het_rates(40) &
       )

         mat(90) = -( het_rates(52) )

         mat(117) = -( rxt(238) + het_rates(53) )
         mat(226) = rxt(53) + rxt(65)
         mat(70) = rxt(231)*y(51)

         mat(55) = -( het_rates(54) )
         mat(143) = rxt(52)

         mat(69) = -( rxt(229) + rxt(231)*y(51) + het_rates(55) )
         mat(497) = rxt(49) + rxt(50) + rxt(51) + rxt(62) + rxt(63) + rxt(64)
         mat(223) = rxt(55) + rxt(56) + rxt(57) + rxt(67) + rxt(68) + rxt(69)

         mat(126) = -( het_rates(56) )
         mat(307) = rxt(7)
         mat(71) = rxt(229)
         mat(118) = rxt(238)

         mat(75) = -( het_rates(58) )

         mat(137) = -( het_rates(57) )
         mat(308) = rxt(7)
         mat(506) = rxt(49) + rxt(50) + rxt(51) + rxt(62) + rxt(63) + rxt(64)
         mat(147) = rxt(52)
         mat(228) = rxt(53) + rxt(55) + rxt(56) + rxt(57) + rxt(65) + rxt(67) + rxt(68) &
                      + rxt(69)

         mat(22) = -( het_rates(59) )

         mat(38) = -( het_rates(60) )

         mat(1) = -( het_rates(61) )

         mat(2) = -( het_rates(62) )

         mat(3) = -( het_rates(63) )

         mat(4) = -( het_rates(64) )

         mat(5) = -( het_rates(65) )

         mat(6) = -( het_rates(66) )

         mat(7) = -( het_rates(67) )

         mat(8) = -( het_rates(68) )

         mat(9) = -( het_rates(69) )

         mat(10) = -( het_rates(70) )

         mat(11) = -( het_rates(71) )

         mat(12) = -( het_rates(72) )

         mat(13) = -( het_rates(73) )

         mat(14) = -( het_rates(74) )

         mat(15) = -( het_rates(75) )

         mat(16) = -( het_rates(76) )

         mat(17) = -( het_rates(77) )


      end subroutine linmat01

      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(in)    ::  het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) ::  mat(nzcnt)

      call linmat01( mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix
