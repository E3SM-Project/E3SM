





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

         mat(496) = -( rxt(3) + rxt(4) + het_rates(1) )

         mat(279) = -( rxt(49) + rxt(50) + rxt(51) + rxt(62) + rxt(63) + rxt(64) &
                 + het_rates(2) )
         mat(206) = rxt(1) + 2.000_r8*rxt(2) + rxt(55) + rxt(56) + rxt(57) &
                      + 2.000_r8*rxt(60) + rxt(67) + rxt(68) + rxt(69) + 2.000_r8*rxt(72)
         mat(486) = rxt(4)
         mat(301) = rxt(6)
         mat(342) = rxt(8)
         mat(19) = rxt(10)
         mat(426) = rxt(12)
         mat(169) = rxt(21)
         mat(245) = rxt(24)
         mat(11) = rxt(25)
         mat(223) = rxt(32)
         mat(360) = rxt(88)

         mat(364) = -( rxt(88) + rxt(92)*y(7) + rxt(93)*y(7) + rxt(95)*y(43) &
                      + rxt(96)*y(44) + rxt(97)*y(45) + rxt(98)*y(46) + rxt(99)*y(47) &
                      + rxt(100)*y(42) + rxt(101)*y(50) + rxt(102)*y(49) + rxt(103)*y(15) &
                      + rxt(104)*y(15) + rxt(105)*y(15) + rxt(106)*y(20) + het_rates(3) )
         mat(209) = rxt(1)
         mat(490) = rxt(3)
         mat(170) = rxt(20)

         mat(205) = -( rxt(1) + rxt(2) + rxt(53) + rxt(55) + rxt(56) + rxt(57) + rxt(60) &
                      + rxt(65) + rxt(67) + rxt(68) + rxt(69) + rxt(72) + het_rates(4) )
         mat(483) = rxt(4)
         mat(425) = rxt(13)
         mat(8) = rxt(83)
         mat(5) = rxt(86) + rxt(87)
         mat(359) = rxt(93)*y(7)

         mat(7) = -( rxt(80) + rxt(83) + rxt(82)*y(51) + het_rates(5) )

         mat(4) = -( rxt(86) + rxt(87) + het_rates(6) )
         mat(480) = rxt(3)
         mat(6) = rxt(80) + rxt(82)*y(51)

         mat(131) = -( rxt(52) + het_rates(8) )
         mat(293) = rxt(6)
         mat(46) = rxt(221)

         mat(302) = -( rxt(6) + rxt(7) + het_rates(9) )
         mat(343) = rxt(8)
         mat(20) = rxt(10)
         mat(427) = rxt(13)
         mat(106) = rxt(230)
         mat(361) = 2.000_r8*rxt(92)*y(7)

         mat(345) = -( rxt(8) + het_rates(10) )
         mat(21) = rxt(9) + rxt(121)
         mat(88) = rxt(11)
         mat(429) = rxt(12)
         mat(32) = rxt(15) + rxt(130)
         mat(162) = rxt(30)
         mat(60) = rxt(35)

         mat(415) = -( rxt(131)*y(15) + rxt(138)*y(19) + rxt(139)*y(19) + rxt(150)*y(20) &
                      + rxt(199)*y(41) + rxt(200)*y(48) + rxt(201)*y(46) + rxt(202)*y(42) &
                 + het_rates(22) )
         mat(89) = rxt(11)
         mat(34) = rxt(14)
         mat(28) = rxt(16)
         mat(171) = rxt(19)
         mat(66) = 2.000_r8*rxt(22)
         mat(152) = rxt(27)
         mat(97) = rxt(33)
         mat(366) = rxt(103)*y(15) + rxt(106)*y(20)

         mat(433) = -( rxt(12) + rxt(13) + het_rates(11) )
         mat(22) = rxt(9) + rxt(10) + rxt(121)
         mat(35) = rxt(14)
         mat(164) = rxt(29)
         mat(61) = rxt(34)

         mat(86) = -( rxt(11) + het_rates(12) )
         mat(18) = 2.000_r8*rxt(203) + 2.000_r8*rxt(209) + 2.000_r8*rxt(214)
         mat(156) = rxt(204) + rxt(210) + rxt(215)
         mat(55) = rxt(205) + rxt(213) + rxt(216)

         mat(29) = -( rxt(14) + rxt(15) + rxt(130) + het_rates(13) )

         mat(17) = -( rxt(9) + rxt(10) + rxt(121) + rxt(203) + rxt(209) + rxt(214) &
                 + het_rates(14) )

         mat(138) = -( het_rates(16) )
         mat(355) = rxt(103)*y(15)
         mat(400) = rxt(131)*y(15)
         mat(439) = rxt(162)*y(15)

         mat(23) = -( rxt(16) + het_rates(17) )

         mat(175) = -( rxt(17) + rxt(18) + het_rates(18) )
         mat(26) = rxt(16)
         mat(357) = rxt(104)*y(15) + rxt(105)*y(15)

         mat(187) = -( het_rates(21) )
         mat(27) = rxt(16)
         mat(176) = 2.000_r8*rxt(17)
         mat(168) = rxt(19) + 2.000_r8*rxt(21)
         mat(465) = rxt(28)
         mat(358) = rxt(104)*y(15) + rxt(106)*y(20)
         mat(405) = rxt(139)*y(19) + rxt(150)*y(20)
         mat(444) = rxt(157)*y(20)

         mat(389) = -( het_rates(23) )
         mat(33) = rxt(15) + rxt(130)
         mat(365) = rxt(104)*y(15)
         mat(414) = rxt(138)*y(19) + rxt(199)*y(41) + rxt(202)*y(42)
         mat(452) = rxt(198)*y(41)

         mat(62) = -( rxt(22) + het_rates(24) )

         mat(167) = -( rxt(19) + rxt(20) + rxt(21) + het_rates(59) )
         mat(403) = rxt(131)*y(15) + rxt(150)*y(20) + rxt(199)*y(41) + rxt(200)*y(48) &
                      + rxt(201)*y(46) + rxt(202)*y(42)

         mat(455) = -( rxt(157)*y(20) + rxt(162)*y(15) + rxt(198)*y(41) + het_rates(27) )
         mat(13) = 2.000_r8*rxt(23)
         mat(253) = rxt(24)
         mat(3) = 2.000_r8*rxt(26)
         mat(153) = rxt(27)
         mat(477) = rxt(28)
         mat(165) = rxt(29)
         mat(16) = rxt(31)
         mat(368) = 3.000_r8*rxt(95)*y(43) + 2.000_r8*rxt(96)*y(44) &
                      + 3.000_r8*rxt(97)*y(45) + rxt(98)*y(46) + 4.000_r8*rxt(99)*y(47)
         mat(417) = rxt(199)*y(41) + 3.000_r8*rxt(200)*y(48) + rxt(201)*y(46)

         mat(12) = -( rxt(23) + het_rates(28) )

         mat(244) = -( rxt(24) + het_rates(29) )
         mat(10) = rxt(25)
         mat(160) = rxt(30)
         mat(2) = 2.000_r8*rxt(173)

         mat(9) = -( rxt(25) + het_rates(30) )

         mat(1) = -( rxt(26) + rxt(173) + het_rates(31) )

         mat(478) = -( rxt(28) + het_rates(32) )
         mat(456) = rxt(157)*y(20) + rxt(162)*y(15) + 2.000_r8*rxt(198)*y(41)

         mat(148) = -( rxt(27) + het_rates(33) )
         mat(157) = rxt(204) + rxt(210) + rxt(215)

         mat(158) = -( rxt(29) + rxt(30) + rxt(204) + rxt(210) + rxt(215) + het_rates(34) &
       )

         mat(14) = -( rxt(31) + het_rates(35) )

         mat(321) = -( het_rates(36) )
         mat(15) = rxt(31)
         mat(225) = rxt(32)
         mat(96) = rxt(33)
         mat(59) = rxt(34)
         mat(362) = rxt(100)*y(42) + rxt(101)*y(50) + rxt(102)*y(49)
         mat(411) = rxt(202)*y(42)

         mat(221) = -( rxt(32) + het_rates(37) )
         mat(57) = rxt(35)

         mat(80) = -( het_rates(38) )

         mat(92) = -( rxt(33) + het_rates(39) )
         mat(56) = rxt(205) + rxt(213) + rxt(216)

         mat(54) = -( rxt(34) + rxt(35) + rxt(205) + rxt(213) + rxt(216) + het_rates(40) &
       )

         mat(71) = -( het_rates(52) )

         mat(100) = -( rxt(230) + het_rates(53) )
         mat(200) = rxt(53) + rxt(65)
         mat(44) = rxt(223)*y(51)

         mat(36) = -( het_rates(54) )
         mat(126) = rxt(52)

         mat(43) = -( rxt(221) + rxt(223)*y(51) + het_rates(55) )
         mat(259) = rxt(49) + rxt(50) + rxt(51) + rxt(62) + rxt(63) + rxt(64)
         mat(197) = rxt(55) + rxt(56) + rxt(57) + rxt(67) + rxt(68) + rxt(69)

         mat(109) = -( het_rates(56) )
         mat(291) = rxt(7)
         mat(45) = rxt(221)
         mat(101) = rxt(230)

         mat(49) = -( het_rates(58) )

         mat(120) = -( het_rates(57) )
         mat(292) = rxt(7)
         mat(269) = rxt(49) + rxt(50) + rxt(51) + rxt(62) + rxt(63) + rxt(64)
         mat(130) = rxt(52)
         mat(202) = rxt(53) + rxt(55) + rxt(56) + rxt(57) + rxt(65) + rxt(67) + rxt(68) &
                      + rxt(69)


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
