











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

         mat(212) = -( rxt(23) + rxt(24) + rxt(25)*y(34) + het_rates(1) )
         mat(186) = rxt(8)
         mat(197) = rxt(9)
         mat(42) = rxt(12)

         mat(235) = -( rxt(28) + rxt(26)*y(15) + rxt(27)*y(15) + rxt(29)*y(34) &
                      + rxt(30)*y(16) + rxt(31)*y(17) + rxt(75)*y(20) + rxt(76)*y(20) &
                      + rxt(88)*y(36) + rxt(89)*y(36) + rxt(90)*y(37) + het_rates(2) )
         mat(37) = 2.000_r8*rxt(3)
         mat(58) = rxt(6)
         mat(73) = rxt(7)
         mat(63) = rxt(13)
         mat(68) = rxt(14)
         mat(213) = 2.000_r8*rxt(23) + rxt(24) + rxt(25)*y(34)

         mat(155) = -( het_rates(3) )
         mat(75) = 2.000_r8*rxt(4)
         mat(57) = rxt(6)
         mat(72) = rxt(7)
         mat(65) = rxt(15) + rxt(57)
         mat(130) = rxt(16)
         mat(121) = .500_r8*rxt(19)
         mat(208) = rxt(24)
         mat(230) = rxt(28) + rxt(27)*y(15) + rxt(31)*y(17) + .500_r8*rxt(89)*y(36)

         mat(35) = -( rxt(3) + het_rates(4) )

         mat(74) = -( rxt(4) + rxt(5) + het_rates(5) )
         mat(56) = rxt(6)
         mat(115) = 2.000_r8*rxt(19)
         mat(222) = rxt(75)*y(20) + rxt(76)*y(20)

         mat(251) = -( het_rates(6) )
         mat(134) = rxt(16)
         mat(126) = .500_r8*rxt(19)
         mat(236) = rxt(29)*y(34)

         mat(55) = -( rxt(6) + het_rates(7) )

         mat(171) = -( het_rates(8) )
         mat(183) = rxt(8)
         mat(194) = rxt(10)
         mat(39) = rxt(12)

         mat(184) = -( rxt(8) + het_rates(9) )
         mat(195) = rxt(9)
         mat(40) = rxt(11) + rxt(58)
         mat(61) = rxt(13)
         mat(66) = rxt(15) + rxt(57)
         mat(34) = rxt(17) + rxt(59)

         mat(196) = -( rxt(9) + rxt(10) + rxt(91)*y(36) + het_rates(10) )
         mat(41) = rxt(11) + rxt(12) + rxt(58)
         mat(67) = rxt(14)

         mat(38) = -( rxt(11) + rxt(12) + rxt(58) + het_rates(11) )

         mat(60) = -( rxt(13) + het_rates(12) )
         mat(190) = rxt(91)*y(36)

         mat(64) = -( rxt(14) + rxt(15) + rxt(57) + het_rates(13) )

         mat(32) = -( rxt(17) + rxt(59) + het_rates(14) )

         mat(97) = -( het_rates(21) )
         mat(225) = rxt(30)*y(16)

         mat(69) = -( rxt(7) + het_rates(22) )

         mat(128) = -( rxt(16) + het_rates(23) )
         mat(71) = rxt(7)
         mat(228) = .200_r8*rxt(31)*y(17)

         mat(137) = -( het_rates(24) )
         mat(33) = rxt(17) + rxt(59)
         mat(229) = rxt(75)*y(20) + rxt(76)*y(20)

         mat(43) = -( het_rates(18) )

         mat(88) = -( het_rates(19) )

         mat(49) = -( het_rates(25) )

         mat(1) = -( rxt(95) + rxt(97) + rxt(99) + het_rates(26) )

         mat(107) = -( het_rates(27) )

         mat(5) = -( rxt(96) + rxt(98) + rxt(100) + het_rates(28) )

         mat(118) = -( rxt(19) + het_rates(29) )

         mat(79) = -( het_rates(30) )

         mat(10) = -( rxt(20) + het_rates(54) )

         mat(11) = -( rxt(21) + het_rates(55) )

         mat(12) = -( rxt(22) + het_rates(56) )

         mat(13) = -( rxt(92) + het_rates(39) )

         mat(15) = -( rxt(93) + het_rates(40) )
         mat(14) = 1.150_r8*rxt(92)

         mat(17) = -( rxt(94) + het_rates(41) )
         mat(16) = 1.150_r8*rxt(93)

         mat(20) = -( het_rates(42) )
         mat(6) = .0436005_r8*rxt(98)
         mat(21) = .460_r8*rxt(104)

         mat(22) = -( rxt(104) + het_rates(43) )
         mat(2) = .0027205_r8*rxt(95) + .0002705_r8*rxt(99)
         mat(7) = .1498605_r8*rxt(96) + .0103505_r8*rxt(98)
         mat(24) = .460_r8*rxt(103)

         mat(25) = -( rxt(103) + het_rates(44) )
         mat(18) = .460_r8*rxt(94)
         mat(3) = .0098105_r8*rxt(95) + .0032705_r8*rxt(97) + .0406005_r8*rxt(99)
         mat(8) = .0545005_r8*rxt(96) + .0980905_r8*rxt(98) + .1749305_r8*rxt(100)
         mat(27) = .460_r8*rxt(102)

         mat(28) = -( rxt(102) + het_rates(45) )
         mat(4) = .0021805_r8*rxt(95) + .0645805_r8*rxt(99)
         mat(9) = .0926405_r8*rxt(96) + .0163505_r8*rxt(98) + .5901905_r8*rxt(100)
         mat(30) = .460_r8*rxt(101)

         mat(31) = -( rxt(101) + het_rates(46) )
         mat(19) = .500_r8*rxt(94)
         mat(31) = mat(31) + .500_r8*rxt(101)
         mat(29) = .500_r8*rxt(102)
         mat(26) = .500_r8*rxt(103)
         mat(23) = .500_r8*rxt(104)


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
