













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

         mat(261) = -( rxt(20) + rxt(21) + rxt(22) + het_rates(1) )
         mat(187) = rxt(8)
         mat(175) = rxt(9)
         mat(41) = rxt(12)

         mat(245) = -( rxt(25) + rxt(26) + rxt(23)*y(15) + rxt(24)*y(15) + rxt(27)*y(16) &
                      + rxt(28)*y(17) + rxt(72)*y(20) + rxt(73)*y(20) + het_rates(2) )
         mat(33) = 2.000_r8*rxt(3)
         mat(73) = rxt(6)
         mat(78) = rxt(7)
         mat(63) = rxt(13)
         mat(68) = rxt(14)
         mat(260) = 2.000_r8*rxt(20) + rxt(21) + rxt(22)

         mat(203) = -( het_rates(3) )
         mat(80) = 2.000_r8*rxt(4)
         mat(71) = rxt(6)
         mat(77) = rxt(7)
         mat(67) = rxt(15) + rxt(54)
         mat(137) = rxt(16)
         mat(128) = .500_r8*rxt(19)
         mat(258) = rxt(21)
         mat(243) = rxt(25) + rxt(24)*y(15) + rxt(28)*y(17)

         mat(31) = -( rxt(3) + het_rates(4) )

         mat(79) = -( rxt(4) + rxt(5) + het_rates(5) )
         mat(70) = rxt(6)
         mat(120) = 2.000_r8*rxt(19)
         mat(232) = rxt(72)*y(20) + rxt(73)*y(20)

         mat(219) = -( het_rates(6) )
         mat(138) = rxt(16)
         mat(129) = .500_r8*rxt(19)
         mat(244) = rxt(26)

         mat(69) = -( rxt(6) + het_rates(7) )

         mat(156) = -( het_rates(8) )
         mat(181) = rxt(8)
         mat(169) = rxt(10)
         mat(38) = rxt(12)

         mat(183) = -( rxt(8) + het_rates(9) )
         mat(171) = rxt(9)
         mat(40) = rxt(11) + rxt(55)
         mat(62) = rxt(13)
         mat(66) = rxt(15) + rxt(54)
         mat(36) = rxt(17) + rxt(56)

         mat(170) = -( rxt(9) + rxt(10) + het_rates(10) )
         mat(39) = rxt(11) + rxt(12) + rxt(55)
         mat(65) = rxt(14)

         mat(37) = -( rxt(11) + rxt(12) + rxt(55) + het_rates(11) )

         mat(60) = -( rxt(13) + het_rates(12) )

         mat(64) = -( rxt(14) + rxt(15) + rxt(54) + het_rates(13) )

         mat(34) = -( rxt(17) + rxt(56) + het_rates(14) )

         mat(102) = -( het_rates(21) )
         mat(235) = rxt(27)*y(16)

         mat(74) = -( rxt(7) + het_rates(22) )

         mat(133) = -( rxt(16) + het_rates(23) )
         mat(76) = rxt(7)
         mat(238) = .200_r8*rxt(28)*y(17)

         mat(142) = -( het_rates(24) )
         mat(35) = rxt(17) + rxt(56)
         mat(239) = rxt(72)*y(20) + rxt(73)*y(20)

         mat(48) = -( het_rates(18) )

         mat(93) = -( het_rates(19) )

         mat(54) = -( het_rates(25) )

         mat(112) = -( het_rates(26) )

         mat(123) = -( rxt(19) + het_rates(27) )

         mat(84) = -( het_rates(28) )

         mat(43) = -( het_rates(31) )

         mat(29) = -( het_rates(32) )

         mat(1) = -( het_rates(33) )

         mat(2) = -( het_rates(34) )

         mat(3) = -( het_rates(35) )

         mat(4) = -( het_rates(36) )

         mat(5) = -( het_rates(37) )

         mat(6) = -( het_rates(38) )

         mat(7) = -( het_rates(39) )

         mat(8) = -( het_rates(40) )

         mat(9) = -( het_rates(41) )

         mat(10) = -( het_rates(42) )

         mat(11) = -( het_rates(43) )

         mat(12) = -( het_rates(44) )

         mat(13) = -( het_rates(45) )

         mat(14) = -( het_rates(46) )

         mat(15) = -( het_rates(47) )

         mat(16) = -( het_rates(48) )

         mat(17) = -( het_rates(49) )

         mat(18) = -( het_rates(50) )

         mat(19) = -( het_rates(51) )

         mat(20) = -( het_rates(52) )

         mat(21) = -( het_rates(53) )

         mat(22) = -( het_rates(54) )

         mat(23) = -( het_rates(55) )

         mat(24) = -( het_rates(56) )

         mat(25) = -( het_rates(57) )

         mat(26) = -( het_rates(58) )

         mat(27) = -( het_rates(59) )


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
