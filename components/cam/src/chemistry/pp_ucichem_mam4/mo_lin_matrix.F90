













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

         mat(256) = -( rxt(17) + rxt(18) + rxt(19)*y(16) + het_rates(1) )
         mat(186) = rxt(8)
         mat(197) = .886_r8*rxt(9)

         mat(240) = -( rxt(22) + rxt(20)*y(15) + rxt(21)*y(15) + rxt(23)*y(16) &
                      + rxt(24)*y(17) + rxt(25)*y(18) + rxt(68)*y(21) + rxt(69)*y(21) &
                 + het_rates(2) )
         mat(37) = 2.000_r8*rxt(3)
         mat(70) = rxt(6)
         mat(75) = rxt(7)
         mat(60) = rxt(11)
         mat(65) = .330_r8*rxt(12)
         mat(255) = 2.000_r8*rxt(17) + rxt(18) + rxt(19)*y(16)

         mat(214) = -( rxt(82) + het_rates(3) )
         mat(77) = 2.000_r8*rxt(4)
         mat(69) = rxt(6)
         mat(74) = rxt(7)
         mat(64) = .670_r8*rxt(12) + rxt(50)
         mat(135) = rxt(13)
         mat(126) = .500_r8*rxt(16)
         mat(254) = rxt(18)
         mat(239) = rxt(22) + rxt(21)*y(15) + rxt(25)*y(18)

         mat(35) = -( rxt(3) + het_rates(4) )

         mat(76) = -( rxt(4) + rxt(5) + het_rates(5) )
         mat(67) = rxt(6)
         mat(117) = 2.000_r8*rxt(16)
         mat(227) = rxt(68)*y(21) + rxt(69)*y(21)

         mat(169) = -( het_rates(6) )
         mat(132) = rxt(13)
         mat(124) = .500_r8*rxt(16)
         mat(236) = rxt(23)*y(16)

         mat(66) = -( rxt(6) + het_rates(7) )

         mat(153) = -( het_rates(8) )
         mat(180) = rxt(8)
         mat(191) = .114_r8*rxt(9)

         mat(182) = -( rxt(8) + het_rates(9) )
         mat(193) = .886_r8*rxt(9)
         mat(43) = rxt(10) + rxt(51)
         mat(58) = rxt(11)
         mat(62) = .670_r8*rxt(12) + rxt(50)
         mat(40) = rxt(14) + rxt(52)

         mat(194) = -( rxt(9) + rxt(81) + het_rates(10) )
         mat(44) = rxt(10) + rxt(51)
         mat(63) = .330_r8*rxt(12)

         mat(41) = -( rxt(10) + rxt(51) + rxt(80) + het_rates(11) )

         mat(57) = -( rxt(11) + het_rates(12) )
         mat(42) = 2.000_r8*rxt(80)
         mat(188) = rxt(81)

         mat(61) = -( rxt(12) + rxt(50) + het_rates(13) )

         mat(38) = -( rxt(14) + rxt(52) + het_rates(14) )

         mat(99) = -( het_rates(22) )
         mat(230) = rxt(24)*y(17)

         mat(71) = -( rxt(7) + het_rates(23) )

         mat(130) = -( rxt(13) + het_rates(24) )
         mat(73) = rxt(7)
         mat(233) = .200_r8*rxt(25)*y(18)

         mat(139) = -( het_rates(25) )
         mat(39) = rxt(14) + rxt(52)
         mat(234) = rxt(68)*y(21) + rxt(69)*y(21)

         mat(45) = -( het_rates(19) )

         mat(90) = -( het_rates(20) )

         mat(51) = -( het_rates(26) )

         mat(109) = -( het_rates(27) )

         mat(120) = -( rxt(16) + het_rates(28) )

         mat(81) = -( het_rates(29) )

         mat(32) = -( het_rates(32) )

         mat(29) = -( het_rates(33) )

         mat(1) = -( het_rates(34) )

         mat(2) = -( het_rates(35) )

         mat(3) = -( het_rates(36) )

         mat(4) = -( het_rates(37) )

         mat(5) = -( het_rates(38) )

         mat(6) = -( het_rates(39) )

         mat(7) = -( het_rates(40) )

         mat(8) = -( het_rates(41) )

         mat(9) = -( het_rates(42) )

         mat(10) = -( het_rates(43) )

         mat(11) = -( het_rates(44) )

         mat(12) = -( het_rates(45) )

         mat(13) = -( het_rates(46) )

         mat(14) = -( het_rates(47) )

         mat(15) = -( het_rates(48) )

         mat(16) = -( het_rates(49) )

         mat(17) = -( het_rates(50) )

         mat(18) = -( het_rates(51) )

         mat(19) = -( het_rates(52) )

         mat(20) = -( het_rates(53) )

         mat(21) = -( het_rates(54) )

         mat(22) = -( het_rates(55) )

         mat(23) = -( het_rates(56) )

         mat(24) = -( het_rates(57) )

         mat(25) = -( het_rates(58) )

         mat(26) = -( het_rates(59) )

         mat(27) = -( het_rates(60) )


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
