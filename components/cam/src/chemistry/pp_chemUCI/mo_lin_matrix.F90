













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

         mat(132) = -( rxt(17) + rxt(18) + rxt(19)*y(16) + het_rates(1) )
         mat(117) = rxt(8)
         mat(158) = .886_r8*rxt(9)

         mat(183) = -( rxt(22) + rxt(20)*y(15) + rxt(21)*y(15) + rxt(23)*y(16) &
                      + rxt(24)*y(17) + rxt(25)*y(18) + rxt(68)*y(21) + rxt(69)*y(21) &
                 + het_rates(2) )
         mat(2) = 2.000_r8*rxt(3)
         mat(26) = rxt(6)
         mat(40) = rxt(7)
         mat(31) = rxt(11)
         mat(35) = .330_r8*rxt(12)
         mat(135) = 2.000_r8*rxt(17) + rxt(18) + rxt(19)*y(16)

         mat(219) = -( rxt(82) + het_rates(3) )
         mat(44) = 2.000_r8*rxt(4)
         mat(27) = rxt(6)
         mat(41) = rxt(7)
         mat(36) = .670_r8*rxt(12) + rxt(50)
         mat(102) = rxt(13)
         mat(94) = .500_r8*rxt(16)
         mat(137) = rxt(18)
         mat(185) = rxt(22) + rxt(21)*y(15) + rxt(25)*y(18)

         mat(1) = -( rxt(3) + het_rates(4) )

         mat(42) = -( rxt(4) + rxt(5) + het_rates(5) )
         mat(24) = rxt(6)
         mat(83) = 2.000_r8*rxt(16)
         mat(171) = rxt(68)*y(21) + rxt(69)*y(21)

         mat(148) = -( het_rates(6) )
         mat(99) = rxt(13)
         mat(91) = .500_r8*rxt(16)
         mat(181) = rxt(23)*y(16)

         mat(23) = -( rxt(6) + het_rates(7) )

         mat(199) = -( het_rates(8) )
         mat(121) = rxt(8)
         mat(162) = .114_r8*rxt(9)

         mat(116) = -( rxt(8) + het_rates(9) )
         mat(157) = .886_r8*rxt(9)
         mat(9) = rxt(10) + rxt(51)
         mat(29) = rxt(11)
         mat(33) = .670_r8*rxt(12) + rxt(50)
         mat(6) = rxt(14) + rxt(52)

         mat(160) = -( rxt(9) + rxt(81) + het_rates(10) )
         mat(10) = rxt(10) + rxt(51)
         mat(34) = .330_r8*rxt(12)

         mat(7) = -( rxt(10) + rxt(51) + rxt(80) + het_rates(11) )

         mat(28) = -( rxt(11) + het_rates(12) )
         mat(8) = 2.000_r8*rxt(80)
         mat(154) = rxt(81)

         mat(32) = -( rxt(12) + rxt(50) + het_rates(13) )

         mat(4) = -( rxt(14) + rxt(52) + het_rates(14) )

         mat(65) = -( het_rates(22) )
         mat(174) = rxt(24)*y(17)

         mat(37) = -( rxt(7) + het_rates(23) )

         mat(96) = -( rxt(13) + het_rates(24) )
         mat(39) = rxt(7)
         mat(177) = .200_r8*rxt(25)*y(18)

         mat(105) = -( het_rates(25) )
         mat(5) = rxt(14) + rxt(52)
         mat(178) = rxt(68)*y(21) + rxt(69)*y(21)

         mat(11) = -( het_rates(19) )

         mat(56) = -( het_rates(20) )

         mat(17) = -( het_rates(26) )

         mat(75) = -( het_rates(27) )

         mat(86) = -( rxt(16) + het_rates(28) )

         mat(47) = -( het_rates(29) )


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
