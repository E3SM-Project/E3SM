













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

         mat(171) = -( rxt(17) + rxt(18) + rxt(19) + het_rates(1) )
         mat(242) = rxt(7)
         mat(156) = .886_r8*rxt(8)

         mat(233) = -( rxt(46) + rxt(47)*y(16) + rxt(48)*y(15) + rxt(49)*y(15) &
                      + rxt(50)*y(17) + rxt(51)*y(18) + rxt(72)*y(20) + rxt(73)*y(20) &
                 + het_rates(2) )
         mat(27) = 2.000_r8*rxt(2)
         mat(51) = rxt(5)
         mat(65) = rxt(6)
         mat(54) = rxt(10)
         mat(59) = .330_r8*rxt(11)
         mat(174) = 2.000_r8*rxt(17) + rxt(18) + 2.000_r8*rxt(19)

         mat(191) = -( rxt(42)*y(34) + het_rates(3) )
         mat(67) = 2.000_r8*rxt(3)
         mat(49) = rxt(5)
         mat(64) = rxt(6)
         mat(58) = .670_r8*rxt(11) + rxt(37)
         mat(123) = rxt(12)
         mat(115) = .500_r8*rxt(16)
         mat(172) = rxt(18)
         mat(231) = rxt(46) + rxt(48)*y(15) + rxt(49)*y(15) + rxt(51)*y(18)

         mat(25) = -( rxt(2) + het_rates(4) )

         mat(66) = -( rxt(3) + rxt(4) + het_rates(5) )
         mat(48) = rxt(5)
         mat(107) = 2.000_r8*rxt(16)
         mat(220) = rxt(72)*y(20) + rxt(73)*y(20)

         mat(207) = -( het_rates(6) )
         mat(124) = rxt(12)
         mat(116) = .500_r8*rxt(16)
         mat(232) = rxt(47)*y(16)

         mat(47) = -( rxt(5) + het_rates(7) )

         mat(143) = -( het_rates(8) )
         mat(240) = rxt(7)
         mat(154) = .114_r8*rxt(8)

         mat(246) = -( rxt(7) + het_rates(9) )
         mat(160) = .886_r8*rxt(8)
         mat(34) = rxt(9) + rxt(38)
         mat(55) = rxt(10)
         mat(60) = .670_r8*rxt(11) + rxt(37)
         mat(30) = rxt(13) + rxt(39)

         mat(155) = -( rxt(8) + rxt(41)*y(34) + het_rates(10) )
         mat(33) = rxt(9) + rxt(38)
         mat(57) = .330_r8*rxt(11)

         mat(31) = -( rxt(9) + rxt(38) + rxt(40)*y(34) + het_rates(11) )

         mat(52) = -( rxt(10) + het_rates(12) )
         mat(32) = 2.000_r8*rxt(40)*y(34)
         mat(151) = rxt(41)*y(34)

         mat(56) = -( rxt(11) + rxt(37) + het_rates(13) )

         mat(28) = -( rxt(13) + rxt(39) + het_rates(14) )

         mat(89) = -( het_rates(21) )
         mat(223) = rxt(50)*y(17)

         mat(61) = -( rxt(6) + het_rates(22) )

         mat(120) = -( rxt(12) + het_rates(23) )
         mat(63) = rxt(6)
         mat(226) = .200_r8*rxt(51)*y(18)

         mat(129) = -( het_rates(24) )
         mat(29) = rxt(13) + rxt(39)
         mat(227) = rxt(72)*y(20) + rxt(73)*y(20)

         mat(35) = -( het_rates(19) )

         mat(80) = -( het_rates(29) )

         mat(41) = -( het_rates(25) )

         mat(99) = -( het_rates(26) )

         mat(110) = -( rxt(16) + het_rates(27) )

         mat(71) = -( het_rates(28) )

         mat(22) = -( het_rates(36) )

         mat(19) = -( het_rates(35) )

         mat(1) = -( het_rates(37) )

         mat(2) = -( het_rates(38) )

         mat(3) = -( het_rates(39) )

         mat(4) = -( het_rates(40) )

         mat(5) = -( het_rates(41) )

         mat(6) = -( het_rates(42) )

         mat(7) = -( het_rates(43) )

         mat(8) = -( het_rates(44) )

         mat(9) = -( het_rates(45) )

         mat(10) = -( het_rates(46) )

         mat(11) = -( het_rates(47) )

         mat(12) = -( het_rates(48) )

         mat(13) = -( het_rates(49) )

         mat(14) = -( het_rates(50) )

         mat(15) = -( het_rates(51) )

         mat(16) = -( het_rates(52) )

         mat(17) = -( het_rates(53) )


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
