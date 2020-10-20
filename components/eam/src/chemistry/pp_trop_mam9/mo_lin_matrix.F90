






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

         mat(1) = -( rxt(1) + rxt(3) + het_rates(1) )

         mat(2) = -( het_rates(2) )
         mat(3) = rxt(4)

         mat(4) = -( rxt(4) + het_rates(3) )
         mat(5) = rxt(5) + .500_r8*rxt(6) + rxt(7)

         mat(6) = -( rxt(5) + rxt(6) + rxt(7) + het_rates(4) )

         mat(7) = -( rxt(8) + het_rates(5) )

         mat(8) = -( het_rates(6) )

         mat(9) = -( het_rates(7) )

         mat(10) = -( het_rates(8) )

         mat(11) = -( het_rates(9) )

         mat(12) = -( het_rates(10) )

         mat(13) = -( het_rates(11) )

         mat(14) = -( het_rates(12) )

         mat(15) = -( het_rates(13) )

         mat(16) = -( het_rates(14) )

         mat(17) = -( het_rates(15) )

         mat(18) = -( het_rates(16) )

         mat(19) = -( het_rates(17) )

         mat(20) = -( het_rates(18) )

         mat(21) = -( het_rates(19) )

         mat(22) = -( het_rates(20) )

         mat(23) = -( het_rates(21) )

         mat(24) = -( het_rates(22) )

         mat(25) = -( het_rates(23) )

         mat(26) = -( het_rates(24) )

         mat(27) = -( het_rates(25) )

         mat(28) = -( het_rates(26) )

         mat(29) = -( het_rates(27) )

         mat(30) = -( het_rates(28) )

         mat(31) = -( het_rates(29) )

         mat(32) = -( het_rates(30) )

         mat(33) = -( het_rates(31) )

         mat(34) = -( het_rates(32) )

         mat(35) = -( het_rates(33) )

         mat(36) = -( het_rates(34) )

         mat(37) = -( het_rates(35) )

         mat(38) = -( het_rates(36) )

         mat(39) = -( het_rates(37) )

         mat(40) = -( het_rates(38) )

         mat(41) = -( het_rates(39) )

         mat(42) = -( het_rates(40) )

         mat(43) = -( het_rates(41) )

         mat(44) = -( het_rates(42) )

         mat(45) = -( het_rates(43) )

         mat(46) = -( het_rates(44) )

         mat(47) = -( het_rates(45) )

         mat(48) = -( het_rates(46) )

         mat(49) = -( het_rates(47) )

         mat(50) = -( het_rates(48) )

         mat(51) = -( het_rates(49) )

         mat(52) = -( het_rates(50) )

         mat(53) = -( het_rates(51) )

         mat(54) = -( het_rates(52) )

         mat(55) = -( het_rates(53) )

         mat(56) = -( het_rates(54) )


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
