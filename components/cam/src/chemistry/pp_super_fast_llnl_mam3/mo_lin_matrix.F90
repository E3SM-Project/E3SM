




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

         mat(92) = -( rxt(1) + het_rates(1) )
         mat(47) = rxt(3)

         mat(57) = -( rxt(15) + rxt(16)*y(8) + het_rates(2) )
         mat(88) = 2.000_r8*rxt(1)
         mat(27) = 2.000_r8*rxt(2)
         mat(40) = rxt(6)

         mat(66) = -( het_rates(3) )
         mat(37) = 2.000_r8*rxt(4)
         mat(41) = rxt(6)
         mat(58) = rxt(16)*y(8)

         mat(26) = -( rxt(2) + het_rates(4) )

         mat(83) = -( het_rates(5) )
         mat(46) = rxt(3)

         mat(44) = -( rxt(3) + rxt(23) + het_rates(6) )

         mat(1) = -( het_rates(7) )
         mat(43) = .500_r8*rxt(23)

         mat(35) = -( rxt(4) + rxt(5) + het_rates(9) )
         mat(38) = rxt(6)

         mat(75) = -( het_rates(10) )
         mat(59) = rxt(15)

         mat(39) = -( rxt(6) + het_rates(11) )

         mat(23) = -( het_rates(12) )

         mat(20) = -( het_rates(13) )

         mat(29) = -( het_rates(14) )

         mat(2) = -( het_rates(15) )

         mat(3) = -( het_rates(16) )

         mat(4) = -( het_rates(17) )

         mat(5) = -( het_rates(18) )

         mat(6) = -( het_rates(19) )

         mat(7) = -( het_rates(20) )

         mat(8) = -( het_rates(21) )

         mat(9) = -( het_rates(22) )

         mat(10) = -( het_rates(23) )

         mat(11) = -( het_rates(24) )

         mat(12) = -( het_rates(25) )

         mat(13) = -( het_rates(26) )

         mat(14) = -( het_rates(27) )

         mat(15) = -( het_rates(28) )

         mat(16) = -( het_rates(29) )

         mat(17) = -( het_rates(30) )

         mat(18) = -( het_rates(31) )


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
