





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

         mat(41) = -( rxt(1) + het_rates(1) )
         mat(33) = rxt(3)

         mat(74) = -( rxt(15) + rxt(16)*y(8) + het_rates(2) )
         mat(44) = 2.000_r8*rxt(1)
         mat(10) = 2.000_r8*rxt(2)
         mat(29) = rxt(6)

         mat(51) = -( het_rates(3) )
         mat(24) = 2.000_r8*rxt(4)
         mat(28) = rxt(6)
         mat(72) = rxt(16)*y(8)

         mat(7) = -( rxt(2) + het_rates(4) )

         mat(59) = -( het_rates(5) )
         mat(34) = rxt(3)

         mat(32) = -( rxt(3) + rxt(23) + het_rates(6) )

         mat(1) = -( het_rates(7) )
         mat(31) = .500_r8*rxt(23)

         mat(23) = -( rxt(4) + rxt(5) + het_rates(9) )
         mat(26) = rxt(6)

         mat(83) = -( het_rates(10) )
         mat(75) = rxt(15)

         mat(27) = -( rxt(6) + het_rates(11) )

         mat(3) = -( het_rates(12) )

         mat(13) = -( het_rates(13) )

         mat(2) = -( het_rates(14) )

         mat(17) = -( het_rates(15) )


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
