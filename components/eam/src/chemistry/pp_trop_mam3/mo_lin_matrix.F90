





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

         mat(7) = -( het_rates(5) )

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
