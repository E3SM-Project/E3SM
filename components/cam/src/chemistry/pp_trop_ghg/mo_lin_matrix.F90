





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

         mat(1) = -( rxt(1) + het_rates(1) )

         mat(3) = -( rxt(2) + het_rates(2) )

         mat(4) = -( rxt(3) + het_rates(3) )

         mat(5) = -( rxt(4) + het_rates(4) )

         mat(6) = -( rxt(5) + het_rates(5) )
         mat(2) = 2.000_r8*rxt(1)


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
