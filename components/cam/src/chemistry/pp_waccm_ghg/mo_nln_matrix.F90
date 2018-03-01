      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: nlnmat
      contains
      subroutine nlnmat( mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      subroutine nlnmat_finit( mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = lmat( 3)
         mat( 4) = lmat( 4)
         mat( 5) = lmat( 5)
         mat( 6) = lmat( 6)
         mat( 1) = mat( 1) - dti
         mat( 3) = mat( 3) - dti
         mat( 4) = mat( 4) - dti
         mat( 5) = mat( 5) - dti
         mat( 6) = mat( 6) - dti
      end subroutine nlnmat_finit
      end module mo_nln_matrix
