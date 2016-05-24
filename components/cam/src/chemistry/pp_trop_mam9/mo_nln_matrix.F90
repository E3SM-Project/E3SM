






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
         mat( 7) = lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = lmat( 12)
         mat( 13) = lmat( 13)
         mat( 14) = lmat( 14)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = lmat( 17)
         mat( 18) = lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = lmat( 24)
         mat( 25) = lmat( 25)
         mat( 26) = lmat( 26)
         mat( 27) = lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = lmat( 29)
         mat( 30) = lmat( 30)
         mat( 31) = lmat( 31)
         mat( 32) = lmat( 32)
         mat( 33) = lmat( 33)
         mat( 34) = lmat( 34)
         mat( 35) = lmat( 35)
         mat( 36) = lmat( 36)
         mat( 37) = lmat( 37)
         mat( 38) = lmat( 38)
         mat( 39) = lmat( 39)
         mat( 40) = lmat( 40)
         mat( 41) = lmat( 41)
         mat( 42) = lmat( 42)
         mat( 43) = lmat( 43)
         mat( 44) = lmat( 44)
         mat( 45) = lmat( 45)
         mat( 46) = lmat( 46)
         mat( 47) = lmat( 47)
         mat( 48) = lmat( 48)
         mat( 49) = lmat( 49)
         mat( 50) = lmat( 50)
         mat( 51) = lmat( 51)
         mat( 52) = lmat( 52)
         mat( 53) = lmat( 53)
         mat( 54) = lmat( 54)
         mat( 55) = lmat( 55)
         mat( 56) = lmat( 56)
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 4) = mat( 4) - dti
         mat( 6) = mat( 6) - dti
         mat( 7) = mat( 7) - dti
         mat( 8) = mat( 8) - dti
         mat( 9) = mat( 9) - dti
         mat( 10) = mat( 10) - dti
         mat( 11) = mat( 11) - dti
         mat( 12) = mat( 12) - dti
         mat( 13) = mat( 13) - dti
         mat( 14) = mat( 14) - dti
         mat( 15) = mat( 15) - dti
         mat( 16) = mat( 16) - dti
         mat( 17) = mat( 17) - dti
         mat( 18) = mat( 18) - dti
         mat( 19) = mat( 19) - dti
         mat( 20) = mat( 20) - dti
         mat( 21) = mat( 21) - dti
         mat( 22) = mat( 22) - dti
         mat( 23) = mat( 23) - dti
         mat( 24) = mat( 24) - dti
         mat( 25) = mat( 25) - dti
         mat( 26) = mat( 26) - dti
         mat( 27) = mat( 27) - dti
         mat( 28) = mat( 28) - dti
         mat( 29) = mat( 29) - dti
         mat( 30) = mat( 30) - dti
         mat( 31) = mat( 31) - dti
         mat( 32) = mat( 32) - dti
         mat( 33) = mat( 33) - dti
         mat( 34) = mat( 34) - dti
         mat( 35) = mat( 35) - dti
         mat( 36) = mat( 36) - dti
         mat( 37) = mat( 37) - dti
         mat( 38) = mat( 38) - dti
         mat( 39) = mat( 39) - dti
         mat( 40) = mat( 40) - dti
         mat( 41) = mat( 41) - dti
         mat( 42) = mat( 42) - dti
         mat( 43) = mat( 43) - dti
         mat( 44) = mat( 44) - dti
         mat( 45) = mat( 45) - dti
         mat( 46) = mat( 46) - dti
         mat( 47) = mat( 47) - dti
         mat( 48) = mat( 48) - dti
         mat( 49) = mat( 49) - dti
         mat( 50) = mat( 50) - dti
         mat( 51) = mat( 51) - dti
         mat( 52) = mat( 52) - dti
         mat( 53) = mat( 53) - dti
         mat( 54) = mat( 54) - dti
         mat( 55) = mat( 55) - dti
         mat( 56) = mat( 56) - dti
      end subroutine nlnmat_finit
      end module mo_nln_matrix
