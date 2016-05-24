      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(41) = -(rxt(7)*y(2) + rxt(8)*y(3) + rxt(12)*y(5) + rxt(28)*y(13) + rxt(30) &
                      *y(15))
         mat(71) = -rxt(7)*y(1)
         mat(50) = -rxt(8)*y(1)
         mat(57) = -rxt(12)*y(1)
         mat(14) = -rxt(28)*y(1)
         mat(19) = -rxt(30)*y(1)
         mat(74) = -(rxt(7)*y(1) + rxt(9)*y(3) + rxt(11)*y(4) + rxt(14)*y(6) + rxt(17) &
                      *y(9) + rxt(19)*y(11) + (rxt(24) + rxt(25)) * y(12) + rxt(26) &
                      *y(13) + rxt(29)*y(15))
         mat(44) = -rxt(7)*y(2)
         mat(53) = -rxt(9)*y(2)
         mat(10) = -rxt(11)*y(2)
         mat(35) = -rxt(14)*y(2)
         mat(25) = -rxt(17)*y(2)
         mat(29) = -rxt(19)*y(2)
         mat(5) = -(rxt(24) + rxt(25)) * y(2)
         mat(16) = -rxt(26)*y(2)
         mat(21) = -rxt(29)*y(2)
         mat(44) = mat(44) + rxt(8)*y(3)
         mat(53) = mat(53) + rxt(8)*y(1) + rxt(13)*y(5)
         mat(60) = rxt(13)*y(3)
         mat(51) = -(rxt(8)*y(1) + rxt(9)*y(2) + 4._r8*rxt(10)*y(3) + rxt(13)*y(5) &
                      + rxt(18)*y(10))
         mat(42) = -rxt(8)*y(3)
         mat(72) = -rxt(9)*y(3)
         mat(58) = -rxt(13)*y(3)
         mat(80) = -rxt(18)*y(3)
         mat(42) = mat(42) + rxt(7)*y(2) + .060_r8*rxt(30)*y(15)
         mat(72) = mat(72) + rxt(7)*y(1) + rxt(11)*y(4) + rxt(17)*y(9)
         mat(9) = rxt(11)*y(2)
         mat(58) = mat(58) + rxt(21)*y(10)
         mat(24) = rxt(17)*y(2)
         mat(80) = mat(80) + rxt(21)*y(5) + 1.600_r8*rxt(22)*y(10)
         mat(20) = .060_r8*rxt(30)*y(1)
         mat(7) = -(rxt(11)*y(2) + rxt(27)*y(13))
         mat(65) = -rxt(11)*y(4)
         mat(12) = -rxt(27)*y(4)
         mat(46) = 2.000_r8*rxt(10)*y(3)
         mat(59) = -(rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10))
         mat(43) = -rxt(12)*y(5)
         mat(52) = -rxt(13)*y(5)
         mat(81) = -rxt(21)*y(5)
         mat(32) = -(rxt(14)*y(2))
         mat(70) = -rxt(14)*y(6)
         mat(40) = rxt(12)*y(5)
         mat(49) = rxt(13)*y(5)
         mat(56) = rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10)
         mat(78) = rxt(21)*y(5)
         mat(62) = rxt(14)*y(6)
         mat(31) = rxt(14)*y(2)
         mat(23) = -(rxt(17)*y(2))
         mat(68) = -rxt(17)*y(9)
         mat(39) = .870_r8*rxt(30)*y(15)
         mat(68) = mat(68) + rxt(20)*y(11)
         mat(55) = rxt(21)*y(10)
         mat(76) = rxt(21)*y(5) + 4.000_r8*rxt(22)*y(10)
         mat(26) = rxt(20)*y(2)
         mat(18) = .870_r8*rxt(30)*y(1)
         mat(83) = -(rxt(18)*y(3) + rxt(21)*y(5) + 4._r8*rxt(22)*y(10))
         mat(54) = -rxt(18)*y(10)
         mat(61) = -rxt(21)*y(10)
         mat(45) = 1.860_r8*rxt(30)*y(15)
         mat(75) = rxt(19)*y(11)
         mat(30) = rxt(19)*y(2)
         mat(22) = 1.860_r8*rxt(30)*y(1)
         mat(27) = -((rxt(19) + rxt(20)) * y(2))
         mat(69) = -(rxt(19) + rxt(20)) * y(11)
         mat(48) = rxt(18)*y(10)
         mat(77) = rxt(18)*y(3)
         mat(3) = -((rxt(24) + rxt(25)) * y(2))
         mat(64) = -(rxt(24) + rxt(25)) * y(12)
         mat(13) = -(rxt(26)*y(2) + rxt(27)*y(4) + rxt(28)*y(1))
         mat(66) = -rxt(26)*y(13)
         mat(8) = -rxt(27)*y(13)
         mat(37) = -rxt(28)*y(13)
         mat(66) = mat(66) + (rxt(24)+.750_r8*rxt(25))*y(12)
         mat(4) = (rxt(24)+.750_r8*rxt(25))*y(2)
         mat(36) = rxt(28)*y(13)
         mat(63) = rxt(26)*y(13)
         mat(6) = rxt(27)*y(13)
         mat(11) = rxt(28)*y(1) + rxt(26)*y(2) + rxt(27)*y(4)
         mat(17) = -(rxt(29)*y(2) + rxt(30)*y(1))
         mat(67) = -rxt(29)*y(15)
         mat(38) = -rxt(30)*y(15)
      end subroutine nlnmat01
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
         mat( 3) = mat( 3) + lmat( 3)
         mat( 7) = mat( 7) + lmat( 7)
         mat( 10) = mat( 10) + lmat( 10)
         mat( 13) = mat( 13) + lmat( 13)
         mat( 17) = mat( 17) + lmat( 17)
         mat( 23) = mat( 23) + lmat( 23)
         mat( 24) = mat( 24) + lmat( 24)
         mat( 26) = mat( 26) + lmat( 26)
         mat( 27) = mat( 27) + lmat( 27)
         mat( 28) = lmat( 28)
         mat( 29) = mat( 29) + lmat( 29)
         mat( 31) = mat( 31) + lmat( 31)
         mat( 32) = mat( 32) + lmat( 32)
         mat( 33) = lmat( 33)
         mat( 34) = lmat( 34)
         mat( 41) = mat( 41) + lmat( 41)
         mat( 44) = mat( 44) + lmat( 44)
         mat( 51) = mat( 51) + lmat( 51)
         mat( 59) = mat( 59) + lmat( 59)
         mat( 72) = mat( 72) + lmat( 72)
         mat( 74) = mat( 74) + lmat( 74)
         mat( 75) = mat( 75) + lmat( 75)
         mat( 83) = mat( 83) + lmat( 83)
         mat( 15) = 0._r8
         mat( 47) = 0._r8
         mat( 73) = 0._r8
         mat( 79) = 0._r8
         mat( 82) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 3) = mat( 3) - dti
         mat( 7) = mat( 7) - dti
         mat( 13) = mat( 13) - dti
         mat( 17) = mat( 17) - dti
         mat( 23) = mat( 23) - dti
         mat( 27) = mat( 27) - dti
         mat( 32) = mat( 32) - dti
         mat( 41) = mat( 41) - dti
         mat( 51) = mat( 51) - dti
         mat( 59) = mat( 59) - dti
         mat( 74) = mat( 74) - dti
         mat( 83) = mat( 83) - dti
      end subroutine nlnmat_finit
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
      call nlnmat01( mat, y, rxt )
      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
