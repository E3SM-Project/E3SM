




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


         mat(92) = -(rxt(7)*y(2) + rxt(8)*y(3) + rxt(12)*y(5) + rxt(28)*y(14))
         mat(61) = -rxt(7)*y(1)
         mat(69) = -rxt(8)*y(1)
         mat(84) = -rxt(12)*y(1)
         mat(34) = -rxt(28)*y(1)

         mat(57) = -(rxt(7)*y(1) + rxt(9)*y(3) + rxt(11)*y(4) + rxt(14)*y(6) + rxt(17) &
                      *y(9) + rxt(19)*y(11) + (rxt(24) + rxt(25)) * y(12) + rxt(26) &
                      *y(13) + rxt(27)*y(14))
         mat(88) = -rxt(7)*y(2)
         mat(65) = -rxt(9)*y(2)
         mat(27) = -rxt(11)*y(2)
         mat(45) = -rxt(14)*y(2)
         mat(36) = -rxt(17)*y(2)
         mat(40) = -rxt(19)*y(2)
         mat(24) = -(rxt(24) + rxt(25)) * y(2)
         mat(21) = -rxt(26)*y(2)
         mat(31) = -rxt(27)*y(2)

         mat(88) = mat(88) + rxt(8)*y(3)
         mat(65) = mat(65) + rxt(8)*y(1) + rxt(13)*y(5)
         mat(80) = rxt(13)*y(3)

         mat(66) = -(rxt(8)*y(1) + rxt(9)*y(2) + 4._r8*rxt(10)*y(3) + rxt(13)*y(5) &
                      + rxt(18)*y(10))
         mat(89) = -rxt(8)*y(3)
         mat(58) = -rxt(9)*y(3)
         mat(81) = -rxt(13)*y(3)
         mat(74) = -rxt(18)*y(3)

         mat(89) = mat(89) + rxt(7)*y(2) + .060_r8*rxt(28)*y(14)
         mat(58) = mat(58) + rxt(7)*y(1) + rxt(11)*y(4) + rxt(17)*y(9) &
                      + .500_r8*rxt(25)*y(12)
         mat(28) = rxt(11)*y(2)
         mat(81) = mat(81) + rxt(21)*y(10)
         mat(37) = rxt(17)*y(2)
         mat(74) = mat(74) + rxt(21)*y(5) + 1.600_r8*rxt(22)*y(10)
         mat(25) = .500_r8*rxt(25)*y(2)
         mat(32) = .060_r8*rxt(28)*y(1)

         mat(26) = -(rxt(11)*y(2))
         mat(52) = -rxt(11)*y(4)

         mat(62) = 2.000_r8*rxt(10)*y(3)

         mat(83) = -(rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10))
         mat(91) = -rxt(12)*y(5)
         mat(68) = -rxt(13)*y(5)
         mat(76) = -rxt(21)*y(5)

         mat(44) = -(rxt(14)*y(2))
         mat(56) = -rxt(14)*y(6)

         mat(87) = rxt(12)*y(5)
         mat(64) = rxt(13)*y(5)
         mat(79) = rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10)
         mat(72) = rxt(21)*y(5)


         mat(48) = rxt(14)*y(6)
         mat(43) = rxt(14)*y(2)

         mat(35) = -(rxt(17)*y(2))
         mat(54) = -rxt(17)*y(9)

         mat(86) = .870_r8*rxt(28)*y(14)
         mat(54) = mat(54) + rxt(20)*y(11)
         mat(78) = rxt(21)*y(10)
         mat(70) = rxt(21)*y(5) + 4.000_r8*rxt(22)*y(10)
         mat(38) = rxt(20)*y(2)
         mat(30) = .870_r8*rxt(28)*y(1)

         mat(75) = -(rxt(18)*y(3) + rxt(21)*y(5) + 4._r8*rxt(22)*y(10))
         mat(67) = -rxt(18)*y(10)
         mat(82) = -rxt(21)*y(10)

         mat(90) = 1.860_r8*rxt(28)*y(14)
         mat(59) = rxt(19)*y(11)
         mat(42) = rxt(19)*y(2)
         mat(33) = 1.860_r8*rxt(28)*y(1)

         mat(39) = -((rxt(19) + rxt(20)) * y(2))
         mat(55) = -(rxt(19) + rxt(20)) * y(11)

         mat(63) = rxt(18)*y(10)
         mat(71) = rxt(18)*y(3)

         mat(23) = -((rxt(24) + rxt(25)) * y(2))
         mat(51) = -(rxt(24) + rxt(25)) * y(12)

         mat(20) = -(rxt(26)*y(2))
         mat(50) = -rxt(26)*y(13)

         mat(50) = mat(50) + (rxt(24)+.500_r8*rxt(25))*y(12)
         mat(22) = (rxt(24)+.500_r8*rxt(25))*y(2)

         mat(29) = -(rxt(27)*y(2) + rxt(28)*y(1))
         mat(53) = -rxt(27)*y(14)
         mat(85) = -rxt(28)*y(14)


         mat(49) = rxt(26)*y(13)
         mat(19) = rxt(26)*y(2)
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
         mat( 20) = mat( 20) + lmat( 20)
         mat( 23) = mat( 23) + lmat( 23)
         mat( 26) = mat( 26) + lmat( 26)
         mat( 27) = mat( 27) + lmat( 27)
         mat( 29) = mat( 29) + lmat( 29)
         mat( 35) = mat( 35) + lmat( 35)
         mat( 37) = mat( 37) + lmat( 37)
         mat( 38) = mat( 38) + lmat( 38)
         mat( 39) = mat( 39) + lmat( 39)
         mat( 40) = mat( 40) + lmat( 40)
         mat( 41) = lmat( 41)
         mat( 43) = mat( 43) + lmat( 43)
         mat( 44) = mat( 44) + lmat( 44)
         mat( 46) = lmat( 46)
         mat( 47) = lmat( 47)
         mat( 57) = mat( 57) + lmat( 57)
         mat( 58) = mat( 58) + lmat( 58)
         mat( 59) = mat( 59) + lmat( 59)
         mat( 66) = mat( 66) + lmat( 66)
         mat( 75) = mat( 75) + lmat( 75)
         mat( 83) = mat( 83) + lmat( 83)
         mat( 88) = mat( 88) + lmat( 88)
         mat( 92) = mat( 92) + lmat( 92)
         mat( 60) = 0._r8
         mat( 73) = 0._r8
         mat( 77) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 3) = mat( 3) - dti
         mat( 4) = mat( 4) - dti
         mat( 5) = mat( 5) - dti
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
         mat( 20) = mat( 20) - dti
         mat( 23) = mat( 23) - dti
         mat( 26) = mat( 26) - dti
         mat( 29) = mat( 29) - dti
         mat( 35) = mat( 35) - dti
         mat( 39) = mat( 39) - dti
         mat( 44) = mat( 44) - dti
         mat( 57) = mat( 57) - dti
         mat( 66) = mat( 66) - dti
         mat( 75) = mat( 75) - dti
         mat( 83) = mat( 83) - dti
         mat( 92) = mat( 92) - dti
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
