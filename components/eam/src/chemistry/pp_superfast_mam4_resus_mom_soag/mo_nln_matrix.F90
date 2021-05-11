













      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: nlnmat

      contains

      subroutine     nlnmat01( mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(inout) ::  mat(nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------


         mat(103) = -(rxt(7)*y(2) + rxt(8)*y(3) + rxt(12)*y(5) + rxt(28)*y(14))
         mat(87) = -rxt(7)*y(1)
         mat(95) = -rxt(8)*y(1)
         mat(65) = -rxt(12)*y(1)
         mat(45) = -rxt(28)*y(1)

         mat(85) = -(rxt(7)*y(1) + rxt(9)*y(3) + rxt(11)*y(4) + rxt(14)*y(6) + rxt(17) &
                      *y(9) + rxt(19)*y(11) + (rxt(24) + rxt(25)) * y(12) + rxt(26) &
                      *y(13) + rxt(27)*y(14))
         mat(101) = -rxt(7)*y(2)
         mat(93) = -rxt(9)*y(2)
         mat(38) = -rxt(11)*y(2)
         mat(57) = -rxt(14)*y(2)
         mat(47) = -rxt(17)*y(2)
         mat(52) = -rxt(19)*y(2)
         mat(35) = -(rxt(24) + rxt(25)) * y(2)
         mat(31) = -rxt(26)*y(2)
         mat(43) = -rxt(27)*y(2)

         mat(101) = mat(101) + rxt(8)*y(3)
         mat(93) = mat(93) + rxt(8)*y(1) + rxt(13)*y(5)
         mat(63) = rxt(13)*y(3)

         mat(94) = -(rxt(8)*y(1) + rxt(9)*y(2) + 4._r8*rxt(10)*y(3) + rxt(13)*y(5) &
                      + rxt(18)*y(10))
         mat(102) = -rxt(8)*y(3)
         mat(86) = -rxt(9)*y(3)
         mat(64) = -rxt(13)*y(3)
         mat(72) = -rxt(18)*y(3)

         mat(102) = mat(102) + rxt(7)*y(2) + .060_r8*rxt(28)*y(14)
         mat(86) = mat(86) + rxt(7)*y(1) + rxt(11)*y(4) + rxt(17)*y(9)  &
                      + .500_r8*rxt(25)*y(12)
         mat(39) = rxt(11)*y(2)
         mat(64) = mat(64) + rxt(21)*y(10)
         mat(48) = rxt(17)*y(2)
         mat(72) = mat(72) + rxt(21)*y(5) + 1.600_r8*rxt(22)*y(10)
         mat(36) = .500_r8*rxt(25)*y(2)
         mat(44) = .060_r8*rxt(28)*y(1)

         mat(37) = -(rxt(11)*y(2))
         mat(78) = -rxt(11)*y(4)

         mat(88) = 2.000_r8*rxt(10)*y(3)

         mat(61) = -(rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10))
         mat(99) = -rxt(12)*y(5)
         mat(91) = -rxt(13)*y(5)
         mat(69) = -rxt(21)*y(5)

         mat(55) = -(rxt(14)*y(2))
         mat(82) = -rxt(14)*y(6)

         mat(98) = rxt(12)*y(5)
         mat(90) = rxt(13)*y(5)
         mat(60) = rxt(12)*y(1) + rxt(13)*y(3) + rxt(21)*y(10)
         mat(68) = rxt(21)*y(5)


         mat(74) = rxt(14)*y(6)
         mat(54) = rxt(14)*y(2)

         mat(46) = -(rxt(17)*y(2))
         mat(80) = -rxt(17)*y(9)

         mat(97) = .870_r8*rxt(28)*y(14)
         mat(80) = mat(80) + rxt(20)*y(11)
         mat(59) = rxt(21)*y(10)
         mat(66) = rxt(21)*y(5) + 4.000_r8*rxt(22)*y(10)
         mat(49) = rxt(20)*y(2)
         mat(41) = .870_r8*rxt(28)*y(1)

         mat(70) = -(rxt(18)*y(3) + rxt(21)*y(5) + 4._r8*rxt(22)*y(10))
         mat(92) = -rxt(18)*y(10)
         mat(62) = -rxt(21)*y(10)

         mat(100) = 1.860_r8*rxt(28)*y(14)
         mat(84) = rxt(19)*y(11)
         mat(51) = rxt(19)*y(2)
         mat(42) = 1.860_r8*rxt(28)*y(1)

         mat(50) = -((rxt(19) + rxt(20)) * y(2))
         mat(81) = -(rxt(19) + rxt(20)) * y(11)

         mat(89) = rxt(18)*y(10)
         mat(67) = rxt(18)*y(3)

         mat(34) = -((rxt(24) + rxt(25)) * y(2))
         mat(77) = -(rxt(24) + rxt(25)) * y(12)

         mat(30) = -(rxt(26)*y(2))
         mat(76) = -rxt(26)*y(13)

         mat(76) = mat(76) + (rxt(24)+.500_r8*rxt(25))*y(12)
         mat(33) = (rxt(24)+.500_r8*rxt(25))*y(2)

         mat(40) = -(rxt(27)*y(2) + rxt(28)*y(1))
         mat(79) = -rxt(27)*y(14)
         mat(96) = -rxt(28)*y(14)


         mat(75) = rxt(26)*y(13)
         mat(29) = rxt(26)*y(2)




























      end subroutine     nlnmat01

      subroutine     nlnmat_finit( mat, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  dti
      real(r8), intent(in)    ::  lmat(nzcnt)
      real(r8), intent(inout) ::  mat(nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------


         mat(   1) = lmat(   1)
         mat(   2) = lmat(   2)
         mat(   3) = lmat(   3)
         mat(   4) = lmat(   4)
         mat(   5) = lmat(   5)
         mat(   6) = lmat(   6)
         mat(   7) = lmat(   7)
         mat(   8) = lmat(   8)
         mat(   9) = lmat(   9)
         mat(  10) = lmat(  10)
         mat(  11) = lmat(  11)
         mat(  12) = lmat(  12)
         mat(  13) = lmat(  13)
         mat(  14) = lmat(  14)
         mat(  15) = lmat(  15)
         mat(  16) = lmat(  16)
         mat(  17) = lmat(  17)
         mat(  18) = lmat(  18)
         mat(  19) = lmat(  19)
         mat(  20) = lmat(  20)
         mat(  21) = lmat(  21)
         mat(  22) = lmat(  22)
         mat(  23) = lmat(  23)
         mat(  24) = lmat(  24)
         mat(  25) = lmat(  25)
         mat(  26) = lmat(  26)
         mat(  27) = lmat(  27)
         mat(  28) = lmat(  28)
         mat(  30) = mat(  30) + lmat(  30)
         mat(  32) = lmat(  32)
         mat(  33) = mat(  33) + lmat(  33)
         mat(  34) = mat(  34) + lmat(  34)
         mat(  37) = mat(  37) + lmat(  37)
         mat(  38) = mat(  38) + lmat(  38)
         mat(  40) = mat(  40) + lmat(  40)
         mat(  46) = mat(  46) + lmat(  46)
         mat(  48) = mat(  48) + lmat(  48)
         mat(  49) = mat(  49) + lmat(  49)
         mat(  50) = mat(  50) + lmat(  50)
         mat(  52) = mat(  52) + lmat(  52)
         mat(  53) = lmat(  53)
         mat(  54) = mat(  54) + lmat(  54)
         mat(  55) = mat(  55) + lmat(  55)
         mat(  56) = lmat(  56)
         mat(  58) = lmat(  58)
         mat(  61) = mat(  61) + lmat(  61)
         mat(  70) = mat(  70) + lmat(  70)
         mat(  84) = mat(  84) + lmat(  84)
         mat(  85) = mat(  85) + lmat(  85)
         mat(  86) = mat(  86) + lmat(  86)
         mat(  94) = mat(  94) + lmat(  94)
         mat( 101) = mat( 101) + lmat( 101)
         mat( 103) = mat( 103) + lmat( 103)
         mat(  71) = 0._r8
         mat(  73) = 0._r8
         mat(  83) = 0._r8
         mat(   1) = mat(   1) - dti
         mat(   2) = mat(   2) - dti
         mat(   3) = mat(   3) - dti
         mat(   4) = mat(   4) - dti
         mat(   5) = mat(   5) - dti
         mat(   6) = mat(   6) - dti
         mat(   7) = mat(   7) - dti
         mat(   8) = mat(   8) - dti
         mat(   9) = mat(   9) - dti
         mat(  10) = mat(  10) - dti
         mat(  11) = mat(  11) - dti
         mat(  12) = mat(  12) - dti
         mat(  13) = mat(  13) - dti
         mat(  14) = mat(  14) - dti
         mat(  15) = mat(  15) - dti
         mat(  16) = mat(  16) - dti
         mat(  17) = mat(  17) - dti
         mat(  18) = mat(  18) - dti
         mat(  19) = mat(  19) - dti
         mat(  20) = mat(  20) - dti
         mat(  21) = mat(  21) - dti
         mat(  22) = mat(  22) - dti
         mat(  23) = mat(  23) - dti
         mat(  24) = mat(  24) - dti
         mat(  25) = mat(  25) - dti
         mat(  26) = mat(  26) - dti
         mat(  27) = mat(  27) - dti
         mat(  28) = mat(  28) - dti
         mat(  30) = mat(  30) - dti
         mat(  34) = mat(  34) - dti
         mat(  37) = mat(  37) - dti
         mat(  40) = mat(  40) - dti
         mat(  46) = mat(  46) - dti
         mat(  50) = mat(  50) - dti
         mat(  55) = mat(  55) - dti
         mat(  61) = mat(  61) - dti
         mat(  70) = mat(  70) - dti
         mat(  85) = mat(  85) - dti
         mat(  94) = mat(  94) - dti
         mat( 103) = mat( 103) - dti

      end subroutine nlnmat_finit

      subroutine     nlnmat( mat, y, rxt, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot,     nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      real(r8), intent(in)    ::  dti
      real(r8), intent(in)    ::  lmat(nzcnt)
      real(r8), intent(in)    ::  y(gas_pcnst)
      real(r8), intent(in)    ::  rxt(rxntot)
      real(r8), intent(inout) ::  mat(nzcnt)

      call     nlnmat01( mat, y, rxt )
      call     nlnmat_finit( mat, lmat, dti )

      end subroutine nlnmat

      end module mo_nln_matrix
