




      module mo_lu_factor

      private
      public :: lu_fac

      contains

      subroutine lu_fac01( lu )


      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: lu(:)

         lu(1) = 1._r8 / lu(1)

         lu(2) = 1._r8 / lu(2)

         lu(3) = 1._r8 / lu(3)

         lu(4) = 1._r8 / lu(4)

         lu(5) = 1._r8 / lu(5)

         lu(6) = 1._r8 / lu(6)

         lu(7) = 1._r8 / lu(7)

         lu(8) = 1._r8 / lu(8)

         lu(9) = 1._r8 / lu(9)

         lu(10) = 1._r8 / lu(10)

         lu(11) = 1._r8 / lu(11)

         lu(12) = 1._r8 / lu(12)

         lu(13) = 1._r8 / lu(13)

         lu(14) = 1._r8 / lu(14)

         lu(15) = 1._r8 / lu(15)

         lu(16) = 1._r8 / lu(16)

         lu(17) = 1._r8 / lu(17)

         lu(18) = 1._r8 / lu(18)

         lu(20) = 1._r8 / lu(20)
         lu(21) = lu(21) * lu(20)
         lu(24) = lu(24) - lu(21) * lu(22)
         lu(57) = lu(57) - lu(21) * lu(50)

         lu(23) = 1._r8 / lu(23)
         lu(24) = lu(24) * lu(23)
         lu(25) = lu(25) * lu(23)
         lu(57) = lu(57) - lu(24) * lu(51)
         lu(58) = lu(58) - lu(25) * lu(51)

         lu(26) = 1._r8 / lu(26)
         lu(27) = lu(27) * lu(26)
         lu(28) = lu(28) * lu(26)
         lu(57) = lu(57) - lu(27) * lu(52)
         lu(58) = lu(58) - lu(28) * lu(52)
         lu(65) = lu(65) - lu(27) * lu(62)
         lu(66) = lu(66) - lu(28) * lu(62)

         lu(29) = 1._r8 / lu(29)
         lu(30) = lu(30) * lu(29)
         lu(31) = lu(31) * lu(29)
         lu(32) = lu(32) * lu(29)
         lu(33) = lu(33) * lu(29)
         lu(34) = lu(34) * lu(29)
         lu(54) = lu(54) - lu(30) * lu(53)
         lu(57) = lu(57) - lu(31) * lu(53)
         lu(58) = lu(58) - lu(32) * lu(53)
         lu(59) = lu(59) - lu(33) * lu(53)
         lu(61) = lu(61) - lu(34) * lu(53)
         lu(86) = lu(86) - lu(30) * lu(85)
         lu(88) = lu(88) - lu(31) * lu(85)
         lu(89) = lu(89) - lu(32) * lu(85)
         lu(90) = lu(90) - lu(33) * lu(85)
         lu(92) = lu(92) - lu(34) * lu(85)

         lu(35) = 1._r8 / lu(35)
         lu(36) = lu(36) * lu(35)
         lu(37) = lu(37) * lu(35)
         lu(40) = lu(40) - lu(36) * lu(38)
         lu(41) = lu(41) - lu(37) * lu(38)
         lu(57) = lu(57) - lu(36) * lu(54)
         lu(58) = lu(58) - lu(37) * lu(54)
         lu(73) = - lu(36) * lu(70)
         lu(74) = lu(74) - lu(37) * lu(70)
         lu(80) = lu(80) - lu(36) * lu(78)
         lu(81) = lu(81) - lu(37) * lu(78)
         lu(88) = lu(88) - lu(36) * lu(86)
         lu(89) = lu(89) - lu(37) * lu(86)

         lu(39) = 1._r8 / lu(39)
         lu(40) = lu(40) * lu(39)
         lu(41) = lu(41) * lu(39)
         lu(42) = lu(42) * lu(39)
         lu(57) = lu(57) - lu(40) * lu(55)
         lu(58) = lu(58) - lu(41) * lu(55)
         lu(59) = lu(59) - lu(42) * lu(55)
         lu(65) = lu(65) - lu(40) * lu(63)
         lu(66) = lu(66) - lu(41) * lu(63)
         lu(67) = lu(67) - lu(42) * lu(63)
         lu(73) = lu(73) - lu(40) * lu(71)
         lu(74) = lu(74) - lu(41) * lu(71)
         lu(75) = lu(75) - lu(42) * lu(71)

         lu(44) = 1._r8 / lu(44)
         lu(45) = lu(45) * lu(44)
         lu(46) = lu(46) * lu(44)
         lu(47) = lu(47) * lu(44)
         lu(57) = lu(57) - lu(45) * lu(56)
         lu(60) = - lu(46) * lu(56)
         lu(61) = lu(61) - lu(47) * lu(56)
         lu(65) = lu(65) - lu(45) * lu(64)
         lu(68) = lu(68) - lu(46) * lu(64)
         lu(69) = lu(69) - lu(47) * lu(64)
         lu(73) = lu(73) - lu(45) * lu(72)
         lu(76) = lu(76) - lu(46) * lu(72)
         lu(77) = - lu(47) * lu(72)
         lu(80) = lu(80) - lu(45) * lu(79)
         lu(83) = lu(83) - lu(46) * lu(79)
         lu(84) = lu(84) - lu(47) * lu(79)
         lu(88) = lu(88) - lu(45) * lu(87)
         lu(91) = lu(91) - lu(46) * lu(87)
         lu(92) = lu(92) - lu(47) * lu(87)

         lu(57) = 1._r8 / lu(57)
         lu(58) = lu(58) * lu(57)
         lu(59) = lu(59) * lu(57)
         lu(60) = lu(60) * lu(57)
         lu(61) = lu(61) * lu(57)
         lu(66) = lu(66) - lu(58) * lu(65)
         lu(67) = lu(67) - lu(59) * lu(65)
         lu(68) = lu(68) - lu(60) * lu(65)
         lu(69) = lu(69) - lu(61) * lu(65)
         lu(74) = lu(74) - lu(58) * lu(73)
         lu(75) = lu(75) - lu(59) * lu(73)
         lu(76) = lu(76) - lu(60) * lu(73)
         lu(77) = lu(77) - lu(61) * lu(73)
         lu(81) = lu(81) - lu(58) * lu(80)
         lu(82) = lu(82) - lu(59) * lu(80)
         lu(83) = lu(83) - lu(60) * lu(80)
         lu(84) = lu(84) - lu(61) * lu(80)
         lu(89) = lu(89) - lu(58) * lu(88)
         lu(90) = lu(90) - lu(59) * lu(88)
         lu(91) = lu(91) - lu(60) * lu(88)
         lu(92) = lu(92) - lu(61) * lu(88)


      end subroutine lu_fac01

      subroutine lu_fac02( lu )


      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: lu(:)

         lu(66) = 1._r8 / lu(66)
         lu(67) = lu(67) * lu(66)
         lu(68) = lu(68) * lu(66)
         lu(69) = lu(69) * lu(66)
         lu(75) = lu(75) - lu(67) * lu(74)
         lu(76) = lu(76) - lu(68) * lu(74)
         lu(77) = lu(77) - lu(69) * lu(74)
         lu(82) = lu(82) - lu(67) * lu(81)
         lu(83) = lu(83) - lu(68) * lu(81)
         lu(84) = lu(84) - lu(69) * lu(81)
         lu(90) = lu(90) - lu(67) * lu(89)
         lu(91) = lu(91) - lu(68) * lu(89)
         lu(92) = lu(92) - lu(69) * lu(89)

         lu(75) = 1._r8 / lu(75)
         lu(76) = lu(76) * lu(75)
         lu(77) = lu(77) * lu(75)
         lu(83) = lu(83) - lu(76) * lu(82)
         lu(84) = lu(84) - lu(77) * lu(82)
         lu(91) = lu(91) - lu(76) * lu(90)
         lu(92) = lu(92) - lu(77) * lu(90)

         lu(83) = 1._r8 / lu(83)
         lu(84) = lu(84) * lu(83)
         lu(92) = lu(92) - lu(84) * lu(91)

         lu(92) = 1._r8 / lu(92)


      end subroutine lu_fac02

      subroutine lu_fac( lu )


      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: lu(:)

      call lu_fac01( lu )
      call lu_fac02( lu )

      end subroutine lu_fac

      end module mo_lu_factor
