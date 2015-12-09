






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

         lu(4) = 1._r8 / lu(4)

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

         lu(19) = 1._r8 / lu(19)

         lu(20) = 1._r8 / lu(20)

         lu(21) = 1._r8 / lu(21)

         lu(22) = 1._r8 / lu(22)

         lu(23) = 1._r8 / lu(23)

         lu(24) = 1._r8 / lu(24)

         lu(25) = 1._r8 / lu(25)

         lu(26) = 1._r8 / lu(26)

         lu(27) = 1._r8 / lu(27)

         lu(28) = 1._r8 / lu(28)

         lu(29) = 1._r8 / lu(29)

         lu(30) = 1._r8 / lu(30)

         lu(31) = 1._r8 / lu(31)

         lu(32) = 1._r8 / lu(32)

         lu(33) = 1._r8 / lu(33)

         lu(34) = 1._r8 / lu(34)

         lu(35) = 1._r8 / lu(35)

         lu(36) = 1._r8 / lu(36)

         lu(37) = 1._r8 / lu(37)

         lu(38) = 1._r8 / lu(38)

         lu(39) = 1._r8 / lu(39)

         lu(40) = 1._r8 / lu(40)

         lu(41) = 1._r8 / lu(41)

         lu(42) = 1._r8 / lu(42)

         lu(43) = 1._r8 / lu(43)

         lu(44) = 1._r8 / lu(44)

         lu(45) = 1._r8 / lu(45)

         lu(46) = 1._r8 / lu(46)

         lu(47) = 1._r8 / lu(47)

         lu(48) = 1._r8 / lu(48)

         lu(49) = 1._r8 / lu(49)

         lu(50) = 1._r8 / lu(50)

         lu(51) = 1._r8 / lu(51)

         lu(52) = 1._r8 / lu(52)

         lu(53) = 1._r8 / lu(53)


      end subroutine lu_fac01

      subroutine lu_fac02( lu )


      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: lu(:)

         lu(54) = 1._r8 / lu(54)

         lu(55) = 1._r8 / lu(55)

         lu(56) = 1._r8 / lu(56)


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
