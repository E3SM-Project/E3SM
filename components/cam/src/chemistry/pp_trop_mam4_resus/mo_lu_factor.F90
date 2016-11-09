




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


      end subroutine lu_fac01

      subroutine lu_fac( lu )


      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: lu(:)

      call lu_fac01( lu )

      end subroutine lu_fac

      end module mo_lu_factor
