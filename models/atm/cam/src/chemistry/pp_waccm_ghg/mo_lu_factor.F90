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
         lu(2) = lu(2) * lu(1)
         lu(3) = 1._r8 / lu(3)
         lu(4) = 1._r8 / lu(4)
         lu(5) = 1._r8 / lu(5)
         lu(6) = 1._r8 / lu(6)
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
