      module mo_lu_factor
      private
      public :: lu_fac
      contains
      subroutine lu_fac( lu )
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: lu(:)
      end subroutine lu_fac
      end module mo_lu_factor
