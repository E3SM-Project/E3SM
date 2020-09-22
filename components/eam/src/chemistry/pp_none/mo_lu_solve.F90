      module mo_lu_solve
      private
      public :: lu_slv
      contains
      subroutine lu_slv( lu, b )
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!-----------------------------------------------------------------------
! ... Dummy args
!-----------------------------------------------------------------------
      real(r8), intent(in) :: lu(:)
      real(r8), intent(inout) :: b(:)
      end subroutine lu_slv
      end module mo_lu_solve
