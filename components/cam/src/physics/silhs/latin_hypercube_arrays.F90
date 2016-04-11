!-------------------------------------------------------------------------------
! $Id: latin_hypercube_arrays.F90 7655 2015-04-29 19:22:26Z raut@uwm.edu $
!===============================================================================
module latin_hypercube_arrays

  use clubb_precision, only: &
    core_rknd

  implicit none

  public :: cleanup_latin_hypercube_arrays

  private

  integer, allocatable, dimension(:,:), public :: & 
    one_height_time_matrix ! matrix of rand ints

!$omp threadprivate(one_height_time_matrix)

  contains

  !-----------------------------------------------------------------------------
  subroutine cleanup_latin_hypercube_arrays( )

    ! Description:
    !   De-allocate latin hypercube arrays
    ! References:
    !   None
    !---------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( one_height_time_matrix ) ) then
      deallocate( one_height_time_matrix )
    end if

    return
  end subroutine cleanup_latin_hypercube_arrays

end module latin_hypercube_arrays
