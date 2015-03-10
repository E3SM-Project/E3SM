module mpi_wrapper_mod
  ! This module wraps some mpi stuff for the sake of unit testing.
  !
  ! NOTE(wjs, 2015-01-16) Currently this has been set up with serial unit tests in mind -
  ! using a serial build of pFUnit and the mpi-serial library (although this probably also
  ! works with a true mpi library). More thought will be needed about how to make this
  ! work with either serial or parallel pFUnit tests. In particular: I think parallel
  ! pFUnit tests would do their own mpi_init and mpi_finalize calls, so it would be an
  ! error for unit tests to make these calls themselves. That explains the motivation for
  ! setting up this module: When the time comes to introduce that generality, I didn't
  ! want to have to add logic in every test module regarding whether to call mpi_init and
  ! mpi_finalize: instead, they can safely call these wrappers, which will do the right
  ! thing.

  implicit none
  private

#include <mpif.h>
  
  public :: mpi_init_wrapper
  public :: mpi_finalize_wrapper

contains

  !-----------------------------------------------------------------------
  subroutine mpi_init_wrapper()
    !
    ! !DESCRIPTION:
    ! Wrapper to mpi_init. Unit tests should call this rather than calling mpi_init
    ! directly.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: ier
    
    character(len=*), parameter :: subname = 'mpi_init_wrapper'
    !-----------------------------------------------------------------------

    call mpi_init(ier)
    if (ier /= MPI_SUCCESS) then
       print *, 'ERROR in mpi_init'
       stop
    end if
    
  end subroutine mpi_init_wrapper

  !-----------------------------------------------------------------------
  subroutine mpi_finalize_wrapper()
    !
    ! !DESCRIPTION:
    ! Wrapper to mpi_finalize. Unit tests should call this rather than calling
    ! mpi_finalize directly.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: ier
    
    character(len=*), parameter :: subname = 'mpi_finalize_wrapper'
    !-----------------------------------------------------------------------

    call mpi_finalize(ier)
    if (ier /= MPI_SUCCESS) then
       print *, 'ERROR in mpi_finalize'
       stop
    end if
    
  end subroutine mpi_finalize_wrapper

end module mpi_wrapper_mod
