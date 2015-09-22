module mct_wrapper_mod
  ! This module provides some variables and convenience functions for the sake of unit
  ! tests that use mct.
  !
  ! Any test that uses mct should call mct_init in its initialization, and mct_clean in
  ! its teardown.
  
  implicit none
  private

#include <mpif.h>

  public :: mct_init  ! initialize data structures needed to use mct 
  public :: mct_clean ! clean up mct data structures that were set up by mct_init

  ! MPI communicator that can be used wherever mct routines expect a communicator
  integer, parameter, public :: mct_communicator = MPI_COMM_WORLD

  ! value that can be used wherever mct routines expect a component ID
  integer, parameter, public :: mct_compid = 1
  
contains
  
  !-----------------------------------------------------------------------
  subroutine mct_init()
    !
    ! !DESCRIPTION:
    ! Initializes data structures needed to use mct.
    !
    ! Expects that mpi_init has already been called.
    !
    ! !USES:
    use seq_comm_mct, only : seq_comm_init
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'mct_init'
    !-----------------------------------------------------------------------

    call seq_comm_init(Comm_in = mct_communicator, nmlfile = ' ')
  end subroutine mct_init

  !-----------------------------------------------------------------------
  subroutine mct_clean()
    !
    ! !DESCRIPTION:
    ! Cleans up mct data structures that were set up by mct_init.
    !
    ! !USES:
    use seq_comm_mct, only : seq_comm_clean
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'mct_clean'
    !-----------------------------------------------------------------------

    call seq_comm_clean()
    
  end subroutine mct_clean


end module mct_wrapper_mod
