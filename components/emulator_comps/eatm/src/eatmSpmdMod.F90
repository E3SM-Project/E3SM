module eatmSpmdMod

  use shr_sys_mod,  only: shr_sys_abort

  implicit none
  save
  private

  logical, public :: masterproc ! proc 0 logical for printing msgs
  integer, public :: iam        ! processor rank
  integer, public :: npes       ! number of processors
  integer, public :: mpicom_atm ! communicator group
  integer, public :: COMPID     ! mct compid

  ! Public methods
  public eatmSpmdInit

  !
  ! Values from mpif.h that can be used
  !
  public :: MPI_INTEGER
  public :: MPI_REAL8
  public :: MPI_LOGICAL
  public :: MPI_SUM
  public :: MPI_MIN
  public :: MPI_MAX
  public :: MPI_LOR
  public :: MPI_STATUS_SIZE
  public :: MPI_ANY_SOURCE
  public :: MPI_CHARACTER
  public :: MPI_COMM_WORLD
  public :: MPI_MAX_PROCESSOR_NAME

#include <mpif.h>

contains

  subroutine eatmSpmdInit(mpicom_local)
    !
    implicit none
    !
    integer, intent (in) :: mpicom_local
    !
    character(*), parameter :: subName = "(eatmSpmdInit) "
    integer :: ierr

    mpicom_atm = mpicom_local

    ! Determine rank
    call mpi_comm_rank(mpicom_atm, iam, ierr)

    ! Determine my rank
    if (iam == 0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors
    call mpi_comm_size(mpicom_atm, npes, ierr)

    ! EATM only runs as a serial model currently, so abort if more than one
    ! pe in mpi communicator
    if (npes.gt.1) then
      call shr_sys_abort(trim(subname)//' ERROR EATM is serial but multiple pes in mpicomm')
    endif

  end subroutine eatmSpmdInit

end module eatmSpmdMod
