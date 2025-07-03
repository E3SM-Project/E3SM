module RDycoreSpmdMod

  implicit none
  save
  private

  logical, public :: masterproc ! proc 0 logical for printing msgs
  integer, public :: iam        ! processor rank
  integer, public :: npes       ! number of processors
  integer, public :: mpicom_rof ! communicator group
  integer, public :: ROFID      ! mct compid

  ! Public methods
  public RDycoreSpmdInit

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

  subroutine RDycoreSpmdInit(mpicom_local)
    !
    implicit none
    !
    integer, intent (in) :: mpicom_local
    !
    integer :: ierr

    mpicom_rof = mpicom_local

    ! Determine rank
    call mpi_comm_rank(mpicom_rof, iam, ierr)

    ! Determine my rank
    if (iam == 0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors
    call mpi_comm_size(mpicom_rof, npes, ierr)

  end subroutine RDycoreSpmdInit

end module RDycoreSpmdMod
