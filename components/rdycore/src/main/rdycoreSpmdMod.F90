module RDycoreSpmdMod

  implicit none
  save
  private

  logical, public :: masterproc ! proc 0 logical for printing msgs
  integer, public :: iam        ! processor rank
  logical, public :: npes       ! number of processors
  integer, public :: mpicom     ! communicator group

  ! Public methods
  public RDycoreSpmdInit

contains

  subroutine RDycoreSpmdInit(mpicom_local)
    !
    implicit none
    !
    integer, intent (in) :: mpicom_local
    !
    integer :: ierr

    mpicom = mpicom_local

    ! Determine rank
    call mpi_comm_rank(mpicom, iam, ierr)

    ! Determine my rank
    if (iam == 0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors
    call mpi_comm_size(mpicom, npes, ierr)

  end subroutine RDycoreSpmdInit

end module RDycoreSpmdMod
