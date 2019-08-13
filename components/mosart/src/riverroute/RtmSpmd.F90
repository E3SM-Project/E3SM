
module RtmSpmd

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmSpmd
!
! !DESCRIPTION:
! RTM SPMD initialization
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------
  implicit none
  save
  private

  ! Default settings valid even if there is no spmd 

  logical, public :: masterproc      ! proc 0 logical for printing msgs
  integer, public :: iam             ! processor number
  integer, public :: npes            ! number of processors for rtm
  integer, public :: mpicom_rof      ! communicator group for rtm
  integer, public :: ROFID           ! mct compid
  integer, public, parameter :: MASTERTASK=0 ! the value of iam which is assigned 
                                             ! the masterproc duties

  !
  ! Public methods
  !
  public :: RtmSpmdInit                ! Initialization

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

!-----------------------------------------------------------------------

  subroutine RtmSpmdInit(mpicom)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! MPI initialization (number of processes, etc)
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: mpicom
    !
    ! !LOCAL VARIABLES:
    integer :: ier  ! return error status
    !-----------------------------------------------------------------------

    ! Initialize mpi communicator group

    mpicom_rof = mpicom

    ! Get my processor id

    call mpi_comm_rank(mpicom_rof, iam, ier)
    if (iam == MASTERTASK) then 
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors

    call mpi_comm_size(mpicom_rof, npes, ier)

  end subroutine RtmSpmdInit

end module RtmSpmd
