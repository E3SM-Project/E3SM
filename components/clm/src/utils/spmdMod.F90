
module spmdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdMod
!
! !DESCRIPTION:
! SPMD initialization
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

  use shr_kind_mod   , only: r8 => shr_kind_r8
  implicit none

  private

#include <mpif.h>

  save

  ! Default settings valid even if there is no spmd 

  logical, public :: masterproc      ! proc 0 logical for printing msgs
  integer, public :: iam             ! processor number
  integer, public :: npes            ! number of processors for clm
  integer, public :: mpicom          ! communicator group for clm
  integer, public :: comp_id         ! component id

  !
  ! Public methods
  !
  public :: spmd_init                ! Initialization

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

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init( clm_mpicom )
!
! !INTERFACE:
  subroutine spmd_init( clm_mpicom, LNDID )
!
! !DESCRIPTION:
! MPI initialization (number of cpus, processes, tids, etc)
!
! !USES
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: clm_mpicom
    integer, intent(in) :: LNDID
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: i,j         ! indices
    integer :: ier         ! return error status
!-----------------------------------------------------------------------

    ! Initialize mpi communicator group

    mpicom = clm_mpicom

    comp_id = LNDID

    ! Get my processor id

    call mpi_comm_rank(mpicom, iam, ier)
    if (iam==0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    ! Get number of processors

    call mpi_comm_size(mpicom, npes, ier)

  end subroutine spmd_init

end module spmdMod
