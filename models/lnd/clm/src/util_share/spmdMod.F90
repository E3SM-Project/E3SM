
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

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only: iulog
  implicit none
  save
  private

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

#include <mpif.h>  

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
    integer :: mylength    ! my processor length
    logical :: mpi_running ! temporary
    integer, allocatable :: length(:)
    integer, allocatable :: displ(:)
    character*(MPI_MAX_PROCESSOR_NAME), allocatable :: procname(:)
    character*(MPI_MAX_PROCESSOR_NAME)              :: myprocname
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

    ! Get my processor names

    allocate (length(0:npes-1), displ(0:npes-1), procname(0:npes-1))

    call mpi_get_processor_name (myprocname, mylength, ier)
    call mpi_allgather(mylength,1,MPI_INTEGER,length,1,MPI_INTEGER,mpicom,ier)

    do i = 0,npes-1
       displ(i)=i*MPI_MAX_PROCESSOR_NAME
    end do
    call mpi_gatherv (myprocname,mylength,MPI_CHARACTER, &
                      procname,length,displ,MPI_CHARACTER,0,mpicom,ier)
    if (masterproc) then
       write(iulog,100)npes
       write(iulog,200)
       write(iulog,220)
       do i=0,min(100,npes-1)
          write(iulog,250)i,(procname((i))(j:j),j=1,length(i))
       end do
    endif

    deallocate (length, displ, procname)

100 format(//,i3," pes participating in computation for CLM")
200 format(/,35('-'))
220 format(/,"NODE#",2x,"NAME")
250 format("(",i5,")",2x,100a1,//)

  end subroutine spmd_init

end module spmdMod
