module parallel_mod
  ! ---------------------------
  use shr_kind_mod,   only: r8=>shr_kind_r8
  ! ---------------------------
  use dimensions_mod, only : nmpi_per_node, nlev, qsize_d
  ! ---------------------------
  use spmd_utils,     only: MPI_STATUS_SIZE, MPI_MAX_ERROR_STRING, MPI_TAG_UB

  implicit none
  private

  integer,  public, parameter   :: ORDERED         = 1
  integer,  public, parameter   :: FAST            = 2
  integer,  public, parameter   :: BNDRY_TAG_BASE  = 0
  integer,  public, parameter   :: THREAD_TAG_BITS = 9
  integer,  public, parameter   :: MAX_ACTIVE_MSG = (MPI_TAG_UB/2**THREAD_TAG_BITS) - 1
  integer,  public, parameter   :: HME_status_size = MPI_STATUS_SIZE

  integer,  public, parameter   :: HME_BNDRY_P2P   = 1
  integer,  public, parameter   :: HME_BNDRY_MASHM = 2
  integer,  public, parameter   :: HME_BNDRY_A2A   = 3
  integer,  public, parameter   :: HME_BNDRY_A2AO  = 4

  integer,  public, parameter   :: nrepro_vars = MAX(10, nlev*qsize_d)

  integer,  public              :: MaxNumberFrames
  integer,  public              :: numframes
  integer,  public              :: useframes
  logical,  public              :: PartitionForNodes
  logical,  public              :: PartitionForFrames

  ! Namelist-selectable type of boundary comms (AUTO,P2P,A2A,MASHM)
  integer,  public              :: boundaryCommMethod

  integer,  public, allocatable :: status(:,:)
  integer,  public, allocatable :: Rrequest(:)
  integer,  public, allocatable :: Srequest(:)

  real(r8), public, allocatable :: FrameWeight(:)
  integer,  public, allocatable :: FrameIndex(:)
  integer,  public, allocatable :: FrameCount(:)
  integer,  public              :: nComPoints
  integer,  public              :: nPackPoints

  real(r8), public, allocatable :: global_shared_buf(:,:)
  real(r8), public              :: global_shared_sum(nrepro_vars)

  ! ==================================================
  ! Define type parallel_t for distributed memory info
  ! ==================================================
  type, public :: parallel_t
    integer :: rank                       ! local rank
    integer :: root                       ! local root
    integer :: nprocs                     ! number of processes in group
    integer :: comm                       ! communicator
    integer :: intracomm                  ! Intra-node communicator
    integer :: commGraphFull              ! distributed graph topo communicator for all neighbors
    integer :: commGraphInter             ! distributed graph topo communicator for off-node neighbors
    integer :: commGraphIntra             ! distributed graph topo communicator for on-node neighbors
    integer :: groupGraphFull
    logical :: masterproc
  end type

  type (parallel_t), public :: par ! info for distributed memory programming

  ! ===================================================
  ! Module Interfaces
  ! ===================================================

  public :: initmpi
  public :: syncmp
  public :: copy_par

  interface assignment ( = )
    module procedure copy_par
  end interface

CONTAINS

! ================================================
!   copy_par: copy constructor for parallel_t type
!
!
!   Overload assignment operator for parallel_t
! ================================================

  subroutine copy_par(par2,par1)
    type(parallel_t), intent(out) :: par2
    type(parallel_t), intent(in)  :: par1

    par2%rank       = par1%rank
    par2%root       = par1%root
    par2%nprocs     = par1%nprocs
    par2%comm       = par1%comm
    par2%intracomm  = par1%intracomm
    par2%commGraphFull   = par1%commGraphFull
    par2%commGraphInter  = par1%commGraphInter
    par2%commGraphIntra  = par1%commGraphIntra
    par2%groupGraphFull  = par1%groupGraphFull
    par2%masterproc = par1%masterproc

  end subroutine copy_par

! ================================================
!  initmpi:
!  Initializes the parallel (message passing)
!  environment, returns a parallel_t structure..
! ================================================

  function initmpi(npes_homme) result(par)
    use cam_logfile,    only: iulog
    use cam_abortutils, only: endrun
    use spmd_utils,     only: mpicom, MPI_COMM_NULL, MPI_MAX_PROCESSOR_NAME
    use spmd_utils,     only: MPI_CHARACTER, MPI_INTEGER, MPI_BAND, iam, npes

    integer, intent(in) :: npes_homme

    type(parallel_t)     :: par

    integer              :: ierr,tmp
    integer              :: FrameNumber
    logical :: running   ! state of MPI at beginning of initmpi call
    character(len=MPI_MAX_PROCESSOR_NAME)               :: my_name
    character(len=MPI_MAX_PROCESSOR_NAME), allocatable  :: the_names(:)

    integer, allocatable :: tarray(:)
    integer              :: namelen, i
    integer              :: color

    !================================================
    !     Basic MPI initialization
    ! ================================================

    call MPI_initialized(running, ierr)

    if (.not.running) then
      call endrun('initmpi: MPI not initialized for SE dycore')
    end if

    par%root          = 0
    par%masterproc    = .FALSE.
    nmpi_per_node     = 2
    PartitionForNodes = .TRUE.

    ! The SE dycore needs to split from CAM communicator for npes > par%nprocs
    color = iam / npes_homme
    call mpi_comm_split(mpicom, color, iam, par%comm, ierr)
    if (iam < npes_homme) then
      call MPI_comm_size(par%comm, par%nprocs, ierr)
      call MPI_comm_rank(par%comm, par%rank,  ierr)
      if ( par%nprocs /= npes_homme) then
        call endrun('INITMPI: SE communicator count mismatch')
      end if

      if(par%rank == par%root) then
        par%masterproc = .TRUE.
      end if
    else
      par%rank   = 0
      par%nprocs = 0
      par%comm   = MPI_COMM_NULL
    end if

    if (par%masterproc) then
      write(iulog, '(a,i0)')'initmpi: Number of MPI processes: ', par%nprocs
    end if

    if (iam < npes_homme) then
      ! ================================================
      !  Determine where this MPI process is running
      !   then use this information to determined the
      !   number of MPI processes per node
      ! ================================================
      my_name(:) = ''
      call MPI_Get_Processor_Name(my_name, namelen, ierr)

      allocate(the_names(par%nprocs))
      do i = 1, par%nprocs
        the_names(i)(:) =  ''
      end do

      ! ================================================
      !   Collect all the machine names
      ! ================================================
      call MPI_Allgather(my_name, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
           the_names,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,par%comm,ierr)

      ! ======================================================================
      !   Calculate how many other MPI processes are on my node
      ! ======================================================================
      nmpi_per_node = 0
      do i = 1, par%nprocs
        if(TRIM(ADJUSTL(my_name)) .eq. TRIM(ADJUSTL(the_names(i)))) then
          nmpi_per_node = nmpi_per_node + 1
        end if
      end do

      ! =======================================================================
      !  Verify that everybody agrees on this number otherwise do not do
      !  the multi-level partitioning
      ! =======================================================================
      call MPI_Allreduce(nmpi_per_node,tmp,1,MPI_INTEGER,MPI_BAND,par%comm,ierr)
      if(tmp /= nmpi_per_node) then
        if (par%masterproc) then
          write(iulog,*)'initmpi:  disagrement accross nodes for nmpi_per_node'
        end if
        nmpi_per_node = 1
        PartitionForNodes = .FALSE.
      else
        PartitionForNodes = .TRUE.
      end if

      if(PartitionForFrames .and. par%masterproc) then
        write(iulog,*)'initmpi: FrameWeight: ', FrameWeight
      end if

      deallocate(the_names)
    end if

  end function initmpi

  ! =====================================
  ! syncmp:
  !
  ! sychronize message passing domains
  !
  ! =====================================
  subroutine syncmp(par)
    use cam_abortutils, only: endrun
    use spmd_utils,     only: MPI_MAX_ERROR_STRING, MPI_ERROR

    type (parallel_t), intent(in)       :: par

    integer                             :: errorcode, errorlen, ierr
    character(len=MPI_MAX_ERROR_STRING) :: errorstring

    call MPI_barrier(par%comm, ierr)

    if(ierr == MPI_ERROR) then
      errorcode = ierr
      call MPI_Error_String(errorcode, errorstring, errorlen, ierr)
      call endrun(errorstring)
    end if
  end subroutine syncmp

end module parallel_mod
