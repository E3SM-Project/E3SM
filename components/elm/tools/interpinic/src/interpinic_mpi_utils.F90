module interpinic_mpi_utils

  !-----------------------------------------------------------------------
  ! MPI plumbing for interpinic: rank state, collective wrappers and the
  ! output-mesh decomposition.  This is a support module, NOT an MPI variant
  ! of interpinic.F90 -- both files are compiled and linked in every build,
  ! serial and parallel alike, and interpinic.F90 depends on this one.
  !
  ! The decomposition it supports, for the three nearest-neighbour searches,
  ! is: replicate the INPUT mesh on every rank, distribute
  ! the OUTPUT mesh, gather the resulting indices back to rank 0.  All netCDF
  ! I/O stays on rank 0.  Because each output point is independent of every
  ! other output point, and the inner loop order over the input mesh is
  ! untouched, results are bit-for-bit identical to a serial run.
  !
  ! Everything MPI is behind HAVE_MPI.  Without it this module still compiles
  ! and supplies serial equivalents -- iam=0, npes=1, masterproc=.true., and
  ! scatter/gather degenerate to copies -- so the calling code has exactly one
  ! control path and the serial build stays byte-identical to the old tool.
  !
  ! "use mpi" rather than "include 'mpif.h'" is deliberate.  gfortran 10 and
  ! later reject a file that passes two different types to the same external
  ! procedure, which is precisely what the typed Bcast/Scatterv/Gatherv
  ! wrappers below do.  The mpi module declares its buffer arguments
  ! NO_ARG_CHECK, so this compiles without needing -fallow-argument-mismatch
  ! in the Makefile.
  !
  ! counts/displs convention: both are dimension(npes), and element p+1
  ! describes rank p.  displs holds 0-based offsets, as MPI requires.
  !-----------------------------------------------------------------------

  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
#ifdef HAVE_MPI
  use mpi
#endif

  implicit none

  private

  public :: mpi_initialize
  public :: mpi_cleanup
  public :: interpinic_abort
  public :: compute_bounds
  public :: bcast_int
  public :: bcast_logical
  public :: bcast_int_1d
  public :: bcast_r8_1d
  public :: scatterv_int_1d
  public :: scatterv_r8_1d
  public :: gatherv_int_1d

  integer, public, save :: iam        = 0        ! this rank
  integer, public, save :: npes       = 1        ! number of ranks
  logical, public, save :: masterproc = .true.   ! iam == 0

#ifdef HAVE_MPI
  integer, save :: mpicom = MPI_COMM_NULL
#endif

contains

  !=======================================================================

  subroutine mpi_initialize()

    ! Start MPI and record the rank layout.  Plain MPI_Init is sufficient:
    ! the tool is MPI-only, with no threading, so no MPI_THREAD_* level has
    ! to be negotiated.

    implicit none
#ifdef HAVE_MPI
    integer :: ier

    call MPI_Init(ier)
    if (ier /= MPI_SUCCESS) then
       write(6,*) 'interpinic: MPI_Init failed'
       call shr_sys_flush(6)
       stop 1
    end if

    mpicom = MPI_COMM_WORLD
    call MPI_Comm_rank(mpicom, iam,  ier)
    call MPI_Comm_size(mpicom, npes, ier)
#endif

    masterproc = (iam == 0)

    ! Say plainly what this binary is.  Without it there is no way to tell from
    ! a log whether a run is using every rank it was given, or whether a serial
    ! binary was launched under srun (in which case every rank redundantly
    ! computes the whole problem and they all write the same output file).
    if (masterproc) then
#ifdef HAVE_MPI
       write(6,'(a,i0,a)') ' interpinic: MPI build, running on ', npes, ' rank(s)'
#else
       write(6,'(a)') ' interpinic: SERIAL build (compiled without -DHAVE_MPI)'
#endif
       call shr_sys_flush(6)
    end if

  end subroutine mpi_initialize

  !=======================================================================

  subroutine mpi_cleanup()

    implicit none
#ifdef HAVE_MPI
    integer :: ier

    call MPI_Finalize(ier)
#endif

  end subroutine mpi_cleanup

  !=======================================================================

  subroutine interpinic_abort(msg)

    ! Collective-safe abort.  A bare "stop" on one rank would leave the
    ! others blocked in the following Gatherv, so every error path reachable
    ! from inside a search must come through here.

    implicit none
    character(len=*), intent(in) :: msg
#ifdef HAVE_MPI
    integer :: ier
#endif

    write(6,'(a,i0,a,a)') 'interpinic: ABORT (rank ', iam, '): ', trim(msg)
    call shr_sys_flush(6)

#ifdef HAVE_MPI
    call MPI_Abort(mpicom, 1, ier)
#endif

    stop 1

  end subroutine interpinic_abort

  !=======================================================================

  subroutine compute_bounds(ntot, active, counts, displs)

    ! Split 1..ntot into npes contiguous, ascending blocks holding roughly
    ! equal numbers of ACTIVE points.
    !
    ! Balancing on active count rather than on raw index count matters
    ! because the searches skip output points with zero weight, and restart
    ! files do cluster those.  Blocks stay contiguous and ascending so that
    ! Scatterv/Gatherv remain one-liners.
    !
    ! "active" is read on rank 0 only; the resulting layout is broadcast, so
    ! every rank leaves with the same counts/displs.

    implicit none
    integer, intent(in)  :: ntot          ! total number of output points
    logical, intent(in)  :: active(:)     ! rank 0: .true. where work is needed
    integer, intent(out) :: counts(:)     ! dimension(npes)
    integer, intent(out) :: displs(:)     ! dimension(npes), 0-based

    integer :: n, p, nact, acc, tgt, nmin, nmax

    if (masterproc) then

       counts(:) = 0
       displs(:) = ntot          ! ranks that get no block sit at the end

       nact = 0
       do n = 1, ntot
          if (active(n)) nact = nact + 1
       end do

       p         = 0
       acc       = 0
       displs(1) = 0

       do n = 1, ntot
          if (active(n)) acc = acc + 1
          if (p < npes-1 .and. n < ntot) then
             ! Cumulative target for the end of rank p, computed from the
             ! ideal split rather than incrementally, so rounding does not
             ! accumulate across ranks.
             tgt = int((int(p+1,8) * int(nact,8)) / int(npes,8))
             if (acc >= tgt) then
                counts(p+1) = n - displs(p+1)
                p           = p + 1
                displs(p+1) = n
             end if
          end if
       end do

       counts(p+1) = ntot - displs(p+1)

       ! Report the balance actually achieved.  counts is the block size; what
       ! costs time is the number of ACTIVE points in each block, so report
       ! that -- a large max/min spread here is what would justify revisiting
       ! this decomposition.
       nmin = huge(0)
       nmax = 0
       do p = 1, npes
          acc = 0
          do n = displs(p)+1, displs(p)+counts(p)
             if (active(n)) acc = acc + 1
          end do
          nmin = min(nmin, acc)
          nmax = max(nmax, acc)
       end do
       write(6,'(a,i0,a,i0,a,i0,a)') '   decomposition: ', ntot, ' points (', &
            nact, ' active) over ', npes, ' rank(s)'
       write(6,'(a,i0,a,i0)') '   active points per rank: min ', nmin, ' max ', nmax
       call shr_sys_flush(6)

    end if

    call bcast_int_1d(counts)
    call bcast_int_1d(displs)

  end subroutine compute_bounds

  !=======================================================================

  subroutine bcast_int(val)

    implicit none
    integer, intent(inout) :: val
#ifdef HAVE_MPI
    integer :: ier

    call MPI_Bcast(val, 1, MPI_INTEGER, 0, mpicom, ier)
#endif

  end subroutine bcast_int

  !=======================================================================

  subroutine bcast_logical(val)

    implicit none
    logical, intent(inout) :: val
#ifdef HAVE_MPI
    integer :: ier

    call MPI_Bcast(val, 1, MPI_LOGICAL, 0, mpicom, ier)
#endif

  end subroutine bcast_logical

  !=======================================================================

  subroutine bcast_int_1d(arr)

    implicit none
    integer, intent(inout) :: arr(:)
#ifdef HAVE_MPI
    integer :: ier

    if (size(arr) > 0) call MPI_Bcast(arr, size(arr), MPI_INTEGER, 0, mpicom, ier)
#endif

  end subroutine bcast_int_1d

  !=======================================================================

  subroutine bcast_r8_1d(arr)

    implicit none
    real(r8), intent(inout) :: arr(:)
#ifdef HAVE_MPI
    integer :: ier

    if (size(arr) > 0) call MPI_Bcast(arr, size(arr), MPI_REAL8, 0, mpicom, ier)
#endif

  end subroutine bcast_r8_1d

  !=======================================================================

  subroutine scatterv_int_1d(src, counts, displs, dst)

    ! src is meaningful on rank 0 only; dst receives counts(iam+1) elements.

    implicit none
    integer, intent(in)  :: src(:)
    integer, intent(in)  :: counts(:)
    integer, intent(in)  :: displs(:)
    integer, intent(out) :: dst(:)
#ifdef HAVE_MPI
    integer :: ier

    call MPI_Scatterv(src, counts, displs, MPI_INTEGER, &
                      dst, counts(iam+1), MPI_INTEGER, 0, mpicom, ier)
#else
    dst(1:counts(1)) = src(displs(1)+1 : displs(1)+counts(1))
#endif

  end subroutine scatterv_int_1d

  !=======================================================================

  subroutine scatterv_r8_1d(src, counts, displs, dst)

    implicit none
    real(r8), intent(in)  :: src(:)
    integer,  intent(in)  :: counts(:)
    integer,  intent(in)  :: displs(:)
    real(r8), intent(out) :: dst(:)
#ifdef HAVE_MPI
    integer :: ier

    call MPI_Scatterv(src, counts, displs, MPI_REAL8, &
                      dst, counts(iam+1), MPI_REAL8, 0, mpicom, ier)
#else
    dst(1:counts(1)) = src(displs(1)+1 : displs(1)+counts(1))
#endif

  end subroutine scatterv_r8_1d

  !=======================================================================

  subroutine gatherv_int_1d(src, counts, displs, dst)

    ! src holds this rank's counts(iam+1) results; dst is the full-length
    ! array, meaningful on rank 0 afterwards.

    implicit none
    integer, intent(in)    :: src(:)
    integer, intent(in)    :: counts(:)
    integer, intent(in)    :: displs(:)
    integer, intent(inout) :: dst(:)
#ifdef HAVE_MPI
    integer :: ier

    call MPI_Gatherv(src, counts(iam+1), MPI_INTEGER, &
                     dst, counts, displs, MPI_INTEGER, 0, mpicom, ier)
#else
    dst(displs(1)+1 : displs(1)+counts(1)) = src(1:counts(1))
#endif

  end subroutine gatherv_int_1d

  !=======================================================================

end module interpinic_mpi_utils
