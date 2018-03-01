module decompMod

  ! !PUBLIC TYPES:
  implicit none

  integer,public :: nclumps     ! total number of clumps across all processors
  integer,public :: numg        ! total number of gridcells on all procs
  integer,public :: numl        ! total number of landunits on all procs
  integer,public :: numc        ! total number of columns on all procs
  integer,public :: nump        ! total number of pfts on all procs
  integer,public :: numCohort   ! total number of ED cohorts on all procs

  ! Define possible bounds levels
  integer, parameter, public :: BOUNDS_LEVEL_PROC  = 1
  integer, parameter, public :: BOUNDS_LEVEL_CLUMP = 2

  type, public :: bounds_type
    integer :: begg, endg       ! beginning and ending gridcell index
    integer :: begl, endl       ! beginning and ending landunit index
    integer :: begc, endc       ! beginning and ending column index
    integer :: begp, endp       ! beginning and ending pft index
    integer :: lbj, ubj
    integer :: begCohort, endCohort ! beginning and ending cohort indices

    integer :: level            ! whether defined on the proc or clump level
    integer :: clump_index      ! if defined on the clump level, this gives the clump index
  end type bounds_type

  !---global information on each pe
  type processor_type
     integer :: nclumps          ! number of clumps for processor_type iam
     integer,pointer :: cid(:)   ! clump indices
     integer :: ncells           ! number of gridcells in proc
     integer :: nlunits          ! number of landunits in proc
     integer :: ncols            ! number of columns in proc
     integer :: npfts            ! number of pfts in proc
     integer :: nCohorts          ! number of cohorts in proc
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending pft index
     integer :: begCohort, endCohort ! beginning and ending cohort indices
  end type processor_type
  public processor_type
  type(processor_type),public :: procinfo

  interface get_proc_bounds
     module procedure get_proc_bounds_old
     module procedure get_proc_bounds_new
  end interface
  public get_proc_bounds    ! this processor beg and end gridcell,landunit,column,pft
  public get_proc_global    ! total gridcells, landunits, columns, pfts across all processors

  contains

  !------------------------------------------------------------------------------
  subroutine get_proc_bounds_new (bounds)
    !
    ! !DESCRIPTION:
    ! Retrieve processor bounds
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(out) :: bounds ! processor bounds bounds

    bounds%begp = procinfo%begp
    bounds%endp = procinfo%endp
    bounds%begc = procinfo%begc
    bounds%endc = procinfo%endc
    bounds%begl = procinfo%begl
    bounds%endl = procinfo%endl
    bounds%begg = procinfo%begg
    bounds%endg = procinfo%endg
    bounds%begCohort = procinfo%begCohort
    bounds%endCohort = procinfo%endCohort

    bounds%level = BOUNDS_LEVEL_PROC
    bounds%clump_index = -1           ! irrelevant for proc, so assigned a bogus value

  end subroutine get_proc_bounds_new


  !------------------------------------------------------------------------------
  subroutine get_proc_bounds_old (begg, endg, begl, endl, begc, endc, begp, endp, &
       begCohort, endCohort)

    integer, optional, intent(out) :: begp, endp  ! proc beg and end pft indices
    integer, optional, intent(out) :: begc, endc  ! proc beg and end column indices
    integer, optional, intent(out) :: begl, endl  ! proc beg and end landunit indices
    integer, optional, intent(out) :: begg, endg  ! proc beg and end gridcell indices
    integer, optional, intent(out) :: begCohort, endCohort  ! cohort beg and end gridcell indices
    !------------------------------------------------------------------------------

    if (present(begp)) begp = procinfo%begp
    if (present(endp)) endp = procinfo%endp
    if (present(begc)) begc = procinfo%begc
    if (present(endc)) endc = procinfo%endc
    if (present(begl)) begl = procinfo%begl
    if (present(endl)) endl = procinfo%endl
    if (present(begg)) begg = procinfo%begg
    if (present(endg)) endg = procinfo%endg
    if (present(begCohort)) begCohort = procinfo%begCohort
    if (present(endCohort)) endCohort = procinfo%endCohort
  end subroutine get_proc_bounds_old


  !------------------------------------------------------------------------------
  subroutine get_proc_global(ng, nl, nc, np, nCohorts)
    !
    ! !DESCRIPTION:
    ! Return number of gridcells, landunits, columns, and pfts across all processes.
    !
    ! !ARGUMENTS:
    integer, optional, intent(out) :: ng        ! total number of gridcells across all processors
    integer, optional, intent(out) :: nl        ! total number of landunits across all processors
    integer, optional, intent(out) :: nc        ! total number of columns across all processors
    integer, optional, intent(out) :: np        ! total number of pfts across all processors
    integer, optional, intent(out) :: nCohorts  ! total number ED cohorts
    !------------------------------------------------------------------------------

    if (present(np)) np             = nump
    if (present(nc)) nc             = numc
    if (present(nl)) nl             = numl
    if (present(ng)) ng             = numg
    if (present(nCohorts)) nCohorts = numCohort

  end subroutine get_proc_global
end module decompMod
