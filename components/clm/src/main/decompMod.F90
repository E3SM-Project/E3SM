module decompMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  ! Must use shr_sys_abort rather than endrun here to avoid circular dependency
  use shr_sys_mod , only : shr_sys_abort 
  use clm_varctl  , only : iulog
  use clm_varcon  , only : grlnd, nameg, namel, namec, namep, nameCohort
  use mct_mod     , only : mct_gsMap
  !
  ! !PUBLIC TYPES:
  implicit none
  integer, public :: clump_pproc ! number of clumps per MPI process

  ! Define possible bounds levels
  integer, parameter, public :: BOUNDS_LEVEL_PROC  = 1
  integer, parameter, public :: BOUNDS_LEVEL_CLUMP = 2
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public get_proc_clumps    ! number of clumps for this processor
  public get_proc_total     ! total no. of gridcells, landunits, columns and pfts for any processor
  public get_proc_global    ! total gridcells, landunits, columns, pfts across all processors
  public get_clmlevel_gsize ! get global size associated with clmlevel
  public get_clmlevel_gsmap ! get gsmap associated with clmlevel

  interface get_clump_bounds
     module procedure get_clump_bounds_old
     module procedure get_clump_bounds_new
  end interface
  public get_clump_bounds   ! clump beg and end gridcell,landunit,column,pft

  interface get_proc_bounds
     module procedure get_proc_bounds_old
     module procedure get_proc_bounds_new
  end interface
  public get_proc_bounds    ! this processor beg and end gridcell,landunit,column,pft

  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! !PRIVATE TYPES:
  private  ! (now mostly public for decompinitmod)

  integer,public :: nclumps     ! total number of clumps across all processors
  integer,public :: numg        ! total number of gridcells on all procs
  integer,public :: numt        ! total number of topographic units on all procs
  integer,public :: numl        ! total number of landunits on all procs
  integer,public :: numc        ! total number of columns on all procs
  integer,public :: nump        ! total number of pfts on all procs
  integer,public :: numCohort   ! total number of ED cohorts on all procs

  type bounds_type
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begt, endt       ! beginning and ending topographic unit index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending pft index
     integer :: begCohort, endCohort ! beginning and ending cohort indices

     integer :: level            ! whether defined on the proc or clump level
     integer :: clump_index      ! if defined on the clump level, this gives the clump index
  end type bounds_type
  public bounds_type

  !---global information on each pe
  type processor_type
     integer :: nclumps          ! number of clumps for processor_type iam
     integer,pointer :: cid(:)   ! clump indices
     integer :: ncells           ! number of gridcells in proc
     integer :: ntopounits       ! number of topographic units in proc
     integer :: nlunits          ! number of landunits in proc
     integer :: ncols            ! number of columns in proc
     integer :: npfts            ! number of pfts in proc
     integer :: nCohorts          ! number of cohorts in proc
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begt, endt       ! beginning and ending topographic unit index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending pft index
     integer :: begCohort, endCohort ! beginning and ending cohort indices
  end type processor_type
  public processor_type
  type(processor_type),public :: procinfo

  !---global information on each pe
  type clump_type
     integer :: owner            ! process id owning clump
     integer :: ncells           ! number of gridcells in clump
     integer :: ntopounits       ! number of topographic units in clump
     integer :: nlunits          ! number of landunits in clump
     integer :: ncols            ! number of columns in clump
     integer :: npfts            ! number of pfts in clump
     integer :: nCohorts          ! number of cohorts in clump
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begt, endt       ! beginning and ending topographic unit index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending pft index
     integer :: begCohort, endCohort ! beginning and ending cohort indices
  end type clump_type
  public clump_type
  type(clump_type),public, allocatable :: clumps(:)

  !---global information on each pe
  !--- glo = 1d global sn ordered
  !--- gdc = 1d global dc ordered compressed
  type decomp_type
     integer,pointer :: gdc2glo(:)    ! 1d gdc to 1d glo
  end type decomp_type
  public decomp_type
  type(decomp_type),public,target :: ldecomp

  type(mct_gsMap)  ,public,target :: gsMap_lnd_gdc2glo
  type(mct_gsMap)  ,public,target :: gsMap_gce_gdc2glo
  type(mct_gsMap)  ,public,target :: gsMap_lun_gdc2glo
  type(mct_gsMap)  ,public,target :: gsMap_col_gdc2glo
  type(mct_gsMap)  ,public,target :: gsMap_patch_gdc2glo
  type(mct_gsMap)  ,public,target :: gsMap_cohort_gdc2glo
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
   subroutine get_clump_bounds_new (n, bounds)
     !
     ! !DESCRIPTION:
     ! Determine clump bounds
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: n                ! processor clump index
     type(bounds_type), intent(out) :: bounds ! clump bounds
     !
     ! !LOCAL VARIABLES:
     character(len=32), parameter :: subname = 'get_clump_bounds'  ! Subroutine name
     integer :: cid                                                ! clump id
#ifdef _OPENMP
     integer, external :: OMP_GET_MAX_THREADS
     integer, external :: OMP_GET_NUM_THREADS
     integer, external :: OMP_GET_THREAD_NUM
#endif
     !------------------------------------------------------------------------------
     !    Make sure this IS being called from a threaded region
#ifdef _OPENMP
     ! FIX(SPM, 090314) - for debugging ED and openMP
     !write(iulog,*) 'SPM omp debug decompMod 1 ', &
          !OMP_GET_NUM_THREADS(),OMP_GET_MAX_THREADS(),OMP_GET_THREAD_NUM()

     if ( OMP_GET_NUM_THREADS() == 1 .and. OMP_GET_MAX_THREADS() > 1 )then
        call shr_sys_abort( trim(subname)//' ERROR: Calling from inside a non-threaded region)')
     end if
#endif

     cid  = procinfo%cid(n)
     bounds%begp = clumps(cid)%begp
     bounds%endp = clumps(cid)%endp
     bounds%begc = clumps(cid)%begc
     bounds%endc = clumps(cid)%endc
     bounds%begl = clumps(cid)%begl
     bounds%endl = clumps(cid)%endl
     bounds%begt = clumps(cid)%begt
     bounds%endt = clumps(cid)%endt
     bounds%begg = clumps(cid)%begg
     bounds%endg = clumps(cid)%endg
     bounds%begCohort = clumps(cid)%begCohort
     bounds%endCohort = clumps(cid)%endCohort
     
     bounds%level = BOUNDS_LEVEL_CLUMP
     bounds%clump_index = n

   end subroutine get_clump_bounds_new

   !------------------------------------------------------------------------------
   subroutine get_clump_bounds_old (n, begg, endg, begl, endl, begc, endc, begp, endp, &
        begCohort, endCohort)
     integer, intent(in)  :: n           ! proc clump index
     integer, intent(out) :: begp, endp  ! clump beg and end pft indices
     integer, intent(out) :: begc, endc  ! clump beg and end column indices
     integer, intent(out) :: begl, endl  ! clump beg and end landunit indices
     integer, intent(out) :: begg, endg  ! clump beg and end gridcell indices
     integer, intent(out) :: begCohort, endCohort  ! cohort beg and end gridcell indices
     integer :: cid                                                ! clump id
     !------------------------------------------------------------------------------

     cid  = procinfo%cid(n)
     begp = clumps(cid)%begp
     endp = clumps(cid)%endp
     begc = clumps(cid)%begc
     endc = clumps(cid)%endc
     begl = clumps(cid)%begl
     endl = clumps(cid)%endl
     begg = clumps(cid)%begg
     endg = clumps(cid)%endg
     begCohort = clumps(cid)%begCohort
     endCohort = clumps(cid)%endCohort
   end subroutine get_clump_bounds_old

   !------------------------------------------------------------------------------
   subroutine get_proc_bounds_new (bounds)
     !
     ! !DESCRIPTION:
     ! Retrieve processor bounds
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(out) :: bounds ! processor bounds bounds
     !
     ! !LOCAL VARIABLES:
#ifdef _OPENMP
     integer, external :: OMP_GET_NUM_THREADS
     integer, external :: OMP_GET_MAX_THREADS
     integer, external :: OMP_GET_THREAD_NUM
#endif
     character(len=32), parameter :: subname = 'get_proc_bounds'  ! Subroutine name
     !------------------------------------------------------------------------------
     !    Make sure this is NOT being called from a threaded region
#ifdef _OPENMP
     ! FIX(SPM, 090314) - for debugging ED and openMP
     !write(*,*) 'SPM omp debug decompMod 2 ', &
          !OMP_GET_NUM_THREADS(),OMP_GET_MAX_THREADS(),OMP_GET_THREAD_NUM()

     if ( OMP_GET_NUM_THREADS() > 1 )then
        call shr_sys_abort( trim(subname)//' ERROR: Calling from inside  a threaded region')
     end if
#endif

     bounds%begp = procinfo%begp
     bounds%endp = procinfo%endp
     bounds%begc = procinfo%begc
     bounds%endc = procinfo%endc
     bounds%begl = procinfo%begl
     bounds%endl = procinfo%endl
     bounds%begt = procinfo%begt
     bounds%endt = procinfo%endt
     bounds%begg = procinfo%begg
     bounds%endg = procinfo%endg
     bounds%begCohort = procinfo%begCohort
     bounds%endCohort = procinfo%endCohort

     bounds%level = BOUNDS_LEVEL_PROC
     bounds%clump_index = -1           ! irrelevant for proc, so assigned a bogus value

   end subroutine get_proc_bounds_new

   !------------------------------------------------------------------------------
   subroutine get_proc_bounds_old (begg, endg, begt, endt, begl, endl, begc, endc, begp, endp, &
        begCohort, endCohort)

     integer, optional, intent(out) :: begp, endp  ! proc beg and end pft indices
     integer, optional, intent(out) :: begc, endc  ! proc beg and end column indices
     integer, optional, intent(out) :: begl, endl  ! proc beg and end landunit indices
     integer, optional, intent(out) :: begt, endt  ! proc beg and end topographic unit indices
     integer, optional, intent(out) :: begg, endg  ! proc beg and end gridcell indices
     integer, optional, intent(out) :: begCohort, endCohort  ! cohort beg and end gridcell indices
     !------------------------------------------------------------------------------

     if (present(begp)) begp = procinfo%begp
     if (present(endp)) endp = procinfo%endp
     if (present(begc)) begc = procinfo%begc
     if (present(endc)) endc = procinfo%endc
     if (present(begl)) begl = procinfo%begl
     if (present(endl)) endl = procinfo%endl
     if (present(begt)) begt = procinfo%begt
     if (present(endt)) endt = procinfo%endt
     if (present(begg)) begg = procinfo%begg
     if (present(endg)) endg = procinfo%endg
     if (present(begCohort)) begCohort = procinfo%begCohort
     if (present(endCohort)) endCohort = procinfo%endCohort
   end subroutine get_proc_bounds_old

   !------------------------------------------------------------------------------
   subroutine get_proc_total(pid, ncells, ntopounits, nlunits, ncols, npfts, nCohorts)
     !
     ! !DESCRIPTION:
     ! Count up gridcells, topographic units, landunits, columns, and pfts on process.
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: pid     ! proc id
     integer, intent(out) :: ncells  ! total number of gridcells on the processor
     integer, intent(out) :: ntopounits  ! total number of topographic units on the processor
     integer, intent(out) :: nlunits ! total number of landunits on the processor
     integer, intent(out) :: ncols   ! total number of columns on the processor
     integer, intent(out) :: npfts   ! total number of pfts on the processor
     integer, intent(out) :: nCohorts ! total number of cohorts on the processor
     !
     ! !LOCAL VARIABLES:
     integer :: cid       ! clump index
     !------------------------------------------------------------------------------

     nCohorts = 0
     npfts   = 0
     ncols   = 0
     nlunits = 0
     ntopounits = 0
     ncells  = 0
     do cid = 1,nclumps
        if (clumps(cid)%owner == pid) then
           ncells  = ncells  + clumps(cid)%ncells
           ntopounits = ntopounits + clumps(cid)%ntopounits
           nlunits = nlunits + clumps(cid)%nlunits
           ncols   = ncols   + clumps(cid)%ncols
           npfts   = npfts   + clumps(cid)%npfts
           nCohorts = nCohorts + clumps(cid)%nCohorts
        end if
     end do
   end subroutine get_proc_total

   !------------------------------------------------------------------------------
   subroutine get_proc_global(ng, nt, nl, nc, np, nCohorts)
     !
     ! !DESCRIPTION:
     ! Return number of gridcells, landunits, columns, and pfts across all processes.
     !
     ! !ARGUMENTS:
     integer, optional, intent(out) :: ng        ! total number of gridcells across all processors
     integer, optional, intent(out) :: nt        ! total number of gridcells across all processors
     integer, optional, intent(out) :: nl        ! total number of landunits across all processors
     integer, optional, intent(out) :: nc        ! total number of columns across all processors
     integer, optional, intent(out) :: np        ! total number of pfts across all processors
     integer, optional, intent(out) :: nCohorts  ! total number ED cohorts
     !------------------------------------------------------------------------------

     if (present(np)) np             = nump
     if (present(nc)) nc             = numc
     if (present(nl)) nl             = numl
     if (present(nt)) nt             = numt
     if (present(ng)) ng             = numg
     if (present(nCohorts)) nCohorts = numCohort

   end subroutine get_proc_global

   !------------------------------------------------------------------------------
   integer function get_proc_clumps()
     !
     ! !DESCRIPTION:
     ! Return the number of clumps.
     !------------------------------------------------------------------------------

     get_proc_clumps = procinfo%nclumps

   end function get_proc_clumps

   !-----------------------------------------------------------------------
   integer function get_clmlevel_gsize (clmlevel)
     !
     ! !DESCRIPTION:
     ! Determine 1d size from clmlevel
     !
     ! !USES:
     use domainMod , only : ldomain
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: clmlevel    !type of clm 1d array
     !-----------------------------------------------------------------------

     select case (clmlevel)
     case(grlnd)
        get_clmlevel_gsize = ldomain%ns
     case(nameg)
        get_clmlevel_gsize = numg
     case(namel)
        get_clmlevel_gsize = numl
     case(namec)
        get_clmlevel_gsize = numc
     case(namep)
        get_clmlevel_gsize = nump
     case(nameCohort)
        get_clmlevel_gsize = numCohort
     case default
        write(iulog,*) 'get_clmlevel_gsize does not match clmlevel type: ', trim(clmlevel)
        call shr_sys_abort()
     end select

   end function get_clmlevel_gsize

   !-----------------------------------------------------------------------
   subroutine get_clmlevel_gsmap (clmlevel, gsmap)
     !
     ! !DESCRIPTION:
     ! Compute arguments for gatherv, scatterv for vectors
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: clmlevel     ! type of input data
     type(mct_gsmap) , pointer    :: gsmap
     !----------------------------------------------------------------------

    select case (clmlevel)
    case(grlnd)
       gsmap => gsmap_lnd_gdc2glo
    case(nameg)
       gsmap => gsmap_gce_gdc2glo
    case(namel)
       gsmap => gsmap_lun_gdc2glo
    case(namec)
       gsmap => gsmap_col_gdc2glo
    case(namep)
       gsmap => gsmap_patch_gdc2glo
    case(nameCohort)
       gsmap => gsMap_cohort_gdc2glo
    case default
       write(iulog,*) 'get_clmlevel_gsmap: Invalid expansion character: ',trim(clmlevel)
       call shr_sys_abort()
    end select

  end subroutine get_clmlevel_gsmap

end module decompMod
