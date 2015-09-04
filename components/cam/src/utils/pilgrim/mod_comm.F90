#define MOD_ASSUMED_SIZE
!#define MOD_SPECIFIED_SHAPE
#if !defined( STAND_ALONE )
#define NOR4 ! Activate to effectively eliminate real*4 window
#endif
!BOP
!
! !MODULE: mod_comm --- SPMD parallel decompostion/communication module
      module mod_comm
!
! !DESCRIPTION:
!
!  \paragraph{Overview}
!
!    This module contains SPMD parallelism decomposition and
!    communication routines.  This library was originally written by
!    W. Putman and S.-J. Lin for simple gridded communications in the
!    Finite-Volume General Circulation Model (FVGCM).  Most of the
!    member functions are specific to the type of gridded data, 
!    ghost communication and decompositions used in FVGCM (which 
!    are, however, very common in atmospheric models).
!
!    The module was extended for irregular communication
!    by W. Sawyer and A. Mirin. It is now
!    a more general tool and has been incorporated into the Parallel
!    Library for Grid Manipulations (PILGRIM) which is used in the
!    Community Atmospheric Model (CAM) and the Physical-space
!    Statistical Analysis System (PSAS).
!
!    **********************************************************
!    The storage associated with the irregular communications
!    is based on CAM requirements. It runs the risk of functioning
!    improperly when used in another code.
!    **********************************************************
!
!    Irregular communication is based on the {\tt blockdescriptor}
!    derived type, which defines a set of parcels which are to be
!    send to (or received from) another PE.  The irregular 
!    communication routines operate on arrays of block descriptors
!    whose length is equal to number of PEs involved in the
!    communication.  This means the irregular communication primitives
!    are merely non-blocking (potentially) all-to-all primitives.
! 
!    This package is based on standard MPI-1 communications, and OpenMP
!    may be implemented. MPI-2 and SHMEM support have been removed.
!
!
!  \paragraph{Use of Global Arrays}
!
!    The module uses the concept of global arrays (coined from former
!    usage of shared memory arenas in the "multi-level parallelism"
!    (MLP) paradigm).  Global arrays are merely buffers into which
!    data are packed for the transfer to other PEs and are not
!    necessarily of global extent. All such arrays are
!    1-dimensional; they are accessed as needed with offset vars.
!
!  \paragraph{Use of Windows}
!
!       All implementations use real*8, real*4, and integer*4 windows
!       which are used with global arrays as follows:
!
!       \begin{itemize}
!         \item   r8\_win -> ga\_r8 - for use with real*8 types
!         \item   r4\_win -> ga\_r4 - for use with real*4 types
!         \item   i4\_win -> ga\_i4 - for use with integer*4 types
!       \end{itemize}
!
!       note: MPI routines need 2 buffers per GA, ga\_<type>\_s & ga\_<type>\_r
!             ga\_<type>\_r is used for the windows
!
!  \paragraph{Compilation}
!
!    This module contains several precompile options:
!
!    \begin{itemize}
!      \item {\tt STAND_ALONE}:  Use as stand-alone library (if
!                                defined) or as part of CAM (if 
!                                undefined)
!      \item {\tt MODCM_TIMING}: Turn on CAM timing routines (only
!                                available if compiled in CAM framework)
!      \item {\tt _OPENMP}:      Implicit token (controlled by
!                                compiler) to enable OpenMP
!    \end{itemize}
!
!    
!  \paragraph{Usage}
!
!    NOTE - must call PILGRIM routine parinit to initialize before
!    making any other calls.
!
!    The public members of this module are:
!
!      \begin{itemize}
!         \item {\tt mp\_init}:          Initialize module
!         \item {\tt mp\_exit}:          Exit module
!         \item {\tt mp\_send4d\_ns}:    Ghost 4D array on north/south
!         \item {\tt mp\_recv4d\_ns}:    Complete 4D N/S ghost operation
!         \item {\tt mp\_send2\_ns}:     Ghost 2 3D arrays on north/south
!         \item {\tt mp\_recv2\_ns}:     Complete 2x3D N/S ghost operation
!         \item {\tt mp\_send3d}:        Send 3D general ghost region
!         \item {\tt mp\_recv3d}:        Complete 3D general ghost operation
!         \item {\tt mp\_send3d\_2}:     Send 2x3D general ghost regions
!         \item {\tt mp\_recv3d\_2}:     Complete 2x3D general ghost operation
!         \item {\tt get\_partneroffset}:Offset for remote write
!         \item {\tt mp\_sendirr}:       Initiate all-to-all send of parcels
!         \item {\tt mp\_recvirr}:       Complete all-to-all chunk commun.
!       \end{itemize}
!
!     There are variants of some of these routines for r4 and i4 data types.
!     There are other public routines, but these are only used internally
!     in PILGRIM, and they should not be called by user applications.
!
! !REVISION HISTORY:
!    2001.09.01   Lin
!    2002.04.16   Putman  Modified for Global Array code
!    2002.04.16   Putman  Added ProTeX documentation
!    2002.05.28   Putman  Added use of precision module
!    2003.06.24   Sawyer  Minor additions for use with mod_irreg
!    2004.01.08   Sawyer  Removed older functionality, no longer needed
!    2004.02.10   Mirin   Major restructuring and simplification. Documentation
!    2004.03.06   Sawyer  Additional documentation; cosmetics
!    2005.03.20   Sawyer  Added extensive support for real*4
!    2005.10.12   Worley  Improved vectorization of buffer copies and general clean-up
!    2006.05.15   Mirin   Make dynamic allocation the default; general clean-up.
! !USES:
#if defined( STAND_ALONE )
# define iulog 6
#else
      use cam_logfile, only: iulog
#endif

!
! Performance bug work around for Gemini interconnect
!
#ifdef _NO_MPI_RSEND
#define MPI_RSEND MPI_SEND
#define mpi_rsend mpi_send
#define MPI_IRSEND MPI_ISEND
#define mpi_irsend mpi_isend
#endif 

!
! Mod_comm has option for stand-alone use as well as within CAM
!

#if defined ( SPMD )

#if defined( STAND_ALONE )
# define r8 selected_real_kind(12)
# define r4 selected_real_kind( 6)
# define i8 selected_int_kind(13)
# define i4 selected_int_kind( 6)
# define PLON        144
# define PLAT         91
# define PLEV         26
# define PCNST         1
#else
      use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4,  &
                               i8 => shr_kind_i8, i4 => shr_kind_i4
#endif
#if defined( MODCM_TIMING )
      use perf_mod
#endif

      implicit none

#include "mpif.h"

! !PUBLIC MEMBER FUNCTIONS:
      public mp_init, mp_exit,                                             &
             mp_send4d_ns, mp_recv4d_ns, mp_send4d_ns_r4, mp_recv4d_ns_r4, &
             mp_send2_ns, mp_recv2_ns, mp_send3d_2, mp_recv3d_2,           &
             mp_send3d, mp_recv3d, mp_sendirr, mp_recvirr,                 &
             mp_sendirr_r4, mp_recvirr_r4, mp_sendirr_i4, mp_recvirr_i4,   &
             mp_swapirr, mp_swapirr_i4, mp_barrier,                        &
             get_partneroffset, mp_r8, mp_r4, mp_i4,                       &
             mp_sendtrirr, mp_recvtrirr, mp_swaptrirr
      public modcam_method, modcam_geopk, modcam_gatscat, modcam_npryz, modcam_maxirr

! !PRIVATE MEMBER FUNCTIONS:
      private ceil2    ! copy of routine in atm/cam/src/utils/spmdutils
      private pair     ! copy of routine in atm/cam/src/utils/spmdutils

!------------------------------------------------------------------------------
!  type declaration for describing an arbitrary number of contiguous parcels
!  this is for irregular communications
!------------------------------------------------------------------------------
      type blockdescriptor
         integer              :: method             ! transpose method
         integer              :: type               ! Ptr to MPI derived type
         integer, pointer     :: displacements(:)   ! Offsets in local segment
         integer, pointer     :: blocksizes(:)      ! Block sizes to transfer
         integer              :: partneroffset      ! Aggregated partner offset
         integer              :: partnertype        ! Ptr to partner's MPI derived type
         integer              :: Nparcels           ! size( displacements )
         integer              :: Tot_Size           ! sum ( blocksizes )
      end type blockdescriptor

! Transpose methods (method)
!      0 for contiguous temporary buffer
!      1 for direct communication (derived types)

! The variables immediately below refer specifically to mpi derived types
      INTEGER, ALLOCATABLE, SAVE :: InHandle(:, :)
      INTEGER, ALLOCATABLE, SAVE :: OutHandle(:, :)
      INTEGER, SAVE :: BegTrf = 0  ! Ongoing overlapped begintransfer #
      INTEGER, SAVE :: EndTrf = 0  ! Ongoing overlapped endtransfer #
      INTEGER, SAVE :: MaxTrf = 0  ! Max no. active Mp_sendirr derived type messages

! !PUBLIC DATA MEMBERS:
      integer, SAVE:: gid                         ! PE id
      integer(i4), SAVE:: masterpro = 0           ! Master process id 
      integer(i4), SAVE:: numpro                  ! Permanent No. of PEs
      integer(i4), SAVE:: numcomm                 ! Local No. of PEs
      integer(i4), SAVE:: numcpu                  ! No. of threads
      integer, SAVE:: commglobal                  ! Global Communicator
      integer, SAVE:: Max_Nparcels = 0            ! Maximum number of parcels in
                                                  !  single blockdescriptor

!------------------------------------------------------------------------------
!  Local parameters
!------------------------------------------------------------------------------
      integer, parameter:: nbuf = 2               ! Max No. of sends per call
! mp_send4d_ns has two sends per call (full border regions to north and south)
! mp_send2_ns has four sends per call (2 directions and 2 variables); however,
!  only one ghost latitude is sent, so nbuf=2 suffices as long as nghost
!  is greater than 1.
! mp_send3d has one send per call (border region in one direction).
! mp_send3d_2 has two sends per call (2 variables, border region in one direction).
      integer, parameter:: nghost = 3             ! No. of ghost indices
      integer, parameter:: max_nq = 1             ! No. of tracers simultaneously
                                                  !  border communicated; can be
                                                  !  overridden with dynamic storage
      integer, parameter:: max_trac = PCNST       ! No. of tracers
      integer, parameter:: max_call = 2           ! Max No. of back-to-back...
                                                  ! ...mp_send calls
! Currently, CAM has at most two overlapping border communication calls
! The above variable is relevant for contiguous irregular communications

      integer, parameter:: idimsize = PLON*nghost*(PLEV+1)*max_nq
                                                  ! Size of MPI buffer region
                                                  ! in mp_send/mp_recv calls, used
                                                  ! to determine offset in GA
      integer, parameter:: platg = PLAT + 2*nghost
      integer, parameter :: mp_r4 = MPI_REAL
      integer, parameter :: mp_r8 = MPI_DOUBLE_PRECISION
      integer, parameter :: mp_i4 = MPI_INTEGER

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------

      integer, SAVE:: max_irr = 0       ! Max No. active Mp_sendirr calls with window
      integer ierror
      integer, SAVE:: sizet1, sizer8, sizer4, sizei4

! CAM-specific variables
      integer, SAVE:: tracmax, tracbmax, dpvarmax, totvar
      integer, SAVE:: phys_transpose_mod
      integer, SAVE:: idimsizz
      integer, SAVE:: modcam_method, modcam_geopk, modcam_gatscat
      integer, SAVE:: modcam_npryz(4), modcam_tagoffset, modcam_maxirr
      integer, parameter :: phys_transpose_modmin = 11
      integer, parameter :: phys_transpose_vars = 7
      data phys_transpose_mod / -1 /
      data modcam_method / -1 /
      data modcam_geopk / -1 /
      data modcam_gatscat / -1 /
      data modcam_npryz / -1, -1, -1, -1 /
      data modcam_tagoffset / 0 /
      data modcam_maxirr / -1 /
!
! tracmax is the maximum number of tracers simultaneously transposed within dynamics (set to 1)
!    (except in dynamics-physics transposes)
! tracbmax is the maximum number of tracers simultaneously border communicated
! dpvarmax is the number of variables communicated in dynamics-physics transposes 
! totvar is the maximum number of variables simultaneously transposed
! phys_transpose_mod is the communication method for dynamics/physics transposes; admissable values
!     are >= phys_transpose_modmin; it is communicated from CAM when such transposes
!     are requested.
! phys_transpose_vars is the number of non-tracer variables transposed between dynamics and
!     physics instantiations in CAM.
! modcam_method, modcam_geopk and modcam_gatscat correspond to mod_method, mod_geopk and
!     mod_gatscat in CAM.
! modcam_npryz corresponds to npr_yz in CAM.
! modcam_maxirr corresonds to mod_maxirr in CAM.

!------------------------------------------------------------------------------
!  Variables to control global array locations and window synchronization
!------------------------------------------------------------------------------
      integer win_count                 ! Counts No. of windows in use
      integer igosouth, igonorth        ! Index of latitudinal send direction
      integer ifromsouth, ifromnorth    ! Index of latitudinal recv direction

!------------------------------------------------------------------------------
!  Local type declaration for mp_windows
!------------------------------------------------------------------------------
      type window
         integer :: id            ! Window id
         integer :: size          ! Size of global window (point based)
         integer :: ncall_s       ! Count send calls on window
         integer :: ncall_r       ! Count recv calls on window
         integer :: offset_s      ! Starting position in GA send
         integer :: offset_r      ! Starting position in GA recv
         integer :: dest          ! For use with send calls
         integer :: src           ! For use with recv calls
         integer :: size_r        ! Size of incoming message
         integer :: nsend         ! Send counter
         integer :: nrecv         ! Receive post counter
         integer :: nread         ! Receive confirm counter
         integer, pointer :: sqest(:) ! Send handle
         integer, pointer :: rqest(:) ! Receive handle
     end type window

!------------------------------------------------------------------------------
! Beginning Global Array variable declaration:
!------------------------------------------------------------------------------

      type (window) :: r8_win
      type (window) :: r4_win
      type (window) :: i4_win
      type (window) :: t1_win

! Upper bound on ratio of local to average storage over subdomains.
! This takes into account different sized subdomains.

      real*8, parameter :: alloc_slack_factor = 1.2_r8

!
!   window variable declarations
!
      real(r8), allocatable,    SAVE:: ga_t1_r(:)
      real(r8), allocatable,    SAVE:: ga_t1_s(:)
      real(r8), allocatable,    SAVE:: ga_r8_r(:)
      real(r8), allocatable,    SAVE:: ga_r8_s(:)
      real(r4), allocatable,    SAVE:: ga_r4_r(:)
      real(r4), allocatable,    SAVE:: ga_r4_s(:)
      integer(i4), allocatable, SAVE:: ga_i4_r(:)
      integer(i4), allocatable, SAVE:: ga_i4_s(:)
!
!   auxiliary variable declarations
!
      integer, SAVE:: Status(MPI_STATUS_SIZE)
      integer, allocatable, SAVE:: Stats(:)
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_init --- Initialize SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_init( comm, npryzxy, mod_method, mod_geopk, mod_gatscat, mod_maxirr )
!
! !INPUT PARAMETERS:
      integer, optional :: comm                        ! communicator
      integer, optional, intent(in) :: npryzxy(4)      ! 2D decomposition
      integer, optional, intent(in) :: mod_method      ! CAM optimization
      integer, optional, intent(in) :: mod_geopk       ! CAM optimization
      integer, optional, intent(in) :: mod_gatscat     ! CAM optimization
      integer, optional, intent(in) :: mod_maxirr      ! CAM optimization
! !DESCRIPTION:
!
!     Initialize SPMD parallel communication.  It is recommended that
!     COMM (main communicator) and NPRYZXY (2D decomposition) be set.
!
!     Set the mod* variables only if you are acquainted with their 
!     meaning (default is 0).
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman        Modified for Global Array code
!    2002.04.09   Putman        Added ProTeX documentation
!    2002.08.06   Sawyer        Added optional communicator input argument
!    2006.06.15   Sawyer        Added CAM-dependent optional arguments
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer mysize
      integer using_window, vertical_lines, latitude_lines
      integer local_dynamic_storage, numpro_use
      real*8 geopkrat, one, ghostrat

! Initialize MPI; allow for general communicator
      if ( present(comm) ) then
        call mpi_start( comm )
      else
        call mpi_start( MPI_COMM_WORLD )
      endif
! Initialize OpenMP
      call omp_start
!
! Adopt 2D decomposition if provided.
!
      modcam_npryz = (/ 1,1,1,1 /)    ! Default value (sequential)
      if ( present( npryzxy ) ) then
          modcam_npryz(1:4) = npryzxy(1:4)
          modcam_tagoffset = modcam_npryz(3) * modcam_npryz(4)
      endif
      if (gid .eq. 0) then
        write (iulog,*) 'MOD_COMM - modcam_npryz = ', &
               modcam_npryz(1), modcam_npryz(2),     &
               modcam_npryz(3), modcam_npryz(4)
        write (iulog,*) 'MOD_COMM - modcam_tagoffset = ', modcam_tagoffset
      endif

!
! Set CAM optimization variables
!
! modcam_method refers to irregular communications for transposes
! modcam_geopk refers to irregular communications for the geopotential
! modcam_gatscat refers to irregular communications for gather/scatters
! For any of these, a value of 0 means source data will be gathered into a contiguous
!  buffer (window), communicated to a contiguous buffer (window) in the target, and
!  then scattered to its final destination; a value of 1 means MPI derived types will
!  be used (hence not requiring window storage).
! modcam_maxirr refers to maximum number of irregular communications to be active at once
      modcam_method  = 0   ! Default value
      modcam_geopk   = 0   ! Default value
      modcam_gatscat = 0   ! Default value
      modcam_maxirr = 1   ! Default value
      if ( present( mod_method ) )  modcam_method  = mod_method
      if ( present( mod_geopk ) )   modcam_geopk   = mod_geopk
      if ( present( mod_gatscat ) ) modcam_gatscat = mod_gatscat
      if ( present( mod_maxirr ) ) modcam_maxirr = mod_maxirr

      if (gid .eq. 0) then
        write(iulog,*) 'MOD_COMM - modcam_method modcam_geopk modcam_gatscat modcam_maxirr = ',    &
        modcam_method, modcam_geopk, modcam_gatscat, modcam_maxirr
      endif

!
! End CAM optimizations
!

      MaxTrf = modcam_maxirr
      max_irr = modcam_maxirr

      win_count = 0
!
!*************************************************************************
! local_dynamic_storage is set to 1 when window storage is based on locally dimensioned
!  arrays, 0 otherwise; this occurs when modcam_gatscat equals 1, as it is only the
!  gather/scatters that require global storage.
!*************************************************************************
!
      local_dynamic_storage = 0
      if (modcam_gatscat .eq. 1) local_dynamic_storage = 1

!*************************************************************************
! Override original strategy, as only single 2D lat-lon variables (rather
!  than 3D with multiple tracers) are used for gather/scatters;
!  set local_dynamic_storage to 1 always, and then allow for gather
!  of 2D lat-lon variable in inidat.
!*************************************************************************

      local_dynamic_storage = 1

      allocate( Stats(MAX(nbuf,numpro)*MAX(max_call,max_irr)*MPI_STATUS_SIZE) )
      allocate( InHandle(numpro,MaxTrf) )
      allocate( OutHandle(numpro,MaxTrf) )

      idimsizz = idimsize
      if (local_dynamic_storage .eq. 1) then
         if (gid .eq. 0) write(iulog,*) 'Using local dynamic storage for mod_comm window'
      else
         if (gid .eq. 0) write(iulog,*) 'Using global dynamic storage for mod_comm window'
      endif
!
! Dynamically allocate target global arrays
!
!*************************************************************************
! Compute additional storage due to ghost latitudes being included in some
!   transposes. Allow 3 ghost points on each side. The required storage
!   could be (6+L)/L times the original storage, where L is the number of
!   latitude lines in the subdomain. Ghost points can also occur in the
!   vertical due to edge quantities, but this would not occur simultaneously
!   with ghost points in latitude; the extra storage due to vertical ghost
!   points is not nearly as great as with latitude.
!*************************************************************************
      using_window = 1   !  This is a local variable
      if (modcam_method .eq. 1) using_window = 0
      one = real(1,r8)
      ghostrat = one
      if (using_window .eq. 1 .and. local_dynamic_storage .eq. 1) then
         latitude_lines = real(PLAT,r8)/real(modcam_npryz(1),r8)
         ghostrat = real(6+latitude_lines,r8)/real(latitude_lines,r8)
      endif
      if (gid .eq. 0) write(iulog,*) 'Mod_comm - ghostrat = ', ghostrat

!*************************************************************************
! Compute extent to which required window storage for geopotential computation
!   exceeds that of transpose - relevant only for local dynamic storage,
!   since with global storage there will be enough space anyway; also,
!   this applies only when using window; further, this applies only when
!   the CAM variable geopktrans equals 1, though we do not test for that here.
! The geopotential calculation sends a latitude line to every other process
!   either vertically above or below the given process; there can be
!   at most modcam_npryz(2)-1 such target processes; compared to transposes
!   (which send all vertical lines), the amount of data sent is expressed
!   as the ratio geopkrat; our concern is making the window (whose size
!   is computed based on transposes) large enough, so we must multiply its
!   size by geopkrat; we never shrink the window, so geopkrat >= 1.
!*************************************************************************
      using_window = 1   !  This is a local variable
      if (modcam_geopk .eq. 1) using_window = 0
      one = real(1,r8)
      geopkrat = one
      if (using_window .eq. 1 .and. local_dynamic_storage .eq. 1) then
         vertical_lines = ceiling(real(PLEV,r8)/real(modcam_npryz(2),r8))
         geopkrat = real(modcam_npryz(2)-1,r8)/real(vertical_lines,r8)
         geopkrat = max(geopkrat,one)
      endif
      if (gid .eq. 0) write(iulog,*) 'Mod_comm - geopkrat = ', geopkrat

!*************************************************************************
! beginning of CAM totvar computation
!*************************************************************************

! CAM contains two kinds of transposes. The most commonly referred to transposes
!  refer to those which connect the xy and yz decompositions. Depending on
!  the physics decomposition, CAM might additionally compute transposes between
!  the dynamics and physics; this depends on the variable phys_loadbalance.
!  Furthermore, these transposes might or might not be computed using mod_comm.
!  The former transposes are generally performed one variable at a time; the
!  latter transposes combine all variables to be transposed, including the
!  full complement of tracers. The maximum number of variables to be 
!  simultaneously subject to irregular communications is dependent on
!  whether or not mod_comm is used to compute dynamics-physics transposes
!  and could depend on the number of tracers.

! Compute maximum number of variables to be simultaneously subject
!  to irregular communications (e.g., transposed variables based on CAM)
!  and store in the variable 'totvar'.

! Tracmax is the number of tracers simultaneously transposed within dynamics;
! Tracbmax is the number of tracers simultaneously border comunicated within trac2d;
!  both of these are currently hardwired to 1.
      tracmax = 1
      tracbmax = 1
      totvar = tracmax

! Now consider dynamics-physics transposes in CAM dp_coupling (dpvarmax)
!  If phys_transpose_mod is still -1, that means it has not been updated
!  by CAM and hence mod_comm will not be used for dynamics-physics transposes.
! (NOTE: phys_transpose_mod is computed in phys_grid_setopts in phys_grid.F90.)

! Also note that the logic involving phys_transpose_mod and phys_transpose_modmin
!  must remain consistent with the coding in phys_grid.F90. Additionally,
!  phys_transpose_vars must remain consistent with the coding in dp_coupling.F90.
!  (See above declaration and initialization for CAM-specific variables.)

! (begin dpvarmax calculation)

      if (phys_transpose_mod .eq. -1) then
         if (gid .eq. 0) write(iulog,*)       &
           '(MOD_COMM) - mod_comm not being used for dynamcis-physics transposes'
         dpvarmax = 0
!
! If phys_transpose_mod is >= phys_transpose_modmin, that is a signal that mod_comm is to be used
!  for dynamics/physics transposes in CAM. In that case, one must allocate enough window
!  storage for those transposes. Presently, the number of such simultaneously transposed
!  variables equals phys_transpose_vars plus the number of constituents.
!
      elseif (phys_transpose_mod .ge. phys_transpose_modmin) then
         dpvarmax = phys_transpose_vars + max_trac
      else
         dpvarmax = 0
      endif

! (end dpvarmax calculation)

! totvar is the maximum of (1) the number of tracers to be simultaneously transposed
!  within the dynamics, and (2) the number of variables to be transposed between
!  dynamics and physics instantiations in CAM

      totvar = max(totvar, dpvarmax)

!*************************************************************************
! end of CAM totvar computation
!*************************************************************************

      if (gid .eq. 0) write(iulog,*) 'Mod_comm - tracmax dpvarmax totvar tracbmax = ',     &
          tracmax, dpvarmax, totvar, tracbmax 

      idimsizz = (idimsize/max_nq)*tracbmax
      sizet1 = idimsizz*nbuf*max_call
! Adjust window sizes for geopotential and/or ghost points
      sizer8 = PLON*platg*(PLEV+1)*totvar*max(geopkrat,ghostrat)*max_irr
      sizer4 = PLON*platg*(PLEV+1)*totvar*max(geopkrat,ghostrat)*max_irr
      sizei4 = PLON*PLAT*PLEV*max_irr

! Compute local storage requirement for irregular communications by dividing
!    global requirement by the number of tasks. Allow slack factor to account
!    for nonuniformity of decomposition and ghost zones. Not valid for global
!    operations such as gathers and scatters when local windows are used.
      if (local_dynamic_storage .eq. 1) then
         numpro_use = modcam_npryz(1) * modcam_npryz(2)
         sizer8 = ceiling( alloc_slack_factor*real(sizer8,r8)/real(numpro_use,r8) )

! Allow for gather of single 2D lat-lon variable in inidat.
         if (modcam_gatscat .eq. 0) sizer8 = max( sizer8, PLON*PLAT*max_irr )  

         sizer4 = ceiling( alloc_slack_factor*real(sizer4,r8)/real(numpro_use,r8) )
! The only i4 irregular communications in CAM occur in io_dist.
         sizei4 = 1
      endif

# if defined ( NOR4 )
      sizer4 = 1
      if (gid .eq. 0) write(iulog,*) 'Mod_comm - r4 windows disabled'
# endif

      using_window = 1   !  This is a local variable
      if (modcam_method .eq. 1 .and. modcam_geopk .eq. 1) using_window = 0
      if (using_window .eq. 0) then
         if (gid .eq. 0) write(iulog,*) 'Mod_comm - r8 and r4 windows set to trivial size'
         sizer8 = 1
         sizer4 = 1
      endif

! Allocate global storage

      allocate( ga_t1_r(sizet1) )
      allocate( ga_t1_s(sizet1) )
      allocate( ga_r8_r(sizer8) )
      allocate( ga_r8_s(sizer8) )
      allocate( ga_r4_r(sizer4) )
      allocate( ga_r4_s(sizer4) )
      allocate( ga_i4_r(sizei4) )
      allocate( ga_i4_s(sizei4) )

! Initialize windows

        mysize = sizet1
        call win_init_r8(comm, t1_win, ga_t1_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm t1_win window size = ', mysize

        mysize = sizer8
        call win_init_r8(comm, r8_win, ga_r8_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm r8_win window size = ', mysize

        mysize = sizer4
        call win_init_r4(comm, r4_win, ga_r4_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm r4_win window size = ', mysize

        mysize = sizei4
        call win_init_i4(comm, i4_win, ga_i4_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm i4_win window size = ', mysize

        igosouth   = 0
        igonorth   = 1
        ifromsouth = 1
        ifromnorth = 0

!EOC
      end subroutine mp_init
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_exit --- End SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_exit( comm )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !DESCRIPTION:
!
!     End SPMD parallel communication
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman        Modified for Global Array code
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
        call MPI_FINALIZE (ierror)
        return
!EOC
      end subroutine mp_exit
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: omp_start --- Start openMP parallelism
!
! !INTERFACE:
      subroutine omp_start
! !DESCRIPTION:
!
!     Start openMP parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer ios, n, nowpro, nowcpu

! Compute number of OpenMP threads

#if defined(_OPENMP)

        integer omp_get_num_threads
!$omp parallel
        numcpu = omp_get_num_threads()
!$omp end parallel

#else
        numcpu = 1
#endif

!EOC
      end subroutine omp_start
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mpi_start --- Start MPI parallelism
!
! !INTERFACE:
      subroutine mpi_start( comm )
! !INPUT PARAMETERS:
      integer :: comm      !  communicator
! !DESCRIPTION:
!
!     Start MPI parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    02.08.06   Sawyer  Added communicator input arguments
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        logical flag
        integer npthreads

        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT( ierror )
          comm = MPI_COMM_WORLD
        endif

        call MPI_COMM_RANK (comm, gid, ierror)
        call MPI_COMM_SIZE (comm, numpro, ierror)
        call MPI_COMM_DUP  (comm, commglobal, ierror)
!EOC
      end subroutine mpi_start
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_r8 --- Initialize real*8 communication window
!
! !INTERFACE:
      subroutine win_init_r8(comm, win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: comm      !  communicator
        integer, intent(in) :: isize
        real(r8), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize real*8 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

        win_count = win_count + 1
        win%id = win_count
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
        win%nsend = 0
        win%nrecv = 0
        win%nread = 0
        allocate( win%sqest(MAX(nbuf,numpro)*MAX(max_call,max_irr)) )
        allocate( win%rqest(MAX(nbuf,numpro)*MAX(max_call,max_irr)) )
!EOC
      end subroutine win_init_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_r4 --- Initialize real*4 communication window
!
! !INTERFACE:
      subroutine win_init_r4(comm, win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: comm      !  communicator
        integer, intent(in) :: isize
        real(r4), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize real*4 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

        win_count = win_count + 1
        win%id = win_count
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
        win%nsend = 0
        win%nrecv = 0
        win%nread = 0
        allocate( win%sqest(MAX(nbuf,numpro)*MAX(max_call,max_irr)) )
        allocate( win%rqest(MAX(nbuf,numpro)*MAX(max_call,max_irr)) )
!EOC
      end subroutine win_init_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_i4 --- Initialize integer*4 communication window
!
! !INTERFACE:
      subroutine win_init_i4(comm, win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: comm      !  communicator
        integer, intent(in) :: isize
        integer(i4), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize integer*4 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

        win_count = win_count + 1
        win%id = win_count
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
        win%nsend = 0
        win%nrecv = 0
        win%nread = 0
        allocate( win%sqest(MAX(nbuf,numpro)*MAX(max_call,max_irr)) )
        allocate( win%rqest(MAX(nbuf,numpro)*MAX(max_call,max_irr)) )
!EOC
      end subroutine win_init_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send4d_ns --- Send 4d north/south ghost latitudes (real*8)
!
! !INTERFACE:
      subroutine mp_send4d_ns(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
      real(r8), intent(in):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Send 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Open(comm, t1_win)

! Send to south
      if ( jfirst > 1 ) then
        t1_win%src = gidu - 1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        t1_win%size_r = im*ng_s*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%dest = gidu - 1
        t1_win%offset_s = igosouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst, jfirst+ng_n-1, kfirst, klast, 1, nq,   &
                         ga_t1_s, ga_t1_r )
      endif
! Send to north
      if ( jlast < jm ) then
        t1_win%src = gidu + 1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        t1_win%size_r = im*ng_n*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%dest = gidu + 1
        t1_win%offset_s = igonorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast-ng_s+1, jlast, kfirst, klast, 1, nq,     &
                         ga_t1_s, ga_t1_r )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv4d_ns --- Receive 4d north/south ghost latitudes (real*8)
!
! !INTERFACE:
      subroutine mp_recv4d_ns(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Receive 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Close(comm, t1_win)

! Recv from south
      if ( jfirst > 1 ) then
        t1_win%src  = gidu-1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst-ng_s, jfirst-1,   kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        t1_win%src  = gidu+1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast+1,     jlast+ng_n, kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send4d_ns_r4 --- Send 4d north/south ghost latitudes (real*4)
!
! !INTERFACE:
      subroutine mp_send4d_ns_r4(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                                 ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
      real(r4), intent(in):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Send 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2005.03.20   Sawyer        Creation from mp_send4d_ns
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: mp_send4d_ns_r4 - r4 windows disabled - exiting'
        stop
#endif

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Open(comm, r4_win)

! Send to south
      if ( jfirst > 1 ) then
        r4_win%src = gidu - 1
        r4_win%offset_r = ifromsouth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        r4_win%size_r = im*ng_s*(klast-kfirst+1)*nq
        call Ga_RecvInit_r4(comm, r4_win, ga_r4_r)
        r4_win%dest = gidu - 1
        r4_win%offset_s = igosouth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst, jfirst+ng_n-1, kfirst, klast, 1, nq,   &
                         ga_r4_s, ga_r4_r )
      endif
! Send to north
      if ( jlast < jm ) then
        r4_win%src = gidu + 1
        r4_win%offset_r = ifromnorth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        r4_win%size_r = im*ng_n*(klast-kfirst+1)*nq
        call Ga_RecvInit_r4(comm, r4_win, ga_r4_r)
        r4_win%dest = gidu + 1
        r4_win%offset_s = igonorth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast-ng_s+1, jlast, kfirst, klast, 1, nq,     &
                         ga_r4_s, ga_r4_r )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send4d_ns_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv4d_ns_r4 --- Receive 4d north/south ghost latitudes (real*4)
!
! !INTERFACE:
      subroutine mp_recv4d_ns_r4(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Receive 4d north/south ghost latitudes (real*4)
!
! !REVISION HISTORY: 
!    2005.03.20   Sawyer        Creation from mp_recv4d_ns
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: mp_recv4d_ns_r4 - r4 windows disabled - exiting'
        stop
#endif

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Close(comm, r4_win)

! Recv from south
      if ( jfirst > 1 ) then
        r4_win%src  = gidu-1
        r4_win%offset_r = ifromsouth*idimsizz + (r4_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst-ng_s, jfirst-1,   kfirst, klast, 1, nq, &
                         ga_r4_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        r4_win%src  = gidu+1
        r4_win%offset_r = ifromnorth*idimsizz + (r4_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast+1,     jlast+ng_n, kfirst, klast, 1, nq, &
                         ga_r4_r  )
      endif

      call Win_Finalize(comm, r4_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv4d_ns_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send2_ns --- Send 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_send2_ns(comm, im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
      real(r8), intent(in):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(in):: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 
!
! !DESCRIPTION:
!
!     Send 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Open(comm, t1_win)

! Send to south
      if ( jfirst > 1 ) then
        t1_win%src  = gidu - 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%dest = gidu - 1
        t1_win%offset_s = igosouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif
! Send to north
      if ( jlast < jm ) then
        t1_win%src  = gidu + 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%dest = gidu + 1
        t1_win%offset_s = igonorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jlast,     jlast,    kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jlast,     jlast,    kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv2_ns --- Receive 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_recv2_ns(comm, im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(inout):: q2(im,jfirst-nd:jlast+nd,kfirst:klast)
!
! !DESCRIPTION:
!
!     Receive 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer j
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Close(comm, t1_win)

! Recv from south
      if ( jfirst > 1 ) then
        j = jfirst - 1
        t1_win%src  = gidu - 1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        j = jlast + 1
        t1_win%src  = gidu + 1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, q1, t1_win, im, jm, km, 2, & 
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d --- Send ghost region
!
! !INTERFACE:
      subroutine mp_send3d(comm, dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                  i1, i2, j1, j2, k1, k2, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Open(comm, t1_win)

! Init Recv src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
      endif
! Send ghost region
      if ( dest >= 0 .and. dest < numcomm ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d --- Recv ghost region
!
! !INTERFACE:
      subroutine mp_recv3d(comm, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                i1, i2, j1, j2, k1, k2, qout)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Close(comm, t1_win)

! Recv from src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, qout, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d_2 --- Send 2 ghost regions
!
! !INTERFACE:
      subroutine mp_send3d_2(comm, dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                        i1, i2, j1, j2, k1, k2, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q1(if:il, jf:jl, kf:kl)
      real(r8), intent(in):: q2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send two general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Open(comm, t1_win)

! Init Recv src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + t1_win%size_r 
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
      endif
! Send ghost region
      if ( dest >= 0 .and. dest < numcomm ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Put4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,  &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d_2 --- Recv 2 ghost regions
!
! !INTERFACE:
      subroutine mp_recv3d_2(comm, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                  i1, i2, j1, j2, k1, k2, qout1, qout2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout1(if:il, jf:jl, kf:kl)
      real(r8), intent(inout):: qout2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv two general 3d real*8 ghost regions
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Close(comm, t1_win)

! Recv from src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, qout1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Get4d_r8( comm, qout2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,   &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_barrier --- Synchronize all SPMD processes
!
! !INTERFACE:
      subroutine mp_barrier (comm)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !DESCRIPTION:
!
!     Synchronize all SPMD processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

        call MPI_BARRIER(comm, ierror)

!EOC
      end subroutine mp_barrier
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Open --- Open a communication window
!
! !INTERFACE:
      subroutine Win_Open(comm, win)
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Begin a communication epoch, by opening a comm window.
!     Update number of send calls on the window (win%ncall_s).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      win%ncall_s = win%ncall_s + 1

!EOC
      end subroutine Win_Open
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Close --- Close a communication window
!
! !INTERFACE:
      subroutine Win_Close(comm, win)
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     End a communication epoch, by closing a comm window.
!     Update number of receive calls on the window (win%ncall_r).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      win%ncall_r = win%ncall_r + 1

!EOC
      end subroutine Win_Close
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Finalize --- Reset a communication window after a comm epoch.
!
! !INTERFACE:
      subroutine Win_Finalize(comm, win)
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Complete a communication epoch and reset a comm window.
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY:
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (win%ncall_s == win%ncall_r) then
        call MPI_WAITALL(win%nsend, win%sqest, Stats, ierror)
        win%nsend = 0
        win%nrecv = 0
        win%nread = 0
        win%ncall_s = 0
        win%ncall_r = 0
      endif

!EOC
      end subroutine Win_Finalize
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_r8 --- Write to real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_r8 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga_s(win%size)
      real(r8), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
      integer :: gidu

      call MPI_COMM_RANK (comm, gidu, ierror)

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length 
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
      send_tag = gidu
      win%nsend = win%nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_r8, win%dest, &
                     send_tag, comm, win%sqest(win%nsend), ierror)

!EOC
      end subroutine Ga_Put4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_r8 --- Initiate real*8 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_r8( comm, win, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate real*8 Non-Blocking receive
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    03.06.06   Sawyer        Added else clause
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

      if (win%size >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        win%nrecv    = win%nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_r8, win%src, &
                       recv_tag, comm, win%rqest(win%nrecv), ierror)
      else
        write(iulog,*) "Fatal ga_recvinit_r8: receive window out of space - exiting"
        write(iulog,*) 'gid win%size win%offset_r win%size_r = ', gid,  &
                  win%size, win%offset_r, win%size_r
        stop
      endif

!EOC
      end subroutine Ga_RecvInit_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_r8 --- Read from real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_r8 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
      win%nread = win%nread + 1
      call MPI_WAIT(win%rqest(win%nread), Status, ierror)

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_r4 --- Write to real*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_r4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r4), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: ga_s(win%size)
      real(r4), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to real*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
      integer :: gidu

      call MPI_COMM_RANK (comm, gidu, ierror)

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_Put4d_r4 - r4 windows disabled - exiting'
        stop
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length 
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
      send_tag = gidu
      win%nsend = win%nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_r4, win%dest, &
                     send_tag, comm, win%sqest(win%nsend), ierror)

!EOC
      end subroutine Ga_Put4d_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_r4 --- Initiate real*4 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_r4( comm, win, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate real*8 Non-Blocking receive
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    03.06.06   Sawyer        Added else clause
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_RecvInit_r4 - r4 windows disabled - exiting'
        stop
#endif

      if (win%size >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        win%nrecv    = win%nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_r4, win%src, &
                       recv_tag, comm, win%rqest(win%nrecv), ierror)
      else
        write(iulog,*) "Fatal ga_recvinit_r4: receive window out of space - exiting"
        write(iulog,*) 'gid win%size win%offset_r win%size_r = ', gid,  &
                  win%size, win%offset_r, win%size_r
        stop
      endif

!EOC
      end subroutine Ga_RecvInit_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_r4 --- Read from real*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_r4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r4), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_Get4d_r4 - r4 windows disabled - exiting'
        stop
#endif

      win%nread = win%nread + 1
      call MPI_WAIT(win%rqest(win%nread), Status, ierror)

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_i4 --- Write to integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_i4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga_s(win%size)
      integer(i4), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
      integer :: gidu

      call MPI_COMM_RANK (comm, gidu, ierror)

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
      send_tag = gidu
      win%nsend = win%nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_i4, win%dest, &
                     send_tag, comm, win%sqest(win%nsend), ierror)

!EOC
      end subroutine Ga_Put4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_i4 --- Initiate integer*4 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_i4( comm, win, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate integer*4 Non-Blocking receive
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    06.05.21   Mirin         Added else clause
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

      if (win%size >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        win%nrecv    = win%nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_i4, win%src, &
                       recv_tag, comm, win%rqest(win%nrecv), ierror)
      else
        write(iulog,*) "Fatal ga_recvinit_i4: receive window out of space - exiting"
        write(iulog,*) 'gid win%size win%offset_r win%size_r = ', gid,  &
                  win%size, win%offset_r, win%size_r
        stop
      endif
!EOC
      end subroutine Ga_RecvInit_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_i4 --- Read from integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_i4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(inout)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      win%nread = win%nread + 1
      call MPI_WAIT(win%rqest(win%nread), Status, ierror)

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_r8 --- Broadcast an real*8 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_r8 ( comm, q, isize )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!     Broadcast an real*8 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      call MPI_BCAST(q, isize, mp_r8, 0, comm, ierror)

!EOC
      end subroutine Ga_Broadcast_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_r4 --- Broadcast an real*4 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_r4 ( comm, q, isize )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!     Broadcast an real*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_Broadcast_r4 - r4 windows disabled - exiting'
        stop
#endif

      call MPI_BCAST(q, isize, mp_r4, 0, comm, ierror)

!EOC
      end subroutine Ga_Broadcast_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_i4 --- Broadcast an integer*4 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_i4 ( comm, q, isize )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!     Broadcast an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      call MPI_BCAST(q, isize, mp_i4, 0, comm, ierror)

!EOC
      end subroutine Ga_Broadcast_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_r8 --- All to All of an real*8 1d global array
!
! !INTERFACE:
      subroutine Ga_AllToAll_r8 ( comm, q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: Gsize    ! Global size of array
      integer, intent(in)  :: Lsize    ! size of Local portion
      integer, intent(in)  :: istart   ! starting point
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!     All to All of a real*8 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      call MPI_ALLGATHER(q(istart), Lsize, mp_r8, q, Lsize, mp_r8, comm, ierror)

!EOC
      end subroutine Ga_AllToAll_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_r4 --- All to All of an real*4 1d global array
!
! !INTERFACE:
      subroutine Ga_AllToAll_r4 ( comm, q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      integer, intent(in)  :: Gsize   ! Global size of array
      integer, intent(in)  :: Lsize   ! size of Local portion
      integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!     All to All of an real*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_AllToAll_r4 - r4 windows disabled - exiting'
        stop
#endif

      call MPI_ALLGATHER(q(istart), Lsize, mp_r4, q, Lsize, mp_r4, comm, ierror)

!EOC
      end subroutine Ga_AllToAll_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_i4 --- All to All of an integer*4 1d global array
!
! !INTERFACE:
      subroutine Ga_AllToAll_i4 ( comm, q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      integer, intent(in)  :: Gsize   ! Global size of array
      integer, intent(in)  :: Lsize   ! size of Local portion
      integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!     All to All of an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      call MPI_ALLGATHER(q(istart), Lsize, mp_i4, q, Lsize, mp_i4, comm, ierror)

!EOC
      end subroutine Ga_AllToAll_i4
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: get_partneroffset --- Computes partneroffset/type from descriptor
!
! !INTERFACE:
      subroutine get_partneroffset ( comm, send_bl, recv_bl )

! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
! !INPUT/OUTPUT PARAMETERS:
      type(blockdescriptor), intent(inout)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(inout)  :: recv_bl(:) ! receive blocks

!
! !DESCRIPTION:
!     Compute partneroffsets/types from other blockdescriptor
!     information.  Used exclusively for irregular communication 
!     in PILGRIM.
!
! !REVISION HISTORY: 
!    03.10.31   Mirin       Creation
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer :: i, j, k, ns, pos, por, numpsq, ierror
      integer :: ami(numpro,numpro), am(numpro,numpro)
      integer mod_method, num_s, num_r

      num_s = size(send_bl)
      num_r = size(recv_bl)

      do j = 1, num_s
         send_bl(j)%partneroffset = 0
         send_bl(j)%partnertype = MPI_DATATYPE_NULL
      enddo
      do j = 1, num_r
         recv_bl(j)%partneroffset = 0
         recv_bl(j)%partnertype = MPI_DATATYPE_NULL
      enddo

      end subroutine get_partneroffset
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr --- Initiate communication of contiguous parcels
!
! !INTERFACE:
      subroutine mp_sendirr ( comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,      &
                              modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r8), intent(in) :: q1in(*)                  ! input array
      real(r8), optional, intent(in) :: q2in(*)        ! second input array
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: q1out(*)                ! output array
      real(r8), optional, intent(out) :: q2out(*)      ! second output array
!
! !DESCRIPTION:
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications here and placing
!     wait points in mp_recvirr; if 1, call swap routine with p2p messages; if 2, call swap
!     routine with a2a messages. 
!     Modc(2): if 1, then apply handshaking (don't send until corresponding receive is posted)
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!     Modc(4): maximum number of outstanding requests (applies to swap routines only)
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!    09.10.07   Worley      eliminated mpi_recv from handshake logic
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nparcels, offset_s, offset_r, ierr, mod_method
      integer p, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize, unitsize, offset_0
      integer i, j, send_tag, recv_tag, num_s, num_r
      integer :: offset_v (Max_Nparcels)
      integer :: hs_snd, hs_rcv(numpro), hs_rcvids(numpro)
      integer ipe2, ceil2num
      integer onetwo
      logical twovar
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_pid


#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .true.
         maxreq_local = -1
      endif

! Do not call mp_swapirr unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

      onetwo = 1
      twovar = .false.
      if (present(q2in)) then
         onetwo = 2
         twovar = .true.
      endif

    if (sw_local .gt. 0) then
         sw_alltoall = (sw_local .eq. 2)
         if (present(q2in)) then
            call mp_swapirr(comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,   &
                            sw_handshake=hs_local, sw_maxreq=maxreq_local,      &
                            sw_alltoall=sw_alltoall, sw_send=send_local)
         else
            call mp_swapirr(comm, send_bl, recv_bl, q1in, q1out,                &
                            sw_handshake=hs_local, sw_maxreq=maxreq_local,      &
                            sw_alltoall=sw_alltoall, sw_send=send_local)
         endif
    else

      call MPI_COMM_RANK (comm, comm_pid, ierr)

      hs_snd = 1
      ceil2num = ceil2(numpro)

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      r8_win%ncall_s = r8_win%ncall_s + 1
     if (mod_method .gt. 0) then
!
! mpi derived types
      if (r8_win%ncall_s .gt. MaxTrf-onetwo+1) then
         write(iulog,*) "mp_sendirr: derived type handle count exceeded - exiting"
         write(iulog,*) "r8_win%ncall_s MaxTrf = ", r8_win%ncall_s, MaxTrf
         stop
      endif
!
! MPI: Irecv over all processes
!
      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               if (ipe-1 /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, ipe-1, comm_pid, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      OutHandle(:,r8_win%ncall_s) = MPI_REQUEST_NULL
      if (twovar) OutHandle(:,r8_win%ncall_s+1) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_r) cycle
!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          recv_tag = ipe-1 + modcam_tagoffset
          call mpi_irecv( q1out, 1, recv_bl(ipe)%type, ipe-1, recv_tag,     &
                          comm, OutHandle(ipe,r8_win%ncall_s), ierr )
          if (twovar) then
             call mpi_irecv( q2out, 1, recv_bl(ipe)%type, ipe-1, recv_tag,     &
                             comm, OutHandle(ipe,r8_win%ncall_s+1), ierr )
          endif
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
               call MPI_SEND ( hs_snd, 1, mp_i4, ipe-1, ipe-1, comm, ierr )
          endif
        endif
      enddo

!
! MPI: Isend/Send over all processes; use risend/rsend with hs
!
      InHandle(:,r8_win%ncall_s) = MPI_REQUEST_NULL
      if (twovar) InHandle(:,r8_win%ncall_s+1) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_s) cycle

!
! Send the individual buffers with non-blocking sends
!
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          send_tag = comm_pid + modcam_tagoffset
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
                call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
             if (send_local) then
                call mpi_rsend( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                comm, ierr )
             else
                call mpi_irsend( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                 comm, InHandle(ipe,r8_win%ncall_s), ierr )
             endif
             if (twovar) then
                if (send_local) then
                   call mpi_rsend( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, ierr )
                else
                   call mpi_irsend( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                    comm, InHandle(ipe,r8_win%ncall_s+1), ierr )
                endif
             endif
          else
             if (send_local) then
                call mpi_send( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                               comm, ierr )
             else
                call mpi_isend( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                comm, InHandle(ipe,r8_win%ncall_s), ierr )
             endif
             if (twovar) then
                if (send_local) then
                   call mpi_send( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                  comm, ierr )
                else
                   call mpi_isend( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, InHandle(ipe,r8_win%ncall_s+1), ierr )
                endif
             endif
          endif
        endif
      enddo
     else

! temporary contiguous buffers

      if (r8_win%ncall_s .gt. max_irr-onetwo+1) then
         write(iulog,*) "mp_sendirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_s max_irr = ", r8_win%ncall_s, max_irr
         stop
      endif
      unitsize = r8_win%size/max_irr

! issue call to receive data in global receive buffer
      offset_0 = (r8_win%ncall_s-1)*unitsize
      offset_s = offset_0
      offset_r = offset_0

      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            qsize = onetwo*send_bl(ipe)%Tot_Size
            if (qsize .ne. 0) then
               r8_win%dest = ipe-1
               send_tag = comm_pid + modcam_tagoffset
               if (r8_win%dest /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, r8_win%dest, send_tag, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         r8_win%size_r = onetwo*recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            r8_win%src = ipe-1
            if (onetwo*unitsize >= offset_r-offset_0) then
              recv_tag = r8_win%src + modcam_tagoffset
              qsize    = r8_win%size_r
              r8_win%nrecv    = r8_win%nrecv + 1
              call MPI_IRECV(ga_r8_r(r8_win%offset_r+1), qsize, mp_r8, r8_win%src, &
                             recv_tag, comm, r8_win%rqest(r8_win%nrecv), ierror)
              if (hs_local) then
                 if (r8_win%src /= comm_pid) &
                   call MPI_SEND ( hs_snd, 1, mp_i4, r8_win%src, recv_tag, comm, ierror)
              endif
            else
              write(iulog,*) "Fatal mp_sendirr: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif
         endif
      enddo

! gather data into global send buffer
      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_s) cycle
         qsize = onetwo*send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
            r8_win%dest = ipe-1
            r8_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s-offset_0 .gt. onetwo*unitsize) then
              write(iulog,*) "Fatal mp_sendirr: send window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_s, offset_0
              stop
            endif

            offset_v(1) = r8_win%offset_s
            do j = 2, send_bl(ipe)%nparcels
               offset_v(j) = offset_v(j-1) + send_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, send_bl(ipe)%nparcels
               do i = 1, send_bl(ipe)%blocksizes(j)
                  ga_r8_s(offset_v(j)+i) = q1in(send_bl(ipe)%displacements(j)+i)
               enddo
            enddo
            if (twovar) then
               do j = 1, send_bl(ipe)%nparcels
                  do i = 1, send_bl(ipe)%blocksizes(j)
                     ga_r8_s(send_bl(ipe)%Tot_Size+offset_v(j)+i) = q2in(send_bl(ipe)%displacements(j)+i)
                  enddo
               enddo
            endif

! nonblocking send
            send_tag = comm_pid + modcam_tagoffset
            r8_win%nsend = r8_win%nsend + 1
            if (hs_local) then
               if (r8_win%dest /= comm_pid) &
                  call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
               if (send_local) then
                  call MPI_RSEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_IRSEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, r8_win%sqest(r8_win%nsend), ierr)
               endif
            else
               if (send_local) then
                  call MPI_SEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_ISEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, r8_win%sqest(r8_win%nsend), ierr)
               endif
            endif
         endif
      enddo

     endif   !  mod_method

      if (twovar) r8_win%ncall_s = r8_win%ncall_s + 1

    endif   !  sw_local

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr --- Finalize communication of contiguous parcels
!
! !INTERFACE:
      subroutine mp_recvirr ( comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,      &
                              modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r8), intent(in) :: q1in(*)                  ! input array
      real(r8), optional, intent(in) :: q2in(*)        ! second input array
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q1out(*)                ! output array
      real(r8), optional, intent(inout) :: q2out(*)      ! second output array
!
! !DESCRIPTION:
!     Complete transfer of a generalized region initiated by {\tt mp\_sendirr}.
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications in mp_sendirr and
!     placing wait points here; otherwise don't do anything - mp_swapirr is called from mp_sendirr.
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, offset_r, mod_method
      integer unitsize, offset_0
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer i, j, num_r, num_s
      integer :: offset_v (Max_Nparcels)
      integer ipe2, ceil2num
      integer onetwo
      logical twovar
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_size, comm_pid

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .true.
         maxreq_local = -1
      endif

! Do not call mp_swapirr (hence return) unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

! Return if swap_irr
      if (sw_local .gt. 0) return

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      onetwo = 1
      twovar = .false.
      if (present(q2in)) then
         onetwo = 2
         twovar = .true.
      endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

      ceil2num = ceil2(numpro)

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      r8_win%ncall_r = r8_win%ncall_r + 1

    if (mod_method .gt. 0) then

! mpi derived types
      if (r8_win%ncall_r .gt. MaxTrf-onetwo+1) then
         write(iulog,*) "mp_recvirr: derived type handle count exceeded - exiting"
         write(iulog,*) "r8_win%ncall_r MaxTrf = ", r8_win%ncall_r, MaxTrf
         stop
      endif

      if (num_s .gt. 0 .and. (.not. send_local)) then
         CALL MPI_WAITALL( comm_size, InHandle(:,r8_win%ncall_r), InStats, Ierr )
         if (twovar) then
            CALL MPI_WAITALL( comm_size, InHandle(:,r8_win%ncall_r+1), InStats, Ierr )
         endif
      endif
      if (num_r .gt. 0) then
         CALL MPI_WAITALL( comm_size, OutHandle(:,r8_win%ncall_r), OutStats, Ierr )
         if (twovar) then
            CALL MPI_WAITALL( comm_size, OutHandle(:,r8_win%ncall_r+1), OutStats, Ierr )
         endif
      endif

    else

! temporary contiguous buffer / global window

      if (r8_win%ncall_r .gt. max_irr-onetwo+1) then
         write(iulog,*) "mp_recvirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_r max_irr = ", r8_win%ncall_r, max_irr
         stop
      endif
      unitsize = r8_win%size/max_irr

! scatter data from global receive buffer to final destination
      offset_0 = (r8_win%ncall_r-1)*unitsize
      offset_r = offset_0

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         r8_win%size_r = onetwo*recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r-offset_0 .gt. onetwo*unitsize) then
              write(iulog,*) "Fatal mp_recvirr: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif

            r8_win%nread = r8_win%nread + 1
            call MPI_WAIT(r8_win%rqest(r8_win%nread), Status, ierr)

            offset_v(1) = r8_win%offset_r
            do j = 2, recv_bl(ipe)%Nparcels
               offset_v(j) = offset_v(j-1) + recv_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, recv_bl(ipe)%Nparcels
               do i = 1, recv_bl(ipe)%blocksizes(j)
                  q1out(recv_bl(ipe)%displacements(j)+i) = ga_r8_r(offset_v(j)+i)
               enddo
            enddo
            if (twovar) then
            do j = 1, recv_bl(ipe)%Nparcels
               do i = 1, recv_bl(ipe)%blocksizes(j)
                  q2out(recv_bl(ipe)%displacements(j)+i) = ga_r8_r(recv_bl(ipe)%Tot_Size+offset_v(j)+i)
               enddo
            enddo
            endif

         endif
      enddo

      if ((r8_win%ncall_s == r8_win%ncall_r + onetwo - 1) .and. (.not. send_local)) then
         call MPI_WAITALL(r8_win%nsend, r8_win%sqest, Stats, ierror)
      endif

    endif    !    mod_method .gt. 0

    if (twovar) r8_win%ncall_r = r8_win%ncall_r + 1

    if (r8_win%ncall_s == r8_win%ncall_r) then
       r8_win%nsend = 0
       r8_win%nrecv = 0
       r8_win%nread = 0
       r8_win%ncall_s = 0
       r8_win%ncall_r = 0
    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr_r4 --- Initiate communication of contiguous parcels - r4
!
! !INTERFACE:
      subroutine mp_sendirr_r4 ( comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,      &
                                 modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r4), intent(in) :: q1in(*)                  ! input array
      real(r4), optional, intent(in) :: q2in(*)        ! second input array
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests

! !OUTPUT PARAMETERS:
      real(r4), intent(out) :: q1out(*)                ! output array
      real(r4), optional, intent(out) :: q2out(*)      ! second output array
!
! !DESCRIPTION:
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications here and placing
!     wait points in mp_recvirr; if 1, call swap routine with p2p messages; if 2, call swap
!     routine with a2a messages. 
!     Modc(2): if 1, then apply handshaking (don't send until corresponding receive is posted)
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!     Modc(4): maximum number of outstanding requests (applies to swap routines only)
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       No-op version
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      write(iulog,*) 'Mod_comm: mp_sendirr_r4 - r4 no longer supported - exiting'
      stop

!EOC
      end subroutine mp_sendirr_r4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr_r4 --- Finalize communication of contiguous parcels - r4
!
! !INTERFACE:
      subroutine mp_recvirr_r4 ( comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,      &
                                 modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r4), intent(in) :: q1in(*)                  ! input array
      real(r4), optional, intent(in) :: q2in(*)        ! second input array
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests
! !INPUT/OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q1out(*)                ! output array
      real(r4), optional, intent(inout) :: q2out(*)      ! second output array
!
! !DESCRIPTION:
!     Complete transfer of a generalized region initiated by {\tt mp\_sendirr}.
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications in mp_sendirr and
!     placing wait points here; otherwise don't do anything - mp_swapirr is called from mp_sendirr.
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       No-op version
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      write(iulog,*) 'Mod_comm: mp_recvirr_r4 - r4 no longer supported - exiting'
      stop

!EOC
      end subroutine mp_recvirr_r4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr_i4 --- Initiate communication of contiguous parcels - i4
!
! !INTERFACE:
      subroutine mp_sendirr_i4 ( comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,      &
                                 modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer(i4), intent(in) :: q1in(*)                  ! input array
      integer(i4), optional, intent(in) :: q2in(*)        ! second input array
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests

! !OUTPUT PARAMETERS:
      integer(i4), intent(out) :: q1out(*)                ! output array
      integer(i4), optional, intent(out) :: q2out(*)      ! second output array
!
! !DESCRIPTION:
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications here and placing
!     wait points in mp_recvirr; if 1, call swap routine with p2p messages; if 2, call swap
!     routine with a2a messages. 
!     Modc(2): if 1, then apply handshaking (don't send until corresponding receive is posted)
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!     Modc(4): maximum number of outstanding requests (applies to swap routines only)
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!    09.10.07   Worley      eliminated mpi_recv from handshake logic
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nparcels, offset_s, offset_r, ierr, mod_method
      integer p, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize, unitsize, offset_0
      integer i, j, send_tag, recv_tag, num_s, num_r
      integer :: offset_v (Max_Nparcels)
      integer :: hs_snd, hs_rcv(numpro), hs_rcvids(numpro)
      integer ipe2, ceil2num
      integer onetwo
      logical twovar
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_pid

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .false.
         maxreq_local = -1
      endif

! Do not call mp_swapirr_i4 unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

      onetwo = 1
      twovar = .false.
      if (present(q2in)) then
         onetwo = 2
         twovar = .true.
      endif

    if (sw_local .gt. 0) then
         sw_alltoall = (sw_local .eq. 2)
         if (present(q2in)) then
            call mp_swapirr_i4(comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,   &
                               sw_handshake=hs_local, sw_maxreq=maxreq_local,      &
                               sw_alltoall=sw_alltoall, sw_send=send_local)
         else
            call mp_swapirr_i4(comm, send_bl, recv_bl, q1in, q1out,                &
                               sw_handshake=hs_local, sw_maxreq=maxreq_local,      &
                               sw_alltoall=sw_alltoall, sw_send=send_local)
         endif
    else

      call MPI_COMM_RANK (comm, comm_pid, ierr)

      hs_snd = 1
      ceil2num = ceil2(numpro)

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      mod_method = recv_bl(1)%method

      i4_win%ncall_s = i4_win%ncall_s + 1
     if (mod_method .gt. 0) then
!
! mpi derived types
      if (i4_win%ncall_s .gt. MaxTrf-onetwo+1) then
         write(iulog,*) "mp_sendirr_i4: derived type handle count exceeded - exiting"
         write(iulog,*) "i4_win%ncall_s MaxTrf = ", i4_win%ncall_s, MaxTrf
         stop
      endif
!
! MPI: Irecv over all processes
!
      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               if (ipe-1 /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, ipe-1, comm_pid, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      OutHandle(:,i4_win%ncall_s) = MPI_REQUEST_NULL
      if (twovar) OutHandle(:,i4_win%ncall_s+1) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_r) cycle
!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          recv_tag = ipe-1 + modcam_tagoffset
          call mpi_irecv( q1out, 1, recv_bl(ipe)%type, ipe-1, recv_tag,     &
                          comm, OutHandle(ipe,i4_win%ncall_s), ierr )
          if (twovar) then
             call mpi_irecv( q2out, 1, recv_bl(ipe)%type, ipe-1, recv_tag,     &
                             comm, OutHandle(ipe,i4_win%ncall_s+1), ierr )
          endif
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
               call MPI_SEND ( hs_snd, 1, mp_i4, ipe-1, ipe-1, comm, ierr )
          endif
        endif
      enddo

!
! MPI: Isend/Send over all processes; use risend/rsend with hs
!
      InHandle(:,i4_win%ncall_s) = MPI_REQUEST_NULL
      if (twovar) InHandle(:,i4_win%ncall_s+1) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_s) cycle

!
! Send the individual buffers with non-blocking sends
!
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          send_tag = comm_pid + modcam_tagoffset
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
                call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
             if (send_local) then
                call mpi_rsend( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                comm, ierr )
             else
                call mpi_irsend( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                 comm, InHandle(ipe,i4_win%ncall_s), ierr )
             endif
             if (twovar) then
                if (send_local) then
                   call mpi_rsend( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, ierr )
                else
                   call mpi_irsend( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                    comm, InHandle(ipe,i4_win%ncall_s+1), ierr )
                endif
             endif
          else
             if (send_local) then
                call mpi_send( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                               comm, ierr )
             else
                call mpi_isend( q1in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                comm, InHandle(ipe,i4_win%ncall_s), ierr )
             endif
             if (twovar) then
                if (send_local) then
                   call mpi_send( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                  comm, ierr )
                else
                   call mpi_isend( q2in, 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, InHandle(ipe,i4_win%ncall_s+1), ierr )
                endif
             endif
          endif
        endif
      enddo
     else

! temporary contiguous buffers

      if (i4_win%ncall_s .gt. max_irr-onetwo+1) then
         write(iulog,*) "mp_sendirr_i4: insufficient window storage - exiting"
         write(iulog,*) "i4_win%ncall_s max_irr = ", i4_win%ncall_s, max_irr
         stop
      endif
      unitsize = i4_win%size/max_irr

! issue call to receive data in global receive buffer
      offset_0 = (i4_win%ncall_s-1)*unitsize
      offset_s = offset_0
      offset_r = offset_0

      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            qsize = onetwo*send_bl(ipe)%Tot_Size
            if (qsize .ne. 0) then
               i4_win%dest = ipe-1
               send_tag = comm_pid + modcam_tagoffset
               if (i4_win%dest /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, i4_win%dest, send_tag, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         i4_win%size_r = onetwo*recv_bl(ipe)%Tot_Size
         if (i4_win%size_r .ne. 0) then
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            i4_win%src = ipe-1
            if (onetwo*unitsize >= offset_r-offset_0) then
              recv_tag = i4_win%src + modcam_tagoffset
              qsize    = i4_win%size_r
              i4_win%nrecv    = i4_win%nrecv + 1
              call MPI_IRECV(ga_i4_r(i4_win%offset_r+1), qsize, mp_i4, i4_win%src, &
                             recv_tag, comm, i4_win%rqest(i4_win%nrecv), ierror)
              if (hs_local) then
                 if (i4_win%src /= comm_pid) &
                   call MPI_SEND ( hs_snd, 1, mp_i4, i4_win%src, recv_tag, comm, ierror)
              endif
            else
              write(iulog,*) "Fatal mp_sendirr_i4: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif
         endif
      enddo
! gather data into global send buffer
      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_s) cycle
         qsize = onetwo*send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
            i4_win%dest = ipe-1
            i4_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s-offset_0 .gt. onetwo*unitsize) then
              write(iulog,*) "Fatal mp_sendirr_i4: send window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_s, offset_0
              stop
            endif

            offset_v(1) = i4_win%offset_s
            do j = 2, send_bl(ipe)%nparcels
               offset_v(j) = offset_v(j-1) + send_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, send_bl(ipe)%nparcels
               do i = 1, send_bl(ipe)%blocksizes(j)
                  ga_i4_s(offset_v(j)+i) = q1in(send_bl(ipe)%displacements(j)+i)
               enddo
            enddo
            if (twovar) then
               do j = 1, send_bl(ipe)%nparcels
                  do i = 1, send_bl(ipe)%blocksizes(j)
                     ga_i4_s(send_bl(ipe)%Tot_Size+offset_v(j)+i) = q2in(send_bl(ipe)%displacements(j)+i)
                  enddo
               enddo
            endif

! nonblocking send
            send_tag = comm_pid + modcam_tagoffset
            i4_win%nsend = i4_win%nsend + 1
            if (hs_local) then
               if (i4_win%dest /= comm_pid) &
                  call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
               if (send_local) then
                  call MPI_RSEND(ga_i4_s(i4_win%offset_s+1), qsize, mp_i4, i4_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_IRSEND(ga_i4_s(i4_win%offset_s+1), qsize, mp_i4, i4_win%dest, &
                                 send_tag, comm, i4_win%sqest(i4_win%nsend), ierr)
               endif
            else
               if (send_local) then
                  call MPI_SEND(ga_i4_s(i4_win%offset_s+1), qsize, mp_i4, i4_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_ISEND(ga_i4_s(i4_win%offset_s+1), qsize, mp_i4, i4_win%dest, &
                                 send_tag, comm, i4_win%sqest(i4_win%nsend), ierr)
               endif
            endif
         endif
      enddo

     endif   !  mod_method

      if (twovar) i4_win%ncall_s = i4_win%ncall_s + 1

    endif   !  sw_local

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendirr_i4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr_i4 --- Finalize communication of contiguous parcels - i4
!
! !INTERFACE:
      subroutine mp_recvirr_i4 ( comm, send_bl, recv_bl, q1in, q1out, q2in, q2out,      &
                                 modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer(i4), intent(in) :: q1in(*)                  ! input array
      integer(i4), optional, intent(in) :: q2in(*)        ! second input array
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests

! !INPUT/OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q1out(*)                ! output array
      integer(i4), optional, intent(inout) :: q2out(*)      ! second output array
!
! !DESCRIPTION:
!     Complete transfer of a generalized region initiated by {\tt mp\_sendirr}.
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications in mp_sendirr and
!     placing wait points here; otherwise don't do anything - mp_swapirr is called from mp_sendirr.
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, offset_r, mod_method
      integer unitsize, offset_0
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer i, j, num_r, num_s
      integer :: offset_v (Max_Nparcels)
      integer ipe2, ceil2num
      integer onetwo
      logical twovar
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_size, comm_pid

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .false.
         maxreq_local = -1
      endif

! Do not call mp_swapirr_i4 (hence return) unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

! Return if swap_irr
      if (sw_local .gt. 0) return

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      onetwo = 1
      twovar = .false.
      if (present(q2in)) then
         onetwo = 2
         twovar = .true.
      endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

      ceil2num = ceil2(numpro)

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      mod_method = recv_bl(1)%method

      i4_win%ncall_r = i4_win%ncall_r + 1

    if (mod_method .gt. 0) then

! mpi derived types
      if (i4_win%ncall_r .gt. MaxTrf-onetwo+1) then
         write(iulog,*) "mp_recvirr_i4: derived type handle count exceeded - exiting"
         write(iulog,*) "i4_win%ncall_r MaxTrf = ", i4_win%ncall_r, MaxTrf
         stop
      endif

      if (num_s .gt. 0 .and. (.not. send_local)) then
         CALL MPI_WAITALL( comm_size, InHandle(:,i4_win%ncall_r), InStats, Ierr )
         if (twovar) then
            CALL MPI_WAITALL( comm_size, InHandle(:,i4_win%ncall_r+1), InStats, Ierr )
         endif
      endif
      if (num_r .gt. 0) then
         CALL MPI_WAITALL( comm_size, OutHandle(:,i4_win%ncall_r), OutStats, Ierr )
         if (twovar) then
            CALL MPI_WAITALL( comm_size, OutHandle(:,i4_win%ncall_r+1), OutStats, Ierr )
         endif
      endif

    else

! temporary contiguous buffer / global window

      if (i4_win%ncall_r .gt. max_irr-onetwo+1) then
         write(iulog,*) "mp_recvirr_i4: insufficient window storage - exiting"
         write(iulog,*) "i4_win%ncall_r max_irr = ", i4_win%ncall_r, max_irr
         stop
      endif
      unitsize = i4_win%size/max_irr

! scatter data from global receive buffer to final destination
      offset_0 = (i4_win%ncall_r-1)*unitsize
      offset_r = offset_0

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         i4_win%size_r = onetwo*recv_bl(ipe)%Tot_Size
         if (i4_win%size_r .ne. 0) then
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            if (offset_r-offset_0 .gt. onetwo*unitsize) then
              write(iulog,*) "Fatal mp_recvirr_i4: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif

            i4_win%nread = i4_win%nread + 1
            call MPI_WAIT(i4_win%rqest(i4_win%nread), Status, ierr)

            offset_v(1) = i4_win%offset_r
            do j = 2, recv_bl(ipe)%Nparcels
               offset_v(j) = offset_v(j-1) + recv_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, recv_bl(ipe)%Nparcels
               do i = 1, recv_bl(ipe)%blocksizes(j)
                  q1out(recv_bl(ipe)%displacements(j)+i) = ga_i4_r(offset_v(j)+i)
               enddo
            enddo
            if (twovar) then
            do j = 1, recv_bl(ipe)%Nparcels
               do i = 1, recv_bl(ipe)%blocksizes(j)
                  q2out(recv_bl(ipe)%displacements(j)+i) = ga_i4_r(recv_bl(ipe)%Tot_Size+offset_v(j)+i)
               enddo
            enddo
            endif

         endif
      enddo

      if ((i4_win%ncall_s == i4_win%ncall_r + onetwo - 1) .and. (.not. send_local)) then
         call MPI_WAITALL(i4_win%nsend, i4_win%sqest, Stats, ierror)
      endif

    endif    !    mod_method .gt. 0

    if (twovar) i4_win%ncall_r = i4_win%ncall_r + 1

    if (i4_win%ncall_s == i4_win%ncall_r) then
       i4_win%nsend = 0
       i4_win%nrecv = 0
       i4_win%nread = 0
       i4_win%ncall_s = 0
       i4_win%ncall_r = 0
    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvirr_i4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_swapirr --- Write r8 contiguous parcels to global array
!                           using XOR swap ordering
!
! !INTERFACE:
      subroutine mp_swapirr ( comm, send_bl, recv_bl, a1in, a1out, &
                               a2in, a2out, sw_handshake, sw_maxreq, &
                               sw_alltoall, sw_send )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm                     ! communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r8), intent(in) :: a1in(*)                  ! local data segment
      real(r8), optional, intent(in) :: a2in(*)        ! local data segment
      logical, optional, intent(in) :: sw_handshake    ! use flow control and 
                                                       !  ready send
      integer, optional, intent(in) :: sw_maxreq       ! maximum number of outstanding
                                                       !  MPI requests
      logical, optional, intent(in) :: sw_alltoall     ! use mpi_alltoall
      logical, optional, intent(in) :: sw_send         ! use mpi_send instead of isend

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: a1out(*)                ! local output segment
      real(r8), optional, intent(out) :: a2out(*)      ! local output segment
!
! !DESCRIPTION:
!     XOR-ordered version of all-to-all communication
!
! WARNING: mod_comm parameter max_irr might need to be set larger than expected
!          when swapping two variables; specifically, max_irr must be at least
!          as large as the incoming r8_win%ncall_s + the number of variables to
!          be swapped
!
! !REVISION HISTORY: 
!    08.06.30   Worley      original: derived from mp_sendirr, but using 
!                            swapm logic and XOR swap order 
!    08.08.22   Worley      removed swapm; reimplemented with native MPI,
!                            added flow control/ready send option and maxreq
!                            throttling, added alltoall option
!    09.10.07   Worley      eliminated mpi_recv from handshake logic
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: i, j, p, istep, num_s, num_r
      integer :: comm_pid, comm_size, steps, ierr
      integer :: ipe, offset_s, offset_r, offset_0, unitsize, onetwo

      integer :: arr_sndlths(0:numpro-1), arr_rcvlths(0:numpro-1)
      integer :: sndlths(0:numpro-1), sdispls(0:numpro-1)
      integer :: rcvlths(0:numpro-1), rdispls(0:numpro-1)
      integer :: swapids(numpro) 
      integer :: sndids(numpro)  ! nonblocking MPI send request ids
      integer :: rcvids(numpro)  ! nonblocking MPI recv request ids
      integer :: hs_snd, hs_rcv(numpro)! handshake variables (send/receive)
      integer :: hs_rcvids(numpro) ! nonblocking MPI handshake recv request ids
      integer :: InStats(numpro*MPI_STATUS_SIZE)
      integer :: OutStats(numpro*MPI_STATUS_SIZE)

      integer :: offset_v

      integer :: rstep

      integer :: maxreq, maxreqh
      logical :: handshake, alltoall, sendd

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

!     num_s = 0 if this process is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this process is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      if ( present(a2in) .and. (.not. present(a2out)) ) then
         write(iulog,*) "Fatal mp_swapirr: a2in specified, but a2out missing - exiting"
         stop
      endif

      if ( (.not. present(a2in)) .and. present(a2out)) then
         write(iulog,*) "Fatal mp_swapirr: a2out specified, but a2in missing - exiting"
         stop
      endif

      if ( present(sw_handshake) ) then
         handshake = sw_handshake
         hs_snd = 1
      else
         handshake = .false.
      endif

      if ( present(sw_alltoall) ) then
         alltoall = sw_alltoall
      else
         alltoall = .false.
      endif

      if ( present(sw_send) ) then
         sendd = sw_send
      else
         sendd = .false.
      endif

      onetwo = 1
      if (present(a2in)) onetwo = 2
      unitsize = r8_win%size/max_irr

! advance to unused portion of storage window
      r8_win%ncall_s = r8_win%ncall_s + 1

      if (r8_win%ncall_s .gt. max_irr-onetwo+1) then
         write(iulog,*) "mp_swapirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_s max_irr = ", r8_win%ncall_s, max_irr
         stop
      endif

! calculate send lengths and displacements
      offset_0 = (r8_win%ncall_s-1)*unitsize
      offset_s = offset_0
      sndlths(:) = 0
      sdispls(:) = 0
      arr_sndlths(:) = 0
      do ipe=1, num_s
         sndlths(ipe-1) = send_bl(ipe)%Tot_Size
         sdispls(ipe-1) = offset_s
         if (sndlths(ipe-1) .ne. 0) then

            ! pack first array
            offset_s = offset_s + sndlths(ipe-1)
            if (offset_s-offset_0 .gt. onetwo*unitsize) then
              write(iulog,*) "Fatal mp_swapirr: send window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                             ipe, unitsize, offset_s, offset_0
              stop
            endif

            arr_sndlths(ipe-1) = sndlths(ipe-1)

            ! calculate for second array (if it exists)
            if ( present(a2in) ) then

               offset_s = offset_s + sndlths(ipe-1)
               if (offset_s-offset_0 .gt. onetwo*unitsize) then
                 write(iulog,*) "Fatal mp_swapirr: send window out of space - exiting"
                 write(iulog,*) '2 comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                                ipe, unitsize, offset_s, offset_0
                 stop
               endif

               sndlths(ipe-1) = sndlths(ipe-1) + arr_sndlths(ipe-1)

            endif

         endif
      enddo

! calculate receive lengths and displacements
      offset_r = offset_0
      rcvlths(:) = 0
      rdispls(:) = 0
      arr_rcvlths(:) = 0
      do ipe=1, num_r
         rcvlths(ipe-1) = recv_bl(ipe)%Tot_Size
         rdispls(ipe-1) = offset_r
         if (rcvlths(ipe-1) .ne. 0) then

            offset_r = offset_r + rcvlths(ipe-1)
            if (onetwo*unitsize < offset_r-offset_0) then
              write(iulog,*) "Fatal mp_swapirr: receive window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif

            arr_rcvlths(ipe-1) = rcvlths(ipe-1)

            ! compute for second array (if it exists)
            if ( present(a2out) ) then

               offset_r = offset_r + rcvlths(ipe-1)
               if (onetwo*unitsize < offset_r-offset_0) then
                 write(iulog,*) "Fatal mp_swapirr: receive window out of space - exiting"
                 write(iulog,*) '2 comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                          ipe, unitsize, offset_r, offset_0
                 stop
               endif

               rcvlths(ipe-1) = rcvlths(ipe-1) + arr_rcvlths(ipe-1)

            endif

         endif
      enddo

! Calculate swap partners and number of steps in point-to-point
! implementations of alltoall algorithm.
      steps = 0
      do ipe=1,ceil2(comm_size)-1
         p = pair(comm_size,ipe,comm_pid)
         if (p >= 0) then
            if (sndlths(p) > 0 .or. rcvlths(p) > 0) then
               steps = steps + 1
               swapids(steps) = p
            end if
         end if
      end do

      if (.not. alltoall) then

         sndids(1:steps) = MPI_REQUEST_NULL
         rcvids(1:steps) = MPI_REQUEST_NULL

         if (steps .eq. 0) then
            maxreq  = 0
            maxreqh = 0
         elseif (steps .eq. 1) then
            maxreq  = 1
            maxreqh = 1
         else
            if ( present(sw_maxreq) ) then
               if ((sw_maxreq .le. steps) .and. (sw_maxreq .ge. 0)) then
                  maxreq  = sw_maxreq
                  if (maxreq > 1) then
                     maxreqh = maxreq/2
                  else
                     maxreq  = 2
                     maxreqh = 1
                  endif
               else
                  maxreq  = steps
                  maxreqh = steps
               endif
            else
               maxreq  = steps
               maxreqh = steps
            endif
         endif

! Post initial handshake receive requests
         if (handshake) then
            do istep=1,maxreq
               p = swapids(istep)
               if (sndlths(p) > 0) then
                  call mpi_irecv  ( hs_rcv(istep), 1, mp_i4, p, comm_pid, comm, &
                                    hs_rcvids(istep), ierr )
               endif
            enddo
         endif

! Post initial receive requests
         do istep=1,maxreq
            p = swapids(istep)
            if (rcvlths(p) > 0) then
               offset_r = rdispls(p)+1
               call mpi_irecv ( ga_r8_r(offset_r), rcvlths(p), mp_r8, &
                                p, p, comm, rcvids(istep), ierr )
               if (handshake) then
                  call mpi_send( hs_snd, 1, mp_i4, p, p, comm, &
                                 ierr )
               endif
            endif
         enddo
         rstep = maxreq
!
      endif

! gather data into global send buffer
      do istep=1,steps
         p = swapids(istep)

         if (sndlths(p) .ne. 0) then
            offset_v = sdispls(p)
            do j = 1, send_bl(p+1)%nparcels
               do i = 1, send_bl(p+1)%blocksizes(j)
                  ga_r8_s(offset_v+i) = a1in(send_bl(p+1)%displacements(j)+i)
               enddo
               offset_v = offset_v + send_bl(p+1)%blocksizes(j)
            enddo

            ! pack second array (if it exists)
            if ( present(a2in) ) then
               offset_v = sdispls(p) + arr_sndlths(p)
               do j = 1, send_bl(p+1)%nparcels
                  do i = 1, send_bl(p+1)%blocksizes(j)
                     ga_r8_s(offset_v+i) = a2in(send_bl(p+1)%displacements(j)+i)
                  enddo
                  offset_v = offset_v + send_bl(p+1)%blocksizes(j)
               enddo
            endif

         endif

         if (.not. alltoall) then

! Submit new i(r)send request
            offset_s = sdispls(p)+1
            if (sndlths(p) > 0) then
               if (handshake) then
                  call mpi_wait( hs_rcvids(istep), MPI_STATUS_IGNORE, ierr )
                  if (sendd) then
                     call mpi_rsend( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_irsend( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               else
                  if (sendd) then
                     call mpi_send ( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_isend ( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               endif
            endif

            if (istep > maxreqh) then
! Wait for oldest irecv request to complete
               call mpi_wait( rcvids(istep-maxreqh), OutStats, ierr )

               if (rstep < steps) then
                  rstep = rstep + 1
                  p = swapids(rstep)

! Submit a new handshake irecv request
                  if (handshake) then
                     if (sndlths(p) > 0) then
                        call mpi_irecv( hs_rcv(rstep), 1, mp_i4, p, comm_pid, comm, &
                                        hs_rcvids(rstep), ierr )
                     endif
                  endif

! Submit a new irecv request
                  if (rcvlths(p) > 0) then
                     offset_r = rdispls(p)+1
                     call mpi_irecv( ga_r8_r(offset_r), rcvlths(p), mp_r8, &
                                     p, p, comm, rcvids(rstep), ierr )
                     if (handshake) then
                        call mpi_send ( hs_snd, 1, mp_i4, p, p, comm, &
                                        ierr )
                     endif
                  endif
               endif

! Wait for outstanding i(r)send request to complete
               if (.not. sendd) then
                  call mpi_wait( sndids(istep-maxreqh), InStats, ierr )
               endif
            endif
!
         endif
!
      enddo

! local copy to send buffer
      if (sndlths(comm_pid) .ne. 0) then

         offset_v = sdispls(comm_pid)
         do j = 1, send_bl(comm_pid+1)%nparcels
            do i = 1, send_bl(comm_pid+1)%blocksizes(j)
               ga_r8_s(offset_v+i) = a1in(send_bl(comm_pid+1)%displacements(j)+i)
            enddo
            offset_v = offset_v + send_bl(comm_pid+1)%blocksizes(j)
         enddo

         ! pack second array (if it exists)
         if ( present(a2in) ) then
            offset_v = sdispls(comm_pid) + arr_sndlths(comm_pid)
            do j = 1, send_bl(comm_pid+1)%nparcels
               do i = 1, send_bl(comm_pid+1)%blocksizes(j)
                  ga_r8_s(offset_v+i) = a2in(send_bl(comm_pid+1)%displacements(j)+i)
               enddo
               offset_v = offset_v + send_bl(comm_pid+1)%blocksizes(j)
            enddo
         endif

         if (.not. alltoall) then
            ga_r8_r(rdispls(comm_pid)+1:rdispls(comm_pid)+rcvlths(comm_pid)) = &
               ga_r8_s(sdispls(comm_pid)+1:sdispls(comm_pid)+sndlths(comm_pid))
         endif

      endif

      if (alltoall) then
         call mpi_alltoallv (ga_r8_s, sndlths, sdispls, mp_r8, &
                             ga_r8_r, rcvlths, rdispls, mp_r8, &
                             comm, ierror)
      endif

! local copy from receive buffer
      if (rcvlths(comm_pid) .ne. 0) then

         offset_v = rdispls(comm_pid)
         do j = 1, recv_bl(comm_pid+1)%Nparcels
            do i = 1, recv_bl(comm_pid+1)%blocksizes(j)
               a1out(recv_bl(comm_pid+1)%displacements(j)+i) = ga_r8_r(offset_v+i)
            enddo
            offset_v = offset_v + recv_bl(comm_pid+1)%blocksizes(j)
         enddo

         ! scatter data for second array (if it exists)
         if ( present(a2out) ) then
            offset_v = rdispls(comm_pid) + arr_rcvlths(comm_pid)
            do j = 1, recv_bl(comm_pid+1)%Nparcels
               do i = 1, recv_bl(comm_pid+1)%blocksizes(j)
                  a2out(recv_bl(comm_pid+1)%displacements(j)+i) = ga_r8_r(offset_v+i)
               enddo
               offset_v = offset_v + recv_bl(comm_pid+1)%blocksizes(j)
            enddo
         endif

      endif

! scatter data from global receive buffer to final destination
      do istep=1,steps
         p = swapids(istep)

         if (.not. alltoall) then
            if (istep > steps-maxreqh) then
               call mpi_wait( rcvids(istep), OutStats, ierr )
            endif
         endif

         if (rcvlths(p) .ne. 0) then

            offset_v = rdispls(p)
            do j = 1, recv_bl(p+1)%Nparcels
               do i = 1, recv_bl(p+1)%blocksizes(j)
                  a1out(recv_bl(p+1)%displacements(j)+i) = ga_r8_r(offset_v+i)
               enddo
               offset_v = offset_v + recv_bl(p+1)%blocksizes(j)
            enddo

            ! scatter data for second array (if it exists)
            if ( present(a2out) ) then

               offset_v = rdispls(p) + arr_rcvlths(p)
               do j = 1, recv_bl(p+1)%Nparcels
                  do i = 1, recv_bl(p+1)%blocksizes(j)
                     a2out(recv_bl(p+1)%displacements(j)+i) = ga_r8_r(offset_v+i)
                  enddo
                  offset_v = offset_v + recv_bl(p+1)%blocksizes(j)
               enddo

            endif

         endif
      enddo

! Wait for any outstanding send requests to complete.
      if (.not. alltoall .and. .not. sendd) then
         call mpi_waitall( maxreqh, sndids(steps-maxreqh+1), InStats, ierr )
      endif

! clean-up
! make used portion of storage window available for reuse
      r8_win%ncall_s = r8_win%ncall_s - 1

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_swapirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_swapirr_i4 --- Write i4 contiguous parcels to global array
!                             using XOR swap ordering
!
! !INTERFACE:
      subroutine mp_swapirr_i4 ( comm, send_bl, recv_bl, a1in, a1out, &
                                 a2in, a2out, sw_handshake, sw_maxreq, &
                                 sw_alltoall, sw_send )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm                     ! communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer(i4), intent(in) :: a1in(*)               ! input array
      integer(i4), optional, intent(in) :: a2in(*)     ! second input array
      logical, optional, intent(in) :: sw_handshake    ! use flow control and 
                                                       !  ready send
      integer, optional, intent(in) :: sw_maxreq       ! maximum number of outstanding
                                                       !  MPI requests
      logical, optional, intent(in) :: sw_alltoall     ! use mpi_alltoall
      logical, optional, intent(in) :: sw_send         ! use mpi_send instead of isend

! !OUTPUT PARAMETERS:
      integer(i4), intent(out) :: a1out(*)             ! output array
      integer(i4), optional, intent(out) :: a2out(*)   ! second output array
!
! !DESCRIPTION:
!     XOR-ordered version of all-to-all communication
!
! WARNING: mod_comm parameter max_irr might need to be set larger than expected
!          when swapping two variables; specifically, max_irr must be at least
!          as large as the incoming i4_win%ncall_s + the number of variables to
!          be swapped
!
! !REVISION HISTORY: 
!    08.06.30   Worley      original: derived from mp_sendirr, but using 
!                            swapm logic and XOR swap order 
!    08.08.22   Worley      removed swapm; reimplemented with native MPI,
!                            added flow control/ready send option and maxreq
!                            throttling, added alltoall option
!    09.10.07   Worley      eliminated mpi_recv from handshake logic
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: i, j, p, istep, num_s, num_r
      integer :: comm_pid, comm_size, steps, ierr
      integer :: ipe, offset_s, offset_r, offset_0, unitsize, onetwo

      integer :: arr_sndlths(0:numpro-1), arr_rcvlths(0:numpro-1)
      integer :: sndlths(0:numpro-1), sdispls(0:numpro-1)
      integer :: rcvlths(0:numpro-1), rdispls(0:numpro-1)
      integer :: swapids(numpro) 
      integer :: sndids(numpro)  ! nonblocking MPI send request ids
      integer :: rcvids(numpro)  ! nonblocking MPI recv request ids
      integer :: hs_snd, hs_rcv(numpro)! handshake variables (send/receive)
      integer :: hs_rcvids(numpro) ! nonblocking MPI handshake recv request ids
      integer :: InStats(numpro*MPI_STATUS_SIZE)
      integer :: OutStats(numpro*MPI_STATUS_SIZE)

      integer :: offset_v

      integer :: rstep

      integer :: maxreq, maxreqh
      logical :: handshake, alltoall, sendd

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

!     num_s = 0 if this process is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this process is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      if ( present(a2in) .and. (.not. present(a2out)) ) then
         write(iulog,*) "Fatal mp_swapirr_i4: a2in specified, but a2out missing - exiting"
         stop
      endif

      if ( (.not. present(a2in)) .and. present(a2out)) then
         write(iulog,*) "Fatal mp_swapirr_i4: a2out specified, but a2in missing - exiting"
         stop
      endif

      if ( present(sw_handshake) ) then
         handshake = sw_handshake
         hs_snd = 1
      else
         handshake = .false.
      endif

      if ( present(sw_alltoall) ) then
         alltoall = sw_alltoall
      else
         alltoall = .false.
      endif

      if ( present(sw_send) ) then
         sendd = sw_send
      else
         sendd = .false.
      endif

      onetwo = 1
      if (present(a2in)) onetwo = 2
      unitsize = i4_win%size/max_irr

! advance to unused portion of storage window
      i4_win%ncall_s = i4_win%ncall_s + 1

      if (i4_win%ncall_s .gt. max_irr-onetwo+1) then
         write(iulog,*) "mp_swapirr_i4: insufficient window storage - exiting"
         write(iulog,*) "i4_win%ncall_s max_irr = ", i4_win%ncall_s, max_irr
         stop
      endif

! calculate send lengths and displacements
      offset_0 = (i4_win%ncall_s-1)*unitsize
      offset_s = offset_0
      sndlths(:) = 0
      sdispls(:) = 0
      arr_sndlths(:) = 0
      do ipe=1, num_s
         sndlths(ipe-1) = send_bl(ipe)%Tot_Size
         sdispls(ipe-1) = offset_s
         if (sndlths(ipe-1) .ne. 0) then

            ! pack first array
            offset_s = offset_s + sndlths(ipe-1)
            if (offset_s-offset_0 .gt. onetwo*unitsize) then
              write(iulog,*) "Fatal mp_swapirr_i4: send window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                             ipe, unitsize, offset_s, offset_0
              stop
            endif

            arr_sndlths(ipe-1) = sndlths(ipe-1)

            ! calculate for second array (if it exists)
            if ( present(a2in) ) then

               offset_s = offset_s + sndlths(ipe-1)
               if (offset_s-offset_0 .gt. onetwo*unitsize) then
                 write(iulog,*) "Fatal mp_swapirr_i4: send window out of space - exiting"
                 write(iulog,*) '2 comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                                ipe, unitsize, offset_s, offset_0
                 stop
               endif

               sndlths(ipe-1) = sndlths(ipe-1) + arr_sndlths(ipe-1)

            endif

         endif
      enddo

! calculate receive lengths and displacements
      offset_r = offset_0
      rcvlths(:) = 0
      rdispls(:) = 0
      arr_rcvlths(:) = 0
      do ipe=1, num_r
         rcvlths(ipe-1) = recv_bl(ipe)%Tot_Size
         rdispls(ipe-1) = offset_r
         if (rcvlths(ipe-1) .ne. 0) then

            offset_r = offset_r + rcvlths(ipe-1)
            if (onetwo*unitsize < offset_r-offset_0) then
              write(iulog,*) "Fatal mp_swapirr_i4: receive window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif

            arr_rcvlths(ipe-1) = rcvlths(ipe-1)

            ! compute for second array (if it exists)
            if ( present(a2out) ) then

               offset_r = offset_r + rcvlths(ipe-1)
               if (onetwo*unitsize < offset_r-offset_0) then
                 write(iulog,*) "Fatal mp_swapirr_i4: receive window out of space - exiting"
                 write(iulog,*) '2 comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                          ipe, unitsize, offset_r, offset_0
                 stop
               endif

               rcvlths(ipe-1) = rcvlths(ipe-1) + arr_rcvlths(ipe-1)

            endif

         endif
      enddo

! Calculate swap partners and number of steps in point-to-point
! implementations of alltoall algorithm.
      steps = 0
      do ipe=1,ceil2(comm_size)-1
         p = pair(comm_size,ipe,comm_pid)
         if (p >= 0) then
            if (sndlths(p) > 0 .or. rcvlths(p) > 0) then
               steps = steps + 1
               swapids(steps) = p
            end if
         end if
      end do

      if (.not. alltoall) then

         sndids(1:steps) = MPI_REQUEST_NULL
         rcvids(1:steps) = MPI_REQUEST_NULL

         if (steps .eq. 0) then
            maxreq  = 0
            maxreqh = 0
         elseif (steps .eq. 1) then
            maxreq  = 1
            maxreqh = 1
         else
            if ( present(sw_maxreq) ) then
               if ((sw_maxreq .le. steps) .and. (sw_maxreq .ge. 0)) then
                  maxreq  = sw_maxreq
                  if (maxreq > 1) then
                     maxreqh = maxreq/2
                  else
                     maxreq  = 2
                     maxreqh = 1
                  endif
               else
                  maxreq  = steps
                  maxreqh = steps
               endif
            else
               maxreq  = steps
               maxreqh = steps
            endif
         endif

! Post initial handshake receive requests
         if (handshake) then
            do istep=1,maxreq
               p = swapids(istep)
               if (sndlths(p) > 0) then
                  call mpi_irecv  ( hs_rcv(istep), 1, mp_i4, p, comm_pid, comm, &
                                    hs_rcvids(istep), ierr )
               endif
            enddo
         endif

! Post initial receive requests
         do istep=1,maxreq
            p = swapids(istep)
            if (rcvlths(p) > 0) then
               offset_r = rdispls(p)+1
               call mpi_irecv ( ga_i4_r(offset_r), rcvlths(p), mp_i4, &
                                p, p, comm, rcvids(istep), ierr )
               if (handshake) then
                  call mpi_send( hs_snd, 1, mp_i4, p, p, comm, &
                                 ierr )
               endif
            endif
         enddo
         rstep = maxreq
!
      endif

! gather data into global send buffer
      do istep=1,steps
         p = swapids(istep)

         if (sndlths(p) .ne. 0) then
            offset_v = sdispls(p)
            do j = 1, send_bl(p+1)%nparcels
               do i = 1, send_bl(p+1)%blocksizes(j)
                  ga_i4_s(offset_v+i) = a1in(send_bl(p+1)%displacements(j)+i)
               enddo
               offset_v = offset_v + send_bl(p+1)%blocksizes(j)
            enddo

            ! pack second array (if it exists)
            if ( present(a2in) ) then
               offset_v = sdispls(p) + arr_sndlths(p)
               do j = 1, send_bl(p+1)%nparcels
                  do i = 1, send_bl(p+1)%blocksizes(j)
                     ga_i4_s(offset_v+i) = a2in(send_bl(p+1)%displacements(j)+i)
                  enddo
                  offset_v = offset_v + send_bl(p+1)%blocksizes(j)
               enddo
            endif

         endif

         if (.not. alltoall) then

! Submit new i(r)send request
            offset_s = sdispls(p)+1
            if (sndlths(p) > 0) then
               if (handshake) then
                  call mpi_wait( hs_rcvids(istep), MPI_STATUS_IGNORE, ierr )
                  if (sendd) then
                     call mpi_rsend( ga_i4_s(offset_s), sndlths(p), mp_i4, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_irsend( ga_i4_s(offset_s), sndlths(p), mp_i4, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               else
                  if (sendd) then
                     call mpi_send ( ga_i4_s(offset_s), sndlths(p), mp_i4, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_isend ( ga_i4_s(offset_s), sndlths(p), mp_i4, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               endif
            endif

            if (istep > maxreqh) then
! Wait for oldest irecv request to complete
               call mpi_wait( rcvids(istep-maxreqh), OutStats, ierr )

               if (rstep < steps) then
                  rstep = rstep + 1
                  p = swapids(rstep)

! Submit a new handshake irecv request
                  if (handshake) then
                     if (sndlths(p) > 0) then
                        call mpi_irecv( hs_rcv(rstep), 1, mp_i4, p, comm_pid, comm, &
                                        hs_rcvids(rstep), ierr )
                     endif
                  endif

! Submit a new irecv request
                  if (rcvlths(p) > 0) then
                     offset_r = rdispls(p)+1
                     call mpi_irecv( ga_i4_r(offset_r), rcvlths(p), mp_i4, &
                                     p, p, comm, rcvids(rstep), ierr )
                     if (handshake) then
                        call mpi_send ( hs_snd, 1, mp_i4, p, p, comm, &
                                        ierr )
                     endif
                  endif
               endif

! Wait for outstanding i(r)send request to complete
               if (.not. sendd) then
                  call mpi_wait( sndids(istep-maxreqh), InStats, ierr )
               endif
            endif
!
         endif
!
      enddo

! local copy to send buffer
      if (sndlths(comm_pid) .ne. 0) then

         offset_v = sdispls(comm_pid)
         do j = 1, send_bl(comm_pid+1)%nparcels
            do i = 1, send_bl(comm_pid+1)%blocksizes(j)
               ga_i4_s(offset_v+i) = a1in(send_bl(comm_pid+1)%displacements(j)+i)
            enddo
            offset_v = offset_v + send_bl(comm_pid+1)%blocksizes(j)
         enddo

         ! pack second array (if it exists)
         if ( present(a2in) ) then
            offset_v = sdispls(comm_pid) + arr_sndlths(comm_pid)
            do j = 1, send_bl(comm_pid+1)%nparcels
               do i = 1, send_bl(comm_pid+1)%blocksizes(j)
                  ga_i4_s(offset_v+i) = a2in(send_bl(comm_pid+1)%displacements(j)+i)
               enddo
               offset_v = offset_v + send_bl(comm_pid+1)%blocksizes(j)
            enddo
         endif

         if (.not. alltoall) then
            ga_i4_r(rdispls(comm_pid)+1:rdispls(comm_pid)+rcvlths(comm_pid)) = &
               ga_i4_s(sdispls(comm_pid)+1:sdispls(comm_pid)+sndlths(comm_pid))
         endif

      endif

      if (alltoall) then
         call mpi_alltoallv (ga_i4_s, sndlths, sdispls, mp_i4, &
                             ga_i4_r, rcvlths, rdispls, mp_i4, &
                             comm, ierror)
      endif

! local copy from receive buffer
      if (rcvlths(comm_pid) .ne. 0) then

         offset_v = rdispls(comm_pid)
         do j = 1, recv_bl(comm_pid+1)%Nparcels
            do i = 1, recv_bl(comm_pid+1)%blocksizes(j)
               a1out(recv_bl(comm_pid+1)%displacements(j)+i) = ga_i4_r(offset_v+i)
            enddo
            offset_v = offset_v + recv_bl(comm_pid+1)%blocksizes(j)
         enddo

         ! scatter data for second array (if it exists)
         if ( present(a2out) ) then
            offset_v = rdispls(comm_pid) + arr_rcvlths(comm_pid)
            do j = 1, recv_bl(comm_pid+1)%Nparcels
               do i = 1, recv_bl(comm_pid+1)%blocksizes(j)
                  a2out(recv_bl(comm_pid+1)%displacements(j)+i) = ga_i4_r(offset_v+i)
               enddo
               offset_v = offset_v + recv_bl(comm_pid+1)%blocksizes(j)
            enddo
         endif

      endif

! scatter data from global receive buffer to final destination
      do istep=1,steps
         p = swapids(istep)

         if (.not. alltoall) then
            if (istep > steps-maxreqh) then
               call mpi_wait( rcvids(istep), OutStats, ierr )
            endif
         endif

         if (rcvlths(p) .ne. 0) then

            offset_v = rdispls(p)
            do j = 1, recv_bl(p+1)%Nparcels
               do i = 1, recv_bl(p+1)%blocksizes(j)
                  a1out(recv_bl(p+1)%displacements(j)+i) = ga_i4_r(offset_v+i)
               enddo
               offset_v = offset_v + recv_bl(p+1)%blocksizes(j)
            enddo

            ! scatter data for second array (if it exists)
            if ( present(a2out) ) then

               offset_v = rdispls(p) + arr_rcvlths(p)
               do j = 1, recv_bl(p+1)%Nparcels
                  do i = 1, recv_bl(p+1)%blocksizes(j)
                     a2out(recv_bl(p+1)%displacements(j)+i) = ga_i4_r(offset_v+i)
                  enddo
                  offset_v = offset_v + recv_bl(p+1)%blocksizes(j)
               enddo

            endif

         endif
      enddo

! Wait for any outstanding send requests to complete.
      if (.not. alltoall .and. .not. sendd) then
         call mpi_waitall( maxreqh, sndids(steps-maxreqh+1), InStats, ierr )
      endif

! clean-up
! make used portion of storage window available for reuse
      i4_win%ncall_s = i4_win%ncall_s - 1

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_swapirr_i4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: pair 
!
! !INTERFACE:
      integer function pair(np,p,k)
!
! !INPUT PARAMETERS:
      integer :: np
      integer :: p
      integer :: k
! !DESCRIPTION:
!
!     Bitwise XOR of arguments p and k, if less than upper bound np
!
! !REVISION HISTORY: 
!    2008.08.21   Worley         Imported from spmdutils
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer q
!
      q = ieor(p,k)
      if ( q > np-1 ) then
         pair = -1
      else
         pair = q
      endif

      return

!EOC
      end function pair
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: ceil2
!
! !INTERFACE:
     integer function ceil2(n)
!
! !INPUT PARAMETERS:
     integer :: n
! !DESCRIPTION:
!
!     Smallest power of 2 greater than or equal to the argument
!
! !REVISION HISTORY: 
!    2008.08.21   Worley         Imported from spmdutils
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
     integer p

     p=1
     do while ( p < n )
        p=p*2
     enddo
     ceil2=p

     return
!EOC
     end function ceil2
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
# if defined( MOD_ASSUMED_SIZE )
!BOP
! !ROUTINE: mp_sendtrirr --- Initiate communication of contiguous tracer parcels
!
! !INTERFACE:
      subroutine mp_sendtrirr ( comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq,  &
                                ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts,  &
                                ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr,  &
                                modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer, intent(in)           :: mbeg            ! initial tracer index
      integer, intent(in)           :: mend            ! final tracer index
      integer, intent(in)           :: mq              ! total tracer indices
      integer, intent(in)           :: ifirsts         ! first I index of source
      integer, intent(in)           :: ilasts          ! last I index of source
      integer, intent(in)           :: jfirsts         ! first j index of source
      integer, intent(in)           :: jlasts          ! last j index of source
      integer, intent(in)           :: kfirsts         ! first k index of source
      integer, intent(in)           :: klasts          ! last k index of source
      integer, intent(in)           :: ifirstr         ! first I index of target
      integer, intent(in)           :: ilastr          ! last I index of target
      integer, intent(in)           :: jfirstr         ! first j index of target
      integer, intent(in)           :: jlastr          ! last j index of target
      integer, intent(in)           :: kfirstr         ! first k index of target
      integer, intent(in)           :: klastr          ! last k index of target
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests
      real(r8), intent(in) :: qin(*) ! input tracer array

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(*) ! output tracer array
!
! !DESCRIPTION:
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications here and placing
!     wait points in mp_recvtrirr; if 1, call swap routine with p2p messages; if 2, call swap
!     routine with a2a messages. 
!     Modc(2): if 1, then apply handshaking (don't send until corresponding receive is posted)
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!     Modc(4): maximum number of outstanding requests (applies to swap routines only)
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!    09.10.07   Worley      eliminated mpi_recv from handshake logic
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nparcels, offset_s, offset_r, ierr, mod_method
      integer p, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize, unitsize, offset_0
      integer i, j, send_tag, recv_tag, num_s, num_r, m
      integer :: offset_v (Max_Nparcels)
      integer :: hs_snd, hs_rcv(numpro), hs_rcvids(numpro)
      integer ipe2, ceil2num
      integer numtr, numtrm
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_pid
      integer ijks, ijkr, ij


#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .true.
         maxreq_local = -1
      endif

! Do not call mp_swaptrirr unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

    if (sw_local .gt. 0) then
         sw_alltoall = (sw_local .eq. 2)
         call mp_swaptrirr(comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq, &
                           ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts, &
                           ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr, &
                           sw_handshake=hs_local, sw_maxreq=maxreq_local,     &
                           sw_alltoall=sw_alltoall, sw_send=send_local)
    else

      call MPI_COMM_RANK (comm, comm_pid, ierr)

      hs_snd = 1
      ceil2num = ceil2(numpro)

      numtrm = mend - mbeg
      numtr = numtrm + 1

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      r8_win%ncall_s = r8_win%ncall_s + 1

      ijks =(klasts-kfirsts+1)*(jlasts-jfirsts+1)*(ilasts-ifirsts+1)
      ijkr =(klastr-kfirstr+1)*(jlastr-jfirstr+1)*(ilastr-ifirstr+1)

     if (mod_method .gt. 0) then
!
! mpi derived types
      if (r8_win%ncall_s .gt. MaxTrf-numtrm) then
         write(iulog,*) "mp_sendtrirr: derived type handle count exceeded - exiting"
         write(iulog,*) "r8_win%ncall_s MaxTrf = ", r8_win%ncall_s, MaxTrf
         stop
      endif
!
! MPI: Irecv over all processes
!
      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               if (ipe-1 /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, ipe-1, comm_pid, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      OutHandle(:,r8_win%ncall_s:r8_win%ncall_s+numtrm) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_r) cycle
!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          recv_tag = ipe-1 + modcam_tagoffset
          do m = mbeg, mend
             call mpi_irecv( qout((m-1)*ijkr+1), 1, recv_bl(ipe)%type, ipe-1, recv_tag,   &
                             comm, OutHandle(ipe,r8_win%ncall_s+m-mbeg), ierr )
          enddo
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
               call MPI_SEND ( hs_snd, 1, mp_i4, ipe-1, ipe-1, comm, ierr )
          endif
        endif
      enddo

!
! MPI: Isend/Send over all processes; use risend/rsend with hs
!
      InHandle(:,r8_win%ncall_s:r8_win%ncall_s+numtrm) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_s) cycle

!
! Send the individual buffers with non-blocking sends
!
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          send_tag = comm_pid + modcam_tagoffset
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
                call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
             if (send_local) then
                do m = mbeg, mend
                   call mpi_rsend( qin((m-1)*ijks+1), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, ierr )
                enddo
             else
                do m = mbeg, mend
                   call mpi_irsend( qin((m-1)*ijks+1), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                    comm, InHandle(ipe,r8_win%ncall_s), ierr )
                enddo
             endif
          else
             if (send_local) then
                do m = mbeg, mend
                   call mpi_send( qin((m-1)*ijks+1), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                  comm, ierr )
                enddo
             else
                do m = mbeg, mend
                   call mpi_isend( qin((m-1)*ijks+1), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, InHandle(ipe,r8_win%ncall_s), ierr )
                enddo
             endif
          endif
        endif
      enddo
     else

! temporary contiguous buffers

      if (r8_win%ncall_s .gt. max_irr-numtrm) then
         write(iulog,*) "mp_sendtrirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_s max_irr = ", r8_win%ncall_s, max_irr
         stop
      endif
      unitsize = r8_win%size/max_irr

! issue call to receive data in global receive buffer
      offset_0 = (r8_win%ncall_s-1)*unitsize
      offset_s = offset_0
      offset_r = offset_0

      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            qsize = numtr*send_bl(ipe)%Tot_Size
            if (qsize .ne. 0) then
               r8_win%dest = ipe-1
               send_tag = comm_pid + modcam_tagoffset
               if (r8_win%dest /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, r8_win%dest, send_tag, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         r8_win%size_r = numtr*recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            r8_win%src = ipe-1
            if (numtr*unitsize >= offset_r-offset_0) then
              recv_tag = r8_win%src + modcam_tagoffset
              qsize    = r8_win%size_r
              r8_win%nrecv    = r8_win%nrecv + 1
              call MPI_IRECV(ga_r8_r(r8_win%offset_r+1), qsize, mp_r8, r8_win%src, &
                             recv_tag, comm, r8_win%rqest(r8_win%nrecv), ierror)
              if (hs_local) then
                 if (r8_win%src /= comm_pid) &
                   call MPI_SEND ( hs_snd, 1, mp_i4, r8_win%src, recv_tag, comm, ierror)
              endif
            else
              write(iulog,*) "Fatal mp_sendtrirr: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif
         endif
      enddo
! gather data into global send buffer
      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_s) cycle
         qsize = numtr*send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
            r8_win%dest = ipe-1
            r8_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s-offset_0 .gt. numtr*unitsize) then
              write(iulog,*) "Fatal mp_sendtrirr: send window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_s, offset_0
              stop
            endif

            offset_v(1) = r8_win%offset_s
            do j = 2, send_bl(ipe)%nparcels
               offset_v(j) = offset_v(j-1) + send_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, send_bl(ipe)%nparcels
               do m = mbeg, mend
                  do i = 1, send_bl(ipe)%blocksizes(j)
                     ij = send_bl(ipe)%displacements(j)+i
                     ga_r8_s(send_bl(ipe)%Tot_Size*(m-mbeg)+offset_v(j)+i) = qin((m-1)*ijks+ij)
                  enddo
               enddo
            enddo

! nonblocking send
            send_tag = comm_pid + modcam_tagoffset
            r8_win%nsend = r8_win%nsend + 1
            if (hs_local) then
               if (r8_win%dest /= comm_pid) &
                  call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
               if (send_local) then
                  call MPI_RSEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_IRSEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, r8_win%sqest(r8_win%nsend), ierr)
               endif
            else
               if (send_local) then
                  call MPI_SEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_ISEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, r8_win%sqest(r8_win%nsend), ierr)
               endif
            endif
         endif
      enddo

     endif   !  mod_method

      r8_win%ncall_s = r8_win%ncall_s + numtrm

    endif   !  sw_local

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendtrirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvtrirr --- Finalize communication of contiguous tracer parcels
!
! !INTERFACE:
      subroutine mp_recvtrirr ( comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq,  &
                                ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts,  &
                                ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr,  &
                                modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer, intent(in)           :: mbeg            ! initial tracer index
      integer, intent(in)           :: mend            ! final tracer index
      integer, intent(in)           :: mq              ! total tracer indices
      integer, intent(in)           :: ifirsts         ! first I index of source
      integer, intent(in)           :: ilasts          ! last I index of source
      integer, intent(in)           :: jfirsts         ! first j index of source
      integer, intent(in)           :: jlasts          ! last j index of source
      integer, intent(in)           :: kfirsts         ! first k index of source
      integer, intent(in)           :: klasts          ! last k index of source
      integer, intent(in)           :: ifirstr         ! first I index of target
      integer, intent(in)           :: ilastr          ! last I index of target
      integer, intent(in)           :: jfirstr         ! first j index of target
      integer, intent(in)           :: jlastr          ! last j index of target
      integer, intent(in)           :: kfirstr         ! first k index of target
      integer, intent(in)           :: klastr          ! last k index of target
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests
      real(r8), intent(in) :: qin(*) ! input tracer array
! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(*) ! output tracer array
!
! !DESCRIPTION:
!     Complete transfer of a generalized region initiated by {\tt mp\_sendtrirr}.
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications in mp_sendtrirr and
!     placing wait points here; otherwise don't do anything - mp_swaptrirr is called from mp_sendirr.
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, offset_r, mod_method
      integer unitsize, offset_0
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer i, j, num_r, num_s, m
      integer :: offset_v (Max_Nparcels)
      integer ipe2, ceil2num
      integer numtr, numtrm
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_size, comm_pid
      integer ijks, ijkr, ij

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .true.
         maxreq_local = -1
      endif

! Do not call mp_swaptrirr (hence return) unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

! Return if swap_irr
      if (sw_local .gt. 0) return

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

      ceil2num = ceil2(numpro)

      numtrm = mend - mbeg
      numtr = numtrm + 1

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      r8_win%ncall_r = r8_win%ncall_r + 1

      ijks =(klasts-kfirsts+1)*(jlasts-jfirsts+1)*(ilasts-ifirsts+1)
      ijkr =(klastr-kfirstr+1)*(jlastr-jfirstr+1)*(ilastr-ifirstr+1)

    if (mod_method .gt. 0) then

! mpi derived types
      if (r8_win%ncall_r .gt. MaxTrf-numtrm) then
         write(iulog,*) "mp_recvtrirr: derived type handle count exceeded - exiting"
         write(iulog,*) "r8_win%ncall_r MaxTrf = ", r8_win%ncall_r, MaxTrf
         stop
      endif

      if (num_s .gt. 0 .and. (.not. send_local)) then
         do m = mbeg, mend
            CALL MPI_WAITALL( comm_size, InHandle(:,r8_win%ncall_r+m-mbeg), InStats, Ierr )
         enddo
      endif
      if (num_r .gt. 0) then
         do m = mbeg, mend
            CALL MPI_WAITALL( comm_size, OutHandle(:,r8_win%ncall_r+m-mbeg), OutStats, Ierr )
         enddo
      endif

    else

! temporary contiguous buffer / global window

      if (r8_win%ncall_r .gt. max_irr-numtrm) then
         write(iulog,*) "mp_recvtrirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_r max_irr = ", r8_win%ncall_r, max_irr
         stop
      endif
      unitsize = r8_win%size/max_irr

! scatter data from global receive buffer to final destination
      offset_0 = (r8_win%ncall_r-1)*unitsize
      offset_r = offset_0

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         r8_win%size_r = numtr*recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r-offset_0 .gt. numtr*unitsize) then
              write(iulog,*) "Fatal mp_recvtrirr: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif

            r8_win%nread = r8_win%nread + 1
            call MPI_WAIT(r8_win%rqest(r8_win%nread), Status, ierr)

            offset_v(1) = r8_win%offset_r
            do j = 2, recv_bl(ipe)%Nparcels
               offset_v(j) = offset_v(j-1) + recv_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, recv_bl(ipe)%Nparcels
               do m = mbeg, mend
                  do i = 1, recv_bl(ipe)%blocksizes(j)
                     ij = recv_bl(ipe)%displacements(j)+i
                     qout((m-1)*ijkr+ij) = ga_r8_r(recv_bl(ipe)%Tot_Size*(m-mbeg)+offset_v(j)+i)
                  enddo
               enddo
            enddo

         endif
      enddo

      if ((r8_win%ncall_s == r8_win%ncall_r + numtrm) .and. (.not. send_local)) then
         call MPI_WAITALL(r8_win%nsend, r8_win%sqest, Stats, ierror)
      endif

    endif    !    mod_method .gt. 0

    r8_win%ncall_r = r8_win%ncall_r + numtrm

    if (r8_win%ncall_s == r8_win%ncall_r) then
       r8_win%nsend = 0
       r8_win%nrecv = 0
       r8_win%nread = 0
       r8_win%ncall_s = 0
       r8_win%ncall_r = 0
    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvtrirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_swaptrirr --- Write r8 contiguous parcels to global array
!                            using XOR swap ordering - for multiple tracers
!
! !INTERFACE:
      subroutine mp_swaptrirr ( comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq,  &
                                ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts,  &
                                ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr,  &
                                sw_handshake, sw_maxreq, sw_alltoall, sw_send   )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm                     ! communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer, intent(in)           :: mbeg            ! initial tracer index
      integer, intent(in)           :: mend            ! final tracer index
      integer, intent(in)           :: mq              ! total tracer indices
      integer, intent(in)           :: ifirsts         ! first I index of source
      integer, intent(in)           :: ilasts          ! last I index of source
      integer, intent(in)           :: jfirsts         ! first j index of source
      integer, intent(in)           :: jlasts          ! last j index of source
      integer, intent(in)           :: kfirsts         ! first k index of source
      integer, intent(in)           :: klasts          ! last k index of source
      integer, intent(in)           :: ifirstr         ! first I index of target
      integer, intent(in)           :: ilastr          ! last I index of target
      integer, intent(in)           :: jfirstr         ! first j index of target
      integer, intent(in)           :: jlastr          ! last j index of target
      integer, intent(in)           :: kfirstr         ! first k index of target
      integer, intent(in)           :: klastr          ! last k index of target
      logical, optional, intent(in) :: sw_handshake    ! use flow control and 
                                                       !  ready send
      integer, optional, intent(in) :: sw_maxreq       ! maximum number of outstanding
                                                       !  MPI requests
      logical, optional, intent(in) :: sw_alltoall     ! use mpi_alltoall
      logical, optional, intent(in) :: sw_send         ! use mpi_send instead of isend
      real(r8), intent(in) :: qin(*) ! input tracer array

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(*) ! output tracer array
!
! !DESCRIPTION:
!
!     XOR-ordered version of all-to-all communication
!
! WARNING: mod_comm parameter max_irr might need to be set larger than expected
!          when swapping multiple variables; specifically, max_irr must be at least
!          as large as the incoming r8_win%ncall_s + the number of variables to
!          be swapped
!
! !REVISION HISTORY: 
!    08.06.30   Worley      original: derived from mp_sendirr, but using 
!                            swapm logic and XOR swap order 
!    08.08.22   Worley      removed swapm; reimplemented with native MPI,
!                            added flow control/ready send option and maxreq
!                            throttling, added alltoall option
!    09.10.07   Worley      eliminated mpi_recv from handshake logic
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: i, j, p, istep, num_s, num_r
      integer :: comm_pid, comm_size, steps, ierr
      integer :: ipe, offset_s, offset_r, offset_0, unitsize

      integer :: sndlths(0:numpro-1), sdispls(0:numpro-1)
      integer :: rcvlths(0:numpro-1), rdispls(0:numpro-1)
      integer :: swapids(numpro) 
      integer :: sndids(numpro)  ! nonblocking MPI send request ids
      integer :: rcvids(numpro)  ! nonblocking MPI recv request ids
      integer :: hs_snd, hs_rcv(numpro)! handshake variables (send/receive)
      integer :: hs_rcvids(numpro) ! nonblocking MPI handshake recv request ids
      integer :: InStats(numpro*MPI_STATUS_SIZE)
      integer :: OutStats(numpro*MPI_STATUS_SIZE)

      integer :: offset_v

      integer :: rstep

      integer :: maxreq, maxreqh
      logical :: handshake, alltoall, sendd
      integer ::  numtr, numtrm, m
      integer ijks, ijkr, ij

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

!     num_s = 0 if this process is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this process is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      if ( present(sw_handshake) ) then
         handshake = sw_handshake
         hs_snd = 1
      else
         handshake = .false.
      endif

      if ( present(sw_alltoall) ) then
         alltoall = sw_alltoall
      else
         alltoall = .false.
      endif

      if ( present(sw_send) ) then
         sendd = sw_send
      else
         sendd = .false.
      endif

      numtrm = mend - mbeg
      numtr = numtrm + 1

      ijks =(klasts-kfirsts+1)*(jlasts-jfirsts+1)*(ilasts-ifirsts+1)
      ijkr =(klastr-kfirstr+1)*(jlastr-jfirstr+1)*(ilastr-ifirstr+1)

      unitsize = r8_win%size/max_irr

! advance to unused portion of storage window
      r8_win%ncall_s = r8_win%ncall_s + 1

      if (r8_win%ncall_s .gt. max_irr-numtrm) then
         write(iulog,*) "mp_swaptrirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_s max_irr = ", r8_win%ncall_s, max_irr
         stop
      endif

! calculate send lengths and displacements
      offset_0 = (r8_win%ncall_s-1)*unitsize
      offset_s = offset_0
      sndlths(:) = 0
      sdispls(:) = 0
      do ipe=1, num_s
         sndlths(ipe-1) = numtr*send_bl(ipe)%Tot_Size
         sdispls(ipe-1) = offset_s
         if (sndlths(ipe-1) .ne. 0) then

            offset_s = offset_s + sndlths(ipe-1)
            if (offset_s-offset_0 .gt. numtr*unitsize) then
              write(iulog,*) "Fatal mp_swaptrirr: send window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                             ipe, unitsize, offset_s, offset_0
              stop
            endif
         endif
      enddo

! calculate receive lengths and displacements
      offset_r = offset_0
      rcvlths(:) = 0
      rdispls(:) = 0
      do ipe=1, num_r
         rcvlths(ipe-1) = numtr*recv_bl(ipe)%Tot_Size
         rdispls(ipe-1) = offset_r
         if (rcvlths(ipe-1) .ne. 0) then

            offset_r = offset_r + rcvlths(ipe-1)
            if (numtr*unitsize < offset_r-offset_0) then
              write(iulog,*) "Fatal mp_swaptrirr: receive window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif
         endif
      enddo

! Calculate swap partners and number of steps in point-to-point
! implementations of alltoall algorithm.
      steps = 0
      do ipe=1,ceil2(comm_size)-1
         p = pair(comm_size,ipe,comm_pid)
         if (p >= 0) then
            if (sndlths(p) > 0 .or. rcvlths(p) > 0) then
               steps = steps + 1
               swapids(steps) = p
            end if
         end if
      end do

      if (.not. alltoall) then

         sndids(1:steps) = MPI_REQUEST_NULL
         rcvids(1:steps) = MPI_REQUEST_NULL

         if (steps .eq. 0) then
            maxreq  = 0
            maxreqh = 0
         elseif (steps .eq. 1) then
            maxreq  = 1
            maxreqh = 1
         else
            if ( present(sw_maxreq) ) then
               if ((sw_maxreq .le. steps) .and. (sw_maxreq .ge. 0)) then
                  maxreq  = sw_maxreq
                  if (maxreq > 1) then
                     maxreqh = maxreq/2
                  else
                     maxreq  = 2
                     maxreqh = 1
                  endif
               else
                  maxreq  = steps
                  maxreqh = steps
               endif
            else
               maxreq  = steps
               maxreqh = steps
            endif
         endif

! Post initial handshake receive requests
         if (handshake) then
            do istep=1,maxreq
               p = swapids(istep)
               if (sndlths(p) > 0) then
                  call mpi_irecv  ( hs_rcv(istep), 1, mp_i4, p, comm_pid, comm, &
                                    hs_rcvids(istep), ierr )
               endif
            enddo
         endif

! Post initial receive requests
         do istep=1,maxreq
            p = swapids(istep)
            if (rcvlths(p) > 0) then
               offset_r = rdispls(p)+1
               call mpi_irecv ( ga_r8_r(offset_r), rcvlths(p), mp_r8, &
                                p, p, comm, rcvids(istep), ierr )
               if (handshake) then
                  call mpi_send( hs_snd, 1, mp_i4, p, p, comm, &
                                 ierr )
               endif
            endif
         enddo
         rstep = maxreq
!
      endif

! gather data into global send buffer
      do istep=1,steps
         p = swapids(istep)

         if (sndlths(p) .ne. 0) then
            offset_v = sdispls(p)
            do j = 1, send_bl(p+1)%nparcels
               do m = mbeg, mend
                  do i = 1, send_bl(p+1)%blocksizes(j)
                     ij = send_bl(p+1)%displacements(j)+i
                     ga_r8_s(send_bl(p+1)%Tot_Size*(m-mbeg)+offset_v+i) = qin((m-1)*ijks+ij)
                  enddo
               enddo
               offset_v = offset_v + send_bl(p+1)%blocksizes(j)
            enddo
         endif

         if (.not. alltoall) then

! Submit new i(r)send request
            offset_s = sdispls(p)+1
            if (sndlths(p) > 0) then
               if (handshake) then
                  call mpi_wait( hs_rcvids(istep), MPI_STATUS_IGNORE, ierr )
                  if (sendd) then
                     call mpi_rsend( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_irsend( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               else
                  if (sendd) then
                     call mpi_send ( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_isend ( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               endif
            endif

            if (istep > maxreqh) then
! Wait for oldest irecv request to complete
               call mpi_wait( rcvids(istep-maxreqh), OutStats, ierr )

               if (rstep < steps) then
                  rstep = rstep + 1
                  p = swapids(rstep)

! Submit a new handshake irecv request
                  if (handshake) then
                     if (sndlths(p) > 0) then
                        call mpi_irecv( hs_rcv(rstep), 1, mp_i4, p, comm_pid, comm, &
                                        hs_rcvids(rstep), ierr )
                     endif
                  endif

! Submit a new irecv request
                  if (rcvlths(p) > 0) then
                     offset_r = rdispls(p)+1
                     call mpi_irecv( ga_r8_r(offset_r), rcvlths(p), mp_r8, &
                                     p, p, comm, rcvids(rstep), ierr )
                     if (handshake) then
                        call mpi_send ( hs_snd, 1, mp_i4, p, p, comm, &
                                        ierr )
                     endif
                  endif
               endif

! Wait for outstanding i(r)send request to complete
               if (.not. sendd) then
                  call mpi_wait( sndids(istep-maxreqh), InStats, ierr )
               endif
            endif
!
         endif
!
      enddo

! local copy to send buffer
      if (sndlths(comm_pid) .ne. 0) then

         offset_v = sdispls(comm_pid)
         do j = 1, send_bl(comm_pid+1)%nparcels
            do m = mbeg, mend
               do i = 1, send_bl(comm_pid+1)%blocksizes(j)
                  ij = send_bl(comm_pid+1)%displacements(j)+i
                  ga_r8_s(send_bl(comm_pid+1)%Tot_Size*(m-mbeg)+offset_v+i) = qin((m-1)*ijks+ij)
               enddo
            enddo
            offset_v = offset_v + send_bl(comm_pid+1)%blocksizes(j)
         enddo

         if (.not. alltoall) then
            ga_r8_r(rdispls(comm_pid)+1:rdispls(comm_pid)+rcvlths(comm_pid)) = &
               ga_r8_s(sdispls(comm_pid)+1:sdispls(comm_pid)+sndlths(comm_pid))
         endif

      endif

      if (alltoall) then
         call mpi_alltoallv (ga_r8_s, sndlths, sdispls, mp_r8, &
                             ga_r8_r, rcvlths, rdispls, mp_r8, &
                             comm, ierror)
      endif

! local copy from receive buffer
      if (rcvlths(comm_pid) .ne. 0) then

         offset_v = rdispls(comm_pid)
         do j = 1, recv_bl(comm_pid+1)%Nparcels
            do m = mbeg, mend
               do i = 1, recv_bl(comm_pid+1)%blocksizes(j)
                  ij = recv_bl(comm_pid+1)%displacements(j)+i
                  qout((m-1)*ijkr+ij) = ga_r8_r(recv_bl(comm_pid+1)%Tot_Size*(m-mbeg)+offset_v+i)
               enddo
            enddo
            offset_v = offset_v + recv_bl(comm_pid+1)%blocksizes(j)
         enddo

      endif

! scatter data from global receive buffer to final destination
      do istep=1,steps
         p = swapids(istep)

         if (.not. alltoall) then
            if (istep > steps-maxreqh) then
               call mpi_wait( rcvids(istep), OutStats, ierr )
            endif
         endif

         if (rcvlths(p) .ne. 0) then

            offset_v = rdispls(p)
            do j = 1, recv_bl(p+1)%Nparcels
               do m = mbeg, mend
                  do i = 1, recv_bl(p+1)%blocksizes(j)
                     ij = recv_bl(p+1)%displacements(j)+i
                     qout((m-1)*ijkr+ij) = ga_r8_r(recv_bl(p+1)%Tot_Size*(m-mbeg)+offset_v+i)
                  enddo
               enddo
               offset_v = offset_v + recv_bl(p+1)%blocksizes(j)
            enddo

         endif
      enddo

! Wait for any outstanding send requests to complete.
      if (.not. alltoall .and. .not. sendd) then
         call mpi_waitall( maxreqh, sndids(steps-maxreqh+1), InStats, ierr )
      endif

! clean-up
! make used portion of storage window available for reuse
      r8_win%ncall_s = r8_win%ncall_s - 1

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_swaptrirr
# endif
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
# if defined( MOD_SPECIFIED_SHAPE )
!BOP
! !ROUTINE: mp_sendtrirr --- Initiate communication of contiguous tracer parcels
!
! !INTERFACE:
      subroutine mp_sendtrirr ( comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq,  &
                                ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts,  &
                                ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr,  &
                                modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer, intent(in)           :: mbeg            ! initial tracer index
      integer, intent(in)           :: mend            ! final tracer index
      integer, intent(in)           :: mq              ! total tracer indices
      integer, intent(in)           :: ifirsts         ! first I index of source
      integer, intent(in)           :: ilasts          ! last I index of source
      integer, intent(in)           :: jfirsts         ! first j index of source
      integer, intent(in)           :: jlasts          ! last j index of source
      integer, intent(in)           :: kfirsts         ! first k index of source
      integer, intent(in)           :: klasts          ! last k index of source
      integer, intent(in)           :: ifirstr         ! first I index of target
      integer, intent(in)           :: ilastr          ! last I index of target
      integer, intent(in)           :: jfirstr         ! first j index of target
      integer, intent(in)           :: jlastr          ! last j index of target
      integer, intent(in)           :: kfirstr         ! first k index of target
      integer, intent(in)           :: klastr          ! last k index of target
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests
      real(r8), intent(in) :: qin(ifirsts:ilasts,jfirsts:jlasts,kfirsts:klasts,1:mq) ! input tracer array

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(ifirstr:ilastr,jfirstr:jlastr,kfirstr:klastr,1:mq) ! output tracer array
!
! !DESCRIPTION:
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications here and placing
!     wait points in mp_recvtrirr; if 1, call swap routine with p2p messages; if 2, call swap
!     routine with a2a messages. 
!     Modc(2): if 1, then apply handshaking (don't send until corresponding receive is posted)
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!     Modc(4): maximum number of outstanding requests (applies to swap routines only)
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!    09.10.07   Worley      eliminated mpi_recv from handshake logic
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nparcels, offset_s, offset_r, ierr, mod_method
      integer p, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize, unitsize, offset_0
      integer i, j, send_tag, recv_tag, num_s, num_r, m
      integer :: offset_v (Max_Nparcels)
      integer :: hs_snd, hs_rcv(numpro), hs_rcvids(numpro)
      integer ipe2, ceil2num
      integer numtr, numtrm
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_pid
      integer ip, jp, kp, mp, ir, jr, jir, mt


#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .true.
         maxreq_local = -1
      endif

! Do not call mp_swaptrirr unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

    if (sw_local .gt. 0) then
         sw_alltoall = (sw_local .eq. 2)
         call mp_swaptrirr(comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq, &
                           ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts, &
                           ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr, &
                           sw_handshake=hs_local, sw_maxreq=maxreq_local,     &
                           sw_alltoall=sw_alltoall, sw_send=send_local)
    else

      call MPI_COMM_RANK (comm, comm_pid, ierr)

      hs_snd = 1
      ceil2num = ceil2(numpro)

      numtrm = mend - mbeg
      numtr = numtrm + 1

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      r8_win%ncall_s = r8_win%ncall_s + 1
     if (mod_method .gt. 0) then
!
! mpi derived types
      if (r8_win%ncall_s .gt. MaxTrf-numtrm) then
         write(iulog,*) "mp_sendtrirr: derived type handle count exceeded - exiting"
         write(iulog,*) "r8_win%ncall_s MaxTrf = ", r8_win%ncall_s, MaxTrf
         stop
      endif
!
! MPI: Irecv over all processes
!
      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               if (ipe-1 /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, ipe-1, comm_pid, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      OutHandle(:,r8_win%ncall_s:r8_win%ncall_s+numtrm) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_r) cycle
!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          recv_tag = ipe-1 + modcam_tagoffset
          do m = mbeg, mend
             call mpi_irecv( qout(:,:,:,m), 1, recv_bl(ipe)%type, ipe-1, recv_tag,   &
                             comm, OutHandle(ipe,r8_win%ncall_s+m-mbeg), ierr )
          enddo
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
               call MPI_SEND ( hs_snd, 1, mp_i4, ipe-1, ipe-1, comm, ierr )
          endif
        endif
      enddo

!
! MPI: Isend/Send over all processes; use risend/rsend with hs
!
      InHandle(:,r8_win%ncall_s:r8_win%ncall_s+numtrm) = MPI_REQUEST_NULL
      do ipe2=1, ceil2num
        ipe = ieor(ipe2-1,comm_pid) + 1
        if (ipe .gt. num_s) cycle

!
! Send the individual buffers with non-blocking sends
!
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          send_tag = comm_pid + modcam_tagoffset
          if (hs_local) then
             if (ipe-1 /= comm_pid) &
                call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
             if (send_local) then
                do m = mbeg, mend
                   call mpi_rsend( qin(:,:,:,m), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, ierr )
                enddo
             else
                do m = mbeg, mend
                   call mpi_irsend( qin(:,:,:,m), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                    comm, InHandle(ipe,r8_win%ncall_s), ierr )
                enddo
             endif
          else
             if (send_local) then
                do m = mbeg, mend
                   call mpi_send( qin(:,:,:,m), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                  comm, ierr )
                enddo
             else
                do m = mbeg, mend
                   call mpi_isend( qin(:,:,:,m), 1, send_bl(ipe)%type, ipe-1, send_tag,        &
                                   comm, InHandle(ipe,r8_win%ncall_s), ierr )
                enddo
             endif
          endif
        endif
      enddo
     else

! temporary contiguous buffers

      jr = jlasts - jfirsts + 1
      ir = ilasts - ifirsts + 1
      jir = jr * ir
      if (r8_win%ncall_s .gt. max_irr-numtrm) then
         write(iulog,*) "mp_sendtrirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_s max_irr = ", r8_win%ncall_s, max_irr
         stop
      endif
      unitsize = r8_win%size/max_irr

! issue call to receive data in global receive buffer
      offset_0 = (r8_win%ncall_s-1)*unitsize
      offset_s = offset_0
      offset_r = offset_0

      if (hs_local) then
         hs_rcvids(:) = MPI_REQUEST_NULL
         do ipe2=1, ceil2num
            ipe = ieor(ipe2-1,comm_pid) + 1
            if (ipe .gt. num_s) cycle
            qsize = numtr*send_bl(ipe)%Tot_Size
            if (qsize .ne. 0) then
               r8_win%dest = ipe-1
               send_tag = comm_pid + modcam_tagoffset
               if (r8_win%dest /= comm_pid) &
                  call MPI_IRECV ( hs_rcv(ipe), 1, mp_i4, r8_win%dest, send_tag, comm, &
                                   hs_rcvids(ipe), ierr )
            endif
         enddo
      endif

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         r8_win%size_r = numtr*recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            r8_win%src = ipe-1
            if (numtr*unitsize >= offset_r-offset_0) then
              recv_tag = r8_win%src + modcam_tagoffset
              qsize    = r8_win%size_r
              r8_win%nrecv    = r8_win%nrecv + 1
              call MPI_IRECV(ga_r8_r(r8_win%offset_r+1), qsize, mp_r8, r8_win%src, &
                             recv_tag, comm, r8_win%rqest(r8_win%nrecv), ierror)
              if (hs_local) then
                 if (r8_win%src /= comm_pid) &
                   call MPI_SEND ( hs_snd, 1, mp_i4, r8_win%src, recv_tag, comm, ierror)
              endif
            else
              write(iulog,*) "Fatal mp_sendtrirr: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif
         endif
      enddo
! gather data into global send buffer
      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_s) cycle
         qsize = numtr*send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
            r8_win%dest = ipe-1
            r8_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s-offset_0 .gt. numtr*unitsize) then
              write(iulog,*) "Fatal mp_sendtrirr: send window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_s, offset_0
              stop
            endif

            offset_v(1) = r8_win%offset_s
            do j = 2, send_bl(ipe)%nparcels
               offset_v(j) = offset_v(j-1) + send_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, send_bl(ipe)%nparcels
               do m = mbeg, mend
                  do i = 1, send_bl(ipe)%blocksizes(j)
                     mp = send_bl(ipe)%displacements(j)+i
                     kp = kfirsts + (mp-1)/jir
                     mt = (kp-kfirsts)*jir
                     jp = jfirsts + (mp-mt-1)/ir
                     ip = mp-mt - (jp-jfirsts)*ir + ifirsts-1
                     ga_r8_s(send_bl(ipe)%Tot_Size*(m-mbeg)+offset_v(j)+i) = qin(ip,jp,kp,m)
                  enddo
               enddo
            enddo

! nonblocking send
            send_tag = comm_pid + modcam_tagoffset
            r8_win%nsend = r8_win%nsend + 1
            if (hs_local) then
               if (r8_win%dest /= comm_pid) &
                  call MPI_WAIT ( hs_rcvids(ipe), MPI_STATUS_IGNORE, ierr )
               if (send_local) then
                  call MPI_RSEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_IRSEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, r8_win%sqest(r8_win%nsend), ierr)
               endif
            else
               if (send_local) then
                  call MPI_SEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, ierr)
               else
                  call MPI_ISEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                                 send_tag, comm, r8_win%sqest(r8_win%nsend), ierr)
               endif
            endif
         endif
      enddo

     endif   !  mod_method

      r8_win%ncall_s = r8_win%ncall_s + numtrm

    endif   !  sw_local

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendtrirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvtrirr --- Finalize communication of contiguous tracer parcels
!
! !INTERFACE:
      subroutine mp_recvtrirr ( comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq,  &
                                ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts,  &
                                ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr,  &
                                modc )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer, intent(in)           :: mbeg            ! initial tracer index
      integer, intent(in)           :: mend            ! final tracer index
      integer, intent(in)           :: mq              ! total tracer indices
      integer, intent(in)           :: ifirsts         ! first I index of source
      integer, intent(in)           :: ilasts          ! last I index of source
      integer, intent(in)           :: jfirsts         ! first j index of source
      integer, intent(in)           :: jlasts          ! last j index of source
      integer, intent(in)           :: kfirsts         ! first k index of source
      integer, intent(in)           :: klasts          ! last k index of source
      integer, intent(in)           :: ifirstr         ! first I index of target
      integer, intent(in)           :: ilastr          ! last I index of target
      integer, intent(in)           :: jfirstr         ! first j index of target
      integer, intent(in)           :: jlastr          ! last j index of target
      integer, intent(in)           :: kfirstr         ! first k index of target
      integer, intent(in)           :: klastr          ! last k index of target
      integer, optional, intent(in) :: modc(4)         ! 1: classical, swap p2p, swap a2a
                                                       ! 2: handshake
                                                       ! 3: send vs isend
                                                       ! 4: max number of outstanding requests
      real(r8), intent(in) :: qin(ifirsts:ilasts,jfirsts:jlasts,kfirsts:klasts,1:mq) ! input tracer array
! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(ifirstr:ilastr,jfirstr:jlastr,kfirstr:klastr,1:mq) ! output tracer array
!
! !DESCRIPTION:
!     Complete transfer of a generalized region initiated by {\tt mp\_sendtrirr}.
!     Communicate a number of contiguous parcels to/from arbitrary set of PEs.
!     Modc(1): if 0, use original approach of posting all communications in mp_sendtrirr and
!     placing wait points here; otherwise don't do anything - mp_swaptrirr is called from mp_sendirr.
!     Modc(3): if 1, then use blocking send; otherwise use nonblocking send
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!    08.09.18   Mirin       Major overhaul, to include approaches from Mirin and Worley
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, offset_r, mod_method
      integer unitsize, offset_0
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer i, j, num_r, num_s, m
      integer :: offset_v (Max_Nparcels)
      integer ipe2, ceil2num
      integer numtr, numtrm
      integer sw_local, maxreq_local
      logical hs_local, send_local
      logical sw_alltoall
      integer comm_size, comm_pid
      integer ip, jp, kp, mp, ir, jr, jir, mt

      if (present(modc)) then
         sw_local   = modc(1)
         hs_local   = (modc(2) .eq. 1)
         send_local = (modc(3) .eq. 1)
         maxreq_local = modc(4)
      else
         sw_local = 0
         hs_local = .true.
         send_local = .true.
         maxreq_local = -1
      endif

! Do not call mp_swaptrirr (hence return) unless mod_method equals 0
      mod_method = recv_bl(1)%method
      if (mod_method .gt. 0) sw_local = 0

! Return if swap_irr
      if (sw_local .gt. 0) return

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

      ceil2num = ceil2(numpro)

      numtrm = mend - mbeg
      numtr = numtrm + 1

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      r8_win%ncall_r = r8_win%ncall_r + 1

    if (mod_method .gt. 0) then

! mpi derived types
      if (r8_win%ncall_r .gt. MaxTrf-numtrm) then
         write(iulog,*) "mp_recvtrirr: derived type handle count exceeded - exiting"
         write(iulog,*) "r8_win%ncall_r MaxTrf = ", r8_win%ncall_r, MaxTrf
         stop
      endif

      if (num_s .gt. 0 .and. (.not. send_local)) then
         do m = mbeg, mend
            CALL MPI_WAITALL( comm_size, InHandle(:,r8_win%ncall_r+m-mbeg), InStats, Ierr )
         enddo
      endif
      if (num_r .gt. 0) then
         do m = mbeg, mend
            CALL MPI_WAITALL( comm_size, OutHandle(:,r8_win%ncall_r+m-mbeg), OutStats, Ierr )
         enddo
      endif

    else

! temporary contiguous buffer / global window

      jr = jlastr - jfirstr + 1
      ir = ilastr - ifirstr + 1
      jir = jr * ir
      if (r8_win%ncall_r .gt. max_irr-numtrm) then
         write(iulog,*) "mp_recvtrirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_r max_irr = ", r8_win%ncall_r, max_irr
         stop
      endif
      unitsize = r8_win%size/max_irr

! scatter data from global receive buffer to final destination
      offset_0 = (r8_win%ncall_r-1)*unitsize
      offset_r = offset_0

      do ipe2=1, ceil2num
         ipe = ieor(ipe2-1,comm_pid) + 1
         if (ipe .gt. num_r) cycle
         r8_win%size_r = numtr*recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r-offset_0 .gt. numtr*unitsize) then
              write(iulog,*) "Fatal mp_recvtrirr: receive window out of space - exiting"
              write(iulog,*) 'comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif

            r8_win%nread = r8_win%nread + 1
            call MPI_WAIT(r8_win%rqest(r8_win%nread), Status, ierr)

            offset_v(1) = r8_win%offset_r
            do j = 2, recv_bl(ipe)%Nparcels
               offset_v(j) = offset_v(j-1) + recv_bl(ipe)%blocksizes(j-1)
            enddo

            do j = 1, recv_bl(ipe)%Nparcels
               do m = mbeg, mend
                  do i = 1, recv_bl(ipe)%blocksizes(j)
                     mp = recv_bl(ipe)%displacements(j)+i
                     kp = kfirstr + (mp-1)/jir
                     mt = (kp-kfirstr)*jir
                     jp = jfirstr + (mp-mt-1)/ir
                     ip = mp-mt - (jp-jfirstr)*ir + ifirstr-1
                     qout(ip,jp,kp,m) = ga_r8_r(recv_bl(ipe)%Tot_Size*(m-mbeg)+offset_v(j)+i)
                  enddo
               enddo
            enddo

         endif
      enddo

      if ((r8_win%ncall_s == r8_win%ncall_r + numtrm) .and. (.not. send_local)) then
         call MPI_WAITALL(r8_win%nsend, r8_win%sqest, Stats, ierror)
      endif

    endif    !    mod_method .gt. 0

    r8_win%ncall_r = r8_win%ncall_r + numtrm

    if (r8_win%ncall_s == r8_win%ncall_r) then
       r8_win%nsend = 0
       r8_win%nrecv = 0
       r8_win%nread = 0
       r8_win%ncall_s = 0
       r8_win%ncall_r = 0
    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvtrirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_swaptrirr --- Write r8 contiguous parcels to global array
!                            using XOR swap ordering - for multiple tracers
!
! !INTERFACE:
      subroutine mp_swaptrirr ( comm, send_bl, recv_bl, qin, qout, mbeg, mend, mq,  &
                                ifirsts, ilasts, jfirsts, jlasts, kfirsts, klasts,  &
                                ifirstr, ilastr, jfirstr, jlastr, kfirstr, klastr,  &
                                sw_handshake, sw_maxreq, sw_alltoall, sw_send   )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm                     ! communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer, intent(in)           :: mbeg            ! initial tracer index
      integer, intent(in)           :: mend            ! final tracer index
      integer, intent(in)           :: mq              ! total tracer indices
      integer, intent(in)           :: ifirsts         ! first I index of source
      integer, intent(in)           :: ilasts          ! last I index of source
      integer, intent(in)           :: jfirsts         ! first j index of source
      integer, intent(in)           :: jlasts          ! last j index of source
      integer, intent(in)           :: kfirsts         ! first k index of source
      integer, intent(in)           :: klasts          ! last k index of source
      integer, intent(in)           :: ifirstr         ! first I index of target
      integer, intent(in)           :: ilastr          ! last I index of target
      integer, intent(in)           :: jfirstr         ! first j index of target
      integer, intent(in)           :: jlastr          ! last j index of target
      integer, intent(in)           :: kfirstr         ! first k index of target
      integer, intent(in)           :: klastr          ! last k index of target
      logical, optional, intent(in) :: sw_handshake    ! use flow control and 
                                                       !  ready send
      integer, optional, intent(in) :: sw_maxreq       ! maximum number of outstanding
                                                       !  MPI requests
      logical, optional, intent(in) :: sw_alltoall     ! use mpi_alltoall
      logical, optional, intent(in) :: sw_send         ! use mpi_send instead of isend
      real(r8), intent(in) :: qin(ifirsts:ilasts,jfirsts:jlasts,kfirsts:klasts,1:mq) ! input tracer array

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(ifirstr:ilastr,jfirstr:jlastr,kfirstr:klastr,1:mq) ! output tracer array
!
! !DESCRIPTION:
!
!     XOR-ordered version of all-to-all communication
!
! WARNING: mod_comm parameter max_irr might need to be set larger than expected
!          when swapping multiple variables; specifically, max_irr must be at least
!          as large as the incoming r8_win%ncall_s + the number of variables to
!          be swapped
!
! !REVISION HISTORY: 
!    08.06.30   Worley      original: derived from mp_sendirr, but using 
!                            swapm logic and XOR swap order 
!    08.08.22   Worley      removed swapm; reimplemented with native MPI,
!                            added flow control/ready send option and maxreq
!                            throttling, added alltoall option
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: i, j, p, istep, num_s, num_r
      integer :: comm_pid, comm_size, steps, ierr
      integer :: ipe, offset_s, offset_r, offset_0, unitsize

      integer :: sndlths(0:numpro-1), sdispls(0:numpro-1)
      integer :: rcvlths(0:numpro-1), rdispls(0:numpro-1)
      integer :: swapids(numpro) 
      integer :: sndids(numpro)  ! nonblocking MPI send request ids
      integer :: rcvids(numpro)  ! nonblocking MPI recv request ids
      integer :: hs_snd, hs_rcv(numpro)! handshake variables (send/receive)
      integer :: hs_rcvids(numpro) ! nonblocking MPI handshake recv request ids
      integer :: InStats(numpro*MPI_STATUS_SIZE)
      integer :: OutStats(numpro*MPI_STATUS_SIZE)

      integer :: offset_v

      integer :: rstep

      integer :: maxreq, maxreqh
      logical :: handshake, alltoall, sendd
      integer :: ip, jp, kp, mp, irs, jrs, jirs, mt
      integer :: numtr, numtrm, irr, jrr, jirr, m

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, comm_size, ierr)
      call MPI_COMM_RANK (comm, comm_pid, ierr)

!     num_s = 0 if this process is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this process is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      if ( present(sw_handshake) ) then
         handshake = sw_handshake
         hs_snd = 1
      else
         handshake = .false.
      endif

      if ( present(sw_alltoall) ) then
         alltoall = sw_alltoall
      else
         alltoall = .false.
      endif

      if ( present(sw_send) ) then
         sendd = sw_send
      else
         sendd = .false.
      endif

      numtrm = mend - mbeg
      numtr = numtrm + 1
      jrs = jlasts - jfirsts + 1
      irs = ilasts - ifirsts + 1
      jirs = jrs * irs
      jrr = jlastr - jfirstr + 1
      irr = ilastr - ifirstr + 1
      jirr = jrr * irr

      unitsize = r8_win%size/max_irr

! advance to unused portion of storage window
      r8_win%ncall_s = r8_win%ncall_s + 1

      if (r8_win%ncall_s .gt. max_irr-numtrm) then
         write(iulog,*) "mp_swaptrirr: insufficient window storage - exiting"
         write(iulog,*) "r8_win%ncall_s max_irr = ", r8_win%ncall_s, max_irr
         stop
      endif

! calculate send lengths and displacements
      offset_0 = (r8_win%ncall_s-1)*unitsize
      offset_s = offset_0
      sndlths(:) = 0
      sdispls(:) = 0
      do ipe=1, num_s
         sndlths(ipe-1) = numtr*send_bl(ipe)%Tot_Size
         sdispls(ipe-1) = offset_s
         if (sndlths(ipe-1) .ne. 0) then

            offset_s = offset_s + sndlths(ipe-1)
            if (offset_s-offset_0 .gt. numtr*unitsize) then
              write(iulog,*) "Fatal mp_swaptrirr: send window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_s offset_0 = ', comm_pid,  &
                             ipe, unitsize, offset_s, offset_0
              stop
            endif
         endif
      enddo

! calculate receive lengths and displacements
      offset_r = offset_0
      rcvlths(:) = 0
      rdispls(:) = 0
      do ipe=1, num_r
         rcvlths(ipe-1) = numtr*recv_bl(ipe)%Tot_Size
         rdispls(ipe-1) = offset_r
         if (rcvlths(ipe-1) .ne. 0) then

            offset_r = offset_r + rcvlths(ipe-1)
            if (numtr*unitsize < offset_r-offset_0) then
              write(iulog,*) "Fatal mp_swaptrirr: receive window out of space - exiting"
              write(iulog,*) '1 comm_pid ipe unitsize offset_r offset_0 = ', comm_pid,  &
                        ipe, unitsize, offset_r, offset_0
              stop
            endif
         endif
      enddo

! Calculate swap partners and number of steps in point-to-point
! implementations of alltoall algorithm.
      steps = 0
      do ipe=1,ceil2(comm_size)-1
         p = pair(comm_size,ipe,comm_pid)
         if (p >= 0) then
            if (sndlths(p) > 0 .or. rcvlths(p) > 0) then
               steps = steps + 1
               swapids(steps) = p
            end if
         end if
      end do

      if (.not. alltoall) then

         sndids(1:steps) = MPI_REQUEST_NULL
         rcvids(1:steps) = MPI_REQUEST_NULL

         if (steps .eq. 0) then
            maxreq  = 0
            maxreqh = 0
         elseif (steps .eq. 1) then
            maxreq  = 1
            maxreqh = 1
         else
            if ( present(sw_maxreq) ) then
               if ((sw_maxreq .le. steps) .and. (sw_maxreq .ge. 0)) then
                  maxreq  = sw_maxreq
                  if (maxreq > 1) then
                     maxreqh = maxreq/2
                  else
                     maxreq  = 2
                     maxreqh = 1
                  endif
               else
                  maxreq  = steps
                  maxreqh = steps
               endif
            else
               maxreq  = steps
               maxreqh = steps
            endif
         endif

! Post initial handshake receive requests
         if (handshake) then
            do istep=1,maxreq
               p = swapids(istep)
               if (sndlths(p) > 0) then
                  call mpi_irecv  ( hs_rcv(istep), 1, mp_i4, p, comm_pid, comm, &
                                    hs_rcvids(istep), ierr )
               endif
            enddo
         endif

! Post initial receive requests
         do istep=1,maxreq
            p = swapids(istep)
            if (rcvlths(p) > 0) then
               offset_r = rdispls(p)+1
               call mpi_irecv ( ga_r8_r(offset_r), rcvlths(p), mp_r8, &
                                p, p, comm, rcvids(istep), ierr )
               if (handshake) then
                  call mpi_send( hs_snd, 1, mp_i4, p, p, comm, &
                                 ierr )
               endif
            endif
         enddo
         rstep = maxreq
!
      endif

! gather data into global send buffer
      do istep=1,steps
         p = swapids(istep)

         if (sndlths(p) .ne. 0) then
            offset_v = sdispls(p)
            do j = 1, send_bl(p+1)%nparcels
               do m = mbeg, mend
                  do i = 1, send_bl(p+1)%blocksizes(j)
                     mp = send_bl(p+1)%displacements(j)+i
                     kp = kfirsts + (mp-1)/jirs
                     mt = (kp-kfirsts)*jirs
                     jp = jfirsts + (mp-mt-1)/irs
                     ip = mp-mt - (jp-jfirsts)*irs + ifirsts-1
                     ga_r8_s(send_bl(p+1)%Tot_Size*(m-mbeg)+offset_v+i) = qin(ip,jp,kp,m)
                  enddo
               enddo
               offset_v = offset_v + send_bl(p+1)%blocksizes(j)
            enddo
         endif

         if (.not. alltoall) then

! Submit new i(r)send request
            offset_s = sdispls(p)+1
            if (sndlths(p) > 0) then
               if (handshake) then
                  call mpi_wait( hs_rcvids(istep), MPI_STATUS_IGNORE, ierr )
                  if (sendd) then
                     call mpi_rsend( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_irsend( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               else
                  if (sendd) then
                     call mpi_send ( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, ierr )
                  else
                     call mpi_isend ( ga_r8_s(offset_s), sndlths(p), mp_r8, &
                                      p, comm_pid, comm, sndids(istep), ierr )
                  endif
               endif
            endif

            if (istep > maxreqh) then
! Wait for oldest irecv request to complete
               call mpi_wait( rcvids(istep-maxreqh), OutStats, ierr )

               if (rstep < steps) then
                  rstep = rstep + 1
                  p = swapids(rstep)

! Submit a new handshake irecv request
                  if (handshake) then
                     if (sndlths(p) > 0) then
                        call mpi_irecv( hs_rcv(rstep), 1, mp_i4, p, comm_pid, comm, &
                                        hs_rcvids(rstep), ierr )
                     endif
                  endif

! Submit a new irecv request
                  if (rcvlths(p) > 0) then
                     offset_r = rdispls(p)+1
                     call mpi_irecv( ga_r8_r(offset_r), rcvlths(p), mp_r8, &
                                     p, p, comm, rcvids(rstep), ierr )
                     if (handshake) then
                        call mpi_send ( hs_snd, 1, mp_i4, p, p, comm, &
                                        ierr )
                     endif
                  endif
               endif

! Wait for outstanding i(r)send request to complete
               if (.not. sendd) then
                  call mpi_wait( sndids(istep-maxreqh), InStats, ierr )
               endif
            endif
!
         endif
!
      enddo

! local copy to send buffer
      if (sndlths(comm_pid) .ne. 0) then

         offset_v = sdispls(comm_pid)
         do j = 1, send_bl(comm_pid+1)%nparcels
            do m = mbeg, mend
               do i = 1, send_bl(comm_pid+1)%blocksizes(j)
                  mp = send_bl(comm_pid+1)%displacements(j)+i
                  kp = kfirsts + (mp-1)/jirs
                  mt = (kp-kfirsts)*jirs
                  jp = jfirsts + (mp-mt-1)/irs
                  ip = mp-mt - (jp-jfirsts)*irs + ifirsts-1
                  ga_r8_s(send_bl(comm_pid+1)%Tot_Size*(m-mbeg)+offset_v+i) = qin(ip,jp,kp,m)
               enddo
            enddo
            offset_v = offset_v + send_bl(comm_pid+1)%blocksizes(j)
         enddo

         if (.not. alltoall) then
            ga_r8_r(rdispls(comm_pid)+1:rdispls(comm_pid)+rcvlths(comm_pid)) = &
               ga_r8_s(sdispls(comm_pid)+1:sdispls(comm_pid)+sndlths(comm_pid))
         endif

      endif

      if (alltoall) then
         call mpi_alltoallv (ga_r8_s, sndlths, sdispls, mp_r8, &
                             ga_r8_r, rcvlths, rdispls, mp_r8, &
                             comm, ierror)
      endif

! local copy from receive buffer
      if (rcvlths(comm_pid) .ne. 0) then

         offset_v = rdispls(comm_pid)
         do j = 1, recv_bl(comm_pid+1)%Nparcels
            do m = mbeg, mend
               do i = 1, recv_bl(comm_pid+1)%blocksizes(j)
                  mp = recv_bl(comm_pid+1)%displacements(j)+i
                  kp = kfirstr + (mp-1)/jirr
                  mt = (kp-kfirstr)*jirr
                  jp = jfirstr + (mp-mt-1)/irr
                  ip = mp-mt - (jp-jfirstr)*irr + ifirstr-1
                  qout(ip,jp,kp,m) = ga_r8_r(recv_bl(comm_pid+1)%Tot_Size*(m-mbeg)+offset_v+i)
               enddo
            enddo
            offset_v = offset_v + recv_bl(comm_pid+1)%blocksizes(j)
         enddo

      endif

! scatter data from global receive buffer to final destination
      do istep=1,steps
         p = swapids(istep)

         if (.not. alltoall) then
            if (istep > steps-maxreqh) then
               call mpi_wait( rcvids(istep), OutStats, ierr )
            endif
         endif

         if (rcvlths(p) .ne. 0) then

            offset_v = rdispls(p)
            do j = 1, recv_bl(p+1)%Nparcels
               do m = mbeg, mend
                  do i = 1, recv_bl(p+1)%blocksizes(j)
                     mp = recv_bl(p+1)%displacements(j)+i
                     kp = kfirstr + (mp-1)/jirr
                     mt = (kp-kfirstr)*jirr
                     jp = jfirstr + (mp-mt-1)/irr
                     ip = mp-mt - (jp-jfirstr)*irr + ifirstr-1
                     qout(ip,jp,kp,m) = ga_r8_r(recv_bl(p+1)%Tot_Size*(m-mbeg)+offset_v+i)
                  enddo
               enddo
               offset_v = offset_v + recv_bl(p+1)%blocksizes(j)
            enddo

         endif
      enddo

! Wait for any outstanding send requests to complete.
      if (.not. alltoall .and. .not. sendd) then
         call mpi_waitall( maxreqh, sndids(steps-maxreqh+1), InStats, ierr )
      endif

! clean-up
! make used portion of storage window available for reuse
      r8_win%ncall_s = r8_win%ncall_s - 1

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_swaptrirr
# endif
!------------------------------------------------------------------------------
#endif
      end module mod_comm

