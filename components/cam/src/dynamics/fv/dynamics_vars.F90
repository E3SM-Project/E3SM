module dynamics_vars
!BOP
!
! !MODULE: dynamics_vars --- GEOS5/CAM fvcore internal variables
!
! !USES:
   use shr_kind_mod,       only: r8 => shr_kind_r8, r4 => shr_kind_r4

!
   use decompmodule,       only: decomptype
   use ghostmodule,        only: ghosttype
   use cam_logfile,        only: iulog
#if defined(SPMD)
   use parutilitiesmodule, only: parpatterntype, REAL4, INT4
#endif

! !PUBLIC MEMBER FUNCTIONS:
   public dynamics_init, dynamics_clean
   public a2d3d, d2a3d, b2d3d, d2b3d, c2a3d

! !PUBLIC DATA MEMBERS:

  public T_FVDYCORE_VARS, T_FVDYCORE_GRID, T_FVDYCORE_STATE

! T_FVDYCORE_VARS contains the prognostic variables for FVdycore
  type T_FVDYCORE_VARS
       real(r8), dimension(:,:,:  ), pointer     :: U      ! U winds (D-grid)
       real(r8), dimension(:,:,:  ), pointer     :: V      ! V winds (D-grid)
       real(r8), dimension(:,:,:  ), pointer     :: PT     ! scaled virtual pot. temp.
       real(r8), dimension(:,:,:  ), pointer     :: PE     ! Pressure at layer edges
       real(r8), dimension(:,:,:  ), pointer     :: PKZ    ! P^kappa mean
       real(r8), dimension(:,:,:,:), pointer     :: tracer ! Tracers
  end type T_FVDYCORE_VARS

! T_FVDYCORE_GRID contains information about the horizontal and vertical
! discretization, unlike in ARIES where these data are split into HORZ_GRID
! and VERT_GRID.  The reason for this: currently all of this information is
! initialized in one call to FVCAM dynamics_init.
 
  type T_FVDYCORE_GRID

!
! PILGRIM communication information (was in spmd_dyn)
!
    integer :: twod_decomp = 0  ! 1 for multi-2D decompositions, 0 otherwise

    integer :: npes_xy= 1    ! number of PEs for XY decomposition
    integer :: npes_yz= 1    ! number of PEs for YZ decomposition
    integer :: myid_y = 0    ! subdomain index (0-based) in latitude (y)
    integer :: myid_z = 0    ! subdomain index (0 based) in level (z)
    integer :: npr_y  = 1    ! number of subdomains in y
    integer :: npr_z  = 1    ! number of subdomains in z

    integer :: myidxy_x = 0  ! subdomain index (0-based) in longitude (x) (second. decomp.)
    integer :: myidxy_y = 0  ! subdomain index (0 based) in latitude (y) (second. decomp.)
    integer :: nprxy_x = 1   ! number of subdomains in x (second. decomp.)
    integer :: nprxy_y = 1   ! number of subdomains in y (second. decomp.)
    integer :: iam = 0       ! 

    integer :: mod_method = 0  ! 1 for mpi derived types with transposes, 0 for contiguous buffers
    integer :: mod_geopk = 0   ! 1 for mpi derived types with transposes, 0 for contiguous buffers
    integer :: mod_gatscat = 0 ! 1 for mpi derived types with transposes, 0 for contiguous buffers

    type(decomptype) :: strip2d, strip2dx, strip3dxyz, strip3dxzy,          &
                        strip3dxyzp, strip3zaty, strip3dxzyp,               &
                        strip3yatz, strip3yatzp, strip3zatypt,              &
                        strip3kxyz, strip3kxzy, strip3kxyzp, strip3kxzyp,   &
                        strip3dyz, checker3kxy

    integer :: commdyn           ! communicator for all dynamics
    integer :: commxy            ! communicator for XY decomposition
    integer :: commyz            ! communicator for YZ decomposition
    integer :: commnyz           ! communicator for multiple YZ decomposition

    integer :: comm_y            ! communicator in latitude
    integer :: comm_z            ! communicator in vertical
    integer :: commxy_x          ! communicator in longitude (xy second. decomp.)
    integer :: commxy_y          ! communicator in latitude (xy second. decomp.)
    logical :: geopkdist         ! use distributed method for geopotential calculation 
                                 !  with 2D decomp.
    logical :: geopk16byte       ! use Z-parallel distributed method for geopotential 
                                 !  calculation with 2D decomp.; otherwise use Z-serial
                                 !  pipeline algorithm when using distributed algoritm
    integer :: geopkblocks       ! number of stages to use in Z-serial pipeline 
                                 !  (non-transpose) geopotential algorithm
    integer :: modc_dynrun(4)    ! 1: mod_comm irregular underlying communication method for dyn_run/misc
                                 ! 2: mod_comm irregular communication handshaking for dyn_run/misc
                                 ! 3: mod_comm irregular communication send protocol for dyn_run/misc
                                 ! 4: mod_comm irregular communication nonblocking request throttle for dyn_run/misc
    integer :: modc_cdcore(4)    ! 1: mod_comm irregular underlying communication method for cd_core/geopk
                                 ! 2: mod_comm irregular communication handshaking for cd_core/geopk
                                 ! 3: geopk_d and mod_comm irregular communication send protocol for cd_core/geopk
                                 ! 4: mod_comm irregular communication nonblocking request throttle for cd_core/geopk
    integer :: modc_gather(4)    ! 1: mod_comm irregular underlying communication method for gather
                                 ! 2: mod_comm irregular communication handshaking for gather
                                 ! 3: mod_comm irregular communication send protocol for gather
                                 ! 4: mod_comm irregular communication nonblocking request throttle for gather
    integer :: modc_scatter(4)   ! 1: mod_comm irregular underlying communication method for scatter
                                 ! 2: mod_comm irregular communication handshaking for scatter
                                 ! 3: mod_comm irregular communication send protocol for scatter
                                 ! 4: mod_comm irregular communication nonblocking request throttle for scatter 
    integer :: modc_tracer(4)    ! 1: mod_comm irregular underlying communication method for multiple tracers
                                 ! 2: mod_comm irregular communication handshaking for multiple tracers
                                 ! 3: mod_comm irregular communication send protocol for multiple tracers
                                 ! 4: mod_comm irregular communication nonblocking request throttle for multiple tracers 
    integer :: modc_onetwo       ! one or two simultaneous mod_comm irregular communications (excl. tracers)
    integer :: modc_tracers      ! max number of tracers for simultaneous mod_comm irregular communications

#if defined(SPMD)
    type (ghosttype)        :: ghostu_yz, ghostv_yz, ghostpt_yz,                   &
                               ghostpe_yz, ghostpkc_yz
    type (parpatterntype)   :: u_to_uxy, uxy_to_u, v_to_vxy, vxy_to_v,             &
                               ikj_yz_to_xy, ikj_xy_to_yz,                         &
                               ijk_yz_to_xy, ijk_xy_to_yz,                         &
                               pe_to_pexy, pexy_to_pe,                             &
                               pt_to_ptxy, ptxy_to_pt, pkxy_to_pkc,                &
                               r4_xy_to_yz, r4_yz_to_xy, q3_to_qxy3, qxy3_to_q3,   &
                               xy2d_to_yz2d, yz2d_to_xy2d, scatter_3d, gather_3d,  &
                               g_2dxy_r8, g_2dxy_r4, g_2dxy_i4,                    &
                               s_2dxy_r8, s_2dxy_r4, s_2dxy_i4,                    &
                               g_3dxyz_r8, g_3dxyz_r4, g_3dxyzp_r8, g_3dxyzp_r4,   &
                               s_3dxyz_r8, s_3dxyz_r4, s_3dxyzp_r8, s_3dxyzp_r4
#endif

!
! END PILGRIM communication information
!

    integer                         :: JFIRST           ! Start latitude (exclusive)
    integer                         :: JLAST            ! End latitude (exclusive)
 
!
    integer                         :: NG_C             ! Ccore ghosting
    integer                         :: NG_D             ! Dcore ghosting
    integer                         :: NG_S             ! Staggered grid ghosting for
                                                        ! certain arrays, max(ng_c+1,ng_d)
!
! For 2D decomposition (currently not used)
!
    integer                         :: IFIRSTXY         ! Start longitude (exclusive)
    integer                         :: ILASTXY          ! End longitude (exclusive)
    integer                         :: JFIRSTXY         ! Start latitude (exclusive)
    integer                         :: JLASTXY          ! End latitude (exclusive)
!
    integer                         :: IM               ! Full longitude dim
    integer                         :: JM               ! Full latitude dim (including poles)
!
    real(r8)                        :: DL
    real(r8)                        :: DP
    real(r8)                        :: ACAP
    real(r8)                        :: RCAP
!
    real(r8), dimension(:), pointer :: COSP             ! Cosine of lat angle -- volume mean
    real(r8), dimension(:), pointer :: SINP             ! Sine of lat angle -- volume mean
    real(r8), dimension(:), pointer :: COSE             ! Cosine at finite volume edge
    real(r8), dimension(:), pointer :: SINE             ! Sine at finite volume edge
    real(r8), dimension(:), pointer :: ACOSP            ! Reciprocal of cosine of lat angle
!
    real(r8), dimension(:), pointer :: ACOSU            ! Reciprocal of cosine of lat angle (staggered)
!
    real(r8), dimension(:), pointer :: COSLON           ! Cosine of longitudes - volume center
    real(r8), dimension(:), pointer :: SINLON           ! Sine of longitudes - volume center
    real(r8), dimension(:), pointer :: COSL5            ! Cosine of longitudes - volume center
    real(r8), dimension(:), pointer :: SINL5            ! Sine of longitudes - volume center
 
!
!   Variables which are used repeatedly in CD_CORE
!

      integer ::       js2g0
      integer ::       jn2g0
      integer ::       jn1g1

      real(r8), pointer :: trigs(:)
      real(r8), pointer :: fc(:), f0(:)
      real(r8), pointer :: dc(:,:), de(:,:), sc(:), se(:)
      real(r8), pointer :: cdx(:,:), cdy(:,:)
      real(r8), pointer :: cdx4(:,:), cdy4(:,:) !for div4 damping
      real(r8), pointer :: cdxde(:,:), cdxdp(:,:),cdyde(:,:),cdydp(:,:)      !for del2 damping
      real(r8), pointer :: cdxdiv(:,:), cdydiv(:,:), cdtau4(:,:) !for del2 damping

      real(r8), pointer :: dcdiv4(:,:), dediv4(:,:), scdiv4(:), sediv4(:)    !for div4 damping

      real(r8), pointer :: dtdx(:), dtdxe(:), txe5(:), dtxe5(:)
      real(r8), pointer :: dyce(:),   dx(:) ,  rdx(:),    cy(:)
      real(r8), pointer :: dtdx2(:), dtdx4(:),  dxdt(:), dxe(:)
      real(r8), pointer :: cye(:),    dycp(:),  rdxe(:)

      real(r8) :: rdy, dtdy, dydt, dtdy5, tdy5
      real(r8) :: dt0 = 0

      integer  :: ifax(13)

      real(r8) ::  zt_c
      real(r8) ::  zt_d

!
! This part refers to the vertical grid
!
    integer                         :: KM              ! Numer of levels
    integer                         :: KMAX            ! KM+1 (?)
!
! For 2D decomposition (currently not used)
!
    integer                         :: KFIRST          ! Start level (exclusive)
    integer                         :: KLAST           ! End level (exclusive)
    integer                         :: KLASTP          ! klast+1, except km+1 when klastp=km+1
!
!
    integer                         :: KORD            ! monotonicity order for mapping (te_map)
    integer                         :: KS              ! Number of true pressure levels (out of KM+1)
    real(r8)                        :: PTOP            ! pressure at top (ak(1))
    real(r8)                        :: PINT            ! initial pressure (ak(km+1))
    real(r8), dimension(:), pointer :: AK              ! Sigma mapping
    real(r8), dimension(:), pointer :: BK              ! Sigma mapping

!
! Tracers
!
    integer                         :: NQ              ! Number of advected tracers
    integer                         :: NTOTQ           ! Total number of tracers (NQ <= NC)

! Extra subdomain bounds for cd_core/trac2d overlap and trac2d decomposition
! Relevant for secondary yz decomposition only; refers back to primary yz decomposition
    integer                         :: JFIRSTCT          ! jfirst
    integer                         :: JLASTCT           ! jlast
    integer                         :: KFIRSTCT          ! kfirst
    integer                         :: KLASTCT           ! klast

! Bounds for tracer decomposition
    integer, dimension(:), pointer  :: ktloa             ! lower tracer index (global map)
    integer, dimension(:), pointer  :: kthia             ! upper tracer index (global map)
    integer                         :: ktlo              ! lower tracer index (local)
    integer                         :: kthi              ! upper tracer index (local)

  end type T_FVDYCORE_GRID

! Constants used by fvcore
  type T_FVDYCORE_CONSTANTS
    real(r8)                             :: pi
    real(r8)                             :: omega    ! angular velocity of earth's rotation  
    real(r8)                             :: cp       ! heat capacity of air at constant pressure
    real(r8)                             :: ae       ! radius of the earth (m)
    real(r8)                             :: rair     ! Gas constant of the air
    real(r8)                             :: cappa    ! Cappa?
    real(r8)                             :: zvir     ! RWV/RAIR-1
  end type T_FVDYCORE_CONSTANTS

  integer, parameter :: NUM_FVDYCORE_ALARMS        = 3
  integer, parameter :: NUM_TIMES      = 8

  type T_FVDYCORE_STATE
!!!    private
    type (T_FVDYCORE_VARS)               :: VARS
    type (T_FVDYCORE_GRID )              :: GRID
    type (T_FVDYCORE_CONSTANTS)          :: CONSTANTS
#if defined( GEOS_MODE )
    type (ESMF_Clock), pointer           :: CLOCK
    type (ESMF_Alarm)                    :: ALARMS(NUM_FVDYCORE_ALARMS)
#endif
    integer(kind=8)                      :: RUN_TIMES(4,NUM_TIMES)
    logical                              :: DOTIME, DODYN
    real(r8)                             :: DT          ! Large time step
    real(r8)                             :: CHECK_DT    ! Time step to check maxmin
    integer                              :: ICD, JCD    ! Algorithm orders (C Grid)
    integer                              :: IORD, JORD  ! Algorithm orders (D Grid)
    integer                              :: KORD        ! Vertical order
    integer                              :: TE_METHOD   ! method for total energy mapping (te_map)
    logical                              :: CONSV       ! dycore conserves tot. en.
    integer                              :: NSPLIT
    integer                              :: NSPLTRAC
    integer                              :: NSPLTVRM
    integer                              :: NUM_CALLS
    integer                              :: FILTCW      ! filter c-grid winds if positive
  end type T_FVDYCORE_STATE

!
! !DESCRIPTION:
!
!      This module provides variables which are specific to the Lin-Rood
!      dynamical core.  Most of them were previously SAVE variables in 
!      different routines and were set with an "if (first)" statement.
!
!      \begin{tabular}{|l|l|} \hline \hline
!        lr\_init    &  Initialize the Lin-Rood variables  \\ \hline
!        lr\_clean   &  Deallocate all internal data structures \\ \hline 
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.06.06   Sawyer     Consolidated from various code snippets
!   01.07.12   Sawyer     Removed CCM common blocks comtim.h and commap.h
!   03.06.25   Sawyer     Cleaned up, used ParPatternCopy (Create)
!   03.07.23   Sawyer     Removed dependencies on params.h, constituents
!   03.08.05   Sawyer     Removed rayf_init and hswf_init, related vars
!   03.09.17   Sawyer     Removed unneeded ghost definitions
!   03.10.22   Sawyer     pmgrid removed (now spmd_dyn)
!   03.11.18   Sawyer     Removed set_eta (ak, bk, now read from restart)
!   03.12.04   Sawyer     Moved T_FVDYCORE_GRID here (removed some vars)
!   04.08.25   Sawyer     Removed all module data members, now GRID only
!   04.10.06   Sawyer     Added spmd_dyn vars here; ESMF transpose vars
!   05.04.12   Sawyer     Added support for r4/r8 tracers
!   05.05.24   Sawyer     CAM/GEOS5 merge (removed GEOS_mod dependencies)
!   05.06.10   Sawyer     Scaled down version for CAM (no ESMF)
!   05.11.10   Sawyer     Removed dyn_interface (now in dyn_comp)
!   06.03.01   Sawyer     Removed m_ttrans, q_to_qxy, qxy_to_q, etc.
!   06.05.09   Sawyer     Added CONSV to dyn_state (conserve energy)
!   06.08.27   Sawyer     Removed unused ESMF code for RouteHandle
!
!EOP
!-----------------------------------------------------------------------
   real(r8), parameter ::  D0_0                    =   0.0_r8
   real(r8), parameter ::  D0_5                    =   0.5_r8
   real(r8), parameter ::  D1_0                    =   1.0_r8
   real(r8), parameter ::  D2_0                    =   2.0_r8
   real(r8), parameter ::  D4_0                    =   4.0_r8
   real(r8), parameter ::  D180_0                  = 180.0_r8
   real(r8), parameter ::  ratmax                  =  0.81_r8


contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dynamics_init --- initialize the lin-rood dynamical core
!
! !INTERFACE: 
   subroutine dynamics_init( dt, jord, im, jm, km,              &
                             pi, ae, om, nq, ntotq,             &
                             ks, ifirstxy, ilastxy,             &
                             jfirstxy, jlastxy,                 &
                             jfirst, jlast, kfirst, klast,      &
                             npes_xy, npes_yz,                  &
                             commdyn, commxy, commyz, commnyz,  &
                             nprxy_x, nprxy_y, npryz_y, npryz_z,&
                             imxy, jmxy, jmyz, kmyz,            &
                             ak, bk, unit,                      &
                             grid )
! !USES:
      use fv_control_mod, only: trac_decomp
      implicit none

! !INPUT PARAMETERS:
      real(r8), intent(in)  :: dt                    !  Initial time step
      integer, intent(in)   :: jord                  !  Horz. scheme #
      integer, intent(in)   :: im, jm, km            !  Global dims
      real(r8), intent(in)  :: pi                    !  Pi
      real(r8), intent(in)  :: ae                    !  Earth radius
      real(r8), intent(in)  :: om                    !  Earth angular velocity
      integer, intent(in)   :: nq                    !  No. adv. tracers
      integer, intent(in)   :: ntotq                 !  No. total tracers
      integer, intent(in)   :: ks                    !  True # pressure levels
      integer, intent(in)   :: ifirstxy, ilastxy     !  Interval
      integer, intent(in)   :: jfirstxy, jlastxy     !  Interval
      integer, intent(in)   :: jfirst, jlast         !  Interval
      integer, intent(in)   :: kfirst, klast         !  Interval
      integer, intent(in)   :: nprxy_x       ! XY decomp - Nr in X
      integer, intent(in)   :: nprxy_y       ! XY decomp - Nr in Y
      integer, intent(in)   :: npryz_y       ! YZ decomp - Nr in Y
      integer, intent(in)   :: npryz_z       ! YZ decomp - Nr in Z
      integer, intent(in)   :: npes_xy       ! XY decomp - Total nr.
      integer, intent(in)   :: npes_yz       ! YZ decomp - Total nr.
      integer, intent(in)   :: commdyn       ! Communicator for dynamics
      integer, intent(in)   :: commxy        ! Communicator for XY decomp 
      integer, intent(in)   :: commyz        ! Communicator for YZ decomp 
      integer, intent(in)   :: commnyz       ! Communicator for multiple YZ decomp 

      integer, dimension(:), intent(in)   :: imxy
      integer, dimension(:), intent(in)   :: jmxy
      integer, dimension(:), intent(in)   :: jmyz
      integer, dimension(:), intent(in)   :: kmyz

      real(r8), dimension(:), intent(in)  :: ak
      real(r8), dimension(:), intent(in)  :: bk

      integer, intent(in)   :: unit

! !INPUT/OUTPUT PARAMETERS:
      type(T_FVDYCORE_GRID), intent(inout)  :: grid     ! Resulting grid

! !DESCRIPTION:
!
!   Initialize Lin-Rood specific variables
!
! !REVISION HISTORY:
!
!   01.06.06   Sawyer     Create
!   03.07.31   Sawyer     Added the 'layout' arguments
!   03.08.05   Sawyer     Removed hswf_init and rayf_init
!   04.08.25   Sawyer     Added GRID, contains all information
!   04.10.04   Sawyer     Added init_spmd here
!   06.03.01   Sawyer     Removed argument m_ttrans_in
!   06.11.27   Sawyer     Removed argument layout (no longer used)
!   06.11.29   Sawyer     Constant PI now passed as argument
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
!!!      real(r8) :: pi !WS  29.11.2006 -- uncomment this for zero diffs
      integer  :: rc

!  Set the basic grid variables

      grid%im       = im
      grid%jm       = jm
      grid%km       = km
      grid%kmax     = km + 1
      grid%nq       = nq  
      grid%ntotq    = ntotq
      grid%ks       = ks  
      grid%ifirstxy = ifirstxy
      grid%ilastxy  = ilastxy
      grid%jfirstxy = jfirstxy
      grid%jlastxy  = jlastxy
      grid%jfirst   = jfirst
      grid%jlast    = jlast
      grid%kfirst   = kfirst
      grid%klast    = klast
      if ( klast == km ) then
        grid%klastp = km+1
      else
        grid%klastp = klast
      endif

!WS  29.11.2006 -- uncomment this for zero diffs
!!!      pi  = D4_0 * atan(D1_0)
!!!      call dynpkg_init( pi, ae_in, om_in, dt_in, im_in, &
      call dynpkg_init( pi, ae, om, dt, im, &
                        jm, jord, grid )

!
! Level-dependent variables  (was in vert_init, now removed)
!
      ALLOCATE(GRID%AK(km+1))
      ALLOCATE(GRID%BK(km+1))

      GRID%AK   = AK
      GRID%BK   = BK
      GRID%PTOP = GRID%AK(1)
      GRID%PINT = GRID%AK(ks+1)

!
! Tracer decomposition limits
!
      allocate(grid%ktloa(trac_decomp))
      allocate(grid%kthia(trac_decomp))
      grid%ktloa(:) = 1
      grid%kthia(:) = nq
      grid%ktlo = 1
      grid%kthi = nq

#if defined( SPMD )
      call spmd_vars_init( nprxy_x, nprxy_y, npryz_y, npryz_z,          &
                           npes_xy, npes_yz, commdyn, commxy, commyz,   &
                           commnyz, imxy, jmxy, jmyz, kmyz, nq,         &
                           grid )
#else
      grid%npes_yz = 1
#endif
      return

CONTAINS

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  dynpkg_init --- Initialization for dynamics package
!
! !INTERFACE:
subroutine dynpkg_init( pi, ae, om, dt, im, jm, jord, grid )

! !USES:
   use pft_module, only : pftinit, pft2d, pft_cf
   implicit none

! !INPUT PARAMETERS:
      real(r8) , intent(in) :: pi
      real(r8) , intent(in) :: ae
      real(r8) , intent(in) :: om
      real(r8) , intent(in) :: dt
      integer, intent(in)   :: im
      integer, intent(in)   :: jm
      integer, intent(in)   :: jord

! !INPUT/OUTPUT PARAMETERS:
      type( T_FVDYCORE_GRID ), intent(inout) :: grid


! !DESCRIPTION:
! 
!   {\bf Purpose:} Initialization of the FV specific GRID vars
! 
! !REVISION HISTORY: 
!   00.01.10    Grant        Creation using code from SJ Lin
!   01.03.26    Sawyer       Added ProTeX documentation
!   01.06.06    Sawyer       Modified for dynamics_vars
!   04.08.25    Sawyer       Now updates GRID
!   05.06.30    Sawyer       Added initializations from cd_core
!   06.09.15    Sawyer       PI now passed as argument
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer  :: i, j, imh, js2g0, jn2g0, jn1g1, js2gc, jn1gc
      integer  :: js2gs, jn2gd, jn1gs
      real(r8) :: zam5, zamda
      real(r8) :: ph5      ! This is to ensure 64-bit for any choice of r8

      real(r8), pointer :: coslon(:), sinlon(:), cosl5(:), sinl5(:)
      real(r8), pointer :: cosp(:), sinp(:), cose(:), sine(:), acosp(:), acosu(:)

!
! Local variables from cd_core
!

      integer  :: icffta
      real(r8) :: rcffta

      real(r8) :: rat, ycrit, dt5

!
! Start initialization
!
      grid%dl  = (pi+pi)/im
      grid%dp  = pi/(jm-1)

      allocate(grid%cosp(jm))
      allocate(grid%sinp(jm))
      allocate(grid%cose(jm))
      allocate(grid%sine(jm))
      allocate(grid%acosp(jm))
      allocate(grid%acosu(jm))

      allocate(grid%coslon(im))
      allocate(grid%sinlon(im))
      allocate(grid%cosl5(im))
      allocate(grid%sinl5(im))

      cosp => grid%cosp
      sinp => grid%sinp
      cose => grid%cose
      sine => grid%sine
      acosp  => grid%acosp
      acosu  => grid%acosu

      coslon => grid%coslon
      sinlon => grid%sinlon
      cosl5  => grid%cosl5
      sinl5  => grid%sinl5

      do j=2,jm
         ph5  = -D0_5*pi + ((j-1)-D0_5)*(pi/(jm-1))
         sine(j) = sin(ph5)
      enddo

      cosp( 1) =  D0_0
      cosp(jm) =  D0_0
      !
      ! cos(theta) at cell center distretized as
      !
      ! cos(theta) = d(sin(theta))/d(theta)
      !
      do j=2,jm-1
         cosp(j) = (sine(j+1)-sine(j)) / grid%dp
      enddo

! Define cosine at edges..

      do j=2,jm
         cose(j) = D0_5 * (cosp(j-1) + cosp(j))
      enddo
         cose(1) = cose(2)

      do j=2,jm-1
         acosu(j) = D2_0 / (cose(j) + cose(j+1))
      enddo

         sinp( 1) = -D1_0
         sinp(jm) =  D1_0

      do j=2,jm-1
         sinp(j) = D0_5 * (sine(j) + sine(j+1))
      enddo


!
! Pole cap area and inverse
      grid%acap = im*(D1_0+sine(2)) / grid%dp
      grid%rcap = D1_0 / grid%acap
 
      imh = im/2
      if(im .ne. 2*imh) then
         write(iulog,*) 'im must be an even integer'
         stop
      endif
 
! Define logitude at the center of the volume
! i=1, Zamda = -pi
 
      do i=1,imh
         zam5          = ((i-1)-D0_5) * grid%dl
         cosl5(i)      =  cos(zam5)
         cosl5(i+imh)  = -cosl5(i)
         sinl5(i)      =  sin(zam5)
         sinl5(i+imh)  = -sinl5(i)
         zamda         = (i-1)*grid%dl
         coslon(i)     =  cos(zamda)
         coslon(i+imh) = -coslon(i)
         sinlon(i)     =  sin(zamda)
         sinlon(i+imh) = -sinlon(i)
      enddo

      do j=2,jm-1
         acosp(j) = D1_0 / cosp(j)
      enddo
      acosp( 1) = grid%rcap * im
      acosp(jm) = grid%rcap * im

#if defined( SPMD )
!
! Calculate the ghost region sizes for the SPMD version (tricky stuff)
!
      grid%ng_c = 2                    ! Avoid the case where ng_c = 1
      grid%ng_d = min( abs(jord), 3)   ! SJL: number of max ghost latitudes
      grid%ng_d = max( grid%ng_d, 2)
      grid%ng_s = max( grid%ng_c+1, grid%ng_d )
#else
      grid%ng_c = 0
      grid%ng_d = 0                   ! No ghosting necessary for pure SMP runs
      grid%ng_s = 0
#endif

!
! cd_core initializations
!

      allocate(grid%dtdx(jm))
      allocate(grid%dtdx2(jm))
      allocate(grid%dtdx4(jm))
      allocate(grid%dtdxe(jm))
      allocate(grid%dxdt(jm))
      allocate(grid%dxe(jm))
      allocate(grid%cye(jm))
      allocate(grid%dycp(jm))
      allocate(grid%rdxe(jm))
      allocate(grid%txe5(jm))
      allocate(grid%dtxe5(jm))
      allocate(grid%dyce(jm))
      allocate(grid%dx(jm))
      allocate(grid%rdx(jm))
      allocate(grid%cy(jm))

      js2g0  = max(2,grid%jfirst)
      jn2g0  = min(jm-1,grid%jlast)
      jn1g1  = min(jm,grid%jlast+1)
      js2gc  = max(2,grid%jfirst-grid%ng_c) ! NG lats on S (starting at 2)
      jn1gc  = min(jm,grid%jlast+grid%ng_c) ! ng_c lats on N (ending at jm)

      grid%js2g0  = js2g0
      grid%jn2g0  = jn2g0
      grid%jn1g1  = jn1g1

      js2gs = max(2,grid%jfirst-grid%ng_s)
      jn2gd = min(jm-1,grid%jlast+grid%ng_d)
      jn1gs = min(jm,grid%jlast+grid%ng_s)

      allocate(grid%sc(js2g0:jn2g0))
      allocate(grid%se(js2g0:jn1g1))
      allocate(grid%dc(im,js2g0:jn2g0))
      allocate(grid%de(im,js2g0:jn1g1))

      allocate(grid%scdiv4(js2gs:jn2gd))   !for filtering of u and v in div4 damping 
      allocate(grid%sediv4(js2gs:jn1gs))   !for filtering of u and v in div4 damping 
      allocate(grid%dcdiv4(im,js2gs:jn2gd))!for filtering of u and v in div4 damping 
      allocate(grid%dediv4(im,js2gs:jn1gs))!for filtering of u and v in div4 damping 

      call pftinit(im)

! Determine ycrit such that effective DX >= DY
      rat = real(im,r8)/real(2*(jm-1),r8)
      ycrit = acos( min(ratmax, rat) ) * (D180_0/pi)

      call pft_cf(im, jm, js2g0, jn2g0, jn1g1, &
                  grid%sc, grid%se, grid%dc, grid%de,  &
                  grid%cosp, grid%cose, ycrit)

      !for filtering of u and v in div4 damping 
      !(needs larger halo than cam3.5 code)
      call pft_cf(im, jm, js2gs, jn2gd, jn1gs,                             & 
                  grid%scdiv4, grid%sediv4, grid%dcdiv4, grid%dediv4,      & 
                  grid%cosp, grid%cose, ycrit)                               


      allocate( grid%cdx   (js2g0:jn1g1,grid%kfirst:grid%klast) )
      allocate( grid%cdy   (js2g0:jn1g1,grid%kfirst:grid%klast) )

      allocate( grid%cdx4  (js2g0:jn1g1,grid%kfirst:grid%klast) )!for div4 damping
      allocate( grid%cdy4  (js2g0:jn1g1,grid%kfirst:grid%klast) )!for div4 damping

      allocate( grid%cdxde (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping
      allocate( grid%cdxdp (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping
      allocate( grid%cdyde (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping
      allocate( grid%cdydp (js2g0:jn1g1,grid%kfirst:grid%klast) )!for del2 damping

      allocate( grid%cdxdiv(jm,grid%kfirst:grid%klast) )!for div4 damping
      allocate( grid%cdydiv(jm,grid%kfirst:grid%klast) )!for div4 damping
      allocate( grid%cdtau4(js2g0:jn1g1,grid%kfirst:grid%klast) )!for div4 damping

! 000304 bug fix: ng_s not ng_d
      allocate( grid%f0(grid%jfirst-grid%ng_s-1:grid%jlast+grid%ng_d) )
      allocate( grid%fc(js2gc:jn1gc) )

! 000304 bug fix
      do j=max(1,grid%jfirst-grid%ng_s-1),min(jm,grid%jlast+grid%ng_d)
        grid%f0(j) = (om+om)*grid%sinp(j)
      enddo

! Compute coriolis parameter at cell corners.
      do j=js2gc, jn1gc                    ! Not the issue with ng_c = ng_d
         grid%fc(j) = D0_5*(grid%f0(j) + grid%f0(j-1))
      enddo

!!!        grid%dt0 = dt
        grid%dt0 = D0_0
        dt5 = D0_5*dt

        grid%rdy   = D1_0/(ae*grid%dp)
        grid%dtdy  = dt *grid%rdy
        grid%dtdy5 = dt5*grid%rdy
        grid%dydt  = (ae*grid%dp) / dt
        grid%tdy5  = D0_5/grid%dtdy

      return
!EOC
end subroutine dynpkg_init
!-----------------------------------------------------------------------

#if defined(SPMD)
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  spmd_vars_init --- Initialization of SPMD-related variables
!
! !INTERFACE:
subroutine spmd_vars_init( nprxy_x, nprxy_y, npryz_y, npryz_z,          &
                           npes_xy, npes_yz, commdyn, commxy, commyz,   &
                           commnyz, imxy, jmxy, jmyz, kmyz, nq,         &
                           grid )

! !USES:
   use decompmodule, only: decompcreate, decompfree
   use ghostmodule, only : ghostcreate, ghostfree
   use parutilitiesmodule, only : gid, parpatterncreate, parsplit
   use mpishorthand, only: mpiint
   use fv_control_mod, only: ct_overlap, trac_decomp
   implicit none

! !INPUT PARAMETERS:
   integer, intent(in)   :: nprxy_x    ! XY decomp - Nr in X
   integer, intent(in)   :: nprxy_y    ! XY decomp - Nr in Y
   integer, intent(in)   :: npryz_y    ! YZ decomp - Nr in Y
   integer, intent(in)   :: npryz_z    ! YZ decomp - Nr in Z
   integer, intent(in)   :: npes_xy    ! XY decomp - Total nr
   integer, intent(in)   :: npes_yz    ! YZ decomp - Total nr
   integer, intent(in)   :: commdyn    ! Communicator for all dyn
   integer, intent(in)   :: commxy     ! Communicator for XY decomp
   integer, intent(in)   :: commyz     ! Communicator for YZ decomp
   integer, intent(in)   :: commnyz    ! Communicator for multiple YZ decomp

   integer, dimension(:), intent(in) :: imxy
   integer, dimension(:), intent(in) :: jmxy
   integer, dimension(:), intent(in) :: jmyz
   integer, dimension(:), intent(in) :: kmyz
   integer, intent(in)               :: nq

! !INPUT/OUTPUT PARAMETERS:
      type( T_FVDYCORE_GRID ), intent(inout) :: grid

! !DESCRIPTION:
! 
!   {\bf Purpose:} Initialization of the SPMD related variables.
!   This has to be done in this module since certain variables
!   (in particular the ghost sizes {\tt ng\_d, ng\_s} are first
!   defined here.
! 
! !REVISION HISTORY: 
!   02.11.08    Sawyer       Creation
!   03.05.07    Sawyer       Use ParPatternCopy for q_to_qxy, etc.
!   03.07.23    Sawyer       Removed dependency on constituents module
!   03.09.10    Sawyer       Reactivated u_to_uxy, etc, redefined pe2pexy
!   03.11.19    Sawyer       Merged in CAM code with mod_method
!   04.08.25    Sawyer       Added GRID as argument
!   04.09.30    Sawyer       Initial ESMF routehandlers
!   04.10.04    Sawyer       Added INIT_SPMD functionality
!   06.08.27    Sawyer       Removed ESMF routehandles -- non-current
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
#if defined( GEOS_MODE )
      character(len=ESMF_MAXSTR), parameter :: IAm='spmd_vars_init'
#endif

! !LOCAL VARIABLES:
      type(decomptype) :: global2d, local2d

      integer   :: im, jm, km         !  Global dims
      integer   :: ifirstxy, ilastxy  !  Interval
      integer   :: jfirstxy, jlastxy  !  Interval
      integer   :: jfirst, jlast      !  Interval
      integer   :: kfirst, klast      !  Interval
      integer   :: ng_s, ng_c, ng_d   !  Ghost widths
      integer   :: rc                 !  return code
      integer   :: rank_y, rank_z, rankxy_x, rankxy_y  ! Currently not used
      integer   :: size_y, size_z, sizexy_x, sizexy_y  ! Currently not used

      integer :: xdist(1), ydistk(1), zdist1(1), zdistxy(1) ! non-distributed dims
      integer, allocatable :: xdist_global(:), ydist_global(:) 
      integer, allocatable :: zdist(:) ! number of levels per subdomain
      integer :: ier       ! error flag
      integer :: ig1, ig2, jg1, jg2, jg1d, jg2d, jg1s, jg2s, kg1, kg2, kg2p
      integer :: ktmod, ml
      integer :: myidmod
      integer :: ictstuff(4)
      integer :: kquot, krem, krun, kt, mlt

!
! Grab crucial variables from Grid
!
      im = grid%im
      jm = grid%jm
      km = grid%km

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      jfirst = grid%jfirst
      jlast  = grid%jlast
      kfirst = grid%kfirst
      klast  = grid%klast

      ng_s   = grid%ng_s
      ng_c   = grid%ng_c
      ng_d   = grid%ng_d


!
! This section of code used to be in INIT_SPMD (FVdycore_GridCompMod)
!
      grid%iam      = gid

      grid%npr_y    = npryz_y
      grid%npr_z    = npryz_z
      grid%nprxy_x  = nprxy_x
      grid%nprxy_y  = nprxy_y

      grid%npes_xy  = npes_xy
      grid%npes_yz  = npes_yz

      grid%myid_z   = gid/grid%npr_y
      grid%myid_y   = gid - grid%myid_z*grid%npr_y
      grid%myidxy_y = gid/grid%nprxy_x
      grid%myidxy_x = gid - grid%myidxy_y*grid%nprxy_x

      grid%commdyn  = commdyn
      grid%commxy   = commxy
      grid%commyz   = commyz
      grid%commnyz  = commnyz

! Split communicators

      call parsplit(commyz, grid%myid_z, gid, grid%comm_y, rank_y, size_y)
      call parsplit(commyz, grid%myid_y, gid, grid%comm_z, rank_z, size_z)
      call parsplit(commxy, grid%myidxy_y, gid, grid%commxy_x, rankxy_x, sizexy_x)
      call parsplit(commxy, grid%myidxy_x, gid, grid%commxy_y, rankxy_y, sizexy_y)

!
! WS: create decompositions for NCAR data structures
!
      allocate (xdist_global(nprxy_x))
      allocate (ydist_global(nprxy_y))
      allocate (zdist (npryz_z))
      xdist(1) = im
!
! Create PILGRIM decompositions (see decompmodule)
!
      if (gid .lt. npes_xy) then
        xdist_global = 0
        ydist_global = 0
        xdist_global(1) = im
        ydist_global(1) = jm
        call decompcreate( nprxy_x, nprxy_y, xdist_global,            &
                           ydist_global, global2d )
        call decompcreate( nprxy_x, nprxy_y, imxy, jmxy, local2d )
      endif

! Decompositions needed on xy decomposition for parpatterncreate
      if (gid .lt. npes_xy) then
        call decompcreate( 1, npryz_y, xdist, jmyz, grid%strip2d )
        call decompcreate( 1, npryz_y, npryz_z, xdist,                &
                           jmyz, kmyz, grid%strip3dxyz )
        call decompcreate( "xzy", 1, npryz_z, grid%npr_y, xdist,         &
                           kmyz, jmyz, grid%strip3dxzy )

! For y communication within z subdomain (klast version)
! Use myidmod to have valid index for inactive processes
!  for smaller yz decomposition
        myidmod = mod(grid%myid_z, grid%npr_z)  ! = myid_z for active yz process
        zdist1(1) = kmyz(myidmod+1)
        call decompcreate( 1, npryz_y, 1, xdist, jmyz, zdist1,        &
                           grid%strip3yatz )

! For z communication within y subdomain

        ydistk(1) = jmyz(grid%myid_y+1)
        call decompcreate( 1, 1, npryz_z, xdist, ydistk, kmyz,        &
                           grid%strip3zaty )

! Arrays dimensioned plev+1

        zdist(:) = kmyz(:)
        zdist(npryz_z) = kmyz(npryz_z) + 1
        call decompcreate( 1, npryz_y, npryz_z, xdist, jmyz, zdist,&
                           grid%strip3dxyzp )
        call decompcreate( "xzy", 1, npryz_z, npryz_y,                &
                           xdist, zdist, jmyz, grid%strip3dxzyp )

! Arrays dimensioned plev+1, within y subdomain

        ydistk(1) = jmyz(grid%myid_y+1)
        call decompcreate( "xzy", 1, npryz_z, 1, xdist, zdist, ydistk,   &
                         grid%strip3zatypt )

! For y communication within z subdomain (klast+1 version)
! Use myidmod to have valid index for inactive processes
!  for smaller yz decomposition
        myidmod = mod(grid%myid_z, grid%npr_z)  ! = myid_z for active yz process
        zdist1(1) = kmyz(myidmod+1)
        call decompcreate( 1, npryz_y, 1, xdist, jmyz, zdist1,        &
                           grid%strip3yatzp )

! For the 2D XY-YZ data transfer, we need a short 3D array
        zdist(:) = 1    ! One copy on each z PE set
        call decompcreate( 1, npryz_y, npryz_z,                       &
                           xdist, jmyz, zdist, grid%strip3dyz )
      endif

! Secondary xy decomposition
!
      if (grid%twod_decomp == 1) then
       if (gid .lt. npes_xy) then
        zdistxy(1) = npryz_z     ! All npr_z copies on 1 PE
        call decompcreate( nprxy_x, nprxy_y, 1,                     &
                           imxy, jmxy, zdistxy, grid%checker3kxy )
        zdistxy(1) = km
        call decompcreate( nprxy_x, nprxy_y, 1,                     &
                           imxy, jmxy, zdistxy, grid%strip3kxyz )
        call decompcreate( "xzy", nprxy_x, 1, nprxy_y,              &
                           imxy, zdistxy, jmxy, grid%strip3kxzy )

        zdistxy(1) = zdistxy(1) + 1
        call decompcreate( nprxy_x, nprxy_y, 1,                     &
                           imxy, jmxy, zdistxy, grid%strip3kxyzp )
        call decompcreate( "xzy", nprxy_x, 1, nprxy_y,              &
                           imxy, zdistxy, jmxy, grid%strip3kxzyp )
        zdistxy(1) = jlastxy - jfirstxy + 1
        call decompcreate( nprxy_x, 1, imxy, zdistxy, grid%strip2dx )
       endif
      endif

      deallocate(zdist)
      deallocate(ydist_global)
      deallocate(xdist_global)
!
! End of section imported from INIT_SPMD (FVdycore_GridCompMod)
!

      if ( grid%twod_decomp == 1 ) then
! Initialize ghost regions
!
   !!!     call t_startf('ghost_creation')


! Set limits for ghostcreate
       ig1 = 1
       ig2 = im
       jg1 = jfirst
       jg2 = jlast
       jg1d = jfirst-ng_d
       jg1s = jfirst-ng_s
       jg2d = jlast+ng_d
       jg2s = jlast+ng_s
       kg1 = kfirst
       kg2 = klast
       kg2p = klast+1

! Call ghostcreate with null ranges for non-yz processes
       if (gid .ge. npes_yz) then
          ig1 = im/2
          ig2 = ig1 - 1
          jg1 = (jfirst+jlast)/2
          jg2 = jg1 - 1
          jg1d = jg1
          jg1s = jg1
          jg2d = jg2
          jg2s = jg2
          kg1 = (kfirst+klast)/2
          kg2 = kg1 - 1
          kg2p = kg2
       endif

! Ghosted decompositions needed on xy decomposition for parpatterncreate
       if (gid .lt. npes_xy) then
        call ghostcreate( grid%strip3dxyz, gid, im, ig1, ig2, .true., &
                          jm, jg1d, jg2s, .false., &
                          km, kg1, kg2, .false., grid%ghostu_yz )
        call ghostcreate( grid%strip3dxyz, gid, im, ig1, ig2, .true., &
                          jm, jg1s, jg2d, .false., &
                          km, kg1, kg2, .false., grid%ghostv_yz )
        call ghostcreate( grid%strip3dxyz, gid, im, ig1, ig2, .true., &
                          jm, jg1d, jg2d, .false., &
                          km, kg1, kg2, .false., grid%ghostpt_yz )
        call ghostcreate( grid%strip3dxzyp, gid, im, ig1, ig2, .true., &
                          km+1, kg1, kg2p, .false., &
                          jm, jg1, jg2, .false., grid%ghostpe_yz)
        call ghostcreate( grid%strip3dxyzp, gid, im, ig1, ig2, .true., &
                          jm, jg1, jg2, .false.,       &
                          km+1, kg1, kg2p, .false., grid%ghostpkc_yz)
       endif
   !!!     call t_stopf('ghost_creation')

! Initialize transposes
!
   !!!     call t_startf('transpose_creation')

       if (gid .lt. npes_xy) then
        call parpatterncreate(commxy, grid%ghostu_yz, grid%strip3kxyz, &
                              grid%u_to_uxy, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3kxyz,grid%ghostu_yz, &
                              grid%uxy_to_u, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%ghostv_yz, grid%strip3kxyz, &
                              grid%v_to_vxy, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3kxyz, grid%ghostv_yz, &
                              grid%vxy_to_v, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3dxyz, grid%strip3kxyz,&
                              grid%ijk_yz_to_xy, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3kxyz, grid%strip3dxyz,&
                              grid%ijk_xy_to_yz, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3dxzy, grid%strip3kxzy,&
                              grid%ikj_yz_to_xy, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3kxzy, grid%strip3dxzy,&
                              grid%ikj_xy_to_yz, mod_method=grid%mod_method)
!
! Note PE <-> PEXY has been redefined for PEXY ijk, but PE ikj
!
        call parpatterncreate(commxy, grid%ghostpe_yz, grid%strip3kxzyp, &
                              grid%pe_to_pexy, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3kxzyp, grid%ghostpe_yz, &
                              grid%pexy_to_pe, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%ghostpt_yz, grid%strip3kxyz,  &
                              grid%pt_to_ptxy, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3kxyz, grid%ghostpt_yz,  &
                              grid%ptxy_to_pt, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3dxyz, grid%strip3kxyz,  &
                              grid%r4_yz_to_xy, mod_method=grid%mod_method,  &
                              T = REAL4 )
        call parpatterncreate(commxy, grid%strip3kxyz, grid%strip3dxyz,  &
                              grid%r4_xy_to_yz, mod_method=grid%mod_method,  &
                              T = REAL4 )
        call parpatterncreate(commxy, grid%strip3kxyzp, grid%ghostpkc_yz, &
                              grid%pkxy_to_pkc, mod_method=grid%mod_method)
!
! These are for 'transposing' 2D arrays from XY YZ
        call parpatterncreate(commxy, grid%checker3kxy, grid%strip3dyz, &
                              grid%xy2d_to_yz2d, mod_method=grid%mod_method)
        call parpatterncreate(commxy, grid%strip3dyz, grid%checker3kxy, &
                              grid%yz2d_to_xy2d, mod_method=grid%mod_method)
       endif
   !!!     call t_stopf('transpose_creation')

!
! Free unneeded decompositions
!
!         call decompfree(grid%strip2dx)     ! apparently added after 3_3_47
          call decompfree(grid%strip3dxzyp)
          call decompfree(grid%strip3dyz)
          call decompfree(grid%strip3yatz)
          call decompfree(grid%strip3yatzp)
          call decompfree(grid%strip3zaty)
          call decompfree(grid%strip3zatypt)
          call decompfree(grid%strip3kxyz)
          call decompfree(grid%strip3kxzy)
          call decompfree(grid%strip3kxyzp)
          call decompfree(grid%strip3kxzyp)
          call decompfree(grid%checker3kxy)

          call ghostfree(grid%ghostu_yz)
          call ghostfree(grid%ghostv_yz)
          call ghostfree(grid%ghostpt_yz)
          call ghostfree(grid%ghostpe_yz)
          call ghostfree(grid%ghostpkc_yz)

      endif

#if !defined( GEOS_MODE )
!
! Define scatter and gather patterns for 2D and 3D unghosted arrays
!

     if (gid .lt. npes_xy) then
      call parpatterncreate( commxy, global2d, local2d, grid%s_2dxy_r8, &
                             mod_method=grid%mod_gatscat )
      call parpatterncreate( commxy, local2d, global2d, grid%g_2dxy_r8,  &
                             mod_method=grid%mod_gatscat )

      call parpatterncreate( commxy, global2d, local2d, grid%s_2dxy_r4, &
                             mod_method=grid%mod_gatscat, T = REAL4 )
      call parpatterncreate( commxy, local2d, global2d, grid%g_2dxy_r4,  &
                             mod_method=grid%mod_gatscat, T = REAL4 )

      call parpatterncreate( commxy, global2d, local2d, grid%s_2dxy_i4, &
                             mod_method=grid%mod_gatscat, T = INT4 )
      call parpatterncreate( commxy, local2d, global2d, grid%g_2dxy_i4,  &
                             mod_method=grid%mod_gatscat, T = INT4 )

!
! 3D XYZ patterns, will replace XZY patterns eventually
!
      call parpatterncreate( commxy, grid%s_2dxy_r8, grid%s_3dxyz_r8, km )
      call parpatterncreate( commxy, grid%g_2dxy_r8, grid%g_3dxyz_r8, km )
      call parpatterncreate( commxy, grid%s_2dxy_r8, grid%s_3dxyzp_r8, km+1 )
      call parpatterncreate( commxy, grid%g_2dxy_r8, grid%g_3dxyzp_r8, km+1 )

      call parpatterncreate( commxy, grid%s_2dxy_r4, grid%s_3dxyz_r4, km )
      call parpatterncreate( commxy, grid%g_2dxy_r4, grid%g_3dxyz_r4, km )
      call parpatterncreate( commxy, grid%s_2dxy_r4, grid%s_3dxyzp_r4, km+1 )
      call parpatterncreate( commxy, grid%g_2dxy_r4, grid%g_3dxyzp_r4, km+1 )

#endif

      call decompfree( global2d )
      call decompfree( local2d )
     endif

! Secondary subdomain limits for cd_core/trac2d overlap and trac2d decomposition

     grid%jfirstct = grid%jfirst
     grid%jlastct  = grid%jlast
     grid%kfirstct = grid%kfirst
     grid%klastct  = grid%klast
     if (ct_overlap .gt. 0) then
        mlt = 2
     elseif (trac_decomp .gt. 1) then
        mlt = trac_decomp
     else
        mlt = 1
     endif

     if (mlt .gt. 1) then
        if (gid .lt. npes_yz) then
           ictstuff(1) = grid%jfirstct
           ictstuff(2) = grid%jlastct
           ictstuff(3) = grid%kfirstct
           ictstuff(4) = grid%klastct
           do ml = 2, mlt
              call mpisend(ictstuff, 4, mpiint, gid+(ml-1)*npes_yz, gid+(ml-1)*npes_yz, commnyz)
           enddo
        elseif (gid .lt. mlt*npes_yz) then
           ktmod = gid/npes_yz
           call mpirecv(ictstuff, 4, mpiint, gid-ktmod*npes_yz, gid, commnyz)
           grid%jfirstct = ictstuff(1)
           grid%jlastct  = ictstuff(2)
           grid%kfirstct = ictstuff(3)
           grid%klastct  = ictstuff(4)
        endif
     endif

     if (trac_decomp .gt. 1) then
        kquot = nq / trac_decomp
        krem = nq - kquot * trac_decomp
        krun = 0
        do kt = 1, krem
           grid%ktloa(kt) = krun + 1
           krun = krun + kquot + 1
           grid%kthia(kt) = krun
        enddo
        do kt = krem+1, trac_decomp
           grid%ktloa(kt) = krun + 1
           krun = krun + kquot
           grid%kthia(kt) = krun
        enddo
        ktmod = gid/npes_yz + 1
        ktmod = min(ktmod, trac_decomp)
        grid%ktlo = grid%ktloa(ktmod)
        grid%kthi = grid%kthia(ktmod)
     endif

   return
!EOC
end subroutine spmd_vars_init
!-----------------------------------------------------------------------
#endif

!EOC
   end subroutine dynamics_init
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dynamics_clean -- clean up Lin-Rood-specific variables
!
! !INTERFACE: 
   subroutine dynamics_clean(grid)

! !USES:
      implicit none

! !INPUT/OUTPUT PARAMETERS:
      type(T_FVDYCORE_GRID), intent(inout)  :: grid     ! Resulting grid


! !DESCRIPTION:
!
! Clean up (deallocate) Lin-Rood-specific variables
!
! !REVISION HISTORY:
!
!   01.06.06   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! Temporary data structures

    if(associated(GRID%SINLON          )) deallocate(GRID%SINLON)
    if(associated(GRID%COSLON          )) deallocate(GRID%COSLON)
    if(associated(GRID%SINL5           )) deallocate(GRID%SINL5)
    if(associated(GRID%COSL5           )) deallocate(GRID%COSL5)

    if(associated(GRID%ACOSP           )) deallocate(GRID%ACOSP)
    if(associated(GRID%ACOSU           )) deallocate(GRID%ACOSU)
    if(associated(GRID%SINP            )) deallocate(GRID%SINP)
    if(associated(GRID%COSP            )) deallocate(GRID%COSP)
    if(associated(GRID%SINE            )) deallocate(GRID%SINE)
    if(associated(GRID%COSE            )) deallocate(GRID%COSE)
    if(associated(GRID%AK              )) deallocate(GRID%AK)
    if(associated(GRID%BK              )) deallocate(GRID%BK)

!
! cd_core variables
!
    if(associated( grid%dtdx  )) deallocate(grid%dtdx)
    if(associated( grid%dtdx2 )) deallocate(grid%dtdx2)
    if(associated( grid%dtdx4 )) deallocate(grid%dtdx4)
    if(associated( grid%dtdxe )) deallocate(grid%dtdxe)
    if(associated( grid%dxdt  )) deallocate(grid%dxdt)
    if(associated( grid%dxe   )) deallocate(grid%dxe)
    if(associated( grid%cye   )) deallocate(grid%cye)
    if(associated( grid%dycp  )) deallocate(grid%dycp)
    if(associated( grid%rdxe  )) deallocate(grid%rdxe)
    if(associated( grid%txe5  )) deallocate(grid%txe5)
    if(associated( grid%dtxe5 )) deallocate(grid%dtxe5)
    if(associated( grid%dyce  )) deallocate(grid%dyce)
    if(associated( grid%dx    )) deallocate(grid%dx)
    if(associated( grid%rdx   )) deallocate(grid%rdx)
    if(associated( grid%cy    )) deallocate(grid%cy)

    if(associated( grid%sc    )) deallocate(grid%sc)
    if(associated( grid%se    )) deallocate(grid%se)
    if(associated( grid%dc    )) deallocate(grid%dc)
    if(associated( grid%de    )) deallocate(grid%de)

    if(associated( grid%cdx   )) deallocate(grid%cdx)
    if(associated( grid%cdy   )) deallocate(grid%cdy)
    if(associated( grid%cdx4  )) deallocate(grid%cdx4)  
    if(associated( grid%cdy4  )) deallocate(grid%cdy4)  
    if(associated( grid%cdxde )) deallocate(grid%cdxde) 
    if(associated( grid%cdxdp )) deallocate(grid%cdxdp) 
    if(associated( grid%cdydp )) deallocate(grid%cdydp) 
    if(associated( grid%cdyde )) deallocate(grid%cdyde) 
    if(associated( grid%cdxdiv)) deallocate(grid%cdxdiv)
    if(associated( grid%cdydiv)) deallocate(grid%cdydiv)
    if(associated( grid%cdtau4)) deallocate(grid%cdtau4)

    if(associated( grid%scdiv4)) deallocate(grid%scdiv4)
    if(associated( grid%sediv4)) deallocate(grid%sediv4)
    if(associated( grid%dcdiv4)) deallocate(grid%dcdiv4)
    if(associated( grid%dediv4)) deallocate(grid%dediv4)


    if(associated( grid%f0    )) deallocate(grid%f0)
    if(associated( grid%fc    )) deallocate(grid%fc)

#if defined(SPMD)
   call spmd_vars_clean(grid)
#endif
   return
!EOC
   end subroutine dynamics_clean
!-----------------------------------------------------------------------

#if defined(SPMD)
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  spmd_vars_clean --- Clean the SPMD-related variables
!
! !INTERFACE:
subroutine spmd_vars_clean(grid)

! !USES:
   use parutilitiesmodule, only : parpatternfree
   implicit none

!------------------------------Commons----------------------------------

! !INPUT PARAMETERS:
    type (T_FVDYCORE_GRID), intent(inout) :: grid

! 
! !DESCRIPTION:
! 
!   {\bf Purpose:} Clean the SPMD related variables.
! 
! !REVISION HISTORY: 
!   02.11.08    Sawyer       Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!

      if ( grid%twod_decomp == 1 ) then
! Clean transposes
!
        call parpatternfree(grid%commxy, grid%u_to_uxy)
        call parpatternfree(grid%commxy, grid%uxy_to_u)
        call parpatternfree(grid%commxy, grid%v_to_vxy)
        call parpatternfree(grid%commxy, grid%vxy_to_v)
        call parpatternfree(grid%commxy, grid%ijk_yz_to_xy)
        call parpatternfree(grid%commxy, grid%ijk_xy_to_yz)
        call parpatternfree(grid%commxy, grid%ikj_xy_to_yz)
        call parpatternfree(grid%commxy, grid%ikj_yz_to_xy)
        call parpatternfree(grid%commxy, grid%pe_to_pexy)
        call parpatternfree(grid%commxy, grid%pexy_to_pe)
        call parpatternfree(grid%commxy, grid%pt_to_ptxy)
        call parpatternfree(grid%commxy, grid%ptxy_to_pt)
        call parpatternfree(grid%commxy, grid%r4_xy_to_yz)
        call parpatternfree(grid%commxy, grid%r4_yz_to_xy)
        call parpatternfree(grid%commxy, grid%pkxy_to_pkc)
        call parpatternfree(grid%commxy, grid%xy2d_to_yz2d)
        call parpatternfree(grid%commxy, grid%yz2d_to_xy2d)
      endif
   return
!EOC
end subroutine spmd_vars_clean
#endif

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: a2d3d -- 2nd order A-to-D grid transform (3D) XY decomp
!                      INOUT array is i,j,k, and is modified in place
!
! !INTERFACE:

      subroutine a2d3d( grid, u, v )

! !USES:

#if defined( SPMD )
      use mod_comm, only: mp_send3d, mp_recv3d
#endif

      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: u(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(inout) :: v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !DESCRIPTION:
!
!     This routine performs a second order interpolation of 
!     three-dimensional wind fields on a A grid to an D grid.
!     In place calculation!
!
! !REVISION HISTORY:
!     WS  03.08.27 : Creation from d2a3d
!     WS  03.10.22 : pmgrid removed (now spmd_dyn)
!     WS  04.08.25 : simplified interfaces with grid (only for XY!!!)
!     WS  04.10.06 : removed spmd_dyn, all those vars. now from grid
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer  :: im      ! Dimensions longitude (total)
      integer  :: jm      ! Dimensions latitude (total)
      integer  :: km      ! Dimensions vertical (total)
      integer  :: ifirst  ! longitude strip start
      integer  :: ilast   ! longitude strip finish
      integer  :: jfirst  ! latitude strip start
      integer  :: jlast   ! latitude strip finish
      integer  :: iam     ! process identifier
      integer  :: myidxy_y, myidxy_x, nprxy_x
      integer  :: commxy

      real(r8), parameter :: UNDEFINED = 1.0D15

      integer  :: i, j, k, itot, jtot
      real(r8) :: vwest(grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) :: usouth(grid%ifirstxy:grid%ilastxy,grid%km)

#if defined( SPMD )
      integer dest, src
#endif

      im     = grid%im
      jm     = grid%jm
      km     = grid%km
      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      jfirst = grid%jfirstxy
      jlast  = grid%jlastxy

      iam      = grid%iam
      myidxy_x = grid%myidxy_x
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x

      commxy   = grid%commxy

      itot = ilast-ifirst+1
      jtot = jlast-jfirst+1

#if defined( SPMD )
! Send one latitude to the north
      call mp_send3d( commxy, iam+nprxy_x, iam-nprxy_x, im, jm, km,    &
                      ifirst, ilast, jfirst, jlast, 1, km,             &
                      ifirst, ilast, jlast, jlast, 1, km, u )
      call mp_recv3d( commxy, iam-nprxy_x, im, jm, km,                 &
                      ifirst, ilast, jfirst-1, jfirst-1, 1, km,        &
                      ifirst, ilast, jfirst-1, jfirst-1, 1, km, usouth )
#endif

!$omp  parallel do private(i,j,k)
      do k=1,km
        do j=jlast, jfirst+1, -1
          do i=ifirst,ilast
            u(i,j,k) = D0_5*(u(i,j-1,k) + u(i,j,k))
          enddo
        enddo
      enddo

#if defined( SPMD )
      if ( jfirst > 1 ) then
!$omp  parallel do private(i, k)
         do k=1,km
            do i=ifirst,ilast
              u(i,jfirst,k) = D0_5 * ( u(i,jfirst,k) + usouth(i,k) )
            enddo
         enddo
      endif
#endif

      if ( jfirst ==  1 ) then
!$omp  parallel do private(i,k)
         do k=1,km
            do i=ifirst,ilast
               u(i,1,k) = UNDEFINED
            enddo
         enddo
      endif

!
! V-winds
!

! Pack vwest with wrap-around condition

!$omp  parallel do private(j,k)
      do k = 1,km
         do j=jfirst,jlast
            vwest(j,k) = v(ilast,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (itot /= im) then
         dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         call mp_send3d( commxy, dest, src, im, jm, km,                   &
                         ifirst, ilast, jfirst, jlast, 1, km,             &
                         ilast, ilast, jfirst, jlast, 1, km, v )
         call mp_recv3d( commxy, src, im, jm, km,                         & 
                         ifirst-1, ifirst-1, jfirst, jlast, 1, km,        &
                         ifirst-1, ifirst-1, jfirst, jlast, 1, km, vwest )
      endif
#endif

!
! Beware: ilast is en route, don't alter its value
!

!$omp  parallel do private(i,j,k)
      do k=1,km
         do j=jfirst, jlast
            do i=ilast,ifirst+1,-1
               v(i,j,k) = D0_5*(v(i-1,j,k) + v(i,j,k))
            enddo
         enddo
      enddo
!
! Clean up shop
!

!$omp  parallel do private(i,j,k)
      do k=1,km
         do j=jfirst, jlast
            v(ifirst,j,k)= D0_5*(vwest(j,k) + v(ifirst,j,k))
         enddo
      enddo

      return
!EOC
      end subroutine a2d3d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: d2a3d -- 2nd order D-to-A grid transform (3D) XY decomp.
!                    Output array is i,j,k
!
! !INTERFACE:

      subroutine d2a3d( grid, u, v, ua, va )

! !USES:

#if defined( SPMD )
      use parutilitiesmodule, only : parcollective3d, sumop
      use mod_comm, only: mp_send3d, mp_recv3d
#endif

      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      real(r8), intent(in) :: u(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(in) :: v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ua(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(inout) :: va(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind


! !DESCRIPTION:
!
!     This routine performs a second order 
!     interpolation of three-dimensional wind
!     fields on a D grid to an A grid.  Only for an XY decomposition!
!
! !REVISION HISTORY:
!     WS  00.12.22 : Creation from d2a3d
!     AAM 01.06.13 : Generalized to 2D decomposition
!     WS  02.04.25 : Newest mod_comm interfaces
!     WS  03.08.27 : Minimal alterations to interface, renamed d2a3d
!     WS  03.10.22 : pmgrid removed (now spmd_dyn)
!     WS  04.08.25 : simplified interfaces with grid (only for XY!!!)
!     WS  04.10.06 : removed spmd_dyn, all those vars. now from grid
!
!EOP
!-----------------------------------------------------------------------
!BOC
      integer  :: im        ! Dimensions longitude (total)
      integer  :: jm        ! Dimensions latitude (total)
      integer  :: km        ! Dimensions level (total)
      integer  :: ifirst    ! longitude strip start
      integer  :: ilast     ! longitude strip finish
      integer  :: jfirst    ! latitude strip start
      integer  :: jlast     ! latitude strip finish
      integer  :: iam, myidxy_y, nprxy_x, commxy, commxy_x

      real(r8), pointer :: coslon(:) ! Cosine in longitude
      real(r8), pointer :: sinlon(:) ! Sine in longitude

      integer imh, i, j, k, itot, jtot, ltot, lbegin, lend, ik

      real(r8) :: un(grid%km), vn(grid%km), us(grid%km), vs(grid%km)
      real(r8) :: veast(grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) :: unorth(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) :: uvaglob(grid%im,grid%km,4)
      real(r8) :: uvaloc(grid%ifirstxy:grid%ilastxy,grid%km,4)
      real(r8) :: uaglob(grid%im),vaglob(grid%im)

#if defined( SPMD )
      integer dest, src, incount, outcount
#endif

!
! Retrieve values from grid
!
      im     = grid%im
      jm     = grid%jm
      km     = grid%km
      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      jfirst = grid%jfirstxy
      jlast  = grid%jlastxy

      iam      = grid%iam
      nprxy_x  = grid%nprxy_x
      commxy   = grid%commxy
      commxy_x = grid%commxy_x
      myidxy_y = grid%myidxy_y

      coslon =>grid%coslon
      sinlon =>grid%sinlon

      itot = ilast-ifirst+1
      jtot = jlast-jfirst+1

      imh = im/2

#if defined( SPMD )
! Set ua on A-grid
      call mp_send3d( commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,          &
                      ifirst, ilast, jfirst, jlast, 1, km,                   &
                      ifirst, ilast, jfirst, jfirst, 1, km, u )
      call mp_recv3d( commxy, iam+nprxy_x, im, jm, km,                       &
                      ifirst, ilast, jlast+1, jlast+1, 1, km,                &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, unorth )

      if ( jlast .lt. jm ) then
!$omp  parallel do private(i, k)

         do k=1,km
            do i=ifirst,ilast
               ua(i,jlast,k) = D0_5 * ( u(i,jlast,k) + unorth(i,k) )
            enddo
         enddo
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
        do j=jfirst, jlast-1
          do i=ifirst,ilast
            ua(i,j,k) = D0_5*(u(i,j,k) + u(i,j+1,k))
          enddo
        enddo
      enddo

! Set va on A-grid

!$omp  parallel do private(j,k)

      do k = 1,km
         do j=jfirst,jlast
            veast(j,k) = v(ifirst,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (itot .ne. im) then
         dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         call mp_send3d( commxy, dest, src, im, jm, km,                           &
                         ifirst, ilast, jfirst, jlast, 1, km,                     &
                         ifirst, ifirst, jfirst, jlast, 1, km, v )
         call mp_recv3d( commxy, src, im, jm, km,                                 & 
                         ilast+1, ilast+1, jfirst, jlast, 1, km,                  &
                         ilast+1, ilast+1, jfirst, jlast, 1, km, veast )
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
         do j=jfirst, jlast
            do i=ifirst,ilast-1
               va(i,j,k) = D0_5*(v(i,j,k) + v(i+1,j,k))
            enddo
            va(ilast,j,k) = D0_5*(v(ilast,j,k) + veast(j,k))
         enddo
      enddo

!$omp  parallel do private(i,ik,k)

      do ik=1,4
         do k=1,km
            do i=1,im
               uvaglob(i,k,ik) = D0_0
            enddo
         enddo
      enddo

      if (jfirst .eq. 1) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirst,ilast
               uvaloc(i,k,1) = ua(i,2,k)
               uvaloc(i,k,2) = va(i,2,k)
               uvaglob(i,k,1) = ua(i,2,k)
               uvaglob(i,k,2) = va(i,2,k)
            enddo
         enddo
         lbegin = 1
         lend = 2
      endif

      if (jlast .eq. jm) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirst,ilast
               uvaloc(i,k,3) = ua(i,jm-1,k)
               uvaloc(i,k,4) = va(i,jm-1,k)
               uvaglob(i,k,3) = ua(i,jm-1,k)
               uvaglob(i,k,4) = va(i,jm-1,k)
            enddo
         enddo
         lbegin = 3
         lend = 4
      endif
      if (jtot .eq. jm) lbegin=1

#if defined( SPMD )
      if (itot .ne. im) then
         ltot = lend-lbegin+1
         if (jfirst .eq. 1 .or. jlast .eq. jm) then
            call parcollective3d(commxy_x, sumop, im, km, ltot, uvaglob(1,1,lbegin))
         endif
      endif
#endif

      if ( jfirst .eq. 1 ) then
! Projection at SP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,km
            us(k) = D0_0
            vs(k) = D0_0
            do i=1,imh
               us(k) = us(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*sinlon(i)  &
                     + (uvaglob(i,k,2)-uvaglob(i+imh,k,2))*coslon(i)
               vs(k) = vs(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*coslon(i)  & 
                     + (uvaglob(i+imh,k,2)-uvaglob(i,k,2))*sinlon(i)
            enddo

            us(k) = us(k)/im
            vs(k) = vs(k)/im
            do i=1,imh
               uaglob(i)   = -us(k)*sinlon(i) - vs(k)*coslon(i)
               vaglob(i)   =  us(k)*coslon(i) - vs(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirst,ilast
               ua(i,1,k) = uaglob(i)
               va(i,1,k) = vaglob(i)
            enddo
         enddo
      endif

      if ( jlast .eq. jm ) then
! Projection at NP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,km
            un(k) = D0_0
            vn(k) = D0_0
            do i=1,imh
               un(k) = un(k) + (uvaglob(i+imh,k,3)-uvaglob(i,k,3))*sinlon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*coslon(i)
               vn(k) = vn(k) + (uvaglob(i,k,3)-uvaglob(i+imh,k,3))*coslon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*sinlon(i)
            enddo

            un(k) = un(k)/im
            vn(k) = vn(k)/im
            do i=1,imh
               uaglob(i) = -un(k)*sinlon(i) + vn(k)*coslon(i)
               vaglob(i) = -un(k)*coslon(i) - vn(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirst,ilast
               ua(i,jm,k) = uaglob(i)
               va(i,jm,k) = vaglob(i)
            enddo
         enddo
      endif

      return
!EOC
      end subroutine d2a3d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: d2b3d -- 2nd order D-to-B grid transform (3D) XY decomp.
!                    Output array is i,j,k
!
! !INTERFACE:

      subroutine d2b3d( grid, u, v, ub, vb )

! !USES:
#if defined( SPMD )
      use mod_comm, only: mp_send3d, mp_recv3d
#endif

      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      real(r8), intent(in) :: u(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(in) :: v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ub(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(inout) :: vb(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind


! !DESCRIPTION:
!
!     This routine performs a second order
!     interpolation of three-dimensional wind
!     fields on a D grid to an B grid.  Only for an XY decomposition!
!
! !REVISION HISTORY:
!     BP  05.02.22 : Creation from d2a3d
!
!EOP
!-----------------------------------------------------------------------
!BOC
      integer  :: im        ! Dimensions longitude (total)
      integer  :: jm        ! Dimensions latitude (total)
      integer  :: km        ! Dimensions level (total)
      integer  :: ifirst    ! longitude strip start
      integer  :: ilast     ! longitude strip finish
      integer  :: jfirst    ! latitude strip start
      integer  :: jlast     ! latitude strip finish
      integer  :: iam, myidxy_y, nprxy_x, commxy

      real(r8), parameter :: UNDEFINED = 1.0D15


      real(r8), pointer :: coslon(:) ! Cosine in longitude
      real(r8), pointer :: sinlon(:) ! Sine in longitude

      integer imh, i, j, k, itot, jtot, ltot, lbegin, lend, ik

      real(r8) :: ueast(grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) :: vsouth(grid%ifirstxy:grid%ilastxy,grid%km)

#if defined( SPMD )
      integer dest, src, incount, outcount
#endif

!
! Retrieve values from grid
!
      im     = grid%im
      jm     = grid%jm
      km     = grid%km
      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      jfirst = grid%jfirstxy
      jlast  = grid%jlastxy

      iam      = grid%iam
      nprxy_x  = grid%nprxy_x
      commxy   = grid%commxy
      myidxy_y = grid%myidxy_y

      coslon =>grid%coslon
      sinlon =>grid%sinlon

      itot = ilast-ifirst+1
      jtot = jlast-jfirst+1

      imh = im/2

#if defined( SPMD )
! Set vb on B-grid
      call mp_send3d( commxy, iam+nprxy_x, iam-nprxy_x, im, jm, km,          &
                      ifirst, ilast, jfirst, jlast, 1, km, &
                      ifirst, ilast, jfirst, jfirst, 1, km, v )
      call mp_recv3d( commxy, iam-nprxy_x, im, jm, km,                       &
                      ifirst, ilast, jlast+1, jlast+1, 1, km,        &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, vsouth )

      if ( jfirst .gt. 1 ) then
!$omp  parallel do private(i, k)

         do k=1,km
            do i=ifirst,ilast
               vb(i,jfirst,k) = D0_5 * ( v(i,jfirst,k) + vsouth(i,k) )
            enddo
         enddo
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
        do j=jfirst+1, jlast
          do i=ifirst,ilast
            vb(i,j,k) = D0_5*(v(i,j,k) + v(i,j-1,k))
          enddo
        enddo
      enddo

! Set ub on B-grid

!$omp  parallel do private(j,k)

      do k = 1,km
         do j=jfirst,jlast
            ueast(j,k) = u(ifirst,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (itot .ne. im) then
         dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         call mp_send3d( commxy, dest, src, im, jm, km,                   &
                         ifirst, ilast, jfirst, jlast, 1, km,             &
                         ifirst, ifirst, jfirst, jlast, 1, km, u )
         call mp_recv3d( commxy, src, im, jm, km,                         &
                         ilast+1, ilast+1, jfirst, jlast, 1, km,          &
                         ilast+1, ilast+1, jfirst, jlast, 1, km, ueast )
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
         do j=jfirst, jlast
            do i=ifirst,ilast-1
               ub(i,j,k) = D0_5*(u(i,j,k) + u(i+1,j,k))
            enddo
            ub(ilast,j,k) = D0_5*(u(ilast,j,k) + ueast(j,k))
         enddo
      enddo

      if ( jfirst == 1 ) then
!$omp  parallel do private(i,k)
         do k=1,km
            do i=ifirst,ilast
               ub(i,1,k) = UNDEFINED
               vb(i,1,k) = UNDEFINED
            enddo
         enddo
      endif

      return
!EOC
      end subroutine d2b3d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: b2d3d -- 2nd order B-to-D grid transform (3D) XY decomp
!                      INOUT array is i,j,k, and is modified in place
!
! !INTERFACE:

      subroutine b2d3d( grid, u, v )

! !USES:
#if defined( SPMD )
      use parutilitiesmodule, only : parcollective3d, sumop, gid
      use mod_comm, only: mp_send3d, mp_recv3d
#endif

      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: u(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(inout) :: v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !DESCRIPTION:
!
!     This routine performs a second order interpolation of
!     three-dimensional wind fields on a B grid to an D grid.
!     In place calculation!
!
! !REVISION HISTORY:
!     BP  05.02.22 : Creation from a2d3d
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer  :: im      ! Dimensions longitude (total)
      integer  :: jm      ! Dimensions latitude (total)
      integer  :: km      ! Dimensions vertical (total)
      integer  :: ifirst  ! longitude strip start
      integer  :: ilast   ! longitude strip finish
      integer  :: jfirst  ! latitude strip start
      integer  :: jlast   ! latitude strip finish
      integer  :: iam     ! process identifier
      integer  :: myidxy_y, myidxy_x, nprxy_x
      integer  :: commxy, commxy_x

      real(r8), parameter :: UNDEFINED = 1.0D15

      real(r8) :: uwest(grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) :: vsouth(grid%ifirstxy:grid%ilastxy,grid%km)

      integer  :: imh, i, j, k, itot, jtot, ltot, lbegin, lend, ik
      real(r8) :: un(grid%km), vn(grid%km), us(grid%km), vs(grid%km)
      real(r8) :: uvbglob(grid%im,grid%km,4)
      real(r8) :: uvbloc(grid%ifirstxy:grid%ilastxy,grid%km,4)
      real(r8) :: ubglob(grid%im),vbglob(grid%im)

      real(r8), pointer :: coslon(:) ! Cosine in longitude
      real(r8), pointer :: sinlon(:) ! Sine in longitude
      real(r8), pointer :: cosl5(:)   ! Cosine in longitude
      real(r8), pointer :: sinl5(:)   ! Sine in longitude

#if defined( SPMD )
      integer dest, src
#endif

      im     = grid%im
      jm     = grid%jm
      km     = grid%km
      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      jfirst = grid%jfirstxy
      jlast  = grid%jlastxy

      iam      = grid%iam
      myidxy_x = grid%myidxy_x
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x

      commxy   = grid%commxy
      commxy_x = grid%commxy_x

      itot = ilast-ifirst+1
      jtot = jlast-jfirst+1

      imh = im/2
      coslon => grid%coslon
      sinlon => grid%sinlon
      cosl5 => grid%cosl5
      sinl5 => grid%sinl5

!
! Initial Preparation for Projection at Poles
!
!$omp  parallel do private(i,ik,k)

      do ik=1,4
         do k=1,km
            do i=1,im
               uvbglob(i,k,ik) = D0_0
            enddo
         enddo
      enddo

      if (jfirst .eq. 1) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirst,ilast
               uvbloc(i,k,1) = u(i,2,k)
               uvbloc(i,k,2) = v(i,2,k)
               uvbglob(i,k,1) = u(i,2,k)
               uvbglob(i,k,2) = v(i,2,k)
            enddo
         enddo
         lbegin = 1
         lend = 2
      endif

      if (jlast .eq. jm) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirst,ilast
               uvbloc(i,k,3) = u(i,jm,k)
               uvbloc(i,k,4) = v(i,jm,k)
               uvbglob(i,k,3) = u(i,jm,k)
               uvbglob(i,k,4) = v(i,jm,k)
            enddo
         enddo
         lbegin = 3
         lend = 4
      endif
      if (jtot .eq. jm) lbegin=1

#if defined( SPMD )
      if (itot .ne. im) then
         ltot = lend-lbegin+1
         if (jfirst .eq. 1 .or. jlast .eq. jm) then
            call parcollective3d(commxy_x, sumop, im, km, ltot, uvbglob(1,1,lbegin))
         endif
      endif
#endif

!
! V-Winds
!

#if defined( SPMD )
! Send one latitude to the north
      call mp_send3d( commxy, iam+nprxy_x, iam-nprxy_x, im, jm, km,          &
                      ifirst, ilast, jfirst, jlast, 1, km, &
                      ifirst, ilast, jlast, jlast, 1, km, v )
      call mp_recv3d( commxy, iam-nprxy_x, im, jm, km,                       &
                      ifirst, ilast, jfirst-1, jfirst-1, 1, km,        &
                      ifirst, ilast, jfirst-1, jfirst-1, 1, km, vsouth )
#endif

!$omp  parallel do private(i,j,k)
      do k=1,km
        do j=jlast, jfirst+1, -1
          do i=ifirst,ilast
            v(i,j,k) = D0_5*(v(i,j-1,k) + v(i,j,k))
          enddo
        enddo
      enddo

#if defined( SPMD )
      if ( jfirst > 1 ) then
!$omp  parallel do private(i, k)
         do k=1,km
            do i=ifirst,ilast
              v(i,jfirst,k) = D0_5 * ( v(i,jfirst,k) + vsouth(i,k) )
            enddo
         enddo
      endif
#endif

!
! U-winds
!

! Pack uwest with wrap-around condition

!$omp  parallel do private(j,k)
      do k = 1,km
         do j=jfirst,jlast
            uwest(j,k) = v(ilast,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (itot /= im) then
         dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         call mp_send3d( commxy, dest, src, im, jm, km,                   &
                         ifirst, ilast, jfirst, jlast, 1, km,             &
                         ilast, ilast, jfirst, jlast, 1, km, u )
         call mp_recv3d( commxy, src, im, jm, km,                         &
                         ifirst-1, ifirst-1, jfirst, jlast, 1, km,        &
                         ifirst-1, ifirst-1, jfirst, jlast, 1, km, uwest )
      endif
#endif

!$omp  parallel do private(i,j,k)
      do k=1,km
         do j=jfirst, jlast
            do i=ilast,ifirst+1,-1
               u(i,j,k) = D0_5*(u(i-1,j,k) + u(i,j,k))
            enddo
         enddo
      enddo

!$omp  parallel do private(i,j,k)
      do k=1,km
         do j=jfirst, jlast
            u(ifirst,j,k)= D0_5*(uwest(j,k) + u(ifirst,j,k))
         enddo
      enddo

      if ( jfirst ==  1 ) then
!$omp  parallel do private(i,k)
         do k=1,km
            do i=ifirst,ilast
               u(i,1,k) = UNDEFINED
            enddo
         enddo
      endif

!
! Project V-Winds to the Poles
!
      if ( jfirst == 1 ) then
! Projection at SP
!$omp  parallel do private(i,k,ubglob,vbglob)
         do k=1,km
            us(k) = D0_0
            vs(k) = D0_0
            do i=1,imh
               us(k) = us(k) + (uvbglob(i+imh,k,1)-uvbglob(i,k,1))*sinlon(i)  &
                     + (uvbglob(i,k,2)-uvbglob(i+imh,k,2))*coslon(i)
               vs(k) = vs(k) + (uvbglob(i+imh,k,1)-uvbglob(i,k,1))*coslon(i)  &
                     + (uvbglob(i+imh,k,2)-uvbglob(i,k,2))*sinlon(i)
            enddo

            us(k) = us(k)/im
            vs(k) = vs(k)/im
            do i=1,imh
               vbglob(i)   =  us(k)*cosl5(i) - vs(k)*sinl5(i)
               vbglob(i+imh) = -vbglob(i)
            enddo
            do i=ifirst,ilast
               v(i,1,k) = vbglob(i)
            enddo
         enddo
      endif

      if ( jlast == jm ) then
! Projection at NP
!$omp  parallel do private(i,k,ubglob,vbglob)
         do k=1,km
            un(k) = D0_0
            vn(k) = D0_0
            do i=1,imh
               un(k) = un(k) + (uvbglob(i+imh,k,3)-uvbglob(i,k,3))*sinlon(i) &
                     + (uvbglob(i+imh,k,4)-uvbglob(i,k,4))*coslon(i)
               vn(k) = vn(k) + (uvbglob(i,k,3)-uvbglob(i+imh,k,3))*coslon(i) &
                     + (uvbglob(i+imh,k,4)-uvbglob(i,k,4))*sinlon(i)
            enddo

            un(k) = un(k)/im
            vn(k) = vn(k)/im
            do i=1,imh
               vbglob(i) = -un(k)*cosl5(i) - vn(k)*sinl5(i)
               vbglob(i+imh) = -vbglob(i)
            enddo
            do i=ifirst,ilast
               v(i,jm,k) = vbglob(i)
            enddo
         enddo
      endif

      return
!EOC
      end subroutine b2d3d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: c2a3d -- 2nd order C-to-A grid transform (3D) XY decomp.
!                    Output array is i,j,k
!
! !INTERFACE:

      subroutine c2a3d( grid, u, v, ua, va ) 
                             
! !USES:                     
                             
#if defined( SPMD )          
      use parutilitiesmodule, only : parcollective3d, sumop          
      use mod_comm, only: mp_send3d, mp_recv3d
#endif                       
                             
      implicit none          
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid                     
      real(r8), intent(in) :: u(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(in) :: v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ua(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(inout) :: va(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !DESCRIPTION:
!
!     This routine performs a second order              
!     interpolation of three-dimensional wind           
!     fields on a C grid to an A grid.  Only for an XY decomposition!
!     
! !REVISION HISTORY:
!     WMP  06.11.03 : Creation from d2a3d
!
!EOP
!-----------------------------------------------------------------------
!BOC
      integer  :: im        ! Dimensions longitude (total)
      integer  :: jm        ! Dimensions latitude (total)
      integer  :: km        ! Dimensions level (total)
      integer  :: ifirst    ! longitude strip start
      integer  :: ilast     ! longitude strip finish
      integer  :: jfirst    ! latitude strip start
      integer  :: jlast     ! latitude strip finish
      integer  :: iam, myidxy_y, nprxy_x, commxy_x, commxy

      real(r8), pointer :: coslon(:) ! Cosine in longitude
      real(r8), pointer :: sinlon(:) ! Sine in longitude

      integer imh, i, j, k, itot, jtot, ltot, lbegin, lend, ik

      real(r8) :: un(grid%km), vn(grid%km), us(grid%km), vs(grid%km)
      real(r8) :: ueast(grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) :: vnorth(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) :: uvaglob(grid%im,grid%km,4)
      real(r8) :: uvaloc(grid%ifirstxy:grid%ilastxy,grid%km,4)
      real(r8) :: uaglob(grid%im),vaglob(grid%im)

#if defined( SPMD )
      integer dest, src, incount, outcount
#endif

!
! Retrieve values from grid
!
      im     = grid%im
      jm     = grid%jm
      km     = grid%km
      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      jfirst = grid%jfirstxy
      jlast  = grid%jlastxy

      iam      = grid%iam
      nprxy_x  = grid%nprxy_x
      commxy   = grid%commxy
      commxy_x = grid%commxy_x
      myidxy_y = grid%myidxy_y

      coslon =>grid%coslon
      sinlon =>grid%sinlon

      itot = ilast-ifirst+1
      jtot = jlast-jfirst+1

      imh = im/2

#if defined( SPMD )
! Set va on A-grid
      call mp_send3d( commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,          &
                      ifirst, ilast, jfirst, jlast, 1, km, &
                      ifirst, ilast, jfirst, jfirst, 1, km, v )
      call mp_recv3d( commxy, iam+nprxy_x, im, jm, km,                       &
                      ifirst, ilast, jlast+1, jlast+1, 1, km,        &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, vnorth )

      if ( jlast .lt. jm ) then
!$omp  parallel do private(i, k)

         do k=1,km
            do i=ifirst,ilast
               va(i,jlast,k) = D0_5 * ( v(i,jlast,k) + vnorth(i,k) )
            enddo
         enddo
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
        do j=jfirst, jlast-1
          do i=ifirst,ilast
            va(i,j,k) = D0_5*(v(i,j,k) + v(i,j+1,k))
          enddo
        enddo
      enddo

! Set ua on A-grid

!$omp  parallel do private(j,k)

      do k = 1,km
         do j=jfirst,jlast
            ueast(j,k) = u(ifirst,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (itot .ne. im) then
         dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         call mp_send3d( commxy, dest, src, im, jm, km,                           &
                         ifirst, ilast, jfirst, jlast, 1, km,   &
                         ifirst, ifirst, jfirst, jlast, 1, km, u )
         call mp_recv3d( commxy, src, im, jm, km,                                 &
                         ilast+1, ilast+1, jfirst, jlast, 1, km,          &
                         ilast+1, ilast+1, jfirst, jlast, 1, km, ueast )
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
         do j=jfirst, jlast
            do i=ifirst,ilast-1
               ua(i,j,k) = D0_5*(u(i,j,k) + u(i+1,j,k))
            enddo
            ua(ilast,j,k) = D0_5*(u(ilast,j,k) + ueast(j,k))
         enddo
      enddo

!$omp  parallel do private(i,ik,k)

      do ik=1,4
         do k=1,km
            do i=1,im
               uvaglob(i,k,ik) = D0_0
            enddo
         enddo
      enddo

      if (jfirst .eq. 1) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirst,ilast
               uvaloc(i,k,1) = ua(i,2,k)
               uvaloc(i,k,2) = va(i,2,k)
               uvaglob(i,k,1) = ua(i,2,k)
               uvaglob(i,k,2) = va(i,2,k)
            enddo
         enddo
         lbegin = 1
         lend = 2
      endif

      if (jlast .eq. jm) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirst,ilast
               uvaloc(i,k,3) = ua(i,jm-1,k)
               uvaloc(i,k,4) = va(i,jm-1,k)
               uvaglob(i,k,3) = ua(i,jm-1,k)
               uvaglob(i,k,4) = va(i,jm-1,k)
            enddo
         enddo
         lbegin = 3
         lend = 4
      endif
      if (jtot .eq. jm) lbegin=1

#if defined( SPMD )
      if (itot .ne. im) then
         ltot = lend-lbegin+1
         if (jfirst .eq. 1 .or. jlast .eq. jm) then
            call parcollective3d(commxy_x, sumop, im, km, ltot, uvaglob(1,1,lbegin))
         endif
      endif
#endif

      if ( jfirst .eq. 1 ) then
! Projection at SP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,km
            us(k) = D0_0
            vs(k) = D0_0
            do i=1,imh
               us(k) = us(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*sinlon(i)  &
                     + (uvaglob(i,k,2)-uvaglob(i+imh,k,2))*coslon(i)
               vs(k) = vs(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*coslon(i)  &
                     + (uvaglob(i+imh,k,2)-uvaglob(i,k,2))*sinlon(i)
            enddo

            us(k) = us(k)/im
            vs(k) = vs(k)/im
            do i=1,imh
               uaglob(i)   = -us(k)*sinlon(i) - vs(k)*coslon(i)
               vaglob(i)   =  us(k)*coslon(i) - vs(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirst,ilast
               ua(i,1,k) = uaglob(i)
               va(i,1,k) = vaglob(i)
            enddo
         enddo
      endif

      if ( jlast .eq. jm ) then
! Projection at NP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,km
            un(k) = D0_0
            vn(k) = D0_0
            do i=1,imh
               un(k) = un(k) + (uvaglob(i+imh,k,3)-uvaglob(i,k,3))*sinlon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*coslon(i)
               vn(k) = vn(k) + (uvaglob(i,k,3)-uvaglob(i+imh,k,3))*coslon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*sinlon(i)
            enddo

            un(k) = un(k)/im
            vn(k) = vn(k)/im
            do i=1,imh
               uaglob(i) = -un(k)*sinlon(i) + vn(k)*coslon(i)
               vaglob(i) = -un(k)*coslon(i) - vn(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirst,ilast
               ua(i,jm,k) = uaglob(i)
               va(i,jm,k) = vaglob(i)
            enddo
         enddo
      endif

      return
!EOC
      end subroutine c2a3d
!-----------------------------------------------------------------------

end module dynamics_vars

