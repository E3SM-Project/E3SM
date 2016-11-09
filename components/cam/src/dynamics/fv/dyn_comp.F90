module dyn_comp

use shr_kind_mod,           only: r8 => shr_kind_r8, r4 => shr_kind_r4
use dynamics_vars,          only: T_FVDYCORE_GRID,            &
                                  T_FVDYCORE_STATE, T_FVDYCORE_CONSTANTS
use cam_abortutils,         only: endrun

#if defined(SPMD)
use mpishorthand,           only: mpicom, mpir8
#endif

use perf_mod
use cam_logfile,            only: iulog
use hycoef,                 only: hycoef_init, hyai, hybi
use pio,                    only: file_desc_t

implicit none
private

public dyn_init, dyn_run, dyn_final

public dyn_import_t, dyn_export_t, dyn_state

type (T_FVDYCORE_STATE), save, target :: dyn_state ! to be moved up later

type dyn_import_t
     real(r8), dimension(:,: ),    pointer     :: phis   ! Surface geopotential
     real(r8), dimension(:,: ),    pointer     :: ps     ! Surface pressure
     real(r8), dimension(:,:,:  ), pointer     :: u3s    ! U-winds (staggered)
     real(r8), dimension(:,:,:  ), pointer     :: v3s    ! V-winds (staggered)
     real(r8), dimension(:,:,:  ), pointer     :: pe     ! Pressure
     real(r8), dimension(:,:,:  ), pointer     :: pt     ! Potential temperature
     real(r8), dimension(:,:,:  ), pointer     :: t3     ! Temperatures
     real(r8), dimension(:,:,:  ), pointer     :: pk     ! Pressure to the kappa
     real(r8), dimension(:,:,:  ), pointer     :: pkz    ! Pressure to the kappa offset
     real(r8), dimension(:,:,:  ), pointer     :: delp   ! Delta pressure
     real(r8), dimension(:,:,:,:), pointer     :: tracer ! Tracers
end type dyn_import_t

type dyn_export_t
     real(r8), dimension(:,: ),    pointer     :: phis   ! Surface geopotential
     real(r8), dimension(:,: ),    pointer     :: ps     ! Surface pressure
     real(r8), dimension(:,:,:  ), pointer     :: u3s    ! U-winds (staggered)
     real(r8), dimension(:,:,:  ), pointer     :: v3s    ! V-winds (staggered)
     real(r8), dimension(:,:,:  ), pointer     :: pe     ! Pressure
     real(r8), dimension(:,:,:  ), pointer     :: pt     ! Potential temperature
     real(r8), dimension(:,:,:  ), pointer     :: t3     ! Temperatures
     real(r8), dimension(:,:,:  ), pointer     :: pk     ! Pressure to the kappa
     real(r8), dimension(:,:,:  ), pointer     :: pkz    ! Pressure to the kappa offset
     real(r8), dimension(:,:,:  ), pointer     :: delp   ! Delta pressure
     real(r8), dimension(:,:,:,:), pointer     :: tracer ! Tracers
     real(r8), dimension(:,:,:  ), pointer     :: peln   !
     real(r8), dimension(:,:,:  ), pointer     :: omga   ! Vertical velocity
     real(r8), dimension(:,:,:  ), pointer     :: mfx    ! Mass flux in X
     real(r8), dimension(:,:,:  ), pointer     :: mfy    ! Mass flux in Y
end type dyn_export_t

! !DESCRIPTION: This module implements the FVCAM Dynamical Core as
!               an ESMF gridded component.  It is specific to FVCAM
!               and does not use ESMF.
!
! \paragraph{Overview}
!
!   This module contains an ESMF wrapper for the Finite-Volume
!   Dynamical Core used in the Community Atmospheric Model
!   (FVCAM). This component will hereafter be referred
!   to as the ``FVdycore'' ESMF gridded component.  FVdycore
!   consists of four sub-components,
!
!   \begin{itemize}
!      \item {\tt cd\_core:}  The C/D-grid dycore component
!      \item {\tt te\_map:}   Vertical remapping algorithm
!      \item {\tt trac2d:}    Tracer advection
!      \item {\tt benergy:}   Energy balance
!   \end{itemize}
!
!   Subsequently the ESMF component design for FV dycore
!   will be described.   
!
! \paragraph{Internal State}
!
!  FVdycore maintains an internal state consisting of the
!  following fields:  control variables
!
!   \begin{itemize}
!     \item {\tt U}:    U winds on a D-grid (m/s)
!     \item {\tt V}:    V winds on a D-grid (m/s)
!     \item {\tt PT}:   Scaled Virtual Potential Temperature (T_v/PKZ)
!     \item {\tt PE}:   Edge pressures
!     \item {\tt Q}:    Tracers
!     \item {\tt PKZ}:  Consistent mean for p^kappa
!   \end{itemize}
!
!  as well as a GRID (to be mentioned later) 
!  and same additional run-specific variables 
!  (dt, iord, jord, nsplit, nspltrac, nspltvrm -- to be mentioned later)
!
! Note: {\tt PT} is not updated if the flag {\tt CONVT} is true.
!
! The internal state is updated each time FVdycore is called.
!
! !REVISION HISTORY:
!
!   WS  05.06.10:  Adapted from FVdycore_GridCompMod
!   WS  05.09.20:  Renamed dyn_comp
!   WS  05.11.10:  Now using dyn_import/export_t containers
!   WS  06.03.01:  Removed tracertrans-related variables
!   WS  06.04.13:  dyn_state moved here from prognostics (temporary?)
!   CC  07.01.29:  Corrected calculation of OMGA
!   AM  07.10.31:  Supports overlap of trac2d and cd_core subcycles
!
!----------------------------------------------------------------------


logical, parameter         :: DEBUG = .true.

! The FV core is always called in its "full physics" mode.  We don't want
! the dycore to know what physics package is responsible for the forcing.
logical, parameter         :: convt = .true.

real(r8), parameter        :: ZERO                 = 0.0_r8
real(r8), parameter        :: HALF                 = 0.5_r8
real(r8), parameter        :: THREE_QUARTERS       = 0.75_r8
real(r8), parameter        :: ONE                  = 1.0_r8
real(r4), parameter        :: ONE_R4               = 1.0_r4
real(r8), parameter        :: SECS_IN_SIX_HOURS    = 21600.0_r8

! Frontogenesis indices
integer, public :: frontgf_idx = -1
integer, public :: frontga_idx = -1

! Index into physics buffer for zonal mean zonal wind.
integer, public :: uzm_idx = -1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dyn_init(file, dyn_state, dyn_in, dyn_out, NLFileName )

! DESCRIPTION: Initialize the FV dynamical core
!
! REVISION HISTORY:
!   05.06.18   Sawyer  Creation
!   06.03.03   Sawyer  Added dyn_state as argument (for reentrancy)
!   06.05.09   Sawyer  Added dyn_conservative to conserve total energy
!   
!==================================================================================

   use constituents,  only : pcnst
   use pmgrid,        only : plon, plat, plev, plevp,                         &
                                    beglonxy, endlonxy, beglatxy, endlatxy,   &
                                    beglat,   endlat,   beglev,   endlev,     &
                                    npr_y, npr_z, nprxy_x, nprxy_y,           &
                                    twod_decomp
   use time_manager,  only : get_step_size
   use pmgrid,        only : dyndecomp_set
   use dynamics_vars, only : dynamics_init
   use dycore,        only : get_resolution
   use dyn_grid,      only : define_cam_grids, initgrid
#if ( defined OFFLINE_DYN )
  use metdata, only:  metdata_dyn_init
#endif
  use spmd_utils, only: npes, masterproc
#if defined(SPMD)
   use parutilitiesmodule, only : gid, parcollective, maxop
   use spmd_dyn, only : geopkdist, geopkblocks, geopk16byte,                 &
                        npes_xy, npes_yz, mpicom_xy, mpicom_yz, mpicom_nyz,  &
                        modcomm_transpose, modcomm_geopk, modcomm_gatscat,   &
                        modc_sw_dynrun, modc_hs_dynrun,                      &
                        modc_send_dynrun, modc_mxreq_dynrun,                 &
                        modc_sw_cdcore, modc_hs_cdcore,                      &
                        modc_send_cdcore, modc_mxreq_cdcore,                 &
                        modc_sw_gather, modc_hs_gather,                      &
                        modc_send_gather, modc_mxreq_gather,                 &
                        modc_sw_scatter, modc_hs_scatter,                    &
                        modc_send_scatter, modc_mxreq_scatter,               &
                        modc_sw_tracer, modc_hs_tracer,                      &
                        modc_send_tracer, modc_mxreq_tracer,                 &
                        modc_onetwo, modc_tracers
   use spmd_dyn, only: spmd_readnl,spmdinit_dyn
   use mpishorthand, only: mpicom
#endif
   use ctem,            only : ctem_init
   use fv_control_mod,  only : dyn_readnl,        & 
                               kord,              & 
                               jord,              & 
                               iord,              &
	                       nsplit,            & 
                               nspltrac,          & 
                               nspltvrm,          & 
                               dyn_conservative,  & 
                               filtcw
   use physconst,       only : omega,             & 
                               rearth,            &
                               rair,              &
                               cpair,             &
                               zvir,              &
                               pi
   use phys_control,    only : use_gw_front
   use qbo,             only : qbo_use_forcing
   use physics_buffer,  only : pbuf_add_field, dtype_r8
   use ppgrid,          only : pcols, pver

   ! ARGUMENTS:
   type(file_desc_t),       intent(in)  :: file       ! PIO file handle for initial or restart file
   type (T_FVDYCORE_STATE), target      :: dyn_state
   type (dyn_import_t),     intent(OUT) :: dyn_in
   type (dyn_export_t),     intent(OUT) :: dyn_out
   character(len=*),        intent(in)  :: NLFileName ! namelist file




! Local variables

  integer, parameter    :: MAXPES = 256
  real(r8), parameter ::  D0_0                  =   0.0_r8
  real(r8), parameter ::  D1E5                  =   1.0e5_r8

  type (T_FVDYCORE_GRID)      , pointer :: GRID      ! For convenience
  type (T_FVDYCORE_CONSTANTS) , pointer :: CONSTANTS ! For convenience
  integer              :: unit

#if defined(SPMD)
  integer :: tmp(npes)
#endif
  integer, allocatable :: jmyz(:),kmyz(:),imxy(:),jmxy(:) ! used for nonblocking receive

  integer :: nstep, nymd, nhms
  integer :: yr, mm, dd, h, m, s, itmp
  integer :: INT_PACK(6)

  integer             :: k
  integer             :: TE_METHOD = 0
  integer             :: NTOTQ                !  Method for total energy remapping
  integer             :: NQ                   !  No. advected tracers
  integer             :: IFIRSTXY             !  No. total tracers
  integer             :: ILASTXY
  integer             :: JFIRSTXY
  integer             :: JLASTXY
  integer             :: JFIRST
  integer             :: JLAST
  integer             :: KFIRST
  integer             :: KLAST
  integer             :: im, jm, km
  real(r8)            :: dt
  real(r8)            :: cp       ! heat capacity of air at constant pressure
  real(r8)            :: ae       ! radius of the earth (m)
  real(r8), allocatable :: ak(:), bk(:)     !  Vertical coordinates
  integer             :: ks               !  True # press. levs

#if !defined( SPMD )
  integer :: npes_xy=1
  integer :: npes_yz=1
  integer :: mpicom=0
  integer :: mpicom_xy=0
  integer :: mpicom_yz=0
  integer :: mpicom_nyz=0
#endif

  !----------------------------------------------------------------------

  if (use_gw_front) then
     call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
          frontgf_idx)
     call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
          frontga_idx)
  end if

  if (qbo_use_forcing) then
     call pbuf_add_field("UZM", "global", dtype_r8, (/pcols,pver/), &
          uzm_idx)
  end if

  ! Initialize hybrid coordinate arrays
  call hycoef_init(file)

  allocate( ak(plev+1) )
  allocate( bk(plev+1) )
  do k = 1, plev+1
    ak(k) = hyai(k) * D1E5
    bk(k) = hybi(k)
    if( bk(k) == D0_0 ) ks = k-1
  end do
!
! Get the layout and store directly in the GRID data structure
!
  GRID => DYN_STATE%GRID     ! For convenience
  CONSTANTS => DYN_STATE%CONSTANTS

  dt = get_step_size()

  call dyn_readnl(nlfilename)

#if defined(SPMD)
  call spmd_readnl(nlfilename)
  call spmdinit_dyn()
#endif

  IFIRSTXY = beglonxy
  ILASTXY  = endlonxy
  JFIRSTXY = beglatxy
  JLASTXY  = endlatxy
  JFIRST   = beglat
  JLAST    = endlat
  KFIRST   = beglev
  KLAST    = endlev
  NTOTQ  = pcnst
  NQ     = pcnst
  IM     = plon
  JM     = plat
  KM     = plev
  cp     = cpair
  ae     = rearth

! Set constants
  constants%pi    = pi
  constants%omega = omega
  constants%ae    = ae
  constants%rair  = rair
  constants%cp    = cp
  constants%cappa = rair/cpair
  constants%zvir  = zvir

  allocate (jmyz(npr_y))
  allocate (kmyz(npr_z))
  allocate (imxy(nprxy_x))
  allocate (jmxy(nprxy_y))

!
! SPMD-related stuff
!
#if defined(SPMD)
    grid%twod_decomp = twod_decomp
    grid%geopkdist= geopkdist
    grid%geopk16byte = geopk16byte
    grid%geopkblocks = geopkblocks
    grid%mod_method = modcomm_transpose
    grid%mod_geopk  = modcomm_geopk
    grid%mod_gatscat  = modcomm_gatscat

  grid%modc_dynrun(1)   = modc_sw_dynrun
  if (modc_hs_dynrun) then
    grid%modc_dynrun(2) = 1
  else
    grid%modc_dynrun(2) = 0
  end if
  if (modc_send_dynrun) then
    grid%modc_dynrun(3) = 1
  else
    grid%modc_dynrun(3) = 0
  end if
  grid%modc_dynrun(4) = modc_mxreq_dynrun

  grid%modc_cdcore(1)   = modc_sw_cdcore
  if (modc_hs_cdcore) then
    grid%modc_cdcore(2) = 1
  else
    grid%modc_cdcore(2) = 0
  end if
  if (modc_send_cdcore) then
    grid%modc_cdcore(3) = 1
  else
    grid%modc_cdcore(3) = 0
  end if
  grid%modc_cdcore(4) = modc_mxreq_cdcore

  grid%modc_gather(1)   = modc_sw_gather
  if (modc_hs_gather) then
    grid%modc_gather(2) = 1
  else
    grid%modc_gather(2) = 0
  end if
  if (modc_send_gather) then
    grid%modc_gather(3) = 1
  else
    grid%modc_gather(3) = 0
  end if
  grid%modc_gather(4) = modc_mxreq_gather

  grid%modc_scatter(1)   = modc_sw_scatter
  if (modc_hs_scatter) then
    grid%modc_scatter(2) = 1
  else
    grid%modc_scatter(2) = 0
  end if
  if (modc_send_scatter) then
    grid%modc_scatter(3) = 1
  else
    grid%modc_scatter(3) = 0
  end if
  grid%modc_scatter(4) = modc_mxreq_scatter

  grid%modc_tracer(1)   = modc_sw_tracer
  if (modc_hs_tracer) then
    grid%modc_tracer(2) = 1
  else
    grid%modc_tracer(2) = 0
  endif
  if (modc_send_tracer) then
    grid%modc_tracer(3) = 1
  else
    grid%modc_tracer(3) = 0
  end if
  grid%modc_tracer(4) = modc_mxreq_tracer

  grid%modc_onetwo       = modc_onetwo
  grid%modc_tracers      = modc_tracers

!
!  Define imxy, jmxy, jmyz, kmyz from ifirstxy, ilastxy, etc.
!
  tmp = 0
  tmp(gid+1) = ilastxy-ifirstxy+1
  call parcollective( mpicom, maxop, npes, tmp )
  imxy(1:nprxy_x) = tmp(1:nprxy_x)

  tmp = 0
  tmp(gid+1) = jlastxy-jfirstxy+1
  call parcollective( mpicom, maxop, npes, tmp )
  do k=1,nprxy_y
    jmxy(k) = tmp((k-1)*nprxy_x+1)
  end do

  tmp = 0
  tmp(gid+1)   = jlast-jfirst+1
  call parcollective( mpicom, maxop, npes, tmp )
  jmyz(1:npr_y) = tmp(1:npr_y)

  tmp = 0
  tmp(gid+1)   = klast-kfirst+1
  call parcollective( mpicom, maxop, npes, tmp )
  do k=1,npr_z
    kmyz(k) = tmp((k-1)*npr_y+1)
  end do

#else
!
! Sensible initializations for OMP-only  (hopefully none of these variables are used...)
!
  grid%twod_decomp = 0
  grid%geopkdist   = .false.
  grid%geopk16byte = .false.
  grid%geopkblocks = 1
  grid%mod_method  = 0
  grid%mod_geopk   = 0
  grid%mod_gatscat = 0

  grid%modc_dynrun(1)  = 0
  grid%modc_dynrun(2)  = 1
  grid%modc_dynrun(3)  = 1
  grid%modc_dynrun(4)  = -1

  grid%modc_cdcore(1)  = 0
  grid%modc_cdcore(2)  = 1
  grid%modc_cdcore(3)  = 1
  grid%modc_cdcore(4)  = -1
    
  grid%modc_gather(1)  = 0
  grid%modc_gather(2)  = 1
  grid%modc_gather(3)  = 1
  grid%modc_gather(4)  = -1

  grid%modc_scatter(1) = 0
  grid%modc_scatter(2) = 1
  grid%modc_scatter(3) = 1
  grid%modc_scatter(4) = -1

  grid%modc_tracer(1)  = 0
  grid%modc_tracer(2)  = 1
  grid%modc_tracer(3)  = 1
  grid%modc_tracer(4)  = -1

  grid%modc_onetwo     = 1
  grid%modc_tracers    = 0

#endif

! These are run-specific variables:  
!     DT              Time step
!     IORD            Order (mode) of X interpolation (1,..,6)
!     JORD            Order (mode) of Y interpolation (1,..,6)
!     NSPLIT          Ratio of big to small timestep (set to zero if in doubt)
!     NSPLTRAC        Ratio of big to tracer timestep
!     NSPLTVRM        Ratio of big to vertical re-mapping timestep
!

  DYN_STATE%DOTIME    = .TRUE.
  DYN_STATE%CHECK_DT  = SECS_IN_SIX_HOURS ! Check max and min every 6 hours.
  DYN_STATE%DT        = DT         ! Should this be part of state??
  DYN_STATE%NSPLIT    = NSPLIT
  DYN_STATE%NSPLTRAC  = NSPLTRAC
  DYN_STATE%NSPLTVRM  = NSPLTVRM
  DYN_STATE%IORD      = IORD
  DYN_STATE%JORD      = JORD
  DYN_STATE%KORD      = KORD
  DYN_STATE%TE_METHOD = TE_METHOD
  DYN_STATE%CONSV     = dyn_conservative
  DYN_STATE%FILTCW    = filtcw
  if (filtcw .gt. 0) then
    if (masterproc) then
      write (iulog,*) ' '
      write (iulog,*) 'Filtering of c-grid winds turned on'
      write (iulog,*) ' '
    end if
  end if

!
! Calculation of orders for the C grid is fixed by D-grid IORD, JORD
  if( iord <= 2 ) then
    DYN_STATE%ICD =  1
  else
    DYN_STATE%ICD = -2
  end if

  if( jord <= 2 ) then
    DYN_STATE%JCD =  1
  else
    DYN_STATE%JCD =  -2
  end if

!
! Calculate NSPLIT if it was specified as 0
  if ( NSPLIT <= 0 ) DYN_STATE%NSPLIT= INIT_NSPLIT(DYN_STATE%DT,IM,JM)

! Calculate NSPLTRAC if it was specified as 0
  if (NSPLTRAC <= 0) then
    if (get_resolution() == '0.23x0.31') then
      DYN_STATE%NSPLTRAC = max ( 1, DYN_STATE%NSPLIT/2 )
    else
      DYN_STATE%NSPLTRAC = max ( 1, DYN_STATE%NSPLIT/4 )
    end if
  end if


! Set NSPLTVRM to 1 if it was specified as 0
  if (NSPLTVRM <= 0) then
    DYN_STATE%NSPLTVRM = 1
  end if
!
!
! Create the dynamics interface
!
  call dyn_create_interface( ifirstxy, ilastxy, jfirstxy, jlastxy, &
                             1, km, ntotq, dyn_in, dyn_out )

!
! Now there is sufficient information to perform the dynamics initialization
! from FVCAM.  Gradually this will be removed, and all the initialization will
! be performed in this module.
!

!
! Initialize the FVDYCORE variables, which are now all in the GRID
!
  call dynamics_init( dt, dyn_state%jord, im, jm, km,         &
                      pi, ae, omega, nq, ntotq, ks,           &
                      ifirstxy, ilastxy, jfirstxy, jlastxy,   &
                      jfirst, jlast, kfirst, klast,           &
                      npes_xy, npes_yz, mpicom, mpicom_xy,    &
                      mpicom_yz, mpicom_nyz,                  &
                      nprxy_x, nprxy_y, npr_y, npr_z,         &
                      imxy,    jmxy,    jmyz,    kmyz,        &
                      ak,      bk,      unit, grid )

!  Clear wall clock time clocks and global budgets

  DYN_STATE%RUN_TIMES = 0
  DYN_STATE%NUM_CALLS = 0

#if ( defined OFFLINE_DYN )
  call metdata_dyn_init(grid)
#endif
  dyndecomp_set=.true.

  deallocate (jmyz)
  deallocate (kmyz)
  deallocate (imxy)
  deallocate (jmxy)
  deallocate (ak)
  deallocate (bk)

  ! Call initgrid (nee initcom, initializes some coordinate info)
  call initgrid()
  ! Define the CAM grids (this has to be after dycore spinup).
  call define_cam_grids()

  ! Setup circulation diagnostics (has to be after define_cam_grids)
  call ctem_init( NLFileName )

  ! Set history defaults (has to be after grids are defined)
  call history_defaults()

  return

contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  dyn_create_interface --- create the dynamics import and export
!
! !INTERFACE:
subroutine dyn_create_interface ( I1, IN, J1, JN, K1, KN, LM, &
                                  dyn_in, dyn_out )
   use infnan, only : inf, assignment(=)
!
! !USES:
  implicit none

! !PARAMETERS:
   integer, intent(in)                 :: I1, IN, J1, JN, K1, KN, LM
   type (dyn_import_t), intent(out)    :: dyn_in
   type (dyn_export_t), intent(out)    :: dyn_out

!EOP
!-----------------------------------------------------------------------

   integer :: l
   integer :: ierror

   allocate( dyn_in%phis( I1:IN, J1:JN  ), stat=ierror  )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PHIS')
   allocate( dyn_in%ps(   I1:IN, J1:JN  ), stat=ierror  )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PS')
   allocate( dyn_in%u3s(  I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array U3S')
   allocate( dyn_in%v3s(  I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array V3S')
   allocate( dyn_in%pe(   I1:IN,K1:KN+1,J1:JN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PE')
   allocate( dyn_in%pt(   I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PT')
   allocate( dyn_in%t3(   I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array T3')
   allocate( dyn_in%pk(   I1:IN,J1:JN,K1:KN+1  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PK')
   allocate( dyn_in%pkz(  I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PKZ')
   allocate( dyn_in%delp( I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array DELP')
!
! allocate tracer contents
!
   allocate( dyn_in%tracer(I1:IN,J1:JN,K1:KN,LM), stat=ierror )
   if ( ierror /= 0 ) then
      write(iulog,*) "Allocation error", ierror, "for tracer"
      call endrun('DYN_COMP ALLOC error: array TRACER')
   endif

   dyn_in%tracer = inf     
   dyn_in%phis = inf
   dyn_in%ps = inf
   dyn_in%u3s = inf
   dyn_in%v3s = inf
   dyn_in%pe = inf
   dyn_in%pt = inf
   dyn_in%t3 = inf
   dyn_in%pk = inf
   dyn_in%pkz = inf
   dyn_in%delp = inf

!
! Output has all of these except phis
!
   dyn_out%phis => dyn_in%phis
   dyn_out%ps   => dyn_in%ps
   dyn_out%u3s  => dyn_in%u3s
   dyn_out%v3s  => dyn_in%v3s
   dyn_out%pe   => dyn_in%pe
   dyn_out%pt   => dyn_in%pt
   dyn_out%t3   => dyn_in%t3
   dyn_out%pk   => dyn_in%pk
   dyn_out%pkz  => dyn_in%pkz
   dyn_out%delp => dyn_in%delp
   dyn_out%tracer => dyn_in%tracer

!
! And several more which are not in the import container
!
   allocate( dyn_out%peln( I1:IN,K1:KN+1,J1:JN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PELN')
   allocate( dyn_out%omga( I1:IN,K1:KN,J1:JN    ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array OMGA')
   allocate( dyn_out%mfx( I1:IN,J1:JN,K1:KN    ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array MFX')
   allocate( dyn_out%mfy( I1:IN,J1:JN,K1:KN    ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array MFY')

   dyn_out%peln = inf
   dyn_out%omga = inf
   dyn_out%mfx  = inf
   dyn_out%mfy  = inf



end subroutine dyn_create_interface
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  init_nsplit --- find proper value for nsplit if not specified
!
! !INTERFACE:
      integer function INIT_NSPLIT(dtime,im,jm) 
!
! !USES:
      implicit none

! !INPUT PARAMETERS:
      real (r8), intent(in) :: dtime      !  time step
      integer, intent(in)   :: im, jm     !  Global horizontal resolution

! !DESCRIPTION:
!
!    If nsplit=0 (module variable) then determine a good value 
!    for ns (used in dynpkg) based on resolution and the large-time-step 
!    (pdt). The user may have to set this manually if instability occurs.
!
! !REVISION HISTORY:
!   00.10.19   Lin     Creation
!   01.03.26   Sawyer  ProTeX documentation
!   01.06.10   Sawyer  Modified for dynamics_init framework
!   03.12.04   Sawyer  Moved here from dynamics_vars.  Now a function
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      real (r8)   pdt                       ! Time-step in seconds
                                            ! Negative dt (backward in time
                                            ! integration) is allowed
      real (r8)   dim
      real (r8)   dim0                      ! base dimension
      real (r8)   dt0                       ! base time step
      real (r8)   ns0                       ! base nsplit for base dimension
      real (r8)   ns                        ! final value to be returned
      real (r8)   one                       ! equal to unity

      parameter ( dim0 = 191._r8  )
      parameter ( dt0  = 1800._r8 )
      parameter ( ns0  = 4._r8    )
      parameter ( one  = 1.0_r8   )

      pdt = int(dtime)   ! dtime is a variable internal to this module
      dim = max ( im, 2*(jm-1) )
      ns  = int ( ns0*abs(pdt)*dim/(dt0*dim0) + THREE_QUARTERS )
      ns  = max ( one, ns )   ! for cases in which dt or dim is too small

      init_nsplit = ns

      return
!EOC
      end function INIT_NSPLIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine history_defaults
        !----------------------------------------------------------------------- 
        ! 
        ! Purpose: 
        !
        ! Build Master Field List of all possible fields in a history file.  Each field has 
        ! associated with it a "long_name" netcdf attribute that describes what the field is, 
        ! and a "units" attribute.
        ! 
        ! Method: Call a subroutine to add each field
        ! 
        ! Author: CCM Core Group
        ! 
        !-----------------------------------------------------------------------

        use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
        use constituents, only: pcnst, cnst_name, cnst_longname, tottnam, cnst_get_ind
        use ppgrid,       only: pver, pverp
        use pmgrid,       only: plev, plevp
        use cam_history,  only: addfld, add_default, horiz_only
        use phys_control, only: phys_getopts

        implicit none

        !-----------------------------------------------------------------------
        !
        ! Local workspace
        !
        integer m                      ! Index
        integer :: ixcldice, ixcldliq  ! constituent indices for cloud liquid and ice water.
        logical :: history_budget      ! output tendencies and state variables for CAM4
                                       ! temperature, water vapor, cloud ice and cloud
                                       ! liquid budgets.
        integer :: history_budget_histfile_num  ! output history file number for budget fields

        !
        ! Call addfld to add each field to the Master Field List.
        !

        !----------------------------------------------------------------------------
        ! Dynamics variables which belong in dynamics specific initialization modules
        !----------------------------------------------------------------------------
        call addfld ('US',    (/ 'lev' /),'A','m/s','Zonal wind, staggered', gridname='fv_u_stagger' )
        call addfld ('VS',    (/ 'lev' /),'A','m/s','Meridional wind, staggered', gridname='fv_v_stagger' )
        call addfld ('US&IC', (/ 'lev' /),'I','m/s','Zonal wind, staggered',gridname='fv_u_stagger' )
        call addfld ('VS&IC', (/ 'lev' /),'I','m/s','Meridional wind, staggered',gridname='fv_v_stagger' )
        call addfld ('PS&IC', horiz_only, 'I','Pa', 'Surface pressure',gridname='fv_centers' )
        call addfld ('T&IC',  (/ 'lev' /),'I','K',  'Temperature',gridname='fv_centers' )
        do m = 1,pcnst
           call addfld (trim(cnst_name(m))//'&IC',(/ 'lev' /),'I','kg/kg',cnst_longname(m),gridname='fv_centers' )
        end do

        do m=1,pcnst
           call addfld (tottnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horz + vert + fixer tendency ',  &
                gridname='fv_centers')
        end do

        call addfld ('DUH',      (/ 'lev' /), 'A','K/s','U horizontal diffusive heating',                 gridname='fv_centers')
        call addfld ('DVH',      (/ 'lev' /), 'A','K/s','V horizontal diffusive heating',                 gridname='fv_centers')
        call addfld ('ENGYCORR', (/ 'lev' /), 'A','W/m2','Energy correction for over-all conservation',   gridname='fv_centers')

        call addfld ('FU',       (/ 'lev' /), 'A','m/s2','Zonal wind forcing term',                       gridname='fv_centers')
        call addfld ('FV',       (/ 'lev' /), 'A','m/s2','Meridional wind forcing term',                  gridname='fv_centers')
        call addfld ('FU_S',     (/ 'lev' /), 'A','m/s2','Zonal wind forcing term on staggered grid',     gridname='fv_u_stagger')
        call addfld ('FV_S',     (/ 'lev' /), 'A','m/s2','Meridional wind forcing term on staggered grid',gridname='fv_v_stagger')
        call addfld ('UTEND',    (/ 'lev' /), 'A','m/s2','U tendency',                                    gridname='fv_centers')
        call addfld ('VTEND',    (/ 'lev' /), 'A','m/s2','V tendency',                                    gridname='fv_centers')
        call addfld ('TTEND',    (/ 'lev' /), 'A','K/s','Total T tendency (all processes)',               gridname='fv_centers')
        call addfld ('LPSTEN',   horiz_only,  'A','Pa/s','Surface pressure tendency',                     gridname='fv_centers')
        call addfld ('VAT',      (/ 'lev' /), 'A','K/s','Vertical advective tendency of T',               gridname='fv_centers')
        call addfld ('KTOOP',    (/ 'lev' /), 'A','K/s','(Kappa*T)*(omega/P)',                            gridname='fv_centers')

        !----------------------------------------------------------------------------
        ! Determine defaults variables
        !----------------------------------------------------------------------------
        call phys_getopts( history_budget_out = history_budget, history_budget_histfile_num_out = history_budget_histfile_num)
        if ( history_budget ) then
           call cnst_get_ind('CLDLIQ', ixcldliq)
           call cnst_get_ind('CLDICE', ixcldice)
           call add_default(tottnam(       1), history_budget_histfile_num, ' ')
           call add_default(tottnam(ixcldliq), history_budget_histfile_num, ' ')
           call add_default(tottnam(ixcldice), history_budget_histfile_num, ' ')
           call add_default('TTEND   '       , history_budget_histfile_num, ' ')
        end if

        call add_default ('US&IC   ', 0, 'I')
        call add_default ('VS&IC   ', 0, 'I')   
        call add_default ('PS&IC      ',0, 'I')
        call add_default ('T&IC       ',0, 'I')

        do m = 1,pcnst
           call add_default(trim(cnst_name(m))//'&IC',0, 'I')
        end do

        !-----------------------------------------------------------------------
        ! End of dynamics variables
        !-----------------------------------------------------------------------

      end subroutine history_defaults

end subroutine dyn_init
!---------------------------------------------------------------------
!BOP
! !ROUTINE:  RUN --- Driver for the NASA finite-volume dynamical core
!
! !INTERFACE:
subroutine dyn_run(ptop, ndt, te0, dyn_state, dyn_in, dyn_out, rc)

! !USES:
   use shr_kind_mod, only  : r8 => shr_kind_r8

   use diag_module, only   : compute_vdot_gradp
   use spmd_utils,   only  : masterproc
   use fv_control_mod, only: ct_overlap, trac_decomp

#if defined( SPMD )
   use mod_comm, only : mp_sendirr,                        &
                        mp_recvirr, mp_sendirr_r4,         &
                        mp_recvirr_r4, mp_send4d_ns,       &
                        mp_recv4d_ns, mp_sendtrirr,        &
                        mp_recvtrirr
#endif
#if ( defined OFFLINE_DYN )
   use metdata,     only: get_met_fields, advance_met, get_us_vs, met_fix_mass, met_rlx
   use pfixer,      only: adjust_press
#endif

   implicit none

#if defined( SPMD )
#include "mpif.h"
#endif

! !PARAMETERS:
   integer, intent(in):: ndt       ! the large time step in seconds
                                   ! Also the mapping time step in this setup
   real(r8), intent(in):: ptop     ! Pressure at model top (interface pres)

   real(r8), intent(out):: te0     ! Total energy before dynamics
   type (T_FVDYCORE_STATE), target :: dyn_state ! Internal state
   type (dyn_import_t)             :: dyn_in    ! Import container
   type (dyn_export_t)             :: dyn_out   ! Export container

   integer, intent(out)               :: rc      ! Return code

! !DESCRIPTION:
!
! Developer: Shian-Jiann Lin, NASA/GSFC; email: lin@dao.gsfc.nasa.gov
!
! Top view of D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!               u(i,j+1)
!                 |
!      v(i,j)---delp(i,j)---v(i+1,j)
!                 |
!               u(i,j)
!
! External routine required: the user needs to supply a subroutine to set up
!                            "Eulerian vertical coordinate" for remapping purpose.
!                             Currently this routine is named as set_eta()
!                             In principle any terrian following vertical
!                             coordinate can be used. The input to fvcore
!                             need not be on the same vertical coordinate
!                             as the output.
!                             If SPMD is defined the Pilgrim communication
!                             library developed by Will Sawyer will be needed.
!
! Remarks: values at poles for both u and v need not be defined; but values for
!          all other scalars needed to be defined at both poles (as polar cap mean
!          quantities). Tracer advection is done "off-line" using the
!          large time step. Consistency is maintained by using the time accumulated
!          Courant numbers and horizontal mass fluxes for the FFSL algorithm.
!          The input "pt" can be either dry potential temperature
!          defined as T/pkz (adiabatic case) or virtual potential temperature
!          defined as T*/pkz (full phys case). IF convt is true, pt is not updated.
!          Instead, virtual temperature is ouput.
!          ipt is updated if convt is false.
!          The user may set the value of nx to optimize the SMP performance
!          The optimal valuse of nx depends on the total number of available
!          shared memory CPUs per node (NS). Assuming the maximm MPI 
!          decomposition is used in the y-direction, set nx=1 if the
!          NS <=4; nx=4 if NS=16.
!
! This version supports overlap of trac2d and cd_core subcycles (Art Mirin, November 2007).
!   This refers to the subcycles described by the "do 2000 n=1,n2" loop and has nothing to
!   do with the "do it=1,nsplit" lower-level subcycling. Each trac2d call (n), other than the last,
!   is overlapped with the subsequent cd_core 'series' (n+1). The controlling namelist variable
!   is ct_overlap. The overlapping trac2d calls are carried out on the second set of
!   npes_yz processes (npes_yz <= iam < 2*npes_yz). The tracer arrays are sent to the
!   auxiliary processes prior to the do-2000 loop. During each subcycle (other than the last),
!   the dp0 array is sent prior to the cd_core series; arrays cx, cy, mfx, mfy are sent directly
!   from cd_core during the last call in the series (it=nsplit). At the completion of the last
!   auxiliary trac2d subcycle (n=n2-1), the updated tracer values are returned to the
!   primary processes; the last tracer subcycle (n=n2) is carried out on the primary processes.
!   Communication calls are nonblocking, with attempt to overlap computation to the extent
!   possible. The CCSM mpi layer (wrap_mpi) is used. Tags with values greater than npes_xy
!   are chosen to avoid possible interference between the messages sent from cd_core and
!   the geopk-related transpose messages called from cd_core thereafter. The auxiliary
!   processes must use values of jfirst, jlast, kfirst, klast corresponding to their primary
!   process antecedents, whereas by design those values are (1,0,1,0), resp. (set in spmdinit_dyn).
!   We therefore add auxiliary subdomain limits to the grid datatype: jfirstct, jlastct,
!   kfirstct, klastct. For the primary processes, these are identical to the actual subdomain
!   limits; for the secondary processes, these correspond to the subdomain limits of the
!   antecedent primary process. These values are communicated to the auxiliary processes
!   during initialization (spmd_vars_init). During the auxiliary calculations (and allocations)
!   we temporarily set jfirst equal to jfirstct (etc.) and when done, restore to the original
!   values. Other information needed by the auxiliary processes is obtained through the grid
!   datatype.
!
! This version supports tracer decomposition with trac2d (Art Mirin, January 2008).
!   This option is mutually exclusive with ct_overlap. Variable "trac_decomp" is the size of the
!   decomposition. The tracers are divided into trac_decomp groups, and the kth group is solved
!   on the kth set of npes_yz processes. Much of the methodology is similar to that for ct_overlap.
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Initial SMP version delivered to Will Sawyer
!   WS  99.10.03:  1D MPI completed and tested; 
!   WS  99.10.11:  Additional documentation
!   WS  99.10.28:  benergy and te_map added; arrays pruned
!   SJL 00.01.01:  SMP and MPI enhancements; documentation
!   WS  00.07.13:  Changed PILGRIM API
!   WS  00.08.28:  SPMD instead of MPI_ON
!   AAM 00.08.10:  Add kfirst:klast
!   WS  00.12.19:  phis now distr., LLNL2DModule initialized here
!   WS  01.02.02:  bug fix: parsplit only called for FIRST time
!   WS  01.04.09:  Added initialization of ghost regions
!   WS  01.06.10:  Removed if(first) section; use module
!   AAM 01.06.27:  Extract te_map call into separate routine
!   AAM 01.07.13:  Get rid of dynpkg2; recombine te_map;
!                  perform forward transposes for 2D decomposition
!   WS  01.12.10:  Ghosted PT (changes benergy, cd_core, te_map, hswf)
!   WS  03.08.05:  removed vars dcaf, rayf, ideal, call to hswf
!                  (idealized physics is now in physics package)
!   WS  03.08.13:  Removed ghost region from UXY
!   WS  05.06.11:  Inserted into FVCAM_GridCompMod
!   WS  06.03.03:  Added dyn_state as argument (for reentrancy)
!   WS  06.06.28:  Using new version of benergy
!
!EOP
!-----------------------------------------------------------------------
!BOC

   integer, parameter  ::  DYN_RUN_SUCCESS           = 0
   integer, parameter  ::  DYN_RUN_FAILURE           = -1
   integer, parameter  ::  DYN_RUN_R4_NOT_SUPPORTED  = -10
   integer, parameter  ::  DYN_RUN_MUST_BE_2D_DECOMP = -20
   real(r8), parameter ::  D1_0                      = 1.0_r8

! Variables from the dynamics interface (import or export)

   real(r8), pointer :: phisxy(:,:)   ! surface geopotential (grav*zs)
   real(r8), pointer :: psxy(:,:)     ! Surface pressure (pa) 
   real(r8), pointer :: t3xy(:,:,:)   ! temperature (K)
   real(r8), pointer :: ptxy(:,:,:)   ! scaled (virtual) potential temperature
   real(r8), pointer :: delpxy(:,:,:) ! Pressure thickness
   real(r8), pointer :: tracer(:,:,:,:) ! Tracers
   real(r8), pointer :: uxy(:,:,:)    ! u wind velocities, staggered grid
   real(r8), pointer :: vxy(:,:,:)    ! v wind velocities, staggered grid

!--------------------------------------------------------------------------------------
! The arrays pexy, pkxy, pkzxy must be pre-computed as input to benergy(). 
! They are NOT needed if dyn_state%consv=.F.; updated on output (to be used 
! by physdrv) Please refer to routine pkez on the algorithm for computing pkz
! from pe and pk
!--------------------------------------------------------------------------------------

   real(r8), pointer :: pexy(:,:,:)   ! Pres at layer edges 
   real(r8), pointer :: pkxy(:,:,:)   ! pe**cappa
   real(r8), pointer :: pkzxy(:,:,:)  ! finite-volume mean of pk

! Export state only variables
   real(r8), pointer :: pelnxy(:,:,:) ! Natural logarithm of pe
   real(r8), pointer :: omgaxy(:,:,:) ! vertical pressure velocity (pa/sec)
   real(r8), pointer :: mfxxy(:,:,:)  ! mass flux in X (Pa m^\2 / s)
   real(r8), pointer :: mfyxy(:,:,:)  ! mass flux in Y (Pa m^\2 / s)

! Other pointers (for convenience)
   type (T_FVDYCORE_GRID)      , pointer :: GRID      ! For convenience
   type (T_FVDYCORE_CONSTANTS) , pointer :: CONSTANTS ! For convenience

!  YZ variables currently allocated on stack... should they be on the heap?

   real(r8) :: ps(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast)
   real(r8) :: phis(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast)
   real(r8) :: pe(dyn_state%grid%im,  &
                  dyn_state%grid%kfirst:dyn_state%grid%klast+1,&
                  dyn_state%grid%jfirst:dyn_state%grid%jlast)
   real(r8) :: delp(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast,&
                    dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: pk(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast,&
                  dyn_state%grid%kfirst:dyn_state%grid%klast+1)
   real(r8) :: pkz(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast, &
                   dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: u(dyn_state%grid%im,   &
                 dyn_state%grid%jfirst-dyn_state%grid%ng_d:dyn_state%grid%jlast+dyn_state%grid%ng_s,&
                 dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: v(dyn_state%grid%im,   &
                 dyn_state%grid%jfirst-dyn_state%grid%ng_s:dyn_state%grid%jlast+dyn_state%grid%ng_d,&
                 dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: pt(dyn_state%grid%im,  &
                  dyn_state%grid%jfirst-dyn_state%grid%ng_d:dyn_state%grid%jlast+dyn_state%grid%ng_d,&
                  dyn_state%grid%kfirst:dyn_state%grid%klast) 

   real(r8) :: pi
   real(r8) :: om       ! angular velocity of earth's rotation  
   real(r8) :: cp       ! heat capacity of air at constant pressure
   real(r8) :: ae       ! radius of the earth (m)

   real(r8) :: rair     ! Gas constant of the air
   real(r8) :: cappa    ! R/Cp
   real(r8) :: zvir     ! Virtual effect constant ( = rwv/rair-1 )

   logical :: consv     ! Energy conserved?

   integer :: im        ! dimension in east-west
   integer :: jm        ! dimension in North-South
   integer :: km        ! number of Lagrangian layers
   integer :: jfirst    ! starting latitude index for MPI
   integer :: jlast     ! ending latitude index for MPI
   integer :: kfirst    ! starting vertical index for MPI
   integer :: klast     ! ending vertical index for MPI
   integer :: klastp    ! klast, except km+1 when klast=km
   integer :: ntotq     ! total # of tracers to be advected
   integer :: iord      ! parameter controlling monotonicity in E-W
                                   ! recommendation: iord=4
   integer :: jord      ! parameter controlling monotonicity in N-S 
                                   ! recommendation: jord=4
   integer :: kord      ! parameter controlling monotonicity in mapping
                                   ! recommendation: kord=4
   integer :: te_method ! parameter controlling total energy mapping
                                   ! recommendation: te_method=0 (PPM)
                                   ! GEOS5 uses te_method=1 (Cubic Interp.)
   integer :: icd       ! X algorithm order on C-grid
   integer :: jcd       ! Y algorithm order on C-grid
   integer :: ng_c      ! Ghosting width on C-grid
   integer :: ng_d      ! Ghosting width on D-grid
   integer :: ng_s      ! Ghosting width (staggered, for winds)
   integer :: ns        ! overall split

   integer :: ifirstxy, ilastxy, jfirstxy, jlastxy  ! xy decomposition
   integer :: npr_z

   logical :: cd_penul

   real(r8), allocatable, target    :: q_internal(:,:,:,:)    ! Pointers to tracers
   integer i, j, k, iq          ! Loop indicies
   real(r8) umax                ! Maximum winds, m/s
   parameter (umax = 300.0_r8)

   integer    nx          ! # of split pieces in x-direction; for performance, the
#if defined( UNICOSMP )
   parameter (nx = 1)
#else
   parameter (nx = 4)     ! user may set nx=1 if there is NO shared memory multitasking
#endif
   integer ipe, it, iv
   integer nsplit, nspltrac, n, n2, nv
   integer incount, outcount
   integer iqa, iqb, iqc, iqd, mq  ! used for tracer transpose grouping
#if (! defined SPMD)
   integer :: mpicom = 0
#endif

! Geometric arrays

! Move the following 3D arrays to an initialization routine?
   real(r8), allocatable :: worka(:,:,:),workb(:,:,:),dp0(:,:,:),cx(:,:,:),cy(:,:,:)
   real(r8), allocatable :: mfx(:,:,:), mfy(:,:,:)
   real(r8), allocatable :: delpf(:,:,:), uc(:,:,:), vc(:,:,:)
   real(r8), allocatable :: dwz(:,:,:), pkc(:,:,:), wz(:,:,:)
   real(r8), allocatable :: dpt(:,:,:), peln(:,:,:)
   real(r8), allocatable :: pkcc(:,:,:), wzc(:,:,:)
! The following variables are work arrays for xy=>yz transpose
   real(r8), allocatable :: pkkp(:,:,:), wzkp(:,:,:)
! The following variables are xy instantiations
   real(r8), allocatable :: tempxy(:,:,:), dp0xy(:,:,:), wzxy(:,:,:)
! psxy3 is dummy 3d variant of psxy
   real(r8), allocatable :: psxy3(:,:,:)
! phisxy3 is dummy 3d variant of phisxy
   real(r8), allocatable :: phisxy3(:,:,:)
   real(r8), pointer     :: q3xypt(:,:,:)
   real(r8), pointer     :: q3yzpt(:,:,:)
   real(r8)              :: tte(dyn_state%grid%jm)
   real(r8)              :: XXX(dyn_state%grid%km)


#if ( defined OFFLINE_DYN )
   real(r8), allocatable :: ps_obs(:,:)
   real(r8), allocatable :: ps_mod(:,:)
   real(r8), allocatable :: u_tmp(:,:,:)
   real(r8), allocatable :: v_tmp(:,:,:)
#endif

   double precision zamda, zam5
   logical fill

   integer imh
   real(r8) dt
   real(r8) bdt
   integer filtcw
   integer modc_tracers, mlast

! cd_core / trac2d overlap and tracer decomposition data (AAM)
   integer :: commnyz                                         ! n*npes_yz communicator
   integer :: jfirstct, jlastct, kfirstct, klastct            ! primary subdomain limits
   integer :: jkstore(4)                                      ! storage for subdomain limits
   integer :: iamlocal                                        ! task number (global indexing)
   integer :: iremotea(trac_decomp)                           ! source/target; id array
   integer :: iremote                                         ! source/target; working id
   integer :: ndp0, ncx, ncy, nmfx, nmfy, ntrac               ! message sizes
   integer :: dp0tag, cxtag, cytag, mfxtag, mfytag, tractag   ! message tags
   integer :: cxtaga(trac_decomp), cytaga(trac_decomp)        ! tag arrays for cd_core
   integer :: mfxtaga(trac_decomp), mfytaga(trac_decomp)      ! tag arrays for cd_core
   logical :: ct_aux                                          ! true if auxiliary process
   logical :: s_trac                                          ! true for cd_core posting tracer-related sends
   integer, allocatable :: ctreq(:,:)                         ! used for nonblocking receive
   integer, allocatable :: ctstat(:,:,:)                      ! used for nonblocking receive
   integer, allocatable :: ctreqs(:,:)                        ! used for nonblocking send
   integer, allocatable :: ctstats(:,:,:)                     ! used for nonblocking send
   integer, allocatable :: cdcreqs(:,:)                       ! used for nonblocking send in cd_core
   integer, pointer :: ktloa(:)                               ! lower limit of tracer decomposition (global)
   integer, pointer :: kthia(:)                               ! upper limit of tracer decomposition (global)
   integer ktlo                                               ! lower limit of tracer decomposition (local)
   integer kthi                                               ! upper limit of tracer decomposition (local)
   integer kt, tagu, naux, kaux, ntg0

   logical :: print_subcycling = .true.
   logical :: c_dotrac, t_dotrac
   logical :: convt_local

   data fill  /.true./              ! perform a simple filling algorithm
                                    ! in case negatives were found

! C.-C. Chen, omega calculation
  real(r8) ::   &
    cx_om(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast, &
          dyn_state%grid%kfirst:dyn_state%grid%klast)! Courant no. in X
  real(r8) ::  &
    cy_om(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast+1, &
          dyn_state%grid%kfirst:dyn_state%grid%klast)        ! Courant no. in Y
  real(r8) :: &
    pexy_om(dyn_state%grid%ifirstxy:dyn_state%grid%ilastxy,dyn_state%grid%km+1, &
            dyn_state%grid%jfirstxy:dyn_state%grid%jlastxy)

   rc       =  DYN_RUN_FAILURE      ! Set initially to fail

   phisxy   => dyn_in%phis
   psxy     => dyn_in%ps
   uxy      => dyn_in%u3s
   vxy      => dyn_in%v3s
   t3xy     => dyn_in%t3
   ptxy     => dyn_in%pt
   delpxy   => dyn_in%delp
   tracer   => dyn_in%tracer
   pexy     => dyn_in%pe
   pkxy     => dyn_in%pk
   pkzxy    => dyn_in%pkz

   pelnxy   => dyn_out%peln
   omgaxy   => dyn_out%omga
   mfxxy    => dyn_out%mfx
   mfyxy    => dyn_out%mfy

   grid => dyn_state%grid    ! For convenience
   constants => DYN_STATE%CONSTANTS

   ns   = dyn_state%nsplit   ! large split (will be subdivided later)
   n2   = dyn_state%nspltrac ! tracer split(will be subdivided later)
   nv   = dyn_state%nspltvrm ! vertical re-mapping split
   icd  = dyn_state%icd
   jcd  = dyn_state%jcd
   iord = dyn_state%iord
   jord = dyn_state%jord
   kord = dyn_state%kord
   filtcw = dyn_state%filtcw

   consv     = dyn_state%consv
   te_method = dyn_state%te_method

   pi   =  constants%pi
   om   =  constants%omega
   ae   =  constants%ae
   rair =  constants%rair
   cp   =  constants%cp
   cappa=  constants%cappa
   zvir =  constants%zvir

   im = grid%im
   jm = grid%jm
   km = grid%km

   ng_c  = grid%ng_c
   ng_d  = grid%ng_d
   ng_s  = grid%ng_s

   ifirstxy = grid%ifirstxy
   ilastxy  = grid%ilastxy
   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy

   jfirst   = grid%jfirst
   jlast    = grid%jlast
   kfirst   = grid%kfirst
   klast    = grid%klast
   klastp   = grid%klastp

   ntotq    = grid%ntotq
   modc_tracers = grid%modc_tracers

   npr_z    = grid%npr_z

! cd_core/trac2d overlap and tracer decomposition
   jfirstct = grid%jfirstct
   jlastct  = grid%jlastct
   kfirstct = grid%kfirstct
   klastct  = grid%klastct
   commnyz = grid%commnyz
   iamlocal = grid%iam
! kaux is an index describing the set of npes_yz processes; 0 for first set, 1 for second set, etc.
   kaux = iamlocal/grid%npes_yz
! ct_aux is true if current process is auxiliary, false otherwise
   ct_aux = ((ct_overlap .gt. 0 .and. kaux .eq. 1) .or.      &
      (trac_decomp .gt. 1 .and. kaux .ge. 1 .and. kaux .lt. trac_decomp))
! define message tags to exceed npes_xy so as not to interfere with geopotential transpose tags
! tags below correspond to communicated variables with ct_overlap and trac_decomp
   dp0tag = grid%npes_xy + 5
   cxtag = dp0tag + 1
   cytag = dp0tag + 2
   mfxtag = dp0tag + 3
   mfytag = dp0tag + 4
   tractag = dp0tag + 5
! ntg0 is upper bound on number of needed tags beyond tracer tags for ct_overlap and trac_decomp
   ntg0 = 10

#if ( defined OFFLINE_DYN )
!
! advance the meteorology data
!
    call advance_met(grid)
!
! set the staggered winds (verticity winds) to offline meteorological data
!
    call get_us_vs( grid, u, v )
#endif

   if ( km > 1 ) then         ! not shallow water equations

      if( consv ) then
       if (grid%iam .lt. grid%npes_xy) then
! Compute globally integrated Total Energy (te0)

         call t_startf ('benergy')
!
! Tests indicate that t3 does not have consistent
! pole values, e.g. t3(:,1,k) are not all the same.
! Not clear why this is not the case: it may be that the pole
! values are not consistent on the restart file.  For the time being, 
! perform a parallel sum over t3 and correct the pole values
!
         if ( jfirstxy == 1 ) then
           call par_xsum(grid, t3xy(:,1,:), km, XXX )
           do k=1, km
             do i=ifirstxy, ilastxy
               t3xy(i,1,k)             = XXX(k) / real(im,r8)
             enddo
           enddo
         endif
         if ( jlastxy == jm ) then
           call par_xsum(grid, t3xy(:,jm,:), km, XXX )
           do k=1, km
             do i=ifirstxy, ilastxy
               t3xy(i,jm,k)             = XXX(k) / real(im,r8)
             enddo
           enddo
         endif
         call benergy(grid, uxy, vxy, t3xy, delpxy,         &
                      tracer(:,:,:,1), pexy, pelnxy, phisxy,           &
                      zvir, cp,  rair, tte, te0 )
         call t_stopf  ('benergy')
       endif  ! (grid%iam .lt. grid%npes_xy)
      endif

   endif


! Allocate temporary work arrays
! Change later to use pointers for SMP performance???
! (prime candidates: uc, vc, delpf)

   call t_startf ('dyn_run_alloc')

      if (ct_aux) then
! Temporarily set subdomain limits in auxiliary process equal to those of antecedent
!   to allow following arrays to have proper size
! (Normally, sizes of unneeded arrays for auxiliary processes will be deliberately small.)
         jkstore(1) = jfirst
         jkstore(2) = jlast
         jkstore(3) = kfirst
         jkstore(4) = klast
         jfirst = jfirstct
         jlast  = jlastct
         kfirst = kfirstct
         klast  = klastct
      endif
      allocate( worka(im,jfirst:     jlast,     kfirst:klast) )
      allocate( workb(im,jfirst:     jlast,     kfirst:klast) )
      allocate(   dp0(im,jfirst-1:   jlast,     kfirst:klast) )
      allocate(   mfx(im,jfirst:     jlast,     kfirst:klast) )
      allocate(   mfy(im,jfirst:     jlast+1,   kfirst:klast) )
      allocate(    cx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(    cy(im,jfirst:     jlast+1,   kfirst:klast) )
      dp0(:,:,:) = 0._r8
      mfx(:,:,:) = 0._r8
      mfy(:,:,:) = 0._r8
      cx(:,:,:) = 0._r8
      cy(:,:,:) = 0._r8
      if (ct_aux) then
! Restore subdomain limits in auxiliary process
         jfirst = jkstore(1)
         jlast  = jkstore(2)
         kfirst = jkstore(3)
         klast  = jkstore(4)
      endif
      allocate( delpf(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(    uc(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(    vc(im,jfirst-2:   jlast+2,   kfirst:klast) )
      allocate(   dpt(im,jfirst-1:   jlast+1,   kfirst:klast) )
      allocate(   dwz(im,jfirst-1:    jlast,    kfirst:klast+1) )
      allocate(   pkc(im,jfirst-1:   jlast+1,   kfirst:klast+1) ) 
      allocate(    wz(im,jfirst-1:   jlast+1,   kfirst:klast+1) )
      allocate(  pkcc(im,jfirst  :   jlast  ,   kfirst:klast+1) ) 
      allocate(   wzc(im,jfirst  :   jlast  ,   kfirst:klast+1) ) 
      allocate(  peln(im,kfirst  :   klast+1,   jfirst:jlast) )    ! For consv = .true.
      allocate(pkkp(im,jfirst:jlast,kfirst:klast+1))
      allocate(wzkp(im,jfirst:jlast,kfirst:klast+1))
      allocate(wzxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1))
      allocate(tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
      allocate(dp0xy(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
      allocate(psxy3(ifirstxy:ilastxy,jfirstxy:jlastxy,npr_z))
      allocate(phisxy3(ifirstxy:ilastxy,jfirstxy:jlastxy,npr_z))

#if ( defined OFFLINE_DYN )
      allocate( ps_obs(im,jfirst:jlast) )
      allocate( ps_mod(im,jfirst:jlast) )
      allocate( u_tmp(im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
      allocate( v_tmp(im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
#endif

!
! Allocation of tracers
!
   if (ct_aux) then
! Temporarily set subdomain limits in auxiliary process equal to those of antecedent
!   to allow trac2d temporary storage to have proper size
      jfirst = jfirstct
      jlast  = jlastct
      kfirst = kfirstct
      klast  = klastct
   endif
   allocate ( q_internal(im, jfirst:jlast, kfirst:klast, ntotq) )
! Trac2d-related mpi quantities for ct_overlap and tracer decomposition
   allocate (ctreq(ntotq+ntg0,trac_decomp))
   allocate (ctreqs(ntotq+ntg0,trac_decomp))
   allocate (cdcreqs(trac_decomp,4))
   cdcreqs(:,:) = 0
#if defined(SPMD)
   allocate (ctstat(MPI_STATUS_SIZE,ntotq+ntg0,trac_decomp))
   allocate (ctstats(MPI_STATUS_SIZE,ntotq+ntg0,trac_decomp))
#endif
! Compute i.d.'s of remote processes for ct_overlap or trac_decomp
   naux = 0
   if ((ct_overlap .gt. 0 .and. kaux .lt. 2) .or.      &
      (trac_decomp .gt. 1 .and. kaux .lt. trac_decomp)) then
! Identify involved processes
      iremotea(:) = -1
      naux = max(1,trac_decomp-1)
      if (kaux .eq. 0) then
! Primary process - identify corresponding auxiliary process(es)
          do kt = 1, naux
             iremotea(kt) = iamlocal + kt*grid%npes_yz
             cxtaga(kt) = cxtag + (kt-1)*(ntotq+ntg0)
             cytaga(kt) = cytag + (kt-1)*(ntotq+ntg0)
             mfxtaga(kt) = mfxtag + (kt-1)*(ntotq+ntg0)
             mfytaga(kt) = mfytag + (kt-1)*(ntotq+ntg0)
          enddo
      else
! Auxiliary process - identify corresponding primary process
          iremotea(1) = iamlocal - kaux*grid%npes_yz
      endif
      iremote = iremotea(1)
! Message sizes
      ndp0  = im*(jlast-jfirst+2       )*(klast-kfirst+1)
      ncx   = im*(jlast-jfirst+2*ng_d+1)*(klast-kfirst+1)
      ncy   = im*(jlast-jfirst+2       )*(klast-kfirst+1)
      nmfx  = im*(jlast-jfirst+1       )*(klast-kfirst+1)
      nmfy  = im*(jlast-jfirst+2       )*(klast-kfirst+1)
      ntrac = im*(jlast-jfirst+1       )*(klast-kfirst+1)
   endif
   if (ct_aux) then
! Restore subdomain limits in auxiliary process
      jfirst = jkstore(1)
      jlast  = jkstore(2)
      kfirst = jkstore(3)
      klast  = jkstore(4)
   endif

! Set tracer limits to be supplied to trac2d (needed even without tracer decomposition)
   ktloa => grid%ktloa
   kthia => grid%kthia
   ktlo = grid%ktlo
   kthi = grid%kthi

   call t_stopf  ('dyn_run_alloc')

! Determine splitting
   bdt = ndt

! Second/third level splitting (nsplit and n2 variables overloaded)
   n2     = (n2+nv   -1) / nv
   nsplit = (ns+n2*nv-1) / (n2*nv)
   dt     = bdt / real(nsplit*n2*nv,r8)

   if (print_subcycling) then
      print_subcycling = .false.
      if (masterproc) then
         write(iulog,*) 'FV subcycling - nv, n2, nsplit, dt = ', nv, n2, nsplit, dt
         if( (nsplit*n2*nv /= dyn_state%nsplit) .or. (n2*nv /= dyn_state%nspltrac) ) then
            write(iulog,*) "ERROR:  Because of loop nesting, FV dycore can't use the specified namelist settings for subcycling"
            write(iulog,*) '  The original namelist settings were:'
            write(iulog,*) '  NSPLIT   = ', dyn_state%nsplit
            write(iulog,*) '  NSPLTRAC = ', dyn_state%nspltrac
            if( dyn_state%nspltvrm /= 1 ) write(iulog,*) '  NSPLTVRM = ', dyn_state%nspltvrm
            write(iulog,*)
            write(iulog,*) '  NSPLIT needs to be a multiple of NSPLTRAC'
            if( dyn_state%nspltvrm /= 1 ) write(iulog,*) '    which in turn needs to be a multiple of NSPLTVRM.'
            write(iulog,*) '  Suggested settings would be:'
            write(iulog,*) '  NSPLIT   = ', nsplit*n2*nv
            write(iulog,*) '  NSPLTRAC = ', n2*nv
            if( dyn_state%nspltvrm /= 1 ) write(iulog,*) '  NSPLTVRM = ', nv
            call endrun("Bad namelist settings for FV subcycling.")
         endif
      endif
   endif

!
! IF convt_local is false, pt is updated for the next iteration of the 3000-loop
! On the last iteration, convt_local is set to convt
!
   convt_local = .false.
!
! Begin vertical re-mapping sub-cycle loop
!
  do 3000 iv = 1, nv

     if(iv == nv) convt_local = convt
!
! Transpose XY arrays to YZ
!
     call t_barrierf('sync_xy_to_yz_1', grid%commdyn)
     call t_startf ('xy_to_yz')

     if (grid%iam .lt. grid%npes_xy) then

        if (grid%twod_decomp .eq. 1) then

#if defined( SPMD )

! Embed psxy and phisxy in 3D array since transpose machinery cannot handle 2D arrays

!$omp parallel do private(i,j,k)
           do k=1,npr_z
              do j=jfirstxy,jlastxy
                 do i=ifirstxy,ilastxy
                    psxy3(i,j,k) = psxy(i,j)
                    phisxy3(i,j,k) = phisxy(i,j)
                 enddo
              enddo
           enddo

           if (grid%modc_onetwo .eq. 1) then
              call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                              grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                              modc=grid%modc_dynrun )
              call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                              grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                              modc=grid%modc_dynrun )

              call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                              grid%xy2d_to_yz2d%RecvDesc, phisxy3, phis,                 &
                              modc=grid%modc_dynrun )
              call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                              grid%xy2d_to_yz2d%RecvDesc, phisxy3, phis,                 &
                              modc=grid%modc_dynrun )
           else
              call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                              grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                              phisxy3, phis,                                             &
                              modc=grid%modc_dynrun )
              call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                              grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                              phisxy3, phis,                                             &
                              modc=grid%modc_dynrun )
           endif

!
! if OFFLINE_DYN is defined, u and v are filled at this point
!
#if defined( OFFLINE_DYN )
           call mp_sendirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                            grid%uxy_to_u%RecvDesc, uxy, u_tmp,                         &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                            grid%uxy_to_u%RecvDesc, uxy, u_tmp,                         &
                            modc=grid%modc_dynrun )

           call mp_sendirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                            grid%vxy_to_v%RecvDesc, vxy, v_tmp,                         &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                            grid%vxy_to_v%RecvDesc, vxy, v_tmp,                         &
                            modc=grid%modc_dynrun )
!$omp parallel do private(i,j,k)
           do k = kfirst,klast
              do j = jfirst,jlast
                 do i = 1,im
                    u(i,j,k)    = (1._r8-met_rlx(k))*u_tmp(i,j,k) + met_rlx(k)*u(i,j,k)
                    v(i,j,k)    = (1._r8-met_rlx(k))*v_tmp(i,j,k) + met_rlx(k)*v(i,j,k)
                 enddo
              enddo
           enddo
#else
           call mp_sendirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                            grid%uxy_to_u%RecvDesc, uxy, u,                             &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                            grid%uxy_to_u%RecvDesc, uxy, u,                             &
                            modc=grid%modc_dynrun )

           call mp_sendirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                            grid%vxy_to_v%RecvDesc, vxy, v,                             &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                            grid%vxy_to_v%RecvDesc, vxy, v,                             &
                            modc=grid%modc_dynrun )
#endif

           call mp_sendirr( grid%commxy, grid%pexy_to_pe%SendDesc,                      &
                            grid%pexy_to_pe%RecvDesc, pexy, pe,                         &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%pexy_to_pe%SendDesc,                      &
                            grid%pexy_to_pe%RecvDesc, pexy, pe,                         &
                            modc=grid%modc_dynrun )

           call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                            grid%ijk_xy_to_yz%RecvDesc, delpxy, delp,                   &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                            grid%ijk_xy_to_yz%RecvDesc, delpxy, delp,                   &
                            modc=grid%modc_dynrun )

           call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,                     &
                            grid%pkxy_to_pkc%RecvDesc, pkxy, pk,                        &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,                     &
                            grid%pkxy_to_pkc%RecvDesc, pkxy, pk,                        &
                            modc=grid%modc_dynrun )

           call mp_sendirr( grid%commxy, grid%ptxy_to_pt%SendDesc,                      &
                            grid%ptxy_to_pt%RecvDesc, ptxy, pt,                         &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%ptxy_to_pt%SendDesc,                      &
                            grid%ptxy_to_pt%RecvDesc, ptxy, pt,                         &
                            modc=grid%modc_dynrun )

           if (modc_tracers .eq. 0) then
              do mq = 1, ntotq
                 q3xypt => tracer(:,:,:,mq)
                 q3yzpt => q_internal(:,:,:,mq)
                 call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                  &
                                  grid%ijk_xy_to_yz%RecvDesc, q3xypt, q3yzpt,               &
                                  modc=grid%modc_dynrun )
                 call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                  &
                                  grid%ijk_xy_to_yz%RecvDesc, q3xypt, q3yzpt,               &
                                  modc=grid%modc_dynrun )
              enddo
           else
              do mq = 1, ntotq, modc_tracers
                 mlast = min(mq+modc_tracers-1,ntotq)
                 call mp_sendtrirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                      &
                                    grid%ijk_xy_to_yz%RecvDesc, tracer, q_internal, mq, mlast, ntotq,    &
                                    grid%ifirstxy, grid%ilastxy, grid%jfirstxy, grid%jlastxy,     &
                                    1, grid%km,                                                   &
                                    1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, grid%klast, &
                                    modc=grid%modc_tracer )
                 call mp_recvtrirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                      &
                                    grid%ijk_xy_to_yz%RecvDesc, tracer, q_internal, mq, mlast, ntotq,    &
                                    grid%ifirstxy, grid%ilastxy, grid%jfirstxy, grid%jlastxy,     &
                                    1, grid%km,                                                   &
                                    1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, grid%klast, &
                                    modc=grid%modc_tracer )
              enddo
           endif

#else
           write(iulog,*)'DYN_COMP:dyn_run -- SPMD must be defined for 2D decomp -- returning'
           rc = DYN_RUN_MUST_BE_2D_DECOMP
           return    ! Not possible to have 2D decomposition with SPMD undefined
#endif
        else ! twod_decomp
           do j=jfirst,jlast
              do i=1,im
                 ps(i,j) = psxy(i,j)
                 phis(i,j) = phisxy(i,j)
              enddo
           enddo

!$omp parallel do private(i,j,k)
           do j = jfirst,jlast
              do k = 1,km+1
                 do i = 1,im
                    pe(i,k,j) = pexy(i,k,j)
                 enddo
              enddo
           enddo

!$omp parallel do private(i,j,k)
           do k = 1,km+1
              do j = jfirst,jlast
                 do i = 1,im
                    pk(i,j,k) = pkxy(i,j,k)
                 enddo
              enddo
           enddo

!$omp parallel do private(i,j,k)
           do k = 1,km
              do j = jfirst,jlast
                 do i = 1,im
#if defined( OFFLINE_DYN )
                    u(i,j,k)    = (1._r8-met_rlx(k))*uxy(i,j,k) + met_rlx(k)*u(i,j,k)
                    v(i,j,k)    = (1._r8-met_rlx(k))*vxy(i,j,k) + met_rlx(k)*v(i,j,k)
#else 
                    u(i,j,k)    = uxy(i,j,k)
                    v(i,j,k)    = vxy(i,j,k)
#endif
                    delp(i,j,k) = delpxy(i,j,k)
                    pt(i,j,k)   = ptxy(i,j,k)
                 enddo
              enddo
           enddo

           do mq = 1, ntotq
!
!  For now just copy in the contents of tracer; later, use pointers
!
! TODO:  q_internal(mq) => tracer(mq)    ! Make sure not to allocate q_internal in this case
!
              q_internal(1:im,jfirst:jlast,kfirst:klast,mq) = &
                  tracer(1:im,jfirst:jlast,kfirst:klast,mq)
           enddo

        endif  !  (grid%twod_decomp .eq. 1)

     endif  !  (grid%iam .lt. grid%npes_xy)

#if defined(SPMD)
! Send tracers to auxiliary processes when overlapping
     if (ct_overlap .gt. 0 .and. n2 .gt. 1 .and. kaux .eq. 0) then
        do iq = 1, ntotq
           call mpiisend(q_internal(:,:,:,iq), ntrac, mpir8, iremote, tractag+iq-1, commnyz, ctreqs(5+iq,1))
        enddo
     endif
! Send tracers to auxiliary processes when decomposing
     if (trac_decomp .gt. 1 .and. kaux .eq. 0) then
        do kt = 2, trac_decomp
           do iq = ktloa(kt), kthia(kt)
              tagu = tractag+iq-1 + (kt-2)*(ntotq+ntg0)
              call mpiisend(q_internal(:,:,:,iq), ntrac, mpir8, iremotea(kt-1), tagu, commnyz, ctreqs(5+iq,kt-1))
           enddo
        enddo
     endif
#endif

     call t_stopf  ('xy_to_yz')

! C.-C. Chen
     omgaxy(:,:,:) = ZERO
!
! Begin tracer sub-cycle loop
!
     do 2000 n=1, n2

        if( ntotq > 0 ) then

           call t_barrierf('sync_small_ts_init', grid%commdyn)
           call t_startf ('small_ts_init')

!$omp parallel do private(i, j, k)
           do k=kfirst,klast
              do j=jfirst,jlast
                 do i=1,im
! Save initial delp field before the small-time-step
! Initialize the CFL number accumulators: (cx, cy)
! Initialize total mass fluxes: (mfx, mfy)
                    dp0(i,j,k) = delp(i,j,k)
                    cx(i,j,k) = ZERO
                    cy(i,j,k) = ZERO
                    mfx(i,j,k) = ZERO
                    mfy(i,j,k) = ZERO
                 enddo
              enddo
           enddo

#if defined( SPMD )
           if (grid%iam .lt. grid%npes_yz) then
              call mp_send4d_ns( grid%commyz, im, jm, km,                     &
                                 1, jfirst, jlast, kfirst, klast, 1, 0, dp0 )
              call mp_recv4d_ns( grid%commyz, im, jm, km,                     &
                                 1, jfirst, jlast, kfirst, klast, 1, 0, dp0 )
           endif  !  (grid%iam .lt. grid%npes_yz)
#endif

           call t_stopf  ('small_ts_init')

        endif

#if defined(SPMD)
! Send dp0 to auxiliary processes when overlapping or tracer decomposition
        if (kaux .eq. 0) then
           if (ct_overlap .gt. 0 .and. n .lt. n2) then
              call mpiisend(dp0, ndp0, mpir8, iremote, dp0tag, commnyz, ctreqs(1,1))
           endif
           if (trac_decomp .gt. 1) then
              do kt = 2, trac_decomp
                 tagu = dp0tag + (kt-2)*(ntotq+ntg0)
                 call mpiisend(dp0, ndp0, mpir8, iremotea(kt-1), tagu, commnyz, ctreqs(1,kt-1))
              enddo
           endif
        endif
#endif
!
! Begin dynamics sub-cycle loop
!
        do it=1, nsplit

           if(it == nsplit .and. n == n2) then
              ipe = 1                     ! end of fvcore; output pe for te_map
           elseif(it == 1 .and. n == 1) then
              ipe = -1                    ! start of cd_core
           else
              ipe = 0
           endif

! determine whether this is the second to last call to cd_core or not 
           cd_penul = .false.
           if ( nsplit > 1 ) then
              if ( (n == n2) .and. (it == nsplit-1) ) cd_penul = .true.
           elseif ( n2 > 1 ) then
              if ( n == n2-1 ) cd_penul = .true.
           endif

           if (cd_penul) then
              if (ipe == -1) then
                 ipe = -2   ! second to last is also the first
              else
                 ipe = 2
              endif
           endif

! s_trac is true if cd_core is to post sends for ct_overlap or trac_decomp
! such sends are posted during last inner cd_core subcycle
           s_trac = ((ct_overlap .gt. 0 .and. it .eq. nsplit .and. n .lt. n2) .or.     &
                    (trac_decomp .gt. 1 .and. it .eq. nsplit))

! C.-C. Chen
           if((it == nsplit).and.(n == n2).and.(iv == nv)) then
!$omp parallel do private(j)
              do j=jfirstxy,jlastxy
                 pexy_om(ifirstxy:ilastxy,1:km+1,j) = pexy(ifirstxy:ilastxy,1:km+1,j)
              end do
           endif


! Call the Lagrangian dynamical core using small tme step

           call t_barrierf('sync_cd_core', grid%commdyn)
           call t_startf ('cd_core')

           if (grid%iam .lt. grid%npes_xy) then
              call cd_core(grid,   nx,     u,   v,   pt,                     &
                            delp,   pe,     pk,  nsplit,  dt,                &
                            ptop,   umax,   pi, ae,  cp,  cappa,             &
                            icd,    jcd, iord, jord,   ipe,                  &
                            om,     phis,     cx  ,  cy, mfx, mfy,           &
                            delpf, uc, vc, pkz, dpt, worka,                  &
                            dwz, pkc, wz,  phisxy, ptxy, pkxy,               &
                            pexy, pkcc, wzc, wzxy, delpxy,                   &
                            pkkp, wzkp, cx_om, cy_om, filtcw, s_trac,        &
                            naux, ncx, ncy, nmfx, nmfy, iremotea,            &
                            cxtaga, cytaga, mfxtaga, mfytaga, cdcreqs(1,1),  &
                            cdcreqs(1,2), cdcreqs(1,3), cdcreqs(1,4))
              ctreqs(2,:) = cdcreqs(:,1)
              ctreqs(3,:) = cdcreqs(:,2)
              ctreqs(4,:) = cdcreqs(:,3)
              ctreqs(5,:) = cdcreqs(:,4)
           endif  !  (grid%iam .lt. grid%npes_yz)

           call t_stopf  ('cd_core')

! C.-C. Chen
           if((it == nsplit).and.(n == n2).and.(iv == nv)) then

!$omp  parallel do     &
!$omp  default(shared) &
!$omp  private(i,j,k)
              do j=jfirstxy,jlastxy
                 do k=1,km
                    do i=ifirstxy,ilastxy
                       omgaxy(i,k,j) = omgaxy(i,k,j)+HALF*(pexy(i,k,j)+pexy(i,k+1,j)- &
                                       pexy_om(i,k,j)-pexy_om(i,k+1,j))/dt
                    end do
                 end do
                 do k=1,km+1
                    do i=ifirstxy,ilastxy
                       pexy_om(i,k,j) = HALF*(pexy_om(i,k,j)+pexy(i,k,j))
                    end do
                 end do
              end do

!-----------------------------------------------------
! Add the v*grad(p) term to omega (dp/dt) for physics
!-----------------------------------------------------
!
!pw           call t_barrief('sync_vdot_gradp', grid%commdyn)
              call t_startf ('vdot_gradp')
              if (grid%iam .lt. grid%npes_xy) then
                 call compute_vdot_gradp( grid, dt, dt/dt, cx_om, cy_om, pexy_om, omgaxy )
              endif  !  (grid%iam .lt. grid%npes_xy)
              call t_stopf  ('vdot_gradp')

           endif

        enddo  !  it = 1, nsplit - dynamics sub-cycle loop

        if( ntotq .ne. 0 ) then
#if ( defined OFFLINE_DYN )
           if (met_fix_mass) then
              ps_mod(:,:) = ps(:,:)
              ! get the observed PS interpolated to current substep
              call get_met_fields( grid, ps_obs, n2, n )

              ! adjust mass fluxes and edge pressures to be consistent with observed PS 
              call adjust_press( grid, ps_mod, ps_obs, mfx, mfy, pexy )

              ! make pkxy consistent with the adjusted pexy
!$omp parallel do private(i,j,k)
              do i=ifirstxy,ilastxy
                 do j=jfirstxy,jlastxy
                    do k=1,km+1
                       pkxy(i,j,k) = pexy(i,k,j)**cappa
                    enddo
                 enddo
              enddo

              ! adjust courant numbers to be consistent with the adjusted mass fluxes
!$omp parallel do private(i,j,k)
              do i=1,im
                 do j=jfirst,jlast
                    do k=kfirst,klast
                       if (i .ne. 1) cx(i,j,k) = mfx(i,j,k)/(HALF*(dp0(i-1,j,k)+dp0(i,j,k)))
                       if (i .eq. 1) cx(i,j,k) = mfx(i,j,k)/(HALF*(dp0(1,j,k)+dp0(im,j,k)))
                    enddo
                 enddo
              enddo
!$omp parallel do private(i,j,k)
              do i=1,im
                 do j=jfirst,jlast
                    do k=kfirst,klast
                       if ((j .gt. 1) .and. (j .lt. jm)) cy(i,j,k) = &
                             mfy(i,j,k)/(HALF*(dp0(i,j-1,k)+dp0(i,j,k)))/grid%cose(j)
                    enddo
                 enddo
              enddo
           endif
#endif

! WS 2006-12-04 : this seems like the safest place to preprocess and
!                 transpose the C-grid mass-flux and later the
!                 Courant numbers for potential output
!        

! Horizontal mass fluxes

           if (grid%iam .lt. grid%npes_xy) then

              if (grid%twod_decomp .eq. 1) then
#if defined( SPMD )
!$omp parallel do private(i,j,k)
                 do k = kfirst,klast
                    do j = jfirst,jlast
                       do i = 1,im
                          worka(i,j,k) = mfx(i,j,k)*(ae*grid%dp)*(grid%dl*ae*grid%cosp(j))/(ndt) ! Pa m^2/s
                          workb(i,j,k) = mfy(i,j,k)*(grid%dl*ae*grid%cosp(j))*(ae*grid%dp)/(ndt*grid%cose(j)) ! Pa m^2 / s
                       enddo
                    enddo
                 enddo
                 if (grid%modc_onetwo .eq. 1) then
                    call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                     grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                     modc=grid%modc_dynrun )
                    call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                     grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                     modc=grid%modc_dynrun )

                    call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                     grid%ijk_yz_to_xy%RecvDesc, workb, mfyxy,                 &
                                     modc=grid%modc_dynrun )
                    call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                     grid%ijk_yz_to_xy%RecvDesc, workb, mfyxy,                 &
                                     modc=grid%modc_dynrun )
                 else
                    call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                     grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                     workb, mfyxy,                                             &
                                     modc=grid%modc_dynrun )
                    call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                     grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                     workb, mfyxy,                                             &
                                     modc=grid%modc_dynrun )
                 endif

#else
                 write(iulog,*)'DYN_COMP:dyn_run -- SPMD must be defined for 2D decomp -- returning'
                 rc = DYN_RUN_MUST_BE_2D_DECOMP
                 return    ! Not possible to have 2D decomposition with SPMD undefined
#endif
              else   ! if not twod_decomp   (1D or sequential)
!$omp parallel do private(i,j,k)
                 do k = kfirst,klast
                    do j = jfirst,jlast
                       do i = 1,im
                          mfxxy(i,j,k) = mfy(i,j,k)*(grid%dl*ae*grid%cosp(j))*(ae*grid%dp)/(ndt*grid%cose(j)) ! Pa m^2 / s
                          mfyxy(i,j,k) = mfy(i,j,k)*(grid%dl*ae*grid%cosp(j))*(ae*grid%dp)/(ndt*grid%cose(j)) ! Pa m^2 / s
                       enddo
                    enddo
                 enddo
              endif

           endif  !  (grid%iam .lt. grid%npes_xy)


! Perform large-tme-step scalar transport using the accumulated CFL and
! mass fluxes
 
           call t_barrierf('sync_trac2d', grid%commdyn)
           call t_startf ('trac2d')

! Overlap trac2d with subsequent cd_core set, or decompose over tracers

           if ((ct_overlap .gt. 0 .and. n .lt. n2 .and. kaux .lt. 2) .or.   &
              (trac_decomp .gt. 1 .and. kaux .lt. trac_decomp)) then

              if (kaux .eq. 0) then

! Primary process

! Send data to auxiliary yz decomposition
! Communicate tracers on first subcycle only
! Also post receive of new tracer values from aux processes

#if defined(SPMD)
                 if (n .eq. 1) then
! Block on send of tracers to aux
                    if (ct_overlap .gt. 0) then
                       do iq = 1, ntotq
                          call mpiwait(ctreqs(5+iq,1), ctstats(1,5+iq,1))
                       enddo
                    endif
                    if (trac_decomp .gt. 1) then
                       do kt = 2, trac_decomp
                          do iq = ktloa(kt), kthia(kt)
                             call mpiwait(ctreqs(5+iq,kt-1), ctstats(1,5+iq,kt-1))
                          enddo
                       enddo
                    endif
! Post receive for updated tracers from aux
                    if (ct_overlap .gt. 0) then
                       do iq = 1, ntotq
                          call mpiirecv(q_internal(:,:,:,iq), ntrac, mpir8, iremote,    &
                               tractag+iq-1, commnyz, ctreq(iq,1))
                       enddo
                    endif
                    if (trac_decomp .gt. 1) then
                       do kt = 2, trac_decomp
                          do iq = ktloa(kt), kthia(kt)
                             tagu = tractag+iq-1 + (kt-2)*(ntotq+ntg0)
                             call mpiirecv(q_internal(:,:,:,iq), ntrac, mpir8, iremotea(kt-1),   &
                                  tagu, commnyz, ctreq(iq,kt-1))
                          enddo
                       enddo
                    endif
                 endif  !  (n .eq. 1)

                 if (ct_overlap .gt. 0) then
! Block on send of dp0 to aux
                    call mpiwait(ctreqs(1,1), ctstats(1,1,1))
! Block on sends from cd_core to aux
                    call mpiwait(ctreqs(2,1), ctstats(1,2,1))
                    call mpiwait(ctreqs(3,1), ctstats(1,3,1))
                    call mpiwait(ctreqs(4,1), ctstats(1,4,1))
                    call mpiwait(ctreqs(5,1), ctstats(1,5,1))
                 endif
                 if (trac_decomp .gt. 1) then
                    do kt = 2, trac_decomp
! Block on send of dp0 to aux
                       call mpiwait(ctreqs(1,kt-1), ctstats(1,1,kt-1))
! Block on sends from cd_core to aux
                       call mpiwait(ctreqs(2,kt-1), ctstats(1,2,kt-1))
                       call mpiwait(ctreqs(3,kt-1), ctstats(1,3,kt-1))
                       call mpiwait(ctreqs(4,kt-1), ctstats(1,4,kt-1))
                       call mpiwait(ctreqs(5,kt-1), ctstats(1,5,kt-1))
                    enddo
                 endif
#endif

              else

! Auxiliary process

! Temporarily set subdomain limits and process index in auxiliary process equal to those of antecedent
                 jfirst = jfirstct
                 jlast  = jlastct
                 kfirst = kfirstct
                 klast  = klastct
                 grid%jfirst = jfirstct
                 grid%jlast  = jlastct
                 grid%kfirst = kfirstct
                 grid%klast  = klastct
! Translate process index to frame of auxiliary yz decomposition for use with auxiliary
!    communication in trac2d
                 grid%iam = iremote

! Receive data from primary yz decomposition
! Include tracers first subcycle only

#if defined(SPMD)
                 if (n .eq. 1) then
                    do iq = ktlo, kthi
                       tagu = tractag+iq-1 + (kaux-1)*(ntotq+ntg0)
                       call mpiirecv(q_internal(:,:,:,iq), ntrac, mpir8, iremote, tagu, commnyz, ctreq(5+iq,1))
                       call mpiwait(ctreq(5+iq,1), ctstat(1,5+iq,1))
                    enddo
                 endif
                 tagu = dp0tag + (kaux-1)*(ntotq+ntg0)
                 call mpiirecv(dp0, ndp0, mpir8, iremote, tagu, commnyz, ctreq(1,1))
                 tagu = cxtag + (kaux-1)*(ntotq+ntg0)
                 call mpiirecv(cx, ncx, mpir8, iremote, tagu, commnyz, ctreq(2,1))
                 tagu = cytag + (kaux-1)*(ntotq+ntg0)
                 call mpiirecv(cy, ncy, mpir8, iremote, tagu, commnyz, ctreq(3,1))
                 tagu = mfxtag + (kaux-1)*(ntotq+ntg0)
                 call mpiirecv(mfx, nmfx, mpir8, iremote, tagu, commnyz, ctreq(4,1))
                 tagu = mfytag + (kaux-1)*(ntotq+ntg0)
                 call mpiirecv(mfy, nmfy, mpir8, iremote, tagu, commnyz, ctreq(5,1))
                 call mpiwait(ctreq(1,1), ctstat(1,1,1))
                 call mpiwait(ctreq(2,1), ctstat(1,2,1))
                 call mpiwait(ctreq(3,1), ctstat(1,3,1))
                 call mpiwait(ctreq(4,1), ctstat(1,4,1))
                 call mpiwait(ctreq(5,1), ctstat(1,5,1))
#endif

              endif  !  (kaux .eq. 0)

           else

! Block on receive of updated tracers from aux (last subcycle)
#if defined(SPMD)
              if (ct_overlap .gt. 0 .and. n .eq. n2 .and. n2 .gt. 1 .and. kaux .eq. 0) then
                 do iq = 1, ntotq
                    call mpiwait(ctreq(iq,1), ctstat(1,iq,1))
                 enddo
              endif  !  (ct_overlap .gt. 0 .and. n .eq. n2 .and. n2 .gt. 1 .and. kaux .eq. 0)
#endif

           endif  !  (ct_overlap .gt. 0 .and. n .lt. n2 .and. kaux .lt. 2)
                  ! or (trac_decomp .gt. 1 .and. kaux .lt. trac_decomp)

! Call tracer advection

           c_dotrac = ct_overlap .gt. 0 .and.    &                       
                ((n .lt. n2 .and. kaux .eq. 1) .or. (n .eq. n2 .and. kaux .eq. 0))
           t_dotrac = ct_overlap .eq. 0 .and. kaux .lt. trac_decomp

           if (c_dotrac .or. t_dotrac) then
              call trac2d( grid, dp0(:,jfirst:jlast,:),    q_internal,         &  
                           cx,    cy,     mfx,    mfy,    iord,   jord,        &
                           fill,  ktlo,   kthi,   workb,  worka  )
           endif

! Return data to primary yz decomposition
! For overlap, next-to-last subcycle only; for tracer decomp, last subcycle only
#if defined(SPMD)
           if (ct_aux .and. ((ct_overlap .gt. 0 .and. n .eq. n2-1) .or.   &
              (trac_decomp .gt. 1 .and. n .eq. n2)))  then
              do iq = ktlo, kthi
                 tagu = tractag+iq-1 + (kaux-1)*(ntotq+ntg0)
                 call mpiisend(q_internal(:,:,:,iq), ntrac, mpir8, iremote, tagu, commnyz, ctreqs(5+iq,1))
                 call mpiwait(ctreqs(5+iq,1), ctstats(1,5+iq,1))
              enddo
           endif
#endif

! For tracer decomposition, block on receive of updated tracers from aux (last subcycle)
#if defined(SPMD)
           if (trac_decomp .gt. 1 .and. n .eq. n2 .and. kaux .eq. 0) then
              do kt = 2, trac_decomp
                 do iq = ktloa(kt), kthia(kt)
                    call mpiwait(ctreq(iq,kt-1), ctstat(1,iq,kt-1))
                 enddo
              enddo
           endif  !  (trac_decomp .gt. 1 .and. n .eq. n2)
#endif

! Restore subdomain limits and process index in auxiliary process
           if (ct_aux) then
              jfirst = jkstore(1)
              jlast  = jkstore(2)
              kfirst = jkstore(3)
              klast  = jkstore(4)
              grid%jfirst = jkstore(1)
              grid%jlast  = jkstore(2)
              grid%kfirst = jkstore(3)
              grid%klast  = jkstore(4)
              grid%iam = iamlocal
           endif

! NOTE: for cd_core / trac2d overlap, tracer data is returned to primary processes
!   prior to n=n2 call to trac2d

           call t_stopf  ('trac2d')

#if ( defined OFFLINE_DYN )
           if (met_fix_mass) then

              if (grid%twod_decomp .eq. 1) then

#if defined( SPMD )
                 call mp_sendirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                                  grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                                  modc=grid%modc_dynrun )
                 call mp_recvirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                                  grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                                  modc=grid%modc_dynrun )
#endif

              else
!$omp parallel do private(i,j,k)
                 do j = jfirst,jlast
                    do k = kfirst,klast+1
                       do i = 1,im
                          pe(i,k,j) = pexy(i,k,j)
                       enddo
                    enddo
                 enddo

              endif

              do j = jfirst,jlast

                 if (klast .eq. km) ps_mod(:,j) = pe(:,km+1,j)

              enddo

           endif
#endif

        endif

2000 continue !   do 2000 n=1, n2 - tracer sub-cycle loop

     call t_barrierf('sync_yz_to_xy_1', grid%commdyn)

     if (grid%iam .lt. grid%npes_xy) then

        if (grid%twod_decomp .eq. 1) then
!
! Transpose ps, u, v, and tracer from yz to xy decomposition
!
! Note: pt, pe and pk will have already been transposed through
! call to geopk in cd_core. geopk does not actually require
! secondary xy decomposition; direct 16-byte technique works just
! as well, perhaps better. However, transpose method is used on last
! call to avoid having to compute these three transposes now.
!
#if defined (SPMD)

           call t_startf ('yz_to_xy_psuv')

! Transpose ps
! Embed in 3D array since transpose machinery cannot handle 2D arrays

           call mp_sendirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                            grid%yz2d_to_xy2d%RecvDesc, ps, psxy3,                     &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                            grid%yz2d_to_xy2d%RecvDesc, ps, psxy3,                     &
                            modc=grid%modc_dynrun )

!$omp parallel do private(i,j)
           do j = jfirstxy,jlastxy
              do i = ifirstxy,ilastxy
                 psxy(i,j) = psxy3(i,j,1)
              enddo
           enddo

! Transpose u
           call mp_sendirr( grid%commxy, grid%u_to_uxy%SendDesc,                       &
                            grid%u_to_uxy%RecvDesc, u, uxy,                            &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%u_to_uxy%SendDesc,                       &
                            grid%u_to_uxy%RecvDesc, u, uxy,                            &
                            modc=grid%modc_dynrun )

! Transpose v
           call mp_sendirr( grid%commxy, grid%v_to_vxy%SendDesc,                       &
                            grid%v_to_vxy%RecvDesc, v, vxy,                            &
                            modc=grid%modc_dynrun )
           call mp_recvirr( grid%commxy, grid%v_to_vxy%SendDesc,                       &
                            grid%v_to_vxy%RecvDesc, v, vxy,                            &
                            modc=grid%modc_dynrun )

           call t_stopf  ('yz_to_xy_psuv')

           call t_startf ('yz_to_xy_q')

           if (modc_tracers .eq. 0) then
              do mq = 1, ntotq

!
! Transpose
!
                 q3yzpt => q_internal(:,:,:,mq)
                 q3xypt => dyn_out%tracer(:,:,:,mq)
                 call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                                  grid%ijk_yz_to_xy%RecvDesc, q3yzpt, q3xypt,                 &
                                  modc=grid%modc_dynrun )
                 call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                                  grid%ijk_yz_to_xy%RecvDesc, q3yzpt, q3xypt,                 &
                                  modc=grid%modc_dynrun )
              enddo
           else
              do mq = 1, ntotq, modc_tracers
                 mlast = min(mq+modc_tracers-1,ntotq)
                 call mp_sendtrirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                     &
                                    grid%ijk_yz_to_xy%RecvDesc, q_internal, dyn_out%tracer,      &
                                    mq, mlast, ntotq, 1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, &
                                    grid%klast, grid%ifirstxy, grid%ilastxy, grid%jfirstxy,      &
                                    grid%jlastxy, 1, grid%km,                                    &
                                    modc=grid%modc_tracer )
                 call mp_recvtrirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                     &
                                    grid%ijk_yz_to_xy%RecvDesc, q_internal, dyn_out%tracer,      &
                                    mq, mlast, ntotq, 1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, &
                                    grid%klast, grid%ifirstxy, grid%ilastxy, grid%jfirstxy,      &
                                    grid%jlastxy, 1, grid%km,                                    &
                                    modc=grid%modc_tracer )
              enddo
           endif

           call t_stopf  ('yz_to_xy_q')

#endif

        else

           call t_startf ('yz_to_xy_psuv')

           do j = jfirst,jlast
              do i = 1,im
                 psxy(i,j) = ps(i,j)
              enddo
           enddo

!$omp parallel do private(i,j,k)
           do k = kfirst,klast
              do j = jfirst,jlast
                 do i = 1,im
                    uxy(i,j,k) = u(i,j,k)
                    vxy(i,j,k) = v(i,j,k)
                 enddo
              enddo
           enddo

           call t_stopf  ('yz_to_xy_psuv')

           call t_startf ('yz_to_xy_q')

!
! TODO:  does the array have to be copied?  Copying pointers sufficient?


!$omp parallel do private(i,j,k,mq)
           do mq = 1,ntotq
!
! Temporary -- here the pointers will ultimately be set, not the contents copied
!
              do k = 1,km
                 do j = jfirst,jlast
                    do i = 1,im
                       dyn_out%tracer(i,j,k,mq) = q_internal(i,j,k,mq)
                    enddo
                 enddo
              enddo
           enddo

           call t_stopf  ('yz_to_xy_q')

        endif  !  (grid%twod_decomp .eq. 1)

     endif  !  (grid%iam .lt. grid%npes_xy)

     if ( km > 1 ) then           ! not shallow water equations

! Perform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.

! 
! te_map requires uxy, vxy, psxy, pexy, pkxy, phisxy, q3xy, and ptxy
!
        call t_barrierf('sync_te_map', grid%commdyn)
        call t_startf ('te_map')

        if (grid%iam .lt. grid%npes_xy) then
           call te_map(grid,     consv,   convt_local, psxy, omgaxy,       &
                       pexy,     delpxy,  pkzxy,  pkxy,   ndt,             &
                       nx,       uxy,     vxy,    ptxy,   dyn_out%tracer,  & 
                       phisxy,   cp,      cappa,  kord,   pelnxy,          &
                       te0,      tempxy,  dp0xy,  mfxxy,  mfyxy,           &
                       te_method )

           if( .not. convt_local ) then
!$omp parallel do private(i,j,k)
              do j=jfirstxy,jlastxy
                 do k=1,km
                    do i=ifirstxy,ilastxy
                       t3xy(i,j,k) = ptxy(i,j,k)*pkzxy(i,j,k)/ &
                             (D1_0+zvir*dyn_out%tracer(i,j,k,1))
                    end do
                 end do
              end do
           end if

        endif  !  (grid%iam .lt. grid%npes_xy)

        call t_stopf ('te_map')

     endif
!
! te_map computes uxy, vxy, psxy, delpxy, pexy, pkxy, pkzxy,
! pelnxy, omgaxy, tracer, ptxy, mfxxy and mfyxy
!
3000 continue  !    do 3000 iv = 1, nv - vertical re-mapping sub-cycle loop

  call t_startf ('dyn_run_dealloc')

  deallocate( worka )
  deallocate( workb )
  deallocate( dp0 )
  deallocate( mfx )
  deallocate( mfy )
  deallocate(  cx )
  deallocate(  cy )
  deallocate( delpf )
  deallocate( uc    )
  deallocate( vc    )
  deallocate( dpt   )
  deallocate( dwz   )
  deallocate( pkc   )
  deallocate(  wz   )
  deallocate( pkcc )
  deallocate( wzc )
  deallocate( peln )
  deallocate( pkkp )
  deallocate( wzkp )
  deallocate( wzxy )
  deallocate( tempxy )
  deallocate( dp0xy )
  deallocate( psxy3 )
  deallocate( phisxy3 )
  deallocate( q_internal )
  deallocate (ctreq)
  deallocate (ctreqs)
  deallocate (cdcreqs)
#if defined(SPMD)
  deallocate (ctstat)
  deallocate (ctstats)
#endif
#if ( defined OFFLINE_DYN )
  deallocate( ps_obs )
  deallocate( ps_mod )
  deallocate( u_tmp )
  deallocate( v_tmp )
#endif

  call t_stopf  ('dyn_run_dealloc')

  rc = DYN_RUN_SUCCESS

!----------------------------------------------------------
! WS 03.08.05: removed idealized physics (Held-Suarez) 
! from here (is now in Physics package).
!----------------------------------------------------------

!EOC
end subroutine dyn_run
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine dyn_final(restart_file, dyn_state, dyn_in, dyn_out)

use dynamics_vars, only : dynamics_clean

  character(LEN=*)             , intent(IN   ) :: restart_file
  type (T_FVDYCORE_STATE), target              :: dyn_state
  type (dyn_import_t), intent(inout)           :: dyn_in
  type (dyn_export_t), intent(inout)           :: dyn_out

! BEGIN

   call dynamics_clean    ( dyn_state%grid  )
   call dyn_free_interface( dyn_in, dyn_out )

  contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  dyn_free_interface --- free the dynamics import and export
!
! !INTERFACE:
subroutine dyn_free_interface ( dyn_in, dyn_out )

! !USES:
  implicit none

! !PARAMETERS:
   type (dyn_import_t), intent(inout) :: dyn_in
   type (dyn_export_t), intent(inout) :: dyn_out
!EOP
!-----------------------------------------------------------------------
   integer :: l

   if ( associated(dyn_in%phis) ) deallocate( dyn_in%phis )
   if ( associated(dyn_in%ps) )   deallocate( dyn_in%ps )
   if ( associated(dyn_in%u3s) )  deallocate( dyn_in%u3s )
   if ( associated(dyn_in%v3s) )  deallocate( dyn_in%v3s )
   if ( associated(dyn_in%pe) )   deallocate( dyn_in%pe )
   if ( associated(dyn_in%pt) )   deallocate( dyn_in%pt )
   if ( associated(dyn_in%t3) )   deallocate( dyn_in%t3 )
   if ( associated(dyn_in%pk) )   deallocate( dyn_in%pk )
   if ( associated(dyn_in%pkz) )  deallocate( dyn_in%pkz )
   if ( associated(dyn_in%delp) ) deallocate( dyn_in%delp )
   if ( associated(dyn_in%tracer) ) deallocate( dyn_in%tracer)

   if ( associated(dyn_out%ps) )   nullify( dyn_out%ps )
   if ( associated(dyn_out%u3s) )  nullify( dyn_out%u3s )
   if ( associated(dyn_out%v3s) )  nullify( dyn_out%v3s )
   if ( associated(dyn_out%pe) )   nullify( dyn_out%pe )
   if ( associated(dyn_out%pt) )   nullify( dyn_out%pt )
   if ( associated(dyn_out%t3) )   nullify( dyn_out%t3 )
   if ( associated(dyn_out%pk) )   nullify( dyn_out%pk )
   if ( associated(dyn_out%pkz) )  nullify( dyn_out%pkz )
   if ( associated(dyn_out%delp) ) nullify( dyn_out%delp )
   if ( associated(dyn_out%tracer) )   nullify( dyn_out%tracer )

   if ( associated(dyn_out%omga) ) deallocate( dyn_out%omga )
   if ( associated(dyn_out%peln) ) deallocate( dyn_out%peln )
   if ( associated(dyn_out%mfx) )  deallocate( dyn_out%mfx )
   if ( associated(dyn_out%mfy) )  deallocate( dyn_out%mfy )

end subroutine dyn_free_interface
!-----------------------------------------------------------------------

end subroutine dyn_final
!-----------------------------------------------------------------------


end module dyn_comp
