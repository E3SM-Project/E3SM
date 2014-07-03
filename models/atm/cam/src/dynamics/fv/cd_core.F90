
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: cd_core --- Dynamical core for both C- and D-grid Lagrangian
!                       dynamics
!
! !INTERFACE:
 subroutine cd_core(grid,   nx,     u,   v,   pt,                  &
                    delp,   pe,     pk,  ns,  dt,                  &
                    ptopin, umax,   pi, ae,  cp,  akap,            &
                    iord_c, jord_c, iord_d, jord_d,   ipe,         &
                    om,     hs,     cx3  ,  cy3, mfx, mfy,         &
                    delpf, uc, vc, ptc, dpt, ptk,                  &
                    wz3, pkc, wz,  hsxy, ptxy, pkxy,               &
                    pexy, pkcc, wzc, wzxy, delpxy,                 &
                    pkkp, wzkp, cx_om, cy_om, filtcw, s_trac,      &
                    mlt, ncx, ncy, nmfx, nmfy, iremote,            &
                    cxtag, cytag, mfxtag, mfytag,                  &
                    cxreqs, cyreqs, mfxreqs, mfyreqs)

! !USES:
   use shr_kind_mod,  only : r8 => shr_kind_r8
   use sw_core,       only : d2a2c_winds, c_sw, d_sw
   use pft_module,    only : pft2d
   use dynamics_vars, only : T_FVDYCORE_GRID
   use FVperf_module, only : FVstartclock, FVstopclock, FVbarrierclock
   use cam_logfile,   only : iulog
   use fv_control_mod, only: div24del2flag, del2coef
   use spmd_utils,     only: masterproc
   use abortutils,     only: endrun

#if defined( SPMD )
   use mod_comm,      only : mp_send4d_ns, mp_recv4d_ns,     &
                             mp_send2_ns, mp_recv2_ns,                   &
                             mp_send3d_2, mp_recv3d_2,                   &
                             mp_send3d, mp_recv3d, mp_sendirr,           &
                             mp_recvirr
   use mpishorthand
#endif

#if defined( OFFLINE_DYN )
   use metdata,       only : get_met_fields, met_winds_on_walls
#endif
   use metdata,       only : met_rlx

   implicit none

! !INPUT PARAMETERS:

  type (T_FVDYCORE_GRID), intent(inout) :: grid! grid (for YZ decomp)
  integer, intent(in) :: nx                 ! # of split pieces in longitude direction
  integer, intent(in) :: ipe                ! ipe=1:  end of cd_core()
                                            ! ipe=-1,-2: start of cd_core()
                                            ! ipe=-2,2: second to last call to cd_core()
                                            ! ipe=0 :
  integer, intent(in) :: ns                 ! Number of internal time steps (splitting)
  integer, intent(in) :: iord_c, jord_c     ! scheme order on C grid in X and Y dir.
  integer, intent(in) :: iord_d, jord_d     ! scheme order on D grid in X and Y dir.
  integer, intent(in) :: filtcw             ! flag for filtering C-grid winds

! ct_overlap data
  logical, intent(in) :: s_trac                        ! true to post send for ct_overlap or
                                                       ! tracer decomposition information
  integer, intent(in) :: mlt                           ! multiplicity of sends
  integer, intent(in) :: ncx, ncy, nmfx, nmfy          ! array sizes
  integer, intent(in) :: cxtag(mlt), cytag(mlt)        ! tags
  integer, intent(in) :: mfxtag(mlt), mfytag(mlt)      ! tags
  integer, intent(in) :: iremote(mlt)                  ! target tasks
  integer, intent(in) :: cxreqs(mlt), cyreqs(mlt)      ! mpi requests
  integer, intent(in) :: mfxreqs(mlt), mfyreqs(mlt)    ! mpi requests


  real(r8), intent(in) :: pi
  real(r8), intent(in) :: ae                ! Radius of the Earth (m)
  real(r8), intent(in) :: om                ! rotation rate
  real(r8), intent(in) :: ptopin
  real(r8), intent(in) :: umax
  real(r8), intent(in) :: dt                !small time step in seconds
  real(r8), intent(in) :: cp
  real(r8), intent(in) :: akap

! Input time independent arrays:
  real(r8), intent(in) ::        &
    hs(grid%im,grid%jfirst:grid%jlast)         !surface geopotential
  real(r8), intent(in) ::        &
    hsxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy) !surface geopotential XY-decomp.

! !INPUT/OUTPUT PARAMETERS:

  real(r8), intent(inout) ::   &   
    u(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_s,grid%kfirst:grid%klast) ! u-Wind (m/s)
  real(r8), intent(inout) ::   &   
    v(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d,grid%kfirst:grid%klast) ! v-Wind (m/s)

  real(r8), intent(inout) ::   & 
    delp(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)          ! Delta pressure (pascal)
  real(r8), intent(inout) ::   &
    pt(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)! Scaled-Pot. temp.

! Input/output: accumulated winds & mass fluxes on c-grid for large-
!               time-step transport
  real(r8), intent(inout) ::   &
    cx3(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)! Accum. Courant no. in X
  real(r8), intent(inout) ::   &
    cy3(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)        ! Accumulated Courant no. in Y
  real(r8), intent(inout) ::   &
    mfx(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)          ! Mass flux in X  (unghosted)
  real(r8), intent(inout) ::   &
    mfy(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast)        ! Mass flux in Y

! Input/output work arrays:
  real(r8), intent(inout) ::   &
    delpf(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)  ! filtered delp
  real(r8), intent(inout) ::   &
    uc(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)   ! u-Winds on C-grid
  real(r8), intent(inout) ::   &
    vc(grid%im,grid%jfirst-2:   grid%jlast+2,   grid%kfirst:grid%klast)   ! v-Winds on C-grid

  real(r8), intent(inout) ::   &
    dpt(grid%im,grid%jfirst-1:grid%jlast+1,grid%kfirst:grid%klast)
  real(r8), intent(inout) ::   &
    wz3(grid%im,grid%jfirst-1:grid%jlast  ,grid%kfirst:grid%klast+1)
  real(r8), intent(inout) ::   &
    pkc(grid%im,grid%jfirst-1:grid%jlast+1,grid%kfirst:grid%klast+1) 
  real(r8), intent(inout) ::   &
    wz(grid%im,grid%jfirst-1:grid%jlast+1,grid%kfirst:grid%klast+1)
  real(r8), intent(inout) ::   &
    pkcc(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1) 
  real(r8), intent(inout) ::   &
    wzc(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1) 
  real(r8), intent(inout) ::   &
    wzxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1)
  real(r8), intent(inout) ::   &
    delpxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
  real(r8), intent(inout) ::   &
    pkkp(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)
  real(r8), intent(inout) ::   &
    wzkp(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)

! !OUTPUT PARAMETERS:
  real(r8), intent(out) ::   & 
    pe(grid%im,grid%kfirst:grid%klast+1,grid%jfirst:grid%jlast)         ! Edge pressure (pascal)
  real(r8), intent(out) ::   &
    pk(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)         ! Pressure to the kappa
  real(r8), intent(out) ::   & 
    ptxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) ! Potential temperature XY decomp
  real(r8), intent(out) ::   & 
    pkxy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1) ! P-to-the-kappa XY decomp
  real(r8), intent(out) ::   &
    pexy(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy) ! Edge pressure XY decomp
  real(r8), intent(out) ::   &
    ptc(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
  real(r8), intent(out) ::   &
    ptk(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
! Work arrays

! ! !DESCRIPTION:
!    Perform a dynamical update for one small time step; the small
!    time step is limitted by the fastest wave within the Lagrangian control-
!    volume 
!
! !REVISION HISTORY:
!     SJL  99.01.01:   Original SMP version
!     WS   99.04.13:   Added jfirst:jlast concept
!     SJL  99.07.15:   Merged c_core and d_core to this routine
!     WS   99.09.07:   Restructuring, cleaning, documentation
!     WS   99.10.18:   Walkthrough corrections; frozen for 1.0.7
!     WS   99.11.23:   Pruning of some 2-D arrays
!     SJL  99.12.23:   More comments; general optimization; reduction
!                      of redundant computation & communication
!     WS   00.05.14:   Modified ghost indices per Kevin's definition
!     WS   00.07.13:   Changed PILGRIM API
!     WS   00.08.28:   Cosmetic changes: removed old loop limit comments
!     AAM  00.08.30:   Introduced kfirst,klast
!     WS   00.12.01:   Replaced MPI_ON with SPMD; hs now distributed
!     WS   01.04.11:   PILGRIM optimizations for begin/endtransfer
!     WS   01.05.08:   Optimizations in the call of c_sw and d_sw
!     AAM  01.06.27:   Reinstituted 2D decomposition for use in ccm
!     WS   01.12.10:   Ghosted PT, code now uses mod_comm primitives
!     WS   01.12.31:   Removed vorticity damping, ghosted U,V,PT
!     WS   02.01.15:   Completed transition to mod_comm
!     WS   02.07.04:   Fixed 2D decomposition bug dest/src for mp_send3d
!     WS   02.09.04:   Integrated fvgcm-1_3_71 zero diff. changes by Lin
!     WS   03.07.22:   Removed HIGH_P option; this is outdated
!     WS   03.10.15:   Fixed hack of 00.04.13 for JORD>1 JCD=1, in clean way
!     WS   03.12.03:   Added grid as argument, some dynamics_vars removed
!     WS   04.08.25:   Interface simplified with GRID argument
!     WS   04.10.07:   Removed dependency on spmd_dyn; info now in GRID 
!     WS   05.05.24:   Incorporated OFFLINE_DYN; merge of CAM/GEOS5
!     PW   05.07.26:   Changes for Cray X1
!     PW   05.10.12:   More changes for Cray X1(E), avoiding array segment copying
!     WS   06.09.08:   Isolated magic numbers as F90 parameters
!     WS   06.09.15:   PI now passed as argument
!     CC   07.01.29:   Corrected calculation of OMEGA
!     PW   08.06.29:   Added options to call geopk_d and swap-based transposes
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local 2D arrays:
      real(r8) ::  wk(grid%im+2,grid%jfirst:  grid%jlast+2)
      real(r8) ::  wk1(grid%im,grid%jfirst-1:grid%jlast+1)
      real(r8) ::  wk2(grid%im+1,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d)
      real(r8) ::  wk3(grid%im,grid%jfirst-1:grid%jlast+1)

      real(r8) :: p1d(grid%im)

! fvitt    cell centered u- and v-Winds (m/s)
      real(r8) ::  u_cen(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)
      real(r8) ::  v_cen(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)
      real(r8) ::  ua(grid%im,grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)
      real(r8) ::  va(grid%im,grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d,grid%kfirst:grid%klast)


! Local scalars

      real(r8), parameter ::  D0_0                    =   0.0_r8
      real(r8), parameter ::  D0_1                    =   0.1_r8
      real(r8), parameter ::  D0_5                    =   0.5_r8
      real(r8), parameter ::  D1_0                    =   1.0_r8
      real(r8), parameter ::  D4_0                    =   4.0_r8
      real(r8), parameter ::  D8_0                    =   8.0_r8
      real(r8), parameter ::  D10_0                   =  10.0_r8
      real(r8), parameter ::  D128_0                  = 128.0_r8
      real(r8), parameter ::  D180_0                  = 180.0_r8
      real(r8), parameter ::  D1E5                    = 1.0e5_r8

      real(r8), parameter ::  ratmax                  = 0.81_r8
      real(r8), parameter ::  tiny                    = 1.0e-10_r8

      real(r8) ::  press
      real(r8) ::  rat, ycrit
      real(r8) ::  dt5

      integer :: msgtag             ! MPI message tag

      integer :: im, jm, km         ! problem dimensions
      integer :: nq                 ! # of tracers to be advected by trac2d
      integer :: ifirstxy,ilastxy   ! xy-decomp. longitude ranges
      integer :: jfirstxy,jlastxy   ! xy-decomp. latitude ranges
      integer :: ng_c               ! ghost latitudes on C grid
      integer :: ng_d               ! ghost lats on D (Max NS dependencies, ng_d >= ng_c)
      integer :: ng_s               ! max(ng_c+1,ng_d) significant if ng_c = ng_d

      integer :: jfirst
      integer :: jlast
      integer :: kfirst
      integer :: klast
      integer :: klastp            ! klast, except km+1 when klast=km

      integer :: iam
      integer :: npr_y
      integer :: npes_xy
      integer :: npes_yz

      integer i, j, k, ml
      integer js1g1, js2g0, js2g1, jn2g1
      integer jn2g0, jn1g1
      integer iord , jord
      integer ktot, ktotp

      real(r8) ::  tau, fac, pk4
      real(r8) ::  tau4 ! coefficient for 4th-order divergence damping

#if defined( SPMD )
      integer dest, src
#endif

  logical :: reset_winds = .false.
  logical :: everytime = .false.
  !
  ! set damping options:
  !
  ! - ldel2: 2nd-order velocity-component damping targetted to top layers,
  !          with coefficient del2coef (default 3E5)
  !
  ! - ldiv2: 2nd-order divergence damping everywhere and increasing in top layers
  !          (default cam3.5 setting)
  !
  ! - ldiv4: 4th-order divergence damping everywhere and increasing in top layers
  !
  ! - div24del2flag: 2 for ldiv2 (default), 4 for ldiv4, 42 for ldiv4 + ldel2
  ! - ldiv2 and ldel2 cannot coexist
  !
  logical :: ldiv2 = .true.  
  logical :: ldiv4 = .false.  
  logical :: ldel2 = .false.   


! C.-C. Chen, omega calculation
  real(r8), intent(out) ::   &
    cx_om(grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)   ! Courant in X
  real(r8), intent(out) ::  &
    cy_om(grid%im,grid%jfirst:grid%jlast+1,grid%kfirst:grid%klast) ! Courant in Y

!******************************************************************
!******************************************************************
!
! IMPORTANT CODE OPTIONS - SEE BELOW
!
!******************************************************************
!******************************************************************

! Option for which version of geopk to use with yz decomposition.
! If geopkdist=false, variables are transposed to/from xy decomposition
!   for use in geopk.
! If geopkdist=true, either geopk_d or geopk16 is used. Both 
!   compute local partial sums in z and then communicate those 
!   sums to combine them. geopk_d does not try to parallelize in the
!   z-direction except in a pipeline fashion controlled by the
!   parameter geopkblocks, and is bit-for-bit the same as the
!   transpose-based algorithm. geopk16 exploits z-direction 
!   parallelism and requires 16-byte arithmetic (DSIZE=16) 
!   to reproduce the same numerics (and to be reproducible with
!   respect to process count). The geopk16 default is to use 
!   8-byte arithmetic (DSIZE=8). This is faster than
!   16-byte, but also gives up reproducibility. On many systems
!   performance of geopk_d is comparable to geopk16 even with 
!   8-byte numerics.
! On the last two small timesteps (ipe=1,2 or 1,-2) for D-grid, 
!   the version of geopk that uses transposes is called regardless, 
!   as some transposed quantities are required for the te_map phase
!   and for the calculation of omega.
! For non-SPMD mode, geopk_[cd]dist are set to false.

      logical geopk_cdist, geopk_ddist

      geopk_cdist = .false.
      geopk_ddist = .false.
#if defined( SPMD )
      if (grid%geopkdist) then
        geopk_cdist = .true.
        if ((ipe == -1) .or. (ipe == 0)) geopk_ddist = .true.
      endif
#endif

!******************************************************************

      npes_xy  = grid%npes_xy
      npes_yz  = grid%npes_yz

      im       = grid%im
      jm       = grid%jm
      km       = grid%km
      nq       = grid%nq

      ng_c     = grid%ng_c
      ng_d     = grid%ng_d
      ng_s     = grid%ng_s

      jfirst   = grid%jfirst
      jlast    = grid%jlast
      kfirst   = grid%kfirst
      klast    = grid%klast
      klastp   = grid%klastp

      iam      = grid%iam
      npr_y    = grid%npr_y

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      ktot = klast - kfirst + 1
      ktotp = ktot + 1

    if (iam .lt. npes_yz) then

      call FVstartclock(grid,'---PRE_C_CORE')

#if defined( SPMD )
      call FVstartclock(grid,'---PRE_C_CORE_COMM')
      call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,  &
                         kfirst, klast, ng_d, ng_s, u )
      call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,  &
                         kfirst, klast, ng_s, ng_d, v )
      call FVstopclock(grid,'---PRE_C_CORE_COMM')
#endif

! Set general loop limits
! jfirst >= 1; jlast <= jm
      js1g1  = max(1,jfirst-1)
      js2g0  = max(2,jfirst)
      js2g1  = max(2,jfirst-1)
      jn2g0  = min(jm-1,jlast)
      jn1g1  = min(jm,jlast+1)
      jn2g1 = min(jm-1,jlast+1)

      if( abs(grid%dt0-dt) > D0_1 ) then

        grid%dt0 = dt
        dt5 = D0_5*dt

        grid%rdy   = D1_0/(ae*grid%dp)
        grid%dtdy  = dt *grid%rdy
        grid%dtdy5 = dt5*grid%rdy
        grid%dydt  = (ae*grid%dp) / dt
        grid%tdy5  = D0_5/grid%dtdy

        do j=2,jm-1
          grid%dx(j)    = grid%dl*ae*grid%cosp(j)
          grid%rdx(j)   = D1_0 / grid%dx(j)
          grid%dtdx(j)  = dt /grid% dx(j)
          grid%dxdt(j)  = grid%dx(j) / dt
          grid%dtdx2(j) = D0_5*grid%dtdx(j)
          grid%dtdx4(j) = D0_5*grid%dtdx2(j)
          grid%dycp(j)  = ae*grid%dp/grid%cosp(j)
          grid%cy(j)    = grid%rdy * grid%acosp(j)
        enddo

        do j=2,jm
          grid%dxe(j)   = ae*grid%dl*grid%cose(j)
          grid%rdxe(j)  = D1_0 / grid%dxe(j)
          grid%dtdxe(j) = dt / grid%dxe(j)
          grid%dtxe5(j) = D0_5*grid%dtdxe(j)
          grid%txe5(j)  = D0_5/grid%dtdxe(j)
          grid%cye(j)   =  D1_0 / (ae*grid%cose(j)*grid%dp)
          grid%dyce(j)  = ae*grid%dp/grid%cose(j)
        enddo

! C-grid
#ifndef WACCM_MOZART
        grid%zt_c = abs(umax*dt5) / (grid%dl*ae)
#else
        grid%zt_c = cos( D10_0 * pi / D180_0 )
#endif

! D-grid
#ifndef WACCM_MOZART
        grid%zt_d = abs(umax*dt) / (grid%dl*ae)
#else
        grid%zt_d = cos( D10_0 * pi / D180_0 )
#endif

        if ( ptopin /= grid%ptop) then
             write(iulog,*) 'PTOP as input to cd_core != ptop from T_FVDYCORE_GRID'
             stop
        endif

        !
        ! damping code
        !
        if (div24del2flag == 2) then
          !
          ! cam3.5 default damping setting
          !
          ldiv2 = .true.
          ldiv4 = .false.
          ldel2 = .false.
          if (masterproc)  write(iulog,*) 'Divergence damping: use 2nd order damping'
        elseif (div24del2flag == 4) then
          !
          ! fourth order divergence damping and no velocity diffusion
          !
          ldiv2 = .false.
          ldiv4 = .true.
          ldel2 = .false.
          if (masterproc) write(iulog,*) 'Divergence damping: use 4th order damping'
        elseif (div24del2flag == 42) then
          !
          ! fourth order divergence damping with velocity diffusion
          !
          ldiv2 = .false.
          ldiv4 = .true.
          ldel2 = .true.
          if (masterproc) write(iulog,*) 'Divergence damping: use 4th order damping'
          if (masterproc) write(iulog,*) 'Velocity del2 damping with coefficient ', del2coef
        else
          ldiv2 = .true.
          ldiv4 = .false.
          ldel2 = .false.
          if (masterproc) write(iulog,*) 'Inadmissable velocity smoothing option - div24del2flag = ', div24del2flag
          call endrun('Inadmissable value of div24del2flag')
        endif

        do k=kfirst,klast
         
          if (ldel2) then
            !
            !***********************************
            !
            ! Laplacian on velocity components
            !
            !***********************************
            !
            press = D0_5 * ( grid%ak(k)+grid%ak(k+1) + &
                 (grid%bk(k)+grid%bk(k+1))*D1E5 )
            tau = D8_0 * (D1_0+ tanh(D1_0*log(grid%ptop/press)) )
            !
            ! tau is strength of damping
            !
            if (tau < 0.3_r8) then
              !
              ! no del2 damping at lower levels
              !
              tau = 0.0_r8
            end if

            do j=js2g0,jn1g1
              !
              ! fac must include dt for the momentum equation
              ! i.e. diffusion coefficient is fac/dt
              !
              ! del2 diffusion coefficient in spectral core is 2.5e5             
              !
              fac =  tau * dt * del2coef
              !
              ! all these coefficients are necessary because of the staggering of the
              ! wind components
              !
              grid%cdxde(j,k) = fac/(ae*ae*grid%cose(j)*grid%cose(j)*grid%dl*grid%dl) 
              grid%cdyde(j,k) = fac/(ae*ae*grid%cose(j)*grid%dp*grid%dp)
            end do
            do j=js2g0,jn2g1
              fac =  tau * dt * del2coef
              grid%cdxdp(j,k) = fac/(ae*ae*grid%cosp(j)*grid%cosp(j)*grid%dl*grid%dl) 
              grid%cdydp(j,k) = fac/(ae*ae*grid%cosp(j)*grid%dp*grid%dp)
            end do
          end if

          if (ldiv2) then
            !
            !***********************************************
            !
            ! cam3 default second-order divergence damping
            !
            !***********************************************
            !
            press = D0_5 * ( grid%ak(k)+grid%ak(k+1) + &
                 (grid%bk(k)+grid%bk(k+1))*D1E5 )
            tau = D8_0 * (D1_0+ tanh(D1_0*log(grid%ptop/press)) )
            tau = max(D1_0, tau) / (D128_0*abs(dt))
            do j=js2g0,jn1g1
              !-----------------------------------------
              ! Explanation of divergence damping coeff.
              ! ========================================
              !
              ! Divergence damping is added to the momentum 
              ! equations through a term tau*div where
              !
              !       tau = C*L**2/dt 
              !
              ! where L is the length scale given by
              ! 
              !       L**2 = a**2*dl*dp
              !
              ! and divergence is given by
              !
              !       div = divx + divy
              !
              ! where
              !
              !       divx = (1/(a*cos(p)))*du/dl
              !       divy = (1/(a*cos(p)))*(d(cos(theta)*v)/dp))
              ! 
              ! du and (d(cos(theta*v)/dp)) are computed in sw_core
              !
              ! The constant terms in divx*tau and divy*tau are
              !
              ! cdx = (1/(a*cos(p)))* (1/dl) * C * a**2 * dl * dp / dt = C * (a*dp/(cos(p)))/dt
              ! cdy = (1/(a*cos(p)))* (1/dp) * C * a**2 * dl * dp / dt = C * (a*dl/(cos(p)))/dt
              !
              !-----------------------------------------
              fac = tau * ae / grid%cose(j) !default
              grid%cdx(j,k) = fac*grid%dp   !default
              grid%cdy(j,k) = fac*grid%dl   !default
            end do
          end if

          if (ldiv4) then
            !
            ! 4th-order divergence damping
            !            
            tau4 = 0.01_r8 / (abs(dt)) 
            !
            !**************************************
            !
            ! fourth order divergence damping
            !
            !**************************************
            !
            do j=1,jm
              !
              ! divergence computation coefficients
              !
              grid%cdxdiv  (j,k) = D1_0/(grid%cose(j)*grid%dl)
              grid%cdydiv  (j,k) = D1_0/(grid%cose(j)*grid%dp)
            end do
            do j=js2g0,jn1g1
              !
              ! div4 coefficients
              !
              fac = grid%dl*grid%cose(j)!*ae
              grid%cdx4   (j,k) = D1_0/(fac*fac)
              fac = grid%dp*grid%dp*grid%cose(j)!*ae*ae
              grid%cdy4   (j,k) = D1_0/fac
              fac = grid%cose(j)*grid%dp*grid%dl
              grid%cdtau4(j,k)  = -ae*tau4*fac*fac
            end do
          endif
        end do
      end if


      if ( ipe < 0 .or. ns == 1 ) then          ! starting cd_core
         call FVstartclock(grid,'---C_DELP_LOOP')
!$omp parallel do private(i, j, k, wk, wk2)
#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (I, J, K, WK, WK2)
#endif
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  delpf(i,j,k) = delp(i,j,k)
               enddo
            enddo
            call pft2d( delpf(1,js2g0,k), grid%sc, &
                        grid%dc, im, jn2g0-js2g0+1,    &
                        wk, wk2 )
         enddo
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif
         call FVstopclock(grid,'---C_DELP_LOOP')

      endif

#if defined( SPMD )
      call FVstartclock(grid,'---PRE_C_CORE_COMM')
      call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,    &
                         kfirst, klast, ng_d, ng_s, u )
      call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,    &
                         kfirst, klast, ng_s, ng_d, v )

      call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,    &
                         kfirst, klast, ng_d, ng_d, pt )
      if ( ipe < 0 .or. ns == 1 ) then          ! starting cd_core
         call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast, &
                            kfirst, klast, ng_d, ng_d, delpf )
      endif                         ! end if ipe < 0 check
      call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,    &
                         kfirst, klast, ng_d, ng_d, pt )
      if ( ipe < 0 .or. ns == 1 ) then          ! starting cd_core
         call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast, &
                            kfirst, klast, ng_d, ng_d, delpf )
      endif                          ! end if ipe < 0 check
      call FVstopclock(grid,'---PRE_C_CORE_COMM')
#endif

!
! Get the cell centered winds if needed for the sub-step
!
#if ( defined OFFLINE_DYN )
      if ( ( (ipe < 0) .or. (everytime) ) .and. (.not. met_winds_on_walls()) ) then
         call get_met_fields( grid, u_cen, v_cen )
         reset_winds = .true.
      else
         reset_winds = .false.
      endif
#endif


! Get D-grid V-wind at the poles and interpolate winds to A- and C-grids;
! This calculation was formerly done in subroutine c_sw but is being done here to
! avoid communication in OpenMP loops

!$omp parallel do private(k, wk, wk2)

#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (K, WK, WK2)
#endif
      do  k=kfirst,klast
         call d2a2c_winds(grid, u(1,jfirst-ng_d,k),         v(1,jfirst-ng_s,k),     &
                                ua(1,jfirst-ng_d,k),       va(1,jfirst-ng_s,k),     &
                                uc(1,jfirst-ng_d,k),       vc(1,jfirst-2,k),        &
                                u_cen(1,jfirst-ng_d,k), v_cen(1,jfirst-ng_s,k),     &
                                reset_winds, met_rlx(k)  )

! Optionally filter advecting C-grid winds
         if (filtcw .gt. 0) then
            call pft2d(uc(1,js2g0,k), grid%sc, grid%dc, im, jn2g0-js2g0+1, wk, wk2 )
            call pft2d(vc(1,js2g0,k), grid%se, grid%de, im, jlast-js2g0+1, wk, wk2 )
         endif

      enddo
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif

! Fill C-grid advecting winds Halo regions
! vc only needs to be ghosted at jlast+1
#if defined( SPMD )
      call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,               &
                         kfirst, klast, ng_d, ng_d, uc )
      call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,               &
                         kfirst, klast, 2, 2, vc )
      call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,               &
                         kfirst, klast, ng_d, ng_d, uc )
      call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,               &
                         kfirst, klast, 2, 2, vc )
#endif

      call FVstopclock(grid,'---PRE_C_CORE')

      call FVbarrierclock(grid,'sync_c_core', grid%commyz)
      call FVstartclock(grid,'---C_CORE')

#if !defined(INNER_OMP)
!$omp parallel do private(i, j, k, iord, jord)    
#endif

#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (K, IORD, JORD)
#endif

      do  k=kfirst,klast     ! This is the main parallel loop.

         if ( k <= km/8 ) then
             iord = 1
             jord = 1
         else
             iord = iord_c
             jord = jord_c
         endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the C-grid
!-----------------------------------------------------------------

         call c_sw( grid, u(1,jfirst-ng_d,k),   v(1,jfirst-ng_s,k),   &
                    pt(1,jfirst-ng_d,k),  delp(1,jfirst,k),           &
                    ua(1,jfirst-ng_d,k),  va(1,jfirst-ng_s,k),        &
                    uc(1,jfirst-ng_d,k),  vc(1,jfirst-2,k),           &
                    ptc(1,jfirst,k),      delpf(1,jfirst-ng_d,k),     &
                    ptk(1,jfirst,k),  tiny, iord,   jord)
      enddo
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif

      call FVstopclock(grid,'---C_CORE')

! MPI note: uc, vc, ptk, and ptc computed within the above k-look from jfirst to jlast
! Needed by D-core: uc(jfirst-ng_d:jlast+ng_d), vc(jfirst:jlast+1)

      call FVbarrierclock(grid,'sync_c_geop', grid%commyz)

    end if  !  (iam .lt. npes_yz)

      if (geopk_cdist) then

    if (iam .lt. npes_yz) then

!
! Stay in yz space and use z communications
!

      if (grid%geopk16byte) then
         call FVstartclock(grid,'---C_GEOP16')
         call geopk16(grid, pe, ptk, pkcc, wzc, hs, ptc, &
                    0, cp, akap)
      else
         call FVstartclock(grid,'---C_GEOP_D')
         call geopk_d(grid, pe, ptk, pkcc, wzc, hs, ptc, &
                      0, cp, akap)
      endif

!
! Geopk does not need j ghost zones of pkc and wz
!

!$omp parallel do private(i, j, k)
         do k = kfirst, klast+1
           do j = jfirst, jlast
              do i = 1, im
               pkc(i,j,k) = pkcc(i,j,k)
               wz(i,j,k) = wzc(i,j,k)
             enddo
           enddo
         enddo

      if (grid%geopk16byte) then
         call FVstopclock(grid,'---C_GEOP16')
      else
         call FVstopclock(grid,'---C_GEOP_D')
      endif

    end if  !  (iam .lt. npes_yz)

      else

! Begin xy geopotential section

         call FVstartclock(grid,'---C_GEOP')

         if (grid%twod_decomp == 1) then

!
! Transpose to xy decomposition
!

#if defined( SPMD )
            call FVstartclock(grid,'YZ_TO_XY_C_GEOP')
            if (grid%modc_onetwo .eq. 1) then
               call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptk, delpxy,     &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptk, delpxy,     &
                                modc=grid%modc_cdcore )
               call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptc, ptxy,       &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptc, ptxy,       &
                                modc=grid%modc_cdcore )
            else
               call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptk, delpxy,     &
                                ptc, ptxy,                                   &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptk, delpxy,     &
                                ptc, ptxy,                                   &
                                modc=grid%modc_cdcore )
            endif
            call FVstopclock(grid,'YZ_TO_XY_C_GEOP')
#endif

         else

!$omp parallel do private(i, j, k)
            do k = kfirst, klast
              do j = jfirst, jlast
                do i = 1, im
                  delpxy(i,j,k) = ptk(i,j,k)
                  ptxy(i,j,k) = ptc(i,j,k)
                enddo
              enddo
            enddo

         endif

         call geopk(grid, pexy, delpxy, pkxy, wzxy, hsxy, ptxy, &
                    cp, akap, nx)

         if (grid%twod_decomp == 1) then
!
! Transpose back to yz decomposition.
! pexy is not output quantity on this call.
! pkkp and wzkp are holding arrays, whose specific z-dimensions
!    are required by Pilgrim.
! Z edge ghost points (klast+1) are automatically filled in
!

#if defined( SPMD )

            call FVstartclock(grid,'XY_TO_YZ_C_GEOP')
            if (grid%modc_onetwo .eq. 1) then
               call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                modc=grid%modc_cdcore )
               call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, wzxy, wzkp,       &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, wzxy, wzkp,       &
                                modc=grid%modc_cdcore )
            else
               call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                wzxy, wzkp,                                  &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                wzxy, wzkp,                                  &
                                modc=grid%modc_cdcore )
            endif
            call FVstopclock(grid,'XY_TO_YZ_C_GEOP')

!$omp parallel do private(i, j, k)
            do k = kfirst, klast+1
               do j = jfirst, jlast
                  do i = 1, im
                     pkc(i,j,k) = pkkp(i,j,k)
                  enddo
               enddo
            enddo

!$omp parallel do private(i, j, k)
            do k = kfirst, klast+1
               do j = jfirst, jlast
                  do i = 1, im
                     wz(i,j,k) = wzkp(i,j,k)
                  enddo
               enddo
            enddo

#endif

         else

!$omp parallel do private(i, j, k)
            do k = kfirst, klast+1
               do j = jfirst, jlast
                  do i = 1, im
                     pkc(i,j,k) = pkxy(i,j,k)
                     wz(i,j,k) = wzxy(i,j,k)
                  enddo
               enddo
            enddo

         endif

         call FVstopclock(grid,'---C_GEOP')

! End xy geopotential section

      endif       ! geopk_cdist

    if (iam .lt. npes_yz) then

      call FVbarrierclock(grid,'sync_pre_d_core', grid%commyz)
      call FVstartclock(grid,'---PRE_D_CORE')

! Upon exit from geopk, the quantities pe, pkc and wz will have been
! updated at klast+1


#if defined( SPMD )
!
! pkc & wz need to be ghosted only at jfirst-1
!
      call FVstartclock(grid,'---PRE_D_CORE_COMM')
      dest = iam+1
      src  = iam-1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      if ( mod(iam,npr_y) == 0 ) src = -1
      call mp_send3d_2( grid%commyz, dest, src, im, jm, km+1,          &
                        1, im, jfirst-1, jlast+1, kfirst, klast+1,    &
                        1, im, jlast, jlast, kfirst, klast+1, pkc, wz)
      call FVstopclock(grid,'---PRE_D_CORE_COMM')
#endif


      call FVstartclock(grid,'---C_U_LOOP')
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k, p1d, wk, wk2)

#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (I, J, K, P1D, WK, WK2)
#endif
      do k=kfirst,klast
         do j=js2g0,jn2g0
            do i=1,im
               p1d(i) = pkc(i,j,k+1) - pkc(i,j,k)
            enddo

            uc(1,j,k) = uc(1,j,k) + grid%dtdx2(j) * (                   &
                   (wz(im,j,k+1)-wz(1,j,k))*(pkc(1,j,k+1)-pkc(im,j,k))  &
                 + (wz(im,j,k)-wz(1,j,k+1))*(pkc(im,j,k+1)-pkc(1,j,k))) &
                                   / (p1d(1)+p1d(im))
            do i=2,im
               uc(i,j,k) = uc(i,j,k) + grid%dtdx2(j) * (                &
                 (wz(i-1,j,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i-1,j,k))  &
               + (wz(i-1,j,k)-wz(i,j,k+1))*(pkc(i-1,j,k+1)-pkc(i,j,k))) &
                                    / (p1d(i)+p1d(i-1))
            enddo

! C.-C. Chen
            do i=1,im
               cx_om(i,j,k) = grid%dtdx(j)*uc(i,j,k)
            enddo
         enddo
         call pft2d(uc(1,js2g0,k), grid%sc,       &
                    grid%dc, im, jn2g0-js2g0+1,       &
                    wk, wk2 )
         if ( jfirst == 1 ) then   ! Clean up
            do i=1,im
               uc(i,1,k) = D0_0
               cx_om(i,1,k) = D0_0
            enddo
         endif
         if ( jlast == jm ) then   ! Clean up
            do i=1,im
               uc(i,jm,k) = D0_0
               cx_om(i,jm,k) = D0_0
            enddo
         endif

      enddo 
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif
      call FVstopclock(grid,'---C_U_LOOP')

#if defined( SPMD )
      call FVstartclock(grid,'---PRE_D_CORE_COMM')
      call mp_recv3d_2( grid%commyz, src, im, jm, km+1,                &
                        1, im, jfirst-1, jlast+1, kfirst, klast+1,    &
                        1, im, jfirst-1, jfirst-1, kfirst, klast+1, pkc, wz)

      call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,    &
                         kfirst, klast, ng_d, ng_d, uc )
      call FVstopclock(grid,'---PRE_D_CORE_COMM')
#endif

      call FVstartclock(grid,'---C_V_PGRAD')
!
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k, wk, wk1 )

! pkc and wz need only to be ghosted jfirst-1

#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (I, J, K, WK, WK1 )
#endif
      do k=kfirst,klast
         do j=js1g1,jlast
            do i=1,im
               wk1(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
            enddo
         enddo

         do j=js2g0,jlast
            do i=1,im
               vc(i,j,k) = vc(i,j,k) + grid%dtdy5/(wk1(i,j)+wk1(i,j-1)) *  &
               ( (wz(i,j-1,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j-1,k))     &
               +  (wz(i,j-1,k)-wz(i,j,k+1))*(pkc(i,j-1,k+1)-pkc(i,j,k)) )

! C.-C. Chen
               cy_om(i,j,k) = grid%dtdy*vc(i,j,k)
            enddo
         enddo

         call pft2d(vc(1,js2g0,k), grid%se,          &
                    grid%de, im, jlast-js2g0+1, wk, wk1 )
      enddo
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif

      call FVstopclock(grid,'---C_V_PGRAD')

#if defined( SPMD )
      call FVstartclock(grid,'---PRE_D_CORE_COMM')
      call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,   &
                         kfirst, klast, ng_d, ng_d, uc )

! vc only needs to be ghosted at jlast+1
      dest = iam-1
      src  = iam+1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( grid%commyz, dest, src, im, jm, km,             &
                      1, im, jfirst-2, jlast+2, kfirst, klast,       &
                      1, im, jfirst, jfirst, kfirst, klast, vc )            
      call mp_recv3d( grid%commyz, src, im, jm, km,                   &
                      1, im, jfirst-2, jlast+2, kfirst, klast,       &
                      1, im, jlast+1, jlast+1, kfirst, klast, vc )
      call FVstopclock(grid,'---PRE_D_CORE_COMM')

! C.-C. Chen
      call mp_send3d( grid%commyz, dest, src, im, jm, km,                      &
                      1, im, jfirst, jlast+1, kfirst, klast,       &
                      1, im, jfirst, jfirst, kfirst, klast, cy_om )            
      call mp_recv3d( grid%commyz, src, im, jm, km,                             &
                      1, im, jfirst, jlast+1, kfirst, klast,       &
                      1, im, jlast+1, jlast+1, kfirst, klast, cy_om )
#endif

      call FVstopclock(grid,'---PRE_D_CORE')

      call FVbarrierclock(grid,'sync_d_core', grid%commyz)
      call FVstartclock(grid,'---D_CORE')

#if !defined(INNER_OMP)
!$omp parallel do private(i, j, k, iord, jord) 
#endif   
#if !defined(USE_OMP)
!CSD$ PARALLEL DO PRIVATE (K, IORD, JORD)
#endif

      do k=kfirst,klast

         if( k <= km/8 ) then
            if( k == 1 ) then
               iord = 1
               jord = 1
            else
               iord = min(2, iord_d)
               jord = min(2, jord_d)
            endif
         else
            iord = iord_d
            jord = jord_d
         endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the D-grid
!-----------------------------------------------------------------

         call d_sw( grid, u(1,jfirst-ng_d,k),      v(1,jfirst-ng_s,k),     &
                    uc(1,jfirst-ng_d,k),    vc(1,jfirst-2,k),        &
                    pt(1,jfirst-ng_d,k),   delp(1,jfirst,k),         &
                    delpf(1,jfirst-ng_d,k), cx3(1,jfirst-ng_d,k),    &
                    cy3(1,jfirst,k),        mfx(1,jfirst,k),         &
                    mfy(1,jfirst,k),              &
                    grid%cdx  (js2g0:,k),grid%cdy (js2g0:,k),        &
                    grid%cdxde (js2g0:,k),grid%cdxdp (js2g0:,k),     & 
                    grid%cdyde(js2g0:,k) ,grid%cdydp(js2g0:,k),      & 
                    grid%cdxdiv(:,k),grid%cdydiv(:,k) ,              & 
                    grid%cdx4 (js2g0:,k),grid%cdy4(js2g0:,k) ,       & 
                    grid%cdtau4(js2g0:,k), ldiv2, ldiv4, ldel2,      & 
                    iord, jord, tiny )

      enddo
#if !defined(USE_OMP)
!CSD$ END PARALLEL DO
#endif

      call FVstopclock(grid,'---D_CORE')

      call FVbarrierclock(grid,'sync_d_geop', grid%commyz)

#if defined( SPMD )
      if (s_trac) then
! post sends for ct_overlap or tracer decomposition information
         do ml = 1, mlt
            call mpiisend(cx3, ncx, mpir8, iremote(ml), cxtag(ml), grid%commnyz, cxreqs(ml))
            call mpiisend(cy3, ncy, mpir8, iremote(ml), cytag(ml), grid%commnyz, cyreqs(ml))
            call mpiisend(mfx, nmfx, mpir8, iremote(ml), mfxtag(ml), grid%commnyz, mfxreqs(ml))
            call mpiisend(mfy, nmfy, mpir8, iremote(ml), mfytag(ml), grid%commnyz, mfyreqs(ml))
         enddo
      endif
#endif

    end if  !  (iam .lt. npes_yz)

      if (geopk_ddist) then

    if (iam .lt. npes_yz) then

!
! Stay in yz space and use z communications

      if (grid%geopk16byte) then
         call FVstartclock(grid,'---D_GEOP16')
         call geopk16(grid, pe, delp, pkcc, wzc, hs, pt, &
                    ng_d, cp, akap)
      else
         call FVstartclock(grid,'---D_GEOP_D')
         call geopk_d(grid, pe, delp, pkcc, wzc, hs, pt, &
                      ng_d, cp, akap)
      endif

!
! Geopk does not need j ghost zones of pkc and wz
!

!$omp parallel do private(i, j, k)
      do k = kfirst, klast+1
         do j = jfirst, jlast
            do i = 1, im
               pkc(i,j,k) = pkcc(i,j,k)
               wz(i,j,k) = wzc(i,j,k)
            enddo
         enddo
      enddo

      if (grid%geopk16byte) then
         call FVstopclock(grid,'---D_GEOP16')
      else
         call FVstopclock(grid,'---D_GEOP_D')
      endif

    end if  !  (iam .lt. npes_yz)

      else

! Begin xy geopotential section

         call FVstartclock(grid,'---D_GEOP')

         if (grid%twod_decomp == 1) then
!
! Transpose to xy decomposition
!

#if defined( SPMD )

!$omp parallel do private(i,j,k)
            do k=kfirst,klast
               do j=jfirst,jlast
                  do i=1,im
                     ptc(i,j,k) = pt(i,j,k)
                  enddo
               enddo
            enddo

            call FVstartclock(grid,'YZ_TO_XY_D_GEOP')
            if (grid%modc_onetwo .eq. 1) then
               call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, delp, delpxy,    &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, delp, delpxy,    &
                                modc=grid%modc_cdcore )
               call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptc, ptxy,       &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, ptc, ptxy,       &
                                modc=grid%modc_cdcore )
            else
               call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, delp, delpxy,    &
                                ptc, ptxy,                                   &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,     &
                                grid%ijk_yz_to_xy%RecvDesc, delp, delpxy,    &
                                ptc, ptxy,                                   &
                                modc=grid%modc_cdcore )
            endif
            call FVstopclock(grid,'YZ_TO_XY_D_GEOP')
#endif

         else

!$omp parallel do private(i,j,k)
            do k=kfirst,klast
               do j=jfirst,jlast
                  do i=1,im
                     delpxy(i,j,k) = delp(i,j,k)
                     ptxy(i,j,k) = pt(i,j,k)
                  enddo
               enddo
            enddo
 
         endif

         call geopk(grid, pexy, delpxy, pkxy, wzxy, hsxy, ptxy, &
                    cp, akap, nx)

         if (grid%twod_decomp == 1) then
!
! Transpose back to yz decomposition
! Z edge ghost points (klast+1) are automatically filled in
! pexy is output quantity on last small timestep
!

#if defined( SPMD )

            call FVstartclock(grid,'XY_TO_YZ_D_GEOP')
            if (grid%modc_onetwo .eq. 1) then
               call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                modc=grid%modc_cdcore )
               call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, wzxy, wzkp,       &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, wzxy, wzkp,       &
                                modc=grid%modc_cdcore )
            else
               call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                wzxy, wzkp,                                  &
                                modc=grid%modc_cdcore )
               call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,      &
                                grid%pkxy_to_pkc%RecvDesc, pkxy, pkkp,       &
                                wzxy, wzkp,                                  &
                                modc=grid%modc_cdcore )
            endif
            call FVstopclock(grid,'XY_TO_YZ_D_GEOP')

!$omp parallel do private(i, j, k)
            do k = kfirst, klast+1
               do j = jfirst, jlast
                  do i = 1, im
                     pkc(i,j,k) = pkkp(i,j,k)
                  enddo
               enddo
            enddo

!$omp parallel do private(i, j, k)
            do k = kfirst, klast+1
               do j = jfirst, jlast
                  do i = 1, im
                     wz(i,j,k) = wzkp(i,j,k)
                  enddo
               enddo
            enddo
#endif

         else

!$omp parallel do private(i, j, k)
            do k = kfirst, klast+1
               do j = jfirst, jlast
                  do i = 1, im
                     pkc(i,j,k) = pkxy(i,j,k)
                     wz(i,j,k) = wzxy(i,j,k)
                  enddo
               enddo
            enddo

         endif

         call FVstopclock(grid,'---D_GEOP')

! End xy geopotential section

      endif       ! geopk_ddist

    if (iam .lt. npes_yz) then

      call FVbarrierclock(grid,'sync_pre_d_pgrad', grid%commyz)

!
! Upon exit from geopk, the quantities pe, pkc and wz will have been
!      updated at klast+1

      call FVstartclock(grid,'---PRE_D_PGRAD')

#if defined( SPMD )
      call FVstartclock(grid,'---PRE_D_PGRAD_COMM_1')
! Exchange boundary regions on north and south for pkc and wz
      call mp_send2_ns( grid%commyz, im, jm, km+1, jfirst, jlast,      &
                        kfirst, klast+1, 1, pkc, wz)
      call FVstopclock(grid,'---PRE_D_PGRAD_COMM_1')
#endif

      if ( ipe /= 1 ) then          !  not the last call
!
! Perform some work while sending data on the way
!

         call FVstartclock(grid,'---D_DELP_LOOP')

!$omp parallel do private(i, j, k, wk, wk2)

#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (I, J, K, WK, WK2)
#endif
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  delpf(i,j,k) = delp(i,j,k)
               enddo
            enddo
            call pft2d( delpf(1,js2g0,k), grid%sc, &
                        grid%dc, im, jn2g0-js2g0+1,    &
                        wk, wk2 )
         enddo
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif
         call FVstopclock(grid,'---D_DELP_LOOP')

      else
! Last call
!$omp parallel do private(i, j, k)
         do k=kfirst,klast+1
            do j=jfirst,jlast
               do i=1,im
                  pk(i,j,k) = pkc(i,j,k)
               enddo
            enddo
         enddo
      endif

#if defined( SPMD )
      call FVstartclock(grid,'---PRE_D_PGRAD_COMM_1')
      call mp_recv2_ns( grid%commyz, im, jm, km+1, jfirst, jlast,          &
                        kfirst, klast+1, 1, pkc, wz)
      if ( ipe /= 1 ) then          !  not the last call
         call mp_send4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,     &
                            kfirst, klast, ng_d, ng_d, delpf )
      endif
      call FVstopclock(grid,'---PRE_D_PGRAD_COMM_1')
#endif


!
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k)

      do k=kfirst,klast
         do j=js1g1,jn1g1                  ! dpt needed NS
            do i=1,im                       ! wz, pkc ghosted NS
               dpt(i,j,k)=(wz(i,j,k+1)+wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j,k))
            enddo
         enddo
      enddo

!  GHOSTING:   wz (input) NS ; pkc (input) NS

      call FVstopclock(grid,'---PRE_D_PGRAD')
      call FVstartclock(grid,'---D_PGRAD_1')

!$omp parallel do private(i, j, k, wk3, wk1)
#if !defined(USE_OMP)
!CSD$ PARALLEL DO PRIVATE (I, J, K, WK3, WK1)
#endif
      do k=kfirst,klast+1

         if (k == 1) then
            do j=js2g0,jlast
               do i=1,im
                  wz3(i,j,1) = D0_0
                  wz(i,j,1) = D0_0
               enddo
            enddo
            pk4 = D4_0*grid%ptop**akap
            do j=js2g0,jn1g1
               do i=1,im
                  pkc(i,j,1) = pk4
               enddo
            enddo
            go to 4500
         endif

         do j=js2g1,jn2g0                             ! wk3 needed S
            wk3(1,j) = (wz(1,j,k)+wz(im,j,k)) *       &
                       (pkc(1,j,k)-pkc(im,j,k))
            do i=2,im
               wk3(i,j) = (wz(i,j,k)+wz(i-1,j,k)) *      & 
                          (pkc(i,j,k)-pkc(i-1,j,k))
            enddo
         enddo

         do j=js2g1,jn2g0                               
            do i=1,im-1
               wk1(i,j) = wk3(i,j) + wk3(i+1,j)        
            enddo
            wk1(im,j) = wk3(im,j) + wk3(1,j)      ! wk3 ghosted S
         enddo

         if ( jfirst == 1 ) then
            do i=1,im
               wk1(i, 1) = D0_0
           enddo
         endif

         if ( jlast == jm ) then
            do i=1,im
               wk1(i,jm) = D0_0
            enddo
         endif

         do j=js2g0,jlast                          ! wk1 ghosted S
            do i=1,im
               wz3(i,j,k) = wk1(i,j) + wk1(i,j-1)
            enddo
         enddo

! N-S walls

         do j=js2g0,jn1g1                     ! wk1 needed N
            do i=1,im                         ! wz, pkc ghosted NS
               wk1(i,j) = (wz(i,j,k)+wz(i,j-1,k))*(pkc(i,j,k)-pkc(i,j-1,k))
            enddo
         enddo

         do j=js2g0,jn1g1                    ! wk3 needed N
            wk3(1,j) = wk1(1,j) + wk1(im,j)      ! wk1 ghosted N
            do i=2,im
               wk3(i,j) = wk1(i,j) + wk1(i-1,j)   ! wk1 ghosted N
            enddo
         enddo

         do j=js2g0,jn2g0
            do i=1,im
               wz(i,j,k) = wk3(i,j) + wk3(i,j+1)  ! wk3 ghosted N
            enddo
         enddo

         do j=js1g1,jn1g1
            wk1(1,j) = pkc(1,j,k) + pkc(im,j,k)
            do i=2,im
               wk1(i,j) = pkc(i,j,k) + pkc(i-1,j,k)
            enddo
         enddo
 
         do j=js2g0,jn1g1
            do i=1,im
               pkc(i,j,k) = wk1(i,j) + wk1(i,j-1)
            enddo
         enddo

4500     continue
      enddo

#if !defined(USE_OMP)
!CSD$ END PARALLEL DO
#endif
      call FVstopclock(grid,'---D_PGRAD_1')
      call FVstartclock(grid,'---D_PGRAD_2')

!     GHOSTING:   dpt (loop 4000) NS ; pkc (loop 4500) N
!
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k, wk, wk1, wk2, wk3)
#if (!defined USE_OMP) 
!CSD$ PARALLEL DO PRIVATE (i, j, k, wk, wk1, wk2, wk3)
#endif
      do 6000 k=kfirst,klast

         do j=js1g1,jn1g1
            wk1(1,j) = dpt(1,j,k) + dpt(im,j,k)
            do i=2,im
               wk1(i,j) = dpt(i,j,k) + dpt(i-1,j,k)
            enddo
         enddo
 
         do j=js2g0,jn1g1
            do i=1,im
               wk2(i,j) = wk1(i,j) + wk1(i,j-1)
               wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
            enddo
         enddo

         do j=js2g0,jlast
            do i=1,im-1
               wk3(i,j) = uc(i,j,k)  +  grid%dtdxe(j)/(wk(i,j) + wk(i+1,j))      &
                          * (wk2(i,j)-wk2(i+1,j)+wz3(i,j,k+1)-wz3(i,j,k))
            enddo
            wk3(im,j) = uc(im,j,k)  +  grid%dtdxe(j)/(wk(im,j) + wk(1,j))       &
                        * (wk2(im,j)-wk2(1,j)+wz3(im,j,k+1)-wz3(im,j,k))
         enddo

         do j=js2g0,jn2g0                  ! Assumes wk2 ghosted on N
            do i=1,im
               wk1(i,j) = vc(i,j,k)  +  grid%dtdy/(wk(i,j)+wk(i,j+1)) *          &
                          (wk2(i,j)-wk2(i,j+1)+wz(i,j,k+1)-wz(i,j,k))
            enddo
         enddo

         call pft2d( wk3(1,js2g0), grid%se,        &
                     grid%de, im, jlast-js2g0+1,       &
                     wk, wk2 )
         call pft2d( wk1(1,js2g0), grid%sc,        &
                     grid%dc, im, jn2g0-js2g0+1,       &
                     wk, wk2 )

         do j=js2g0,jn2g0
            do i=1,im
               v(i,j,k) = v(i,j,k) + wk1(i,j)
               u(i,j,k) = u(i,j,k) + wk3(i,j)
            enddo
         enddo
 
         if ( jlast == jm ) then
            do i=1,im
               u(i,jlast,k) = u(i,jlast,k) + wk3(i,jlast)
            enddo
         endif

6000  continue
#if (!defined USE_OMP) 
!CSD$ END PARALLEL DO
#endif
      call FVstopclock(grid,'---D_PGRAD_2')

#if defined( SPMD )
      if ( ipe /= 1 ) then
         call FVstartclock(grid,'---PRE_D_PGRAD_COMM_2')
         call mp_recv4d_ns( grid%commyz, im, jm, km, 1, jfirst, jlast,   &
                            kfirst, klast, ng_d, ng_d, delpf )
         call FVstopclock(grid,'---PRE_D_PGRAD_COMM_2')
      endif
#endif

    end if  !  (iam .lt. npes_yz)

      return
!EOC
      end subroutine cd_core
!-----------------------------------------------------------------------
