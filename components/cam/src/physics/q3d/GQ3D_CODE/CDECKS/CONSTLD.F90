MODULE constld
!  Contains CRM parameters.

   USE shr_kind_mod, only: r8 => shr_kind_r8
   USE parmsld,      only: nk2,nk3,nk2_ext

   USE physconst, only: piv=>pi         ! = 3.14159265358979323846
   USE physconst, only: rearth          ! radius of the earth (m) = 6.37122e6
   USE physconst, only: omega           ! earth rotation rate (rad/s) = 7.292e-5
   USE physconst, only: epsilo
   USE physconst, only: cappa           ! R/cp
   USE physconst, only: vk=>karman      ! Von Karman constant = 0.4
   USE physconst, only: hlm=>latice
   USE physconst, only: hlf=>latvap     ! latent heat of evaporation (J/kg) = 2.501e6
   USE physconst, only: grav=>gravit    ! acceleration of gravity (m/s^2) = 9.80616
   USE physconst, only: cp=>cpair       ! specific heat of dry air (J/kg/K) = 1.00464e3
   USE physconst, only: gas_const=>rair ! Dry air gas constant (J/K/kg) = 287.042
   USE physconst, only: wvg_const=>rh2o ! Water vapor gas constant (J/K/kg) = 461.5
   USE physconst, only: delta=>zvir     ! = (Rv/Ra) - 1, TV = (1+delta*qv)T = 0.608
   USE physconst, only: h2otrip         ! Triple point temperature of water (K) = 287.16
   USE physconst, only: t_melt=>tmelt   ! Freezing point of water (K) = 273.15
   USE physconst, only: cp_h2o=>cpliq   ! specific heat of liquid water (J/K/kg) = 4.188e3
   USE physconst, only: cp_ice=>cpice   ! specific heat of ice (J/K/kg) = 2.11727e3
   USE physconst, only: cpwv            ! specific heat of water vapor (J/K/kg) = 1.810e3
   USE physconst, only: sec_cday=>cday  ! sec in calendar day (sec) = 86400.0
   USE physconst, only: sec_sday=>sday  ! sec in siderial day (sec) = 86164.0
   USE physconst, only: pstd            ! Standard pressure (Pascals) = 101325.0
   USE physconst, only: rho_h2o=>rhoh2o ! Density of liquid water (STP) = 1.0e3

IMPLICIT NONE

!-------------------------------------------------------------------------------------------
!   Constant Parameters
!-------------------------------------------------------------------------------------------
    real (kind=r8), parameter :: val_missing = -1.0E30

    real (kind=r8), parameter :: pzero = 100000.0_r8   ! Reference pressure (Pascals)

    real (kind=r8) :: hls = hlf+hlm                    ! latent heat of sublimation (J/kg)
                                                       ! = 2.835e6

    real (kind=r8) :: rad2deg = 180.0_r8/piv           ! deg. of 1 rad

!-------------------------------------------------------------------------------------------
!   Values of the below are specified in "Subroutine INIT_common" (q3d_init_module.F90)
!-------------------------------------------------------------------------------------------
    ! Logical variables (some will be removed later)
    logical :: nosfx                  ! T, no surface flux
    logical :: notopo                 ! T, no topography
    logical :: nomap                  ! T, no (horizontal) mapping
    logical :: zconst                 ! T, use of constant-depth vertical grids
    logical :: physics                ! T, include physical processes
    logical :: turbulence             ! T, include turbulence
    logical :: radcode                ! T, include radiation (RRTMG)
    logical :: microcode              ! T, include microphysics
    logical :: buoy                   ! T, include buoyancy
    logical :: notherm                ! T, no thermodynamical processes
    logical :: fconst                 ! T, no Coriolis effect (f=0)
    
    logical :: localp                 ! T, use local p (& pi) in physics

    ! Integer variables (some will be removed later)
    integer :: itt_max                ! # of CRM time steps per one GCM time step (DT_GCM/DT_CRM)
    integer :: itinit                 ! starting time for random perturbation
    integer :: itstop                 ! ending time for random perturbation
    integer :: nsflux                 ! frequency of surface flux calculation
    integer :: nrad                   ! frequency of radiation calculation
    integer :: niterw                 ! # of iterations for w
    integer :: niterxy                ! # of iterations for psi and chi
    integer :: niterxy_init           ! # of iterations for psi and chi (initialization)
    integer :: klow_gwd               ! index of the lowest layer where GWD is applied
    integer :: iver_wind              ! Version number of Subroutine wind_3d (0, 1, 2) 

    ! Real variables (some will be removed later)
    real (kind=r8) :: a         ! coefficient for time scheme
    real (kind=r8) :: b         ! coefficient for time scheme
    real (kind=r8) :: dt        ! integration timestep (s)
    real (kind=r8) :: dt_gcm    ! dt of GCM
    real (kind=r8) :: dz        ! grid size in z-(vertical) direction (m)
    real (kind=r8) :: dzsq      ! dz**2   (m^2)
    real (kind=r8) :: dz_ext    ! grid size of the extended vertical layers (m)
    real (kind=r8) :: dx        ! grid size in x-direction (m)
    real (kind=r8) :: dy        ! grid size in y-direction (m)
    real (kind=r8) :: dxsq      ! dx**2   (m^2)
    real (kind=r8) :: dysq      ! dy**2   (m^2)
    real (kind=r8) :: dxdy      ! dx x dy (m^2)
    real (kind=r8) :: dx_gcm    ! gcm grid size in x-direction (m)
    real (kind=r8) :: dy_gcm    ! gcm grid size in y-direction (m)
    real (kind=r8) :: zrsea     ! surface roughness over sea (m)
    real (kind=r8) :: zrland    ! surface roughness over land (m)
    real (kind=r8) :: sst       ! sea surface temperature (K)
    real (kind=r8) :: psfc      ! surface pressure (Pa)
    real (kind=r8) :: pisfc     ! Exner function at surface
    real (kind=r8) :: aladv     ! alpha-value in advection
    real (kind=r8) :: rxtau_q   ! relaxation timescale for thermodynamic variables (s)
    real (kind=r8) :: rxtau_d   ! relaxation timescale for dynamical variables (s)
    real (kind=r8) :: wrxmu     ! (mu/dt) in relaxed method
    real (kind=r8) :: dtrad     ! time interval between radiation calls (s)
    real (kind=r8) :: crad      ! coefficient for gravity wave damping
    real (kind=r8) :: gwdht     ! apply GWD above this height (m)
    real (kind=r8) :: dthmax    ! coefficient for random perturbation (K)
    real (kind=r8) :: z1pert    ! z_low for random perturbation (m)
    real (kind=r8) :: z2pert    ! z_high for random perturbation (m)

    ! Vertical Profiles
    ! Note that any variable used as a hist_coord (via add_hist_coord) must
    !   be a target or a pointer.
    real (kind=r8), target :: ZZ(nk2_ext)     ! height at the level position (m)
    real (kind=r8), target :: ZT(nk2_ext)     ! height at the layer position (m)
    real (kind=r8)         :: FNZ(nk2_ext)    ! vertical map factor at the level position
    real (kind=r8)         :: FNT(nk2_ext)    ! vertical map factor at the layer position

    real (kind=r8)         :: THBAR(nk2_ext)  ! mean profile of potential temp. (K)

    real (kind=r8), target :: RHO(nk2_ext)    ! air density at the layer position (kg/m**3)
    real (kind=r8), target :: PBAR(nk2_ext)   ! pressure at the layer position (Pa)
    real (kind=r8)         :: PIBAR(nk2_ext)  ! Exner function at the layer position, T=TH3D*PIBAR

    real (kind=r8), target :: RHOZ(nk2_ext)   ! air density at the level position (kg/m**3)
    real (kind=r8), target :: PBARZ(nk2_ext)  ! pressure at the level position (Pa)
    real (kind=r8)         :: PIBARZ(nk2_ext) ! Exner function at the level position

    real (kind=r8)         :: FN1(nk2)        ! vertical map factor ratio (FNZ(K)/FNT(K+1))
    real (kind=r8)         :: FN2(nk2)        ! vertical map factor ratio (FNZ(K)/FNT(K))

END MODULE constld
