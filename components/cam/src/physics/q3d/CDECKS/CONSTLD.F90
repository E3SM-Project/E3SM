MODULE constld
!  Contains CRM parameters.

   USE shr_kind_mod, only: dbl_kind => shr_kind_r8
   USE parmsld,      only: nk2,nk3
   
IMPLICIT NONE

    REAL (KIND=dbl_kind), PARAMETER :: D0_0  = 0.0_dbl_kind
    REAL (KIND=dbl_kind), PARAMETER :: D1_0  = 1.0_dbl_kind
    REAL (KIND=dbl_kind), PARAMETER :: D2_0  = 2.0_dbl_kind
    REAL (KIND=dbl_kind), PARAMETER :: D3_0  = 3.0_dbl_kind
    REAL (KIND=dbl_kind), PARAMETER :: D4_0  = 4.0_dbl_kind
    REAL (KIND=dbl_kind), PARAMETER :: D6_0  = 6.0_dbl_kind
    REAL (KIND=dbl_kind), PARAMETER :: D8_0  = 8.0_dbl_kind
    
    REAL (KIND=dbl_kind), PARAMETER :: D0_5  = 0.5_dbl_kind
    REAL (KIND=dbl_kind), PARAMETER :: D1_5  = 1.5_dbl_kind 
    REAL (KIND=dbl_kind), PARAMETER :: D0_25 = 0.25_dbl_kind
    
    ! Logical Parameters (some will be removed later)
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
    
    ! Integer Parameters (some will be removed later)
    integer :: itt_max                ! # of CRM time steps per one GCM time step (DT_GCM/DT_CRM) 
    integer :: nrestart               ! restart writing interval
    integer :: ittadd                 ! time addition to restart time
    integer :: nxs                    ! accumulation frequency for averaging
    integer :: nxsavg                 ! averaging period for output
    integer :: nflprt                 ! frequency of writing the full output
    integer :: nout                   ! # of data outputs
    integer :: itinit                 ! starting time for random perturbation 
    integer :: itstop                 ! ending time for random perturbation
    integer :: nsflux                 ! frequency of surface flux calculation
    integer :: nrad                   ! frequency of radiation calculation
    integer :: niterw                 ! # of iterations for w
    integer :: niterxy                ! # of iterations for psi and chi
    integer :: klow_gwd               ! index of the lowest layer where GWD is applied

    ! Real Parameters (some will be removed later)
    real (kind=dbl_kind) :: a         ! coefficient for time scheme
    real (kind=dbl_kind) :: b         ! coefficient for time scheme
    real (kind=dbl_kind) :: dt        ! integration timestep (s)
    real (kind=dbl_kind) :: dt_gcm    ! dt of GCM 
    real (kind=dbl_kind) :: dz        ! grid size in z-(vertical) direction (m)
    real (kind=dbl_kind) :: dzsq      ! dz**2   (m^2)
    real (kind=dbl_kind) :: dx        ! grid size in x-direction (m)
    real (kind=dbl_kind) :: dy        ! grid size in y-direction (m)    
    real (kind=dbl_kind) :: dxsq      ! dx**2   (m^2)
    real (kind=dbl_kind) :: dysq      ! dy**2   (m^2)
    real (kind=dbl_kind) :: dxdy      ! dx x dy (m^2)    
    real (kind=dbl_kind) :: dx_gcm    ! gcm grid size in x-direction (m)
    real (kind=dbl_kind) :: dy_gcm    ! gcm grid size in y-direction (m)
    real (kind=dbl_kind) :: dz1       ! the lowest layer depth in z (m)
    real (kind=dbl_kind) :: domain    ! vertical domain for z'-selection
    real (kind=dbl_kind) :: zb        ! surface height (m)
    real (kind=dbl_kind) :: cz1       ! coef. for vertical coordinate
    real (kind=dbl_kind) :: cz2       ! coef. for vertical coordinate
    real (kind=dbl_kind) :: pi        ! 3.14159...        
    real (kind=dbl_kind) :: grav      ! gravitational constant (m/s^2)
    real (kind=dbl_kind) :: hlf       ! latent heat of vaporization (J/kg)
    real (kind=dbl_kind) :: hlm       ! latent heat of fusion (J/kg)
    real (kind=dbl_kind) :: cp        ! specific heat of dry air at constant pressure (J/kg/K)
    real (kind=dbl_kind) :: gas_const ! gas constant
    real (kind=dbl_kind) :: delta     ! 0.608
    real (kind=dbl_kind) :: vk        ! von Karman's constant
    real (kind=dbl_kind) :: zrsea     ! surface roughness over sea (m)
    real (kind=dbl_kind) :: zrland    ! surface roughness over land (m)
    real (kind=dbl_kind) :: sst       ! sea surface temperature (K)
    real (kind=dbl_kind) :: psfc      ! surface pressure (Pa)
    real (kind=dbl_kind) :: pisfc     ! Exner function at surface
    real (kind=dbl_kind) :: rearth    ! radius of the earth (m)
    real (kind=dbl_kind) :: omega     ! earth rotation rate (s^-1)    
    real (kind=dbl_kind) :: rad2deg   ! deg. of 1 rad (i.e., 180/pi)
    real (kind=dbl_kind) :: aladv     ! alpha in advection
    real (kind=dbl_kind) :: wrxmu     ! (mu/dt) in relaxed method 
    real (kind=dbl_kind) :: dtrad     ! time interval between radiation calls (s)       
    real (kind=dbl_kind) :: crad      ! coefficient for gravity wave damping
    real (kind=dbl_kind) :: gwdht     ! apply GWD above this height (m)
    real (kind=dbl_kind) :: dthmax    ! coefficient for random perturbation (K)
    real (kind=dbl_kind) :: z1pert    ! z_low for random perturbation (m)
    real (kind=dbl_kind) :: z2pert    ! z_high for random perturbation (m)

    ! Vertical Profiles  
    real (kind=dbl_kind) :: ZZ(nk3)    ! height at the level position (m)
    real (kind=dbl_kind) :: ZT(nk3)    ! height at the layer position (m)
    real (kind=dbl_kind) :: FNZ(nk3)   ! vertical map factor at the level position
    real (kind=dbl_kind) :: FNT(nk3)   ! vertical map factor at the layer position
    real (kind=dbl_kind) :: RHO(nk3)   ! air density at the layer position (kg/m**3)
    real (kind=dbl_kind) :: THBAR(nk3) ! mean profile of potential temp. (K)
    real (kind=dbl_kind) :: QVBAR(nk3) ! mean profile of water vapor mixing ratio (kg/kg)
    real (kind=dbl_kind) :: PBAR(nk3)  ! pressure at the layer position (Pa)
    real (kind=dbl_kind) :: PIBAR(nk3) ! Exner function at the layer position, T=TH3D*PIBAR
    
    real (kind=dbl_kind) :: FN1(nk2)   ! vertical map factor ratio (FNZ(K)/FNT(K+1))
    real (kind=dbl_kind) :: FN2(nk2)   ! vertical map factor ratio (FNZ(K)/FNT(K))
    real (kind=dbl_kind) :: RHOZ(nk2)  ! air density at the level position (kg/m**3)

END MODULE constld