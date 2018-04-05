MODULE constld
!  Contains CRM parameters.

   USE shr_kind_mod, only: r8 => shr_kind_r8
   USE parmsld,      only: nk2,nk3
   
IMPLICIT NONE
    
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
    real (kind=r8) :: a         ! coefficient for time scheme
    real (kind=r8) :: b         ! coefficient for time scheme
    real (kind=r8) :: dt        ! integration timestep (s)
    real (kind=r8) :: dt_gcm    ! dt of GCM 
    real (kind=r8) :: dz        ! grid size in z-(vertical) direction (m)
    real (kind=r8) :: dzsq      ! dz**2   (m^2)
    real (kind=r8) :: dx        ! grid size in x-direction (m)
    real (kind=r8) :: dy        ! grid size in y-direction (m)    
    real (kind=r8) :: dxsq      ! dx**2   (m^2)
    real (kind=r8) :: dysq      ! dy**2   (m^2)
    real (kind=r8) :: dxdy      ! dx x dy (m^2)    
    real (kind=r8) :: dx_gcm    ! gcm grid size in x-direction (m)
    real (kind=r8) :: dy_gcm    ! gcm grid size in y-direction (m)
    real (kind=r8) :: dz1       ! the lowest layer depth in z (m)
    real (kind=r8) :: domain    ! vertical domain for stretched grids
    real (kind=r8) :: zb        ! surface height (m)
    real (kind=r8) :: cz1       ! coef. for vertical coordinate
    real (kind=r8) :: cz2       ! coef. for vertical coordinate
    real (kind=r8) :: pi        ! 3.14159...        
    real (kind=r8) :: grav      ! gravitational constant (m/s^2)
    real (kind=r8) :: hlf       ! latent heat of vaporization (J/kg)
    real (kind=r8) :: hlm       ! latent heat of fusion (J/kg)
    real (kind=r8) :: cp        ! specific heat of dry air at constant pressure (J/kg/K)
    real (kind=r8) :: gas_const ! gas constant
    real (kind=r8) :: delta     ! 0.608
    real (kind=r8) :: vk        ! von Karman constant
    real (kind=r8) :: zrsea     ! surface roughness over sea (m)
    real (kind=r8) :: zrland    ! surface roughness over land (m)
    real (kind=r8) :: sst       ! sea surface temperature (K)
    real (kind=r8) :: psfc      ! surface pressure (Pa)
    real (kind=r8) :: pisfc     ! Exner function at surface
    real (kind=r8) :: rearth    ! radius of the earth (m)
    real (kind=r8) :: omega     ! earth rotation rate (s^-1)    
    real (kind=r8) :: rad2deg   ! deg. of 1 rad (i.e., 180/pi)
    real (kind=r8) :: aladv     ! alpha in advection
    real (kind=r8) :: wrxmu     ! (mu/dt) in relaxed method 
    real (kind=r8) :: dtrad     ! time interval between radiation calls (s)       
    real (kind=r8) :: crad      ! coefficient for gravity wave damping
    real (kind=r8) :: gwdht     ! apply GWD above this height (m)
    real (kind=r8) :: dthmax    ! coefficient for random perturbation (K)
    real (kind=r8) :: z1pert    ! z_low for random perturbation (m)
    real (kind=r8) :: z2pert    ! z_high for random perturbation (m)

    ! Vertical Profiles  
    real (kind=r8) :: ZZ(nk3)    ! height at the level position (m)
    real (kind=r8) :: ZT(nk3)    ! height at the layer position (m)
    real (kind=r8) :: FNZ(nk3)   ! vertical map factor at the level position
    real (kind=r8) :: FNT(nk3)   ! vertical map factor at the layer position
    real (kind=r8) :: RHO(nk3)   ! air density at the layer position (kg/m**3)
    real (kind=r8) :: THBAR(nk3) ! mean profile of potential temp. (K)
    real (kind=r8) :: QVBAR(nk3) ! mean profile of water vapor mixing ratio (kg/kg)
    real (kind=r8) :: PBAR(nk3)  ! pressure at the layer position (Pa)
    real (kind=r8) :: PIBAR(nk3) ! Exner function at the layer position, T=TH3D*PIBAR
    
    real (kind=r8) :: FN1(nk2)   ! vertical map factor ratio (FNZ(K)/FNT(K+1))
    real (kind=r8) :: FN2(nk2)   ! vertical map factor ratio (FNZ(K)/FNT(K))
    real (kind=r8) :: RHOZ(nk2)  ! air density at the level position (kg/m**3)
    real (kind=r8) :: PBARZ(nk2) ! pressure at the level position (Pa)
    real (kind=r8) :: PIBARZ(nk2)! Exner function at the level position

END MODULE constld