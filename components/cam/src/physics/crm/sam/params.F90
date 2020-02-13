module params

  use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
  shr_const_latice, shr_const_latsub, shr_const_rgas, &
  shr_const_mwwv, shr_const_stebol, shr_const_tkfrz, &
  shr_const_mwdair, shr_const_g, shr_const_karman, &
  shr_const_rhofw

  implicit none

#ifdef CRM_SINGLE_PRECISION
  integer, parameter :: crm_rknd = selected_real_kind( 6) ! 4 byte real
#else
  ! default precision of real - kind(1.d0)
  integer, parameter :: crm_rknd = selected_real_kind(12) ! 8 byte real
#endif

  !   Constants:

  real(crm_rknd), parameter :: cp    = real( shr_const_cpdair ,crm_rknd)
  real(crm_rknd), parameter :: ggr   = real( shr_const_g      ,crm_rknd)
  real(crm_rknd), parameter :: lcond = real( shr_const_latvap ,crm_rknd)
  real(crm_rknd), parameter :: lfus  = real( shr_const_latice ,crm_rknd)
  real(crm_rknd), parameter :: lsub  = real( lcond + lfus     ,crm_rknd)
  real(crm_rknd), parameter :: rgas  = real( shr_const_rdair  ,crm_rknd)
  real(crm_rknd), parameter :: rv    = real( shr_const_rgas/shr_const_mwwv ,crm_rknd)

  real(crm_rknd), parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
  real(crm_rknd), parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
  real(crm_rknd), parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

  real(crm_rknd), parameter :: fac_cond = lcond/cp
  real(crm_rknd), parameter :: fac_fus  = lfus/cp
  real(crm_rknd), parameter :: fac_sub  = lsub/cp

  real(crm_rknd), parameter ::  pi = 3.141592653589793

  !
  ! internally set parameters:

  real(crm_rknd)   epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
  logical:: dosubsidence = .false.
  real(crm_rknd), allocatable :: fcorz(:)      ! Vertical Coriolis parameter

  !----------------------------------------------
  ! Parameters set by PARAMETERS namelist:
  ! Initialized to default values.
  !----------------------------------------------

  real(crm_rknd):: ug = 0.        ! Velocity of the Domain's drift in x direction
  real(crm_rknd):: vg = 0.        ! Velocity of the Domain's drift in y direction
  real(crm_rknd), allocatable :: fcor(:)  ! Coriolis parameter
  real(crm_rknd), allocatable :: longitude0(:)    ! latitude of the domain's center
  real(crm_rknd), allocatable :: latitude0 (:)    ! longitude of the domain's center

  real(crm_rknd), allocatable :: z0(:)            ! roughness length
  logical :: les =.false.    ! flag for Large-Eddy Simulation
  logical, allocatable :: ocean(:)           ! flag indicating that surface is water
  logical, allocatable :: land(:)            ! flag indicating that surface is land

  logical:: docloud       = .false.   ! 
  logical:: doprecip      = .false.   ! 
  logical:: dodamping     = .false.   ! Newtonian damping for upper levels
  logical:: dosgs         = .false.   ! sub-grid turbulence scheme
  logical:: dosurface     = .false.   ! surface scheme to calculate friction within CRM
  logical:: docoriolis    = .false.   ! not normally used for MMF
  logical:: dowallx       = .false.   ! not normally used for MMF
  logical:: dowally       = .false.   ! not normally used for MMF
  logical:: docolumn      = .false.   ! not normally used for MMF
  logical:: dotracers     = .false.   ! not normally used for MMF
  logical:: dosmoke       = .false.   ! not normally used for MMF

  integer, parameter :: asyncid = 1

  real(crm_rknd), allocatable :: uhl(:)      ! current large-scale velocity in x near sfc
  real(crm_rknd), allocatable :: vhl(:)      ! current large-scale velocity in y near sfc
  real(crm_rknd), allocatable :: taux0(:)    ! surface stress in x, m2/s2
  real(crm_rknd), allocatable :: tauy0(:)    ! surface stress in y, m2/s2


contains

  
  subroutine allocate_params(ncrms)
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    allocate(fcor (ncrms))
    allocate(fcorz(ncrms))
    allocate(longitude0(ncrms))
    allocate(latitude0 (ncrms))
    allocate(z0        (ncrms))
    allocate(ocean     (ncrms))
    allocate(land      (ncrms))
    allocate(uhl       (ncrms))
    allocate(vhl       (ncrms))
    allocate(taux0     (ncrms))
    allocate(tauy0     (ncrms))

    call prefetch(fcor )
    call prefetch(fcorz)
    call prefetch(longitude0)
    call prefetch(latitude0 )
    call prefetch(z0)
    call prefetch(ocean)
    call prefetch(land)
    call prefetch(uhl)
    call prefetch(vhl)
    call prefetch(taux0)
    call prefetch(tauy0)

    fcor  = 0
    fcorz = 0
    longitude0 = 0
    latitude0  = 0
    z0 = 0.035
    ocean = .false.
    land = .false.
    uhl = 0
    vhl = 0
    taux0 = 0
    tauy0 = 0
  end subroutine allocate_params

  
  subroutine deallocate_params()
    implicit none
    deallocate(fcor )
    deallocate(fcorz)
    deallocate(longitude0)
    deallocate(latitude0 )
    deallocate(z0)
    deallocate(ocean)
    deallocate(land)
    deallocate(uhl)
    deallocate(vhl)
    deallocate(taux0)
    deallocate(tauy0)
  end subroutine deallocate_params


end module params
