module CNDecompCascadeConType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Decomposition Cascade Type
  !
  ! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : ndecomp_cascade_transitions, ndecomp_pools
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_decomp_cascade_constants
  !
  type, public :: decomp_cascade_type
     !-- properties of each pathway along decomposition cascade 
     character(len=8)  , pointer :: cascade_step_name(:)               ! name of transition
     integer           , pointer :: cascade_donor_pool(:)              ! which pool is C taken from for a given decomposition step
     integer           , pointer :: cascade_receiver_pool(:)           ! which pool is C added to for a given decomposition step

     !-- properties of each decomposing pool
     logical           , pointer  :: floating_cn_ratio_decomp_pools(:) ! TRUE => pool has fixed C:N ratio
     character(len=8)  , pointer  :: decomp_pool_name_restart(:)       ! name of pool for restart files
     character(len=8)  , pointer  :: decomp_pool_name_history(:)       ! name of pool for history files
     character(len=20) , pointer  :: decomp_pool_name_long(:)          ! name of pool for netcdf long names
     character(len=8)  , pointer  :: decomp_pool_name_short(:)         ! name of pool for netcdf short names
     logical           , pointer  :: is_litter(:)                      ! TRUE => pool is a litter pool
     logical           , pointer  :: is_soil(:)                        ! TRUE => pool is a soil pool
     logical           , pointer  :: is_cwd(:)                         ! TRUE => pool is a cwd pool
     real(r8)          , pointer  :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_stock(:)                  ! initial concentration for seeding at spinup
     logical           , pointer  :: is_metabolic(:)                   ! TRUE => pool is metabolic material
     logical           , pointer  :: is_cellulose(:)                   ! TRUE => pool is cellulose
     logical           , pointer  :: is_lignin(:)                      ! TRUE => pool is lignin
     real(r8)          , pointer  :: spinup_factor(:)                  ! factor by which to scale AD and relevant processes by
  end type decomp_cascade_type

  type(decomp_cascade_type), public :: decomp_cascade_con
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_decomp_cascade_constants()
    !
    ! !DESCRIPTION:
    ! Initialize decomposition cascade state
    !------------------------------------------------------------------------

    !-- properties of each pathway along decomposition cascade 
    allocate(decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions))

    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_litter(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_soil(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cwd(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_stock(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_metabolic(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cellulose(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_lignin(0:ndecomp_pools))
    allocate(decomp_cascade_con%spinup_factor(0:ndecomp_pools))

    !-- properties of each pathway along decomposition cascade 
    decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions) = ''
    decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions) = 0
    decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions) = 0

    !-- first initialization of properties of each decomposing pool
    decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools)          = ''
    decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools)         = ''
    decomp_cascade_con%is_litter(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%is_soil(0:ndecomp_pools)                        = .false.
    decomp_cascade_con%is_cwd(0:ndecomp_pools)                         = .false.
    decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools)               = nan
    decomp_cascade_con%initial_stock(0:ndecomp_pools)                  = nan
    decomp_cascade_con%is_metabolic(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_cellulose(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_lignin(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%spinup_factor(0:ndecomp_pools)                  = nan

  end subroutine init_decomp_cascade_constants

end module CNDecompCascadeConType
