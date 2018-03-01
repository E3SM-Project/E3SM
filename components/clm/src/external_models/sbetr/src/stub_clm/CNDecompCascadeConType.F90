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

     !-- properties of each decomposing pool
     logical           , pointer  :: floating_cn_ratio_decomp_pools(:) => null() ! TRUE => pool has fixed C:N ratio
     logical           , pointer  :: floating_cp_ratio_decomp_pools(:) => null() ! TRUE => pool has fixed C:N ratio
     character(len=8)  , pointer  :: decomp_pool_name_restart(:)    => null()   ! name of pool for restart files
     character(len=8)  , pointer  :: decomp_pool_name_history(:)   => null()    ! name of pool for history files
     character(len=20) , pointer  :: decomp_pool_name_long(:)    => null()      ! name of pool for netcdf long names
     character(len=8)  , pointer  :: decomp_pool_name_short(:)   => null()      ! name of pool for netcdf short names
     logical           , pointer  :: is_litter(:)          => null()            ! TRUE => pool is a litter pool
     logical           , pointer  :: is_soil(:)          => null()              ! TRUE => pool is a soil pool
     logical           , pointer  :: is_cwd(:)             => null()            ! TRUE => pool is a cwd pool
     real(r8)          , pointer  :: initial_cn_ratio(:)   => null()            ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_cp_ratio(:)  => null()             ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_stock(:)    => null()              ! initial concentration for seeding at spinup
     logical           , pointer  :: is_metabolic(:)     => null()              ! TRUE => pool is metabolic material
     logical           , pointer  :: is_cellulose(:)     => null()              ! TRUE => pool is cellulose
     logical           , pointer  :: is_lignin(:)       => null()               ! TRUE => pool is lignin
     real(r8)          , pointer  :: spinup_factor(:)    => null()              ! factor by which to scale AD and relevant processes by
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

    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%floating_cp_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_litter(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_soil(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cwd(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cp_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_stock(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_metabolic(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cellulose(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_lignin(0:ndecomp_pools))
    allocate(decomp_cascade_con%spinup_factor(0:ndecomp_pools))

    !-- properties of each pathway along decomposition cascade

    !-- first initialization of properties of each decomposing pool
    decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%floating_cp_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools)          = ''
    decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools)         = ''
    decomp_cascade_con%is_litter(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%is_soil(0:ndecomp_pools)                        = .false.
    decomp_cascade_con%is_cwd(0:ndecomp_pools)                         = .false.
    decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools)               = nan
    decomp_cascade_con%initial_cp_ratio(0:ndecomp_pools)               = nan
    decomp_cascade_con%initial_stock(0:ndecomp_pools)                  = nan
    decomp_cascade_con%is_metabolic(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_cellulose(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_lignin(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%spinup_factor(0:ndecomp_pools)                  = nan

  end subroutine init_decomp_cascade_constants

end module CNDecompCascadeConType
