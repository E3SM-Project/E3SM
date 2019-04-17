module PhenologyFLuxLimitMod

!DESCRIPTION
! limit the allocation fluxes resulting from pheonology
! calculation to avoid negative state varaibles for
! carbon, nitrogen and phosphours.
! The current implementation does not use sparse
! matrix, future version should consider this possibility.

  use LSparseMatMod               , only : sparseMat_type, flux_correction,spm_list_type
  use shr_kind_mod                , only : r8 => shr_kind_r8
  use VegetationType              , only : veg_pp
  use VegetationPropertiesType    , only : veg_vp
  use clm_time_manager            , only : get_step_size
  use pftvarcon                   , only : npcropmin
  use clm_varctl                  , only : iulog
  use abortutils                  , only : endrun
implicit none
  private
#define mov(a,b) a=b
#define imov(a,b) b=max(a,0._r8)
#define ascal(a,b) a=a*b
  integer :: s_cpool
  integer :: s_leafc
  integer :: s_leafc_xfer
  integer :: s_leafc_storage
  integer :: s_frootc
  integer :: s_frootc_xfer
  integer :: s_frootc_sotrage
  integer :: s_livestemc
  integer :: s_livestemc_xfer
  integer :: s_livestemc_storage
  integer :: s_livecrootc
  integer :: s_livecrootc_xfer
  integer :: s_livecrootc_storage
  integer :: s_deadstemc
  integer :: s_deadstemc_xfer
  integer :: s_deadstemc_storage
  integer :: s_deadcrootc
  integer :: s_deadcrootc_xfer
  integer :: s_deadcrootc_storage
  integer :: s_grainc
  integer :: s_grainc_xfer
  integer :: s_grainc_storage
  integer :: s_gresp_xfer
  integer :: s_gresp_storage

  integer :: f_cpool_to_leafc
  integer :: f_cpool_to_leafc_storage
  integer :: f_cpool_to_frootc
  integer :: f_cpool_to_frootc_storage
  integer :: f_cpool_to_xsmrpool
  integer :: f_cpool_to_gresp_storage
  integer :: f_cpool_to_ar
  integer :: f_cpool_to_livestemc
  integer :: f_cpool_to_livestemc_storage
  integer :: f_cpool_to_deadstemc
  integer :: f_cpool_to_deadstemc_storage
  integer :: f_cpool_to_livecrootc
  integer :: f_cpool_to_livecrootc_storage
  integer :: f_cpool_to_deadcrootc
  integer :: f_cpool_to_deadcrootc_storage
  integer :: f_cpool_to_grainc
  integer :: f_cpool_to_grainc_storage
  integer :: f_leafc_to_litter
  integer :: f_leafc_xfer_to_leafc
  integer :: f_leafc_storage_to_xfer
  integer :: f_frootc_to_litter
  integer :: f_frootc_xfer_to_frootc
  integer :: f_frootc_storage_to_xfer
  integer :: f_livestemc_to_deadstemc
  integer :: f_livestemc_xfer_to_livestemc
  integer :: f_livestemc_storage_to_xfer
  integer :: f_livestemc_to_litter
  integer :: f_deadstemc_xfer_to_deadstemc
  integer :: f_deadstemc_storage_to_xfer
  integer :: f_livecrootc_xfer_to_livecrootc
  integer :: f_livecrootc_storage_to_xfer
  integer :: f_livecrootc_to_deadcrootc
  integer :: f_deadcrootc_xfer_to_deadcrootc
  integer :: f_deadcrootc_storage_to_xfer
  integer :: f_grainc_to_food
  integer :: f_grainc_xfer_to_grainc
  integer :: f_grainc_storage_to_xfer
  integer :: f_gresp_storage_to_xfer


  integer :: s_npool
  integer :: s_leafn
  integer :: s_leafn_xfer
  integer :: s_leafn_storage
  integer :: s_frootn
  integer :: s_frootn_xfer
  integer :: s_frootn_sotrage
  integer :: s_livestemn
  integer :: s_livestemn_xfer
  integer :: s_livestemn_storage
  integer :: s_deadstemn
  integer :: s_deadstemn_xfer
  integer :: s_deadstemn_storage
  integer :: s_livecrootn
  integer :: s_livecrootn_xfer
  integer :: s_livecrootn_storage
  integer :: s_deadcrootn
  integer :: s_deadcrootn_xfer
  integer :: s_deadcrootn_storage
  integer :: s_retransn
  integer :: s_grainn
  integer :: s_grainn_xfer
  integer :: s_grainn_storage

  integer :: f_npool_to_leafn
  integer :: f_npool_to_leafn_storage
  integer :: f_npool_to_frootn
  integer :: f_npool_to_frootn_storage
  integer :: f_npool_to_livestemn_storage
  integer :: f_npool_to_liverootn_storage
  integer :: f_npool_to_livestemn
  integer :: f_npool_to_livecrootn
  integer :: f_npool_to_livecrootn_storage
  integer :: f_npool_to_deadstemn
  integer :: f_npool_to_deadcrootn
  integer :: f_npool_to_deadstemn_storage
  integer :: f_npool_to_deadcrootn_storage
  integer :: f_npool_to_grainn
  integer :: f_npool_to_grainn_storage
  integer :: f_leafn_to_retransn
  integer :: f_leafn_to_litter
  integer :: f_leafn_xfer_to_leafn
  integer :: f_leafn_storage_to_xfer
  integer :: f_frootn_xfer_to_frootn
  integer :: f_frootn_storage_to_xfer
  integer :: f_frootn_to_retransn
  integer :: f_frootn_to_litter
  integer :: f_livestemn_to_litter
  integer :: f_livestemn_storage_to_xfer
  integer :: f_livestemn_xfer_to_livestemn
  integer :: f_livestemn_to_retransn
  integer :: f_livestemn_to_deadstemn
  integer :: f_livecrootn_storage_to_xfer
  integer :: f_livecrootn_xfer_to_livecrootn
  integer :: f_livecrootn_to_deadcrootn
  integer :: f_livecrootn_to_retransn
  integer :: f_deadstemn_storage_to_xfer
  integer :: f_deadstemn_xfer_to_deadstem
  integer :: f_deadcrootn_storage_to_xfer
  integer :: f_deadcrootn_xfer_to_deadcrootn
  integer :: f_grainn_xfer_to_grainn
  integer :: f_grainn_to_food
  integer :: f_retransn_to_npool
  integer :: f_supplement_to_plantn
  integer :: n_carbon_fluxes
  integer :: n_nutrient_fluxes
  integer :: n_carbon_states
  integer :: n_nutrient_states

  class(sparseMat_type), pointer :: spm_carbon_p,spm_carbon_d
  class(sparseMat_type), pointer :: spm_nutrient_p,spm_nutrient_d
  type(spm_list_type), pointer :: spm_list
  public :: InitPhenoFluxLimiter
  public :: phenology_flux_limiter
contains


  function ic_next(id)result(ans)

  implicit none
  integer, intent(inout) :: id
  integer :: ans
  id=id+1
  ans=id
  end function ic_next
  !-------------------------------------------------------------------------------
  subroutine init_phenofluxl_counters()

  implicit none
  integer :: id

  id = 0
  s_cpool              = ic_next(id)
  s_leafc              = ic_next(id)
  s_leafc_xfer         = ic_next(id)
  s_leafc_storage      = ic_next(id)
  s_frootc             = ic_next(id)
  s_frootc_xfer        = ic_next(id)
  s_frootc_sotrage     = ic_next(id)
  s_livestemc          = ic_next(id)
  s_livestemc_xfer     = ic_next(id)
  s_livestemc_storage  = ic_next(id)
  s_livecrootc         = ic_next(id)
  s_livecrootc_xfer    = ic_next(id)
  s_livecrootc_storage = ic_next(id)
  s_deadstemc          = ic_next(id)
  s_deadstemc_xfer     = ic_next(id)
  s_deadstemc_storage  = ic_next(id)
  s_deadcrootc         = ic_next(id)
  s_deadcrootc_xfer    = ic_next(id)
  s_deadcrootc_storage = ic_next(id)
  s_grainc             = ic_next(id)
  s_grainc_xfer        = ic_next(id)
  s_grainc_storage     = ic_next(id)
  s_gresp_xfer         = ic_next(id)
  s_gresp_storage      = ic_next(id)
  n_carbon_states = id
  id=0
  f_cpool_to_leafc               = ic_next(id)
  f_cpool_to_leafc_storage       = ic_next(id)
  f_cpool_to_frootc              = ic_next(id)
  f_cpool_to_frootc_storage      = ic_next(id)
  f_cpool_to_xsmrpool            = ic_next(id)
  f_cpool_to_gresp_storage       = ic_next(id)
  f_cpool_to_ar                  = ic_next(id)
  f_cpool_to_livestemc           = ic_next(id)
  f_cpool_to_livestemc_storage   = ic_next(id)
  f_cpool_to_deadstemc           = ic_next(id)
  f_cpool_to_deadstemc_storage   = ic_next(id)
  f_cpool_to_livecrootc          = ic_next(id)
  f_cpool_to_livecrootc_storage  = ic_next(id)
  f_cpool_to_deadcrootc          = ic_next(id)
  f_cpool_to_deadcrootc_storage  = ic_next(id)
  f_cpool_to_grainc              = ic_next(id)
  f_cpool_to_grainc_storage      = ic_next(id)
  f_leafc_to_litter              = ic_next(id)
  f_leafc_xfer_to_leafc          = ic_next(id)
  f_leafc_storage_to_xfer        = ic_next(id)
  f_frootc_to_litter             = ic_next(id)
  f_frootc_xfer_to_frootc        = ic_next(id)
  f_frootc_storage_to_xfer       = ic_next(id)
  f_livestemc_to_deadstemc       = ic_next(id)
  f_livestemc_xfer_to_livestemc  = ic_next(id)
  f_livestemc_storage_to_xfer    = ic_next(id)
  f_livestemc_to_litter          = ic_next(id)
  f_deadstemc_xfer_to_deadstemc  = ic_next(id)
  f_deadstemc_storage_to_xfer    = ic_next(id)
  f_livecrootc_xfer_to_livecrootc= ic_next(id)
  f_livecrootc_storage_to_xfer   = ic_next(id)
  f_livecrootc_to_deadcrootc     = ic_next(id)
  f_deadcrootc_xfer_to_deadcrootc= ic_next(id)
  f_deadcrootc_storage_to_xfer   = ic_next(id)
  f_grainc_to_food               = ic_next(id)
  f_grainc_xfer_to_grainc        = ic_next(id)
  f_grainc_storage_to_xfer       = ic_next(id)
  f_gresp_storage_to_xfer        = ic_next(id)
  n_carbon_fluxes = id
  id=0
  s_npool               = ic_next(id)
  s_leafn               = ic_next(id)
  s_leafn_xfer          = ic_next(id)
  s_leafn_storage       = ic_next(id)
  s_frootn              = ic_next(id)
  s_frootn_xfer         = ic_next(id)
  s_frootn_sotrage      = ic_next(id)
  s_livestemn           = ic_next(id)
  s_livestemn_xfer      = ic_next(id)
  s_livestemn_storage   = ic_next(id)
  s_deadstemn           = ic_next(id)
  s_deadstemn_xfer      = ic_next(id)
  s_deadstemn_storage   = ic_next(id)
  s_livecrootn          = ic_next(id)
  s_livecrootn_xfer     = ic_next(id)
  s_livecrootn_storage  = ic_next(id)
  s_deadcrootn          = ic_next(id)
  s_deadcrootn_xfer     = ic_next(id)
  s_deadcrootn_storage  = ic_next(id)
  s_grainn              = ic_next(id)
  s_grainn_xfer         = ic_next(id)
  s_grainn_storage      = ic_next(id)
  s_retransn            = ic_next(id)
  n_nutrient_states = id
  id = 0
  f_npool_to_leafn                = ic_next(id)
  f_npool_to_leafn_storage        = ic_next(id)
  f_npool_to_frootn               = ic_next(id)
  f_npool_to_frootn_storage       = ic_next(id)
  f_npool_to_livestemn            = ic_next(id)
  f_npool_to_livestemn_storage    = ic_next(id)
  f_npool_to_livecrootn           = ic_next(id)
  f_npool_to_livecrootn_storage   = ic_next(id)
  f_npool_to_deadstemn            = ic_next(id)
  f_npool_to_deadcrootn           = ic_next(id)
  f_npool_to_deadstemn_storage    = ic_next(id)
  f_npool_to_deadcrootn_storage   = ic_next(id)
  f_npool_to_grainn               = ic_next(id)
  f_npool_to_grainn_storage       = ic_next(id)
  f_leafn_to_retransn             = ic_next(id)
  f_leafn_to_litter               = ic_next(id)
  f_leafn_xfer_to_leafn           = ic_next(id)
  f_leafn_storage_to_xfer         = ic_next(id)
  f_frootn_xfer_to_frootn         = ic_next(id)
  f_frootn_storage_to_xfer        = ic_next(id)
  f_frootn_to_retransn            = ic_next(id)
  f_frootn_to_litter              = ic_next(id)
  f_livestemn_to_litter           = ic_next(id)
  f_livestemn_storage_to_xfer     = ic_next(id)
  f_livestemn_xfer_to_livestemn   = ic_next(id)
  f_livestemn_to_retransn         = ic_next(id)
  f_livestemn_to_deadstemn        = ic_next(id)
  f_livecrootn_to_deadcrootn      = ic_next(id)
  f_livecrootn_to_retransn        = ic_next(id)
  f_livecrootn_storage_to_xfer    = ic_next(id)
  f_livecrootn_xfer_to_livecrootn = ic_next(id)
  f_deadstemn_storage_to_xfer     = ic_next(id)
  f_deadstemn_xfer_to_deadstem    = ic_next(id)
  f_deadcrootn_storage_to_xfer    = ic_next(id)
  f_deadcrootn_xfer_to_deadcrootn = ic_next(id)
  f_grainn_xfer_to_grainn         = ic_next(id)
  f_grainn_to_food                = ic_next(id)
  f_retransn_to_npool             = ic_next(id)
  f_supplement_to_plantn          = ic_next(id)
  n_nutrient_fluxes = id
  end subroutine init_phenofluxl_counters
  !-------------------------------------------------------------------------------
  subroutine InitPhenoFluxLimiter()

  use LSparseMatMod, only : spm_list_type,  spm_list_to_mat
  use LSparseMatMod, only : spm_list_init, spm_list_insert
  implicit none
  integer :: nelms
  type(spm_list_type), pointer :: spm_list

  call init_phenofluxl_counters

  call spm_list_init(spm_list, -1._r8, f_cpool_to_leafc                , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_leafc_storage      , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_frootc              , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_frootc_storage      , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_xsmrpool            , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_gresp_storage       , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_ar                  , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_livestemc           , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_livestemc_storage   , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_deadstemc           , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_deadstemc_storage   , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_livecrootc          , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_livecrootc_storage  , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_deadcrootc          , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_deadcrootc_storage  , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_grainc              , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_cpool_to_grainc_storage      , s_cpool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_leafc_to_litter              , s_leafc, nelms)
  call spm_list_insert(spm_list, -1._r8, f_leafc_xfer_to_leafc          , s_leafc_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_leafc_storage_to_xfer        , s_leafc_storage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_frootc_to_litter             , s_frootc, nelms)
  call spm_list_insert(spm_list, -1._r8, f_frootc_xfer_to_frootc        , s_frootc_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_frootc_storage_to_xfer       , s_frootc_sotrage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemc_to_deadstemc       , s_livestemc,nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemc_xfer_to_livestemc  , s_livestemc_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemc_storage_to_xfer    , s_livestemc_storage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemc_to_litter          , s_livestemc, nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadstemc_xfer_to_deadstemc  , s_deadstemc_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadstemc_storage_to_xfer    , s_deadstemc_storage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livecrootc_to_deadcrootc     , s_livecrootc, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livecrootc_xfer_to_livecrootc, s_livecrootc_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livecrootc_storage_to_xfer   , s_livecrootc_storage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadcrootc_xfer_to_deadcrootc, s_deadcrootc_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadcrootc_storage_to_xfer   , s_deadcrootc_storage,nelms)
  call spm_list_insert(spm_list, -1._r8, f_grainc_to_food               , s_grainc, nelms)
  call spm_list_insert(spm_list, -1._r8, f_grainc_xfer_to_grainc        , s_grainc_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_grainc_storage_to_xfer       , s_grainc_storage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_gresp_storage_to_xfer        , s_gresp_storage,nelms)

  call spm_list_to_mat(spm_list, spm_carbon_d, nelms, f_gresp_storage_to_xfer)

  call spm_list_init(spm_list, 0._r8, f_cpool_to_leafc               , s_cpool, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_leafc             , s_leafc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafc_xfer_to_leafc        , s_leafc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafc_storage_to_xfer      , s_leafc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_leafc_storage     , s_leafc_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_frootc            , s_frootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootc_xfer_to_frootc      , s_frootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootc_storage_to_xfer     , s_frootc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_frootc_storage    , s_frootc_sotrage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_livestemc         , s_livestemc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livestemc_xfer_to_livestemc, s_livestemc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livestemc_storage_to_xfer  , s_livestemc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_livestemc_storage , s_livestemc_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadstemc_xfer_to_deadstemc, s_deadstemc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_deadstemc         , s_deadstemc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livestemc_to_deadstemc     , s_deadstemc,nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadstemc_storage_to_xfer  , s_deadstemc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_deadstemc_storage , s_deadstemc_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_livecrootc        , s_livecrootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livecrootc_xfer_to_livecrootc, s_livecrootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livecrootc_storage_to_xfer , s_livecrootc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_livecrootc_storage, s_livecrootc_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_deadcrootc        , s_deadcrootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livecrootc_to_deadcrootc   , s_deadcrootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadcrootc_xfer_to_deadcrootc, s_deadcrootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadcrootc_storage_to_xfer , s_deadcrootc_xfer,nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_deadcrootc_storage, s_deadcrootc_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_grainc_xfer_to_grainc      , s_grainc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_grainc            , s_grainc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_grainc_storage_to_xfer     , s_grainc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_grainc_storage    , s_grainc_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_gresp_storage_to_xfer      , s_gresp_xfer,nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_gresp_storage     , s_gresp_storage, nelms)

  call spm_list_to_mat(spm_list, spm_carbon_p, nelms, f_gresp_storage_to_xfer)

  call spm_list_init(spm_list, -1._r8, f_npool_to_leafn                 , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_leafn_storage       , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_frootn              , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_frootn_storage      , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_livestemn           , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_livestemn_storage   , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_deadstemn           , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_deadstemn_storage   , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_livecrootn          , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_livecrootn_storage  , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_deadcrootn          , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_deadcrootn_storage  , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_grainn              , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_npool_to_grainn_storage      , s_npool, nelms)
  call spm_list_insert(spm_list, -1._r8, f_leafn_to_litter              , s_leafn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_leafn_to_retransn            , s_leafn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_leafn_xfer_to_leafn          , s_leafn_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_leafn_storage_to_xfer        , s_leafn_storage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_frootn_to_litter             , s_frootn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_frootn_to_retransn           , s_frootn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_frootn_xfer_to_frootn        , s_frootn_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_frootn_storage_to_xfer       , s_frootn_sotrage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemn_to_litter          , s_livestemn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemn_to_retransn        , s_livestemn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemn_to_deadstemn       , s_livestemn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemn_xfer_to_livestemn  , s_livestemn_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livestemn_storage_to_xfer    , s_livestemn_storage,nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadstemn_xfer_to_deadstem   , s_deadstemn_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadstemn_storage_to_xfer    , s_deadstemn_storage, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livecrootn_to_deadcrootn     , s_livecrootn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livecrootn_to_retransn       , s_livecrootn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livecrootn_xfer_to_livecrootn, s_livecrootn_xfer, nelms)
  call spm_list_insert(spm_list, -1._r8, f_livecrootn_storage_to_xfer   , s_livecrootn_storage,nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadcrootn_xfer_to_deadcrootn, s_deadcrootn_xfer,nelms)
  call spm_list_insert(spm_list, -1._r8, f_deadcrootn_storage_to_xfer   , s_deadcrootn_storage,nelms)
  call spm_list_insert(spm_list, -1._r8, f_grainn_to_food               , s_grainn, nelms)
  call spm_list_insert(spm_list, -1._r8, f_grainn_xfer_to_grainn        , s_grainn_xfer,nelms)
  call spm_list_insert(spm_list, -1._r8, f_retransn_to_npool            , s_retransn, nelms)

  call spm_list_to_mat(spm_list, spm_nutrient_d, nelms, f_supplement_to_plantn)

  call spm_list_init(spm_list, 1._r8, f_retransn_to_npool                , s_npool, nelms)
  call spm_list_insert(spm_list, 1._r8, f_supplement_to_plantn           , s_npool, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_leafn                 , s_leafn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafn_xfer_to_leafn            , s_leafn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafn_storage_to_xfer          , s_leafn_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_leafn_storage         , s_leafn_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_frootn                , s_frootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootn_xfer_to_frootn          , s_frootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootn_storage_to_xfer         , s_frootn_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_frootn_storage        , s_frootn_sotrage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_livestemn             , s_livestemn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livestemn_xfer_to_livestemn    , s_livestemn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livestemn_storage_to_xfer      , s_livestemn_xfer,nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_livestemn_storage     , s_livestemn_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_deadstemn             , s_deadstemn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livestemn_to_deadstemn         , s_deadstemn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadstemn_xfer_to_deadstem     , s_deadstemn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadstemn_storage_to_xfer      , s_deadstemn_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_deadstemn_storage     , s_deadstemn_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_livecrootn            , s_livecrootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livecrootn_xfer_to_livecrootn  , s_livecrootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livecrootn_storage_to_xfer     , s_livecrootn_xfer,nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_livecrootn_storage    , s_livecrootn_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_deadcrootn            , s_deadcrootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livecrootn_to_deadcrootn       , s_deadcrootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadcrootn_xfer_to_deadcrootn  , s_deadcrootn,nelms)
  call spm_list_insert(spm_list, 1._r8, f_deadcrootn_storage_to_xfer     , s_deadcrootn_xfer,nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_deadcrootn_storage    , s_deadcrootn_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_grainn_xfer_to_grainn          , s_grainn,nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_grainn                , s_grainn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_grainn_storage        , s_grainn_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafn_to_retransn              , s_retransn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootn_to_retransn             , s_retransn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livestemn_to_retransn          , s_retransn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_livecrootn_to_retransn         , s_retransn, nelms)

  call spm_list_to_mat(spm_list, spm_nutrient_p, nelms, f_supplement_to_plantn)
  end subroutine InitPhenoFluxLimiter
!---------------------------------------------------------------------------
  subroutine phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
      carbonflux_vars, carbonstate_vars, &
      c13_carbonflux_vars, c13_carbonstate_vars, &
      c14_carbonflux_vars, c14_carbonstate_vars, &
      nitrogenflux_vars, nitrogenstate_vars, &
      phosphorusflux_vars, phosphorusstate_vars)

    use decompMod           , only : bounds_type
    use CropType                  , only : crop_type
    use CNStateType               , only : cnstate_type
    use CNCarbonFluxType          , only : carbonflux_type
    use CNCarbonStateType         , only : carbonstate_type
    use CNNitrogenFluxType        , only : nitrogenflux_type
    use CNNitrogenStateType       , only : nitrogenstate_type
    use PhosphorusFluxType        , only : phosphorusflux_type
    use PhosphorusStateType       , only : phosphorusstate_type
    use clm_varctl          , only : use_c13, use_c14
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(crop_type)        , intent(inout) :: crop_vars
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(carbonflux_type)  , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type) , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)  , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type) , intent(inout) :: c14_carbonstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars

  call carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      carbonflux_vars, carbonstate_vars)

  if (use_c13)then
    call carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      c13_carbonflux_vars, c13_carbonstate_vars)
  endif

  if (use_c14)then
    call carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      c14_carbonflux_vars, c14_carbonstate_vars)
  endif

  call nitrogen_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      nitrogenflux_vars, nitrogenstate_vars)

  call phosphorus_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      phosphorusflux_vars, phosphorusstate_vars)


  end subroutine phenology_flux_limiter
  !-------------------------------------------------------------------------------

  subroutine carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      carbonflux_vars, carbonstate_vars)

    use decompMod           , only : bounds_type
    use CropType                  , only : crop_type
    use CNCarbonFluxType          , only : carbonflux_type
    use CNCarbonStateType         , only : carbonstate_type
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(crop_type)        , intent(inout) :: crop_vars
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars

  integer :: fp, p
  real(r8) :: ystates(n_carbon_states)
  real(r8) :: rfluxes(n_carbon_fluxes)
  real(r8) :: rscal
  real(r8) :: dt
  real(r8) :: ar_p
  associate(                                                    &
         ivt                   =>    veg_pp%itype             , & ! Input:  [integer  (:)     ]  pft vegetation type
         woody                 =>    veg_vp%woody             , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
         harvdate              =>    crop_vars%harvdate_patch , & ! Input:  [integer  (:)     ]  harvest date
         cf                    => carbonflux_vars             , &
         cs                    => carbonstate_vars              &
         )
  ! set time steps
  dt = real( get_step_size(), r8 )

  do fp = 1,num_soilp
    p = filter_soilp(fp)
    !assemble state variables
    ystates(:)=0._r8
    mov(ystates(s_cpool)              , cs%cpool_patch(p))
    mov(ystates(s_leafc)              , cs%leafc_patch(p))
    mov(ystates(s_leafc_xfer)         , cs%leafc_xfer_patch(p))
    mov(ystates(s_leafc_storage)      , cs%leafc_storage_patch(p))
    mov(ystates(s_frootc)             , cs%frootc_patch(p))
    mov(ystates(s_frootc_xfer)        , cs%frootc_xfer_patch(p))
    mov(ystates(s_frootc_sotrage)     , cs%frootc_storage_patch(p))
    if (woody(ivt(p)) == 1._r8) then
      mov(ystates(s_livestemc)          , cs%livestemc_patch(p))
      mov(ystates(s_livestemc_xfer)     , cs%livestemc_xfer_patch(p))
      mov(ystates(s_livestemc_storage)  , cs%livestemc_storage_patch(p))
      mov(ystates(s_livecrootc)         , cs%livecrootc_patch(p))
      mov(ystates(s_livecrootc_xfer)    , cs%livecrootc_xfer_patch(p))
      mov(ystates(s_livecrootc_storage) , cs%livecrootc_storage_patch(p))
      mov(ystates(s_deadstemc)          , cs%deadstemc_patch(p))
      mov(ystates(s_deadstemc_xfer)     , cs%deadstemc_xfer_patch(p))
      mov(ystates(s_deadstemc_storage)  , cs%deadstemc_storage_patch(p))
      mov(ystates(s_deadcrootc)         , cs%deadcrootc_patch(p))
      mov(ystates(s_deadcrootc_xfer)    , cs%deadcrootc_xfer_patch(p))
      mov(ystates(s_deadcrootc_storage) , cs%deadcrootc_storage_patch(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(ystates(s_livestemc)          , cs%livestemc_patch(p))
      mov(ystates(s_livestemc_xfer)     , cs%livestemc_xfer_patch(p))
      mov(ystates(s_livestemc_storage)  , cs%livestemc_storage_patch(p))
      mov(ystates(s_grainc)             , cs%grainc_patch(p))
      mov(ystates(s_grainc_xfer)        , cs%grainc_xfer_patch(p))
      mov(ystates(s_grainc_storage)     , cs%grainc_storage_patch(p))
    endif
    mov(ystates(s_gresp_xfer)         , cs%gresp_xfer_patch(p))
    mov(ystates(s_gresp_storage)      , cs%gresp_storage_patch(p))

    rfluxes(:) = 0._r8
    ar_p = cf%leaf_curmr_patch(p)             &
         + cf%froot_curmr_patch(p)            &
         + cf%cpool_leaf_gr_patch(p)          &
         + cf%cpool_froot_gr_patch(p)         &
         + cf%cpool_leaf_storage_gr_patch(p)  &
         + cf%cpool_froot_storage_gr_patch(p)
    if (woody(ivt(p)) == 1._r8) then
      ar_p = ar_p                                 &
         + cf%livestem_curmr_patch(p)             &
         + cf%livecroot_curmr_patch(p)            &
         + cf%cpool_livestem_gr_patch(p)          &
         + cf%cpool_deadstem_gr_patch(p)          &
         + cf%cpool_livecroot_gr_patch(p)         &
         + cf%cpool_deadcroot_gr_patch(p)         &
         + cf%cpool_livestem_storage_gr_patch(p)  &
         + cf%cpool_deadstem_storage_gr_patch(p)  &
         + cf%cpool_livecroot_storage_gr_patch(p) &
         + cf%cpool_deadcroot_storage_gr_patch(p)
    endif
    if (ivt(p) >= npcropmin) then
      ar_p= ar_p                                 &
        + cf%livestem_curmr_patch(p)             &
        + cf%grain_curmr_patch(p)                &
        + cf%cpool_livestem_gr_patch(p)          &
        + cf%cpool_grain_gr_patch(p)             &
        + cf%cpool_livestem_storage_gr_patch(p)  &
        + cf%cpool_grain_storage_gr_patch(p)
    endif
    !assemble reactive fluxes
    mov(rfluxes(f_cpool_to_leafc)               , cf%cpool_to_leafc_patch(p))
    mov(rfluxes(f_cpool_to_leafc_storage)       , cf%cpool_to_leafc_storage_patch(p))
    mov(rfluxes(f_cpool_to_frootc)              , cf%cpool_to_frootc_patch(p))
    mov(rfluxes(f_cpool_to_frootc_storage)      , cf%cpool_to_frootc_storage_patch(p))
    mov(rfluxes(f_cpool_to_xsmrpool)            , cf%cpool_to_xsmrpool_patch(p))
    mov(rfluxes(f_cpool_to_gresp_storage)       , cf%cpool_to_gresp_storage_patch(p))
    mov(rfluxes(f_cpool_to_ar)                  , ar_p)
    if (woody(ivt(p)) == 1._r8) then
      mov(rfluxes(f_cpool_to_livestemc)           , cf%cpool_to_livestemc_patch(p))
      mov(rfluxes(f_cpool_to_livestemc_storage)   , cf%cpool_to_livestemc_storage_patch(p))
      mov(rfluxes(f_cpool_to_deadstemc)           , cf%cpool_to_deadstemc_patch(p))
      mov(rfluxes(f_cpool_to_deadstemc_storage)   , cf%cpool_to_deadstemc_storage_patch(p))
      mov(rfluxes(f_cpool_to_livecrootc)          , cf%cpool_to_livecrootc_patch(p))
      mov(rfluxes(f_cpool_to_livecrootc_storage)  , cf%cpool_to_livecrootc_storage_patch(p))
      mov(rfluxes(f_cpool_to_deadcrootc)          , cf%cpool_to_deadcrootc_patch(p))
      mov(rfluxes(f_cpool_to_deadcrootc_storage)  , cf%cpool_to_deadcrootc_storage_patch(p))
      mov(rfluxes(f_livestemc_to_deadstemc)       , cf%livestemc_to_deadstemc_patch(p))
      mov(rfluxes(f_livestemc_xfer_to_livestemc)  , cf%livestemc_xfer_to_livestemc_patch(p))
      mov(rfluxes(f_livestemc_storage_to_xfer)    , cf%livestemc_storage_to_xfer_patch(p))
      mov(rfluxes(f_deadstemc_xfer_to_deadstemc)  , cf%deadstemc_xfer_to_deadstemc_patch(p))
      mov(rfluxes(f_deadstemc_storage_to_xfer)    , cf%deadstemc_storage_to_xfer_patch(p))
      mov(rfluxes(f_livecrootc_xfer_to_livecrootc), cf%livecrootc_xfer_to_livecrootc_patch(p))
      mov(rfluxes(f_livecrootc_storage_to_xfer)   , cf%livecrootc_storage_to_xfer_patch(p))
      mov(rfluxes(f_livecrootc_to_deadcrootc)     , cf%livecrootc_to_deadcrootc_patch(p))
      mov(rfluxes(f_deadcrootc_xfer_to_deadcrootc), cf%deadcrootc_xfer_to_deadcrootc_patch(p))
      mov(rfluxes(f_deadcrootc_storage_to_xfer)   , cf%deadcrootc_storage_to_xfer_patch(p))
      mov(rfluxes(f_gresp_storage_to_xfer)        , cf%gresp_storage_to_xfer_patch(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(rfluxes(f_cpool_to_livestemc)           , cf%cpool_to_livestemc_patch(p))
      mov(rfluxes(f_cpool_to_livestemc_storage)   , cf%cpool_to_livestemc_storage_patch(p))
      mov(rfluxes(f_cpool_to_grainc)              , cf%cpool_to_grainc_patch(p))
      mov(rfluxes(f_cpool_to_grainc_storage)      , cf%cpool_to_grainc_storage_patch(p))
      mov(rfluxes(f_livestemc_xfer_to_livestemc)  , cf%livestemc_xfer_to_livestemc_patch(p))
      mov(rfluxes(f_livestemc_storage_to_xfer)    , cf%livestemc_storage_to_xfer_patch(p))
      mov(rfluxes(f_livestemc_to_litter)          , cf%livestemc_to_litter_patch(p))
      mov(rfluxes(f_grainc_to_food)               , cf%grainc_to_food_patch(p))
      mov(rfluxes(f_grainc_xfer_to_grainc)        , cf%grainc_xfer_to_grainc_patch(p))
      mov(rfluxes(f_grainc_storage_to_xfer)       , cf%grainc_storage_to_xfer_patch(p))
    endif
    mov(rfluxes(f_leafc_to_litter)              , cf%leafc_to_litter_patch(p))
    mov(rfluxes(f_leafc_xfer_to_leafc)          , cf%leafc_xfer_to_leafc_patch(p))
    mov(rfluxes(f_leafc_storage_to_xfer)        , cf%leafc_storage_to_xfer_patch(p))
    mov(rfluxes(f_frootc_to_litter)             , cf%frootc_to_litter_patch(p))
    mov(rfluxes(f_frootc_xfer_to_frootc)        , cf%frootc_xfer_to_frootc_patch(p))
    mov(rfluxes(f_frootc_storage_to_xfer)       , cf%frootc_storage_to_xfer_patch(p))

    !obtain the limiting factor

    call flux_correction(spm_carbon_d%szrow, spm_carbon_d%szcol, spm_carbon_p, &
       spm_carbon_d, dt, ystates, rfluxes)

    !correct the fluxes

    imov(rfluxes(f_cpool_to_leafc)               , cf%cpool_to_leafc_patch(p))
    imov(rfluxes(f_cpool_to_leafc_storage)       , cf%cpool_to_leafc_storage_patch(p))
    imov(rfluxes(f_cpool_to_frootc)              , cf%cpool_to_frootc_patch(p))
    imov(rfluxes(f_cpool_to_frootc_storage)      , cf%cpool_to_frootc_storage_patch(p))
    imov(rfluxes(f_cpool_to_xsmrpool)            , cf%cpool_to_xsmrpool_patch(p))
    imov(rfluxes(f_cpool_to_gresp_storage)       , cf%cpool_to_gresp_storage_patch(p))

    if (woody(ivt(p)) == 1._r8) then
      imov(rfluxes(f_cpool_to_livestemc)           , cf%cpool_to_livestemc_patch(p))
      imov(rfluxes(f_cpool_to_livestemc_storage)   , cf%cpool_to_livestemc_storage_patch(p))
      imov(rfluxes(f_cpool_to_deadstemc)           , cf%cpool_to_deadstemc_patch(p))
      imov(rfluxes(f_cpool_to_deadstemc_storage)   , cf%cpool_to_deadstemc_storage_patch(p))
      imov(rfluxes(f_cpool_to_livecrootc)          , cf%cpool_to_livecrootc_patch(p))
      imov(rfluxes(f_cpool_to_livecrootc_storage)  , cf%cpool_to_livecrootc_storage_patch(p))
      imov(rfluxes(f_cpool_to_deadcrootc)          , cf%cpool_to_deadcrootc_patch(p))
      imov(rfluxes(f_cpool_to_deadcrootc_storage)  , cf%cpool_to_deadcrootc_storage_patch(p))
      imov(rfluxes(f_livestemc_to_deadstemc)       , cf%livestemc_to_deadstemc_patch(p))
      imov(rfluxes(f_livestemc_xfer_to_livestemc)  , cf%livestemc_xfer_to_livestemc_patch(p))
      imov(rfluxes(f_livestemc_storage_to_xfer)    , cf%livestemc_storage_to_xfer_patch(p))
      imov(rfluxes(f_deadstemc_xfer_to_deadstemc)  , cf%deadstemc_xfer_to_deadstemc_patch(p))
      imov(rfluxes(f_deadstemc_storage_to_xfer)    , cf%deadstemc_storage_to_xfer_patch(p))
      imov(rfluxes(f_livecrootc_xfer_to_livecrootc), cf%livecrootc_xfer_to_livecrootc_patch(p))
      imov(rfluxes(f_livecrootc_storage_to_xfer)   , cf%livecrootc_storage_to_xfer_patch(p))
      imov(rfluxes(f_livecrootc_to_deadcrootc)     , cf%livecrootc_to_deadcrootc_patch(p))
      imov(rfluxes(f_deadcrootc_xfer_to_deadcrootc), cf%deadcrootc_xfer_to_deadcrootc_patch(p))
      imov(rfluxes(f_deadcrootc_storage_to_xfer)   , cf%deadcrootc_storage_to_xfer_patch(p))
      imov(rfluxes(f_gresp_storage_to_xfer)        , cf%gresp_storage_to_xfer_patch(p))
    endif
    if (ivt(p) >= npcropmin) then
      imov(rfluxes(f_cpool_to_livestemc)           , cf%cpool_to_livestemc_patch(p))
      imov(rfluxes(f_cpool_to_livestemc_storage)   , cf%cpool_to_livestemc_storage_patch(p))
      imov(rfluxes(f_cpool_to_grainc)              , cf%cpool_to_grainc_patch(p))
      imov(rfluxes(f_cpool_to_grainc_storage)      , cf%cpool_to_grainc_storage_patch(p))
      imov(rfluxes(f_livestemc_xfer_to_livestemc)  , cf%livestemc_xfer_to_livestemc_patch(p))
      imov(rfluxes(f_livestemc_storage_to_xfer)    , cf%livestemc_storage_to_xfer_patch(p))
      imov(rfluxes(f_livestemc_to_litter)          , cf%livestemc_to_litter_patch(p))
      imov(rfluxes(f_grainc_to_food)               , cf%grainc_to_food_patch(p))
      imov(rfluxes(f_grainc_xfer_to_grainc)        , cf%grainc_xfer_to_grainc_patch(p))
      imov(rfluxes(f_grainc_storage_to_xfer)       , cf%grainc_storage_to_xfer_patch(p))
    endif
    imov(rfluxes(f_leafc_to_litter)              , cf%leafc_to_litter_patch(p))
    imov(rfluxes(f_leafc_xfer_to_leafc)          , cf%leafc_xfer_to_leafc_patch(p))
    imov(rfluxes(f_leafc_storage_to_xfer)        , cf%leafc_storage_to_xfer_patch(p))
    imov(rfluxes(f_frootc_to_litter)             , cf%frootc_to_litter_patch(p))
    imov(rfluxes(f_frootc_xfer_to_frootc)        , cf%frootc_xfer_to_frootc_patch(p))
    imov(rfluxes(f_frootc_storage_to_xfer)       , cf%frootc_storage_to_xfer_patch(p))

    if(rfluxes(f_cpool_to_ar)  <  ar_p)then
      rscal= rfluxes(f_cpool_to_ar)/ar_p
      ascal(cf%leaf_curmr_patch(p)                , rscal)
      ascal(cf%froot_curmr_patch(p)               , rscal)
      ascal(cf%cpool_leaf_gr_patch(p)             , rscal)
      ascal(cf%cpool_froot_gr_patch(p)            , rscal)
      ascal(cf%cpool_leaf_storage_gr_patch(p)     , rscal)
      ascal(cf%cpool_froot_storage_gr_patch(p)    , rscal)
      if (woody(ivt(p)) == 1._r8) then
        ascal(cf%livestem_curmr_patch(p)            , rscal)
        ascal(cf%livecroot_curmr_patch(p)           , rscal)
        ascal(cf%cpool_livestem_gr_patch(p)         , rscal)
        ascal(cf%cpool_deadstem_gr_patch(p)         , rscal)
        ascal(cf%cpool_livecroot_gr_patch(p)        , rscal)
        ascal(cf%cpool_deadcroot_gr_patch(p)        , rscal)
        ascal(cf%cpool_livestem_storage_gr_patch(p) , rscal)
        ascal(cf%cpool_deadstem_storage_gr_patch(p) , rscal)
        ascal(cf%cpool_livecroot_storage_gr_patch(p), rscal)
        ascal(cf%cpool_deadcroot_storage_gr_patch(p), rscal)
      endif
      if (ivt(p) >= npcropmin) then
        ascal(cf%livestem_curmr_patch(p)            , rscal)
        ascal(cf%grain_curmr_patch(p)               , rscal)
        ascal(cf%cpool_livestem_gr_patch(p)         , rscal)
        ascal(cf%cpool_grain_gr_patch(p)            , rscal)
        ascal(cf%cpool_livestem_storage_gr_patch(p) , rscal)
        ascal(cf%cpool_grain_storage_gr_patch(p)    , rscal)
      endif
    endif
  enddo
  end associate
  end subroutine carbon_flux_limiter
!---------------------------------------------------------------------------
  subroutine nitrogen_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      nitrogenflux_vars, nitrogenstate_vars)

    use decompMod           , only : bounds_type
    use CNStateType               , only : cnstate_type
    use CNNitrogenFluxType        , only : nitrogenflux_type
    use CNNitrogenStateType       , only : nitrogenstate_type
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars

  integer :: fp, p
  real(r8) :: ystates(n_nutrient_states)
  real(r8) :: rfluxes(n_nutrient_fluxes)
  real(r8) :: dt
  associate(                                                       &
         ivt                   => veg_pp%itype                   , & ! Input:  [integer  (:)     ]  pft vegetation type
         woody                 => veg_vp%woody                   , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
         nf                    => nitrogenflux_vars              , &
         ns                    => nitrogenstate_vars               &
  )
  ! set time steps
  dt = real( get_step_size(), r8 )
  do fp = 1,num_soilp
    p = filter_soilp(fp)
    ystates(:) = 0._r8
    !assemble state variables
    mov(ystates(s_npool)               , ns%npool_patch(p))
    mov(ystates(s_leafn)               , ns%leafn_patch(p))
    mov(ystates(s_leafn_xfer)          , ns%leafn_xfer_patch(p))
    mov(ystates(s_leafn_storage)       , nf%npool_to_leafn_storage_patch(p))
    mov(ystates(s_frootn)              , ns%frootn_patch(p))
    mov(ystates(s_frootn_xfer)         , ns%frootn_xfer_patch(p))
    mov(ystates(s_frootn_sotrage)      , ns%frootn_storage_patch(p))
    if (woody(ivt(p)) == 1.0_r8) then
      mov(ystates(s_livestemn)           , ns%livestemn_patch(p))
      mov(ystates(s_livestemn_xfer)      , ns%livestemn_xfer_patch(p))
      mov(ystates(s_livestemn_storage)   , ns%livestemn_storage_patch(p))
      mov(ystates(s_deadstemn)           , ns%deadstemn_patch(p))
      mov(ystates(s_deadstemn_xfer)      , ns%deadstemn_xfer_patch(p))
      mov(ystates(s_deadstemn_storage)   , ns%deadstemn_storage_patch(p))
      mov(ystates(s_livecrootn)          , ns%livecrootn_patch(p))
      mov(ystates(s_livecrootn_xfer)     , ns%livecrootn_xfer_patch(p))
      mov(ystates(s_livecrootn_storage)  , ns%livecrootn_storage_patch(p))
      mov(ystates(s_deadcrootn)          , ns%deadcrootn_patch(p))
      mov(ystates(s_deadcrootn_xfer)     , ns%deadcrootn_xfer_patch(p))
      mov(ystates(s_deadcrootn_storage)  , ns%deadcrootn_storage_patch(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(ystates(s_grainn)              , ns%grainn_patch(p))
      mov(ystates(s_grainn_xfer)         , ns%grainn_xfer_patch(p))
      mov(ystates(s_grainn_storage)      , ns%grainn_storage_patch(p))
      mov(ystates(s_livestemn)           , ns%livestemn_patch(p))
      mov(ystates(s_livestemn_xfer)      , ns%livestemn_xfer_patch(p))
      mov(ystates(s_livestemn_storage)   , ns%livestemn_storage_patch(p))
    endif
    mov(ystates(s_retransn)            , ns%retransn_patch(p))

    rfluxes(:) = 0._r8
    !assemble reactive fluxes
    mov(rfluxes(f_npool_to_leafn)                , nf%npool_to_leafn_patch(p))
    mov(rfluxes(f_npool_to_leafn_storage)        , nf%npool_to_leafn_storage_patch(p))
    mov(rfluxes(f_npool_to_frootn)               , nf%npool_to_frootn_patch(p))
    mov(rfluxes(f_npool_to_frootn_storage)       , nf%npool_to_frootn_storage_patch(p))
    if (woody(ivt(p)) == 1.0_r8) then
      mov(rfluxes(f_npool_to_livestemn)            , nf%npool_to_livestemn_patch(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , nf%npool_to_livestemn_storage_patch(p))
      mov(rfluxes(f_npool_to_livecrootn)           , nf%npool_to_livecrootn_patch(p))
      mov(rfluxes(f_npool_to_livecrootn_storage)   , nf%npool_to_livecrootn_storage_patch(p))
      mov(rfluxes(f_npool_to_deadstemn)            , nf%npool_to_deadstemn_patch(p))
      mov(rfluxes(f_npool_to_deadcrootn)           , nf%npool_to_deadcrootn_patch(p))
      mov(rfluxes(f_npool_to_deadstemn_storage)    , nf%npool_to_deadstemn_storage_patch(p))
      mov(rfluxes(f_npool_to_deadcrootn_storage)   , nf%npool_to_deadcrootn_storage_patch(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , nf%livestemn_storage_to_xfer_patch(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   , nf%livestemn_xfer_to_livestemn_patch(p))
      mov(rfluxes(f_livestemn_to_retransn)         , nf%livestemn_to_retransn_patch(p))
      mov(rfluxes(f_livestemn_to_deadstemn)        , nf%livestemn_to_deadstemn_patch(p))
      mov(rfluxes(f_livecrootn_to_deadcrootn)      , nf%livecrootn_to_deadcrootn_patch(p))
      mov(rfluxes(f_livecrootn_to_retransn)        , nf%livecrootn_to_retransn_patch(p))
      mov(rfluxes(f_livecrootn_storage_to_xfer)    , nf%livecrootn_storage_to_xfer_patch(p))
      mov(rfluxes(f_livecrootn_xfer_to_livecrootn) , nf%livecrootn_xfer_to_livecrootn_patch(p))
      mov(rfluxes(f_deadstemn_storage_to_xfer)     , nf%deadstemn_storage_to_xfer_patch(p))
      mov(rfluxes(f_deadstemn_xfer_to_deadstem)    , nf%deadstemn_xfer_to_deadstemn_patch(p))
      mov(rfluxes(f_deadcrootn_storage_to_xfer)    , nf%deadcrootn_storage_to_xfer_patch(p))
      mov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , nf%deadcrootn_xfer_to_deadcrootn_patch(p))
    endif

    if (ivt(p) >= npcropmin) then
      mov(rfluxes(f_npool_to_livestemn)            , nf%npool_to_livestemn_patch(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , nf%npool_to_livestemn_storage_patch(p))
      mov(rfluxes(f_npool_to_grainn)               , nf%npool_to_grainn_patch(p))
      mov(rfluxes(f_npool_to_grainn_storage)       , nf%npool_to_grainn_storage_patch(p))
      mov(rfluxes(f_frootn_to_retransn)            , nf%frootn_to_retransn_patch(p))
      mov(rfluxes(f_livestemn_to_litter)           , nf%livestemn_to_litter_patch(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , nf%livestemn_storage_to_xfer_patch(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   ,  nf%livestemn_xfer_to_livestemn_patch(p))
      mov(rfluxes(f_livestemn_to_retransn)         , nf%livestemn_to_retransn_patch(p))
      mov(rfluxes(f_grainn_xfer_to_grainn)         , nf%grainn_xfer_to_grainn_patch(p))
      mov(rfluxes(f_grainn_to_food)                , nf%grainn_to_food_patch(p))
    endif

    mov(rfluxes(f_leafn_to_litter)               , nf%leafn_to_litter_patch(p))
    mov(rfluxes(f_leafn_xfer_to_leafn)           , nf%leafn_xfer_to_leafn_patch(p))
    mov(rfluxes(f_leafn_storage_to_xfer)         , nf%leafn_storage_to_xfer_patch(p))
    mov(rfluxes(f_frootn_xfer_to_frootn)         , nf%frootn_xfer_to_frootn_patch(p))
    mov(rfluxes(f_frootn_storage_to_xfer)        , nf%frootn_storage_to_xfer_patch(p))
    mov(rfluxes(f_leafn_to_retransn)             , nf%leafn_to_retransn_patch(p))
    mov(rfluxes(f_frootn_to_litter)              , nf%frootn_to_litter_patch(p))

    mov(rfluxes(f_retransn_to_npool)             , nf%retransn_to_npool_patch(p))
    mov(rfluxes(f_supplement_to_plantn)          , nf%supplement_to_plantn(p))

    call flux_correction(spm_nutrient_d%szrow, spm_nutrient_d%szcol, spm_nutrient_p, &
      spm_nutrient_d, dt, ystates, rfluxes)

    !correct the fluxes
    imov(rfluxes(f_npool_to_leafn)                , nf%npool_to_leafn_patch(p))
    imov(rfluxes(f_npool_to_leafn_storage)        , nf%npool_to_leafn_storage_patch(p))
    imov(rfluxes(f_npool_to_frootn)               , nf%npool_to_frootn_patch(p))
    imov(rfluxes(f_npool_to_frootn_storage)       , nf%npool_to_frootn_storage_patch(p))
    if (woody(ivt(p)) == 1.0_r8) then
      imov(rfluxes(f_npool_to_livestemn)            , nf%npool_to_livestemn_patch(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , nf%npool_to_livestemn_storage_patch(p))
      imov(rfluxes(f_npool_to_livecrootn)           , nf%npool_to_livecrootn_patch(p))
      imov(rfluxes(f_npool_to_livecrootn_storage)   , nf%npool_to_livecrootn_storage_patch(p))
      imov(rfluxes(f_npool_to_deadstemn)            , nf%npool_to_deadstemn_patch(p))
      imov(rfluxes(f_npool_to_deadcrootn)           , nf%npool_to_deadcrootn_patch(p))
      imov(rfluxes(f_npool_to_deadstemn_storage)    , nf%npool_to_deadstemn_storage_patch(p))
      imov(rfluxes(f_npool_to_deadcrootn_storage)   , nf%npool_to_deadcrootn_storage_patch(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , nf%livestemn_storage_to_xfer_patch(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   , nf%livestemn_xfer_to_livestemn_patch(p))
      imov(rfluxes(f_livestemn_to_retransn)         , nf%livestemn_to_retransn_patch(p))
      imov(rfluxes(f_livestemn_to_deadstemn)        , nf%livestemn_to_deadstemn_patch(p))
      imov(rfluxes(f_livecrootn_to_deadcrootn)      , nf%livecrootn_to_deadcrootn_patch(p))
      imov(rfluxes(f_livecrootn_to_retransn)        , nf%livecrootn_to_retransn_patch(p))
      imov(rfluxes(f_livecrootn_storage_to_xfer)    , nf%livecrootn_storage_to_xfer_patch(p))
      imov(rfluxes(f_livecrootn_xfer_to_livecrootn) , nf%livecrootn_xfer_to_livecrootn_patch(p))
      imov(rfluxes(f_deadstemn_storage_to_xfer)     , nf%deadstemn_storage_to_xfer_patch(p))
      imov(rfluxes(f_deadstemn_xfer_to_deadstem)    , nf%deadstemn_xfer_to_deadstemn_patch(p))
      imov(rfluxes(f_deadcrootn_storage_to_xfer)    , nf%deadcrootn_storage_to_xfer_patch(p))
      imov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , nf%deadcrootn_xfer_to_deadcrootn_patch(p))
    endif

    if (ivt(p) >= npcropmin) then
      imov(rfluxes(f_npool_to_livestemn)            , nf%npool_to_livestemn_patch(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , nf%npool_to_livestemn_storage_patch(p))
      imov(rfluxes(f_npool_to_grainn)               , nf%npool_to_grainn_patch(p))
      imov(rfluxes(f_npool_to_grainn_storage)       , nf%npool_to_grainn_storage_patch(p))
      imov(rfluxes(f_frootn_to_retransn)            , nf%frootn_to_retransn_patch(p))
      imov(rfluxes(f_livestemn_to_litter)           , nf%livestemn_to_litter_patch(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , nf%livestemn_storage_to_xfer_patch(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   ,  nf%livestemn_xfer_to_livestemn_patch(p))
      imov(rfluxes(f_livestemn_to_retransn)         , nf%livestemn_to_retransn_patch(p))
      imov(rfluxes(f_grainn_xfer_to_grainn)         , nf%grainn_xfer_to_grainn_patch(p))
      imov(rfluxes(f_grainn_to_food)                , nf%grainn_to_food_patch(p))
    endif

    mov(rfluxes(f_leafn_to_litter)               , nf%leafn_to_litter_patch(p))
    mov(rfluxes(f_leafn_xfer_to_leafn)           , nf%leafn_xfer_to_leafn_patch(p))
    mov(rfluxes(f_leafn_storage_to_xfer)         , nf%leafn_storage_to_xfer_patch(p))
    mov(rfluxes(f_frootn_xfer_to_frootn)         , nf%frootn_xfer_to_frootn_patch(p))
    mov(rfluxes(f_frootn_storage_to_xfer)        , nf%frootn_storage_to_xfer_patch(p))
    mov(rfluxes(f_leafn_to_retransn)             , nf%leafn_to_retransn_patch(p))
    mov(rfluxes(f_frootn_to_litter)              , nf%frootn_to_litter_patch(p))

    mov(rfluxes(f_retransn_to_npool)             , nf%retransn_to_npool_patch(p))
    mov(rfluxes(f_supplement_to_plantn)          , nf%supplement_to_plantn(p))

  enddo
  end associate
  end subroutine nitrogen_flux_limiter

  !-------------------------------------------------------------------------------
  subroutine phosphorus_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      phosphorusflux_vars, phosphorusstate_vars)

    use decompMod           , only : bounds_type
    use CNStateType               , only : cnstate_type
    use PhosphorusFluxType        , only : phosphorusflux_type
    use PhosphorusStateType       , only : phosphorusstate_type
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars

  integer :: fp, p
  real(r8) :: ystates(n_nutrient_states)
  real(r8) :: rfluxes(n_nutrient_fluxes)

  real(r8) :: dt

  associate(                                                     &
    ivt                   => veg_pp%itype                      , & ! Input:  [integer  (:)     ]  pft vegetation type
    woody                 => veg_vp%woody                      , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
    pf                    => phosphorusflux_vars               , &
    ps                    => phosphorusstate_vars                &
  )

  ! set time steps
  dt = real( get_step_size(), r8 )
  do fp = 1,num_soilp
    p = filter_soilp(fp)
    ystates(:) = 0._r8
    !assemble state variables
    mov(ystates(s_npool)               , ps%ppool_patch(p))
    mov(ystates(s_leafn)               , ps%leafp_patch(p))
    mov(ystates(s_leafn_xfer)          , ps%leafp_xfer_patch(p))
    mov(ystates(s_leafn_storage)       , ps%leafp_storage_patch(p))
    mov(ystates(s_frootn)              , ps%frootp_patch(p))
    mov(ystates(s_frootn_xfer)         , ps%frootp_xfer_patch(p))
    mov(ystates(s_frootn_sotrage)      , pf%ppool_to_frootp_storage_patch(p))
    if (woody(ivt(p)) == 1.0_r8) then
      mov(ystates(s_livestemn)           , ps%livestemp_patch(p))
      mov(ystates(s_livestemn_xfer)      , ps%livestemp_xfer_patch(p))
      mov(ystates(s_livestemn_storage)   , ps%livestemp_storage_patch(p))
      mov(ystates(s_deadstemn)           , ps%deadstemp_patch(p))
      mov(ystates(s_deadstemn_xfer)      , ps%deadstemp_xfer_patch(p))
      mov(ystates(s_deadstemn_storage)   , ps%deadstemp_storage_patch(p))
      mov(ystates(s_livecrootn)          , ps%livecrootp_patch(p))
      mov(ystates(s_livecrootn_xfer)     , ps%livecrootp_xfer_patch(p))
      mov(ystates(s_livecrootn_storage)  , ps%livecrootp_storage_patch(p))
      mov(ystates(s_deadcrootn)          , ps%deadcrootp_patch(p))
      mov(ystates(s_deadcrootn_xfer)     , ps%deadcrootp_xfer_patch(p))
      mov(ystates(s_deadcrootn_storage)  , ps%deadcrootp_storage_patch(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(ystates(s_grainn)              , ps%grainp_patch(p))
      mov(ystates(s_grainn_xfer)         , ps%grainp_xfer_patch(p))
      mov(ystates(s_grainn_storage)      , ps%grainp_storage_patch(p))
      mov(ystates(s_livestemn)           , ps%livestemp_patch(p))
      mov(ystates(s_livestemn_xfer)      , ps%livestemp_xfer_patch(p))
      mov(ystates(s_livestemn_storage)   , ps%livestemp_storage_patch(p))
    endif
    mov(ystates(s_retransn)              , ps%retransp_patch(p))

    rfluxes(:) = 0._r8
    !assemble reactive fluxes
    mov(rfluxes(f_npool_to_leafn)                , pf%ppool_to_leafp_patch(p))
    mov(rfluxes(f_npool_to_leafn_storage)        , pf%ppool_to_leafp_storage_patch(p))
    mov(rfluxes(f_npool_to_frootn)               , pf%ppool_to_frootp_patch(p))
    mov(rfluxes(f_npool_to_frootn_storage)       , pf%ppool_to_frootp_storage_patch(p))
    if (woody(ivt(p)) == 1._r8) then
      mov(rfluxes(f_npool_to_livestemn)            , pf%ppool_to_livestemp_patch(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , pf%ppool_to_livestemp_storage_patch(p))
      mov(rfluxes(f_npool_to_livecrootn)           , pf%ppool_to_livecrootp_patch(p))
      mov(rfluxes(f_npool_to_livecrootn_storage)   , pf%ppool_to_livecrootp_storage_patch(p))
      mov(rfluxes(f_npool_to_deadstemn)            , pf%ppool_to_deadstemp_patch(p))
      mov(rfluxes(f_npool_to_deadcrootn)           , pf%ppool_to_deadcrootp_patch(p))
      mov(rfluxes(f_npool_to_deadstemn_storage)    , pf%ppool_to_deadstemp_storage_patch(p))
      mov(rfluxes(f_npool_to_deadcrootn_storage)   , pf%ppool_to_deadcrootp_storage_patch(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , pf%livestemp_storage_to_xfer_patch(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   , pf%livestemp_xfer_to_livestemp_patch(p))
      mov(rfluxes(f_livestemn_to_retransn)         , pf%livestemp_to_retransp_patch(p))
      mov(rfluxes(f_livestemn_to_deadstemn)        , pf%livestemp_to_deadstemp_patch(p))
      mov(rfluxes(f_livecrootn_to_deadcrootn)      , pf%livecrootp_to_deadcrootp_patch(p))
      mov(rfluxes(f_livecrootn_to_retransn)        , pf%livecrootp_to_retransp_patch(p))
      mov(rfluxes(f_livecrootn_storage_to_xfer)    , pf%livecrootp_storage_to_xfer_patch(p))
      mov(rfluxes(f_livecrootn_xfer_to_livecrootn) , pf%livecrootp_xfer_to_livecrootp_patch(p))
      mov(rfluxes(f_deadstemn_storage_to_xfer)     , pf%deadstemp_storage_to_xfer_patch(p))
      mov(rfluxes(f_deadstemn_xfer_to_deadstem)    , pf%deadstemp_xfer_to_deadstemp_patch(p))
      mov(rfluxes(f_deadcrootn_storage_to_xfer)    , pf%deadcrootp_storage_to_xfer_patch(p))
      mov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , pf%deadcrootp_xfer_to_deadcrootp_patch(p))
    endif

    if (ivt(p) >= npcropmin) then
      mov(rfluxes(f_npool_to_livestemn)            , pf%ppool_to_livestemp_patch(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , pf%ppool_to_livestemp_storage_patch(p))
      mov(rfluxes(f_npool_to_grainn)               , pf%ppool_to_grainp_patch(p))
      mov(rfluxes(f_npool_to_grainn_storage)       , pf%ppool_to_grainp_storage_patch(p))
      mov(rfluxes(f_frootn_to_retransn)            , pf%frootp_to_retransp_patch(p))
      mov(rfluxes(f_livestemn_to_litter)           , pf%livestemp_to_litter_patch(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , pf%livestemp_storage_to_xfer_patch(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   , pf%livestemp_xfer_to_livestemp_patch(p))
      mov(rfluxes(f_livestemn_to_retransn)         , pf%livestemp_to_retransp_patch(p))
      mov(rfluxes(f_grainn_xfer_to_grainn)         , pf%grainp_xfer_to_grainp_patch(p))
      mov(rfluxes(f_grainn_to_food)                , pf%grainp_to_food_patch(p))
    endif

    mov(rfluxes(f_leafn_to_litter)               , pf%leafp_to_litter_patch(p))
    mov(rfluxes(f_leafn_xfer_to_leafn)           , pf%leafp_xfer_to_leafp_patch(p))
    mov(rfluxes(f_leafn_storage_to_xfer)         , pf%ppool_to_leafp_storage_patch(p))
    mov(rfluxes(f_frootn_xfer_to_frootn)         , pf%frootp_xfer_to_frootp_patch(p))
    mov(rfluxes(f_frootn_storage_to_xfer)        , pf%frootp_storage_to_xfer_patch(p))
    mov(rfluxes(f_leafn_to_retransn)             , pf%leafp_to_retransp_patch(p))
    mov(rfluxes(f_frootn_to_litter)              , pf%frootp_to_litter_patch(p))
    mov(rfluxes(f_retransn_to_npool)             , pf%retransp_to_ppool_patch(p))
    mov(rfluxes(f_supplement_to_plantn)          , pf%supplement_to_plantp(p))
    !assemble stoichiometry matrix

    !obtain the limiting factor
    call flux_correction(spm_nutrient_d%szrow, spm_nutrient_d%szcol, spm_nutrient_p, &
      spm_nutrient_d, dt, ystates, rfluxes)
    !correct the fluxes
    imov(rfluxes(f_npool_to_leafn)                , pf%ppool_to_leafp_patch(p))
    imov(rfluxes(f_npool_to_leafn_storage)        , pf%ppool_to_leafp_storage_patch(p))
    imov(rfluxes(f_npool_to_frootn)               , pf%ppool_to_frootp_patch(p))
    imov(rfluxes(f_npool_to_frootn_storage)       , pf%ppool_to_frootp_storage_patch(p))
    if (woody(ivt(p)) == 1._r8) then
      imov(rfluxes(f_npool_to_livestemn)            , pf%ppool_to_livestemp_patch(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , pf%ppool_to_livestemp_storage_patch(p))
      imov(rfluxes(f_npool_to_livecrootn)           , pf%ppool_to_livecrootp_patch(p))
      imov(rfluxes(f_npool_to_livecrootn_storage)   , pf%ppool_to_livecrootp_storage_patch(p))
      imov(rfluxes(f_npool_to_deadstemn)            , pf%ppool_to_deadstemp_patch(p))
      imov(rfluxes(f_npool_to_deadcrootn)           , pf%ppool_to_deadcrootp_patch(p))
      imov(rfluxes(f_npool_to_deadstemn_storage)    , pf%ppool_to_deadstemp_storage_patch(p))
      imov(rfluxes(f_npool_to_deadcrootn_storage)   , pf%ppool_to_deadcrootp_storage_patch(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , pf%livestemp_storage_to_xfer_patch(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   , pf%livestemp_xfer_to_livestemp_patch(p))
      imov(rfluxes(f_livestemn_to_retransn)         , pf%livestemp_to_retransp_patch(p))
      imov(rfluxes(f_livestemn_to_deadstemn)        , pf%livestemp_to_deadstemp_patch(p))
      imov(rfluxes(f_livecrootn_to_deadcrootn)      , pf%livecrootp_to_deadcrootp_patch(p))
      imov(rfluxes(f_livecrootn_to_retransn)        , pf%livecrootp_to_retransp_patch(p))
      imov(rfluxes(f_livecrootn_storage_to_xfer)    , pf%livecrootp_storage_to_xfer_patch(p))
      imov(rfluxes(f_livecrootn_xfer_to_livecrootn) , pf%livecrootp_xfer_to_livecrootp_patch(p))
      imov(rfluxes(f_deadstemn_storage_to_xfer)     , pf%deadstemp_storage_to_xfer_patch(p))
      imov(rfluxes(f_deadstemn_xfer_to_deadstem)    , pf%deadstemp_xfer_to_deadstemp_patch(p))
      imov(rfluxes(f_deadcrootn_storage_to_xfer)    , pf%deadcrootp_storage_to_xfer_patch(p))
      imov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , pf%deadcrootp_xfer_to_deadcrootp_patch(p))
    endif

    if (ivt(p) >= npcropmin) then
      imov(rfluxes(f_npool_to_livestemn)            , pf%ppool_to_livestemp_patch(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , pf%ppool_to_livestemp_storage_patch(p))
      imov(rfluxes(f_npool_to_grainn)               , pf%ppool_to_grainp_patch(p))
      imov(rfluxes(f_npool_to_grainn_storage)       , pf%ppool_to_grainp_storage_patch(p))
      imov(rfluxes(f_frootn_to_retransn)            , pf%frootp_to_retransp_patch(p))
      imov(rfluxes(f_livestemn_to_litter)           , pf%livestemp_to_litter_patch(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , pf%livestemp_storage_to_xfer_patch(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   , pf%livestemp_xfer_to_livestemp_patch(p))
      imov(rfluxes(f_livestemn_to_retransn)         , pf%livestemp_to_retransp_patch(p))
      imov(rfluxes(f_grainn_xfer_to_grainn)         , pf%grainp_xfer_to_grainp_patch(p))
      imov(rfluxes(f_grainn_to_food)                , pf%grainp_to_food_patch(p))
    endif

    imov(rfluxes(f_leafn_to_litter)               , pf%leafp_to_litter_patch(p))
    imov(rfluxes(f_leafn_xfer_to_leafn)           , pf%leafp_xfer_to_leafp_patch(p))
    imov(rfluxes(f_leafn_storage_to_xfer)         , pf%ppool_to_leafp_storage_patch(p))
    imov(rfluxes(f_frootn_xfer_to_frootn)         , pf%frootp_xfer_to_frootp_patch(p))
    imov(rfluxes(f_frootn_storage_to_xfer)        , pf%frootp_storage_to_xfer_patch(p))
    imov(rfluxes(f_leafn_to_retransn)             , pf%leafp_to_retransp_patch(p))
    imov(rfluxes(f_frootn_to_litter)              , pf%frootp_to_litter_patch(p))
    imov(rfluxes(f_retransn_to_npool)             , pf%retransp_to_ppool_patch(p))
    imov(rfluxes(f_supplement_to_plantn)          , pf%supplement_to_plantp(p))

  enddo
  end associate
  end subroutine phosphorus_flux_limiter

end module PhenologyFLuxLimitMod
