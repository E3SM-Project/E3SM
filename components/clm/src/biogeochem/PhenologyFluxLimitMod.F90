module PhenologyFLuxLimitMod

!DESCRIPTION
! limit the allocation fluxes resulting from pheonology
! calculation to avoid negative state varaibles for
! carbon, nitrogen and phosphours.

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
      veg_cf, veg_cs , c13_veg_cf, c13_veg_cs , c14_veg_cf, c14_veg_cs , &
      veg_nf, veg_ns, veg_pf, veg_ps)

    use decompMod           , only : bounds_type
    use CropType                  , only : crop_type
    use CNStateType               , only : cnstate_type
    use VegetationDataType, only : vegetation_carbon_flux
    use VegetationDataType, only : vegetation_carbon_state
    use VegetationDataType, only : vegetation_nitrogen_flux
    use VegetationDataType, only : vegetation_nitrogen_state
    use VegetationDataType, only : vegetation_phosphorus_flux
    use VegetationDataType, only : vegetation_phosphorus_state
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
    type(vegetation_carbon_flux) , intent(inout) :: veg_cf
    type(vegetation_carbon_state), intent(inout) :: veg_cs
    type(vegetation_carbon_flux)  , intent(inout) :: c13_veg_cf
    type(vegetation_carbon_state) , intent(inout) :: c13_veg_cs
    type(vegetation_carbon_flux)  , intent(inout) :: c14_veg_cf
    type(vegetation_carbon_state) , intent(inout) :: c14_veg_cs
    type(vegetation_nitrogen_flux)  , intent(inout) :: veg_nf
    type(vegetation_nitrogen_state) , intent(inout) :: veg_ns
    type(vegetation_phosphorus_flux)  , intent(inout) :: veg_pf
    type(vegetation_phosphorus_state) , intent(inout) :: veg_ps

  call carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      veg_cf, veg_cs )

  if (use_c13)then
    call carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      c13_veg_cf, c13_veg_cs )
  endif

  if (use_c14)then
    call carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      c14_veg_cf, c14_veg_cs )
  endif

  call nitrogen_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_nf, veg_ns)

  call phosphorus_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_pf, veg_ps)


  end subroutine phenology_flux_limiter
  !-------------------------------------------------------------------------------

  subroutine carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      veg_cf, veg_cs )

    use decompMod           , only : bounds_type
    use CropType                  , only : crop_type
    use VegetationDataType, only : vegetation_carbon_flux
    use VegetationDataType, only : vegetation_carbon_state
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(crop_type)        , intent(inout) :: crop_vars
    type(vegetation_carbon_flux) , intent(inout) :: veg_cf
    type(vegetation_carbon_state), intent(inout) :: veg_cs

  integer :: fp, p
  real(r8) :: ystates(n_carbon_states)
  real(r8) :: rfluxes(n_carbon_fluxes)
  real(r8) :: rscal
  real(r8) :: dt
  real(r8) :: ar_p
  associate(                                                    &
         ivt                   =>    veg_pp%itype             , & ! Input:  [integer  (:)     ]  pft vegetation type
         woody                 =>    veg_vp%woody             , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
         harvdate              =>    crop_vars%harvdate_patch   & ! Input:  [integer  (:)     ]  harvest date
         )
  ! set time steps
  dt = real( get_step_size(), r8 )

  do fp = 1,num_soilp
    p = filter_soilp(fp)
    !assemble state variables
    ystates(:)=0._r8
    mov(ystates(s_cpool)              , veg_cs%cpool(p))
    mov(ystates(s_leafc)              , veg_cs%leafc(p))
    mov(ystates(s_leafc_xfer)         , veg_cs%leafc_xfer(p))
    mov(ystates(s_leafc_storage)      , veg_cs%leafc_storage(p))
    mov(ystates(s_frootc)             , veg_cs%frootc(p))
    mov(ystates(s_frootc_xfer)        , veg_cs%frootc_xfer(p))
    mov(ystates(s_frootc_sotrage)     , veg_cs%frootc_storage(p))
    if (woody(ivt(p)) == 1._r8) then
      mov(ystates(s_livestemc)          , veg_cs%livestemc(p))
      mov(ystates(s_livestemc_xfer)     , veg_cs%livestemc_xfer(p))
      mov(ystates(s_livestemc_storage)  , veg_cs%livestemc_storage(p))
      mov(ystates(s_livecrootc)         , veg_cs%livecrootc(p))
      mov(ystates(s_livecrootc_xfer)    , veg_cs%livecrootc_xfer(p))
      mov(ystates(s_livecrootc_storage) , veg_cs%livecrootc_storage(p))
      mov(ystates(s_deadstemc)          , veg_cs%deadstemc(p))
      mov(ystates(s_deadstemc_xfer)     , veg_cs%deadstemc_xfer(p))
      mov(ystates(s_deadstemc_storage)  , veg_cs%deadstemc_storage(p))
      mov(ystates(s_deadcrootc)         , veg_cs%deadcrootc(p))
      mov(ystates(s_deadcrootc_xfer)    , veg_cs%deadcrootc_xfer(p))
      mov(ystates(s_deadcrootc_storage) , veg_cs%deadcrootc_storage(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(ystates(s_livestemc)          , veg_cs%livestemc(p))
      mov(ystates(s_livestemc_xfer)     , veg_cs%livestemc_xfer(p))
      mov(ystates(s_livestemc_storage)  , veg_cs%livestemc_storage(p))
      mov(ystates(s_grainc)             , veg_cs%grainc(p))
      mov(ystates(s_grainc_xfer)        , veg_cs%grainc_xfer(p))
      mov(ystates(s_grainc_storage)     , veg_cs%grainc_storage(p))
    endif
    mov(ystates(s_gresp_xfer)         , veg_cs%gresp_xfer(p))
    mov(ystates(s_gresp_storage)      , veg_cs%gresp_storage(p))

    rfluxes(:) = 0._r8
    ar_p = veg_cf%leaf_curmr(p)             &
         + veg_cf%froot_curmr(p)            &
         + veg_cf%cpool_leaf_gr(p)          &
         + veg_cf%cpool_froot_gr(p)         &
         + veg_cf%cpool_leaf_storage_gr(p)  &
         + veg_cf%cpool_froot_storage_gr(p)
    if (woody(ivt(p)) == 1._r8) then
      ar_p = ar_p                                 &
         + veg_cf%livestem_curmr(p)             &
         + veg_cf%livecroot_curmr(p)            &
         + veg_cf%cpool_livestem_gr(p)          &
         + veg_cf%cpool_deadstem_gr(p)          &
         + veg_cf%cpool_livecroot_gr(p)         &
         + veg_cf%cpool_deadcroot_gr(p)         &
         + veg_cf%cpool_livestem_storage_gr(p)  &
         + veg_cf%cpool_deadstem_storage_gr(p)  &
         + veg_cf%cpool_livecroot_storage_gr(p) &
         + veg_cf%cpool_deadcroot_storage_gr(p)
    endif
    if (ivt(p) >= npcropmin) then
      ar_p= ar_p                                 &
        + veg_cf%livestem_curmr(p)             &
        + veg_cf%grain_curmr(p)                &
        + veg_cf%cpool_livestem_gr(p)          &
        + veg_cf%cpool_grain_gr(p)             &
        + veg_cf%cpool_livestem_storage_gr(p)  &
        + veg_cf%cpool_grain_storage_gr(p)
    endif
    !assemble reactive fluxes
    mov(rfluxes(f_cpool_to_leafc)               , veg_cf%cpool_to_leafc(p))
    mov(rfluxes(f_cpool_to_leafc_storage)       , veg_cf%cpool_to_leafc_storage(p))
    mov(rfluxes(f_cpool_to_frootc)              , veg_cf%cpool_to_frootc(p))
    mov(rfluxes(f_cpool_to_frootc_storage)      , veg_cf%cpool_to_frootc_storage(p))
    mov(rfluxes(f_cpool_to_xsmrpool)            , veg_cf%cpool_to_xsmrpool(p))
    mov(rfluxes(f_cpool_to_gresp_storage)       , veg_cf%cpool_to_gresp_storage(p))
    mov(rfluxes(f_cpool_to_ar)                  , ar_p)
    if (woody(ivt(p)) == 1._r8) then
      mov(rfluxes(f_cpool_to_livestemc)           , veg_cf%cpool_to_livestemc(p))
      mov(rfluxes(f_cpool_to_livestemc_storage)   , veg_cf%cpool_to_livestemc_storage(p))
      mov(rfluxes(f_cpool_to_deadstemc)           , veg_cf%cpool_to_deadstemc(p))
      mov(rfluxes(f_cpool_to_deadstemc_storage)   , veg_cf%cpool_to_deadstemc_storage(p))
      mov(rfluxes(f_cpool_to_livecrootc)          , veg_cf%cpool_to_livecrootc(p))
      mov(rfluxes(f_cpool_to_livecrootc_storage)  , veg_cf%cpool_to_livecrootc_storage(p))
      mov(rfluxes(f_cpool_to_deadcrootc)          , veg_cf%cpool_to_deadcrootc(p))
      mov(rfluxes(f_cpool_to_deadcrootc_storage)  , veg_cf%cpool_to_deadcrootc_storage(p))
      mov(rfluxes(f_livestemc_to_deadstemc)       , veg_cf%livestemc_to_deadstemc(p))
      mov(rfluxes(f_livestemc_xfer_to_livestemc)  , veg_cf%livestemc_xfer_to_livestemc(p))
      mov(rfluxes(f_livestemc_storage_to_xfer)    , veg_cf%livestemc_storage_to_xfer(p))
      mov(rfluxes(f_deadstemc_xfer_to_deadstemc)  , veg_cf%deadstemc_xfer_to_deadstemc(p))
      mov(rfluxes(f_deadstemc_storage_to_xfer)    , veg_cf%deadstemc_storage_to_xfer(p))
      mov(rfluxes(f_livecrootc_xfer_to_livecrootc), veg_cf%livecrootc_xfer_to_livecrootc(p))
      mov(rfluxes(f_livecrootc_storage_to_xfer)   , veg_cf%livecrootc_storage_to_xfer(p))
      mov(rfluxes(f_livecrootc_to_deadcrootc)     , veg_cf%livecrootc_to_deadcrootc(p))
      mov(rfluxes(f_deadcrootc_xfer_to_deadcrootc), veg_cf%deadcrootc_xfer_to_deadcrootc(p))
      mov(rfluxes(f_deadcrootc_storage_to_xfer)   , veg_cf%deadcrootc_storage_to_xfer(p))
      mov(rfluxes(f_gresp_storage_to_xfer)        , veg_cf%gresp_storage_to_xfer(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(rfluxes(f_cpool_to_livestemc)           , veg_cf%cpool_to_livestemc(p))
      mov(rfluxes(f_cpool_to_livestemc_storage)   , veg_cf%cpool_to_livestemc_storage(p))
      mov(rfluxes(f_cpool_to_grainc)              , veg_cf%cpool_to_grainc(p))
      mov(rfluxes(f_cpool_to_grainc_storage)      , veg_cf%cpool_to_grainc_storage(p))
      mov(rfluxes(f_livestemc_xfer_to_livestemc)  , veg_cf%livestemc_xfer_to_livestemc(p))
      mov(rfluxes(f_livestemc_storage_to_xfer)    , veg_cf%livestemc_storage_to_xfer(p))
      mov(rfluxes(f_livestemc_to_litter)          , veg_cf%livestemc_to_litter(p))
      mov(rfluxes(f_grainc_to_food)               , veg_cf%grainc_to_food(p))
      mov(rfluxes(f_grainc_xfer_to_grainc)        , veg_cf%grainc_xfer_to_grainc(p))
      mov(rfluxes(f_grainc_storage_to_xfer)       , veg_cf%grainc_storage_to_xfer(p))
    endif
    mov(rfluxes(f_leafc_to_litter)              , veg_cf%leafc_to_litter(p))
    mov(rfluxes(f_leafc_xfer_to_leafc)          , veg_cf%leafc_xfer_to_leafc(p))
    mov(rfluxes(f_leafc_storage_to_xfer)        , veg_cf%leafc_storage_to_xfer(p))
    mov(rfluxes(f_frootc_to_litter)             , veg_cf%frootc_to_litter(p))
    mov(rfluxes(f_frootc_xfer_to_frootc)        , veg_cf%frootc_xfer_to_frootc(p))
    mov(rfluxes(f_frootc_storage_to_xfer)       , veg_cf%frootc_storage_to_xfer(p))

    !obtain the limiting factor

    call flux_correction(spm_carbon_d%szrow, spm_carbon_d%szcol, spm_carbon_p, &
       spm_carbon_d, dt, ystates, rfluxes)

    !correct the fluxes

    imov(rfluxes(f_cpool_to_leafc)               , veg_cf%cpool_to_leafc(p))
    imov(rfluxes(f_cpool_to_leafc_storage)       , veg_cf%cpool_to_leafc_storage(p))
    imov(rfluxes(f_cpool_to_frootc)              , veg_cf%cpool_to_frootc(p))
    imov(rfluxes(f_cpool_to_frootc_storage)      , veg_cf%cpool_to_frootc_storage(p))
    imov(rfluxes(f_cpool_to_xsmrpool)            , veg_cf%cpool_to_xsmrpool(p))
    imov(rfluxes(f_cpool_to_gresp_storage)       , veg_cf%cpool_to_gresp_storage(p))

    if (woody(ivt(p)) == 1._r8) then
      imov(rfluxes(f_cpool_to_livestemc)           , veg_cf%cpool_to_livestemc(p))
      imov(rfluxes(f_cpool_to_livestemc_storage)   , veg_cf%cpool_to_livestemc_storage(p))
      imov(rfluxes(f_cpool_to_deadstemc)           , veg_cf%cpool_to_deadstemc(p))
      imov(rfluxes(f_cpool_to_deadstemc_storage)   , veg_cf%cpool_to_deadstemc_storage(p))
      imov(rfluxes(f_cpool_to_livecrootc)          , veg_cf%cpool_to_livecrootc(p))
      imov(rfluxes(f_cpool_to_livecrootc_storage)  , veg_cf%cpool_to_livecrootc_storage(p))
      imov(rfluxes(f_cpool_to_deadcrootc)          , veg_cf%cpool_to_deadcrootc(p))
      imov(rfluxes(f_cpool_to_deadcrootc_storage)  , veg_cf%cpool_to_deadcrootc_storage(p))
      imov(rfluxes(f_livestemc_to_deadstemc)       , veg_cf%livestemc_to_deadstemc(p))
      imov(rfluxes(f_livestemc_xfer_to_livestemc)  , veg_cf%livestemc_xfer_to_livestemc(p))
      imov(rfluxes(f_livestemc_storage_to_xfer)    , veg_cf%livestemc_storage_to_xfer(p))
      imov(rfluxes(f_deadstemc_xfer_to_deadstemc)  , veg_cf%deadstemc_xfer_to_deadstemc(p))
      imov(rfluxes(f_deadstemc_storage_to_xfer)    , veg_cf%deadstemc_storage_to_xfer(p))
      imov(rfluxes(f_livecrootc_xfer_to_livecrootc), veg_cf%livecrootc_xfer_to_livecrootc(p))
      imov(rfluxes(f_livecrootc_storage_to_xfer)   , veg_cf%livecrootc_storage_to_xfer(p))
      imov(rfluxes(f_livecrootc_to_deadcrootc)     , veg_cf%livecrootc_to_deadcrootc(p))
      imov(rfluxes(f_deadcrootc_xfer_to_deadcrootc), veg_cf%deadcrootc_xfer_to_deadcrootc(p))
      imov(rfluxes(f_deadcrootc_storage_to_xfer)   , veg_cf%deadcrootc_storage_to_xfer(p))
      imov(rfluxes(f_gresp_storage_to_xfer)        , veg_cf%gresp_storage_to_xfer(p))
    endif
    if (ivt(p) >= npcropmin) then
      imov(rfluxes(f_cpool_to_livestemc)           , veg_cf%cpool_to_livestemc(p))
      imov(rfluxes(f_cpool_to_livestemc_storage)   , veg_cf%cpool_to_livestemc_storage(p))
      imov(rfluxes(f_cpool_to_grainc)              , veg_cf%cpool_to_grainc(p))
      imov(rfluxes(f_cpool_to_grainc_storage)      , veg_cf%cpool_to_grainc_storage(p))
      imov(rfluxes(f_livestemc_xfer_to_livestemc)  , veg_cf%livestemc_xfer_to_livestemc(p))
      imov(rfluxes(f_livestemc_storage_to_xfer)    , veg_cf%livestemc_storage_to_xfer(p))
      imov(rfluxes(f_livestemc_to_litter)          , veg_cf%livestemc_to_litter(p))
      imov(rfluxes(f_grainc_to_food)               , veg_cf%grainc_to_food(p))
      imov(rfluxes(f_grainc_xfer_to_grainc)        , veg_cf%grainc_xfer_to_grainc(p))
      imov(rfluxes(f_grainc_storage_to_xfer)       , veg_cf%grainc_storage_to_xfer(p))
    endif
    imov(rfluxes(f_leafc_to_litter)              , veg_cf%leafc_to_litter(p))
    imov(rfluxes(f_leafc_xfer_to_leafc)          , veg_cf%leafc_xfer_to_leafc(p))
    imov(rfluxes(f_leafc_storage_to_xfer)        , veg_cf%leafc_storage_to_xfer(p))
    imov(rfluxes(f_frootc_to_litter)             , veg_cf%frootc_to_litter(p))
    imov(rfluxes(f_frootc_xfer_to_frootc)        , veg_cf%frootc_xfer_to_frootc(p))
    imov(rfluxes(f_frootc_storage_to_xfer)       , veg_cf%frootc_storage_to_xfer(p))

    if(rfluxes(f_cpool_to_ar)  <  ar_p)then
      rscal= rfluxes(f_cpool_to_ar)/ar_p
      ascal(veg_cf%leaf_curmr(p)                , rscal)
      ascal(veg_cf%froot_curmr(p)               , rscal)
      ascal(veg_cf%cpool_leaf_gr(p)             , rscal)
      ascal(veg_cf%cpool_froot_gr(p)            , rscal)
      ascal(veg_cf%cpool_leaf_storage_gr(p)     , rscal)
      ascal(veg_cf%cpool_froot_storage_gr(p)    , rscal)
      if (woody(ivt(p)) == 1._r8) then
        ascal(veg_cf%livestem_curmr(p)            , rscal)
        ascal(veg_cf%livecroot_curmr(p)           , rscal)
        ascal(veg_cf%cpool_livestem_gr(p)         , rscal)
        ascal(veg_cf%cpool_deadstem_gr(p)         , rscal)
        ascal(veg_cf%cpool_livecroot_gr(p)        , rscal)
        ascal(veg_cf%cpool_deadcroot_gr(p)        , rscal)
        ascal(veg_cf%cpool_livestem_storage_gr(p) , rscal)
        ascal(veg_cf%cpool_deadstem_storage_gr(p) , rscal)
        ascal(veg_cf%cpool_livecroot_storage_gr(p), rscal)
        ascal(veg_cf%cpool_deadcroot_storage_gr(p), rscal)
      endif
      if (ivt(p) >= npcropmin) then
        ascal(veg_cf%livestem_curmr(p)            , rscal)
        ascal(veg_cf%grain_curmr(p)               , rscal)
        ascal(veg_cf%cpool_livestem_gr(p)         , rscal)
        ascal(veg_cf%cpool_grain_gr(p)            , rscal)
        ascal(veg_cf%cpool_livestem_storage_gr(p) , rscal)
        ascal(veg_cf%cpool_grain_storage_gr(p)    , rscal)
      endif
    endif
  enddo
  end associate
  end subroutine carbon_flux_limiter
!---------------------------------------------------------------------------
  subroutine nitrogen_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_nf, veg_ns)

    use decompMod           , only : bounds_type
    use CNStateType               , only : cnstate_type
    use VegetationDataType, only : vegetation_nitrogen_flux
    use VegetationDataType, only : vegetation_nitrogen_state
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(vegetation_nitrogen_flux)  , intent(inout) :: veg_nf
    type(vegetation_nitrogen_state) , intent(inout) :: veg_ns

  integer :: fp, p
  real(r8) :: ystates(n_nutrient_states)
  real(r8) :: rfluxes(n_nutrient_fluxes)
  real(r8) :: dt
  associate(                                                       &
         ivt                   => veg_pp%itype                   , & ! Input:  [integer  (:)     ]  pft vegetation type
         woody                 => veg_vp%woody                   , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
         nf                    => veg_nf              , &
         ns                    => veg_ns               &
  )
  ! set time steps
  dt = real( get_step_size(), r8 )
  do fp = 1,num_soilp
    p = filter_soilp(fp)
    ystates(:) = 0._r8
    !assemble state variables
    mov(ystates(s_npool)               , veg_ns%npool(p))
    mov(ystates(s_leafn)               , veg_ns%leafn(p))
    mov(ystates(s_leafn_xfer)          , veg_ns%leafn_xfer(p))
    mov(ystates(s_leafn_storage)       , veg_nf%npool_to_leafn_storage(p))
    mov(ystates(s_frootn)              , veg_ns%frootn(p))
    mov(ystates(s_frootn_xfer)         , veg_ns%frootn_xfer(p))
    mov(ystates(s_frootn_sotrage)      , veg_ns%frootn_storage(p))
    if (woody(ivt(p)) == 1.0_r8) then
      mov(ystates(s_livestemn)           , veg_ns%livestemn(p))
      mov(ystates(s_livestemn_xfer)      , veg_ns%livestemn_xfer(p))
      mov(ystates(s_livestemn_storage)   , veg_ns%livestemn_storage(p))
      mov(ystates(s_deadstemn)           , veg_ns%deadstemn(p))
      mov(ystates(s_deadstemn_xfer)      , veg_ns%deadstemn_xfer(p))
      mov(ystates(s_deadstemn_storage)   , veg_ns%deadstemn_storage(p))
      mov(ystates(s_livecrootn)          , veg_ns%livecrootn(p))
      mov(ystates(s_livecrootn_xfer)     , veg_ns%livecrootn_xfer(p))
      mov(ystates(s_livecrootn_storage)  , veg_ns%livecrootn_storage(p))
      mov(ystates(s_deadcrootn)          , veg_ns%deadcrootn(p))
      mov(ystates(s_deadcrootn_xfer)     , veg_ns%deadcrootn_xfer(p))
      mov(ystates(s_deadcrootn_storage)  , veg_ns%deadcrootn_storage(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(ystates(s_grainn)              , veg_ns%grainn(p))
      mov(ystates(s_grainn_xfer)         , veg_ns%grainn_xfer(p))
      mov(ystates(s_grainn_storage)      , veg_ns%grainn_storage(p))
      mov(ystates(s_livestemn)           , veg_ns%livestemn(p))
      mov(ystates(s_livestemn_xfer)      , veg_ns%livestemn_xfer(p))
      mov(ystates(s_livestemn_storage)   , veg_ns%livestemn_storage(p))
    endif
    mov(ystates(s_retransn)            , veg_ns%retransn(p))

    rfluxes(:) = 0._r8
    !assemble reactive fluxes
    mov(rfluxes(f_npool_to_leafn)                , veg_nf%npool_to_leafn(p))
    mov(rfluxes(f_npool_to_leafn_storage)        , veg_nf%npool_to_leafn_storage(p))
    mov(rfluxes(f_npool_to_frootn)               , veg_nf%npool_to_frootn(p))
    mov(rfluxes(f_npool_to_frootn_storage)       , veg_nf%npool_to_frootn_storage(p))
    if (woody(ivt(p)) == 1.0_r8) then
      mov(rfluxes(f_npool_to_livestemn)            , veg_nf%npool_to_livestemn(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , veg_nf%npool_to_livestemn_storage(p))
      mov(rfluxes(f_npool_to_livecrootn)           , veg_nf%npool_to_livecrootn(p))
      mov(rfluxes(f_npool_to_livecrootn_storage)   , veg_nf%npool_to_livecrootn_storage(p))
      mov(rfluxes(f_npool_to_deadstemn)            , veg_nf%npool_to_deadstemn(p))
      mov(rfluxes(f_npool_to_deadcrootn)           , veg_nf%npool_to_deadcrootn(p))
      mov(rfluxes(f_npool_to_deadstemn_storage)    , veg_nf%npool_to_deadstemn_storage(p))
      mov(rfluxes(f_npool_to_deadcrootn_storage)   , veg_nf%npool_to_deadcrootn_storage(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , veg_nf%livestemn_storage_to_xfer(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_nf%livestemn_xfer_to_livestemn(p))
      mov(rfluxes(f_livestemn_to_retransn)         , veg_nf%livestemn_to_retransn(p))
      mov(rfluxes(f_livestemn_to_deadstemn)        , veg_nf%livestemn_to_deadstemn(p))
      mov(rfluxes(f_livecrootn_to_deadcrootn)      , veg_nf%livecrootn_to_deadcrootn(p))
      mov(rfluxes(f_livecrootn_to_retransn)        , veg_nf%livecrootn_to_retransn(p))
      mov(rfluxes(f_livecrootn_storage_to_xfer)    , veg_nf%livecrootn_storage_to_xfer(p))
      mov(rfluxes(f_livecrootn_xfer_to_livecrootn) , veg_nf%livecrootn_xfer_to_livecrootn(p))
      mov(rfluxes(f_deadstemn_storage_to_xfer)     , veg_nf%deadstemn_storage_to_xfer(p))
      mov(rfluxes(f_deadstemn_xfer_to_deadstem)    , veg_nf%deadstemn_xfer_to_deadstemn(p))
      mov(rfluxes(f_deadcrootn_storage_to_xfer)    , veg_nf%deadcrootn_storage_to_xfer(p))
      mov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , veg_nf%deadcrootn_xfer_to_deadcrootn(p))
    endif

    if (ivt(p) >= npcropmin) then
      mov(rfluxes(f_npool_to_livestemn)            , veg_nf%npool_to_livestemn(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , veg_nf%npool_to_livestemn_storage(p))
      mov(rfluxes(f_npool_to_grainn)               , veg_nf%npool_to_grainn(p))
      mov(rfluxes(f_npool_to_grainn_storage)       , veg_nf%npool_to_grainn_storage(p))
      mov(rfluxes(f_frootn_to_retransn)            , veg_nf%frootn_to_retransn(p))
      mov(rfluxes(f_livestemn_to_litter)           , veg_nf%livestemn_to_litter(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , veg_nf%livestemn_storage_to_xfer(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   ,  veg_nf%livestemn_xfer_to_livestemn(p))
      mov(rfluxes(f_livestemn_to_retransn)         , veg_nf%livestemn_to_retransn(p))
      mov(rfluxes(f_grainn_xfer_to_grainn)         , veg_nf%grainn_xfer_to_grainn(p))
      mov(rfluxes(f_grainn_to_food)                , veg_nf%grainn_to_food(p))
    endif

    mov(rfluxes(f_leafn_to_litter)               , veg_nf%leafn_to_litter(p))
    mov(rfluxes(f_leafn_xfer_to_leafn)           , veg_nf%leafn_xfer_to_leafn(p))
    mov(rfluxes(f_leafn_storage_to_xfer)         , veg_nf%leafn_storage_to_xfer(p))
    mov(rfluxes(f_frootn_xfer_to_frootn)         , veg_nf%frootn_xfer_to_frootn(p))
    mov(rfluxes(f_frootn_storage_to_xfer)        , veg_nf%frootn_storage_to_xfer(p))
    mov(rfluxes(f_leafn_to_retransn)             , veg_nf%leafn_to_retransn(p))
    mov(rfluxes(f_frootn_to_litter)              , veg_nf%frootn_to_litter(p))

    mov(rfluxes(f_retransn_to_npool)             , veg_nf%retransn_to_npool(p))
    mov(rfluxes(f_supplement_to_plantn)          , veg_nf%supplement_to_plantn(p))

    call flux_correction(spm_nutrient_d%szrow, spm_nutrient_d%szcol, spm_nutrient_p, &
      spm_nutrient_d, dt, ystates, rfluxes)

    !correct the fluxes
    imov(rfluxes(f_npool_to_leafn)                , veg_nf%npool_to_leafn(p))
    imov(rfluxes(f_npool_to_leafn_storage)        , veg_nf%npool_to_leafn_storage(p))
    imov(rfluxes(f_npool_to_frootn)               , veg_nf%npool_to_frootn(p))
    imov(rfluxes(f_npool_to_frootn_storage)       , veg_nf%npool_to_frootn_storage(p))
    if (woody(ivt(p)) == 1.0_r8) then
      imov(rfluxes(f_npool_to_livestemn)            , veg_nf%npool_to_livestemn(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , veg_nf%npool_to_livestemn_storage(p))
      imov(rfluxes(f_npool_to_livecrootn)           , veg_nf%npool_to_livecrootn(p))
      imov(rfluxes(f_npool_to_livecrootn_storage)   , veg_nf%npool_to_livecrootn_storage(p))
      imov(rfluxes(f_npool_to_deadstemn)            , veg_nf%npool_to_deadstemn(p))
      imov(rfluxes(f_npool_to_deadcrootn)           , veg_nf%npool_to_deadcrootn(p))
      imov(rfluxes(f_npool_to_deadstemn_storage)    , veg_nf%npool_to_deadstemn_storage(p))
      imov(rfluxes(f_npool_to_deadcrootn_storage)   , veg_nf%npool_to_deadcrootn_storage(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , veg_nf%livestemn_storage_to_xfer(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_nf%livestemn_xfer_to_livestemn(p))
      imov(rfluxes(f_livestemn_to_retransn)         , veg_nf%livestemn_to_retransn(p))
      imov(rfluxes(f_livestemn_to_deadstemn)        , veg_nf%livestemn_to_deadstemn(p))
      imov(rfluxes(f_livecrootn_to_deadcrootn)      , veg_nf%livecrootn_to_deadcrootn(p))
      imov(rfluxes(f_livecrootn_to_retransn)        , veg_nf%livecrootn_to_retransn(p))
      imov(rfluxes(f_livecrootn_storage_to_xfer)    , veg_nf%livecrootn_storage_to_xfer(p))
      imov(rfluxes(f_livecrootn_xfer_to_livecrootn) , veg_nf%livecrootn_xfer_to_livecrootn(p))
      imov(rfluxes(f_deadstemn_storage_to_xfer)     , veg_nf%deadstemn_storage_to_xfer(p))
      imov(rfluxes(f_deadstemn_xfer_to_deadstem)    , veg_nf%deadstemn_xfer_to_deadstemn(p))
      imov(rfluxes(f_deadcrootn_storage_to_xfer)    , veg_nf%deadcrootn_storage_to_xfer(p))
      imov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , veg_nf%deadcrootn_xfer_to_deadcrootn(p))
    endif

    if (ivt(p) >= npcropmin) then
      imov(rfluxes(f_npool_to_livestemn)            , veg_nf%npool_to_livestemn(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , veg_nf%npool_to_livestemn_storage(p))
      imov(rfluxes(f_npool_to_grainn)               , veg_nf%npool_to_grainn(p))
      imov(rfluxes(f_npool_to_grainn_storage)       , veg_nf%npool_to_grainn_storage(p))
      imov(rfluxes(f_frootn_to_retransn)            , veg_nf%frootn_to_retransn(p))
      imov(rfluxes(f_livestemn_to_litter)           , veg_nf%livestemn_to_litter(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , veg_nf%livestemn_storage_to_xfer(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   ,  veg_nf%livestemn_xfer_to_livestemn(p))
      imov(rfluxes(f_livestemn_to_retransn)         , veg_nf%livestemn_to_retransn(p))
      imov(rfluxes(f_grainn_xfer_to_grainn)         , veg_nf%grainn_xfer_to_grainn(p))
      imov(rfluxes(f_grainn_to_food)                , veg_nf%grainn_to_food(p))
    endif

    mov(rfluxes(f_leafn_to_litter)               , veg_nf%leafn_to_litter(p))
    mov(rfluxes(f_leafn_xfer_to_leafn)           , veg_nf%leafn_xfer_to_leafn(p))
    mov(rfluxes(f_leafn_storage_to_xfer)         , veg_nf%leafn_storage_to_xfer(p))
    mov(rfluxes(f_frootn_xfer_to_frootn)         , veg_nf%frootn_xfer_to_frootn(p))
    mov(rfluxes(f_frootn_storage_to_xfer)        , veg_nf%frootn_storage_to_xfer(p))
    mov(rfluxes(f_leafn_to_retransn)             , veg_nf%leafn_to_retransn(p))
    mov(rfluxes(f_frootn_to_litter)              , veg_nf%frootn_to_litter(p))

    mov(rfluxes(f_retransn_to_npool)             , veg_nf%retransn_to_npool(p))
    mov(rfluxes(f_supplement_to_plantn)          , veg_nf%supplement_to_plantn(p))

  enddo
  end associate
  end subroutine nitrogen_flux_limiter

  !-------------------------------------------------------------------------------
  subroutine phosphorus_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_pf, veg_ps)

    use decompMod           , only : bounds_type
    use CNStateType               , only : cnstate_type
    use VegetationDataType, only : vegetation_phosphorus_flux
    use VegetationDataType, only : vegetation_phosphorus_state
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(vegetation_phosphorus_flux)  , intent(inout) :: veg_pf
    type(vegetation_phosphorus_state) , intent(inout) :: veg_ps

  integer :: fp, p
  real(r8) :: ystates(n_nutrient_states)
  real(r8) :: rfluxes(n_nutrient_fluxes)

  real(r8) :: dt

  associate(                                                     &
    ivt                   => veg_pp%itype                      , & ! Input:  [integer  (:)     ]  pft vegetation type
    woody                 => veg_vp%woody                      , & ! Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
    pf                    => veg_pf               , &
    ps                    => veg_ps                &
  )

  ! set time steps
  dt = real( get_step_size(), r8 )
  do fp = 1,num_soilp
    p = filter_soilp(fp)
    ystates(:) = 0._r8
    !assemble state variables
    mov(ystates(s_npool)               , veg_ps%ppool(p))
    mov(ystates(s_leafn)               , veg_ps%leafp(p))
    mov(ystates(s_leafn_xfer)          , veg_ps%leafp_xfer(p))
    mov(ystates(s_leafn_storage)       , veg_ps%leafp_storage(p))
    mov(ystates(s_frootn)              , veg_ps%frootp(p))
    mov(ystates(s_frootn_xfer)         , veg_ps%frootp_xfer(p))
    mov(ystates(s_frootn_sotrage)      , veg_pf%ppool_to_frootp_storage(p))
    if (woody(ivt(p)) == 1.0_r8) then
      mov(ystates(s_livestemn)           , veg_ps%livestemp(p))
      mov(ystates(s_livestemn_xfer)      , veg_ps%livestemp_xfer(p))
      mov(ystates(s_livestemn_storage)   , veg_ps%livestemp_storage(p))
      mov(ystates(s_deadstemn)           , veg_ps%deadstemp(p))
      mov(ystates(s_deadstemn_xfer)      , veg_ps%deadstemp_xfer(p))
      mov(ystates(s_deadstemn_storage)   , veg_ps%deadstemp_storage(p))
      mov(ystates(s_livecrootn)          , veg_ps%livecrootp(p))
      mov(ystates(s_livecrootn_xfer)     , veg_ps%livecrootp_xfer(p))
      mov(ystates(s_livecrootn_storage)  , veg_ps%livecrootp_storage(p))
      mov(ystates(s_deadcrootn)          , veg_ps%deadcrootp(p))
      mov(ystates(s_deadcrootn_xfer)     , veg_ps%deadcrootp_xfer(p))
      mov(ystates(s_deadcrootn_storage)  , veg_ps%deadcrootp_storage(p))
    endif
    if (ivt(p) >= npcropmin) then
      mov(ystates(s_grainn)              , veg_ps%grainp(p))
      mov(ystates(s_grainn_xfer)         , veg_ps%grainp_xfer(p))
      mov(ystates(s_grainn_storage)      , veg_ps%grainp_storage(p))
      mov(ystates(s_livestemn)           , veg_ps%livestemp(p))
      mov(ystates(s_livestemn_xfer)      , veg_ps%livestemp_xfer(p))
      mov(ystates(s_livestemn_storage)   , veg_ps%livestemp_storage(p))
    endif
    mov(ystates(s_retransn)              , veg_ps%retransp(p))

    rfluxes(:) = 0._r8
    !assemble reactive fluxes
    mov(rfluxes(f_npool_to_leafn)                , veg_pf%ppool_to_leafp(p))
    mov(rfluxes(f_npool_to_leafn_storage)        , veg_pf%ppool_to_leafp_storage(p))
    mov(rfluxes(f_npool_to_frootn)               , veg_pf%ppool_to_frootp(p))
    mov(rfluxes(f_npool_to_frootn_storage)       , veg_pf%ppool_to_frootp_storage(p))
    if (woody(ivt(p)) == 1._r8) then
      mov(rfluxes(f_npool_to_livestemn)            , veg_pf%ppool_to_livestemp(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , veg_pf%ppool_to_livestemp_storage(p))
      mov(rfluxes(f_npool_to_livecrootn)           , veg_pf%ppool_to_livecrootp(p))
      mov(rfluxes(f_npool_to_livecrootn_storage)   , veg_pf%ppool_to_livecrootp_storage(p))
      mov(rfluxes(f_npool_to_deadstemn)            , veg_pf%ppool_to_deadstemp(p))
      mov(rfluxes(f_npool_to_deadcrootn)           , veg_pf%ppool_to_deadcrootp(p))
      mov(rfluxes(f_npool_to_deadstemn_storage)    , veg_pf%ppool_to_deadstemp_storage(p))
      mov(rfluxes(f_npool_to_deadcrootn_storage)   , veg_pf%ppool_to_deadcrootp_storage(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , veg_pf%livestemp_storage_to_xfer(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_pf%livestemp_xfer_to_livestemp(p))
      mov(rfluxes(f_livestemn_to_retransn)         , veg_pf%livestemp_to_retransp(p))
      mov(rfluxes(f_livestemn_to_deadstemn)        , veg_pf%livestemp_to_deadstemp(p))
      mov(rfluxes(f_livecrootn_to_deadcrootn)      , veg_pf%livecrootp_to_deadcrootp(p))
      mov(rfluxes(f_livecrootn_to_retransn)        , veg_pf%livecrootp_to_retransp(p))
      mov(rfluxes(f_livecrootn_storage_to_xfer)    , veg_pf%livecrootp_storage_to_xfer(p))
      mov(rfluxes(f_livecrootn_xfer_to_livecrootn) , veg_pf%livecrootp_xfer_to_livecrootp(p))
      mov(rfluxes(f_deadstemn_storage_to_xfer)     , veg_pf%deadstemp_storage_to_xfer(p))
      mov(rfluxes(f_deadstemn_xfer_to_deadstem)    , veg_pf%deadstemp_xfer_to_deadstemp(p))
      mov(rfluxes(f_deadcrootn_storage_to_xfer)    , veg_pf%deadcrootp_storage_to_xfer(p))
      mov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , veg_pf%deadcrootp_xfer_to_deadcrootp(p))
    endif

    if (ivt(p) >= npcropmin) then
      mov(rfluxes(f_npool_to_livestemn)            , veg_pf%ppool_to_livestemp(p))
      mov(rfluxes(f_npool_to_livestemn_storage)    , veg_pf%ppool_to_livestemp_storage(p))
      mov(rfluxes(f_npool_to_grainn)               , veg_pf%ppool_to_grainp(p))
      mov(rfluxes(f_npool_to_grainn_storage)       , veg_pf%ppool_to_grainp_storage(p))
      mov(rfluxes(f_frootn_to_retransn)            , veg_pf%frootp_to_retransp(p))
      mov(rfluxes(f_livestemn_to_litter)           , veg_pf%livestemp_to_litter(p))
      mov(rfluxes(f_livestemn_storage_to_xfer)     , veg_pf%livestemp_storage_to_xfer(p))
      mov(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_pf%livestemp_xfer_to_livestemp(p))
      mov(rfluxes(f_livestemn_to_retransn)         , veg_pf%livestemp_to_retransp(p))
      mov(rfluxes(f_grainn_xfer_to_grainn)         , veg_pf%grainp_xfer_to_grainp(p))
      mov(rfluxes(f_grainn_to_food)                , veg_pf%grainp_to_food(p))
    endif

    mov(rfluxes(f_leafn_to_litter)               , veg_pf%leafp_to_litter(p))
    mov(rfluxes(f_leafn_xfer_to_leafn)           , veg_pf%leafp_xfer_to_leafp(p))
    mov(rfluxes(f_leafn_storage_to_xfer)         , veg_pf%ppool_to_leafp_storage(p))
    mov(rfluxes(f_frootn_xfer_to_frootn)         , veg_pf%frootp_xfer_to_frootp(p))
    mov(rfluxes(f_frootn_storage_to_xfer)        , veg_pf%frootp_storage_to_xfer(p))
    mov(rfluxes(f_leafn_to_retransn)             , veg_pf%leafp_to_retransp(p))
    mov(rfluxes(f_frootn_to_litter)              , veg_pf%frootp_to_litter(p))
    mov(rfluxes(f_retransn_to_npool)             , veg_pf%retransp_to_ppool(p))
    mov(rfluxes(f_supplement_to_plantn)          , veg_pf%supplement_to_plantp(p))
    !assemble stoichiometry matrix

    !obtain the limiting factor
    call flux_correction(spm_nutrient_d%szrow, spm_nutrient_d%szcol, spm_nutrient_p, &
      spm_nutrient_d, dt, ystates, rfluxes)
    !correct the fluxes
    imov(rfluxes(f_npool_to_leafn)                , veg_pf%ppool_to_leafp(p))
    imov(rfluxes(f_npool_to_leafn_storage)        , veg_pf%ppool_to_leafp_storage(p))
    imov(rfluxes(f_npool_to_frootn)               , veg_pf%ppool_to_frootp(p))
    imov(rfluxes(f_npool_to_frootn_storage)       , veg_pf%ppool_to_frootp_storage(p))
    if (woody(ivt(p)) == 1._r8) then
      imov(rfluxes(f_npool_to_livestemn)            , veg_pf%ppool_to_livestemp(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , veg_pf%ppool_to_livestemp_storage(p))
      imov(rfluxes(f_npool_to_livecrootn)           , veg_pf%ppool_to_livecrootp(p))
      imov(rfluxes(f_npool_to_livecrootn_storage)   , veg_pf%ppool_to_livecrootp_storage(p))
      imov(rfluxes(f_npool_to_deadstemn)            , veg_pf%ppool_to_deadstemp(p))
      imov(rfluxes(f_npool_to_deadcrootn)           , veg_pf%ppool_to_deadcrootp(p))
      imov(rfluxes(f_npool_to_deadstemn_storage)    , veg_pf%ppool_to_deadstemp_storage(p))
      imov(rfluxes(f_npool_to_deadcrootn_storage)   , veg_pf%ppool_to_deadcrootp_storage(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , veg_pf%livestemp_storage_to_xfer(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_pf%livestemp_xfer_to_livestemp(p))
      imov(rfluxes(f_livestemn_to_retransn)         , veg_pf%livestemp_to_retransp(p))
      imov(rfluxes(f_livestemn_to_deadstemn)        , veg_pf%livestemp_to_deadstemp(p))
      imov(rfluxes(f_livecrootn_to_deadcrootn)      , veg_pf%livecrootp_to_deadcrootp(p))
      imov(rfluxes(f_livecrootn_to_retransn)        , veg_pf%livecrootp_to_retransp(p))
      imov(rfluxes(f_livecrootn_storage_to_xfer)    , veg_pf%livecrootp_storage_to_xfer(p))
      imov(rfluxes(f_livecrootn_xfer_to_livecrootn) , veg_pf%livecrootp_xfer_to_livecrootp(p))
      imov(rfluxes(f_deadstemn_storage_to_xfer)     , veg_pf%deadstemp_storage_to_xfer(p))
      imov(rfluxes(f_deadstemn_xfer_to_deadstem)    , veg_pf%deadstemp_xfer_to_deadstemp(p))
      imov(rfluxes(f_deadcrootn_storage_to_xfer)    , veg_pf%deadcrootp_storage_to_xfer(p))
      imov(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , veg_pf%deadcrootp_xfer_to_deadcrootp(p))
    endif

    if (ivt(p) >= npcropmin) then
      imov(rfluxes(f_npool_to_livestemn)            , veg_pf%ppool_to_livestemp(p))
      imov(rfluxes(f_npool_to_livestemn_storage)    , veg_pf%ppool_to_livestemp_storage(p))
      imov(rfluxes(f_npool_to_grainn)               , veg_pf%ppool_to_grainp(p))
      imov(rfluxes(f_npool_to_grainn_storage)       , veg_pf%ppool_to_grainp_storage(p))
      imov(rfluxes(f_frootn_to_retransn)            , veg_pf%frootp_to_retransp(p))
      imov(rfluxes(f_livestemn_to_litter)           , veg_pf%livestemp_to_litter(p))
      imov(rfluxes(f_livestemn_storage_to_xfer)     , veg_pf%livestemp_storage_to_xfer(p))
      imov(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_pf%livestemp_xfer_to_livestemp(p))
      imov(rfluxes(f_livestemn_to_retransn)         , veg_pf%livestemp_to_retransp(p))
      imov(rfluxes(f_grainn_xfer_to_grainn)         , veg_pf%grainp_xfer_to_grainp(p))
      imov(rfluxes(f_grainn_to_food)                , veg_pf%grainp_to_food(p))
    endif

    imov(rfluxes(f_leafn_to_litter)               , veg_pf%leafp_to_litter(p))
    imov(rfluxes(f_leafn_xfer_to_leafn)           , veg_pf%leafp_xfer_to_leafp(p))
    imov(rfluxes(f_leafn_storage_to_xfer)         , veg_pf%ppool_to_leafp_storage(p))
    imov(rfluxes(f_frootn_xfer_to_frootn)         , veg_pf%frootp_xfer_to_frootp(p))
    imov(rfluxes(f_frootn_storage_to_xfer)        , veg_pf%frootp_storage_to_xfer(p))
    imov(rfluxes(f_leafn_to_retransn)             , veg_pf%leafp_to_retransp(p))
    imov(rfluxes(f_frootn_to_litter)              , veg_pf%frootp_to_litter(p))
    imov(rfluxes(f_retransn_to_npool)             , veg_pf%retransp_to_ppool(p))
    imov(rfluxes(f_supplement_to_plantn)          , veg_pf%supplement_to_plantp(p))

  enddo
  end associate
  end subroutine phosphorus_flux_limiter

end module PhenologyFLuxLimitMod
