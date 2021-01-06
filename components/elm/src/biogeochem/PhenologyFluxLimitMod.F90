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
  use elm_varctl                  , only : iulog
  use abortutils                  , only : endrun
implicit none
  private

  !carbon state indices
  integer :: s_cpool
  integer :: s_leafc
  integer :: s_leafc_xfer
  integer :: s_leafc_storage
  integer :: s_frootc
  integer :: s_frootc_xfer
  integer :: s_frootc_storage
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
!---------------------------------------------
  !carbon flux indices
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
!---------------------------------------------
  !nutrient state indices (N & P use the same indices)
  integer :: s_npool
  integer :: s_leafn
  integer :: s_leafn_xfer
  integer :: s_leafn_storage
  integer :: s_frootn
  integer :: s_frootn_xfer
  integer :: s_frootn_storage
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
!---------------------------------------------
  !nutrient flux indices (N & P use the same indices)
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
  integer :: num_carbon_fluxes
  integer :: num_nutrient_fluxes
  integer :: num_carbon_states
  integer :: num_nutrient_states

  class(sparseMat_type), pointer :: spm_carbon_p,spm_carbon_d
  class(sparseMat_type), pointer :: spm_nutrient_p,spm_nutrient_d
  type(spm_list_type), pointer :: spm_list
  public :: InitPhenoFluxLimiter
  public :: phenology_flux_limiter
contains

  subroutine ascal(a, b)

  !DESCRIPTION
  !a = a * b
  implicit none
  real(r8), intent(inout) :: a
  real(r8), intent(in) :: b

  a = a * b
  end subroutine ascal
!------------------------------------------------------------
  subroutine fpmax(a,b)
  !
  ! DESCRIPTION
  ! b=max(a,0._r8)
  implicit none
  real(r8), intent(in) :: a
  real(r8), intent(out):: b

  b=max(a,0._r8)

  end subroutine fpmax
!------------------------------------------------------------
  function ic_next(id)result(ans)
  !
  ! DESCRIPTION
  ! add one to id
  implicit none
  integer, intent(inout) :: id
  integer :: ans
  id=id+1
  ans=id
  end function ic_next
  !-------------------------------------------------------------------------------
  subroutine init_phenofluxl_counters()

  !
  ! DESCRIPTION
  ! initialize the array member location counters
  implicit none
  integer :: id

  !carbon state variables
  id = 0
  s_cpool              = ic_next(id)
  s_leafc              = ic_next(id)
  s_leafc_xfer         = ic_next(id)
  s_leafc_storage      = ic_next(id)
  s_frootc             = ic_next(id)
  s_frootc_xfer        = ic_next(id)
  s_frootc_storage     = ic_next(id)
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
  num_carbon_states = id

  !carbon flux variables
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
  num_carbon_fluxes = id

  !nutrient state variables
  id=0
  s_npool               = ic_next(id)
  s_leafn               = ic_next(id)
  s_leafn_xfer          = ic_next(id)
  s_leafn_storage       = ic_next(id)
  s_frootn              = ic_next(id)
  s_frootn_xfer         = ic_next(id)
  s_frootn_storage      = ic_next(id)
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
  num_nutrient_states = id

  !nutrient flux variables
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
  num_nutrient_fluxes = id
  end subroutine init_phenofluxl_counters
  !-------------------------------------------------------------------------------
  subroutine InitPhenoFluxLimiter()

  !
  ! DESCRIPTION
  ! initialize the phenology flux limiter
  use LSparseMatMod, only : spm_list_type,  spm_list_to_mat
  use LSparseMatMod, only : spm_list_init, spm_list_insert
  implicit none
  integer :: nelms
  type(spm_list_type), pointer :: spm_list

  call init_phenofluxl_counters

  !initialize stoichiometric relationship between carbon consumption flux and state varaibles
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
  call spm_list_insert(spm_list, -1._r8, f_frootc_storage_to_xfer       , s_frootc_storage, nelms)
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
  !turn the list into sparse matrix form
  call spm_list_to_mat(spm_list, spm_carbon_d, nelms, f_gresp_storage_to_xfer)
  
  !initialize stoichiometric relationship between carbon production flux and corresponding state varaibles
  call spm_list_init(spm_list, 1._r8, f_cpool_to_leafc             , s_leafc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafc_xfer_to_leafc        , s_leafc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafc_storage_to_xfer      , s_leafc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_leafc_storage     , s_leafc_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_frootc            , s_frootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootc_xfer_to_frootc      , s_frootc, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootc_storage_to_xfer     , s_frootc_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_cpool_to_frootc_storage    , s_frootc_storage, nelms)
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

   
  !turn the list into sparse matrix form
  call spm_list_to_mat(spm_list, spm_carbon_p, nelms, f_gresp_storage_to_xfer)
  !initialize stoichiometry relationship between nutrient consumption and corresponding state variables
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
  call spm_list_insert(spm_list, -1._r8, f_frootn_storage_to_xfer       , s_frootn_storage, nelms)
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
  !turn the list into sparse matrix form  
  call spm_list_to_mat(spm_list, spm_nutrient_d, nelms, f_supplement_to_plantn)
  !initialize stoichiometry relationship between nutrient production and corresponding state variables
  call spm_list_init(spm_list, 1._r8, f_retransn_to_npool                , s_npool, nelms)
  call spm_list_insert(spm_list, 1._r8, f_supplement_to_plantn           , s_npool, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_leafn                 , s_leafn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafn_xfer_to_leafn            , s_leafn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_leafn_storage_to_xfer          , s_leafn_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_leafn_storage         , s_leafn_storage, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_frootn                , s_frootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootn_xfer_to_frootn          , s_frootn, nelms)
  call spm_list_insert(spm_list, 1._r8, f_frootn_storage_to_xfer         , s_frootn_xfer, nelms)
  call spm_list_insert(spm_list, 1._r8, f_npool_to_frootn_storage        , s_frootn_storage, nelms)
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
  !turn the list into sparse matrix form  
  call spm_list_to_mat(spm_list, spm_nutrient_p, nelms, f_supplement_to_plantn)
  end subroutine InitPhenoFluxLimiter
!---------------------------------------------------------------------------
  subroutine phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
      veg_cf, veg_cs , c13_veg_cf, c13_veg_cs , c14_veg_cf, c14_veg_cs , &
      veg_nf, veg_ns, veg_pf, veg_ps)

    !
    ! DESCRIPTION
    !  apply the phenology flux limiter to avoid potential negative fluxes.
    use decompMod           , only : bounds_type
    use CropType                  , only : crop_type
    use CNStateType               , only : cnstate_type
    use VegetationDataType, only : vegetation_carbon_flux
    use VegetationDataType, only : vegetation_carbon_state
    use VegetationDataType, only : vegetation_nitrogen_flux
    use VegetationDataType, only : vegetation_nitrogen_state
    use VegetationDataType, only : vegetation_phosphorus_flux
    use VegetationDataType, only : vegetation_phosphorus_state
    use elm_varctl          , only : use_c13, use_c14
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

  !apply the carbon flux limiter
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

  !apply the nitrogen flux limiter
  call nitrogen_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_nf, veg_ns)

  !apply the phoshporus flux limiter
  call phosphorus_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_pf, veg_ps)


  end subroutine phenology_flux_limiter
  !-------------------------------------------------------------------------------

  subroutine carbon_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, crop_vars,  &
      veg_cf, veg_cs )
    !
    ! DESCRIPTION
    ! the flux limiter for phenology carbon fluxes
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
  real(r8) :: ystates(num_carbon_states)
  real(r8) :: rfluxes(num_carbon_fluxes)
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
    ystates(s_cpool)              = veg_cs%cpool(p)
    ystates(s_leafc)              = veg_cs%leafc(p)
    ystates(s_leafc_xfer)         = veg_cs%leafc_xfer(p)
    ystates(s_leafc_storage)      = veg_cs%leafc_storage(p)
    ystates(s_frootc)             = veg_cs%frootc(p)
    ystates(s_frootc_xfer)        = veg_cs%frootc_xfer(p)
    ystates(s_frootc_storage)     = veg_cs%frootc_storage(p)
    if (woody(ivt(p)) == 1._r8) then
      ystates(s_livestemc)          = veg_cs%livestemc(p)
      ystates(s_livestemc_xfer)     = veg_cs%livestemc_xfer(p)
      ystates(s_livestemc_storage)  = veg_cs%livestemc_storage(p)
      ystates(s_livecrootc)         = veg_cs%livecrootc(p)
      ystates(s_livecrootc_xfer)    = veg_cs%livecrootc_xfer(p)
      ystates(s_livecrootc_storage) = veg_cs%livecrootc_storage(p)
      ystates(s_deadstemc)          = veg_cs%deadstemc(p)
      ystates(s_deadstemc_xfer)     = veg_cs%deadstemc_xfer(p)
      ystates(s_deadstemc_storage)  = veg_cs%deadstemc_storage(p)
      ystates(s_deadcrootc)         = veg_cs%deadcrootc(p)
      ystates(s_deadcrootc_xfer)    = veg_cs%deadcrootc_xfer(p)
      ystates(s_deadcrootc_storage) = veg_cs%deadcrootc_storage(p)
    endif
    if (ivt(p) >= npcropmin) then
      ystates(s_livestemc)          = veg_cs%livestemc(p)
      ystates(s_livestemc_xfer)     = veg_cs%livestemc_xfer(p)
      ystates(s_livestemc_storage)  = veg_cs%livestemc_storage(p)
      ystates(s_grainc)             = veg_cs%grainc(p)
      ystates(s_grainc_xfer)        = veg_cs%grainc_xfer(p)
      ystates(s_grainc_storage)     = veg_cs%grainc_storage(p)
    endif
    ystates(s_gresp_xfer)         = veg_cs%gresp_xfer(p)
    ystates(s_gresp_storage)      = veg_cs%gresp_storage(p)

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
    rfluxes(f_cpool_to_leafc)               = veg_cf%cpool_to_leafc(p)
    rfluxes(f_cpool_to_leafc_storage)       = veg_cf%cpool_to_leafc_storage(p)
    rfluxes(f_cpool_to_frootc)              = veg_cf%cpool_to_frootc(p)
    rfluxes(f_cpool_to_frootc_storage)      = veg_cf%cpool_to_frootc_storage(p)
    rfluxes(f_cpool_to_xsmrpool)            = veg_cf%cpool_to_xsmrpool(p)
    rfluxes(f_cpool_to_gresp_storage)       = veg_cf%cpool_to_gresp_storage(p)
    rfluxes(f_cpool_to_ar)                  = ar_p
    if (woody(ivt(p)) == 1._r8) then
      rfluxes(f_cpool_to_livestemc)           = veg_cf%cpool_to_livestemc(p)
      rfluxes(f_cpool_to_livestemc_storage)   = veg_cf%cpool_to_livestemc_storage(p)
      rfluxes(f_cpool_to_deadstemc)           = veg_cf%cpool_to_deadstemc(p)
      rfluxes(f_cpool_to_deadstemc_storage)   = veg_cf%cpool_to_deadstemc_storage(p)
      rfluxes(f_cpool_to_livecrootc)          = veg_cf%cpool_to_livecrootc(p)
      rfluxes(f_cpool_to_livecrootc_storage)  = veg_cf%cpool_to_livecrootc_storage(p)
      rfluxes(f_cpool_to_deadcrootc)          = veg_cf%cpool_to_deadcrootc(p)
      rfluxes(f_cpool_to_deadcrootc_storage)  = veg_cf%cpool_to_deadcrootc_storage(p)
      rfluxes(f_livestemc_to_deadstemc)       = veg_cf%livestemc_to_deadstemc(p)
      rfluxes(f_livestemc_xfer_to_livestemc)  = veg_cf%livestemc_xfer_to_livestemc(p)
      rfluxes(f_livestemc_storage_to_xfer)    = veg_cf%livestemc_storage_to_xfer(p)
      rfluxes(f_deadstemc_xfer_to_deadstemc)  = veg_cf%deadstemc_xfer_to_deadstemc(p)
      rfluxes(f_deadstemc_storage_to_xfer)    = veg_cf%deadstemc_storage_to_xfer(p)
      rfluxes(f_livecrootc_xfer_to_livecrootc)= veg_cf%livecrootc_xfer_to_livecrootc(p)
      rfluxes(f_livecrootc_storage_to_xfer)   = veg_cf%livecrootc_storage_to_xfer(p)
      rfluxes(f_livecrootc_to_deadcrootc)     = veg_cf%livecrootc_to_deadcrootc(p)
      rfluxes(f_deadcrootc_xfer_to_deadcrootc)= veg_cf%deadcrootc_xfer_to_deadcrootc(p)
      rfluxes(f_deadcrootc_storage_to_xfer)   = veg_cf%deadcrootc_storage_to_xfer(p)
      rfluxes(f_gresp_storage_to_xfer)        = veg_cf%gresp_storage_to_xfer(p)
    endif
    if (ivt(p) >= npcropmin) then
      rfluxes(f_cpool_to_livestemc)           = veg_cf%cpool_to_livestemc(p)
      rfluxes(f_cpool_to_livestemc_storage)   = veg_cf%cpool_to_livestemc_storage(p)
      rfluxes(f_cpool_to_grainc)              = veg_cf%cpool_to_grainc(p)
      rfluxes(f_cpool_to_grainc_storage)      = veg_cf%cpool_to_grainc_storage(p)
      rfluxes(f_livestemc_xfer_to_livestemc)  = veg_cf%livestemc_xfer_to_livestemc(p)
      rfluxes(f_livestemc_storage_to_xfer)    = veg_cf%livestemc_storage_to_xfer(p)
      rfluxes(f_livestemc_to_litter)          = veg_cf%livestemc_to_litter(p)
      rfluxes(f_grainc_to_food)               = veg_cf%grainc_to_food(p)
      rfluxes(f_grainc_xfer_to_grainc)        = veg_cf%grainc_xfer_to_grainc(p)
      rfluxes(f_grainc_storage_to_xfer)       = veg_cf%grainc_storage_to_xfer(p)
    endif
    rfluxes(f_leafc_to_litter)              = veg_cf%leafc_to_litter(p)
    rfluxes(f_leafc_xfer_to_leafc)          = veg_cf%leafc_xfer_to_leafc(p)
    rfluxes(f_leafc_storage_to_xfer)        = veg_cf%leafc_storage_to_xfer(p)
    rfluxes(f_frootc_to_litter)             = veg_cf%frootc_to_litter(p)
    rfluxes(f_frootc_xfer_to_frootc)        = veg_cf%frootc_xfer_to_frootc(p)
    rfluxes(f_frootc_storage_to_xfer)       = veg_cf%frootc_storage_to_xfer(p)

    !obtain the limiting factor

    call flux_correction(spm_carbon_d%szrow, spm_carbon_d%szcol, spm_carbon_p, &
       spm_carbon_d, dt, ystates, rfluxes)

    !correct the fluxes

    call fpmax(rfluxes(f_cpool_to_leafc)               , veg_cf%cpool_to_leafc(p))
    call fpmax(rfluxes(f_cpool_to_leafc_storage)       , veg_cf%cpool_to_leafc_storage(p))
    call fpmax(rfluxes(f_cpool_to_frootc)              , veg_cf%cpool_to_frootc(p))
    call fpmax(rfluxes(f_cpool_to_frootc_storage)      , veg_cf%cpool_to_frootc_storage(p))
    call fpmax(rfluxes(f_cpool_to_xsmrpool)            , veg_cf%cpool_to_xsmrpool(p))
    call fpmax(rfluxes(f_cpool_to_gresp_storage)       , veg_cf%cpool_to_gresp_storage(p))

    if (woody(ivt(p)) == 1._r8) then
      call fpmax(rfluxes(f_cpool_to_livestemc)           , veg_cf%cpool_to_livestemc(p))
      call fpmax(rfluxes(f_cpool_to_livestemc_storage)   , veg_cf%cpool_to_livestemc_storage(p))
      call fpmax(rfluxes(f_cpool_to_deadstemc)           , veg_cf%cpool_to_deadstemc(p))
      call fpmax(rfluxes(f_cpool_to_deadstemc_storage)   , veg_cf%cpool_to_deadstemc_storage(p))
      call fpmax(rfluxes(f_cpool_to_livecrootc)          , veg_cf%cpool_to_livecrootc(p))
      call fpmax(rfluxes(f_cpool_to_livecrootc_storage)  , veg_cf%cpool_to_livecrootc_storage(p))
      call fpmax(rfluxes(f_cpool_to_deadcrootc)          , veg_cf%cpool_to_deadcrootc(p))
      call fpmax(rfluxes(f_cpool_to_deadcrootc_storage)  , veg_cf%cpool_to_deadcrootc_storage(p))
      call fpmax(rfluxes(f_livestemc_to_deadstemc)       , veg_cf%livestemc_to_deadstemc(p))
      call fpmax(rfluxes(f_livestemc_xfer_to_livestemc)  , veg_cf%livestemc_xfer_to_livestemc(p))
      call fpmax(rfluxes(f_livestemc_storage_to_xfer)    , veg_cf%livestemc_storage_to_xfer(p))
      call fpmax(rfluxes(f_deadstemc_xfer_to_deadstemc)  , veg_cf%deadstemc_xfer_to_deadstemc(p))
      call fpmax(rfluxes(f_deadstemc_storage_to_xfer)    , veg_cf%deadstemc_storage_to_xfer(p))
      call fpmax(rfluxes(f_livecrootc_xfer_to_livecrootc), veg_cf%livecrootc_xfer_to_livecrootc(p))
      call fpmax(rfluxes(f_livecrootc_storage_to_xfer)   , veg_cf%livecrootc_storage_to_xfer(p))
      call fpmax(rfluxes(f_livecrootc_to_deadcrootc)     , veg_cf%livecrootc_to_deadcrootc(p))
      call fpmax(rfluxes(f_deadcrootc_xfer_to_deadcrootc), veg_cf%deadcrootc_xfer_to_deadcrootc(p))
      call fpmax(rfluxes(f_deadcrootc_storage_to_xfer)   , veg_cf%deadcrootc_storage_to_xfer(p))
      call fpmax(rfluxes(f_gresp_storage_to_xfer)        , veg_cf%gresp_storage_to_xfer(p))
    endif
    if (ivt(p) >= npcropmin) then
      call fpmax(rfluxes(f_cpool_to_livestemc)           , veg_cf%cpool_to_livestemc(p))
      call fpmax(rfluxes(f_cpool_to_livestemc_storage)   , veg_cf%cpool_to_livestemc_storage(p))
      call fpmax(rfluxes(f_cpool_to_grainc)              , veg_cf%cpool_to_grainc(p))
      call fpmax(rfluxes(f_cpool_to_grainc_storage)      , veg_cf%cpool_to_grainc_storage(p))
      call fpmax(rfluxes(f_livestemc_xfer_to_livestemc)  , veg_cf%livestemc_xfer_to_livestemc(p))
      call fpmax(rfluxes(f_livestemc_storage_to_xfer)    , veg_cf%livestemc_storage_to_xfer(p))
      call fpmax(rfluxes(f_livestemc_to_litter)          , veg_cf%livestemc_to_litter(p))
      call fpmax(rfluxes(f_grainc_to_food)               , veg_cf%grainc_to_food(p))
      call fpmax(rfluxes(f_grainc_xfer_to_grainc)        , veg_cf%grainc_xfer_to_grainc(p))
      call fpmax(rfluxes(f_grainc_storage_to_xfer)       , veg_cf%grainc_storage_to_xfer(p))
    endif
    call fpmax(rfluxes(f_leafc_to_litter)              , veg_cf%leafc_to_litter(p))
    call fpmax(rfluxes(f_leafc_xfer_to_leafc)          , veg_cf%leafc_xfer_to_leafc(p))
    call fpmax(rfluxes(f_leafc_storage_to_xfer)        , veg_cf%leafc_storage_to_xfer(p))
    call fpmax(rfluxes(f_frootc_to_litter)             , veg_cf%frootc_to_litter(p))
    call fpmax(rfluxes(f_frootc_xfer_to_frootc)        , veg_cf%frootc_xfer_to_frootc(p))
    call fpmax(rfluxes(f_frootc_storage_to_xfer)       , veg_cf%frootc_storage_to_xfer(p))

    if(rfluxes(f_cpool_to_ar)  <  ar_p)then
      rscal= rfluxes(f_cpool_to_ar)/ar_p
      call ascal(veg_cf%leaf_curmr(p)                , rscal)
      call ascal(veg_cf%froot_curmr(p)               , rscal)
      call ascal(veg_cf%cpool_leaf_gr(p)             , rscal)
      call ascal(veg_cf%cpool_froot_gr(p)            , rscal)
      call ascal(veg_cf%cpool_leaf_storage_gr(p)     , rscal)
      call ascal(veg_cf%cpool_froot_storage_gr(p)    , rscal)
      if (woody(ivt(p)) == 1._r8) then
        call ascal(veg_cf%livestem_curmr(p)            , rscal)
        call ascal(veg_cf%livecroot_curmr(p)           , rscal)
        call ascal(veg_cf%cpool_livestem_gr(p)         , rscal)
        call ascal(veg_cf%cpool_deadstem_gr(p)         , rscal)
        call ascal(veg_cf%cpool_livecroot_gr(p)        , rscal)
        call ascal(veg_cf%cpool_deadcroot_gr(p)        , rscal)
        call ascal(veg_cf%cpool_livestem_storage_gr(p) , rscal)
        call ascal(veg_cf%cpool_deadstem_storage_gr(p) , rscal)
        call ascal(veg_cf%cpool_livecroot_storage_gr(p), rscal)
        call ascal(veg_cf%cpool_deadcroot_storage_gr(p), rscal)
      endif
      if (ivt(p) >= npcropmin) then
        call ascal(veg_cf%livestem_curmr(p)            , rscal)
        call ascal(veg_cf%grain_curmr(p)               , rscal)
        call ascal(veg_cf%cpool_livestem_gr(p)         , rscal)
        call ascal(veg_cf%cpool_grain_gr(p)            , rscal)
        call ascal(veg_cf%cpool_livestem_storage_gr(p) , rscal)
        call ascal(veg_cf%cpool_grain_storage_gr(p)    , rscal)
      endif
    endif
  enddo
  end associate
  end subroutine carbon_flux_limiter
!---------------------------------------------------------------------------
  subroutine nitrogen_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_nf, veg_ns)

    !
    ! DESCRIPTION
    ! the flux limiter for phenology nitrogen fluxes
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
  real(r8) :: ystates(num_nutrient_states)
  real(r8) :: rfluxes(num_nutrient_fluxes)
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
    ystates(s_npool)               = veg_ns%npool(p)
    ystates(s_leafn)               = veg_ns%leafn(p)
    ystates(s_leafn_xfer)          = veg_ns%leafn_xfer(p)
    ystates(s_leafn_storage)       = veg_ns%leafn_storage(p)
    ystates(s_frootn)              = veg_ns%frootn(p)
    ystates(s_frootn_xfer)         = veg_ns%frootn_xfer(p)
    ystates(s_frootn_storage)      = veg_ns%frootn_storage(p)
    if (woody(ivt(p)) == 1.0_r8) then
      ystates(s_livestemn)           = veg_ns%livestemn(p)
      ystates(s_livestemn_xfer)      = veg_ns%livestemn_xfer(p)
      ystates(s_livestemn_storage)   = veg_ns%livestemn_storage(p)
      ystates(s_deadstemn)           = veg_ns%deadstemn(p)
      ystates(s_deadstemn_xfer)      = veg_ns%deadstemn_xfer(p)
      ystates(s_deadstemn_storage)   = veg_ns%deadstemn_storage(p)
      ystates(s_livecrootn)          = veg_ns%livecrootn(p)
      ystates(s_livecrootn_xfer)     = veg_ns%livecrootn_xfer(p)
      ystates(s_livecrootn_storage)  = veg_ns%livecrootn_storage(p)
      ystates(s_deadcrootn)          = veg_ns%deadcrootn(p)
      ystates(s_deadcrootn_xfer)     = veg_ns%deadcrootn_xfer(p)
      ystates(s_deadcrootn_storage)  = veg_ns%deadcrootn_storage(p)
    endif
    if (ivt(p) >= npcropmin) then
      ystates(s_grainn)              = veg_ns%grainn(p)
      ystates(s_grainn_xfer)         = veg_ns%grainn_xfer(p)
      ystates(s_grainn_storage)      = veg_ns%grainn_storage(p)
      ystates(s_livestemn)           = veg_ns%livestemn(p)
      ystates(s_livestemn_xfer)      = veg_ns%livestemn_xfer(p)
      ystates(s_livestemn_storage)   = veg_ns%livestemn_storage(p)
    endif
    ystates(s_retransn)            = veg_ns%retransn(p)

    rfluxes(:) = 0._r8
    !assemble reactive fluxes
    rfluxes(f_npool_to_leafn)                = veg_nf%npool_to_leafn(p)
    rfluxes(f_npool_to_leafn_storage)        = veg_nf%npool_to_leafn_storage(p)
    rfluxes(f_npool_to_frootn)               = veg_nf%npool_to_frootn(p)
    rfluxes(f_npool_to_frootn_storage)       = veg_nf%npool_to_frootn_storage(p)
    if (woody(ivt(p)) == 1.0_r8) then
      rfluxes(f_npool_to_livestemn)            = veg_nf%npool_to_livestemn(p)
      rfluxes(f_npool_to_livestemn_storage)    = veg_nf%npool_to_livestemn_storage(p)
      rfluxes(f_npool_to_livecrootn)           = veg_nf%npool_to_livecrootn(p)
      rfluxes(f_npool_to_livecrootn_storage)   = veg_nf%npool_to_livecrootn_storage(p)
      rfluxes(f_npool_to_deadstemn)            = veg_nf%npool_to_deadstemn(p)
      rfluxes(f_npool_to_deadcrootn)           = veg_nf%npool_to_deadcrootn(p)
      rfluxes(f_npool_to_deadstemn_storage)    = veg_nf%npool_to_deadstemn_storage(p)
      rfluxes(f_npool_to_deadcrootn_storage)   = veg_nf%npool_to_deadcrootn_storage(p)
      rfluxes(f_livestemn_storage_to_xfer)     = veg_nf%livestemn_storage_to_xfer(p)
      rfluxes(f_livestemn_xfer_to_livestemn)   = veg_nf%livestemn_xfer_to_livestemn(p)
      rfluxes(f_livestemn_to_retransn)         = veg_nf%livestemn_to_retransn(p)
      rfluxes(f_livestemn_to_deadstemn)        = veg_nf%livestemn_to_deadstemn(p)
      rfluxes(f_livecrootn_to_deadcrootn)      = veg_nf%livecrootn_to_deadcrootn(p)
      rfluxes(f_livecrootn_to_retransn)        = veg_nf%livecrootn_to_retransn(p)
      rfluxes(f_livecrootn_storage_to_xfer)    = veg_nf%livecrootn_storage_to_xfer(p)
      rfluxes(f_livecrootn_xfer_to_livecrootn) = veg_nf%livecrootn_xfer_to_livecrootn(p)
      rfluxes(f_deadstemn_storage_to_xfer)     = veg_nf%deadstemn_storage_to_xfer(p)
      rfluxes(f_deadstemn_xfer_to_deadstem)    = veg_nf%deadstemn_xfer_to_deadstemn(p)
      rfluxes(f_deadcrootn_storage_to_xfer)    = veg_nf%deadcrootn_storage_to_xfer(p)
      rfluxes(f_deadcrootn_xfer_to_deadcrootn) = veg_nf%deadcrootn_xfer_to_deadcrootn(p)
    endif

    if (ivt(p) >= npcropmin) then
      rfluxes(f_npool_to_livestemn)            = veg_nf%npool_to_livestemn(p)
      rfluxes(f_npool_to_livestemn_storage)    = veg_nf%npool_to_livestemn_storage(p)
      rfluxes(f_npool_to_grainn)               = veg_nf%npool_to_grainn(p)
      rfluxes(f_npool_to_grainn_storage)       = veg_nf%npool_to_grainn_storage(p)
      rfluxes(f_frootn_to_retransn)            = veg_nf%frootn_to_retransn(p)
      rfluxes(f_livestemn_to_litter)           = veg_nf%livestemn_to_litter(p)
      rfluxes(f_livestemn_storage_to_xfer)     = veg_nf%livestemn_storage_to_xfer(p)
      rfluxes(f_livestemn_xfer_to_livestemn)   =  veg_nf%livestemn_xfer_to_livestemn(p)
      rfluxes(f_livestemn_to_retransn)         = veg_nf%livestemn_to_retransn(p)
      rfluxes(f_grainn_xfer_to_grainn)         = veg_nf%grainn_xfer_to_grainn(p)
      rfluxes(f_grainn_to_food)                = veg_nf%grainn_to_food(p)
    endif

    rfluxes(f_leafn_to_litter)               = veg_nf%leafn_to_litter(p)
    rfluxes(f_leafn_xfer_to_leafn)           = veg_nf%leafn_xfer_to_leafn(p)
    rfluxes(f_leafn_storage_to_xfer)         = veg_nf%leafn_storage_to_xfer(p)
    rfluxes(f_frootn_xfer_to_frootn)         = veg_nf%frootn_xfer_to_frootn(p)
    rfluxes(f_frootn_storage_to_xfer)        = veg_nf%frootn_storage_to_xfer(p)
    rfluxes(f_leafn_to_retransn)             = veg_nf%leafn_to_retransn(p)
    rfluxes(f_frootn_to_litter)              = veg_nf%frootn_to_litter(p)

    rfluxes(f_retransn_to_npool)             = veg_nf%retransn_to_npool(p)
    rfluxes(f_supplement_to_plantn)          = veg_nf%supplement_to_plantn(p)

    call flux_correction(spm_nutrient_d%szrow, spm_nutrient_d%szcol, spm_nutrient_p, &
      spm_nutrient_d, dt, ystates, rfluxes)

    !correct the fluxes
    call fpmax(rfluxes(f_npool_to_leafn)                , veg_nf%npool_to_leafn(p))
    call fpmax(rfluxes(f_npool_to_leafn_storage)        , veg_nf%npool_to_leafn_storage(p))
    call fpmax(rfluxes(f_npool_to_frootn)               , veg_nf%npool_to_frootn(p))
    call fpmax(rfluxes(f_npool_to_frootn_storage)       , veg_nf%npool_to_frootn_storage(p))
    if (woody(ivt(p)) == 1.0_r8) then
      call fpmax(rfluxes(f_npool_to_livestemn)            , veg_nf%npool_to_livestemn(p))
      call fpmax(rfluxes(f_npool_to_livestemn_storage)    , veg_nf%npool_to_livestemn_storage(p))
      call fpmax(rfluxes(f_npool_to_livecrootn)           , veg_nf%npool_to_livecrootn(p))
      call fpmax(rfluxes(f_npool_to_livecrootn_storage)   , veg_nf%npool_to_livecrootn_storage(p))
      call fpmax(rfluxes(f_npool_to_deadstemn)            , veg_nf%npool_to_deadstemn(p))
      call fpmax(rfluxes(f_npool_to_deadcrootn)           , veg_nf%npool_to_deadcrootn(p))
      call fpmax(rfluxes(f_npool_to_deadstemn_storage)    , veg_nf%npool_to_deadstemn_storage(p))
      call fpmax(rfluxes(f_npool_to_deadcrootn_storage)   , veg_nf%npool_to_deadcrootn_storage(p))
      call fpmax(rfluxes(f_livestemn_storage_to_xfer)     , veg_nf%livestemn_storage_to_xfer(p))
      call fpmax(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_nf%livestemn_xfer_to_livestemn(p))
      call fpmax(rfluxes(f_livestemn_to_retransn)         , veg_nf%livestemn_to_retransn(p))
      call fpmax(rfluxes(f_livestemn_to_deadstemn)        , veg_nf%livestemn_to_deadstemn(p))
      call fpmax(rfluxes(f_livecrootn_to_deadcrootn)      , veg_nf%livecrootn_to_deadcrootn(p))
      call fpmax(rfluxes(f_livecrootn_to_retransn)        , veg_nf%livecrootn_to_retransn(p))
      call fpmax(rfluxes(f_livecrootn_storage_to_xfer)    , veg_nf%livecrootn_storage_to_xfer(p))
      call fpmax(rfluxes(f_livecrootn_xfer_to_livecrootn) , veg_nf%livecrootn_xfer_to_livecrootn(p))
      call fpmax(rfluxes(f_deadstemn_storage_to_xfer)     , veg_nf%deadstemn_storage_to_xfer(p))
      call fpmax(rfluxes(f_deadstemn_xfer_to_deadstem)    , veg_nf%deadstemn_xfer_to_deadstemn(p))
      call fpmax(rfluxes(f_deadcrootn_storage_to_xfer)    , veg_nf%deadcrootn_storage_to_xfer(p))
      call fpmax(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , veg_nf%deadcrootn_xfer_to_deadcrootn(p))
    endif

    if (ivt(p) >= npcropmin) then
      call fpmax(rfluxes(f_npool_to_livestemn)            , veg_nf%npool_to_livestemn(p))
      call fpmax(rfluxes(f_npool_to_livestemn_storage)    , veg_nf%npool_to_livestemn_storage(p))
      call fpmax(rfluxes(f_npool_to_grainn)               , veg_nf%npool_to_grainn(p))
      call fpmax(rfluxes(f_npool_to_grainn_storage)       , veg_nf%npool_to_grainn_storage(p))
      call fpmax(rfluxes(f_frootn_to_retransn)            , veg_nf%frootn_to_retransn(p))
      call fpmax(rfluxes(f_livestemn_to_litter)           , veg_nf%livestemn_to_litter(p))
      call fpmax(rfluxes(f_livestemn_storage_to_xfer)     , veg_nf%livestemn_storage_to_xfer(p))
      call fpmax(rfluxes(f_livestemn_xfer_to_livestemn)   ,  veg_nf%livestemn_xfer_to_livestemn(p))
      call fpmax(rfluxes(f_livestemn_to_retransn)         , veg_nf%livestemn_to_retransn(p))
      call fpmax(rfluxes(f_grainn_xfer_to_grainn)         , veg_nf%grainn_xfer_to_grainn(p))
      call fpmax(rfluxes(f_grainn_to_food)                , veg_nf%grainn_to_food(p))
    endif

    rfluxes(f_leafn_to_litter)               = veg_nf%leafn_to_litter(p)
    rfluxes(f_leafn_xfer_to_leafn)           = veg_nf%leafn_xfer_to_leafn(p)
    rfluxes(f_leafn_storage_to_xfer)         = veg_nf%leafn_storage_to_xfer(p)
    rfluxes(f_frootn_xfer_to_frootn)         = veg_nf%frootn_xfer_to_frootn(p)
    rfluxes(f_frootn_storage_to_xfer)        = veg_nf%frootn_storage_to_xfer(p)
    rfluxes(f_leafn_to_retransn)             = veg_nf%leafn_to_retransn(p)
    rfluxes(f_frootn_to_litter)              = veg_nf%frootn_to_litter(p)

    rfluxes(f_retransn_to_npool)             = veg_nf%retransn_to_npool(p)
    rfluxes(f_supplement_to_plantn)          = veg_nf%supplement_to_plantn(p)

  enddo
  end associate
  end subroutine nitrogen_flux_limiter

  !-------------------------------------------------------------------------------
  subroutine phosphorus_flux_limiter(bounds, num_soilc, filter_soilc,&
      num_soilp, filter_soilp, cnstate_vars,  &
      veg_pf, veg_ps)
    ! DESCRIPTION
    ! the flux limiter for phenology phosphorus  fluxes
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
  real(r8) :: ystates(num_nutrient_states)
  real(r8) :: rfluxes(num_nutrient_fluxes)

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
    ystates(s_npool)               = veg_ps%ppool(p)
    ystates(s_leafn)               = veg_ps%leafp(p)
    ystates(s_leafn_xfer)          = veg_ps%leafp_xfer(p)
    ystates(s_leafn_storage)       = veg_ps%leafp_storage(p)
    ystates(s_frootn)              = veg_ps%frootp(p)
    ystates(s_frootn_xfer)         = veg_ps%frootp_xfer(p)
    ystates(s_frootn_storage)      = veg_ps%frootp_storage(p)
    if (woody(ivt(p)) == 1.0_r8) then
      ystates(s_livestemn)           = veg_ps%livestemp(p)
      ystates(s_livestemn_xfer)      = veg_ps%livestemp_xfer(p)
      ystates(s_livestemn_storage)   = veg_ps%livestemp_storage(p)
      ystates(s_deadstemn)           = veg_ps%deadstemp(p)
      ystates(s_deadstemn_xfer)      = veg_ps%deadstemp_xfer(p)
      ystates(s_deadstemn_storage)   = veg_ps%deadstemp_storage(p)
      ystates(s_livecrootn)          = veg_ps%livecrootp(p)
      ystates(s_livecrootn_xfer)     = veg_ps%livecrootp_xfer(p)
      ystates(s_livecrootn_storage)  = veg_ps%livecrootp_storage(p)
      ystates(s_deadcrootn)          = veg_ps%deadcrootp(p)
      ystates(s_deadcrootn_xfer)     = veg_ps%deadcrootp_xfer(p)
      ystates(s_deadcrootn_storage)  = veg_ps%deadcrootp_storage(p)
    endif
    if (ivt(p) >= npcropmin) then
      ystates(s_grainn)              = veg_ps%grainp(p)
      ystates(s_grainn_xfer)         = veg_ps%grainp_xfer(p)
      ystates(s_grainn_storage)      = veg_ps%grainp_storage(p)
      ystates(s_livestemn)           = veg_ps%livestemp(p)
      ystates(s_livestemn_xfer)      = veg_ps%livestemp_xfer(p)
      ystates(s_livestemn_storage)   = veg_ps%livestemp_storage(p)
    endif
    ystates(s_retransn)              = veg_ps%retransp(p)

    rfluxes(:) = 0._r8
    !assemble reactive fluxes
    rfluxes(f_npool_to_leafn)                = veg_pf%ppool_to_leafp(p)
    rfluxes(f_npool_to_leafn_storage)        = veg_pf%ppool_to_leafp_storage(p)
    rfluxes(f_npool_to_frootn)               = veg_pf%ppool_to_frootp(p)
    rfluxes(f_npool_to_frootn_storage)       = veg_pf%ppool_to_frootp_storage(p)
    if (woody(ivt(p)) == 1._r8) then
      rfluxes(f_npool_to_livestemn)            = veg_pf%ppool_to_livestemp(p)
      rfluxes(f_npool_to_livestemn_storage)    = veg_pf%ppool_to_livestemp_storage(p)
      rfluxes(f_npool_to_livecrootn)           = veg_pf%ppool_to_livecrootp(p)
      rfluxes(f_npool_to_livecrootn_storage)   = veg_pf%ppool_to_livecrootp_storage(p)
      rfluxes(f_npool_to_deadstemn)            = veg_pf%ppool_to_deadstemp(p)
      rfluxes(f_npool_to_deadcrootn)           = veg_pf%ppool_to_deadcrootp(p)
      rfluxes(f_npool_to_deadstemn_storage)    = veg_pf%ppool_to_deadstemp_storage(p)
      rfluxes(f_npool_to_deadcrootn_storage)   = veg_pf%ppool_to_deadcrootp_storage(p)
      rfluxes(f_livestemn_storage_to_xfer)     = veg_pf%livestemp_storage_to_xfer(p)
      rfluxes(f_livestemn_xfer_to_livestemn)   = veg_pf%livestemp_xfer_to_livestemp(p)
      rfluxes(f_livestemn_to_retransn)         = veg_pf%livestemp_to_retransp(p)
      rfluxes(f_livestemn_to_deadstemn)        = veg_pf%livestemp_to_deadstemp(p)
      rfluxes(f_livecrootn_to_deadcrootn)      = veg_pf%livecrootp_to_deadcrootp(p)
      rfluxes(f_livecrootn_to_retransn)        = veg_pf%livecrootp_to_retransp(p)
      rfluxes(f_livecrootn_storage_to_xfer)    = veg_pf%livecrootp_storage_to_xfer(p)
      rfluxes(f_livecrootn_xfer_to_livecrootn) = veg_pf%livecrootp_xfer_to_livecrootp(p)
      rfluxes(f_deadstemn_storage_to_xfer)     = veg_pf%deadstemp_storage_to_xfer(p)
      rfluxes(f_deadstemn_xfer_to_deadstem)    = veg_pf%deadstemp_xfer_to_deadstemp(p)
      rfluxes(f_deadcrootn_storage_to_xfer)    = veg_pf%deadcrootp_storage_to_xfer(p)
      rfluxes(f_deadcrootn_xfer_to_deadcrootn) = veg_pf%deadcrootp_xfer_to_deadcrootp(p)
    endif

    if (ivt(p) >= npcropmin) then
      rfluxes(f_npool_to_livestemn)            = veg_pf%ppool_to_livestemp(p)
      rfluxes(f_npool_to_livestemn_storage)    = veg_pf%ppool_to_livestemp_storage(p)
      rfluxes(f_npool_to_grainn)               = veg_pf%ppool_to_grainp(p)
      rfluxes(f_npool_to_grainn_storage)       = veg_pf%ppool_to_grainp_storage(p)
      rfluxes(f_frootn_to_retransn)            = veg_pf%frootp_to_retransp(p)
      rfluxes(f_livestemn_to_litter)           = veg_pf%livestemp_to_litter(p)
      rfluxes(f_livestemn_storage_to_xfer)     = veg_pf%livestemp_storage_to_xfer(p)
      rfluxes(f_livestemn_xfer_to_livestemn)   = veg_pf%livestemp_xfer_to_livestemp(p)
      rfluxes(f_livestemn_to_retransn)         = veg_pf%livestemp_to_retransp(p)
      rfluxes(f_grainn_xfer_to_grainn)         = veg_pf%grainp_xfer_to_grainp(p)
      rfluxes(f_grainn_to_food)                = veg_pf%grainp_to_food(p)
    endif

    rfluxes(f_leafn_to_litter)               = veg_pf%leafp_to_litter(p)
    rfluxes(f_leafn_xfer_to_leafn)           = veg_pf%leafp_xfer_to_leafp(p)
    rfluxes(f_leafn_storage_to_xfer)         = veg_pf%ppool_to_leafp_storage(p)
    rfluxes(f_frootn_xfer_to_frootn)         = veg_pf%frootp_xfer_to_frootp(p)
    rfluxes(f_frootn_storage_to_xfer)        = veg_pf%frootp_storage_to_xfer(p)
    rfluxes(f_leafn_to_retransn)             = veg_pf%leafp_to_retransp(p)
    rfluxes(f_frootn_to_litter)              = veg_pf%frootp_to_litter(p)
    rfluxes(f_retransn_to_npool)             = veg_pf%retransp_to_ppool(p)
    rfluxes(f_supplement_to_plantn)          = veg_pf%supplement_to_plantp(p)
    !assemble stoichiometry matrix

    !obtain the limiting factor
    call flux_correction(spm_nutrient_d%szrow, spm_nutrient_d%szcol, spm_nutrient_p, &
      spm_nutrient_d, dt, ystates, rfluxes)
    !correct the fluxes
    call fpmax(rfluxes(f_npool_to_leafn)                , veg_pf%ppool_to_leafp(p))
    call fpmax(rfluxes(f_npool_to_leafn_storage)        , veg_pf%ppool_to_leafp_storage(p))
    call fpmax(rfluxes(f_npool_to_frootn)               , veg_pf%ppool_to_frootp(p))
    call fpmax(rfluxes(f_npool_to_frootn_storage)       , veg_pf%ppool_to_frootp_storage(p))
    if (woody(ivt(p)) == 1._r8) then
      call fpmax(rfluxes(f_npool_to_livestemn)            , veg_pf%ppool_to_livestemp(p))
      call fpmax(rfluxes(f_npool_to_livestemn_storage)    , veg_pf%ppool_to_livestemp_storage(p))
      call fpmax(rfluxes(f_npool_to_livecrootn)           , veg_pf%ppool_to_livecrootp(p))
      call fpmax(rfluxes(f_npool_to_livecrootn_storage)   , veg_pf%ppool_to_livecrootp_storage(p))
      call fpmax(rfluxes(f_npool_to_deadstemn)            , veg_pf%ppool_to_deadstemp(p))
      call fpmax(rfluxes(f_npool_to_deadcrootn)           , veg_pf%ppool_to_deadcrootp(p))
      call fpmax(rfluxes(f_npool_to_deadstemn_storage)    , veg_pf%ppool_to_deadstemp_storage(p))
      call fpmax(rfluxes(f_npool_to_deadcrootn_storage)   , veg_pf%ppool_to_deadcrootp_storage(p))
      call fpmax(rfluxes(f_livestemn_storage_to_xfer)     , veg_pf%livestemp_storage_to_xfer(p))
      call fpmax(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_pf%livestemp_xfer_to_livestemp(p))
      call fpmax(rfluxes(f_livestemn_to_retransn)         , veg_pf%livestemp_to_retransp(p))
      call fpmax(rfluxes(f_livestemn_to_deadstemn)        , veg_pf%livestemp_to_deadstemp(p))
      call fpmax(rfluxes(f_livecrootn_to_deadcrootn)      , veg_pf%livecrootp_to_deadcrootp(p))
      call fpmax(rfluxes(f_livecrootn_to_retransn)        , veg_pf%livecrootp_to_retransp(p))
      call fpmax(rfluxes(f_livecrootn_storage_to_xfer)    , veg_pf%livecrootp_storage_to_xfer(p))
      call fpmax(rfluxes(f_livecrootn_xfer_to_livecrootn) , veg_pf%livecrootp_xfer_to_livecrootp(p))
      call fpmax(rfluxes(f_deadstemn_storage_to_xfer)     , veg_pf%deadstemp_storage_to_xfer(p))
      call fpmax(rfluxes(f_deadstemn_xfer_to_deadstem)    , veg_pf%deadstemp_xfer_to_deadstemp(p))
      call fpmax(rfluxes(f_deadcrootn_storage_to_xfer)    , veg_pf%deadcrootp_storage_to_xfer(p))
      call fpmax(rfluxes(f_deadcrootn_xfer_to_deadcrootn) , veg_pf%deadcrootp_xfer_to_deadcrootp(p))
    endif

    if (ivt(p) >= npcropmin) then
      call fpmax(rfluxes(f_npool_to_livestemn)            , veg_pf%ppool_to_livestemp(p))
      call fpmax(rfluxes(f_npool_to_livestemn_storage)    , veg_pf%ppool_to_livestemp_storage(p))
      call fpmax(rfluxes(f_npool_to_grainn)               , veg_pf%ppool_to_grainp(p))
      call fpmax(rfluxes(f_npool_to_grainn_storage)       , veg_pf%ppool_to_grainp_storage(p))
      call fpmax(rfluxes(f_frootn_to_retransn)            , veg_pf%frootp_to_retransp(p))
      call fpmax(rfluxes(f_livestemn_to_litter)           , veg_pf%livestemp_to_litter(p))
      call fpmax(rfluxes(f_livestemn_storage_to_xfer)     , veg_pf%livestemp_storage_to_xfer(p))
      call fpmax(rfluxes(f_livestemn_xfer_to_livestemn)   , veg_pf%livestemp_xfer_to_livestemp(p))
      call fpmax(rfluxes(f_livestemn_to_retransn)         , veg_pf%livestemp_to_retransp(p))
      call fpmax(rfluxes(f_grainn_xfer_to_grainn)         , veg_pf%grainp_xfer_to_grainp(p))
      call fpmax(rfluxes(f_grainn_to_food)                , veg_pf%grainp_to_food(p))
    endif

    call fpmax(rfluxes(f_leafn_to_litter)               , veg_pf%leafp_to_litter(p))
    call fpmax(rfluxes(f_leafn_xfer_to_leafn)           , veg_pf%leafp_xfer_to_leafp(p))
    call fpmax(rfluxes(f_leafn_storage_to_xfer)         , veg_pf%ppool_to_leafp_storage(p))
    call fpmax(rfluxes(f_frootn_xfer_to_frootn)         , veg_pf%frootp_xfer_to_frootp(p))
    call fpmax(rfluxes(f_frootn_storage_to_xfer)        , veg_pf%frootp_storage_to_xfer(p))
    call fpmax(rfluxes(f_leafn_to_retransn)             , veg_pf%leafp_to_retransp(p))
    call fpmax(rfluxes(f_frootn_to_litter)              , veg_pf%frootp_to_litter(p))
    call fpmax(rfluxes(f_retransn_to_npool)             , veg_pf%retransp_to_ppool(p))
    call fpmax(rfluxes(f_supplement_to_plantn)          , veg_pf%supplement_to_plantp(p))

  enddo
  end associate
  end subroutine phosphorus_flux_limiter

end module PhenologyFLuxLimitMod
