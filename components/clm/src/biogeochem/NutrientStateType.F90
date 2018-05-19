module NutrientStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use decompMod              , only : bounds_type
  use clm_varpar             , only : nlevdecomp_full, ndecomp_pools, crop_prog
  use clm_varctl             , only : use_fates, use_crop, use_vertsoilc
  use clm_varcon             , only : spval, ispval
  use CNDecompCascadeConType , only : decomp_cascade_con
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use abortutils             , only : endrun
  use VegetationType         , only : veg_pp
  use pftvarcon              , only : npcropmin
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  type, public :: nutrientstate_type

     integer           :: species                         ! C12, C13, C14, N, P

     character(len=3)  :: name
     character(len=4)  :: history_name_prefix
     character(len=3)  :: restart_name

     ! Vegetation pools
     real(r8), pointer :: deadcroot_patch         (:)     ! (mass/m2) dead coarse root
     real(r8), pointer :: deadcroot_storage_patch (:)     ! (mass/m2) dead coarse root storage
     real(r8), pointer :: deadcroot_xfer_patch    (:)     ! (mass/m2) dead coarse root transfer
     real(r8), pointer :: deadstem_patch          (:)     ! (mass/m2) dead stem
     real(r8), pointer :: deadstem_storage_patch  (:)     ! (mass/m2) dead stem storage
     real(r8), pointer :: deadstem_xfer_patch     (:)     ! (mass/m2) dead stem transfer
     real(r8), pointer :: froot_patch             (:)     ! (mass/m2) fine root
     real(r8), pointer :: froot_storage_patch     (:)     ! (mass/m2) fine root storage
     real(r8), pointer :: froot_xfer_patch        (:)     ! (mass/m2) fine root transfer
     real(r8), pointer :: grain_patch             (:)     ! (mass/m2) grain
     real(r8), pointer :: grain_storage_patch     (:)     ! (mass/m2) grain storage
     real(r8), pointer :: grain_xfer_patch        (:)     ! (mass/m2) grain transfer
     real(r8), pointer :: leaf_patch              (:)     ! (mass/m2) leaf
     real(r8), pointer :: leaf_storage_patch      (:)     ! (mass/m2) leaf storage
     real(r8), pointer :: leaf_xfer_patch         (:)     ! (mass/m2) leaf transfer
     real(r8), pointer :: livecroot_patch         (:)     ! (mass/m2) live coarse root
     real(r8), pointer :: livecroot_storage_patch (:)     ! (mass/m2) live coarse root storage
     real(r8), pointer :: livecroot_xfer_patch    (:)     ! (mass/m2) live coarse root transfer
     real(r8), pointer :: livestem_patch          (:)     ! (mass/m2) live stem
     real(r8), pointer :: livestem_storage_patch  (:)     ! (mass/m2) live stem storage
     real(r8), pointer :: livestem_xfer_patch     (:)     ! (mass/m2) live stem transfer

     ! Balance check
     real(r8), pointer :: beg_bal_patch           (:)     ! (mass/m2) beginning mass balance at patch
     real(r8), pointer :: beg_bal_col             (:)     ! (mass/m2) beginning mass balance at column
     real(r8), pointer :: beg_bal_grc             (:)     ! (mass/m2) beginning mass balance at grid
     real(r8), pointer :: end_bal_patch           (:)     ! (mass/m2) ending mass balance at patch
     real(r8), pointer :: end_bal_col             (:)     ! (mass/m2) ending mass balance at column
     real(r8), pointer :: end_bal_grc             (:)     ! (mass/m2) ending mass balance at grid
     real(r8), pointer :: err_bal_patch           (:)     ! (mass/m2) mass balance error at patch
     real(r8), pointer :: err_bal_col             (:)     ! (mass/m2) mass balance error at column
     real(r8), pointer :: err_bal_grc             (:)     ! (mass/m2) mass balance error at grid

     ! Soil pools
     real(r8), pointer :: decomp_pools_vr_col     (:,:,:) ! (mass/m2) vertically-resolved soil decomposing pools

     ! (mass/m2) Truncation
     real(r8), pointer :: veg_trunc_col           (:)     ! (mass/m2) patch-level vegetation nutrient sink for truncation
     real(r8), pointer :: veg_trunc_patch         (:)     ! (mass/m2) column-level vegetation nutrient sink for truncation
     real(r8), pointer :: soil_trunc_vr_col       (:,:)   ! (mass/m2) column-level soil nutrient sink for truncation

     ! Pools for dynamic landcover
     real(r8), pointer :: cropseed_deficit_patch  (:)     ! (mass/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
     real(r8), pointer :: seed_col                (:)     ! (mass/m2) column-level pool for seeding new patches
     real(r8), pointer :: seed_grc                (:)     ! (mass/m2) gridcell-level pool for seeding new PFTs via dynamic landcover
     real(r8), pointer :: prod1_col               (:)     ! (mass/m2) wood product pool, 1-year lifespan
     real(r8), pointer :: prod10_col              (:)     ! (mass/m2) wood product pool, 10-year lifespan
     real(r8), pointer :: prod100_col             (:)     ! (mass/m2) wood product pool, 100-year lifespan

     ! Summary (diagnostic) state variables
     real(r8), pointer :: cwd_col                 (:)     ! (mass/m2) coarse woody debris
     real(r8), pointer :: decomp_pools_1m_col     (:,:)   ! (mass/m2) decomposing (litter, cwd, soil) pools to 1 meter
     real(r8), pointer :: decomp_pools_col        (:,:)   ! (mass/m2) decomposing (litter, cwd, soil) pools
     real(r8), pointer :: dispveg_patch           (:)     ! (mass/m2) displayed veg carbon
     real(r8), pointer :: dyn_bal_adjustments_col (:)     ! (mass/m2) adjustments to each column made in this timestep via dynamic column area adjustments
     real(r8), pointer :: pool_patch              (:)     ! (mass/m2) temporary photosynthate pool
     real(r8), pointer :: storveg_patch           (:)     ! (mass/m2) stored vegetation nutrient
     real(r8), pointer :: totabg_col              (:)     ! (mass/m2) total column above ground nutrient
     real(r8), pointer :: totblg_col              (:)     ! (mass/m2) total column below ground nutrient
     real(r8), pointer :: totcol_col              (:)     ! (mass/m2) total column nutrient
     real(r8), pointer :: totecosys_col           (:)     ! (mass/m2) total ecosystem nutrient
     real(r8), pointer :: totlit_1m_col           (:)     ! (mass/m2) total litter upto 1 meter
     real(r8), pointer :: totlit_col              (:)     ! (mass/m2) total litter
     real(r8), pointer :: totpft_col              (:)     ! (mass/m2) total PFT nutrient at column-level
     real(r8), pointer :: totpft_patch            (:)     ! (mass/m2) total PFT nutrient at patch-level
     real(r8), pointer :: totprod_col             (:)     ! (mass/m2) total product nutrient at column-level
     real(r8), pointer :: totsom_1m_col           (:)     ! (mass/m2) total soil organic matter upto 1 meter
     real(r8), pointer :: totsom_col              (:)     ! (mass/m2) total soil organic matter
     real(r8), pointer :: totveg_patch            (:)     ! (mass/m2) total vegetation nutrient at patch-level
     real(r8), pointer :: totveg_col              (:)     ! (mass/m2) total vegetation nutrient at column-level
     real(r8), pointer :: totveg_abg_patch        (:)     ! (mass/m2) total above ground vegetation nutrient at patch-level
     real(r8), pointer :: totveg_abg_col          (:)     ! (mass/m2) total above ground vegetation nutrient at column-level

  end type nutrientstate_type

  public :: NutrientStateInitAllocate, &
            NutrientStateInitHistory, &
            NutrientStateDynamicPatchAdjustments, &
            NutrientStateRestart, &
            NutrientStatePatchSummary, &
            NutrientStateColumnSummary

contains

  !------------------------------------------------------------------------
  subroutine NutrientStateInitAllocate(nutrient_state, bounds)
    !
    implicit none
    !
    class (nutrientstate_type)    :: nutrient_state
    type(bounds_type), intent(in) :: bounds
    !
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    if (.not. use_fates) then
       allocate(nutrient_state%deadcroot_patch         (begp :endp)) ;  nutrient_state%deadcroot_patch         (:)   = nan
       allocate(nutrient_state%deadcroot_storage_patch (begp :endp)) ;  nutrient_state%deadcroot_storage_patch (:)   = nan
       allocate(nutrient_state%deadcroot_xfer_patch    (begp :endp)) ;  nutrient_state%deadcroot_xfer_patch    (:)   = nan
       allocate(nutrient_state%deadstem_patch          (begp :endp)) ;  nutrient_state%deadstem_patch          (:)   = nan
       allocate(nutrient_state%deadstem_storage_patch  (begp :endp)) ;  nutrient_state%deadstem_storage_patch  (:)   = nan
       allocate(nutrient_state%deadstem_xfer_patch     (begp :endp)) ;  nutrient_state%deadstem_xfer_patch     (:)   = nan
       allocate(nutrient_state%froot_patch             (begp :endp)) ;  nutrient_state%froot_patch             (:)   = nan
       allocate(nutrient_state%froot_storage_patch     (begp :endp)) ;  nutrient_state%froot_storage_patch     (:)   = nan

       allocate(nutrient_state%froot_xfer_patch        (begp :endp)) ;  nutrient_state%froot_xfer_patch        (:)   = nan
       allocate(nutrient_state%grain_patch             (begp :endp)) ;  nutrient_state%grain_patch             (:)   = nan
       allocate(nutrient_state%grain_storage_patch     (begp :endp)) ;  nutrient_state%grain_storage_patch     (:)   = nan
       allocate(nutrient_state%grain_xfer_patch        (begp :endp)) ;  nutrient_state%grain_xfer_patch        (:)   = nan
       allocate(nutrient_state%leaf_patch              (begp :endp)) ;  nutrient_state%leaf_patch              (:)   = nan
       allocate(nutrient_state%leaf_storage_patch      (begp :endp)) ;  nutrient_state%leaf_storage_patch      (:)   = nan
       allocate(nutrient_state%leaf_xfer_patch         (begp :endp)) ;  nutrient_state%leaf_xfer_patch         (:)   = nan
       allocate(nutrient_state%livecroot_patch         (begp :endp)) ;  nutrient_state%livecroot_patch         (:)   = nan
       allocate(nutrient_state%livecroot_storage_patch (begp :endp)) ;  nutrient_state%livecroot_storage_patch (:)   = nan
       allocate(nutrient_state%livecroot_xfer_patch    (begp :endp)) ;  nutrient_state%livecroot_xfer_patch    (:)   = nan
       allocate(nutrient_state%livestem_patch          (begp :endp)) ;  nutrient_state%livestem_patch          (:)   = nan
       allocate(nutrient_state%livestem_storage_patch  (begp :endp)) ;  nutrient_state%livestem_storage_patch  (:)   = nan
       allocate(nutrient_state%livestem_xfer_patch     (begp :endp)) ;  nutrient_state%livestem_xfer_patch     (:)   = nan

       allocate(nutrient_state%veg_trunc_patch         (begp :endp)) ;  nutrient_state%veg_trunc_patch         (:)   = nan

       allocate(nutrient_state%dispveg_patch           (begp :endp)) ;  nutrient_state%dispveg_patch           (:)   = nan
       allocate(nutrient_state%pool_patch              (begp :endp)) ;  nutrient_state%pool_patch              (:)   = nan
       allocate(nutrient_state%storveg_patch           (begp :endp)) ;  nutrient_state%storveg_patch           (:)   = nan
       allocate(nutrient_state%totveg_patch            (begp :endp)) ;  nutrient_state%totveg_patch            (:)   = nan
       allocate(nutrient_state%totpft_patch            (begp :endp)) ;  nutrient_state%totpft_patch            (:)   = nan
       allocate(nutrient_state%totveg_abg_patch        (begp :endp)) ;  nutrient_state%totveg_abg_patch        (:)   = nan
    end if

    allocate(nutrient_state%beg_bal_col             (begc:endc))                   ; nutrient_state%beg_bal_col             (:)     = nan
    allocate(nutrient_state%beg_bal_grc             (begg:endg))                   ; nutrient_state%beg_bal_grc             (:)     = nan
    allocate(nutrient_state%end_bal_patch           (begp:endp))                   ; nutrient_state%end_bal_patch           (:)     = nan
    allocate(nutrient_state%end_bal_col             (begc:endc))                   ; nutrient_state%end_bal_col             (:)     = nan
    allocate(nutrient_state%end_bal_grc             (begg:endg))                   ; nutrient_state%end_bal_grc             (:)     = nan
    allocate(nutrient_state%err_bal_patch           (begp:endp))                   ; nutrient_state%err_bal_patch           (:)     = nan
    allocate(nutrient_state%err_bal_col             (begc:endc))                   ; nutrient_state%err_bal_col             (:)     = nan
    allocate(nutrient_state%err_bal_grc             (begg:endg))                   ; nutrient_state%err_bal_grc             (:)     = nan

    allocate(nutrient_state%cropseed_deficit_patch  (begp:endp))                   ; nutrient_state%cropseed_deficit_patch  (:)     = nan
    allocate(nutrient_state%seed_grc                (begg:endg))                   ; nutrient_state%seed_grc                (:)     = nan
    allocate(nutrient_state%seed_col                (begc:endc))                   ; nutrient_state%seed_col                (:)     = nan
    allocate(nutrient_state%prod1_col               (begc:endc))                   ; nutrient_state%prod1_col               (:)     = nan
    allocate(nutrient_state%prod10_col              (begc:endc))                   ; nutrient_state%prod10_col              (:)     = nan
    allocate(nutrient_state%prod100_col             (begc:endc))                   ; nutrient_state%prod100_col             (:)     = nan
    
    allocate(nutrient_state%decomp_pools_vr_col     (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    nutrient_state%decomp_pools_vr_col(:,:,:) = nan

    allocate(nutrient_state%veg_trunc_col           (begc:endc))                   ; nutrient_state%veg_trunc_col           (:)     = nan
    allocate(nutrient_state%soil_trunc_vr_col       (begc:endc,1:nlevdecomp_full)) ; nutrient_state%soil_trunc_vr_col       (:,:)   = nan

    allocate(nutrient_state%cwd_col                 (begc:endc))                   ; nutrient_state%cwd_col                 (:)     = nan
    allocate(nutrient_state%decomp_pools_col        (begc:endc,1:ndecomp_pools))   ; nutrient_state%decomp_pools_col        (:,:)   = nan
    allocate(nutrient_state%decomp_pools_1m_col     (begc:endc,1:ndecomp_pools))   ; nutrient_state%decomp_pools_1m_col     (:,:)   = nan
    allocate(nutrient_state%dyn_bal_adjustments_col (begc:endc))                   ; nutrient_state%dyn_bal_adjustments_col (:)     = nan

    allocate(nutrient_state%totabg_col              (begc:endc))                   ; nutrient_state%totabg_col              (:)     = nan
    allocate(nutrient_state%totblg_col              (begc:endc))                   ; nutrient_state%totblg_col              (:)     = nan
    allocate(nutrient_state%totcol_col              (begc:endc))                   ; nutrient_state%totcol_col              (:)     = nan
    allocate(nutrient_state%totecosys_col           (begc:endc))                   ; nutrient_state%totecosys_col           (:)     = nan
    allocate(nutrient_state%totlit_1m_col           (begc:endc))                   ; nutrient_state%totlit_1m_col           (:)     = nan
    allocate(nutrient_state%totlit_col              (begc:endc))                   ; nutrient_state%totlit_col              (:)     = nan
    allocate(nutrient_state%totpft_col              (begc:endc))                   ; nutrient_state%totpft_col              (:)     = nan
    allocate(nutrient_state%totprod_col             (begc:endc))                   ; nutrient_state%totprod_col             (:)     = nan
    allocate(nutrient_state%totsom_col              (begc:endc))                   ; nutrient_state%totsom_col              (:)     = nan
    allocate(nutrient_state%totsom_1m_col           (begc:endc))                   ; nutrient_state%totsom_1m_col           (:)     = nan
    allocate(nutrient_state%totveg_col              (begc:endc))                   ; nutrient_state%totveg_col              (:)     = nan
    allocate(nutrient_state%totveg_col              (begc:endc))                   ; nutrient_state%totveg_col              (:)     = nan
    allocate(nutrient_state%totveg_abg_col          (begc:endc))                   ; nutrient_state%totveg_abg_col          (:)     = nan

  end subroutine NutrientStateInitAllocate

  !------------------------------------------------------------------------
  subroutine NutrientStateInitHistory(nutrient_state, bounds)
    !
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    !
    implicit none
    !
    ! !ARGUMENTS:
    class (nutrientstate_type)    :: nutrient_state
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: l
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    character(len=10) :: unit
    character(len=3)  :: name
    character(len=4)  :: name_prefix
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    write(unit,*) 'g' // trim(nutrient_state%name) // '/m^2'
    unit = trim(unit)

    name = nutrient_state%name
    name_prefix = nutrient_state%history_name_prefix

    nutrient_state%deadcroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'DEADCROOT' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead coarse root', &
         ptr_patch=nutrient_state%deadcroot_patch)

    nutrient_state%deadcroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'DEADCROOT' // trim(name) // '_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead coarse root storage', &
         ptr_patch=nutrient_state%deadcroot_storage_patch,  default='inactive')

    nutrient_state%deadcroot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'DEADCROOT'// trim(name) // '_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead coarse root transfer', &
         ptr_patch=nutrient_state%deadcroot_xfer_patch, default='inactive')

    nutrient_state%deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'DEADSTEM' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead stem', &
         ptr_patch=nutrient_state%deadstem_patch)

    nutrient_state%deadstem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'DEADSTEM' // trim(name) // '_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead stem storage', &
         ptr_patch=nutrient_state%deadstem_storage_patch, default='inactive')

    nutrient_state%deadstem_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'DEADSTEM' // trim(name) // '_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead stem transfer', &
         ptr_patch=nutrient_state%deadstem_xfer_patch, default='inactive')

    nutrient_state%froot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'FROOT' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' fine root', &
         ptr_patch=nutrient_state%froot_patch)

    nutrient_state%froot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'FROOT' // trim(name) // '_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' fine root storage', &
         ptr_patch=nutrient_state%froot_storage_patch, default='inactive')

    nutrient_state%froot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'FROOT' // trim(name) // '_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' fine root transfer', &
         ptr_patch=nutrient_state%froot_xfer_patch, default='inactive')

    if (use_crop) then
       nutrient_state%grain_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name_prefix) // 'GRAIN' // trim(name), units=trim(unit), &
            avgflag='A', long_name = trim(name) // ' grain', &
            ptr_patch=nutrient_state%grain_patch, default='inactive')

       nutrient_state%grain_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name_prefix) // 'GRAIN' // trim(name) // '_STORAGE', units=trim(unit), &
            avgflag='A', long_name = trim(name) // ' grain storage', &
            ptr_patch=nutrient_state%grain_storage_patch, default='inactive')

       nutrient_state%grain_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name_prefix) // 'GRAIN' // trim(name) //'_XFER', units=trim(unit), &
            avgflag='A', long_name = trim(name) // ' grain transfer', &
            ptr_patch=nutrient_state%grain_xfer_patch, default='inactive')

       nutrient_state%cropseed_deficit_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name_prefix) // 'CROPSEED' // trim(name) // '_DEFICIT', units='gN/m^2', &
            avgflag='A', long_name=trim(name) // ' used for crop seed that needs to be repaid', &
            ptr_patch=nutrient_state%cropseed_deficit_patch)
    end if

    nutrient_state%leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'LEAF' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' leaf', &
         ptr_patch=nutrient_state%leaf_patch)

    nutrient_state%leaf_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'LEAF' // trim(name) //'_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' leaf storage', &
         ptr_patch=nutrient_state%leaf_storage_patch, default='inactive')

    nutrient_state%leaf_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'LEAF' // trim(name) //'_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' leaf transfer', &
         ptr_patch=nutrient_state%leaf_xfer_patch, default='inactive')

    nutrient_state%livestem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'LIVESTEM' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live stem', &
         ptr_patch=nutrient_state%livestem_patch)

    nutrient_state%livestem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'LIVESTEM' // trim(name) //'_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live stem storage', &
         ptr_patch=nutrient_state%livestem_storage_patch, default='inactive')

    nutrient_state%livestem_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'LIVESTEMC' // trim(name) //'_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live stem transfer', &
         ptr_patch=nutrient_state%livestem_xfer_patch, default='inactive')

    nutrient_state%livecroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'LIVECROOT' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live coarse root', &
         ptr_patch=nutrient_state%livecroot_patch)

    nutrient_state%livecroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // trim(name) // '_LIVECROOT' // trim(name) //'_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live coarse root storage', &
         ptr_patch=nutrient_state%livecroot_storage_patch, default='inactive')

    nutrient_state%livecroot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // trim(name) // '_LIVECROOT' // trim(name) //'_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live coarse root transfer', &
         ptr_patch=nutrient_state%livecroot_xfer_patch, default='inactive')

    nutrient_state%decomp_pools_vr_col(begc:endc,:,:) = spval
    do l = 1, ndecomp_pools
       if ( nlevdecomp_full > 1 ) then
          data2dptr => nutrient_state%decomp_pools_vr_col(:,:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))// &
               trim(name) // '_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' ' &
               // trim(name) // ' (vertically resolved)'
          call hist_addfld2d (fname=trim(name_prefix) // fieldname, units='g' // trim(name) // '/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif

    end do

    nutrient_state%veg_trunc_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'COL_' // trim(name) // 'TRUNC', units=trim(unit),  &
         avgflag='A', long_name = trim(name) // ' column-level sink for truncation', &
         ptr_col=nutrient_state%veg_trunc_col, default='inactive')

    nutrient_state%veg_trunc_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'PFT_' // trim(name) // 'TRUNC', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' patch-level sink for truncation', &
         ptr_patch=nutrient_state%veg_trunc_patch, default='inactive')

    nutrient_state%seed_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'SEED' // trim(name) // '_COL', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' pool for seeding new Patches', &
         ptr_col=nutrient_state%seed_col, default='inactive')

    nutrient_state%seed_grc(begg:endg) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'SEED' // trim(name) // '_GRC', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' pool for seeding new PFTs via dynamic landcover', &
         ptr_gcell=nutrient_state%seed_grc)

    nutrient_state%prod1_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'PROD1' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' 1-yr crop product', &
         ptr_col=nutrient_state%prod1_col, default='inactive')

    nutrient_state%prod10_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'PROD10' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' 10-yr wood product', &
         ptr_col=nutrient_state%prod10_col, default='inactive')

    nutrient_state%prod100_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'PROD100' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' 100-yr wood product', &
         ptr_col=nutrient_state%prod100_col, default='inactive')

    do l = 1, ndecomp_pools

       if ( nlevdecomp_full > 1 ) then
          data1dptr => nutrient_state%decomp_pools_1m_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))// trim(name) // '_1m'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))// trim(name) // ' to 1 meter'
          call hist_addfld1d (fname=trim(name_prefix) // fieldname, units=trim(unit), &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default = 'inactive')
       endif

       data1dptr => nutrient_state%decomp_pools_col(:,l)
       fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l)) // trim(name)
       longname =  trim(decomp_cascade_con%decomp_pool_name_history(l)) // ' ' // trim(name)
       call hist_addfld1d (fname=trim(name_prefix) // fieldname, units=trim(unit), &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

    end do

    nutrient_state%dispveg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'DISPVEG' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' displayed veg carbon, excluding storage and cpool', &
         ptr_patch=nutrient_state%dispveg_patch)

    nutrient_state%pool_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // trim(name) // 'POOL', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' temporary photosynthate pool', &
         ptr_patch=nutrient_state%pool_patch)

    nutrient_state%storveg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'STORVEG' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' stored vegetation, excluding cpool', &
         ptr_patch=nutrient_state%storveg_patch)

    nutrient_state%totcol_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'TOTCOL' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total column incl veg and pool but excl product pools', &
         ptr_col=nutrient_state%totcol_col)

    nutrient_state%totecosys_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'TOTECOSYS' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total ecosystem, incl veg but excl cpool but excl product pools', &
         ptr_col=nutrient_state%totecosys_col)

    if ( nlevdecomp_full > 1 ) then
       nutrient_state%totlit_1m_col(begc:endc) = spval
       call hist_addfld1d (fname=trim(name_prefix) // 'TOTLIT' // trim(name) // '_1m', units=trim(unit), &
            avgflag='A', long_name=trim(name) // ' total litter to 1 meter depth', &
            ptr_col=nutrient_state%totlit_1m_col)
    end if

    nutrient_state%totlit_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'TOTLIT' // trim(name), units=trim(unit), &
         avgflag='A', long_name=trim(name) // ' total litter', &
         ptr_col=nutrient_state%totlit_col)

    nutrient_state%totpft_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'TOTPFT' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total patch-level, including pool', &
         ptr_patch=nutrient_state%totpft_patch)

    nutrient_state%totprod_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'TOTPROD' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total wood product', &
         ptr_col=nutrient_state%totprod_col, default='inactive')

    if ( nlevdecomp_full > 1 ) then
       nutrient_state%totsom_1m_col(begc:endc) = spval
       call hist_addfld1d (fname=trim(name_prefix) // 'TOTSOM' // trim(name) // '_1m', units=trim(unit), &
            avgflag='A', long_name=trim(name) // ' total soil organic matter to 1 meter depth', &
            ptr_col=nutrient_state%totsom_1m_col)
    end if

    nutrient_state%totsom_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'TOTSOM' // trim(name), units=trim(unit), &
         avgflag='A', long_name=trim(name) // ' total soil organic matter', &
         ptr_col=nutrient_state%totsom_col)

    nutrient_state%totveg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) //  'TOTVEG' // trim(name), units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total vegetation, excluding pool', &
         ptr_patch=nutrient_state%totveg_patch)

    nutrient_state%totveg_abg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix) // 'TOTVEG' // trim(name) // '_ABG', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total aboveground vegetation, excluding pool', &
         ptr_patch=nutrient_state%totveg_abg_patch)

  end subroutine NutrientStateInitHistory

  !-----------------------------------------------------------------------
  subroutine NutrientStateDynamicPatchAdjustments( &
       this,                      &
       bounds,                    &
       num_filterp_with_inactive, &
       filterp_with_inactive,     &
       prior_weights,             &
       patch_state_updater,       &
       species_type,              &
       dwt_leaf_seed,             &
       dwt_deadstem_seed,         &
       conv_flux,                 &
       dwt_froot_to_litter,       &
       dwt_livecroot_to_litter,   &
       dwt_deadcroot_to_litter,   &
       prod10_flux,               &
       prod100_flux,              &
       crop_product_flux,         &
       dwt_pool_seed,             &
       pool_seed_param,           &
       pool_seed_patch            &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use CNComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    class(nutrientstate_type)      , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    integer                        , intent(in)    :: species_type
    real(r8)                       , intent(inout) :: dwt_leaf_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstem_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_flux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_froot_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecroot_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcroot_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_flux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_flux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_flux       (bounds%begp:)
    real(r8), optional             , intent(inout) :: dwt_pool_seed           (bounds%begp:)
    real(r8), optional             , intent(in)    :: pool_seed_param
    real(r8), optional             , intent(inout) :: pool_seed_patch         (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leaf_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leaf_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leaf_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstem_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: wood_product_flux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstem_patch_temp     (bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'DynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leaf_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstem_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_pool_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_flux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_froot_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecroot_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcroot_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_flux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_flux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_flux       ) == (/endp/)), errMsg(__FILE__, __LINE__))

    if (present(pool_seed_patch)) then
       SHR_ASSERT_ALL((ubound(pool_seed_patch       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    endif
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = species_type                        , &
         leaf_patch                 = this%leaf_patch(begp:endp)          , &
         leaf_storage_patch         = this%leaf_storage_patch(begp:endp)  , &
         leaf_xfer_patch            = this%leaf_xfer_patch(begp:endp)     , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leaf_patch(begp:endp)          , &
         seed_leaf_storage_patch    = seed_leaf_storage_patch(begp:endp)  , &
         seed_leaf_xfer_patch       = seed_leaf_xfer_patch(begp:endp)     , &
         seed_deadstem_patch        = seed_deadstem_patch(begp:endp)     , &
         pool_seed_param            = pool_seed_param                     , &
         pool_seed_patch            = pool_seed_patch(begp:endp))

    ! 1) LEAF_PATCH
    call patch_state_updater%update_patch_state(            &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = this%leaf_patch   (begp:endp) , &
         flux_out_grc_area = conv_flux       (begp:endp) , &
         seed              = seed_leaf_patch (begp:endp) , &
         seed_addition     = dwt_leaf_seed   (begp:endp))

    ! 2) LEAF_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                    &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = this%leaf_storage_patch   (begp:endp) , &
         flux_out_grc_area = conv_flux               (begp:endp) , &
         seed              = seed_leaf_storage_patch (begp:endp) , &
         seed_addition     = dwt_leaf_seed           (begp:endp))

    ! 3) LEAF_XFER_PATCH
    call patch_state_updater%update_patch_state( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = this%leaf_xfer_patch   (begp:endp), &
         flux_out_grc_area = conv_flux            (begp:endp), &
         seed              = seed_leaf_xfer_patch (begp:endp), &
         seed_addition     = dwt_leaf_seed        (begp:endp))

    ! 4) FROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_patch(begp:endp)             , &
         flux_out_col_area = dwt_froot_to_litter(begp:endp))

    ! 5) FROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 6) FROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 7) LIVESTEMN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_patch(begp:endp)          , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 8) LIVESTEMN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 9) LIVESTEMN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 10) PROD10_FLUX
    wood_product_flux(begp:endp)      = 0._r8
    deadstem_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstem_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_flux            (begp:endp) , &
         flux2_out                  = wood_product_flux      (begp:endp) , &
         seed                       = seed_deadstem_patch    (begp:endp) )

    ! 11) PROD100_FLUX
    wood_product_flux(begp:endp)      = 0._r8
    deadstem_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstem_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_flux           (begp:endp) , &
         flux2_out                  = wood_product_flux      (begp:endp) , &
         seed                       = seed_deadstem_patch    (begp:endp))

    ! 12) DEADSTEM_PATCH
    wood_product_flux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = this%deadstem_patch   (begp:endp)    , &
         flux1_out                  = conv_flux           (begp:endp)    , &
         flux2_out                  = wood_product_flux   (begp:endp)    , &
         seed                       = seed_deadstem_patch (begp:endp)    , &
         seed_addition              = dwt_deadstem_seed   (begp:endp))

    ! 13) DEADSTEM_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 14) DEADSTEM_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 15) LIVECROOTN_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_livecroot_to_litter(begp:endp))

    ! 16) LIVECROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 17) LIVECROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 18) DEADCROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_deadcroot_to_litter(begp:endp))

    ! 19) DEADCROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 21) NTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%veg_trunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 22) CPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%pool_patch(begp:endp)              , &
         flux_out_grc_area = conv_flux(begp:endp))

    ! 23) DISPVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%dispveg_patch(begp:endp))

    ! 24) STORVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%storveg_patch(begp:endp))

    ! 25) TOTVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totveg_patch(begp:endp))

    ! 26) TOTPFTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totpft_patch(begp:endp))

    ! 27) CROPSEED_DEFICIT
    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call patch_state_updater%update_patch_state(         &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = this%cropseed_deficit_patch(begp:endp) , &
            flux_out_grc_area = conv_flux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_froot_to_litter(p)     = -1._r8 * dwt_froot_to_litter(p)
       dwt_livecroot_to_litter(p) = -1._r8 * dwt_livecroot_to_litter(p)
       dwt_deadcroot_to_litter(p) = -1._r8 * dwt_deadcroot_to_litter(p)
    end do

  end subroutine NutrientStateDynamicPatchAdjustments

  !-----------------------------------------------------------------------
  subroutine NutrientStateRestart ( this,  bounds, ncid, flag)
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan
    use shr_const_mod    , only : SHR_CONST_PDB
    use clm_varctl       , only : spinup_mortality_factor, spinup_state
    use tracer_varcon    , only : is_active_betr_bgc
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (nutrientstate_type)           :: this
    type(bounds_type)    , intent(in)    :: bounds 
    type(file_desc_t)    , intent(inout) :: ncid       ! netcdf id
    character(len=*)     , intent(in)    :: flag       !'read' or 'write'
                                                       !
                                                       ! !LOCAL VARIABLES:
    logical                              :: readvar
    integer                              :: k
    character(len=128)                   :: varname    ! temporary
    real(r8), pointer                    :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer                    :: ptr1d(:)   ! temp. pointers for slicing larger arrays

    call restartvar(ncid=ncid, flag=flag, varname='leaf' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leaf' // trim(this%restart_name) // '_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leaf' // trim(this%restart_name) // '_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='froot' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='froot' // trim(this%restart_name) // '_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='froot' // trim(this%restart_name) // '_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestem' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestem' // trim(this%restart_name) // '_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestem' // trim(this%restart_name) // '_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstem' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstem' // trim(this%restart_name) // '_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstem' // trim(this%restart_name) // '_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecroot' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecroot' // trim(this%restart_name) // '_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecroot' // trim(this%restart_name) // '_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcroot' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcroot' // trim(this%restart_name) // '_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcroot' // trim(this%restart_name) // '_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname=trim(this%restart_name) // 'pool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_' //trim(this%restart_name) // 'trunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='totveg' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totveg_patch)

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grain' // trim(this%restart_name), xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grain' // trim(this%restart_name) // '_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grain' // trim(this%restart_name) // '_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_xfer_patch)

       call restartvar(ncid=ncid, flag=flag, varname='cropseed' // trim(this%restart_name) // '_deficit', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cropseed_deficit_patch)
    end if

    do k = 1, ndecomp_pools
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k)) // trim(this%restart_name)
       if (use_vertsoilc) then
          ptr2d => this%decomp_pools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%decomp_pools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end do

    if (is_active_betr_bgc)then
       call restartvar(ncid=ncid, flag=flag, varname='totblg' // trim(this%restart_name), xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totblg_col)

       call restartvar(ncid=ncid, flag=flag, varname='cwd' // trim(this%restart_name), xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cwd_col)
    endif

    if (use_vertsoilc) then
       ptr2d => this%soil_trunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname='col_' // trim(this%restart_name) // 'trunc_vr', xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%soil_trunc_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname='col_' // trim(this%restart_name) // 'trunc', xtype=ncd_double,  &
            dim1name='column', long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
            errMsg(__FILE__, __LINE__))
    end if

    call restartvar(ncid=ncid, flag=flag, varname='seed' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%seed_col) 

    call restartvar(ncid=ncid, flag=flag, varname='totlit' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totlit_col) 

    call restartvar(ncid=ncid, flag=flag, varname='totcol' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod10' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod100' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod1' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod1_col)

    call restartvar(ncid=ncid, flag=flag, varname='totsom' // trim(this%restart_name), xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totsom_col)

  end subroutine NutrientStateRestart

  !-----------------------------------------------------------------------
  subroutine NutrientStatePatchSummary(this, bounds, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform patch-level summary calculations
    !
    ! !USES:
    use clm_varctl       , only: iulog
    use clm_time_manager , only: get_step_size
    use clm_varcon       , only: secspday
    use clm_varpar       , only: nlevdecomp, ndecomp_pools, nlevdecomp_full
    !
    ! !ARGUMENTS:
    class(nutrientstate_type)              :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    integer  :: p, fp
    !-----------------------------------------------------------------------

    ! calculate patch -level summary of carbon state

    if (use_fates) return

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation pools
       this%dispveg_patch(p) =        &
            this%leaf_patch(p)      + &
            this%froot_patch(p)     + &
            this%livestem_patch(p)  + &
            this%deadstem_patch(p)  + &
            this%livecroot_patch(p) + &
            this%deadcroot_patch(p)

       ! stored vegetation pools
       this%storveg_patch(p) =                &
            this%pool_patch(p)              + &
            this%leaf_storage_patch(p)      + &
            this%froot_storage_patch(p)     + &
            this%livestem_storage_patch(p)  + &
            this%deadstem_storage_patch(p)  + &
            this%livecroot_storage_patch(p) + &
            this%deadcroot_storage_patch(p) + &
            this%leaf_xfer_patch(p)         + &
            this%froot_xfer_patch(p)        + &
            this%livestem_xfer_patch(p)     + &
            this%deadstem_xfer_patch(p)     + &
            this%livecroot_xfer_patch(p)    + &
            this%deadcroot_xfer_patch(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%storveg_patch(p) =            &
               this%storveg_patch(p)       + &
               this%grain_storage_patch(p) + &
               this%grain_xfer_patch(p)

          this%dispveg_patch(p) =            &
               this%dispveg_patch(p)       + &
               this%grain_patch(p)
       end if

       this%totveg_abg_patch(p) = &
               this%leaf_patch(p)              + &
               this%leaf_storage_patch(p)      + &
               this%leaf_xfer_patch(p)         + &
               this%livestem_patch(p)          + &
               this%livestem_storage_patch(p)  + &
               this%livestem_xfer_patch(p)     + &
               this%deadstem_patch(p)          + &
               this%deadstem_storage_patch(p)  + &
               this%deadstem_xfer_patch(p)

    end do

  end subroutine NutrientStatePatchSummary

  !-----------------------------------------------------------------------
  subroutine NutrientStateColumnSummary(this, bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform patch and column-level carbon summary calculations
    !
    ! !USES:
    use clm_varctl       , only: iulog
    use clm_time_manager , only: get_step_size
    use clm_varcon       , only: secspday
    use clm_varpar       , only: nlevdecomp, ndecomp_pools, nlevdecomp_full
    use clm_varcon       , only: dzsoi_decomp, zisoi, zsoi
    use clm_varctl       , only: use_pflotran, pf_cmode
    !
    ! !ARGUMENTS:
    class(nutrientstate_type)      :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer                        :: c,p,j,k,l       ! indices
    integer                        :: fp,fc           ! lake filter indices
    real(r8)                       :: maxdepth        ! depth to integrate soil variables
    integer                        :: nlev
    !-----------------------------------------------------------------------

    ! calculate patch -level summary of carbon state

    nlev = nlevdecomp
    if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

    ! initialize
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_pools_col(c,l) = 0._r8
       end do
    end do
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totlit_col(c)    = 0._r8
       this%totsom_col(c)    = 0._r8
       this%cwd_col(c)       = 0._r8
       this%veg_trunc_col(c) = 0._r8
    end do

    ! vertically integrate each of the decomposing pools
    do l = 1, ndecomp_pools
       do j = 1, nlev
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_pools_col(c,l) = &
                  this%decomp_pools_col(c,l) + &
                  this%decomp_pools_vr_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
    end do

    if ( nlevdecomp > 1) then

       ! vertically integrate each of the decomposing pools to 1 meter
       maxdepth = 1._r8
       do l = 1, ndecomp_pools
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_pools_1m_col(c,l) = 0._r8
          end do
       end do
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp
             if ( zisoi(j) <= maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   this%decomp_pools_1m_col(c,l) = &
                        this%decomp_pools_1m_col(c,l) + &
                        this%decomp_pools_vr_col(c,j,l) * dzsoi_decomp(j)
                end do
             elseif ( zisoi(j-1) < maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   this%decomp_pools_1m_col(c,l) = &
                        this%decomp_pools_1m_col(c,l) + &
                        this%decomp_pools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
                end do
             endif
          end do
       end do

       ! total litter in the top 1 meter
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%totlit_1m_col(c) = 0._r8
          this%totsom_1m_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%totlit_1m_col(c) = &
                     this%totlit_1m_col(c) + &
                     this%decomp_pools_1m_col(c,l)
             end do
          endif
       end do

       ! total soil organic matter in the top 1 meter
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_soil(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%totsom_1m_col(c) = &
                     this%totsom_1m_col(c) + &
                     this%decomp_pools_1m_col(c,l)
             end do
          end if
       end do

    endif

    do l = 1, ndecomp_pools
       ! total litter
       if ( decomp_cascade_con%is_litter(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totlit_col(c) = &
                  this%totlit_col(c) + &
                  this%decomp_pools_col(c,l)
          end do
       endif

       ! total soil organic matter
       if ( decomp_cascade_con%is_soil(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totsom_col(c) = &
                  this%totsom_col(c) + &
                  this%decomp_pools_col(c,l)
          end do
       end if

       ! coarse woody debris
       if ( decomp_cascade_con%is_cwd(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%cwd_col(c) = &
                  this%cwd_col(c) + &
                  this%decomp_pools_col(c,l)
          end do
       end if
    end do

    ! truncation
    do j = 1, nlev
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%veg_trunc_col(c) = &
               this%veg_trunc_col(c) + &
               this%soil_trunc_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

  end subroutine NutrientStateColumnSummary

end module NutrientStateType
