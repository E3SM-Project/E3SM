module NutrientStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use decompMod              , only : bounds_type
  use clm_varpar             , only : nlevdecomp_full, ndecomp_pools
  use clm_varctl             , only : use_fates, use_crop
  use clm_varcon             , only : spval, ispval
  use CNDecompCascadeConType , only : decomp_cascade_con
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  type, public :: nutrientstate_type

     integer           :: species                         ! C12, C13, C14, N, P

     character(len=3)  :: name

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
            NutrientStateInitHistory

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

    nutrient_state%deadcroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_DEADCROOT', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead coarse root', &
         ptr_patch=nutrient_state%deadcroot_patch)

    nutrient_state%deadcroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_DEADCROOT_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead coarse root storage', &
         ptr_patch=nutrient_state%deadcroot_storage_patch,  default='inactive')

    nutrient_state%deadcroot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_DEADCROOT_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead coarse root transfer', &
         ptr_patch=nutrient_state%deadcroot_xfer_patch, default='inactive')

    nutrient_state%deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_DEADSTEM', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead stem', &
         ptr_patch=nutrient_state%deadstem_patch)

    nutrient_state%deadstem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_DEADSTEM_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead stem storage', &
         ptr_patch=nutrient_state%deadstem_storage_patch, default='inactive')

    nutrient_state%deadstem_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_DEADSTEM_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' dead stem transfer', &
         ptr_patch=nutrient_state%deadstem_xfer_patch, default='inactive')

    nutrient_state%froot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_FROOT', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' fine root', &
         ptr_patch=nutrient_state%froot_patch)

    nutrient_state%froot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_FROOT_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' fine root storage', &
         ptr_patch=nutrient_state%froot_storage_patch, default='inactive')

    nutrient_state%froot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_FROOT_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' fine root transfer', &
         ptr_patch=nutrient_state%froot_xfer_patch, default='inactive')

    if (use_crop) then
       nutrient_state%grain_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name) // '_GRAIN', units=trim(unit), &
            avgflag='A', long_name = trim(name) // ' grain', &
            ptr_patch=nutrient_state%grain_patch, default='inactive')

       nutrient_state%grain_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name) // '_GRAIN_STORAGE', units=trim(unit), &
            avgflag='A', long_name = trim(name) // ' grain storage', &
            ptr_patch=nutrient_state%grain_storage_patch, default='inactive')

       nutrient_state%grain_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name) // '_GRAIN_XFER', units=trim(unit), &
            avgflag='A', long_name = trim(name) // ' grain transfer', &
            ptr_patch=nutrient_state%grain_xfer_patch, default='inactive')
    end if

    nutrient_state%leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LEAF', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' leaf', &
         ptr_patch=nutrient_state%leaf_patch)

    nutrient_state%leaf_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LEAF_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' leaf storage', &
         ptr_patch=nutrient_state%leaf_storage_patch, default='inactive')

    nutrient_state%leaf_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LEAF_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' leaf transfer', &
         ptr_patch=nutrient_state%leaf_xfer_patch, default='inactive')

    nutrient_state%livestem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LIVESTEM', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live stem', &
         ptr_patch=nutrient_state%livestem_patch)

    nutrient_state%livestem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LIVESTEM_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live stem storage', &
         ptr_patch=nutrient_state%livestem_storage_patch, default='inactive')

    nutrient_state%livestem_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LIVESTEMC_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live stem transfer', &
         ptr_patch=nutrient_state%livestem_xfer_patch, default='inactive')

    nutrient_state%livecroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LIVECROOT', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live coarse root', &
         ptr_patch=nutrient_state%livecroot_patch)

    nutrient_state%livecroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LIVECROOT_STORAGE', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live coarse root storage', &
         ptr_patch=nutrient_state%livecroot_storage_patch, default='inactive')

    nutrient_state%livecroot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_LIVECROOT_XFER', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' live coarse root transfer', &
         ptr_patch=nutrient_state%livecroot_xfer_patch, default='inactive')

    nutrient_state%decomp_pools_vr_col(begc:endc,:,:) = spval
    do l = 1, ndecomp_pools
       if ( nlevdecomp_full > 1 ) then
          data2dptr => nutrient_state%decomp_pools_vr_col(:,:,l)
          fieldname = trim(name) // '_'//trim(decomp_cascade_con%decomp_pool_name_history(l))// &
               trim(name) // '_vr'
          longname =  trim(name) // ' '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' ' &
               // trim(name) // ' (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='g' // trim(name) // '/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif

    end do

    nutrient_state%veg_trunc_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_COL_CTRUN', units=trim(unit),  &
         avgflag='A', long_name = trim(name) // ' column-level sink for truncation', &
         ptr_col=nutrient_state%veg_trunc_col, default='inactive')

    nutrient_state%veg_trunc_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_PFT_CTRUN', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' patch-level sink for truncation', &
         ptr_patch=nutrient_state%veg_trunc_patch, default='inactive')

    nutrient_state%seed_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_SEED_COL', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' pool for seeding new Patches', &
         ptr_col=nutrient_state%seed_col, default='inactive')

    nutrient_state%seed_grc(begg:endg) = spval
    call hist_addfld1d (fname=trim(name) // '_SEED_GRC', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' pool for seeding new PFTs via dynamic landcover', &
         ptr_gcell=nutrient_state%seed_grc)

    nutrient_state%prod1_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_PROD1', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' 1-yr crop product', &
         ptr_col=nutrient_state%prod1_col, default='inactive')

    nutrient_state%prod10_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_PROD10', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' 10-yr wood product', &
         ptr_col=nutrient_state%prod10_col, default='inactive')

    nutrient_state%prod100_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_PROD100', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' 100-yr wood product', &
         ptr_col=nutrient_state%prod100_col, default='inactive')

    do l = 1, ndecomp_pools

       if ( nlevdecomp_full > 1 ) then
          data1dptr => nutrient_state%decomp_pools_1m_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))// trim(name) // '_1m'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))// trim(name) // ' to 1 meter'
          call hist_addfld1d (fname=fieldname, units=trim(unit), &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default = 'inactive')
       endif

       data1dptr => nutrient_state%decomp_pools_col(:,l)
       fieldname = trim(name) // '_' // trim(decomp_cascade_con%decomp_pool_name_history(l)) // trim(name)
       longname =  trim(name) // ' ' // trim(decomp_cascade_con%decomp_pool_name_history(l)) // ' ' // trim(name)
       call hist_addfld1d (fname=fieldname, units=trim(unit), &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

    end do

    nutrient_state%dispveg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_DISPVEG', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' displayed veg carbon, excluding storage and cpool', &
         ptr_patch=nutrient_state%dispveg_patch)

    nutrient_state%pool_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_POOL', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' temporary photosynthate pool', &
         ptr_patch=nutrient_state%pool_patch)

    nutrient_state%storveg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_STORVEG', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' stored vegetation, excluding cpool', &
         ptr_patch=nutrient_state%storveg_patch)

    nutrient_state%totcol_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTCOL', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total column incl veg and pool but excl product pools', &
         ptr_col=nutrient_state%totcol_col)

    nutrient_state%totecosys_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTECOSYS', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total ecosystem, incl veg but excl cpool but excl product pools', &
         ptr_col=nutrient_state%totecosys_col)

    if ( nlevdecomp_full > 1 ) then
       nutrient_state%totlit_1m_col(begc:endc) = spval
       call hist_addfld1d (fname=trim(name) // '_TOTLIT_1m', units=trim(unit), &
            avgflag='A', long_name=trim(name) // ' total litter to 1 meter depth', &
            ptr_col=nutrient_state%totlit_1m_col)
    end if

    nutrient_state%totlit_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTLIT', units=trim(unit), &
         avgflag='A', long_name=trim(name) // ' total litter', &
         ptr_col=nutrient_state%totlit_col)

    nutrient_state%totpft_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTPFT', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total patch-level, including pool', &
         ptr_patch=nutrient_state%totpft_patch)

    nutrient_state%totprod_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTPROD', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total wood product', &
         ptr_col=nutrient_state%totprod_col, default='inactive')

    if ( nlevdecomp_full > 1 ) then
       nutrient_state%totsom_1m_col(begc:endc) = spval
       call hist_addfld1d (fname=trim(name) // '_TOTSOM_1m', units=trim(unit), &
            avgflag='A', long_name=trim(name) // ' total soil organic matter to 1 meter depth', &
            ptr_col=nutrient_state%totsom_1m_col)
    end if

    nutrient_state%totsom_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTSOM', units=trim(unit), &
         avgflag='A', long_name=trim(name) // ' total soil organic matter', &
         ptr_col=nutrient_state%totsom_col)

    nutrient_state%totveg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTVEG', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total vegetation, excluding pool', &
         ptr_patch=nutrient_state%totveg_patch)

    nutrient_state%totveg_abg_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name) // '_TOTVEG_ABG', units=trim(unit), &
         avgflag='A', long_name = trim(name) // ' total aboveground vegetation, excluding pool', &
         ptr_patch=nutrient_state%totveg_abg_patch)

  end subroutine NutrientStateInitHistory

end module NutrientStateType
