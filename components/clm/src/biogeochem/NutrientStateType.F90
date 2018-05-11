module NutrientStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use decompMod              , only : bounds_type
  use clm_varpar             , only : nlevdecomp_full, ndecomp_pools
  use clm_varctl             , only : use_fates
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

  public :: NutrientStateInitAllocate

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

end module NutrientStateType
