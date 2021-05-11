module CNCarbonStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, crop_prog, nlevdecomp
  !use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi, zsoi
  use landunit_varcon        , only : istcrop 
  use clm_varctl             , only : iulog, use_vertsoilc, use_cndv, spinup_state 
  use decompMod              , only : bounds_type
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use clm_varctl             , only : nu_com, use_fates, use_crop

  ! bgc interface & pflotran
  use clm_varctl             , only : use_elm_interface, use_pflotran, pf_cmode
  
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: carbonstate_type
     
     integer :: species  ! c12, c13, c14

     real(r8), pointer :: grainc_patch             (:)     ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainc_storage_patch     (:)     ! (gC/m2) grain C storage (crop model)
     real(r8), pointer :: grainc_xfer_patch        (:)     ! (gC/m2) grain C transfer (crop model)
     real(r8), pointer :: leafc_patch              (:)     ! (gC/m2) leaf C
     real(r8), pointer :: leafc_storage_patch      (:)     ! (gC/m2) leaf C storage
     real(r8), pointer :: leafc_xfer_patch         (:)     ! (gC/m2) leaf C transfer
     real(r8), pointer :: frootc_patch             (:)     ! (gC/m2) fine root C
     real(r8), pointer :: frootc_storage_patch     (:)     ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_patch        (:)     ! (gC/m2) fine root C transfer
     real(r8), pointer :: livestemc_patch          (:)     ! (gC/m2) live stem C
     real(r8), pointer :: livestemc_storage_patch  (:)     ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemc_xfer_patch     (:)     ! (gC/m2) live stem C transfer
     real(r8), pointer :: deadstemc_patch          (:)     ! (gC/m2) dead stem C
     real(r8), pointer :: deadstemc_storage_patch  (:)     ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemc_xfer_patch     (:)     ! (gC/m2) dead stem C transfer
     real(r8), pointer :: livecrootc_patch         (:)     ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootc_storage_patch (:)     ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootc_xfer_patch    (:)     ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: deadcrootc_patch         (:)     ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootc_storage_patch (:)     ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootc_xfer_patch    (:)     ! (gC/m2) dead coarse root C transfer
     real(r8), pointer :: gresp_storage_patch      (:)     ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch         (:)     ! (gC/m2) growth respiration transfer
     real(r8), pointer :: cpool_patch              (:)     ! (gC/m2) temporary photosynthate C pool
     real(r8), pointer :: xsmrpool_patch           (:)     ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: ctrunc_patch             (:)     ! (gC/m2) patch-level sink for C truncation
     real(r8), pointer :: woodc_patch              (:)     ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch           (:)     ! (gC/m2) ann max leaf C
     real(r8), pointer :: rootc_col                (:)     ! col (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: totvegc_col              (:)     ! col (gC/m2) column-level totvegc (fire)
     real(r8), pointer :: leafc_col                (:)     ! col (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col            (:)     ! col (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                (:)     ! col fuel avalability factor for Reg.C (0-1)
     real(r8), pointer :: fuelc_crop_col           (:)     ! col fuel avalability factor for Reg.A (0-1)

     ! all c pools involved in decomposition
     real(r8), pointer :: decomp_cpools_vr_col    (:,:,:)  ! col (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: ctrunc_vr_col           (:,:)    ! col (gC/m3) vertically-resolved column-level sink for C truncation
     

     ! pools for dynamic landcover
     real(r8), pointer :: cropseedc_deficit_patch  (:) ! (gC/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
     real(r8), pointer :: seedc_grc                (:)     ! (gC/m2) gridcell-level pool for seeding new PFTs via dynamic landcover
     real(r8), pointer :: frootc_col               (:)     ! col (gC/m2) column-level C pool for fine root
     real(r8), pointer :: seedc_col                (:)     ! col (gC/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod10c_col              (:)     ! col (gC/m2) wood product C pool, 10-year lifespan
     real(r8), pointer :: prod100c_col             (:)     ! col (gC/m2) wood product C pool, 100-year lifespan
     real(r8), pointer :: totprodc_col             (:)     ! col (gC/m2) total wood product C
     real(r8), pointer :: dyn_cbal_adjustments_col (:)     ! (gC/m2) adjustments to each column made in this timestep via dynamic column area adjustments
     ! pools for crop harvest
     real(r8), pointer :: prod1c_col               (:)     ! col (gC/m2) crop product C pool, 1-year lifespan

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegc_patch           (:)     ! (gC/m2) displayed veg carbon, excluding storage and cpool
     real(r8), pointer :: storvegc_patch           (:)     ! (gC/m2) stored vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_patch            (:)     ! (gC/m2) total vegetation carbon, excluding cpool
     real(r8), pointer :: totvegcc_col             (:)     ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)
     real(r8), pointer :: totpftc_patch            (:)     ! (gC/m2) total patch-level carbon, including cpool
     real(r8), pointer :: totpftc_col              (:)     ! (gC/m2) total patch-level carbon, including cpool averaged to column (p2c)
     real(r8), pointer :: decomp_cpools_1m_col     (:,:)   ! col (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
     real(r8), pointer :: decomp_cpools_col        (:,:)   ! col (gC/m2)  decomposing (litter, cwd, soil) c pools

     real(r8), pointer :: cwdc_col                 (:)     ! col (gC/m2) Diagnostic: coarse woody debris C
     real(r8), pointer :: ctrunc_col               (:)     ! col (gC/m2) column-level sink for C truncation
     real(r8), pointer :: totlitc_col              (:)     ! col (gC/m2) total litter carbon
     real(r8), pointer :: totsomc_col              (:)     ! col (gC/m2) total soil organic matter carbon
     real(r8), pointer :: totlitc_1m_col           (:)     ! col (gC/m2) total litter carbon to 1 meter
     real(r8), pointer :: totsomc_1m_col           (:)     ! col (gC/m2) total soil organic matter carbon to 1 meter
     real(r8), pointer :: totecosysc_col           (:)     ! col (gC/m2) total ecosystem carbon, incl veg but excl cpool
     real(r8), pointer :: totcolc_col              (:)     ! col (gC/m2) total column carbon, incl veg and cpool
     real(r8), pointer :: totabgc_col              (:)     ! col (gC/m2) total column above ground carbon, excluding som 
     real(r8), pointer :: totblgc_col              (:)     ! col (gc/m2) total column non veg carbon

     ! variables for above ground vegetation biomass
     real(r8), pointer :: totvegc_abg_patch            (:)     ! (gC/m2) total above vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_abg_col             (:)     ! (gC/m2) total above vegetation carbon, excluding cpool averaged to column (p2c)


     ! Balance checks
     real(r8), pointer :: begcb_patch              (:)     ! patch carbon mass, beginning of time step (gC/m**2)
     real(r8), pointer :: begcb_col                (:)     ! column carbon mass, beginning of time step (gC/m**2)
     real(r8), pointer :: begcb_grc                (:)     ! grid cell carbon mass, beginning of time step (gC/m**2)
     real(r8), pointer :: endcb_patch              (:)     ! patch carbon mass, end of time step (gC/m**2)
     real(r8), pointer :: endcb_col                (:)     ! column carbon mass, end of time step (gC/m**2)
     real(r8), pointer :: endcb_grc                (:)     ! grid cell carbon mass, end of time step (gC/m**2)
     real(r8), pointer :: errcb_patch              (:)     ! patch carbon balance error for the timestep (gC/m**2)
     real(r8), pointer :: errcb_col                (:)     ! column carbon balance error for the timestep (gC/m**2)
     real(r8), pointer :: errcb_grc                (:)     ! grid cell carbon balance error for the timestep (gC/m**2)
     
     real(r8), pointer :: totpftc_beg_col(:)
     real(r8), pointer :: cwdc_beg_col(:)
     real(r8), pointer :: totlitc_beg_col(:)
     real(r8), pointer :: totsomc_beg_col(:)
     
     real(r8), pointer :: totpftc_end_col(:)
     real(r8), pointer :: cwdc_end_col(:)
     real(r8), pointer :: totlitc_end_col(:)
     real(r8), pointer :: totsomc_end_col(:)
     real(r8), pointer :: decomp_som2c_vr_col(:,:)

   contains

     procedure , public  :: Init   
     procedure , private :: InitAllocate 

  end type carbonstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(carbonstate_type)                       :: this
    type(bounds_type)      , intent(in)           :: bounds  

    call this%InitAllocate ( bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    if ( .not. use_fates ) then
       allocate(this%leafc_patch              (begp :endp))                   ;     this%leafc_patch              (:)   = nan
       allocate(this%leafc_storage_patch      (begp :endp))                   ;     this%leafc_storage_patch      (:)   = nan
       allocate(this%leafc_xfer_patch         (begp :endp))                   ;     this%leafc_xfer_patch         (:)   = nan
       allocate(this%frootc_patch             (begp :endp))                   ;     this%frootc_patch             (:)   = nan
       allocate(this%frootc_storage_patch     (begp :endp))                   ;     this%frootc_storage_patch     (:)   = nan
       allocate(this%frootc_xfer_patch        (begp :endp))                   ;     this%frootc_xfer_patch        (:)   = nan
       allocate(this%livestemc_patch          (begp :endp))                   ;     this%livestemc_patch          (:)   = nan
       allocate(this%livestemc_storage_patch  (begp :endp))                   ;     this%livestemc_storage_patch  (:)   = nan
       allocate(this%livestemc_xfer_patch     (begp :endp))                   ;     this%livestemc_xfer_patch     (:)   = nan
       allocate(this%deadstemc_patch          (begp :endp))                   ;     this%deadstemc_patch          (:)   = nan
       allocate(this%deadstemc_storage_patch  (begp :endp))                   ;     this%deadstemc_storage_patch  (:)   = nan
       allocate(this%deadstemc_xfer_patch     (begp :endp))                   ;     this%deadstemc_xfer_patch     (:)   = nan
       allocate(this%livecrootc_patch         (begp :endp))                   ;     this%livecrootc_patch         (:)   = nan
       allocate(this%livecrootc_storage_patch (begp :endp))                   ;     this%livecrootc_storage_patch (:)   = nan
       allocate(this%livecrootc_xfer_patch    (begp :endp))                   ;     this%livecrootc_xfer_patch    (:)   = nan
       allocate(this%deadcrootc_patch         (begp :endp))                   ;     this%deadcrootc_patch         (:)   = nan
       allocate(this%deadcrootc_storage_patch (begp :endp))                   ;     this%deadcrootc_storage_patch (:)   = nan
       allocate(this%deadcrootc_xfer_patch    (begp :endp))                   ;     this%deadcrootc_xfer_patch    (:)   = nan
       allocate(this%gresp_storage_patch      (begp :endp))                   ;     this%gresp_storage_patch      (:)   = nan
       allocate(this%gresp_xfer_patch         (begp :endp))                   ;     this%gresp_xfer_patch         (:)   = nan
       allocate(this%cpool_patch              (begp :endp))                   ;     this%cpool_patch              (:)   = nan
       allocate(this%xsmrpool_patch           (begp :endp))                   ;     this%xsmrpool_patch           (:)   = nan
       allocate(this%ctrunc_patch             (begp :endp))                   ;     this%ctrunc_patch             (:)   = nan
       allocate(this%dispvegc_patch           (begp :endp))                   ;     this%dispvegc_patch           (:)   = nan
       allocate(this%storvegc_patch           (begp :endp))                   ;     this%storvegc_patch           (:)   = nan
       allocate(this%totvegc_patch            (begp :endp))                   ;     this%totvegc_patch            (:)   = nan
       allocate(this%totpftc_patch            (begp :endp))                   ;     this%totpftc_patch            (:)   = nan
       allocate(this%leafcmax_patch           (begp :endp))                   ;     this%leafcmax_patch           (:)   = nan
       allocate(this%grainc_patch             (begp :endp))                   ;     this%grainc_patch             (:)   = nan
       allocate(this%grainc_storage_patch     (begp :endp))                   ;     this%grainc_storage_patch     (:)   = nan
       allocate(this%grainc_xfer_patch        (begp :endp))                   ;     this%grainc_xfer_patch        (:)   = nan
       allocate(this%woodc_patch              (begp :endp))                   ;     this%woodc_patch              (:)   = nan     
       allocate(this%totvegc_abg_patch        (begp :endp))                   ;     this%totvegc_abg_patch            (:)   = nan

    endif
    allocate(this%cwdc_col                 (begc :endc))                   ;     this%cwdc_col                 (:)   = nan
    allocate(this%ctrunc_col               (begc :endc))                   ;     this%ctrunc_col               (:)   = nan
    allocate(this%ctrunc_vr_col            (begc :endc,1:nlevdecomp_full)) ;     this%ctrunc_vr_col            (:,:) = nan
    allocate(this%seedc_col                (begc :endc))                   ;     this%seedc_col                (:)   = nan
    allocate(this%prod10c_col              (begc :endc))                   ;     this%prod10c_col              (:)   = nan
    allocate(this%prod100c_col             (begc :endc))                   ;     this%prod100c_col             (:)   = nan
    allocate(this%prod1c_col               (begc :endc))                   ;     this%prod1c_col               (:)   = nan
    allocate(this%totprodc_col             (begc :endc))                   ;     this%totprodc_col             (:)   = nan
    allocate(this%dyn_cbal_adjustments_col (begc :endc))                   ;     this%dyn_cbal_adjustments_col (:)   = nan
    allocate(this%totlitc_col              (begc :endc))                   ;     this%totlitc_col              (:)   = nan
    allocate(this%totsomc_col              (begc :endc))                   ;     this%totsomc_col              (:)   = nan
    allocate(this%totlitc_1m_col           (begc :endc))                   ;     this%totlitc_1m_col           (:)   = nan
    allocate(this%totsomc_1m_col           (begc :endc))                   ;     this%totsomc_1m_col           (:)   = nan
    allocate(this%totecosysc_col           (begc :endc))                   ;     this%totecosysc_col           (:)   = nan
    allocate(this%totcolc_col              (begc :endc))                   ;     this%totcolc_col              (:)   = nan
    allocate(this%rootc_col                (begc :endc))                   ;     this%rootc_col                (:)   = nan
    allocate(this%totvegc_col              (begc :endc))                   ;     this%totvegc_col              (:)   = nan
    allocate(this%leafc_col                (begc :endc))                   ;     this%leafc_col                (:)   = nan
    allocate(this%deadstemc_col            (begc :endc))                   ;     this%deadstemc_col            (:)   = nan
    allocate(this%fuelc_col                (begc :endc))                   ;     this%fuelc_col                (:)   = nan
    allocate(this%fuelc_crop_col           (begc :endc))                   ;     this%fuelc_crop_col           (:)   = nan
    allocate(this%decomp_cpools_col        (begc :endc,1:ndecomp_pools))   ;     this%decomp_cpools_col        (:,:) = nan
    allocate(this%decomp_cpools_1m_col     (begc :endc,1:ndecomp_pools))   ;     this%decomp_cpools_1m_col     (:,:) = nan
    allocate(this%totpftc_col              (begc :endc))                   ;     this%totpftc_col              (:)   = nan
    allocate(this%totvegc_col              (begc :endc))                   ;     this%totvegc_col              (:)   = nan

    allocate(this%totvegc_abg_col          (begc :endc))                   ;     this%totvegc_abg_col              (:)   = nan


    allocate(this%totabgc_col              (begc :endc))                   ;     this%totabgc_col              (:)   = nan
    allocate(this%totblgc_col              (begc:endc))                    ;     this%totblgc_col              (:)   = nan
    allocate(this%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
    this%decomp_cpools_vr_col(:,:,:)= nan

    allocate(this%decomp_som2c_vr_col(begc:endc,1:nlevdecomp_full)); this%decomp_som2c_vr_col(:,:)= nan
    allocate(this%begcb_patch (begp:endp));     this%begcb_patch (:) = nan
    allocate(this%begcb_col   (begc:endc));     this%begcb_col   (:) = nan
    allocate(this%begcb_grc   (begg:endg));     this%begcb_grc   (:) = nan
    allocate(this%endcb_patch (begp:endp));     this%endcb_patch (:) = nan
    allocate(this%endcb_col   (begc:endc));     this%endcb_col   (:) = nan
    allocate(this%endcb_grc   (begg:endg));     this%endcb_grc   (:) = nan
    allocate(this%errcb_patch (begp:endp));     this%errcb_patch (:) = nan
    allocate(this%errcb_col   (begc:endc));     this%errcb_col   (:) = nan
    allocate(this%errcb_grc   (begg:endg));     this%errcb_grc   (:) = nan

    allocate(this%cropseedc_deficit_patch  (begp:endp)) ; this%cropseedc_deficit_patch  (:) = nan
    allocate(this%seedc_grc                (begg:endg)) ; this%seedc_grc                (:) = nan
    allocate(this%totpftc_beg_col(begc:endc));  this%totpftc_beg_col (:) = nan
    allocate(this%cwdc_beg_col   (begc:endc));  this%cwdc_beg_col    (:) = nan
    allocate(this%totlitc_beg_col(begc:endc));  this%totlitc_beg_col (:) = nan
    allocate(this%totsomc_beg_col(begc:endc));  this%totsomc_beg_col (:) = nan
    
    allocate(this%totpftc_end_col(begc:endc));  this%totpftc_end_col (:) = nan
    allocate(this%cwdc_end_col   (begc:endc));  this%cwdc_end_col    (:) = nan
    allocate(this%totlitc_end_col(begc:endc));  this%totlitc_end_col (:) = nan
    allocate(this%totsomc_end_col(begc:endc));  this%totsomc_end_col (:) = nan

  end subroutine InitAllocate


end module CNCarbonStateType
