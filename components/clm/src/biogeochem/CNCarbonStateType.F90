module CNCarbonStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, crop_prog, nlevdecomp
  use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi, zsoi
  use landunit_varcon        , only : istcrop 
  use clm_varctl             , only : iulog, use_vertsoilc, use_cndv, spinup_state 
  use decompMod              , only : bounds_type
  use CNStateType            , only : cnstate_type
  use pftvarcon              , only : npcropmin
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use subgridAveMod          , only : p2c
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use clm_varctl             , only : nu_com, use_fates, use_crop
  use VegetationType         , only : veg_pp
  use SpeciesMod           , only : species_from_string
  use dynPatchStateUpdaterMod, only : patch_state_updater_type

  ! bgc interface & pflotran
  use clm_varctl             , only : use_clm_interface, use_pflotran, pf_cmode
  
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
     procedure , public  :: SetValues 
     procedure , public  :: ZeroDWT
     procedure , public  :: Restart
     procedure , public  :: Summary
     procedure , public  :: DynamicPatchAdjustments
     procedure , public  :: DynamicColumnAdjustments
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type carbonstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type, ratio, c12_carbonstate_vars)

    class(carbonstate_type)                       :: this
    type(bounds_type)      , intent(in)           :: bounds  
    character(len=3)       , intent(in)           :: carbon_type
    real(r8)               , intent(in)           :: ratio
    type(carbonstate_type) , intent(in), optional :: c12_carbonstate_vars

    this%species = species_from_string(carbon_type)

    call this%InitAllocate ( bounds)
    call this%InitHistory ( bounds, carbon_type)
    if (present(c12_carbonstate_vars)) then
       call this%InitCold  ( bounds, ratio, c12_carbonstate_vars)
    else
       call this%InitCold  ( bounds, ratio)
    end if

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

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use clm_varpar , only : nlevdecomp, nlevdecomp_full, nlevgrnd
    use clm_varctl , only : use_c13, use_c14
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    character(len=3)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    ! Note (RGK 04-2017) - I am taking the ultra-conservative
    ! approach to identifying which carbon state variables should 
    ! be output when FATES is turned on.  Over time we should relax
    ! this restriction and identify which output state variables
    ! have all the required information.

    if ( use_fates ) then
       if (carbon_type == 'c12') then
          if ( nlevdecomp_full > 1 ) then
             this%totlitc_1m_col(begc:endc) = spval
             call hist_addfld1d (fname='TOTLITC_1m', units='gC/m^2', &
                  avgflag='A', long_name='total litter carbon to 1 meter depth', &
                  ptr_col=this%totlitc_1m_col)
             
             this%totsomc_1m_col(begc:endc) = spval
             call hist_addfld1d (fname='TOTSOMC_1m', units='gC/m^2', &
                  avgflag='A', long_name='total soil organic matter carbon to 1 meter depth', &
                  ptr_col=this%totsomc_1m_col)
          end if

          this%totlitc_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
               avgflag='A', long_name='total litter carbon', &
               ptr_col=this%totlitc_col)
          
          this%totsomc_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
               avgflag='A', long_name='total soil organic matter carbon', &
               ptr_col=this%totsomc_col)

       end if
       return
    end if
    !-------------------------------
    ! C12 state variables
    !-------------------------------

    if (carbon_type == 'c12') then

       if (crop_prog) then
          this%grainc_patch(begp:endp) = spval
          call hist_addfld1d (fname='GRAINC', units='gC/m^2', &
                avgflag='A', long_name='grain C', &
                ptr_patch=this%grainc_patch, default='inactive')

          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch)
       end if

       this%woodc_patch(begp:endp) = spval
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
             avgflag='A', long_name='wood C', &
             ptr_patch=this%woodc_patch)

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
             avgflag='A', long_name='leaf C', &
             ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='leaf C storage', &
             ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
             avgflag='A', long_name='leaf C transfer', &
             ptr_patch=this%leafc_xfer_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
             avgflag='A', long_name='fine root C', &
             ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='fine root C storage', &
             ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='fine root C transfer', &
             ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
             avgflag='A', long_name='live stem C', &
             ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='live stem C storage', &
             ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
             avgflag='A', long_name='live stem C transfer', &
             ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
             avgflag='A', long_name='dead stem C', &
             ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='dead stem C storage', &
             ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
             avgflag='A', long_name='dead stem C transfer', &
             ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C', &
             ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C storage', &
             ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C transfer', &
             ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C', &
             ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C storage', &
             ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C transfer', &
             ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='growth respiration storage', &
             ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
             avgflag='A', long_name='growth respiration transfer', &
             ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
             avgflag='A', long_name='temporary photosynthate C pool', &
             ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
             avgflag='A', long_name='temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool_patch, default='active')

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
             avgflag='A', long_name='patch-level sink for C truncation', &
             ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
             avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
             avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
             avgflag='A', long_name='total vegetation carbon, excluding cpool', &
             ptr_patch=this%totvegc_patch)

       this%totpftc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
             avgflag='A', long_name='total patch-level carbon, including cpool', &
             ptr_patch=this%totpftc_patch)

       this%totvegc_abg_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC_ABG', units='gC/m^2', &
            avgflag='A', long_name='total aboveground vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_abg_patch)

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='SEEDC_GRC', units='gC/m^2', &
            avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc)

    end if

    !-------------------------------
    ! C13 state variables 
    !-------------------------------

    if ( carbon_type == 'c13' ) then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C', &
             ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C storage', &
             ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C transfer', &
             ptr_patch=this%leafc_xfer_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C', &
             ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C storage', &
             ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C transfer', &
             ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C', &
             ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C storage', &
             ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C transfer', &
             ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C', &
             ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C storage', &
             ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C transfer', &
             ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C', &
             ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C storage', &
             ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C transfer', &
             ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C', &
             ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C storage', &
             ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C transfer', &
             ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 growth respiration storage', &
             ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 growth respiration transfer', &
             ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
             avgflag='A', long_name='C13 temporary photosynthate C pool', &
             ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
             avgflag='A', long_name='C13 temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool_patch)

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
             avgflag='A', long_name='C13 patch-level sink for C truncation', &
             ptr_patch=this%ctrunc_patch)

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
             ptr_patch=this%totvegc_patch)

       this%totpftc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total patch-level carbon, including cpool', &
             ptr_patch=this%totpftc_patch)

    endif

    !-------------------------------
    ! C14 state variables 
    !-------------------------------

    if ( carbon_type == 'c14') then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C', &
             ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C storage', &
             ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C transfer', &
             ptr_patch=this%leafc_xfer_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C', &
             ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C storage', &
             ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C transfer', &
             ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C', &
             ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C storage', &
             ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C transfer', &
             ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C', &
             ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C storage', &
             ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C transfer', &
             ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C', &
             ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C storage', &
             ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C transfer', &
             ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C', &
             ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C storage', &
             ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C transfer', &
             ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 growth respiration storage', &
             ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 growth respiration transfer', &
             ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL', units='gC14/m^2', &
             avgflag='A', long_name='C14 temporary photosynthate C pool', &
             ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_XSMRPOOL', units='gC14/m^2', &
             avgflag='A', long_name='C14 temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool_patch)

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PFT_CTRUNC', units='gC14/m^2', &
             avgflag='A', long_name='C14 patch-level sink for C truncation', &
             ptr_patch=this%ctrunc_patch)

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DISPVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_STORVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total vegetation carbon, excluding cpool', &
             ptr_patch=this%totvegc_patch)

       this%totpftc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTPFTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total patch-level carbon, including cpool', &
             ptr_patch=this%totpftc_patch)
    endif

    !-------------------------------
    ! C state variables - column
    !-------------------------------

    ! add history fields for all CLAMP CN variables



    if (carbon_type == 'c12') then


       !those variables are now ouput in betr
       this%decomp_cpools_col(begc:endc,:) = spval
       do l  = 1, ndecomp_pools
          if(trim(decomp_cascade_con%decomp_pool_name_history(l))=='')exit
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_cpools_vr_col(:,:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'

             call hist_addfld2d (fname=fieldname, units='gC/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                avgflag='A', long_name=longname, &
                ptr_col=data1dptr)

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_cpools_1m_col(:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr, default = 'inactive')
          endif
       end do

       if ( nlevdecomp_full > 1 ) then

          this%totlitc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTLITC_1m', units='gC/m^2', &
                avgflag='A', long_name='total litter carbon to 1 meter depth', &
                ptr_col=this%totlitc_1m_col)

          this%totsomc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTSOMC_1m', units='gC/m^2', &
                avgflag='A', long_name='total soil organic matter carbon to 1 meter depth', &
                ptr_col=this%totsomc_1m_col)
       end if


       this%ctrunc_col(begc:endc) = spval
       call hist_addfld1d (fname='COL_CTRUNC', units='gC/m^2',  &
             avgflag='A', long_name='column-level sink for C truncation', &
             ptr_col=this%ctrunc_col, default='inactive')

       this%seedc_col(begc:endc) = spval
       call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
             avgflag='A', long_name='pool for seeding new Patches', &
             ptr_col=this%seedc_col, default='inactive')

       this%totlitc_col(begc:endc) = spval
       call hist_addfld1d (fname='LITTERC', units='gC/m^2', &
             avgflag='A', long_name='litter C', &
             ptr_col=this%totlitc_col)
       call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
             avgflag='A', long_name='total litter carbon', &
             ptr_col=this%totlitc_col)

       this%totsomc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
             avgflag='A', long_name='total soil organic matter carbon', &
             ptr_col=this%totsomc_col)

       call hist_addfld1d (fname='SOILC', units='gC/m^2', &
             avgflag='A', long_name='soil C', &
             ptr_col=this%totsomc_col)

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
             avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool but excl product pools', &
             ptr_col=this%totecosysc_col)

       this%totcolc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
             avgflag='A', long_name='total column carbon, incl veg and cpool but excl product pools', &
             ptr_col=this%totcolc_col)

       this%prod10c_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD10C', units='gC/m^2', &
             avgflag='A', long_name='10-yr wood product C', &
             ptr_col=this%prod10c_col, default='inactive')

       this%prod100c_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD100C', units='gC/m^2', &
             avgflag='A', long_name='100-yr wood product C', &
             ptr_col=this%prod100c_col, default='inactive')

       this%prod1c_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD1C', units='gC/m^2', &
             avgflag='A', long_name='1-yr crop product C', &
             ptr_col=this%prod1c_col, default='inactive')

       this%totprodc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTPRODC', units='gC/m^2', &
             avgflag='A', long_name='total wood product C', &
             ptr_col=this%totprodc_col, default='inactive')

       this%fuelc_col(begc:endc) = spval
       call hist_addfld1d (fname='FUELC', units='gC/m^2', &
             avgflag='A', long_name='fuel load', &
             ptr_col=this%fuelc_col, default='inactive')

    end if

    !-------------------------------
    ! C13 state variables - column
    !-------------------------------

    if ( carbon_type == 'c13' ) then


       this%decomp_cpools_vr_col(begc:endc,:,:) = spval
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_cpools_vr_col(:,:,l)
             fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC13/m^2', &
                avgflag='A', long_name=longname, &
                ptr_col=data1dptr)
       end do

       this%seedc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
             avgflag='A', long_name='C13 pool for seeding new Patches', &
             ptr_col=this%seedc_col)

       this%ctrunc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_COL_CTRUNC', units='gC13/m^2',  &
             avgflag='A', long_name='C13 column-level sink for C truncation', &
             ptr_col=this%ctrunc_col)

       this%totlitc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTLITC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total litter carbon', &
             ptr_col=this%totlitc_col)

       this%totsomc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTSOMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total soil organic matter carbon', &
             ptr_col=this%totsomc_col)

       if ( nlevdecomp_full > 1 ) then
          this%totlitc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C13_TOTLITC_1m', units='gC13/m^2', &
                avgflag='A', long_name='C13 total litter carbon to 1 meter', &
                ptr_col=this%totlitc_1m_col)

          this%totsomc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C13_TOTSOMC_1m', units='gC13/m^2', &
                avgflag='A', long_name='C13 total soil organic matter carbon to 1 meter', &
                ptr_col=this%totsomc_1m_col)
       endif

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool but excl product pools', &
             ptr_col=this%totecosysc_col)

       this%totcolc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total column carbon, incl veg and cpool but excl product pools', &
             ptr_col=this%totcolc_col)

       this%prod10c_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD10C', units='gC13/m^2', &
             avgflag='A', long_name='C13 10-yr wood product C', &
             ptr_col=this%prod10c_col)

       this%prod100c_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD100C', units='gC13/m^2', &
             avgflag='A', long_name='C13 100-yr wood product C', &
             ptr_col=this%prod100c_col)

       this%prod1c_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD1C', units='gC13/m^2', &
             avgflag='A', long_name='C13 1-yr crop product C', &
             ptr_col=this%prod1c_col)

       this%totprodc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTPRODC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total wood product C', &
             ptr_col=this%totprodc_col)

       if (use_crop) then
          this%grainc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C13 grain C (does not equal yield)', &
               ptr_patch=this%grainc_patch)
          this%cropseedc_deficit_patch(begp:endp) = spval

          call hist_addfld1d (fname='C13_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C13 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch)
       end if
       
    endif

    !-------------------------------
    ! C14 state variables - column
    !-------------------------------

    if ( carbon_type == 'c14' ) then

       this%decomp_cpools_vr_col(begc:endc,:,:) = spval
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_cpools_vr_col(:,:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_cpools_col(:,l)
          fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC14/m^2', &
                avgflag='A', long_name=longname, ptr_col=data1dptr)

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_cpools_1m_col(:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
             longname =  'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                   avgflag='A', long_name=longname, ptr_col=data1dptr, default='inactive')
          endif
       end do

       this%seedc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_SEEDC', units='gC14/m^2', &
             avgflag='A', long_name='C14 pool for seeding new Patches', &
             ptr_col=this%seedc_col)

       this%ctrunc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_COL_CTRUNC', units='gC14/m^2', &
             avgflag='A', long_name='C14 column-level sink for C truncation', &
             ptr_col=this%ctrunc_col)

       this%totlitc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTLITC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total litter carbon', &
             ptr_col=this%totlitc_col)

       this%totsomc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTSOMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total soil organic matter carbon', &
             ptr_col=this%totsomc_col)

       if ( nlevdecomp_full > 1 ) then       
          this%totlitc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C14_TOTLITC_1m', units='gC14/m^2', &
                avgflag='A', long_name='C14 total litter carbon to 1 meter', &
                ptr_col=this%totlitc_1m_col)

          this%totsomc_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C14_TOTSOMC_1m', units='gC14/m^2', &
                avgflag='A', long_name='C14 total soil organic matter carbon to 1 meter', &
                ptr_col=this%totsomc_1m_col)
       endif

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTECOSYSC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total ecosystem carbon, incl veg but excl cpool but excl product pools', &
             ptr_col=this%totecosysc_col)

       this%totcolc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTCOLC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total column carbon, incl veg and cpool but excl product pools', &
             ptr_col=this%totcolc_col)

       this%prod10c_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD10C', units='gC14/m^2', &
             avgflag='A', long_name='C14 10-yr wood product C', &
             ptr_col=this%prod10c_col)

       this%prod100c_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD100C', units='gC14/m^2', &
             avgflag='A', long_name='C14 100-yr wood product C', &
             ptr_col=this%prod100c_col)

       this%prod1c_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD1C', units='gC14/m^2', &
             avgflag='A', long_name='C14 1-yr crop product C', &
             ptr_col=this%prod1c_col)

       this%totprodc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTPRODC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total wood product C', &
             ptr_col=this%totprodc_col)

       if (use_crop) then
          this%grainc_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C14 grain C (does not equal yield)', &
               ptr_patch=this%grainc_patch)
          this%cropseedc_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C14 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseedc_deficit_patch)
       end if

    endif

 end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, ratio, c12_carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon , only: istsoil
    use pftvarcon       , only: noveg, npcropmin
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this 
    type(bounds_type), intent(in) :: bounds  
    real(r8), intent(in) :: ratio
    type(carbonstate_type), optional, intent(in) :: c12_carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,g,j,k
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !-----------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do


    if ( .not. use_fates ) then
       !-----------------------------------------------
       ! initialize patch-level carbon state variables
       !-----------------------------------------------

       do p = bounds%begp,bounds%endp

          this%leafcmax_patch(p) = 0._r8

          l = veg_pp%landunit(p)
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

             if (veg_pp%itype(p) == noveg) then
                this%leafc_patch(p)         = 0._r8
                this%leafc_storage_patch(p) = 0._r8
             else
                if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                   this%leafc_patch(p)         = 1._r8 * ratio
                   this%leafc_storage_patch(p) = 0._r8
                else if (veg_pp%itype(p) >= npcropmin) then ! prognostic crop types
                   this%leafc_patch(p) = 0._r8
                   this%leafc_storage_patch(p) = 0._r8
                else
                   this%leafc_patch(p) = 0._r8
                   this%leafc_storage_patch(p) = 1._r8 * ratio
                end if
             end if
             this%leafc_xfer_patch(p) = 0._r8

             this%frootc_patch(p)            = 0._r8 
             this%frootc_storage_patch(p)    = 0._r8 
             this%frootc_xfer_patch(p)       = 0._r8 

             this%livestemc_patch(p)         = 0._r8 
             this%livestemc_storage_patch(p) = 0._r8 
             this%livestemc_xfer_patch(p)    = 0._r8 

             if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
                this%deadstemc_patch(p) = 0.1_r8 * ratio
             else
                this%deadstemc_patch(p) = 0._r8 
             end if
             this%deadstemc_storage_patch(p)  = 0._r8 
             this%deadstemc_xfer_patch(p)     = 0._r8

             if (nu_com .ne. 'RD') then
                ! ECA competition calculate root NP uptake as a function of fine root biomass
                ! better to initialize root CNP pools with a non-zero value
                if (veg_pp%itype(p) .ne. noveg) then
                   if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                      this%leafc_patch(p) = 20._r8 * ratio
                      this%leafc_storage_patch(p) = 0._r8
                      this%frootc_patch(p) = 20._r8 * ratio
                      this%frootc_storage_patch(p) = 0._r8
                   else
                      this%leafc_patch(p) = 0._r8 
                      this%leafc_storage_patch(p) = 20._r8 * ratio
                      this%frootc_patch(p) = 0._r8
                      this%frootc_storage_patch(p) = 20._r8 * ratio
                   end if
                end if
             end if

             this%livecrootc_patch(p)         = 0._r8 
             this%livecrootc_storage_patch(p) = 0._r8 
             this%livecrootc_xfer_patch(p)    = 0._r8 

             this%deadcrootc_patch(p)         = 0._r8 
             this%deadcrootc_storage_patch(p) = 0._r8 
             this%deadcrootc_xfer_patch(p)    = 0._r8 

             this%gresp_storage_patch(p)      = 0._r8 
             this%gresp_xfer_patch(p)         = 0._r8 

             this%cpool_patch(p)              = 0._r8 
             this%xsmrpool_patch(p)           = 0._r8 
             this%ctrunc_patch(p)             = 0._r8 
             this%dispvegc_patch(p)           = 0._r8 
             this%storvegc_patch(p)           = 0._r8 
             this%totpftc_patch(p)            = 0._r8 
             this%woodc_patch(p)              = 0._r8

             if ( crop_prog )then
                this%grainc_patch(p)            = 0._r8 
                this%grainc_storage_patch(p)    = 0._r8 
                this%grainc_xfer_patch(p)       = 0._r8 
                this%cropseedc_deficit_patch(p) = 0._r8
             end if

             ! calculate totvegc explicitly so that it is available for the isotope 
             ! code on the first time step.

             this%totvegc_patch(p) = &
                  this%leafc_patch(p)              + &
                  this%leafc_storage_patch(p)      + &
                  this%leafc_xfer_patch(p)         + &
                  this%frootc_patch(p)             + &
                  this%frootc_storage_patch(p)     + &
                  this%frootc_xfer_patch(p)        + &
                  this%livestemc_patch(p)          + &
                  this%livestemc_storage_patch(p)  + &
                  this%livestemc_xfer_patch(p)     + &
                  this%deadstemc_patch(p)          + &
                  this%deadstemc_storage_patch(p)  + &
                  this%deadstemc_xfer_patch(p)     + &
                  this%livecrootc_patch(p)         + &
                  this%livecrootc_storage_patch(p) + &
                  this%livecrootc_xfer_patch(p)    + &
                  this%deadcrootc_patch(p)         + &
                  this%deadcrootc_storage_patch(p) + &
                  this%deadcrootc_xfer_patch(p)    + &
                  this%gresp_storage_patch(p)      + &
                  this%gresp_xfer_patch(p)         + &
                  this%cpool_patch(p)

             if ( crop_prog )then
                this%totvegc_patch(p) =  this%totvegc_patch(p) + &
                     this%grainc_patch(p)                            + &
                     this%grainc_storage_patch(p)                    + &
                     this%grainc_xfer_patch(p)
             end if
          endif

       end do
    endif ! .not. use_fates
    
    ! initialize column-level variables
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          if (.not. present(c12_carbonstate_vars)) then !c12

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   if (zsoi(j) < 0.3 ) then  !! only initialize upper soil column
                      this%decomp_cpools_vr_col(c,j,k) = decomp_cascade_con%initial_stock(k)
                   else
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                   endif
                end do
                this%ctrunc_vr_col(c,j) = 0._r8
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                   end do
                   this%ctrunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%decomp_cpools_col(c,1:ndecomp_pools)    = decomp_cascade_con%initial_stock(1:ndecomp_pools)
             this%decomp_cpools_1m_col(c,1:ndecomp_pools) = decomp_cascade_con%initial_stock(1:ndecomp_pools)

          else

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   this%decomp_cpools_vr_col(c,j,k) = c12_carbonstate_vars%decomp_cpools_vr_col(c,j,k) * ratio
                end do
                this%ctrunc_vr_col(c,j) = c12_carbonstate_vars%ctrunc_vr_col(c,j) * ratio
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_cpools_vr_col(c,j,k) = 0._r8
                   end do
                   this%ctrunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%cwdc_col(c) = c12_carbonstate_vars%cwdc_col(c) * ratio
             do k = 1, ndecomp_pools
                this%decomp_cpools_col(c,k)    = c12_carbonstate_vars%decomp_cpools_col(c,k) * ratio
                this%decomp_cpools_1m_col(c,k) = c12_carbonstate_vars%decomp_cpools_1m_col(c,k) * ratio
             end do

          endif

          this%cwdc_col(c)       = 0._r8
          this%ctrunc_col(c)     = 0._r8
          this%totlitc_col(c)    = 0._r8
          this%totsomc_col(c)    = 0._r8
          this%totlitc_1m_col(c) = 0._r8
          this%totsomc_1m_col(c) = 0._r8
          this%totecosysc_col(c) = 0._r8
          this%totcolc_col(c)    = 0._r8

          ! dynamic landcover state variables
          this%seedc_col(c)      = 0._r8
          this%prod10c_col(c)    = 0._r8
          this%prod100c_col(c)   = 0._r8
          this%prod1c_col(c)     = 0._r8
          this%totprodc_col(c)   = 0._r8

       end if

    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics
    
    do fc = 1,num_special_col
       c = special_col(fc)

       this%seedc_col(c)      = 0._r8
       this%prod10c_col(c)    = 0._r8
       this%prod100c_col(c)   = 0._r8
       this%prod1c_col(c)     = 0._r8
       this%totprodc_col(c)   = 0._r8
    end do

    do g = bounds%begg, bounds%endg
       this%seedc_grc(g) = 0._r8
    end do

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag, carbon_type, c12_carbonstate_vars, cnstate_vars)
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use shr_const_mod    , only : SHR_CONST_PDB
    use clm_time_manager , only : is_restart, get_nstep
    use clm_varcon       , only : c13ratio, c14ratio
    use clm_varctl       , only : spinup_mortality_factor, spinup_state

    use restUtilMod
    use ncdio_pio
    use tracer_varcon  , only : is_active_betr_bgc
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    type(bounds_type)         , intent(in)           :: bounds 
    type(file_desc_t)         , intent(inout)        :: ncid   ! netcdf id
    character(len=*)          , intent(in)           :: flag   !'read' or 'write'
    character(len=3)          , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    type (carbonstate_type)   , intent(in), optional :: c12_carbonstate_vars 
    type (cnstate_type)       , intent(in)           :: cnstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: i,j,k,l,c
    real(r8) :: vwc,psi             ! for calculating soilpsi
    real(r8) :: c3_del13c           ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c           ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1               ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1               ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2               ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2               ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8) :: m, m_veg            ! multiplier for the exit_spinup code
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname   ! temporary
    logical  :: readvar
    integer  :: idata
    logical  :: exit_spinup  = .false.
    logical  :: enter_spinup = .false.
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer  :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer  :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_carbonstate_vars)) then
          call endrun(msg=' ERROR: for C14 must pass in c12_carbonstate_vars as argument' //&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    if ( .not. use_fates ) then
       
       !--------------------------------
       ! patch carbon state variables (c12)
       !--------------------------------
       
       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='leafc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='cpool', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
       end if

       if (carbon_type == 'c12' .and. use_cndv) then
          call restartvar(ncid=ncid, flag=flag, varname='leafcmax', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafcmax_patch)
       end if

       !--------------------------------
       ! C13 pft carbon state variables 
       !--------------------------------

       if ( carbon_type == 'c13')  then
          if (.not. present(c12_carbonstate_vars)) then
             call endrun(msg=' ERROR: for C13 must pass in c12_carbonstate_vars as argument' // errMsg(__FILE__, __LINE__))
          end if

          if ( .not. is_restart() .and. get_nstep() == 1 ) then
             c3_del13c = -28._r8
             c4_del13c = -13._r8
             c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
             c3_r2 = c3_r1/(1._r8 + c3_r1)
             c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
             c4_r2 = c4_r1/(1._r8 + c4_r1)

             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%grainc_patch(i)            = c12_carbonstate_vars%grainc_patch(i)         * c3_r2
                   this%grainc_storage_patch(i)    = c12_carbonstate_vars%grainc_storage_patch(i) * c3_r2
                   this%grainc_xfer_patch(i)       = c12_carbonstate_vars%grainc_xfer_patch(i)    * c3_r2
                   this%dispvegc_patch(i)          = c12_carbonstate_vars%dispvegc_patch(i)       * c3_r2
                   this%storvegc_patch(i)          = c12_carbonstate_vars%storvegc_patch(i)       * c3_r2
                   this%totvegc_patch(i)           = c12_carbonstate_vars%totvegc_patch(i)        * c3_r2
                   this%totpftc_patch(i)           = c12_carbonstate_vars%totpftc_patch(i)        * c3_r2
                   this%woodc_patch(i)             = c12_carbonstate_vars%woodc_patch(i)          * c3_r2
                else
                   this%grainc_patch(i)            = c12_carbonstate_vars%grainc_patch(i)         * c4_r2
                   this%grainc_storage_patch(i)    = c12_carbonstate_vars%grainc_storage_patch(i) * c4_r2
                   this%grainc_xfer_patch(i)       = c12_carbonstate_vars%grainc_xfer_patch(i)    * c4_r2
                   this%dispvegc_patch(i)          = c12_carbonstate_vars%dispvegc_patch(i)       * c4_r2
                   this%storvegc_patch(i)          = c12_carbonstate_vars%storvegc_patch(i)       * c4_r2
                   this%totvegc_patch(i)           = c12_carbonstate_vars%totvegc_patch(i)        * c4_r2
                   this%totpftc_patch(i)           = c12_carbonstate_vars%totpftc_patch(i)        * c4_r2
                   this%woodc_patch(i)             = c12_carbonstate_vars%woodc_patch(i)          * c4_r2
                end if
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_patch)
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leafc_patch(i) = c12_carbonstate_vars%leafc_patch(i) * c3_r2
                else
                   this%leafc_patch(i) = c12_carbonstate_vars%leafc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leafc_storage_patch(i) = c12_carbonstate_vars%leafc_storage_patch(i) * c3_r2
                else
                   this%leafc_storage_patch(i) = c12_carbonstate_vars%leafc_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leafc_xfer_patch(i) = c12_carbonstate_vars%leafc_xfer_patch(i) * c3_r2
                else
                   this%leafc_xfer_patch(i) = c12_carbonstate_vars%leafc_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%frootc_patch(i) = c12_carbonstate_vars%frootc_patch(i) * c3_r2
                else
                   this%frootc_patch(i) = c12_carbonstate_vars%frootc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%frootc_storage_patch(i) = c12_carbonstate_vars%frootc_storage_patch(i) * c3_r2
                else
                   this%frootc_storage_patch(i) = c12_carbonstate_vars%frootc_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%frootc_xfer_patch(i) = c12_carbonstate_vars%frootc_xfer_patch(i) * c3_r2
                else
                   this%frootc_xfer_patch(i) = c12_carbonstate_vars%frootc_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestemc_patch(i) = c12_carbonstate_vars%livestemc_patch(i) * c3_r2
                else
                   this%livestemc_patch(i) = c12_carbonstate_vars%livestemc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestemc_storage_patch(i) = c12_carbonstate_vars%livestemc_storage_patch(i) * c3_r2
                else
                   this%livestemc_storage_patch(i) = c12_carbonstate_vars%livestemc_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestemc_xfer_patch(i) = c12_carbonstate_vars%livestemc_xfer_patch(i) * c3_r2
                else
                   this%livestemc_xfer_patch(i) = c12_carbonstate_vars%livestemc_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstemc_patch(i) = c12_carbonstate_vars%deadstemc_patch(i) * c3_r2
                else
                   this%deadstemc_patch(i) = c12_carbonstate_vars%deadstemc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstemc_storage_patch(i) = c12_carbonstate_vars%deadstemc_storage_patch(i) * c3_r2
                else
                   this%deadstemc_storage_patch(i) = c12_carbonstate_vars%deadstemc_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstemc_xfer_patch(i) = c12_carbonstate_vars%deadstemc_xfer_patch(i) * c3_r2
                else
                   this%deadstemc_xfer_patch(i) = c12_carbonstate_vars%deadstemc_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecrootc_patch(i) = c12_carbonstate_vars%livecrootc_patch(i) * c3_r2
                else
                   this%livecrootc_patch(i) = c12_carbonstate_vars%livecrootc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecrootc_storage_patch(i) = c12_carbonstate_vars%livecrootc_storage_patch(i) * c3_r2
                else
                   this%livecrootc_storage_patch(i) = c12_carbonstate_vars%livecrootc_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecrootc_xfer_patch(i) = c12_carbonstate_vars%livecrootc_xfer_patch(i) * c3_r2
                else
                   this%livecrootc_xfer_patch(i) = c12_carbonstate_vars%livecrootc_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcrootc_patch(i) = c12_carbonstate_vars%deadcrootc_patch(i) * c3_r2
                else
                   this%deadcrootc_patch(i) = c12_carbonstate_vars%deadcrootc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcrootc_storage_patch(i) = c12_carbonstate_vars%deadcrootc_storage_patch(i) * c3_r2
                else
                   this%deadcrootc_storage_patch(i) = c12_carbonstate_vars%deadcrootc_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcrootc_xfer_patch(i) = c12_carbonstate_vars%deadcrootc_xfer_patch(i) * c3_r2
                else
                   this%deadcrootc_xfer_patch(i) = c12_carbonstate_vars%deadcrootc_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%gresp_storage_patch(i) = c12_carbonstate_vars%gresp_storage_patch(i) * c3_r2
                else
                   this%gresp_storage_patch(i) = c12_carbonstate_vars%gresp_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%gresp_xfer_patch(i) = c12_carbonstate_vars%gresp_xfer_patch(i) * c3_r2
                else
                   this%gresp_xfer_patch(i) = c12_carbonstate_vars%gresp_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='cpool_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%cpool with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%cpool_patch(i) = c12_carbonstate_vars%cpool_patch(i) * c3_r2
                else
                   this%cpool_patch(i) = c12_carbonstate_vars%cpool_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%xsmrpool with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%xsmrpool_patch(i) = c12_carbonstate_vars%xsmrpool_patch(i) * c3_r2
                else
                   this%xsmrpool_patch(i) = c12_carbonstate_vars%xsmrpool_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%ctrunc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%ctrunc_patch(i) = c12_carbonstate_vars%ctrunc_patch(i) * c3_r2
                else
                   this%ctrunc_patch(i) = c12_carbonstate_vars%ctrunc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing carbonstate_vars %totvegc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%totvegc_patch(i) = c12_carbonstate_vars%totvegc_patch(i) * c3_r2
                else
                   this%totvegc_patch(i) = c12_carbonstate_vars%totvegc_patch(i) * c4_r2
                endif
             end do
          end if
       endif

       !--------------------------------
       ! C14 pft carbon state variables 
       !--------------------------------

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leafc_patch(i) /= spval .and. &
                     .not. isnan(this%leafc_patch(i)) ) then
                   this%leafc_patch(i) = c12_carbonstate_vars%leafc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leafc_storage_patch(i) /= spval .and. &
                     .not. isnan(this%leafc_storage_patch(i)) ) then
                   this%leafc_storage_patch(i) = c12_carbonstate_vars%leafc_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leafc_xfer_patch(i) /= spval .and. .not. isnan(this%leafc_xfer_patch(i)) ) then
                   this%leafc_xfer_patch(i) = c12_carbonstate_vars%leafc_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%frootc_patch(i) /= spval .and. &
                     .not. isnan(this%frootc_patch(i)) ) then
                   this%frootc_patch(i) = c12_carbonstate_vars%frootc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%frootc_storage_patch(i) /= spval .and. &
                     .not. isnan(this%frootc_storage_patch(i)) ) then
                   this%frootc_storage_patch(i) = c12_carbonstate_vars%frootc_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%frootc_xfer_patch(i) /= spval .and. &
                     .not. isnan(this%frootc_xfer_patch(i)) ) then
                   this%frootc_xfer_patch(i) = c12_carbonstate_vars%frootc_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestemc_patch(i) /= spval .and. .not. isnan(this%livestemc_patch(i)) ) then
                   this%livestemc_patch(i) = c12_carbonstate_vars%livestemc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestemc_storage_patch(i) /= spval .and. .not. isnan(this%livestemc_storage_patch(i)) ) then
                   this%livestemc_storage_patch(i) = c12_carbonstate_vars%livestemc_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestemc_xfer_patch(i) /= spval .and. .not. isnan(this%livestemc_xfer_patch(i)) ) then
                   this%livestemc_xfer_patch(i) = c12_carbonstate_vars%livestemc_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstemc_patch(i) /= spval .and. .not. isnan(this%deadstemc_patch(i)) ) then
                   this%deadstemc_patch(i) = c12_carbonstate_vars%deadstemc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstemc_storage_patch(i) /= spval .and. .not. isnan(this%deadstemc_storage_patch(i)) ) then
                   this%deadstemc_storage_patch(i) = c12_carbonstate_vars%deadstemc_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstemc_xfer_patch(i) /= spval .and. .not. isnan(this%deadstemc_xfer_patch(i)) ) then
                   this%deadstemc_xfer_patch(i) = c12_carbonstate_vars%deadstemc_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecrootc_patch(i) /= spval .and. .not. isnan(this%livecrootc_patch(i)) ) then
                   this%livecrootc_patch(i) = c12_carbonstate_vars%livecrootc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecrootc_storage_patch(i) /= spval .and. .not. isnan(this%livecrootc_storage_patch(i)) ) then
                   this%livecrootc_storage_patch(i) = c12_carbonstate_vars%livecrootc_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecrootc_xfer_patch(i) /= spval .and. .not. isnan(this%livecrootc_xfer_patch(i)) ) then
                   this%livecrootc_xfer_patch(i) = c12_carbonstate_vars%livecrootc_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcrootc_patch(i) /= spval .and. .not. isnan(this%deadcrootc_patch(i)) ) then
                   this%deadcrootc_patch(i) = c12_carbonstate_vars%deadcrootc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcrootc_storage_patch(i) /= spval .and. .not. isnan(this%deadcrootc_storage_patch(i)) ) then
                   this%deadcrootc_storage_patch(i) = c12_carbonstate_vars%deadcrootc_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog) 'initializing this%deadcrootc_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcrootc_xfer_patch(i) /= spval .and. .not. isnan(this%deadcrootc_xfer_patch(i)) ) then
                   this%deadcrootc_xfer_patch(i) = c12_carbonstate_vars%deadcrootc_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%gresp_storage_patch(i) /= spval .and. .not. isnan(this%gresp_storage_patch(i)) ) then
                   this%gresp_storage_patch(i) = c12_carbonstate_vars%gresp_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%gresp_xfer_patch(i) /= spval .and. .not. isnan(this%gresp_xfer_patch(i)) ) then
                   this%gresp_xfer_patch(i) = c12_carbonstate_vars%gresp_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='cpool_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%cpool_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%cpool_patch(i) /= spval .and. .not. isnan(this%cpool_patch(i)) ) then
                   this%cpool_patch(i) = c12_carbonstate_vars%cpool_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%xsmrpool_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%xsmrpool_patch(i) /= spval .and. .not. isnan(this%xsmrpool_patch(i)) ) then
                   this%xsmrpool_patch(i) = c12_carbonstate_vars%xsmrpool_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%ctrunc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%ctrunc_patch(i) /= spval .and. .not. isnan(this%ctrunc_patch(i)) ) then
                   this%ctrunc_patch(i) = c12_carbonstate_vars%ctrunc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%totvegc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%totvegc_patch(i) /= spval .and. .not. isnan(this%totvegc_patch(i)) ) then
                   this%totvegc_patch(i) = c12_carbonstate_vars%totvegc_patch(i) * c14ratio
                endif
             end do
          end if
       endif

       !--------------------------------
       ! patch prognostic crop variables
       !--------------------------------

       if (crop_prog) then
          if (carbon_type == 'c12') then
             call restartvar(ncid=ncid, flag=flag,  varname='grainc', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grainc_patch)

             call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C storage', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grainc_storage_patch)

             call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C transfer', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer_patch)

             call restartvar(ncid=ncid, flag=flag, varname='cropseedc_deficit', xtype=ncd_double,  &
                  dim1name='pft', long_name='pool for seeding new crop growth', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%cropseedc_deficit_patch)
          end if

       end if

    endif  ! .not. use_fates
    
    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    if (carbon_type == 'c12') then
       do k = 1, ndecomp_pools
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
          if (use_vertsoilc) then
             ptr2d => this%decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_cpools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
                  errMsg(__FILE__, __LINE__))
          end if
       end do
    end if
    if(is_active_betr_bgc)then
      if (carbon_type == 'c12') then
        call restartvar(ncid=ncid, flag=flag, varname='totblgc', xtype=ncd_double,  &
           dim1name='column', long_name='', units='', &
           interpinic_flag='interp', readvar=readvar, data=this%totblgc_col)

        call restartvar(ncid=ncid, flag=flag, varname='cwdc', xtype=ncd_double,  &
           dim1name='column', long_name='', units='', &
           interpinic_flag='interp', readvar=readvar, data=this%cwdc_col)
      endif
    endif
    if (carbon_type == 'c12') then
       if (use_vertsoilc) then
          ptr2d => this%ctrunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc_vr', xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%ctrunc_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
          call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc', xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlitc_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcolc_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10c_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100c_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1c_col)
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='totsomc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totsomc_col) 
    end if

    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c13' ) then
       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
          if (use_vertsoilc) then
             ptr2d => this%decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_cpools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%decomp_cpools_vr_col with atmospheric c13 value for: '//varname
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (c12_carbonstate_vars%decomp_cpools_vr_col(i,j,k) /= spval .and. &
                        .not. isnan(this%decomp_cpools_vr_col(i,j,k)) ) then
                         this%decomp_cpools_vr_col(i,j,k) = c12_carbonstate_vars%decomp_cpools_vr_col(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%seedc_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%seedc_col(i)) ) then
             this%seedc_col(i) = c12_carbonstate_vars%seedc_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       if (use_vertsoilc) then
          ptr2d => this%ctrunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%ctrunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlitc_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%totlitc_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%totlitc_col(i) ) ) then
             this%totlitc_col(i) = c12_carbonstate_vars%totlitc_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcolc_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%totcolc_col(i) /= spval .and. &
               .not. isnan (c12_carbonstate_vars%totcolc_col(i) ) ) then
             this%totcolc_col(i) = c12_carbonstate_vars%totcolc_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10c_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod10c_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod10c_col(i) ) ) then
             this%prod10c_col(i) = c12_carbonstate_vars%prod10c_col(i) * c3_r2
          endif
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100c_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod100c_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod100c_col(i) ) ) then
             this%prod100c_col(i) = c12_carbonstate_vars%prod100c_col(i) * c3_r2
          endif
       end if
    endif

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1c_col)
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod1c_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod1c_col(i) ) ) then
             this%prod1c_col(i) = c12_carbonstate_vars%prod1c_col(i) * c3_r2
          endif
       end if
    end if

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c14' ) then
       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
          if (use_vertsoilc) then
             ptr2d => this%decomp_cpools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_cpools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%decomp_cpools_vr_col with atmospheric c14 value for: '//trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (c12_carbonstate_vars%decomp_cpools_vr_col(i,j,k) /= spval .and. &
                        .not. isnan(c12_carbonstate_vars%decomp_cpools_vr_col(i,j,k)) ) then
                         this%decomp_cpools_vr_col(i,j,k) = c12_carbonstate_vars%decomp_cpools_vr_col(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%seedc_col with atmospheric c14 value'
          do i = bounds%begc,bounds%endc
             if (c12_carbonstate_vars%seedc_col(i) /= spval .and. &
                  .not. isnan(c12_carbonstate_vars%seedc_col(i)) ) then
                this%seedc_col(i) = c12_carbonstate_vars%seedc_col(i) * c14ratio
             endif
          end do
       end if
    end if

    if ( carbon_type == 'c14' ) then
       if (use_vertsoilc) then 
          ptr2d => this%ctrunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%ctrunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlitc_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%totlitc_col with atmospheric c14 value'
          if (c12_carbonstate_vars%totlitc_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%totlitc_col(i)) ) then
             this%totlitc_col(i) = c12_carbonstate_vars%totlitc_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcolc_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%totcolc_col with atmospheric c14 value'
          if (c12_carbonstate_vars%totcolc_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%totcolc_col(i)) ) then
             this%totcolc_col(i) = c12_carbonstate_vars%totcolc_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10c_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod10c_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod10c_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod10c_col(i)) ) then
             this%prod10c_col(i) = c12_carbonstate_vars%prod10c_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100c_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod100c_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod100c_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod100c_col(i)) ) then
             this%prod100c_col(i) = c12_carbonstate_vars%prod100c_col(i) * c14ratio
          endif
       end if
    endif

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1c_col)
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod1c_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod1c_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod1c_col(i)) ) then
             this%prod1c_col(i) = c12_carbonstate_vars%prod1c_col(i) * c14ratio
          endif
       end if
    end if

    !--------------------------------
    ! Spinup state
    !--------------------------------

    if (carbon_type == 'c12'  .or. carbon_type == 'c14') then
        if (flag == 'write') then
           idata = spinup_state
        end if
        if (carbon_type == 'c12' .or. (carbon_type == 'c14' .and. flag == 'read')) then
        call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
             long_name='Spinup state of the model that wrote this restart file: ' &
             // ' 0 = normal model mode, 1 = AD spinup', units='', &
             interpinic_flag='copy', readvar=readvar,  data=idata)
        end if

        if (flag == 'read') then
           if (readvar) then
              restart_file_spinup_state = idata
           else
              ! assume, for sake of backwards compatibility, that if spinup_state is not in 
              ! the restart file then current model state is the same as prior model state
              restart_file_spinup_state = spinup_state
              if ( masterproc ) then
                 write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                      // ' on spinup state used to generate the restart file. '
                 write(iulog,*) '   Assuming the same as current setting: ', spinup_state
              end if
           end if
        end if

        ! now compare the model and restart file spinup states, and either take the 
        ! model into spinup mode or out of it if they are not identical
        ! taking model out of spinup mode requires multiplying each decomposing pool 
        ! by the associated AD factor.
        ! putting model into spinup mode requires dividing each decomposing pool 
        ! by the associated AD factor.
        ! only allow this to occur on first timestep of model run.
        
        if (flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
           if (spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
              if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools out of AD spinup mode'
              exit_spinup = .true.
           else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
              if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools into AD spinup mode'
              enter_spinup = .true.
           else
              call endrun(msg=' CNRest: error in entering/exiting spinup.  spinup_state ' &
                   // ' != restart_file_spinup_state, but do not know what to do'//&
                   errMsg(__FILE__, __LINE__))
           end if
           if (get_nstep() >= 2) then
              call endrun(msg=' CNRest: error in entering/exiting spinup - should occur only when nstep = 1'//&
                   errMsg(__FILE__, __LINE__))
           endif
           do k = 1, ndecomp_pools
              do c = bounds%begc, bounds%endc
                 do j = 1, nlevdecomp
		    if ( exit_spinup ) then
		       m = decomp_cascade_con%spinup_factor(k)
                       if (decomp_cascade_con%spinup_factor(k) > 1) m = m / cnstate_vars%scalaravg_col(c,j)
                    else if ( enter_spinup ) then 
		       m = 1. / decomp_cascade_con%spinup_factor(k)
		       if (decomp_cascade_con%spinup_factor(k) > 1) m = m * cnstate_vars%scalaravg_col(c,j)
		    end if
                    this%decomp_cpools_vr_col(c,j,k) = this%decomp_cpools_vr_col(c,j,k) * m
                 end do
              end do
           end do
           do i = bounds%begp, bounds%endp
              if (exit_spinup) then 
                 m_veg = spinup_mortality_factor
              else if (enter_spinup) then 
                 m_veg = 1._r8 / spinup_mortality_factor
              end if
              this%deadstemc_patch(i)  = this%deadstemc_patch(i) * m_veg
              this%deadcrootc_patch(i) = this%deadcrootc_patch(i) * m_veg
           end do
        end if
     end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state variables
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    if ( .not. use_fates ) then
       do fi = 1,num_patch
          i = filter_patch(fi)

          this%leafc_patch(i)              = value_patch
          this%leafc_storage_patch(i)      = value_patch
          this%leafc_xfer_patch(i)         = value_patch
          this%frootc_patch(i)             = value_patch
          this%frootc_storage_patch(i)     = value_patch
          this%frootc_xfer_patch(i)        = value_patch
          this%livestemc_patch(i)          = value_patch
          this%livestemc_storage_patch(i)  = value_patch
          this%livestemc_xfer_patch(i)     = value_patch
          this%deadstemc_patch(i)          = value_patch
          this%deadstemc_storage_patch(i)  = value_patch
          this%deadstemc_xfer_patch(i)     = value_patch
          this%livecrootc_patch(i)         = value_patch
          this%livecrootc_storage_patch(i) = value_patch
          this%livecrootc_xfer_patch(i)    = value_patch
          this%deadcrootc_patch(i)         = value_patch
          this%deadcrootc_storage_patch(i) = value_patch
          this%deadcrootc_xfer_patch(i)    = value_patch
          this%gresp_storage_patch(i)      = value_patch
          this%gresp_xfer_patch(i)         = value_patch
          this%cpool_patch(i)              = value_patch
          this%xsmrpool_patch(i)           = value_patch
          this%ctrunc_patch(i)             = value_patch
          this%dispvegc_patch(i)           = value_patch
          this%storvegc_patch(i)           = value_patch
          this%totvegc_patch(i)            = value_patch
          this%totpftc_patch(i)            = value_patch
          this%woodc_patch(i)              = value_patch
          this%totvegc_abg_patch(i)        = value_patch
       end do
       
       if ( crop_prog ) then
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%grainc_patch(i)            = value_patch
             this%grainc_storage_patch(i)    = value_patch
             this%grainc_xfer_patch(i)       = value_patch
             this%cropseedc_deficit_patch(i) = value_patch
          end do
       endif
    endif ! .not. use_fates
    
    do fi = 1,num_column
       i = filter_column(fi)
       this%cwdc_col(i)       = value_column
       this%ctrunc_col(i)     = value_column
       this%totlitc_col(i)    = value_column
       this%totsomc_col(i)    = value_column
       this%totecosysc_col(i) = value_column
       this%totcolc_col(i)    = value_column
       this%rootc_col(i)      = value_column
       this%totvegc_col(i)    = value_column
       this%leafc_col(i)      = value_column
       this%deadstemc_col(i)  = value_column
       this%fuelc_col(i)      = value_column
       this%fuelc_crop_col(i) = value_column
       this%totlitc_1m_col(i) = value_column
       this%totsomc_1m_col(i) = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%ctrunc_vr_col(i,j) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cpools_col(i,k) = value_column
          this%decomp_cpools_1m_col(i,k) = value_column
       end do
    end do

    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cpools_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegc_patch(p)   = 0._r8
       this%storvegc_patch(p)   = 0._r8
       this%totpftc_patch(p)    = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform patch and column-level carbon summary calculations
    !
    ! !USES:
    use clm_varctl       , only: iulog
    use clm_time_manager , only: get_step_size
    use clm_varcon       , only: secspday
    use clm_varpar       , only: nlevdecomp, ndecomp_pools, nlevdecomp_full
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    real(r8) :: maxdepth        ! depth to integrate soil variables
    integer  :: nlev
    real(r8) :: cropseedc_deficit_col(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------

    ! calculate patch -level summary of carbon state

    if (use_fates) return

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
       this%dispvegc_patch(p) =        &
            this%leafc_patch(p)      + &
            this%frootc_patch(p)     + &
            this%livestemc_patch(p)  + &
            this%deadstemc_patch(p)  + &
            this%livecrootc_patch(p) + &
            this%deadcrootc_patch(p)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
       this%storvegc_patch(p) =                &
            this%cpool_patch(p)              + &
            this%leafc_storage_patch(p)      + &
            this%frootc_storage_patch(p)     + &
            this%livestemc_storage_patch(p)  + &
            this%deadstemc_storage_patch(p)  + &
            this%livecrootc_storage_patch(p) + &
            this%deadcrootc_storage_patch(p) + &
            this%leafc_xfer_patch(p)         + &
            this%frootc_xfer_patch(p)        + &
            this%livestemc_xfer_patch(p)     + &
            this%deadstemc_xfer_patch(p)     + &
            this%livecrootc_xfer_patch(p)    + &
            this%deadcrootc_xfer_patch(p)    + &
            this%gresp_storage_patch(p)      + &
            this%gresp_xfer_patch(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%storvegc_patch(p) =            &
               this%storvegc_patch(p)       + &
               this%grainc_storage_patch(p) + &
               this%grainc_xfer_patch(p)

          this%dispvegc_patch(p) =            &
               this%dispvegc_patch(p)       + &
               this%grainc_patch(p)
       end if

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totvegc_patch(p) = &
            this%dispvegc_patch(p) + &
            this%storvegc_patch(p)

       ! total pft-level carbon, including xsmrpool, ctrunc
       this%totpftc_patch(p) = &
            this%totvegc_patch(p) + &
            this%xsmrpool_patch(p) + &
            this%ctrunc_patch(p)

       ! (WOODC) - wood C
       this%woodc_patch(p) = &
            this%deadstemc_patch(p)    + &
            this%livestemc_patch(p)    + &
            this%deadcrootc_patch(p)   + &
            this%livecrootc_patch(p)

       this%totvegc_abg_patch(p) = &
               this%leafc_patch(p)              + &
               this%leafc_storage_patch(p)      + &
               this%leafc_xfer_patch(p)         + &
               this%livestemc_patch(p)          + &
               this%livestemc_storage_patch(p)  + &
               this%livestemc_xfer_patch(p)     + &
               this%deadstemc_patch(p)          + &
               this%deadstemc_storage_patch(p)  + &
               this%deadstemc_xfer_patch(p)


    end do


    call p2c(bounds, num_soilc, filter_soilc, &
         this%totpftc_patch(bounds%begp:bounds%endp), &
         this%totpftc_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegc_patch(bounds%begp:bounds%endp), &
         this%totvegc_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegc_abg_patch(bounds%begp:bounds%endp), &
         this%totvegc_abg_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%cropseedc_deficit_patch(bounds%begp:bounds%endp), &
         cropseedc_deficit_col(bounds%begc:bounds%endc))

    ! column level summary

     nlev = nlevdecomp
     if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full


      ! vertically integrate each of the decomposing C pools
      do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_cpools_col(c,l) = 0._r8
       end do
      end do
      do l = 1, ndecomp_pools
       do j = 1, nlev
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cpools_col(c,l) = &
                  this%decomp_cpools_col(c,l) + &
                  this%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
      end do

      if ( nlevdecomp > 1) then

       ! vertically integrate each of the decomposing C pools to 1 meter
       maxdepth = 1._r8
       do l = 1, ndecomp_pools
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cpools_1m_col(c,l) = 0._r8
          end do
       end do
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp
             if ( zisoi(j) <= maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   this%decomp_cpools_1m_col(c,l) = &
                        this%decomp_cpools_1m_col(c,l) + &
                        this%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
                end do
             elseif ( zisoi(j-1) < maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   this%decomp_cpools_1m_col(c,l) = &
                        this%decomp_cpools_1m_col(c,l) + &
                        this%decomp_cpools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
                end do
             endif
          end do
       end do

       ! total litter carbon in the top meter (TOTLITC_1m)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%totlitc_1m_col(c)         = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%totlitc_1m_col(c) = &
                     this%totlitc_1m_col(c) + &
                     this%decomp_cpools_1m_col(c,l)
             end do
          endif
       end do

       ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%totsomc_1m_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_soil(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%totsomc_1m_col(c) = &
                     this%totsomc_1m_col(c) + &
                     this%decomp_cpools_1m_col(c,l)
             end do
          end if
       end do

      endif
    
      ! total litter carbon (TOTLITC)
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totlitc_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totlitc_col(c) = &
                  this%totlitc_col(c) + &
                  this%decomp_cpools_col(c,l)
          end do
       endif
      end do

      ! total soil organic matter carbon (TOTSOMC)
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totsomc_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_soil(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totsomc_col(c) = &
                  this%totsomc_col(c) + &
                  this%decomp_cpools_col(c,l)
          end do
       end if
      end do


      ! coarse woody debris carbon
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%cwdc_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_cwd(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%cwdc_col(c) = &
                  this%cwdc_col(c) + &
                  this%decomp_cpools_col(c,l)
          end do
       end if
      end do

    ! truncation carbon
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%ctrunc_col(c) = 0._r8
    end do
    do j = 1, nlev
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%ctrunc_col(c) = &
               this%ctrunc_col(c) + &
               this%ctrunc_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! total product carbon
       this%totprodc_col(c) = &
            this%prod10c_col(c)  + &
            this%prod100c_col(c) + &
            this%prod1c_col(c) 

       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosysc_col(c) = &
            this%cwdc_col(c)     + &
            this%totlitc_col(c)  + &
            this%totsomc_col(c)  + &
            this%totprodc_col(c) + &
            this%totvegc_col(c)

       ! total column carbon, including veg and cpool (TOTCOLC)
       ! adding col_ctrunc, seedc
       this%totcolc_col(c) = &
            this%totpftc_col(c)  + &
            this%cwdc_col(c)     + &
            this%totlitc_col(c)  + &
            this%totsomc_col(c)  + &
            this%prod1c_col(c)   + &
            this%ctrunc_col(c)   + &
            cropseedc_deficit_col(c)
            
       this%totabgc_col(c) = &
            this%totpftc_col(c)  + &
            this%totprodc_col(c) + &
            this%seedc_col(c)    + &
            this%ctrunc_col(c)    
    end do

  end subroutine Summary

  !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments( this, &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafc_seed,                      &
       dwt_deadstemc_seed,                  &
       conv_cflux,                          &
       dwt_frootc_to_litter,                &
       dwt_livecrootc_to_litter,            &
       dwt_deadcrootc_to_litter,            &
       prod10_cflux,                        &
       prod100_cflux,                       &
       crop_product_cflux                   &
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
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    class(carbonstate_type)        , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:) ! patch filter that includes inactive points
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafc_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemc_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_cflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootc_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_cflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_cflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_cflux       (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical                     :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafc_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_storage_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_xfer_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemc_patch(bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_cflux(bounds%begp:bounds%endp)
    real(r8)                    :: deadstemc_patch_temp(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'CStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafc_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemc_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_cflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootc_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_cflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_cflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_cflux       ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = this%species                          , &
         leaf_patch                 = this%leafc_patch(begp:endp)           , &
         leaf_storage_patch         = this%leafc_storage_patch(begp:endp)   , &
         leaf_xfer_patch            = this%leafc_xfer_patch(begp:endp)      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafc_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafc_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafc_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemc_patch(begp:endp))

    ! 1) LEAFC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%leafc_patch   (begp:endp)           , &
         flux_out_grc_area = conv_cflux       (begp:endp)           , &
         seed              = seed_leafc_patch (begp:endp)           , &
         seed_addition     = dwt_leafc_seed   (begp:endp))


    ! 2) LEAFC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%leafc_storage_patch   (begp:endp)   , &
         flux_out_grc_area = conv_cflux               (begp:endp)   , &
         seed              = seed_leafc_storage_patch (begp:endp)   , &
         seed_addition     = dwt_leafc_seed           (begp:endp))

    ! 3) LEAF_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%leafc_xfer_patch   (begp:endp)      , &
         flux_out_grc_area = conv_cflux            (begp:endp)      , &
         seed              = seed_leafc_xfer_patch (begp:endp)      , &
         seed_addition     = dwt_leafc_seed        (begp:endp))

    ! 4) FROOTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%frootc_patch(begp:endp)             , &
         flux_out_col_area = dwt_frootc_to_litter(begp:endp))
    
    ! 5) FROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%frootc_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 6) FROOTC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%frootc_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 7) LIVESTEMC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestemc_patch(begp:endp)          , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 8) LIVESTEMC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestemc_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 9) LIVESTEMC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestemc_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 10) PROD10_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = this%deadstemc_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod10                          , &
         var                        = deadstemc_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp) )

    ! 11) PROD100_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = this%deadstemc_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod100                         , &
         var                        = deadstemc_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp))

    ! 12) DEADSTEMC_PATCH
    wood_product_cflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pconv                          , &
         var                        = this%deadstemc_patch(begp:endp) , &
         flux1_out                  = conv_cflux              (begp:endp) , &
         flux2_out                  = wood_product_cflux       (begp:endp) , &
         seed                       = seed_deadstemc_patch    (begp:endp) , &
         seed_addition              = dwt_deadstemc_seed     (begp:endp))

    ! 13) DEADSTEMC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstemc_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 14) DEADSTEMC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstemc_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 15) LIVECROOTC_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecrootc_patch(begp:endp)         , &
         flux_out_col_area = dwt_livecrootc_to_litter(begp:endp))

    ! 16) LIVECROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecrootc_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 17) LIVECROOTC_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecrootc_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 18) DEADCROOTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcrootc_patch(begp:endp)         , &
         flux_out_col_area = dwt_deadcrootc_to_litter(begp:endp))

    ! 19) DEADCROOTC_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcrootc_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcrootc_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 21) GRESP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%gresp_storage_patch(begp:endp)      , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 22) GRESP_XFER_STORAGE
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%gresp_xfer_patch(begp:endp)         , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 23) CPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%cpool_patch(begp:endp)              , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 24) XSMRPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%xsmrpool_patch(begp:endp)           , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 25) CTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%ctrunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 26) DISPVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%dispvegc_patch(begp:endp))

    ! 27) STORVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%storvegc_patch(begp:endp))

    ! 28) TOTVEGC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totvegc_patch(begp:endp))

    ! 29) TOTPFTC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totpftc_patch(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call patch_state_updater%update_patch_state(         &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = this%cropseedc_deficit_patch(begp:endp) , &
            flux_out_grc_area = conv_cflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp, endp
       dwt_frootc_to_litter(p)     = -1._r8 * dwt_frootc_to_litter(p)
       dwt_livecrootc_to_litter(p) = -1._r8 * dwt_livecrootc_to_litter(p)
       dwt_deadcrootc_to_litter(p) = -1._r8 * dwt_deadcrootc_to_litter(p)
    end do

  end subroutine DynamicPatchAdjustments
  
  !-----------------------------------------------------------------------
  subroutine DynamicColumnAdjustments( this, &
       bounds, clump_index, column_state_updater)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use dynColumnStateUpdaterMod, only : column_state_updater_type
    !
    ! !ARGUMENTS:
    class(carbonstate_type)         , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'CStateDynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    this%dyn_cbal_adjustments_col(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call column_state_updater%update_column_state_no_special_handling( &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = this%decomp_cpools_vr_col(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))

          this%dyn_cbal_adjustments_col(begc:endc) = &
               this%dyn_cbal_adjustments_col(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%ctrunc_vr_col(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_cbal_adjustments_col(begc:endc) = &
            this%dyn_cbal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do

  end subroutine DynamicColumnAdjustments

end module CNCarbonStateType
