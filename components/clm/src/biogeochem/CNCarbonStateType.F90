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
  use CNSpeciesMod           , only : species_from_string
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use NutrientStateType      , only : nutrientstate_type

  ! bgc interface & pflotran
  use clm_varctl             , only : use_clm_interface, use_pflotran, pf_cmode
  
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public, extends(nutrientstate_type) :: carbonstate_type
     
     real(r8), pointer :: gresp_storage_patch      (:)     ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch         (:)     ! (gC/m2) growth respiration transfer
     real(r8), pointer :: xsmrpool_patch           (:)     ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: woodc_patch              (:)     ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch           (:)     ! (gC/m2) ann max leaf C
     real(r8), pointer :: rootc_col                (:)     ! (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: leafc_col                (:)     ! (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col            (:)     ! (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                (:)     ! fuel avalability factor for Reg.C (0-1)
     real(r8), pointer :: fuelc_crop_col           (:)     ! fuel avalability factor for Reg.A (0-1)
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
       allocate(this%leaf_patch              (begp :endp))                   ;     this%leaf_patch              (:)   = nan
       allocate(this%leaf_storage_patch      (begp :endp))                   ;     this%leaf_storage_patch      (:)   = nan
       allocate(this%leaf_xfer_patch         (begp :endp))                   ;     this%leaf_xfer_patch         (:)   = nan
       allocate(this%froot_patch             (begp :endp))                   ;     this%froot_patch             (:)   = nan
       allocate(this%froot_storage_patch     (begp :endp))                   ;     this%froot_storage_patch     (:)   = nan
       allocate(this%froot_xfer_patch        (begp :endp))                   ;     this%froot_xfer_patch        (:)   = nan
       allocate(this%livestem_patch          (begp :endp))                   ;     this%livestem_patch          (:)   = nan
       allocate(this%livestem_storage_patch  (begp :endp))                   ;     this%livestem_storage_patch  (:)   = nan
       allocate(this%livestem_xfer_patch     (begp :endp))                   ;     this%livestem_xfer_patch     (:)   = nan
       allocate(this%deadstem_patch          (begp :endp))                   ;     this%deadstem_patch          (:)   = nan
       allocate(this%deadstem_storage_patch  (begp :endp))                   ;     this%deadstem_storage_patch  (:)   = nan
       allocate(this%deadstem_xfer_patch     (begp :endp))                   ;     this%deadstem_xfer_patch     (:)   = nan
       allocate(this%livecroot_patch         (begp :endp))                   ;     this%livecroot_patch         (:)   = nan
       allocate(this%livecroot_storage_patch (begp :endp))                   ;     this%livecroot_storage_patch (:)   = nan
       allocate(this%livecroot_xfer_patch    (begp :endp))                   ;     this%livecroot_xfer_patch    (:)   = nan
       allocate(this%deadcroot_patch         (begp :endp))                   ;     this%deadcroot_patch         (:)   = nan
       allocate(this%deadcroot_storage_patch (begp :endp))                   ;     this%deadcroot_storage_patch (:)   = nan
       allocate(this%deadcroot_xfer_patch    (begp :endp))                   ;     this%deadcroot_xfer_patch    (:)   = nan
       allocate(this%gresp_storage_patch      (begp :endp))                   ;     this%gresp_storage_patch      (:)   = nan
       allocate(this%gresp_xfer_patch         (begp :endp))                   ;     this%gresp_xfer_patch         (:)   = nan
       allocate(this%pool_patch              (begp :endp))                   ;     this%pool_patch              (:)   = nan
       allocate(this%xsmrpool_patch           (begp :endp))                   ;     this%xsmrpool_patch           (:)   = nan
       allocate(this%veg_trunc_patch             (begp :endp))                   ;     this%veg_trunc_patch             (:)   = nan
       allocate(this%dispveg_patch           (begp :endp))                   ;     this%dispveg_patch           (:)   = nan
       allocate(this%storveg_patch           (begp :endp))                   ;     this%storveg_patch           (:)   = nan
       allocate(this%totveg_patch            (begp :endp))                   ;     this%totveg_patch            (:)   = nan
       allocate(this%totpft_patch            (begp :endp))                   ;     this%totpft_patch            (:)   = nan
       allocate(this%leafcmax_patch           (begp :endp))                   ;     this%leafcmax_patch           (:)   = nan
       allocate(this%grain_patch             (begp :endp))                   ;     this%grain_patch             (:)   = nan
       allocate(this%grain_storage_patch     (begp :endp))                   ;     this%grain_storage_patch     (:)   = nan
       allocate(this%grain_xfer_patch        (begp :endp))                   ;     this%grain_xfer_patch        (:)   = nan
       allocate(this%woodc_patch              (begp :endp))                   ;     this%woodc_patch              (:)   = nan     
       allocate(this%totveg_abg_patch        (begp :endp))                   ;     this%totveg_abg_patch            (:)   = nan

    endif
    allocate(this%cwd_col                 (begc :endc))                   ;     this%cwd_col                 (:)   = nan
    allocate(this%veg_trunc_col               (begc :endc))                   ;     this%veg_trunc_col               (:)   = nan
    allocate(this%soil_trunc_vr_col            (begc :endc,1:nlevdecomp_full)) ;     this%soil_trunc_vr_col            (:,:) = nan
    allocate(this%seed_col                (begc :endc))                   ;     this%seed_col                (:)   = nan
    allocate(this%prod10_col              (begc :endc))                   ;     this%prod10_col              (:)   = nan
    allocate(this%prod100_col             (begc :endc))                   ;     this%prod100_col             (:)   = nan
    allocate(this%prod1_col               (begc :endc))                   ;     this%prod1_col               (:)   = nan
    allocate(this%totprod_col             (begc :endc))                   ;     this%totprod_col             (:)   = nan
    allocate(this%dyn_bal_adjustments_col (begc :endc))                   ;     this%dyn_bal_adjustments_col (:)   = nan
    allocate(this%totlit_col              (begc :endc))                   ;     this%totlit_col              (:)   = nan
    allocate(this%totsom_col              (begc :endc))                   ;     this%totsom_col              (:)   = nan
    allocate(this%totlit_1m_col           (begc :endc))                   ;     this%totlit_1m_col           (:)   = nan
    allocate(this%totsom_1m_col           (begc :endc))                   ;     this%totsom_1m_col           (:)   = nan
    allocate(this%totecosys_col           (begc :endc))                   ;     this%totecosys_col           (:)   = nan
    allocate(this%totcol_col              (begc :endc))                   ;     this%totcol_col              (:)   = nan
    allocate(this%rootc_col                (begc :endc))                   ;     this%rootc_col                (:)   = nan
    allocate(this%totveg_col              (begc :endc))                   ;     this%totveg_col              (:)   = nan
    allocate(this%leafc_col                (begc :endc))                   ;     this%leafc_col                (:)   = nan
    allocate(this%deadstemc_col            (begc :endc))                   ;     this%deadstemc_col            (:)   = nan
    allocate(this%fuelc_col                (begc :endc))                   ;     this%fuelc_col                (:)   = nan
    allocate(this%fuelc_crop_col           (begc :endc))                   ;     this%fuelc_crop_col           (:)   = nan
    allocate(this%decomp_pools_col        (begc :endc,1:ndecomp_pools))   ;     this%decomp_pools_col        (:,:) = nan
    allocate(this%decomp_pools_1m_col     (begc :endc,1:ndecomp_pools))   ;     this%decomp_pools_1m_col     (:,:) = nan
    allocate(this%totpft_col              (begc :endc))                   ;     this%totpft_col              (:)   = nan
    allocate(this%totveg_col              (begc :endc))                   ;     this%totveg_col              (:)   = nan

    allocate(this%totveg_abg_col          (begc :endc))                   ;     this%totveg_abg_col              (:)   = nan


    allocate(this%totabg_col              (begc :endc))                   ;     this%totabg_col              (:)   = nan
    allocate(this%totblg_col              (begc:endc))                    ;     this%totblg_col              (:)   = nan
    allocate(this%decomp_pools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
    this%decomp_pools_vr_col(:,:,:)= nan

    allocate(this%decomp_som2c_vr_col(begc:endc,1:nlevdecomp_full)); this%decomp_som2c_vr_col(:,:)= nan
    allocate(this%beg_bal_col   (begc:endc));     this%beg_bal_col   (:) = nan
    allocate(this%beg_bal_grc   (begg:endg));     this%beg_bal_grc   (:) = nan
    allocate(this%end_bal_patch (begp:endp));     this%end_bal_patch (:) = nan
    allocate(this%end_bal_col   (begc:endc));     this%end_bal_col   (:) = nan
    allocate(this%end_bal_grc   (begg:endg));     this%end_bal_grc   (:) = nan
    allocate(this%err_bal_patch (begp:endp));     this%err_bal_patch (:) = nan
    allocate(this%err_bal_col   (begc:endc));     this%err_bal_col   (:) = nan
    allocate(this%err_bal_grc   (begg:endg));     this%err_bal_grc   (:) = nan

    allocate(this%cropseed_deficit_patch  (begp:endp)) ; this%cropseed_deficit_patch  (:) = nan
    allocate(this%seed_grc                (begg:endg)) ; this%seed_grc                (:) = nan

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
             this%totlit_1m_col(begc:endc) = spval
             call hist_addfld1d (fname='TOTLITC_1m', units='gC/m^2', &
                  avgflag='A', long_name='total litter carbon to 1 meter depth', &
                  ptr_col=this%totlit_1m_col)
             
             this%totsom_1m_col(begc:endc) = spval
             call hist_addfld1d (fname='TOTSOMC_1m', units='gC/m^2', &
                  avgflag='A', long_name='total soil organic matter carbon to 1 meter depth', &
                  ptr_col=this%totsom_1m_col)
          end if

          this%totlit_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
               avgflag='A', long_name='total litter carbon', &
               ptr_col=this%totlit_col)
          
          this%totsom_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
               avgflag='A', long_name='total soil organic matter carbon', &
               ptr_col=this%totsom_col)

       end if
       return
    end if
    !-------------------------------
    ! C12 state variables
    !-------------------------------

    if (carbon_type == 'c12') then

       if (crop_prog) then
          this%grain_patch(begp:endp) = spval
          call hist_addfld1d (fname='GRAINC', units='gC/m^2', &
                avgflag='A', long_name='grain C', &
                ptr_patch=this%grain_patch, default='inactive')

          this%cropseed_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseed_deficit_patch)
       end if

       this%woodc_patch(begp:endp) = spval
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
             avgflag='A', long_name='wood C', &
             ptr_patch=this%woodc_patch)

       this%leaf_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
             avgflag='A', long_name='leaf C', &
             ptr_patch=this%leaf_patch)

       this%leaf_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='leaf C storage', &
             ptr_patch=this%leaf_storage_patch, default='inactive')

       this%leaf_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
             avgflag='A', long_name='leaf C transfer', &
             ptr_patch=this%leaf_xfer_patch, default='inactive')

       this%froot_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
             avgflag='A', long_name='fine root C', &
             ptr_patch=this%froot_patch)

       this%froot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='fine root C storage', &
             ptr_patch=this%froot_storage_patch, default='inactive')

       this%froot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='fine root C transfer', &
             ptr_patch=this%froot_xfer_patch, default='inactive')

       this%livestem_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
             avgflag='A', long_name='live stem C', &
             ptr_patch=this%livestem_patch)

       this%livestem_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='live stem C storage', &
             ptr_patch=this%livestem_storage_patch, default='inactive')

       this%livestem_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
             avgflag='A', long_name='live stem C transfer', &
             ptr_patch=this%livestem_xfer_patch, default='inactive')

       this%deadstem_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
             avgflag='A', long_name='dead stem C', &
             ptr_patch=this%deadstem_patch)

       this%deadstem_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='dead stem C storage', &
             ptr_patch=this%deadstem_storage_patch, default='inactive')

       this%deadstem_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
             avgflag='A', long_name='dead stem C transfer', &
             ptr_patch=this%deadstem_xfer_patch, default='inactive')

       this%livecroot_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C', &
             ptr_patch=this%livecroot_patch)

       this%livecroot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C storage', &
             ptr_patch=this%livecroot_storage_patch, default='inactive')

       this%livecroot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='live coarse root C transfer', &
             ptr_patch=this%livecroot_xfer_patch, default='inactive')

       this%deadcroot_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C', &
             ptr_patch=this%deadcroot_patch)

       this%deadcroot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C storage', &
             ptr_patch=this%deadcroot_storage_patch,  default='inactive')

       this%deadcroot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
             avgflag='A', long_name='dead coarse root C transfer', &
             ptr_patch=this%deadcroot_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
             avgflag='A', long_name='growth respiration storage', &
             ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
             avgflag='A', long_name='growth respiration transfer', &
             ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%pool_patch(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
             avgflag='A', long_name='temporary photosynthate C pool', &
             ptr_patch=this%pool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
             avgflag='A', long_name='temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool_patch, default='active')

       this%veg_trunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
             avgflag='A', long_name='patch-level sink for C truncation', &
             ptr_patch=this%veg_trunc_patch, default='inactive')

       this%dispveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
             avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispveg_patch)

       this%storveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
             avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storveg_patch)

       this%totveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
             avgflag='A', long_name='total vegetation carbon, excluding cpool', &
             ptr_patch=this%totveg_patch)

       this%totpft_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
             avgflag='A', long_name='total patch-level carbon, including cpool', &
             ptr_patch=this%totpft_patch)

       this%totveg_abg_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC_ABG', units='gC/m^2', &
            avgflag='A', long_name='total aboveground vegetation carbon, excluding cpool', &
            ptr_patch=this%totveg_abg_patch)

       this%seed_grc(begg:endg) = spval
       call hist_addfld1d (fname='seed_grc', units='gC/m^2', &
            avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seed_grc)

    end if

    !-------------------------------
    ! C13 state variables 
    !-------------------------------

    if ( carbon_type == 'c13' ) then

       this%leaf_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C', &
             ptr_patch=this%leaf_patch)

       this%leaf_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C storage', &
             ptr_patch=this%leaf_storage_patch, default='inactive')

       this%leaf_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 leaf C transfer', &
             ptr_patch=this%leaf_xfer_patch, default='inactive')

       this%froot_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C', &
             ptr_patch=this%froot_patch)

       this%froot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C storage', &
             ptr_patch=this%froot_storage_patch, default='inactive')

       this%froot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 fine root C transfer', &
             ptr_patch=this%froot_xfer_patch, default='inactive')

       this%livestem_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C', &
             ptr_patch=this%livestem_patch)

       this%livestem_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C storage', &
             ptr_patch=this%livestem_storage_patch, default='inactive')

       this%livestem_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 live stem C transfer', &
             ptr_patch=this%livestem_xfer_patch, default='inactive')

       this%deadstem_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C', &
             ptr_patch=this%deadstem_patch)

       this%deadstem_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C storage', &
             ptr_patch=this%deadstem_storage_patch, default='inactive')

       this%deadstem_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead stem C transfer', &
             ptr_patch=this%deadstem_xfer_patch, default='inactive')

       this%livecroot_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C', &
             ptr_patch=this%livecroot_patch)

       this%livecroot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C storage', &
             ptr_patch=this%livecroot_storage_patch, default='inactive')

       this%livecroot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 live coarse root C transfer', &
             ptr_patch=this%livecroot_xfer_patch, default='inactive')

       this%deadcroot_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C', &
             ptr_patch=this%deadcroot_patch)

       this%deadcroot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C storage', &
             ptr_patch=this%deadcroot_storage_patch,  default='inactive')

       this%deadcroot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 dead coarse root C transfer', &
             ptr_patch=this%deadcroot_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
             avgflag='A', long_name='C13 growth respiration storage', &
             ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
             avgflag='A', long_name='C13 growth respiration transfer', &
             ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%pool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
             avgflag='A', long_name='C13 temporary photosynthate C pool', &
             ptr_patch=this%pool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
             avgflag='A', long_name='C13 temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool_patch)

       this%veg_trunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
             avgflag='A', long_name='C13 patch-level sink for C truncation', &
             ptr_patch=this%veg_trunc_patch)

       this%dispveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispveg_patch)

       this%storveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storveg_patch)

       this%totveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
             ptr_patch=this%totveg_patch)

       this%totpft_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total patch-level carbon, including cpool', &
             ptr_patch=this%totpft_patch)

    endif

    !-------------------------------
    ! C14 state variables 
    !-------------------------------

    if ( carbon_type == 'c14') then

       this%leaf_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C', &
             ptr_patch=this%leaf_patch)

       this%leaf_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C storage', &
             ptr_patch=this%leaf_storage_patch, default='inactive')

       this%leaf_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 leaf C transfer', &
             ptr_patch=this%leaf_xfer_patch, default='inactive')

       this%froot_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C', &
             ptr_patch=this%froot_patch)

       this%froot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C storage', &
             ptr_patch=this%froot_storage_patch, default='inactive')

       this%froot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 fine root C transfer', &
             ptr_patch=this%froot_xfer_patch, default='inactive')

       this%livestem_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C', &
             ptr_patch=this%livestem_patch)

       this%livestem_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C storage', &
             ptr_patch=this%livestem_storage_patch, default='inactive')

       this%livestem_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 live stem C transfer', &
             ptr_patch=this%livestem_xfer_patch, default='inactive')

       this%deadstem_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C', &
             ptr_patch=this%deadstem_patch)

       this%deadstem_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C storage', &
             ptr_patch=this%deadstem_storage_patch, default='inactive')

       this%deadstem_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead stem C transfer', &
             ptr_patch=this%deadstem_xfer_patch, default='inactive')

       this%livecroot_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C', &
             ptr_patch=this%livecroot_patch)

       this%livecroot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C storage', &
             ptr_patch=this%livecroot_storage_patch, default='inactive')

       this%livecroot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 live coarse root C transfer', &
             ptr_patch=this%livecroot_xfer_patch, default='inactive')

       this%deadcroot_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C', &
             ptr_patch=this%deadcroot_patch)

       this%deadcroot_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C storage', &
             ptr_patch=this%deadcroot_storage_patch,  default='inactive')

       this%deadcroot_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 dead coarse root C transfer', &
             ptr_patch=this%deadcroot_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_STORAGE', units='gC14/m^2', &
             avgflag='A', long_name='C14 growth respiration storage', &
             ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_XFER', units='gC14/m^2', &
             avgflag='A', long_name='C14 growth respiration transfer', &
             ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%pool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL', units='gC14/m^2', &
             avgflag='A', long_name='C14 temporary photosynthate C pool', &
             ptr_patch=this%pool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_XSMRPOOL', units='gC14/m^2', &
             avgflag='A', long_name='C14 temporary photosynthate C pool', &
             ptr_patch=this%xsmrpool_patch)

       this%veg_trunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PFT_CTRUNC', units='gC14/m^2', &
             avgflag='A', long_name='C14 patch-level sink for C truncation', &
             ptr_patch=this%veg_trunc_patch)

       this%dispveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DISPVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 displayed veg carbon, excluding storage and cpool', &
             ptr_patch=this%dispveg_patch)

       this%storveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_STORVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 stored vegetation carbon, excluding cpool', &
             ptr_patch=this%storveg_patch)

       this%totveg_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTVEGC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total vegetation carbon, excluding cpool', &
             ptr_patch=this%totveg_patch)

       this%totpft_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTPFTC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total patch-level carbon, including cpool', &
             ptr_patch=this%totpft_patch)
    endif

    !-------------------------------
    ! C state variables - column
    !-------------------------------

    ! add history fields for all CLAMP CN variables



    if (carbon_type == 'c12') then


       !those variables are now ouput in betr
       this%decomp_pools_col(begc:endc,:) = spval
       do l  = 1, ndecomp_pools
          if(trim(decomp_cascade_con%decomp_pool_name_history(l))=='')exit
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_pools_vr_col(:,:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'

             call hist_addfld2d (fname=fieldname, units='gC/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_pools_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                avgflag='A', long_name=longname, &
                ptr_col=data1dptr)

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_pools_1m_col(:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
             longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr, default = 'inactive')
          endif
       end do

       if ( nlevdecomp_full > 1 ) then

          this%totlit_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTLITC_1m', units='gC/m^2', &
                avgflag='A', long_name='total litter carbon to 1 meter depth', &
                ptr_col=this%totlit_1m_col)

          this%totsom_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='TOTSOMC_1m', units='gC/m^2', &
                avgflag='A', long_name='total soil organic matter carbon to 1 meter depth', &
                ptr_col=this%totsom_1m_col)
       end if


       this%veg_trunc_col(begc:endc) = spval
       call hist_addfld1d (fname='COL_CTRUNC', units='gC/m^2',  &
             avgflag='A', long_name='column-level sink for C truncation', &
             ptr_col=this%veg_trunc_col, default='inactive')

       this%seed_col(begc:endc) = spval
       call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
             avgflag='A', long_name='pool for seeding new Patches', &
             ptr_col=this%seed_col, default='inactive')

       this%totlit_col(begc:endc) = spval
       call hist_addfld1d (fname='LITTERC', units='gC/m^2', &
             avgflag='A', long_name='litter C', &
             ptr_col=this%totlit_col)
       call hist_addfld1d (fname='TOTLITC', units='gC/m^2', &
             avgflag='A', long_name='total litter carbon', &
             ptr_col=this%totlit_col)

       this%totsom_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTSOMC', units='gC/m^2', &
             avgflag='A', long_name='total soil organic matter carbon', &
             ptr_col=this%totsom_col)

       call hist_addfld1d (fname='SOILC', units='gC/m^2', &
             avgflag='A', long_name='soil C', &
             ptr_col=this%totsom_col)

       this%totecosys_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
             avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool but excl product pools', &
             ptr_col=this%totecosys_col)

       this%totcol_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
             avgflag='A', long_name='total column carbon, incl veg and cpool but excl product pools', &
             ptr_col=this%totcol_col)

       this%prod10_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD10C', units='gC/m^2', &
             avgflag='A', long_name='10-yr wood product C', &
             ptr_col=this%prod10_col, default='inactive')

       this%prod100_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD100C', units='gC/m^2', &
             avgflag='A', long_name='100-yr wood product C', &
             ptr_col=this%prod100_col, default='inactive')

       this%prod1_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD1C', units='gC/m^2', &
             avgflag='A', long_name='1-yr crop product C', &
             ptr_col=this%prod1_col, default='inactive')

       this%totprod_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTPRODC', units='gC/m^2', &
             avgflag='A', long_name='total wood product C', &
             ptr_col=this%totprod_col, default='inactive')

       this%fuelc_col(begc:endc) = spval
       call hist_addfld1d (fname='FUELC', units='gC/m^2', &
             avgflag='A', long_name='fuel load', &
             ptr_col=this%fuelc_col, default='inactive')

    end if

    !-------------------------------
    ! C13 state variables - column
    !-------------------------------

    if ( carbon_type == 'c13' ) then


       this%decomp_pools_vr_col(begc:endc,:,:) = spval
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_pools_vr_col(:,:,l)
             fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_pools_col(:,l)
          fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC13/m^2', &
                avgflag='A', long_name=longname, &
                ptr_col=data1dptr)
       end do

       this%seed_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
             avgflag='A', long_name='C13 pool for seeding new Patches', &
             ptr_col=this%seed_col)

       this%veg_trunc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_COL_CTRUNC', units='gC13/m^2',  &
             avgflag='A', long_name='C13 column-level sink for C truncation', &
             ptr_col=this%veg_trunc_col)

       this%totlit_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTLITC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total litter carbon', &
             ptr_col=this%totlit_col)

       this%totsom_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTSOMC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total soil organic matter carbon', &
             ptr_col=this%totsom_col)

       if ( nlevdecomp_full > 1 ) then
          this%totlit_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C13_TOTLITC_1m', units='gC13/m^2', &
                avgflag='A', long_name='C13 total litter carbon to 1 meter', &
                ptr_col=this%totlit_1m_col)

          this%totsom_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C13_TOTSOMC_1m', units='gC13/m^2', &
                avgflag='A', long_name='C13 total soil organic matter carbon to 1 meter', &
                ptr_col=this%totsom_1m_col)
       endif

       this%totecosys_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool but excl product pools', &
             ptr_col=this%totecosys_col)

       this%totcol_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total column carbon, incl veg and cpool but excl product pools', &
             ptr_col=this%totcol_col)

       this%prod10_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD10C', units='gC13/m^2', &
             avgflag='A', long_name='C13 10-yr wood product C', &
             ptr_col=this%prod10_col)

       this%prod100_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD100C', units='gC13/m^2', &
             avgflag='A', long_name='C13 100-yr wood product C', &
             ptr_col=this%prod100_col)

       this%prod1_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD1C', units='gC13/m^2', &
             avgflag='A', long_name='C13 1-yr crop product C', &
             ptr_col=this%prod1_col)

       this%totprod_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTPRODC', units='gC13/m^2', &
             avgflag='A', long_name='C13 total wood product C', &
             ptr_col=this%totprod_col)

       this%seed_grc(begg:endg) = spval
       call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
            avgflag='A', long_name='C13 pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seed_grc)

       if (use_crop) then
          this%grain_patch(begp:endp) = spval
          call hist_addfld1d (fname='C13_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C13 grain C (does not equal yield)', &
               ptr_patch=this%grain_patch)
          this%cropseed_deficit_patch(begp:endp) = spval

          call hist_addfld1d (fname='C13_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C13 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseed_deficit_patch)
       end if
       
    endif

    !-------------------------------
    ! C14 state variables - column
    !-------------------------------

    if ( carbon_type == 'c14' ) then

       this%decomp_pools_vr_col(begc:endc,:,:) = spval
       do l = 1, ndecomp_pools
          if ( nlevdecomp_full > 1 ) then
             data2dptr => this%decomp_pools_vr_col(:,:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_vr'
             longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C (vertically resolved)'
             call hist_addfld2d (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, ptr_col=data2dptr)
          endif

          data1dptr => this%decomp_pools_col(:,l)
          fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C'
          longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C'
          call hist_addfld1d (fname=fieldname, units='gC14/m^2', &
                avgflag='A', long_name=longname, ptr_col=data1dptr)

          if ( nlevdecomp_full > 1 ) then
             data1dptr => this%decomp_pools_1m_col(:,l)
             fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//'C_1m'
             longname =  'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(l))//' C to 1 meter'
             call hist_addfld1d (fname=fieldname, units='gC/m^2', &
                   avgflag='A', long_name=longname, ptr_col=data1dptr, default='inactive')
          endif
       end do

       this%seed_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_SEEDC', units='gC14/m^2', &
             avgflag='A', long_name='C14 pool for seeding new Patches', &
             ptr_col=this%seed_col)

       this%veg_trunc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_COL_CTRUNC', units='gC14/m^2', &
             avgflag='A', long_name='C14 column-level sink for C truncation', &
             ptr_col=this%veg_trunc_col)

       this%totlit_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTLITC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total litter carbon', &
             ptr_col=this%totlit_col)

       this%totsom_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTSOMC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total soil organic matter carbon', &
             ptr_col=this%totsom_col)

       if ( nlevdecomp_full > 1 ) then       
          this%totlit_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C14_TOTLITC_1m', units='gC14/m^2', &
                avgflag='A', long_name='C14 total litter carbon to 1 meter', &
                ptr_col=this%totlit_1m_col)

          this%totsom_1m_col(begc:endc) = spval
          call hist_addfld1d (fname='C14_TOTSOMC_1m', units='gC14/m^2', &
                avgflag='A', long_name='C14 total soil organic matter carbon to 1 meter', &
                ptr_col=this%totsom_1m_col)
       endif

       this%totecosys_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTECOSYSC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total ecosystem carbon, incl veg but excl cpool but excl product pools', &
             ptr_col=this%totecosys_col)

       this%totcol_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTCOLC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total column carbon, incl veg and cpool but excl product pools', &
             ptr_col=this%totcol_col)

       this%prod10_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD10C', units='gC14/m^2', &
             avgflag='A', long_name='C14 10-yr wood product C', &
             ptr_col=this%prod10_col)

       this%prod100_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD100C', units='gC14/m^2', &
             avgflag='A', long_name='C14 100-yr wood product C', &
             ptr_col=this%prod100_col)

       this%prod1_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD1C', units='gC14/m^2', &
             avgflag='A', long_name='C14 1-yr crop product C', &
             ptr_col=this%prod1_col)

       this%totprod_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTPRODC', units='gC14/m^2', &
             avgflag='A', long_name='C14 total wood product C', &
             ptr_col=this%totprod_col)

       this%seed_grc(begg:endg) = spval
       call hist_addfld1d (fname='C14_SEEDC', units='gC14/m^2', &
            avgflag='A', long_name='C14 pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seed_grc)

       if (use_crop) then
          this%grain_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_GRAINC', units='gC/m^2', &
               avgflag='A', long_name='C14 grain C (does not equal yield)', &
               ptr_patch=this%grain_patch)
          this%cropseed_deficit_patch(begp:endp) = spval
          call hist_addfld1d (fname='C14_CROPSEEDC_DEFICIT', units='gC/m^2', &
               avgflag='A', long_name='C14 C used for crop seed that needs to be repaid', &
               ptr_patch=this%cropseed_deficit_patch)
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
                this%leaf_patch(p)         = 0._r8
                this%leaf_storage_patch(p) = 0._r8
             else
                if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                   this%leaf_patch(p)         = 1._r8 * ratio
                   this%leaf_storage_patch(p) = 0._r8
                else if (veg_pp%itype(p) >= npcropmin) then ! prognostic crop types
                   this%leaf_patch(p) = 0._r8
                   this%leaf_storage_patch(p) = 0._r8
                else
                   this%leaf_patch(p) = 0._r8
                   this%leaf_storage_patch(p) = 1._r8 * ratio
                end if
             end if
             this%leaf_xfer_patch(p) = 0._r8

             this%froot_patch(p)            = 0._r8 
             this%froot_storage_patch(p)    = 0._r8 
             this%froot_xfer_patch(p)       = 0._r8 

             this%livestem_patch(p)         = 0._r8 
             this%livestem_storage_patch(p) = 0._r8 
             this%livestem_xfer_patch(p)    = 0._r8 

             if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
                this%deadstem_patch(p) = 0.1_r8 * ratio
             else
                this%deadstem_patch(p) = 0._r8 
             end if
             this%deadstem_storage_patch(p)  = 0._r8 
             this%deadstem_xfer_patch(p)     = 0._r8

             if (nu_com .ne. 'RD') then
                ! ECA competition calculate root NP uptake as a function of fine root biomass
                ! better to initialize root CNP pools with a non-zero value
                if (veg_pp%itype(p) .ne. noveg) then
                   if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                      this%leaf_patch(p) = 20._r8 * ratio
                      this%leaf_storage_patch(p) = 0._r8
                      this%froot_patch(p) = 20._r8 * ratio
                      this%froot_storage_patch(p) = 0._r8
                   else
                      this%leaf_patch(p) = 0._r8 
                      this%leaf_storage_patch(p) = 20._r8 * ratio
                      this%froot_patch(p) = 0._r8
                      this%froot_storage_patch(p) = 20._r8 * ratio
                   end if
                end if
             end if

             this%livecroot_patch(p)         = 0._r8 
             this%livecroot_storage_patch(p) = 0._r8 
             this%livecroot_xfer_patch(p)    = 0._r8 

             this%deadcroot_patch(p)         = 0._r8 
             this%deadcroot_storage_patch(p) = 0._r8 
             this%deadcroot_xfer_patch(p)    = 0._r8 

             this%gresp_storage_patch(p)      = 0._r8 
             this%gresp_xfer_patch(p)         = 0._r8 

             this%pool_patch(p)              = 0._r8 
             this%xsmrpool_patch(p)           = 0._r8 
             this%veg_trunc_patch(p)             = 0._r8 
             this%dispveg_patch(p)           = 0._r8 
             this%storveg_patch(p)           = 0._r8 
             this%totpft_patch(p)            = 0._r8 
             this%woodc_patch(p)              = 0._r8

             if ( crop_prog )then
                this%grain_patch(p)            = 0._r8 
                this%grain_storage_patch(p)    = 0._r8 
                this%grain_xfer_patch(p)       = 0._r8 
                this%cropseed_deficit_patch(p) = 0._r8
             end if

             ! calculate totvegc explicitly so that it is available for the isotope 
             ! code on the first time step.

             this%totveg_patch(p) = &
                  this%leaf_patch(p)              + &
                  this%leaf_storage_patch(p)      + &
                  this%leaf_xfer_patch(p)         + &
                  this%froot_patch(p)             + &
                  this%froot_storage_patch(p)     + &
                  this%froot_xfer_patch(p)        + &
                  this%livestem_patch(p)          + &
                  this%livestem_storage_patch(p)  + &
                  this%livestem_xfer_patch(p)     + &
                  this%deadstem_patch(p)          + &
                  this%deadstem_storage_patch(p)  + &
                  this%deadstem_xfer_patch(p)     + &
                  this%livecroot_patch(p)         + &
                  this%livecroot_storage_patch(p) + &
                  this%livecroot_xfer_patch(p)    + &
                  this%deadcroot_patch(p)         + &
                  this%deadcroot_storage_patch(p) + &
                  this%deadcroot_xfer_patch(p)    + &
                  this%gresp_storage_patch(p)      + &
                  this%gresp_xfer_patch(p)         + &
                  this%pool_patch(p)

             if ( crop_prog )then
                this%totveg_patch(p) =  this%totveg_patch(p) + &
                     this%grain_patch(p)                            + &
                     this%grain_storage_patch(p)                    + &
                     this%grain_xfer_patch(p)
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
                      this%decomp_pools_vr_col(c,j,k) = decomp_cascade_con%initial_stock(k)
                   else
                      this%decomp_pools_vr_col(c,j,k) = 0._r8
                   endif
                end do
                this%soil_trunc_vr_col(c,j) = 0._r8
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_pools_vr_col(c,j,k) = 0._r8
                   end do
                   this%soil_trunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%decomp_pools_col(c,1:ndecomp_pools)    = decomp_cascade_con%initial_stock(1:ndecomp_pools)
             this%decomp_pools_1m_col(c,1:ndecomp_pools) = decomp_cascade_con%initial_stock(1:ndecomp_pools)

          else

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   this%decomp_pools_vr_col(c,j,k) = c12_carbonstate_vars%decomp_pools_vr_col(c,j,k) * ratio
                end do
                this%soil_trunc_vr_col(c,j) = c12_carbonstate_vars%soil_trunc_vr_col(c,j) * ratio
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_pools_vr_col(c,j,k) = 0._r8
                   end do
                   this%soil_trunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%cwd_col(c) = c12_carbonstate_vars%cwd_col(c) * ratio
             do k = 1, ndecomp_pools
                this%decomp_pools_col(c,k)    = c12_carbonstate_vars%decomp_pools_col(c,k) * ratio
                this%decomp_pools_1m_col(c,k) = c12_carbonstate_vars%decomp_pools_1m_col(c,k) * ratio
             end do

          endif

          this%cwd_col(c)       = 0._r8
          this%veg_trunc_col(c)     = 0._r8
          this%totlit_col(c)    = 0._r8
          this%totsom_col(c)    = 0._r8
          this%totlit_1m_col(c) = 0._r8
          this%totsom_1m_col(c) = 0._r8
          this%totecosys_col(c) = 0._r8
          this%totcol_col(c)    = 0._r8

          ! dynamic landcover state variables
          this%seed_col(c)      = 0._r8
          this%prod10_col(c)    = 0._r8
          this%prod100_col(c)   = 0._r8
          this%prod1_col(c)     = 0._r8
          this%totprod_col(c)   = 0._r8

       end if

    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics
    
    do fc = 1,num_special_col
       c = special_col(fc)

       this%seed_col(c)      = 0._r8
       this%prod10_col(c)    = 0._r8
       this%prod100_col(c)   = 0._r8
       this%prod1_col(c)     = 0._r8
       this%totprod_col(c)   = 0._r8
    end do

    do g = bounds%begg, bounds%endg
       this%seed_grc(g) = 0._r8
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
               interpinic_flag='interp', readvar=readvar, data=this%leaf_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 
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
               interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 
       end if

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totveg_patch) 
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
                   this%grain_patch(i)            = c12_carbonstate_vars%grain_patch(i)         * c3_r2
                   this%grain_storage_patch(i)    = c12_carbonstate_vars%grain_storage_patch(i) * c3_r2
                   this%grain_xfer_patch(i)       = c12_carbonstate_vars%grain_xfer_patch(i)    * c3_r2
                   this%dispveg_patch(i)          = c12_carbonstate_vars%dispveg_patch(i)       * c3_r2
                   this%storveg_patch(i)          = c12_carbonstate_vars%storveg_patch(i)       * c3_r2
                   this%totveg_patch(i)           = c12_carbonstate_vars%totveg_patch(i)        * c3_r2
                   this%totpft_patch(i)           = c12_carbonstate_vars%totpft_patch(i)        * c3_r2
                   this%woodc_patch(i)             = c12_carbonstate_vars%woodc_patch(i)          * c3_r2
                else
                   this%grain_patch(i)            = c12_carbonstate_vars%grain_patch(i)         * c4_r2
                   this%grain_storage_patch(i)    = c12_carbonstate_vars%grain_storage_patch(i) * c4_r2
                   this%grain_xfer_patch(i)       = c12_carbonstate_vars%grain_xfer_patch(i)    * c4_r2
                   this%dispveg_patch(i)          = c12_carbonstate_vars%dispveg_patch(i)       * c4_r2
                   this%storveg_patch(i)          = c12_carbonstate_vars%storveg_patch(i)       * c4_r2
                   this%totveg_patch(i)           = c12_carbonstate_vars%totveg_patch(i)        * c4_r2
                   this%totpft_patch(i)           = c12_carbonstate_vars%totpft_patch(i)        * c4_r2
                   this%woodc_patch(i)             = c12_carbonstate_vars%woodc_patch(i)          * c4_r2
                end if
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_patch)
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leaf_patch(i) = c12_carbonstate_vars%leaf_patch(i) * c3_r2
                else
                   this%leaf_patch(i) = c12_carbonstate_vars%leaf_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leaf_storage_patch(i) = c12_carbonstate_vars%leaf_storage_patch(i) * c3_r2
                else
                   this%leaf_storage_patch(i) = c12_carbonstate_vars%leaf_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leaf_xfer_patch(i) = c12_carbonstate_vars%leaf_xfer_patch(i) * c3_r2
                else
                   this%leaf_xfer_patch(i) = c12_carbonstate_vars%leaf_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%froot_patch(i) = c12_carbonstate_vars%froot_patch(i) * c3_r2
                else
                   this%froot_patch(i) = c12_carbonstate_vars%froot_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%froot_storage_patch(i) = c12_carbonstate_vars%froot_storage_patch(i) * c3_r2
                else
                   this%froot_storage_patch(i) = c12_carbonstate_vars%froot_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%froot_xfer_patch(i) = c12_carbonstate_vars%froot_xfer_patch(i) * c3_r2
                else
                   this%froot_xfer_patch(i) = c12_carbonstate_vars%froot_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestem_patch(i) = c12_carbonstate_vars%livestem_patch(i) * c3_r2
                else
                   this%livestem_patch(i) = c12_carbonstate_vars%livestem_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestem_storage_patch(i) = c12_carbonstate_vars%livestem_storage_patch(i) * c3_r2
                else
                   this%livestem_storage_patch(i) = c12_carbonstate_vars%livestem_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestem_xfer_patch(i) = c12_carbonstate_vars%livestem_xfer_patch(i) * c3_r2
                else
                   this%livestem_xfer_patch(i) = c12_carbonstate_vars%livestem_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstem_patch(i) = c12_carbonstate_vars%deadstem_patch(i) * c3_r2
                else
                   this%deadstem_patch(i) = c12_carbonstate_vars%deadstem_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstem_storage_patch(i) = c12_carbonstate_vars%deadstem_storage_patch(i) * c3_r2
                else
                   this%deadstem_storage_patch(i) = c12_carbonstate_vars%deadstem_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstem_xfer_patch(i) = c12_carbonstate_vars%deadstem_xfer_patch(i) * c3_r2
                else
                   this%deadstem_xfer_patch(i) = c12_carbonstate_vars%deadstem_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecroot_patch(i) = c12_carbonstate_vars%livecroot_patch(i) * c3_r2
                else
                   this%livecroot_patch(i) = c12_carbonstate_vars%livecroot_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecroot_storage_patch(i) = c12_carbonstate_vars%livecroot_storage_patch(i) * c3_r2
                else
                   this%livecroot_storage_patch(i) = c12_carbonstate_vars%livecroot_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecroot_xfer_patch(i) = c12_carbonstate_vars%livecroot_xfer_patch(i) * c3_r2
                else
                   this%livecroot_xfer_patch(i) = c12_carbonstate_vars%livecroot_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcroot_patch(i) = c12_carbonstate_vars%deadcroot_patch(i) * c3_r2
                else
                   this%deadcroot_patch(i) = c12_carbonstate_vars%deadcroot_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcroot_storage_patch(i) = c12_carbonstate_vars%deadcroot_storage_patch(i) * c3_r2
                else
                   this%deadcroot_storage_patch(i) = c12_carbonstate_vars%deadcroot_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcroot_xfer_patch(i) = c12_carbonstate_vars%deadcroot_xfer_patch(i) * c3_r2
                else
                   this%deadcroot_xfer_patch(i) = c12_carbonstate_vars%deadcroot_xfer_patch(i) * c4_r2
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
               interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%cpool with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%pool_patch(i) = c12_carbonstate_vars%pool_patch(i) * c3_r2
                else
                   this%pool_patch(i) = c12_carbonstate_vars%pool_patch(i) * c4_r2
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
               interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%ctrunc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%veg_trunc_patch(i) = c12_carbonstate_vars%veg_trunc_patch(i) * c3_r2
                else
                   this%veg_trunc_patch(i) = c12_carbonstate_vars%veg_trunc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totveg_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing carbonstate_vars %totvegc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%totveg_patch(i) = c12_carbonstate_vars%totveg_patch(i) * c3_r2
                else
                   this%totveg_patch(i) = c12_carbonstate_vars%totveg_patch(i) * c4_r2
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
               interpinic_flag='interp', readvar=readvar, data=this%leaf_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leaf_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leaf_patch(i) /= spval .and. &
                     .not. isnan(this%leaf_patch(i)) ) then
                   this%leaf_patch(i) = c12_carbonstate_vars%leaf_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leaf_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leaf_storage_patch(i) /= spval .and. &
                     .not. isnan(this%leaf_storage_patch(i)) ) then
                   this%leaf_storage_patch(i) = c12_carbonstate_vars%leaf_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leaf_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leaf_xfer_patch(i) /= spval .and. .not. isnan(this%leaf_xfer_patch(i)) ) then
                   this%leaf_xfer_patch(i) = c12_carbonstate_vars%leaf_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%froot_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%froot_patch(i) /= spval .and. &
                     .not. isnan(this%froot_patch(i)) ) then
                   this%froot_patch(i) = c12_carbonstate_vars%froot_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%froot_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%froot_storage_patch(i) /= spval .and. &
                     .not. isnan(this%froot_storage_patch(i)) ) then
                   this%froot_storage_patch(i) = c12_carbonstate_vars%froot_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%froot_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%froot_xfer_patch(i) /= spval .and. &
                     .not. isnan(this%froot_xfer_patch(i)) ) then
                   this%froot_xfer_patch(i) = c12_carbonstate_vars%froot_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestem_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestem_patch(i) /= spval .and. .not. isnan(this%livestem_patch(i)) ) then
                   this%livestem_patch(i) = c12_carbonstate_vars%livestem_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestem_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestem_storage_patch(i) /= spval .and. .not. isnan(this%livestem_storage_patch(i)) ) then
                   this%livestem_storage_patch(i) = c12_carbonstate_vars%livestem_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestem_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestem_xfer_patch(i) /= spval .and. .not. isnan(this%livestem_xfer_patch(i)) ) then
                   this%livestem_xfer_patch(i) = c12_carbonstate_vars%livestem_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstem_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstem_patch(i) /= spval .and. .not. isnan(this%deadstem_patch(i)) ) then
                   this%deadstem_patch(i) = c12_carbonstate_vars%deadstem_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstem_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstem_storage_patch(i) /= spval .and. .not. isnan(this%deadstem_storage_patch(i)) ) then
                   this%deadstem_storage_patch(i) = c12_carbonstate_vars%deadstem_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstem_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstem_xfer_patch(i) /= spval .and. .not. isnan(this%deadstem_xfer_patch(i)) ) then
                   this%deadstem_xfer_patch(i) = c12_carbonstate_vars%deadstem_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecroot_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecroot_patch(i) /= spval .and. .not. isnan(this%livecroot_patch(i)) ) then
                   this%livecroot_patch(i) = c12_carbonstate_vars%livecroot_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecroot_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecroot_storage_patch(i) /= spval .and. .not. isnan(this%livecroot_storage_patch(i)) ) then
                   this%livecroot_storage_patch(i) = c12_carbonstate_vars%livecroot_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecroot_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecroot_xfer_patch(i) /= spval .and. .not. isnan(this%livecroot_xfer_patch(i)) ) then
                   this%livecroot_xfer_patch(i) = c12_carbonstate_vars%livecroot_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcroot_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcroot_patch(i) /= spval .and. .not. isnan(this%deadcroot_patch(i)) ) then
                   this%deadcroot_patch(i) = c12_carbonstate_vars%deadcroot_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcroot_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcroot_storage_patch(i) /= spval .and. .not. isnan(this%deadcroot_storage_patch(i)) ) then
                   this%deadcroot_storage_patch(i) = c12_carbonstate_vars%deadcroot_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog) 'initializing this%deadcroot_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcroot_xfer_patch(i) /= spval .and. .not. isnan(this%deadcroot_xfer_patch(i)) ) then
                   this%deadcroot_xfer_patch(i) = c12_carbonstate_vars%deadcroot_xfer_patch(i) * c14ratio
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
               interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%pool_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%pool_patch(i) /= spval .and. .not. isnan(this%pool_patch(i)) ) then
                   this%pool_patch(i) = c12_carbonstate_vars%pool_patch(i) * c14ratio
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
               interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%veg_trunc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%veg_trunc_patch(i) /= spval .and. .not. isnan(this%veg_trunc_patch(i)) ) then
                   this%veg_trunc_patch(i) = c12_carbonstate_vars%veg_trunc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totveg_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%totveg_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%totveg_patch(i) /= spval .and. .not. isnan(this%totveg_patch(i)) ) then
                   this%totveg_patch(i) = c12_carbonstate_vars%totveg_patch(i) * c14ratio
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
                  interpinic_flag='interp', readvar=readvar, data=this%grain_patch)

             call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C storage', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grain_storage_patch)

             call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer', xtype=ncd_double,  &
                  dim1name='pft', long_name='grain C transfer', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%grain_xfer_patch)

             call restartvar(ncid=ncid, flag=flag, varname='cropseedc_deficit', xtype=ncd_double,  &
                  dim1name='pft', long_name='pool for seeding new crop growth', units='gC/m2', &
                  interpinic_flag='interp', readvar=readvar, data=this%cropseed_deficit_patch)
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
    end if
    if(is_active_betr_bgc)then
      if (carbon_type == 'c12') then
        call restartvar(ncid=ncid, flag=flag, varname='totblgc', xtype=ncd_double,  &
           dim1name='column', long_name='', units='', &
           interpinic_flag='interp', readvar=readvar, data=this%totblg_col)

        call restartvar(ncid=ncid, flag=flag, varname='cwdc', xtype=ncd_double,  &
           dim1name='column', long_name='', units='', &
           interpinic_flag='interp', readvar=readvar, data=this%cwd_col)
      endif
    endif
    if (carbon_type == 'c12') then
       if (use_vertsoilc) then
          ptr2d => this%soil_trunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc_vr', xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%soil_trunc_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
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
            interpinic_flag='interp', readvar=readvar, data=this%seed_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlit_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1_col)
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='totsomc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totsom_col) 
    end if

    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c13' ) then
       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
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
             write(iulog,*) 'initializing this%decomp_pools_vr_col with atmospheric c13 value for: '//varname
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) /= spval .and. &
                        .not. isnan(this%decomp_pools_vr_col(i,j,k)) ) then
                         this%decomp_pools_vr_col(i,j,k) = c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seed_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%seed_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%seed_col(i)) ) then
             this%seed_col(i) = c12_carbonstate_vars%seed_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       if (use_vertsoilc) then
          ptr2d => this%soil_trunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%soil_trunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlit_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%totlit_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%totlit_col(i) ) ) then
             this%totlit_col(i) = c12_carbonstate_vars%totlit_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%totcol_col(i) /= spval .and. &
               .not. isnan (c12_carbonstate_vars%totcol_col(i) ) ) then
             this%totcol_col(i) = c12_carbonstate_vars%totcol_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod10_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod10_col(i) ) ) then
             this%prod10_col(i) = c12_carbonstate_vars%prod10_col(i) * c3_r2
          endif
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod100_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod100_col(i) ) ) then
             this%prod100_col(i) = c12_carbonstate_vars%prod100_col(i) * c3_r2
          endif
       end if
    endif

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1_col)
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod1_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod1_col(i) ) ) then
             this%prod1_col(i) = c12_carbonstate_vars%prod1_col(i) * c3_r2
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
             ptr2d => this%decomp_pools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_pools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%decomp_pools_vr_col with atmospheric c14 value for: '//trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) /= spval .and. &
                        .not. isnan(c12_carbonstate_vars%decomp_pools_vr_col(i,j,k)) ) then
                         this%decomp_pools_vr_col(i,j,k) = c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) * c3_r2
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
            interpinic_flag='interp', readvar=readvar, data=this%seed_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%seed_col with atmospheric c14 value'
          do i = bounds%begc,bounds%endc
             if (c12_carbonstate_vars%seed_col(i) /= spval .and. &
                  .not. isnan(c12_carbonstate_vars%seed_col(i)) ) then
                this%seed_col(i) = c12_carbonstate_vars%seed_col(i) * c14ratio
             endif
          end do
       end if
    end if

    if ( carbon_type == 'c14' ) then
       if (use_vertsoilc) then 
          ptr2d => this%soil_trunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%soil_trunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlit_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%totlit_col with atmospheric c14 value'
          if (c12_carbonstate_vars%totlit_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%totlit_col(i)) ) then
             this%totlit_col(i) = c12_carbonstate_vars%totlit_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%totcol_col with atmospheric c14 value'
          if (c12_carbonstate_vars%totcol_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%totcol_col(i)) ) then
             this%totcol_col(i) = c12_carbonstate_vars%totcol_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod10_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod10_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod10_col(i)) ) then
             this%prod10_col(i) = c12_carbonstate_vars%prod10_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod100_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod100_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod100_col(i)) ) then
             this%prod100_col(i) = c12_carbonstate_vars%prod100_col(i) * c14ratio
          endif
       end if
    endif

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1_col)
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod1_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod1_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod1_col(i)) ) then
             this%prod1_col(i) = c12_carbonstate_vars%prod1_col(i) * c14ratio
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
                    this%decomp_pools_vr_col(c,j,k) = this%decomp_pools_vr_col(c,j,k) * m
                 end do
              end do
           end do
           do i = bounds%begp, bounds%endp
              if (exit_spinup) then 
                 m_veg = spinup_mortality_factor
              else if (enter_spinup) then 
                 m_veg = 1._r8 / spinup_mortality_factor
              end if
              this%deadstem_patch(i)  = this%deadstem_patch(i) * m_veg
              this%deadcroot_patch(i) = this%deadcroot_patch(i) * m_veg
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

          this%leaf_patch(i)              = value_patch
          this%leaf_storage_patch(i)      = value_patch
          this%leaf_xfer_patch(i)         = value_patch
          this%froot_patch(i)             = value_patch
          this%froot_storage_patch(i)     = value_patch
          this%froot_xfer_patch(i)        = value_patch
          this%livestem_patch(i)          = value_patch
          this%livestem_storage_patch(i)  = value_patch
          this%livestem_xfer_patch(i)     = value_patch
          this%deadstem_patch(i)          = value_patch
          this%deadstem_storage_patch(i)  = value_patch
          this%deadstem_xfer_patch(i)     = value_patch
          this%livecroot_patch(i)         = value_patch
          this%livecroot_storage_patch(i) = value_patch
          this%livecroot_xfer_patch(i)    = value_patch
          this%deadcroot_patch(i)         = value_patch
          this%deadcroot_storage_patch(i) = value_patch
          this%deadcroot_xfer_patch(i)    = value_patch
          this%gresp_storage_patch(i)      = value_patch
          this%gresp_xfer_patch(i)         = value_patch
          this%pool_patch(i)              = value_patch
          this%xsmrpool_patch(i)           = value_patch
          this%veg_trunc_patch(i)             = value_patch
          this%dispveg_patch(i)           = value_patch
          this%storveg_patch(i)           = value_patch
          this%totveg_patch(i)            = value_patch
          this%totpft_patch(i)            = value_patch
          this%woodc_patch(i)              = value_patch
          this%totveg_abg_patch(i)        = value_patch
       end do
       
       if ( crop_prog ) then
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%grain_patch(i)            = value_patch
             this%grain_storage_patch(i)    = value_patch
             this%grain_xfer_patch(i)       = value_patch
             this%cropseed_deficit_patch(i) = value_patch
          end do
       endif
    endif ! .not. use_fates
    
    do fi = 1,num_column
       i = filter_column(fi)
       this%cwd_col(i)       = value_column
       this%veg_trunc_col(i)     = value_column
       this%totlit_col(i)    = value_column
       this%totsom_col(i)    = value_column
       this%totecosys_col(i) = value_column
       this%totcol_col(i)    = value_column
       this%rootc_col(i)      = value_column
       this%totveg_col(i)    = value_column
       this%leafc_col(i)      = value_column
       this%deadstemc_col(i)  = value_column
       this%fuelc_col(i)      = value_column
       this%fuelc_crop_col(i) = value_column
       this%totlit_1m_col(i) = value_column
       this%totsom_1m_col(i) = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%soil_trunc_vr_col(i,j) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_pools_col(i,k) = value_column
          this%decomp_pools_1m_col(i,k) = value_column
       end do
    end do

    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_pools_vr_col(i,j,k) = value_column
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
       this%dispveg_patch(p)   = 0._r8
       this%storveg_patch(p)   = 0._r8
       this%totpft_patch(p)    = 0._r8
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
       this%dispveg_patch(p) =        &
            this%leaf_patch(p)      + &
            this%froot_patch(p)     + &
            this%livestem_patch(p)  + &
            this%deadstem_patch(p)  + &
            this%livecroot_patch(p) + &
            this%deadcroot_patch(p)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
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
            this%deadcroot_xfer_patch(p)    + &
            this%gresp_storage_patch(p)      + &
            this%gresp_xfer_patch(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%storveg_patch(p) =            &
               this%storveg_patch(p)       + &
               this%grain_storage_patch(p) + &
               this%grain_xfer_patch(p)

          this%dispveg_patch(p) =            &
               this%dispveg_patch(p)       + &
               this%grain_patch(p)
       end if

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totveg_patch(p) = &
            this%dispveg_patch(p) + &
            this%storveg_patch(p)

       ! total pft-level carbon, including xsmrpool, ctrunc
       this%totpft_patch(p) = &
            this%totveg_patch(p) + &
            this%xsmrpool_patch(p) + &
            this%veg_trunc_patch(p)

       ! (WOODC) - wood C
       this%woodc_patch(p) = &
            this%deadstem_patch(p)    + &
            this%livestem_patch(p)    + &
            this%deadcroot_patch(p)   + &
            this%livecroot_patch(p)

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


    call p2c(bounds, num_soilc, filter_soilc, &
         this%totpft_patch(bounds%begp:bounds%endp), &
         this%totpft_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totveg_patch(bounds%begp:bounds%endp), &
         this%totveg_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totveg_abg_patch(bounds%begp:bounds%endp), &
         this%totveg_abg_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%cropseed_deficit_patch(bounds%begp:bounds%endp), &
         cropseedc_deficit_col(bounds%begc:bounds%endc))

    ! column level summary

     nlev = nlevdecomp
     if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full


      ! vertically integrate each of the decomposing C pools
      do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_pools_col(c,l) = 0._r8
       end do
      end do
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

       ! vertically integrate each of the decomposing C pools to 1 meter
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

       ! total litter carbon in the top meter (TOTLITC_1m)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%totlit_1m_col(c)         = 0._r8
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

       ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%totsom_1m_col(c) = 0._r8
       end do
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
    
      ! total litter carbon (TOTLITC)
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totlit_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totlit_col(c) = &
                  this%totlit_col(c) + &
                  this%decomp_pools_col(c,l)
          end do
       endif
      end do

      ! total soil organic matter carbon (TOTSOMC)
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totsom_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_soil(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%totsom_col(c) = &
                  this%totsom_col(c) + &
                  this%decomp_pools_col(c,l)
          end do
       end if
      end do


      ! coarse woody debris carbon
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%cwd_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_cwd(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%cwd_col(c) = &
                  this%cwd_col(c) + &
                  this%decomp_pools_col(c,l)
          end do
       end if
      end do

    ! truncation carbon
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%veg_trunc_col(c) = 0._r8
    end do
    do j = 1, nlev
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%veg_trunc_col(c) = &
               this%veg_trunc_col(c) + &
               this%soil_trunc_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! total product carbon
       this%totprod_col(c) = &
            this%prod10_col(c)  + &
            this%prod100_col(c) + &
            this%prod1_col(c) 

       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosys_col(c) = &
            this%cwd_col(c)     + &
            this%totlit_col(c)  + &
            this%totsom_col(c)  + &
            this%totprod_col(c) + &
            this%totveg_col(c)

       ! total column carbon, including veg and cpool (TOTCOLC)
       ! adding col_ctrunc, seedc
       this%totcol_col(c) = &
            this%totpft_col(c)  + &
            this%cwd_col(c)     + &
            this%totlit_col(c)  + &
            this%totsom_col(c)  + &
            this%prod1_col(c)   + &
            this%veg_trunc_col(c)   + &
            cropseedc_deficit_col(c)
            
       this%totabg_col(c) = &
            this%totpft_col(c)  + &
            this%totprod_col(c) + &
            this%seed_col(c)    + &
            this%veg_trunc_col(c)    
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
    use CNComputeSeedMod   , only : ComputeSeedAmounts
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
    real(r8)                    :: seed_leaf_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leaf_storage_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leaf_xfer_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstem_patch(bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_cflux(bounds%begp:bounds%endp)
    real(r8)                    :: deadstem_patch_temp(bounds%begp:bounds%endp)

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
         leaf_patch                 = this%leaf_patch(begp:endp)           , &
         leaf_storage_patch         = this%leaf_storage_patch(begp:endp)   , &
         leaf_xfer_patch            = this%leaf_xfer_patch(begp:endp)      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leaf_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leaf_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leaf_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstem_patch(begp:endp))

    ! 1) LEAF_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%leaf_patch   (begp:endp)           , &
         flux_out_grc_area = conv_cflux       (begp:endp)           , &
         seed              = seed_leaf_patch (begp:endp)           , &
         seed_addition     = dwt_leafc_seed   (begp:endp))


    ! 2) LEAF_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%leaf_storage_patch   (begp:endp)   , &
         flux_out_grc_area = conv_cflux               (begp:endp)   , &
         seed              = seed_leaf_storage_patch (begp:endp)   , &
         seed_addition     = dwt_leafc_seed           (begp:endp))

    ! 3) LEAF_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%leaf_xfer_patch   (begp:endp)      , &
         flux_out_grc_area = conv_cflux            (begp:endp)      , &
         seed              = seed_leaf_xfer_patch (begp:endp)      , &
         seed_addition     = dwt_leafc_seed        (begp:endp))

    ! 4) FROOT_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_patch(begp:endp)             , &
         flux_out_col_area = dwt_frootc_to_litter(begp:endp))
    
    ! 5) FROOT_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 6) FROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 7) LIVESTEM_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_patch(begp:endp)          , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 8) LIVESTEM_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 9) LIVESTEM_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 10) PROD10_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstem_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod10                          , &
         var                        = deadstem_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstem_patch    (begp:endp) )

    ! 11) PROD100_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstem_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod100                         , &
         var                        = deadstem_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_cflux            (begp:endp) , &
         flux2_out                  = wood_product_cflux      (begp:endp) , &
         seed                       = seed_deadstem_patch    (begp:endp))

    ! 12) DEADSTEM_PATCH
    wood_product_cflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pconv                          , &
         var                        = this%deadstem_patch(begp:endp) , &
         flux1_out                  = conv_cflux              (begp:endp) , &
         flux2_out                  = wood_product_cflux       (begp:endp) , &
         seed                       = seed_deadstem_patch    (begp:endp) , &
         seed_addition              = dwt_deadstemc_seed     (begp:endp))

    ! 13) DEADSTEM_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 14) DEADSTEM_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 15) LIVECROOT_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_livecrootc_to_litter(begp:endp))

    ! 16) LIVECROOT_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 17) LIVECROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 18) DEADCROOT_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_deadcrootc_to_litter(begp:endp))

    ! 19) DEADCROOT_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_xfer_patch(begp:endp)    , &
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
         var               = this%pool_patch(begp:endp)              , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 24) XSMRPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%xsmrpool_patch(begp:endp)           , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 25) veg_trunc_patch
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%veg_trunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 26) dispveg_patch
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%dispveg_patch(begp:endp))

    ! 27) storveg_patch
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%storveg_patch(begp:endp))

    ! 28) totveg_patch
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totveg_patch(begp:endp))

    ! 29) totpft_patch
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totpft_patch(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call patch_state_updater%update_patch_state(         &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = this%cropseed_deficit_patch(begp:endp) , &
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

    this%dyn_bal_adjustments_col(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call column_state_updater%update_column_state_no_special_handling( &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = this%decomp_pools_vr_col(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))

          this%dyn_bal_adjustments_col(begc:endc) = &
               this%dyn_bal_adjustments_col(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%soil_trunc_vr_col(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_bal_adjustments_col(begc:endc) = &
            this%dyn_bal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do

  end subroutine DynamicColumnAdjustments

end module CNCarbonStateType
