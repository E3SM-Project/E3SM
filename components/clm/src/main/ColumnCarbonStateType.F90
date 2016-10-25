module ColumnCarbonStateType
  

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval, spval
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak

  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, crop_prog, nlevdecomp
  use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi, zsoi
  use landunit_varcon        , only : istcrop 
  use clm_varctl             , only : iulog, use_vertsoilc, use_cndv, spinup_state 
  use decompMod              , only : bounds_type
!!!!!*******DW*******  use CNStateType            , only : cnstate_type
  use pftvarcon              , only : npcropmin
  use CNDecompCascadeConType , only : decomp_cascade_con
  use EcophysConType         , only : ecophyscon
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
!  use subgridAveMod          , only : p2c
!!!!!  use LandunitType           , only : lun                
! use ColumnType             , only : col                
  use PatchType              , only : pft
  use clm_varctl             , only : nu_com
  

!  moved from CNStateType.F90

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevsoifl, nlevsoi, crop_prog
  use clm_varpar     , only : ndecomp_cascade_transitions, nlevdecomp, nlevdecomp_full, more_vertlayers  
  use clm_varcon     , only : spval, ispval, c14ratio, grlnd
  use landunit_varcon, only : istsoil, istcrop
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, crop_prog 
  use clm_varctl     , only : use_vertsoilc, use_c14, use_cn 
  use clm_varctl     , only : iulog, fsurdat
!!!!!!  use LandunitType   , only : lun                
!  use ColumnType     , only : col                
!  use PatchType      , only : pft   



  implicit none
  save
  private

!!!!!!!! *********DW********* checkDates should be inculded at patch or veg level?????
  ! !PRIVATE MEMBER FUNCTIONS: 
 !!!!!!!! private :: checkDates


  type, public :: column_carbon_state
!   Moved from CNCarbonStateType.F90
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
     real(r8), pointer :: frootc_col               (:)     ! col (gC/m2) column-level C pool for fine root
     real(r8), pointer :: seedc_col                (:)     ! col (gC/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod10c_col              (:)     ! col (gC/m2) wood product C pool, 10-year lifespan
     real(r8), pointer :: prod100c_col             (:)     ! col (gC/m2) wood product C pool, 100-year lifespan
     real(r8), pointer :: totprodc_col             (:)     ! col (gC/m2) total wood product C
     ! pools for crop harvest
     real(r8), pointer :: prod1c_col               (:)     ! col (gC/m2) crop product C pool, 1-year lifespan

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: totvegcc_col             (:)     ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)
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

     ! Balance checks
     real(r8), pointer :: begcb_col                (:)     ! patch carbon mass, beginning of time step (gC/m**2)
     real(r8), pointer :: endcb_col                (:)     ! patch carbon mass, end of time step (gC/m**2)
     real(r8), pointer :: errcb_col                (:)     ! patch carbon balance error for the timestep (gC/m**2)
     
     real(r8), pointer :: totpftc_beg_col(:)
     real(r8), pointer :: cwdc_beg_col(:)
     real(r8), pointer :: totlitc_beg_col(:)
     real(r8), pointer :: totsomc_beg_col(:)
     
     real(r8), pointer :: totpftc_end_col(:)
     real(r8), pointer :: cwdc_end_col(:)
     real(r8), pointer :: totlitc_end_col(:)
     real(r8), pointer :: totsomc_end_col(:)

!   From CNStateType.F90

     real(r8) , pointer :: gdp_lf_col                  (:)     ! col global real gdp data (k US$/capita)
     real(r8) , pointer :: peatf_lf_col                (:)     ! col global peatland fraction data (0-1)
     integer  , pointer :: abm_lf_col                  (:)     ! col global peak month of crop fire emissions 

     real(r8) , pointer :: lgdp_col                    (:)     ! col gdp limitation factor for fire occurrence (0-1)
     real(r8) , pointer :: lgdp1_col                   (:)     ! col gdp limitation factor for fire spreading (0-1)
     real(r8) , pointer :: lpop_col                    (:)     ! col pop limitation factor for fire spreading (0-1)

     real(r8) , pointer :: fpi_vr_col                  (:,:)   ! col fraction of potential immobilization (no units) 
     real(r8) , pointer :: fpi_col                     (:)     ! col fraction of potential immobilization (no units) 
     real(r8),  pointer :: fpg_col                     (:)     ! col fraction of potential gpp (no units)

     !!! add phosphorus  -X. YANG

     integer  ,pointer  :: isoilorder                  (:)     ! col global soil order data
     real(r8) , pointer :: fpi_p_vr_col                (:,:)   ! col fraction of potential immobilization (no units) 
     real(r8) , pointer :: fpi_p_col                   (:)     ! col fraction of potential immobilization (no units) 
     real(r8),  pointer :: fpg_p_col                   (:)     ! col fraction of potential gpp (no units)

     real(r8) , pointer :: rf_decomp_cascade_col       (:,:,:) ! col respired fraction in decomposition step (frac)
     real(r8) , pointer :: pathfrac_decomp_cascade_col (:,:,:) ! col what fraction of C leaving a given pool passes through a given transition (frac) 
     real(r8) , pointer :: nfixation_prof_col          (:,:)   ! col (1/m) profile for N fixation additions 
     real(r8) , pointer :: ndep_prof_col               (:,:)   ! col (1/m) profile for N fixation additions
     real(r8) , pointer :: pdep_prof_col               (:,:)   ! col (1/m) profile for P deposition additions 
     real(r8) , pointer :: som_adv_coef_col            (:,:)   ! col SOM advective flux (m/s) 
     real(r8) , pointer :: som_diffus_coef_col         (:,:)   ! col SOM diffusivity due to bio/cryo-turbation (m2/s) 

     real(r8) , pointer :: annavg_t2m_col              (:)     ! col annual average of 2m air temperature, averaged from pft-level (K)
     real(r8) , pointer :: scalaravg_col               (:)     ! column average scalar for decompostion (for ad_spinup)
     real(r8) , pointer :: annsum_counter_col          (:)     ! col seconds since last annual accumulator turnover

     ! Fire
     real(r8) , pointer :: nfire_col                   (:)     ! col fire counts (count/km2/sec), valid only in Reg. C
     real(r8) , pointer :: fsr_col                     (:)     ! col fire spread rate at column level (m/s)
     real(r8) , pointer :: fd_col                      (:)     ! col fire duration at column level (hr)
     real(r8) , pointer :: lfc_col                     (:)     ! col conversion area fraction of BET and BDT that haven't burned before (/timestep)
     real(r8) , pointer :: lfc2_col                    (:)     ! col conversion area fraction of BET and BDT that burned (/sec)
     real(r8) , pointer :: dtrotr_col                  (:)     ! col annual decreased fraction coverage of BET on the gridcell (0-1)
     real(r8) , pointer :: trotr1_col                  (:)     ! col patch weight of BET and BDT on the gridcell(0-1)
     real(r8) , pointer :: trotr2_col                  (:)     ! col patch weight of BDT on the gridcell (0-1)
     real(r8) , pointer :: cropf_col                   (:)     ! col crop fraction in veg column (0-1)
     real(r8) , pointer :: baf_crop_col                (:)     ! col baf for cropland(/sec)
     real(r8) , pointer :: baf_peatf_col               (:)     ! col baf for peatland (/sec)
     real(r8) , pointer :: fbac_col                    (:)     ! col total burned area out of conversion (/sec)
     real(r8) , pointer :: fbac1_col                   (:)     ! col burned area out of conversion region due to land use fire (/sec)
     real(r8) , pointer :: wtlf_col                    (:)     ! col fractional coverage of non-crop Patches (0-1)
     real(r8) , pointer :: lfwt_col                    (:)     ! col fractional coverage of non-crop and non-bare-soil Patches (0-1)
     real(r8) , pointer :: farea_burned_col            (:)     ! col fractional area burned (/sec) 

     real(r8), pointer :: frootc_nfix_scalar_col       (:)     ! col scalar for nitrogen fixation
     real(r8), pointer :: decomp_litpool_rcn_col       (:,:,:) ! cn ratios of the decomposition pools
!!!!!!********DW**************
     integer           :: CropRestYear                         ! restart year from initial conditions file - increment as time elapses

     real(r8), pointer :: fpg_nh4_vr_col               (:,:)   ! fraction of plant nh4 demand that is satisfied (no units) BGC mode
     real(r8), pointer :: fpg_no3_vr_col               (:,:)   ! fraction of plant no3 demand that is satisfied (no units) BGC mode
     real(r8), pointer :: fpg_vr_col                   (:,:)   ! fraction of plant N demand that is satisfied (no units) CN mode
     real(r8), pointer :: fpg_p_vr_col                 (:,:)   ! fraction of plant p demand that is satisfied (no units) 


  contains
     procedure, public :: Init => init_col_cs
     procedure, public :: Clean => clean_col_cs
     procedure , public  :: SetValues => setvalues_col_cs
!     procedure , public  :: ZeroDWT => zerodwt_col_cs  !zerodwt only at patch level
     procedure , public  :: Restart => restart_col_cs
     procedure , public  :: Summary => summary_col_cs  !moved summary at patch level
     procedure , private :: InitAllocate => initallocate_col_cs
     procedure , private :: InitHistory => inithistory_col_cs
     procedure , private :: InitCold => initcold_col_cs

!!!!!!!!!!!************DW**************  CropRestIncYear Crop_prog should be at Veg level????
!        procedure, public  :: CropRestIncYear

  end type soilcol_carbon_state

  ! declare the public instances of soilcolumn types

  type(columnl_carbon_state) , public, target :: col_cs  ! column carbon state data structure

  !------------------------------------------------------------------------

  contains 


  !------------------------------------------------------------------------

  subroutine init_col_cs(this, bounds, carbon_type, ratio, c12_carbonstate_vars)
    class(column_carbon_state)               :: this
    type(bounds_type)      , intent(in)           :: bounds  
    character(len=3)       , intent(in)           :: carbon_type
    real(r8)               , intent(in)           :: ratio
    type(column_carbonstate_type) , intent(in), optional :: c12_carbonstate_vars

    call this%InitAllocate ( bounds)
    call this%InitHistory ( bounds, carbon_type)
    if (present(c12_carbonstate_vars)) then
       call this%InitCold  ( bounds, ratio, c12_carbonstate_vars)
    else
       call this%InitCold  ( bounds, ratio)
    end if

  end subroutine init_col_cs


  !------------------------------------------------------------------------
  subroutine initallocate_col_cs(this, bounds)
    !
    ! !ARGUMENTS:
    class (column_carbon_state) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
!    integer           :: begp,endp
    integer           :: begc,endc
    !------------------------------------------------------------------------

!    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%cwdc_col                 (begc :endc))                   ;     this%cwdc_col                 (:)   = nan
    allocate(this%ctrunc_col               (begc :endc))                   ;     this%ctrunc_col               (:)   = nan
    allocate(this%ctrunc_vr_col            (begc :endc,1:nlevdecomp_full)) ;     this%ctrunc_vr_col            (:,:) = nan
    allocate(this%seedc_col                (begc :endc))                   ;     this%seedc_col                (:)   = nan
    allocate(this%prod10c_col              (begc :endc))                   ;     this%prod10c_col              (:)   = nan
    allocate(this%prod100c_col             (begc :endc))                   ;     this%prod100c_col             (:)   = nan
    allocate(this%prod1c_col               (begc :endc))                   ;     this%prod1c_col               (:)   = nan
    allocate(this%totprodc_col             (begc :endc))                   ;     this%totprodc_col             (:)   = nan
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

    allocate(this%totabgc_col              (begc :endc))                   ;     this%totabgc_col              (:)   = nan
    allocate(this%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
    this%decomp_cpools_vr_col(:,:,:)= nan

    allocate(this%begcb_col   (begc:endc));     this%begcb_col   (:) = nan
    allocate(this%endcb_col   (begc:endc));     this%endcb_col   (:) = nan
    allocate(this%errcb_col   (begc:endc));     this%errcb_col   (:) = nan

    allocate(this%totpftc_beg_col(begc:endc));  this%totpftc_beg_col (:) = nan
    allocate(this%cwdc_beg_col   (begc:endc));  this%cwdc_beg_col    (:) = nan
    allocate(this%totlitc_beg_col(begc:endc));  this%totlitc_beg_col (:) = nan
    allocate(this%totsomc_beg_col(begc:endc));  this%totsomc_beg_col (:) = nan
    
    allocate(this%totpftc_end_col(begc:endc));  this%totpftc_end_col (:) = nan
    allocate(this%cwdc_end_col   (begc:endc));  this%cwdc_end_col    (:) = nan
    allocate(this%totlitc_end_col(begc:endc));  this%totlitc_end_col (:) = nan
    allocate(this%totsomc_end_col(begc:endc));  this%totsomc_end_col (:) = nan


! from CNStateType.F90

    allocate(this%gdp_lf_col          (begc:endc))                   ;
    allocate(this%peatf_lf_col        (begc:endc))                   ; 
    allocate(this%abm_lf_col          (begc:endc))                   ; 

    allocate(this%lgdp_col            (begc:endc))                   ; 
    allocate(this%lgdp1_col           (begc:endc))                   ; 
    allocate(this%lpop_col            (begc:endc))                   ;  

    allocate(this%fpi_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%fpi_vr_col          (:,:) = nan
    allocate(this%fpi_col             (begc:endc))                   ; this%fpi_col             (:)   = nan
    allocate(this%fpg_col             (begc:endc))                   ; this%fpg_col             (:)   = nan
    !!! add phosphours related variables
    allocate(this%isoilorder            (begc:endc))                   ; 
    allocate(this%fpi_p_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%fpi_p_vr_col          (:,:) = nan
    allocate(this%fpi_p_col             (begc:endc))                   ; this%fpi_p_col             (:)   = nan
    allocate(this%fpg_p_col             (begc:endc))                   ; this%fpg_p_col             (:)   = nan
    allocate(this%pdep_prof_col         (begc:endc,1:nlevdecomp_full)) ; this%pdep_prof_col       (:,:) = spval

    allocate(this%rf_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); 
    this%rf_decomp_cascade_col(:,:,:) = nan

    allocate(this%pathfrac_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));     
    this%pathfrac_decomp_cascade_col(:,:,:) = nan

    allocate(this%nfixation_prof_col  (begc:endc,1:nlevdecomp_full)) ; this%nfixation_prof_col  (:,:) = spval
    allocate(this%ndep_prof_col       (begc:endc,1:nlevdecomp_full)) ; this%ndep_prof_col       (:,:) = spval
    allocate(this%pdep_prof_col       (begc:endc,1:nlevdecomp_full)) ; this%pdep_prof_col       (:,:) = spval
    allocate(this%som_adv_coef_col    (begc:endc,1:nlevdecomp_full)) ; this%som_adv_coef_col    (:,:) = spval
    allocate(this%som_diffus_coef_col (begc:endc,1:nlevdecomp_full)) ; this%som_diffus_coef_col (:,:) = spval

    allocate(this%annsum_counter_col  (begc:endc))                   ; this%annsum_counter_col  (:)   = nan
    allocate(this%annavg_t2m_col      (begc:endc))                   ; this%annavg_t2m_col      (:)   = nan
    allocate(this%scalaravg_col       (begc:endc))                   ; this%scalaravg_col       (:)   = nan

    allocate(this%nfire_col           (begc:endc))                   ; this%nfire_col           (:)   = spval
    allocate(this%fsr_col             (begc:endc))                   ; this%fsr_col             (:)   = nan
    allocate(this%fd_col              (begc:endc))                   ; this%fd_col              (:)   = nan
    allocate(this%lfc_col             (begc:endc))                   ; this%lfc_col             (:)   = spval
    allocate(this%lfc2_col            (begc:endc))                   ; this%lfc2_col            (:)   = 0._r8
    allocate(this%dtrotr_col          (begc:endc))                   ; this%dtrotr_col          (:)   = 0._r8
    allocate(this%trotr1_col          (begc:endc))                   ; this%trotr1_col          (:)   = 0._r8
    allocate(this%trotr2_col          (begc:endc))                   ; this%trotr2_col          (:)   = 0._r8
    allocate(this%cropf_col           (begc:endc))                   ; this%cropf_col           (:)   = nan
    allocate(this%baf_crop_col        (begc:endc))                   ; this%baf_crop_col        (:)   = nan
    allocate(this%baf_peatf_col       (begc:endc))                   ; this%baf_peatf_col       (:)   = nan
    allocate(this%fbac_col            (begc:endc))                   ; this%fbac_col            (:)   = nan
    allocate(this%fbac1_col           (begc:endc))                   ; this%fbac1_col           (:)   = nan
    allocate(this%wtlf_col            (begc:endc))                   ; this%wtlf_col            (:)   = nan
    allocate(this%lfwt_col            (begc:endc))                   ; this%lfwt_col            (:)   = nan
    allocate(this%farea_burned_col    (begc:endc))                   ; this%farea_burned_col    (:)   = nan
    allocate(this%decomp_litpool_rcn_col (begc:endc, 1:nlevdecomp_full, 4)); this%decomp_litpool_rcn_col (:,:,:) = nan
    allocate(this%frootc_nfix_scalar_col (begc:endc))                ; this%frootc_nfix_scalar_col(:) = nan
 
!!!!!!!************DW****************
    this%CropRestYear = 0

    allocate(this%fpg_nh4_vr_col              (begc:endc,1:nlevdecomp_full)) ; this%fpg_nh4_vr_col(:,:) = nan 
    allocate(this%fpg_no3_vr_col              (begc:endc,1:nlevdecomp_full)) ; this%fpg_no3_vr_col(:,:) = nan
    allocate(this%fpg_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%fpg_vr_col    (:,:) = nan
    allocate(this%fpg_p_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%fpg_p_vr_col  (:,:) = nan


  end subroutine initallocate_col_cs


  !-----------------------------------------------------------------------
!   subroutine SetValues ( this, &
!       num_patch, filter_patch, value_patch, &
!       num_column, filter_column, value_column)

   subroutine setvalues_col_cs ( this, &
       num_column, filter_column, value_column)

    !
    ! !DESCRIPTION:
    ! Set carbon state variables at soil column level
    !
    ! !ARGUMENTS:
    class (column_carbon_state) :: this
!    integer , intent(in) :: num_patch
!    integer , intent(in) :: filter_patch(:)
!    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

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

  end subroutine setvalues_col_cs

  !-----------------------------------------------------------------------
 !  subroutine summary_col_cs(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
  subroutine summary_col_cs(this, bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform column-level carbon summary calculations
    !
    ! !USES:
    use clm_varctl       , only: iulog
    use clm_time_manager , only: get_step_size
    use clm_varcon       , only: secspday
    use clm_varpar       , only: nlevdecomp, ndecomp_pools 
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
!    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
!    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    real(r8) :: maxdepth        ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    ! column level summary

      ! vertically integrate each of the decomposing C pools
      do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_cpools_col(c,l) = 0._r8
       end do
      end do
      do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
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
    do j = 1, nlevdecomp
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
            this%totprodc_col(c) + &
            this%seedc_col(c)    + &
            this%ctrunc_col(c)
            
       this%totabgc_col(c) = &
            this%totpftc_col(c)  + &
            this%totprodc_col(c) + &
            this%seedc_col(c)    + &
            this%ctrunc_col(c)    
    end do

  end subroutine summary_col_cs


  subroutine inithistory_col_cs(this, bounds, carbon_type)
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
    class (column_carbon_state) :: this
    type(bounds_type)         , intent(in) :: bounds 
    character(len=3)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
!    integer           :: begp,endp
    integer           :: begc,endc
!    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays

!   from CNStateType.90
    character(8)      :: vr_suffix
    character(10)     :: active
!    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays


    !---------------------------------------------------------------------

!    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
 !   begg = bounds%begg; endg = bounds%endg
 

    !-------------------------------
    ! C state variables - column
    !-------------------------------

    ! add history fields for all CLAMP CN variables

    if (carbon_type == 'c12') then


         !those variables are now ouput in betr
         this%decomp_cpools_col(begc:endc,:) = spval
         do l  = 1, ndecomp_pools
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
            avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=this%totecosysc_col)

       this%totcolc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
            avgflag='A', long_name='total column carbon, incl veg and cpool', &
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
            avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=this%totecosysc_col)

       this%totcolc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total column carbon, incl veg and cpool', &
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
            avgflag='A', long_name='C14 total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=this%totecosysc_col)

       this%totcolc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTCOLC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total column carbon, incl veg and cpool', &
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
    endif


!  from CNStateType.F90


    this%nfixation_prof_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='NFIXATION_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for biological N fixation', &
         ptr_col=this%nfixation_prof_col, default='inactive')

    this%ndep_prof_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='NDEP_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for atmospheric N  deposition', &
         ptr_col=this%ndep_prof_col, default='inactive')

    this%pdep_prof_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='PDEP_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for atmospheric P  deposition', &
         ptr_col=this%pdep_prof_col, default='inactive')

    this%som_adv_coef_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SOM_ADV_COEF', units='m/s',  type2d='levdcmp', &
         avgflag='A', long_name='advection term for vertical SOM translocation', &
         ptr_col=this%som_adv_coef_col, default='inactive')

    this%som_diffus_coef_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SOM_DIFFUS_COEF', units='m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='diffusion coefficient for vertical SOM translocation', &
         ptr_col=this%som_diffus_coef_col, default='inactive')

    this%lfc2_col(begc:endc) = spval
    call hist_addfld1d (fname='LFC2', units='per sec', &
         avgflag='A', long_name='conversion area fraction of BET and BDT that burned', &
         ptr_col=this%lfc2_col)

    if ( nlevdecomp_full > 1 ) then
       this%fpi_col(begc:endc) = spval
       call hist_addfld1d (fname='FPI', units='proportion', &
            avgflag='A', long_name='fraction of potential immobilization of nitrogen', &
            ptr_col=this%fpi_col)
       this%fpi_p_col(begc:endc) = spval
       call hist_addfld1d (fname='FPI_P', units='proportion', &
            avgflag='A', long_name='fraction of potential immobilization of phosphorus', &
            ptr_col=this%fpi_p_col)
    endif

    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif
    this%fpi_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='FPI'//trim(vr_suffix), units='proportion', type2d='levdcmp', & 
         avgflag='A', long_name='fraction of potential immobilization of nitrogen', &
         ptr_col=this%fpi_vr_col)
    this%fpi_p_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='FPI_P'//trim(vr_suffix), units='proportion', type2d='levdcmp', & 
         avgflag='A', long_name='fraction of potential immobilization of phosphorus', &
         ptr_col=this%fpi_p_vr_col)

    this%fpg_col(begc:endc) = spval
    call hist_addfld1d (fname='FPG', units='proportion', &
         avgflag='A', long_name='fraction of potential gpp due to N limitation', &
         ptr_col=this%fpg_col)

    this%fpg_p_col(begc:endc) = spval
    call hist_addfld1d (fname='FPG_P', units='proportion', &
         avgflag='A', long_name='fraction of potential gpp due to P limitation', &
         ptr_col=this%fpg_p_col)

    this%annsum_counter_col(begc:endc) = spval
    call hist_addfld1d (fname='ANNSUM_COUNTER', units='s', &
         avgflag='A', long_name='seconds since last annual accumulator turnover', &
         ptr_col=this%annsum_counter_col, default='inactive')

    this%annavg_t2m_col(begc:endc) = spval
    call hist_addfld1d (fname='CANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average of 2m air temperature', &
         ptr_col=this%annavg_t2m_col, default='inactive')

    this%scalaravg_col(begc:endc) = spval
    call hist_addfld1d(fname='SCALARAVG', units='fraction', &
         avgflag='A', long_name='average of decomposition scalar', &
         ptr_col=this%scalaravg_col)

    this%nfire_col(begc:endc) = spval
    call hist_addfld1d (fname='NFIRE',  units='counts/km2/sec', &
         avgflag='A', long_name='fire counts valid only in Reg.C', &
         ptr_col=this%nfire_col)

    this%farea_burned_col(begc:endc) = spval
    call hist_addfld1d (fname='FAREA_BURNED',  units='proportion', &
         avgflag='A', long_name='timestep fractional area burned', &
         ptr_col=this%farea_burned_col)

    this%baf_crop_col(begc:endc) = spval
    call hist_addfld1d (fname='BAF_CROP',  units='proportion/sec', &
         avgflag='A', long_name='fractional area burned for crop', &
         ptr_col=this%baf_crop_col)

    this%baf_peatf_col(begc:endc) = spval
    call hist_addfld1d (fname='BAF_PEATF',  units='proportion/sec', &
         avgflag='A', long_name='fractional area burned in peatland', &
         ptr_col=this%baf_peatf_col)
 
  end subroutine inithistory_col_cs

 !-----------------------------------------------------------------------
  subroutine initcold_col_cs(this, bounds, ratio, c12_carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon , only: istsoil
    use pftvarcon       , only: noveg, npcropmin
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this 
    type(bounds_type), intent(in) :: bounds  
    real(r8), intent(in) :: ratio
    type(column_carbon_state), optional, intent(in) :: c12_carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,j,k
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
!    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
!    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches


!    Moved from CNStateType.F90

    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use clm_varctl , only : nsrest, nsrStartup
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(cnstate_type) :: this
    type(bounds_type), intent(in) :: bounds   
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p,n,j,m            ! indices
    real(r8) ,pointer     :: gdp (:)                  ! global gdp data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: peatf (:)                ! global peatf data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: soilorder_rdin (:)       ! global soil order data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: abm (:)                  ! global abm data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: gti (:)                  ! read in - fmax (needs to be a pointer for use in ncdio)
    integer               :: dimid                    ! dimension id
    integer               :: ier                      ! error status
    type(file_desc_t)     :: ncid                     ! netcdf id
    logical               :: readvar 
    character(len=256)    :: locfn                    ! local filename
    integer               :: begc, endc
    integer               :: begg, endg

    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg


    !   End of section from CNStateType.F90

    !-----------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do


    ! initialize column-level variables
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

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

    ! initialize fields for special filters

    call this%SetValues (&
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

!!!!!   Moved from CNStateType.F90

        ! --------------------------------------------------------------------
    ! Open surface dataset
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil color, sand and clay boundary data .....'
    end if

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_inqdlen(ncid,dimid,nlevsoifl,name='nlevsoi')
    if ( .not. more_vertlayers )then
       if ( nlevsoifl /= nlevsoi )then
          call endrun(msg=' ERROR: Number of soil layers on file does NOT match the number being used'//&
               errMsg(__FILE__, __LINE__))
       end if
    else
       ! read in layers, interpolate to high resolution grid later
    end if


    ! --------------------------------------------------------------------
    ! Read in GDP data 
    ! --------------------------------------------------------------------

    allocate(gdp(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='gdp', flag='read', data=gdp, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: gdp NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%gdp_lf_col(c) = gdp(g)
    end do
    deallocate(gdp)

    ! --------------------------------------------------------------------
    ! Read in peatf data 
    ! --------------------------------------------------------------------

    allocate(peatf(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='peatf', flag='read', data=peatf, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: peatf NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%peatf_lf_col(c) = peatf(g)
    end do
    deallocate(peatf)

    ! --------------------------------------------------------------------
    ! Read in soilorder data 
    ! --------------------------------------------------------------------

    allocate(soilorder_rdin(bounds%begg:bounds%endg))
    !call ncd_io(ncid=ncid, varname='SOIL_ORDER', flag='read',data=soilorder_rdin, dim1name=grlnd, readvar=readvar)
    !if (.not. readvar) then
    !   call endrun(msg=' ERROR: SOIL_ORDER NOT on surfdata file'//errMsg(__FILE__, __LINE__))
    !end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
!       this%isoilorder(c) = soilorder_rdin(g)
       this%isoilorder(c) = 12
    end do
    deallocate(soilorder_rdin)


    ! --------------------------------------------------------------------
    ! Read in ABM data 
    ! --------------------------------------------------------------------

    allocate(abm(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='abm', flag='read', data=abm, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: abm NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%abm_lf_col(c) = abm(g)
    end do
    deallocate(abm)

    ! Close file

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) 'Successfully read fmax, soil color, sand and clay boundary data'
       write(iulog,*)
    endif
    
    ! --------------------------------------------------------------------
    ! Initialize terms needed for dust model
    ! TODO - move these terms to DUSTMod module variables 
    ! --------------------------------------------------------------------
       
    do c = bounds%begc, bounds%endc

!!!  ******DW*****
!!!!  need to deal with col and lun here

       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          this%annsum_counter_col (c) = spval
          this%annavg_t2m_col     (c) = spval
          this%scalaravg_col      (c) = spval
          this%nfire_col          (c) = spval
          this%baf_crop_col       (c) = spval
          this%baf_peatf_col      (c) = spval
          this%fbac_col           (c) = spval
          this%fbac1_col          (c) = spval
          this%farea_burned_col   (c) = spval
          this%fpi_col            (c) = spval
          this%fpg_col            (c) = spval
          this%fpi_p_col          (c) = spval
          this%fpg_p_col          (c) = spval
          do j = 1,nlevdecomp_full
             this%fpi_vr_col(c,j) = spval
             this%fpi_p_vr_col(c,j) = spval
          end do
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%annsum_counter_col(c) = 0._r8   
          this%annavg_t2m_col(c)     = 280._r8 

          ! fire related variables 
          this%baf_crop_col(c)       = 0._r8 
          this%baf_peatf_col(c)      = 0._r8 
          this%fbac_col(c)           = 0._r8 
          this%fbac1_col(c)          = 0._r8 
          this%farea_burned_col(c)   = 0._r8 

          if (nsrest == nsrStartup) this%nfire_col(c) = 0._r8

          ! initialize fpi_vr so that levels below nlevsoi are not nans
          this%fpi_vr_col(c,1:nlevdecomp_full)          = 0._r8 
          this%fpi_p_vr_col(c,1:nlevdecomp_full)          = 0._r8 
          this%som_adv_coef_col(c,1:nlevdecomp_full)    = 0._r8 
          this%som_diffus_coef_col(c,1:nlevdecomp_full) = 0._r8 

          ! initialize the profiles for converting to vertically resolved carbon pools
          this%nfixation_prof_col(c,1:nlevdecomp_full)  = 0._r8 
          this%ndep_prof_col(c,1:nlevdecomp_full)       = 0._r8 
          this%pdep_prof_col(c,1:nlevdecomp_full)       = 0._r8 
       end if
    end do

    ! ecophysiological variables
    ! phenology variables
       
    ! ecophysiological variables


    ! fire variables

    do c = bounds%begc,bounds%endc
       this%lfc2_col(c) = 0._r8
    end do

  end subroutine initcold_col_cs

subroutine restart_col_cs ( this,  bounds, ncid, flag, carbon_type, c12_carbonstate_vars, cnstate_vars)
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

!   moved from CNStateType.F90

    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun


    !
    ! !ARGUMENTS:
    class (column_carbon_state) :: this
    type(bounds_type)         , intent(in)           :: bounds 
    type(file_desc_t)         , intent(inout)        :: ncid   ! netcdf id
    character(len=*)          , intent(in)           :: flag   !'read' or 'write'
    character(len=3)          , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    type (column_carbon_state)   , intent(in), optional :: c12_carbonstate_vars 
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


!  moved from CNStateType.F90

    integer, pointer :: temp1d(:) ! temporary
    integer          :: p,j,c,i   ! indices


    !------------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_carbonstate_vars)) then
          call endrun(msg=' ERROR: for C14 must pass in c12_carbonstate_vars as argument' //&
               errMsg(__FILE__, __LINE__))
       end if
    end if
 
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

    if (carbon_type == 'c12') then
        if (flag == 'write') then
           idata = spinup_state
        end if
        call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
             long_name='Spinup state of the model that wrote this restart file: ' &
             // ' 0 = normal model mode, 1 = AD spinup', units='', &
             interpinic_flag='copy', readvar=readvar,  data=idata)
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

!!!!!!!!****************DW****************the if statement inside IF/THEN/ELSE statement
           do k = 1, ndecomp_pools
            do c = bounds%begc, bounds%endc
              do j = 1, nlevdecomp

	    	       if ( exit_spinup ) then
    		        m = decomp_cascade_con%spinup_factor(k)
                if (decomp_cascade_con%spinup_factor(k) > 1)   m = m / cnstate_vars%scalaravg_col(c)
               else if ( enter_spinup ) then 
		            m = 1. / decomp_cascade_con%spinup_factor(k)
		            if (decomp_cascade_con%spinup_factor(k) > 1) m = m * cnstate_vars%scalaravg_col(c)
		           end if
               this%decomp_cpools_vr_col(c,j,k) = this%decomp_cpools_vr_col(c,j,k) * m
 
              end do
            end do
           end do

          end if
    end if

!!!!  moved from CNStateType.F90  crop_prog is at PFT level???

    if (use_vertsoilc) then
       ptr2d => this%fpi_vr_col
       call restartvar(ncid=ncid, flag=flag, varname='fpi_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='fraction of potential immobilization of nitrogen',  units='unitless', &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
       ptr2d => this%fpi_p_vr_col
       call restartvar(ncid=ncid, flag=flag, varname='fpi_p_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='fraction of potential immobilization of phosphorus',  units='unitless', &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%fpi_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
       call restartvar(ncid=ncid, flag=flag, varname='fpi', xtype=ncd_double,  &
            dim1name='column', &
            long_name='fraction of potential immobilization of nitrogen',  units='unitless', &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       ptr1d => this%fpi_p_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
       call restartvar(ncid=ncid, flag=flag, varname='fpi_p', xtype=ncd_double,  &
            dim1name='column', &
            long_name='fraction of potential immobilization of phosphorus',  units='unitless', &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    if (use_vertsoilc) then
       ptr2d => this%som_adv_coef_col
       call restartvar(ncid=ncid, flag=flag, varname='som_adv_coef_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='SOM advective flux', units='m/s', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    end if
    
    if (use_vertsoilc) then
       ptr2d => this%som_diffus_coef_col
       call restartvar(ncid=ncid, flag=flag, varname='som_diffus_coef_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='SOM diffusivity due to bio/cryo-turbation',  units='m^2/s', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='fpg', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpg_col) 

    call restartvar(ncid=ncid, flag=flag, varname='fpg_p', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpg_p_col) 

    call restartvar(ncid=ncid, flag=flag, varname='annsum_counter', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_counter_col) 

    call restartvar(ncid=ncid, flag=flag, varname='burndate', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%burndate_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='lfc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lfc_col) 

    call restartvar(ncid=ncid, flag=flag, varname='cannavg_t2m', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annavg_t2m_col) 

    call restartvar(ncid=ncid, flag=flag, varname='scalaravg_col', xtype=ncd_double, &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag = 'interp', readvar=readvar, data=this%scalaravg_col)


  end subroutine restart_col_cs


end module ColumnCarbonStateType

