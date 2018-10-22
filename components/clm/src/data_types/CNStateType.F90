module CNStateType

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
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  use VegetationType      , only : veg_pp                
  use clm_varctl     , only: forest_fert_exp
  use clm_varctl          , only : nu_com
  use clm_varctl   , only:  use_fates,use_crop

  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PRIVATE MEMBER FUNCTIONS: 
  private :: checkDates
  ! !PUBLIC TYPES:
  integer    , pointer, public :: fert_type         (:)
  integer    , pointer, public :: fert_continue     (:)
  real(r8)   , pointer, public :: fert_dose         (:,:)
  integer    , pointer, public :: fert_start        (:)
  integer    , pointer, public :: fert_end          (:)
  ! 
  ! !PUBLIC TYPES:
  type, public :: cnstate_type

     integer  , pointer :: burndate_patch              (:)     ! patch crop burn date
     real(r8) , pointer :: lfpftd_patch                (:)     ! patch decrease of pft weight (0-1) on the column for the timestep 

     ! Prognostic crop model -  Note that cropplant and harvdate could be 2D to facilitate rotation
     real(r8) , pointer :: hdidx_patch                 (:)     ! patch cold hardening index?
     real(r8) , pointer :: cumvd_patch                 (:)     ! patch cumulative vernalization d?ependence?
     real(r8) , pointer :: gddmaturity_patch           (:)     ! patch growing degree days (gdd) needed to harvest (ddays)
     real(r8) , pointer :: huileaf_patch               (:)     ! patch heat unit index needed from planting to leaf emergence
     real(r8) , pointer :: huigrain_patch              (:)     ! patch heat unit index needed to reach vegetative maturity
     real(r8) , pointer :: aleafi_patch                (:)     ! patch saved leaf allocation coefficient from phase 2
     real(r8) , pointer :: astemi_patch                (:)     ! patch saved stem allocation coefficient from phase 2
     real(r8) , pointer :: aleaf_patch                 (:)     ! patch leaf allocation coefficient
     real(r8) , pointer :: astem_patch                 (:)     ! patch stem allocation coefficient
     real(r8) , pointer :: htmx_patch                  (:)     ! patch max hgt attained by a crop during yr (m)
     integer  , pointer :: peaklai_patch               (:)     ! patch 1: max allowed lai; 0: not at max

     integer  , pointer :: idop_patch                  (:)     ! patch date of planting
     real(r8) , pointer :: leaf_prof_patch             (:,:)   ! patch (1/m) profile of leaves (vertical profiles for calculating fluxes)
     real(r8) , pointer :: froot_prof_patch            (:,:)   ! patch (1/m) profile of fine roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: croot_prof_patch            (:,:)   ! patch (1/m) profile of coarse roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: stem_prof_patch             (:,:)   ! patch (1/m) profile of stems (vertical profiles for calculating fluxes)

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

     real(r8) , pointer :: tempavg_t2m_patch           (:)     ! patch temporary average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_patch            (:)     ! patch annual average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_col              (:)     ! col annual average of 2m air temperature, averaged from pft-level (K)
     real(r8) , pointer :: scalaravg_col               (:,:)   ! column average scalar for decompostion (for ad_spinup)
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

     real(r8), pointer :: dormant_flag_patch           (:)     ! patch dormancy flag
     real(r8), pointer :: days_active_patch            (:)     ! patch number of days since last dormancy
     real(r8), pointer :: onset_flag_patch             (:)     ! patch onset flag
     real(r8), pointer :: onset_counter_patch          (:)     ! patch onset days counter
     real(r8), pointer :: onset_gddflag_patch          (:)     ! patch onset flag for growing degree day sum
     real(r8), pointer :: onset_fdd_patch              (:)     ! patch onset freezing degree days counter
     real(r8), pointer :: onset_gdd_patch              (:)     ! patch onset growing degree days
     real(r8), pointer :: onset_swi_patch              (:)     ! patch onset soil water index
     real(r8), pointer :: offset_flag_patch            (:)     ! patch offset flag
     real(r8), pointer :: offset_counter_patch         (:)     ! patch offset days counter
     real(r8), pointer :: offset_fdd_patch             (:)     ! patch offset freezing degree days counter
     real(r8), pointer :: offset_swi_patch             (:)     ! patch offset soil water index
     real(r8), pointer :: grain_flag_patch             (:)     ! patch 1: grain fill stage; 0: not
     real(r8), pointer :: lgsf_patch                   (:)     ! patch long growing season factor [0-1]
     real(r8), pointer :: bglfr_patch                  (:)     ! patch background litterfall rate (1/s)
     real(r8), pointer :: bglfr_leaf_patch             (:)     ! patch background leaf litterfall rate (1/s)
     real(r8), pointer :: bglfr_froot_patch            (:)     ! patch background fine root litterfall rate (1/s)
     real(r8), pointer :: bgtr_patch                   (:)     ! patch background transfer growth rate (1/s)
     real(r8), pointer :: alloc_pnow_patch             (:)     ! patch fraction of current allocation to display as new growth (DIM)
     real(r8), pointer :: c_allometry_patch            (:)     ! patch C allocation index (DIM)
     real(r8), pointer :: n_allometry_patch            (:)     ! patch N allocation index (DIM)
     real(r8), pointer :: tempsum_potential_gpp_patch  (:)     ! patch temporary annual sum of potential GPP
     real(r8), pointer :: annsum_potential_gpp_patch   (:)     ! patch annual sum of potential GPP
     real(r8), pointer :: tempmax_retransn_patch       (:)     ! patch temporary annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: annmax_retransn_patch        (:)     ! patch annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: downreg_patch                (:)     ! patch fractional reduction in GPP due to N limitation (DIM)
     real(r8), pointer :: rc14_atm_patch               (:)     ! patch C14O2/C12O2 in atmosphere


     !!! add phosphorus  -X. YANG
     real(r8), pointer :: p_allometry_patch            (:)     ! patch P allocation index (DIM)
     real(r8), pointer :: tempmax_retransp_patch       (:)     ! patch temporary annual max of retranslocated P pool (gP/m2)
     real(r8), pointer :: annmax_retransp_patch        (:)     ! patch annual max of retranslocated P pool (gP/m2)

     real(r8), pointer :: frootc_nfix_scalar_col       (:)     ! col scalar for nitrogen fixation
     real(r8), pointer :: decomp_litpool_rcn_col       (:,:,:) ! cn ratios of the decomposition pools

     integer           :: CropRestYear                         ! restart year from initial conditions file - increment as time elapses

     real(r8), pointer :: fpg_nh4_vr_col               (:,:)   ! fraction of plant nh4 demand that is satisfied (no units) BGC mode
     real(r8), pointer :: fpg_no3_vr_col               (:,:)   ! fraction of plant no3 demand that is satisfied (no units) BGC mode
     real(r8), pointer :: fpg_vr_col                   (:,:)   ! fraction of plant N demand that is satisfied (no units) CN mode
     real(r8), pointer :: fpg_p_vr_col                 (:,:)   ! fraction of plant p demand that is satisfied (no units) 
     real(r8), pointer :: cn_scalar                    (:)     ! cn scaling factor for root n uptake kinetics (no units) 
     real(r8), pointer :: cp_scalar                    (:)     ! cp scaling factor for root p uptake kinetics (no units)
     real(r8), pointer :: np_scalar                    (:)     ! np scaling factor for root n/p uptake kinetics (no units)
     real(r8), pointer :: cost_ben_scalar              (:)     ! cost benefit analysis scaling factor for root n uptake kinetics (no units)
     real(r8), pointer :: cn_scalar_runmean            (:)     ! long term average of cn scaling factor for root n uptake kinetics (no units) 
     real(r8), pointer :: cp_scalar_runmean            (:)     ! long term average of cp scaling factor for root p uptake kinetics (no units)

     real(r8), pointer :: frac_loss_lit_to_fire_col        (:)
     real(r8), pointer :: frac_loss_cwd_to_fire_col        (:)

     ! soil phosphorus pools Qing Z. 2017
     real(r8), pointer :: labp_col             (:)   ! labile phosphorus g/m2
     real(r8), pointer :: secp_col             (:)   ! secondary phosphorus g/m2
     real(r8), pointer :: occp_col             (:)   ! occluded phosphorus g/m2
     real(r8), pointer :: prip_col             (:)   ! parent material phosphorus g/m2
     logical           :: pdatasets_present          ! surface dataset has p pools info
     !!! annual mortality rate dynamically calcaulted at patch
     real(r8), pointer :: r_mort_cal_patch                 (:)     ! patch annual mortality rate  

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, public  :: CropRestIncYear
     procedure, private :: InitAllocate
     procedure, private :: InitHistory  
     procedure, private :: InitCold     
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars

  end type cnstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(cnstate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate ( bounds )
    if (use_cn) then
       call this%InitHistory ( bounds )
    end if
    call this%InitCold ( bounds ) 

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(cnstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%burndate_patch      (begp:endp))                   ; this%burndate_patch      (:)   = ispval
    allocate(this%lfpftd_patch        (begp:endp))                   ;

    allocate(this%hdidx_patch         (begp:endp))                   ; this%hdidx_patch         (:)   = nan
    allocate(this%cumvd_patch         (begp:endp))                   ; this%cumvd_patch         (:)   = nan
    allocate(this%gddmaturity_patch   (begp:endp))                   ; this%gddmaturity_patch   (:)   = spval
    allocate(this%huileaf_patch       (begp:endp))                   ; this%huileaf_patch       (:)   = nan
    allocate(this%huigrain_patch      (begp:endp))                   ; this%huigrain_patch      (:)   = 0.0_r8
    allocate(this%aleafi_patch        (begp:endp))                   ; this%aleafi_patch        (:)   = nan
    allocate(this%astemi_patch        (begp:endp))                   ; this%astemi_patch        (:)   = nan
    allocate(this%aleaf_patch         (begp:endp))                   ; this%aleaf_patch         (:)   = nan
    allocate(this%astem_patch         (begp:endp))                   ; this%astem_patch         (:)   = nan
    allocate(this%htmx_patch          (begp:endp))                   ; this%htmx_patch          (:)   = 0.0_r8
    allocate(this%peaklai_patch       (begp:endp))                   ; this%peaklai_patch       (:)   = 0

    allocate(this%idop_patch          (begp:endp))                   ; this%idop_patch          (:)   = huge(1)
    allocate(this%leaf_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%leaf_prof_patch     (:,:) = spval
    allocate(this%froot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%froot_prof_patch    (:,:) = spval
    allocate(this%croot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%croot_prof_patch    (:,:) = spval
    allocate(this%stem_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%stem_prof_patch     (:,:) = spval

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

    allocate(this%tempavg_t2m_patch   (begp:endp))                   ; this%tempavg_t2m_patch   (:)   = nan
    allocate(this%annsum_counter_col  (begc:endc))                   ; this%annsum_counter_col  (:)   = nan
    allocate(this%annavg_t2m_col      (begc:endc))                   ; this%annavg_t2m_col      (:)   = nan
    allocate(this%scalaravg_col       (begc:endc,1:nlevdecomp_full)) ; this%scalaravg_col       (:,:) = spval
    allocate(this%annavg_t2m_patch    (begp:endp))                   ; this%annavg_t2m_patch    (:)   = nan

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
    this%CropRestYear = 0

    allocate(this%dormant_flag_patch          (begp:endp)) ;    this%dormant_flag_patch          (:) = nan
    allocate(this%days_active_patch           (begp:endp)) ;    this%days_active_patch           (:) = nan
    allocate(this%onset_flag_patch            (begp:endp)) ;    this%onset_flag_patch            (:) = nan
    allocate(this%onset_counter_patch         (begp:endp)) ;    this%onset_counter_patch         (:) = nan
    allocate(this%onset_gddflag_patch         (begp:endp)) ;    this%onset_gddflag_patch         (:) = nan
    allocate(this%onset_fdd_patch             (begp:endp)) ;    this%onset_fdd_patch             (:) = nan
    allocate(this%onset_gdd_patch             (begp:endp)) ;    this%onset_gdd_patch             (:) = nan
    allocate(this%onset_swi_patch             (begp:endp)) ;    this%onset_swi_patch             (:) = nan
    allocate(this%offset_flag_patch           (begp:endp)) ;    this%offset_flag_patch           (:) = nan
    allocate(this%offset_counter_patch        (begp:endp)) ;    this%offset_counter_patch        (:) = nan
    allocate(this%offset_fdd_patch            (begp:endp)) ;    this%offset_fdd_patch            (:) = nan
    allocate(this%offset_swi_patch            (begp:endp)) ;    this%offset_swi_patch            (:) = nan
    allocate(this%grain_flag_patch            (begp:endp)) ;    this%grain_flag_patch            (:) = nan
    allocate(this%lgsf_patch                  (begp:endp)) ;    this%lgsf_patch                  (:) = nan
    allocate(this%bglfr_patch                 (begp:endp)) ;    this%bglfr_patch                 (:) = nan
    allocate(this%bglfr_leaf_patch            (begp:endp)) ;    this%bglfr_leaf_patch            (:) = nan
    allocate(this%bglfr_froot_patch           (begp:endp)) ;    this%bglfr_froot_patch           (:) = nan
    allocate(this%bgtr_patch                  (begp:endp)) ;    this%bgtr_patch                  (:) = nan
    allocate(this%alloc_pnow_patch            (begp:endp)) ;    this%alloc_pnow_patch            (:) = nan
    allocate(this%c_allometry_patch           (begp:endp)) ;    this%c_allometry_patch           (:) = nan
    allocate(this%n_allometry_patch           (begp:endp)) ;    this%n_allometry_patch           (:) = nan
    allocate(this%tempsum_potential_gpp_patch (begp:endp)) ;    this%tempsum_potential_gpp_patch (:) = nan
    allocate(this%annsum_potential_gpp_patch  (begp:endp)) ;    this%annsum_potential_gpp_patch  (:) = nan
    allocate(this%tempmax_retransn_patch      (begp:endp)) ;    this%tempmax_retransn_patch      (:) = nan
    allocate(this%annmax_retransn_patch       (begp:endp)) ;    this%annmax_retransn_patch       (:) = nan
    allocate(this%downreg_patch               (begp:endp)) ;    this%downreg_patch               (:) = nan
    allocate(this%rc14_atm_patch              (begp:endp)) ;    this%rc14_atm_patch              (:) = nan    


    !! add phosphorus -X.YANG
    allocate(this%p_allometry_patch           (begp:endp)) ;    this%p_allometry_patch           (:) = nan
    allocate(this%tempmax_retransp_patch      (begp:endp)) ;    this%tempmax_retransp_patch      (:) = nan
    allocate(this%annmax_retransp_patch       (begp:endp)) ;    this%annmax_retransp_patch       (:) = nan

    allocate(this%fpg_nh4_vr_col              (begc:endc,1:nlevdecomp_full)) ; this%fpg_nh4_vr_col(:,:) = nan 
    allocate(this%fpg_no3_vr_col              (begc:endc,1:nlevdecomp_full)) ; this%fpg_no3_vr_col(:,:) = nan
    allocate(this%fpg_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%fpg_vr_col    (:,:) = nan
    allocate(this%fpg_p_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%fpg_p_vr_col  (:,:) = nan
    allocate(this%cn_scalar                   (begp:endp))                   ; this%cn_scalar     (:) = 0.0
    allocate(this%cp_scalar                   (begp:endp))                   ; this%cp_scalar     (:) = 0.0
    allocate(this%np_scalar                   (begp:endp))                   ; this%np_scalar     (:) = 0.0
    allocate(this%cost_ben_scalar             (begp:endp))                   ; this%cost_ben_scalar(:) = 0.0
    allocate(this%cn_scalar_runmean           (begp:endp))                   ; this%cn_scalar_runmean (:) = 0.0
    allocate(this%cp_scalar_runmean           (begp:endp))                   ; this%cp_scalar_runmean (:) = 0.0
    allocate(this%frac_loss_lit_to_fire_col       (begc:endc))               ; this%frac_loss_lit_to_fire_col(:) =0._r8
    allocate(this%frac_loss_cwd_to_fire_col       (begc:endc))               ; this%frac_loss_cwd_to_fire_col(:) =0._r8
    allocate(fert_type                        (begc:endc))                   ; fert_type     (:) = 0
    allocate(fert_continue                    (begc:endc))                   ; fert_continue (:) =0
    allocate(fert_dose                        (begc:endc,1:12))              ; fert_dose     (:,:) = nan

    ! soil phosphorus pools Qing Z. 2017
    allocate(this%labp_col                    (begc:endc))                   ; this%labp_col(:) = nan
    allocate(this%secp_col                    (begc:endc))                   ; this%secp_col(:) = nan
    allocate(this%occp_col                    (begc:endc))                   ; this%occp_col(:) = nan
    allocate(this%prip_col                    (begc:endc))                   ; this%prip_col(:) = nan
    
    allocate(fert_start                       (begc:endc))                   ; fert_start    (:) = 0
    allocate(fert_end                         (begc:endc))                   ; fert_end      (:) = 0
    allocate(this%r_mort_cal_patch                (begp:endp))               ; this%r_mort_cal_patch   (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp, no_snow_normal
    !
    ! !ARGUMENTS:
    class(cnstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    character(8)      :: vr_suffix
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    if ( crop_prog) then
       this%gddmaturity_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDDHARV', units='ddays', &
            avgflag='A', long_name='Growing degree days (gdd) needed to harvest', &
            ptr_patch=this%gddmaturity_patch, default='inactive')
    end if

    this%croot_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='CROOT_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from coarse roots', &
         ptr_patch=this%croot_prof_patch, default='inactive')

    this%froot_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='FROOT_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from fine roots', &
         ptr_patch=this%froot_prof_patch, default='inactive')

    this%leaf_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='LEAF_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from leaves', &
         ptr_patch=this%leaf_prof_patch, default='inactive')

    this%stem_prof_patch(begp:endp,:) = spval
    call hist_addfld_decomp (fname='STEM_PROF', units='1/m',  type2d='levdcmp', &
         avgflag='A', long_name='profile for litter C and N inputs from stems', &
         ptr_patch=this%stem_prof_patch, default='inactive')

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

    this%scalaravg_col(begc:endc,:) = spval
    call hist_addfld_decomp(fname='SCALARAVG'//trim(vr_suffix), units='fraction', &
         type2d='levdcmp', avgflag='A', long_name='average of decomposition scalar',&
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
 
    this%annavg_t2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='ANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average 2m air temperature', &
         ptr_patch=this%annavg_t2m_patch, default='inactive')

    this%tempavg_t2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEMPAVG_T2M', units='K', &
         avgflag='A', long_name='temporary average 2m air temperature', &
         ptr_patch=this%tempavg_t2m_patch, default='inactive')

    this%dormant_flag_patch(begp:endp) = spval
    call hist_addfld1d (fname='DORMANT_FLAG', units='none', &
         avgflag='A', long_name='dormancy flag', &
         ptr_patch=this%dormant_flag_patch, default='inactive')

    this%days_active_patch(begp:endp) = spval
    call hist_addfld1d (fname='DAYS_ACTIVE', units='days', &
         avgflag='A', long_name='number of days since last dormancy', &
         ptr_patch=this%days_active_patch, default='inactive')

    this%onset_flag_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_FLAG', units='none', &
         avgflag='A', long_name='onset flag', &
         ptr_patch=this%onset_flag_patch, default='inactive')

    this%onset_counter_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_COUNTER', units='days', &
         avgflag='A', long_name='onset days counter', &
         ptr_patch=this%onset_counter_patch, default='inactive')

    this%onset_gddflag_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_GDDFLAG', units='none', &
         avgflag='A', long_name='onset flag for growing degree day sum', &
         ptr_patch=this%onset_gddflag_patch, default='inactive')

    this%onset_fdd_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_FDD', units='C degree-days', &
         avgflag='A', long_name='onset freezing degree days counter', &
         ptr_patch=this%onset_fdd_patch, default='inactive')

    this%onset_gdd_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_GDD', units='C degree-days', &
         avgflag='A', long_name='onset growing degree days', &
         ptr_patch=this%onset_gdd_patch, default='inactive')

    this%onset_swi_patch(begp:endp) = spval
    call hist_addfld1d (fname='ONSET_SWI', units='none', &
         avgflag='A', long_name='onset soil water index', &
         ptr_patch=this%onset_swi_patch, default='inactive')

    this%offset_flag_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_FLAG', units='none', &
         avgflag='A', long_name='offset flag', &
         ptr_patch=this%offset_flag_patch, default='inactive')

    this%offset_counter_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_COUNTER', units='days', &
         avgflag='A', long_name='offset days counter', &
         ptr_patch=this%offset_counter_patch, default='inactive')

    this%offset_fdd_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_FDD', units='C degree-days', &
         avgflag='A', long_name='offset freezing degree days counter', &
         ptr_patch=this%offset_fdd_patch, default='inactive')

    this%offset_swi_patch(begp:endp) = spval
    call hist_addfld1d (fname='OFFSET_SWI', units='none', &
         avgflag='A', long_name='offset soil water index', &
         ptr_patch=this%offset_swi_patch, default='inactive')

    this%lgsf_patch(begp:endp) = spval
    call hist_addfld1d (fname='LGSF', units='proportion', &
         avgflag='A', long_name='long growing season factor', &
         ptr_patch=this%lgsf_patch, default='inactive')

    this%bglfr_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname='BGLFR_LEAF', units='1/s', &
         avgflag='A', long_name='background leaf litterfall rate', &
         ptr_patch=this%bglfr_leaf_patch, default='inactive')

    this%bglfr_froot_patch(begp:endp) = spval
    call hist_addfld1d (fname='BGLFR_FROOT', units='1/s', &
         avgflag='A', long_name='background fine root litterfall rate', &
         ptr_patch=this%bglfr_froot_patch, default='inactive')

    this%bgtr_patch(begp:endp) = spval
    call hist_addfld1d (fname='BGTR', units='1/s', &
         avgflag='A', long_name='background transfer growth rate', &
         ptr_patch=this%bgtr_patch, default='inactive')

    this%alloc_pnow_patch(begp:endp) = spval
    call hist_addfld1d (fname='ALLOC_PNOW', units='proportion', &
         avgflag='A', long_name='fraction of current allocation to display as new growth', &
         ptr_patch=this%alloc_pnow_patch, default='inactive')

    this%c_allometry_patch(begp:endp) = spval
    call hist_addfld1d (fname='C_ALLOMETRY', units='none', &
         avgflag='A', long_name='C allocation index', &
         ptr_patch=this%c_allometry_patch, default='inactive')

    this%n_allometry_patch(begp:endp) = spval
    call hist_addfld1d (fname='N_ALLOMETRY', units='none', &
         avgflag='A', long_name='N allocation index', &
         ptr_patch=this%n_allometry_patch, default='inactive')

    this%tempsum_potential_gpp_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEMPSUM_POTENTIAL_GPP', units='gC/m^2/yr', &
         avgflag='A', long_name='temporary annual sum of potential GPP', &
         ptr_patch=this%tempsum_potential_gpp_patch, default='inactive')

    this%annsum_potential_gpp_patch(begp:endp) = spval
    call hist_addfld1d (fname='ANNSUM_POTENTIAL_GPP', units='gN/m^2/yr', &
         avgflag='A', long_name='annual sum of potential GPP', &
         ptr_patch=this%annsum_potential_gpp_patch, default='inactive')

    this%tempmax_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEMPMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='temporary annual max of retranslocated N pool', &
         ptr_patch=this%tempmax_retransn_patch, default='inactive')

    this%annmax_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='ANNMAX_RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='annual max of retranslocated N pool', &
         ptr_patch=this%annmax_retransn_patch, default='inactive')

    this%downreg_patch(begp:endp) = spval
    call hist_addfld1d (fname='DOWNREG', units='proportion', &
         avgflag='A', long_name='fractional reduction in GPP due to N limitation', &
         ptr_patch=this%downreg_patch, default='inactive')


    !! add phosphorus -X.YANG
    this%p_allometry_patch(begp:endp) = spval
    call hist_addfld1d (fname='P_ALLOMETRY', units='none', &
         avgflag='A', long_name='P allocation index', &
         ptr_patch=this%p_allometry_patch, default='inactive')

    this%tempmax_retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='TEMPMAX_RETRANSP', units='gP/m^2', &
         avgflag='A', long_name='temporary annual max of retranslocated P pool', &
         ptr_patch=this%tempmax_retransp_patch, default='inactive')

    this%annmax_retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='ANNMAX_RETRANSP', units='gP/m^2', &
         avgflag='A', long_name='annual max of retranslocated P pool', &
         ptr_patch=this%annmax_retransp_patch, default='inactive')

    this%cn_scalar(begp:endp) = spval
    call hist_addfld1d (fname='cn_scalar', units='', &
       avgflag='A', long_name='N limitation factor', &
       ptr_patch=this%cn_scalar, default='active')
         
    this%cp_scalar(begp:endp) = spval
    call hist_addfld1d (fname='cp_scalar', units='', &
       avgflag='A', long_name='P limitation factor', &
       ptr_patch=this%cp_scalar, default='active')

    this%cn_scalar_runmean(begp:endp) = spval
    call hist_addfld1d (fname='nlim_m', units='', &
       avgflag='A', long_name='runmean N limitation factor', &
       ptr_patch=this%cn_scalar_runmean, default='active')

    this%cp_scalar_runmean(begp:endp) = spval
    call hist_addfld1d (fname='plim_m', units='', &
       avgflag='A', long_name='runmean P limitation factor', &
       ptr_patch=this%cp_scalar_runmean, default='active')

    this%r_mort_cal_patch(begp:endp) = spval
    call hist_addfld1d (fname='R_MORT_CAL', units='none', &
         avgflag='A', long_name='calcualted annual mortality rate', &
         ptr_patch=this%r_mort_cal_patch, default='inactive')


  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
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

    integer     ,pointer     :: fert_type_rdin (:)
    integer     ,pointer     :: fert_continue_rdin (:)
    real(r8)    ,pointer     :: fert_dose_rdin (:,:)
    ! soil phosphorus pool Qing Z. 2017
    real(r8) ,pointer  :: labp_g (:)                       ! read in - LABILE_P
    real(r8) ,pointer  :: secp_g (:)                       ! read in - SECONDARY_P
    real(r8) ,pointer  :: occp_g (:)                       ! read in - OCCLUDED_P
    real(r8) ,pointer  :: prip_g (:)                       ! read in - APATITE_P
    integer     ,pointer     :: fert_start_rdin (:)
    integer     ,pointer     :: fert_end_rdin (:)
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

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
       g = col_pp%gridcell(c)
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
       g = col_pp%gridcell(c)
       this%peatf_lf_col(c) = peatf(g)
    end do
    deallocate(peatf)

    ! --------------------------------------------------------------------
    ! Read in soilorder data 
    ! --------------------------------------------------------------------

    if ( (nu_com .eq. 'RD' .or. nu_com .eq. 'ECA') .and. (use_cn .and. .not. use_fates .and. .not. use_crop) )  then 
       allocate(soilorder_rdin(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncid, varname='SOIL_ORDER', flag='read',data=soilorder_rdin, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: SOIL_ORDER NOT on surfdata file'//errMsg(__FILE__, __LINE__))
       end if
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          this%isoilorder(c) = soilorder_rdin(g)
       end do
       deallocate(soilorder_rdin)

    else
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          this%isoilorder(c) = 12
       end do 
    end if

    ! --------------------------------------------------------------------
    ! forest fertilization experiments info, Q. Z. 2017
    ! --------------------------------------------------------------------
    if (forest_fert_exp) then
       allocate(fert_type_rdin(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncid, varname='fert_type', flag='read',data=fert_type_rdin, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: fert_type NOT on surfdata file'//errMsg(__FILE__, __LINE__))
       end if
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          fert_type(c) = fert_type_rdin(g)
       end do
       deallocate(fert_type_rdin)

       allocate(fert_continue_rdin(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncid, varname='fert_continue', flag='read',data=fert_continue_rdin, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: fert_continue NOT on surfdata file'//errMsg(__FILE__, __LINE__))
       end if
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          fert_continue(c) = fert_continue_rdin(g)
       end do
       deallocate(fert_continue_rdin)

       allocate(fert_dose_rdin(bounds%begg:bounds%endg,1:12))
       call ncd_io(ncid=ncid, varname='fert_dose', flag='read',data=fert_dose_rdin, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: fert_dose NOT on surfdata file'//errMsg(__FILE__, __LINE__))
       end if
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          do j = 1 , 12
             fert_dose(c,j) = fert_dose_rdin(g,j)
          end do
       end do
       deallocate(fert_dose_rdin)

       allocate(fert_start_rdin(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncid, varname='fert_startyr', flag='read',data=fert_start_rdin, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: fert_start NOT on surfdata file'//errMsg(__FILE__, __LINE__))
       end if
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          fert_start(c) = fert_start_rdin(g)
       end do
       deallocate(fert_start_rdin)

       allocate(fert_end_rdin(bounds%begg:bounds%endg))
       call ncd_io(ncid=ncid, varname='fert_endyr', flag='read',data=fert_end_rdin, dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun(msg=' ERROR: fert_endyr NOT on surfdata file'//errMsg(__FILE__, __LINE__))
       end if
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          fert_end(c) = fert_end_rdin(g)
       end do
       deallocate(fert_end_rdin)

    end if

    ! --------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! Read in ABM data 
    ! --------------------------------------------------------------------

    allocate(abm(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='abm', flag='read', data=abm, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: abm NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
       this%abm_lf_col(c) = abm(g)
    end do
    deallocate(abm)

    ! Close file

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) 'Successfully read fmax, soil color, sand and clay boundary data'
       write(iulog,*)
    endif
    
    if (masterproc) then
       write(iulog,*) 'Attempting to read initial phosphorus pools data .....'
    end if

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! Read soil phosphorus pool Qing Z. 2017 
    this%pdatasets_present = .true.
    allocate(labp_g(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='LABILE_P', flag='read', data=labp_g, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       this%pdatasets_present = .false.
    else
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          this%labp_col(c) = labp_g(g)
       end do
    end if
    deallocate(labp_g)

    allocate(secp_g(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='SECONDARY_P', flag='read', data=secp_g, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       this%pdatasets_present = .false.
    else
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          this%secp_col(c) = secp_g(g)
       end do
    end if
    deallocate(secp_g)

    allocate(occp_g(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='OCCLUDED_P', flag='read', data=occp_g, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       this%pdatasets_present = .false.
    else
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          this%occp_col(c) = occp_g(g)
       end do
    end if
    deallocate(occp_g)

    allocate(prip_g(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='APATITE_P', flag='read', data=prip_g, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       this%pdatasets_present = .false.
    else
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          this%prip_col(c) = prip_g(g)
       end do
    end if
    deallocate(prip_g)  

    ! --------------------------------------------------------------------
    ! Initialize terms needed for dust model
    ! TODO - move these terms to DUSTMod module variables 
    ! --------------------------------------------------------------------
       
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%ifspecial(l)) then
          this%annsum_counter_col (c) = spval
          this%annavg_t2m_col     (c) = spval
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
             this%scalaravg_col(c,j) = spval
          end do
       end if

       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
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
          !this%scalaravg_col(c,1:nlevdecomp_full)       = 0._r8

          ! initialize the profiles for converting to vertically resolved carbon pools
          this%nfixation_prof_col(c,1:nlevdecomp_full)  = 0._r8 
          this%ndep_prof_col(c,1:nlevdecomp_full)       = 0._r8 
          this%pdep_prof_col(c,1:nlevdecomp_full)       = 0._r8 
       end if
    end do

    ! ecophysiological variables
    ! phenology variables

    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)
       this%rc14_atm_patch(p)              = c14ratio 

       if (lun_pp%ifspecial(l)) then
          this%annavg_t2m_patch  (p)          = spval
          this%tempavg_t2m_patch (p)          = spval
          this%dormant_flag_patch(p)          = spval
          this%days_active_patch(p)           = spval
          this%onset_flag_patch(p)            = spval
          this%onset_counter_patch(p)         = spval
          this%onset_gddflag_patch(p)         = spval
          this%onset_fdd_patch(p)             = spval
          this%onset_gdd_patch(p)             = spval
          this%onset_swi_patch(p)             = spval
          this%offset_flag_patch(p)           = spval
          this%offset_counter_patch(p)        = spval
          this%offset_fdd_patch(p)            = spval
          this%offset_swi_patch(p)            = spval
          this%grain_flag_patch(p)            = spval
          this%lgsf_patch(p)                  = spval
          this%bglfr_patch(p)                 = spval
          this%bglfr_leaf_patch(p)            = spval
          this%bglfr_froot_patch(p)           = spval
          this%bgtr_patch(p)                  = spval
          this%alloc_pnow_patch(p)            = spval
          this%c_allometry_patch(p)           = spval
          this%n_allometry_patch(p)           = spval
          this%tempsum_potential_gpp_patch(p) = spval
          this%annsum_potential_gpp_patch(p)  = spval
          this%tempmax_retransn_patch(p)      = spval
          this%annmax_retransn_patch(p)       = spval
          this%downreg_patch(p)               = spval

          this%p_allometry_patch(p)           = spval
          this%tempmax_retransp_patch(p)      = spval
          this%annmax_retransp_patch(p)       = spval

          this%r_mort_cal_patch(p)           = spval

 
       end if
    end do
       
    ! ecophysiological variables

    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          this%rc14_atm_patch(p) = c14ratio

          ! phenology variables
          this%dormant_flag_patch(p)   = 1._r8
          this%days_active_patch(p)    = 0._r8
          this%onset_flag_patch(p)     = 0._r8
          this%onset_counter_patch(p)  = 0._r8
          this%onset_gddflag_patch(p)  = 0._r8
          this%onset_fdd_patch(p)      = 0._r8
          this%onset_gdd_patch(p)      = 0._r8
          this%onset_swi_patch(p)      = 0._r8
          this%offset_flag_patch(p)    = 0._r8
          this%offset_counter_patch(p) = 0._r8
          this%offset_fdd_patch(p)     = 0._r8
          this%offset_swi_patch(p)     = 0._r8
          this%lgsf_patch(p)           = 0._r8
          this%bglfr_patch(p)          = 0._r8
          this%bglfr_leaf_patch(p)     = 0._r8
          this%bglfr_froot_patch(p)    = 0._r8
          this%bgtr_patch(p)           = 0._r8
          this%annavg_t2m_patch(p)     = 280._r8
          this%tempavg_t2m_patch(p)    = 0._r8
          this%grain_flag_patch(p)     = 0._r8

          ! non-phenology variables
          this%alloc_pnow_patch(p)            = 1._r8
          this%c_allometry_patch(p)           = 0._r8
          this%n_allometry_patch(p)           = 0._r8
          this%tempsum_potential_gpp_patch(p) = 0._r8
          this%annsum_potential_gpp_patch(p)  = 0._r8
          this%tempmax_retransn_patch(p)      = 0._r8
          this%annmax_retransn_patch(p)       = 0._r8
          this%downreg_patch(p)               = 0._r8

          this%p_allometry_patch(p)           = 0._r8
          this%tempmax_retransp_patch(p)      = 0._r8
          this%annmax_retransp_patch(p)       = 0._r8

          this%r_mort_cal_patch(p)           = 0._r8


       end if
    end do

    ! fire variables

    do c = bounds%begc,bounds%endc
       this%lfc2_col(c) = 0._r8
    end do

  end subroutine initCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use restUtilMod
    use ncdio_pio
    use fileutils  , only : getfil
    !
    ! !ARGUMENTS:
    class(cnstate_type)              :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer, pointer :: temp1d(:) ! temporary
    integer          :: p,j,c,i,g   ! indices
    logical          :: readvar   ! determine if variable is on initial file
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    !-----------------------------------------------------------------------
  
    call restartvar(ncid=ncid, flag=flag, varname='dormant_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='dormancy flag', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%dormant_flag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='days_active', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='number of days since last dormancy', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%days_active_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='flag if critical growing degree-day sum is exceeded', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_flag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_counter_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_gddflag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset flag for growing degree day sum', units='' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_gddflag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_fdd_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_gdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset growing degree days', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_gdd_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset soil water index', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_swi_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset flag', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_flag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_counter_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_fdd_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%offset_swi_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='lgsf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lgsf_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='bglfr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bglfr_patch)

    call restartvar(ncid=ncid, flag=flag, varname='bglfr_leaf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bglfr_leaf_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='bglfr_froot', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bglfr_froot_patch)

    call restartvar(ncid=ncid, flag=flag, varname='bgtr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bgtr_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='annavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annavg_t2m_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempavg_t2m_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='alloc_pnow', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%alloc_pnow_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='c_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%c_allometry_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='n_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%n_allometry_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempsum_potential_gpp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='annsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_potential_gpp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempmax_retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='annmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annmax_retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='downreg', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%downreg_patch) 


    call restartvar(ncid=ncid, flag=flag, varname='p_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%p_allometry_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempmax_retransp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempmax_retransp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='annmax_retransp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annmax_retransp_patch) 

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

    ptr2d => this%scalaravg_col
    call restartvar(ncid=ncid, flag=flag, varname='scalaravg_col', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='fraction of potential immobilization of phosphorus', units='unitless', &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)

    if (crop_prog) then

       call restartvar(ncid=ncid, flag=flag,  varname='restyear', xtype=ncd_int,  &
            long_name='Number of years prognostic crop ran', units="years", &
            interpinic_flag='copy', readvar=readvar, data=this%CropRestYear)
       if (flag=='read' .and. readvar)  then
          call checkDates( )
       end if

       call restartvar(ncid=ncid, flag=flag,  varname='htmx', xtype=ncd_double,  &
            dim1name='pft', long_name='max height attained by a crop during year', units='m', &
            interpinic_flag='interp', readvar=readvar, data=this%htmx_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='peaklai', xtype=ncd_int,  &
            dim1name='pft', long_name='Flag if at max allowed LAI or not', &
            flag_values=(/0,1/), nvalid_range=(/0,1/), &
            flag_meanings=(/'NOT-at-peak', 'AT_peak-LAI' /) , &
            interpinic_flag='interp', readvar=readvar, data=this%peaklai_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='idop', xtype=ncd_int,  &
            dim1name='pft', long_name='Date of planting', units='jday', nvalid_range=(/1,366/), & 
            interpinic_flag='interp', readvar=readvar, data=this%idop_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='aleaf', xtype=ncd_double,  &
            dim1name='pft', long_name='leaf allocation coefficient', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%aleaf_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='aleafi', xtype=ncd_double,  &
            dim1name='pft', long_name='Saved leaf allocation coefficient from phase 2', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%aleafi_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='astem', xtype=ncd_double,  &
            dim1name='pft', long_name='stem allocation coefficient', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%astem_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='astemi', xtype=ncd_double,  &
            dim1name='pft', long_name='Saved stem allocation coefficient from phase 2', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%astemi_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='hdidx', xtype=ncd_double,  &
            dim1name='pft', long_name='cold hardening index', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%hdidx_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cumvd', xtype=ncd_double,  &
            dim1name='pft', long_name='cumulative vernalization d', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cumvd_patch)

      call restartvar(ncid=ncid, flag=flag,  varname='gddmaturity', xtype=ncd_double,  &
            dim1name='pft', long_name='Growing degree days needed to harvest', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gddmaturity_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='huileaf', xtype=ncd_double,  &
            dim1name='pft', long_name='heat unit index needed from planting to leaf emergence', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%huileaf_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='huigrain', xtype=ncd_double,  &
            dim1name='pft', long_name='heat unit index needed to reach vegetative maturity', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%huigrain_patch)

       call restartvar(ncid=ncid, flag=flag, varname='grain_flag', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_flag_patch)
    end if

    if (use_c14) then
       call restartvar(ncid=ncid, flag=flag, varname='rc14_atm', xtype=ncd_double,  &
            dim1name='pft',    long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%rc14_atm_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%rc14_atm with atmospheric c14 value'
          do i = bounds%begp, bounds%endp
             this%rc14_atm_patch(i) = c14ratio
          end do
       end if
    end if

    call restartvar(ncid=ncid, flag=flag, varname='cn_scalar', xtype=ncd_double,  &
            dim1name='pft', long_name='cn_scalar', units='-', &
            interpinic_flag='interp', readvar=readvar, data=this%cn_scalar)
    call restartvar(ncid=ncid, flag=flag, varname='cp_scalar', xtype=ncd_double,  &
            dim1name='pft', long_name='cp_scalar', units='-', &
            interpinic_flag='interp', readvar=readvar, data=this%cp_scalar)

    call restartvar(ncid=ncid, flag=flag, varname='nlim_m', xtype=ncd_double,  &
            dim1name='pft', long_name='cn_scalar_runmean', units='-', &
            interpinic_flag='interp', readvar=readvar, data=this%cn_scalar_runmean)
    call restartvar(ncid=ncid, flag=flag, varname='plim_m', xtype=ncd_double,  &
            dim1name='pft', long_name='cp_scalar_runmean', units='-', &
            interpinic_flag='interp', readvar=readvar, data=this%cp_scalar_runmean)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    ! !USES
    use accumulMod       , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(cnstate_type)           :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    this%cn_scalar_runmean(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='nlim_m', units='-',                                              &
         desc='runing average of N limitation strength',  accum_type='runmean', accum_period=-7300,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%cp_scalar_runmean(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='plim_m', units='-',                                              &
         desc='runing average of P limitation strength',  accum_type='runmean', accum_period=-7300,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)

    ! !USES
    use accumulMod       , only : init_accum_field, extract_accum_field
    use clm_time_manager , only : get_nstep
    use clm_varctl       , only : nsrest
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(cnstate_type)           :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level pft field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="extract_accum_hist allocation error for rbufslp"//&
            errMsg(__FILE__, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('nlim_m', rbufslp, nstep)
    this%cn_scalar_runmean(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('plim_m', rbufslp, nstep)
    this%cp_scalar_runmean(begp:endp) = rbufslp(begp:endp)
    deallocate(rbufslp)
  
  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod       , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(cnstate_type)                    :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: m,g,l,c,p                 ! indices
    integer :: ier                       ! error status
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    logical :: end_cd                    ! temporary for is_end_curr_day() value
    integer :: begp, endp
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract
    do p = begp,endp
       rbufslp(p) = this%cn_scalar(p)
    end do
    call update_accum_field  ('nlim_m' , rbufslp             , nstep)
    call extract_accum_field ('nlim_m' , this%cn_scalar_runmean  , nstep)
    do p = begp,endp
       rbufslp(p) = this%cp_scalar(p)
    end do
    call update_accum_field  ('plim_m' , rbufslp             , nstep)
    call extract_accum_field ('plim_m' , this%cp_scalar_runmean  , nstep)
    deallocate(rbufslp)

  end subroutine UpdateAccVars

  !-----------------------------------------------------------------------
  subroutine CropRestIncYear (this)
    !
    ! !DESCRIPTION: 
    ! Increment the crop restart year, if appropriate
    !
    ! This routine should be called every time step, but only once per clump (to avoid
    ! inadvertently updating nyrs multiple times)
    !
    ! !USES:
    use clm_varpar       , only : crop_prog
    use clm_time_manager , only : get_curr_date, is_first_step
    !
    ! !ARGUMENTS:
    class(cnstate_type) :: this
    !
    ! !LOCAL VARIABLES:
    integer kyr   ! current year
    integer kmo   ! month of year  (1, ..., 12)
    integer kda   ! day of month   (1, ..., 31)
    integer mcsec ! seconds of day (0, ..., seconds/day)
    !-----------------------------------------------------------------------

    ! Update restyear only when running with prognostic crop
    if ( crop_prog )then

       ! Update restyear when it's the start of a new year - but don't do that at the
       ! very start of the run
       call get_curr_date (   kyr, kmo, kda, mcsec)
       if ((kmo == 1 .and. kda == 1 .and. mcsec == 0) .and. .not. is_first_step()) then
          this%CropRestYear = this%CropRestYear + 1
       end if

    end if

  end subroutine CropRestIncYear

  !-----------------------------------------------------------------------
  subroutine checkDates( )
    !
    ! !DESCRIPTION: 
    ! Make sure the dates are compatible. The date given to startup the model
    ! and the date on the restart file must be the same although years can be
    ! different. The dates need to be checked when the restart file is being
    ! read in for a startup or branch case (they are NOT allowed to be different
    ! for a restart case).
    !
    ! For the prognostic crop model the date of planting is tracked and growing
    ! degree days is tracked (with a 20 year mean) -- so shifting the start dates
    ! messes up these bits of saved information.
    !
    ! !ARGUMENTS:
    use clm_time_manager, only : get_driver_start_ymd, get_start_date
    use clm_varctl      , only : iulog
    use clm_varctl      , only : nsrest, nsrBranch, nsrStartup
    !
    ! !LOCAL VARIABLES:
    integer :: stymd       ! Start date YYYYMMDD from driver
    integer :: styr        ! Start year from driver
    integer :: stmon_day   ! Start date MMDD from driver
    integer :: rsmon_day   ! Restart date MMDD from restart file
    integer :: rsyr        ! Restart year from restart file
    integer :: rsmon       ! Restart month from restart file
    integer :: rsday       ! Restart day from restart file
    integer :: tod         ! Restart time of day from restart file
    character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)' ! log output format
    character(len=32) :: subname = 'CropRest::checkDates'
    !-----------------------------------------------------------------------
    !
    ! If branch or startup make sure the startdate is compatible with the date
    ! on the restart file.
    !
    if ( nsrest == nsrBranch .or. nsrest == nsrStartup )then
       stymd       = get_driver_start_ymd()
       styr        = stymd / 10000
       stmon_day   = stymd - styr*10000
       call get_start_date( rsyr, rsmon, rsday, tod )
       rsmon_day = rsmon*100 + rsday
       if ( masterproc ) &
            write(iulog,formDate) 'Date on the restart file is: ', rsyr, rsmon, rsday
       if ( stmon_day /= rsmon_day )then
          write(iulog,formDate) 'Start date is: ', styr, stmon_day/100, &
               (stmon_day - stmon_day/100)
          call endrun(msg=' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
               ' and the date on the restart file needs to match (years can be different)'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if

  end subroutine checkDates

end module CNStateType
