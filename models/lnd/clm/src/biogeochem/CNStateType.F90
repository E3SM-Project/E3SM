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
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : pft                
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PRIVATE MEMBER FUNCTIONS: 
  private :: checkDates
  ! 
  ! !PUBLIC TYPES:
  type, public :: cnstate_type

     integer  , pointer :: burndate_patch              (:)     ! patch crop burn date
     real(r8) , pointer :: lfpftd_patch                (:)     ! patch decrease of pft weight (0-1) on the column for the timestep 

     ! Prognostic crop model -  Note that cropplant and harvdate could be 2D to facilitate rotation
     real(r8) , pointer :: hdidx_patch                 (:)     ! patch cold hardening index?
     real(r8) , pointer :: cumvd_patch                 (:)     ! patch cumulative vernalization d?ependence?
     real(r8) , pointer :: vf_patch                    (:)     ! patch vernalization factor for cereal
     real(r8) , pointer :: gddmaturity_patch           (:)     ! patch growing degree days (gdd) needed to harvest (ddays)
     real(r8) , pointer :: huileaf_patch               (:)     ! patch heat unit index needed from planting to leaf emergence
     real(r8) , pointer :: huigrain_patch              (:)     ! patch heat unit index needed to reach vegetative maturity
     real(r8) , pointer :: aleafi_patch                (:)     ! patch saved leaf allocation coefficient from phase 2
     real(r8) , pointer :: astemi_patch                (:)     ! patch saved stem allocation coefficient from phase 2
     real(r8) , pointer :: aleaf_patch                 (:)     ! patch leaf allocation coefficient
     real(r8) , pointer :: astem_patch                 (:)     ! patch stem allocation coefficient
     logical  , pointer :: croplive_patch              (:)     ! patch Flag, true if planted, not harvested
     logical  , pointer :: cropplant_patch             (:)     ! patch Flag, true if planted
     integer  , pointer :: harvdate_patch              (:)     ! patch harvest date
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

     real(r8) , pointer :: rf_decomp_cascade_col       (:,:,:) ! col respired fraction in decomposition step (frac)
     real(r8) , pointer :: pathfrac_decomp_cascade_col (:,:,:) ! col what fraction of C leaving a given pool passes through a given transition (frac) 
     real(r8) , pointer :: nfixation_prof_col          (:,:)   ! col (1/m) profile for N fixation additions 
     real(r8) , pointer :: ndep_prof_col               (:,:)   ! col (1/m) profile for N fixation additions 
     real(r8) , pointer :: som_adv_coef_col            (:,:)   ! col SOM advective flux (m/s) 
     real(r8) , pointer :: som_diffus_coef_col         (:,:)   ! col SOM diffusivity due to bio/cryo-turbation (m2/s) 

     real(r8) , pointer :: tempavg_t2m_patch           (:)     ! patch temporary average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_patch            (:)     ! patch annual average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_col              (:)     ! col annual average of 2m air temperature, averaged from pft-level (K)
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

     integer           :: CropRestYear                        ! restart year from initial conditions file - increment as time elapses

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, public  :: CropRestIncYear
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

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
    allocate(this%vf_patch            (begp:endp))                   ; this%vf_patch            (:)   = 0.0_r8
    allocate(this%gddmaturity_patch   (begp:endp))                   ; this%gddmaturity_patch   (:)   = spval
    allocate(this%huileaf_patch       (begp:endp))                   ; this%huileaf_patch       (:)   = nan
    allocate(this%huigrain_patch      (begp:endp))                   ; this%huigrain_patch      (:)   = nan
    allocate(this%aleafi_patch        (begp:endp))                   ; this%aleafi_patch        (:)   = nan
    allocate(this%astemi_patch        (begp:endp))                   ; this%astemi_patch        (:)   = nan
    allocate(this%aleaf_patch         (begp:endp))                   ; this%aleaf_patch         (:)   = nan
    allocate(this%astem_patch         (begp:endp))                   ; this%astem_patch         (:)   = nan
    allocate(this%croplive_patch      (begp:endp))                   ; this%croplive_patch      (:)   = .false.
    allocate(this%cropplant_patch     (begp:endp))                   ; this%cropplant_patch     (:)   = .false.
    allocate(this%harvdate_patch      (begp:endp))                   ; this%harvdate_patch      (:)   = huge(1) 
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

    allocate(this%rf_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); 
    this%rf_decomp_cascade_col(:,:,:) = nan

    allocate(this%pathfrac_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));     
    this%pathfrac_decomp_cascade_col(:,:,:) = nan

    allocate(this%nfixation_prof_col  (begc:endc,1:nlevdecomp_full)) ; this%nfixation_prof_col  (:,:) = spval
    allocate(this%ndep_prof_col       (begc:endc,1:nlevdecomp_full)) ; this%ndep_prof_col       (:,:) = spval
    allocate(this%som_adv_coef_col    (begc:endc,1:nlevdecomp_full)) ; this%som_adv_coef_col    (:,:) = spval
    allocate(this%som_diffus_coef_col (begc:endc,1:nlevdecomp_full)) ; this%som_diffus_coef_col (:,:) = spval

    allocate(this%tempavg_t2m_patch   (begp:endp))                   ; this%tempavg_t2m_patch   (:)   = nan
    allocate(this%annsum_counter_col  (begc:endc))                   ; this%annsum_counter_col  (:)   = nan
    allocate(this%annavg_t2m_col      (begc:endc))                   ; this%annavg_t2m_col      (:)   = nan
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
            avgflag='A', long_name='fraction of potential immobilization', &
            ptr_col=this%fpi_col)
    endif

    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif
    this%fpi_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='FPI'//trim(vr_suffix), units='proportion', type2d='levdcmp', & 
         avgflag='A', long_name='fraction of potential immobilization', &
         ptr_col=this%fpi_vr_col)

    this%fpg_col(begc:endc) = spval
    call hist_addfld1d (fname='FPG', units='proportion', &
         avgflag='A', long_name='fraction of potential gpp', &
         ptr_col=this%fpg_col)

    this%annsum_counter_col(begc:endc) = spval
    call hist_addfld1d (fname='ANNSUM_COUNTER', units='s', &
         avgflag='A', long_name='seconds since last annual accumulator turnover', &
         ptr_col=this%annsum_counter_col, default='inactive')

    this%annavg_t2m_col(begc:endc) = spval
    call hist_addfld1d (fname='CANNAVG_T2M', units='K', &
         avgflag='A', long_name='annual average of 2m air temperature', &
         ptr_col=this%annavg_t2m_col, default='inactive')

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

    this%bglfr_patch(begp:endp) = spval
    call hist_addfld1d (fname='BGLFR', units='1/s', &
         avgflag='A', long_name='background litterfall rate', &
         ptr_patch=this%bglfr_patch, default='inactive')

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
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
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
          do j = 1,nlevdecomp_full
             this%fpi_vr_col(c,j) = spval
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
          this%som_adv_coef_col(c,1:nlevdecomp_full)    = 0._r8 
          this%som_diffus_coef_col(c,1:nlevdecomp_full) = 0._r8 

          ! initialize the profiles for converting to vertically resolved carbon pools
          this%nfixation_prof_col(c,1:nlevdecomp_full)  = 0._r8 
          this%ndep_prof_col(c,1:nlevdecomp_full)       = 0._r8 
       end if
    end do

    ! ecophysiological variables
    ! phenology variables

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       this%rc14_atm_patch(p)              = c14ratio 

       if (lun%ifspecial(l)) then
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
          this%bgtr_patch(p)                  = spval
          this%alloc_pnow_patch(p)            = spval
          this%c_allometry_patch(p)           = spval
          this%n_allometry_patch(p)           = spval
          this%tempsum_potential_gpp_patch(p) = spval
          this%annsum_potential_gpp_patch(p)  = spval
          this%tempmax_retransn_patch(p)      = spval
          this%annmax_retransn_patch(p)       = spval
          this%downreg_patch(p)               = spval
       end if
    end do
       
    ! ecophysiological variables

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

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
    !
    ! !ARGUMENTS:
    class(cnstate_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer, pointer :: temp1d(:) ! temporary
    integer          :: p,j,c,i   ! indices
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

    if (use_vertsoilc) then
       ptr2d => this%fpi_vr_col
       call restartvar(ncid=ncid, flag=flag, varname='fpi_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='fraction of potential immobilization',  units='unitless', &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%fpi_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
       call restartvar(ncid=ncid, flag=flag, varname='fpi', xtype=ncd_double,  &
            dim1name='column', &
            long_name='fraction of potential immobilization',  units='unitless', &
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

       call restartvar(ncid=ncid, flag=flag,  varname='vf', xtype=ncd_double,  &
            dim1name='pft', long_name='vernalization factor', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%vf_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cumvd', xtype=ncd_double,  &
            dim1name='pft', long_name='cumulative vernalization d', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cumvd_patch)

       allocate(temp1d(bounds%begp:bounds%endp))
       if (flag == 'write') then 
          do p= bounds%begp,bounds%endp
             if (this%croplive_patch(p)) then
                temp1d(p) = 1
             else
                temp1d(p) = 0
             end if
          end do
       end if
       call restartvar(ncid=ncid, flag=flag,  varname='croplive', xtype=ncd_log,  &
            dim1name='pft', &
            long_name='Flag that crop is alive, but not harvested', &
            interpinic_flag='interp', readvar=readvar, data=temp1d)
       if (flag == 'read') then 
          do p= bounds%begp,bounds%endp
             if (temp1d(p) == 1) then
                this%croplive_patch(p) = .true.
             else
                this%croplive_patch(p) = .false.
             end if
          end do
       end if
       deallocate(temp1d)

       allocate(temp1d(bounds%begp:bounds%endp))
       if (flag == 'write') then 
          do p= bounds%begp,bounds%endp
             if (this%cropplant_patch(p)) then
                temp1d(p) = 1
             else
                temp1d(p) = 0
             end if
          end do
       end if
       call restartvar(ncid=ncid, flag=flag,  varname='cropplant', xtype=ncd_log,  &
            dim1name='pft', &
            long_name='Flag that crop is planted, but not harvested' , &
            interpinic_flag='interp', readvar=readvar, data=temp1d)
       if (flag == 'read') then 
          do p= bounds%begp,bounds%endp
             if (temp1d(p) == 1) then
                this%cropplant_patch(p) = .true.
             else
                this%cropplant_patch(p) = .false.
             end if
          end do
       end if
       deallocate(temp1d)

       call restartvar(ncid=ncid, flag=flag,  varname='harvdate', xtype=ncd_int,  &
            dim1name='pft', long_name='harvest date', units='jday', nvalid_range=(/1,366/), & 
            interpinic_flag='interp', readvar=readvar, data=this%harvdate_patch)

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

  end subroutine Restart

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
