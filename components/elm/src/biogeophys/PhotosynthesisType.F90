module PhotosynthesisType

  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use elm_varpar     , only : nlevcan
  use clm_varctl     , only : iulog, use_cn, use_c13, use_c14
  use elm_varcon     , only : spval
  use LandunitType   , only : lun_pp                
  use VegetationType      , only : veg_pp                
  !
  implicit none
  save
  private
  !
  type, public :: photosyns_type   

     logical , pointer :: c3flag_patch      (:)   ! patch true if C3 and false if C4
     real(r8), pointer :: ac_patch          (:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: aj_patch          (:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: ap_patch          (:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: ag_patch          (:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: an_patch          (:,:) ! patch net leaf photosynthesis (umol CO2/m**2/s)   
     real(r8), pointer :: vcmax_z_patch     (:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer :: cp_patch          (:)   ! patch CO2 compensation point (Pa)
     real(r8), pointer :: kc_patch          (:)   ! patch Michaelis-Menten constant for CO2 (Pa)
     real(r8), pointer :: ko_patch          (:)   ! patch Michaelis-Menten constant for O2 (Pa)
     real(r8), pointer :: qe_patch          (:)   ! patch quantum efficiency, used only for C4 (mol CO2 / mol photons)
     real(r8), pointer :: tpu_z_patch       (:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer :: kp_z_patch        (:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer :: theta_cj_patch    (:)   ! patch empirical curvature parameter for ac, aj photosynthesis co-limitation
     real(r8), pointer :: bbb_patch         (:)   ! patch Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
     real(r8), pointer :: mbb_patch         (:)   ! patch Ball-Berry slope of conductance-photosynthesis relationship
     real(r8), pointer :: gs_mol_patch      (:,:) ! patch leaf stomatal conductance       (umol H2O/m**2/s)
     real(r8), pointer :: gb_mol_patch      (:)   ! patch leaf boundary layer conductance (umol H2O/m**2/s)
     real(r8), pointer :: rh_leaf_patch     (:)   ! patch fractional humidity at leaf surface (dimensionless)

     real(r8), pointer :: alphapsnsun_patch (:)   ! patch sunlit 13c fractionation ([])                      
     real(r8), pointer :: alphapsnsha_patch (:)   ! patch shaded 13c fractionation ([])                      
     real(r8), pointer :: rc13_canair_patch (:)   ! patch C13O2/C12O2 in canopy air                          
     real(r8), pointer :: rc13_psnsun_patch (:)   ! patch C13O2/C12O2 in sunlit canopy psn flux              
     real(r8), pointer :: rc13_psnsha_patch (:)   ! patch C13O2/C12O2 in shaded canopy psn flux              

     real(r8), pointer :: psnsun_patch      (:)   ! patch sunlit leaf photosynthesis     (umol CO2/m**2/s)       
     real(r8), pointer :: psnsha_patch      (:)   ! patch shaded leaf photosynthesis     (umol CO2/m**2/s)       
     real(r8), pointer :: c13_psnsun_patch  (:)   ! patch c13 sunlit leaf photosynthesis (umol 13CO2/m**2/s) 
     real(r8), pointer :: c13_psnsha_patch  (:)   ! patch c13 shaded leaf photosynthesis (umol 13CO2/m**2/s) 
     real(r8), pointer :: c14_psnsun_patch  (:)   ! patch c14 sunlit leaf photosynthesis (umol 14CO2/m**2/s) 
     real(r8), pointer :: c14_psnsha_patch  (:)   ! patch c14 shaded leaf photosynthesis (umol 14CO2/m**2/s) 

     real(r8), pointer :: psnsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer :: psnsha_z_patch    (:,:) ! patch canopy layer: shaded leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer :: psnsun_wc_patch   (:)   ! patch Rubsico-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: psnsha_wc_patch   (:)   ! patch Rubsico-limited shaded leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: psnsun_wj_patch   (:)   ! patch RuBP-limited sunlit leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer :: psnsha_wj_patch   (:)   ! patch RuBP-limited shaded leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer :: psnsun_wp_patch   (:)   ! patch product-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: psnsha_wp_patch   (:)   ! patch product-limited shaded leaf photosynthesis (umol CO2/m**2/s)

     real(r8), pointer :: fpsn_patch        (:)   ! patch photosynthesis                 (umol CO2/m**2/s)
     real(r8), pointer :: fpsn_wc_patch     (:)   ! patch Rubisco-limited photosynthesis (umol CO2/m**2/s)
     real(r8), pointer :: fpsn_wj_patch     (:)   ! patch RuBP-limited photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer :: fpsn_wp_patch     (:)   ! patch product-limited photosynthesis (umol CO2/m**2/s)

     real(r8), pointer :: lmrsun_patch      (:)   ! patch sunlit leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer :: lmrsha_patch      (:)   ! patch shaded leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer :: lmrsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
     real(r8), pointer :: lmrsha_z_patch    (:,:) ! patch canopy layer: shaded leaf maintenance respiration rate (umol CO2/m**2/s)

     real(r8) , pointer :: cisun_z_patch    (:,:) ! patch intracellular sunlit leaf CO2 (Pa)
     real(r8) , pointer :: cisha_z_patch    (:,:) ! patch intracellular shaded leaf CO2 (Pa)

     real(r8) , pointer :: rssun_z_patch    (:,:) ! patch canopy layer: sunlit leaf stomatal resistance (s/m)
     real(r8) , pointer :: rssha_z_patch    (:,:) ! patch canopy layer: shaded leaf stomatal resistance (s/m)
     real(r8) , pointer :: rssun_patch      (:)   ! patch sunlit stomatal resistance (s/m)
     real(r8) , pointer :: rssha_patch      (:)   ! patch shaded stomatal resistance (s/m)

     real(r8) , pointer :: psncanopy_patch  (:)   ! patch sunlit leaf photosynthesis (umol CO2 /m**2/ s) (ED specific)

     ! ED specific variables
     real(r8), pointer :: lmrcanopy_patch   (:)   ! sunlit leaf maintenance respiration rate (umol CO2/m**2/s) (ED specific)
     ! Plant hydraulic stress specific variables
     real(r8), pointer, public :: ac_phs_patch      (:,:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, public :: aj_phs_patch      (:,:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, public :: ap_phs_patch      (:,:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, public :: ag_phs_patch      (:,:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, public :: an_sun_patch      (:,:)   ! patch sunlit net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, public :: an_sha_patch      (:,:)   ! patch shaded net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, public :: vcmax_z_phs_patch (:,:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer, public :: kp_z_phs_patch    (:,:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer, public :: tpu_z_phs_patch   (:,:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer, public :: gs_mol_sun_patch  (:,:) ! patch sunlit leaf stomatal conductance (umol H2O/m**2/s)
     real(r8), pointer, public :: gs_mol_sha_patch  (:,:) ! patch shaded leaf stomatal conductance (umol H2O/m**2/s)

   contains

     procedure, public  :: Init
     procedure, public  :: Restart      
     procedure, public  :: TimeStepInit
     procedure, public  :: NewPatchInit
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type photosyns_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds) 

    class(photosyns_type) :: this 
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds) 
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%c3flag_patch      (begp:endp))           ; this%c3flag_patch      (:)   =.false.
    allocate(this%ac_patch          (begp:endp,1:nlevcan)) ; this%ac_patch          (:,:) = nan
    allocate(this%aj_patch          (begp:endp,1:nlevcan)) ; this%aj_patch          (:,:) = nan
    allocate(this%ap_patch          (begp:endp,1:nlevcan)) ; this%ap_patch          (:,:) = nan
    allocate(this%ag_patch          (begp:endp,1:nlevcan)) ; this%ag_patch          (:,:) = nan
    allocate(this%an_patch          (begp:endp,1:nlevcan)) ; this%an_patch          (:,:) = nan
    allocate(this%vcmax_z_patch     (begp:endp,1:nlevcan)) ; this%vcmax_z_patch     (:,:) = nan
    allocate(this%cp_patch          (begp:endp))           ; this%cp_patch          (:)   = nan
    allocate(this%kc_patch          (begp:endp))           ; this%kc_patch          (:)   = nan
    allocate(this%ko_patch          (begp:endp))           ; this%ko_patch          (:)   = nan
    allocate(this%qe_patch          (begp:endp))           ; this%qe_patch          (:)   = nan
    allocate(this%tpu_z_patch       (begp:endp,1:nlevcan)) ; this%tpu_z_patch       (:,:) = nan
    allocate(this%kp_z_patch        (begp:endp,1:nlevcan)) ; this%kp_z_patch        (:,:) = nan
    allocate(this%theta_cj_patch    (begp:endp))           ; this%theta_cj_patch    (:)   = nan
    allocate(this%bbb_patch         (begp:endp))           ; this%bbb_patch         (:)   = nan
    allocate(this%mbb_patch         (begp:endp))           ; this%mbb_patch         (:)   = nan
    allocate(this%gb_mol_patch      (begp:endp))           ; this%gb_mol_patch      (:)   = nan
    allocate(this%gs_mol_patch      (begp:endp,1:nlevcan)) ; this%gs_mol_patch      (:,:) = nan
    allocate(this%rh_leaf_patch     (begp:endp))           ; this%rh_leaf_patch     (:)   = nan

    allocate(this%psnsun_patch      (begp:endp))           ; this%psnsun_patch      (:)   = nan
    allocate(this%psnsha_patch      (begp:endp))           ; this%psnsha_patch      (:)   = nan
    allocate(this%c13_psnsun_patch  (begp:endp))           ; this%c13_psnsun_patch  (:)   = nan
    allocate(this%c13_psnsha_patch  (begp:endp))           ; this%c13_psnsha_patch  (:)   = nan
    allocate(this%c14_psnsun_patch  (begp:endp))           ; this%c14_psnsun_patch  (:)   = nan
    allocate(this%c14_psnsha_patch  (begp:endp))           ; this%c14_psnsha_patch  (:)   = nan

    allocate(this%psnsun_z_patch    (begp:endp,1:nlevcan)) ; this%psnsun_z_patch    (:,:) = nan
    allocate(this%psnsha_z_patch    (begp:endp,1:nlevcan)) ; this%psnsha_z_patch    (:,:) = nan
    allocate(this%psnsun_wc_patch   (begp:endp))           ; this%psnsun_wc_patch   (:)   = nan
    allocate(this%psnsha_wc_patch   (begp:endp))           ; this%psnsha_wc_patch   (:)   = nan
    allocate(this%psnsun_wj_patch   (begp:endp))           ; this%psnsun_wj_patch   (:)   = nan
    allocate(this%psnsha_wj_patch   (begp:endp))           ; this%psnsha_wj_patch   (:)   = nan
    allocate(this%psnsun_wp_patch   (begp:endp))           ; this%psnsun_wp_patch   (:)   = nan
    allocate(this%psnsha_wp_patch   (begp:endp))           ; this%psnsha_wp_patch   (:)   = nan
    allocate(this%fpsn_patch        (begp:endp))           ; this%fpsn_patch        (:)   = nan
    allocate(this%fpsn_wc_patch     (begp:endp))           ; this%fpsn_wc_patch     (:)   = nan
    allocate(this%fpsn_wj_patch     (begp:endp))           ; this%fpsn_wj_patch     (:)   = nan
    allocate(this%fpsn_wp_patch     (begp:endp))           ; this%fpsn_wp_patch     (:)   = nan

    allocate(this%lmrsun_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsun_z_patch    (:,:) = nan
    allocate(this%lmrsha_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsha_z_patch    (:,:) = nan
    allocate(this%lmrsun_patch      (begp:endp))           ; this%lmrsun_patch      (:)   = nan
    allocate(this%lmrsha_patch      (begp:endp))           ; this%lmrsha_patch      (:)   = nan

    allocate(this%alphapsnsun_patch (begp:endp))           ; this%alphapsnsun_patch (:)   = nan
    allocate(this%alphapsnsha_patch (begp:endp))           ; this%alphapsnsha_patch (:)   = nan
    allocate(this%rc13_canair_patch (begp:endp))           ; this%rc13_canair_patch (:)   = nan
    allocate(this%rc13_psnsun_patch (begp:endp))           ; this%rc13_psnsun_patch (:)   = nan
    allocate(this%rc13_psnsha_patch (begp:endp))           ; this%rc13_psnsha_patch (:)   = nan

    allocate(this%cisun_z_patch     (begp:endp,1:nlevcan)) ; this%cisun_z_patch     (:,:) = nan
    allocate(this%cisha_z_patch     (begp:endp,1:nlevcan)) ; this%cisha_z_patch     (:,:) = nan

    allocate(this%rssun_z_patch     (begp:endp,1:nlevcan)) ; this%rssun_z_patch     (:,:) = nan
    allocate(this%rssha_z_patch     (begp:endp,1:nlevcan)) ; this%rssha_z_patch     (:,:) = nan
    allocate(this%rssun_patch       (begp:endp))           ; this%rssun_patch       (:)   = nan
    allocate(this%rssha_patch       (begp:endp))           ; this%rssha_patch       (:)   = nan

    allocate(this%psncanopy_patch   (begp:endp))           ; this%psncanopy_patch   (:)   = nan

    allocate(this%lmrcanopy_patch   (begp:endp))           ; this%lmrcanopy_patch   (:)   = nan
! plant hydraulics
    allocate(this%ac_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ac_phs_patch (:,:,:) = nan
    allocate(this%aj_phs_patch      (begp:endp,2,1:nlevcan)) ; this%aj_phs_patch (:,:,:) = nan
    allocate(this%ap_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ap_phs_patch (:,:,:) = nan
    allocate(this%ag_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ag_phs_patch (:,:,:) = nan
    allocate(this%an_sun_patch      (begp:endp,1:nlevcan))   ; this%an_sun_patch (:,:)   = nan
    allocate(this%an_sha_patch      (begp:endp,1:nlevcan))   ; this%an_sha_patch (:,:)   = nan
    allocate(this%vcmax_z_phs_patch (begp:endp,2,1:nlevcan)) ; this%vcmax_z_phs_patch (:,:,:) = nan
    allocate(this%tpu_z_phs_patch   (begp:endp,2,1:nlevcan)) ; this%tpu_z_phs_patch   (:,:,:) = nan
    allocate(this%kp_z_phs_patch    (begp:endp,2,1:nlevcan)) ; this%kp_z_phs_patch    (:,:,:) = nan
    allocate(this%gs_mol_sun_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sun_patch  (:,:)   = nan
    allocate(this%gs_mol_sha_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sha_patch  (:,:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%rh_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname='RH_LEAF', units='fraction', &
         avgflag='A', long_name='fractional humidity at leaf surface', &
         ptr_patch=this%rh_leaf_patch, set_spec=spval, default='inactive')

    this%fpsn_patch(begp:endp) = spval
    call hist_addfld1d (fname='FPSN', units='umol/m2s',  &
         avgflag='A', long_name='photosynthesis', &
         ptr_patch=this%fpsn_patch, set_lake=0._r8, set_urb=0._r8)

    this%fpsn_wc_patch(begp:endp) = spval
    call hist_addfld1d (fname='FPSN_WC', units='umol/m2s',  &
         avgflag='A', long_name='Rubisco-limited photosynthesis', &
         ptr_patch=this%fpsn_wc_patch, set_lake=0._r8, set_urb=0._r8)

    this%fpsn_wj_patch(begp:endp) = spval
    call hist_addfld1d (fname='FPSN_WJ', units='umol/m2s',  &
         avgflag='A', long_name='RuBP-limited photosynthesis', &
         ptr_patch=this%fpsn_wj_patch, set_lake=0._r8, set_urb=0._r8)

    this%fpsn_wp_patch(begp:endp) = spval
    call hist_addfld1d (fname='FPSN_WP', units='umol/m2s',  &
         avgflag='A', long_name='Product-limited photosynthesis', &
         ptr_patch=this%fpsn_wp_patch, set_lake=0._r8, set_urb=0._r8)

    if (use_cn) then
       this%psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='sunlit leaf photosynthesis', &
            ptr_patch=this%psnsun_patch)

       this%psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='shaded leaf photosynthesis', &
            ptr_patch=this%psnsha_patch)
    end if

    if ( use_c13 ) then
       this%c13_psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C13 sunlit leaf photosynthesis', &
            ptr_patch=this%c13_psnsun_patch)

       this%c13_psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C13 shaded leaf photosynthesis', &
            ptr_patch=this%c13_psnsha_patch)
    end if

    if ( use_c14 ) then
       this%c14_psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C14 sunlit leaf photosynthesis', &
            ptr_patch=this%c14_psnsun_patch)

       this%c14_psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='C14 shaded leaf photosynthesis', &
            ptr_patch=this%c14_psnsha_patch)
    end if

    if ( use_c13 ) then
       this%rc13_canair_patch(begp:endp) = spval
       call hist_addfld1d (fname='RC13_CANAIR', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for canopy air', &
            ptr_patch=this%rc13_canair_patch)

       this%rc13_psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='RC13_PSNSUN', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for sunlit photosynthesis', &
            ptr_patch=this%rc13_psnsun_patch)

       this%rc13_psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='RC13_PSNSHA', units='proportion', &
            avgflag='A', long_name='C13/C(12+13) for shaded photosynthesis', &
            ptr_patch=this%rc13_psnsha_patch)
    endif

    ! Canopy physiology

    if ( use_c13 ) then
       this%alphapsnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='ALPHAPSNSUN', units='proportion', &
            avgflag='A', long_name='sunlit c13 fractionation', &
            ptr_patch=this%alphapsnsun_patch, default='inactive')
       
       this%alphapsnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='ALPHAPSNSHA', units='proportion', &
            avgflag='A', long_name='shaded c13 fractionation', &
            ptr_patch=this%alphapsnsha_patch, default='inactive')
    endif

    this%rssun_patch(begp:endp) = spval
    call hist_addfld1d (fname='RSSUN', units='s/m',  &
         avgflag='M', long_name='sunlit leaf stomatal resistance', &
         ptr_patch=this%rssun_patch, set_lake=spval, set_urb=spval, default='inactive')

    this%rssha_patch(begp:endp) = spval
    call hist_addfld1d (fname='RSSHA', units='s/m',  &
         avgflag='M', long_name='shaded leaf stomatal resistance', &
         ptr_patch=this%rssha_patch, set_lake=spval, set_urb=spval, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p,l                        ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)

       this%lmrcanopy_patch(p) =  0.0_r8 

       this%alphapsnsun_patch(p) = spval
       this%alphapsnsha_patch(p) = spval

       if (lun_pp%ifspecial(l)) then
          this%psnsun_patch(p) = 0._r8
          this%psnsha_patch(p) = 0._r8
          if ( use_c13 ) then
             this%c13_psnsun_patch(p) = 0._r8
             this%c13_psnsha_patch(p) = 0._r8
          endif
          if ( use_c14 ) then
             this%c14_psnsun_patch(p) = 0._r8
             this%c14_psnsha_patch(p) = 0._r8
          endif
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    if ( use_c13 ) then
       call restartvar(ncid=ncid, flag=flag, varname='rc13_canair', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%rc13_canair_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='rc13_psnsun', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%rc13_psnsun_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='rc13_psnsha', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%rc13_psnsha_patch) 
    endif

  end subroutine Restart

  !------------------------------------------------------------------------------
  subroutine TimeStepInit (this, bounds)
    !
    ! Time step initialization
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice, istice_mec, istwet 
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type) , intent(in)    :: bounds                         
    !
    ! !LOCAL VARIABLES:
    integer :: p,l ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       l = veg_pp%landunit(p)
       if (.not. lun_pp%lakpoi(l)) then
          this%psnsun_patch(p)    = 0._r8
          this%psnsun_wc_patch(p) = 0._r8
          this%psnsun_wj_patch(p) = 0._r8
          this%psnsun_wp_patch(p) = 0._r8

          this%psnsha_patch(p)    = 0._r8
          this%psnsha_wc_patch(p) = 0._r8
          this%psnsha_wj_patch(p) = 0._r8
          this%psnsha_wp_patch(p) = 0._r8

          this%fpsn_patch(p)      = 0._r8
          this%fpsn_wc_patch(p)   = 0._r8
          this%fpsn_wj_patch(p)   = 0._r8
          this%fpsn_wp_patch(p)   = 0._r8

          if ( use_c13 ) then
             this%alphapsnsun_patch(p) = 0._r8
             this%alphapsnsha_patch(p) = 0._r8
             this%c13_psnsun_patch(p)  = 0._r8
             this%c13_psnsha_patch(p)  = 0._r8
          endif
          if ( use_c14 ) then
             this%c14_psnsun_patch(p) = 0._r8
             this%c14_psnsha_patch(p) = 0._r8
          endif
       end if
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop &
            .or. lun_pp%itype(l) == istice .or. lun_pp%itype(l) == istice_mec &
            .or. lun_pp%itype(l) == istwet) then
          if (use_c13) then
             this%rc13_canair_patch(p) = 0._r8
             this%rc13_psnsun_patch(p) = 0._r8
             this%rc13_psnsha_patch(p) = 0._r8
          end if
       end if
    end do

  end subroutine TimeStepInit

  !------------------------------------------------------------------------------
  subroutine NewPatchInit (this, p)
    !
    ! For new run-time pft, modify state and flux variables to maintain 
    ! carbon and nitrogen balance with dynamic pft-weights. 
    ! Called from dyn_cnbal_patch
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    integer, intent(in) :: p
    !-----------------------------------------------------------------------

    if ( use_c13 ) then
       this%alphapsnsun_patch(p) = 0._r8
       this%alphapsnsha_patch(p) = 0._r8
       this%rc13_canair_patch(p) = 0._r8
       this%rc13_psnsun_patch(p) = 0._r8
       this%rc13_psnsha_patch(p) = 0._r8
    endif

    this%psnsun_patch(p) = 0._r8
    this%psnsha_patch(p) = 0._r8

    if (use_c13) then
       this%c13_psnsun_patch(p) = 0._r8
       this%c13_psnsha_patch(p) = 0._r8
    end if
    if ( use_c14 ) then
       this%c14_psnsun_patch(p) = 0._r8
       this%c14_psnsha_patch(p) = 0._r8
    end if
    
  end subroutine NewPatchInit

end module PhotosynthesisType
