Module SoilHydrologyType

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use decompMod             , only : bounds_type
  use abortutils            , only : endrun
  use clm_varpar            , only : nlevgrnd, nlayer, nlayert, nlevsoi 
  use clm_varpar            , only : more_vertlayers, nlevsoifl, toplev_equalspace 
  use clm_varctl            , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  type, public :: soilhydrology_type

     integer :: h2osfcflag              ! true => surface water is active (namelist)       
     integer :: origflag                ! used to control soil hydrology properties (namelist)

     ! NON-VIC
     real(r8), pointer :: frost_table_col   (:)     ! col frost table depth                    
     real(r8), pointer :: zwt_col           (:)     ! col water table depth
     real(r8), pointer :: zwts_col          (:)     ! col water table depth, the shallower of the two water depths     
     real(r8), pointer :: zwt_perched_col   (:)     ! col perched water table depth
     real(r8), pointer :: wa_col            (:)     ! col water in the unconfined aquifer (mm)
     real(r8), pointer :: beg_wa_grc        (:)     ! grid-level water in the unconfined aquifer at beginning of the time step (mm)
     real(r8), pointer :: end_wa_grc        (:)     ! grid-level water in the unconfined aquifer at end of the time step (mm)
     real(r8), pointer :: qflx_bot_col      (:)
     real(r8), pointer :: qcharge_col       (:)     ! col aquifer recharge rate (mm/s) 
     real(r8), pointer :: fracice_col       (:,:)   ! col fractional impermeability (-)
     real(r8), pointer :: icefrac_col       (:,:)   ! col fraction of ice       
     real(r8), pointer :: fcov_col          (:)     ! col fractional impermeable area
     real(r8), pointer :: fsat_col          (:)     ! col fractional area with water table at surface
     real(r8), pointer :: h2osfc_thresh_col (:)     ! col level at which h2osfc "percolates"   (time constant)

     ! VIC 
     real(r8), pointer :: hkdepth_col       (:)     ! col VIC decay factor (m) (time constant)                    
     real(r8), pointer :: b_infil_col       (:)     ! col VIC b infiltration parameter (time constant)                    
     real(r8), pointer :: ds_col            (:)     ! col VIC fracton of Dsmax where non-linear baseflow begins (time constant)                    
     real(r8), pointer :: dsmax_col         (:)     ! col VIC max. velocity of baseflow (mm/day) (time constant)
     real(r8), pointer :: Wsvic_col         (:)     ! col VIC fraction of maximum soil moisutre where non-liear base flow occurs (time constant)
     real(r8), pointer :: porosity_col      (:,:)   ! col VIC porosity (1-bulk_density/soil_density)
     real(r8), pointer :: vic_clm_fract_col (:,:,:) ! col VIC fraction of VIC layers in CLM layers 
     real(r8), pointer :: depth_col         (:,:)   ! col VIC layer depth of upper layer  
     real(r8), pointer :: c_param_col       (:)     ! col VIC baseflow exponent (Qb) 
     real(r8), pointer :: expt_col          (:,:)   ! col VIC pore-size distribution related paramter(Q12) 
     real(r8), pointer :: ksat_col          (:,:)   ! col VIC Saturated hydrologic conductivity 
     real(r8), pointer :: phi_s_col         (:,:)   ! col VIC soil moisture dissusion parameter 
     real(r8), pointer :: moist_col         (:,:)   ! col VIC soil moisture (kg/m2) for VIC soil layers 
     real(r8), pointer :: moist_vol_col     (:,:)   ! col VIC volumetric soil moisture for VIC soil layers 
     real(r8), pointer :: max_moist_col     (:,:)   ! col VIC max layer moist + ice (mm) 
     real(r8), pointer :: max_infil_col     (:)     ! col VIC maximum infiltration rate calculated in VIC
     real(r8), pointer :: i_0_col           (:)     ! col VIC average saturation in top soil layers 
     real(r8), pointer :: ice_col           (:,:)   ! col VIC soil ice (kg/m2) for VIC soil layers

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate

  end type soilhydrology_type
  !-----------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilhydrology_type) :: this
    type(bounds_type), intent(in)    :: bounds  

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%frost_table_col   (begc:endc))                 ; this%frost_table_col   (:)     = nan
    allocate(this%zwt_col           (begc:endc))                 ; this%zwt_col           (:)     = nan
    allocate(this%qflx_bot_col      (begc:endc))                 ; this%qflx_bot_col      (:)     = nan
    allocate(this%zwt_perched_col   (begc:endc))                 ; this%zwt_perched_col   (:)     = nan
    allocate(this%zwts_col          (begc:endc))                 ; this%zwts_col          (:)     = nan

    allocate(this%wa_col            (begc:endc))                 ; this%wa_col            (:)     = nan
    allocate(this%beg_wa_grc        (begg:endg))                 ; this%beg_wa_grc        (:)     = nan
    allocate(this%end_wa_grc        (begg:endg))                 ; this%end_wa_grc        (:)     = nan
    allocate(this%qcharge_col       (begc:endc))                 ; this%qcharge_col       (:)     = nan
    allocate(this%fracice_col       (begc:endc,nlevgrnd))        ; this%fracice_col       (:,:)   = nan
    allocate(this%icefrac_col       (begc:endc,nlevgrnd))        ; this%icefrac_col       (:,:)   = nan
    allocate(this%fcov_col          (begc:endc))                 ; this%fcov_col          (:)     = nan   
    allocate(this%fsat_col          (begc:endc))                 ; this%fsat_col          (:)     = nan
    allocate(this%h2osfc_thresh_col (begc:endc))                 ; this%h2osfc_thresh_col (:)     = nan

    allocate(this%hkdepth_col       (begc:endc))                 ; this%hkdepth_col       (:)     = nan
    allocate(this%b_infil_col       (begc:endc))                 ; this%b_infil_col       (:)     = nan
    allocate(this%ds_col            (begc:endc))                 ; this%ds_col            (:)     = nan
    allocate(this%dsmax_col         (begc:endc))                 ; this%dsmax_col         (:)     = nan
    allocate(this%Wsvic_col         (begc:endc))                 ; this%Wsvic_col         (:)     = nan
    allocate(this%depth_col         (begc:endc,nlayert))         ; this%depth_col         (:,:)   = nan
    allocate(this%porosity_col      (begc:endc,nlayer))          ; this%porosity_col      (:,:)   = nan
    allocate(this%vic_clm_fract_col (begc:endc,nlayer, nlevsoi)) ; this%vic_clm_fract_col (:,:,:) = nan
    allocate(this%c_param_col       (begc:endc))                 ; this%c_param_col       (:)     = nan
    allocate(this%expt_col          (begc:endc,nlayer))          ; this%expt_col          (:,:)   = nan
    allocate(this%ksat_col          (begc:endc,nlayer))          ; this%ksat_col          (:,:)   = nan
    allocate(this%phi_s_col         (begc:endc,nlayer))          ; this%phi_s_col         (:,:)   = nan
    allocate(this%moist_col         (begc:endc,nlayert))         ; this%moist_col         (:,:)   = nan
    allocate(this%moist_vol_col     (begc:endc,nlayert))         ; this%moist_vol_col     (:,:)   = nan
    allocate(this%max_moist_col     (begc:endc,nlayer))          ; this%max_moist_col     (:,:)   = nan
    allocate(this%max_infil_col     (begc:endc))                 ; this%max_infil_col     (:)     = nan
    allocate(this%i_0_col           (begc:endc))                 ; this%i_0_col           (:)     = nan
    allocate(this%ice_col           (begc:endc,nlayert))         ; this%ice_col           (:,:)   = nan

  end subroutine InitAllocate

 end Module SoilHydrologyType
