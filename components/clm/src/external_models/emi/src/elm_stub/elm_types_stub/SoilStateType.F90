module SoilStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varpar      , only : more_vertlayers, numpft, numrad
  use clm_varpar      , only : nlevsoi, nlevgrnd, nlevlak, nlevsoifl, nlayer, nlayert, nlevurb, nlevsno
  use landunit_varcon , only : istice, istdlak, istwet, istsoil, istcrop, istice_mec
  use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv 
  use clm_varcon      , only : secspday, pc, mu, denh2o, denice, grlnd
  use clm_varctl      , only : use_cn, use_lch4,use_dynroot, use_fates
  use clm_varctl      , only : use_var_soil_thick
  use clm_varctl      , only : iulog, fsurdat, hist_wrtch4diag
  use clm_varcon      , only : spval
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: soilstate_type

     ! sand/ clay/ organic matter
     real(r8), pointer :: sandfrac_patch       (:)   ! patch sand fraction
     real(r8), pointer :: clayfrac_patch       (:)   ! patch clay fraction
     real(r8), pointer :: mss_frc_cly_vld_col  (:)   ! col mass fraction clay limited to 0.20
     real(r8), pointer :: cellorg_col          (:,:) ! col organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellsand_col         (:,:) ! sand value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellclay_col         (:,:) ! clay value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: bd_col               (:,:) ! col bulk density of dry soil material [kg/m^3] (CN)

     ! hydraulic properties
     real(r8), pointer :: hksat_col            (:,:) ! col hydraulic conductivity at saturation (mm H2O /s) 
     real(r8), pointer :: hksat_min_col        (:,:) ! col mineral hydraulic conductivity at saturation (hksat) (mm/s)
     real(r8), pointer :: hk_l_col             (:,:) ! col hydraulic conductivity (mm/s)
     real(r8), pointer :: smp_l_col            (:,:) ! col soil matric potential (mm)
     real(r8), pointer :: smpmin_col           (:)   ! col restriction for min of soil potential (mm) 
     real(r8), pointer :: bsw_col              (:,:) ! col Clapp and Hornberger "b" (nlevgrnd)  
     real(r8), pointer :: watsat_col           (:,:) ! col volumetric soil water at saturation (porosity) 
     real(r8), pointer :: watdry_col           (:,:) ! col btran parameter for btran = 0
     real(r8), pointer :: watopt_col           (:,:) ! col btran parameter for btran = 1
     real(r8), pointer :: watfc_col            (:,:) ! col volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: watmin_col           (:,:) ! col minimum volumetric soil water (nlevsoi)
     real(r8), pointer :: sucsat_col           (:,:) ! col minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: sucmin_col           (:,:) ! col minimum allowable soil liquid suction pressure (mm) [Note: sucmin_col is a negative value, while sucsat_col is a positive quantity]
     real(r8), pointer :: soilbeta_col         (:)   ! col factor that reduces ground evaporation L&P1992(-)
     real(r8), pointer :: soilalpha_col        (:)   ! col factor that reduces ground saturated specific humidity (-)
     real(r8), pointer :: soilalpha_u_col      (:)   ! col urban factor that reduces ground saturated specific humidity (-) 
     real(r8), pointer :: soilpsi_col          (:,:) ! col soil water potential in each soil layer (MPa) (CN)
     real(r8), pointer :: wtfact_col           (:)   ! col maximum saturated fraction for a gridcell
     real(r8), pointer :: porosity_col         (:,:) ! col soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity_col     (:,:) ! col effective porosity = porosity - vol_ice (nlevgrnd) 
     real(r8), pointer :: gwc_thr_col          (:)   ! col threshold soil moisture based on clay content

     ! thermal conductivity / heat capacity
     real(r8), pointer :: thk_col              (:,:) ! col thermal conductivity of each layer [W/m-K] 
     real(r8), pointer :: tkmg_col             (:,:) ! col thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: tkdry_col            (:,:) ! col thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd) 
     real(r8), pointer :: tksatu_col           (:,:) ! col thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: csol_col             (:,:) ! col heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd) 

     ! roots
     real(r8), pointer :: rootr_patch          (:,:) ! patch effective fraction of roots in each soil layer (nlevgrnd)
     real(r8), pointer :: rootr_col            (:,:) ! col effective fraction of roots in each soil layer (nlevgrnd)  
     real(r8), pointer :: rootfr_col           (:,:) ! col fraction of roots in each soil layer (nlevgrnd) 
     real(r8), pointer :: rootfr_patch         (:,:) ! patch fraction of roots in each soil layer (nlevgrnd)
     real(r8), pointer :: rootr_road_perv_col  (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: rootfr_road_perv_col (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: root_depth_patch     (:)   ! rooting depth of each PFT (m)

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate

  end type soilstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: begc_all, endc_all
    !------------------------------------------------------------------------

    begp     = bounds%begp    ; endp     = bounds%endp
    begc     = bounds%begc    ; endc     = bounds%endc
    begg     = bounds%begg    ; endg     = bounds%endg
    begc_all = bounds%begc    ; endc_all = bounds%endc

    allocate(this%mss_frc_cly_vld_col  (begc:endc))                     ; this%mss_frc_cly_vld_col  (:)   = nan
    allocate(this%sandfrac_patch       (begp:endp))                     ; this%sandfrac_patch       (:)   = nan
    allocate(this%clayfrac_patch       (begp:endp))                     ; this%clayfrac_patch       (:)   = nan
    allocate(this%cellorg_col          (begc:endc,nlevgrnd))            ; this%cellorg_col          (:,:) = nan 
    allocate(this%cellsand_col         (begc:endc,nlevgrnd))            ; this%cellsand_col         (:,:) = nan 
    allocate(this%cellclay_col         (begc:endc,nlevgrnd))            ; this%cellclay_col         (:,:) = nan 
    allocate(this%bd_col               (begc:endc,nlevgrnd))            ; this%bd_col               (:,:) = nan

    allocate(this%hksat_col            (begc_all:endc_all,nlevgrnd))    ; this%hksat_col            (:,:) = spval
    allocate(this%hksat_min_col        (begc:endc,nlevgrnd))            ; this%hksat_min_col        (:,:) = spval
    allocate(this%hk_l_col             (begc:endc,nlevgrnd))            ; this%hk_l_col             (:,:) = nan   
    allocate(this%smp_l_col            (begc:endc,nlevgrnd))            ; this%smp_l_col            (:,:) = nan   
    allocate(this%smpmin_col           (begc:endc))                     ; this%smpmin_col           (:)   = nan

    allocate(this%bsw_col              (begc_all:endc_all,nlevgrnd))    ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col           (begc_all:endc_all,nlevgrnd))    ; this%watsat_col           (:,:) = nan
    allocate(this%watdry_col           (begc:endc,nlevgrnd))            ; this%watdry_col           (:,:) = spval
    allocate(this%watopt_col           (begc:endc,nlevgrnd))            ; this%watopt_col           (:,:) = spval
    allocate(this%watfc_col            (begc:endc,nlevgrnd))            ; this%watfc_col            (:,:) = nan
    allocate(this%watmin_col           (begc:endc,nlevgrnd))            ; this%watmin_col           (:,:) = nan
    allocate(this%sucsat_col           (begc:endc,nlevgrnd))            ; this%sucsat_col           (:,:) = spval
    allocate(this%sucmin_col           (begc:endc,nlevgrnd))            ; this%sucmin_col           (:,:) = spval
    allocate(this%soilbeta_col         (begc:endc))                     ; this%soilbeta_col         (:)   = nan   
    allocate(this%soilalpha_col        (begc:endc))                     ; this%soilalpha_col        (:)   = nan
    allocate(this%soilalpha_u_col      (begc:endc))                     ; this%soilalpha_u_col      (:)   = nan
    allocate(this%soilpsi_col          (begc:endc,nlevgrnd))            ; this%soilpsi_col          (:,:) = nan
    allocate(this%wtfact_col           (begc:endc))                     ; this%wtfact_col           (:)   = nan
    allocate(this%porosity_col         (begc:endc,nlayer))              ; this%porosity_col         (:,:) = spval
    allocate(this%eff_porosity_col     (begc:endc,nlevgrnd))            ; this%eff_porosity_col     (:,:) = spval
    allocate(this%gwc_thr_col          (begc:endc))                     ; this%gwc_thr_col          (:)   = nan

    allocate(this%thk_col              (begc:endc,-nlevsno+1:nlevgrnd)) ; this%thk_col              (:,:) = nan
    allocate(this%tkmg_col             (begc:endc,nlevgrnd))            ; this%tkmg_col             (:,:) = nan
    allocate(this%tkdry_col            (begc:endc,nlevgrnd))            ; this%tkdry_col            (:,:) = nan
    allocate(this%tksatu_col           (begc:endc,nlevgrnd))            ; this%tksatu_col           (:,:) = nan
    allocate(this%csol_col             (begc:endc,nlevgrnd))            ; this%csol_col             (:,:) = nan

    allocate(this%rootr_patch          (begp:endp,1:nlevgrnd))          ; this%rootr_patch          (:,:) = nan
    allocate(this%rootr_col            (begc:endc,nlevgrnd))            ; this%rootr_col            (:,:) = nan
    allocate(this%rootr_road_perv_col  (begc:endc,1:nlevgrnd))          ; this%rootr_road_perv_col  (:,:) = nan
    allocate(this%rootfr_patch         (begp:endp,1:nlevgrnd))          ; this%rootfr_patch         (:,:) = nan
    allocate(this%rootfr_col           (begc:endc,1:nlevgrnd))          ; this%rootfr_col           (:,:) = nan 
    allocate(this%rootfr_road_perv_col (begc:endc,1:nlevgrnd))          ; this%rootfr_road_perv_col (:,:) = nan
    allocate(this%root_depth_patch     (begp:endp))                     ; this%root_depth_patch     (:)   = spval

  end subroutine InitAllocate


end module SoilStateType
