module ColumnWaterFluxType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use decompMod    , only : bounds_type, get_proc_global
  use clm_varcon   , only : spval
  use LandunitType , only : lun                
  use ColumnType   , only : col                
  use PatchType    , only : pft                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: soilcol_water_flux

      real(r8), pointer :: qflx_prec_grnd_col       (:)   ! col water onto ground including canopy runoff [kg/(m2 s)]
      real(r8), pointer :: qflx_rain_grnd_col       (:)   ! col rain on ground after interception (mm H2O/s) [+]
      real(r8), pointer :: qflx_snow_grnd_col       (:)   ! col snow on ground after interception (mm H2O/s) [+]
      real(r8), pointer :: qflx_sub_snow_col        (:)   ! col sublimation rate from snow pack (mm H2O /s) [+]
      real(r8), pointer :: qflx_evap_soi_col        (:)   ! col soil evaporation (mm H2O/s) (+ = to atm)
      real(r8), pointer :: qflx_evap_veg_col        (:)   ! col vegetation evaporation (mm H2O/s) (+ = to atm)
      real(r8), pointer :: qflx_evap_can_col        (:)   ! col evaporation from leaves and stems (mm H2O/s) (+ = to atm)
      real(r8), pointer :: qflx_evap_tot_col        (:)   ! col col_qflx_evap_soi + col_qflx_evap_veg + qflx_tran_veg
      real(r8), pointer :: qflx_evap_grnd_col       (:)   ! col ground surface evaporation rate (mm H2O/s) [+]
      real(r8), pointer :: qflx_snwcp_liq_col       (:)   ! col excess rainfall due to snow capping (mm H2O /s)
      real(r8), pointer :: qflx_snwcp_ice_col       (:)   ! col excess snowfall due to snow capping (mm H2O /s)
      real(r8), pointer :: qflx_tran_veg_col        (:)   ! col vegetation transpiration (mm H2O/s) (+ = to atm)
      real(r8), pointer :: qflx_dew_snow_col        (:)   ! col surface dew added to snow pack (mm H2O /s) [+]
      real(r8), pointer :: qflx_dew_grnd_col        (:)   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
      real(r8), pointer :: qflx_prec_intr_col       (:)   ! col interception of precipitation [mm/s]

      real(r8), pointer :: qflx_ev_snow_col         (:)   ! col evaporation heat flux from snow         (W/m**2) [+ to atm]
      real(r8), pointer :: qflx_ev_soil_col         (:)   ! col evaporation heat flux from soil         (W/m**2) [+ to atm]
      real(r8), pointer :: qflx_ev_h2osfc_col       (:)   ! col evaporation heat flux from soil         (W/m**2) [+ to atm]

      real(r8), pointer :: qflx_gross_evap_soil_col (:)   ! col gross infiltration from soil, this satisfies the relationship qflx_infl_col = qflx_gross_infl_soil_col-qflx_gross_evap_soil_col
      real(r8), pointer :: qflx_gross_infl_soil_col (:)   ! col gross infiltration, before considering the evaporation
      real(r8), pointer :: qflx_adv_col             (:,:) ! col advective flux across different soil layer interfaces [mm H2O/s] [+ downward]
      real(r8), pointer :: qflx_rootsoi_col         (:,:) ! col root and soil water exchange [mm H2O/s] [+ into root]     

      real(r8), pointer :: qflx_infl_col            (:)   ! col infiltration (mm H2O /s)
      real(r8), pointer :: qflx_surf_col            (:)   ! col surface runoff (mm H2O /s)
      real(r8), pointer :: qflx_drain_col           (:)   ! col sub-surface runoff (mm H2O /s)
      real(r8), pointer :: qflx_top_soil_col        (:)   ! col net water input into soil from top (mm/s)
      real(r8), pointer :: qflx_h2osfc_to_ice_col   (:)   ! col conversion of h2osfc to ice
      real(r8), pointer :: qflx_h2osfc_surf_col     (:)   ! col surface water runoff
      real(r8), pointer :: qflx_snow_h2osfc_col     (:)   ! col snow falling on surface water
      real(r8), pointer :: qflx_drain_perched_col   (:)   ! col sub-surface runoff from perched wt (mm H2O /s)
      real(r8), pointer :: qflx_deficit_col         (:)   ! col water deficit to keep non-negative liquid water content (mm H2O)   
      real(r8), pointer :: qflx_floodc_col          (:)   ! col flood water flux at column level
      real(r8), pointer :: qflx_sl_top_soil_col     (:)   ! col liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
      real(r8), pointer :: qflx_snomelt_col         (:)   ! col snow melt (mm H2O /s)
      real(r8), pointer :: qflx_snow_melt_col       (:)   ! col snow melt (net)
      real(r8), pointer :: qflx_qrgwl_col           (:)   ! col qflx_surf at glaciers, wetlands, lakes
      real(r8), pointer :: qflx_runoff_col          (:)   ! col total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
      real(r8), pointer :: qflx_runoff_r_col        (:)   ! col Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
      real(r8), pointer :: qflx_runoff_u_col        (:)   ! col urban total runoff (qflx_drain+qflx_surf) (mm H2O /s) 
      real(r8), pointer :: qflx_rsub_sat_col        (:)   ! col soil saturation excess [mm/s]
      real(r8), pointer :: qflx_snofrz_lyr_col      (:,:) ! col snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
      real(r8), pointer :: qflx_snofrz_col          (:)   ! col column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]
      real(r8), pointer :: qflx_glcice_col          (:)   ! col net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC
      real(r8), pointer :: qflx_glcice_frz_col      (:)   ! col ice growth (positive definite) (mm H2O/s)
      real(r8), pointer :: qflx_glcice_melt_col     (:)   ! col ice melt (positive definite) (mm H2O/s)
      real(r8), pointer :: qflx_drain_vr_col        (:,:) ! col liquid water losted as drainage (m /time step)
      real(r8), pointer :: qflx_h2osfc2topsoi_col   (:)   ! col liquid water coming from surface standing water top soil (mm H2O/s)
      real(r8), pointer :: qflx_snow2topsoi_col     (:)   ! col liquid water coming from residual snow to topsoil (mm H2O/s)

      real(r8), pointer :: snow_sources_col         (:)   ! col snow sources (mm H2O/s)
      real(r8), pointer :: snow_sinks_col           (:)   ! col snow sinks (mm H2O/s)

      real(r8), pointer :: qflx_irrig_col           (:)   ! col irrigation flux (mm H2O/s)

      real(r8), pointer :: mflx_infl_col_1d         (:)   ! infiltration source in top soil control volume (kg H2O /s)
      real(r8), pointer :: mflx_dew_col_1d          (:)   ! liquid+snow dew source in top soil control volume (kg H2O /s)
      real(r8), pointer :: mflx_et_col_1d           (:)   ! evapotranspiration sink from all soil coontrol volumes (kg H2O /s)
      real(r8), pointer :: mflx_drain_col_1d        (:)   ! drainage from groundwater table (kg H2O /s)
      real(r8), pointer :: mflx_drain_perched_col_1d(:)   ! drainage from perched water table (kg H2O /s)
      real(r8), pointer :: mflx_snowlyr_col_1d      (:)   ! mass flux to top soil layer due to disappearance of snow (kg H2O /s)
      real(r8), pointer :: mflx_sub_snow_col_1d     (:)   ! mass flux from top soil layer due to sublimation of snow (kg H2O /s)
      real(r8), pointer :: mflx_snowlyr_col         (:)   ! mass flux to top soil layer due to disappearance of snow (kg H2O /s). This is for restart
      real(r8), pointer :: mflx_neg_snow_col_1d     (:)   ! mass flux from top soil layer due to negative water content in snow layers (kg H2O /s)

  contains
      procedure, public :: Init => init_col_wf
      procedure, public :: InitAllocate => initallocate_col_wf
      procedure, public :: InitCold => initcold_col_wf
      procedure, public :: InitHistory => inithistory_col_wf
      procedure, public :: Restart => restart_col_wf
      procedure, public :: Reset => reset_col_wf
      procedure, public :: Clean => clean_col_wf

  end type soilcol_water_flux

  ! Declare public API
  type(soilcol_water_flux), public, target :: col_wf

  subroutine init_col_wf(this, bounds)

    class(soilcol_water_flux) :: this
    type(bounds_type), intent(in)    :: bounds  

    call this%InitAllocate(bounds) ! same as "call initAllocate_type(hydro, bounds)"
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init


  subroutine initallocate_col_wf(this, begc, endc)
    class(soilcol_water_flux) :: this
    integer, intent(in) :: begc   ! beginning soil column index
    integer, intent(in) :: endc   ! ending soil column index    

    allocate(this%qflx_prec_intr_col       (begc:endc))              ; this%qflx_prec_intr_col       (:)   = nan
    allocate(this%qflx_prec_grnd_col       (begc:endc))              ; this%qflx_prec_grnd_col       (:)   = nan
    allocate(this%qflx_rain_grnd_col       (begc:endc))              ; this%qflx_rain_grnd_col       (:)   = nan
    allocate(this%qflx_snow_grnd_col       (begc:endc))              ; this%qflx_snow_grnd_col       (:)   = nan
    allocate(this%qflx_sub_snow_col        (begc:endc))              ; this%qflx_sub_snow_col        (:)   = 0.0_r8
    allocate(this%qflx_snwcp_liq_col       (begc:endc))              ; this%qflx_snwcp_liq_col       (:)   = nan
    allocate(this%qflx_snwcp_ice_col       (begc:endc))              ; this%qflx_snwcp_ice_col       (:)   = nan
    allocate(this%qflx_tran_veg_col        (begc:endc))              ; this%qflx_tran_veg_col        (:)   = nan
    allocate(this%qflx_evap_veg_col        (begc:endc))              ; this%qflx_evap_veg_col        (:)   = nan
    allocate(this%qflx_evap_can_col        (begc:endc))              ; this%qflx_evap_can_col        (:)   = nan
    allocate(this%qflx_evap_soi_col        (begc:endc))              ; this%qflx_evap_soi_col        (:)   = nan
    allocate(this%qflx_evap_tot_col        (begc:endc))              ; this%qflx_evap_tot_col        (:)   = nan
    allocate(this%qflx_evap_grnd_col       (begc:endc))              ; this%qflx_evap_grnd_col       (:)   = nan
    allocate(this%qflx_dew_grnd_col        (begc:endc))              ; this%qflx_dew_grnd_col        (:)   = nan
    allocate(this%qflx_dew_snow_col        (begc:endc))              ; this%qflx_dew_snow_col        (:)   = nan

    allocate( this%qflx_ev_snow_col        (begc:endc))              ; this%qflx_ev_snow_col         (:)   = nan
    allocate( this%qflx_ev_soil_col        (begc:endc))              ; this%qflx_ev_soil_col         (:)   = nan
    allocate( this%qflx_ev_h2osfc_col      (begc:endc))              ; this%qflx_ev_h2osfc_col       (:)   = nan

    allocate(this%qflx_gross_evap_soil_col (begc:endc))              ; this%qflx_gross_evap_soil_col (:)   = nan
    allocate(this%qflx_gross_infl_soil_col (begc:endc))              ; this%qflx_gross_infl_soil_col (:)   = nan
    allocate(this%qflx_drain_vr_col        (begc:endc,1:nlevsoi))    ; this%qflx_drain_vr_col        (:,:) = nan
    allocate(this%qflx_adv_col             (begc:endc,0:nlevsoi))    ; this%qflx_adv_col             (:,:) = nan
    allocate(this%qflx_rootsoi_col         (begc:endc,1:nlevsoi))    ; this%qflx_rootsoi_col         (:,:) = nan
    
    allocate(this%qflx_infl_col            (begc:endc))              ; this%qflx_infl_col            (:)   = nan
    allocate(this%qflx_surf_col            (begc:endc))              ; this%qflx_surf_col            (:)   = nan
    allocate(this%qflx_drain_col           (begc:endc))              ; this%qflx_drain_col           (:)   = nan
    allocate(this%qflx_top_soil_col        (begc:endc))              ; this%qflx_top_soil_col        (:)   = nan
    allocate(this%qflx_h2osfc_to_ice_col   (begc:endc))              ; this%qflx_h2osfc_to_ice_col   (:)   = nan
    allocate(this%qflx_h2osfc_surf_col     (begc:endc))              ; this%qflx_h2osfc_surf_col     (:)   = nan
    allocate(this%qflx_snow_h2osfc_col     (begc:endc))              ; this%qflx_snow_h2osfc_col     (:)   = nan
    allocate(this%qflx_snomelt_col         (begc:endc))              ; this%qflx_snomelt_col         (:)   = nan
    allocate(this%qflx_snow_melt_col       (begc:endc))              ; this%qflx_snow_melt_col       (:)   = nan
    allocate(this%qflx_snofrz_col          (begc:endc))              ; this%qflx_snofrz_col          (:)   = nan
    allocate(this%qflx_snofrz_lyr_col      (begc:endc,-nlevsno+1:0)) ; this%qflx_snofrz_lyr_col      (:,:) = nan
    allocate(this%qflx_qrgwl_col           (begc:endc))              ; this%qflx_qrgwl_col           (:)   = nan
    allocate(this%qflx_drain_perched_col   (begc:endc))              ; this%qflx_drain_perched_col   (:)   = nan
    allocate(this%qflx_deficit_col         (begc:endc))              ; this%qflx_deficit_col         (:)   = nan
    allocate(this%qflx_floodc_col          (begc:endc))              ; this%qflx_floodc_col          (:)   = nan
    allocate(this%qflx_sl_top_soil_col     (begc:endc))              ; this%qflx_sl_top_soil_col     (:)   = nan
    allocate(this%qflx_runoff_col          (begc:endc))              ; this%qflx_runoff_col          (:)   = nan
    allocate(this%qflx_runoff_r_col        (begc:endc))              ; this%qflx_runoff_r_col        (:)   = nan
    allocate(this%qflx_runoff_u_col        (begc:endc))              ; this%qflx_runoff_u_col        (:)   = nan
    allocate(this%qflx_rsub_sat_col        (begc:endc))              ; this%qflx_rsub_sat_col        (:)   = nan
    allocate(this%qflx_glcice_col          (begc:endc))              ; this%qflx_glcice_col          (:)   = nan
    allocate(this%qflx_glcice_frz_col      (begc:endc))              ; this%qflx_glcice_frz_col      (:)   = nan
    allocate(this%qflx_glcice_melt_col     (begc:endc))              ; this%qflx_glcice_melt_col     (:)   = nan
    allocate(this%snow_sources_col         (begc:endc))              ; this%snow_sources_col         (:)   = nan   
    allocate(this%snow_sinks_col           (begc:endc))              ; this%snow_sinks_col           (:)   = nan   


    allocate(this%qflx_irrig_col           (begc:endc))              ; this%qflx_irrig_col           (:)   = nan

    allocate(this%qflx_snow2topsoi_col     (begc:endc))              ; this%qflx_snow2topsoi_col     (:)   = nan
    allocate(this%qflx_h2osfc2topsoi_col   (begc:endc))              ; this%qflx_h2osfc2topsoi_col   (:)   = nan
    
    ncells = endc - begc + 1
    allocate(this%mflx_infl_col_1d(            ncells))              ; this%mflx_infl_col_1d         (:)   = nan
    allocate(this%mflx_dew_col_1d(             ncells))              ; this%mflx_dew_col_1d          (:)   = nan
    allocate(this%mflx_snowlyr_col_1d(         ncells))              ; this%mflx_snowlyr_col_1d      (:)   = nan
    allocate(this%mflx_sub_snow_col_1d(        ncells))              ; this%mflx_sub_snow_col_1d     (:)   = nan
    allocate(this%mflx_snowlyr_col(         begc:endc))              ; this%mflx_snowlyr_col         (:)   = 0._r8
    allocate(this%mflx_neg_snow_col_1d(        ncells))              ; this%mflx_neg_snow_col_1d     (:)   = nan

    ncells = (endc - begc + 1)*nlevgrnd
    allocate(this%mflx_et_col_1d(              ncells))              ; this%mflx_et_col_1d           (:)   = nan
    allocate(this%mflx_drain_col_1d(           ncells))              ; this%mflx_drain_col_1d        (:)   = nan
    allocate(this%mflx_drain_perched_col_1d(   ncells))              ; this%mflx_drain_perched_col_1d(:)  = nan

  end subroutine init_col_wf


  !----------------------------------------------------------------------
  subroutine inithistory_col_wf(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : create_glacier_mec_landunit, use_cn, use_lch4
    use clm_varpar     , only : nlevsno, crop_prog, nlevsoi 
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(soilcol_water_flux) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    this%qflx_top_soil_col(begc:endc) = spval
    call hist_addfld1d (fname='QTOPSOIL',  units='mm/s',  &
         avgflag='A', long_name='water input to surface', &
         ptr_col=this%qflx_top_soil_col, c2l_scale_type='urbanf', default='inactive')

    this%qflx_infl_col(begc:endc) = spval
    call hist_addfld1d (fname='QINFL',  units='mm/s',  &
         avgflag='A', long_name='infiltration', &
         ptr_col=this%qflx_infl_col, c2l_scale_type='urbanf')

    this%qflx_surf_col(begc:endc) = spval
    call hist_addfld1d (fname='QOVER',  units='mm/s',  &
         avgflag='A', long_name='surface runoff', &
         ptr_col=this%qflx_surf_col, c2l_scale_type='urbanf')

    this%qflx_qrgwl_col(begc:endc) = spval
    call hist_addfld1d (fname='QRGWL',  units='mm/s',  &
         avgflag='A', long_name='surface runoff at glaciers (liquid only), wetlands, lakes', &
         ptr_col=this%qflx_qrgwl_col, c2l_scale_type='urbanf')

    this%dwb_col(begc:endc) = spval
    call hist_addfld1d (fname='DWB',  units='mm/s',  &
         avgflag='A', long_name='net change in total water mass', &
         ptr_col=this%dwb_col, c2l_scale_type='urbanf')

    this%qflx_drain_col(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI',  units='mm/s',  &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_col=this%qflx_drain_col, c2l_scale_type='urbanf')

    this%qflx_runoff_col(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_NODYNLNDUSE',  units='mm/s',  &
         avgflag='A', &
         long_name='total liquid runoff (does not include QSNWCPICE) not including correction for land use change', &
         ptr_col=this%qflx_runoff_col, c2l_scale_type='urbanf')

    this%qflx_runoff_u_col(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_U', units='mm/s',  &
         avgflag='A', long_name='Urban total runoff', &
         ptr_col=this%qflx_runoff_u_col, set_nourb=spval, c2l_scale_type='urbanf')

    this%qflx_runoff_r_col(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_R', units='mm/s',  &
         avgflag='A', long_name='Rural total runoff', &
         ptr_col=this%qflx_runoff_r_col, set_spec=spval)

    this%qflx_snow_melt_col(begc:endc) = spval
    call hist_addfld1d (fname='QSNOMELT',  units='mm/s',  &
         avgflag='A', long_name='snow melt', &
         ptr_col=this%qflx_snow_melt_col, c2l_scale_type='urbanf')

    this%qflx_snofrz_col(begc:endc) = spval
    call hist_addfld1d (fname='QSNOFRZ', units='kg/m2/s', &
         avgflag='A', long_name='column-integrated snow freezing rate', &
         ptr_col=this%qflx_snofrz_col, set_lake=spval, c2l_scale_type='urbanf', &
         default='inactive')

    if (create_glacier_mec_landunit) then
       this%qflx_glcice_frz_col(begc:endc) = spval
       call hist_addfld1d (fname='QICE_FRZ',  units='mm/s',  &
            avgflag='A', long_name='ice growth', &
            ptr_col=this%qflx_glcice_frz_col, l2g_scale_type='ice')
    end if

    if (create_glacier_mec_landunit) then
       this%qflx_glcice_melt_col(begc:endc) = spval
       call hist_addfld1d (fname='QICE_MELT',  units='mm/s',  &
            avgflag='A', long_name='ice melt', &
            ptr_col=this%qflx_glcice_melt_col, l2g_scale_type='ice')

    ! Use qflx_snwcp_ice_col rather than qflx_snwcp_ice_patch, because the column version 
    ! is the final version, which includes some  additional corrections beyond the patch-level version
    end if

    this%qflx_h2osfc_surf_col(begc:endc) = spval
    call hist_addfld1d (fname='QH2OSFC',  units='mm/s',  &
         avgflag='A', long_name='surface water runoff', &
         ptr_col=this%qflx_h2osfc_surf_col)

    this%qflx_drain_perched_col(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI_PERCH',  units='mm/s',  &
         avgflag='A', long_name='perched wt drainage', &
         ptr_col=this%qflx_drain_perched_col, c2l_scale_type='urbanf')

    this%qflx_rsub_sat_col(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI_XS',  units='mm/s',  &
         avgflag='A', long_name='saturation excess drainage', &
         ptr_col=this%qflx_rsub_sat_col, c2l_scale_type='urbanf')

    ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at any
    ! given time step but only if there is at least one snow layer (for all landunits 
    ! except lakes).  Also note that monthly average files of snow_sources and snow_sinks
    ! sinks must be weighted by number of days in the month to diagnose, for example, an 
    ! annual value of the change in h2osno. 

    this%snow_sources_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOW_SOURCES',  units='mm/s',  &
         avgflag='A', long_name='snow sources (liquid water)', &
         ptr_col=this%snow_sources_col, c2l_scale_type='urbanf')

    this%snow_sinks_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOW_SINKS',  units='mm/s',  &
         avgflag='A', long_name='snow sinks (liquid water)', &
         ptr_col=this%snow_sinks_col, c2l_scale_type='urbanf')

  end subroutine inithistory_col_wf

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  subroutine initcold_col_wf(this, bounds)
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(soilcol_water_flux) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    !-----------------------------------------------------------------------

    this%qflx_evap_grnd_col(bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_grnd_col (bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_snow_col (bounds%begc:bounds%endc) = 0.0_r8

    this%qflx_h2osfc_surf_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_snow_melt_col(bounds%begc:bounds%endc)   = 0._r8

    this%dwb_col(bounds%begc:bounds%endc) = 0._r8
    ! needed for CNNLeaching 
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%qflx_drain_col(c) = 0._r8
          this%qflx_surf_col(c)  = 0._r8
       end if
    end do


  end subroutine initcold_col_wf



  !------------------------------------------------------------------------
  subroutine restart_col_wf(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : denice, denh2o, pondmx, watmin
    use landunit_varcon  , only : istcrop, istdlak, istsoil 
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varpar       , only : nlevgrnd, nlevurb, nlevsno   
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(soilcol_water_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: dimlen       ! dimension length
    integer :: nump_global  ! total number of pfts, globally
    integer :: err_code     ! error code
    logical :: do_io
    !-----------------------------------------------------------------------

    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    ! needed for SNICAR
    call restartvar(ncid=ncid, flag=flag, varname='qflx_snofrz_lyr', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer ice freezing rate', units='kg m-2 s-1', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snofrz_lyr_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snofrz_lyr to zero
       this%qflx_snofrz_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    call restartvar(ncid=ncid, flag=flag, varname='qflx_snow_melt', xtype=ncd_double,  &
         dim1name='column', &
         long_name='net snow melt', units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snow_melt_col)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snow_melt to zero
       this%qflx_snow_melt_col(bounds%begc:bounds%endc) = 0._r8
    endif

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'n_irrig_steps_left', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'irrig_rate', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if

    call restartvar(ncid=ncid, flag=flag, varname='MFLX_SNOW_LYR', xtype=ncd_double,  &
         dim1name='column', &
         long_name='mass flux due to disapperance of last snow layer', units='kg/s', &
         interpinic_flag='interp', readvar=readvar, data=this%mflx_snowlyr_col)

  end subroutine restart_col_wf


  subroutine reset_col_wf(this, bounds, numf, filter)
  !
  ! !DESCRIPTION:
  ! Intitialize SNICAR variables for fresh snow column
  !
  ! !ARGUMENTS:
  class(waterflux_type)              :: this
  type(bounds_type)    , intent(in)  :: bounds
  integer              , intent(in)  :: numf
  integer              , intent(in)  :: filter(:)
  !-----------------------------------------------------------------------
  
  integer :: fc, column
  
  do fc = 1, numf
    column = filter(fc)
    this%qflx_snow2topsoi_col     (column)   = 0._r8
    this%qflx_h2osfc2topsoi_col   (column)   = 0._r8      
  enddo  
end subroutine reset_col_wf


  subroutine clean_col_wf(this)

    class(soilcol_water_flux) :: this

    deallocate(this%qflx_prec_intr_col )       
    deallocate(this%qflx_prec_grnd_col )      
    deallocate(this%qflx_rain_grnd_col )      
    deallocate(this%qflx_snow_grnd_col )      
    deallocate(this%qflx_sub_snow_col  )      
    deallocate(this%qflx_snwcp_liq_col )      
    deallocate(this%qflx_snwcp_ice_col )      
    deallocate(this%qflx_tran_veg_col  )      
    deallocate(this%qflx_evap_veg_col  )      
    deallocate(this%qflx_evap_can_col  )      
    deallocate(this%qflx_evap_soi_col  )      
    deallocate(this%qflx_evap_tot_col  )      
    deallocate(this%qflx_evap_grnd_col )      
    deallocate(this%qflx_dew_grnd_col  )      
    deallocate(this%qflx_dew_snow_col  )      

    deallocate( this%qflx_ev_snow_col  )
    deallocate( this%qflx_ev_soil_col  )
    deallocate( this%qflx_ev_h2osfc_col)

    deallocate(this%qflx_gross_evap_soil_col 
    deallocate(this%qflx_gross_infl_soil_col 
    deallocate(this%qflx_drain_vr_col        
    deallocate(this%qflx_adv_col             
    deallocate(this%qflx_rootsoi_col         
    
    deallocate(this%qflx_surf_col          )
    deallocate(this%qflx_infl_col          )
    deallocate(this%qflx_drain_col         )
    deallocate(this%qflx_top_soil_col      )
    deallocate(this%qflx_h2osfc_to_ice_col )
    deallocate(this%qflx_h2osfc_surf_col   )
    deallocate(this%qflx_snow_h2osfc_col   )
    deallocate(this%qflx_snomelt_col       )
    deallocate(this%qflx_snow_melt_col     )
    deallocate(this%qflx_snofrz_col        )
    deallocate(this%qflx_snofrz_lyr_col    )
    deallocate(this%qflx_qrgwl_col         )
    deallocate(this%qflx_drain_perched_col )
    deallocate(this%qflx_deficit_col       )
    deallocate(this%qflx_floodc_col        )
    deallocate(this%qflx_sl_top_soil_col   )
    deallocate(this%qflx_runoff_col        )
    deallocate(this%qflx_runoff_r_col      )
    deallocate(this%qflx_runoff_u_col      )
    deallocate(this%qflx_rsub_sat_col      )
    deallocate(this%qflx_glcice_col        )
    deallocate(this%qflx_glcice_frz_col    )
    deallocate(this%qflx_glcice_melt_col   )  
    deallocate(this%snow_sources_col       )  
    deallocate(this%snow_sinks_col         )

    deallocate(this%qflx_irrig_col         )

    deallocate(this%qflx_snow2topsoi_col    )
    deallocate(this%qflx_h2osfc2topsoi_col   )
        
    deallocate(this%mflx_infl_col_1d)        
    deallocate(this%mflx_dew_col_1d)          
    deallocate(this%mflx_snowlyr_col_1d)
    deallocate(this%mflx_sub_snow_col_1d)     
    deallocate(this%mflx_snowlyr_col)     
    deallocate(this%mflx_neg_snow_col_1d)
    deallocate(this%mflx_et_col_1d)     
    deallocate(this%mflx_drain_col_1d)        
    deallocate(this%mflx_drain_perched_col_1d)

  end subroutine clean_col_wf

end module ColumnWaterFluxType