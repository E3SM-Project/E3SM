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
      procedure, public :: Clean => clean_col_wf

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