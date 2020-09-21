module WaterstateType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module variables for hydrology
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use clm_varctl     , only : use_vancouver, use_mexicocity, use_cn, iulog, use_fates_planthydro, &
                              use_hydrstress
  use elm_varpar     , only : nlevgrnd, nlevurb, nlevsno 
  use elm_varcon     , only : spval
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterstate_type

     logical , pointer :: do_capsnow_col         (:)   ! col true => do snow capping
     real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)
     real(r8), pointer :: snow_persistence_col   (:)   ! col length of time that ground has had non-zero snow thickness (sec)
     real(r8), pointer :: snowdp_col             (:)   ! col gridcell averaged snow height (m)
     real(r8), pointer :: snowice_col            (:)   ! col average snow ice lens
     real(r8), pointer :: snowliq_col            (:)   ! col average snow liquid water
     real(r8), pointer :: int_snow_col           (:)   ! col integrated snowfall (mm H2O)
     real(r8), pointer :: snow_layer_unity_col   (:,:) ! value 1 for each snow layer, used for history diagnostics
     real(r8), pointer :: bw_col                 (:,:) ! col partial density of water in the snow pack (ice + liquid) [kg/m3] 
     real(r8), pointer :: finundated_col         (:)   ! fraction of column that is inundated, this is for bgc caclulation in betr

     real(r8), pointer :: rhvap_soi_col          (:,:)
     real(r8), pointer :: rho_vap_col            (:,:)
     real(r8), pointer :: smp_l_col              (:,:) ! col liquid phase soil matric potential, mm
     real(r8), pointer :: h2osno_col             (:)   ! col snow water (mm H2O)
     real(r8), pointer :: h2osno_old_col         (:)   ! col snow mass for previous time step (kg/m2) (new)
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_liq_old_col     (:,:)
     real(r8), pointer :: h2osoi_ice_old_col     (:,:)
     real(r8), pointer :: h2osoi_liqice_10cm_col (:)   ! col liquid water + ice lens in top 10cm of soil (kg/m2)
     real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: air_vol_col            (:,:) ! col air filled porosity     
     real(r8), pointer :: h2osoi_liqvol_col      (:,:) ! col volumetric liquid water content (v/v)
     real(r8), pointer :: h2osoi_icevol_col      (:,:) ! col volumetric ice content (v/v)     
     real(r8), pointer :: h2ocan_patch           (:)   ! patch canopy water (mm H2O)
     real(r8), pointer :: h2ocan_col             (:)   ! col canopy water (mm H2O)
     real(r8), pointer :: h2osfc_col             (:)   ! col surface water (mm H2O)
     real(r8), pointer :: swe_old_col            (:,:) ! col initial snow water
     real(r8), pointer :: liq1_grc               (:)   ! grc initial gridcell total h2o liq content
     real(r8), pointer :: liq2_grc               (:)   ! grc post land cover change total liq content
     real(r8), pointer :: ice1_grc               (:)   ! grc initial gridcell total h2o ice content
     real(r8), pointer :: ice2_grc               (:)   ! grc post land cover change total ice content
     real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)
     real(r8), pointer :: tws_month_beg_grc      (:)   ! grc total water storage at the beginning of a month
     real(r8), pointer :: tws_month_end_grc      (:)   ! grc total water storage at the end of a month

     real(r8), pointer :: total_plant_stored_h2o_col(:) ! col water that is bound in plants, including roots, sapwood, leaves, etc
                                                        ! in most cases, the vegetation scheme does not have a dynamic
                                                        ! water storage in plants, and thus 0.0 is a suitable for the trivial case.
                                                        ! When FATES is coupled in with plant hydraulics turned on, this storage
                                                        ! term is set to non-zero. (kg/m2 H2O)

     real(r8), pointer :: snw_rds_col            (:,:) ! col snow grain radius (col,lyr)    [m^-6, microns]
     real(r8), pointer :: snw_rds_top_col        (:)   ! col snow grain radius (top layer)  [m^-6, microns]
     real(r8), pointer :: h2osno_top_col         (:)   ! col top-layer mass of snow  [kg]
     real(r8), pointer :: sno_liq_top_col        (:)   ! col snow liquid water fraction (mass), top layer  [fraction]

     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)
     real(r8), pointer :: rh_ref2m_patch         (:)   ! patch 2 m height surface relative humidity (%)
     real(r8), pointer :: rh_ref2m_r_patch       (:)   ! patch 2 m height surface relative humidity - rural (%)
     real(r8), pointer :: rh_ref2m_u_patch       (:)   ! patch 2 m height surface relative humidity - urban (%)
     real(r8), pointer :: rh_af_patch            (:)   ! patch fractional humidity of canopy air (dimensionless) ! private
     real(r8), pointer :: qg_snow_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_soil_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_h2osfc_col          (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_col                 (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: dqgdT_col              (:)   ! col d(qg)/dT
     real(r8), pointer :: qaf_lun                (:)   ! lun urban canopy air specific humidity (kg/kg)

     ! Fractions
     real(r8), pointer :: frac_sno_col           (:)   ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_sno_eff_col       (:)   ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_iceold_col        (:,:) ! col fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: frac_h2osfc_col        (:)   ! col fractional area with surface water greater than zero
     real(r8), pointer :: wf_col                 (:)   ! col soil water as frac. of whc for top 0.05 m (0-1) 
     real(r8), pointer :: wf2_col                (:)   ! col soil water as frac. of whc for top 0.17 m (0-1) 
     real(r8), pointer :: fwet_patch             (:)   ! patch canopy fraction that is wet (0 to 1)
     real(r8), pointer :: fdry_patch             (:)   ! patch canopy fraction of foliage that is green and dry [-] (new)

     ! Balance Checks

     real(r8), pointer :: begwb_patch            (:)   ! water mass begining of the time step
     real(r8), pointer :: begwb_col              (:)   ! water mass begining of the time step
     real(r8), pointer :: begwb_grc              (:)   ! water mass begining of the time step
     real(r8), pointer :: endwb_patch            (:)   ! water mass end of the time step
     real(r8), pointer :: endwb_col              (:)   ! water mass end of the time step
     real(r8), pointer :: endwb_grc              (:)   ! water mass end of the time step
     real(r8), pointer :: errh2o_patch           (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2o_col             (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2o_grc             (:)   ! water conservation error (mm H2O)
     real(r8), pointer :: errh2osno_col          (:)   ! snow water conservation error(mm H2O)

     real(r8), pointer :: h2osoi_liq_depth_intg_col(:) ! grid-level depth integrated liquid soil water
     real(r8), pointer :: h2osoi_ice_depth_intg_col(:) ! grid-level depth integrated ice soil water

     real(r8), pointer :: beg_h2ocan_grc         (:)   ! grid-level canopy water at begining of the time step (mm H2O)
     real(r8), pointer :: beg_h2osno_grc         (:)   ! grid-level snow water at begining of the time step (mm H2O)
     real(r8), pointer :: beg_h2osfc_grc         (:)   ! grid-level surface water at begining of the time step (mm H2O)
     real(r8), pointer :: beg_h2osoi_liq_grc     (:)   ! grid-level liquid water at begining of the time step (kg/m2)
     real(r8), pointer :: beg_h2osoi_ice_grc     (:)   ! grid-level ice lens at begining of the time step (kg/m2)

     real(r8), pointer :: end_h2ocan_grc         (:)   ! grid-level canopy water at end of the time step (mm H2O)
     real(r8), pointer :: end_h2osno_grc         (:)   ! grid-level snow water at end of the time step (mm H2O)
     real(r8), pointer :: end_h2osfc_grc         (:)   ! grid-level surface water at end of the time step (mm H2O)
     real(r8), pointer :: end_h2osoi_liq_grc     (:)   ! grid-level liquid water at end of the time step (kg/m2)
     real(r8), pointer :: end_h2osoi_ice_grc     (:)   ! grid-level ice lens at end of the time step (kg/m2)

     ! For VSFM
     real(r8), pointer :: vsfm_fliq_col_1d       (:)   ! fraction of liquid saturation for VSFM [-]
     real(r8), pointer :: vsfm_sat_col_1d        (:)   ! liquid saturation from VSFM [-]
     real(r8), pointer :: vsfm_mass_col_1d       (:)   ! liquid mass per unit area from VSFM [kg H2O/m^2]
     real(r8), pointer :: vsfm_smpl_col_1d       (:)   ! 1D soil matrix potential liquid from VSFM [m]
     real(r8), pointer :: vsfm_soilp_col_1d      (:)   ! 1D soil liquid pressure from VSFM [Pa]
     real(r8), pointer :: soilp_col              (:,:) ! col soil pressure (Pa)

   contains

     procedure          :: Init         
     procedure          :: Restart      
     procedure, public  :: Reset 
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     
     procedure, public  :: save_h2osoi_old
  end type waterstate_type
  ! minimum allowed snow effective radius (also "fresh snow" value) [microns]
  real(r8), public, parameter :: snw_rds_min = 54.526_r8    

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, &
       h2osno_input_col, snow_depth_input_col, watsat_col, t_soisno_col)

    class(waterstate_type)            :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(inout) :: h2osno_input_col(bounds%begc:)
    real(r8)          , intent(inout) :: snow_depth_input_col(bounds%begc:)
    !real(r8)          , intent(inout) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    !real(r8)          , intent(inout) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    real(r8)              , intent(inout)    :: watsat_col(bounds%begc:bounds%endc, 1:nlevgrnd)          ! volumetric soil water at saturation (porosity)
    real(r8)              , intent(inout)    :: t_soisno_col(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd) ! col soil temperature (Kelvin)

    call this%InitAllocate(bounds) 

    call this%InitHistory(bounds)

    call this%InitCold(bounds, &
       h2osno_input_col, snow_depth_input_col, watsat_col, t_soisno_col)

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
    class(waterstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    integer :: ncells
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate(this%do_capsnow_col         (begc:endc))                   
    allocate(this%snow_depth_col         (begc:endc))                     ; this%snow_depth_col         (:)   = nan
    allocate(this%snow_persistence_col   (begc:endc))                     ; this%snow_persistence_col   (:)   = nan
    allocate(this%snowdp_col             (begc:endc))                     ; this%snowdp_col             (:)   = nan
    allocate(this%snowice_col            (begc:endc))                     ; this%snowice_col            (:)   = nan   
    allocate(this%snowliq_col            (begc:endc))                     ; this%snowliq_col            (:)   = nan   
    allocate(this%int_snow_col           (begc:endc))                     ; this%int_snow_col           (:)   = nan   
    allocate(this%snow_layer_unity_col   (begc:endc,-nlevsno+1:0))        ; this%snow_layer_unity_col   (:,:) = nan
    allocate(this%bw_col                 (begc:endc,-nlevsno+1:0))        ; this%bw_col                 (:,:) = nan   
    allocate(this%smp_l_col              (begc:endc,-nlevsno+1:nlevgrnd)) ; this%smp_l_col              (:,:) = nan
    allocate(this%finundated_col         (begc:endc))                     ; this%finundated_col         (:)   = nan

    allocate(this%h2osno_col             (begc:endc))                     ; this%h2osno_col             (:)   = nan   
    allocate(this%h2osno_old_col         (begc:endc))                     ; this%h2osno_old_col         (:)   = nan   
    allocate(this%h2osoi_liqice_10cm_col (begc:endc))                     ; this%h2osoi_liqice_10cm_col (:)   = nan
    
    allocate(this%h2osoi_vol_col         (begc:endc, 1:nlevgrnd))         ; this%h2osoi_vol_col         (:,:) = nan
    allocate(this%air_vol_col            (begc:endc, 1:nlevgrnd))         ; this%air_vol_col            (:,:) = nan
    allocate(this%h2osoi_liqvol_col      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liqvol_col      (:,:) = nan
    allocate(this%h2osoi_icevol_col      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_icevol_col      (:,:) = nan    
    allocate(this%rho_vap_col            (begc:endc,-nlevsno+1:nlevgrnd)) ; this%rho_vap_col            (:,:) = nan
    allocate(this%rhvap_soi_col          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%rhvap_soi_col          (:,:) = nan
    allocate(this%h2osoi_liq_old_col     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_old_col     (:,:) = nan
    allocate(this%h2osoi_ice_old_col     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_old_col     (:,:) = nan 
    allocate(this%h2osoi_ice_col         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_col         (:,:) = nan
    allocate(this%h2osoi_liq_col         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_col         (:,:) = nan
    allocate(this%h2ocan_patch           (begp:endp))                     ; this%h2ocan_patch           (:)   = nan  
    allocate(this%h2ocan_col             (begc:endc))                     ; this%h2ocan_col             (:)   = nan  
    allocate(this%h2osfc_col             (begc:endc))                     ; this%h2osfc_col             (:)   = nan   
    allocate(this%swe_old_col            (begc:endc,-nlevsno+1:0))        ; this%swe_old_col            (:,:) = nan   
    allocate(this%liq1_grc               (begg:endg))                     ; this%liq1_grc               (:)   = nan
    allocate(this%liq2_grc               (begg:endg))                     ; this%liq2_grc               (:)   = nan
    allocate(this%ice1_grc               (begg:endg))                     ; this%ice1_grc               (:)   = nan
    allocate(this%ice2_grc               (begg:endg))                     ; this%ice2_grc               (:)   = nan
    allocate(this%tws_grc                (begg:endg))                     ; this%tws_grc                (:)   = nan
    allocate(this%tws_month_beg_grc      (begg:endg))                     ; this%tws_month_beg_grc      (:)   = nan
    allocate(this%tws_month_end_grc      (begg:endg))                     ; this%tws_month_end_grc      (:)   = nan

    allocate(this%total_plant_stored_h2o_col(begc:endc))                  ; this%total_plant_stored_h2o_col(:) = nan

    allocate(this%snw_rds_col            (begc:endc,-nlevsno+1:0))        ; this%snw_rds_col            (:,:) = nan
    allocate(this%snw_rds_top_col        (begc:endc))                     ; this%snw_rds_top_col        (:)   = nan
    allocate(this%h2osno_top_col         (begc:endc))                     ; this%h2osno_top_col         (:)   = nan
    allocate(this%sno_liq_top_col        (begc:endc))                     ; this%sno_liq_top_col        (:)   = nan

    allocate(this%qg_snow_col            (begc:endc))                     ; this%qg_snow_col            (:)   = nan   
    allocate(this%qg_soil_col            (begc:endc))                     ; this%qg_soil_col            (:)   = nan   
    allocate(this%qg_h2osfc_col          (begc:endc))                     ; this%qg_h2osfc_col          (:)   = nan   
    allocate(this%qg_col                 (begc:endc))                     ; this%qg_col                 (:)   = nan   
    allocate(this%dqgdT_col              (begc:endc))                     ; this%dqgdT_col              (:)   = nan   
    allocate(this%qaf_lun                (begl:endl))                     ; this%qaf_lun                (:)   = nan
    allocate(this%q_ref2m_patch          (begp:endp))                     ; this%q_ref2m_patch          (:)   = nan
    allocate(this%rh_ref2m_patch         (begp:endp))                     ; this%rh_ref2m_patch         (:)   = nan
    allocate(this%rh_ref2m_u_patch       (begp:endp))                     ; this%rh_ref2m_u_patch       (:)   = nan
    allocate(this%rh_ref2m_r_patch       (begp:endp))                     ; this%rh_ref2m_r_patch       (:)   = nan
    allocate(this%rh_af_patch            (begp:endp))                     ; this%rh_af_patch            (:)   = nan

    allocate(this%frac_sno_col           (begc:endc))                     ; this%frac_sno_col           (:)   = nan
    allocate(this%frac_sno_eff_col       (begc:endc))                     ; this%frac_sno_eff_col       (:)   = nan
    allocate(this%frac_iceold_col        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%frac_iceold_col        (:,:) = nan
    allocate(this%frac_h2osfc_col        (begc:endc))                     ; this%frac_h2osfc_col        (:)   = nan 
    allocate(this%wf_col                 (begc:endc))                     ; this%wf_col                 (:)   = nan
    allocate(this%wf2_col                (begc:endc))                     ; 
    allocate(this%fwet_patch             (begp:endp))                     ; this%fwet_patch             (:)   = nan
    allocate(this%fdry_patch             (begp:endp))                     ; this%fdry_patch             (:)   = nan

    allocate(this%begwb_patch            (begp:endp))                     ; this%begwb_patch            (:)   = nan
    allocate(this%begwb_col              (begc:endc))                     ; this%begwb_col              (:)   = nan
    allocate(this%begwb_grc              (begg:endg))                     ; this%begwb_grc              (:)   = nan
    allocate(this%endwb_patch            (begp:endp))                     ; this%endwb_patch            (:)   = nan
    allocate(this%endwb_col              (begc:endc))                     ; this%endwb_col              (:)   = nan
    allocate(this%endwb_grc              (begg:endg))                     ; this%endwb_grc              (:)   = nan
    allocate(this%errh2o_patch           (begp:endp))                     ; this%errh2o_patch           (:)   = nan
    allocate(this%errh2o_col             (begc:endc))                     ; this%errh2o_col             (:)   = nan
    allocate(this%errh2o_grc             (begg:endg))                     ; this%errh2o_grc             (:)   = nan
    allocate(this%errh2osno_col          (begc:endc))                     ; this%errh2osno_col          (:)   = nan

    allocate(this%h2osoi_liq_depth_intg_col(begc:endc))                   ; this%h2osoi_liq_depth_intg_col(:) = nan
    allocate(this%h2osoi_ice_depth_intg_col(begc:endc))                   ; this%h2osoi_ice_depth_intg_col(:) = nan

    allocate(this%beg_h2ocan_grc          (begg:endg))                    ; this%beg_h2ocan_grc         (:)   = nan
    allocate(this%beg_h2osno_grc          (begg:endg))                    ; this%beg_h2osno_grc         (:)   = nan
    allocate(this%beg_h2osfc_grc          (begg:endg))                    ; this%beg_h2osfc_grc         (:)   = nan
    allocate(this%beg_h2osoi_liq_grc      (begg:endg))                    ; this%beg_h2osoi_liq_grc     (:)   = nan
    allocate(this%beg_h2osoi_ice_grc      (begg:endg))                    ; this%beg_h2osoi_ice_grc     (:)   = nan

    allocate(this%end_h2ocan_grc          (begg:endg))                    ; this%end_h2ocan_grc         (:)   = nan
    allocate(this%end_h2osno_grc          (begg:endg))                    ; this%end_h2osno_grc         (:)   = nan
    allocate(this%end_h2osfc_grc          (begg:endg))                    ; this%end_h2osfc_grc         (:)   = nan
    allocate(this%end_h2osoi_liq_grc      (begg:endg))                    ; this%end_h2osoi_liq_grc     (:)   = nan
    allocate(this%end_h2osoi_ice_grc      (begg:endg))                    ; this%end_h2osoi_ice_grc     (:)   = nan

    ncells = (endc - begc + 1)*nlevgrnd
    allocate(this%vsfm_fliq_col_1d(          ncells))                     ; this%vsfm_fliq_col_1d       (:)   = nan
    allocate(this%vsfm_sat_col_1d(           ncells))                     ; this%vsfm_sat_col_1d        (:)   = nan
    allocate(this%vsfm_mass_col_1d(          ncells))                     ; this%vsfm_mass_col_1d       (:)   = nan
    allocate(this%vsfm_smpl_col_1d(          ncells))                     ; this%vsfm_smpl_col_1d       (:)   = nan
    allocate(this%vsfm_soilp_col_1d(         ncells))                     ; this%vsfm_soilp_col_1d      (:)   = nan
    allocate(this%soilp_col              (begc:endc,1:nlevgrnd))          ; this%soilp_col              (:,:) = 0._r8

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : create_glacier_mec_landunit, use_cn, use_lch4
    use clm_varctl     , only : hist_wrtch4diag
    use elm_varpar     , only : nlevsno, crop_prog 
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal, no_snow_zero
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
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

    ! h2osno also includes snow that is part of the soil column (an 
    ! initial snow layer is only created if h2osno > 10mm). 

    ! Snow properties - these will be vertically averaged over the snow profile

    ! We determine the fractional time (and fraction of the grid cell) over which each
    ! snow layer existed by running the snow averaging routine on a field whose value is 1
    ! everywhere
    data2dptr => this%snow_layer_unity_col(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_EXISTENCE', units='unitless', type2d='levsno', &
         avgflag='A', long_name='Fraction of averaging period for which each snow layer existed', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_zero, default='inactive')


  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       h2osno_input_col, snow_depth_input_col, watsat_col, t_soisno_col)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use shr_const_mod   , only : shr_const_pi
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use shr_spfn_mod    , only : shr_spfn_erf
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use elm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use landunit_varcon , only : istice, istwet, istsoil, istdlak, istcrop, istice_mec  
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use elm_varcon      , only : denice, denh2o, spval, sb, bdsno 
    use elm_varcon      , only : h2osno_max, zlnd, tfrz, spval, pc
    use clm_varctl      , only : fsurdat, iulog
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(waterstate_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: snow_depth_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)              , intent(in)    :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs,nlevbed 
    real(r8)           :: maxslope, slopemax, minslope
    real(r8)           :: d, fd, dfdd, slope0,slopebeta
    real(r8) ,pointer  :: std (:)     
    logical            :: readvar 
    type(file_desc_t)  :: ncid        
    character(len=256) :: locfn       
    real(r8)           :: snowbd      ! temporary calculation of snow bulk density (kg/m3)
    real(r8)           :: fmelt       ! snowbd/100
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osno_input_col)     == (/bounds%endc/))          , errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(snow_depth_input_col) == (/bounds%endc/))          , errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(watsat_col)           == (/bounds%endc,nlevgrnd/)) , errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col)         == (/bounds%endc,nlevgrnd/)) , errMsg(__FILE__, __LINE__))

    ! The first three arrays are initialized from the input argument
    do c = bounds%begc,bounds%endc
       this%h2osno_col(c)             = h2osno_input_col(c) 
       this%int_snow_col(c)           = h2osno_input_col(c) 
       this%snow_depth_col(c)         = snow_depth_input_col(c)
       this%snow_persistence_col(c)   = 0._r8
       this%snow_layer_unity_col(c,:) = 1._r8
    end do

    do c = bounds%begc,bounds%endc
       this%wf_col(c) = spval
       this%wf2_col(c) = spval
    end do

    do l = bounds%begl, bounds%endl 
       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%qaf_lun(l) = 0.0111_r8
          else if (use_mexicocity) then
             this%qaf_lun(l) = 0.00248_r8
          else
             this%qaf_lun(l) = 1.e-4_r8 ! Arbitrary set since forc_q is not yet available
          end if
       end if
    end do

    ! Water Stored in plants is almost always a static entity, with the exception
    ! of when FATES-hydraulics is used. As such, this is trivially set to 0.0 (rgk 03-2017)
    this%total_plant_stored_h2o_col(bounds%begc:bounds%endc) = 0.0_r8
    
    associate(snl => col_pp%snl, nlev2bed => col_pp%nlevbed) 

      this%h2osfc_col(bounds%begc:bounds%endc) = 0._r8
      this%h2ocan_patch(bounds%begp:bounds%endp) = 0._r8
      this%h2ocan_col(bounds%begc:bounds%endc) = 0._r8

      this%frac_h2osfc_col(bounds%begc:bounds%endc) = 0._r8

      this%fwet_patch(bounds%begp:bounds%endp) = 0._r8
      this%fdry_patch(bounds%begp:bounds%endp) = 0._r8

      !--------------------------------------------
      ! Set snow water
      !--------------------------------------------

      ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
      ! This gives more realistic values of qflx_glcice sooner in the simulation
      ! for columns with net ablation, at the cost of delaying ice formation
      ! in columns with net accumulation.

      do c = bounds%begc, bounds%endc
         l = col_pp%landunit(c)
         if (lun_pp%urbpoi(l)) then
            ! From Bonan 1996 (LSM technical note)
            this%frac_sno_col(c) = min( this%snow_depth_col(c)/0.05_r8, 1._r8)
         else
            this%frac_sno_col(c) = 0._r8
            ! snow cover fraction as in Niu and Yang 2007
            if(this%snow_depth_col(c) > 0.0)  then
               snowbd   = min(400._r8, this%h2osno_col(c)/this%snow_depth_col(c)) !bulk density of snow (kg/m3)
               fmelt    = (snowbd/100.)**1.
               ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
               ! reconsidered, optimal value of 1.5 in Niu et al., 2007
               this%frac_sno_col(c) = tanh( this%snow_depth_col(c) /(2.5 * zlnd * fmelt) )
            endif
         end if
      end do

      do c = bounds%begc,bounds%endc
         if (snl(c) < 0) then
            this%snw_rds_col(c,snl(c)+1:0)        = snw_rds_min
            this%snw_rds_col(c,-nlevsno+1:snl(c)) = 0._r8
            this%snw_rds_top_col(c)               = snw_rds_min
         elseif (this%h2osno_col(c) > 0._r8) then
            this%snw_rds_col(c,0)                 = snw_rds_min
            this%snw_rds_col(c,-nlevsno+1:-1)     = 0._r8
            this%snw_rds_top_col(c)               = spval
            this%sno_liq_top_col(c)               = spval
         else
            this%snw_rds_col(c,:)                 = 0._r8
            this%snw_rds_top_col(c)               = spval
            this%sno_liq_top_col(c)               = spval
         endif
      end do

    end associate

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       watsat_col)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use elm_varcon       , only : denice, denh2o, pondmx, watmin, spval  
    use landunit_varcon  , only : istcrop, istdlak, istsoil  
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    use subgridAveMod    , only : c2g
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    real(r8)         , intent(in)    :: watsat_col (bounds%begc:, 1:)  ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j,nlevs
    logical  :: readvar
    real(r8) :: maxwatsat    ! maximum porosity    
    real(r8) :: excess       ! excess volumetric soil water
    real(r8) :: totwat       ! total soil water (mm)
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col) == (/bounds%endc,nlevgrnd/)) , errMsg(__FILE__, __LINE__))


    call restartvar(ncid=ncid, flag=flag, varname='TWS_MONTH_BEGIN', xtype=ncd_double,  &
         dim1name='gridcell', &
         long_name='surface watertotal water storage at the beginning of a month', units='mm', &
          interpinic_flag='interp', readvar=readvar, data=this%tws_month_beg_grc)

    call restartvar(ncid=ncid, flag=flag, varname='ENDWB_COL', xtype=ncd_double, &
         dim1name='column', long_name='col-level water mass end of the time step', &
         units='mm', interpinic_flag='interp', readvar=readvar, data=this%endwb_col)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Reset(this, column)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    integer , intent(in)   :: column     ! column index
    !-----------------------------------------------------------------------

    this%snw_rds_col(column,0)  = snw_rds_min

  end subroutine Reset

  !-----------------------------------------------------------------------
  subroutine save_h2osoi_old(this,bounds)
    !
    ! !DESCRIPTION:
    ! save current water status to old
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type) , intent(in)    :: bounds  


    integer :: begc, endc

    begc = bounds%begc; endc=bounds%endc

    this%h2osoi_liq_old_col(begc:endc,:) = this%h2osoi_liq_col(begc:endc,:)
    this%h2osoi_ice_old_col(begc:endc,:) = this%h2osoi_ice_col(begc:endc,:)
      
  end subroutine save_h2osoi_old
end module WaterstateType
