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
  use clm_varctl     , only : use_vancouver, use_mexicocity, use_cn, iulog, use_fates_planthydro
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
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
     procedure, private :: InitAllocate
  end type waterstate_type

  ! minimum allowed snow effective radius (also "fresh snow" value) [microns]
  real(r8), public, parameter :: snw_rds_min = 54.526_r8    
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(waterstate_type)            :: this
    type(bounds_type) , intent(in)    :: bounds  

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

end module WaterstateType
