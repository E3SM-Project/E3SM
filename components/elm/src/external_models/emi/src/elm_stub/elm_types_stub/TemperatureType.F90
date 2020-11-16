module TemperatureType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : use_cndv, iulog
  use clm_varpar      , only : nlevsno, nlevgrnd, nlevlak, nlevlak, nlevurb, crop_prog 
  use clm_varcon      , only : spval
  !
  implicit none
  save
  private
  !
  type, public :: temperature_type

     ! Temperatures
     real(r8), pointer :: t_veg_patch              (:)   ! patch vegetation temperature (Kelvin)
     real(r8), pointer :: t_h2osfc_col             (:)   ! col surface water temperature
     real(r8), pointer :: t_h2osfc_bef_col         (:)   ! col surface water temperature from time-step before  
     real(r8), pointer :: t_ssbef_col              (:,:) ! col soil/snow temperature before update (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: t_soisno_col             (:,:) ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: t_soi10cm_col            (:)   ! col soil temperature in top 10cm of soil (Kelvin)
     real(r8), pointer :: t_soi17cm_col            (:)   ! col soil temperature in top 17cm of soil (Kelvin)
     real(r8), pointer :: t_lake_col               (:,:) ! col lake temperature (Kelvin)  (1:nlevlak)          
     real(r8), pointer :: t_grnd_col               (:)   ! col ground temperature (Kelvin)
     real(r8), pointer :: t_grnd_r_col             (:)   ! col rural ground temperature (Kelvin)
     real(r8), pointer :: t_grnd_u_col             (:)   ! col urban ground temperature (Kelvin) (needed by Hydrology2Mod)
     real(r8), pointer :: t_building_lun           (:)   ! lun internal building temperature (K)
     real(r8), pointer :: snot_top_col             (:)   ! col temperature of top snow layer [K]
     real(r8), pointer :: dTdz_top_col             (:)   ! col temperature gradient in top layer  [K m-1]
     real(r8), pointer :: dt_veg_patch             (:)   ! patch change in t_veg, last iteration (Kelvin)

     real(r8), pointer :: dt_grnd_col              (:)   ! col change in t_grnd, last iteration (Kelvin)
     real(r8), pointer :: thv_col                  (:)   ! col virtual potential temperature (kelvin)
     real(r8), pointer :: thm_patch                (:)   ! patch intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
     real(r8), pointer :: t_a10_patch              (:)   ! patch 10-day running mean of the 2 m temperature (K)
     real(r8), pointer :: t_a10min_patch           (:)   ! patch 10-day running mean of min 2-m temperature
     real(r8), pointer :: t_a5min_patch            (:)   ! patch 5-day running mean of min 2-m temperature

     real(r8), pointer :: taf_lun                  (:)   ! lun urban canopy air temperature (K)

     real(r8), pointer :: t_ref2m_patch            (:)   ! patch 2 m height surface air temperature (Kelvin)
     real(r8), pointer :: t_ref2m_r_patch          (:)   ! patch rural 2 m height surface air temperature (Kelvin)
     real(r8), pointer :: t_ref2m_u_patch          (:)   ! patch urban 2 m height surface air temperature (Kelvin)
     real(r8), pointer :: t_ref2m_min_patch        (:)   ! patch daily minimum of average 2 m height surface air temperature (K)
     real(r8), pointer :: t_ref2m_min_r_patch      (:)   ! patch daily minimum of average 2 m height surface air temperature - rural(K)
     real(r8), pointer :: t_ref2m_min_u_patch      (:)   ! patch daily minimum of average 2 m height surface air temperature - urban (K)
     real(r8), pointer :: t_ref2m_max_patch        (:)   ! patch daily maximum of average 2 m height surface air temperature (K)
     real(r8), pointer :: t_ref2m_max_r_patch      (:)   ! patch daily maximum of average 2 m height surface air temperature - rural(K)
     real(r8), pointer :: t_ref2m_max_u_patch      (:)   ! patch daily maximum of average 2 m height surface air temperature - urban (K)
     real(r8), pointer :: t_ref2m_min_inst_patch   (:)   ! patch instantaneous daily min of average 2 m height surface air temp (K)
     real(r8), pointer :: t_ref2m_min_inst_r_patch (:)   ! patch instantaneous daily min of average 2 m height surface air temp - rural (K)
     real(r8), pointer :: t_ref2m_min_inst_u_patch (:)   ! patch instantaneous daily min of average 2 m height surface air temp - urban (K)
     real(r8), pointer :: t_ref2m_max_inst_patch   (:)   ! patch instantaneous daily max of average 2 m height surface air temp (K)
     real(r8), pointer :: t_ref2m_max_inst_r_patch (:)   ! patch instantaneous daily max of average 2 m height surface air temp - rural (K)
     real(r8), pointer :: t_ref2m_max_inst_u_patch (:)   ! patch instantaneous daily max of average 2 m height surface air temp - urban (K)

     ! Accumulated quantities
     !
     ! TODO(wjs, 2014-08-05) Move these to the module(s) where they are used, to improve
     ! modularity. In cases where they are used by two completely different modules,
     ! which only use the same variable out of convenience, introduce a duplicate (point
     ! being: that way one parameterization is free to change the exact meaning of its
     ! accumulator without affecting the other).
     !
     real(r8), pointer :: t_veg24_patch           (:)   ! patch 24hr average vegetation temperature (K)
     real(r8), pointer :: t_veg240_patch          (:)   ! patch 240hr average vegetation temperature (Kelvin)
     real(r8), pointer :: gdd0_patch              (:)   ! patch growing degree-days base  0C from planting  (ddays)
     real(r8), pointer :: gdd8_patch              (:)   ! patch growing degree-days base  8C from planting  (ddays)
     real(r8), pointer :: gdd10_patch             (:)   ! patch growing degree-days base 10C from planting  (ddays)
     real(r8), pointer :: gdd020_patch            (:)   ! patch 20-year average of gdd0                     (ddays)
     real(r8), pointer :: gdd820_patch            (:)   ! patch 20-year average of gdd8                     (ddays)
     real(r8), pointer :: gdd1020_patch           (:)   ! patch 20-year average of gdd10                    (ddays)

     ! Heat content
     real(r8), pointer :: beta_col                 (:)   ! coefficient of convective velocity [-]
     real(r8), pointer :: hc_soi_col               (:)   ! col soil heat content (MJ/m2)
     real(r8), pointer :: hc_soisno_col            (:)   ! col soil plus snow heat content (MJ/m2)
     real(r8), pointer :: heat1_grc                (:)   ! grc initial gridcell total heat content
     real(r8), pointer :: heat2_grc                (:)   ! grc post land cover change total heat content
     real(r8), pointer :: liquid_water_temp1_grc   (:)   ! grc initial weighted average liquid water temperature (K)
     real(r8), pointer :: liquid_water_temp2_grc   (:)   ! grc post land cover change weighted average liquid water temperature (K)

     ! Flags
     integer , pointer :: imelt_col                (:,:) ! flag for melting (=1), freezing (=2), Not=0 (-nlevsno+1:nlevgrnd) 

     ! Emissivities
     real(r8), pointer :: emv_patch                (:)   ! patch vegetation emissivity 
     real(r8), pointer :: emg_col                  (:)   ! col ground emissivity

     ! Misc
     real(r8), pointer    :: xmf_col               (:)   ! total latent heat of phase change of ground water
     real(r8), pointer    :: xmf_h2osfc_col        (:)   ! latent heat of phase change of surface water
     real(r8), pointer    :: fact_col              (:,:) ! used in computing tridiagonal matrix
     real(r8), pointer    :: c_h2osfc_col          (:)   ! heat capacity of surface water

     ! For VSFM model
     real(r8), pointer :: t_soil_col_1d            (:)   ! 1D temperature of soil layers (Kelvin)

     ! For coupling with pflotran
     real(r8), pointer :: t_nearsurf_col           (:)   ! near-surface air temperature averaged over bare-veg as BC  (Kelvin)

   contains

     procedure, public  :: Init         
     procedure, private :: InitAllocate

  end type temperature_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  

    call this%InitAllocate ( bounds )

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
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    ! Temperatures
    allocate(this%t_veg_patch              (begp:endp))                      ; this%t_veg_patch              (:)   = nan
    allocate(this%t_h2osfc_col             (begc:endc))                      ; this%t_h2osfc_col             (:)   = nan
    allocate(this%t_h2osfc_bef_col         (begc:endc))                      ; this%t_h2osfc_bef_col         (:)   = nan
    allocate(this%t_ssbef_col              (begc:endc,-nlevsno+1:nlevgrnd))  ; this%t_ssbef_col              (:,:) = nan
    allocate(this%t_soisno_col             (begc:endc,-nlevsno+1:nlevgrnd))  ; this%t_soisno_col             (:,:) = nan
    allocate(this%t_lake_col               (begc:endc,1:nlevlak))            ; this%t_lake_col               (:,:) = nan
    allocate(this%t_grnd_col               (begc:endc))                      ; this%t_grnd_col               (:)   = nan
    allocate(this%t_grnd_r_col             (begc:endc))                      ; this%t_grnd_r_col             (:)   = nan
    allocate(this%t_grnd_u_col             (begc:endc))                      ; this%t_grnd_u_col             (:)   = nan
    allocate(this%t_building_lun           (begl:endl))                      ; this%t_building_lun           (:)   = nan
    allocate(this%snot_top_col             (begc:endc))                      ; this%snot_top_col             (:)   = nan
    allocate(this%dTdz_top_col             (begc:endc))                      ; this%dTdz_top_col             (:)   = nan
    allocate(this%dt_veg_patch             (begp:endp))                      ; this%dt_veg_patch             (:)   = nan

    allocate(this%t_soi10cm_col            (begc:endc))                      ; this%t_soi10cm_col            (:)   = nan
    allocate(this%t_soi17cm_col            (begc:endc))                      ; this%t_soi17cm_col            (:)   = spval
    allocate(this%dt_grnd_col              (begc:endc))                      ; this%dt_grnd_col              (:)   = nan
    allocate(this%thv_col                  (begc:endc))                      ; this%thv_col                  (:)   = nan
    allocate(this%thm_patch                (begp:endp))                      ; this%thm_patch                (:)   = nan
    allocate(this%t_a10_patch              (begp:endp))                      ; this%t_a10_patch              (:)   = nan
    allocate(this%t_a10min_patch           (begp:endp))                      ; this%t_a10min_patch           (:)   = nan
    allocate(this%t_a5min_patch            (begp:endp))                      ; this%t_a5min_patch            (:)   = nan

    allocate(this%taf_lun                  (begl:endl))                      ; this%taf_lun                  (:)   = nan

    allocate(this%t_ref2m_patch            (begp:endp))                      ; this%t_ref2m_patch            (:)   = nan
    allocate(this%t_ref2m_r_patch          (begp:endp))                      ; this%t_ref2m_r_patch          (:)   = nan
    allocate(this%t_ref2m_u_patch          (begp:endp))                      ; this%t_ref2m_u_patch          (:)   = nan
    allocate(this%t_ref2m_min_patch        (begp:endp))                      ; this%t_ref2m_min_patch        (:)   = nan
    allocate(this%t_ref2m_min_r_patch      (begp:endp))                      ; this%t_ref2m_min_r_patch      (:)   = nan
    allocate(this%t_ref2m_min_u_patch      (begp:endp))                      ; this%t_ref2m_min_u_patch      (:)   = nan
    allocate(this%t_ref2m_max_patch        (begp:endp))                      ; this%t_ref2m_max_patch        (:)   = nan
    allocate(this%t_ref2m_max_r_patch      (begp:endp))                      ; this%t_ref2m_max_r_patch      (:)   = nan
    allocate(this%t_ref2m_max_u_patch      (begp:endp))                      ; this%t_ref2m_max_u_patch      (:)   = nan
    allocate(this%t_ref2m_max_inst_patch   (begp:endp))                      ; this%t_ref2m_max_inst_patch   (:)   = nan
    allocate(this%t_ref2m_max_inst_r_patch (begp:endp))                      ; this%t_ref2m_max_inst_r_patch (:)   = nan
    allocate(this%t_ref2m_max_inst_u_patch (begp:endp))                      ; this%t_ref2m_max_inst_u_patch (:)   = nan
    allocate(this%t_ref2m_min_inst_patch   (begp:endp))                      ; this%t_ref2m_min_inst_patch   (:)   = nan
    allocate(this%t_ref2m_min_inst_r_patch (begp:endp))                      ; this%t_ref2m_min_inst_r_patch (:)   = nan
    allocate(this%t_ref2m_min_inst_u_patch (begp:endp))                      ; this%t_ref2m_min_inst_u_patch (:)   = nan

    ! Accumulated fields
    allocate(this%t_veg24_patch            (begp:endp))                      ; this%t_veg24_patch            (:)   = nan
    allocate(this%t_veg240_patch           (begp:endp))                      ; this%t_veg240_patch           (:)   = nan
    allocate(this%gdd0_patch               (begp:endp))                      ; this%gdd0_patch               (:)   = spval
    allocate(this%gdd8_patch               (begp:endp))                      ; this%gdd8_patch               (:)   = spval
    allocate(this%gdd10_patch              (begp:endp))                      ; this%gdd10_patch              (:)   = spval
    allocate(this%gdd020_patch             (begp:endp))                      ; this%gdd020_patch             (:)   = spval
    allocate(this%gdd820_patch             (begp:endp))                      ; this%gdd820_patch             (:)   = spval
    allocate(this%gdd1020_patch            (begp:endp))                      ; this%gdd1020_patch            (:)   = spval

    ! Heat content
    allocate(this%beta_col                 (begc:endc))                      ; this%beta_col                 (:)   = nan
    allocate(this%hc_soi_col               (begc:endc))                      ; this%hc_soi_col               (:)   = nan
    allocate(this%hc_soisno_col            (begc:endc))                      ; this%hc_soisno_col            (:)   = nan
    allocate(this%heat1_grc                (begg:endg))                      ; this%heat1_grc                (:)   = nan
    allocate(this%heat2_grc                (begg:endg))                      ; this%heat2_grc                (:)   = nan
    allocate(this%liquid_water_temp1_grc   (begg:endg))                      ; this%liquid_water_temp1_grc   (:)   = nan
    allocate(this%liquid_water_temp2_grc   (begg:endg))                      ; this%liquid_water_temp2_grc   (:)   = nan

    ! flags
    allocate(this%imelt_col                (begc:endc,-nlevsno+1:nlevgrnd))  ; this%imelt_col                (:,:) = huge(1)

    ! emissivities
    allocate(this%emv_patch                (begp:endp))                      ; this%emv_patch                (:)   = nan
    allocate(this%emg_col                  (begc:endc))                      ; this%emg_col                  (:)   = nan

    allocate(this%xmf_col                  (begc:endc))                      ; this%xmf_col                  (:)   = nan
    allocate(this%xmf_h2osfc_col           (begc:endc))                      ; this%xmf_h2osfc_col           (:)   = nan
    allocate(this%fact_col                 (begc:endc, -nlevsno+1:nlevgrnd)) ; this%fact_col                 (:,:) = nan
    allocate(this%c_h2osfc_col             (begc:endc))                      ; this%c_h2osfc_col             (:)   = nan

    ! For VSFM model
    allocate(this%t_soil_col_1d            ((endc-begc+1)*nlevgrnd))         ; this%t_soil_col_1d            (:)   = nan

    ! for coupling with pflotran
    allocate(this%t_nearsurf_col           (begc:endc))                      ; this%t_nearsurf_col           (:)   = nan

  end subroutine InitAllocate

end module TemperatureType
