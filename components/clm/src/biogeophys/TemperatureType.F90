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
  use GridcellType    , only : grc_pp
  use LandunitType    , only : lun_pp                
  use ColumnType      , only : col_pp                
  use VegetationType       , only : veg_pp                
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
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars

  end type temperature_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, &
       em_roof_lun,  em_wall_lun, em_improad_lun, em_perroad_lun)

    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: em_roof_lun(bounds%begl:)
    real(r8)          , intent(in) :: em_wall_lun(bounds%begl:)
    real(r8)          , intent(in) :: em_improad_lun(bounds%begl:)
    real(r8)          , intent(in) :: em_perroad_lun(bounds%begl:)

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds,                  &
         em_roof_lun(bounds%begl:bounds%endl),    &
         em_wall_lun(bounds%begl:bounds%endl),    &
         em_improad_lun(bounds%begl:bounds%endl), &
         em_perroad_lun(bounds%begl:bounds%endl))

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

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varctl     , only : use_cn, use_cndv
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begl, endl
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    this%t_h2osfc_col(begc:endc) = spval
    call hist_addfld1d (fname='TH2OSFC',  units='K',  &
         avgflag='A', long_name='surface water temperature', &
         ptr_col=this%t_h2osfc_col)

    this%t_grnd_u_col(begc:endc) = spval
    call hist_addfld1d (fname='TG_U', units='K',  &
         avgflag='A', long_name='Urban ground temperature', &
         ptr_col=this%t_grnd_u_col, set_nourb=spval, c2l_scale_type='urbans')

    this%t_lake_col(begc:endc,:) = spval
    call hist_addfld2d (fname='TLAKE',  units='K', type2d='levlak', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=this%t_lake_col)

    this%t_soisno_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%t_soisno_col(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_T', units='K', type2d='levsno',  &
         avgflag='A', long_name='Snow temperatures', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%t_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_patch=this%t_ref2m_patch)

    this%t_ref2m_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='TSA_R', units='K',  &
         avgflag='A', long_name='Rural 2m air temperature', &
         ptr_patch=this%t_ref2m_r_patch, set_spec=spval)

    this%t_ref2m_min_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_patch)

    this%t_ref2m_max_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_patch)

    this%t_ref2m_min_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV_R', units='K',  &
         avgflag='A', long_name='Rural daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_r_patch, set_spec=spval)

    this%t_ref2m_max_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV_R', units='K',  &
         avgflag='A', long_name='Rural daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_r_patch, set_spec=spval)

    this%t_ref2m_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='TSA_U', units='K',  &
         avgflag='A', long_name='Urban 2m air temperature', &
         ptr_patch=this%t_ref2m_u_patch, set_nourb=spval)

    this%t_ref2m_min_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV_U', units='K',  &
         avgflag='A', long_name='Urban daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_u_patch, set_nourb=spval)

    this%t_ref2m_max_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV_U', units='K',  &
         avgflag='A', long_name='Urban daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_u_patch, set_nourb=spval)

    this%t_veg_patch(begp:endp) = spval
    call hist_addfld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_patch=this%t_veg_patch)

    this%t_grnd_col(begc:endc) = spval
    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=this%t_grnd_col, c2l_scale_type='urbans')

    this%t_grnd_r_col(begc:endc) = spval
    call hist_addfld1d (fname='TG_R', units='K',  &
         avgflag='A', long_name='Rural ground temperature', &
         ptr_col=this%t_grnd_r_col, set_spec=spval)

    this%t_soisno_col(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (vegetated landunits only)', &
         ptr_col=this%t_soisno_col, l2g_scale_type='veg')

    this%t_soisno_col(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI_ICE',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (ice landunits only)', &
         ptr_col=this%t_soisno_col, l2g_scale_type='ice')

    this%t_soi10cm_col(begc:endc) = spval
    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=this%t_soi10cm_col, set_urb=spval)

    if (use_cndv .or. crop_prog) then
       active = "active"
    else
       active = "inactive"
    end if
    this%t_a10_patch(begp:endp) = spval
    call hist_addfld1d (fname='T10', units='K',  &
         avgflag='A', long_name='10-day running mean of 2-m temperature', &
         ptr_patch=this%t_a10_patch, default=active)

    if (use_cn .and.  crop_prog )then
       this%t_a5min_patch(begp:endp) = spval
       call hist_addfld1d (fname='A5TMIN', units='K',  &
            avgflag='A', long_name='5-day running mean of min 2-m temperature', &
            ptr_patch=this%t_a5min_patch, default='inactive')
    end if

    if (use_cn .and. crop_prog )then
       this%t_a10min_patch(begp:endp) = spval
       call hist_addfld1d (fname='A10TMIN', units='K',  &
            avgflag='A', long_name='10-day running mean of min 2-m temperature', &
            ptr_patch=this%t_a10min_patch, default='inactive')
    end if

    this%t_building_lun(begl:endl) = spval
    call hist_addfld1d(fname='TBUILD', units='K',  &
         avgflag='A', long_name='internal urban building temperature', &
         ptr_lunit=this%t_building_lun, set_nourb=spval, l2g_scale_type='unity')

    this%hc_soi_col(begc:endc) = spval
    call hist_addfld1d (fname='HCSOI',  units='MJ/m2',  &
         avgflag='A', long_name='soil heat content', &
         ptr_col=this%hc_soi_col, set_lake=spval, set_urb=spval, l2g_scale_type='veg')

    this%hc_soisno_col(begc:endc) = spval
    call hist_addfld1d (fname='HC',  units='MJ/m2',  &
         avgflag='A', long_name='heat content of soil/snow/lake', &
         ptr_col=this%hc_soisno_col, set_urb=spval)

    this%heat1_grc(begg:endg) = spval
    call hist_addfld1d (fname='GC_HEAT1',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content', &
         ptr_lnd=this%heat1_grc)

    this%heat2_grc(begg:endg) = spval
    call hist_addfld1d (fname='GC_HEAT2',  units='J/m^2',  &
         avgflag='A', long_name='post land cover change total heat content', &
         ptr_lnd=this%heat2_grc, default='inactive')  

    this%liquid_water_temp1_grc(begg:endg) = spval
    call hist_addfld1d (fname='LIQUID_WATER_TEMP1', units='K', &
         avgflag='A', long_name='initial gridcell weighted average liquid water temperature', &
         ptr_lnd=this%liquid_water_temp1_grc, default='inactive')

    this%snot_top_col(begc:endc) = spval 
    call hist_addfld1d (fname='SNOTTOPL', units='K/m', &
         avgflag='A', long_name='snow temperature (top layer)', &
         ptr_col=this%snot_top_col, set_urb=spval, default='inactive')

    this%dTdz_top_col(begc:endc) = spval 
    call hist_addfld1d (fname='SNOdTdzL', units='K/m', &
         avgflag='A', long_name='top snow layer temperature gradient (land)', &
         ptr_col=this%dTdz_top_col, set_urb=spval, default='inactive')

    if (use_cn) then
       this%dt_veg_patch(begp:endp) = spval
       call hist_addfld1d (fname='DT_VEG', units='K', &
            avgflag='A', long_name='change in t_veg, last iteration', &
            ptr_patch=this%dt_veg_patch, default='inactive')
    end if

    if (use_cn ) then
       this%emv_patch(begp:endp) = spval
       call hist_addfld1d (fname='EMV', units='proportion', &
            avgflag='A', long_name='vegetation emissivity', &
            ptr_patch=this%emv_patch, default='inactive')
    end if

    if (use_cn) then
       this%emg_col(begc:endc) = spval
       call hist_addfld1d (fname='EMG', units='proportion', &
            avgflag='A', long_name='ground emissivity', &
            ptr_col=this%emg_col, default='inactive')
    end if

    if (use_cn) then
       this%beta_col(begc:endc) = spval
       call hist_addfld1d (fname='BETA', units='none', &
            avgflag='A', long_name='coefficient of convective velocity', &
            ptr_col=this%beta_col, default='inactive')
    end if

    ! Accumulated quantities

    this%t_veg24_patch(begp:endp) = spval
    call hist_addfld1d (fname='TV24', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 24hrs)', &
         ptr_patch=this%t_veg24_patch, default='inactive')

    this%t_veg240_patch(begp:endp)  = spval
    call hist_addfld1d (fname='TV240', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 240hrs)', &
         ptr_patch=this%t_veg240_patch, default='inactive')

    if (crop_prog) then
       this%gdd0_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDD0', units='ddays', &
            avgflag='A', long_name='Growing degree days base  0C from planting', &
            ptr_patch=this%gdd0_patch, default='inactive')
    end if

    if (crop_prog) then
       this%gdd8_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDD8', units='ddays', &
            avgflag='A', long_name='Growing degree days base  8C from planting', &
            ptr_patch=this%gdd8_patch, default='inactive')

       this%gdd10_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDD10', units='ddays', &
            avgflag='A', long_name='Growing degree days base 10C from planting', &
            ptr_patch=this%gdd10_patch, default='inactive')

       this%gdd020_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDD020', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base  0C from planting', &
            ptr_patch=this%gdd020_patch, default='inactive')

       this%gdd820_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDD820', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base  8C from planting', &
            ptr_patch=this%gdd820_patch, default='inactive')

       this%gdd1020_patch(begp:endp) = spval
       call hist_addfld1d (fname='GDD1020', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base 10C from planting', &
            ptr_patch=this%gdd1020_patch, default='inactive')

    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       em_roof_lun,  em_wall_lun, em_improad_lun, em_perroad_lun)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_const_mod  , only : SHR_CONST_TKFRZ
    use clm_varcon     , only : denice, denh2o, sb
    use landunit_varcon, only : istice, istwet, istsoil, istdlak, istice_mec
    use column_varcon  , only : icol_road_imperv, icol_roof, icol_sunwall
    use column_varcon  , only : icol_shadewall, icol_road_perv
    use clm_varctl     , only : iulog, use_vancouver, use_mexicocity
    !
    ! !ARGUMENTS:
    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: em_roof_lun(bounds%begl:bounds%endl)
    real(r8)          , intent(in) :: em_wall_lun(bounds%begl:bounds%endl)
    real(r8)          , intent(in) :: em_improad_lun(bounds%begl:bounds%endl)
    real(r8)          , intent(in) :: em_perroad_lun(bounds%begl:bounds%endl)
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p ! indices
    integer  :: nlevs   ! number of levels
    real(r8) :: snowbd  ! temporary calculation of snow bulk density (kg/m3)
    real(r8) :: fmelt   ! snowbd/100
    integer  :: lev
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(em_roof_lun)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_wall_lun)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_improad_lun) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_perroad_lun) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))

    associate(snl => col_pp%snl) ! Output: [integer (:)    ]  number of snow layers   

      ! Set snow/soil temperature
      ! t_lake only has valid values over non-lake   
      ! t_soisno, t_grnd and t_veg have valid values over all land 

      do c = bounds%begc,bounds%endc
         l = col_pp%landunit(c)

         this%t_soisno_col(c,-nlevsno+1:nlevgrnd) = spval

         ! Snow level temperatures - all land points
         if (snl(c) < 0) then
            do j = snl(c)+1, 0
               this%t_soisno_col(c,j) = 250._r8
            end do
         end if

         ! Below snow temperatures - nonlake points (lake points are set below)
         if (.not. lun_pp%lakpoi(l)) then 

            if (lun_pp%itype(l)==istice .or. lun_pp%itype(l)==istice_mec) then
               this%t_soisno_col(c,1:nlevgrnd) = 250._r8

            else if (lun_pp%itype(l) == istwet) then
               this%t_soisno_col(c,1:nlevgrnd) = 277._r8

            else if (lun_pp%urbpoi(l)) then
               if (use_vancouver) then
                  if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                     ! Set road top layer to initial air temperature and interpolate other
                     ! layers down to 20C in bottom layer
                     do j = 1, nlevgrnd
                        this%t_soisno_col(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevgrnd-1)) 
                     end do
                     ! Set wall and roof layers to initial air temperature
                  else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                       .or. col_pp%itype(c) == icol_roof) then
                     this%t_soisno_col(c,1:nlevurb) = 297.56
                  else
                     this%t_soisno_col(c,1:nlevgrnd) = 283._r8
                  end if
               else if (use_mexicocity) then
                  if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                     ! Set road top layer to initial air temperature and interpolate other
                     ! layers down to 22C in bottom layer
                     do j = 1, nlevgrnd
                        this%t_soisno_col(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevgrnd-1)) 
                     end do
                  else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                       .or. col_pp%itype(c) == icol_roof) then
                     ! Set wall and roof layers to initial air temperature
                     this%t_soisno_col(c,1:nlevurb) = 289.46
                  else
                     this%t_soisno_col(c,1:nlevgrnd) = 283._r8
                  end if
               else
                  if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                     this%t_soisno_col(c,1:nlevgrnd) = 274._r8
                  else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                       .or. col_pp%itype(c) == icol_roof) then
                     ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
                     ! shock from large heating/air conditioning flux
                     this%t_soisno_col(c,1:nlevurb) = 292._r8
                  end if
               end if
            else
               this%t_soisno_col(c,1:nlevgrnd) = 274._r8
            endif
         endif
      end do

      ! Set Ground temperatures

      do c = bounds%begc,bounds%endc
         l = col_pp%landunit(c)

         if (lun_pp%lakpoi(l)) then 
            this%t_grnd_col(c) = 277._r8
         else
            this%t_grnd_col(c) = this%t_soisno_col(c,snl(c)+1)
         end if
         this%t_soi17cm_col(c) = this%t_grnd_col(c)
      end do

      do c = bounds%begc,bounds%endc
         l = col_pp%landunit(c)
         if (lun_pp%lakpoi(l)) then ! lake
            this%t_lake_col(c,1:nlevlak) = this%t_grnd_col(c)
            this%t_soisno_col(c,1:nlevgrnd) = this%t_grnd_col(c)
         end if
      end do

      ! Set t_h2osfc_col

      this%t_h2osfc_col(bounds%begc:bounds%endc)  = 274._r8

      ! Set t_veg, t_ref2m, t_ref2m_u and tref2m_r 

      do p = bounds%begp, bounds%endp
         c = veg_pp%column(p)
         l = veg_pp%landunit(p)

         if (use_vancouver) then
            this%t_veg_patch(p)   = 297.56
         else if (use_mexicocity) then
            this%t_veg_patch(p)   = 289.46
         else
            this%t_veg_patch(p)   = 283._r8
         end if

         if (use_vancouver) then
            this%t_ref2m_patch(p) = 297.56
         else if (use_mexicocity) then
            this%t_ref2m_patch(p) = 289.46
         else
            this%t_ref2m_patch(p) = 283._r8
         end if

         if (lun_pp%urbpoi(l)) then
            if (use_vancouver) then
               this%t_ref2m_u_patch(p) = 297.56
            else if (use_mexicocity) then
               this%t_ref2m_u_patch(p) = 289.46
            else
               this%t_ref2m_u_patch(p) = 283._r8
            end if
         else 
            if (.not. lun_pp%ifspecial(l)) then 
               if (use_vancouver) then
                  this%t_ref2m_r_patch(p) = 297.56
               else if (use_mexicocity) then
                  this%t_ref2m_r_patch(p) = 289.46
               else
                  this%t_ref2m_r_patch(p) = 283._r8
               end if
            else 
               this%t_ref2m_r_patch(p) = spval
            end if
         end if

      end do

    end associate

    do l = bounds%begl, bounds%endl 
       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%taf_lun(l) = 297.56_r8
          else if (use_mexicocity) then
             this%taf_lun(l) = 289.46_r8
          else
             this%taf_lun(l) = 283._r8
          end if
       end if
    end do

    do c = bounds%begc,bounds%endc
       l = col_pp%landunit(c)

       if (col_pp%itype(c) == icol_roof       ) this%emg_col(c) = em_roof_lun(l)
       if (col_pp%itype(c) == icol_sunwall    ) this%emg_col(c) = em_wall_lun(l)
       if (col_pp%itype(c) == icol_shadewall  ) this%emg_col(c) = em_wall_lun(l)
       if (col_pp%itype(c) == icol_road_imperv) this%emg_col(c) = em_improad_lun(l)
       if (col_pp%itype(c) == icol_road_perv  ) this%emg_col(c) = em_perroad_lun(l)
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer :: j,c       ! indices
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='T_SOISNO', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='soil-snow temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_soisno_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_VEG', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='vegetation temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_veg_patch)

    call restartvar(ncid=ncid, flag=flag, varname='TH2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_h2osfc_col)
    if (flag=='read' .and. .not. readvar) then
       this%t_h2osfc_col(bounds%begc:bounds%endc) = 274.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='T_LAKE', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_lake_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_R', xtype=ncd_double,  &
         dim1name='column', &
         long_name='rural ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_r_col)
         
    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_U', xtype=ncd_double, &
         dim1name='column',                    &
         long_name='urban ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_u_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_patch)
    if (flag=='read' .and. .not. readvar) call endrun(msg=errMsg(__FILE__, __LINE__))

    call restartvar(ncid=ncid, flag=flag, varname="T_REF2M_R", xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Rural 2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_r_patch)

    call restartvar(ncid=ncid, flag=flag, varname="T_REF2M_U", xtype=ncd_double, dim1name='pft',                      &
         long_name='Urban 2m height surface air temperature', units='K',                                              &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_u_patch)


    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_r_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_U', xtype=ncd_double, dim1name='pft',                  &
         long_name='urban daily minimum of average 2 m height surface air temperature (K)', units='K',                &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_u_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_r_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_U', xtype=ncd_double, dim1name='pft',                  &
         long_name='urban daily maximum of average 2 m height surface air temperature (K)', units='K',                &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_u_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_r_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_U', xtype=ncd_double, dim1name='pft',             &
         long_name='urban instantaneous daily min of average 2 m height surface air temp (K)', units='K',             &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_u_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_r_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_U', xtype=ncd_double,  dim1name='pft',            &
         long_name='urban instantaneous daily max of average 2 m height surface air temp (K)', units='K',             &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_u_patch)

    call restartvar(ncid=ncid, flag=flag, varname='taf', xtype=ncd_double, dim1name='landunit',                       &
         long_name='urban canopy air temperature', units='K',                                                         &
         interpinic_flag='interp', readvar=readvar, data=this%taf_lun)

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='gdd1020', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 10C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd1020_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='gdd820', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 8C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd820_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='gdd020', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 0C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd020_patch)
    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    ! Each interval and accumulation type is unique to each field processed.
    ! Routine [initAccBuffer] defines the fields to be processed
    ! and the type of accumulation. 
    ! Routine [updateAccVars] does the actual accumulation for a given field.
    ! Fields are accumulated by calls to subroutine [update_accum_field]. 
    ! To accumulate a field, it must first be defined in subroutine [initAccVars] 
    ! and then accumulated by calls to [updateAccVars].
    ! Four types of accumulations are possible:
    !   o average over time interval
    !   o running mean over time interval
    !   o running accumulation over time interval
    ! Time average fields are only valid at the end of the averaging interval.
    ! Running means are valid once the length of the simulation exceeds the
    ! averaging interval. Accumulated fields are continuously accumulated.
    ! The trigger value "-99999." resets the accumulation to zero.
    !
    ! !USES 
    use accumulMod       , only : init_accum_field
    use clm_time_manager , only : get_step_size
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES: 
    real(r8) :: dtime
    integer, parameter :: not_used = huge(1)
    !---------------------------------------------------------------------

    dtime = get_step_size()

    this%t_veg24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='T_VEG24', units='K',                                              &
         desc='24hr average of vegetation temperature',  accum_type='runmean', accum_period=-1,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%t_veg240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='T_VEG240', units='K',                                             &
         desc='240hr average of vegetation temperature',  accum_type='runmean', accum_period=-10,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='TREFAV_U', units='K', &
         desc='average over an hour of urban 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='TREFAV_R', units='K', &
         desc='average over an hour of rural 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following is a running mean. The accumulation period is set to -10 for a 10-day running mean.
    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20._r8)

    if ( crop_prog )then
       call init_accum_field (name='TDM10', units='K', &
            desc='10-day running mean of min 2-m temperature', accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       call init_accum_field (name='TDM5', units='K', &
            desc='5-day running mean of min 2-m temperature', accum_type='runmean', accum_period=-5, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)
    end if

    if ( crop_prog )then

       ! All GDD summations are relative to the planting date (Kucharik & Brye 2003)
       call init_accum_field (name='GDD0', units='K', &
            desc='growing degree-days base 0C from planting', accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD8', units='K', &
            desc='growing degree-days base 8C from planting', accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD10', units='K', &
            desc='growing degree-days base 10C from planting', accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)

    end if

    if (use_cndv) then
       ! 30-day average of 2m temperature.
       call init_accum_field (name='TDA', units='K', &
            desc='30-day average of 2-m temperature', accum_type='timeavg', accum_period=-30, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

    end if

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : init_accum_field, extract_accum_field
    use clm_time_manager , only : get_nstep
    use clm_varctl       , only : nsrest, nsrStartup
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
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

    call extract_accum_field ('T_VEG24', rbufslp, nstep)
    this%t_veg24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('T_VEG240', rbufslp, nstep)
    this%t_veg240_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('T10', rbufslp, nstep)
    this%t_a10_patch(begp:endp) = rbufslp(begp:endp)

    if (crop_prog) then
       call extract_accum_field ('TDM10', rbufslp, nstep) 
       this%t_a10min_patch(begp:endp)= rbufslp(begp:endp)

       call extract_accum_field ('TDM5', rbufslp, nstep) 
       this%t_a5min_patch(begp:endp) = rbufslp(begp:endp)
    end if

    ! Initialize variables that are to be time accumulated
    ! Initialize 2m ref temperature max and min values

    if (nsrest == nsrStartup) then 
       this%t_ref2m_max_patch(begp:endp)        =  spval
       this%t_ref2m_max_r_patch(begp:endp)      =  spval
       this%t_ref2m_max_u_patch(begp:endp)      =  spval

       this%t_ref2m_min_patch(begp:endp)        =  spval
       this%t_ref2m_min_r_patch(begp:endp)      =  spval
       this%t_ref2m_min_u_patch(begp:endp)      =  spval

       this%t_ref2m_max_inst_patch(begp:endp)   = -spval
       this%t_ref2m_max_inst_r_patch(begp:endp) = -spval
       this%t_ref2m_max_inst_u_patch(begp:endp) = -spval

       this%t_ref2m_min_inst_patch(begp:endp)   =  spval
       this%t_ref2m_min_inst_r_patch(begp:endp) =  spval
       this%t_ref2m_min_inst_u_patch(begp:endp) =  spval
    end if

    if ( crop_prog ) then

       call extract_accum_field ('GDD0', rbufslp, nstep)
       this%gdd0_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD8', rbufslp, nstep) ;
       this%gdd8_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD10', rbufslp, nstep) 
       this%gdd10_patch(begp:endp) = rbufslp(begp:endp)

    end if


    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod       , only : update_accum_field, extract_accum_field, accumResetVal
    !
    ! !ARGUMENTS:
    class(temperature_type)                :: this
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

    ! Accumulate and extract T_VEG24 & T_VEG240 
    do p = begp,endp
       rbufslp(p) = this%t_veg_patch(p)
    end do
    call update_accum_field  ('T_VEG24' , rbufslp             , nstep)
    call extract_accum_field ('T_VEG24' , this%t_veg24_patch  , nstep)
    call update_accum_field  ('T_VEG240', rbufslp             , nstep)
    call extract_accum_field ('T_VEG240', this%t_veg240_patch , nstep)

    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', this%t_ref2m_patch, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)
    end_cd = is_end_curr_day()
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst_patch(p) = max(rbufslp(p), this%t_ref2m_max_inst_patch(p))
          this%t_ref2m_min_inst_patch(p) = min(rbufslp(p), this%t_ref2m_min_inst_patch(p))
       endif
       if (end_cd) then
          this%t_ref2m_max_patch(p) = this%t_ref2m_max_inst_patch(p)
          this%t_ref2m_min_patch(p) = this%t_ref2m_min_inst_patch(p)
          this%t_ref2m_max_inst_patch(p) = -spval
          this%t_ref2m_min_inst_patch(p) =  spval
       else if (secs == int(dtime)) then
          this%t_ref2m_max_patch(p) = spval
          this%t_ref2m_min_patch(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_U - hourly average urban 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_U', this%t_ref2m_u_patch, nstep)
    call extract_accum_field ('TREFAV_U', rbufslp, nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst_u_patch(p) = max(rbufslp(p), this%t_ref2m_max_inst_u_patch(p))
          this%t_ref2m_min_inst_u_patch(p) = min(rbufslp(p), this%t_ref2m_min_inst_u_patch(p))
       endif
       if (end_cd) then
         if (lun_pp%urbpoi(l)) then
          this%t_ref2m_max_u_patch(p) = this%t_ref2m_max_inst_u_patch(p)
          this%t_ref2m_min_u_patch(p) = this%t_ref2m_min_inst_u_patch(p)
          this%t_ref2m_max_inst_u_patch(p) = -spval
          this%t_ref2m_min_inst_u_patch(p) =  spval
         end if
       else if (secs == int(dtime)) then
          this%t_ref2m_max_u_patch(p) = spval
          this%t_ref2m_min_u_patch(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_R - hourly average rural 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_R', this%t_ref2m_r_patch, nstep)
    call extract_accum_field ('TREFAV_R', rbufslp, nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst_r_patch(p) = max(rbufslp(p), this%t_ref2m_max_inst_r_patch(p))
          this%t_ref2m_min_inst_r_patch(p) = min(rbufslp(p), this%t_ref2m_min_inst_r_patch(p))
       endif
       if (end_cd) then
         if (.not.(lun_pp%ifspecial(l))) then
          this%t_ref2m_max_r_patch(p) = this%t_ref2m_max_inst_r_patch(p)
          this%t_ref2m_min_r_patch(p) = this%t_ref2m_min_inst_r_patch(p)
          this%t_ref2m_max_inst_r_patch(p) = -spval
          this%t_ref2m_min_inst_r_patch(p) =  spval
         end if
       else if (secs == int(dtime)) then
          this%t_ref2m_max_r_patch(p) = spval
          this%t_ref2m_min_r_patch(p) = spval
       endif
    end do

    call update_accum_field  ('T10', this%t_ref2m_patch, nstep)
    call extract_accum_field ('T10', this%t_a10_patch, nstep)

    if ( crop_prog )then
       ! Accumulate and extract TDM10

       do p = begp,endp
          rbufslp(p) = min(this%t_ref2m_min_patch(p),this%t_ref2m_min_inst_patch(p)) !slevis: ok choice?
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                                     !'min_inst' not initialized?
       call update_accum_field  ('TDM10', rbufslp, nstep)
       call extract_accum_field ('TDM10', this%t_a10min_patch, nstep)

       ! Accumulate and extract TDM5

       do p = begp,endp
          rbufslp(p) = min(this%t_ref2m_min_patch(p),this%t_ref2m_min_inst_patch(p)) !slevis: ok choice?
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                         !'min_inst' not initialized?
       call update_accum_field  ('TDM5', rbufslp, nstep)
       call extract_accum_field ('TDM5', this%t_a5min_patch, nstep)

    end if

    if ( crop_prog )then

       ! Accumulate and extract GDD0

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(26._r8, this%t_ref2m_patch(p)-SHR_CONST_TKFRZ)) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
!       write(iulog,*) 'SPM before this one line 1258 '
       call update_accum_field  ('GDD0', rbufslp, nstep)
       call extract_accum_field ('GDD0', this%gdd0_patch, nstep)

       ! Accumulate and extract GDD8

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  this%t_ref2m_patch(p)-(SHR_CONST_TKFRZ + 8._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD8', rbufslp, nstep)
       call extract_accum_field ('GDD8', this%gdd8_patch, nstep)

       ! Accumulate and extract GDD10

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  this%t_ref2m_patch(p)-(SHR_CONST_TKFRZ + 10._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD10', rbufslp, nstep)
       call extract_accum_field ('GDD10', this%gdd10_patch, nstep)

    end if

    ! Accumulate and extract T10
    !(acumulates TSA as 10-day running mean)

    deallocate(rbufslp)

  end subroutine UpdateAccVars

end module TemperatureType
