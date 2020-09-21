module TemperatureType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog
  use clm_varpar      , only : nlevsno, nlevgrnd, nlevlak, nlevlak, nlevurb, crop_prog 
  use elm_varcon      , only : spval
  use GridcellType    , only : grc_pp
  use LandunitType    , only : lun_pp                
  use ColumnType      , only : col_pp                
  use VegetationType       , only : veg_pp
  use VegetationDataType   , only : veg_es  
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

    allocate(this%t_soi10cm_col            (begc:endc))                      ; this%t_soi10cm_col            (:)   = nan
    allocate(this%t_soi17cm_col            (begc:endc))                      ; this%t_soi17cm_col            (:)   = spval
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
    use clm_varctl     , only : use_cn
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


    if (crop_prog) then
       active = "active"
    else
       active = "inactive"
    end if


    ! Accumulated quantities


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
    use elm_varcon     , only : denice, denh2o, sb
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


  end subroutine UpdateAccVars

end module TemperatureType
