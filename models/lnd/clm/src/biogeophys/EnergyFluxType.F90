module EnergyFluxType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : pft                
  !
  implicit none
  save
  private
  !
  type, public :: energyflux_type

     ! Fluxes
     real(r8), pointer :: eflx_sh_grnd_patch      (:)   ! patch sensible heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_veg_patch       (:)   ! patch sensible heat flux from leaves (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_snow_patch      (:)   ! patch sensible heat flux from snow (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_soil_patch      (:)   ! patch sensible heat flux from soil  (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_h2osfc_patch    (:)   ! patch sensible heat flux from surface water (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_patch       (:)   ! patch total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_u_patch     (:)   ! patch urban total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_r_patch     (:)   ! patch rural total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_patch       (:)   ! patch total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_u_patch     (:)   ! patch urban total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_r_patch     (:)   ! patch rural total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vegt_patch      (:)   ! patch transpiration heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vege_patch      (:)   ! patch evaporation heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_grnd_patch      (:)   ! patch evaporation heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_soil_grnd_patch    (:)   ! patch soil heat flux (W/m**2) [+ = into soil] 
     real(r8), pointer :: eflx_soil_grnd_u_patch  (:)   ! patch urban soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_soil_grnd_r_patch  (:)   ! patch rural soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_lwrad_net_patch    (:)   ! patch net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_net_r_patch  (:)   ! patch rural net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_net_u_patch  (:)   ! patch urban net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_out_patch    (:)   ! patch emitted infrared (longwave) radiation (W/m**2)
     real(r8), pointer :: eflx_lwrad_out_r_patch  (:)   ! patch rural emitted infrared (longwave) rad (W/m**2)
     real(r8), pointer :: eflx_lwrad_out_u_patch  (:)   ! patch urban emitted infrared (longwave) rad (W/m**2)
     real(r8), pointer :: eflx_snomelt_col        (:)   ! col snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_snomelt_r_col      (:)   ! col rural snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_snomelt_u_col      (:)   ! col urban snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_gnet_patch         (:)   ! patch net heat flux into ground  (W/m**2)
     real(r8), pointer :: eflx_grnd_lake_patch    (:)   ! patch net heat flux into lake / snow surface, excluding light transmission (W/m**2)
     real(r8), pointer :: eflx_dynbal_grc         (:)   ! grc dynamic land cover change conversion energy flux (W/m**2)
     real(r8), pointer :: eflx_bot_col            (:)   ! col heat flux from beneath the soil or ice column (W/m**2)
     real(r8), pointer :: eflx_fgr12_col          (:)   ! col ground heat flux between soil layers 1 and 2 (W/m**2)
     real(r8), pointer :: eflx_fgr_col            (:,:) ! col (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
     real(r8), pointer :: eflx_building_heat_col  (:)   ! col heat flux from urban building interior to urban walls, roof (W/m**2)
     real(r8), pointer :: eflx_urban_ac_col       (:)   ! col urban air conditioning flux (W/m**2)
     real(r8), pointer :: eflx_urban_heat_col     (:)   ! col urban heating flux (W/m**2)
     real(r8), pointer :: eflx_anthro_patch       (:)   ! patch total anthropogenic heat flux (W/m**2)
     real(r8), pointer :: eflx_traffic_patch      (:)   ! patch traffic sensible heat flux (W/m**2)
     real(r8), pointer :: eflx_wasteheat_patch    (:)   ! patch sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
     real(r8), pointer :: eflx_heat_from_ac_patch (:)   ! patch sensible heat flux put back into canyon due to removal by AC (W/m**2)
     real(r8), pointer :: eflx_traffic_lun        (:)   ! lun traffic sensible heat flux (W/m**2)
     real(r8), pointer :: eflx_wasteheat_lun      (:)   ! lun sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
     real(r8), pointer :: eflx_heat_from_ac_lun   (:)   ! lun sensible heat flux to be put back into canyon due to removal by AC (W/m**2)

     ! Derivatives of energy fluxes
     real(r8), pointer :: dgnetdT_patch           (:)   ! patch derivative of net ground heat flux wrt soil temp  (W/m**2 K)
     real(r8), pointer :: netrad_patch            (:)   ! col net radiation (W/m**2) [+ = to sfc]
     real(r8), pointer :: cgrnd_patch             (:)   ! col deriv. of soil energy flux wrt to soil temp [W/m2/k]
     real(r8), pointer :: cgrndl_patch            (:)   ! col deriv. of soil latent heat flux wrt soil temp  [W/m**2/k]
     real(r8), pointer :: cgrnds_patch            (:)   ! col deriv. of soil sensible heat flux wrt soil temp [W/m2/k]

     ! Canopy radiation
     real(r8), pointer :: dlrad_patch             (:)   ! col downward longwave radiation below the canopy [W/m2]
     real(r8), pointer :: ulrad_patch             (:)   ! col upward longwave radiation above the canopy [W/m2]

     ! Wind Stress
     real(r8), pointer :: taux_patch              (:)   ! patch wind (shear) stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy_patch              (:)   ! patch wind (shear) stress: n-s (kg/m/s**2)

     ! Conductance
     real(r8), pointer :: canopy_cond_patch       (:)   ! patch tracer conductance for canopy [m/s] 

     ! Transpiration
     real(r8), pointer :: btran_patch             (:)   ! patch transpiration wetness factor (0 to 1)

     ! Roots
     real(r8), pointer :: btran2_patch            (:)   ! patch root zone soil wetness factor (0 to 1) 
     real(r8), pointer :: rresis_patch            (:,:) ! patch root resistance by layer (0-1)  (nlevgrnd)

     ! Latent heat
     real(r8), pointer :: htvp_col                (:)   ! latent heat of vapor of water (or sublimation) [j/kg]

     ! Balance Checks
     real(r8), pointer :: errsoi_patch            (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errsoi_col              (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errseb_patch            (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errseb_col              (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errsol_patch            (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errsol_col              (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errlon_patch            (:)   ! longwave radiation conservation error (W/m**2)
     real(r8), pointer :: errlon_col              (:)   ! longwave radiation conservation error (W/m**2)

   contains

     procedure, public  :: Init         
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     
     procedure, public  :: Restart      

  end type energyflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, t_grnd_col)

    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds, t_grnd_col ) 

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, crop_prog 
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
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

    allocate( this%eflx_sh_snow_patch      (begp:endp))             ; this%eflx_sh_snow_patch      (:)   = nan
    allocate( this%eflx_sh_soil_patch      (begp:endp))             ; this%eflx_sh_soil_patch      (:)   = nan
    allocate( this%eflx_sh_h2osfc_patch    (begp:endp))             ; this%eflx_sh_h2osfc_patch    (:)   = nan
    allocate( this%eflx_sh_tot_patch       (begp:endp))             ; this%eflx_sh_tot_patch       (:)   = nan
    allocate( this%eflx_sh_tot_u_patch     (begp:endp))             ; this%eflx_sh_tot_u_patch     (:)   = nan
    allocate( this%eflx_sh_tot_r_patch     (begp:endp))             ; this%eflx_sh_tot_r_patch     (:)   = nan
    allocate( this%eflx_sh_grnd_patch      (begp:endp))             ; this%eflx_sh_grnd_patch      (:)   = nan
    allocate( this%eflx_sh_veg_patch       (begp:endp))             ; this%eflx_sh_veg_patch       (:)   = nan
    allocate( this%eflx_lh_tot_u_patch     (begp:endp))             ; this%eflx_lh_tot_u_patch     (:)   = nan
    allocate( this%eflx_lh_tot_patch       (begp:endp))             ; this%eflx_lh_tot_patch       (:)   = nan
    allocate( this%eflx_lh_tot_r_patch     (begp:endp))             ; this%eflx_lh_tot_r_patch     (:)   = nan
    allocate( this%eflx_lh_grnd_patch      (begp:endp))             ; this%eflx_lh_grnd_patch      (:)   = nan
    allocate( this%eflx_lh_vege_patch      (begp:endp))             ; this%eflx_lh_vege_patch      (:)   = nan
    allocate( this%eflx_lh_vegt_patch      (begp:endp))             ; this%eflx_lh_vegt_patch      (:)   = nan
    allocate( this%eflx_soil_grnd_patch    (begp:endp))             ; this%eflx_soil_grnd_patch    (:)   = nan
    allocate( this%eflx_soil_grnd_u_patch  (begp:endp))             ; this%eflx_soil_grnd_u_patch  (:)   = nan
    allocate( this%eflx_soil_grnd_r_patch  (begp:endp))             ; this%eflx_soil_grnd_r_patch  (:)   = nan
    allocate( this%eflx_lwrad_net_patch    (begp:endp))             ; this%eflx_lwrad_net_patch    (:)   = nan
    allocate( this%eflx_lwrad_net_u_patch  (begp:endp))             ; this%eflx_lwrad_net_u_patch  (:)   = nan
    allocate( this%eflx_lwrad_net_r_patch  (begp:endp))             ; this%eflx_lwrad_net_r_patch  (:)   = nan
    allocate( this%eflx_lwrad_out_patch    (begp:endp))             ; this%eflx_lwrad_out_patch    (:)   = nan
    allocate( this%eflx_lwrad_out_u_patch  (begp:endp))             ; this%eflx_lwrad_out_u_patch  (:)   = nan
    allocate( this%eflx_lwrad_out_r_patch  (begp:endp))             ; this%eflx_lwrad_out_r_patch  (:)   = nan
    allocate( this%eflx_gnet_patch         (begp:endp))             ; this%eflx_gnet_patch         (:)   = nan
    allocate( this%eflx_grnd_lake_patch    (begp:endp))             ; this%eflx_grnd_lake_patch    (:)   = nan
    allocate( this%eflx_dynbal_grc         (begg:endg))             ; this%eflx_dynbal_grc         (:)   = nan
    allocate( this%eflx_bot_col            (begc:endc))             ; this%eflx_bot_col            (:)   = nan
    allocate( this%eflx_snomelt_col        (begc:endc))             ; this%eflx_snomelt_col        (:)   = nan
    allocate( this%eflx_snomelt_r_col      (begc:endc))             ; this%eflx_snomelt_r_col      (:)   = nan
    allocate( this%eflx_snomelt_u_col      (begc:endc))             ; this%eflx_snomelt_u_col      (:)   = nan
    allocate( this%eflx_fgr12_col          (begc:endc))             ; this%eflx_fgr12_col          (:)   = nan
    allocate( this%eflx_fgr_col            (begc:endc, 1:nlevgrnd)) ; this%eflx_fgr_col            (:,:) = nan
    allocate( this%eflx_building_heat_col  (begc:endc))             ; this%eflx_building_heat_col  (:)   = nan
    allocate( this%eflx_urban_ac_col       (begc:endc))             ; this%eflx_urban_ac_col       (:)   = nan
    allocate( this%eflx_urban_heat_col     (begc:endc))             ; this%eflx_urban_heat_col     (:)   = nan
    allocate( this%eflx_wasteheat_patch    (begp:endp))             ; this%eflx_wasteheat_patch    (:)   = nan
    allocate( this%eflx_traffic_patch      (begp:endp))             ; this%eflx_traffic_patch      (:)   = nan
    allocate( this%eflx_heat_from_ac_patch (begp:endp))             ; this%eflx_heat_from_ac_patch (:)   = nan
    allocate( this%eflx_heat_from_ac_lun   (begl:endl))             ; this%eflx_heat_from_ac_lun   (:)   = nan
    allocate( this%eflx_traffic_lun        (begl:endl))             ; this%eflx_traffic_lun        (:)   = nan
    allocate( this%eflx_wasteheat_lun      (begl:endl))             ; this%eflx_wasteheat_lun      (:)   = nan
    allocate( this%eflx_anthro_patch       (begp:endp))             ; this%eflx_anthro_patch       (:)   = nan

    allocate( this%dgnetdT_patch           (begp:endp))             ; this%dgnetdT_patch           (:)   = nan
    allocate( this%cgrnd_patch             (begp:endp))             ; this%cgrnd_patch             (:)   = nan
    allocate( this%cgrndl_patch            (begp:endp))             ; this%cgrndl_patch            (:)   = nan
    allocate( this%cgrnds_patch            (begp:endp))             ; this%cgrnds_patch            (:)   = nan
    allocate( this%dlrad_patch             (begp:endp))             ; this%dlrad_patch             (:)   = nan
    allocate( this%ulrad_patch             (begp:endp))             ; this%ulrad_patch             (:)   = nan
    allocate( this%netrad_patch            (begp:endp))             ; this%netrad_patch            (:)   = nan  

    allocate( this%taux_patch              (begp:endp))             ; this%taux_patch              (:)   = nan
    allocate( this%tauy_patch              (begp:endp))             ; this%tauy_patch              (:)   = nan

    allocate( this%canopy_cond_patch       (begp:endp))             ; this%canopy_cond_patch       (:)   = nan

    allocate( this%htvp_col                (begc:endc))             ; this%htvp_col                (:)   = nan

    allocate(this%rresis_patch             (begp:endp,1:nlevgrnd))  ; this%rresis_patch            (:,:) = nan
    allocate(this%btran_patch              (begp:endp))             ; this%btran_patch             (:)   = nan
    allocate(this%btran2_patch             (begp:endp))             ; this%btran2_patch            (:)   = nan

    allocate( this%errsoi_patch            (begp:endp))             ; this%errsoi_patch            (:)   = nan
    allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = nan
    allocate( this%errseb_patch            (begp:endp))             ; this%errseb_patch            (:)   = nan
    allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = nan
    allocate( this%errsol_patch            (begp:endp))             ; this%errsol_patch            (:)   = nan
    allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = nan
    allocate( this%errlon_patch            (begp:endp))             ; this%errlon_patch            (:)   = nan
    allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = nan

  end subroutine InitAllocate
    
  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd, crop_prog 
    use clm_varctl     , only : use_cn
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
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


    this%eflx_dynbal_grc(begg:endg) = spval 
    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=this%eflx_dynbal_grc)

    this%eflx_snomelt_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=this%eflx_snomelt_col, c2l_scale_type='urbanf')

    this%eflx_snomelt_r_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM_R',  units='W/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=this%eflx_snomelt_r_col, set_spec=spval)

    this%eflx_snomelt_u_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM_U',  units='W/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=this%eflx_snomelt_u_col, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_lwrad_net_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRA', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_patch, c2l_scale_type='urbanf')

    this%eflx_lwrad_net_r_patch(begp:endp) = spval 
    call hist_addfld1d (fname='FIRA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_r_patch, set_spec=spval)

    this%eflx_lwrad_out_patch(begp:endp) = spval 
    call hist_addfld1d (fname='FIRE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_patch, c2l_scale_type='urbanf')
    ! Rename of FIRE for Urban intercomparision project
    this%eflx_lwrad_out_patch(begp:endp) = spval 
    call hist_addfld1d (fname='LWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling longwave radiation', &
         ptr_patch=this%eflx_lwrad_out_patch, c2l_scale_type='urbanf', default='inactive')

    this%eflx_lwrad_out_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_R', units='W/m^2',  &
         avgflag='A', long_name='Rural emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_r_patch, set_spec=spval)

    this%eflx_lh_vegt_patch(begp:endp) = spval
    call hist_addfld1d (fname='FCTR', units='W/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_patch=this%eflx_lh_vegt_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_vege_patch(begp:endp) = spval
    call hist_addfld1d (fname='FCEV', units='W/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_patch=this%eflx_lh_vege_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGEV', units='W/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_patch=this%eflx_lh_grnd_patch, c2l_scale_type='urbanf') 

    this%eflx_sh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_NODYNLNDUSE', units='W/m^2',  &
         avgflag='A', long_name='sensible heat not including correction for land use change', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf')

    this%eflx_sh_tot_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_R', units='W/m^2',  &
         avgflag='A', long_name='Rural sensible heat', &
         ptr_patch=this%eflx_sh_tot_r_patch, set_spec=spval)

    this%eflx_sh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qh', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qle', units='W/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
         avgflag='A', long_name='total latent heat flux [+ to atm]', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf')

    this%eflx_lh_tot_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_R', units='W/m^2',  &
         avgflag='A', long_name='Rural total evaporation', &
         ptr_patch=this%eflx_lh_tot_r_patch, set_spec=spval)

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qstor', units='W/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf', &
         default = 'inactive')
    this%eflx_sh_veg_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_V', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_patch=this%eflx_sh_veg_patch, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_sh_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_G', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_patch=this%eflx_sh_grnd_patch, c2l_scale_type='urbanf')

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGR', units='W/m^2',  &
         avgflag='A', long_name='heat flux into soil/snow including snow melt and lake / snow light transmission', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf')

    this%eflx_soil_grnd_r_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGR_R', units='W/m^2',  &
         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt and snow light transmission', &
         ptr_patch=this%eflx_soil_grnd_r_patch, set_spec=spval)

    this%eflx_lwrad_net_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_u_patch, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
         avgflag='A', long_name='soil heat flux [+ into soil]', &
         ptr_patch=this%eflx_soil_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%eflx_lwrad_out_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_U', units='W/m^2',  &
         avgflag='A', long_name='Urban emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_u_patch, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_sh_tot_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_U', units='W/m^2',  &
         avgflag='A', long_name='Urban sensible heat', &
         ptr_patch=this%eflx_sh_tot_u_patch, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_lh_tot_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_U', units='W/m^2',  &
         avgflag='A', long_name='Urban total evaporation', &
         ptr_patch=this%eflx_lh_tot_u_patch, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_soil_grnd_u_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGR_U', units='W/m^2',  &
         avgflag='A', long_name='Urban heat flux into soil/snow including snow melt', &
         ptr_patch=this%eflx_soil_grnd_u_patch, c2l_scale_type='urbanf', set_nourb=spval)

    this%netrad_patch(begp:endp) = spval
    call hist_addfld1d (fname='Rnet', units='W/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_patch=this%netrad_patch, c2l_scale_type='urbanf', &
         default='inactive')

    if (use_cn) then
       this%dlrad_patch(begp:endp) = spval
       call hist_addfld1d (fname='DLRAD', units='W/m^2', &
            avgflag='A', long_name='downward longwave radiation below the canopy', &
            ptr_patch=this%dlrad_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%ulrad_patch(begp:endp) = spval
       call hist_addfld1d (fname='ULRAD', units='W/m^2', &
            avgflag='A', long_name='upward longwave radiation above the canopy', &
            ptr_patch=this%ulrad_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%cgrnd_patch(begp:endp) = spval
       call hist_addfld1d (fname='CGRND', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil energy flux wrt to soil temp', &
            ptr_patch=this%cgrnd_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%cgrndl_patch(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDL', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil latent heat flux wrt soil temp', &
            ptr_patch=this%cgrndl_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%cgrnds_patch(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDS', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil sensible heat flux wrt soil temp', &
            ptr_patch=this%cgrnds_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    if (use_cn) then
       this%eflx_gnet_patch(begp:endp) = spval
       call hist_addfld1d (fname='EFLX_GNET', units='W/m^2', &
            avgflag='A', long_name='net heat flux into ground', &
            ptr_patch=this%eflx_gnet_patch, default='inactive', c2l_scale_type='urbanf')
    end if 

    this%eflx_grnd_lake_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_GRND_LAKE', units='W/m^2', &
         avgflag='A', long_name='net heat flux into lake/snow surface, excluding light transmission', &
         ptr_patch=this%eflx_grnd_lake_patch, set_nolake=spval)

    this%eflx_building_heat_col(begc:endc) = spval
    call hist_addfld1d (fname='BUILDHEAT', units='W/m^2',  &
         avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
         ptr_col=this%eflx_building_heat_col, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_ac_col(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
         avgflag='A', long_name='urban air conditioning flux', &
         ptr_col=this%eflx_urban_ac_col, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_heat_col(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
         avgflag='A', long_name='urban heating flux', &
         ptr_col=this%eflx_urban_heat_col, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%dgnetdT_patch(begp:endp) = spval
    call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
         avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
         ptr_patch=this%dgnetdT_patch, default='inactive', c2l_scale_type='urbanf')

    this%eflx_fgr12_col(begc:endc) = spval
    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=this%eflx_fgr12_col, set_lake=spval)

    this%eflx_fgr_col(begc:endc,:) = spval
    call hist_addfld2d (fname='FGR_SOIL_R', units='watt/m^2', type2d='levgrnd', &
         avgflag='A', long_name='Rural downward heat flux at interface below each soil layer', &
         ptr_col=this%eflx_fgr_col, set_spec=spval, default='inactive')

    this%eflx_traffic_patch(begp:endp) = spval
    call hist_addfld1d (fname='TRAFFICFLUX', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from urban traffic', &
         ptr_patch=this%eflx_traffic_patch, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    this%eflx_wasteheat_patch(begp:endp) = spval
    call hist_addfld1d (fname='WASTEHEAT', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from heating/cooling sources of urban waste heat', &
         ptr_patch=this%eflx_wasteheat_patch, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_heat_from_ac_patch(begp:endp) = spval
    call hist_addfld1d (fname='HEAT_FROM_AC', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux put into canyon due to heat removed from air conditioning', &
         ptr_patch=this%eflx_heat_from_ac_patch, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_anthro_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qanth', units='W/m^2',  &
         avgflag='A', long_name='anthropogenic heat flux', &
         ptr_patch=this%eflx_anthro_patch, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    this%taux_patch(begp:endp) = spval
    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_patch=this%taux_patch)
    ! Rename of TAUX for Urban intercomparision project (when U=V)    
    call hist_addfld1d (fname='Qtau', units='kg/m/s^2',  &  
         avgflag='A', long_name='momentum flux', &
         ptr_patch=this%taux_patch, default='inactive')

    this%tauy_patch(begp:endp) = spval
    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_patch=this%tauy_patch)

    this%btran_patch(begp:endp) = spval
    call hist_addfld1d (fname='BTRAN', units='unitless',  &
         avgflag='A', long_name='transpiration beta factor', &
         ptr_patch=this%btran_patch, set_lake=spval, set_urb=spval)

    if (use_cn) then
       this%rresis_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='RRESIS', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='root resistance in each soil layer', &
            ptr_patch=this%rresis_patch, default='inactive')
    end if

    this%errsoi_col(begc:endc) = spval
    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=this%errsoi_col)

    this%errseb_patch(begp:endp) = spval
    call hist_addfld1d (fname='ERRSEB',  units='W/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_patch=this%errseb_patch)

    this%errsol_patch(begp:endp) = spval
    call hist_addfld1d (fname='ERRSOL',  units='W/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_patch=this%errsol_patch, set_urb=spval)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, t_grnd_col)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varcon      , only : denice, denh2o, sb
    use landunit_varcon , only : istice, istwet, istsoil, istdlak, istice_mec
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use clm_varctl      , only : iulog, use_vancouver, use_mexicocity
    !
    ! !ARGUMENTS:
    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p,levs,lev
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    associate(snl => col%snl) ! Output: [integer (:)    ]  number of snow layers   

      do c = bounds%begc, bounds%endc
         l = col%landunit(c)

         if (lun%urbpoi(l)) then
            this%eflx_building_heat_col(c) = 0._r8
            this%eflx_urban_ac_col(c)      = 0._r8
            this%eflx_urban_heat_col(c)    = 0._r8
         else
            this%eflx_building_heat_col(c) = 0._r8
            this%eflx_urban_ac_col(c)      = 0._r8
            this%eflx_urban_heat_col(c)    = 0._r8
         end if

    end do

      do p = bounds%begp, bounds%endp 
         c = pft%column(p)
         l = pft%landunit(p)

         if (.not. lun%urbpoi(l)) then ! non-urban
            this%eflx_lwrad_net_u_patch(p) = spval
            this%eflx_lwrad_out_u_patch(p) = spval
            this%eflx_lh_tot_u_patch(p)    = spval
            this%eflx_sh_tot_u_patch(p)    = spval
            this%eflx_soil_grnd_u_patch(p) = spval
         end if

         this%eflx_lwrad_out_patch(p) = sb * (t_grnd_col(c))**4
      end do

    end associate

    do p = bounds%begp, bounds%endp 
       l = pft%landunit(p)

       if (.not. lun%urbpoi(l)) then
          this%eflx_traffic_lun(l)        = spval
          this%eflx_wasteheat_lun(l)      = spval

          this%eflx_wasteheat_patch(p)    = 0._r8
          this%eflx_heat_from_ac_patch(p) = 0._r8
          this%eflx_traffic_patch(p)      = 0._r8
          this%eflx_anthro_patch(p)       = 0._r8
       end if
    end do

    ! initialize rresis, for use in ecosystemdyn
    do p = bounds%begp,bounds%endp
       do lev = 1,nlevgrnd
          this%rresis_patch(p,lev) = 0._r8
       end do
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
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  & 
         dim1name='pft', &
         long_name='emitted infrared (longwave) radiation', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_lwrad_out_patch)

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC', xtype=ncd_double,  dim1name='column', &
         long_name='urban air conditioning flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_ac_col)

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT', xtype=ncd_double, dim1name='column', &
         long_name='urban heating flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_heat_col)

    call restartvar(ncid=ncid, flag=flag, varname='btran2', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran2_patch) 

  end subroutine Restart

end module EnergyFluxType
