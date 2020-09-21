module EnergyFluxType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use elm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  use VegetationType      , only : veg_pp                
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  !
  implicit none
  save
  private
  !
  type, public :: energyflux_type

     ! Fluxes
     real(r8), pointer :: eflx_h2osfc_to_snow_col (:)   ! col snow melt to h2osfc heat flux (W/m**2)
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
     real(r8), pointer :: eflx_hs_h2osfc_col      (:)   ! heat flux on standing water [W/m2]
     real(r8), pointer :: eflx_hs_top_snow_col    (:)   ! heat flux on top snow layer [W/m2]
     real(r8), pointer :: eflx_hs_soil_col        (:)   ! heat flux on soil [W/m2
     real(r8), pointer :: eflx_sabg_lyr_col       (:,:) ! absorbed solar radiation (col,lyr) [W/m2]

     ! Derivatives of energy fluxes
     real(r8), pointer :: dgnetdT_patch           (:)   ! patch derivative of net ground heat flux wrt soil temp  (W/m**2 K)
     real(r8), pointer :: netrad_patch            (:)   ! col net radiation (W/m**2) [+ = to sfc]
     real(r8), pointer :: cgrnd_patch             (:)   ! col deriv. of soil energy flux wrt to soil temp [W/m2/k]
     real(r8), pointer :: cgrndl_patch            (:)   ! col deriv. of soil latent heat flux wrt soil temp  [W/m**2/k]
     real(r8), pointer :: cgrnds_patch            (:)   ! col deriv. of soil sensible heat flux wrt soil temp [W/m2/k]
     real(r8), pointer :: eflx_dhsdT_col          (:)   ! col deriv. of energy flux into surface layer wrt temp [W/m2/K]

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
     !plant hydraulics
     real(r8), pointer :: bsun_patch              (:)   ! patch sunlit canopy transpiration wetness factor (0 to 1)
     real(r8), pointer :: bsha_patch              (:)   ! patch shaded canopy transpiration wetness factor (0 to 1)

     ! Roots
     real(r8), pointer :: btran2_patch            (:)   ! patch root zone soil wetness factor (0 to 1) 
     real(r8), pointer :: rresis_patch            (:,:) ! patch root resistance by layer (0-1)  (nlevgrnd)

     ! Latent heat
     real(r8), pointer :: htvp_col                (:)   ! latent heat of vapor of water (or sublimation) [j/kg]

     ! for couplig with pflotran
     real(r8), pointer :: eflx_soil_grnd_col      (:)   ! col integrated soil ground heat flux (W/m2)  [+ = into ground]
     real(r8), pointer :: eflx_rnet_soil_col      (:)   ! col soil net (sw+lw) radiation flux (W/m2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_soil_col      (:)   ! col soil-air heat flux (W/m2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_snow_col      (:)   ! col soil-snow heat flux (W/m2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_h2osfc_col    (:)   ! col soil-surfacewater heat flux (W/m2) [+ = into soil]

     ! Balance Checks
     real(r8), pointer :: errsoi_patch            (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errsoi_col              (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errseb_patch            (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errseb_col              (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errsol_patch            (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errsol_col              (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errlon_patch            (:)   ! longwave radiation conservation error (W/m**2)
     real(r8), pointer :: errlon_col              (:)   ! longwave radiation conservation error (W/m**2)

     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: eflx_dynbal_dribbler

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

    allocate( this%eflx_h2osfc_to_snow_col (begc:endc))             ; this%eflx_h2osfc_to_snow_col (:)   = nan
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

    allocate( this%eflx_hs_top_snow_col    (begc:endc))             ; this%eflx_hs_top_snow_col    (:)   = nan
    allocate( this%eflx_hs_h2osfc_col      (begc:endc))             ; this%eflx_hs_h2osfc_col      (:)   = nan
    allocate( this%eflx_hs_soil_col        (begc:endc))             ; this%eflx_hs_soil_col        (:)   = nan
    allocate( this%eflx_sabg_lyr_col       (begc:endc,-nlevsno+1:1)); this%eflx_sabg_lyr_col       (:,:) = nan
    allocate( this%eflx_dhsdT_col          (begc:endc))             ; this%eflx_dhsdT_col          (:)   = nan

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

    ! for coupling with pflotran
    allocate( this%eflx_soil_grnd_col      (begc:endc))             ; this%eflx_soil_grnd_col      (:)   = nan
    allocate( this%eflx_rnet_soil_col      (begc:endc))             ; this%eflx_rnet_soil_col      (:)   = nan
    allocate( this%eflx_fgr0_soil_col      (begc:endc))             ; this%eflx_fgr0_soil_col      (:)   = nan
    allocate( this%eflx_fgr0_snow_col      (begc:endc))             ; this%eflx_fgr0_snow_col      (:)   = nan
    allocate( this%eflx_fgr0_h2osfc_col    (begc:endc))             ; this%eflx_fgr0_h2osfc_col    (:)   = nan

    allocate(this%rresis_patch             (begp:endp,1:nlevgrnd))  ; this%rresis_patch            (:,:) = nan
    allocate(this%btran_patch              (begp:endp))             ; this%btran_patch             (:)   = nan
    allocate(this%btran2_patch             (begp:endp))             ; this%btran2_patch            (:)   = nan
    allocate( this%bsun_patch              (begp:endp))             ; this%bsun_patch              (:)   = nan
    allocate( this%bsha_patch              (begp:endp))             ; this%bsha_patch              (:)   = nan

    allocate( this%errsoi_patch            (begp:endp))             ; this%errsoi_patch            (:)   = nan
    allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = nan
    allocate( this%errseb_patch            (begp:endp))             ; this%errseb_patch            (:)   = nan
    allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = nan
    allocate( this%errsol_patch            (begp:endp))             ; this%errsol_patch            (:)   = nan
    allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = nan
    allocate( this%errlon_patch            (begp:endp))             ; this%errlon_patch            (:)   = nan
    allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = nan

    this%eflx_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'eflx_dynbal', &
         units = 'J/m**2')

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
    use elm_varcon      , only : denice, denh2o, sb
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

    associate(snl => col_pp%snl) ! Output: [integer (:)    ]  number of snow layers   


    end associate


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

    call restartvar(ncid=ncid, flag=flag, varname='btran2', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran2_patch) 

    call this%eflx_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine Restart

end module EnergyFluxType
