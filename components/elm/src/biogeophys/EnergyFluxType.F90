module EnergyFluxType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use elm_varcon     , only : spval
  use elm_varctl     , only : iulog
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
     real(r8), pointer :: btran_min_patch         (:)   ! patch daily minimum transpiration wetness factor (0 to 1)
     real(r8), pointer :: btran_min_inst_patch    (:)   ! patch instantaneous daily minimum transpiration wetness factor (0 to 1)
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

     procedure, public  :: Init          ! MAin interface for initialisation
     procedure, private :: InitAllocate  ! Allocate and initialise variables
     procedure, private :: InitHistory   ! Set up history output
     procedure, private :: InitCold      ! Initialise for cold start
     procedure, public  :: Restart       ! Set up restart fields
     procedure, public  :: InitAccBuffer ! Initialise accumulation buffer
     procedure, public  :: InitAccVars   ! Initialise time accumulation variables
     procedure, public  :: UpdateAccVars ! Update accumulation variables

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
    use elm_varpar     , only : nlevsno, nlevgrnd, nlevlak, crop_prog 
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

    allocate( this%eflx_h2osfc_to_snow_col (begc:endc))             ; this%eflx_h2osfc_to_snow_col (:)   = spval 
    allocate( this%eflx_sh_snow_patch      (begp:endp))             ; this%eflx_sh_snow_patch      (:)   = spval 
    allocate( this%eflx_sh_soil_patch      (begp:endp))             ; this%eflx_sh_soil_patch      (:)   = spval 
    allocate( this%eflx_sh_h2osfc_patch    (begp:endp))             ; this%eflx_sh_h2osfc_patch    (:)   = spval 
    allocate( this%eflx_sh_tot_patch       (begp:endp))             ; this%eflx_sh_tot_patch       (:)   = spval 
    allocate( this%eflx_sh_tot_u_patch     (begp:endp))             ; this%eflx_sh_tot_u_patch     (:)   = spval 
    allocate( this%eflx_sh_tot_r_patch     (begp:endp))             ; this%eflx_sh_tot_r_patch     (:)   = spval 
    allocate( this%eflx_sh_grnd_patch      (begp:endp))             ; this%eflx_sh_grnd_patch      (:)   = spval 
    allocate( this%eflx_sh_veg_patch       (begp:endp))             ; this%eflx_sh_veg_patch       (:)   = spval 
    allocate( this%eflx_lh_tot_u_patch     (begp:endp))             ; this%eflx_lh_tot_u_patch     (:)   = spval 
    allocate( this%eflx_lh_tot_patch       (begp:endp))             ; this%eflx_lh_tot_patch       (:)   = spval 
    allocate( this%eflx_lh_tot_r_patch     (begp:endp))             ; this%eflx_lh_tot_r_patch     (:)   = spval 
    allocate( this%eflx_lh_grnd_patch      (begp:endp))             ; this%eflx_lh_grnd_patch      (:)   = spval 
    allocate( this%eflx_lh_vege_patch      (begp:endp))             ; this%eflx_lh_vege_patch      (:)   = spval 
    allocate( this%eflx_lh_vegt_patch      (begp:endp))             ; this%eflx_lh_vegt_patch      (:)   = spval 
    allocate( this%eflx_soil_grnd_patch    (begp:endp))             ; this%eflx_soil_grnd_patch    (:)   = spval 
    allocate( this%eflx_soil_grnd_u_patch  (begp:endp))             ; this%eflx_soil_grnd_u_patch  (:)   = spval 
    allocate( this%eflx_soil_grnd_r_patch  (begp:endp))             ; this%eflx_soil_grnd_r_patch  (:)   = spval 
    allocate( this%eflx_lwrad_net_patch    (begp:endp))             ; this%eflx_lwrad_net_patch    (:)   = spval 
    allocate( this%eflx_lwrad_net_u_patch  (begp:endp))             ; this%eflx_lwrad_net_u_patch  (:)   = spval 
    allocate( this%eflx_lwrad_net_r_patch  (begp:endp))             ; this%eflx_lwrad_net_r_patch  (:)   = spval 
    allocate( this%eflx_lwrad_out_patch    (begp:endp))             ; this%eflx_lwrad_out_patch    (:)   = spval 
    allocate( this%eflx_lwrad_out_u_patch  (begp:endp))             ; this%eflx_lwrad_out_u_patch  (:)   = spval 
    allocate( this%eflx_lwrad_out_r_patch  (begp:endp))             ; this%eflx_lwrad_out_r_patch  (:)   = spval 
    allocate( this%eflx_gnet_patch         (begp:endp))             ; this%eflx_gnet_patch         (:)   = spval 
    allocate( this%eflx_grnd_lake_patch    (begp:endp))             ; this%eflx_grnd_lake_patch    (:)   = spval 
    allocate( this%eflx_dynbal_grc         (begg:endg))             ; this%eflx_dynbal_grc         (:)   = spval 
    allocate( this%eflx_bot_col            (begc:endc))             ; this%eflx_bot_col            (:)   = spval 
    allocate( this%eflx_snomelt_col        (begc:endc))             ; this%eflx_snomelt_col        (:)   = spval 
    allocate( this%eflx_snomelt_r_col      (begc:endc))             ; this%eflx_snomelt_r_col      (:)   = spval 
    allocate( this%eflx_snomelt_u_col      (begc:endc))             ; this%eflx_snomelt_u_col      (:)   = spval 
    allocate( this%eflx_fgr12_col          (begc:endc))             ; this%eflx_fgr12_col          (:)   = spval 
    allocate( this%eflx_fgr_col            (begc:endc, 1:nlevgrnd)) ; this%eflx_fgr_col            (:,:) = spval 
    allocate( this%eflx_building_heat_col  (begc:endc))             ; this%eflx_building_heat_col  (:)   = spval 
    allocate( this%eflx_urban_ac_col       (begc:endc))             ; this%eflx_urban_ac_col       (:)   = spval 
    allocate( this%eflx_urban_heat_col     (begc:endc))             ; this%eflx_urban_heat_col     (:)   = spval 
    allocate( this%eflx_wasteheat_patch    (begp:endp))             ; this%eflx_wasteheat_patch    (:)   = spval 
    allocate( this%eflx_traffic_patch      (begp:endp))             ; this%eflx_traffic_patch      (:)   = spval 
    allocate( this%eflx_heat_from_ac_patch (begp:endp))             ; this%eflx_heat_from_ac_patch (:)   = spval 
    allocate( this%eflx_heat_from_ac_lun   (begl:endl))             ; this%eflx_heat_from_ac_lun   (:)   = spval 
    allocate( this%eflx_traffic_lun        (begl:endl))             ; this%eflx_traffic_lun        (:)   = spval 
    allocate( this%eflx_wasteheat_lun      (begl:endl))             ; this%eflx_wasteheat_lun      (:)   = spval 
    allocate( this%eflx_anthro_patch       (begp:endp))             ; this%eflx_anthro_patch       (:)   = spval 

    allocate( this%eflx_hs_top_snow_col    (begc:endc))             ; this%eflx_hs_top_snow_col    (:)   = spval
    allocate( this%eflx_hs_h2osfc_col      (begc:endc))             ; this%eflx_hs_h2osfc_col      (:)   = spval
    allocate( this%eflx_hs_soil_col        (begc:endc))             ; this%eflx_hs_soil_col        (:)   = spval
    allocate( this%eflx_sabg_lyr_col       (begc:endc,-nlevsno+1:1)); this%eflx_sabg_lyr_col       (:,:) = spval
    allocate( this%eflx_dhsdT_col          (begc:endc))             ; this%eflx_dhsdT_col          (:)   = spval

    allocate( this%dgnetdT_patch           (begp:endp))             ; this%dgnetdT_patch           (:)   = spval
    allocate( this%cgrnd_patch             (begp:endp))             ; this%cgrnd_patch             (:)   = spval
    allocate( this%cgrndl_patch            (begp:endp))             ; this%cgrndl_patch            (:)   = spval
    allocate( this%cgrnds_patch            (begp:endp))             ; this%cgrnds_patch            (:)   = spval
    allocate( this%dlrad_patch             (begp:endp))             ; this%dlrad_patch             (:)   = spval
    allocate( this%ulrad_patch             (begp:endp))             ; this%ulrad_patch             (:)   = spval
    allocate( this%netrad_patch            (begp:endp))             ; this%netrad_patch            (:)   = spval  

    allocate( this%taux_patch              (begp:endp))             ; this%taux_patch              (:)   = spval
    allocate( this%tauy_patch              (begp:endp))             ; this%tauy_patch              (:)   = spval

    allocate( this%canopy_cond_patch       (begp:endp))             ; this%canopy_cond_patch       (:)   = spval

    allocate( this%htvp_col                (begc:endc))             ; this%htvp_col                (:)   = spval

    ! for coupling with pflotran
    allocate( this%eflx_soil_grnd_col      (begc:endc))             ; this%eflx_soil_grnd_col      (:)   = spval
    allocate( this%eflx_rnet_soil_col      (begc:endc))             ; this%eflx_rnet_soil_col      (:)   = spval
    allocate( this%eflx_fgr0_soil_col      (begc:endc))             ; this%eflx_fgr0_soil_col      (:)   = spval
    allocate( this%eflx_fgr0_snow_col      (begc:endc))             ; this%eflx_fgr0_snow_col      (:)   = spval
    allocate( this%eflx_fgr0_h2osfc_col    (begc:endc))             ; this%eflx_fgr0_h2osfc_col    (:)   = spval

    allocate(this%rresis_patch             (begp:endp,1:nlevgrnd))  ; this%rresis_patch            (:,:) = spval
    allocate(this%btran_patch              (begp:endp))             ; this%btran_patch             (:)   = spval
    allocate(this%btran2_patch             (begp:endp))             ; this%btran2_patch            (:)   = spval
    allocate(this%btran_min_patch          (begp:endp))             ; this%btran_min_patch         (:)   = spval
    allocate(this%btran_min_inst_patch     (begp:endp))             ; this%btran_min_inst_patch    (:)   = spval
    allocate( this%bsun_patch              (begp:endp))             ; this%bsun_patch              (:)   = spval 
    allocate( this%bsha_patch              (begp:endp))             ; this%bsha_patch              (:)   = spval 

    allocate( this%errsoi_patch            (begp:endp))             ; this%errsoi_patch            (:)   = spval
    allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = spval
    allocate( this%errseb_patch            (begp:endp))             ; this%errseb_patch            (:)   = spval
    allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = spval
    allocate( this%errsol_patch            (begp:endp))             ; this%errsol_patch            (:)   = spval
    allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = spval
    allocate( this%errlon_patch            (begp:endp))             ; this%errlon_patch            (:)   = spval
    allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = spval

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
    use elm_varpar     , only : nlevsno, nlevgrnd, crop_prog 
    use elm_varctl     , only : use_cn
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
    call hist_addfld1d (fname='BTRAN', units='1',  &
         avgflag='A', long_name='transpiration beta factor', &
         ptr_patch=this%btran_patch, set_lake=spval, set_urb=spval)

    this%btran_min_patch(begp:endp) = spval
    call hist_addfld1d (fname='BTRAN_DAILY_MIN', units='unitless',  &
         avgflag='A', long_name='daily minimum of transpiration beta factor', &
         ptr_patch=this%btran_min_patch, l2g_scale_type='veg')

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
    use elm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use elm_varcon      , only : denice, denh2o, sb
    use landunit_varcon , only : istice, istwet, istsoil, istdlak, istice_mec
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use column_varcon   , only : icol_shadewall, icol_road_perv
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

    call restartvar(ncid=ncid, flag=flag, varname='BTRAN_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of transpiration wetness factor', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran_min_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='BTRAN_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily minimum of transpiration wetness factor', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran_min_inst_patch) 

    call this%eflx_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine Restart


  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! Original source of this sub-routine:
    ! https://github.com/ESCOMP/CTSM
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
    use elm_time_manager , only : get_step_size_real
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES: 
    real(r8) :: dtime
    integer, parameter :: not_used = huge(1)
    !---------------------------------------------------------------------

    dtime = get_step_size_real()

    call init_accum_field(name='BTRANAVG', units='-', &
         desc='average over an hour of btran', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! Original source of this sub-routine:
    ! https://github.com/ESCOMP/CTSM
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use elm_time_manager , only : get_nstep
    use elm_varctl       , only : nsrest, nsrStartup
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Initialize variables that are to be time accumulated
    ! Initialize btran min values
    if (nsrest == nsrStartup) then
       this%btran_min_patch(begp:endp)        =  spval

       this%btran_min_inst_patch(begp:endp)   =  spval
    end if

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! Original source of this sub-routine:
    ! https://github.com/ESCOMP/CTSM
    !
    !
    ! USES
    use elm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod       , only : update_accum_field, extract_accum_field
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(energyflux_type)                :: this
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

    ! Accumulate and extract BTRAN_HRLY_AVG - hourly average btran
    ! Used to compute minimum of hourly averaged btran
    ! over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('BTRANAVG', this%btran_patch, nstep)
    call extract_accum_field ('BTRANAVG', rbufslp, nstep)
    end_cd = is_end_curr_day()
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          if (this%btran_min_inst_patch(p) == spval) then
             this%btran_min_inst_patch(p) = rbufslp(p)
          else
             this%btran_min_inst_patch(p) = min(rbufslp(p), this%btran_min_inst_patch(p))
          end if
       endif
       if (end_cd) then
          this%btran_min_patch(p) = this%btran_min_inst_patch(p)
          this%btran_min_inst_patch(p) =  spval
       else if (secs == dtime) then
          this%btran_min_patch(p) = spval
       endif
    end do

    deallocate(rbufslp)

  end subroutine UpdateAccVars
end module EnergyFluxType
