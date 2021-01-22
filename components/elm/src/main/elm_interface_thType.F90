module elm_interface_thType

!=================================================================================================
! ALM Thermal(T)-Hydrology (H) Interface: Data Type (Variables)
! created: 8/25/2015
! updated: 9/16/2016, 2/2/2017, June-2017
!=================================================================================================
  ! USES:
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)

  implicit none

  private

  type, public :: elm_interface_th_datatype

     ! soilstate_vars:
     real(r8), pointer :: soilpsi_col                               (:,:)   ! col soil water potential in each soil layer (MPa) (CN)

     ! waterstate_vars:
     real(r8), pointer :: frac_sno_eff_col                          (:)     ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_h2osfc_col                           (:)     ! col fractional area with surface water greater than zero
     real(r8), pointer :: h2osoi_vol_col                            (:,:)   ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2osoi_liq_col                            (:,:)   ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: h2osoi_ice_col                            (:,:)   ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)

     ! temperature_vars:
     real(r8), pointer :: t_soisno_col                              (:,:)   ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: t_grnd_col                                (:)     ! col ground(-air interface averaged) temperature (Kelvin)
     real(r8), pointer :: t_h2osfc_col                              (:)     ! col surface-water temperature [Kelvin]
     real(r8), pointer :: t_nearsurf_col                            (:)     ! col mixed air/veg. temperature near ground surface (for coupling with PFLOTRAN as BC)

     ! canopystate_vars
     integer , pointer :: alt_indx_col                              (:)     ! col current depth of thaw


     ! waterflux_vars:
     real(r8), pointer :: qflx_top_soil_col                         (:)     ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_subl_snow_col                        (:)     ! col sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_evap_soil_col                        (:)     ! col soil surface evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_snow_col                        (:)     ! col snow surface evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_h2osfc_col                      (:)     ! col evaporation flux from surface water (mm H2O /s) [+ to atm] ! note: the energy unit of orginal definition incorrect
     real(r8), pointer :: qflx_tran_veg_col                         (:)     ! col vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_rootsoil_col                         (:,:)   ! col p-aggregated vertically-resolved vegetation/soil water exchange (m H2O/s) (+ = to atm)

     real(r8), pointer :: qflx_infl_col                             (:)     ! col infiltration (mm H2O/s)
     real(r8), pointer :: qflx_surf_col                             (:)     ! col surface runoff (mm H2O/s)
     real(r8), pointer :: qflx_drain_col                            (:)     ! col sub-surface runoff (mm H2O/s)
     real(r8), pointer :: qflx_drain_vr_col                         (:,:)   ! col liquid water losted as drainage (mm H2O/s)

     ! energyflux_vars:
     real(r8), pointer :: htvp_col                                  (:)     ! latent heat of vapor of water (or sublimation) [j/kg]
     real(r8), pointer :: eflx_bot_col                              (:)     ! col heat flux from beneath the soil or ice column (W/m**2)
     real(r8), pointer :: eflx_soil_grnd_col                        (:)     ! col soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_snow_col                        (:)     ! col heat flux from snow bottom to first soil layer (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_h2osfc_col                      (:)     ! col heat flux from surface water bottom to first soil layer (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_soil_col                        (:)     ! col heat flux from near-surface air to first soil layer (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_rnet_soil_col                        (:)     ! net radiation flux between soil layer 1 and above-air, excluding SH and LE (i.e. radiation form only ) (W/m2) [+ = into soil]

     ! soilhydrology_vars:
     real(r8), pointer :: frost_table_col                           (:)     ! col frost table depth
     real(r8), pointer :: zwt_col                                   (:)     ! col water table depth
     real(r8), pointer :: zwt_perched_col                           (:)     ! col perched water table depth
     real(r8), pointer :: qcharge_col                               (:)     ! col aquifer recharge rate (mm/s)

     ! atm2lnd:
     real(r8), pointer :: forc_pbot_grc                             (:)     ! grid atm pressure (Pa)

  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type elm_interface_th_datatype
!-------------------------------------------------------------------------------------------------
  

contains

!-------------------------------------------------------------------------------------------------
  subroutine Init(this, bounds)
     use decompMod               , only : bounds_type
     class(elm_interface_th_datatype)  :: this
     type(bounds_type), intent(in)     :: bounds

     call this%InitAllocate (bounds)
  end subroutine Init
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
    ! USES
    use elm_varpar            , only : nlevsno, nlevgrnd
    use elm_varcon            , only : spval
    use decompMod             , only : bounds_type

    ! ARGUMENTS:
    real(r8) :: ival  = 0.0_r8  ! initial value
    class(elm_interface_th_datatype)  :: this
    type(bounds_type), intent(in)     :: bounds

    ! LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !------------------------------------------------------------------------
    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    !soilstate_vars:
    allocate(this%soilpsi_col           (begc:endc, 1:nlevgrnd))            ; this%soilpsi_col          (:,:) = nan

    ! waterstate_vars:
    allocate(this%frac_sno_eff_col      (begc:endc))                        ; this%frac_sno_eff_col     (:)   = nan
    allocate(this%frac_h2osfc_col       (begc:endc))                        ; this%frac_h2osfc_col      (:)   = nan
    allocate(this%h2osoi_liq_col        (begc:endc,-nlevsno+1:nlevgrnd))    ; this%h2osoi_liq_col       (:,:) = nan
    allocate(this%h2osoi_ice_col        (begc:endc,-nlevsno+1:nlevgrnd))    ; this%h2osoi_ice_col       (:,:) = nan
    allocate(this%h2osoi_vol_col        (begc:endc, 1:nlevgrnd))            ; this%h2osoi_vol_col       (:,:) = nan

    ! temperature_vars:
    allocate(this%t_soisno_col          (begc:endc,-nlevsno+1:nlevgrnd))    ; this%t_soisno_col         (:,:) = nan
    allocate(this%t_grnd_col            (begc:endc))                        ; this%t_grnd_col           (:)   = nan
    allocate(this%t_h2osfc_col          (begc:endc))                        ; this%t_h2osfc_col         (:)   = nan
    allocate(this%t_nearsurf_col        (begc:endc))                        ; this%t_nearsurf_col       (:)   = nan

    ! canopystate_vars:
    allocate(this%alt_indx_col          (begc:endc))                        ; this%alt_indx_col         (:)   = huge(1)

    !------------------------------------------------------------------------------------------
    ! pflotran variables: BEGIN
    !------------------------------------------------------------------------------------------
    ! waterflux_vars:
    allocate(this%qflx_top_soil_col     (begc:endc))                        ; this%qflx_top_soil_col     (:)   = ival
    allocate(this%qflx_evap_h2osfc_col  (begc:endc))                        ; this%qflx_evap_h2osfc_col  (:)   = ival
    allocate(this%qflx_evap_soil_col    (begc:endc))                        ; this%qflx_evap_soil_col    (:)   = ival
    allocate(this%qflx_evap_snow_col    (begc:endc))                        ; this%qflx_evap_snow_col    (:)   = ival
    allocate(this%qflx_subl_snow_col    (begc:endc))                        ; this%qflx_subl_snow_col    (:)   = ival
    allocate(this%qflx_tran_veg_col     (begc:endc))                        ; this%qflx_tran_veg_col     (:)   = ival
    allocate(this%qflx_rootsoil_col     (begc:endc,1:nlevgrnd))             ; this%qflx_rootsoil_col     (:,:) = ival
    allocate(this%qflx_infl_col         (begc:endc))                        ; this%qflx_infl_col         (:)   = ival
    allocate(this%qflx_surf_col         (begc:endc))                        ; this%qflx_surf_col         (:)   = ival
    allocate(this%qflx_drain_col        (begc:endc))                        ; this%qflx_drain_col        (:)   = ival
    allocate(this%qflx_drain_vr_col     (begc:endc,1:nlevgrnd))             ; this%qflx_drain_vr_col     (:,:) = ival


    ! energyflux_vars:
    allocate( this%htvp_col             (begc:endc))                        ; this%htvp_col             (:)   = ival
    allocate( this%eflx_bot_col         (begc:endc))                        ; this%eflx_bot_col         (:)   = ival
    allocate( this%eflx_soil_grnd_col   (begc:endc))                        ; this%eflx_soil_grnd_col   (:)   = ival
    allocate( this%eflx_fgr0_soil_col   (begc:endc))                        ; this%eflx_fgr0_soil_col   (:)   = ival
    allocate( this%eflx_fgr0_snow_col   (begc:endc))                        ; this%eflx_fgr0_snow_col   (:)   = ival
    allocate( this%eflx_fgr0_h2osfc_col (begc:endc))                        ; this%eflx_fgr0_h2osfc_col (:)   = ival
    allocate( this%eflx_rnet_soil_col   (begc:endc))                        ; this%eflx_rnet_soil_col   (:)   = ival

    ! atm2lnd:
    allocate(this%forc_pbot_grc         (begg:endg))                        ; this%forc_pbot_grc        (:)   = ival

    ! soilhydrology_vars:
    allocate( this%frost_table_col      (begc:endc))                        ; this%frost_table_col      (:)   = ival
    allocate( this%zwt_col              (begc:endc))                        ; this%zwt_col              (:)   = ival
    allocate( this%zwt_perched_col      (begc:endc))                        ; this%zwt_perched_col      (:)   = ival
    allocate( this%qcharge_col          (begc:endc))                        ; this%qcharge_col          (:)   = ival

  end subroutine InitAllocate
!-------------------------------------------------------------------------------------------------

end module elm_interface_thType
