module clm_interface_thType

!=================================================================================================
! ALM Thermal(T)-Hydrology (H) Interface: Data Type (Variables)
! created: 8/25/2015
! updated: 9/16/2016, 2/2/2017, June-2017
! update: 5/13/2019 (all are IMPLICITLY column-based coupling,
!        esp. after ELM v2 data-structure modification @3/12/2019, commit 1bf22e32d)
!=================================================================================================
  ! USES:
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)

  implicit none

  private

  type, public :: clm_interface_th_datatype
     ! coupling level (grid/column/patch):
     integer :: cpl_level = 2                                ! coupling level: 1 - grid, 2 - column (default), 3 - patch

     ! soilstate_vars:
     real(r8), pointer :: soilpsi                               (:,:)   ! col soil water potential in each soil layer (MPa) (CN)

     ! waterstate_vars:
     real(r8), pointer :: frac_sno_eff                          (:)     ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_h2osfc                           (:)     ! col fractional area with surface water greater than zero
     real(r8), pointer :: h2osoi_vol                            (:,:)   ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2osoi_liq                            (:,:)   ! col liquid water (kg/m2) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: h2osoi_ice                            (:,:)   ! col ice lens (kg/m2) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: h2osfc                                (:)     ! col surface water (mm H2O)

     ! temperature_vars:
     real(r8), pointer :: t_soisno                              (:,:)   ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: t_grnd                                (:)     ! col ground(-air interface averaged) temperature (Kelvin)
     real(r8), pointer :: t_h2osfc                              (:)     ! col surface-water temperature [Kelvin]
     real(r8), pointer :: t_nearsurf                            (:)     ! col mixed air/veg. temperature near ground surface (for coupling with PFLOTRAN as BC)

     ! canopystate_vars
     integer , pointer :: alt_indx                              (:)     ! col current depth of thaw


     ! waterflux_vars:
     real(r8), pointer :: qflx_top_soil                         (:)     ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_subl_snow                        (:)     ! col sublimation rate from snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_evap_soil                        (:)     ! col soil surface evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_snow                        (:)     ! col snow surface evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_h2osfc                      (:)     ! col evaporation flux from surface water (mm H2O /s) [+ to atm] ! note: the energy unit of orginal definition incorrect
     real(r8), pointer :: qflx_tran_veg                         (:)     ! col vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_rootsoil                         (:,:)   ! col p-aggregated vertically-resolved vegetation/soil water exchange (m H2O/s) (+ = to atm)
     ! new: when coupling with pflotran, evap+tran from CLM is to be as an input, which may be re-adjusted by pflotran flow module.
     real(r8), pointer :: qflx_et_reduced                       (:)     ! col net water exchange re-adjusted for mass-balance check (m H2O/s) (+ = to atm)

     real(r8), pointer :: qflx_infl                             (:)     ! col infiltration (mm H2O/s)
     real(r8), pointer :: qflx_exfl                             (:)     ! col saturate soil excess outflow upwardly (mm H2O/s)
     real(r8), pointer :: qflx_surf                             (:)     ! col ground surface runoff/flood/overland flow (mm H2O/s)
     real(r8), pointer :: qflx_h2osfc                           (:)     ! col water-covered (pond) surface runoff/flood/overflow (mm H2O/s)
     real(r8), pointer :: qflx_base                             (:)     ! col drainage/upflow at bottom (baseflow) (mm H2O/s)
     real(r8), pointer :: qflx_drain                            (:)     ! col liquid water losted as drainage (mm H2O/s)
     real(r8), pointer :: qflx_drain_vr                         (:,:)   ! col vertically-resolved liquid water losted as drainage (mm H2O/s)
     real(r8), pointer :: qflx_lateral                          (:)     ! col net liquid water lateral flow (mm H2O/s)
     real(r8), pointer :: qflx_lateral_vr                       (:,:)   ! col vertically-resolved net liquid water lateral flow (mm H2O/s)

     ! energyflux_vars:
     real(r8), pointer :: htvp                                  (:)     ! latent heat of vapor of water (or sublimation) [j/kg]
     real(r8), pointer :: eflx_bot                              (:)     ! col heat flux from beneath the soil or ice column (W/m**2)
     real(r8), pointer :: eflx_soil_grnd                        (:)     ! col soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_snow                        (:)     ! col heat flux from snow bottom to first soil layer (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_h2osfc                      (:)     ! col heat flux from surface water bottom to first soil layer (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_fgr0_soil                        (:)     ! col heat flux from near-surface air to first soil layer (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_rnet_soil                        (:)     ! net radiation flux between soil layer 1 and above-air, excluding SH and LE (i.e. radiation form only ) (W/m2) [+ = into soil]

     ! soilhydrology_vars:
     real(r8), pointer :: frost_table                           (:)     ! col frost table depth
     real(r8), pointer :: zwt                                   (:)     ! col water table depth
     real(r8), pointer :: zwt_perched                           (:)     ! col perched water table depth
     real(r8), pointer :: qcharge                               (:)     ! col aquifer recharge rate (mm/s)

     ! atm2lnd:
     real(r8), pointer :: forc_pbot                             (:)     ! col atm pressure (Pa)

  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type clm_interface_th_datatype
!-------------------------------------------------------------------------------------------------
  

contains

!-------------------------------------------------------------------------------------------------
  subroutine Init(this, bounds)
     use decompMod               , only : bounds_type
     class(clm_interface_th_datatype)  :: this
     type(bounds_type), intent(in)     :: bounds

     call this%InitAllocate (bounds)
  end subroutine Init
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
    ! USES
    use clm_varpar            , only : nlevsno, nlevgrnd
    use clm_varcon            , only : spval
    use decompMod             , only : bounds_type

    ! ARGUMENTS:
    real(r8) :: ival  = 0.0_r8  ! initial value
    class(clm_interface_th_datatype)  :: this
    type(bounds_type), intent(in)     :: bounds

    ! LOCAL VARIABLES:
    integer  :: begc, endc
    !------------------------------------------------------------------------
    begc = bounds%begc; endc= bounds%endc
    if (this%cpl_level==1) then
      begc = bounds%begg; endc= bounds%endg
    elseif (this%cpl_level==3) then
      begc = bounds%begp; endc= bounds%endp
    endif

    !soilstate_vars:
    allocate(this%soilpsi           (begc:endc, 1:nlevgrnd))            ; this%soilpsi          (:,:) = nan

    ! waterstate_vars:
    allocate(this%frac_sno_eff      (begc:endc))                        ; this%frac_sno_eff     (:)   = nan
    allocate(this%frac_h2osfc       (begc:endc))                        ; this%frac_h2osfc      (:)   = nan
    allocate(this%h2osoi_liq        (begc:endc,-nlevsno+1:nlevgrnd))    ; this%h2osoi_liq       (:,:) = nan
    allocate(this%h2osoi_ice        (begc:endc,-nlevsno+1:nlevgrnd))    ; this%h2osoi_ice       (:,:) = nan
    allocate(this%h2osoi_vol        (begc:endc, 1:nlevgrnd))            ; this%h2osoi_vol       (:,:) = nan
    allocate(this%h2osfc            (begc:endc))                        ; this%h2osfc           (:)   = nan

    ! temperature_vars:
    allocate(this%t_soisno          (begc:endc,-nlevsno+1:nlevgrnd))    ; this%t_soisno         (:,:) = nan
    allocate(this%t_grnd            (begc:endc))                        ; this%t_grnd           (:)   = nan
    allocate(this%t_h2osfc          (begc:endc))                        ; this%t_h2osfc         (:)   = nan
    allocate(this%t_nearsurf        (begc:endc))                        ; this%t_nearsurf       (:)   = nan

    ! canopystate_vars:
    allocate(this%alt_indx          (begc:endc))                        ; this%alt_indx         (:)   = huge(1)

    !------------------------------------------------------------------------------------------
    ! pflotran variables: BEGIN
    !------------------------------------------------------------------------------------------
    ! waterflux_vars:
    allocate(this%qflx_top_soil     (begc:endc))                        ; this%qflx_top_soil     (:)   = ival
    allocate(this%qflx_evap_h2osfc  (begc:endc))                        ; this%qflx_evap_h2osfc  (:)   = ival
    allocate(this%qflx_evap_soil    (begc:endc))                        ; this%qflx_evap_soil    (:)   = ival
    allocate(this%qflx_evap_snow    (begc:endc))                        ; this%qflx_evap_snow    (:)   = ival
    allocate(this%qflx_subl_snow    (begc:endc))                        ; this%qflx_subl_snow    (:)   = ival
    allocate(this%qflx_tran_veg     (begc:endc))                        ; this%qflx_tran_veg     (:)   = ival
    allocate(this%qflx_rootsoil     (begc:endc,1:nlevgrnd))             ; this%qflx_rootsoil     (:,:) = ival
    allocate(this%qflx_et_reduced   (begc:endc))                        ; this%qflx_et_reduced   (:)   = ival
    allocate(this%qflx_infl         (begc:endc))                        ; this%qflx_infl         (:)   = ival
    allocate(this%qflx_exfl         (begc:endc))                        ; this%qflx_exfl         (:)   = ival
    allocate(this%qflx_surf         (begc:endc))                        ; this%qflx_surf         (:)   = ival
    allocate(this%qflx_base         (begc:endc))                        ; this%qflx_base         (:)   = ival
    allocate(this%qflx_drain        (begc:endc))                        ; this%qflx_drain        (:)   = ival
    allocate(this%qflx_drain_vr     (begc:endc,1:nlevgrnd))             ; this%qflx_drain_vr     (:,:) = ival
    allocate(this%qflx_lateral      (begc:endc))                        ; this%qflx_lateral      (:)   = ival
    allocate(this%qflx_lateral_vr   (begc:endc,1:nlevgrnd))             ; this%qflx_lateral_vr   (:,:) = ival
    allocate(this%qflx_h2osfc       (begc:endc))                        ; this%qflx_h2osfc       (:)   = ival


    ! energyflux_vars:
    allocate( this%htvp             (begc:endc))                        ; this%htvp             (:)   = ival
    allocate( this%eflx_bot         (begc:endc))                        ; this%eflx_bot         (:)   = ival
    allocate( this%eflx_soil_grnd   (begc:endc))                        ; this%eflx_soil_grnd   (:)   = ival
    allocate( this%eflx_fgr0_soil   (begc:endc))                        ; this%eflx_fgr0_soil   (:)   = ival
    allocate( this%eflx_fgr0_snow   (begc:endc))                        ; this%eflx_fgr0_snow   (:)   = ival
    allocate( this%eflx_fgr0_h2osfc (begc:endc))                        ; this%eflx_fgr0_h2osfc (:)   = ival
    allocate( this%eflx_rnet_soil   (begc:endc))                        ; this%eflx_rnet_soil   (:)   = ival

    ! atm2lnd:
    allocate(this%forc_pbot         (begc:endc))                        ; this%forc_pbot        (:)   = ival

    ! soilhydrology_vars:
    allocate( this%frost_table      (begc:endc))                        ; this%frost_table      (:)   = ival
    allocate( this%zwt              (begc:endc))                        ; this%zwt              (:)   = ival
    allocate( this%zwt_perched      (begc:endc))                        ; this%zwt_perched      (:)   = ival
    allocate( this%qcharge          (begc:endc))                        ; this%qcharge          (:)   = ival

  end subroutine InitAllocate
!-------------------------------------------------------------------------------------------------

end module clm_interface_thType
