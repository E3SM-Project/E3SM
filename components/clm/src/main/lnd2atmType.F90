module lnd2atmType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd, lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use shr_megan_mod , only : shr_megan_mechcomps_n
  use clm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon    , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl    , only : iulog, use_c13, use_cn, use_lch4
  use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND
  use decompMod     , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! ----------------------------------------------------
  ! land -> atmosphere variables structure
  !----------------------------------------------------
  type, public :: lnd2atm_type

     ! lnd->atm
     real(r8), pointer :: t_rad_grc          (:)   => null() ! radiative temperature (Kelvin)
     real(r8), pointer :: t_ref2m_grc        (:)   => null() ! 2m surface air temperature (Kelvin)
     real(r8), pointer :: q_ref2m_grc        (:)   => null() ! 2m surface specific humidity (kg/kg)
     real(r8), pointer :: u_ref10m_grc       (:)   => null() ! 10m surface wind speed (m/sec)
     real(r8), pointer :: h2osno_grc         (:)   => null() ! snow water (mm H2O)
     real(r8), pointer :: h2osoi_vol_grc     (:,:) => null() ! volumetric soil water (0~watsat, m3/m3, nlevgrnd) (for dust model)
     real(r8), pointer :: albd_grc           (:,:) => null() ! (numrad) surface albedo (direct)
     real(r8), pointer :: albi_grc           (:,:) => null() ! (numrad) surface albedo (diffuse)
     real(r8), pointer :: taux_grc           (:)   => null() ! wind stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy_grc           (:)   => null() ! wind stress: n-s (kg/m/s**2)
     real(r8), pointer :: eflx_lh_tot_grc    (:)   => null() ! total latent HF (W/m**2)  [+ to atm]
     real(r8), pointer :: eflx_sh_tot_grc    (:)   => null() ! total sensible HF (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lwrad_out_grc (:)   => null() ! IR (longwave) radiation (W/m**2)
     real(r8), pointer :: qflx_evap_tot_grc  (:)   => null() ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
     real(r8), pointer :: fsa_grc            (:)   => null() ! solar rad absorbed (total) (W/m**2)
     real(r8), pointer :: nee_grc            (:)   => null() ! net CO2 flux (kg CO2/m**2/s) [+ to atm]
     real(r8), pointer :: nem_grc            (:)   => null() ! gridcell average net methane correction to CO2 flux (g C/m^2/s)
     real(r8), pointer :: ram1_grc           (:)   => null() ! aerodynamical resistance (s/m)
     real(r8), pointer :: fv_grc             (:)   => null() ! friction velocity (m/s) (for dust model)
     real(r8), pointer :: flxdst_grc         (:,:) => null() ! dust flux (size bins)
     real(r8), pointer :: ddvel_grc          (:,:) => null() ! dry deposition velocities
     real(r8), pointer :: flxvoc_grc         (:,:) => null() ! VOC flux (size bins)
     real(r8), pointer :: flux_ch4_grc       (:)   => null() ! net CH4 flux (kg C/m**2/s) [+ to atm]
     ! lnd->rof
     real(r8), pointer :: qflx_rofliq_grc      (:) => null() ! rof liq forcing
     real(r8), pointer :: qflx_rofliq_qsur_grc (:) => null() ! rof liq -- surface runoff component
     real(r8), pointer :: qflx_rofliq_qsurp_grc(:) => null() ! rof liq -- surface ponding runoff component
     real(r8), pointer :: qflx_rofliq_qsub_grc (:) => null() ! rof liq -- subsurface runoff component
     real(r8), pointer :: qflx_rofliq_qsubp_grc(:) => null() ! rof liq -- perched subsurface runoff component
     real(r8), pointer :: qflx_rofliq_qgwl_grc (:) => null() ! rof liq -- glacier, wetland and lakes water balance residual component
     real(r8), pointer :: qflx_rofice_grc      (:) => null() ! rof ice forcing

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  

  end type lnd2atm_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(lnd2atm_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize lnd2atm derived type
    !
    ! !ARGUMENTS:
    class (lnd2atm_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer  :: begg, endg
    !------------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg

    allocate(this%t_rad_grc            (begg:endg))            ; this%t_rad_grc            (:) =ival
    allocate(this%t_ref2m_grc          (begg:endg))            ; this%t_ref2m_grc          (:) =ival
    allocate(this%q_ref2m_grc          (begg:endg))            ; this%q_ref2m_grc          (:) =ival
    allocate(this%u_ref10m_grc         (begg:endg))            ; this%u_ref10m_grc         (:) =ival
    allocate(this%h2osno_grc           (begg:endg))            ; this%h2osno_grc           (:) =ival
    allocate(this%h2osoi_vol_grc       (begg:endg,1:nlevgrnd)) ; this%h2osoi_vol_grc     (:,:) =ival
    allocate(this%albd_grc             (begg:endg,1:numrad))   ; this%albd_grc           (:,:) =ival
    allocate(this%albi_grc             (begg:endg,1:numrad))   ; this%albi_grc           (:,:) =ival
    allocate(this%taux_grc             (begg:endg))            ; this%taux_grc             (:) =ival
    allocate(this%tauy_grc             (begg:endg))            ; this%tauy_grc             (:) =ival
    allocate(this%eflx_lwrad_out_grc   (begg:endg))            ; this%eflx_lwrad_out_grc   (:) =ival
    allocate(this%eflx_sh_tot_grc      (begg:endg))            ; this%eflx_sh_tot_grc      (:) =ival
    allocate(this%eflx_lh_tot_grc      (begg:endg))            ; this%eflx_lh_tot_grc      (:) =ival
    allocate(this%qflx_evap_tot_grc    (begg:endg))            ; this%qflx_evap_tot_grc    (:) =ival
    allocate(this%fsa_grc              (begg:endg))            ; this%fsa_grc              (:) =ival
    allocate(this%nee_grc              (begg:endg))            ; this%nee_grc              (:) =ival
    allocate(this%nem_grc              (begg:endg))            ; this%nem_grc              (:) =ival
    allocate(this%ram1_grc             (begg:endg))            ; this%ram1_grc             (:) =ival
    allocate(this%fv_grc               (begg:endg))            ; this%fv_grc               (:) =ival
    allocate(this%flxdst_grc           (begg:endg,1:ndst))     ; this%flxdst_grc         (:,:) =ival
    allocate(this%flux_ch4_grc         (begg:endg))            ; this%flux_ch4_grc         (:) =ival
    allocate(this%qflx_rofliq_grc      (begg:endg))            ; this%qflx_rofliq_grc      (:) =ival
    allocate(this%qflx_rofliq_qsur_grc (begg:endg))            ; this%qflx_rofliq_qsur_grc (:) =ival
    allocate(this%qflx_rofliq_qsurp_grc(begg:endg))            ; this%qflx_rofliq_qsurp_grc(:) =ival
    allocate(this%qflx_rofliq_qsub_grc (begg:endg))            ; this%qflx_rofliq_qsub_grc (:) =ival
    allocate(this%qflx_rofliq_qsubp_grc(begg:endg))            ; this%qflx_rofliq_qsubp_grc(:) =ival
    allocate(this%qflx_rofliq_qgwl_grc (begg:endg))            ; this%qflx_rofliq_qgwl_grc (:) =ival
    allocate(this%qflx_rofice_grc      (begg:endg))            ; this%qflx_rofice_grc      (:) =ival

    if (shr_megan_mechcomps_n>0) then
       allocate(this%flxvoc_grc(begg:endg,1:shr_megan_mechcomps_n));  this%flxvoc_grc(:,:)=ival
    endif
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(this%ddvel_grc(begg:endg,1:n_drydep));   this%ddvel_grc(:,:)=ival
    end if

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(lnd2atm_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg

    this%eflx_sh_tot_grc(begg:endg) = 0._r8
    call hist_addfld1d (fname='FSH', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_lnd=this%eflx_sh_tot_grc)
       
    this%qflx_rofliq_grc(begg:endg) = 0._r8
    call hist_addfld1d (fname='QRUNOFF',  units='mm/s',  &
         avgflag='A', long_name='total liquid runoff (does not include QSNWCPICE)', &
         ptr_lnd=this%qflx_rofliq_grc)

    this%qflx_rofice_grc(begg:endg) = 0._r8
    call hist_addfld1d (fname='QSNWCPICE',  units='mm/s',  &
         avgflag='A', long_name='excess snowfall due to snow capping', &
         ptr_lnd=this%qflx_rofice_grc)

    if (use_lch4) then
       this%flux_ch4_grc(begg:endg) = 0._r8
       call hist_addfld1d (fname='FCH4', units='kgC/m2/s', &
            avgflag='A', long_name='Gridcell surface CH4 flux to atmosphere (+ to atm)', &
            ptr_lnd=this%flux_ch4_grc)

       this%nem_grc(begg:endg) = spval
       call hist_addfld1d (fname='NEM', units='gC/m2/s', &
            avgflag='A', long_name='Gridcell net adjustment to NEE passed to atm. for methane production', &
            ptr_lnd=this%nem_grc)
    end if

  end subroutine InitHistory

end module lnd2atmType
