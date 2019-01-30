module VegetationDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Vegetation data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevurb, crop_prog
  use clm_varcon     , only : spval, ispval, sb
  use landunit_varcon, only : istsoil, istcrop
  use clm_varctl     , only : iulog, use_cn 
  use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ncdio_pio      , only : file_desc_t, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
  use decompMod      , only : bounds_type, get_proc_global
  use restUtilMod
  use VegetationType , only : veg_pp
  use LandunitType   , only : lun_pp
  use GridcellType   , only : grc_pp
  use ColumnDataType , only : col_es
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_state
    ! temperature variables
    real(r8), pointer :: t_veg              (:) => null() ! vegetation temperature (K)
    real(r8), pointer :: t_ref2m            (:) => null() ! 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_r          (:) => null() ! rural 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_u          (:) => null() ! urban 2 m height surface air temperature (K)
    ! temperature summary and accumulator variables
    real(r8), pointer :: t_a10              (:) => null() ! 10-day running mean of the 2 m temperature (K)
    real(r8), pointer :: t_a10min           (:) => null() ! 10-day running mean of min 2-m temperature
    real(r8), pointer :: t_a5min            (:) => null() ! 5-day running mean of min 2-m temperature
    real(r8), pointer :: t_ref2m_min        (:) => null() ! daily minimum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_min_r      (:) => null() ! daily minimum of average 2 m height surface air temperature - rural(K)
    real(r8), pointer :: t_ref2m_min_u      (:) => null() ! daily minimum of average 2 m height surface air temperature - urban (K)
    real(r8), pointer :: t_ref2m_max        (:) => null() ! daily maximum of average 2 m height surface air temperature (K)
    real(r8), pointer :: t_ref2m_max_r      (:) => null() ! daily maximum of average 2 m height surface air temperature - rural(K)
    real(r8), pointer :: t_ref2m_max_u      (:) => null() ! daily maximum of average 2 m height surface air temperature - urban (K)
    real(r8), pointer :: t_ref2m_min_inst   (:) => null() ! instantaneous daily min of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_min_inst_r (:) => null() ! instantaneous daily min of average 2 m height surface air temp - rural (K)
    real(r8), pointer :: t_ref2m_min_inst_u (:) => null() ! instantaneous daily min of average 2 m height surface air temp - urban (K)
    real(r8), pointer :: t_ref2m_max_inst   (:) => null() ! instantaneous daily max of average 2 m height surface air temp (K)
    real(r8), pointer :: t_ref2m_max_inst_r (:) => null() ! instantaneous daily max of average 2 m height surface air temp - rural (K)
    real(r8), pointer :: t_ref2m_max_inst_u (:) => null() ! instantaneous daily max of average 2 m height surface air temp - urban (K)
    real(r8), pointer :: t_veg24            (:) => null() ! 24hr average vegetation temperature (K)
    real(r8), pointer :: t_veg240           (:) => null() ! 240hr average vegetation temperature (K)
    real(r8), pointer :: gdd0               (:) => null() ! growing degree-days base  0C from planting  (ddays)
    real(r8), pointer :: gdd8               (:) => null() ! growing degree-days base  8C from planting  (ddays)
    real(r8), pointer :: gdd10              (:) => null() ! growing degree-days base 10C from planting  (ddays)
    real(r8), pointer :: gdd020             (:) => null() ! 20-year average of gdd0                     (ddays)
    real(r8), pointer :: gdd820             (:) => null() ! 20-year average of gdd8                     (ddays)
    real(r8), pointer :: gdd1020            (:) => null() ! 20-year average of gdd10                    (ddays)
    ! temperature-related variables
    real(r8), pointer :: thm                (:) => null() ! intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
    real(r8), pointer :: emv                (:) => null() ! vegetation emissivity (unitless)

    contains
    procedure, public :: Init    => veg_es_init
    procedure, public :: Restart => veg_es_restart
    procedure, public :: Clean   => veg_es_clean
    procedure, public :: InitAccBuffer => init_acc_buffer_veg_es
    procedure, public :: InitAccVars   => init_acc_vars_veg_es
    procedure, public :: UpdateAccVars => update_acc_vars_veg_es
  end type vegetation_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_state
    ! Note: All units given as kg in this data type imply kg of H2O
    real(r8), pointer :: h2ocan       (:) => null() ! canopy water (kg/m2)
    real(r8), pointer :: q_ref2m      (:) => null() ! 2 m height surface specific humidity (kg H2O/kg moist air)
    real(r8), pointer :: rh_ref2m     (:) => null() ! 2 m height surface relative humidity (%)
    real(r8), pointer :: rh_ref2m_r   (:) => null() ! 2 m height surface relative humidity - rural (%)
    real(r8), pointer :: rh_ref2m_u   (:) => null() ! 2 m height surface relative humidity - urban (%)
    real(r8), pointer :: rh_af        (:) => null() ! fractional humidity of canopy air (dimensionless)
    real(r8), pointer :: fwet         (:) => null() ! canopy fraction that is wet (0 to 1)
    real(r8), pointer :: fdry         (:) => null() ! canopy fraction of foliage that is green and dry (0 to 1)
    real(r8), pointer :: begwb        (:) => null() ! water mass begining of the time step
    real(r8), pointer :: endwb        (:) => null() ! water mass end of the time step
    real(r8), pointer :: errh2o       (:) => null() ! water conservation error (mm H2O)
  contains
    procedure, public :: Init    => veg_ws_init
    procedure, public :: Restart => veg_ws_restart
    procedure, public :: Clean   => veg_ws_clean
  end type vegetation_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => veg_cs_init
    procedure, public :: Clean => veg_cs_init
  end type vegetation_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => veg_ns_init
    procedure, public :: Clean => veg_ns_clean
  end type vegetation_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => veg_ps_init
    procedure, public :: Clean => veg_ps_clean
  end type vegetation_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_flux
    real(r8), pointer :: eflx_sh_grnd      (:) => null() ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_veg       (:) => null() ! sensible heat flux from leaves (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_snow      (:) => null() ! sensible heat flux from snow (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_soil      (:) => null() ! sensible heat flux from soil  (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_h2osfc    (:) => null() ! sensible heat flux from surface water (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot       (:) => null() ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_u     (:) => null() ! urban total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_r     (:) => null() ! rural total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot       (:) => null() ! total latent heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot_u     (:) => null() ! urban total latent heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot_r     (:) => null() ! rural total latent heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vegt      (:) => null() ! transpiration heat flux from veg (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vege      (:) => null() ! evaporation heat flux from veg (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_grnd      (:) => null() ! evaporation heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_soil_grnd    (:) => null() ! soil heat flux (W/m**2) [+ = into soil] 
    real(r8), pointer :: eflx_soil_grnd_u  (:) => null() ! urban soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_soil_grnd_r  (:) => null() ! rural soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_lwrad_net    (:) => null() ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_net_r  (:) => null() ! rural net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_net_u  (:) => null() ! urban net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_out    (:) => null() ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_out_r  (:) => null() ! rural emitted infrared (longwave) rad (W/m**2)
    real(r8), pointer :: eflx_lwrad_out_u  (:) => null() ! urban emitted infrared (longwave) rad (W/m**2)
    real(r8), pointer :: eflx_gnet         (:) => null() ! net heat flux into ground  (W/m**2)
    real(r8), pointer :: eflx_grnd_lake    (:) => null() ! net heat flux into lake / snow surface, excluding light transmission (W/m**2)
    real(r8), pointer :: eflx_anthro       (:) => null() ! total anthropogenic heat flux (W/m**2)
    real(r8), pointer :: eflx_traffic      (:) => null() ! traffic sensible heat flux (W/m2)
    real(r8), pointer :: eflx_wasteheat    (:) => null() ! sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac (:) => null() ! sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: dlrad             (:) => null() ! downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad             (:) => null() ! upward longwave radiation above the canopy [W/m2]
    ! Wind Stress
    real(r8), pointer :: taux              (:) => null() ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy              (:) => null() ! wind (shear) stress: n-s (kg/m/s**2)
    ! Derivatives of energy fluxes
    real(r8), pointer :: dgnetdT           (:) => null() ! derivative of net ground heat flux wrt soil temp  (W/m**2 K)
    real(r8), pointer :: netrad            (:) => null() ! col net radiation (W/m**2) [+ = to sfc]
    real(r8), pointer :: cgrnd             (:) => null() ! col deriv. of soil energy flux wrt to soil temp [W/m2/k]
    real(r8), pointer :: cgrndl            (:) => null() ! col deriv. of soil latent heat flux wrt soil temp  [W/m**2/k]
    real(r8), pointer :: cgrnds            (:) => null() ! col deriv. of soil sensible heat flux wrt soil temp [W/m2/k]
    ! Balance Checks
    real(r8), pointer :: errsoi            (:) => null() ! soil/lake energy conservation error   (W/m**2)
    real(r8), pointer :: errseb            (:) => null() ! surface energy conservation error     (W/m**2)
    real(r8), pointer :: errsol            (:) => null() ! solar radiation conservation error    (W/m**2)
    real(r8), pointer :: errlon            (:) => null() ! longwave radiation conservation error (W/m**2)

  contains
    procedure, public :: Init    => veg_ef_init
    procedure, public :: Restart => veg_ef_restart
    procedure, public :: Clean   => veg_ef_clean
  end type vegetation_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_flux
    real(r8), pointer :: qflx_prec_grnd     (:)   => null() ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd     (:)   => null() ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd     (:)   => null() ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow      (:)   => null() ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_soi      (:)   => null() ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_veg      (:)   => null() ! vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can      (:)   => null() ! evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot      (:)   => null() ! pft_qflx_evap_soi + pft_qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd     (:)   => null() ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq     (:)   => null() ! excess rainfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_ice     (:)   => null() ! excess snowfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_tran_veg      (:)   => null() ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_dew_snow      (:)   => null() ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd      (:)   => null() ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_prec_intr     (:)   => null() ! interception of precipitation [mm/s]
    real(r8), pointer :: qflx_ev_snow       (:)   => null() ! evaporation heat flux from snow       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_soil       (:)   => null() ! evaporation heat flux from soil       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_h2osfc     (:)   => null() ! evaporation heat flux from soil       (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_rootsoi_frac  (:,:) => null() !  
    real(r8), pointer :: qflx_irrig         (:)   => null() ! irrigation flux (mm H2O/s)
    real(r8), pointer :: irrig_rate         (:)   => null() ! current irrigation rate [mm/s]
    integer , pointer :: n_irrig_steps_left (:)   => null() ! number of time steps for which we still need to irrigate today (if 0, ignore)

  contains
    procedure, public :: Init    => veg_wf_init
    procedure, public :: Restart => veg_wf_restart
    procedure, public :: Clean   => veg_wf_clean
  end type vegetation_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => veg_cf_init
    procedure, public :: Clean => veg_cf_clean
  end type vegetation_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => veg_nf_init
    procedure, public :: Clean => veg_nf_clean
  end type vegetation_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => veg_pf_init
    procedure, public :: Clean => veg_pf_clean
  end type vegetation_phosphorus_flux

   
  !-----------------------------------------------------------------------
  ! declare the public instances of vegetation-level data types
  !-----------------------------------------------------------------------
  type(vegetation_energy_state)          , public, target :: veg_es    ! vegetation energy state
  type(vegetation_water_state)           , public, target :: veg_ws    ! vegetation water state
  type(vegetation_carbon_state)          , public, target :: veg_cs    ! vegetation carbon state
  type(vegetation_nitrogen_state)        , public, target :: veg_ns    ! vegetation nitrogen state
  type(vegetation_phosphorus_state)      , public, target :: veg_ps    ! vegetation phosphorus state
  type(vegetation_energy_flux)           , public, target :: veg_ef    ! vegetation energy flux
  type(vegetation_water_flux)            , public, target :: veg_wf    ! vegetation water flux
  type(vegetation_carbon_flux)           , public, target :: veg_cf    ! vegetation carbon flux
  type(vegetation_nitrogen_flux)         , public, target :: veg_nf    ! vegetation nitrogen flux
  type(vegetation_phosphorus_flux)       , public, target :: veg_pf    ! vegetation phosphorus flux

  !------------------------------------------------------------------------
  
  contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy state data structure
  !------------------------------------------------------------------------
  subroutine veg_es_init(this, begp, endp)
    !
    ! !USES:
    use clm_varctl     , only : use_vancouver, use_mexicocity
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,j                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_es
    !-----------------------------------------------------------------------
    allocate(this%t_veg              (begp:endp))                   ; this%t_veg              (:)   = nan
    allocate(this%t_ref2m            (begp:endp))                   ; this%t_ref2m            (:)   = nan
    allocate(this%t_ref2m_r          (begp:endp))                   ; this%t_ref2m_r          (:)   = nan
    allocate(this%t_ref2m_u          (begp:endp))                   ; this%t_ref2m_u          (:)   = nan
    allocate(this%t_a10              (begp:endp))                   ; this%t_a10              (:)   = nan
    allocate(this%t_a10min           (begp:endp))                   ; this%t_a10min           (:)   = nan
    allocate(this%t_a5min            (begp:endp))                   ; this%t_a5min            (:)   = nan
    allocate(this%t_ref2m_min        (begp:endp))                   ; this%t_ref2m_min        (:)   = nan
    allocate(this%t_ref2m_min_r      (begp:endp))                   ; this%t_ref2m_min_r      (:)   = nan
    allocate(this%t_ref2m_min_u      (begp:endp))                   ; this%t_ref2m_min_u      (:)   = nan
    allocate(this%t_ref2m_max        (begp:endp))                   ; this%t_ref2m_max        (:)   = nan
    allocate(this%t_ref2m_max_r      (begp:endp))                   ; this%t_ref2m_max_r      (:)   = nan
    allocate(this%t_ref2m_max_u      (begp:endp))                   ; this%t_ref2m_max_u      (:)   = nan
    allocate(this%t_ref2m_min_inst   (begp:endp))                   ; this%t_ref2m_min_inst   (:)   = nan
    allocate(this%t_ref2m_min_inst_r (begp:endp))                   ; this%t_ref2m_min_inst_r (:)   = nan
    allocate(this%t_ref2m_min_inst_u (begp:endp))                   ; this%t_ref2m_min_inst_u (:)   = nan
    allocate(this%t_ref2m_max_inst   (begp:endp))                   ; this%t_ref2m_max_inst   (:)   = nan
    allocate(this%t_ref2m_max_inst_r (begp:endp))                   ; this%t_ref2m_max_inst_r (:)   = nan
    allocate(this%t_ref2m_max_inst_u (begp:endp))                   ; this%t_ref2m_max_inst_u (:)   = nan
    allocate(this%t_veg24            (begp:endp))                   ; this%t_veg24            (:)   = nan
    allocate(this%t_veg240           (begp:endp))                   ; this%t_veg240           (:)   = nan
    allocate(this%gdd0               (begp:endp))                   ; this%gdd0               (:)   = spval
    allocate(this%gdd8               (begp:endp))                   ; this%gdd8               (:)   = spval
    allocate(this%gdd10              (begp:endp))                   ; this%gdd10              (:)   = spval
    allocate(this%gdd020             (begp:endp))                   ; this%gdd020             (:)   = spval
    allocate(this%gdd820             (begp:endp))                   ; this%gdd820             (:)   = spval
    allocate(this%gdd1020            (begp:endp))                   ; this%gdd1020            (:)   = spval
    allocate(this%thm                (begp:endp))                   ; this%thm                (:)   = nan
    allocate(this%emv                (begp:endp))                   ; this%emv                (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_es
    !-----------------------------------------------------------------------
    this%t_veg(begp:endp) = spval
    call hist_addfld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_patch=this%t_veg)

    this%t_ref2m(begp:endp) = spval
    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_patch=this%t_ref2m)

    this%t_ref2m_r(begp:endp) = spval
    call hist_addfld1d (fname='TSA_R', units='K',  &
         avgflag='A', long_name='Rural 2m air temperature', &
         ptr_patch=this%t_ref2m_r, set_spec=spval)

    this%t_ref2m_u(begp:endp) = spval
    call hist_addfld1d (fname='TSA_U', units='K',  &
         avgflag='A', long_name='Urban 2m air temperature', &
         ptr_patch=this%t_ref2m_u, set_nourb=spval)

    this%t_a10(begp:endp) = spval
    call hist_addfld1d (fname='T10', units='K',  &
         avgflag='A', long_name='10-day running mean of 2-m temperature', &
         ptr_patch=this%t_a10, default='inactive')

    if (use_cn .and. crop_prog )then
       this%t_a10min(begp:endp) = spval
       call hist_addfld1d (fname='A10TMIN', units='K',  &
            avgflag='A', long_name='10-day running mean of min 2-m temperature', &
            ptr_patch=this%t_a10min, default='inactive')
    end if

    if (use_cn .and.  crop_prog )then
       this%t_a5min(begp:endp) = spval
       call hist_addfld1d (fname='A5TMIN', units='K',  &
            avgflag='A', long_name='5-day running mean of min 2-m temperature', &
            ptr_patch=this%t_a5min, default='inactive')
    end if

    this%t_ref2m_min(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min)

    this%t_ref2m_max(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max)

    this%t_ref2m_min_r(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV_R', units='K',  &
         avgflag='A', long_name='Rural daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_r, set_spec=spval)

    this%t_ref2m_max_r(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV_R', units='K',  &
         avgflag='A', long_name='Rural daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_r, set_spec=spval)

    this%t_ref2m_min_u(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV_U', units='K',  &
         avgflag='A', long_name='Urban daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_u, set_nourb=spval)

    this%t_ref2m_max_u(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV_U', units='K',  &
         avgflag='A', long_name='Urban daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_u, set_nourb=spval)

    this%t_veg24(begp:endp) = spval
    call hist_addfld1d (fname='TV24', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 24hrs)', &
         ptr_patch=this%t_veg24, default='inactive')

    this%t_veg240(begp:endp)  = spval
    call hist_addfld1d (fname='TV240', units='K',  &
         avgflag='A', long_name='vegetation temperature (last 240hrs)', &
         ptr_patch=this%t_veg240, default='inactive')

    if (crop_prog) then
       this%gdd0(begp:endp) = spval
       call hist_addfld1d (fname='GDD0', units='ddays', &
            avgflag='A', long_name='Growing degree days base  0C from planting', &
            ptr_patch=this%gdd0, default='inactive')

       this%gdd8(begp:endp) = spval
       call hist_addfld1d (fname='GDD8', units='ddays', &
            avgflag='A', long_name='Growing degree days base  8C from planting', &
            ptr_patch=this%gdd8, default='inactive')

       this%gdd10(begp:endp) = spval
       call hist_addfld1d (fname='GDD10', units='ddays', &
            avgflag='A', long_name='Growing degree days base 10C from planting', &
            ptr_patch=this%gdd10, default='inactive')

       this%gdd020(begp:endp) = spval
       call hist_addfld1d (fname='GDD020', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base  0C from planting', &
            ptr_patch=this%gdd020, default='inactive')

       this%gdd820(begp:endp) = spval
       call hist_addfld1d (fname='GDD820', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base  8C from planting', &
            ptr_patch=this%gdd820, default='inactive')

       this%gdd1020(begp:endp) = spval
       call hist_addfld1d (fname='GDD1020', units='ddays', &
            avgflag='A', long_name='Twenty year average of growing degree days base 10C from planting', &
            ptr_patch=this%gdd1020, default='inactive')
    end if
    
    if (use_cn ) then
       this%emv(begp:endp) = spval
       call hist_addfld1d (fname='EMV', units='proportion', &
            avgflag='A', long_name='vegetation emissivity', &
            ptr_patch=this%emv, default='inactive')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_es
    !-----------------------------------------------------------------------
    
    ! Set t_veg, t_ref2m, t_ref2m_u and tref2m_r 

    do p = begp, endp
       c = veg_pp%column(p)
       l = veg_pp%landunit(p)

       if (use_vancouver) then
          this%t_veg(p)   = 297.56
       else if (use_mexicocity) then
          this%t_veg(p)   = 289.46
       else
          this%t_veg(p)   = 283._r8
       end if

       if (use_vancouver) then
          this%t_ref2m(p) = 297.56
       else if (use_mexicocity) then
          this%t_ref2m(p) = 289.46
       else
          this%t_ref2m(p) = 283._r8
       end if

       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%t_ref2m_u(p) = 297.56
          else if (use_mexicocity) then
             this%t_ref2m_u(p) = 289.46
          else
             this%t_ref2m_u(p) = 283._r8
          end if
       else 
          if (.not. lun_pp%ifspecial(l)) then 
             if (use_vancouver) then
                this%t_ref2m_r(p) = 297.56
             else if (use_mexicocity) then
                this%t_ref2m_r(p) = 289.46
             else
                this%t_ref2m_r(p) = 283._r8
             end if
          else 
             this%t_ref2m_r(p) = spval
          end if
       end if

    end do ! veg loop

  end subroutine veg_es_init
    
  !------------------------------------------------------------------------
  subroutine veg_es_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='T_VEG', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='vegetation temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_veg)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Rural 2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_U', xtype=ncd_double, dim1name='pft',                      &
         long_name='Urban 2m height surface air temperature', units='K',                                              &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_U', xtype=ncd_double, dim1name='pft',                  &
         long_name='urban daily minimum of average 2 m height surface air temperature (K)', units='K',                &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_U', xtype=ncd_double, dim1name='pft',                  &
         long_name='urban daily maximum of average 2 m height surface air temperature (K)', units='K',                &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_U', xtype=ncd_double, dim1name='pft',             &
         long_name='urban instantaneous daily min of average 2 m height surface air temp (K)', units='K',             &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_u)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_r)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_U', xtype=ncd_double,  dim1name='pft',            &
         long_name='urban instantaneous daily max of average 2 m height surface air temp (K)', units='K',             &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_u)

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='gdd1020', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 10C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd1020)

       call restartvar(ncid=ncid, flag=flag,  varname='gdd820', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 8C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd820)

       call restartvar(ncid=ncid, flag=flag,  varname='gdd020', xtype=ncd_double,  &
            dim1name='pft', long_name='20 year average of growing degree-days base 0C from planting', units='ddays', &
            interpinic_flag='interp', readvar=readvar, data=this%gdd020)
    end if


  end subroutine veg_es_restart

  !------------------------------------------------------------------------
  subroutine veg_es_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    !------------------------------------------------------------------------
  end subroutine veg_es_clean
  
  !-----------------------------------------------------------------------
  subroutine init_acc_buffer_veg_es (this, bounds)
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for accumulated fields for vegetation energy state
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use accumulMod       , only : init_accum_field
    use clm_time_manager , only : get_step_size
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES: 
    real(r8) :: dtime
    integer, parameter :: not_used = huge(1)
    !---------------------------------------------------------------------

    dtime = get_step_size()

    ! The following is a running mean. The accumulation period is set to -10 for a 10-day running mean.
    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20._r8)

    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='TREFAV_U', units='K', &
         desc='average over an hour of urban 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field(name='TREFAV_R', units='K', &
         desc='average over an hour of rural 2-m temperature', accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%t_veg24(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='T_VEG24', units='K',                                              &
         desc='24hr average of vegetation temperature',  accum_type='runmean', accum_period=-1,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%t_veg240(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='T_VEG240', units='K',                                             &
         desc='240hr average of vegetation temperature',  accum_type='runmean', accum_period=-10,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    if ( crop_prog )then
       call init_accum_field (name='TDM10', units='K', &
            desc='10-day running mean of min 2-m temperature', accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       call init_accum_field (name='TDM5', units='K', &
            desc='5-day running mean of min 2-m temperature', accum_type='runmean', accum_period=-5, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

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
    
  end subroutine init_acc_buffer_veg_es

  !-----------------------------------------------------------------------
  subroutine init_acc_vars_veg_es(this, bounds)
    ! !DESCRIPTION:
    ! Initialize variables associated with vegetation energy state
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    use clm_varctl       , only : nsrest, nsrStartup
    use abortutils       , only : endrun
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds  
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

    call extract_accum_field ('T10', rbufslp, nstep)
    this%t_a10(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('T_VEG24', rbufslp, nstep)
    this%t_veg24(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('T_VEG240', rbufslp, nstep)
    this%t_veg240(begp:endp) = rbufslp(begp:endp)

    if (crop_prog) then
       call extract_accum_field ('TDM10', rbufslp, nstep) 
       this%t_a10min(begp:endp)= rbufslp(begp:endp)

       call extract_accum_field ('TDM5', rbufslp, nstep) 
       this%t_a5min(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD0', rbufslp, nstep)
       this%gdd0(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD8', rbufslp, nstep) ;
       this%gdd8(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('GDD10', rbufslp, nstep) 
       this%gdd10(begp:endp) = rbufslp(begp:endp)

    end if

    ! Initialize variables that are to be time accumulated
    ! Initialize 2m ref temperature max and min values

    if (nsrest == nsrStartup) then 
       this%t_ref2m_max(begp:endp)        =  spval
       this%t_ref2m_max_r(begp:endp)      =  spval
       this%t_ref2m_max_u(begp:endp)      =  spval

       this%t_ref2m_min(begp:endp)        =  spval
       this%t_ref2m_min_r(begp:endp)      =  spval
       this%t_ref2m_min_u(begp:endp)      =  spval

       this%t_ref2m_max_inst(begp:endp)   = -spval
       this%t_ref2m_max_inst_r(begp:endp) = -spval
       this%t_ref2m_max_inst_u(begp:endp) = -spval

       this%t_ref2m_min_inst(begp:endp)   =  spval
       this%t_ref2m_min_inst_r(begp:endp) =  spval
       this%t_ref2m_min_inst_u(begp:endp) =  spval
    end if

    deallocate(rbufslp)
    
  end subroutine init_acc_vars_veg_es

  !-----------------------------------------------------------------------
  subroutine update_acc_vars_veg_es (this, bounds)
    !
    ! USES
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod       , only : update_accum_field, extract_accum_field, accumResetVal
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state)    :: this
    type(bounds_type)      , intent(in) :: bounds  
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
    
    ! fill the temporary variable
    do p = begp,endp
       rbufslp(p) = this%t_veg(p)
    end do
    
    call update_accum_field  ('T10', this%t_ref2m, nstep)
    call extract_accum_field ('T10', this%t_a10  , nstep)
    call update_accum_field  ('T_VEG24' , rbufslp       , nstep)
    call extract_accum_field ('T_VEG24' , this%t_veg24  , nstep)
    call update_accum_field  ('T_VEG240', rbufslp       , nstep)
    call extract_accum_field ('T_VEG240', this%t_veg240 , nstep)


    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', this%t_ref2m, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)
    end_cd = is_end_curr_day()
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst(p) = max(rbufslp(p), this%t_ref2m_max_inst(p))
          this%t_ref2m_min_inst(p) = min(rbufslp(p), this%t_ref2m_min_inst(p))
       endif
       if (end_cd) then
          this%t_ref2m_max(p) = this%t_ref2m_max_inst(p)
          this%t_ref2m_min(p) = this%t_ref2m_min_inst(p)
          this%t_ref2m_max_inst(p) = -spval
          this%t_ref2m_min_inst(p) =  spval
       else if (secs == int(dtime)) then
          this%t_ref2m_max(p) = spval
          this%t_ref2m_min(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_U - hourly average urban 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_U', this%t_ref2m_u, nstep)
    call extract_accum_field ('TREFAV_U', rbufslp, nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst_u(p) = max(rbufslp(p), this%t_ref2m_max_inst_u(p))
          this%t_ref2m_min_inst_u(p) = min(rbufslp(p), this%t_ref2m_min_inst_u(p))
       endif
       if (end_cd) then
         if (lun_pp%urbpoi(l)) then
          this%t_ref2m_max_u(p) = this%t_ref2m_max_inst_u(p)
          this%t_ref2m_min_u(p) = this%t_ref2m_min_inst_u(p)
          this%t_ref2m_max_inst_u(p) = -spval
          this%t_ref2m_min_inst_u(p) =  spval
         end if
       else if (secs == int(dtime)) then
          this%t_ref2m_max_u(p) = spval
          this%t_ref2m_min_u(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_R - hourly average rural 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_R', this%t_ref2m_r, nstep)
    call extract_accum_field ('TREFAV_R', rbufslp, nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          this%t_ref2m_max_inst_r(p) = max(rbufslp(p), this%t_ref2m_max_inst_r(p))
          this%t_ref2m_min_inst_r(p) = min(rbufslp(p), this%t_ref2m_min_inst_r(p))
       endif
       if (end_cd) then
         if (.not.(lun_pp%ifspecial(l))) then
          this%t_ref2m_max_r(p) = this%t_ref2m_max_inst_r(p)
          this%t_ref2m_min_r(p) = this%t_ref2m_min_inst_r(p)
          this%t_ref2m_max_inst_r(p) = -spval
          this%t_ref2m_min_inst_r(p) =  spval
         end if
       else if (secs == int(dtime)) then
          this%t_ref2m_max_r(p) = spval
          this%t_ref2m_min_r(p) = spval
       endif
    end do

    if ( crop_prog )then
       ! Accumulate and extract TDM10

       do p = begp,endp
          rbufslp(p) = min(this%t_ref2m_min(p),this%t_ref2m_min_inst(p)) 
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                                     !'min_inst' not initialized?
       call update_accum_field  ('TDM10', rbufslp, nstep)
       call extract_accum_field ('TDM10', this%t_a10min, nstep)

       ! Accumulate and extract TDM5

       do p = begp,endp
          rbufslp(p) = min(this%t_ref2m_min(p),this%t_ref2m_min_inst(p)) 
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                         !'min_inst' not initialized?
       call update_accum_field  ('TDM5', rbufslp, nstep)
       call extract_accum_field ('TDM5', this%t_a5min, nstep)

       ! Accumulate and extract GDD0

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(26._r8, this%t_ref2m(p)-SHR_CONST_TKFRZ)) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD0', rbufslp, nstep)
       call extract_accum_field ('GDD0', this%gdd0, nstep)

       ! Accumulate and extract GDD8

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  this%t_ref2m(p)-(SHR_CONST_TKFRZ + 8._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD8', rbufslp, nstep)
       call extract_accum_field ('GDD8', this%gdd8, nstep)

       ! Accumulate and extract GDD10

       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  this%t_ref2m(p)-(SHR_CONST_TKFRZ + 10._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD10', rbufslp, nstep)
       call extract_accum_field ('GDD10', this%gdd10, nstep)
    end if

    deallocate(rbufslp)

  end subroutine update_acc_vars_veg_es

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water state data structure
  !------------------------------------------------------------------------
  subroutine veg_ws_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer             :: p
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of veg_ws
    !-----------------------------------------------------------------------
    allocate(this%h2ocan              (begp:endp))          ; this%h2ocan            (:) = nan
    allocate(this%q_ref2m             (begp:endp))          ; this%q_ref2m           (:) = nan
    allocate(this%rh_ref2m            (begp:endp))          ; this%rh_ref2m          (:) = nan
    allocate(this%rh_ref2m_r          (begp:endp))          ; this%rh_ref2m_r        (:) = nan
    allocate(this%rh_ref2m_u          (begp:endp))          ; this%rh_ref2m_u        (:) = nan
    allocate(this%rh_af               (begp:endp))          ; this%rh_af             (:) = nan
    allocate(this%fwet                (begp:endp))          ; this%fwet              (:) = nan
    allocate(this%fdry                (begp:endp))          ; this%fdry              (:) = nan
    allocate(this%begwb               (begp:endp))          ; this%begwb             (:) = nan
    allocate(this%endwb               (begp:endp))          ; this%endwb             (:) = nan
    allocate(this%errh2o              (begp:endp))          ; this%errh2o            (:) = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_ws
    !-----------------------------------------------------------------------
    this%h2ocan(begp:endp) = spval 
    call hist_addfld1d (fname='H2OCAN', units='mm',  &
         avgflag='A', long_name='intercepted water', &
         ptr_patch=this%h2ocan, set_lake=0._r8)

    this%q_ref2m(begp:endp) = spval
    call hist_addfld1d (fname='Q2M', units='kg/kg',  &
         avgflag='A', long_name='2m specific humidity', &
         ptr_patch=this%q_ref2m)

    this%rh_ref2m(begp:endp) = spval
    call hist_addfld1d (fname='RH2M', units='%',  &
         avgflag='A', long_name='2m relative humidity', &
         ptr_patch=this%rh_ref2m)

    this%rh_ref2m_r(begp:endp) = spval
    call hist_addfld1d (fname='RH2M_R', units='%',  &
         avgflag='A', long_name='Rural 2m specific humidity', &
         ptr_patch=this%rh_ref2m_r, set_spec=spval)

    this%rh_ref2m_u(begp:endp) = spval
    call hist_addfld1d (fname='RH2M_U', units='%',  &
         avgflag='A', long_name='Urban 2m relative humidity', &
         ptr_patch=this%rh_ref2m_u, set_nourb=spval)

    this%rh_af(begp:endp) = spval
    call hist_addfld1d (fname='RHAF', units='fraction', &
         avgflag='A', long_name='fractional humidity of canopy air', &
         ptr_patch=this%rh_af, set_spec=spval, default='inactive')

    if (use_cn) then
       this%fwet(begp:endp) = spval
       call hist_addfld1d (fname='FWET', units='proportion', &
            avgflag='A', long_name='fraction of canopy that is wet', &
            ptr_patch=this%fwet, default='inactive')
    end if

    if (use_cn) then
       this%fdry(begp:endp) = spval
       call hist_addfld1d (fname='FDRY', units='proportion', &
            avgflag='A', long_name='fraction of foliage that is green and dry', &
            ptr_patch=this%fdry, default='inactive')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_ws
    !-----------------------------------------------------------------------
    do p = begp,endp
       this%h2ocan(begp:endp) = 0._r8
       this%fwet(begp:endp)   = 0._r8
       this%fdry(begp:endp)   = 0._r8
    end do
    
  end subroutine veg_ws_init
    
  !------------------------------------------------------------------------
  subroutine veg_ws_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_water_state)    :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='H2OCAN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='canopy water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2ocan)

    call restartvar(ncid=ncid, flag=flag, varname='FWET', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='fraction of canopy that is wet (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fwet)
   
  end subroutine veg_ws_restart
  
  !------------------------------------------------------------------------
  subroutine veg_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ws_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon state data structure
  !------------------------------------------------------------------------
  subroutine veg_cs_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,j                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_cs
    !-----------------------------------------------------------------------
    allocate(this%xxx              (begp:endp))                   ; this%xxx              (:)   = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_cs
    !-----------------------------------------------------------------------
  end subroutine veg_cs_init
    
  !------------------------------------------------------------------------
  subroutine veg_cs_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_cs_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine veg_ns_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine veg_ns_init
    
  !------------------------------------------------------------------------
  subroutine veg_ns_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ns_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine veg_ps_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    integer, intent(in) :: begp,endp
    !------------------------------------------------------------------------
    
  end subroutine veg_ps_init
    
  !------------------------------------------------------------------------
  subroutine veg_ps_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ps_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy flux data structure
  !------------------------------------------------------------------------
  subroutine veg_ef_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer  :: l,c,p
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of veg_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_sh_grnd        (begp:endp))   ; this%eflx_sh_grnd       (:)   = nan
    allocate(this%eflx_sh_veg         (begp:endp))   ; this%eflx_sh_veg        (:)   = nan
    allocate(this%eflx_sh_snow        (begp:endp))   ; this%eflx_sh_snow       (:)   = nan
    allocate(this%eflx_sh_soil        (begp:endp))   ; this%eflx_sh_soil       (:)   = nan
    allocate(this%eflx_sh_h2osfc      (begp:endp))   ; this%eflx_sh_h2osfc     (:)   = nan
    allocate(this%eflx_sh_tot         (begp:endp))   ; this%eflx_sh_tot        (:)   = nan
    allocate(this%eflx_sh_tot_u       (begp:endp))   ; this%eflx_sh_tot_u      (:)   = nan
    allocate(this%eflx_sh_tot_r       (begp:endp))   ; this%eflx_sh_tot_r      (:)   = nan
    allocate(this%eflx_lh_tot         (begp:endp))   ; this%eflx_lh_tot        (:)   = nan
    allocate(this%eflx_lh_tot_u       (begp:endp))   ; this%eflx_lh_tot_u      (:)   = nan
    allocate(this%eflx_lh_tot_r       (begp:endp))   ; this%eflx_lh_tot_r      (:)   = nan
    allocate(this%eflx_lh_vegt        (begp:endp))   ; this%eflx_lh_vegt       (:)   = nan
    allocate(this%eflx_lh_vege        (begp:endp))   ; this%eflx_lh_vege       (:)   = nan
    allocate(this%eflx_lh_grnd        (begp:endp))   ; this%eflx_lh_grnd       (:)   = nan
    allocate(this%eflx_soil_grnd      (begp:endp))   ; this%eflx_soil_grnd     (:)   = nan
    allocate(this%eflx_soil_grnd_u    (begp:endp))   ; this%eflx_soil_grnd_u   (:)   = nan
    allocate(this%eflx_soil_grnd_r    (begp:endp))   ; this%eflx_soil_grnd_r   (:)   = nan
    allocate(this%eflx_lwrad_net      (begp:endp))   ; this%eflx_lwrad_net     (:)   = nan
    allocate(this%eflx_lwrad_net_r    (begp:endp))   ; this%eflx_lwrad_net_r   (:)   = nan
    allocate(this%eflx_lwrad_net_u    (begp:endp))   ; this%eflx_lwrad_net_u   (:)   = nan
    allocate(this%eflx_lwrad_out      (begp:endp))   ; this%eflx_lwrad_out     (:)   = nan
    allocate(this%eflx_lwrad_out_r    (begp:endp))   ; this%eflx_lwrad_out_r   (:)   = nan
    allocate(this%eflx_lwrad_out_u    (begp:endp))   ; this%eflx_lwrad_out_u   (:)   = nan
    allocate(this%eflx_gnet           (begp:endp))   ; this%eflx_gnet          (:)   = nan
    allocate(this%eflx_grnd_lake      (begp:endp))   ; this%eflx_grnd_lake     (:)   = nan
    allocate(this%eflx_anthro         (begp:endp))   ; this%eflx_anthro        (:)   = nan
    allocate(this%eflx_traffic        (begp:endp))   ; this%eflx_traffic       (:)   = nan
    allocate(this%eflx_wasteheat      (begp:endp))   ; this%eflx_wasteheat     (:)   = nan
    allocate(this%eflx_heat_from_ac   (begp:endp))   ; this%eflx_heat_from_ac  (:)   = nan
    allocate(this%dlrad               (begp:endp))   ; this%dlrad              (:)   = nan
    allocate(this%ulrad               (begp:endp))   ; this%ulrad              (:)   = nan
    allocate(this%taux                (begp:endp))   ; this%taux               (:)   = nan
    allocate(this%tauy                (begp:endp))   ; this%tauy               (:)   = nan
    allocate(this%dgnetdT             (begp:endp))   ; this%dgnetdT            (:)   = nan
    allocate(this%netrad              (begp:endp))   ; this%netrad             (:)   = nan
    allocate(this%cgrnd               (begp:endp))   ; this%cgrnd              (:)   = nan
    allocate(this%cgrndl              (begp:endp))   ; this%cgrndl             (:)   = nan
    allocate(this%cgrnds              (begp:endp))   ; this%cgrnds             (:)   = nan
    allocate(this%errsoi              (begp:endp))   ; this%errsoi             (:)   = nan
    allocate(this%errseb              (begp:endp))   ; this%errseb             (:)   = nan
    allocate(this%errsol              (begp:endp))   ; this%errsol             (:)   = nan
    allocate(this%errlon              (begp:endp))   ; this%errlon             (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_ef
    !-----------------------------------------------------------------------

    this%eflx_lwrad_net(begp:endp) = spval
    call hist_addfld1d (fname='FIRA', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net, c2l_scale_type='urbanf')

    this%eflx_lwrad_net_r(begp:endp) = spval 
    call hist_addfld1d (fname='FIRA_R', units='W/m^2',  &
         avgflag='A', long_name='Rural net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_r, set_spec=spval)

    this%eflx_lwrad_out(begp:endp) = spval 
    call hist_addfld1d (fname='FIRE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out, c2l_scale_type='urbanf')

    this%eflx_lwrad_out(begp:endp) = spval 
    call hist_addfld1d (fname='LWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling longwave radiation', &
         ptr_patch=this%eflx_lwrad_out, c2l_scale_type='urbanf', default='inactive')

    this%eflx_lwrad_out_r(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_R', units='W/m^2',  &
         avgflag='A', long_name='Rural emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_r, set_spec=spval)

    this%eflx_lh_vegt(begp:endp) = spval
    call hist_addfld1d (fname='FCTR', units='W/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_patch=this%eflx_lh_vegt, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_vege(begp:endp) = spval
    call hist_addfld1d (fname='FCEV', units='W/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_patch=this%eflx_lh_vege, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_lh_grnd(begp:endp) = spval
    call hist_addfld1d (fname='FGEV', units='W/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_patch=this%eflx_lh_grnd, c2l_scale_type='urbanf') 

    this%eflx_sh_tot(begp:endp) = spval
    call hist_addfld1d (fname='FSH_NODYNLNDUSE', units='W/m^2',  &
         avgflag='A', long_name='sensible heat not including correction for land use change', &
         ptr_patch=this%eflx_sh_tot, c2l_scale_type='urbanf')

    this%eflx_sh_tot_r(begp:endp) = spval
    call hist_addfld1d (fname='FSH_R', units='W/m^2',  &
         avgflag='A', long_name='Rural sensible heat', &
         ptr_patch=this%eflx_sh_tot_r, set_spec=spval)

    this%eflx_sh_tot(begp:endp) = spval
    call hist_addfld1d (fname='Qh', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_patch=this%eflx_sh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot(begp:endp) = spval
    call hist_addfld1d (fname='Qle', units='W/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_patch=this%eflx_lh_tot, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
         avgflag='A', long_name='total latent heat flux [+ to atm]', &
         ptr_patch=this%eflx_lh_tot, c2l_scale_type='urbanf')

    this%eflx_lh_tot_r(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_R', units='W/m^2',  &
         avgflag='A', long_name='Rural total evaporation', &
         ptr_patch=this%eflx_lh_tot_r, set_spec=spval)

    this%eflx_soil_grnd(begp:endp) = spval
    call hist_addfld1d (fname='Qstor', units='W/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_patch=this%eflx_soil_grnd, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_sh_veg(begp:endp) = spval
    call hist_addfld1d (fname='FSH_V', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_patch=this%eflx_sh_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    this%eflx_sh_grnd(begp:endp) = spval
    call hist_addfld1d (fname='FSH_G', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_patch=this%eflx_sh_grnd, c2l_scale_type='urbanf')

    this%eflx_soil_grnd(begp:endp) = spval
    call hist_addfld1d (fname='FGR', units='W/m^2',  &
         avgflag='A', long_name='heat flux into soil/snow including snow melt and lake / snow light transmission', &
         ptr_patch=this%eflx_soil_grnd, c2l_scale_type='urbanf')

    this%eflx_soil_grnd_r(begp:endp) = spval
    call hist_addfld1d (fname='FGR_R', units='W/m^2',  &
         avgflag='A', long_name='Rural heat flux into soil/snow including snow melt and snow light transmission', &
         ptr_patch=this%eflx_soil_grnd_r, set_spec=spval)

    this%eflx_lwrad_net_u(begp:endp) = spval
    call hist_addfld1d (fname='FIRA_U', units='W/m^2',  &
         avgflag='A', long_name='Urban net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_soil_grnd(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
         avgflag='A', long_name='soil heat flux [+ into soil]', &
         ptr_patch=this%eflx_soil_grnd, default='inactive', c2l_scale_type='urbanf')

    this%eflx_lwrad_out_u(begp:endp) = spval
    call hist_addfld1d (fname='FIRE_U', units='W/m^2',  &
         avgflag='A', long_name='Urban emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_sh_tot_u(begp:endp) = spval
    call hist_addfld1d (fname='FSH_U', units='W/m^2',  &
         avgflag='A', long_name='Urban sensible heat', &
         ptr_patch=this%eflx_sh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_lh_tot_u(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT_U', units='W/m^2',  &
         avgflag='A', long_name='Urban total evaporation', &
         ptr_patch=this%eflx_lh_tot_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_soil_grnd_u(begp:endp) = spval
    call hist_addfld1d (fname='FGR_U', units='W/m^2',  &
         avgflag='A', long_name='Urban heat flux into soil/snow including snow melt', &
         ptr_patch=this%eflx_soil_grnd_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%netrad(begp:endp) = spval
    call hist_addfld1d (fname='Rnet', units='W/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_patch=this%netrad, c2l_scale_type='urbanf', &
         default='inactive')

    if (use_cn) then
       this%dlrad(begp:endp) = spval
       call hist_addfld1d (fname='DLRAD', units='W/m^2', &
            avgflag='A', long_name='downward longwave radiation below the canopy', &
            ptr_patch=this%dlrad, default='inactive', c2l_scale_type='urbanf')

       this%ulrad(begp:endp) = spval
       call hist_addfld1d (fname='ULRAD', units='W/m^2', &
            avgflag='A', long_name='upward longwave radiation above the canopy', &
            ptr_patch=this%ulrad, default='inactive', c2l_scale_type='urbanf')

       this%cgrnd(begp:endp) = spval
       call hist_addfld1d (fname='CGRND', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil energy flux wrt to soil temp', &
            ptr_patch=this%cgrnd, default='inactive', c2l_scale_type='urbanf')

       this%cgrndl(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDL', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil latent heat flux wrt soil temp', &
            ptr_patch=this%cgrndl, default='inactive', c2l_scale_type='urbanf')

       this%cgrnds(begp:endp) = spval
       call hist_addfld1d (fname='CGRNDS', units='W/m^2/K', &
            avgflag='A', long_name='deriv. of soil sensible heat flux wrt soil temp', &
            ptr_patch=this%cgrnds, default='inactive', c2l_scale_type='urbanf')

       this%eflx_gnet(begp:endp) = spval
       call hist_addfld1d (fname='EFLX_GNET', units='W/m^2', &
            avgflag='A', long_name='net heat flux into ground', &
            ptr_patch=this%eflx_gnet, default='inactive', c2l_scale_type='urbanf')
    end if ! use_cn

    this%eflx_grnd_lake(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_GRND_LAKE', units='W/m^2', &
         avgflag='A', long_name='net heat flux into lake/snow surface, excluding light transmission', &
         ptr_patch=this%eflx_grnd_lake, set_nolake=spval)

    this%dgnetdT(begp:endp) = spval
    call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
         avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
         ptr_patch=this%dgnetdT, default='inactive', c2l_scale_type='urbanf')

    this%eflx_traffic(begp:endp) = spval
    call hist_addfld1d (fname='TRAFFICFLUX', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from urban traffic', &
         ptr_patch=this%eflx_traffic, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    this%eflx_wasteheat(begp:endp) = spval
    call hist_addfld1d (fname='WASTEHEAT', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux from heating/cooling sources of urban waste heat', &
         ptr_patch=this%eflx_wasteheat, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_heat_from_ac(begp:endp) = spval
    call hist_addfld1d (fname='HEAT_FROM_AC', units='W/m^2',  &
         avgflag='A', long_name='sensible heat flux put into canyon due to heat removed from air conditioning', &
         ptr_patch=this%eflx_heat_from_ac, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_anthro(begp:endp) = spval
    call hist_addfld1d (fname='Qanth', units='W/m^2',  &
         avgflag='A', long_name='anthropogenic heat flux', &
         ptr_patch=this%eflx_anthro, set_nourb=0._r8, c2l_scale_type='urbanf', &
         default='inactive')

    this%taux(begp:endp) = spval
    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_patch=this%taux)

    this%tauy(begp:endp) = spval
    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_patch=this%tauy)

    this%errseb(begp:endp) = spval
    call hist_addfld1d (fname='ERRSEB',  units='W/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_patch=this%errseb)

    this%errsol(begp:endp) = spval
    call hist_addfld1d (fname='ERRSOL',  units='W/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_patch=this%errsol, set_urb=spval)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_ef
    !-----------------------------------------------------------------------
    
     do p = begp, endp 
        c = veg_pp%column(p)
        l = veg_pp%landunit(p)

        if (.not. lun_pp%urbpoi(l)) then ! non-urban
           this%eflx_lwrad_net_u(p)  = spval
           this%eflx_lwrad_out_u(p)  = spval
           this%eflx_lh_tot_u(p)     = spval
           this%eflx_sh_tot_u(p)     = spval
           this%eflx_soil_grnd_u(p)  = spval
           this%eflx_wasteheat(p)    = 0._r8
           this%eflx_heat_from_ac(p) = 0._r8
           this%eflx_traffic(p)      = 0._r8
           this%eflx_anthro(p)       = 0._r8
        end if

        this%eflx_lwrad_out(p) = sb * (col_es%t_grnd(c))**4
     end do
    
  end subroutine veg_ef_init
    
  !------------------------------------------------------------------------
  subroutine veg_ef_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation energy flux information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  & 
         dim1name='pft', &
         long_name='emitted infrared (longwave) radiation', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_lwrad_out)

  end subroutine veg_ef_restart

  !------------------------------------------------------------------------
  subroutine veg_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_ef_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water flux data structure
  !------------------------------------------------------------------------
  subroutine veg_wf_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,l                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_wf
    !-----------------------------------------------------------------------
    allocate(this%qflx_prec_grnd         (begp:endp))             ; this%qflx_prec_grnd       (:)   = nan
    allocate(this%qflx_rain_grnd         (begp:endp))             ; this%qflx_rain_grnd       (:)   = nan
    allocate(this%qflx_snow_grnd         (begp:endp))             ; this%qflx_snow_grnd       (:)   = nan
    allocate(this%qflx_sub_snow          (begp:endp))             ; this%qflx_sub_snow        (:)   = 0._r8
    allocate(this%qflx_evap_soi          (begp:endp))             ; this%qflx_evap_soi        (:)   = nan
    allocate(this%qflx_evap_veg          (begp:endp))             ; this%qflx_evap_veg        (:)   = nan
    allocate(this%qflx_evap_can          (begp:endp))             ; this%qflx_evap_can        (:)   = nan
    allocate(this%qflx_evap_tot          (begp:endp))             ; this%qflx_evap_tot        (:)   = nan
    allocate(this%qflx_evap_grnd         (begp:endp))             ; this%qflx_evap_grnd       (:)   = nan
    allocate(this%qflx_snwcp_liq         (begp:endp))             ; this%qflx_snwcp_liq       (:)   = nan
    allocate(this%qflx_snwcp_ice         (begp:endp))             ; this%qflx_snwcp_ice       (:)   = nan
    allocate(this%qflx_tran_veg          (begp:endp))             ; this%qflx_tran_veg        (:)   = nan
    allocate(this%qflx_dew_snow          (begp:endp))             ; this%qflx_dew_snow        (:)   = nan
    allocate(this%qflx_dew_grnd          (begp:endp))             ; this%qflx_dew_grnd        (:)   = nan
    allocate(this%qflx_prec_intr         (begp:endp))             ; this%qflx_prec_intr       (:)   = nan
    allocate(this%qflx_ev_snow           (begp:endp))             ; this%qflx_ev_snow         (:)   = nan
    allocate(this%qflx_ev_soil           (begp:endp))             ; this%qflx_ev_soil         (:)   = nan
    allocate(this%qflx_ev_h2osfc         (begp:endp))             ; this%qflx_ev_h2osfc       (:)   = nan
    allocate(this%qflx_rootsoi_frac      (begp:endp,1:nlevgrnd))  ; this%qflx_rootsoi_frac    (:,:) = nan
    allocate(this%qflx_irrig             (begp:endp))             ; this%qflx_irrig           (:)   = nan
    allocate(this%irrig_rate             (begp:endp))             ; this%irrig_rate           (:)   = nan
    allocate(this%n_irrig_steps_left     (begp:endp))             ; this%n_irrig_steps_left   (:)   = 0
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_wf
    !-----------------------------------------------------------------------
    this%qflx_irrig(begp:endp) = spval
    call hist_addfld1d (fname='QIRRIG', units='mm/s', &
         avgflag='A', long_name='water added through irrigation', &
         ptr_patch=this%qflx_irrig)

    this%qflx_prec_intr(begp:endp) = spval
    call hist_addfld1d (fname='QINTR', units='mm/s',  &
         avgflag='A', long_name='interception', &
         ptr_patch=this%qflx_prec_intr, set_lake=0._r8)

    this%qflx_prec_grnd(begp:endp) = spval
    call hist_addfld1d (fname='QDRIP', units='mm/s',  &
         avgflag='A', long_name='throughfall', &
         ptr_patch=this%qflx_prec_grnd, c2l_scale_type='urbanf')

    this%qflx_evap_soi(begp:endp) = spval
    call hist_addfld1d (fname='QSOIL', units='mm/s',  &
         avgflag='A', long_name= 'Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)', &
         ptr_patch=this%qflx_evap_soi, c2l_scale_type='urbanf')

    this%qflx_evap_can(begp:endp) = spval
    call hist_addfld1d (fname='QVEGE', units='mm/s',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_patch=this%qflx_evap_can, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_tran_veg(begp:endp) = spval
    call hist_addfld1d (fname='QVEGT', units='mm/s',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_patch=this%qflx_tran_veg, set_lake=0._r8, c2l_scale_type='urbanf')

    this%qflx_snwcp_liq(begp:endp) = spval
    call hist_addfld1d (fname='QSNWCPLIQ', units='mm H2O/s', &
         avgflag='A', long_name='excess rainfall due to snow capping', &
         ptr_patch=this%qflx_snwcp_liq, c2l_scale_type='urbanf', default='inactive')

    if (use_cn) then
       this%qflx_rain_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_RAIN_GRND', units='mm H2O/s', &
            avgflag='A', long_name='rain on ground after interception', &
            ptr_patch=this%qflx_rain_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_snow_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_SNOW_GRND', units='mm H2O/s', &
            avgflag='A', long_name='snow on ground after interception', &
            ptr_patch=this%qflx_snow_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_evap_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_EVAP_GRND', units='mm H2O/s', &
            avgflag='A', long_name='ground surface evaporation', &
            ptr_patch=this%qflx_evap_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_evap_veg(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_EVAP_VEG', units='mm H2O/s', &
            avgflag='A', long_name='vegetation evaporation', &
            ptr_patch=this%qflx_evap_veg, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_evap_tot(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_EVAP_TOT', units='mm H2O/s', &
            avgflag='A', long_name='qflx_evap_soi + qflx_evap_can + qflx_tran_veg', &
            ptr_patch=this%qflx_evap_tot, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_dew_grnd(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_DEW_GRND', units='mm H2O/s', &
            avgflag='A', long_name='ground surface dew formation', &
            ptr_patch=this%qflx_dew_grnd, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_sub_snow(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_SUB_SNOW', units='mm H2O/s', &
            avgflag='A', long_name='sublimation rate from snow pack', &
            ptr_patch=this%qflx_sub_snow, default='inactive', c2l_scale_type='urbanf')
    end if

    if (use_cn) then
       this%qflx_dew_snow(begp:endp) = spval
       call hist_addfld1d (fname='QFLX_DEW_SNOW', units='mm H2O/s', &
            avgflag='A', long_name='surface dew added to snow pacK', &
            ptr_patch=this%qflx_dew_snow, default='inactive', c2l_scale_type='urbanf')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_wf
    !-----------------------------------------------------------------------
    this%qflx_evap_grnd(begp:endp) = 0.0_r8
    this%qflx_dew_grnd (begp:endp) = 0.0_r8
    this%qflx_dew_snow (begp:endp) = 0.0_r8

    do p = begp, endp
       l = veg_pp%landunit(p)
       
       if (lun_pp%itype(l)==istsoil) then
          this%n_irrig_steps_left(p) = 0
          this%irrig_rate(p)         = 0.0_r8
       end if
    end do
    
  end subroutine veg_wf_init
    
  !------------------------------------------------------------------------
  subroutine veg_wf_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write vegetation water flux information to/from restart file.
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux)     :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: dimlen       ! dimension length
    integer :: nump_global  ! total number of pfts, globally
    integer :: err_code     ! error code
    logical :: do_io
    !-----------------------------------------------------------------------

    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'n_irrig_steps_left', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='n_irrig_steps_left', xtype=ncd_int,  &
            dim1name='pft', &
            long_name='number of irrigation time steps left', units='#', &
            interpinic_flag='interp', readvar=readvar, data=this%n_irrig_steps_left)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%n_irrig_steps_left = 0
    end if

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'irrig_rate', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='irrig_rate', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='irrigation rate', units='mm/s', &
            interpinic_flag='interp', readvar=readvar, data=this%irrig_rate)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%irrig_rate = 0.0_r8
    end if

  end subroutine veg_wf_restart

  !------------------------------------------------------------------------
  subroutine veg_wf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_wf_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon flux data structure
  !------------------------------------------------------------------------
  subroutine veg_cf_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,j                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_cf
    !-----------------------------------------------------------------------
    allocate(this%xxx              (begp:endp))                   ; this%xxx              (:)   = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_cf
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_cf
    !-----------------------------------------------------------------------
    
  end subroutine veg_cf_init
    
  !------------------------------------------------------------------------
  subroutine veg_cf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_cf_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine veg_nf_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,j                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_nf
    !-----------------------------------------------------------------------
    allocate(this%xxx              (begp:endp))                   ; this%xxx              (:)   = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_nf
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_nf
    !-----------------------------------------------------------------------
    
  end subroutine veg_nf_init
    
  !------------------------------------------------------------------------
  subroutine veg_nf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine veg_nf_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine veg_pf_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    integer, intent(in) :: begp,endp
    !
    ! !LOCAL VARIABLES:
    integer           :: p,c,l,j                        ! indices
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of veg_pf
    !-----------------------------------------------------------------------
    allocate(this%xxx              (begp:endp))                   ; this%xxx              (:)   = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of veg_pf
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of veg_pf
    !------------------------------------------------------------------------
    
  end subroutine veg_pf_init
    
  !------------------------------------------------------------------------
  subroutine veg_pf_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    !-----------------------------------------------------------------------
   
  end subroutine veg_pf_clean
  
    !------------------------------------------------------------------------
    
end module vegetationDataType
