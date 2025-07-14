module TopounitDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Topounit data type allocation and initialization
  ! --------------------------------------------------------
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use elm_varcon     , only : spval, ispval
  use elm_varctl     , only : iulog, use_cn, use_fates, use_lch4
  use elm_varpar     , only : numrad
  use histFileMod    , only : hist_addfld1d, hist_addfld2d
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_topounit
  use TopounitType   , only : top_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  !-----------------------------------------------------------------------
  ! Define the data structure where land model receives atmospheric state information.
  type, public :: topounit_atmospheric_state
    real(r8), pointer :: tbot       (:) => null() ! temperature of air at atmospheric forcing height (K)
    real(r8), pointer :: thbot      (:) => null() ! potential temperature of air at atmospheric forcing height (K)
    real(r8), pointer :: pbot       (:) => null() ! air pressure at atmospheric forcing height (Pa)
    real(r8), pointer :: rhobot     (:) => null() ! air density at atmospheric forcing height (kg/m**3)
    real(r8), pointer :: qbot       (:) => null() ! specific humidity at atmospheric forcing height (kg H2O/kg moist air)
    real(r8), pointer :: rhbot      (:) => null() ! relative humidity at atmospheric forcing height (%)
    real(r8), pointer :: ubot       (:) => null() ! wind speed in U (east) direction at atmospheric forcing height (m/s)
    real(r8), pointer :: vbot       (:) => null() ! wind speed in V (north) direction at atmospheric forcing height (m/s)
    real(r8), pointer :: wsresp     (:) => null() ! first-order response of wind to surface stress (m/s/Pa)
    real(r8), pointer :: tau_est    (:) => null() ! estimate of stress in equilibrium with winds (Pa)
    real(r8), pointer :: ugust      (:) => null() ! extra gustiness from atmosphere (m/s)
    real(r8), pointer :: windbot    (:) => null() ! horizontal component of wind at atmospheric forcing height (m/s)
    real(r8), pointer :: zbot       (:) => null() ! atmospheric forcing height (m)
    real(r8), pointer :: po2bot     (:) => null() ! partial pressure of O2 at atmospheric forcing height (Pa)
    real(r8), pointer :: pco2bot    (:) => null() ! partial pressure of CO2 at atmospheric forcing height (Pa)
    real(r8), pointer :: pc13o2bot  (:) => null() ! partial pressure of C13O2 at atmospheric forcing height (Pa)
    real(r8), pointer :: pch4bot    (:) => null() ! partial pressure of CH4 at atmospheric forcing height (Pa)
    ! Accumulated fields
    real(r8), pointer :: rh24h      (:) => null() ! 24-hour running mean of relative humidity at atmospheric forcing height (%)
    real(r8), pointer :: wind24h    (:) => null() ! 24-hour running mean of horizontal wind at atmospheric forcing height (m/s)
  contains
    procedure, public :: Init  => init_top_as
    procedure, public :: Clean => clean_top_as
    procedure, public :: InitAccBuffer => init_acc_buffer_top_as
    procedure, public :: InitAccVars   => init_acc_vars_top_as
     procedure, public :: UpdateAccVars => update_acc_vars_top_as
  end type topounit_atmospheric_state

  !-----------------------------------------------------------------------
  ! Define the data structure that where land model receives atmospheric flux information.
  type, public :: topounit_atmospheric_flux
    real(r8), pointer :: rain      (:)   => null() ! rain rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: snow      (:)   => null() ! snow rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: solad     (:,:) => null() ! direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll) (W/m**2)
    real(r8), pointer :: solai     (:,:) => null() ! diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld) (W/m**2)
    real(r8), pointer :: solar     (:)   => null() ! incident solar radiation (W/m**2)
    real(r8), pointer :: lwrad     (:)   => null() ! atm downwrd IR longwave radiation (W/m**2)
    ! Accumulated fields
    real(r8), pointer :: prec24h   (:)   => null() ! 24-hour mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: prec10d   (:)   => null() ! 10-day mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: prec60d   (:)   => null() ! 60-day mean precip rate (kg H2O/m**2/s, equivalent to mm liquid H2O/s)
    real(r8), pointer :: fsd24h    (:)   => null() ! 24hr average of direct beam radiation (W/m**2)
    real(r8), pointer :: fsd240h   (:)   => null() ! 240hr average of direct beam radiation (W/m**2)
    real(r8), pointer :: fsi24h    (:)   => null() ! 24hr average of diffuse beam radiation (W/m**2)
    real(r8), pointer :: fsi240h   (:)   => null() ! 240hr average of diffuse beam radiation (W/m**2)

  contains
    procedure, public :: Init  => init_top_af
    procedure, public :: Clean => clean_top_af
    procedure, public :: InitAccBuffer => init_acc_buffer_top_af
    procedure, public :: InitAccVars   => init_acc_vars_top_af
     procedure, public :: UpdateAccVars => update_acc_vars_top_af
  end type topounit_atmospheric_flux

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information for land at the level of topographic unit.
  type, public :: topounit_energy_state
    real(r8), pointer :: t_rad      (:) => null() ! mean radiative temperature of land surface (K)
    real(r8), pointer :: eflx_lwrad_out_topo      (:) => null() ! Topounit level longwave radiation flux to be used in downscaling
    real(r8), pointer :: t_grnd     (:) => null()
	
	! temperature variables
    real(r8), pointer :: heat1_tgu                (:)   ! initial tgu total heat content
    real(r8), pointer :: heat2_tgu                (:)   ! post land cover change total heat content
    real(r8), pointer :: liquid_water_temp1_tgu   (:)   ! initial weighted average liquid water temperature (K)
    real(r8), pointer :: liquid_water_temp2_tgu   (:)   ! post land cover change weighted average liquid water temperature (K)

  contains
    procedure, public :: Init  => init_top_es
    procedure, public :: Restart => restart_top_es
    procedure, public :: Clean => clean_top_es
  end type topounit_energy_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: topounit_energy_flux
    real(r8), pointer :: eflx_dynbal           (:)   ! dynamic land cover change conversion energy flux (W/m**2)

  contains
    procedure, public :: Init    => top_ef_init
    procedure, public :: Restart => top_ef_restart
    procedure, public :: Clean   => top_ef_clean
  end type topounit_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: topounit_water_state
    ! Note: All units given as kg in this data type imply kg of H2O
    real(r8), pointer :: liq1               (:)   ! initial topounit total h2o liq content (kg/m2)
    real(r8), pointer :: liq2               (:)   ! post land cover change total liq content (kg/m2)
    real(r8), pointer :: ice1               (:)   ! initial topounit total h2o ice content (kg/m2)
    real(r8), pointer :: ice2               (:)   ! post land cover change total ice content (kg/m2)
    real(r8), pointer :: tws                (:)   ! total water storage (kg/m2)
    real(r8), pointer :: tws_month_beg      (:)   ! top total water storage at the beginning of a month
    real(r8), pointer :: tws_month_end      (:)   ! top total water storage at the end of a month
    real(r8), pointer :: begwb              (:)   ! water mass begining of the time step
    real(r8), pointer :: endwb              (:)   ! water mass end of the time step
    real(r8), pointer :: errh2o             (:)   ! water conservation error (mm H2O)
    real(r8), pointer :: beg_h2ocan         (:)   ! topounit-level canopy water at begining of the time step (mm H2O)
    real(r8), pointer :: beg_h2osno         (:)   ! topounit-level snow water at begining of the time step (mm H2O)
    real(r8), pointer :: beg_h2osfc         (:)   ! topounit-level surface water at begining of the time step (mm H2O)
    real(r8), pointer :: beg_h2osoi_liq     (:)   ! topounit-level liquid water at begining of the time step (kg/m2)
    real(r8), pointer :: beg_h2osoi_ice     (:)   ! topounit-level ice lens at begining of the time step (kg/m2)
    real(r8), pointer :: end_h2ocan         (:)   ! topounit-level canopy water at end of the time step (mm H2O)
    real(r8), pointer :: end_h2osno         (:)   ! topounit-level snow water at end of the time step (mm H2O)
    real(r8), pointer :: end_h2osfc         (:)   ! topounit-level surface water at end of the time step (mm H2O)
    real(r8), pointer :: end_h2osoi_liq     (:)   ! topounit-level liquid water at end of the time step (kg/m2)
    real(r8), pointer :: end_h2osoi_ice     (:)   ! topounit-level ice lens at end of the time step (kg/m2)

  contains
    procedure, public :: Init    => top_ws_init
    procedure, public :: Restart => top_ws_restart
    procedure, public :: Clean   => top_ws_clean
  end type topounit_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: topounit_water_flux
    ! Dynamic land cover change
    real(r8), pointer :: qflx_liq_dynbal      (:)   ! liq dynamic land cover change conversion runoff flux
    real(r8), pointer :: qflx_ice_dynbal      (:)   ! ice dynamic land cover change conversion runoff flux

    ! Objects that help convert once-per-year dynamic land cover changes into fluxes
    ! that are dribbled throughout the year
    type(annual_flux_dribbler_type) :: qflx_liq_dynbal_dribbler
    type(annual_flux_dribbler_type) :: qflx_ice_dynbal_dribbler

  contains
    procedure, public :: Init    => top_wf_init
    procedure, public :: Restart => top_wf_restart
    procedure, public :: Clean   => top_wf_clean
  end type topounit_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the topounit level.
  !-----------------------------------------------------------------------
  type, public :: topounit_carbon_state
    real(r8), pointer :: seedc                 (:) => null() ! (gC/m2) pool for seeding new PFTs via dynamic landcover
    real(r8), pointer :: tcs_month_beg         (:) => null()  ! top total carbon storage at the beginning of a month
    real(r8), pointer :: tcs_month_end         (:) => null()  ! top total carbon storage at the end of a month
    real(r8), pointer :: begcb                 (:) => null() ! carbon mass, beginning of time step (gC/m**2)
    real(r8), pointer :: endcb                 (:) => null() ! carbon mass, end of time step (gC/m**2)
    real(r8), pointer :: errcb                 (:) => null() ! carbon balance error for the timestep (gC/m**2)
    real(r8), pointer :: beg_totc              (:) => null() ! (gC/m2) total carbon, including veg and cpool at begining of the time step
    real(r8), pointer :: beg_totpftc           (:) => null() ! (gC/m2) total carbon, including cpool at begining of the time step
    real(r8), pointer :: beg_cwdc              (:) => null() ! (gC/m2) Diagnostic: coarse woody debris C at begining of the time step
    real(r8), pointer :: beg_totsomc           (:) => null() ! (gC/m2) total soil organic matter carbon at begining of the time step
    real(r8), pointer :: beg_totlitc           (:) => null() ! (gC/m2) total litter carbon at begining of the time step
    real(r8), pointer :: beg_totprodc          (:) => null() ! (gC/m2) total wood product C at begining of the time step
    real(r8), pointer :: beg_ctrunc            (:) => null() ! (gC/m2) column-level sink for C truncation at begining of the time step
    real(r8), pointer :: beg_cropseedc_deficit (:) => null() ! (gC/m2) crop seed C deficit at begining of the time step
    real(r8), pointer :: end_totc              (:) => null() ! (gC/m2) total carbon, including veg and cpool at end of the time step
    real(r8), pointer :: end_totpftc           (:) => null() ! (gC/m2) total carbon, including cpool at end of the time step
    real(r8), pointer :: end_cwdc              (:) => null() ! (gC/m2) Diagnostic: coarse woody debris C at end of the time step
    real(r8), pointer :: end_totsomc           (:) => null() ! (gC/m2) total soil organic matter carbon at end of the time step
    real(r8), pointer :: end_totlitc           (:) => null() ! (gC/m2) total litter carbon at end of the time step
    real(r8), pointer :: end_totprodc          (:) => null() ! (gC/m2) total wood product C at end of the time step
    real(r8), pointer :: end_ctrunc            (:) => null() ! (gC/m2) column-level sink for C truncation at end of the time step
    real(r8), pointer :: end_cropseedc_deficit (:) => null() ! (gC/m2) crop seed C deficit at end of the time step
  contains
    procedure, public :: Init    => top_cs_init
    procedure, public :: Restart => top_cs_restart
    procedure, public :: Clean   => top_cs_clean
  end type topounit_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the topounit level.
  !-----------------------------------------------------------------------
  type, public :: topounit_carbon_flux
    real(r8), pointer :: dwt_seedc_to_leaf         (:) => null()  ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the topounit-level
    real(r8), pointer :: dwt_seedc_to_deadstem     (:) => null()  ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the topounit-level
    real(r8), pointer :: dwt_conv_cflux            (:) => null()  ! (gC/m2/s) dwt_conv_cflux_patch summed to the topounit-level
    real(r8), pointer :: dwt_conv_cflux_dribbled   (:) => null()  ! (gC/m2/s) dwt_conv_cflux dribbled evenly throughout the year
    real(r8), pointer :: dwt_prod10c_gain          (:) => null()  ! (gC/m2/s) dynamic landcover addition to 10-year wood product pool
    real(r8), pointer :: dwt_prod100c_gain         (:) => null()  ! (gC/m2/s) dynamic landcover addition to 100-year wood product pool
    real(r8), pointer :: hrv_deadstemc_to_prod10c  (:) => null()  ! (gC/m2/s) dead stem harvest to 10-year wood product pool
    real(r8), pointer :: hrv_deadstemc_to_prod100c (:) => null()  ! (gC/m2/s) dead stem harvest to 100-year wood product pool
    real(r8), pointer :: cinputs                   (:) => null()  ! (gC/m2/s) topounit-level C inputs
    real(r8), pointer :: coutputs                  (:) => null()  ! (gC/m2/s) topounit-level C outputs
    real(r8), pointer :: gpp                       (:) => null()  ! (gC/m2/s) topounit-level gross primary production
    real(r8), pointer :: er                        (:) => null()  ! (gC/m2/s) topounit-level total ecosystem respiration
    real(r8), pointer :: fire_closs                (:) => null()  ! (gC/m2/s) topounit-level total fire C loss
    real(r8), pointer :: prod1_loss                (:) => null()  ! (gC/m2/s) topounit-level crop leafc harvested
    real(r8), pointer :: prod10_loss               (:) => null()  ! (gC/m2/s) topounit-level 10-year wood C harvested
    real(r8), pointer :: prod100_loss              (:) => null()  ! (gC/m2/s) topounit-level 100-year wood C harvested
    real(r8), pointer :: hrv_xsmrpool_to_atm       (:) => null()  ! (gC/m2/s) topounit-level excess MR pool harvest mortality
    real(r8), pointer :: som_c_leached             (:) => null()  ! (gC/m2/s) topounit-level total SOM C loss from vertical transport
    real(r8), pointer :: somc_yield                (:) => null()  ! (gC/m2/s) topounit-level total SOM C loss by erosion

  contains
    procedure, public :: Init    => top_cf_init
    procedure, public :: ZeroDWT => top_cf_zerodwt
    procedure, public :: Clean   => top_cf_clean
  end type topounit_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the topounit level.
  !-----------------------------------------------------------------------
  type, public :: topounit_nitrogen_state
    real(r8), pointer :: seedn          (:) => null()   ! (gNm2) nitrogen pool for seeding new PFTs via dynamic landcover
    real(r8), pointer :: begnb          (:) => null()   ! (gNm2) nitrogen mass, beginning of time step
    real(r8), pointer :: endnb          (:) => null()   ! (gNm2) nitrogen mass, end of time step 
    real(r8), pointer :: errnb          (:) => null()   ! (gNm2) nitrogen balance error for the timestep

  contains
    procedure, public :: Init    => top_ns_init
    procedure, public :: Clean   => top_ns_clean
  end type topounit_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the topounit level.
  !-----------------------------------------------------------------------
  type, public :: topounit_nitrogen_flux
    ! dynamic landcover fluxes
    real(r8), pointer :: dwt_seedn_to_leaf      (:) => null()  ! (gN/m2/s) dwt_seedn_to_leaf_patch summed to the topounit-level
    real(r8), pointer :: dwt_seedn_to_deadstem  (:) => null()  ! (gN/m2/s) dwt_seedn_to_deadstem_patch summed to the topounit-level
    real(r8), pointer :: dwt_conv_nflux         (:) => null()  ! (gN/m2/s) dwt_conv_nflux_patch summed to the topounit-level
    real(r8), pointer :: dwt_seedn_to_npool     (:) => null()  ! (gN/m2/s) seed source to PFT level
    real(r8), pointer :: dwt_prod10n_gain       (:) => null()  ! (gN/m2/s) addition to 10-yr wood product pool
    real(r8), pointer :: dwt_prod100n_gain      (:) => null()  ! (gN/m2/s) addition to 100-yr wood product pool
    real(r8), pointer :: ninputs                (:) => null()  ! (gN/m2/s) topounit-level N inputs
    real(r8), pointer :: noutputs               (:) => null()  ! (gN/m2/s) topounit-level N outputs
  contains
    procedure, public :: Init    => top_nf_init
    procedure, public :: ZeroDWT => top_nf_zerodwt
    procedure, public :: Clean   => top_nf_clean
  end type topounit_nitrogen_flux
 
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the topounit level.
  !-----------------------------------------------------------------------
  type, public :: topounit_phosphorus_state
    real(r8), pointer :: seedp                    (:)     ! (gP/m2) topounit-level pool for seeding new PFTs via dynamic landcover
    real(r8), pointer :: begpb                    (:)     ! phosphorus mass, beginning of time step (gP/m**2)
    real(r8), pointer :: endpb                    (:)     ! phosphorus mass, end of time step (gP/m**2)
    real(r8), pointer :: errpb                    (:)     ! phosphorus balance error for the timestep (gP/m**2)
  contains
    procedure, public :: Init    => top_ps_init
    procedure, public :: Clean   => top_ps_clean
  end type topounit_phosphorus_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the topounit level.
  !-----------------------------------------------------------------------
  type, public :: topounit_phosphorus_flux
    real(r8), pointer :: dwt_seedp_to_leaf        (:)  ! (gP/m2/s) dwt_seedn_to_leaf_patch summed to the topounit-level
    real(r8), pointer :: dwt_seedp_to_deadstem    (:)  ! (gP/m2/s) dwt_seedn_to_deadstem_patch summed to the topounit-level
    real(r8), pointer :: dwt_conv_pflux           (:)  ! (gP/m2/s) dwt_conv_nflux_patch summed to the topounit-level
    real(r8), pointer :: dwt_seedp_to_ppool       (:)  ! (gP/m2/s) seed source to PFT-level
    real(r8), pointer :: dwt_prod10p_gain         (:)  ! (gP/m2/s) addition to 10-yr wood product pool
    real(r8), pointer :: dwt_prod100p_gain        (:)  ! (gP/m2/s) addition to 100-yr wood product pool
    real(r8), pointer :: pinputs                  (:)  ! (gP/m2/s) topounit-level P inputs
    real(r8), pointer :: poutputs                 (:)  ! (gP/m2/s) topounit-level P outputs
  contains
    procedure, public :: Init    => top_pf_init
    procedure, public :: ZeroDWT => top_pf_zerodwt
    procedure, public :: Clean   => top_pf_clean
  end type topounit_phosphorus_flux
  
  !-----------------------------------------------------------------------
  ! declare the public instances of topounit data types
  type(topounit_atmospheric_state),    public, target :: top_as
  type(topounit_atmospheric_flux),     public, target :: top_af
  type(topounit_energy_state),         public, target :: top_es
  type(topounit_energy_flux)           , public, target :: top_ef     ! topounit energy flux
  type(topounit_water_state)           , public, target :: top_ws     ! topounit water state
  type(topounit_water_flux)            , public, target :: top_wf     ! topounit water flux
  type(topounit_carbon_state)          , public, target :: top_cs     ! topounit carbon state
  type(topounit_carbon_state)          , public, target :: c13_top_cs ! topounit carbon state (C13)
  type(topounit_carbon_state)          , public, target :: c14_top_cs ! topounit carbon state (C14)
  type(topounit_carbon_flux)           , public, target :: top_cf     ! topounit carbon flux
  type(topounit_carbon_flux)           , public, target :: c13_top_cf ! topounit carbon flux (C13)
  type(topounit_carbon_flux)           , public, target :: c14_top_cf ! topounit carbon flux (C14)
  type(topounit_nitrogen_state)        , public, target :: top_ns     ! topounit nitrogen state
  type(topounit_nitrogen_flux)         , public, target :: top_nf     ! topounit nitrogen flux
  type(topounit_phosphorus_state)      , public, target :: top_ps     ! topounit phosphorus state
  type(topounit_phosphorus_flux)       , public, target :: top_pf     ! topounit phosphorus flux
  !------------------------------------------------------------------------

  !$acc declare create(top_as)
  !$acc declare create(top_af)
  !$acc declare create(top_es)
  !$acc declare create(top_ef)
  !$acc declare create(top_ws)
  !$acc declare create(top_wf)
  !$acc declare create(top_cs)
  !$acc declare create(c13_top_cs)
  !$acc declare create(c14_top_cs)
  !$acc declare create(top_cf)
  !$acc declare create(c13_top_cf)
  !$acc declare create(c14_top_cf)
  !$acc declare create(top_ns)
  !$acc declare create(top_nf)
  !$acc declare create(top_ps)
  !$acc declare create(top_pf)

  contains

  !-----------------------------------------------------------------------
  subroutine init_top_as(this, begt, endt)
    class(topounit_atmospheric_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    ! Allocate for atmospheric state forcing variables, initialize to special value
    allocate(this%tbot     (begt:endt)) ; this%tbot      (:) = spval
    allocate(this%thbot    (begt:endt)) ; this%thbot     (:) = spval
    allocate(this%pbot     (begt:endt)) ; this%pbot      (:) = spval
    allocate(this%rhobot   (begt:endt)) ; this%rhobot    (:) = spval
    allocate(this%qbot     (begt:endt)) ; this%qbot      (:) = spval
    allocate(this%rhbot    (begt:endt)) ; this%rhbot     (:) = spval
    allocate(this%ubot     (begt:endt)) ; this%ubot      (:) = spval
    allocate(this%vbot     (begt:endt)) ; this%vbot      (:) = spval
    allocate(this%wsresp   (begt:endt)) ; this%wsresp    (:) = spval
    allocate(this%tau_est  (begt:endt)) ; this%tau_est   (:) = spval
    allocate(this%ugust    (begt:endt)) ; this%ugust     (:) = spval
    allocate(this%windbot  (begt:endt)) ; this%windbot   (:) = spval
    allocate(this%zbot     (begt:endt)) ; this%zbot      (:) = spval
    allocate(this%po2bot   (begt:endt)) ; this%po2bot    (:) = spval
    allocate(this%pco2bot  (begt:endt)) ; this%pco2bot   (:) = spval
    allocate(this%pc13o2bot(begt:endt)) ; this%pc13o2bot (:) = spval
    allocate(this%pch4bot  (begt:endt)) ; this%pch4bot   (:) = spval
    if (use_fates) then
      allocate(this%rh24h  (begt:endt)) ; this%rh24h     (:) = spval
      allocate(this%wind24h(begt:endt)) ; this%wind24h   (:) = spval
    end if

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_as
    !-----------------------------------------------------------------------
    this%tbot(begt:endt) = spval
    call hist_addfld1d (fname='TBOT', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_topo=this%tbot,t2g_scale_type='topounit')

    this%thbot(begt:endt) = spval
    call hist_addfld1d (fname='THBOT', units='K',  &
         avgflag='A', long_name='atmospheric air potential temperature', &
         ptr_topo=this%thbot,t2g_scale_type='topounit')

    this%pbot(begt:endt) = spval
    call hist_addfld1d (fname='PBOT', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure', &
         ptr_topo=this%pbot,t2g_scale_type='topounit')

    this%qbot(begt:endt) = spval
    call hist_addfld1d (fname='QBOT', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_topo=this%qbot,t2g_scale_type='topounit')

    this%rhbot(begt:endt) = spval
     call hist_addfld1d (fname='RH', units='%',  &
          avgflag='A', long_name='atmospheric relative humidity', &
           ptr_topo=this%rhbot, default='inactive')

    this%windbot(begt:endt) = spval
    call hist_addfld1d (fname='WIND', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_topo=this%windbot,t2g_scale_type='topounit')

    this%zbot(begt:endt) = spval
    call hist_addfld1d (fname='ZBOT', units='m',  &
         avgflag='A', long_name='atmospheric reference height', &
         ptr_topo=this%zbot,t2g_scale_type='topounit')

    this%pco2bot(begt:endt) = spval
    call hist_addfld1d (fname='PCO2', units='Pa',  &
         avgflag='A', long_name='atmospheric partial pressure of CO2', &
         ptr_topo=this%pco2bot,t2g_scale_type='topounit')

    if (use_lch4) then
       this%pch4bot(begt:endt) = spval
       call hist_addfld1d (fname='PCH4', units='Pa',  &
            avgflag='A', long_name='atmospheric partial pressure of CH4', &
            ptr_topo=this%pch4bot,t2g_scale_type='topounit')
    end if
  end subroutine init_top_as

  !-----------------------------------------------------------------------
  subroutine clean_top_as(this, begt, endt)
    class(topounit_atmospheric_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    deallocate(this%tbot)
    deallocate(this%thbot)
    deallocate(this%pbot)
    deallocate(this%rhobot)
    deallocate(this%qbot)
    deallocate(this%rhbot)
    deallocate(this%ubot)
    deallocate(this%vbot)
    deallocate(this%wsresp)
    deallocate(this%tau_est)
    deallocate(this%ugust)
    deallocate(this%windbot)
    deallocate(this%zbot)
    deallocate(this%po2bot)
    deallocate(this%pco2bot)
    deallocate(this%pc13o2bot)
    deallocate(this%pch4bot)
    if (use_fates) then
      deallocate(this%rh24h)
      deallocate(this%wind24h)
    end if
  end subroutine clean_top_as

  !-----------------------------------------------------------------------
  subroutine init_acc_buffer_top_as (this, bounds)
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for accumulated fields for atmospheric state
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    use elm_varcon  , only : spval
     use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_state) :: this
    type(bounds_type), intent(in) :: bounds
    !---------------------------------------------------------------------
    ! Accumulator state variables used by FATES
    if (use_fates) then
        call init_accum_field (name='RH24H', units='%', &
             desc='24hr running mean of relative humidity', accum_type='runmean', accum_period=-1, &
              subgrid_type='topounit', numlev=1, init_value=0._r8)
        call init_accum_field (name='WIND24H', units='m/s', &
             desc='24hr running mean of wind', accum_type='runmean', accum_period=-1, &
              subgrid_type='topounit', numlev=1, init_value=0._r8)
    end if
  end subroutine init_acc_buffer_top_as

  !-----------------------------------------------------------------------
  subroutine init_acc_vars_top_as(this, bounds)
    ! !DESCRIPTION:
    ! Initialize variables associated with atmospheric state
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES
     use accumulMod       , only : extract_accum_field
     use elm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_state) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begt, endt
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslt(:)  ! temporary
    !---------------------------------------------------------------------

    begt = bounds%begt; endt = bounds%endt

    ! Allocate needed dynamic memory for single level topounit field
    allocate(rbufslt(begt:endt), stat=ier)
    if (ier/=0) then
        write(iulog,*)' in '
        call endrun(msg="extract_accum_hist allocation error for rbufslt"//&
              errMsg(__FILE__, __LINE__))
    endif

    ! Determine time step
     nstep = get_nstep()

    if (use_fates) then
        call extract_accum_field ('RH24H', rbufslt, nstep)
       this%rh24h(begt:endt) = rbufslt(begt:endt)

        call extract_accum_field ('WIND24H', rbufslt, nstep)
       this%wind24h(begt:endt) = rbufslt(begt:endt)
    end if

    deallocate(rbufslt)
  end subroutine init_acc_vars_top_as

  !-----------------------------------------------------------------------
   subroutine update_acc_vars_top_as (this, bounds)
     !
     ! USES
      use elm_time_manager, only : get_nstep
     use accumulMod      , only : update_accum_field, extract_accum_field
     !
     ! !ARGUMENTS:
     class(topounit_atmospheric_state)    :: this
     type(bounds_type)      , intent(in) :: bounds
     !
     ! !LOCAL VARIABLES:
     integer :: g,t,c,p                   ! indices
     integer :: dtime                     ! timestep size [seconds]
     integer :: nstep                     ! timestep number
     integer :: ier                       ! error status
     integer :: begt, endt
     real(r8), pointer :: rbufslt(:)      ! temporary single level - topounit level
     !---------------------------------------------------------------------

     begt = bounds%begt; endt = bounds%endt

      nstep = get_nstep()

     ! Allocate needed dynamic memory for single level topounit field

     allocate(rbufslt(begt:endt), stat=ier)
     if (ier/=0) then
         write(iulog,*)'update_accum_hist allocation error for rbufslt'
          call endrun(msg=errMsg(__FILE__, __LINE__))
     endif

     ! Accumulate 24-hour running mean of relative humidity
     do t = begt,endt
        rbufslt(t) = this%rhbot(t)
     end do
     if (use_fates) then
        call update_accum_field  ('RH24H', rbufslt, nstep)
        call extract_accum_field ('RH24H', this%rh24h, nstep)
     end if

     ! Accumulate 24-hour running mean of wind speed
     do t = begt,endt
        rbufslt(t) = this%windbot(t)
     end do
     if (use_fates) then
        call update_accum_field  ('WIND24H', rbufslt, nstep)
        call extract_accum_field ('WIND24H', this%wind24h, nstep)
     end if

     deallocate(rbufslt)
    end subroutine update_acc_vars_top_as

  !-----------------------------------------------------------------------
  subroutine init_top_af(this, begt, endt)
    class(topounit_atmospheric_flux) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    ! Allocate for atmospheric flux forcing variables, initialize to special value
    allocate(this%rain     (begt:endt))          ; this%rain      (:) = spval
    allocate(this%snow     (begt:endt))          ; this%snow      (:) = spval
    allocate(this%solad    (begt:endt, numrad))  ; this%solad     (:,:) = spval
    allocate(this%solai    (begt:endt, numrad))  ; this%solai     (:,:) = spval
    allocate(this%solar    (begt:endt))          ; this%solar     (:) = spval
    allocate(this%lwrad    (begt:endt))          ; this%lwrad     (:) = spval
    if (use_fates) then
      allocate(this%prec24h  (begt:endt)) ; this%prec24h   (:) = spval
    end if
    if (use_cn) then
      allocate(this%prec10d  (begt:endt)) ; this%prec10d   (:) = spval
      allocate(this%prec60d  (begt:endt)) ; this%prec60d   (:) = spval
    end if
    allocate(this%fsd24h   (begt:endt))          ; this%fsd24h    (:) = spval
    allocate(this%fsd240h  (begt:endt))          ; this%fsd240h   (:) = spval
    allocate(this%fsi24h   (begt:endt))          ; this%fsi24h    (:) = spval
    allocate(this%fsi240h  (begt:endt))          ; this%fsi240h   (:) = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_af
    !-----------------------------------------------------------------------
    this%rain(begt:endt) = spval
    call hist_addfld1d (fname='RAIN', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_topo=this%rain,t2g_scale_type='topounit')

    this%snow(begt:endt) = spval
    call hist_addfld1d (fname='SNOW', units='mm/s',  &
         avgflag='A', long_name='atmospheric snow', &
         ptr_topo=this%snow,t2g_scale_type='topounit')

    this%solar(begt:endt) = spval
    call hist_addfld1d (fname='FSDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_topo=this%solar,t2g_scale_type='topounit')

    this%lwrad(begt:endt) = spval
    call hist_addfld1d (fname='FLDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_topo=this%lwrad,t2g_scale_type='topounit')
  
  end subroutine init_top_af

  !-----------------------------------------------------------------------
  subroutine clean_top_af(this, begt, endt)
    class(topounit_atmospheric_flux) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    deallocate(this%rain)
    deallocate(this%snow)
    deallocate(this%solad)
    deallocate(this%solai)
    deallocate(this%solar)
    deallocate(this%lwrad)
    if (use_fates) then
      deallocate(this%prec24h)
    end if
    if (use_cn) then
      deallocate(this%prec10d)
      deallocate(this%prec60d)
    end if
    deallocate(this%fsd24h)
    deallocate(this%fsd240h)
    deallocate(this%fsi24h)
    deallocate(this%fsi240h)
  end subroutine clean_top_af

  !-----------------------------------------------------------------------
  subroutine init_acc_buffer_top_af (this, bounds)
    ! !DESCRIPTION:
    ! Initialize accumulation buffer accumulated fields for atmospheric flux
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    use elm_varcon  , only : spval
     use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_flux) :: this
    type(bounds_type), intent(in) :: bounds
    !---------------------------------------------------------------------
    ! Accumulator flux variables used by CN
    if (use_cn) then
       call init_accum_field (name='PREC10D', units='MM H2O/S', &
          desc='10-day running mean of total precipitation', accum_type='runmean', accum_period=-10, &
           subgrid_type='topounit', numlev=1, init_value=0._r8)
       call init_accum_field (name='PREC60D', units='MM H2O/S', &
          desc='60-day running mean of total precipitation', accum_type='runmean', accum_period=-60, &
           subgrid_type='topounit', numlev=1, init_value=0._r8)
    end if
    ! Accumulator flux variables used by FATES
    if (use_fates) then
        call init_accum_field (name='PREC24H', units='MM H2O/S', &
             desc='24hr running mean of total precipitation', accum_type='runmean', accum_period=-1, &
              subgrid_type='topounit', numlev=1, init_value=0._r8)
    end if
    ! Accumulator variables for radiation fluxes
    call init_accum_field (name='FSD24H', units='W/m2',                                             &
         desc='24hr average of direct solar radiation',  accum_type='runmean', accum_period=-1,    &
         subgrid_type='topounit', numlev=1, init_value=0._r8)

    call init_accum_field (name='FSD240H', units='W/m2',                                            &
         desc='240hr average of direct solar radiation',  accum_type='runmean', accum_period=-10,  &
         subgrid_type='topounit', numlev=1, init_value=0._r8)

    call init_accum_field (name='FSI24H', units='W/m2',                                             &
         desc='24hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-1,   &
         subgrid_type='topounit', numlev=1, init_value=0._r8)

    call init_accum_field (name='FSI240H', units='W/m2',                                            &
         desc='240hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-10, &
         subgrid_type='topounit', numlev=1, init_value=0._r8)
  end subroutine init_acc_buffer_top_af

  !-----------------------------------------------------------------------
  subroutine init_acc_vars_top_af(this, bounds)
    ! !DESCRIPTION:
    ! Initialize variables associated with atmospheric flux
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES
     use accumulMod       , only : extract_accum_field
     use elm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_flux) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begt, endt
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslt(:)  ! temporary
    !---------------------------------------------------------------------

    begt = bounds%begt; endt = bounds%endt

    ! Allocate needed dynamic memory for single level topounit field
    allocate(rbufslt(begt:endt), stat=ier)
    if (ier/=0) then
        write(iulog,*)' in '
        call endrun(msg="extract_accum_hist allocation error for rbufslt"//&
              errMsg(__FILE__, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('FSD24H', rbufslt, nstep)
    this%fsd24h(begt:endt) = rbufslt(begt:endt)

    call extract_accum_field ('FSD240H', rbufslt, nstep)
    this%fsd240h(begt:endt) = rbufslt(begt:endt)

    call extract_accum_field ('FSI24H', rbufslt, nstep)
    this%fsi24h(begt:endt) = rbufslt(begt:endt)

    call extract_accum_field ('FSI240H', rbufslt, nstep)
    this%fsi240h(begt:endt) = rbufslt(begt:endt)

    if (use_cn) then
        call extract_accum_field ('PREC10D', rbufslt, nstep)
       this%prec10d(begt:endt) = rbufslt(begt:endt)

        call extract_accum_field ('PREC60D', rbufslt, nstep)
       this%prec60d(begt:endt) = rbufslt(begt:endt)
    end if

    if (use_fates) then
        call extract_accum_field ('PREC24H', rbufslt, nstep)
       this%prec24h(begt:endt) = rbufslt(begt:endt)
    end if

    deallocate(rbufslt)
  end subroutine init_acc_vars_top_af

  !-----------------------------------------------------------------------
  subroutine update_acc_vars_top_af (this, bounds)
    !
    ! USES
    use elm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(topounit_atmospheric_flux)    :: this
    type(bounds_type)      , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,t,c,p                   ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begt, endt
    real(r8), pointer :: rbufslt(:)      ! temporary single level - topounit level
    !---------------------------------------------------------------------

    begt = bounds%begt; endt = bounds%endt

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level topounit field

    allocate(rbufslt(begt:endt), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbufslt'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract forc_solad24 & forc_solad240
    do t = begt,endt
       rbufslt(t) = this%solad(t,1)
    end do
    call update_accum_field  ('FSD240H', rbufslt              , nstep)
    call extract_accum_field ('FSD240H', this%fsd240h         , nstep)
    call update_accum_field  ('FSD24H' , rbufslt              , nstep)
    call extract_accum_field ('FSD24H' , this%fsd24h          , nstep)

    ! Accumulate and extract forc_solai24 & forc_solai240
    do t = begt,endt
       rbufslt(t) = this%solai(t,1)
    end do
    call update_accum_field  ('FSI24H' , rbufslt              , nstep)
    call extract_accum_field ('FSI24H' , this%fsi24h          , nstep)
    call update_accum_field  ('FSI240H', rbufslt              , nstep)
    call extract_accum_field ('FSI240H', this%fsi240h         , nstep)

    ! Accumulate and extract total precip
    do t = begt,endt
       rbufslt(t) = this%rain(t) + this%snow(t)
    end do
    if (use_cn) then
       ! Accumulate and extract PREC60D (accumulates total precipitation as 60-day running mean)
       call update_accum_field  ('PREC60D', rbufslt, nstep)
       call extract_accum_field ('PREC60D', this%prec60d, nstep)

       ! Accumulate and extract PREC10D (accumulates total precipitation as 10-day running mean)
       call update_accum_field  ('PREC10D', rbufslt, nstep)
       call extract_accum_field ('PREC10D', this%prec10d, nstep)
    end if

    if (use_fates) then
       call update_accum_field  ('PREC24H', rbufslt, nstep)
       call extract_accum_field ('PREC24H', this%prec24h, nstep)
    end if

    deallocate(rbufslt)

  end subroutine update_acc_vars_top_af

  !-----------------------------------------------------------------------
  subroutine init_top_es(this, begt, endt)
    class(topounit_energy_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    allocate(this%t_rad   (begt:endt)) ; this%t_rad   (:) = spval
    allocate(this%eflx_lwrad_out_topo   (begt:endt)) ; this%eflx_lwrad_out_topo   (:) = spval
    allocate(this%t_grnd  (begt:endt)) ; this%t_grnd  (:) = spval
	
	allocate(this%heat1_tgu                (begt:endt))                      ; this%heat1_tgu                (:)   = spval
    allocate(this%heat2_tgu                (begt:endt))                      ; this%heat2_tgu                (:)   = spval
    allocate(this%liquid_water_temp1_tgu   (begt:endt))                      ; this%liquid_water_temp1_tgu   (:)   = spval
    allocate(this%liquid_water_temp2_tgu   (begt:endt))                      ; this%liquid_water_temp2_tgu   (:)   = spval
	
	!this%heat1_tgu(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_HEAT1',  units='J/m^2',  &
    !     avgflag='A', long_name='initial topounit total heat content', &
    !     ptr_topo=this%heat1_tgu)

    !this%heat2_tgu(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_HEAT2',  units='J/m^2',  &
    !     avgflag='A', long_name='post land cover change total heat content', &
    !     ptr_topo=this%heat2_tgu, default='inactive')  

    !this%liquid_water_temp1_tgu(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_LIQUID_WATER_TEMP1', units='K', &
    !     avgflag='A', long_name='initial topounit weighted average liquid water temperature', &
    !     ptr_topo=this%liquid_water_temp1_tgu, default='inactive')
	
  end subroutine init_top_es

  !-----------------------------------------------------------------------
  subroutine restart_top_es(this, bounds, ncid, flag)
    class(topounit_energy_state) :: this
    type(bounds_type), intent(in) :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*), intent(in) :: flag

    logical :: readvar ! determine if variable is on initial file

    call restartvar(ncid=ncid, flag=flag, varname='TS_TOPO', xtype=ncd_double, &
        dim1name='topounit', long_name='surface radiative temperature', &
        units='K', interpinic_flag='copy', readvar=readvar, data=this%t_rad)

  end subroutine restart_top_es
  
  !-----------------------------------------------------------------------
  subroutine clean_top_es(this, begt, endt)
    class(topounit_energy_state) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    deallocate(this%t_rad    )
    deallocate(this%eflx_lwrad_out_topo    )
    deallocate(this%t_grnd   )
	
	deallocate(this%heat1_tgu)
    deallocate(this%heat2_tgu)
    deallocate(this%liquid_water_temp1_tgu)
    deallocate(this%liquid_water_temp2_tgu)
	
  end subroutine clean_top_es

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell energy flux data structure
  !------------------------------------------------------------------------
  subroutine top_ef_init(this, begt, endt)
    !
    ! !ARGUMENTS:
    class(topounit_energy_flux) :: this
    integer, intent(in) :: begt,endt
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of top_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_dynbal              (begt:endt))                  ; this%eflx_dynbal             (:)   = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_ef
    !-----------------------------------------------------------------------
    this%eflx_dynbal(begt:endt) = spval 
    call hist_addfld1d (fname='TOP_EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_topo=this%eflx_dynbal)

  end subroutine top_ef_init

  !------------------------------------------------------------------------
  subroutine top_ef_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write gridcell energy flux information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(topounit_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------
  end subroutine top_ef_restart

  !------------------------------------------------------------------------
  subroutine top_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_energy_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%eflx_dynbal)
  end subroutine top_ef_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell water state data structure
  !------------------------------------------------------------------------
  subroutine top_ws_init(this, begt, endt)
    !
    ! !ARGUMENTS:
    class(topounit_water_state) :: this
    integer, intent(in) :: begt,endt
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of top_ws
    !-----------------------------------------------------------------------
    allocate(this%liq1           (begt:endt))       ; this%liq1           (:)   = spval
    allocate(this%liq2           (begt:endt))       ; this%liq2           (:)   = spval
    allocate(this%ice1           (begt:endt))       ; this%ice1           (:)   = spval
    allocate(this%ice2           (begt:endt))       ; this%ice2           (:)   = spval
    allocate(this%tws            (begt:endt))       ; this%tws            (:)   = spval
    allocate(this%tws_month_beg  (begt:endt))       ; this%tws_month_beg  (:)   = spval
    allocate(this%tws_month_end  (begt:endt))       ; this%tws_month_end  (:)   = spval
    allocate(this%begwb          (begt:endt))       ; this%begwb          (:)   = spval
    allocate(this%endwb          (begt:endt))       ; this%endwb          (:)   = spval
    allocate(this%errh2o         (begt:endt))       ; this%errh2o         (:)   = spval
    allocate(this%beg_h2ocan     (begt:endt))       ; this%beg_h2ocan     (:)   = spval
    allocate(this%beg_h2osno     (begt:endt))       ; this%beg_h2osno     (:)   = spval
    allocate(this%beg_h2osfc     (begt:endt))       ; this%beg_h2osfc     (:)   = spval
    allocate(this%beg_h2osoi_liq (begt:endt))       ; this%beg_h2osoi_liq (:)   = spval
    allocate(this%beg_h2osoi_ice (begt:endt))       ; this%beg_h2osoi_ice (:)   = spval
    allocate(this%end_h2ocan     (begt:endt))       ; this%end_h2ocan     (:)   = spval
    allocate(this%end_h2osno     (begt:endt))       ; this%end_h2osno     (:)   = spval
    allocate(this%end_h2osfc     (begt:endt))       ; this%end_h2osfc     (:)   = spval
    allocate(this%end_h2osoi_liq (begt:endt))       ; this%end_h2osoi_liq (:)   = spval
    allocate(this%end_h2osoi_ice (begt:endt))       ; this%end_h2osoi_ice (:)   = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_ws
    !-----------------------------------------------------------------------
    this%liq1(begt:endt) = spval
    call hist_addfld1d (fname='TOP_LIQ1',  units='mm',  &
         avgflag='A', long_name='initial topounit total liq content', &
         ptr_topo=this%liq1)

    this%liq2(begt:endt) = spval
    call hist_addfld1d (fname='TOP_LIQ2',  units='mm',  &  
         avgflag='A', long_name='post landuse change topounit total liq content', &              
         ptr_topo=this%liq2, default='inactive')     

    this%ice1(begt:endt) = spval
    call hist_addfld1d (fname='TOP_ICE1',  units='mm',  &  
         avgflag='A', long_name='initial topounit total ice content', &              
         ptr_topo=this%ice1)     

    this%ice2(begt:endt) = spval
    call hist_addfld1d (fname='TOP_ICE2',  units='mm',  &  
         avgflag='A', long_name='post land cover change total ice content', &              
         ptr_topo=this%ice2, default='inactive')

    this%tws(begt:endt) = spval
    call hist_addfld1d (fname='TOP_TWS',  units='mm',  &
         avgflag='A', long_name='total water storage', &
         ptr_topo=this%tws)

    this%tws_month_beg(begt:endt) = spval
    call hist_addfld1d (fname='TOP_TWS_MONTH_BEGIN',  units='mm',  &
         avgflag='I', long_name='total water storage at the beginning of a month', &
         ptr_topo=this%tws_month_beg)

    this%tws_month_end(begt:endt) = spval
    call hist_addfld1d (fname='TOP_TWS_MONTH_END',  units='mm',  &
         avgflag='I', long_name='total water storage at the end of a month', &
         ptr_topo=this%tws_month_end)
  end subroutine top_ws_init

  !------------------------------------------------------------------------
  subroutine top_ws_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write topounit water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(topounit_water_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------
    call restartvar(ncid=ncid, flag=flag, varname='TOP_TWS_MONTH_BEGIN', xtype=ncd_double,  &
         dim1name='topounit', &
         long_name='total water storage at the beginning of a month', units='mm', &
          interpinic_flag='interp', readvar=readvar, data=this%tws_month_beg)
  
  end subroutine top_ws_restart
  
  !------------------------------------------------------------------------
  subroutine top_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_water_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%liq1)
    deallocate(this%liq2)
    deallocate(this%ice1)
    deallocate(this%ice2)
    deallocate(this%tws )
  end subroutine top_ws_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell water flux data structure
  !------------------------------------------------------------------------
  subroutine top_wf_init(this, begt, endt, bounds)
    !
    ! !ARGUMENTS:
    class(topounit_water_flux) :: this
    integer, intent(in) :: begt,endt
    type(bounds_type)    , intent(in)  :: bounds
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of top_wf
    !-----------------------------------------------------------------------
    allocate(this%qflx_liq_dynbal      (begt:endt))              ; this%qflx_liq_dynbal      (:)   = spval
    allocate(this%qflx_ice_dynbal      (begt:endt))              ; this%qflx_ice_dynbal      (:)   = spval

    this%qflx_liq_dynbal_dribbler = annual_flux_dribbler_topounit( &
         bounds = bounds, &
         name = 'qflx_liq_dynbal', &
         units = 'mm H2O')

    this%qflx_ice_dynbal_dribbler = annual_flux_dribbler_topounit( &
         bounds = bounds, &
         name = 'qflx_ice_dynbal', &
         units = 'mm H2O')

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_wf
    !-----------------------------------------------------------------------
    this%qflx_liq_dynbal(begt:endt) = spval
    call hist_addfld1d (fname='TOP_QFLX_LIQ_DYNBAL',  units='mm/s',  &  
         avgflag='A', long_name='liq dynamic land cover change conversion runoff flux', &              
         ptr_topo=this%qflx_liq_dynbal)     

    this%qflx_ice_dynbal(begt:endt) = spval
    call hist_addfld1d (fname='TOP_QFLX_ICE_DYNBAL',  units='mm/s',  &
         avgflag='A', long_name='ice dynamic land cover change conversion runoff flux', &                                   
         ptr_topo=this%qflx_ice_dynbal)
  
  end subroutine top_wf_init

  !------------------------------------------------------------------------
  subroutine top_wf_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write gridcell water flux information to/from restart file.
    !
    ! !ARGUMENTS:
    class(topounit_water_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call this%qflx_liq_dynbal_dribbler%Restart(bounds, ncid, flag)
    call this%qflx_ice_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine top_wf_restart

  !------------------------------------------------------------------------
  subroutine top_wf_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_water_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%qflx_liq_dynbal)
    deallocate(this%qflx_ice_dynbal)

  end subroutine top_wf_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean topounit carbon state data structure
  !------------------------------------------------------------------------
  subroutine top_cs_init(this, begt, endt,carbon_type)
    !
    ! !ARGUMENTS:
    class(topounit_carbon_state) :: this
    integer, intent(in) :: begt,endt
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']    

    !
    ! !LOCAL VARIABLES:
    integer :: t
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of top_cs
    !-----------------------------------------------------------------------
    allocate(this%seedc                 (begt:endt));     this%seedc                 (:) = spval
    allocate(this%tcs_month_beg         (begt:endt));     this%tcs_month_beg         (:) = spval
    allocate(this%tcs_month_end         (begt:endt));     this%tcs_month_end         (:) = spval
    allocate(this%begcb                 (begt:endt));     this%begcb                 (:) = spval
    allocate(this%endcb                 (begt:endt));     this%endcb                 (:) = spval
    allocate(this%errcb                 (begt:endt));     this%errcb                 (:) = spval

    allocate(this%beg_totc              (begt:endt));     this%beg_totc              (:) = spval
    allocate(this%beg_totpftc           (begt:endt));     this%beg_totpftc           (:) = spval
    allocate(this%beg_cwdc              (begt:endt));     this%beg_cwdc              (:) = spval
    allocate(this%beg_totsomc           (begt:endt));     this%beg_totsomc           (:) = spval
    allocate(this%beg_totlitc           (begt:endt));     this%beg_totlitc           (:) = spval
    allocate(this%beg_totprodc          (begt:endt));     this%beg_totprodc          (:) = spval
    allocate(this%beg_ctrunc            (begt:endt));     this%beg_ctrunc            (:) = spval
    allocate(this%beg_cropseedc_deficit (begt:endt));     this%beg_cropseedc_deficit (:) = spval

    allocate(this%end_totc              (begt:endt));     this%end_totc              (:) = spval
    allocate(this%end_totpftc           (begt:endt));     this%end_totpftc           (:) = spval
    allocate(this%end_cwdc              (begt:endt));     this%end_cwdc              (:) = spval
    allocate(this%end_totsomc           (begt:endt));     this%end_totsomc           (:) = spval
    allocate(this%end_totlitc           (begt:endt));     this%end_totlitc           (:) = spval
    allocate(this%end_totprodc          (begt:endt));     this%end_totprodc          (:) = spval
    allocate(this%end_ctrunc            (begt:endt));     this%end_ctrunc            (:) = spval
    allocate(this%end_cropseedc_deficit (begt:endt));     this%end_cropseedc_deficit (:) = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_cs
    !-----------------------------------------------------------------------
    if (.not. use_fates) then
       if (carbon_type == 'c12') then
          this%seedc(begt:endt) = spval
          call hist_addfld1d (fname='TOP_SEEDC', units='gC/m^2', &
               avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
               ptr_topo=this%seedc)
       end if 
       if (carbon_type == 'c13') then
          this%seedc(begt:endt) = spval
          call hist_addfld1d (fname='TOP_C13_SEEDC', units='gC/m^2', &
               avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
               ptr_topo=this%seedc)
       end if 
       if (carbon_type == 'c14') then
          this%seedc(begt:endt) = spval
          call hist_addfld1d (fname='TOP_C14_SEEDC', units='gC/m^2', &
               avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
               ptr_topo=this%seedc)
       end if 
    end if
    
    this%tcs_month_beg(begt:endt) = spval
    call hist_addfld1d (fname='TOP_TCS_MONTH_BEGIN',  units='gC/m^2',  &
         avgflag='I', long_name='total carbon storage at the beginning of a month', &
         ptr_topo=this%tcs_month_beg)

    this%tcs_month_end(begt:endt) = spval
    call hist_addfld1d (fname='TOP_TCS_MONTH_END',  units='gC/m^2',  &
         avgflag='I', long_name='total carbon storage at the end of a month', &
         ptr_topo=this%tcs_month_end)

    call hist_addfld1d (fname='TOP_CMASS_BALANCE_ERROR',  units='gC/m^2',  &
         avgflag='A', long_name='Topounit carbon mass balance error', &
         ptr_topo=this%errcb)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of top_cs
    !-----------------------------------------------------------------------
    do t = begt, endt
       this%seedc(t) = 0._r8
    end do
    
  end subroutine top_cs_init

  !------------------------------------------------------------------------
  subroutine top_cs_restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write topounit water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(topounit_carbon_state) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: flag
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------
    call restartvar(ncid=ncid, flag=flag, varname='TOP_TCS_MONTH_BEGIN', xtype=ncd_double,  &
         dim1name='topounit', &
         long_name='total carbon storage at the beginning of a month', units='gC/m^2', &
          interpinic_flag='interp', readvar=readvar, data=this%tcs_month_beg)

    call restartvar(ncid=ncid, flag=flag, varname='TOP_TCS_MONTH_END', xtype=ncd_double,  &
         dim1name='topounit', &
         long_name='total carbon storage at the end of a month', units='gC/m^2', &
          interpinic_flag='interp', readvar=readvar, data=this%tcs_month_end)
  end subroutine top_cs_restart

  !------------------------------------------------------------------------
  subroutine top_cs_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_carbon_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%seedc)
    deallocate(this%begcb)
    deallocate(this%endcb)
    deallocate(this%errcb)

  end subroutine top_cs_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean topounit carbon flux data structure
  !------------------------------------------------------------------------
  subroutine top_cf_init(this, begt, endt, carbon_type)
    !
    ! !ARGUMENTS:
    class(topounit_carbon_flux) :: this
    integer, intent(in) :: begt,endt
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer :: t
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of top_cf
    !-----------------------------------------------------------------------
    allocate(this%dwt_seedc_to_leaf            (begt:endt)) ; this%dwt_seedc_to_leaf         (:) = spval
    allocate(this%dwt_seedc_to_deadstem        (begt:endt)) ; this%dwt_seedc_to_deadstem     (:) = spval
    allocate(this%dwt_conv_cflux               (begt:endt)) ; this%dwt_conv_cflux            (:) = spval
    allocate(this%dwt_conv_cflux_dribbled      (begt:endt)) ; this%dwt_conv_cflux_dribbled   (:) = spval
    allocate(this%dwt_prod10c_gain             (begt:endt)) ; this%dwt_prod10c_gain          (:) = spval
    allocate(this%dwt_prod100c_gain            (begt:endt)) ; this%dwt_prod100c_gain         (:) = spval
    allocate(this%hrv_deadstemc_to_prod10c     (begt:endt)) ; this%hrv_deadstemc_to_prod10c  (:) = spval
    allocate(this%hrv_deadstemc_to_prod100c    (begt:endt)) ; this%hrv_deadstemc_to_prod100c (:) = spval
    allocate(this%cinputs                      (begt:endt)) ; this%cinputs                   (:) = spval
    allocate(this%coutputs                     (begt:endt)) ; this%coutputs                  (:) = spval
    allocate(this%gpp                          (begt:endt)) ; this%gpp                       (:) = spval
    allocate(this%er                           (begt:endt)) ; this%er                        (:) = spval
    allocate(this%fire_closs                   (begt:endt)) ; this%fire_closs                (:) = spval
    allocate(this%prod1_loss                   (begt:endt)) ; this%prod1_loss                (:) = spval
    allocate(this%prod10_loss                  (begt:endt)) ; this%prod10_loss               (:) = spval
    allocate(this%prod100_loss                 (begt:endt)) ; this%prod100_loss              (:) = spval
    allocate(this%hrv_xsmrpool_to_atm          (begt:endt)) ; this%hrv_xsmrpool_to_atm       (:) = spval
    allocate(this%som_c_leached                (begt:endt)) ; this%som_c_leached             (:) = spval
    allocate(this%somc_yield                   (begt:endt)) ; this%somc_yield                (:) = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_cf
    !-----------------------------------------------------------------------
    ! no history fields or cold-start initialization for gridcell carbon flux
    ! if using fates
    if (use_fates) then
       return
    end if

    !if (carbon_type == 'c12') then

    !   this%dwt_seedc_to_leaf(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_DWT_SEEDC_TO_LEAF', units='gC/m^2/s', &
    !        avgflag='A', long_name='seed source to patch-level leaf', &
    !        ptr_topo=this%dwt_seedc_to_leaf, default='inactive')

    !   this%dwt_seedc_to_deadstem(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_DWT_SEEDC_TO_DEADSTEM', units='gC/m^2/s', &
    !        avgflag='A', long_name='seed source to patch-level deadstem', &
    !        ptr_topo=this%dwt_seedc_to_deadstem, default='inactive')

    !   this%dwt_conv_cflux(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_DWT_CONV_CFLUX', units='gC/m^2/s', &
    !        avgflag='A', &
    !        long_name='conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
    !        ptr_topo=this%dwt_conv_cflux)

    !   this%dwt_conv_cflux_dribbled(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_DWT_CONV_CFLUX_DRIBBLED', units='gC/m^2/s', &
    !        avgflag='A', &
    !        long_name='conversion C flux (immediate loss to atm), dribbled throughout the year', &
    !        ptr_topo=this%dwt_conv_cflux_dribbled)

    !   this%dwt_prod10c_gain(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_DWT_PROD10C_GAIN', units='gC/m^2/s', &
    !        avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
    !        ptr_col=this%dwt_prod10c_gain, default='inactive')

    !   this%dwt_prod100c_gain(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_DWT_PROD100C_GAIN', units='gC/m^2/s', &
    !        avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
    !        ptr_col=this%dwt_prod100c_gain, default='inactive')

    !   this%hrv_deadstemc_to_prod10c(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_HRV_DEADSTEM_TO_PROD10C', units='gC/m^2/s', &
    !        avgflag='A', long_name='dead stem harvest to 10-yr wood product pool', &
    !        ptr_col=this%hrv_deadstemc_to_prod10c, default='inactive')

    !   this%hrv_deadstemc_to_prod100c(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_HRV_DEADSTEM_TO_PROD100C', units='gC/m^2/s', &
    !        avgflag='A', long_name='dead stem harvest to 100-yr wood product pool', &
    !        ptr_col=this%hrv_deadstemc_to_prod100c, default='inactive')
    !end if

    !if ( carbon_type == 'c13' ) then

    !   this%dwt_seedc_to_leaf(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_DWT_SEEDC_TO_LEAF', units='gC13/m^2/s', &
    !        avgflag='A', long_name='C13 seed source to patch-level leaf', &
    !        ptr_topo=this%dwt_seedc_to_leaf, default='inactive')

    !   this%dwt_seedc_to_deadstem(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_DWT_SEEDC_TO_DEADSTEM', units='gC13/m^2/s', &
    !        avgflag='A', long_name='C13 seed source to patch-level deadstem', &
    !        ptr_topo=this%dwt_seedc_to_deadstem, default='inactive')

    !   this%dwt_conv_cflux(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_DWT_CONV_CFLUX', units='gC13/m^2/s', &
    !        avgflag='A', &
    !        long_name='C13 conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
    !        ptr_topo=this%dwt_conv_cflux)

    !   this%dwt_conv_cflux_dribbled(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_DWT_CONV_CFLUX_DRIBBLED', units='gC13/m^2/s', &
    !        avgflag='A', &
    !        long_name='C13 conversion C flux (immediate loss to atm), dribbled throughout the year', &
    !        ptr_topo=this%dwt_conv_cflux_dribbled)

    !   this%dwt_seedc_to_deadstem(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_DWT_SEEDC_TO_DEADSTEM', units='gC13/m^2/s', &
    !        avgflag='A', long_name='C13 seed source to patch-level deadstem', &
    !        ptr_topo=this%dwt_seedc_to_deadstem, default='inactive')

    !   this%dwt_prod10c_gain(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_DWT_PROD10C_GAIN', units='gC13/m^2/s', &
    !        avgflag='A', long_name='C13 landcover change-driven addition to 10-yr wood product pool', &
    !        ptr_col=this%dwt_prod10c_gain, default='inactive')

    !   this%dwt_prod100c_gain(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_DWT_PROD100C_GAIN', units='gC13/m^2/s', &
    !        avgflag='A', long_name='C13 landcover change-driven addition to 100-yr wood product pool', &
    !        ptr_col=this%dwt_prod100c_gain, default='inactive')

    !   this%hrv_deadstemc_to_prod10c(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_HRV_DEADSTEM_TO_PROD10C', units='gC13/m^2/s', &
    !        avgflag='A', long_name='C13 dead stem harvest to 10-yr wood product pool', &
    !        ptr_col=this%hrv_deadstemc_to_prod10c, default='inactive')

    !   this%hrv_deadstemc_to_prod100c(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C13_HRV_DEADSTEM_TO_PROD100C', units='gC13/m^2/s', &
    !        avgflag='A', long_name='C13 dead stem harvest to 100-yr wood product pool', &
    !        ptr_col=this%hrv_deadstemc_to_prod100c, default='inactive')
    !endif

    !if (carbon_type == 'c14') then

    !   this%dwt_seedc_to_leaf(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_DWT_SEEDC_TO_LEAF', units='gC14/m^2/s', &
    !        avgflag='A', long_name='C14 seed source to patch-level leaf', &
    !        ptr_topo=this%dwt_seedc_to_leaf, default='inactive')

    !   this%dwt_seedc_to_deadstem(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_DWT_SEEDC_TO_DEADSTEM', units='gC14/m^2/s', &
    !        avgflag='A', long_name='C14 seed source to patch-level deadstem', &
    !        ptr_topo=this%dwt_seedc_to_deadstem, default='inactive')

    !   this%dwt_conv_cflux(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_DWT_CONV_CFLUX', units='gC14/m^2/s', &
    !        avgflag='A', &
    !        long_name='C14 conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
    !        ptr_topo=this%dwt_conv_cflux)

    !   this%dwt_conv_cflux_dribbled(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_DWT_CONV_CFLUX_DRIBBLED', units='gC14/m^2/s', &
    !        avgflag='A', &
    !        long_name='C14 conversion C flux (immediate loss to atm), dribbled throughout the year', &
    !        ptr_topo=this%dwt_conv_cflux_dribbled)

    !   this%dwt_seedc_to_deadstem(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_DWT_SEEDC_TO_DEADSTEM', units='gC14/m^2/s', &
    !        avgflag='A', long_name='C14 seed source to patch-level deadstem', &
    !        ptr_topo=this%dwt_seedc_to_deadstem, default='inactive')

    !   this%dwt_prod10c_gain(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_DWT_PROD10C_GAIN', units='gC14/m^2/s', &
    !        avgflag='A', long_name='C14 landcover change-driven addition to 10-yr wood product pool', &
    !        ptr_col=this%dwt_prod10c_gain, default='inactive')

    !   this%dwt_prod100c_gain(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_DWT_PROD100C_GAIN', units='gC14/m^2/s', &
    !        avgflag='A', long_name='C14 landcover change-driven addition to 100-yr wood product pool', &
    !        ptr_col=this%dwt_prod100c_gain, default='inactive')

    !   this%hrv_deadstemc_to_prod10c(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_HRV_DEADSTEM_TO_PROD10C', units='gC14/m^2/s', &
    !        avgflag='A', long_name='C14 dead stem harvest to 10-yr wood product pool', &
    !        ptr_col=this%hrv_deadstemc_to_prod10c, default='inactive')

    !   this%hrv_deadstemc_to_prod100c(begt:endt) = spval
    !   call hist_addfld1d (fname='TOP_C14_HRV_DEADSTEM_TO_PROD100C', units='gC14/m^2/s', &
    !        avgflag='A', long_name='C14 dead stem harvest to 100-yr wood product pool', &
    !        ptr_col=this%hrv_deadstemc_to_prod100c, default='inactive')
    !endif
    
    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of top_cf
    !-----------------------------------------------------------------------
    do t = begt, endt
       this%dwt_prod10c_gain(t)          = 0._r8
       this%dwt_prod100c_gain(t)         = 0._r8
       this%hrv_deadstemc_to_prod10c(t)  = 0._r8
       this%hrv_deadstemc_to_prod100c(t) = 0._r8
       this%cinputs(t)                   = 0._r8
       this%coutputs(t)                  = 0._r8
    end do
    
  end subroutine top_cf_init
  
  !-----------------------------------------------------------------------
  subroutine top_cf_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(topounit_carbon_flux)      :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: t          ! indices
    !-----------------------------------------------------------------------

    ! set column-level conversion and product pool fluxes
    ! to 0 at the beginning of every timestep

    do t = bounds%begt, bounds%endt
       this%dwt_seedc_to_leaf(t)         = 0._r8
       this%dwt_seedc_to_deadstem(t)     = 0._r8
       this%dwt_conv_cflux(t)            = 0._r8
       this%dwt_prod10c_gain(t)          = 0._r8
       this%dwt_prod100c_gain(t)         = 0._r8
       this%hrv_deadstemc_to_prod10c(t)  = 0._r8
       this%hrv_deadstemc_to_prod100c(t) = 0._r8
    end do

  end subroutine top_cf_zerodwt

  !------------------------------------------------------------------------
  subroutine top_cf_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_carbon_flux) :: this
    !------------------------------------------------------------------------

  end subroutine top_cf_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean topounit nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine top_ns_init(this, begt, endt)
    !
    ! !ARGUMENTS:
    class(topounit_nitrogen_state) :: this
    integer, intent(in) :: begt,endt
    !
    ! !LOCAL VARIABLES:
    integer :: t
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of top_ns
    !-----------------------------------------------------------------------
    allocate(this%seedn   (begt:endt));     this%seedn   (:) = spval
    allocate(this%begnb   (begt:endt));     this%begnb   (:) = spval
    allocate(this%endnb   (begt:endt));     this%endnb   (:) = spval
    allocate(this%errnb   (begt:endt));     this%errnb   (:) = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_cf
    !-----------------------------------------------------------------------
    this%seedn(begt:endt) = spval
    call hist_addfld1d (fname='TOP_SEEDN', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs ', &
         ptr_topo=this%seedn, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of top_cf
    !-----------------------------------------------------------------------
    do t = begt, endt
       this%seedn(t) = 0._r8
    end do


  end subroutine top_ns_init
  
  !------------------------------------------------------------------------
  subroutine top_ns_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_nitrogen_state) :: this
    !------------------------------------------------------------------------

  end subroutine top_ns_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean topounit nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine top_nf_init(this, begt, endt)
    !
    ! !ARGUMENTS:
    class(topounit_nitrogen_flux) :: this
    integer, intent(in) :: begt,endt
    !
    ! !LOCAL VARIABLES:
    integer :: t
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of top_nf
    !-----------------------------------------------------------------------
    allocate(this%dwt_seedn_to_leaf     (begt:endt)) ; this%dwt_seedn_to_leaf     (:) = spval
    allocate(this%dwt_seedn_to_deadstem (begt:endt)) ; this%dwt_seedn_to_deadstem (:) = spval
    allocate(this%dwt_conv_nflux        (begt:endt)) ; this%dwt_conv_nflux        (:) = spval
    allocate(this%dwt_seedn_to_npool    (begt:endt)) ; this%dwt_seedn_to_npool    (:) = spval
    allocate(this%dwt_prod10n_gain      (begt:endt)) ; this%dwt_prod10n_gain      (:) = spval
    allocate(this%dwt_prod100n_gain     (begt:endt)) ; this%dwt_prod100n_gain     (:) = spval
    allocate(this%ninputs               (begt:endt)) ; this%ninputs               (:) = spval
    allocate(this%noutputs              (begt:endt)) ; this%noutputs              (:) = spval
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_nf
    !-----------------------------------------------------------------------
    !this%dwt_seedn_to_leaf(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_DWT_SEEDN_TO_LEAF', units='gN/m^2/s', &
    !     avgflag='A', long_name='seed source to patch-level leaf', &
    !     ptr_topo=this%dwt_seedn_to_leaf, default='inactive')

    !this%dwt_seedn_to_deadstem(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_DWT_SEEDN_TO_DEADSTEM', units='gN/m^2/s', &
    !     avgflag='A', long_name='seed source to patch-level deadstem', &
    !     ptr_topo=this%dwt_seedn_to_deadstem, default='inactive')

    !this%dwt_conv_nflux(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_DWT_CONV_NFLUX', units='gN/m^2/s', &
    !     avgflag='A', &
    !     long_name='conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
    !     ptr_topo=this%dwt_conv_nflux)

    !this%dwt_seedn_to_npool(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_DWT_SEEDN_TO_NPOOL', units='gN/m^2/s', &
    !     avgflag='A', long_name='seed source to PFT-level npool', &
    !     ptr_topo=this%dwt_seedn_to_npool, default='inactive')

    !this%dwt_prod10n_gain(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_DWT_PROD10N_GAIN', units='gN/m^2/s', &
    !     avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
    !     ptr_topo=this%dwt_prod10n_gain, default='inactive')

    !this%dwt_prod100n_gain(begt:endt) = spval
    !call hist_addfld1d (fname='TOP_DWT_PROD100N_GAIN', units='gN/m^2/s', &
    !     avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
    !     ptr_topo=this%dwt_prod100n_gain, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of top_nf
    !-----------------------------------------------------------------------
    do t = begt, endt
       this%dwt_prod10n_gain(t)          = 0._r8
       this%dwt_prod100n_gain(t)         = 0._r8
       this%ninputs(t)                   = 0._r8
       this%noutputs(t)                  = 0._r8
    end do
    
  end subroutine top_nf_init
  
  !-----------------------------------------------------------------------
  subroutine top_nf_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(topounit_nitrogen_flux)  :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: t          ! indices
    !-----------------------------------------------------------------------
  
    do t = bounds%begt, bounds%endt
       this%dwt_seedn_to_leaf(t)     = 0._r8
       this%dwt_seedn_to_deadstem(t) = 0._r8
       this%dwt_conv_nflux(t)        = 0._r8
       this%dwt_seedn_to_npool(t)    = 0._r8
       this%dwt_prod10n_gain(t)      = 0._r8
       this%dwt_prod100n_gain(t)     = 0._r8
    end do

  end subroutine top_nf_zerodwt
  
  !------------------------------------------------------------------------
  subroutine top_nf_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_nitrogen_flux) :: this
    !------------------------------------------------------------------------
  
  end subroutine top_nf_clean
  
  subroutine top_ps_init (this, begt, endt)
    !
    ! !ARGUMENTS:
    class(topounit_phosphorus_state) :: this
    integer, intent(in) :: begt,endt
    !
    ! !LOCAL VARIABLES:
    integer :: t
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of top_ps
    !-----------------------------------------------------------------------
    allocate(this%seedp   (begt:endt)) ; this%seedp   (:) = spval
    allocate(this%begpb   (begt:endt)) ; this%begpb   (:) = spval
    allocate(this%endpb   (begt:endt)) ; this%endpb   (:) = spval
    allocate(this%errpb   (begt:endt)) ; this%errpb   (:) = spval
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_ps
    !-----------------------------------------------------------------------
    this%seedp(begt:endt) = spval
    call hist_addfld1d (fname='TOP_SEEDP', units='gP/m^2', &
         avgflag='A', long_name='P pool for seeding new PFTs ', &
         ptr_topo=this%seedp, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of top_ps
    !------------------------------------------------------------------------
    do t = begt, endt
       this%seedp(t) = 0._r8
    end do

  
  end subroutine top_ps_init

  !------------------------------------------------------------------------
  subroutine top_ps_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_phosphorus_state) :: this
    !------------------------------------------------------------------------
  
  end subroutine top_ps_clean

  subroutine top_pf_init (this, begt, endt)
    !
    ! !ARGUMENTS:
    class(topounit_phosphorus_flux) :: this
    integer, intent(in) :: begt,endt
    !
    ! !LOCAL VARIABLES:
    integer :: t
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of top_pf
    !-----------------------------------------------------------------------
    allocate(this%dwt_seedp_to_leaf      (begt:endt))   ; this%dwt_seedp_to_leaf      (:) = spval
    allocate(this%dwt_seedp_to_deadstem  (begt:endt))   ; this%dwt_seedp_to_deadstem  (:) = spval
    allocate(this%dwt_conv_pflux         (begt:endt))   ; this%dwt_conv_pflux         (:) = spval
    allocate(this%dwt_seedp_to_ppool     (begt:endt))   ; this%dwt_seedp_to_ppool     (:) = spval
    allocate(this%dwt_prod10p_gain       (begt:endt))   ; this%dwt_prod10p_gain       (:) = spval
    allocate(this%dwt_prod100p_gain      (begt:endt))   ; this%dwt_prod100p_gain      (:) = spval
    allocate(this%pinputs                (begt:endt))   ; this%pinputs                (:) = spval
    allocate(this%poutputs               (begt:endt))   ; this%poutputs               (:) = spval
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of top_pf
    !-----------------------------------------------------------------------
    this%dwt_seedp_to_leaf(begt:endt) = spval
    call hist_addfld1d (fname='TOP_DWT_SEEDP_TO_LEAF', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to patch-level leaf', &
         ptr_topo=this%dwt_seedp_to_leaf, default='inactive')

    this%dwt_seedp_to_deadstem(begt:endt) = spval
    call hist_addfld1d (fname='TOP_DWT_SEEDP_TO_DEADSTEM', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to patch-level deadstem', &
         ptr_topo=this%dwt_seedp_to_deadstem, default='inactive')

    this%dwt_conv_pflux(begt:endt) = spval
    call hist_addfld1d (fname='TOP_DWT_CONV_PFLUX', units='gP/m^2/s', &
         avgflag='A', &
         long_name='conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
         ptr_topo=this%dwt_conv_pflux)

    this%dwt_seedp_to_ppool(begt:endt) = spval
    call hist_addfld1d (fname='TOP_DWT_SEEDP_TO_PPOOL', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level', &
         ptr_topo=this%dwt_seedp_to_ppool, default='inactive')

    this%dwt_prod10p_gain(begt:endt) = spval
    call hist_addfld1d (fname='TOP_DWT_PROD10P_GAIN', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_topo=this%dwt_prod10p_gain, default='inactive')

    this%dwt_prod100p_gain(begt:endt) = spval
    call hist_addfld1d (fname='TOP_DWT_PROD100P_GAIN', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_topo=this%dwt_prod100p_gain, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of top_pf
    !------------------------------------------------------------------------
    do t = begt, endt
       this%dwt_prod10p_gain(t)          = 0._r8
       this%dwt_prod100p_gain(t)         = 0._r8
       this%pinputs(t)                   = 0._r8
       this%poutputs(t)                  = 0._r8
    end do
  
  end subroutine top_pf_init

  !-----------------------------------------------------------------------
  subroutine top_pf_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(topounit_phosphorus_flux) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: t          ! indices
    !-----------------------------------------------------------------------
    do t = bounds%begt, bounds%endt
       this%dwt_seedp_to_leaf(t)     = 0._r8
       this%dwt_seedp_to_deadstem(t) = 0._r8
       this%dwt_conv_pflux(t)        = 0._r8
       this%dwt_seedp_to_ppool(t)    = 0._r8
       this%dwt_prod10p_gain(t)      = 0._r8
       this%dwt_prod100p_gain(t)     = 0._r8
    end do
  
  end subroutine top_pf_zerodwt
  
  !------------------------------------------------------------------------
  subroutine top_pf_clean(this)
    !
    ! !ARGUMENTS:
    class(topounit_phosphorus_flux) :: this
    !------------------------------------------------------------------------
  
  end subroutine top_pf_clean
  
end module TopounitDataType
