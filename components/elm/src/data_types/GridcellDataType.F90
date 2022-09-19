module GridcellDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gridcell data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_infnan_mod    , only : nan => shr_infnan_nan, assignment(=)
  use elm_varpar        , only : nlevsno, nlevgrnd, nlevlak, nlevurb
  use elm_varcon        , only : spval, ispval
  use elm_varctl        , only : use_fates
  use histFileMod       , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ncdio_pio         , only : file_desc_t, ncd_double
  use decompMod         , only : bounds_type
  use restUtilMod
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  use ColumnType        , only : col_pp
  use LandunitType      , only : lun_pp
  use GridcellType      , only : grc_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_energy_state
    ! temperature variables
    real(r8), pointer :: heat1                (:)   ! initial gridcell total heat content
    real(r8), pointer :: heat2                (:)   ! post land cover change total heat content
    real(r8), pointer :: liquid_water_temp1   (:)   ! initial weighted average liquid water temperature (K)
    real(r8), pointer :: liquid_water_temp2   (:)   ! post land cover change weighted average liquid water temperature (K)

  contains
    procedure, public :: Init    => grc_es_init
    procedure, public :: Clean   => grc_es_clean
  end type gridcell_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_energy_flux
    real(r8), pointer :: eflx_dynbal           (:)   ! dynamic land cover change conversion energy flux (W/m**2)

  contains
    procedure, public :: Init    => grc_ef_init
    procedure, public :: Restart => grc_ef_restart
    procedure, public :: Clean   => grc_ef_clean
  end type gridcell_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_water_state
    ! Note: All units given as kg in this data type imply kg of H2O
    real(r8), pointer :: liq1               (:)   ! initial gridcell total h2o liq content (kg/m2)
    real(r8), pointer :: liq2               (:)   ! post land cover change total liq content (kg/m2)
    real(r8), pointer :: ice1               (:)   ! initial gridcell total h2o ice content (kg/m2)
    real(r8), pointer :: ice2               (:)   ! post land cover change total ice content (kg/m2)
    real(r8), pointer :: tws                (:)   ! total water storage (kg/m2)
    real(r8), pointer :: tws_month_beg      (:)   ! grc total water storage at the beginning of a month
    real(r8), pointer :: tws_month_end      (:)   ! grc total water storage at the end of a month
    real(r8), pointer :: begwb              (:)   ! water mass begining of the time step
    real(r8), pointer :: endwb              (:)   ! water mass end of the time step
    real(r8), pointer :: errh2o             (:)   ! water conservation error (mm H2O)
    real(r8), pointer :: beg_h2ocan         (:)   ! grid-level canopy water at begining of the time step (mm H2O)
    real(r8), pointer :: beg_h2osno         (:)   ! grid-level snow water at begining of the time step (mm H2O)
    real(r8), pointer :: beg_h2osfc         (:)   ! grid-level surface water at begining of the time step (mm H2O)
    real(r8), pointer :: beg_h2osoi_liq     (:)   ! grid-level liquid water at begining of the time step (kg/m2)
    real(r8), pointer :: beg_h2osoi_ice     (:)   ! grid-level ice lens at begining of the time step (kg/m2)
    real(r8), pointer :: end_h2ocan         (:)   ! grid-level canopy water at end of the time step (mm H2O)
    real(r8), pointer :: end_h2osno         (:)   ! grid-level snow water at end of the time step (mm H2O)
    real(r8), pointer :: end_h2osfc         (:)   ! grid-level surface water at end of the time step (mm H2O)
    real(r8), pointer :: end_h2osoi_liq     (:)   ! grid-level liquid water at end of the time step (kg/m2)
    real(r8), pointer :: end_h2osoi_ice     (:)   ! grid-level ice lens at end of the time step (kg/m2)

  contains
    procedure, public :: Init    => grc_ws_init
    procedure, public :: Restart => grc_ws_restart
    procedure, public :: Clean   => grc_ws_clean
  end type gridcell_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_water_flux
    ! Dynamic land cover change
    real(r8), pointer :: qflx_liq_dynbal      (:)   ! liq dynamic land cover change conversion runoff flux
    real(r8), pointer :: qflx_ice_dynbal      (:)   ! ice dynamic land cover change conversion runoff flux

    ! Objects that help convert once-per-year dynamic land cover changes into fluxes
    ! that are dribbled throughout the year
    type(annual_flux_dribbler_type) :: qflx_liq_dynbal_dribbler
    type(annual_flux_dribbler_type) :: qflx_ice_dynbal_dribbler

  contains
    procedure, public :: Init    => grc_wf_init
    procedure, public :: Restart => grc_wf_restart
    procedure, public :: Clean   => grc_wf_clean
  end type gridcell_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_carbon_state
    real(r8), pointer :: seedc                 (:) => null() ! (gC/m2) pool for seeding new PFTs via dynamic landcover
    real(r8), pointer :: tcs_month_beg         (:) => null()  ! grc total carbon storage at the beginning of a month
    real(r8), pointer :: tcs_month_end         (:) => null()  ! grc total carbon storage at the end of a month
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
    procedure, public :: Init    => grc_cs_init
    procedure, public :: Restart => grc_cs_restart
    procedure, public :: Clean   => grc_cs_clean
  end type gridcell_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_carbon_flux
    real(r8), pointer :: dwt_seedc_to_leaf         (:) => null()  ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the gridcell-level
    real(r8), pointer :: dwt_seedc_to_deadstem     (:) => null()  ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the gridcell-level
    real(r8), pointer :: dwt_conv_cflux            (:) => null()  ! (gC/m2/s) dwt_conv_cflux_patch summed to the gridcell-level
    real(r8), pointer :: dwt_conv_cflux_dribbled   (:) => null()  ! (gC/m2/s) dwt_conv_cflux dribbled evenly throughout the year
    real(r8), pointer :: dwt_prod10c_gain          (:) => null()  ! (gC/m2/s) dynamic landcover addition to 10-year wood product pool
    real(r8), pointer :: dwt_prod100c_gain         (:) => null()  ! (gC/m2/s) dynamic landcover addition to 100-year wood product pool
    real(r8), pointer :: hrv_deadstemc_to_prod10c  (:) => null()  ! (gC/m2/s) dead stem harvest to 10-year wood product pool
    real(r8), pointer :: hrv_deadstemc_to_prod100c (:) => null()  ! (gC/m2/s) dead stem harvest to 100-year wood product pool
    real(r8), pointer :: cinputs                   (:) => null()  ! (gC/m2/s) grid-level C inputs
    real(r8), pointer :: coutputs                  (:) => null()  ! (gC/m2/s) grid-level C outputs
    real(r8), pointer :: gpp                       (:) => null()  ! (gC/m2/s) grid-level gross primary production
    real(r8), pointer :: er                        (:) => null()  ! (gC/m2/s) grid-level total ecosystem respiration
    real(r8), pointer :: fire_closs                (:) => null()  ! (gC/m2/s) grid-level total fire C loss
    real(r8), pointer :: prod1_loss                (:) => null()  ! (gC/m2/s) grid-level crop leafc harvested
    real(r8), pointer :: prod10_loss               (:) => null()  ! (gC/m2/s) grid-level 10-year wood C harvested
    real(r8), pointer :: prod100_loss              (:) => null()  ! (gC/m2/s) grid-level 100-year wood C harvested
    real(r8), pointer :: hrv_xsmrpool_to_atm       (:) => null()  ! (gC/m2/s) grid-level excess MR pool harvest mortality
    real(r8), pointer :: som_c_leached             (:) => null()  ! (gC/m2/s) grid-level total SOM C loss from vertical transport
    real(r8), pointer :: somc_yield                (:) => null()  ! (gC/m2/s) grid-level total SOM C loss by erosion

  contains
    procedure, public :: Init    => grc_cf_init
    procedure, public :: ZeroDWT => grc_cf_zerodwt
    procedure, public :: Clean   => grc_cf_clean
  end type gridcell_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_nitrogen_state
    real(r8), pointer :: seedn          (:) => null()   ! (gNm2) nitrogen pool for seeding new PFTs via dynamic landcover
    real(r8), pointer :: begnb          (:) => null()   ! (gNm2) nitrogen mass, beginning of time step
    real(r8), pointer :: endnb          (:) => null()   ! (gNm2) nitrogen mass, end of time step 
    real(r8), pointer :: errnb          (:) => null()   ! (gNm2) nitrogen balance error for the timestep

  contains
    procedure, public :: Init    => grc_ns_init
    procedure, public :: Clean   => grc_ns_clean
  end type gridcell_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_nitrogen_flux
    ! dynamic landcover fluxes
    real(r8), pointer :: dwt_seedn_to_leaf      (:) => null()  ! (gN/m2/s) dwt_seedn_to_leaf_patch summed to the gridcell-level
    real(r8), pointer :: dwt_seedn_to_deadstem  (:) => null()  ! (gN/m2/s) dwt_seedn_to_deadstem_patch summed to the gridcell-level
    real(r8), pointer :: dwt_conv_nflux         (:) => null()  ! (gN/m2/s) dwt_conv_nflux_patch summed to the gridcell-level
    real(r8), pointer :: dwt_seedn_to_npool     (:) => null()  ! (gN/m2/s) seed source to PFT level
    real(r8), pointer :: dwt_prod10n_gain       (:) => null()  ! (gN/m2/s) addition to 10-yr wood product pool
    real(r8), pointer :: dwt_prod100n_gain      (:) => null()  ! (gN/m2/s) addition to 100-yr wood product pool
    real(r8), pointer :: ninputs                (:) => null()  ! (gN/m2/s) grid-level N inputs
    real(r8), pointer :: noutputs               (:) => null()  ! (gN/m2/s) grid-level N outputs
  contains
    procedure, public :: Init    => grc_nf_init
    procedure, public :: ZeroDWT => grc_nf_zerodwt
    procedure, public :: Clean   => grc_nf_clean
  end type gridcell_nitrogen_flux
 
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_phosphorus_state
    real(r8), pointer :: seedp                    (:)     ! (gP/m2) gridcell-level pool for seeding new PFTs via dynamic landcover
    real(r8), pointer :: begpb                    (:)     ! phosphorus mass, beginning of time step (gP/m**2)
    real(r8), pointer :: endpb                    (:)     ! phosphorus mass, end of time step (gP/m**2)
    real(r8), pointer :: errpb                    (:)     ! phosphorus balance error for the timestep (gP/m**2)
  contains
    procedure, public :: Init    => grc_ps_init
    procedure, public :: Clean   => grc_ps_clean
  end type gridcell_phosphorus_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_phosphorus_flux
    real(r8), pointer :: dwt_seedp_to_leaf        (:)  ! (gP/m2/s) dwt_seedn_to_leaf_patch summed to the gridcell-level
    real(r8), pointer :: dwt_seedp_to_deadstem    (:)  ! (gP/m2/s) dwt_seedn_to_deadstem_patch summed to the gridcell-level
    real(r8), pointer :: dwt_conv_pflux           (:)  ! (gP/m2/s) dwt_conv_nflux_patch summed to the gridcell-level
    real(r8), pointer :: dwt_seedp_to_ppool       (:)  ! (gP/m2/s) seed source to PFT-level
    real(r8), pointer :: dwt_prod10p_gain         (:)  ! (gP/m2/s) addition to 10-yr wood product pool
    real(r8), pointer :: dwt_prod100p_gain        (:)  ! (gP/m2/s) addition to 100-yr wood product pool
    real(r8), pointer :: pinputs                  (:)  ! (gP/m2/s) grid-level P inputs
    real(r8), pointer :: poutputs                 (:)  ! (gP/m2/s) grid-level P outputs
  contains
    procedure, public :: Init    => grc_pf_init
    procedure, public :: ZeroDWT => grc_pf_zerodwt
    procedure, public :: Clean   => grc_pf_clean
  end type gridcell_phosphorus_flux
  
  !-----------------------------------------------------------------------
  ! declare the public instances of gridcell-level data types
  !-----------------------------------------------------------------------
  type(gridcell_energy_state)          , public, target :: grc_es     ! gridcell energy state
  type(gridcell_energy_flux)           , public, target :: grc_ef     ! gridcell energy flux
  type(gridcell_water_state)           , public, target :: grc_ws     ! gridcell water state
  type(gridcell_water_flux)            , public, target :: grc_wf     ! gridcell water flux
  type(gridcell_carbon_state)          , public, target :: grc_cs     ! gridcell carbon state
  type(gridcell_carbon_state)          , public, target :: c13_grc_cs ! gridcell carbon state (C13)
  type(gridcell_carbon_state)          , public, target :: c14_grc_cs ! gridcell carbon state (C14)
  type(gridcell_carbon_flux)           , public, target :: grc_cf     ! gridcell carbon flux
  type(gridcell_carbon_flux)           , public, target :: c13_grc_cf ! gridcell carbon flux (C13)
  type(gridcell_carbon_flux)           , public, target :: c14_grc_cf ! gridcell carbon flux (C14)
  type(gridcell_nitrogen_state)        , public, target :: grc_ns     ! gridcell nitrogen state
  type(gridcell_nitrogen_flux)         , public, target :: grc_nf     ! gridcell nitrogen flux
  type(gridcell_phosphorus_state)      , public, target :: grc_ps     ! gridcell phosphorus state
  type(gridcell_phosphorus_flux)       , public, target :: grc_pf     ! gridcell phosphorus flux
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell energy state data structure
  !------------------------------------------------------------------------
  subroutine grc_es_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_state) :: this
    integer, intent(in) :: begg,endg
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_es
    !-----------------------------------------------------------------------
    allocate(this%heat1                (begg:endg))                      ; this%heat1                (:)   = nan
    allocate(this%heat2                (begg:endg))                      ; this%heat2                (:)   = nan
    allocate(this%liquid_water_temp1   (begg:endg))                      ; this%liquid_water_temp1   (:)   = nan
    allocate(this%liquid_water_temp2   (begg:endg))                      ; this%liquid_water_temp2   (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_es
    !-----------------------------------------------------------------------
    this%heat1(begg:endg) = spval
    call hist_addfld1d (fname='GC_HEAT1',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content', &
         ptr_lnd=this%heat1)

    this%heat2(begg:endg) = spval
    call hist_addfld1d (fname='GC_HEAT2',  units='J/m^2',  &
         avgflag='A', long_name='post land cover change total heat content', &
         ptr_lnd=this%heat2, default='inactive')  

    this%liquid_water_temp1(begg:endg) = spval
    call hist_addfld1d (fname='LIQUID_WATER_TEMP1', units='K', &
         avgflag='A', long_name='initial gridcell weighted average liquid water temperature', &
         ptr_lnd=this%liquid_water_temp1, default='inactive')

  end subroutine grc_es_init

  !------------------------------------------------------------------------
  subroutine grc_es_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%heat1)
    deallocate(this%heat2)
    deallocate(this%liquid_water_temp1)
    deallocate(this%liquid_water_temp2)
  end subroutine grc_es_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell energy flux data structure
  !------------------------------------------------------------------------
  subroutine grc_ef_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_flux) :: this
    integer, intent(in) :: begg,endg
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_dynbal              (begg:endg))                  ; this%eflx_dynbal             (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_ef
    !-----------------------------------------------------------------------
    this%eflx_dynbal(begg:endg) = spval 
    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=this%eflx_dynbal)

  end subroutine grc_ef_init

  !------------------------------------------------------------------------
  subroutine grc_ef_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write gridcell energy flux information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(gridcell_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------
  end subroutine grc_ef_restart

  !------------------------------------------------------------------------
  subroutine grc_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%eflx_dynbal)
  end subroutine grc_ef_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell water state data structure
  !------------------------------------------------------------------------
  subroutine grc_ws_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_water_state) :: this
    integer, intent(in) :: begg,endg
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_ws
    !-----------------------------------------------------------------------
    allocate(this%liq1           (begg:endg))       ; this%liq1           (:)   = nan
    allocate(this%liq2           (begg:endg))       ; this%liq2           (:)   = nan
    allocate(this%ice1           (begg:endg))       ; this%ice1           (:)   = nan
    allocate(this%ice2           (begg:endg))       ; this%ice2           (:)   = nan
    allocate(this%tws            (begg:endg))       ; this%tws            (:)   = nan
    allocate(this%tws_month_beg  (begg:endg))       ; this%tws_month_beg  (:)   = nan
    allocate(this%tws_month_end  (begg:endg))       ; this%tws_month_end  (:)   = nan
    allocate(this%begwb          (begg:endg))       ; this%begwb          (:)   = nan
    allocate(this%endwb          (begg:endg))       ; this%endwb          (:)   = nan
    allocate(this%errh2o         (begg:endg))       ; this%errh2o         (:)   = nan
    allocate(this%beg_h2ocan     (begg:endg))       ; this%beg_h2ocan     (:)   = nan
    allocate(this%beg_h2osno     (begg:endg))       ; this%beg_h2osno     (:)   = nan
    allocate(this%beg_h2osfc     (begg:endg))       ; this%beg_h2osfc     (:)   = nan
    allocate(this%beg_h2osoi_liq (begg:endg))       ; this%beg_h2osoi_liq (:)   = nan
    allocate(this%beg_h2osoi_ice (begg:endg))       ; this%beg_h2osoi_ice (:)   = nan
    allocate(this%end_h2ocan     (begg:endg))       ; this%end_h2ocan     (:)   = nan
    allocate(this%end_h2osno     (begg:endg))       ; this%end_h2osno     (:)   = nan
    allocate(this%end_h2osfc     (begg:endg))       ; this%end_h2osfc     (:)   = nan
    allocate(this%end_h2osoi_liq (begg:endg))       ; this%end_h2osoi_liq (:)   = nan
    allocate(this%end_h2osoi_ice (begg:endg))       ; this%end_h2osoi_ice (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_ws
    !-----------------------------------------------------------------------
    this%liq1(begg:endg) = spval
    call hist_addfld1d (fname='GC_LIQ1',  units='mm',  &
         avgflag='A', long_name='initial gridcell total liq content', &
         ptr_lnd=this%liq1)

    this%liq2(begg:endg) = spval
    call hist_addfld1d (fname='GC_LIQ2',  units='mm',  &  
         avgflag='A', long_name='post landuse change gridcell total liq content', &              
         ptr_lnd=this%liq2, default='inactive')     

    this%ice1(begg:endg) = spval
    call hist_addfld1d (fname='GC_ICE1',  units='mm',  &  
         avgflag='A', long_name='initial gridcell total ice content', &              
         ptr_lnd=this%ice1)     

    this%ice2(begg:endg) = spval
    call hist_addfld1d (fname='GC_ICE2',  units='mm',  &  
         avgflag='A', long_name='post land cover change total ice content', &              
         ptr_lnd=this%ice2, default='inactive')

    this%tws(begg:endg) = spval
    call hist_addfld1d (fname='TWS',  units='mm',  &
         avgflag='A', long_name='total water storage', &
         ptr_lnd=this%tws)

    this%tws_month_beg(begg:endg) = spval
    call hist_addfld1d (fname='TWS_MONTH_BEGIN',  units='mm',  &
         avgflag='I', long_name='total water storage at the beginning of a month', &
         ptr_lnd=this%tws_month_beg)

    this%tws_month_end(begg:endg) = spval
    call hist_addfld1d (fname='TWS_MONTH_END',  units='mm',  &
         avgflag='I', long_name='total water storage at the end of a month', &
         ptr_lnd=this%tws_month_end)
  end subroutine grc_ws_init

  !------------------------------------------------------------------------
  subroutine grc_ws_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write gridcell water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(gridcell_water_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------
    call restartvar(ncid=ncid, flag=flag, varname='TWS_MONTH_BEGIN', xtype=ncd_double,  &
         dim1name='gridcell', &
         long_name='surface watertotal water storage at the beginning of a month', units='mm', &
          interpinic_flag='interp', readvar=readvar, data=this%tws_month_beg)
  
  end subroutine grc_ws_restart
  
  !------------------------------------------------------------------------
  subroutine grc_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_water_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%liq1)
    deallocate(this%liq2)
    deallocate(this%ice1)
    deallocate(this%ice2)
    deallocate(this%tws )
  end subroutine grc_ws_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell water flux data structure
  !------------------------------------------------------------------------
  subroutine grc_wf_init(this, begg, endg, bounds)
    !
    ! !ARGUMENTS:
    class(gridcell_water_flux) :: this
    integer, intent(in) :: begg,endg
    type(bounds_type)    , intent(in)  :: bounds
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_wf
    !-----------------------------------------------------------------------
    allocate(this%qflx_liq_dynbal      (begg:endg))              ; this%qflx_liq_dynbal      (:)   = nan
    allocate(this%qflx_ice_dynbal      (begg:endg))              ; this%qflx_ice_dynbal      (:)   = nan

    this%qflx_liq_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'qflx_liq_dynbal', &
         units = 'mm H2O')

    this%qflx_ice_dynbal_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &
         name = 'qflx_ice_dynbal', &
         units = 'mm H2O')

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_wf
    !-----------------------------------------------------------------------
    this%qflx_liq_dynbal(begg:endg) = spval
    call hist_addfld1d (fname='QFLX_LIQ_DYNBAL',  units='mm/s',  &  
         avgflag='A', long_name='liq dynamic land cover change conversion runoff flux', &              
         ptr_lnd=this%qflx_liq_dynbal)     

    this%qflx_ice_dynbal(begg:endg) = spval
    call hist_addfld1d (fname='QFLX_ICE_DYNBAL',  units='mm/s',  &
         avgflag='A', long_name='ice dynamic land cover change conversion runoff flux', &                                   
         ptr_lnd=this%qflx_ice_dynbal)
  
  end subroutine grc_wf_init

  !------------------------------------------------------------------------
  subroutine grc_wf_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write gridcell water flux information to/from restart file.
    !
    ! !ARGUMENTS:
    class(gridcell_water_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call this%qflx_liq_dynbal_dribbler%Restart(bounds, ncid, flag)
    call this%qflx_ice_dynbal_dribbler%Restart(bounds, ncid, flag)

  end subroutine grc_wf_restart

  !------------------------------------------------------------------------
  subroutine grc_wf_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_water_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%qflx_liq_dynbal)
    deallocate(this%qflx_ice_dynbal)

  end subroutine grc_wf_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell carbon state data structure
  !------------------------------------------------------------------------
  subroutine grc_cs_init(this, begg, endg,carbon_type)
    !
    ! !ARGUMENTS:
    class(gridcell_carbon_state) :: this
    integer, intent(in) :: begg,endg
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']    

    !
    ! !LOCAL VARIABLES:
    integer :: g
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_cs
    !-----------------------------------------------------------------------
    allocate(this%seedc                 (begg:endg));     this%seedc                 (:) = nan
    allocate(this%tcs_month_beg         (begg:endg));     this%tcs_month_beg         (:) = nan
    allocate(this%tcs_month_end         (begg:endg));     this%tcs_month_end         (:) = nan
    allocate(this%begcb                 (begg:endg));     this%begcb                 (:) = nan
    allocate(this%endcb                 (begg:endg));     this%endcb                 (:) = nan
    allocate(this%errcb                 (begg:endg));     this%errcb                 (:) = nan

    allocate(this%beg_totc              (begg:endg));     this%beg_totc              (:) = nan
    allocate(this%beg_totpftc           (begg:endg));     this%beg_totpftc           (:) = nan
    allocate(this%beg_cwdc              (begg:endg));     this%beg_cwdc              (:) = nan
    allocate(this%beg_totsomc           (begg:endg));     this%beg_totsomc           (:) = nan
    allocate(this%beg_totlitc           (begg:endg));     this%beg_totlitc           (:) = nan
    allocate(this%beg_totprodc          (begg:endg));     this%beg_totprodc          (:) = nan
    allocate(this%beg_ctrunc            (begg:endg));     this%beg_ctrunc            (:) = nan
    allocate(this%beg_cropseedc_deficit (begg:endg));     this%beg_cropseedc_deficit (:) = nan

    allocate(this%end_totc              (begg:endg));     this%end_totc              (:) = nan
    allocate(this%end_totpftc           (begg:endg));     this%end_totpftc           (:) = nan
    allocate(this%end_cwdc              (begg:endg));     this%end_cwdc              (:) = nan
    allocate(this%end_totsomc           (begg:endg));     this%end_totsomc           (:) = nan
    allocate(this%end_totlitc           (begg:endg));     this%end_totlitc           (:) = nan
    allocate(this%end_totprodc          (begg:endg));     this%end_totprodc          (:) = nan
    allocate(this%end_ctrunc            (begg:endg));     this%end_ctrunc            (:) = nan
    allocate(this%end_cropseedc_deficit (begg:endg));     this%end_cropseedc_deficit (:) = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_cs
    !-----------------------------------------------------------------------
    if (.not. use_fates) then
       if (carbon_type == 'c12') then
          this%seedc(begg:endg) = spval
          call hist_addfld1d (fname='SEEDC_GRC', units='gC/m^2', &
               avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
               ptr_gcell=this%seedc)
       end if 
       if (carbon_type == 'c13') then
          this%seedc(begg:endg) = spval
          call hist_addfld1d (fname='C13_SEEDC_GRC', units='gC/m^2', &
               avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
               ptr_gcell=this%seedc)
       end if 
       if (carbon_type == 'c14') then
          this%seedc(begg:endg) = spval
          call hist_addfld1d (fname='C14_SEEDC_GRC', units='gC/m^2', &
               avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
               ptr_gcell=this%seedc)
       end if 
    end if
    
    this%tcs_month_beg(begg:endg) = spval
    call hist_addfld1d (fname='TCS_MONTH_BEGIN',  units='mm',  &
         avgflag='I', long_name='total carbon storage at the beginning of a month', &
         ptr_lnd=this%tcs_month_beg)

    this%tcs_month_end(begg:endg) = spval
    call hist_addfld1d (fname='TCS_MONTH_END',  units='mm',  &
         avgflag='I', long_name='total carbon storage at the end of a month', &
         ptr_lnd=this%tcs_month_end)

    call hist_addfld1d (fname='CMASS_BALANCE_ERROR',  units='gC/m^2',  &
         avgflag='A', long_name='Gridcell carbon mass balance error', &
         ptr_lnd=this%errcb)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of grc_cs
    !-----------------------------------------------------------------------
    do g = begg, endg
       this%seedc(g) = 0._r8
    end do
    
  end subroutine grc_cs_init

  !------------------------------------------------------------------------
  subroutine grc_cs_restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write gridcell water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(gridcell_carbon_state) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid
    character(len=*) , intent(in)    :: flag
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------
    call restartvar(ncid=ncid, flag=flag, varname='TCS_MONTH_BEGIN', xtype=ncd_double,  &
         dim1name='gridcell', &
         long_name='surface watertotal carbon storage at the beginning of a month', units='mm', &
          interpinic_flag='interp', readvar=readvar, data=this%tcs_month_beg)

  end subroutine grc_cs_restart

  !------------------------------------------------------------------------
  subroutine grc_cs_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_carbon_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%seedc)
    deallocate(this%begcb)
    deallocate(this%endcb)
    deallocate(this%errcb)

  end subroutine grc_cs_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell carbon flux data structure
  !------------------------------------------------------------------------
  subroutine grc_cf_init(this, begg, endg, carbon_type)
    !
    ! !ARGUMENTS:
    class(gridcell_carbon_flux) :: this
    integer, intent(in) :: begg,endg
    character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer :: g
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_cf
    !-----------------------------------------------------------------------
    allocate(this%dwt_seedc_to_leaf            (begg:endg)) ; this%dwt_seedc_to_leaf         (:) = nan
    allocate(this%dwt_seedc_to_deadstem        (begg:endg)) ; this%dwt_seedc_to_deadstem     (:) = nan
    allocate(this%dwt_conv_cflux               (begg:endg)) ; this%dwt_conv_cflux            (:) = nan
    allocate(this%dwt_conv_cflux_dribbled      (begg:endg)) ; this%dwt_conv_cflux_dribbled   (:) = nan
    allocate(this%dwt_prod10c_gain             (begg:endg)) ; this%dwt_prod10c_gain          (:) = nan
    allocate(this%dwt_prod100c_gain            (begg:endg)) ; this%dwt_prod100c_gain         (:) = nan
    allocate(this%hrv_deadstemc_to_prod10c     (begg:endg)) ; this%hrv_deadstemc_to_prod10c  (:) = nan
    allocate(this%hrv_deadstemc_to_prod100c    (begg:endg)) ; this%hrv_deadstemc_to_prod100c (:) = nan
    allocate(this%cinputs                      (begg:endg)) ; this%cinputs                   (:) = nan
    allocate(this%coutputs                     (begg:endg)) ; this%coutputs                  (:) = nan
    allocate(this%gpp                          (begg:endg)) ; this%gpp                       (:) = nan
    allocate(this%er                           (begg:endg)) ; this%er                        (:) = nan
    allocate(this%fire_closs                   (begg:endg)) ; this%fire_closs                (:) = nan
    allocate(this%prod1_loss                   (begg:endg)) ; this%prod1_loss                (:) = nan
    allocate(this%prod10_loss                  (begg:endg)) ; this%prod10_loss               (:) = nan
    allocate(this%prod100_loss                 (begg:endg)) ; this%prod100_loss              (:) = nan
    allocate(this%hrv_xsmrpool_to_atm          (begg:endg)) ; this%hrv_xsmrpool_to_atm       (:) = nan
    allocate(this%som_c_leached                (begg:endg)) ; this%som_c_leached             (:) = nan
    allocate(this%somc_yield                   (begg:endg)) ; this%somc_yield                (:) = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_cf
    !-----------------------------------------------------------------------
    ! no history fields or cold-start initialization for gridcell carbon flux
    ! if using fates
    if (use_fates) then
       return
    end if

    if (carbon_type == 'c12') then

       this%dwt_seedc_to_leaf(begg:endg) = spval
       call hist_addfld1d (fname='DWT_SEEDC_TO_LEAF_GRC', units='gC/m^2/s', &
            avgflag='A', long_name='seed source to patch-level leaf', &
            ptr_gcell=this%dwt_seedc_to_leaf, default='inactive')

       this%dwt_seedc_to_deadstem(begg:endg) = spval
       call hist_addfld1d (fname='DWT_SEEDC_TO_DEADSTEM_GRC', units='gC/m^2/s', &
            avgflag='A', long_name='seed source to patch-level deadstem', &
            ptr_gcell=this%dwt_seedc_to_deadstem, default='inactive')

       this%dwt_conv_cflux(begg:endg) = spval
       call hist_addfld1d (fname='DWT_CONV_CFLUX_GRC', units='gC/m^2/s', &
            avgflag='A', &
            long_name='conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
            ptr_gcell=this%dwt_conv_cflux)

       this%dwt_conv_cflux_dribbled(begg:endg) = spval
       call hist_addfld1d (fname='DWT_CONV_CFLUX_DRIBBLED', units='gC/m^2/s', &
            avgflag='A', &
            long_name='conversion C flux (immediate loss to atm), dribbled throughout the year', &
            ptr_gcell=this%dwt_conv_cflux_dribbled)

       this%dwt_prod10c_gain(begg:endg) = spval
       call hist_addfld1d (fname='DWT_PROD10C_GAIN_GRC', units='gC/m^2/s', &
            avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
            ptr_col=this%dwt_prod10c_gain, default='inactive')

       this%dwt_prod100c_gain(begg:endg) = spval
       call hist_addfld1d (fname='DWT_PROD100C_GAIN_GRC', units='gC/m^2/s', &
            avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
            ptr_col=this%dwt_prod100c_gain, default='inactive')

       this%hrv_deadstemc_to_prod10c(begg:endg) = spval
       call hist_addfld1d (fname='HRV_DEADSTEM_TO_PROD10C_GRC', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem harvest to 10-yr wood product pool', &
            ptr_col=this%hrv_deadstemc_to_prod10c, default='inactive')

       this%hrv_deadstemc_to_prod100c(begg:endg) = spval
       call hist_addfld1d (fname='HRV_DEADSTEM_TO_PROD100C_GRC', units='gC/m^2/s', &
            avgflag='A', long_name='dead stem harvest to 100-yr wood product pool', &
            ptr_col=this%hrv_deadstemc_to_prod100c, default='inactive')
    end if

    if ( carbon_type == 'c13' ) then

       this%dwt_seedc_to_leaf(begg:endg) = spval
       call hist_addfld1d (fname='C13_DWT_SEEDC_TO_LEAF_GRC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 seed source to patch-level leaf', &
            ptr_gcell=this%dwt_seedc_to_leaf, default='inactive')

       this%dwt_seedc_to_deadstem(begg:endg) = spval
       call hist_addfld1d (fname='C13_DWT_SEEDC_TO_DEADSTEM_GRC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 seed source to patch-level deadstem', &
            ptr_gcell=this%dwt_seedc_to_deadstem, default='inactive')

       this%dwt_conv_cflux(begg:endg) = spval
       call hist_addfld1d (fname='C13_DWT_CONV_CFLUX_GRC', units='gC13/m^2/s', &
            avgflag='A', &
            long_name='C13 conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
            ptr_gcell=this%dwt_conv_cflux)

       this%dwt_conv_cflux_dribbled(begg:endg) = spval
       call hist_addfld1d (fname='C13_DWT_CONV_CFLUX_DRIBBLED', units='gC13/m^2/s', &
            avgflag='A', &
            long_name='C13 conversion C flux (immediate loss to atm), dribbled throughout the year', &
            ptr_gcell=this%dwt_conv_cflux_dribbled)

       this%dwt_seedc_to_deadstem(begg:endg) = spval
       call hist_addfld1d (fname='C13_DWT_SEEDC_TO_DEADSTEM', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 seed source to patch-level deadstem', &
            ptr_gcell=this%dwt_seedc_to_deadstem, default='inactive')

       this%dwt_prod10c_gain(begg:endg) = spval
       call hist_addfld1d (fname='C13_DWT_PROD10C_GAIN_GRC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 landcover change-driven addition to 10-yr wood product pool', &
            ptr_col=this%dwt_prod10c_gain, default='inactive')

       this%dwt_prod100c_gain(begg:endg) = spval
       call hist_addfld1d (fname='C13_DWT_PROD100C_GAIN_GRC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 landcover change-driven addition to 100-yr wood product pool', &
            ptr_col=this%dwt_prod100c_gain, default='inactive')

       this%hrv_deadstemc_to_prod10c(begg:endg) = spval
       call hist_addfld1d (fname='C13_HRV_DEADSTEM_TO_PROD10C_GRC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem harvest to 10-yr wood product pool', &
            ptr_col=this%hrv_deadstemc_to_prod10c, default='inactive')

       this%hrv_deadstemc_to_prod100c(begg:endg) = spval
       call hist_addfld1d (fname='C13_HRV_DEADSTEM_TO_PROD100C_GRC', units='gC13/m^2/s', &
            avgflag='A', long_name='C13 dead stem harvest to 100-yr wood product pool', &
            ptr_col=this%hrv_deadstemc_to_prod100c, default='inactive')
    endif

    if (carbon_type == 'c14') then

       this%dwt_seedc_to_leaf(begg:endg) = spval
       call hist_addfld1d (fname='C14_DWT_SEEDC_TO_LEAF_GRC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 seed source to patch-level leaf', &
            ptr_gcell=this%dwt_seedc_to_leaf, default='inactive')

       this%dwt_seedc_to_deadstem(begg:endg) = spval
       call hist_addfld1d (fname='C14_DWT_SEEDC_TO_DEADSTEM_GRC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 seed source to patch-level deadstem', &
            ptr_gcell=this%dwt_seedc_to_deadstem, default='inactive')

       this%dwt_conv_cflux(begg:endg) = spval
       call hist_addfld1d (fname='C14_DWT_CONV_CFLUX_GRC', units='gC14/m^2/s', &
            avgflag='A', &
            long_name='C14 conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
            ptr_gcell=this%dwt_conv_cflux)

       this%dwt_conv_cflux_dribbled(begg:endg) = spval
       call hist_addfld1d (fname='C14_DWT_CONV_CFLUX_DRIBBLED', units='gC14/m^2/s', &
            avgflag='A', &
            long_name='C14 conversion C flux (immediate loss to atm), dribbled throughout the year', &
            ptr_gcell=this%dwt_conv_cflux_dribbled)

       this%dwt_seedc_to_deadstem(begg:endg) = spval
       call hist_addfld1d (fname='C14_DWT_SEEDC_TO_DEADSTEM', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 seed source to patch-level deadstem', &
            ptr_gcell=this%dwt_seedc_to_deadstem, default='inactive')

       this%dwt_prod10c_gain(begg:endg) = spval
       call hist_addfld1d (fname='C14_DWT_PROD10C_GAIN_GRC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 landcover change-driven addition to 10-yr wood product pool', &
            ptr_col=this%dwt_prod10c_gain, default='inactive')

       this%dwt_prod100c_gain(begg:endg) = spval
       call hist_addfld1d (fname='C14_DWT_PROD100C_GAIN_GRC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 landcover change-driven addition to 100-yr wood product pool', &
            ptr_col=this%dwt_prod100c_gain, default='inactive')

       this%hrv_deadstemc_to_prod10c(begg:endg) = spval
       call hist_addfld1d (fname='C14_HRV_DEADSTEM_TO_PROD10C_GRC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem harvest to 10-yr wood product pool', &
            ptr_col=this%hrv_deadstemc_to_prod10c, default='inactive')

       this%hrv_deadstemc_to_prod100c(begg:endg) = spval
       call hist_addfld1d (fname='C14_HRV_DEADSTEM_TO_PROD100C_GRC', units='gC14/m^2/s', &
            avgflag='A', long_name='C14 dead stem harvest to 100-yr wood product pool', &
            ptr_col=this%hrv_deadstemc_to_prod100c, default='inactive')
    endif
    
    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of grc_cf
    !-----------------------------------------------------------------------
    do g = begg, endg
       this%dwt_prod10c_gain(g)          = 0._r8
       this%dwt_prod100c_gain(g)         = 0._r8
       this%hrv_deadstemc_to_prod10c(g)  = 0._r8
       this%hrv_deadstemc_to_prod100c(g) = 0._r8
       this%cinputs(g)                   = 0._r8
       this%coutputs(g)                  = 0._r8
    end do
    
  end subroutine grc_cf_init
  
  !-----------------------------------------------------------------------
  subroutine grc_cf_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(gridcell_carbon_flux)      :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: g          ! indices
    !-----------------------------------------------------------------------

    ! set column-level conversion and product pool fluxes
    ! to 0 at the beginning of every timestep

    do g = bounds%begg, bounds%endg
       this%dwt_seedc_to_leaf(g)         = 0._r8
       this%dwt_seedc_to_deadstem(g)     = 0._r8
       this%dwt_conv_cflux(g)            = 0._r8
       this%dwt_prod10c_gain(g)          = 0._r8
       this%dwt_prod100c_gain(g)         = 0._r8
       this%hrv_deadstemc_to_prod10c(g)  = 0._r8
       this%hrv_deadstemc_to_prod100c(g) = 0._r8
    end do

  end subroutine grc_cf_zerodwt

  !------------------------------------------------------------------------
  subroutine grc_cf_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_carbon_flux) :: this
    !------------------------------------------------------------------------

  end subroutine grc_cf_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine grc_ns_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_nitrogen_state) :: this
    integer, intent(in) :: begg,endg
    !
    ! !LOCAL VARIABLES:
    integer :: g
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_ns
    !-----------------------------------------------------------------------
    allocate(this%seedn   (begg:endg));     this%seedn   (:) = nan
    allocate(this%begnb   (begg:endg));     this%begnb   (:) = nan
    allocate(this%endnb   (begg:endg));     this%endnb   (:) = nan
    allocate(this%errnb   (begg:endg));     this%errnb   (:) = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_cf
    !-----------------------------------------------------------------------
    this%seedn(begg:endg) = spval
    call hist_addfld1d (fname='SEEDN_GRC', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs ', &
         ptr_gcell=this%seedn, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of grc_cf
    !-----------------------------------------------------------------------
    do g = begg, endg
       this%seedn(g) = 0._r8
    end do


  end subroutine grc_ns_init
  
  !------------------------------------------------------------------------
  subroutine grc_ns_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_nitrogen_state) :: this
    !------------------------------------------------------------------------

  end subroutine grc_ns_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine grc_nf_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_nitrogen_flux) :: this
    integer, intent(in) :: begg,endg
    !
    ! !LOCAL VARIABLES:
    integer :: g
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of grc_nf
    !-----------------------------------------------------------------------
    allocate(this%dwt_seedn_to_leaf     (begg:endg)) ; this%dwt_seedn_to_leaf     (:) = nan
    allocate(this%dwt_seedn_to_deadstem (begg:endg)) ; this%dwt_seedn_to_deadstem (:) = nan
    allocate(this%dwt_conv_nflux        (begg:endg)) ; this%dwt_conv_nflux        (:) = nan
    allocate(this%dwt_seedn_to_npool    (begg:endg)) ; this%dwt_seedn_to_npool    (:) = nan
    allocate(this%dwt_prod10n_gain      (begg:endg)) ; this%dwt_prod10n_gain      (:) = nan
    allocate(this%dwt_prod100n_gain     (begg:endg)) ; this%dwt_prod100n_gain     (:) = nan
    allocate(this%ninputs               (begg:endg)) ; this%ninputs               (:) = nan
    allocate(this%noutputs              (begg:endg)) ; this%noutputs              (:) = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_nf
    !-----------------------------------------------------------------------
    this%dwt_seedn_to_leaf(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF_GRC', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to patch-level leaf', &
         ptr_gcell=this%dwt_seedn_to_leaf, default='inactive')

    this%dwt_seedn_to_deadstem(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM_GRC', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to patch-level deadstem', &
         ptr_gcell=this%dwt_seedn_to_deadstem, default='inactive')

    this%dwt_conv_nflux(begg:endg) = spval
    call hist_addfld1d (fname='DWT_CONV_NFLUX_GRC', units='gN/m^2/s', &
         avgflag='A', &
         long_name='conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
         ptr_gcell=this%dwt_conv_nflux)

    this%dwt_seedn_to_npool(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_NPOOL_GRC', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level npool', &
         ptr_gcell=this%dwt_seedn_to_npool, default='inactive')

    this%dwt_prod10n_gain(begg:endg) = spval
    call hist_addfld1d (fname='DWT_PROD10N_GAIN_GRC', units='gN/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_gcell=this%dwt_prod10n_gain, default='inactive')

    this%dwt_prod100n_gain(begg:endg) = spval
    call hist_addfld1d (fname='DWT_PROD100N_GAIN_GRC', units='gN/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_gcell=this%dwt_prod100n_gain, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of grc_nf
    !-----------------------------------------------------------------------
    do g = begg, endg
       this%dwt_prod10n_gain(g)          = 0._r8
       this%dwt_prod100n_gain(g)         = 0._r8
       this%ninputs(g)                   = 0._r8
       this%noutputs(g)                  = 0._r8
    end do
    
  end subroutine grc_nf_init
  
  !-----------------------------------------------------------------------
  subroutine grc_nf_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(gridcell_nitrogen_flux)  :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: g          ! indices
    !-----------------------------------------------------------------------
  
    do g = bounds%begg, bounds%endg
       this%dwt_seedn_to_leaf(g)     = 0._r8
       this%dwt_seedn_to_deadstem(g) = 0._r8
       this%dwt_conv_nflux(g)        = 0._r8
       this%dwt_seedn_to_npool(g)    = 0._r8
       this%dwt_prod10n_gain(g)      = 0._r8
       this%dwt_prod100n_gain(g)     = 0._r8
    end do

  end subroutine grc_nf_zerodwt
  
  !------------------------------------------------------------------------
  subroutine grc_nf_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_nitrogen_flux) :: this
    !------------------------------------------------------------------------
  
  end subroutine grc_nf_clean
  
  subroutine grc_ps_init (this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_phosphorus_state) :: this
    integer, intent(in) :: begg,endg
    !
    ! !LOCAL VARIABLES:
    integer :: g
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of grc_ps
    !-----------------------------------------------------------------------
    allocate(this%seedp   (begg:endg)) ; this%seedp   (:) = nan
    allocate(this%begpb   (begg:endg)) ; this%begpb   (:) = nan
    allocate(this%endpb   (begg:endg)) ; this%endpb   (:) = nan
    allocate(this%errpb   (begg:endg)) ; this%errpb   (:) = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_ps
    !-----------------------------------------------------------------------
    this%seedp(begg:endg) = spval
    call hist_addfld1d (fname='SEEDP_GRC', units='gP/m^2', &
         avgflag='A', long_name='P pool for seeding new PFTs ', &
         ptr_gcell=this%seedp, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of grc_ps
    !------------------------------------------------------------------------
    do g = begg, endg
       this%seedp(g) = 0._r8
    end do

  
  end subroutine grc_ps_init

  !------------------------------------------------------------------------
  subroutine grc_ps_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_phosphorus_state) :: this
    !------------------------------------------------------------------------
  
  end subroutine grc_ps_clean

  subroutine grc_pf_init (this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_phosphorus_flux) :: this
    integer, intent(in) :: begg,endg
    !
    ! !LOCAL VARIABLES:
    integer :: g
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of grc_pf
    !-----------------------------------------------------------------------
    allocate(this%dwt_seedp_to_leaf      (begg:endg))   ; this%dwt_seedp_to_leaf      (:) = nan
    allocate(this%dwt_seedp_to_deadstem  (begg:endg))   ; this%dwt_seedp_to_deadstem  (:) = nan
    allocate(this%dwt_conv_pflux         (begg:endg))   ; this%dwt_conv_pflux         (:) = nan
    allocate(this%dwt_seedp_to_ppool     (begg:endg))   ; this%dwt_seedp_to_ppool     (:) = nan
    allocate(this%dwt_prod10p_gain       (begg:endg))   ; this%dwt_prod10p_gain       (:) = nan
    allocate(this%dwt_prod100p_gain      (begg:endg))   ; this%dwt_prod100p_gain      (:) = nan
    allocate(this%pinputs                (begg:endg))   ; this%pinputs                (:) = nan
    allocate(this%poutputs               (begg:endg))   ; this%poutputs               (:) = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_pf
    !-----------------------------------------------------------------------
    this%dwt_seedp_to_leaf(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_LEAF_GRC', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to patch-level leaf', &
         ptr_gcell=this%dwt_seedp_to_leaf, default='inactive')

    this%dwt_seedp_to_deadstem(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_DEADSTEM_GRC', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to patch-level deadstem', &
         ptr_gcell=this%dwt_seedp_to_deadstem, default='inactive')

    this%dwt_conv_pflux(begg:endg) = spval
    call hist_addfld1d (fname='DWT_CONV_PFLUX_GRC', units='gP/m^2/s', &
         avgflag='A', &
         long_name='conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)', &
         ptr_gcell=this%dwt_conv_pflux)

    this%dwt_seedp_to_ppool(begg:endg) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_PPOOL_GRC', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level', &
         ptr_gcell=this%dwt_seedp_to_ppool, default='inactive')

    this%dwt_prod10p_gain(begg:endg) = spval
    call hist_addfld1d (fname='DWT_PROD10P_GAIN_GRC', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_gcell=this%dwt_prod10p_gain, default='inactive')

    this%dwt_prod100p_gain(begg:endg) = spval
    call hist_addfld1d (fname='DWT_PROD100P_GAIN_GRC', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_gcell=this%dwt_prod100p_gain, default='inactive')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of grc_pf
    !------------------------------------------------------------------------
    do g = begg, endg
       this%dwt_prod10p_gain(g)          = 0._r8
       this%dwt_prod100p_gain(g)         = 0._r8
       this%pinputs(g)                   = 0._r8
       this%poutputs(g)                  = 0._r8
    end do
  
  end subroutine grc_pf_init

  !-----------------------------------------------------------------------
  subroutine grc_pf_zerodwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(gridcell_phosphorus_flux) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: g          ! indices
    !-----------------------------------------------------------------------
    do g = bounds%begg, bounds%endg
       this%dwt_seedp_to_leaf(g)     = 0._r8
       this%dwt_seedp_to_deadstem(g) = 0._r8
       this%dwt_conv_pflux(g)        = 0._r8
       this%dwt_seedp_to_ppool(g)    = 0._r8
       this%dwt_prod10p_gain(g)      = 0._r8
       this%dwt_prod100p_gain(g)     = 0._r8
    end do
  
  end subroutine grc_pf_zerodwt
  
  !------------------------------------------------------------------------
  subroutine grc_pf_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_phosphorus_flux) :: this
    !------------------------------------------------------------------------
  
  end subroutine grc_pf_clean
  
end module GridcellDataType

  
    
