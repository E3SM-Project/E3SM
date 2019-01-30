module ColumnDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_const_mod  , only : SHR_CONST_TKFRZ
  use clm_varpar     , only : nlevsoi, nlevsno, nlevgrnd, nlevlak, nlevurb
  use clm_varcon     , only : spval, ispval, zlnd, snw_rds_min, denice, denh2o, tfrz, pondmx
  use clm_varcon     , only : watmin, bdsno
  use clm_varctl     , only : use_fates_planthydro, create_glacier_mec_landunit
  use clm_varctl     , only : bound_h2osoi, use_cn, iulog
  use clm_time_manager, only : is_first_step
  use landunit_varcon, only : istice, istwet, istsoil, istdlak, istcrop, istice_mec  
  use column_varcon  , only : icol_road_perv, icol_road_imperv, icol_roof, icol_sunwall, icol_shadewall
  use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use spmdMod         , only : masterproc
  use restUtilMod
  use ColumnType     , only : col_pp
  use LandunitType   , only : lun_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_state
    ! temperature variables
    real(r8), pointer :: t_soisno      (:,:) => null() ! soil temperature (K)  (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: t_ssbef       (:,:) => null() ! soil/snow temperature before update (K) (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: t_h2osfc      (:)   => null() ! surface water temperature (K)
    real(r8), pointer :: t_h2osfc_bef  (:)   => null() ! surface water temperature at start of time step (K)
    real(r8), pointer :: t_soi10cm     (:)   => null() ! soil temperature in top 10cm of soil (K)
    real(r8), pointer :: t_soi17cm     (:)   => null() ! soil temperature in top 17cm of soil (K)
    real(r8), pointer :: t_grnd        (:)   => null() ! ground temperature (K)
    real(r8), pointer :: t_lake        (:,:) => null() ! lake temperature (K)  (1:nlevlak)          
    real(r8), pointer :: t_grnd_r      (:)   => null() ! rural ground temperature (K)
    real(r8), pointer :: t_grnd_u      (:)   => null() ! urban ground temperature (K) (needed by Hydrology2Mod)
    real(r8), pointer :: snot_top      (:)   => null() ! temperature of top snow layer (K)
    real(r8), pointer :: dTdz_top      (:)   => null() ! temperature gradient in top layer  (K m-1)
    real(r8), pointer :: thv           (:)   => null() ! virtual potential temperature (K)
    ! heat content variables (diagnostic)
    real(r8), pointer :: hc_soi        (:)   => null() ! soil heat content (MJ/m2)
    real(r8), pointer :: hc_soisno     (:)   => null() ! soil plus snow heat content (MJ/m2)
    ! other quantities related to energy state
    real(r8), pointer :: emg           (:)   => null() ! ground emissivity (unitless)
    real(r8), pointer :: fact          (:,:) => null() ! factors used in computing tridiagonal matrix
    real(r8), pointer :: c_h2osfc      (:)   => null() ! heat capacity of surface water (J/K)
    ! For coupling with pflotran TH
    real(r8), pointer :: t_nearsurf    (:)   => null() ! near-surface air temperature averaged over bare-veg (K)
    
  contains
    procedure, public :: Init    => col_es_init
    procedure, public :: Restart => col_es_restart
    procedure, public :: Clean   => col_es_clean
  end type column_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_state
    ! Note: All units given as kg in this data type imply kg of H2O
    ! Primary water and ice state variables for soil/snow column
    real(r8), pointer :: h2osoi_liq         (:,:) => null() ! liquid water (-nlevsno+1:nlevgrnd) (kg/m2)     
    real(r8), pointer :: h2osoi_ice         (:,:) => null() ! ice lens (-nlevsno+1:nlevgrnd) (kg/m2)    
    real(r8), pointer :: h2osoi_vol         (:,:) => null() ! volumetric soil water (0<=h2osoi_vol<=watsat) (1:nlevgrnd) (m3/m3)  
    real(r8), pointer :: h2osfc             (:)   => null() ! surface water (kg/m2)
    real(r8), pointer :: h2ocan             (:)   => null() ! canopy water integrated to column (kg/m2)
    real(r8), pointer :: total_plant_stored_h2o(:)=> null() ! total water in plants (used??)
    ! Derived water and ice state variables for soil/snow column, depth varying
    real(r8), pointer :: h2osoi_liqvol      (:,:) => null() ! volumetric liquid water content (-nlevsno+1:nlevgrnd) (m3/m3)
    real(r8), pointer :: h2osoi_icevol      (:,:) => null() ! volumetric ice content (-nlevsno+1:nlevgrnd) (m3/m3)     
    real(r8), pointer :: h2osoi_liq_old     (:,:) => null() ! values from previous time step
    real(r8), pointer :: h2osoi_ice_old     (:,:) => null() ! values from previous time step
    real(r8), pointer :: bw                 (:,:) => null() ! partial density of water in the snow pack (ice + liquid) [kg/m3] 
    real(r8), pointer :: smp_l              (:,:) => null() ! liquid phase soil matric potential (-nlevsno+1:nlevgrnd) (mm h2o)
    real(r8), pointer :: soilp              (:,:) => null() ! soil pressure (1:nlevgrnd) (Pa)
    real(r8), pointer :: swe_old            (:,:) => null() ! initial snow water content (-nlevsno+1:0) (kg/m2)
    real(r8), pointer :: snw_rds            (:,:) => null() ! col snow grain radius (-nlevsno+1:0) (m^-6, or microns)
    real(r8), pointer :: air_vol            (:,:) => null() ! air filled porosity (m3/m3)    
    ! Derived water, ice, and snow variables, column aggregate
    real(r8), pointer :: qg_snow            (:)   => null() ! specific humidity over snow (kg H2O/kg moist air)
    real(r8), pointer :: qg_soil            (:)   => null() ! specific humidity over soil (kg H2O/kg moist air)
    real(r8), pointer :: qg_h2osfc          (:)   => null() ! specific humidity over surface water (kg H2O/kg moist air)
    real(r8), pointer :: qg                 (:)   => null() ! average surface specific humidity (kg H2O/kg moist air)
    real(r8), pointer :: dqgdT              (:)   => null() ! rate of change in specific humidity with temperature (kg H2O/kg moist air/K)
    real(r8), pointer :: h2osoi_liqice_10cm (:)   => null() ! liquid water + ice in top 10cm of soil (kg/m2)
    real(r8), pointer :: h2osno             (:)   => null() ! snow mass (kg/m2)
    real(r8), pointer :: h2osno_old         (:)   => null() ! snow mass for previous time step (kg/m2)
    real(r8), pointer :: h2osno_top         (:)   => null() ! top-layer mass of snow  (kg/m2)
    real(r8), pointer :: sno_liq_top        (:)   => null() ! snow liquid water fraction, by mass, top layer (kg/kg)
    real(r8), pointer :: snowice            (:)   => null() ! total snow ice (kg/m2)
    real(r8), pointer :: snowliq            (:)   => null() ! total snow liquid water (kg/m2)
    real(r8), pointer :: int_snow           (:)   => null() ! integrated snowfall (kg/m2)
    real(r8), pointer :: snow_depth         (:)   => null() ! snow height of snow covered area (m)
    real(r8), pointer :: snowdp             (:)   => null() ! snow height averaged for area with and without snow cover(m)
    real(r8), pointer :: snow_persistence   (:)   => null() ! length of time that ground has had non-zero snow thickness (sec)
    real(r8), pointer :: snw_rds_top        (:)   => null() ! snow grain radius (top layer)  (m^-6, microns)
    logical , pointer :: do_capsnow         (:)   => null() ! true => do snow capping
    ! Area fractions
    real(r8), pointer :: frac_sno           (:)   => null() ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_sno_eff       (:)   => null() ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_iceold        (:,:) => null() ! fraction of ice relative to the tot water (-nlevsno+1:nlevgrnd) 
    real(r8), pointer :: frac_h2osfc        (:)   => null() ! fractional area with surface water greater than zero
    real(r8), pointer :: wf                 (:)   => null() ! soil water as frac. of whc for top 0.05 m (0-1) 
    real(r8), pointer :: wf2                (:)   => null() ! soil water as frac. of whc for top 0.17 m (0-1) 
    real(r8), pointer :: finundated         (:)   => null() ! fraction of column inundated, for bgc caclulation (0-1)
    ! Balance checks
    real(r8), pointer :: begwb              (:)   => null() ! water mass begining of the time step (kg/m2)
    real(r8), pointer :: endwb              (:)   => null() ! water mass end of the time step (kg/m2)
    real(r8), pointer :: errh2o             (:)   => null() ! water conservation error (kg/m2)
    real(r8), pointer :: errh2osno          (:)   => null() ! snow water conservation error(kg/m2)
   
  contains
    procedure, public :: Init    => col_ws_init
    procedure, public :: Restart => col_ws_restart
    procedure, public :: Clean   => col_ws_clean
  end type column_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_cs
    procedure, public :: Clean => clean_col_cs
  end type column_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ns
    procedure, public :: Clean => clean_col_ns
  end type column_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ps
    procedure, public :: Clean => clean_col_ps
  end type column_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_flux
    real(r8), pointer :: eflx_h2osfc_to_snow     (:)   => null() ! snow melt to h2osfc heat flux (W/m2)
    real(r8), pointer :: eflx_snomelt            (:)   => null() ! snow melt heat flux (W/m**2)
    real(r8), pointer :: eflx_snomelt_r          (:)   => null() ! rural snow melt heat flux (W/m2)
    real(r8), pointer :: eflx_snomelt_u          (:)   => null() ! urban snow melt heat flux (W/m2)
    real(r8), pointer :: eflx_bot                (:)   => null() ! heat flux from beneath the soil or ice column (W/m2)
    real(r8), pointer :: eflx_fgr12              (:)   => null() ! ground heat flux between soil layers 1 and 2 (W/m2)
    real(r8), pointer :: eflx_fgr                (:,:) => null() ! (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
    real(r8), pointer :: eflx_building_heat      (:)   => null() ! heat flux from urban building interior to urban walls, roof (W/m2)
    real(r8), pointer :: eflx_urban_ac           (:)   => null() ! urban air conditioning flux (W/m2)
    real(r8), pointer :: eflx_urban_heat         (:)   => null() ! urban heating flux (W/m**2)
    real(r8), pointer :: eflx_hs_h2osfc          (:)   => null() ! heat flux on standing water (W/m2)
    real(r8), pointer :: eflx_hs_top_snow        (:)   => null() ! heat flux on top snow layer (W/m2)
    real(r8), pointer :: eflx_hs_soil            (:)   => null() ! heat flux on soil [W/m2
    real(r8), pointer :: eflx_sabg_lyr           (:,:) => null() ! absorbed solar radiation (col,lyr) (W/m2)
    ! Derivatives of energy fluxes                      
    real(r8), pointer :: eflx_dhsdT              (:)   => null() ! deriv. of energy flux into surface layer wrt temp (W/m2/K)
    ! Latent heat terms
    real(r8), pointer :: htvp                    (:)   => null() ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: xmf                     (:)   => null() ! total latent heat of phase change of ground water
    real(r8), pointer :: xmf_h2osfc              (:)   => null() ! latent heat of phase change of surface water
    integer , pointer :: imelt                   (:,:) ! flag for melting (=1), freezing (=2), Not=0 (-nlevsno+1:nlevgrnd) 
    ! for couplig with pflotran                         
    real(r8), pointer :: eflx_soil_grnd          (:)   => null() ! integrated soil ground heat flux (W/m2)  [+ = into ground]
    real(r8), pointer :: eflx_rnet_soil          (:)   => null() ! soil net (sw+lw) radiation flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_soil          (:)   => null() ! soil-air heat flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_snow          (:)   => null() ! soil-snow heat flux (W/m2) [+ = into soil]
    real(r8), pointer :: eflx_fgr0_h2osfc        (:)   => null() ! soil-surfacewater heat flux (W/m2) [+ = into soil]
    ! Balance Checks                                    
    real(r8), pointer :: errsoi                  (:)   => null() ! soil/lake energy conservation error   (W/m2)
    real(r8), pointer :: errseb                  (:)   => null() ! surface energy conservation error     (W/m2)
    real(r8), pointer :: errsol                  (:)   => null() ! solar radiation conservation error    (W/m2)
    real(r8), pointer :: errlon                  (:)   => null() ! longwave radiation conservation error (W/m2)

  contains
    procedure, public :: Init    => col_ef_init
    procedure, public :: Restart => col_ef_restart
    procedure, public :: Clean   => col_ef_clean
  end type column_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_flux
    real(r8), pointer :: qflx_prec_grnd       (:)   => null() ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_rain_grnd       (:)   => null() ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd       (:)   => null() ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow        (:)   => null() ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow_vol    (:)   => null() !
    real(r8), pointer :: qflx_evap_soi        (:)   => null() ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_veg        (:)   => null() ! vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can        (:)   => null() ! evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot        (:)   => null() ! col_qflx_evap_soi + col_qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_evap_grnd       (:)   => null() ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_liq       (:)   => null() ! excess rainfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_ice       (:)   => null() ! excess snowfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_tran_veg        (:)   => null() ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_dew_snow        (:)   => null() ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd        (:)   => null() ! ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
    real(r8), pointer :: qflx_prec_intr       (:)   => null() ! interception of precipitation [mm/s]
    real(r8), pointer :: qflx_ev_snow         (:)   => null() ! evaporation heat flux from snow         (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_soil         (:)   => null() ! evaporation heat flux from soil         (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_ev_h2osfc       (:)   => null() ! evaporation heat flux from soil         (W/m**2) [+ to atm] ! NOTE: unit shall be mm H2O/s for water NOT heat
    real(r8), pointer :: qflx_gross_evap_soil (:)   => null() ! gross infiltration from soil, this satisfies the relationship qflx_infl = qflx_gross_infl_soil-qflx_gross_evap_soil
    real(r8), pointer :: qflx_gross_infl_soil (:)   => null() ! gross infiltration, before considering the evaporation
    real(r8), pointer :: qflx_adv             (:,:) => null() ! advective flux across different soil layer interfaces [mm H2O/s] [+ downward]
    real(r8), pointer :: qflx_rootsoi         (:,:) => null() ! root and soil water exchange [mm H2O/s] [+ into root]     
    real(r8), pointer :: dwb                  (:)   => null() !  water mass change [+ increase] [mm H2O/s] 
    real(r8), pointer :: qflx_infl            (:)   => null() ! infiltration (mm H2O /s)
    real(r8), pointer :: qflx_surf            (:)   => null() ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_drain           (:)   => null() ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_totdrain        (:)   => null() 
    real(r8), pointer :: qflx_top_soil        (:)   => null() ! net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_h2osfc_to_ice   (:)   => null() ! conversion of h2osfc to ice
    real(r8), pointer :: qflx_h2osfc_surf     (:)   => null() ! surface water runoff
    real(r8), pointer :: qflx_snow_h2osfc     (:)   => null() ! snow falling on surface water
    real(r8), pointer :: qflx_drain_perched   (:)   => null() ! sub-surface runoff from perched wt (mm H2O /s)
    real(r8), pointer :: qflx_deficit         (:)   => null() ! water deficit to keep non-negative liquid water content (mm H2O)   
    real(r8), pointer :: qflx_floodc          (:)   => null() ! flood water flux at column level
    real(r8), pointer :: qflx_sl_top_soil     (:)   => null() ! liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
    real(r8), pointer :: qflx_snomelt         (:)   => null() ! snow melt (mm H2O /s)
    real(r8), pointer :: qflx_snow_melt       (:)   => null() ! snow melt (net)
    real(r8), pointer :: qflx_qrgwl           (:)   => null() ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_runoff          (:)   => null() ! total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: qflx_runoff_r        (:)   => null() ! Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
    real(r8), pointer :: qflx_runoff_u        (:)   => null() ! urban total runoff (qflx_drain+qflx_surf) (mm H2O /s) 
    real(r8), pointer :: qflx_rsub_sat        (:)   => null() ! soil saturation excess [mm/s]
    real(r8), pointer :: qflx_snofrz_lyr      (:,:) => null() ! snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
    real(r8), pointer :: qflx_snofrz          (:)   => null() ! column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]
    real(r8), pointer :: qflx_glcice          (:)   => null() ! net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC
    real(r8), pointer :: qflx_glcice_frz      (:)   => null() ! ice growth (positive definite) (mm H2O/s)
    real(r8), pointer :: qflx_glcice_melt     (:)   => null() ! ice melt (positive definite) (mm H2O/s)
    real(r8), pointer :: qflx_drain_vr        (:,:) => null() ! liquid water lost as drainage (m /time step)
    real(r8), pointer :: qflx_h2osfc2topsoi   (:)   => null() ! liquid water coming from surface standing water top soil (mm H2O/s)
    real(r8), pointer :: qflx_snow2topsoi     (:)   => null() ! liquid water coming from residual snow to topsoil (mm H2O/s)
    real(r8), pointer :: qflx_lateral         (:)   => null() ! lateral subsurface flux (mm H2O /s)
    real(r8), pointer :: snow_sources         (:)   => null() ! snow sources (mm H2O/s)
    real(r8), pointer :: snow_sinks           (:)   => null() ! snow sinks (mm H2O/s)
    real(r8), pointer :: qflx_irrig           (:)   => null() ! irrigation flux (mm H2O/s)

  contains
    procedure, public :: Init    => col_wf_init
    procedure, public :: Restart => col_wf_restart
    procedure, public :: Clean   => col_wf_clean
  end type column_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_cf
    procedure, public :: Clean => clean_col_cf
  end type column_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_nf
    procedure, public :: Clean => clean_col_nf
  end type column_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_pf
    procedure, public :: Clean => clean_col_pf
  end type column_phosphorus_flux

   
  !-----------------------------------------------------------------------
  ! declare the public instances of column-level data types
  !-----------------------------------------------------------------------
  type(column_energy_state)          , public, target :: col_es    ! column energy state
  type(column_water_state)           , public, target :: col_ws    ! column water state
  type(column_carbon_state)          , public, target :: col_cs    ! column carbon state
  type(column_nitrogen_state)        , public, target :: col_ns    ! column nitrogen state
  type(column_phosphorus_state)      , public, target :: col_ps    ! column phosphorus state
  type(column_energy_flux)           , public, target :: col_ef    ! column energy flux
  type(column_water_flux)            , public, target :: col_wf    ! column water flux
  type(column_carbon_flux)           , public, target :: col_cf    ! column carbon flux
  type(column_nitrogen_flux)         , public, target :: col_nf    ! column nitrogen flux
  type(column_phosphorus_flux)       , public, target :: col_pf    ! column phosphorus flux

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy state data structure
  !------------------------------------------------------------------------
  subroutine col_es_init(this, begc, endc)
    !
    ! !USES:
    use landunit_varcon, only : istice, istwet, istsoil, istdlak, istice_mec
    use clm_varctl     , only : iulog, use_cn, use_vancouver, use_mexicocity
    use column_varcon  , only : icol_road_perv, icol_road_imperv, icol_roof, icol_sunwall, icol_shadewall
    use UrbanParamsType, only : urbanparams_vars
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer           :: c,l,j                        ! indices
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays

    !------------------------------------------------------------------------------
    associate(snl => col_pp%snl) ! Output: [integer (:)    ]  number of snow layers   

    !-----------------------------------------------------------------------
    ! allocate for each member of col_es
    !-----------------------------------------------------------------------
    allocate(this%t_soisno         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_soisno           (:,:) = nan
    allocate(this%t_ssbef          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_ssbef            (:,:) = nan
    allocate(this%t_h2osfc         (begc:endc))                     ; this%t_h2osfc           (:)   = nan
    allocate(this%t_h2osfc_bef     (begc:endc))                     ; this%t_h2osfc_bef       (:)   = nan
    allocate(this%t_soi10cm        (begc:endc))                     ; this%t_soi10cm          (:)   = nan
    allocate(this%t_soi17cm        (begc:endc))                     ; this%t_soi17cm          (:)   = spval
    allocate(this%t_grnd           (begc:endc))                     ; this%t_grnd             (:)   = nan
    allocate(this%t_lake           (begc:endc,1:nlevlak))           ; this%t_lake             (:,:) = nan
    allocate(this%t_grnd_r         (begc:endc))                     ; this%t_grnd_r           (:)   = nan
    allocate(this%t_grnd_u         (begc:endc))                     ; this%t_grnd_u           (:)   = nan
    allocate(this%snot_top         (begc:endc))                     ; this%snot_top           (:)   = nan
    allocate(this%dTdz_top         (begc:endc))                     ; this%dTdz_top           (:)   = nan
    allocate(this%thv              (begc:endc))                     ; this%thv                (:)   = nan
    allocate(this%hc_soi           (begc:endc))                     ; this%hc_soi             (:)   = nan
    allocate(this%hc_soisno        (begc:endc))                     ; this%hc_soisno          (:)   = nan
    allocate(this%emg              (begc:endc))                     ; this%emg                (:)   = nan
    allocate(this%fact             (begc:endc, -nlevsno+1:nlevgrnd)); this%fact               (:,:) = nan
    allocate(this%c_h2osfc         (begc:endc))                     ; this%c_h2osfc           (:)   = nan
    allocate(this%t_nearsurf       (begc:endc))                     ; this%t_nearsurf         (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of col_es
    !-----------------------------------------------------------------------
    this%t_soisno(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%t_soisno(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_T', units='K', type2d='levsno',  &
         avgflag='A', long_name='Snow temperatures', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%t_soisno(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (vegetated landunits only)', &
         ptr_col=this%t_soisno, l2g_scale_type='veg')

    this%t_soisno(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI_ICE',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (ice landunits only)', &
         ptr_col=this%t_soisno, l2g_scale_type='ice')

    this%t_h2osfc(begc:endc) = spval
    call hist_addfld1d (fname='TH2OSFC',  units='K',  &
         avgflag='A', long_name='surface water temperature', &
         ptr_col=this%t_h2osfc)

    this%t_soi10cm(begc:endc) = spval
    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=this%t_soi10cm, set_urb=spval)

    this%t_grnd(begc:endc) = spval
    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=this%t_grnd, c2l_scale_type='urbans')

    this%t_lake(begc:endc,:) = spval
    call hist_addfld2d (fname='TLAKE',  units='K', type2d='levlak', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=this%t_lake)

    this%t_grnd_r(begc:endc) = spval
    call hist_addfld1d (fname='TG_R', units='K',  &
         avgflag='A', long_name='Rural ground temperature', &
         ptr_col=this%t_grnd_r, set_spec=spval)

    this%t_grnd_u(begc:endc) = spval
    call hist_addfld1d (fname='TG_U', units='K',  &
         avgflag='A', long_name='Urban ground temperature', &
         ptr_col=this%t_grnd_u, set_nourb=spval, c2l_scale_type='urbans')

    this%snot_top(begc:endc) = spval 
    call hist_addfld1d (fname='SNOTTOPL', units='K', &
         avgflag='A', long_name='snow temperature (top layer)', &
         ptr_col=this%snot_top, set_urb=spval, default='inactive')

    this%dTdz_top(begc:endc) = spval 
    call hist_addfld1d (fname='SNOdTdzL', units='K/m', &
         avgflag='A', long_name='top snow layer temperature gradient (land)', &
         ptr_col=this%dTdz_top, set_urb=spval, default='inactive')

    this%hc_soi(begc:endc) = spval
    call hist_addfld1d (fname='HCSOI',  units='MJ/m2',  &
         avgflag='A', long_name='soil heat content', &
         ptr_col=this%hc_soi, set_lake=spval, set_urb=spval, l2g_scale_type='veg')

    this%hc_soisno(begc:endc) = spval
    call hist_addfld1d (fname='HC',  units='MJ/m2',  &
         avgflag='A', long_name='heat content of soil/snow/lake', &
         ptr_col=this%hc_soisno, set_urb=spval)

    if (use_cn) then
       this%emg(begc:endc) = spval
       call hist_addfld1d (fname='EMG', units='proportion', &
            avgflag='A', long_name='ground emissivity', &
            ptr_col=this%emg, default='inactive')
    end if

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of col_es
    !-----------------------------------------------------------------------
    
    ! Initialize soil+snow temperatures
    do c = begc,endc
       l = col_pp%landunit(c)

       ! Snow level temperatures - all land points
       if (snl(c) < 0) then
          do j = snl(c)+1, 0
             this%t_soisno(c,j) = 250._r8
          end do
       end if

       ! Below snow temperatures - nonlake points (lake points are set below)
       if (.not. lun_pp%lakpoi(l)) then 

          if (lun_pp%itype(l)==istice .or. lun_pp%itype(l)==istice_mec) then
             this%t_soisno(c,1:nlevgrnd) = 250._r8

          else if (lun_pp%itype(l) == istwet) then
             this%t_soisno(c,1:nlevgrnd) = 277._r8

          else if (lun_pp%urbpoi(l)) then
             if (use_vancouver) then
                if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 20C in bottom layer
                   do j = 1, nlevgrnd
                      this%t_soisno(c,j) = 297.56 - (j-1) * ((297.56-293.16)/(nlevgrnd-1)) 
                   end do
                   ! Set wall and roof layers to initial air temperature
                else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                     .or. col_pp%itype(c) == icol_roof) then
                   this%t_soisno(c,1:nlevurb) = 297.56
                else
                   this%t_soisno(c,1:nlevgrnd) = 283._r8
                end if
             else if (use_mexicocity) then
                if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                   ! Set road top layer to initial air temperature and interpolate other
                   ! layers down to 22C in bottom layer
                   do j = 1, nlevgrnd
                      this%t_soisno(c,j) = 289.46 - (j-1) * ((289.46-295.16)/(nlevgrnd-1)) 
                   end do
                else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                     .or. col_pp%itype(c) == icol_roof) then
                   ! Set wall and roof layers to initial air temperature
                   this%t_soisno(c,1:nlevurb) = 289.46
                else
                   this%t_soisno(c,1:nlevgrnd) = 283._r8
                end if
             else
                if (col_pp%itype(c) == icol_road_perv .or. col_pp%itype(c) == icol_road_imperv) then 
                   this%t_soisno(c,1:nlevgrnd) = 274._r8
                else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                     .or. col_pp%itype(c) == icol_roof) then
                   ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
                   ! shock from large heating/air conditioning flux
                   this%t_soisno(c,1:nlevurb) = 292._r8
                end if
             end if
          else
             this%t_soisno(c,1:nlevgrnd) = 274._r8
          endif
          this%t_grnd(c) = this%t_soisno(c,snl(c)+1)
       endif
       
       if (lun_pp%lakpoi(l)) then ! special handling for lake points
          this%t_soisno(c,1:nlevgrnd) = 277._r8
          this%t_lake(c,1:nlevlak) = 277._r8
          this%t_grnd(c) = 277._r8
       end if

       this%t_soi17cm(c) = this%t_grnd(c)
       
       ! set emissivity for urban columns
       if (col_pp%itype(c) == icol_roof       ) this%emg(c) = urbanparams_vars%em_roof(l)
       if (col_pp%itype(c) == icol_sunwall    ) this%emg(c) = urbanparams_vars%em_wall(l)
       if (col_pp%itype(c) == icol_shadewall  ) this%emg(c) = urbanparams_vars%em_wall(l)
       if (col_pp%itype(c) == icol_road_imperv) this%emg(c) = urbanparams_vars%em_improad(l)
       if (col_pp%itype(c) == icol_road_perv  ) this%emg(c) = urbanparams_vars%em_perroad(l)

    end do ! columns loop

    ! Initialize surface water temperatures
    this%t_h2osfc = 274._r8
    
    end associate

  end subroutine col_es_init
    
  !------------------------------------------------------------------------
  subroutine col_es_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write column energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='T_SOISNO', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='soil-snow temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_soisno)

    call restartvar(ncid=ncid, flag=flag, varname='TH2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_h2osfc)
    if (flag=='read' .and. .not. readvar) then
       this%t_h2osfc(bounds%begc:bounds%endc) = 274.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd)

    call restartvar(ncid=ncid, flag=flag, varname='T_LAKE', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_lake)

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_R', xtype=ncd_double,  &
         dim1name='column', &
         long_name='rural ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_r)
         
    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_U', xtype=ncd_double, &
         dim1name='column',                    &
         long_name='urban ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_u)

  end subroutine col_es_restart

  !------------------------------------------------------------------------
  subroutine col_es_clean(this)
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%t_h2osfc)
  end subroutine col_es_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water state data structure
  !------------------------------------------------------------------------
  subroutine col_ws_init(this, begc, endc, h2osno_input, snow_depth_input, watsat_input)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    integer , intent(in)      :: begc,endc
    real(r8), intent(in)      :: h2osno_input(begc:)
    real(r8), intent(in)      :: snow_depth_input(begc:)
    real(r8), intent(in)      :: watsat_input(begc:, 1:)          ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer  :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    real(r8)           :: snowbd      ! temporary calculation of snow bulk density (kg/m3)
    real(r8)           :: fmelt       ! snowbd/100
    integer            :: c,l,j,nlevs,nlevbed
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_ws
    !-----------------------------------------------------------------------
    allocate(this%h2osoi_liq         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq         (:,:) = nan
    allocate(this%h2osoi_ice         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice         (:,:) = nan
    allocate(this%h2osoi_vol         (begc:endc, 1:nlevgrnd))         ; this%h2osoi_vol         (:,:) = nan
    allocate(this%h2osfc             (begc:endc))                     ; this%h2osfc             (:)   = nan   
    allocate(this%h2ocan             (begc:endc))                     ; this%h2ocan             (:)   = nan  
    allocate(this%total_plant_stored_h2o(begc:endc))                  ; this%total_plant_stored_h2o(:)= nan  
    allocate(this%h2osoi_liqvol      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liqvol      (:,:) = nan
    allocate(this%h2osoi_icevol      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_icevol      (:,:) = nan    
    allocate(this%h2osoi_liq_old     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_old     (:,:) = nan
    allocate(this%h2osoi_ice_old     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_old     (:,:) = nan 
    allocate(this%bw                 (begc:endc,-nlevsno+1:0))        ; this%bw                 (:,:) = nan   
    allocate(this%smp_l              (begc:endc,-nlevsno+1:nlevgrnd)) ; this%smp_l              (:,:) = nan
    allocate(this%soilp              (begc:endc,1:nlevgrnd))          ; this%soilp              (:,:) = 0._r8
    allocate(this%swe_old            (begc:endc,-nlevsno+1:0))        ; this%swe_old            (:,:) = nan   
    allocate(this%snw_rds            (begc:endc,-nlevsno+1:0))        ; this%snw_rds            (:,:) = nan
    allocate(this%air_vol            (begc:endc, 1:nlevgrnd))         ; this%air_vol            (:,:) = nan
    allocate(this%qg_snow            (begc:endc))                     ; this%qg_snow            (:)   = nan   
    allocate(this%qg_soil            (begc:endc))                     ; this%qg_soil            (:)   = nan   
    allocate(this%qg_h2osfc          (begc:endc))                     ; this%qg_h2osfc          (:)   = nan   
    allocate(this%qg                 (begc:endc))                     ; this%qg                 (:)   = nan   
    allocate(this%dqgdT              (begc:endc))                     ; this%dqgdT              (:)   = nan   
    allocate(this%h2osoi_liqice_10cm (begc:endc))                     ; this%h2osoi_liqice_10cm (:)   = nan
    allocate(this%h2osno             (begc:endc))                     ; this%h2osno             (:)   = nan   
    allocate(this%h2osno_old         (begc:endc))                     ; this%h2osno_old         (:)   = nan   
    allocate(this%h2osno_top         (begc:endc))                     ; this%h2osno_top         (:)   = nan
    allocate(this%sno_liq_top        (begc:endc))                     ; this%sno_liq_top        (:)   = nan
    allocate(this%snowice            (begc:endc))                     ; this%snowice            (:)   = nan   
    allocate(this%snowliq            (begc:endc))                     ; this%snowliq            (:)   = nan   
    allocate(this%int_snow           (begc:endc))                     ; this%int_snow           (:)   = nan   
    allocate(this%snow_depth         (begc:endc))                     ; this%snow_depth         (:)   = nan
    allocate(this%snowdp             (begc:endc))                     ; this%snowdp             (:)   = nan
    allocate(this%snow_persistence   (begc:endc))                     ; this%snow_persistence   (:)   = nan
    allocate(this%snw_rds_top        (begc:endc))                     ; this%snw_rds_top        (:)   = nan
    allocate(this%do_capsnow         (begc:endc))                   
    allocate(this%frac_sno           (begc:endc))                     ; this%frac_sno           (:)   = nan
    allocate(this%frac_sno_eff       (begc:endc))                     ; this%frac_sno_eff       (:)   = nan
    allocate(this%frac_iceold        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%frac_iceold        (:,:) = nan
    allocate(this%frac_h2osfc        (begc:endc))                     ; this%frac_h2osfc        (:)   = nan 
    allocate(this%wf                 (begc:endc))                     ; this%wf                 (:)   = nan
    allocate(this%wf2                (begc:endc))                     ; this%wf2                (:)   = nan
    allocate(this%finundated         (begc:endc))                     ; this%finundated         (:)   = nan
    allocate(this%begwb              (begc:endc))                     ; this%begwb              (:)   = nan
    allocate(this%endwb              (begc:endc))                     ; this%endwb              (:)   = nan
    allocate(this%errh2o             (begc:endc))                     ; this%errh2o             (:)   = nan
    allocate(this%errh2osno          (begc:endc))                     ; this%errh2osno          (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of col_ws
    !-----------------------------------------------------------------------
    data2dptr => this%h2osoi_liq(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_LIQH2O', units='kg/m2', type2d='levsno',  &
         avgflag='A', long_name='Snow liquid water content', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    data2dptr => this%h2osoi_ice(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_ICE', units='kg/m2', type2d='levsno',  &
         avgflag='A', long_name='Snow ice content', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%h2osoi_liq(begc:endc,:) = spval
    call hist_addfld2d (fname='SOILLIQ',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil liquid water (vegetated landunits only)', &
         ptr_col=this%h2osoi_liq, l2g_scale_type='veg')

    this%h2osoi_ice(begc:endc,:) = spval
    call hist_addfld2d (fname='SOILICE',  units='kg/m2', type2d='levgrnd', &
         avgflag='A', long_name='soil ice (vegetated landunits only)', &
         ptr_col=this%h2osoi_ice, l2g_scale_type='veg')

    this%h2osfc(begc:endc) = spval
    call hist_addfld1d (fname='H2OSFC',  units='mm',  &
         avgflag='A', long_name='surface water depth', &
         ptr_col=this%h2osfc)

    this%h2osoi_vol(begc:endc,:) = spval
    call hist_addfld2d (fname='H2OSOI',  units='mm3/mm3', type2d='levgrnd', &
         avgflag='A', long_name='volumetric soil water (vegetated landunits only)', &
         ptr_col=this%h2osoi_vol, l2g_scale_type='veg')

    this%bw(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%bw(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_BW', units='kg/m3', type2d='levsno', &
         avgflag='A', long_name='Partial density of water in the snow pack (ice + liquid)', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%snw_rds(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%snw_rds(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_GS', units='Microns', type2d='levsno',  &
         avgflag='A', long_name='Mean snow grain size', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%snw_rds_top(begc:endc) = spval 
    call hist_addfld1d (fname='SNORDSL', units='m^-6', &
         avgflag='A', long_name='top snow layer effective grain radius', &
         ptr_col=this%snw_rds_top, set_urb=spval, default='inactive')

    this%sno_liq_top(begc:endc) = spval 
    call hist_addfld1d (fname='SNOLIQFL', units='fraction', &
         avgflag='A', long_name='top snow layer liquid water fraction (land)', &
         ptr_col=this%sno_liq_top, set_urb=spval, default='inactive')

    this%h2osoi_liqice_10cm(begc:endc) = spval
    call hist_addfld1d (fname='SOILWATER_10CM',  units='kg/m2', &
         avgflag='A', long_name='soil liquid water + ice in top 10cm of soil (veg landunits only)', &
         ptr_col=this%h2osoi_liqice_10cm, set_urb=spval, set_lake=spval, l2g_scale_type='veg')

    call hist_addfld1d (fname='H2OSNO',  units='mm',  &
         avgflag='A', long_name='snow depth (liquid water)', &
         ptr_col=this%h2osno, c2l_scale_type='urbanf')

    this%h2osno_top(begc:endc) = spval
    call hist_addfld1d (fname='H2OSNO_TOP', units='kg/m2', &
         avgflag='A', long_name='mass of snow in top snow layer', &
         ptr_col=this%h2osno_top, set_urb=spval)

    this%snowice(begc:endc) = spval
    call hist_addfld1d (fname='SNOWICE',  units='kg/m2', &
         avgflag='A', long_name='snow ice', &
         ptr_col=this%snowice)

    this%snowliq(begc:endc) = spval
    call hist_addfld1d (fname='SNOWLIQ',  units='kg/m2',  &
         avgflag='A', long_name='snow liquid water', &
         ptr_col=this%snowliq)

    this%int_snow(begc:endc) = spval
    call hist_addfld1d (fname='INT_SNOW',  units='mm',  &
         avgflag='A', long_name='accumulated swe (vegetated landunits only)', &
         ptr_col=this%int_snow, l2g_scale_type='veg')

    this%snow_depth(begc:endc) = spval
    call hist_addfld1d (fname='SNOW_DEPTH',  units='m',  &
         avgflag='A', long_name='snow height of snow covered area', &
         ptr_col=this%snow_depth, c2l_scale_type='urbanf')!, default='inactive')

    this%snowdp(begc:endc) = spval
    call hist_addfld1d (fname='SNOWDP',  units='m',  &
         avgflag='A', long_name='gridcell mean snow height', &
         ptr_col=this%snowdp, c2l_scale_type='urbanf')

    if (create_glacier_mec_landunit) then
       this%snow_persistence(begc:endc) = spval
       call hist_addfld1d (fname='SNOW_PERSISTENCE',  units='seconds',  &
            avgflag='I', long_name='Length of time of continuous snow cover (nat. veg. landunits only)', &
            ptr_col=this%snow_persistence, l2g_scale_type='natveg') 
    end if

    this%frac_sno(begc:endc) = spval
    call hist_addfld1d (fname='FSNO',  units='unitless',  &
         avgflag='A', long_name='fraction of ground covered by snow', &
         ptr_col=this%frac_sno, c2l_scale_type='urbanf')

    this%frac_sno_eff(begc:endc) = spval
    call hist_addfld1d (fname='FSNO_EFF',  units='unitless',  &
         avgflag='A', long_name='effective fraction of ground covered by snow', &
         ptr_col=this%frac_sno_eff, c2l_scale_type='urbanf')!, default='inactive')

    if (use_cn)then
       this%frac_iceold(begc:endc,:) = spval
       call hist_addfld2d (fname='FRAC_ICEOLD', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='fraction of ice relative to the tot water', &
            ptr_col=this%frac_iceold, default='inactive')
    end if

    this%frac_h2osfc(begc:endc) = spval
    call hist_addfld1d (fname='FH2OSFC',  units='unitless',  &
         avgflag='A', long_name='fraction of ground covered by surface water', &
         ptr_col=this%frac_h2osfc)

    if (use_cn) then
       this%wf(begc:endc) = spval
       call hist_addfld1d (fname='WF', units='proportion', &
            avgflag='A', long_name='soil water as frac. of whc for top 0.05 m', &
            ptr_col=this%wf)
    end if

    this%errh2o(begc:endc) = spval
    call hist_addfld1d (fname='ERRH2O', units='mm',  &
         avgflag='A', long_name='total water conservation error', &
         ptr_col=this%errh2o)

    this%errh2osno(begc:endc) = spval
    call hist_addfld1d (fname='ERRH2OSNO',  units='mm',  &
         avgflag='A', long_name='imbalance in snow depth (liquid water)', &
         ptr_col=this%errh2osno, c2l_scale_type='urbanf')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of col_ws
    !-----------------------------------------------------------------------
    
    ! Arrays that are initialized from input arguments
    do c = begc,endc
       l = col_pp%landunit(c)
       this%h2osno(c)                 = h2osno_input(c) 
       this%int_snow(c)               = h2osno_input(c) 
       this%snow_depth(c)             = snow_depth_input(c)
       this%snow_persistence(c)       = 0._r8
       this%wf(c)                     = spval
       this%wf2(c)                    = spval
       this%total_plant_stored_h2o(c) = 0._r8
       this%h2osfc(c)                 = 0._r8
       this%h2ocan(c)                 = 0._r8
       this%frac_h2osfc(c)            = 0._r8

       if (lun_pp%urbpoi(l)) then
          ! From Bonan 1996 (LSM technical note)
          this%frac_sno(c) = min( this%snow_depth(c)/0.05_r8, 1._r8)
       else
          this%frac_sno(c) = 0._r8
          ! snow cover fraction as in Niu and Yang 2007
          if(this%snow_depth(c) > 0.0)  then
             snowbd   = min(400._r8, this%h2osno(c)/this%snow_depth(c)) !bulk density of snow (kg/m3)
             fmelt    = (snowbd/100.)**1.
             ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
             ! reconsidered, optimal value of 1.5 in Niu et al., 2007
             this%frac_sno(c) = tanh( this%snow_depth(c) /(2.5 * zlnd * fmelt) )
          endif
       end if
       
       if (col_pp%snl(c) < 0) then
          this%snw_rds(c,col_pp%snl(c)+1:0)          = snw_rds_min
          this%snw_rds(c,-nlevsno+1:col_pp%snl(c))   = 0._r8
          this%snw_rds_top(c)                 = snw_rds_min
       elseif (this%h2osno(c) > 0._r8) then   
          this%snw_rds(c,0)                   = snw_rds_min
          this%snw_rds(c,-nlevsno+1:-1)       = 0._r8
          this%snw_rds_top(c)                 = spval
          this%sno_liq_top(c)                 = spval
       else                                   
          this%snw_rds(c,:)                   = 0._r8
          this%snw_rds_top(c)                 = spval
          this%sno_liq_top(c)                 = spval
       endif
       
       !--------------------------------------------
       ! Set soil water
       !--------------------------------------------
       
       ! volumetric water is set first and liquid content and ice lens are obtained
       ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
       ! and urban pervious road (other urban columns have zero soil water)
       
       this%h2osoi_vol(c,         1:) = spval
       this%h2osoi_liq(c,-nlevsno+1:) = spval
       this%h2osoi_ice(c,-nlevsno+1:) = spval
       
       if (.not. lun_pp%lakpoi(l)) then  !not lake
	       nlevbed = col_pp%nlevbed(c)
          ! volumetric water
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevbed) then
                   this%h2osoi_vol(c,j) = 0.0_r8
                else
		             if (use_fates_planthydro) then
                      this%h2osoi_vol(c,j) = 0.70_r8*watsat_input(c,j) !0.15_r8 to avoid very dry conditions that cause errors in FATES HYDRO
                   else
                      this%h2osoi_vol(c,j) = 0.15_r8
                   endif
                endif
             end do
          else if (lun_pp%urbpoi(l)) then
             if (col_pp%itype(c) == icol_road_perv) then
                nlevs = nlevgrnd
                do j = 1, nlevs
                   if (j <= nlevbed) then
                      this%h2osoi_vol(c,j) = 0.3_r8
                   else
                      this%h2osoi_vol(c,j) = 0.0_r8
                   end if
                end do
             else if (col_pp%itype(c) == icol_road_imperv) then
                nlevs = nlevgrnd
                do j = 1, nlevs
                   this%h2osoi_vol(c,j) = 0.0_r8
                end do
             else
                nlevs = nlevurb
                do j = 1, nlevs
                   this%h2osoi_vol(c,j) = 0.0_r8
                end do
             end if
          else if (lun_pp%itype(l) == istwet) then
             nlevs = nlevgrnd
             do j = 1, nlevs
                if (j > nlevbed) then
                   this%h2osoi_vol(c,j) = 0.0_r8
                else
                   this%h2osoi_vol(c,j) = 1.0_r8
                endif
             end do
          else if (lun_pp%itype(l) == istice .or. lun_pp%itype(l) == istice_mec) then
             nlevs = nlevgrnd 
             do j = 1, nlevs
                this%h2osoi_vol(c,j) = 1.0_r8
             end do
          endif
          do j = 1, nlevs
             this%h2osoi_vol(c,j) = min(this%h2osoi_vol(c,j), watsat_input(c,j))

             if (col_es%t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                this%h2osoi_ice(c,j) = col_pp%dz(c,j)*denice*this%h2osoi_vol(c,j)
                this%h2osoi_liq(c,j) = 0._r8
             else
                this%h2osoi_ice(c,j) = 0._r8
                this%h2osoi_liq(c,j) = col_pp%dz(c,j)*denh2o*this%h2osoi_vol(c,j)
             endif
          end do
          do j = -nlevsno+1, 0
             if (j > col_pp%snl(c)) then
                this%h2osoi_ice(c,j) = col_pp%dz(c,j)*250._r8
                this%h2osoi_liq(c,j) = 0._r8
             end if
          end do
       end if
       !--------------------------------------------
       ! Set Lake water
       !--------------------------------------------
       if (lun_pp%lakpoi(l)) then
          do j = -nlevsno+1, 0
             if (j > col_pp%snl(c)) then
                this%h2osoi_ice(c,j) = col_pp%dz(c,j)*bdsno
                this%h2osoi_liq(c,j) = 0._r8
             end if
          end do
          do j = 1,nlevgrnd
             if (j <= nlevsoi) then ! soil
                this%h2osoi_vol(c,j) = watsat_input(c,j)
                this%h2osoi_liq(c,j) = spval
                this%h2osoi_ice(c,j) = spval
             else                  ! bedrock
                this%h2osoi_vol(c,j) = 0._r8
             end if
          end do
       end if

       !--------------------------------------------
       ! For frozen layers !TODO - does the following make sense ???? it seems to overwrite everything
       !--------------------------------------------
       do j = 1,nlevgrnd
          if (col_es%t_soisno(c,j) <= tfrz) then
             this%h2osoi_ice(c,j) = col_pp%dz(c,j)*denice*this%h2osoi_vol(c,j)
             this%h2osoi_liq(c,j) = 0._r8
          else
             this%h2osoi_ice(c,j) = 0._r8
             this%h2osoi_liq(c,j) = col_pp%dz(c,j)*denh2o*this%h2osoi_vol(c,j)
          endif
       end do
       
       this%h2osoi_liq_old(c,:) = this%h2osoi_liq(c,:)
       this%h2osoi_ice_old(c,:) = this%h2osoi_ice(c,:)
    end do
    
  end subroutine col_ws_init
    
  !------------------------------------------------------------------------
  subroutine col_ws_restart(this, bounds, ncid, flag, watsat_input)
    ! 
    ! !DESCRIPTION:
    ! Read/Write column water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    real(r8)         , intent(in)    :: watsat_input (bounds%begc:, 1:)  ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: c,l,j,nlevs  ! indices
    real(r8):: maxwatsat    ! maximum porosity    
    real(r8):: excess       ! excess volumetric soil water
    real(r8):: totwat       ! total soil water (mm)
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='H2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osfc)
    if (flag=='read' .and. .not. readvar) then
       this%h2osfc(bounds%begc:bounds%endc) = 0.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='H2OSOI_LIQ', xtype=ncd_double,  &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='liquid water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_liq)

    call restartvar(ncid=ncid, flag=flag, varname='H2OSOI_ICE', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='ice lens', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_ice)
         
    call restartvar(ncid=ncid, flag=flag, varname='SOILP', xtype=ncd_double,  &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='soil pressure ', units='Pa', &
         interpinic_flag='interp', readvar=readvar, data=this%soilp)

    call restartvar(ncid=ncid, flag=flag, varname='H2OSNO', xtype=ncd_double,  &
         dim1name='column', &
         long_name='snow water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osno)

    call restartvar(ncid=ncid, flag=flag, varname='snw_rds', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer effective radius', units='um', &
         interpinic_flag='interp', readvar=readvar, data=this%snw_rds)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize snw_rds
       if (masterproc) then
          write(iulog,*) "SNICAR: This is an initial run (not a restart), and grain size/aerosol " // &
               "mass data are not defined in initial condition file. Initialize snow " // &
               "effective radius to fresh snow value, and snow/aerosol masses to zero."
       endif

       do c= bounds%begc, bounds%endc
          if (col_pp%snl(c) < 0) then
             this%snw_rds(c,col_pp%snl(c)+1:0) = snw_rds_min
             this%snw_rds(c,-nlevsno+1:col_pp%snl(c)) = 0._r8
             this%snw_rds_top(c) = snw_rds_min
             this%sno_liq_top(c) = this%h2osoi_liq(c,col_pp%snl(c)+1) / &
                                      (this%h2osoi_liq(c,col_pp%snl(c)+1)+this%h2osoi_ice(c,col_pp%snl(c)+1))
          elseif (this%h2osno(c) > 0._r8) then
             this%snw_rds(c,0) = snw_rds_min
             this%snw_rds(c,-nlevsno+1:-1) = 0._r8
             this%snw_rds_top(c) = spval
             this%sno_liq_top(c) = spval
          else
             this%snw_rds(c,:) = 0._r8
             this%snw_rds_top(c) = spval
             this%sno_liq_top(c) = spval
          endif
       enddo
    endif
 
    call restartvar(ncid=ncid, flag=flag, varname='INT_SNOW', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='accuumulated snow', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%int_snow)
    if (flag=='read' .and. .not. readvar) then
       this%int_snow(:) = 0.0_r8
    end if
 
    call restartvar(ncid=ncid, flag=flag, varname='SNOW_DEPTH', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='snow depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%snow_depth) 

    call restartvar(ncid=ncid, flag=flag, varname='SNOW_PERS', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='continuous snow cover time', units='sec', &
         interpinic_flag='interp', readvar=readvar, data=this%snow_persistence)    
    if (flag=='read' .and. .not. readvar) then
         this%snow_persistence(:) = 0.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='frac_sno', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='fraction of ground covered by snow (0 to 1)',units='unitless',&
         interpinic_flag='interp', readvar=readvar, data=this%frac_sno)

    call restartvar(ncid=ncid, flag=flag, varname='frac_sno_eff', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='fraction of ground covered by snow (0 to 1)',units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%frac_sno_eff)
    if (flag == 'read' .and. .not. readvar) then
       this%frac_sno_eff(bounds%begc:bounds%endc) = 0.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='FH2OSFC', xtype=ncd_double,  &
         dim1name='column',&
         long_name='fraction of ground covered by h2osfc (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frac_h2osfc)
    if (flag == 'read' .and. .not. readvar) then
       this%frac_h2osfc(bounds%begc:bounds%endc) = 0.0_r8
    end if

    if (use_cn) then
       call restartvar(ncid=ncid, flag=flag, varname='wf', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%wf) 
    end if

    ! Determine volumetric soil water (for read only)
    if (flag == 'read' ) then
       do c = bounds%begc, bounds%endc
          l = col_pp%landunit(c)
          if ( col_pp%itype(c) == icol_sunwall   .or. &
               col_pp%itype(c) == icol_shadewall .or. &
               col_pp%itype(c) == icol_roof )then
             nlevs = nlevurb
          else
             nlevs = nlevgrnd
          end if
          if ( lun_pp%itype(l) /= istdlak ) then ! This calculation is now done for lakes in initLake.
             do j = 1,nlevs
                this%h2osoi_vol(c,j) = this%h2osoi_liq(c,j)/(col_pp%dz(c,j)*denh2o) &
                                         + this%h2osoi_ice(c,j)/(col_pp%dz(c,j)*denice)
             end do
          end if
       end do
    end if

    ! If initial run -- ensure that water is properly bounded (read only)
    if (flag == 'read' ) then
       if ( is_first_step() .and. bound_h2osoi) then
          do c = bounds%begc, bounds%endc
             l = col_pp%landunit(c)
             if ( col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall .or. &
                  col_pp%itype(c) == icol_roof )then
                nlevs = nlevurb
             else
                nlevs = nlevgrnd
             end if
             do j = 1,nlevs
                l = col_pp%landunit(c)
                if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
                   this%h2osoi_liq(c,j) = max(0._r8,this%h2osoi_liq(c,j))
                   this%h2osoi_ice(c,j) = max(0._r8,this%h2osoi_ice(c,j))
                   this%h2osoi_vol(c,j) = this%h2osoi_liq(c,j)/(col_pp%dz(c,j)*denh2o) &
                                       + this%h2osoi_ice(c,j)/(col_pp%dz(c,j)*denice)
                   if (j == 1) then
                      maxwatsat = (watsat_input(c,j)*col_pp%dz(c,j)*1000.0_r8 + pondmx) / (col_pp%dz(c,j)*1000.0_r8)
                   else
                      maxwatsat =  watsat_input(c,j)
                   end if
                   if (this%h2osoi_vol(c,j) > maxwatsat) then 
                      excess = (this%h2osoi_vol(c,j) - maxwatsat)*col_pp%dz(c,j)*1000.0_r8
                      totwat = this%h2osoi_liq(c,j) + this%h2osoi_ice(c,j)
                      this%h2osoi_liq(c,j) = this%h2osoi_liq(c,j) - &
                                           (this%h2osoi_liq(c,j)/totwat) * excess
                      this%h2osoi_ice(c,j) = this%h2osoi_ice(c,j) - &
                                           (this%h2osoi_ice(c,j)/totwat) * excess
                   end if
                   this%h2osoi_liq(c,j) = max(watmin,this%h2osoi_liq(c,j))
                   this%h2osoi_ice(c,j) = max(watmin,this%h2osoi_ice(c,j))
                   this%h2osoi_vol(c,j) = this%h2osoi_liq(c,j)/(col_pp%dz(c,j)*denh2o) &
                                             + this%h2osoi_ice(c,j)/(col_pp%dz(c,j)*denice)
                end if
             end do
          end do
       end if

    endif   ! end if if-read flag

    end subroutine col_ws_restart


  !------------------------------------------------------------------------
  subroutine col_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_ws_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon state data structure
  !------------------------------------------------------------------------
  subroutine init_col_cs(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_cs
    
  !------------------------------------------------------------------------
  subroutine clean_col_cs(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_cs
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ns(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ns
    
  !------------------------------------------------------------------------
  subroutine clean_col_ns(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ns
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ps(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ps
    
  !------------------------------------------------------------------------
  subroutine clean_col_ps(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ps

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy flux data structure
  !------------------------------------------------------------------------
  subroutine col_ef_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    integer, intent(in) :: begc,endc
    ! !LOCAL VARIABLES:
    integer  :: l,c
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_h2osfc_to_snow  (begc:endc))              ; this%eflx_h2osfc_to_snow  (:)   = nan
    allocate(this%eflx_snomelt         (begc:endc))              ; this%eflx_snomelt         (:)   = nan
    allocate(this%eflx_snomelt_r       (begc:endc))              ; this%eflx_snomelt_r       (:)   = nan
    allocate(this%eflx_snomelt_u       (begc:endc))              ; this%eflx_snomelt_u       (:)   = nan
    allocate(this%eflx_bot             (begc:endc))              ; this%eflx_bot             (:)   = nan
    allocate(this%eflx_fgr12           (begc:endc))              ; this%eflx_fgr12           (:)   = nan
    allocate(this%eflx_fgr             (begc:endc, 1:nlevgrnd))  ; this%eflx_fgr             (:,:) = nan
    allocate(this%eflx_building_heat   (begc:endc))              ; this%eflx_building_heat   (:)   = nan
    allocate(this%eflx_urban_ac        (begc:endc))              ; this%eflx_urban_ac        (:)   = nan
    allocate(this%eflx_urban_heat      (begc:endc))              ; this%eflx_urban_heat      (:)   = nan
    allocate(this%eflx_hs_h2osfc       (begc:endc))              ; this%eflx_hs_h2osfc       (:)   = nan
    allocate(this%eflx_hs_top_snow     (begc:endc))              ; this%eflx_hs_top_snow     (:)   = nan
    allocate(this%eflx_hs_soil         (begc:endc))              ; this%eflx_hs_soil         (:)   = nan
    allocate(this%eflx_sabg_lyr        (begc:endc, -nlevsno+1:1)); this%eflx_sabg_lyr        (:,:) = nan
    allocate(this%eflx_dhsdT           (begc:endc))              ; this%eflx_dhsdT           (:)   = nan
    allocate(this%htvp                 (begc:endc))              ; this%htvp                 (:)   = nan
    allocate(this%xmf                  (begc:endc))              ; this%xmf                  (:)   = nan
    allocate(this%xmf_h2osfc           (begc:endc))              ; this%xmf_h2osfc           (:)   = nan
    allocate(this%imelt                (begc:endc,-nlevsno+1:nlevgrnd))  ; this%imelt        (:,:) = huge(1)
    allocate(this%eflx_soil_grnd       (begc:endc))              ; this%eflx_soil_grnd       (:)   = nan
    allocate(this%eflx_rnet_soil       (begc:endc))              ; this%eflx_rnet_soil       (:)   = nan
    allocate(this%eflx_fgr0_soil       (begc:endc))              ; this%eflx_fgr0_soil       (:)   = nan
    allocate(this%eflx_fgr0_snow       (begc:endc))              ; this%eflx_fgr0_snow       (:)   = nan
    allocate(this%eflx_fgr0_h2osfc     (begc:endc))              ; this%eflx_fgr0_h2osfc     (:)   = nan
    allocate(this%errsoi               (begc:endc))              ; this%errsoi               (:)   = nan
    allocate(this%errseb               (begc:endc))              ; this%errseb               (:)   = nan
    allocate(this%errsol               (begc:endc))              ; this%errsol               (:)   = nan
    allocate(this%errlon               (begc:endc))              ; this%errlon               (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of col_ef
    !-----------------------------------------------------------------------
    this%eflx_snomelt(begc:endc) = spval
    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=this%eflx_snomelt, c2l_scale_type='urbanf')

    this%eflx_snomelt_r(begc:endc) = spval
    call hist_addfld1d (fname='FSM_R',  units='W/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=this%eflx_snomelt_r, set_spec=spval)

    this%eflx_snomelt_u(begc:endc) = spval
    call hist_addfld1d (fname='FSM_U',  units='W/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=this%eflx_snomelt_u, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_building_heat(begc:endc) = spval
    call hist_addfld1d (fname='BUILDHEAT', units='W/m^2',  &
         avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
         ptr_col=this%eflx_building_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_ac(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
         avgflag='A', long_name='urban air conditioning flux', &
         ptr_col=this%eflx_urban_ac, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_heat(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
         avgflag='A', long_name='urban heating flux', &
         ptr_col=this%eflx_urban_heat, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_fgr12(begc:endc) = spval
    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=this%eflx_fgr12, set_lake=spval)

    this%eflx_fgr(begc:endc,:) = spval
    call hist_addfld2d (fname='FGR_SOIL_R', units='watt/m^2', type2d='levgrnd', &
         avgflag='A', long_name='Rural downward heat flux at interface below each soil layer', &
         ptr_col=this%eflx_fgr, set_spec=spval, default='inactive')

    this%errsoi(begc:endc) = spval
    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=this%errsoi)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of col_ef
    !-----------------------------------------------------------------------
    do c = begc, endc
       l = col_pp%landunit(c)
       if (lun_pp%urbpoi(l)) then
          this%eflx_building_heat(c) = 0._r8
          this%eflx_urban_ac(c)      = 0._r8
          this%eflx_urban_heat(c)    = 0._r8
       else
          this%eflx_building_heat(c) = 0._r8
          this%eflx_urban_ac(c)      = 0._r8
          this%eflx_urban_heat(c)    = 0._r8
       end if
    end do
    
  end subroutine col_ef_init
    
  !------------------------------------------------------------------------
  subroutine col_ef_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write column energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC', xtype=ncd_double,  dim1name='column', &
         long_name='urban air conditioning flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_ac)

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT', xtype=ncd_double, dim1name='column', &
         long_name='urban heating flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_heat)

  end subroutine col_ef_restart
  
  !------------------------------------------------------------------------
  subroutine col_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_ef_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water flux data structure
  !------------------------------------------------------------------------
  subroutine col_wf_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    integer, intent(in) :: begc,endc
    ! !LOCAL VARIABLES:
    integer  :: l,c
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of col_wf
    !-----------------------------------------------------------------------
    allocate(this%qflx_prec_grnd         (begc:endc))             ; this%qflx_prec_grnd       (:)   = nan
    allocate(this%qflx_rain_grnd         (begc:endc))             ; this%qflx_rain_grnd       (:)   = nan
    allocate(this%qflx_snow_grnd         (begc:endc))             ; this%qflx_snow_grnd       (:)   = nan
    allocate(this%qflx_sub_snow          (begc:endc))             ; this%qflx_sub_snow        (:)   = nan
    allocate(this%qflx_sub_snow_vol      (begc:endc))             ; this%qflx_sub_snow_vol    (:)   = nan
    allocate(this%qflx_evap_soi          (begc:endc))             ; this%qflx_evap_soi        (:)   = nan
    allocate(this%qflx_evap_veg          (begc:endc))             ; this%qflx_evap_veg        (:)   = nan
    allocate(this%qflx_evap_can          (begc:endc))             ; this%qflx_evap_can        (:)   = nan
    allocate(this%qflx_evap_tot          (begc:endc))             ; this%qflx_evap_tot        (:)   = nan
    allocate(this%qflx_evap_grnd         (begc:endc))             ; this%qflx_evap_grnd       (:)   = nan
    allocate(this%qflx_snwcp_liq         (begc:endc))             ; this%qflx_snwcp_liq       (:)   = nan
    allocate(this%qflx_snwcp_ice         (begc:endc))             ; this%qflx_snwcp_ice       (:)   = nan
    allocate(this%qflx_tran_veg          (begc:endc))             ; this%qflx_tran_veg        (:)   = nan
    allocate(this%qflx_dew_snow          (begc:endc))             ; this%qflx_dew_snow        (:)   = nan
    allocate(this%qflx_dew_grnd          (begc:endc))             ; this%qflx_dew_grnd        (:)   = nan
    allocate(this%qflx_prec_intr         (begc:endc))             ; this%qflx_prec_intr       (:)   = nan
    allocate(this%qflx_ev_snow           (begc:endc))             ; this%qflx_ev_snow         (:)   = nan
    allocate(this%qflx_ev_soil           (begc:endc))             ; this%qflx_ev_soil         (:)   = nan
    allocate(this%qflx_ev_h2osfc         (begc:endc))             ; this%qflx_ev_h2osfc       (:)   = nan
    allocate(this%qflx_gross_evap_soil   (begc:endc))             ; this%qflx_gross_evap_soil (:)   = nan
    allocate(this%qflx_gross_infl_soil   (begc:endc))             ; this%qflx_gross_infl_soil (:)   = nan
    allocate(this%qflx_adv               (begc:endc,0:nlevgrnd))  ; this%qflx_adv             (:,:) = nan
    allocate(this%qflx_rootsoi           (begc:endc,1:nlevgrnd))  ; this%qflx_rootsoi         (:,:) = nan
    allocate(this%dwb                    (begc:endc))             ; this%dwb                  (:)   = nan
    allocate(this%qflx_infl              (begc:endc))             ; this%qflx_infl            (:)   = nan
    allocate(this%qflx_surf              (begc:endc))             ; this%qflx_surf            (:)   = nan
    allocate(this%qflx_drain             (begc:endc))             ; this%qflx_drain           (:)   = nan
    allocate(this%qflx_totdrain          (begc:endc))             ; this%qflx_totdrain        (:)   = nan
    allocate(this%qflx_top_soil          (begc:endc))             ; this%qflx_top_soil        (:)   = nan
    allocate(this%qflx_h2osfc_to_ice     (begc:endc))             ; this%qflx_h2osfc_to_ice   (:)   = nan
    allocate(this%qflx_h2osfc_surf       (begc:endc))             ; this%qflx_h2osfc_surf     (:)   = nan
    allocate(this%qflx_snow_h2osfc       (begc:endc))             ; this%qflx_snow_h2osfc     (:)   = nan
    allocate(this%qflx_drain_perched     (begc:endc))             ; this%qflx_drain_perched   (:)   = nan
    allocate(this%qflx_deficit           (begc:endc))             ; this%qflx_deficit         (:)   = nan
    allocate(this%qflx_floodc            (begc:endc))             ; this%qflx_floodc          (:)   = nan
    allocate(this%qflx_sl_top_soil       (begc:endc))             ; this%qflx_sl_top_soil     (:)   = nan
    allocate(this%qflx_snomelt           (begc:endc))             ; this%qflx_snomelt         (:)   = nan
    allocate(this%qflx_snow_melt         (begc:endc))             ; this%qflx_snow_melt       (:)   = nan
    allocate(this%qflx_qrgwl             (begc:endc))             ; this%qflx_qrgwl           (:)   = nan
    allocate(this%qflx_runoff            (begc:endc))             ; this%qflx_runoff          (:)   = nan
    allocate(this%qflx_runoff_r          (begc:endc))             ; this%qflx_runoff_r        (:)   = nan
    allocate(this%qflx_runoff_u          (begc:endc))             ; this%qflx_runoff_u        (:)   = nan
    allocate(this%qflx_rsub_sat          (begc:endc))             ; this%qflx_rsub_sat        (:)   = nan
    allocate(this%qflx_snofrz_lyr        (begc:endc,-nlevsno+1:0)); this%qflx_snofrz_lyr      (:,:) = nan
    allocate(this%qflx_snofrz            (begc:endc))             ; this%qflx_snofrz          (:)   = nan
    allocate(this%qflx_glcice            (begc:endc))             ; this%qflx_glcice          (:)   = nan
    allocate(this%qflx_glcice_frz        (begc:endc))             ; this%qflx_glcice_frz      (:)   = nan
    allocate(this%qflx_glcice_melt       (begc:endc))             ; this%qflx_glcice_melt     (:)   = nan
    allocate(this%qflx_drain_vr          (begc:endc,1:nlevgrnd))  ; this%qflx_drain_vr        (:,:) = nan
    allocate(this%qflx_h2osfc2topsoi     (begc:endc))             ; this%qflx_h2osfc2topsoi   (:)   = nan
    allocate(this%qflx_snow2topsoi       (begc:endc))             ; this%qflx_snow2topsoi     (:)   = nan
    allocate(this%qflx_lateral           (begc:endc))             ; this%qflx_lateral         (:)   = nan
    allocate(this%snow_sources           (begc:endc))             ; this%snow_sources         (:)   = nan
    allocate(this%snow_sinks             (begc:endc))             ; this%snow_sinks           (:)   = nan
    allocate(this%qflx_irrig             (begc:endc))             ; this%qflx_irrig           (:)   = nan
    
    !-----------------------------------------------------------------------
    ! initialize history fields for select members of col_wf
    !-----------------------------------------------------------------------
    this%qflx_snwcp_ice(begc:endc) = spval
    call hist_addfld1d (fname='QSNWCPICE_NODYNLNDUSE', units='mm H2O/s', &
         avgflag='A', long_name='excess snowfall due to snow capping not including correction for land use change', &
         ptr_col=this%qflx_snwcp_ice, c2l_scale_type='urbanf')

    this%dwb(begc:endc) = spval
    call hist_addfld1d (fname='DWB',  units='mm/s',  &
         avgflag='A', long_name='net change in total water mass', &
         ptr_col=this%dwb, c2l_scale_type='urbanf')

    this%qflx_infl(begc:endc) = spval
    call hist_addfld1d (fname='QINFL',  units='mm/s',  &
         avgflag='A', long_name='infiltration', &
         ptr_col=this%qflx_infl, c2l_scale_type='urbanf')

    this%qflx_surf(begc:endc) = spval
    call hist_addfld1d (fname='QOVER',  units='mm/s',  &
         avgflag='A', long_name='surface runoff', &
         ptr_col=this%qflx_surf, c2l_scale_type='urbanf')

    this%qflx_drain(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI',  units='mm/s',  &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_col=this%qflx_drain, c2l_scale_type='urbanf')

    this%qflx_top_soil(begc:endc) = spval
    call hist_addfld1d (fname='QTOPSOIL',  units='mm/s',  &
         avgflag='A', long_name='water input to surface', &
         ptr_col=this%qflx_top_soil, c2l_scale_type='urbanf', default='inactive')

    this%qflx_h2osfc_surf(begc:endc) = spval
    call hist_addfld1d (fname='QH2OSFC',  units='mm/s',  &
         avgflag='A', long_name='surface water runoff', &
         ptr_col=this%qflx_h2osfc_surf)

    this%qflx_drain_perched(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI_PERCH',  units='mm/s',  &
         avgflag='A', long_name='perched wt drainage', &
         ptr_col=this%qflx_drain_perched, c2l_scale_type='urbanf')

    this%qflx_snow_melt(begc:endc) = spval
    call hist_addfld1d (fname='QSNOMELT',  units='mm/s',  &
         avgflag='A', long_name='snow melt', &
         ptr_col=this%qflx_snow_melt, c2l_scale_type='urbanf')

    this%qflx_qrgwl(begc:endc) = spval
    call hist_addfld1d (fname='QRGWL',  units='mm/s',  &
         avgflag='A', long_name='surface runoff at glaciers (liquid only), wetlands, lakes', &
         ptr_col=this%qflx_qrgwl, c2l_scale_type='urbanf')

    this%qflx_runoff(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_NODYNLNDUSE',  units='mm/s',  &
         avgflag='A', &
         long_name='total liquid runoff (does not include QSNWCPICE) not including correction for land use change', &
         ptr_col=this%qflx_runoff, c2l_scale_type='urbanf')

    this%qflx_runoff_r(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_R', units='mm/s',  &
         avgflag='A', long_name='Rural total runoff', &
         ptr_col=this%qflx_runoff_r, set_spec=spval)

    this%qflx_runoff_u(begc:endc) = spval
    call hist_addfld1d (fname='QRUNOFF_U', units='mm/s',  &
         avgflag='A', long_name='Urban total runoff', &
         ptr_col=this%qflx_runoff_u, set_nourb=spval, c2l_scale_type='urbanf')

    this%qflx_rsub_sat(begc:endc) = spval
    call hist_addfld1d (fname='QDRAI_XS',  units='mm/s',  &
         avgflag='A', long_name='saturation excess drainage', &
         ptr_col=this%qflx_rsub_sat, c2l_scale_type='urbanf')

    this%qflx_snofrz(begc:endc) = spval
    call hist_addfld1d (fname='QSNOFRZ', units='kg/m2/s', &
         avgflag='A', long_name='column-integrated snow freezing rate', &
         ptr_col=this%qflx_snofrz, set_lake=spval, c2l_scale_type='urbanf', default='inactive')
         
    if (create_glacier_mec_landunit) then
       this%qflx_glcice(begc:endc) = spval
       call hist_addfld1d (fname='QICE',  units='mm/s',  &
            avgflag='A', long_name='ice growth/melt', &
            ptr_col=this%qflx_glcice, l2g_scale_type='ice')
    end if

    if (create_glacier_mec_landunit) then
       this%qflx_glcice_frz(begc:endc) = spval
       call hist_addfld1d (fname='QICE_FRZ',  units='mm/s',  &
            avgflag='A', long_name='ice growth', &
            ptr_col=this%qflx_glcice_frz, l2g_scale_type='ice')
    end if

    if (create_glacier_mec_landunit) then
       this%qflx_glcice_melt(begc:endc) = spval
       call hist_addfld1d (fname='QICE_MELT',  units='mm/s',  &
            avgflag='A', long_name='ice melt', &
            ptr_col=this%qflx_glcice_melt, l2g_scale_type='ice')
    end if

    ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at any
    ! given time step but only if there is at least one snow layer (for all landunits 
    ! except lakes).  Also note that monthly average files of snow_sources and snow_sinks
    ! sinks must be weighted by number of days in the month to diagnose, for example, an 
    ! annual value of the change in h2osno. 
    this%snow_sources(begc:endc) = spval
    call hist_addfld1d (fname='SNOW_SOURCES',  units='mm/s',  &
         avgflag='A', long_name='snow sources (liquid water)', &
         ptr_col=this%snow_sources, c2l_scale_type='urbanf')

    this%snow_sinks(begc:endc) = spval
    call hist_addfld1d (fname='SNOW_SINKS',  units='mm/s',  &
         avgflag='A', long_name='snow sinks (liquid water)', &
         ptr_col=this%snow_sinks, c2l_scale_type='urbanf')

    !-----------------------------------------------------------------------
    ! set cold-start initial values for select members of col_wf
    !-----------------------------------------------------------------------
    this%qflx_evap_grnd(begc:endc) = 0.0_r8
    this%qflx_dew_grnd (begc:endc) = 0.0_r8
    this%qflx_dew_snow (begc:endc) = 0.0_r8

    this%qflx_h2osfc_surf(begc:endc) = 0._r8
    this%qflx_snow_melt  (begc:endc)   = 0._r8

    this%dwb(begc:endc) = 0._r8
    ! needed for CNNLeaching 
    do c = begc, endc
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%qflx_drain(c) = 0._r8
          this%qflx_surf(c)  = 0._r8
       end if
    end do

  end subroutine col_wf_init
    
  !------------------------------------------------------------------------
  subroutine col_wf_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write column water state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='qflx_snofrz_lyr', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer ice freezing rate', units='kg m-2 s-1', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snofrz_lyr)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snofrz_lyr to zero
       this%qflx_snofrz_lyr(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    call restartvar(ncid=ncid, flag=flag, varname='qflx_snow_melt', xtype=ncd_double,  &
         dim1name='column', &
         long_name='net snow melt', units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%qflx_snow_melt)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snow_melt to zero
       this%qflx_snow_melt(bounds%begc:bounds%endc) = 0._r8
    endif

  end subroutine col_wf_restart
  
  !------------------------------------------------------------------------
  subroutine col_wf_clean(this)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine col_wf_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_cf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_cf
    
  !------------------------------------------------------------------------
  subroutine clean_col_cf(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_cf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_nf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_nf
    
  !------------------------------------------------------------------------
  subroutine clean_col_nf(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_nf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_pf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_pf
    
  !------------------------------------------------------------------------
  subroutine clean_col_pf(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_pf
  
    !------------------------------------------------------------------------
    
end module ColumnDataType

