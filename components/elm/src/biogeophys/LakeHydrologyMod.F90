module LakeHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of Lake Hydrology. Full hydrology, aerosol deposition, etc. of snow layers is
  ! done. However, there is no infiltration, and the water budget is balanced with 
  ! qflx_qrgwl. Lake water mass is kept constant. The soil is simply maintained at
  ! volumetric saturation if ice melting frees up pore space. Likewise, if the water
  ! portion alone at some point exceeds pore capacity, it is reduced. This is consistent
  ! with the possibility of initializing the soil layer with excess ice.
  ! 
  ! If snow layers are present over an unfrozen lake, and the top layer of the lake
  ! is capable of absorbing the latent heat without going below freezing, 
  ! the snow-water is runoff and the latent heat is subtracted from the lake.
  !
  ! Minimum snow layer thickness for lakes has been increased to avoid instabilities with 30 min timestep.
  ! Also frost / dew is prevented from being added to top snow layers that have already melted during the phase change step.
  !
  ! ! USES
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use decompMod            , only : bounds_type
  use TopounitDataType     , only : top_as, top_af    ! atmospheric state and flux variables
  use ColumnType           , only : col_pp
  use ColumnDataType       , only : col_es, col_ef, col_ws, col_wf  
  use VegetationType       , only : veg_pp                
  use VegetationDataType   , only : veg_ef, veg_wf
  use atm2lndType          , only : atm2lnd_type
  use AerosolType          , only : aerosol_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use LakeStateType        , only : lakestate_type
  use SoilStateType        , only : soilstate_type
  use TemperatureType      , only : temperature_type
  use WaterfluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use elm_varcon           , only : snw_rds_min  
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LakeHydrology              ! Calculates soil/snow hydrology
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LakeHydrology(bounds, &
       num_lakec, filter_lakec, num_lakep, filter_lakep, &
       num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc, &
       atm2lnd_vars, temperature_vars, soilstate_vars, waterstate_vars, waterflux_vars, &
       energyflux_vars, aerosol_vars, lakestate_vars)
    !
    ! !DESCRIPTION:
    ! WARNING: This subroutine assumes lake columns have one and only one pft.
    !
    ! Sequence is:
    !  LakeHydrology:
    !    Do needed tasks from CanopyHydrology, Biogeophysics2, & top of SoilHydrology.
    !    -> SnowWater:             change of snow mass and snow water onto soil
    !    -> SnowCompaction:        compaction of snow layers
    !    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
    !    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
    !
    !    Add water to soil if melting has left it with open pore space.
    !    If snow layers are found above a lake with unfrozen top layer, whose top
    !    layer has enough heat to melt all the snow ice without freezing, do so
    !    and eliminate the snow layers.
    !    Cleanup and do water balance.
    !
    ! !USES:
    use elm_varcon      , only : denh2o, denice, spval, hfus, tfrz, cpliq, cpice
    use elm_varpar      , only : nlevsno, nlevgrnd, nlevsoi
    use clm_varctl      , only : iulog
    use clm_time_manager, only : get_step_size
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, SnowWater, BuildSnowFilter
    use SnowHydrologyMod, only : DivideSnowLayers
    use LakeCon         , only : lsadz
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_lakec               ! number of column lake points in column filter
    integer                , intent(in)    :: filter_lakec(:)         ! column filter for lake points
    integer                , intent(in)    :: num_lakep               ! number of pft lake points in column filter
    integer                , intent(in)    :: filter_lakep(:)         ! patch filter for lake points
    integer                , intent(out)   :: num_shlakesnowc         ! number of column snow points
    integer                , intent(out)   :: filter_shlakesnowc(:)   ! column filter for snow points
    integer                , intent(out)   :: num_shlakenosnowc       ! number of column non-snow points
    integer                , intent(out)   :: filter_shlakenosnowc(:) ! column filter for non-snow points
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(temperature_type) , intent(inout) :: temperature_vars
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(waterstate_type)  , intent(inout) :: waterstate_vars
    type(waterflux_type)   , intent(inout) :: waterflux_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(aerosol_type)     , intent(inout) :: aerosol_vars
    type(lakestate_type)   , intent(inout) :: lakestate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: p,fp,g,t,l,c,j,fc,jtop                          ! indices
    real(r8) :: dtime                                           ! land model time step (sec)
    integer  :: newnode                                         ! flag when new snow node is set, (1=yes, 0=no)
    real(r8) :: dz_snowf                                        ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: bifall                                          ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: fracsnow(bounds%begp:bounds%endp)               ! frac of precipitation that is snow
    real(r8) :: fracrain(bounds%begp:bounds%endp)               ! frac of precipitation that is rain
    real(r8) :: qflx_prec_grnd_snow(bounds%begp:bounds%endp)    ! snow precipitation incident on ground [mm/s]
    real(r8) :: qflx_prec_grnd_rain(bounds%begp:bounds%endp)    ! rain precipitation incident on ground [mm/s]
    real(r8) :: qflx_evap_soi_lim                               ! temporary evap_soi limited by top snow layer content [mm/s]
    real(r8) :: h2osno_temp                                     ! temporary h2osno [kg/m^2]
    real(r8) :: sumsnowice(bounds%begc:bounds%endc)             ! sum of snow ice if snow layers found above unfrozen lake [kg/m&2]
    logical  :: unfrozen(bounds%begc:bounds%endc)               ! true if top lake layer is unfrozen with snow layers above
    real(r8) :: heatrem                                         ! used in case above [J/m^2]
    real(r8) :: heatsum(bounds%begc:bounds%endc)                ! used in case above [J/m^2]
    real(r8) :: snowmass                                        ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct                                 ! temporary factor used to correct for snow capping
    real(r8), parameter :: snow_bd = 250._r8                    ! assumed snow bulk density (for lakes w/out resolved snow layers) [kg/m^3]
                                                                ! Should only be used for frost below.
    !-----------------------------------------------------------------------

    associate(                                                            & 
         pcolumn              =>  veg_pp%column                            , & ! Input:  [integer  (:)   ]  pft's column index                       
         ptopounit            =>  veg_pp%topounit                          , & ! Input:  [integer  (:)   ]  pft's topounit index                     
         pgridcell            =>  veg_pp%gridcell                          , & ! Input:  [integer  (:)   ]  pft's gridcell index                     
         cgridcell            =>  col_pp%gridcell                          , & ! Input:  [integer  (:)   ]  column's gridcell                        
         clandunit            =>  col_pp%landunit                          , & ! Input:  [integer  (:)   ]  column's landunit                        
         dz_lake              =>  col_pp%dz_lake                           , & ! Input:  [real(r8) (:,:) ]  layer thickness for lake (m)          
         z                    =>  col_pp%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth  (m)                      
         dz                   =>  col_pp%dz                                , & ! Input:  [real(r8) (:,:) ]  layer thickness depth (m)             
         zi                   =>  col_pp%zi                                , & ! Input:  [real(r8) (:,:) ]  interface depth (m)                   
         snl                  =>  col_pp%snl                               , & ! Input:  [integer  (:)   ]  number of snow layers                    

         forc_rain            =>  top_af%rain                           , & ! Input:  [real(r8) (:)   ]  rain rate (kg H2O/m**2/s, or mm liquid H2O/s)                        
         forc_snow            =>  top_af%snow                           , & ! Input:  [real(r8) (:)   ]  snow rate (kg H2O/m**2/s, or mm liquid H2O/s)                        
         forc_t               =>  top_as%tbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)        
         qflx_floodg          =>  atm2lnd_vars%forc_flood_grc           , & ! Input:  [real(r8) (:)   ]  gridcell flux of flood water from RTM   

         watsat               =>  soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)

         t_lake               =>  col_es%t_lake           , & ! Input:  [real(r8) (:,:) ]  lake temperature (Kelvin)             
         t_grnd               =>  col_es%t_grnd           , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)             
         t_soisno             =>  col_es%t_soisno         , & ! Output: [real(r8) (:,:) ]  snow temperature (Kelvin)             
         dTdz_top             =>  col_es%dTdz_top         , & ! Output: [real(r8) (:)   ]  temperature gradient in top layer K m-1] !TOD 
         snot_top             =>  col_es%snot_top         , & ! Output: [real(r8) (:)   ]  snow temperature in top layer [K]  !TODO

         do_capsnow           =>  col_ws%do_capsnow        , & ! Input:  [logical  (:)   ]  true => do snow capping                  
         begwb                =>  col_ws%begwb             , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step    
         endwb                =>  col_ws%endwb             , & ! Output: [real(r8) (:)   ]  water mass end of the time step         
         h2osoi_liq_depth_intg=>  col_ws%h2osoi_liq_depth_intg, & ! Output: [real(r8) (:)   ]  grid-level depth integrated liquid soil water
         h2osoi_ice_depth_intg=>  col_ws%h2osoi_ice_depth_intg, & ! Output: [real(r8) (:)   ]  grid-level depth integrated ice soil water
         snw_rds              =>  col_ws%snw_rds           , & ! Output: [real(r8) (:,:) ]  effective snow grain radius (col,lyr) [microns, m^-6] 
         snw_rds_top          =>  col_ws%snw_rds_top       , & ! Output: [real(r8) (:)   ]  effective snow grain size, top layer [microns] 
         h2osno_top           =>  col_ws%h2osno_top        , & ! Output: [real(r8) (:)   ]  mass of snow in top layer [kg]    
         sno_liq_top          =>  col_ws%sno_liq_top       , & ! Output: [real(r8) (:)   ]  liquid water fraction in top snow layer [frc] 
         frac_sno_eff         =>  col_ws%frac_sno_eff      , & ! Output: [real(r8) (:)   ]  needed for snicar code                  
         frac_iceold          =>  col_ws%frac_iceold       , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water
         snow_depth           =>  col_ws%snow_depth        , & ! Output: [real(r8) (:)   ]  snow height (m)                         
         h2osno               =>  col_ws%h2osno            , & ! Output: [real(r8) (:)   ]  snow water (mm H2O)                     
         snowice              =>  col_ws%snowice           , & ! Output: [real(r8) (:)   ]  average snow ice lens                   
         snowliq              =>  col_ws%snowliq           , & ! Output: [real(r8) (:)   ]  average snow liquid water               
         h2osoi_ice           =>  col_ws%h2osoi_ice        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq           =>  col_ws%h2osoi_liq        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                  
         h2osoi_vol           =>  col_ws%h2osoi_vol        , & ! Output: [real(r8) (:,:) ]  volumetric soil water [m3/m3]         

         qflx_floodc          =>  col_wf%qflx_floodc        , & ! Output: [real(r8) (:)   ]  column flux of flood water from RTM     
         qflx_prec_grnd       =>  veg_wf%qflx_prec_grnd   , & ! Output: [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
         qflx_snow_grnd_patch =>  veg_wf%qflx_snow_grnd   , & ! Output: [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
         qflx_rain_grnd       =>  veg_wf%qflx_rain_grnd   , & ! Output: [real(r8) (:)   ]  rain on ground after interception (mm H2O/s) [+]
         qflx_rain_grnd_col   =>  col_wf%qflx_rain_grnd     , & ! Output: [real(r8) (:)   ]  rain on ground after interception (mm H2O/s) [+]
         qflx_evap_tot        =>  veg_wf%qflx_evap_tot    , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
         qflx_evap_soi        =>  veg_wf%qflx_evap_soi    , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_sub_snow        =>  veg_wf%qflx_sub_snow    , & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]
         qflx_evap_grnd       =>  veg_wf%qflx_evap_grnd   , & ! Output: [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
         qflx_dew_snow        =>  veg_wf%qflx_dew_snow    , & ! Output: [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
         qflx_dew_grnd        =>  veg_wf%qflx_dew_grnd    , & ! Output: [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_snwcp_ice       =>  veg_wf%qflx_snwcp_ice   , & ! Output: [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]
         qflx_snwcp_liq       =>  veg_wf%qflx_snwcp_liq   , & ! Output: [real(r8) (:)   ]  excess rainfall due to snow capping (mm H2O /s) [+]
         qflx_dirct_rain      =>  veg_wf%qflx_dirct_rain  , & ! Output: [real(r8) (:)   ]  direct rain throughfall (mm H2O/s)
         qflx_leafdrip        =>  veg_wf%qflx_leafdrip    , & ! Output: [real(r8) (:)   ]  leaf rain drip (mm H2O/s)
         qflx_snomelt         =>  col_wf%qflx_snomelt       , & ! Output: [real(r8) (:)   ]  snow melt (mm H2O /s)                   
         qflx_dirct_rain_col  =>  col_wf%qflx_dirct_rain    , & ! Output: [real(r8) (:)   ]  direct rain throughfall (mm H2O/s)
         qflx_leafdrip_col    =>  col_wf%qflx_leafdrip      , & ! Output: [real(r8) (:)   ]  leaf rain drip (mm H2O/s)
         qflx_prec_grnd_col   =>  col_wf%qflx_prec_grnd     , & ! Output: [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
         qflx_evap_grnd_col   =>  col_wf%qflx_evap_grnd     , & ! Output: [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
         qflx_dew_grnd_col    =>  col_wf%qflx_dew_grnd      , & ! Output: [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_dew_snow_col    =>  col_wf%qflx_dew_snow      , & ! Output: [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
         qflx_sub_snow_col    =>  col_wf%qflx_sub_snow      , & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]
         qflx_snow_grnd_col   =>  col_wf%qflx_snow_grnd     , & ! Output: [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
         qflx_evap_tot_col    =>  col_wf%qflx_evap_tot      , & ! Output: [real(r8) (:)   ]  pft quantity averaged to the column (assuming one pft)
         qflx_snwcp_ice_col   =>  col_wf%qflx_snwcp_ice     , & ! Output: [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]
         qflx_snwcp_liq_col   =>  col_wf%qflx_snwcp_liq     , & ! Output: [real(r8) (:)   ]  excess rainfall due to snow capping (mm H2O /s) [+]
         qflx_drain_perched   =>  col_wf%qflx_drain_perched , & ! Output: [real(r8) (:)   ]  perched wt sub-surface runoff (mm H2O /s) !TODO - move this to somewhere else
         qflx_h2osfc_surf     =>  col_wf%qflx_h2osfc_surf   , & ! Output: [real(r8) (:)   ]  surface water runoff (mm H2O /s)        
         qflx_snow_melt       =>  col_wf%qflx_snow_melt     , & ! Output: [real(r8) (:)   ]  net snow melt                           
         qflx_rsub_sat        =>  col_wf%qflx_rsub_sat      , & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]        
         qflx_surf            =>  col_wf%qflx_surf          , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)              
         qflx_drain           =>  col_wf%qflx_drain         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
         qflx_infl            =>  col_wf%qflx_infl          , & ! Output: [real(r8) (:)   ]  infiltration (mm H2O /s)                
         qflx_qrgwl           =>  col_wf%qflx_qrgwl         , & ! Output: [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes  
         qflx_runoff          =>  col_wf%qflx_runoff        , & ! Output: [real(r8) (:)   ]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_irrig           =>  veg_wf%qflx_irrig_patch   , & ! Output: [real(r8) (:)   ]  irrigation flux (mm H2O /s)             
         qflx_irrig_col       =>  col_wf%qflx_irrig         , & ! Output: [real(r8) (:)   ]  irrigation flux (mm H2O /s)             
         qflx_top_soil        =>  col_wf%qflx_top_soil      , & ! Output: [real(r8) (:)   ]  net water input into soil from top (mm/s)
         qflx_sl_top_soil     =>  col_wf%qflx_sl_top_soil   , & ! Output: [real(r8) (:)   ]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)

         eflx_snomelt         =>  col_ef%eflx_snomelt      , & ! Output: [real(r8) (:)   ]  snow melt heat flux (W/m**2)
         eflx_sh_tot          =>  veg_ef%eflx_sh_tot     , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_grnd         =>  veg_ef%eflx_sh_grnd    , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_soil_grnd       =>  veg_ef%eflx_soil_grnd  , & ! Output: [real(r8) (:)   ]  heat flux into snow / lake (W/m**2) [+ = into soil]
         eflx_gnet            =>  veg_ef%eflx_gnet       , & ! Output: [reay(r8) (:)   ]  net heat flux into ground (W/m**2)      
         eflx_grnd_lake       =>  veg_ef%eflx_grnd_lake  , & ! Output: [real(r8) (:)   ]  net heat flux into lake / snow surface, excluding light transmission (W/m**2)

         lake_icefrac         =>  lakestate_vars%lake_icefrac_col       , & ! Output: [real(r8) (:,:) ]  mass fraction of lake layer that is frozen

         begc => bounds%begc, &
         endc => bounds%endc  &
         )

      ! Determine step size

      dtime = get_step_size()

      ! Add soil water to water balance.
      do j = 1, nlevgrnd
         do fc = 1, num_lakec
            c = filter_lakec(fc)
            begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
            h2osoi_liq_depth_intg(c) = h2osoi_liq_depth_intg(c) + h2osoi_liq(c,j)
            h2osoi_ice_depth_intg(c) = h2osoi_ice_depth_intg(c) + h2osoi_ice(c,j)
         end do
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do precipitation onto ground, etc., from CanopyHydrology

      do fp = 1, num_lakep
         p = filter_lakep(fp)
         c = pcolumn(p)
         t = ptopounit(p)

         qflx_prec_grnd_snow(p) = forc_snow(t)
         qflx_prec_grnd_rain(p) = forc_rain(t)
         qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

         qflx_dirct_rain(p) = 0._r8
         qflx_leafdrip(p) = 0._r8

         if (do_capsnow(c)) then
            qflx_snwcp_ice(p) = qflx_prec_grnd_snow(p)
            qflx_snwcp_liq(p) = qflx_prec_grnd_rain(p)
            qflx_snow_grnd_patch(p) = 0._r8
            qflx_rain_grnd(p) = 0._r8
         else
            qflx_snwcp_ice(p) = 0._r8
            qflx_snwcp_liq(p) = 0._r8
            qflx_snow_grnd_patch(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
            qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
         end if
         ! Assuming one PFT; needed for below
         qflx_snow_grnd_col(c) = qflx_snow_grnd_patch(p)
         qflx_rain_grnd_col(c) = qflx_rain_grnd(p)

      end do ! (end pft loop)

      ! Determine snow height and snow water

      do fc = 1, num_lakec
         c = filter_lakec(fc)
         t = col_pp%topounit(c)

         ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
         ! U.S.Department of Agriculture Forest Service, Project F,
         ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

         if (do_capsnow(c)) then
            dz_snowf = 0._r8
         else
            if (forc_t(t) > tfrz + 2._r8) then
               bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
            else if (forc_t(t) > tfrz - 15._r8) then
               bifall=50._r8 + 1.7_r8*(forc_t(t) - tfrz + 15._r8)**1.5_r8
            else
               bifall=50._r8
            end if
            dz_snowf = qflx_snow_grnd_col(c)/bifall
            snow_depth(c) = snow_depth(c) + dz_snowf*dtime
            h2osno(c) = h2osno(c) + qflx_snow_grnd_col(c)*dtime  ! snow water equivalent (mm)
         end if

         ! When the snow accumulation exceeds 40 mm, initialize snow layer
         ! Currently, the water temperature for the precipitation is simply set
         ! as the surface air temperature

         newnode = 0    ! flag for when snow node will be initialized
         if (snl(c) == 0 .and. qflx_snow_grnd_col(c) > 0.0_r8 .and. snow_depth(c) >= 0.01_r8 + lsadz) then
            newnode = 1
            snl(c) = -1
            dz(c,0) = snow_depth(c)                       ! meter
            z(c,0) = -0.5_r8*dz(c,0)
            zi(c,-1) = -dz(c,0)
            t_soisno(c,0) = min(tfrz, forc_t(t))      ! K
            h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
            h2osoi_liq(c,0) = 0._r8                   ! kg/m2
            frac_iceold(c,0) = 1._r8

             ! intitialize SNICAR variables for fresh snow:
             call aerosol_vars%Reset(column=c)
             ! call waterstate_vars%Reset(column=c)
             col_ws%snw_rds(c,0) = snw_rds_min

         end if

         ! The change of ice partial density of surface node due to precipitation.
         ! Only ice part of snowfall is added here, the liquid part will be added
         ! later.

         if (snl(c) < 0 .and. newnode == 0) then
            h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+dtime*qflx_snow_grnd_col(c)
            dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
         end if

      end do

      ! Calculate sublimation and dew, adapted from HydrologyLake and Biogeophysics2.

      do fp = 1,num_lakep
         p = filter_lakep(fp)
         c = pcolumn(p)
         jtop = snl(c)+1

         qflx_evap_grnd(p) = 0._r8
         qflx_sub_snow(p) = 0._r8
         qflx_dew_snow(p) = 0._r8
         qflx_dew_grnd(p) = 0._r8

         if (jtop <= 0) then ! snow layers
            j = jtop
            ! Assign ground evaporation to sublimation from soil ice or to dew
            ! on snow or ground

            if (qflx_evap_soi(p) >= 0._r8) then
               ! for evaporation partitioning between liquid evap and ice sublimation, 
               ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
               ! Since we're not limiting evap over lakes, but still can't remove more from top
               ! snow layer than there is there, create temp. limited evap_soi.
               qflx_evap_soi_lim = min(qflx_evap_soi(p), (h2osoi_liq(c,j)+h2osoi_ice(c,j))/dtime)
               if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0._r8) then
                  qflx_evap_grnd(p) = max(qflx_evap_soi_lim*(h2osoi_liq(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
               else
                  qflx_evap_grnd(p) = 0._r8
               end if
               qflx_sub_snow(p) = qflx_evap_soi_lim - qflx_evap_grnd(p)     
            else
               ! if (t_grnd(c) < tfrz) then
               ! Causes rare blowup when thin snow layer should completely melt and has a high temp after thermal physics,
               ! but then is not eliminated in SnowHydrology because of this added frost. Also see below removal of
               ! completely melted single snow layer.
               if (t_grnd(c) < tfrz .and. t_soisno(c,j) < tfrz) then
                  qflx_dew_snow(p) = abs(qflx_evap_soi(p))
                  ! If top layer is only snow layer, SnowHydrology won't eliminate it if dew is added.
               else if (j < 0 .or. (t_grnd(c) == tfrz .and. t_soisno(c,j) == tfrz)) then
                  qflx_dew_grnd(p) = abs(qflx_evap_soi(p))
               end if
            end if
            ! Update the pft-level qflx_snowcap
            ! This was moved in from Hydrology2 to keep all pft-level
            ! calculations out of Hydrology2
            if (do_capsnow(c)) then
               qflx_snwcp_ice(p) = qflx_snwcp_ice(p) + qflx_dew_snow(p) 
               qflx_snwcp_liq(p) = qflx_snwcp_liq(p) + qflx_dew_grnd(p)
            end if

         else ! No snow layers
            if (qflx_evap_soi(p) >= 0._r8) then
               ! Sublimation: do not allow for more sublimation than there is snow
               ! after melt.  Remaining surface evaporation used for infiltration.
               qflx_sub_snow(p) = min(qflx_evap_soi(p), h2osno(c)/dtime)
               qflx_evap_grnd(p) = qflx_evap_soi(p) - qflx_sub_snow(p)
            else
               if (t_grnd(c) < tfrz-0.1_r8) then
                  qflx_dew_snow(p) = abs(qflx_evap_soi(p))
               else
                  qflx_dew_grnd(p) = abs(qflx_evap_soi(p))
               end if
            end if

            ! Update snow pack for dew & sub.

            h2osno_temp = h2osno(c)
            if (do_capsnow(c)) then
               h2osno(c) = h2osno(c) - qflx_sub_snow(p)*dtime
               qflx_snwcp_ice(p) = qflx_snwcp_ice(p) + qflx_dew_snow(p) 
               qflx_snwcp_liq(p) = qflx_snwcp_liq(p) + qflx_dew_grnd(p)
            else
               h2osno(c) = h2osno(c) + (-qflx_sub_snow(p)+qflx_dew_snow(p))*dtime
            end if
            if (h2osno_temp > 0._r8) then
               snow_depth(c) = snow_depth(c) * h2osno(c) / h2osno_temp
            else
               snow_depth(c) = h2osno(c)/snow_bd !Assume a constant snow bulk density = 250.
            end if

            h2osno(c) = max(h2osno(c), 0._r8)
         end if

         qflx_snwcp_ice_col(c) = qflx_snwcp_ice(p)
         qflx_snwcp_liq_col(c) = qflx_snwcp_liq(p)


      end do

      ! patch averages must be done here -- BEFORE SNOW CALCULATIONS AS THEY USE IT.
      ! for output to history tape and other uses
      ! (note that pft2col is called before LakeHydrology, so we can't use that routine
      ! to do these column -> pft averages)
      do fp = 1,num_lakep
         p = filter_lakep(fp)
         c = pcolumn(p)

         qflx_evap_tot_col(c)  = qflx_evap_tot(p)
         qflx_prec_grnd_col(c) = qflx_prec_grnd(p)
         qflx_evap_grnd_col(c) = qflx_evap_grnd(p)
         qflx_dew_grnd_col(c)  = qflx_dew_grnd(p)
         qflx_dew_snow_col(c)  = qflx_dew_snow(p)
         qflx_sub_snow_col(c)  = qflx_sub_snow(p)
         qflx_dirct_rain_col(c) = qflx_dirct_rain(p)
         qflx_leafdrip_col(c) = qflx_leafdrip(p)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Determine initial snow/no-snow filters (will be modified possibly by
      ! routines CombineSnowLayers and DivideSnowLayers below)

      call BuildSnowFilter(bounds, num_lakec, filter_lakec, &
           num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

      ! specify snow fraction
      do fc = 1, num_lakec
         c = filter_lakec(fc)
         if (h2osno(c) > 0.0_r8) then
            frac_sno_eff(c)     = 1._r8 
         else
            frac_sno_eff(c)     = 0._r8 
         endif
      enddo

      ! Determine the change of snow mass and the snow water onto soil

      call SnowWater(bounds, &
           num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc, &
           atm2lnd_vars, waterflux_vars, waterstate_vars, aerosol_vars)

      ! Determine soil hydrology
      ! Here this consists only of making sure that soil is saturated even as it melts and
      ! pore space opens up. Conversely, if excess ice is melting and the liquid water exceeds the
      ! saturation value, then remove water.

      do j = 1,nlevsoi  !nlevgrnd
         ! changed to nlevsoi on 8/11/10 to make consistent with non-lake bedrock
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
            ! Could have changed during phase change! (Added 8/11/10)

            if (h2osoi_vol(c,j) < watsat(c,j)) then
               h2osoi_liq(c,j) = (watsat(c,j)*dz(c,j) - h2osoi_ice(c,j)/denice)*denh2o
               ! h2osoi_vol will be updated below, and this water addition will come from qflx_qrgwl
            else if (h2osoi_liq(c,j) > watsat(c,j)*denh2o*dz(c,j)) then
               h2osoi_liq(c,j) = watsat(c,j)*denh2o*dz(c,j)
               ! Another way to do this would be: if h2osoi_vol > watsat then remove min(h2osoi_liq,
               !(h2osoi_vol-watsat)*dz*denh2o) from h2osoi_liq.  The question is whether the excess ice
               ! melts first or last (or simultaneously) to the pore ice.  Because excess ice is often in chunks,
               ! requiring greater convergence of heat to melt, assume it melts last.
               ! This will also improve the initialization behavior or in an occasionally warm year, the excess ice
               ! won't start going away if a layer is briefly at freezing.

               ! Allow up to 10% excess ice over watsat in refreezing soil,
               ! e.g. heaving soil.  (As with > 10% excess ice modeling, and for the lake water,
               ! the thermal conductivity will be adjusted down to compensate for the fact that the nominal dz is smaller
               ! than the real soil volume.)  The current solution is consistent but perhaps unrealistic in real soils,
               ! where slow drainage may occur during freezing; drainage is only assumed to occur here when >10% excess
               ! ice melts. The latter is more likely to be permanent rather than seasonal anyway. Attempting to remove the
               ! ice volume after some has already frozen during the timestep would not conserve energy unless this were
               ! incorporated into the ice stream.

            end if

         end do
      end do
      !!!!!!!!!!

      ! Natural compaction and metamorphosis.

      call SnowCompaction(bounds, num_shlakesnowc, filter_shlakesnowc, &
           temperature_vars, waterstate_vars)

      ! Combine thin snow elements

      call CombineSnowLayers(bounds, num_shlakesnowc, filter_shlakesnowc, &
           aerosol_vars, temperature_vars, waterflux_vars, waterstate_vars)

      ! Divide thick snow elements

      call DivideSnowLayers(bounds, num_shlakesnowc, filter_shlakesnowc, &
           aerosol_vars, temperature_vars, waterstate_vars, is_lake=.true.)

      ! Check for single completely unfrozen snow layer over lake.  Modeling this ponding is unnecessary and
      ! can cause instability after the timestep when melt is completed, as the temperature after melt can be
      ! excessive because the fluxes were calculated with a fixed ground temperature of freezing, but the 
      ! phase change was unable to restore the temperature to freezing.

      do fp = 1, num_lakep
         p = filter_lakep(fp)
         c = pcolumn(p)

         j = 0

         if (snl(c) == -1 .and. h2osoi_ice(c,j) == 0._r8) then
            ! Remove layer
            ! Take extra heat of layer and release to sensible heat in order to maintain energy conservation.
            heatrem             = cpliq*h2osoi_liq(c,j)*(t_soisno(c,j) - tfrz)
            eflx_sh_tot(p)      = eflx_sh_tot(p)    + heatrem/dtime
            eflx_sh_grnd(p)     = eflx_sh_grnd(p)   + heatrem/dtime  ! Added this line 7/22/11 for consistency.
            eflx_soil_grnd(p)   = eflx_soil_grnd(p) - heatrem/dtime
            eflx_gnet(p)        = eflx_gnet(p)      - heatrem/dtime

            eflx_grnd_lake(p)   = eflx_gnet(p) - heatrem/dtime
            qflx_sl_top_soil(c) = qflx_sl_top_soil(c) + h2osno(c)
            snl(c)              = 0
            h2osno(c)           = 0._r8
            snow_depth(c)       = 0._r8
            ! Rest of snow layer book-keeping will be done below.
         else
            eflx_grnd_lake(p) = eflx_gnet(p)
         end if
      end do

      ! Check for snow layers above lake with unfrozen top layer.  Mechanically,
      ! the snow will fall into the lake and melt or turn to ice.  If the top layer has
      ! sufficient heat to melt the snow without freezing, then that will be done.
      ! Otherwise, the top layer will undergo freezing, but only if the top layer will
      ! not freeze completely.  Otherwise, let the snow layers persist and melt by diffusion.

      do fc = 1, num_lakec
         c = filter_lakec(fc)

         if (t_lake(c,1) > tfrz .and. lake_icefrac(c,1) == 0._r8 .and. snl(c) < 0) then
            unfrozen(c) = .true.
         else
            unfrozen(c) = .false.
         end if
      end do

      do j = -nlevsno+1,0
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            if (unfrozen(c)) then
               if (j == -nlevsno+1) then
                  sumsnowice(c) = 0._r8
                  heatsum(c) = 0._r8
               end if
               if (j >= snl(c)+1) then
                  sumsnowice(c) = sumsnowice(c) + h2osoi_ice(c,j)
                  heatsum(c) = heatsum(c) + h2osoi_ice(c,j)*cpice*(tfrz - t_soisno(c,j)) &
                       + h2osoi_liq(c,j)*cpliq*(tfrz - t_soisno(c,j))
               end if
            end if
         end do
      end do

      do fc = 1, num_lakec
         c = filter_lakec(fc)

         if (unfrozen(c)) then
            heatsum(c) = heatsum(c) + sumsnowice(c)*hfus
            heatrem = (t_lake(c,1) - tfrz)*cpliq*denh2o*dz_lake(c,1) - heatsum(c)

            if (heatrem + denh2o*dz_lake(c,1)*hfus > 0._r8) then            
               ! Remove snow and subtract the latent heat from the top layer.
               qflx_snomelt(c) = qflx_snomelt(c) + h2osno(c)/dtime

               eflx_snomelt(c) = eflx_snomelt(c) + h2osno(c)*hfus/dtime 

               ! update snow melt for this case
               qflx_snow_melt(c)     = qflx_snow_melt(c)  + qflx_snomelt(c)

               qflx_sl_top_soil(c) = qflx_sl_top_soil(c) + h2osno(c)

               h2osno(c) = 0._r8
               snow_depth(c) = 0._r8
               snl(c) = 0
               ! The rest of the bookkeeping for the removed snow will be done below.
               if (heatrem > 0._r8) then ! simply subtract the heat from the layer
                  t_lake(c,1) = t_lake(c,1) - heatrem/(cpliq*denh2o*dz_lake(c,1))
               else !freeze part of the layer
                  t_lake(c,1) = tfrz
                  lake_icefrac(c,1) = -heatrem/(denh2o*dz_lake(c,1)*hfus)
               end if
            end if
         end if
      end do

      ! Set empty snow layers to zero

      do j = -nlevsno+1,0
         do fc = 1, num_shlakesnowc
            c = filter_shlakesnowc(fc)
            if (j <= snl(c) .and. snl(c) > -nlevsno) then
               h2osoi_ice(c,j) = 0._r8
               h2osoi_liq(c,j) = 0._r8
               t_soisno(c,j) = 0._r8
               dz(c,j) = 0._r8
               z(c,j) = 0._r8
               zi(c,j-1) = 0._r8
            end if
         end do
      end do

      ! Build new snow filter

      call BuildSnowFilter(bounds, num_lakec, filter_lakec, &
           num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

      ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
      ! over all snow layers for history output

      do fc = 1, num_lakec
         c = filter_lakec(fc)
         snowice(c) = 0._r8
         snowliq(c) = 0._r8
      end do

      do j = -nlevsno+1, 0
         do fc = 1, num_shlakesnowc
            c = filter_shlakesnowc(fc)
            if (j >= snl(c)+1) then
               snowice(c) = snowice(c) + h2osoi_ice(c,j)
               snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
            end if
         end do
      end do

      ! Determine ending water balance and volumetric soil water

      do fc = 1, num_lakec
         c = filter_lakec(fc)
         endwb(c) = h2osno(c)
      end do

      do j = 1, nlevgrnd
         do fc = 1, num_lakec
            c = filter_lakec(fc)
            endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
            h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
         end do
      end do

      do fp = 1,num_lakep
         p = filter_lakep(fp)
         c = pcolumn(p)
         t = ptopounit(p)
         g = pgridcell(p)

         qflx_drain_perched(c) = 0._r8
         qflx_h2osfc_surf(c)   = 0._r8
         qflx_rsub_sat(c)      = 0._r8
         qflx_infl(c)          = 0._r8
         qflx_surf(c)          = 0._r8
         qflx_drain(c)         = 0._r8
         qflx_irrig(p)         = 0._r8
         qflx_irrig_col(c)     = 0._r8

         ! Insure water balance using qflx_qrgwl
         qflx_qrgwl(c)     = forc_rain(t) + forc_snow(t) - qflx_evap_tot(p) - qflx_snwcp_ice(p) - &
              (endwb(c)-begwb(c))/dtime + qflx_floodg(g)
         qflx_floodc(c)    = qflx_floodg(g)
         qflx_runoff(c)    = qflx_drain(c) + qflx_qrgwl(c)
         qflx_top_soil(c)  = qflx_prec_grnd_rain(p) + qflx_snomelt(c)

      enddo

      ! top-layer diagnostics
      do fc = 1, num_shlakesnowc
         c = filter_shlakesnowc(fc)
         h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
      end do

      ! Zero variables in columns without snow
      do fc = 1, num_shlakenosnowc
         c = filter_shlakenosnowc(fc)
            
         h2osno_top(c)      = 0._r8
         snw_rds(c,:)       = 0._r8

         ! top-layer diagnostics (spval is not averaged when computing history fields)
         snot_top(c)        = spval
         dTdz_top(c)        = spval
         snw_rds_top(c)     = spval
         sno_liq_top(c)     = spval
      end do

    end associate

  end subroutine LakeHydrology

end module LakeHydrologyMod
