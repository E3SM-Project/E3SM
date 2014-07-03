module SLakeHydrologyMod

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
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SLakeHydrology        ! Calculates soil/snow hydrology
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SLakeHydrology(bounds, num_lakec, filter_lakec, num_lakep, filter_lakep)
    ! !DESCRIPTION:
    ! Snow filter for lakes is not returned to driver.  That's okay, because it looks like it is only
    ! needed for the call to SnowAge_grain, which will be done at the bottom of this module.
    !
    ! WARNING: This subroutine assumes lake columns have one and only one pft.
    !
    ! Sequence is:
    !  SLakeHydrology:
    !    Do needed tasks from Hydrology1, Biogeophysics2, & top of Hydrology2.
    !    -> SnowWater:             change of snow mass and snow water onto soil
    !    -> SnowCompaction:        compaction of snow layers
    !    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
    !    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
    !    Add water to soil if melting has left it with open pore space.
    !    If snow layers are found above a lake with unfrozen top layer, whose top
    !    layer has enough heat to melt all the snow ice without freezing, do so
    !    and eliminate the snow layers.
    !    Cleanup and do water balance.
    !    Do SNICAR stuff and diagnostics.
    !    Call SnowAge_grain (it must be done here because the snow filters at the driver level are non-lakec only.
    !
    ! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd      , only : clm_a2l, a2l_downscaled_col
    use clm_varcon      , only : denh2o, denice, spval, hfus, tfrz, cpliq, cpice
    use SLakeCon        , only : lsadz
    use clm_varpar      , only : nlevsno, nlevgrnd, nlevsoi
    use clm_varctl      , only : iulog
    use SnowHydrologyMod, only : SnowCompaction, CombineSnowLayers, SnowWater, BuildSnowFilter
    use SnowHydrologyMod, only : DivideSnowLayers_Lake
    use clm_time_manager, only : get_step_size
    use SNICARMod       , only : SnowAge_grain, snw_rds_min
    use decompMod       , only: bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_lakec         ! number of column lake points in column filter
    integer, intent(in) :: filter_lakec(:)   ! column filter for lake points
    integer, intent(in) :: num_lakep         ! number of pft lake points in column filter
    integer, intent(in) :: filter_lakep(:)   ! pft filter for lake points
    !
    ! !LOCAL VARIABLES:
    integer  :: p,fp,g,l,c,j,fc,jtop             ! indices
    integer  :: num_shlakesnowc                  ! number of column snow points
    integer  :: filter_shlakesnowc(bounds%endc-bounds%begc+1)    ! column filter for snow points
    integer  :: num_shlakenosnowc                ! number of column non-snow points
    integer  :: filter_shlakenosnowc(bounds%endc-bounds%begc+1)  ! column filter for non-snow points
    real(r8) :: dtime                        ! land model time step (sec)
    integer  :: newnode                      ! flag when new snow node is set, (1=yes, 0=no)
    real(r8) :: dz_snowf                     ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: bifall                       ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: fracsnow(bounds%begp:bounds%endp)            ! frac of precipitation that is snow
    real(r8) :: fracrain(bounds%begp:bounds%endp)            ! frac of precipitation that is rain
    real(r8) :: qflx_prec_grnd_snow(bounds%begp:bounds%endp) ! snow precipitation incident on ground [mm/s]
    real(r8) :: qflx_prec_grnd_rain(bounds%begp:bounds%endp) ! rain precipitation incident on ground [mm/s]
    real(r8) :: qflx_evap_soi_lim            ! temporary evap_soi limited by top snow layer content [mm/s]
    real(r8) :: h2osno_temp                  ! temporary h2osno [kg/m^2]
    real(r8) :: sumsnowice(bounds%begc:bounds%endc)          ! sum of snow ice if snow layers found above unfrozen lake [kg/m&2]
    logical  :: unfrozen(bounds%begc:bounds%endc)            ! true if top lake layer is unfrozen with snow layers above
    real(r8) :: heatrem                      ! used in case above [J/m^2]
    real(r8) :: heatsum(bounds%begc:bounds%endc)             ! used in case above [J/m^2]
    real(r8) :: snowmass                     ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct              ! temporary factor used to correct for snow capping
    real(r8), parameter :: snow_bd = 250._r8 ! assumed snow bulk density (for lakes w/out resolved snow layers) [kg/m^3]
    ! Should only be used for frost below.
    !-----------------------------------------------------------------------


   associate(& 
   forc_rain            =>  a2l_downscaled_col%forc_rain , & ! Input:  [real(r8) (:)]  rain rate [mm/s]                        
   forc_snow            =>  a2l_downscaled_col%forc_snow , & ! Input:  [real(r8) (:)]  snow rate [mm/s]                        
   forc_t               =>  a2l_downscaled_col%forc_t    , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)        
   frac_sno_eff         =>  cps%frac_sno_eff             , & ! Input:  [real(r8) (:)]  needed for snicar code                  
   pcolumn              =>  pft%column                   , & ! Input:  [integer (:)]  pft's column index                       
   pgridcell            =>  pft%gridcell                 , & ! Input:  [integer (:)]  pft's gridcell index                     
   cgridcell            =>  col%gridcell                 , & ! Input:  [integer (:)]  column's gridcell                        
   clandunit            =>  col%landunit                 , & ! Input:  [integer (:)]  column's landunit                        
   snl                  =>  cps%snl                      , & ! Input:  [integer (:)]  number of snow layers                    
   t_grnd               =>  ces%t_grnd                   , & ! Input:  [real(r8) (:)]  ground temperature (Kelvin)             
   h2osno               =>  cws%h2osno                   , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   snowice              =>  cws%snowice                  , & ! Output: [real(r8) (:)]  average snow ice lens                   
   snowliq              =>  cws%snowliq                  , & ! Output: [real(r8) (:)]  average snow liquid water               
   zwt                  =>  cws%zwt                      , & ! Output: [real(r8) (:)]  water table depth                       
   fcov                 =>  cws%fcov                     , & ! Output: [real(r8) (:)]  fractional area with water table at surface
   fsat                 =>  cws%fsat                     , & ! Output: [real(r8) (:)]  fractional area with water table at surface
   qcharge              =>  cws%qcharge                  , & ! Output: [real(r8) (:)]  aquifer recharge rate (mm/s)            
   qflx_prec_grnd_col   =>  pwf_a%qflx_prec_grnd         , & ! Output: [real(r8) (:)]  water onto ground including canopy runoff [kg/(m2 s)]
   qflx_evap_grnd_col   =>  pwf_a%qflx_evap_grnd         , & ! Output: [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]
   qflx_dew_grnd_col    =>  pwf_a%qflx_dew_grnd          , & ! Output: [real(r8) (:)]  ground surface dew formation (mm H2O /s) [+]
   qflx_dew_snow_col    =>  pwf_a%qflx_dew_snow          , & ! Output: [real(r8) (:)]  surface dew added to snow pack (mm H2O /s) [+]
   qflx_sub_snow_col    =>  pwf_a%qflx_sub_snow          , & ! Output: [real(r8) (:)]  sublimation rate from snow pack (mm H2O /s) [+]
   watsat               =>  cps%watsat                   , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)
   z                    =>  cps%z                        , & ! Input:  [real(r8) (:,:)]  layer depth  (m)                      
   dz                   =>  cps%dz                       , & ! Input:  [real(r8) (:,:)]  layer thickness depth (m)             
   zi                   =>  cps%zi                       , & ! Input:  [real(r8) (:,:)]  interface depth (m)                   
   t_soisno             =>  ces%t_soisno                 , & ! Output: [real(r8) (:,:)]  snow temperature (Kelvin)             
   h2osoi_ice           =>  cws%h2osoi_ice               , & ! Output: [real(r8) (:,:)]  ice lens (kg/m2)                      
   h2osoi_liq           =>  cws%h2osoi_liq               , & ! Output: [real(r8) (:,:)]  liquid water (kg/m2)                  
   h2osoi_vol           =>  cws%h2osoi_vol               , & ! Output: [real(r8) (:,:)]  volumetric soil water [m3/m3]         
   qflx_drain           =>  cwf%qflx_drain               , & ! Output: [real(r8) (:)]  sub-surface runoff (mm H2O /s)          
   qflx_surf            =>  cwf%qflx_surf                , & ! Output: [real(r8) (:)]  surface runoff (mm H2O /s)              
   qflx_infl            =>  cwf%qflx_infl                , & ! Output: [real(r8) (:)]  infiltration (mm H2O /s)                
   qflx_qrgwl           =>  cwf%qflx_qrgwl               , & ! Output: [real(r8) (:)]  qflx_surf at glaciers, wetlands, lakes  
   qflx_runoff          =>  cwf%qflx_runoff              , & ! Output: [real(r8) (:)]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
   qflx_irrig           =>  pwf%qflx_irrig               , & ! Output: [real(r8) (:)]  irrigation flux (mm H2O /s)             
   qflx_irrig_col       =>  pwf_a%qflx_irrig             , & ! Output: [real(r8) (:)]  irrigation flux (mm H2O /s)             
   endwb                =>  cwbal%endwb                  , & ! Output: [real(r8) (:)]  water mass end of the time step         
   begwb                =>  cwbal%begwb                  , & ! Input:  [real(r8) (:)]  water mass begining of the time step    
   dz_lake              =>  cps%dz_lake                  , & ! Input:  [real(r8) (:,:)]  layer thickness for lake (m)          
   t_lake               =>  ces%t_lake                   , & ! Input:  [real(r8) (:,:)]  lake temperature (Kelvin)             
   lake_icefrac         =>  cws%lake_icefrac             , & ! Input:  [real(r8) (:,:)]  mass fraction of lake layer that is frozen
   do_capsnow           =>  cps%do_capsnow               , & ! Input:  [logical (:)]  true => do snow capping                  
   snow_depth           =>  cps%snow_depth               , & ! Input:  [real(r8) (:)]  snow height (m)                         
   qflx_snow_grnd_col   =>  pwf_a%qflx_snow_grnd         , & ! Output: [real(r8) (:)]  snow on ground after interception (mm H2O/s) [+]
   frac_iceold          =>  cps%frac_iceold              , & ! Output: [real(r8) (:,:)]  fraction of ice relative to the tot water
   qflx_evap_tot_col    =>  pwf_a%qflx_evap_tot          , & ! Output: [real(r8) (:)]  pft quantity averaged to the column (assuming one pft)
   soilalpha            =>  cws%soilalpha                , & ! Output: [real(r8) (:)]  factor that reduces ground saturated specific humidity (-)
   rootr_column         =>  cps%rootr_column             , & ! Output: [real(r8) (:,:)]  effective fraction of roots in each soil layer
   qflx_rain_grnd_col   =>  pwf_a%qflx_rain_grnd         , & ! Output: [real(r8) (:)]  rain on ground after interception (mm H2O/s) [+]
   qflx_snomelt         =>  cwf%qflx_snomelt             , & ! Input:  [real(r8) (:)]  snow melt (mm H2O /s)                   
   eflx_snomelt         =>  cef%eflx_snomelt             , & ! Input:  [real(r8) (:)]  snow melt heat flux (W/m**2)            
   ! Use column variables here.
   qflx_snwcp_ice_col   =>  pwf_a%qflx_snwcp_ice         , & ! Output: [real(r8) (:)]  excess snowfall due to snow capping (mm H2O /s) [+]
   qflx_snwcp_liq_col   =>  pwf_a%qflx_snwcp_liq         , & ! Output: [real(r8) (:)]  excess rainfall due to snow capping (mm H2O /s) [+]
   !SNICAR variables from Hydrology1
   snw_rds              =>  cps%snw_rds                  , & ! Output: [real(r8) (:,:)]  effective snow grain radius (col,lyr) [microns, m^-6]
   mss_bcpho            =>  cps%mss_bcpho                , & ! Output: [real(r8) (:,:)]  mass of hydrophobic BC in snow (col,lyr) [kg]
   mss_bcphi            =>  cps%mss_bcphi                , & ! Output: [real(r8) (:,:)]  mass of hydrophilic BC in snow (col,lyr) [kg]
   mss_bctot            =>  cps%mss_bctot                , & ! Output: [real(r8) (:,:)]  total mass of BC in snow (col,lyr) [kg]
   mss_bc_col           =>  cps%mss_bc_col               , & ! Output: [real(r8) (:)]  total column mass of BC in snow (col,lyr) [kg]
   mss_bc_top           =>  cps%mss_bc_top               , & ! Output: [real(r8) (:)]  total top-layer mass of BC (col,lyr) [kg]
   mss_ocpho            =>  cps%mss_ocpho                , & ! Output: [real(r8) (:,:)]  mass of hydrophobic OC in snow (col,lyr) [kg]
   mss_ocphi            =>  cps%mss_ocphi                , & ! Output: [real(r8) (:,:)]  mass of hydrophilic OC in snow (col,lyr) [kg]
   mss_octot            =>  cps%mss_octot                , & ! Output: [real(r8) (:,:)]  total mass of OC in snow (col,lyr) [kg]
   mss_oc_col           =>  cps%mss_oc_col               , & ! Output: [real(r8) (:)]  total column mass of OC in snow (col,lyr) [kg]
   mss_oc_top           =>  cps%mss_oc_top               , & ! Output: [real(r8) (:)]  total top-layer mass of OC (col,lyr) [kg]
   mss_dst1             =>  cps%mss_dst1                 , & ! Output: [real(r8) (:,:)]  mass of dust species 1 in snow (col,lyr) [kg]
   mss_dst2             =>  cps%mss_dst2                 , & ! Output: [real(r8) (:,:)]  mass of dust species 2 in snow (col,lyr) [kg]
   mss_dst3             =>  cps%mss_dst3                 , & ! Output: [real(r8) (:,:)]  mass of dust species 3 in snow (col,lyr) [kg]
   mss_dst4             =>  cps%mss_dst4                 , & ! Output: [real(r8) (:,:)]  mass of dust species 4 in snow (col,lyr) [kg]
   mss_dsttot           =>  cps%mss_dsttot               , & ! Output: [real(r8) (:,:)]  total mass of dust in snow (col,lyr) [kg]
   mss_dst_col          =>  cps%mss_dst_col              , & ! Output: [real(r8) (:)]  total column mass of dust in snow (col,lyr) [kg]
   mss_dst_top          =>  cps%mss_dst_top              , & ! Output: [real(r8) (:)]  total top-layer mass of dust in snow (col,lyr) [kg]
   ! Diagnostics
   snot_top             =>  cps%snot_top                 , & ! Output: [real(r8) (:)]  snow temperature in top layer (col) [K] 
   dTdz_top             =>  cps%dTdz_top                 , & ! Output: [real(r8) (:)]  temperature gradient in top layer (col) [K m-1]
   snw_rds_top          =>  cps%snw_rds_top              , & ! Output: [real(r8) (:)]  effective snow grain size, top layer(col) [microns]
   sno_liq_top          =>  cps%sno_liq_top              , & ! Output: [real(r8) (:)]  liquid water fraction in top snow layer (col) [frc]
   h2osno_top           =>  cps%h2osno_top               , & ! Output: [real(r8) (:)]  mass of snow in top layer (col) [kg]    
   ! SNICAR variables from Hydrology2
   mss_cnc_bcphi        =>  cps%mss_cnc_bcphi            , & ! Output: [real(r8) (:,:)]  mass concentration of BC species 1 (col,lyr) [kg/kg]
   mss_cnc_bcpho        =>  cps%mss_cnc_bcpho            , & ! Output: [real(r8) (:,:)]  mass concentration of BC species 2 (col,lyr) [kg/kg]
   mss_cnc_ocphi        =>  cps%mss_cnc_ocphi            , & ! Output: [real(r8) (:,:)]  mass concentration of OC species 1 (col,lyr) [kg/kg]
   mss_cnc_ocpho        =>  cps%mss_cnc_ocpho            , & ! Output: [real(r8) (:,:)]  mass concentration of OC species 2 (col,lyr) [kg/kg]
   mss_cnc_dst1         =>  cps%mss_cnc_dst1             , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 1 (col,lyr) [kg/kg]
   mss_cnc_dst2         =>  cps%mss_cnc_dst2             , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 2 (col,lyr) [kg/kg]
   mss_cnc_dst3         =>  cps%mss_cnc_dst3             , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 3 (col,lyr) [kg/kg]
   mss_cnc_dst4         =>  cps%mss_cnc_dst4             , & ! Output: [real(r8) (:,:)]  mass concentration of dust species 4 (col,lyr) [kg/kg]
   ! Flooding terms
   qflx_floodg          =>  clm_a2l%forc_flood           , & ! Input:  [real(r8) (:)]  gridcell flux of flood water from RTM   
   qflx_floodc          =>  cwf%qflx_floodc              , & ! Input:  [real(r8) (:)]  column flux of flood water from RTM     
   frost_table          =>  cws%frost_table              , & ! Input:  [real(r8) (:)]  frost table depth (m)                   
   zwt_perched          =>  cws%zwt_perched              , & ! Input:  [real(r8) (:)]  perched water table depth (m)           
   qflx_drain_perched   =>  cwf%qflx_drain_perched       , & ! Input:  [real(r8) (:)]  perched wt sub-surface runoff (mm H2O /s)
   qflx_h2osfc_surf     =>  cwf%qflx_h2osfc_surf         , & ! Input:  [real(r8) (:)]  surface water runoff (mm H2O /s)        
   qflx_snow_melt       =>  cwf%qflx_snow_melt           , & ! Input:  [real(r8) (:)]  net snow melt                           
   qflx_rsub_sat        =>  cwf%qflx_rsub_sat            , & ! Input:  [real(r8) (:)] soil saturation excess [mm h2o/s]        
   qflx_top_soil        =>  cwf%qflx_top_soil            , & ! Output: [real(r8) (:)]  net water input into soil from top (mm/s)
   qflx_sl_top_soil     =>  cwf%qflx_sl_top_soil         , & ! Output: [real(r8) (:)]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
   qflx_sub_snow        =>  pwf%qflx_sub_snow            , & ! Output: [real(r8) (:)]  sublimation rate from snow pack (mm H2O /s) [+]
   qflx_evap_grnd       =>  pwf%qflx_evap_grnd           , & ! Output: [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]
   qflx_dew_snow        =>  pwf%qflx_dew_snow            , & ! Output: [real(r8) (:)]  surface dew added to snow pack (mm H2O /s) [+]
   qflx_dew_grnd        =>  pwf%qflx_dew_grnd            , & ! Output: [real(r8) (:)]  ground surface dew formation (mm H2O /s) [+]
   qflx_prec_grnd       =>  pwf%qflx_prec_grnd           , & ! Output: [real(r8) (:)]  water onto ground including canopy runoff [kg/(m2 s)]
   qflx_snow_grnd_pft   =>  pwf%qflx_snow_grnd           , & ! Output: [real(r8) (:)]  snow on ground after interception (mm H2O/s) [+]
   qflx_rain_grnd       =>  pwf%qflx_rain_grnd           , & ! Output: [real(r8) (:)]  rain on ground after interception (mm H2O/s) [+]
   qflx_evap_tot        =>  pwf%qflx_evap_tot            , & ! Input:  [real(r8) (:)]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
   qflx_evap_soi        =>  pwf%qflx_evap_soi            , & ! Input:  [real(r8) (:)]  soil evaporation (mm H2O/s) (+ = to atm)
   qflx_snwcp_ice       =>  pwf%qflx_snwcp_ice           , & ! Output: [real(r8) (:)]  excess snowfall due to snow capping (mm H2O /s) [+]
   qflx_snwcp_liq       =>  pwf%qflx_snwcp_liq           , & ! Output: [real(r8) (:)]  excess rainfall due to snow capping (mm H2O /s) [+]
   eflx_sh_tot          =>  pef%eflx_sh_tot              , & ! Input:  [real(r8) (:)]  total sensible heat flux (W/m**2) [+ to atm]
   eflx_sh_grnd         =>  pef%eflx_sh_grnd             , & ! Input:  [real(r8) (:)]  sensible heat flux from ground (W/m**2) [+ to atm]
   eflx_soil_grnd       =>  pef%eflx_soil_grnd           , & ! Input:  [real(r8) (:)]  heat flux into snow / lake (W/m**2) [+ = into soil]
   eflx_gnet            =>  pef%eflx_gnet                , & ! Input:  [real(r8) (:)]  net heat flux into ground (W/m**2)      
   eflx_grnd_lake       =>  pef%eflx_grnd_lake             & ! Input:  [real(r8) (:)]  net heat flux into lake / snow surface, excluding light transmission (W/m**2)
   )


    ! Determine step size

    dtime = get_step_size()

    ! Add soil water to water balance.
    do j = 1, nlevgrnd
      do fc = 1, num_lakec
         c = filter_lakec(fc)
         begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
      end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Do precipitation onto ground, etc., from Hydrology1.

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)

       qflx_prec_grnd_snow(p) = forc_snow(c)
       qflx_prec_grnd_rain(p) = forc_rain(c)
       qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

       if (do_capsnow(c)) then
          qflx_snwcp_ice(p) = qflx_prec_grnd_snow(p)
          qflx_snwcp_liq(p) = qflx_prec_grnd_rain(p)
          qflx_snow_grnd_pft(p) = 0._r8
          qflx_rain_grnd(p) = 0._r8
       else
          qflx_snwcp_ice(p) = 0._r8
          qflx_snwcp_liq(p) = 0._r8
          qflx_snow_grnd_pft(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
          qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
       end if
       ! Assuming one PFT; needed for below
       qflx_snow_grnd_col(c) = qflx_snow_grnd_pft(p)
       qflx_rain_grnd_col(c) = qflx_rain_grnd(p)

    end do ! (end pft loop)

    ! Determine snow height and snow water

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
       ! U.S.Department of Agriculture Forest Service, Project F,
       ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

       if (do_capsnow(c)) then
          dz_snowf = 0._r8
       else
          if (forc_t(c) > tfrz + 2._r8) then
             bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
          else if (forc_t(c) > tfrz - 15._r8) then
             bifall=50._r8 + 1.7_r8*(forc_t(c) - tfrz + 15._r8)**1.5_r8
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
          t_soisno(c,0) = min(tfrz, forc_t(c))      ! K
          h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
          h2osoi_liq(c,0) = 0._r8                   ! kg/m2
          frac_iceold(c,0) = 1._r8

          ! intitialize SNICAR variables for fresh snow:
          snw_rds(c,0)    = snw_rds_min

          mss_bcpho(c,:)  = 0._r8
          mss_bcphi(c,:)  = 0._r8
          mss_bctot(c,:)  = 0._r8
          mss_bc_col(c)   = 0._r8
          mss_bc_top(c)   = 0._r8

          mss_ocpho(c,:)  = 0._r8
          mss_ocphi(c,:)  = 0._r8
          mss_octot(c,:)  = 0._r8
          mss_oc_col(c)   = 0._r8
          mss_oc_top(c)   = 0._r8

          mss_dst1(c,:)   = 0._r8
          mss_dst2(c,:)   = 0._r8
          mss_dst3(c,:)   = 0._r8
          mss_dst4(c,:)   = 0._r8
          mss_dsttot(c,:) = 0._r8
          mss_dst_col(c)  = 0._r8
          mss_dst_top(c)  = 0._r8
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

       else ! No snow layers: do as in HydrologyLake but with actual clmtype variables
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

    ! pft averages must be done here -- BEFORE SNOW CALCULATIONS AS THEY USE IT.
    ! for output to history tape and other uses
    ! (note that pft2col is called before SLakeHydrology, so we can't use that routine
    ! to do these column -> pft averages)
    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       qflx_evap_tot_col(c) = qflx_evap_tot(p)
       qflx_prec_grnd_col(c) = qflx_prec_grnd(p)
       qflx_evap_grnd_col(c) = qflx_evap_grnd(p)
       qflx_dew_grnd_col(c) = qflx_dew_grnd(p)
       qflx_dew_snow_col(c) = qflx_dew_snow(p)
       qflx_sub_snow_col(c) = qflx_sub_snow(p)

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

    call SnowWater(bounds, num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)


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
    
    call SnowCompaction(bounds, num_shlakesnowc, filter_shlakesnowc)
    
    ! Combine thin snow elements
    
    call CombineSnowLayers(bounds, num_shlakesnowc, filter_shlakesnowc)
    
    ! Divide thick snow elements
    
    call DivideSnowLayers_Lake(bounds, num_shlakesnowc, filter_shlakesnowc)

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
          heatrem = cpliq*h2osoi_liq(c,j)*(t_soisno(c,j) - tfrz)
          eflx_sh_tot(p) = eflx_sh_tot(p) + heatrem/dtime
          eflx_sh_grnd(p) = eflx_sh_grnd(p) + heatrem/dtime  ! Added this line 7/22/11 for consistency.
          eflx_soil_grnd(p) = eflx_soil_grnd(p) - heatrem/dtime
          eflx_gnet(p) = eflx_gnet(p) - heatrem/dtime
          eflx_grnd_lake(p) = eflx_grnd_lake(p) - heatrem/dtime
          qflx_sl_top_soil(c) = qflx_sl_top_soil(c) + h2osno(c)
          snl(c) = 0
          h2osno(c) = 0._r8
          snow_depth(c) = 0._r8
          ! Rest of snow layer book-keeping will be done below.
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
!!!!!!!!!!!!

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

!!!!!!!!!!!!!
    ! Do history variables and set special landunit runoff (adapted from end of HydrologyLake)
    do fp = 1,num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       zwt_perched(c)    = spval
       frost_table(c)    = spval
       qflx_drain_perched(c)= 0._r8
       qflx_h2osfc_surf(c)  = 0._r8
       qflx_rsub_sat(c)     = 0._r8
       qflx_infl(c)      = 0._r8
       qflx_surf(c)      = 0._r8
       qflx_drain(c)     = 0._r8
       qflx_irrig(p)     = 0._r8
       qflx_irrig_col(c) = 0._r8
       rootr_column(c,:) = spval
       soilalpha(c)      = spval
       zwt(c)            = spval
       fcov(c)           = spval
       fsat(c)           = spval
       qcharge(c)        = spval

       ! Insure water balance using qflx_qrgwl
       qflx_qrgwl(c)     = forc_rain(c) + forc_snow(c) - qflx_evap_tot(p) - qflx_snwcp_ice(p) &
                         - (endwb(c)-begwb(c))/dtime + qflx_floodg(g)
       qflx_floodc(c)    = qflx_floodg(g)
       qflx_runoff(c)    = qflx_drain(c) + qflx_surf(c) + qflx_qrgwl(c)
       qflx_top_soil(c)  = qflx_prec_grnd_rain(p) + qflx_snomelt(c)

    enddo

    !  SNICAR Code and diagnostics

    !  Calculate column-integrated aerosol masses, and
    !  mass concentrations for radiative calculations and output
    !  (based on new snow level state, after SnowFilter is rebuilt.
    !  NEEDS TO BE AFTER SnowFiler is rebuilt, otherwise there 
    !  can be zero snow layers but an active column in filter)

    do fc = 1, num_shlakesnowc
       c = filter_shlakesnowc(fc)

       ! Zero column-integrated aerosol mass before summation
       mss_bc_col(c)  = 0._r8
       mss_oc_col(c)  = 0._r8
       mss_dst_col(c) = 0._r8

       do j = -nlevsno+1, 0

          ! layer mass of snow:
          snowmass = h2osoi_ice(c,j)+h2osoi_liq(c,j)

          ! Correct the top layer aerosol mass to account for snow capping. 
          ! This approach conserves the aerosol mass concentration
          ! (but not the aerosol amss) when snow-capping is invoked

          if (j == snl(c)+1) then
             if (do_capsnow(c)) then
                snowcap_scl_fct = snowmass / (snowmass+(qflx_snwcp_ice_col(c)*dtime))
                                                        ! Make sure column variable here
                mss_bcpho(c,j) = mss_bcpho(c,j)*snowcap_scl_fct
                mss_bcphi(c,j) = mss_bcphi(c,j)*snowcap_scl_fct
                mss_ocpho(c,j) = mss_ocpho(c,j)*snowcap_scl_fct
                mss_ocphi(c,j) = mss_ocphi(c,j)*snowcap_scl_fct
                
                mss_dst1(c,j)  = mss_dst1(c,j)*snowcap_scl_fct
                mss_dst2(c,j)  = mss_dst2(c,j)*snowcap_scl_fct
                mss_dst3(c,j)  = mss_dst3(c,j)*snowcap_scl_fct
                mss_dst4(c,j)  = mss_dst4(c,j)*snowcap_scl_fct 
             endif
          endif

          if (j >= snl(c)+1) then
             mss_bctot(c,j)     = mss_bcpho(c,j) + mss_bcphi(c,j)
             mss_bc_col(c)      = mss_bc_col(c)  + mss_bctot(c,j)
             mss_cnc_bcphi(c,j) = mss_bcphi(c,j) / snowmass
             mss_cnc_bcpho(c,j) = mss_bcpho(c,j) / snowmass

             mss_octot(c,j)     = mss_ocpho(c,j) + mss_ocphi(c,j)
             mss_oc_col(c)      = mss_oc_col(c)  + mss_octot(c,j)
             mss_cnc_ocphi(c,j) = mss_ocphi(c,j) / snowmass
             mss_cnc_ocpho(c,j) = mss_ocpho(c,j) / snowmass
             
             mss_dsttot(c,j)    = mss_dst1(c,j)  + mss_dst2(c,j) + mss_dst3(c,j) + mss_dst4(c,j)
             mss_dst_col(c)     = mss_dst_col(c) + mss_dsttot(c,j)
             mss_cnc_dst1(c,j)  = mss_dst1(c,j)  / snowmass
             mss_cnc_dst2(c,j)  = mss_dst2(c,j)  / snowmass
             mss_cnc_dst3(c,j)  = mss_dst3(c,j)  / snowmass
             mss_cnc_dst4(c,j)  = mss_dst4(c,j)  / snowmass
         
          else
             !set variables of empty snow layers to zero
             snw_rds(c,j)       = 0._r8

             mss_bcpho(c,j)     = 0._r8
             mss_bcphi(c,j)     = 0._r8
             mss_bctot(c,j)     = 0._r8
             mss_cnc_bcphi(c,j) = 0._r8
             mss_cnc_bcpho(c,j) = 0._r8

             mss_ocpho(c,j)     = 0._r8
             mss_ocphi(c,j)     = 0._r8
             mss_octot(c,j)     = 0._r8
             mss_cnc_ocphi(c,j) = 0._r8
             mss_cnc_ocpho(c,j) = 0._r8

             mss_dst1(c,j)      = 0._r8
             mss_dst2(c,j)      = 0._r8
             mss_dst3(c,j)      = 0._r8
             mss_dst4(c,j)      = 0._r8
             mss_dsttot(c,j)    = 0._r8
             mss_cnc_dst1(c,j)  = 0._r8
             mss_cnc_dst2(c,j)  = 0._r8
             mss_cnc_dst3(c,j)  = 0._r8
             mss_cnc_dst4(c,j)  = 0._r8
          endif
       enddo
       
       ! top-layer diagnostics
       h2osno_top(c)  = h2osoi_ice(c,snl(c)+1) + h2osoi_liq(c,snl(c)+1)
       mss_bc_top(c)  = mss_bctot(c,snl(c)+1)
       mss_oc_top(c)  = mss_octot(c,snl(c)+1)
       mss_dst_top(c) = mss_dsttot(c,snl(c)+1)
    enddo
    
    ! Zero mass variables in columns without snow
    do fc = 1, num_shlakenosnowc
       c = filter_shlakenosnowc(fc)
            
       h2osno_top(c)      = 0._r8
       snw_rds(c,:)       = 0._r8

       mss_bc_top(c)      = 0._r8
       mss_bc_col(c)      = 0._r8    
       mss_bcpho(c,:)     = 0._r8
       mss_bcphi(c,:)     = 0._r8
       mss_bctot(c,:)     = 0._r8
       mss_cnc_bcphi(c,:) = 0._r8
       mss_cnc_bcpho(c,:) = 0._r8

       mss_oc_top(c)      = 0._r8
       mss_oc_col(c)      = 0._r8    
       mss_ocpho(c,:)     = 0._r8
       mss_ocphi(c,:)     = 0._r8
       mss_octot(c,:)     = 0._r8
       mss_cnc_ocphi(c,:) = 0._r8
       mss_cnc_ocpho(c,:) = 0._r8

       mss_dst_top(c)     = 0._r8
       mss_dst_col(c)     = 0._r8
       mss_dst1(c,:)      = 0._r8
       mss_dst2(c,:)      = 0._r8
       mss_dst3(c,:)      = 0._r8
       mss_dst4(c,:)      = 0._r8
       mss_dsttot(c,:)    = 0._r8
       mss_cnc_dst1(c,:)  = 0._r8
       mss_cnc_dst2(c,:)  = 0._r8
       mss_cnc_dst3(c,:)  = 0._r8
       mss_cnc_dst4(c,:)  = 0._r8

       ! top-layer diagnostics (spval is not averaged when computing history fields)
       snot_top(c)        = spval
       dTdz_top(c)        = spval
       snw_rds_top(c)     = spval
       sno_liq_top(c)     = spval
    enddo

    !Must be done here because the snow filter used in Hydrology2 & the Driver are for non-lake columns.
    call SnowAge_grain(bounds, num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc)

    end associate 
   end subroutine SLakeHydrology

end module SLakeHydrologyMod
