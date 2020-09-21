module BalanceCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water and energy balance check.
  !
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type
  use abortutils         , only : endrun
  use elm_varctl         , only : iulog, use_var_soil_thick
  use elm_varcon         , only : namep, namec
  use GetGlobalValuesMod , only : GetGlobalIndex
  use atm2lndType        , only : atm2lnd_type
  use glc2lndMod         , only : glc2lnd_type
  use EnergyFluxType     , only : energyflux_type
  use SolarAbsorbedType  , only : solarabs_type
  use SoilHydrologyType  , only : soilhydrology_type  
  use WaterstateType     , only : waterstate_type
  use WaterfluxType      , only : waterflux_type
  use GridcellType       , only : grc_pp
  use GridcellDataType   , only : grc_ef, grc_wf, grc_ws
  use TopounitDataType   , only : top_af ! atmospheric flux variables  
  use LandunitType       , only : lun_pp                
  use ColumnType         , only : col_pp
  use ColumnDataType     , only : col_ef, col_ws, col_wf  
  use VegetationType     , only : veg_pp
  use VegetationDataType , only : veg_ef, veg_ws
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginColWaterBalance  ! Initialize water balance check
  public :: ColWaterBalanceCheck  ! Water and energy balance check
  public :: BeginGridWaterBalance
  public :: GridBalanceCheck
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BeginColWaterBalance(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Initialize column-level water balance at beginning of time step
    !
    ! !USES:
    use subgridAveMod    , only : p2c, c2g
    use elm_varpar       , only : nlevgrnd, nlevsoi, nlevurb
    use elm_varcon       , only : spval
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall 
    use column_varcon    , only : icol_road_perv, icol_road_imperv
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds     
    integer                   , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                   , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                   , intent(in)    :: num_lakec            ! number of column non-lake points in column filter
    integer                   , intent(in)    :: filter_lakec(:)      ! column filter for non-lake points
    integer                   , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                   , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type)  , intent(inout) :: soilhydrology_vars
    type(waterstate_type)     , intent(inout) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c, p, f, j, fc                  ! indices
    real(r8):: h2osoi_vol
    !-----------------------------------------------------------------------

    associate(                                                         & 
         zi                     =>    col_pp%zi                                  , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m) 
         h2ocan_patch           =>    veg_ws%h2ocan               , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O) (pft-level)       
         h2osfc                 =>    col_ws%h2osfc                 , & ! Input:  [real(r8) (:)   ]  surface water (mm)                      
         h2osno                 =>    col_ws%h2osno                 , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
         h2osoi_ice             =>    col_ws%h2osoi_ice             , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq             =>    col_ws%h2osoi_liq             , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                  
         total_plant_stored_h2o =>    col_ws%total_plant_stored_h2o , & ! Input: [real(r8) (:) dynamic water stored in plants
         zwt                    =>    soilhydrology_vars%zwt_col                 , & ! Input:  [real(r8) (:)   ]  water table depth (m)                   
         wa                     =>    soilhydrology_vars%wa_col                  , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)    
         h2ocan_col             =>    col_ws%h2ocan                 , & ! Output: [real(r8) (:)   ]  canopy water (mm H2O) (column level)    
         begwb                  =>    col_ws%begwb                    & ! Output: [real(r8) (:)   ]  water mass begining of the time step    
         )

      ! Determine beginning water balance for time step
      ! pft-level canopy water averaged to column

      call p2c(bounds, num_nolakec, filter_nolakec, &
            h2ocan_patch(bounds%begp:bounds%endp), &
            h2ocan_col(bounds%begc:bounds%endc))
      
      if (use_var_soil_thick) then
	 do f = 1, num_hydrologyc
            c = filter_hydrologyc(f)
      	    wa(c) = 0._r8                ! Made 0 for variable soil thickness
	 end do
      end if
      
      do f = 1, num_nolakec
         c = filter_nolakec(f)
         if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_sunwall &
              .or. col_pp%itype(c) == icol_shadewall .or. col_pp%itype(c) == icol_road_imperv) then
            begwb(c) = h2ocan_col(c) + h2osno(c)
         else
            begwb(c) = h2ocan_col(c) + h2osno(c) + h2osfc(c) + wa(c)
         end if
         
      end do

      do j = 1, nlevgrnd
         do f = 1, num_nolakec
            c = filter_nolakec(f)
            if ((col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) .and. j > nlevurb) then
            else
               begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
            end if
         end do
      end do
      
      ! ---------------------------------------------------------------------------------
      ! Add stored plant water to the column water balance
      ! currently, stored plant water is only dynamic when FATES is turned on.
      ! Other orthogonal modules should not need to worry about this term,
      ! and it should be zero in all other cases and all other columns.
      ! (rgk 02-02-2017)
      ! ---------------------------------------------------------------------------------
      do f = 1, num_nolakec
         c = filter_nolakec(f)
         begwb(c) = begwb(c) + total_plant_stored_h2o(c)
      end do

      do f = 1, num_lakec
         c = filter_lakec(f)
         begwb(c) = h2osno(c)
      end do

    end associate

  end subroutine BeginColWaterBalance

   !-----------------------------------------------------------------------
   subroutine ColWaterBalanceCheck( bounds, num_do_smb_c, filter_do_smb_c, &
        atm2lnd_vars, glc2lnd_vars, solarabs_vars, waterflux_vars, waterstate_vars, &
        energyflux_vars, canopystate_vars)
     !
     ! !DESCRIPTION:
     ! This subroutine accumulates the numerical truncation errors of the water
     ! and energy balance calculation. It is helpful to see the performance of
     ! the process of integration.
     !
     ! The error for energy balance:
     !
     ! error = abs(Net radiation - change of internal energy - Sensible heat
     !             - Latent heat)
     !
     ! The error for water balance:
     !
     ! error = abs(precipitation - change of water storage - evaporation - runoff)
     !
     ! !USES:
     use elm_varcon        , only : spval
     use column_varcon     , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon     , only : icol_road_perv, icol_road_imperv
     use landunit_varcon   , only : istice_mec, istdlak, istsoil,istcrop,istwet
     use elm_varctl        , only : create_glacier_mec_landunit
     use clm_time_manager  , only : get_step_size, get_nstep
     use elm_initializeMod , only : surfalb_vars
     use domainMod         , only : ldomain
     use CanopyStateType   , only : canopystate_type
     use subgridAveMod
     use clm_time_manager  , only : get_curr_date, get_nstep
     !
     ! !ARGUMENTS:
     type(bounds_type)     , intent(in)    :: bounds  
     integer               , intent(in)    :: num_do_smb_c        ! number of columns in filter_do_smb_c
     integer               , intent(in)    :: filter_do_smb_c (:) ! column filter for points where SMB calculations are done
     type(atm2lnd_type)    , intent(in)    :: atm2lnd_vars
     type(glc2lnd_type)    , intent(in)    :: glc2lnd_vars
     type(solarabs_type)   , intent(in)    :: solarabs_vars
     type(waterflux_type)  , intent(inout) :: waterflux_vars
     type(waterstate_type) , intent(inout) :: waterstate_vars
     type(energyflux_type) , intent(inout) :: energyflux_vars
     type(canopystate_type), intent(inout) :: canopystate_vars
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c,l,t,g,fc                           ! indices
     real(r8) :: dtime                                  ! land model time step (sec)
     integer  :: nstep                                  ! time step number
     logical  :: found                                  ! flag in search loop
     integer  :: indexp,indexc,indexl,indext,indexg     ! index of first found in search loop
     real(r8) :: forc_rain_col(bounds%begc:bounds%endc) ! column level rain rate [mm/s]
     real(r8) :: forc_snow_col(bounds%begc:bounds%endc) ! column level snow rate [mm/s]
     !-----------------------------------------------------------------------

     associate(                                                                         & 
          volr                       =>    atm2lnd_vars%volr_grc                      , & ! Input:  [real(r8) (:)   ]  river water storage (m3)                 
          forc_solad                 =>    top_af%solad                               , & ! Input:  [real(r8) (:,:) ]  direct beam radiation (vis=forc_sols , nir=forc_soll) (W/m**2)
          forc_solai                 =>    top_af%solai                               , & ! Input:  [real(r8) (:,:) ]  diffuse radiation     (vis=forc_solsd, nir=forc_solld) (W/m**2)
          forc_rain                  =>    top_af%rain                                , & ! Input:  [real(r8) (:)   ]  rain rate (kg H2O/m**2/s, or mm liquid H2O/s)
          forc_snow                  =>    top_af%snow                                , & ! Input:  [real(r8) (:)   ]  snow rate (kg H2O/m**2/s, or mm liquid H2O/s)
          forc_lwrad                 =>    top_af%lwrad                               , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)
          glc_dyn_runoff_routing     =>    glc2lnd_vars%glc_dyn_runoff_routing_grc    , & ! Input:  [real(r8) (:)   ]  whether we're doing runoff routing appropriate for having a dynamic icesheet

          do_capsnow                 =>    col_ws%do_capsnow             , & ! Input:  [logical (:)    ]  true => do snow capping                  
          h2osno                     =>    col_ws%h2osno                 , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
          h2osno_old                 =>    col_ws%h2osno_old             , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O) at previous time step
          frac_sno_eff               =>    col_ws%frac_sno_eff           , & ! Input:  [real(r8) (:)   ]  effective snow fraction                 
          frac_sno                   =>    col_ws%frac_sno               , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          begwb                      =>    col_ws%begwb                  , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step    
          errh2o                     =>    col_ws%errh2o                 , & ! Output: [real(r8) (:)   ]  water conservation error (mm H2O)       
          errh2osno                  =>    col_ws%errh2osno              , & ! Output: [real(r8) (:)   ]  error in h2osno (kg m-2)                
          endwb                      =>    col_ws%endwb                  , & ! Output: [real(r8) (:)   ]  water mass end of the time step         
          total_plant_stored_h2o_col =>    col_ws%total_plant_stored_h2o , & ! Input: [real(r8) (:)   ]  water mass in plant tissues (kg m-2)
          dwb                        =>    col_wf%dwb                     , & ! Output: [real(r8) (:)   ]  change of water mass within the time step [kg/m2/s]
          qflx_rain_grnd_col         =>    col_wf%qflx_rain_grnd          , & ! Input:  [real(r8) (:)   ]  rain on ground after interception (mm H2O/s) [+]
          qflx_snow_grnd_col         =>    col_wf%qflx_snow_grnd          , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          qflx_evap_soi              =>    col_wf%qflx_evap_soi           , & ! Input:  [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
          qflx_irrig                 =>    col_wf%qflx_irrig              , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)
          qflx_surf_irrig_col        =>    col_wf%qflx_surf_irrig         , & ! Input:  [real(r8) (:)   ]  real surface irrigation flux (mm H2O /s)     
          qflx_over_supply_col       =>    col_wf%qflx_over_supply        , & ! Input:  [real(r8) (:)   ]  over supply irrigation flux (mm H2O /s)
          qflx_snwcp_ice             =>    col_wf%qflx_snwcp_ice          , & ! Input:  [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]`
          qflx_evap_tot              =>    col_wf%qflx_evap_tot           , & ! Input:  [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
          qflx_dew_snow              =>    col_wf%qflx_dew_snow           , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
          qflx_sub_snow              =>    col_wf%qflx_sub_snow           , & ! Input:  [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]
          qflx_evap_grnd             =>    col_wf%qflx_evap_grnd          , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
          qflx_dew_grnd              =>    col_wf%qflx_dew_grnd           , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
          qflx_prec_grnd             =>    col_wf%qflx_prec_grnd          , & ! Input:  [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
          qflx_snwcp_liq             =>    col_wf%qflx_snwcp_liq          , & ! Input:  [real(r8) (:)   ]  excess liquid water due to snow capping (mm H2O /s) [+]`
          qflx_snow_h2osfc           =>    col_wf%qflx_snow_h2osfc        , & ! Input:  [real(r8) (:)   ]  snow falling on surface water (mm/s)    
          qflx_h2osfc_to_ice         =>    col_wf%qflx_h2osfc_to_ice      , & ! Input:  [real(r8) (:)   ]  conversion of h2osfc to ice             
          qflx_drain_perched         =>    col_wf%qflx_drain_perched      , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_floodc                =>    col_wf%qflx_floodc             , & ! Input:  [real(r8) (:)   ]  total runoff due to flooding            
          qflx_h2osfc_surf           =>    col_wf%qflx_h2osfc_surf        , & ! Input:  [real(r8) (:)   ]  surface water runoff (mm/s)              
          qflx_snow_melt             =>    col_wf%qflx_snow_melt          , & ! Input:  [real(r8) (:)   ]  snow melt (net)                         
          qflx_surf                  =>    col_wf%qflx_surf               , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)              
          qflx_qrgwl                 =>    col_wf%qflx_qrgwl              , & ! Input:  [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes  
          qflx_drain                 =>    col_wf%qflx_drain              , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_runoff                =>    col_wf%qflx_runoff             , & ! Input:  [real(r8) (:)   ]  total runoff (mm H2O /s)                
          qflx_glcice                =>    col_wf%qflx_glcice             , & ! Input:  [real(r8) (:)   ]  flux of new glacier ice (mm H2O /s) [+ if ice grows]
          qflx_glcice_melt           =>    col_wf%qflx_glcice_melt        , & ! Input:  [real(r8) (:)   ]  ice melt (mm H2O/s)              
          qflx_glcice_frz            =>    col_wf%qflx_glcice_frz         , & ! Input:  [real(r8) (:)   ]  ice growth (mm H2O/s) [+]               
          qflx_top_soil              =>    col_wf%qflx_top_soil           , & ! Input:  [real(r8) (:)   ]  net water input into soil from top (mm/s)
          qflx_sl_top_soil           =>    col_wf%qflx_sl_top_soil        , & ! Input:  [real(r8) (:)   ]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
          qflx_liq_dynbal            =>    grc_wf%qflx_liq_dynbal         , & ! Input:  [real(r8) (:)   ]  liq runoff due to dynamic land cover change (mm H2O /s)
          qflx_ice_dynbal            =>    grc_wf%qflx_ice_dynbal         , & ! Input:  [real(r8) (:)   ]  ice runoff due to dynamic land cover change (mm H2O /s)
          snow_sources               =>    col_wf%snow_sources            , & ! Output: [real(r8) (:)   ]  snow sources (mm H2O /s)  
          snow_sinks                 =>    col_wf%snow_sinks              , & ! Output: [real(r8) (:)   ]  snow sinks (mm H2O /s)    
          qflx_lateral               =>    col_wf%qflx_lateral            , & ! Input:  [real(r8) (:)   ]  lateral flux of water to neighboring column (mm H2O /s)

          eflx_lwrad_out             =>    veg_ef%eflx_lwrad_out       , & ! Input:  [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)
          eflx_lwrad_net             =>    veg_ef%eflx_lwrad_net       , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
          eflx_sh_tot                =>    veg_ef%eflx_sh_tot          , & ! Input:  [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
          eflx_lh_tot                =>    veg_ef%eflx_lh_tot          , & ! Input:  [real(r8) (:)   ]  total latent heat flux (W/m8*2)  [+ to atm]
          eflx_soil_grnd             =>    veg_ef%eflx_soil_grnd       , & ! Input:  [real(r8) (:)   ]  soil heat flux (W/m**2) [+ = into soil]
          eflx_wasteheat_patch       =>    veg_ef%eflx_wasteheat       , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
          eflx_heat_from_ac_patch    =>    veg_ef%eflx_heat_from_ac    , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
          eflx_traffic_patch         =>    veg_ef%eflx_traffic         , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)
          eflx_dynbal                =>    grc_ef%eflx_dynbal          , & ! Input:  [real(r8) (:)   ]  energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]

          sabg_soil                  =>    solarabs_vars%sabg_soil_patch              , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
          sabg_snow                  =>    solarabs_vars%sabg_snow_patch              , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
          sabg_chk                   =>    solarabs_vars%sabg_chk_patch               , & ! Input:  [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
          fsa                        =>    solarabs_vars%fsa_patch                    , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)
          fsr                        =>    solarabs_vars%fsr_patch                    , & ! Input:  [real(r8) (:)   ]  solar radiation reflected (W/m**2)      
          sabv                       =>    solarabs_vars%sabv_patch                   , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)
          sabg                       =>    solarabs_vars%sabg_patch                   , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
          
          errsoi_col                 =>    col_ef%errsoi                 , & ! Output: [real(r8) (:)   ]  column-level soil/lake energy conservation error (W/m**2)
          errsol                     =>    veg_ef%errsol               , & ! Output: [real(r8) (:)   ]  solar radiation conservation error (W/m**2)
          errseb                     =>    veg_ef%errseb               , & ! Output: [real(r8) (:)   ]  surface energy conservation error (W/m**2)
          errlon                     =>    veg_ef%errlon               , & ! Output: [real(r8) (:)   ]  longwave radiation conservation error (W/m**2)

          fabd                       =>    surfalb_vars%fabd_patch                    , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit direct flux
          fabi                       =>    surfalb_vars%fabi_patch                    , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit indirect flux
          elai                       =>    canopystate_vars%elai_patch                , & ! Input:  [real(r8) (:,:)]  
          esai                       =>    canopystate_vars%esai_patch                , & ! Input:  [real(r8) (:,:)]  

          albd                       =>    surfalb_vars%albd_patch                    , & ! Output: [real(r8) (:,:)]  surface albedo (direct)
          albi                       =>    surfalb_vars%albi_patch                    , & ! Output: [real(r8) (:,:)]  surface albedo (diffuse)
          ftdd                       =>    surfalb_vars%ftdd_patch                    , & ! Input:  [real(r8) (:,:)]  down direct flux below canopy per unit direct flux
          ftid                       =>    surfalb_vars%ftid_patch                    , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit direct flux
          ftii                       =>    surfalb_vars%ftii_patch                    , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit diffuse flux

          netrad                     =>    veg_ef%netrad                 & ! Output: [real(r8) (:)   ]  net radiation (positive downward) (W/m**2)
          )

       ! Get step size and time step

       nstep = get_nstep()
       dtime = get_step_size()

       ! Determine column level incoming snow and rain.
       ! Assume that all columns on a topounit have the same atmospheric forcing.
       ! Assume no incident precipitation on urban wall columns (as in CanopyHydrologyMod.F90).

       do c = bounds%begc,bounds%endc
          g = col_pp%gridcell(c)
          t = col_pp%topounit(c)
          l = col_pp%landunit(c)       

          if (col_pp%itype(c) == icol_sunwall .or.  col_pp%itype(c) == icol_shadewall) then
             forc_rain_col(c) = 0.
             forc_snow_col(c) = 0.
          else
             forc_rain_col(c) = forc_rain(t)
             forc_snow_col(c) = forc_snow(t)
          end if
       end do

       ! Water balance check

       do c = bounds%begc, bounds%endc

          ! add qflx_drain_perched and qflx_flood
          if (col_pp%active(c)) then
             errh2o(c) = endwb(c) - begwb(c) &
                  - (forc_rain_col(c) + forc_snow_col(c)  + qflx_floodc(c) + qflx_surf_irrig_col(c) + qflx_over_supply_col(c) &
                  - qflx_evap_tot(c) - qflx_surf(c)  - qflx_h2osfc_surf(c) &  
                  - qflx_qrgwl(c) - qflx_drain(c) - qflx_drain_perched(c) - qflx_snwcp_ice(c) &
                  - qflx_lateral(c)) * dtime
             dwb(c) = (endwb(c)-begwb(c))/dtime

          else

             errh2o(c) = 0.0_r8
             dwb(c)    = 0.0_r8

          end if

       end do

       ! Suppose glc_dyn_runoff_routing = T:   
       ! (1) We have qflx_snwcp_ice = 0, and excess snow has been incorporated in qflx_glcice_frz.
       !     This flux must be included here to complete the water balance, because it is a
       !     sink of water as far as CLM is concerned (this water will now be owned by CISM).
       ! (2) Meltwater from ice (qflx_glcice_melt) is allowed to run off and is included in qflx_qrgwl,
       !     but the water content of the ice column has not changed (at least for now) because
       !     an equivalent ice mass has been "borrowed" from the base of the column.  So this mass
       !     has to be added back to the column, as far as the error correction is concerned, by
       !     adding back the equivalent flux*timestep.
       
       do fc = 1,num_do_smb_c
          c = filter_do_smb_c(fc)
          g = col_pp%gridcell(c)
          if (glc_dyn_runoff_routing(g)) then
             errh2o(c) = errh2o(c) + qflx_glcice_frz(c)*dtime
             errh2o(c) = errh2o(c) - qflx_glcice_melt(c)*dtime
          endif
       end do

       found = .false.
       do c = bounds%begc, bounds%endc
          if (abs(errh2o(c)) > 1.e-7_r8) then
             found = .true.
             indexc = c
          end if
       end do

       if ( found ) then
          write(iulog,*)'WARNING:  water balance error ',&
               ' nstep= ',nstep, &
               ' local indexc= ',indexc,&
               !' global indexc= ',GetGlobalIndex(decomp_index=indexc, clmlevel=namec), &
               ' errh2o= ',errh2o(indexc)

          if ((col_pp%itype(indexc) == icol_roof .or. &
               col_pp%itype(indexc) == icol_road_imperv .or. &
               col_pp%itype(indexc) == icol_road_perv) .and. &
               abs(errh2o(indexc)) > 1.e-4_r8 .and. (nstep > 2) ) then

             write(iulog,*)'clm urban model is stopping - error is greater than 1e-4 (mm)'
             write(iulog,*)'nstep                      = ',nstep
             write(iulog,*)'errh2o                     = ',errh2o(indexc)
             write(iulog,*)'forc_rain                  = ',forc_rain_col(indexc)
             write(iulog,*)'forc_snow                  = ',forc_snow_col(indexc)
             write(iulog,*)'endwb                      = ',endwb(indexc)
             write(iulog,*)'begwb                      = ',begwb(indexc)
             write(iulog,*)'qflx_evap_tot              = ',qflx_evap_tot(indexc)
             write(iulog,*)'qflx_irrig                 = ',qflx_irrig(indexc)
             write(iulog,*)'qflx_supply                = ',atm2lnd_vars%supply_grc(g)
             write(iulog,*)'f_grd                      = ',ldomain%f_grd(g)
             write(iulog,*)'qflx_surf                  = ',qflx_surf(indexc)
             write(iulog,*)'qflx_qrgwl                 = ',qflx_qrgwl(indexc)
             write(iulog,*)'qflx_drain                 = ',qflx_drain(indexc)
             write(iulog,*)'qflx_snwcp_ice             = ',qflx_snwcp_ice(indexc)
             write(iulog,*)'qflx_lateral               = ',qflx_lateral(indexc)
             write(iulog,*)'total_plant_stored_h2o_col = ',total_plant_stored_h2o_col(indexc)
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))

          else if (abs(errh2o(indexc)) > 1.e-4_r8 .and. (nstep > 2) ) then

             write(iulog,*)'clm model is stopping - error is greater than 1e-4 (mm)'
             write(iulog,*)'colum number               = ',col_pp%gridcell(indexc)
             write(iulog,*)'nstep                      = ',nstep
             write(iulog,*)'errh2o                     = ',errh2o(indexc)
             write(iulog,*)'forc_rain                  = ',forc_rain_col(indexc)
             write(iulog,*)'forc_snow                  = ',forc_snow_col(indexc)
             write(iulog,*)'endwb                      = ',endwb(indexc)
             write(iulog,*)'begwb                      = ',begwb(indexc)
             write(iulog,*)'qflx_evap_tot              = ',qflx_evap_tot(indexc)
             write(iulog,*)'qflx_irrig                 = ',qflx_irrig(indexc)
             write(iulog,*)'qflx_surf_irrig_col        = ',qflx_surf_irrig_col(indexc)
             write(iulog,*)'qflx_over_supply_col       = ',qflx_over_supply_col(indexc)
             write(iulog,*)'qflx_supply                = ',atm2lnd_vars%supply_grc(g)
             write(iulog,*)'f_grd                      = ',ldomain%f_grd(g)
             write(iulog,*)'qflx_surf                  = ',qflx_surf(indexc)
             write(iulog,*)'qflx_h2osfc_surf           = ',qflx_h2osfc_surf(indexc)
             write(iulog,*)'qflx_qrgwl                 = ',qflx_qrgwl(indexc)
             write(iulog,*)'qflx_drain                 = ',qflx_drain(indexc)
             write(iulog,*)'qflx_drain_perched         = ',qflx_drain_perched(indexc)
             write(iulog,*)'qflx_flood                 = ',qflx_floodc(indexc)
             write(iulog,*)'qflx_snwcp_ice             = ',qflx_snwcp_ice(indexc)
             write(iulog,*)'qflx_glcice_melt           = ',qflx_glcice_melt(indexc)
             write(iulog,*)'qflx_glcice_frz            = ',qflx_glcice_frz(indexc) 
             write(iulog,*)'qflx_lateral               = ',qflx_lateral(indexc)
             write(iulog,*)'total_plant_stored_h2o_col = ',total_plant_stored_h2o_col(indexc)
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

       ! Snow balance check

       do c = bounds%begc,bounds%endc
          if (col_pp%active(c)) then
             g = col_pp%gridcell(c)
             l = col_pp%landunit(c)

             ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
             ! any given time step but only if there is at least one snow layer.  h2osno 
             ! also includes snow that is part of the soil column (an initial snow layer is 
             ! only created if h2osno > 10mm).

             if (col_pp%snl(c) < 0) then
                snow_sources(c) = qflx_prec_grnd(c) + qflx_dew_snow(c) + qflx_dew_grnd(c)
                snow_sinks(c)  = qflx_sub_snow(c) + qflx_evap_grnd(c) + qflx_snow_melt(c) &
                     + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) + qflx_sl_top_soil(c)

                if (lun_pp%itype(l) == istdlak) then 
                   if ( do_capsnow(c) ) then
                      snow_sources(c) = qflx_snow_grnd_col(c) &
                           + frac_sno_eff(c) * (qflx_dew_snow(c) + qflx_dew_grnd(c) ) 

                      snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c) ) &
                           + (qflx_snwcp_ice(c) + qflx_snwcp_liq(c) - qflx_prec_grnd(c))  &
                           + qflx_snow_melt(c)  + qflx_sl_top_soil(c)
                   else
                      snow_sources(c) = qflx_snow_grnd_col(c) &
                           + frac_sno_eff(c) * (qflx_rain_grnd_col(c) &
                           +  qflx_dew_snow(c) + qflx_dew_grnd(c) ) 

                      snow_sinks(c)  = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c) ) &
                           + qflx_snow_melt(c)  + qflx_sl_top_soil(c)
                   endif
                endif

                if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop .or. lun_pp%itype(l) == istwet ) then
                   if ( do_capsnow(c) ) then
                      snow_sources(c) = frac_sno_eff(c) * (qflx_dew_snow(c) + qflx_dew_grnd(c) ) &
                           + qflx_h2osfc_to_ice(c) + qflx_prec_grnd(c)

                      snow_sinks(c) = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c)) &
                           + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                           + qflx_snow_melt(c) + qflx_sl_top_soil(c)
                   else
                      snow_sources(c) = (qflx_snow_grnd_col(c) - qflx_snow_h2osfc(c) ) &
                           + frac_sno_eff(c) * (qflx_rain_grnd_col(c) &
                           +  qflx_dew_snow(c) + qflx_dew_grnd(c) ) + qflx_h2osfc_to_ice(c)

                      snow_sinks(c) = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c)) &
                           + qflx_snow_melt(c) + qflx_sl_top_soil(c)
                   endif
                endif

                if (glc_dyn_runoff_routing(g)) then
                   ! Need to add qflx_glcice_frz to snow_sinks for the same reason as it is
                   ! added to errh2o above - see the comment above for details.
                   snow_sinks(c) = snow_sinks(c) + qflx_glcice_frz(c)
                end if

                errh2osno(c) = (h2osno(c) - h2osno_old(c)) - (snow_sources(c) - snow_sinks(c)) * dtime
             else
                snow_sources(c) = 0._r8
                snow_sinks(c) = 0._r8
                errh2osno(c) = 0._r8
             end if

          end if
       end do

       found = .false.
       do c = bounds%begc,bounds%endc
          if (col_pp%active(c)) then
             if (abs(errh2osno(c)) > 1.0e-7_r8) then
                found = .true.
                indexc = c
             end if
          end if
       end do
       if ( found ) then
          write(iulog,*)'WARNING:  snow balance error '
          write(iulog,*)'nstep= ',nstep, &
               ' local indexc= ',indexc, &
               !' global indexc= ',GetGlobalIndex(decomp_index=indexc, clmlevel=namec), &
               ' col_pp%itype= ',col_pp%itype(indexc), &
               ' lun_pp%itype= ',lun_pp%itype(col_pp%landunit(indexc)), &
               ' errh2osno= ',errh2osno(indexc)

          if (abs(errh2osno(indexc)) > 1.e-4_r8 .and. (nstep > 2) ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-4 (mm)'
             write(iulog,*)'nstep            = ',nstep
             write(iulog,*)'errh2osno        = ',errh2osno(indexc)
             write(iulog,*)'snl              = ',col_pp%snl(indexc)
             write(iulog,*)'h2osno           = ',h2osno(indexc)
             write(iulog,*)'h2osno_old       = ',h2osno_old(indexc)
             write(iulog,*)'snow_sources     = ',snow_sources(indexc)
             write(iulog,*)'snow_sinks       = ',snow_sinks(indexc)
             write(iulog,*)'qflx_prec_grnd   = ',qflx_prec_grnd(indexc)*dtime
             write(iulog,*)'qflx_sub_snow    = ',qflx_sub_snow(indexc)*dtime
             write(iulog,*)'qflx_evap_grnd   = ',qflx_evap_grnd(indexc)*dtime
             write(iulog,*)'qflx_top_soil    = ',qflx_top_soil(indexc)*dtime
             write(iulog,*)'qflx_dew_snow    = ',qflx_dew_snow(indexc)*dtime
             write(iulog,*)'qflx_dew_grnd    = ',qflx_dew_grnd(indexc)*dtime
             write(iulog,*)'qflx_snwcp_ice   = ',qflx_snwcp_ice(indexc)*dtime
             write(iulog,*)'qflx_snwcp_liq   = ',qflx_snwcp_liq(indexc)*dtime
             write(iulog,*)'qflx_sl_top_soil = ',qflx_sl_top_soil(indexc)*dtime
             if (create_glacier_mec_landunit) then
                write(iulog,*)'qflx_glcice_frz  = ',qflx_glcice_frz(indexc)*dtime
             end if
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

       ! Energy balance checks

       do p = bounds%begp, bounds%endp
          if (veg_pp%active(p)) then
             c = veg_pp%column(p)
             l = veg_pp%landunit(p)
             t = veg_pp%topounit(p)
             g = veg_pp%gridcell(p)

             ! Solar radiation energy balance
             ! Do not do this check for an urban pft since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun_pp%urbpoi(l)) then
                errsol(p) = fsa(p) + fsr(p) &
                     - (forc_solad(t,1) + forc_solad(t,2) + forc_solai(t,1) + forc_solai(t,2))
             else
                errsol(p) = spval
             end if

             ! Longwave radiation energy balance
             ! Do not do this check for an urban pft since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun_pp%urbpoi(l)) then
                errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(t)
             else
                errlon(p) = spval
             end if

             ! Surface energy balance
             ! Changed to using (eflx_lwrad_net) here instead of (forc_lwrad - eflx_lwrad_out) because
             ! there are longwave interactions between urban columns (and therefore patches). 
             ! For surfaces other than urban, (eflx_lwrad_net) equals (forc_lwrad - eflx_lwrad_out),
             ! and a separate check is done above for these terms.

             if (.not. lun_pp%urbpoi(l)) then
                errseb(p) = sabv(p) + sabg_chk(p) + forc_lwrad(t) - eflx_lwrad_out(p) &
                     - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p)
             else
                errseb(p) = sabv(p) + sabg(p) &
                     - eflx_lwrad_net(p) &
                     - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) &
                     + eflx_wasteheat_patch(p) + eflx_heat_from_ac_patch(p) + eflx_traffic_patch(p)
             end if
             !TODO MV - move this calculation to a better place - does not belong in BalanceCheck 
             netrad(p) = fsa(p) - eflx_lwrad_net(p) 
          end if
       end do

       ! Solar radiation energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (veg_pp%active(p)) then
             if ( (errsol(p) /= spval) .and. (abs(errsol(p)) > 1.e-7_r8) ) then
                found = .true.
                indexp = p
                indext = veg_pp%topounit(indexp)
                indexg = veg_pp%gridcell(indexp)
             end if
          end if
       end do
       if ( found  .and. (nstep > 2) ) then
          write(iulog,*)'WARNING:: BalanceCheck, solar radiation balance error (W/m2)'
          write(iulog,*)'nstep         = ',nstep
          write(iulog,*)'errsol        = ',errsol(indexp)
          if (abs(errsol(indexp)) > 1.e-5_r8 ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (W/m2)'
             write(iulog,*)'fsa           = ',fsa(indexp)
             write(iulog,*)'fsr           = ',fsr(indexp)
             write(iulog,*)'forc_solad(1) = ',forc_solad(indext,1)
             write(iulog,*)'forc_solad(2) = ',forc_solad(indext,2)
             write(iulog,*)'forc_solai(1) = ',forc_solai(indext,1)
             write(iulog,*)'forc_solai(2) = ',forc_solai(indext,2)
             write(iulog,*)'forc_tot      = ',forc_solad(indext,1)+forc_solad(indext,2) &
               +forc_solai(indext,1)+forc_solai(indext,2)
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

       ! Longwave radiation energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (veg_pp%active(p)) then
             if ( (errlon(p) /= spval) .and. (abs(errlon(p)) > 1.e-7_r8) ) then
                found = .true.
                indexp = p
             end if
          end if
       end do
       if ( found  .and. (nstep > 2) ) then
          write(iulog,*)'WARNING: BalanceCheck: longwave energy balance error (W/m2)' 
          write(iulog,*)'nstep        = ',nstep 
          write(iulog,*)'errlon       = ',errlon(indexp)
          if (abs(errlon(indexp)) > 1.e-5_r8 ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (W/m2)'
             call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

       ! Surface energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (veg_pp%active(p)) then
             if (abs(errseb(p)) > 1.e-7_r8 ) then
                found = .true.
                indexp = p
                indexc = veg_pp%column(indexp)
                indext = veg_pp%topounit(indexp)
             end if
          end if
       end do
       if ( found  .and. (nstep > 2) ) then
          write(iulog,*)'WARNING: BalanceCheck: surface flux energy balance error (W/m2)'
          write(iulog,*)'nstep          = ' ,nstep
          write(iulog,*)'errseb         = ' ,errseb(indexp)
          if (abs(errseb(indexp)) > 1.e-5_r8 ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-5 (W/m2)'
             write(iulog,*)'sabv           = ' ,sabv(indexp)

             write(iulog,*)'sabg           = ' ,sabg(indexp), ((1._r8- frac_sno(indexc))*sabg_soil(indexp) + &
                  frac_sno(indexc)*sabg_snow(indexp)),sabg_chk(indexp)

             write(iulog,*)'forc_tot      = '  ,forc_solad(indext,1) + forc_solad(indext,2) + &
                  forc_solai(indext,1) + forc_solai(indext,2)

             write(iulog,*)'eflx_lwrad_net = ' ,eflx_lwrad_net(indexp)
             write(iulog,*)'eflx_sh_tot    = ' ,eflx_sh_tot(indexp)
             write(iulog,*)'eflx_lh_tot    = ' ,eflx_lh_tot(indexp)
             write(iulog,*)'eflx_soil_grnd = ' ,eflx_soil_grnd(indexp)
             write(iulog,*)'fsa fsr = '        ,fsa(indexp),    fsr(indexp)
             write(iulog,*)'fabd fabi = '      ,fabd(indexp,:), fabi(indexp,:)
             write(iulog,*)'albd albi = '      ,albd(indexp,:), albi(indexp,:)
             write(iulog,*)'ftii ftdd ftid = ' ,ftii(indexp,:), ftdd(indexp,:),ftid(indexp,:)
             write(iulog,*)'elai esai = '      ,elai(indexp),   esai(indexp)      
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

       ! Soil energy balance check

       found = .false.
       do c = bounds%begc,bounds%endc
          if (col_pp%active(c)) then
             if (abs(errsoi_col(c)) > 1.0e-5_r8 ) then
                found = .true.
                indexc = c
             end if
          end if
       end do
       if ( found ) then
          write(iulog,*)'WARNING: BalanceCheck: soil balance error (W/m2)'
          write(iulog,*)'nstep         = ',nstep
          write(iulog,*)'errsoi_col    = ',errsoi_col(indexc)
          write(iulog,*)'colum number  = ',col_pp%gridcell(indexc)
          if (abs(errsoi_col(indexc)) > 1.e-4_r8 .and. (nstep > 2) ) then
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

     end associate

   end subroutine ColWaterBalanceCheck

  !-----------------------------------------------------------------------
  subroutine BeginGridWaterBalance(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Initialize column-level water balance at beginning of time step
    !
    ! !USES:
    use subgridAveMod , only : p2c,c2g
    use elm_varpar    , only : nlevgrnd, nlevsoi, nlevurb
    use column_varcon , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon , only : icol_road_perv, icol_road_imperv
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                   , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                   , intent(in)    :: num_lakec            ! number of column non-lake points in column filter
    integer                   , intent(in)    :: filter_lakec(:)      ! column filter for non-lake points
    integer                   , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                   , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type)  , intent(inout) :: soilhydrology_vars
    type(waterstate_type)     , intent(inout) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c, p, f, j, fc,g                  ! indices
    real(r8) :: h2osoi_vol
    real(r8) :: h2ocan_col(bounds%begc:bounds%endc)
    real(r8) :: begwb_col (bounds%begc:bounds%endc)
    real(r8) :: h2osoi_liq_depth_intg(bounds%begc:bounds%endc)
    real(r8) :: h2osoi_ice_depth_intg(bounds%begc:bounds%endc)
    real(r8) :: wa_local_col(bounds%begc:bounds%endc)

    associate(                                                                        &
         zi                        =>    col_pp%zi                                  , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         h2ocan_patch              =>    veg_ws%h2ocan               , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O) (pft-level)
         h2osfc                    =>    col_ws%h2osfc                 , & ! Input:  [real(r8) (:)   ]  surface water (mm)
         h2osno                    =>    col_ws%h2osno                 , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
         h2osoi_ice                =>    col_ws%h2osoi_ice             , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq                =>    col_ws%h2osoi_liq             , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         total_plant_stored_h2o    =>    col_ws%total_plant_stored_h2o , & ! Input:  [real(r8) (:)   ]  dynamic water stored in plants
         zwt                       =>    soilhydrology_vars%zwt_col                 , & ! Input:  [real(r8) (:)   ]  water table depth (m)
         wa                        =>    soilhydrology_vars%wa_col                  , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)
         beg_wa_grc                =>    soilhydrology_vars%beg_wa_grc              , & ! Output: [real(r8) (:)   ]  grid-level water in the unconfined aquifer at begining of the time step (mm)
         begwb_grc                 =>    grc_ws%begwb                  , & ! Output: [real(r8) (:)   ]  grid-level water mass at begining of the time step (mm)
         beg_h2ocan_grc            =>    grc_ws%beg_h2ocan             , & ! Output: [real(r8) (:)   ]  grid-level canopy water at begining of the time step (mm)
         beg_h2osno_grc            =>    grc_ws%beg_h2osno             , & ! Output: [real(r8) (:)   ]  grid-level snow at begining of the time step (mm)
         beg_h2osfc_grc            =>    grc_ws%beg_h2osfc             , & ! Output: [real(r8) (:)   ]  grid-level surface water at begining of the time step (mm)
         beg_h2osoi_liq_grc        =>    grc_ws%beg_h2osoi_liq         , & ! Output: [real(r8) (:)   ]  grid-level depth integrated liquid soil water at begining of the time step (mm)
         beg_h2osoi_ice_grc        =>    grc_ws%beg_h2osoi_ice         , & ! Output: [real(r8) (:)   ]  grid-level depth integrated ice soil water at begining of the time step (mm)
         h2osoi_liq_depth_intg_col =>    col_ws%h2osoi_liq_depth_intg  , &
         h2osoi_ice_depth_intg_col =>    col_ws%h2osoi_ice_depth_intg    &
         )

      ! Set to zero
      begwb_col (bounds%begc:bounds%endc) = 0._r8
      h2ocan_col(bounds%begc:bounds%endc) = 0._r8
      h2osoi_liq_depth_intg(bounds%begc:bounds%endc) = 0._r8
      h2osoi_ice_depth_intg(bounds%begc:bounds%endc) = 0._r8
      h2osoi_liq_depth_intg_col(bounds%begc:bounds%endc) = 0._r8
      h2osoi_ice_depth_intg_col(bounds%begc:bounds%endc) = 0._r8

      ! Determine beginning water balance for time step
      ! pft-level canopy water averaged to column

      call p2c(bounds, num_nolakec, filter_nolakec, &
            h2ocan_patch(bounds%begp:bounds%endp), &
            h2ocan_col(bounds%begc:bounds%endc))

      wa_local_col(bounds%begc:bounds%endc) = wa(bounds%begc:bounds%endc)

      do f = 1, num_nolakec
         c = filter_nolakec(f)
         g = col_pp%gridcell(c)
         if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_sunwall &
              .or. col_pp%itype(c) == icol_shadewall .or. col_pp%itype(c) == icol_road_imperv) then
            begwb_col(c) = h2ocan_col(c) + h2osno(c)
            wa_local_col(c) = 0._r8
         else
            begwb_col(c) = h2ocan_col(c) + h2osno(c) + h2osfc(c) + wa(c)
         end if
         begwb_col(c) = begwb_col(c) + total_plant_stored_h2o(c)
      end do

      do j = 1, nlevgrnd
         do f = 1, num_nolakec
            c = filter_nolakec(f)
            if ((col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall &
                 .or. col_pp%itype(c) == icol_roof) .and. j > nlevurb) then
            else
               begwb_col(c) = begwb_col(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
               h2osoi_liq_depth_intg(c) = h2osoi_liq_depth_intg(c) + h2osoi_liq(c,j)
               h2osoi_ice_depth_intg(c) = h2osoi_ice_depth_intg(c) + h2osoi_ice(c,j)
            end if
         end do
      end do

      do f = 1, num_lakec
         c = filter_lakec(f)
         begwb_col(c) = h2osno(c)
         do j = 1, nlevgrnd
            begwb_col(c) = begwb_col(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
            h2osoi_liq_depth_intg(c) = h2osoi_liq_depth_intg(c) + h2osoi_liq(c,j)
            h2osoi_ice_depth_intg(c) = h2osoi_ice_depth_intg(c) + h2osoi_ice(c,j)
         enddo
      end do

      call c2g(bounds, begwb_col(bounds%begc:bounds%endc), &
           begwb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, wa_local_col(bounds%begc:bounds%endc), &
           beg_wa_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2ocan_col(bounds%begc:bounds%endc), &
           beg_h2ocan_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osno(bounds%begc:bounds%endc), &
           beg_h2osno_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osfc(bounds%begc:bounds%endc), &
           beg_h2osfc_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osoi_liq_depth_intg(bounds%begc:bounds%endc), &
           beg_h2osoi_liq_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osoi_ice_depth_intg(bounds%begc:bounds%endc), &
           beg_h2osoi_ice_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
      
    end associate

  end subroutine BeginGridWaterBalance

   !-----------------------------------------------------------------------
   subroutine GridBalanceCheck( bounds, num_do_smb_c, filter_do_smb_c, &
        atm2lnd_vars, glc2lnd_vars, solarabs_vars, waterflux_vars, &
        waterstate_vars, energyflux_vars, canopystate_vars, soilhydrology_vars)
     !
     ! !DESCRIPTION:
     !
     ! !USES:
     use elm_varcon        , only : spval
     use column_varcon     , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon     , only : icol_road_perv, icol_road_imperv
     use landunit_varcon   , only : istice_mec, istdlak, istsoil,istcrop,istwet
     use elm_varctl        , only : create_glacier_mec_landunit
     use clm_time_manager  , only : get_step_size, get_nstep
     use elm_initializeMod , only : surfalb_vars
     use CanopyStateType   , only : canopystate_type
     use subgridAveMod
     !
     ! !ARGUMENTS:
     type(bounds_type)     , intent(in)    :: bounds  
     integer               , intent(in)    :: num_do_smb_c        ! number of columns in filter_do_smb_c
     integer               , intent(in)    :: filter_do_smb_c (:) ! column filter for points where SMB calculations are done
     type(atm2lnd_type)    , intent(in)    :: atm2lnd_vars
     type(glc2lnd_type)    , intent(in)    :: glc2lnd_vars
     type(solarabs_type)   , intent(in)    :: solarabs_vars
     type(waterflux_type)  , intent(inout) :: waterflux_vars
     type(waterstate_type) , intent(inout) :: waterstate_vars
     type(energyflux_type) , intent(inout) :: energyflux_vars
     type(canopystate_type), intent(inout) :: canopystate_vars
     type(soilhydrology_type), intent(inout) :: soilhydrology_vars
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c,l,g,fc                             ! indices
     real(r8) :: dtime                                  ! land model time step (sec)
     real(r8) :: qflx_net_col (bounds%begc:bounds%endc)
     real(r8) :: forc_rain_col(bounds%begc:bounds%endc) ! column level rain rate [mm/s]
     real(r8) :: forc_snow_col(bounds%begc:bounds%endc) ! column level snow rate [mm/s]
     real(r8) :: wa_local_col(bounds%begc:bounds%endc)

     associate(                                                                         & 
          volr                       =>    atm2lnd_vars%volr_grc                      , & ! Input:  [real(r8) (:)   ]  river water storage (m3)                 
          forc_solad                 =>    atm2lnd_vars%forc_solad_grc                , & ! Input:  [real(r8) (:,:) ]  direct beam radiation (vis=forc_sols , nir=forc_soll )
          forc_solai                 =>    atm2lnd_vars%forc_solai_grc                , & ! Input:  [real(r8) (:,:) ]  diffuse radiation     (vis=forc_solsd, nir=forc_solld)
          forc_rain                  =>    atm2lnd_vars%forc_rain_downscaled_col      , & ! Input:  [real(r8) (:)   ]  rain rate [mm/s]
          forc_snow                  =>    atm2lnd_vars%forc_snow_downscaled_col      , & ! Input:  [real(r8) (:)   ]  snow rate [mm/s]
          forc_lwrad                 =>    atm2lnd_vars%forc_lwrad_downscaled_col     , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)
          glc_dyn_runoff_routing     =>    glc2lnd_vars%glc_dyn_runoff_routing_grc    , & ! Input:  [real(r8) (:)   ]  whether we're doing runoff routing appropriate for having a dynamic icesheet

          do_capsnow                 =>    col_ws%do_capsnow             , & ! Input:  [logical (:)    ]  true => do snow capping                  
          h2osno_col                 =>    col_ws%h2osno                 , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
          h2osno_old                 =>    col_ws%h2osno_old             , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O) at previous time step
          frac_sno_eff               =>    col_ws%frac_sno_eff           , & ! Input:  [real(r8) (:)   ]  effective snow fraction                 
          frac_sno                   =>    col_ws%frac_sno               , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          begwb_col                  =>    col_ws%begwb                  , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step    
          begwb_grc                  =>    grc_ws%begwb                  , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step    
          errh2o                     =>    col_ws%errh2o                 , & ! Output: [real(r8) (:)   ]  water conservation error (mm H2O)       
          errh2o_grc                 =>    grc_ws%errh2o                 , & ! Output: [real(r8) (:)   ]  water conservation error (mm H2O)
          errh2osno                  =>    col_ws%errh2osno              , & ! Output: [real(r8) (:)   ]  error in h2osno (kg m-2)                
          endwb_col                  =>    col_ws%endwb                  , & ! Output: [real(r8) (:)   ]  water mass end of the time step         
          endwb_grc                  =>    grc_ws%endwb                  , & ! Output: [real(r8) (:)   ]  water mass end of the time step         
          total_plant_stored_h2o_col =>    col_ws%total_plant_stored_h2o , & ! Input: [real(r8) (:)   ]  water mass in plant tissues (kg m-2)
          dwb                        =>    col_wf%dwb                     , & ! Output: [real(r8) (:)   ]  change of water mass within the time step [kg/m2/s]
          qflx_rain_grnd_col         =>    col_wf%qflx_rain_grnd          , & ! Input:  [real(r8) (:)   ]  rain on ground after interception (mm H2O/s) [+]
          qflx_snow_grnd_col         =>    col_wf%qflx_snow_grnd          , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          qflx_evap_soi              =>    col_wf%qflx_evap_soi           , & ! Input:  [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
          qflx_irrig                 =>    col_wf%qflx_irrig              , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)             
          qflx_snwcp_ice             =>    col_wf%qflx_snwcp_ice          , & ! Input:  [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]`
          qflx_evap_tot              =>    col_wf%qflx_evap_tot           , & ! Input:  [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
          qflx_dew_snow              =>    col_wf%qflx_dew_snow           , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
          qflx_sub_snow              =>    col_wf%qflx_sub_snow           , & ! Input:  [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]
          qflx_evap_grnd             =>    col_wf%qflx_evap_grnd          , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
          qflx_dew_grnd              =>    col_wf%qflx_dew_grnd           , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
          qflx_prec_grnd             =>    col_wf%qflx_prec_grnd          , & ! Input:  [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
          qflx_snwcp_liq             =>    col_wf%qflx_snwcp_liq          , & ! Input:  [real(r8) (:)   ]  excess liquid water due to snow capping (mm H2O /s) [+]`
          qflx_snow_h2osfc           =>    col_wf%qflx_snow_h2osfc        , & ! Input:  [real(r8) (:)   ]  snow falling on surface water (mm/s)    
          qflx_h2osfc_to_ice         =>    col_wf%qflx_h2osfc_to_ice      , & ! Input:  [real(r8) (:)   ]  conversion of h2osfc to ice             
          qflx_drain_perched         =>    col_wf%qflx_drain_perched      , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_floodc                =>    col_wf%qflx_floodc             , & ! Input:  [real(r8) (:)   ]  total runoff due to flooding            
          qflx_h2osfc_surf           =>    col_wf%qflx_h2osfc_surf        , & ! Input:  [real(r8) (:)   ]  surface water runoff (mm/s)              
          qflx_snow_melt             =>    col_wf%qflx_snow_melt          , & ! Input:  [real(r8) (:)   ]  snow melt (net)                         
          qflx_surf                  =>    col_wf%qflx_surf               , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)              
          qflx_qrgwl                 =>    col_wf%qflx_qrgwl              , & ! Input:  [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes  
          qflx_drain                 =>    col_wf%qflx_drain              , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_runoff                =>    col_wf%qflx_runoff             , & ! Input:  [real(r8) (:)   ]  total runoff (mm H2O /s)                
          qflx_glcice                =>    col_wf%qflx_glcice             , & ! Input:  [real(r8) (:)   ]  flux of new glacier ice (mm H2O /s) [+ if ice grows]
          qflx_glcice_melt           =>    col_wf%qflx_glcice_melt        , & ! Input:  [real(r8) (:)   ]  ice melt (mm H2O/s)              
          qflx_glcice_frz            =>    col_wf%qflx_glcice_frz         , & ! Input:  [real(r8) (:)   ]  ice growth (mm H2O/s) [+]               
          qflx_top_soil              =>    col_wf%qflx_top_soil           , & ! Input:  [real(r8) (:)   ]  net water input into soil from top (mm/s)
          qflx_sl_top_soil           =>    col_wf%qflx_sl_top_soil        , & ! Input:  [real(r8) (:)   ]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
          qflx_liq_dynbal            =>    grc_wf%qflx_liq_dynbal         , & ! Input:  [real(r8) (:)   ]  liq runoff due to dynamic land cover change (mm H2O /s)
          qflx_ice_dynbal            =>    grc_wf%qflx_ice_dynbal         , & ! Input:  [real(r8) (:)   ]  ice runoff due to dynamic land cover change (mm H2O /s)
          snow_sources               =>    col_wf%snow_sources            , & ! Output: [real(r8) (:)   ]  snow sources (mm H2O /s)  
          snow_sinks                 =>    col_wf%snow_sinks              , & ! Output: [real(r8) (:)   ]  snow sinks (mm H2O /s)    
          qflx_lateral               =>    col_wf%qflx_lateral            , & ! Input:  [real(r8) (:)   ]  lateral flux of water to neighboring column (mm H2O /s)

          eflx_lwrad_out             =>    veg_ef%eflx_lwrad_out       , & ! Input:  [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)
          eflx_lwrad_net             =>    veg_ef%eflx_lwrad_net       , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
          eflx_sh_tot                =>    veg_ef%eflx_sh_tot          , & ! Input:  [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
          eflx_lh_tot                =>    veg_ef%eflx_lh_tot          , & ! Input:  [real(r8) (:)   ]  total latent heat flux (W/m8*2)  [+ to atm]
          eflx_soil_grnd             =>    veg_ef%eflx_soil_grnd       , & ! Input:  [real(r8) (:)   ]  soil heat flux (W/m**2) [+ = into soil] 
          eflx_wasteheat_patch       =>    veg_ef%eflx_wasteheat       , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
          eflx_heat_from_ac_patch    =>    veg_ef%eflx_heat_from_ac    , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
          eflx_traffic_patch         =>    veg_ef%eflx_traffic         , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
          eflx_dynbal                =>    grc_ef%eflx_dynbal          , & ! Input:  [real(r8) (:)   ]  energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]

          sabg_soil                  =>    solarabs_vars%sabg_soil_patch              , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
          sabg_snow                  =>    solarabs_vars%sabg_snow_patch              , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
          sabg_chk                   =>    solarabs_vars%sabg_chk_patch               , & ! Input:  [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
          fsa                        =>    solarabs_vars%fsa_patch                    , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)
          fsr                        =>    solarabs_vars%fsr_patch                    , & ! Input:  [real(r8) (:)   ]  solar radiation reflected (W/m**2)      
          sabv                       =>    solarabs_vars%sabv_patch                   , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)
          sabg                       =>    solarabs_vars%sabg_patch                   , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
          
          errsoi_col                 =>    col_ef%errsoi                 , & ! Output: [real(r8) (:)   ]  column-level soil/lake energy conservation error (W/m**2)
          errsol                     =>    veg_ef%errsol               , & ! Output: [real(r8) (:)   ]  solar radiation conservation error (W/m**2)
          errseb                     =>    veg_ef%errseb               , & ! Output: [real(r8) (:)   ]  surface energy conservation error (W/m**2)
          errlon                     =>    veg_ef%errlon               , & ! Output: [real(r8) (:)   ]  longwave radiation conservation error (W/m**2)

          fabd                       =>    surfalb_vars%fabd_patch                    , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit direct flux
          fabi                       =>    surfalb_vars%fabi_patch                    , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit indirect flux
          elai                       =>    canopystate_vars%elai_patch                , & ! Input:  [real(r8) (:,:)]  
          esai                       =>    canopystate_vars%esai_patch                , & ! Input:  [real(r8) (:,:)]  

          albd                       =>    surfalb_vars%albd_patch                    , & ! Output: [real(r8) (:,:)]  surface albedo (direct)
          albi                       =>    surfalb_vars%albi_patch                    , & ! Output: [real(r8) (:,:)]  surface albedo (diffuse)
          ftdd                       =>    surfalb_vars%ftdd_patch                    , & ! Input:  [real(r8) (:,:)]  down direct flux below canopy per unit direct flux
          ftid                       =>    surfalb_vars%ftid_patch                    , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit direct flux
          ftii                       =>    surfalb_vars%ftii_patch                    , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit diffuse flux

          netrad                     =>    veg_ef%netrad               , & ! Output: [real(r8) (:)   ]  net radiation (positive downward) (W/m**2)
          h2ocan_patch               =>    veg_ws%h2ocan               , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O) (pft-level)
          wa                         =>    soilhydrology_vars%wa_col                  , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)
          h2ocan_col                 =>    col_ws%h2ocan                 , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O)
          h2osfc_col                 =>    col_ws%h2osfc                 , & ! Input:  [real(r8) (:)   ]  surface water (mm)
          h2osoi_liq_depth_intg      =>    col_ws%h2osoi_liq_depth_intg  , & ! Input:  [real(r8) (:)   ]  depth integrated liquid soil water (kg/m**2)
          h2osoi_ice_depth_intg      =>    col_ws%h2osoi_ice_depth_intg  , & ! Input:  [real(r8) (:)   ]  depth integrated ice soil water (kg/m**2)
          end_wa_grc                 =>    soilhydrology_vars%end_wa_grc              , & ! Output: [real(r8) (:)   ]  grid-level water in the unconfined aquifer at end of the time step (mm)
          end_h2ocan_grc             =>    grc_ws%end_h2ocan             , & ! Output: [real(r8) (:)   ]  grid-level canopy water at end of the time step (mm)
          end_h2osno_grc             =>    grc_ws%end_h2osno             , & ! Output: [real(r8) (:)   ]  grid-level snow at end of the time step (mm)
          end_h2osfc_grc             =>    grc_ws%end_h2osfc             , & ! Output: [real(r8) (:)   ]  grid-level surface water at end of the time step (mm)
          end_h2osoi_liq_grc         =>    grc_ws%end_h2osoi_liq         , & ! Output: [real(r8) (:)   ]  grid-level depth integrated liquid soil water at end of the time step (mm)
          end_h2osoi_ice_grc         =>    grc_ws%end_h2osoi_ice           & ! Output: [real(r8) (:)   ]  grid-level depth integrated liquid soil water at end of the time step (mm)
          )

       dtime = get_step_size()

       wa_local_col(bounds%begc:bounds%endc) = wa(bounds%begc:bounds%endc)

       do c = bounds%begc,bounds%endc
          g = col_pp%gridcell(c)
          l = col_pp%landunit(c)       

          if (col_pp%itype(c) == icol_sunwall .or.  col_pp%itype(c) == icol_shadewall) then
             forc_rain_col(c) = 0.
             forc_snow_col(c) = 0.
          else
             forc_rain_col(c) = forc_rain(c)
             forc_snow_col(c) = forc_snow(c)
          end if
          if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_sunwall &
               .or. col_pp%itype(c) == icol_shadewall .or. col_pp%itype(c) == icol_road_imperv) then
             wa_local_col(c) = 0._r8
          end if
       end do

      do c = bounds%begc, bounds%endc

         ! add qflx_drain_perched and qflx_flood
         if (col_pp%active(c)) then

            qflx_net_col(c) = &
                 - forc_rain_col(c) - forc_snow_col(c)  - qflx_floodc(c) - qflx_irrig(c) &
                 + qflx_evap_tot(c) + qflx_surf(c)  + qflx_h2osfc_surf(c) &
                 + qflx_qrgwl(c) + qflx_drain(c) + qflx_drain_perched(c) + qflx_snwcp_ice(c) &
                 + qflx_lateral(c)

         else

            qflx_net_col(c) = 0.0_r8

         end if

      end do

      call c2g(bounds, endwb_col(bounds%begc:bounds%endc)             , &
           endwb_grc(bounds%begg:bounds%endg)                         , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, wa_local_col(bounds%begc:bounds%endc)          , &
           end_wa_grc(bounds%begg:bounds%endg)                        , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2ocan_col(bounds%begc:bounds%endc)            , &
           end_h2ocan_grc(bounds%begg:bounds%endg)                    , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osno_col(bounds%begc:bounds%endc)            , &
           end_h2osno_grc(bounds%begg:bounds%endg)                    , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osfc_col(bounds%begc:bounds%endc)            , &
           end_h2osfc_grc(bounds%begg:bounds%endg)                    , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osoi_liq_depth_intg(bounds%begc:bounds%endc) , &
           end_h2osoi_liq_grc(bounds%begg:bounds%endg)                , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, h2osoi_ice_depth_intg(bounds%begc:bounds%endc) , &
           end_h2osoi_ice_grc(bounds%begg:bounds%endg)                , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

      call c2g(bounds, errh2o(bounds%begc:bounds%endc)                , &
           errh2o_grc(bounds%begg:bounds%endg)                        , &
           c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    end associate

    end subroutine GridBalanceCheck

 end module BalanceCheckMod
