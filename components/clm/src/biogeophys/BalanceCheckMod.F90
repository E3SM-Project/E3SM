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
  use clm_varctl         , only : iulog, do_varsoil
  use clm_varcon         , only : namep, namec
  use GetGlobalValuesMod , only : GetGlobalIndex
  use atm2lndType        , only : atm2lnd_type
  use glc2lndMod         , only : glc2lnd_type
  use EnergyFluxType     , only : energyflux_type
  use SolarAbsorbedType  , only : solarabs_type
  use SoilHydrologyType  , only : soilhydrology_type  
  use WaterstateType     , only : waterstate_type
  use WaterfluxType      , only : waterflux_type
  use GridcellType       , only : grc                
  use LandunitType       , only : lun                
  use ColumnType         , only : col                
  use PatchType          , only : pft                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginWaterBalance  ! Initialize water balance check
  public :: BalanceCheck       ! Water and energy balance check
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BeginWaterBalance(bounds, &
       num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
       num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Initialize column-level water balance at beginning of time step
    !
    ! !USES:
    use subgridAveMod , only : p2c
    use clm_varpar    , only : nlevgrnd, nlevsoi, nlevurb
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
    integer :: c, p, f, j, fc                  ! indices
    real(r8):: h2osoi_vol
    !-----------------------------------------------------------------------

    associate(                                               & 
         zi           =>    col%zi                         , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m) 

         h2ocan_patch =>    waterstate_vars%h2ocan_patch   , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O) (pft-level)       
         h2osfc       =>    waterstate_vars%h2osfc_col     , & ! Input:  [real(r8) (:)   ]  surface water (mm)                      
         h2osno       =>    waterstate_vars%h2osno_col     , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
         h2osoi_ice   =>    waterstate_vars%h2osoi_ice_col , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                      
         h2osoi_liq   =>    waterstate_vars%h2osoi_liq_col , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                  

         zwt          =>    soilhydrology_vars%zwt_col     , & ! Input:  [real(r8) (:)   ]  water table depth (m)                   
         wa           =>    soilhydrology_vars%wa_col      , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)    
         
         h2ocan_col   =>    waterstate_vars%h2ocan_col     , & ! Output: [real(r8) (:)   ]  canopy water (mm H2O) (column level)    
         begwb        =>    waterstate_vars%begwb_col        & ! Output: [real(r8) (:)   ]  water mass begining of the time step    
         )

      ! Determine beginning water balance for time step
      ! pft-level canopy water averaged to column

      call p2c(bounds, num_nolakec, filter_nolakec, &
           h2ocan_patch(bounds%begp:bounds%endp), &
           h2ocan_col(bounds%begc:bounds%endc))

      do f = 1, num_hydrologyc
         c = filter_hydrologyc(f)
         if(do_varsoil) then
	   wa(c) = 0._r8               ! Made 0 for variable soil thickness
	 else
           if(zwt(c) <= zi(c,nlevsoi)) then
             wa(c) = 5000._r8
	   end if
         end if
      end do

      do f = 1, num_nolakec
         c = filter_nolakec(f)
         if (col%itype(c) == icol_roof .or. col%itype(c) == icol_sunwall &
              .or. col%itype(c) == icol_shadewall .or. col%itype(c) == icol_road_imperv) then
            begwb(c) = h2ocan_col(c) + h2osno(c)
         else
            begwb(c) = h2ocan_col(c) + h2osno(c) + h2osfc(c) + wa(c)
         end if

      end do
      do j = 1, nlevgrnd
         do f = 1, num_nolakec
            c = filter_nolakec(f)
            if ((col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                 .or. col%itype(c) == icol_roof) .and. j > nlevurb) then
            else
               begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
            end if
         end do
      end do

      do f = 1, num_lakec
         c = filter_lakec(f)
         begwb(c) = h2osno(c)
      end do

    end associate 

   end subroutine BeginWaterBalance

   !-----------------------------------------------------------------------
   subroutine BalanceCheck( bounds, num_do_smb_c, filter_do_smb_c, &
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
     use clm_varcon        , only : spval
     use column_varcon     , only : icol_roof, icol_sunwall, icol_shadewall
     use column_varcon     , only : icol_road_perv, icol_road_imperv
     use landunit_varcon   , only : istice_mec, istdlak, istsoil,istcrop,istwet
     use clm_varctl        , only : create_glacier_mec_landunit
     use clm_time_manager  , only : get_step_size, get_nstep
     use clm_initializeMod , only : surfalb_vars
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
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c,l,g,fc                             ! indices
     real(r8) :: dtime                                  ! land model time step (sec)
     integer  :: nstep                                  ! time step number
     logical  :: found                                  ! flag in search loop
     integer  :: indexp,indexc,indexl,indexg            ! index of first found in search loop
     real(r8) :: forc_rain_col(bounds%begc:bounds%endc) ! column level rain rate [mm/s]
     real(r8) :: forc_snow_col(bounds%begc:bounds%endc) ! column level snow rate [mm/s]
     !-----------------------------------------------------------------------

     associate(                                                                   & 
          volr                    =>    atm2lnd_vars%volr_grc                   , & ! Input:  [real(r8) (:)   ]  river water storage (m3)                 
          forc_solad              =>    atm2lnd_vars%forc_solad_grc             , & ! Input:  [real(r8) (:,:) ]  direct beam radiation (vis=forc_sols , nir=forc_soll )
          forc_solai              =>    atm2lnd_vars%forc_solai_grc             , & ! Input:  [real(r8) (:,:) ]  diffuse radiation     (vis=forc_solsd, nir=forc_solld)
          forc_rain               =>    atm2lnd_vars%forc_rain_downscaled_col   , & ! Input:  [real(r8) (:)   ]  rain rate [mm/s]
          forc_snow               =>    atm2lnd_vars%forc_snow_downscaled_col   , & ! Input:  [real(r8) (:)   ]  snow rate [mm/s]
          forc_lwrad              =>    atm2lnd_vars%forc_lwrad_downscaled_col  , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)

          glc_dyn_runoff_routing  =>    glc2lnd_vars%glc_dyn_runoff_routing_grc , & ! Input:  [real(r8) (:)   ]  whether we're doing runoff routing appropriate for having a dynamic icesheet

          do_capsnow              =>    waterstate_vars%do_capsnow_col          , & ! Input:  [logical (:)    ]  true => do snow capping                  
          h2osno                  =>    waterstate_vars%h2osno_col              , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
          h2osno_old              =>    waterstate_vars%h2osno_old_col          , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O) at previous time step
          frac_sno_eff            =>    waterstate_vars%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  effective snow fraction                 
          frac_sno                =>    waterstate_vars%frac_sno_col            , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          begwb                   =>    waterstate_vars%begwb_col               , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step    
          errh2o                  =>    waterstate_vars%errh2o_col              , & ! Output: [real(r8) (:)   ]  water conservation error (mm H2O)       
          errh2osno               =>    waterstate_vars%errh2osno_col           , & ! Output: [real(r8) (:)   ]  error in h2osno (kg m-2)                
          endwb                   =>    waterstate_vars%endwb_col               , & ! Output: [real(r8) (:)   ]  water mass end of the time step         
          
          dwb                     =>    waterflux_vars%dwb_col                  , & ! Output: [real(r8) (:)   ]  change of water mass within the time step [kg/m2/s]
          qflx_rain_grnd_col      =>    waterflux_vars%qflx_rain_grnd_col       , & ! Input:  [real(r8) (:)   ]  rain on ground after interception (mm H2O/s) [+]
          qflx_snow_grnd_col      =>    waterflux_vars%qflx_snow_grnd_col       , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          qflx_evap_soi           =>    waterflux_vars%qflx_evap_soi_col        , & ! Input:  [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
          qflx_irrig              =>    waterflux_vars%qflx_irrig_col           , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)             
          qflx_snwcp_ice          =>    waterflux_vars%qflx_snwcp_ice_col       , & ! Input:  [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]`
          qflx_evap_tot           =>    waterflux_vars%qflx_evap_tot_col        , & ! Input:  [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
          qflx_dew_snow           =>    waterflux_vars%qflx_dew_snow_col        , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
          qflx_sub_snow           =>    waterflux_vars%qflx_sub_snow_col        , & ! Input:  [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]
          qflx_evap_grnd          =>    waterflux_vars%qflx_evap_grnd_col       , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
          qflx_dew_grnd           =>    waterflux_vars%qflx_dew_grnd_col        , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
          qflx_prec_grnd          =>    waterflux_vars%qflx_prec_grnd_col       , & ! Input:  [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
          qflx_snwcp_liq          =>    waterflux_vars%qflx_snwcp_liq_col       , & ! Input:  [real(r8) (:)   ]  excess liquid water due to snow capping (mm H2O /s) [+]`
          qflx_snow_h2osfc        =>    waterflux_vars%qflx_snow_h2osfc_col     , & ! Input:  [real(r8) (:)   ]  snow falling on surface water (mm/s)    
          qflx_h2osfc_to_ice      =>    waterflux_vars%qflx_h2osfc_to_ice_col   , & ! Input:  [real(r8) (:)   ]  conversion of h2osfc to ice             
          qflx_drain_perched      =>    waterflux_vars%qflx_drain_perched_col   , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_floodc             =>    waterflux_vars%qflx_floodc_col          , & ! Input:  [real(r8) (:)   ]  total runoff due to flooding            
          qflx_h2osfc_surf        =>    waterflux_vars%qflx_h2osfc_surf_col     , & ! Input:  [real(r8) (:)   ]  surface water runoff (mm/s)              
          qflx_snow_melt          =>    waterflux_vars%qflx_snow_melt_col       , & ! Input:  [real(r8) (:)   ]  snow melt (net)                         
          qflx_surf               =>    waterflux_vars%qflx_surf_col            , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)              
          qflx_qrgwl              =>    waterflux_vars%qflx_qrgwl_col           , & ! Input:  [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes  
          qflx_drain              =>    waterflux_vars%qflx_drain_col           , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)          
          qflx_runoff             =>    waterflux_vars%qflx_runoff_col          , & ! Input:  [real(r8) (:)   ]  total runoff (mm H2O /s)                
          qflx_glcice             =>    waterflux_vars%qflx_glcice_col          , & ! Input:  [real(r8) (:)   ]  flux of new glacier ice (mm H2O /s) [+ if ice grows]
          qflx_glcice_melt        =>    waterflux_vars%qflx_glcice_melt_col     , & ! Input:  [real(r8) (:)   ]  ice melt (mm H2O/s)              
          qflx_glcice_frz         =>    waterflux_vars%qflx_glcice_frz_col      , & ! Input:  [real(r8) (:)   ]  ice growth (mm H2O/s) [+]               
          qflx_top_soil           =>    waterflux_vars%qflx_top_soil_col        , & ! Input:  [real(r8) (:)   ]  net water input into soil from top (mm/s)
          qflx_sl_top_soil        =>    waterflux_vars%qflx_sl_top_soil_col     , & ! Input:  [real(r8) (:)   ]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
          qflx_liq_dynbal         =>    waterflux_vars%qflx_liq_dynbal_grc      , & ! Input:  [real(r8) (:)   ]  liq runoff due to dynamic land cover change (mm H2O /s)
          qflx_ice_dynbal         =>    waterflux_vars%qflx_ice_dynbal_grc      , & ! Input:  [real(r8) (:)   ]  ice runoff due to dynamic land cover change (mm H2O /s)
          snow_sources            =>    waterflux_vars%snow_sources_col         , & ! Output: [real(r8) (:)   ]  snow sources (mm H2O /s)  
          snow_sinks              =>    waterflux_vars%snow_sinks_col           , & ! Output: [real(r8) (:)   ]  snow sinks (mm H2O /s)    

          eflx_lwrad_out          =>    energyflux_vars%eflx_lwrad_out_patch    , & ! Input:  [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)
          eflx_lwrad_net          =>    energyflux_vars%eflx_lwrad_net_patch    , & ! Input:  [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm]
          eflx_sh_tot             =>    energyflux_vars%eflx_sh_tot_patch       , & ! Input:  [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
          eflx_lh_tot             =>    energyflux_vars%eflx_lh_tot_patch       , & ! Input:  [real(r8) (:)   ]  total latent heat flux (W/m8*2)  [+ to atm]
          eflx_soil_grnd          =>    energyflux_vars%eflx_soil_grnd_patch    , & ! Input:  [real(r8) (:)   ]  soil heat flux (W/m**2) [+ = into soil] 
          eflx_wasteheat_patch    =>    energyflux_vars%eflx_wasteheat_patch    , & ! Input:  [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
          eflx_heat_from_ac_patch =>    energyflux_vars%eflx_heat_from_ac_patch , & ! Input:  [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
          eflx_traffic_patch      =>    energyflux_vars%eflx_traffic_patch      , & ! Input:  [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)     
          eflx_dynbal             =>    energyflux_vars%eflx_dynbal_grc         , & ! Input:  [real(r8) (:)   ]  energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]

          sabg_soil               =>    solarabs_vars%sabg_soil_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by soil (W/m**2)
          sabg_snow               =>    solarabs_vars%sabg_snow_patch           , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by snow (W/m**2)
          sabg_chk                =>    solarabs_vars%sabg_chk_patch            , & ! Input:  [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
          fsa                     =>    solarabs_vars%fsa_patch                 , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)
          fsr                     =>    solarabs_vars%fsr_patch                 , & ! Input:  [real(r8) (:)   ]  solar radiation reflected (W/m**2)      
          sabv                    =>    solarabs_vars%sabv_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)
          sabg                    =>    solarabs_vars%sabg_patch                , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
          
          errsoi_col              =>    energyflux_vars%errsoi_col              , & ! Output: [real(r8) (:)   ]  column-level soil/lake energy conservation error (W/m**2)
          errsol                  =>    energyflux_vars%errsol_patch            , & ! Output: [real(r8) (:)   ]  solar radiation conservation error (W/m**2)
          errseb                  =>    energyflux_vars%errseb_patch            , & ! Output: [real(r8) (:)   ]  surface energy conservation error (W/m**2)
          errlon                  =>    energyflux_vars%errlon_patch            , & ! Output: [real(r8) (:)   ]  longwave radiation conservation error (W/m**2)

          fabd                    =>    surfalb_vars%fabd_patch                 , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit direct flux
          fabi                    =>    surfalb_vars%fabi_patch                 , & ! Input:  [real(r8) (:,:)]  flux absorbed by canopy per unit indirect flux
          elai                    =>    canopystate_vars%elai_patch             , & ! Input:  [real(r8) (:,:)]  
          esai                    =>    canopystate_vars%esai_patch             , & ! Input:  [real(r8) (:,:)]  

          albd                    =>    surfalb_vars%albd_patch                 , & ! Output: [real(r8) (:,:)]  surface albedo (direct)
          albi                    =>    surfalb_vars%albi_patch                 , & ! Output: [real(r8) (:,:)]  surface albedo (diffuse)
          ftdd                    =>    surfalb_vars%ftdd_patch                , & ! Input:  [real(r8) (:,:)]  down direct flux below canopy per unit direct flux
          ftid                    =>    surfalb_vars%ftid_patch                , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit direct flux
          ftii                    =>    surfalb_vars%ftii_patch                 , & ! Input:  [real(r8) (:,:)]  down diffuse flux below canopy per unit diffuse flux

          netrad                  =>    energyflux_vars%netrad_patch              & ! Output: [real(r8) (:)   ]  net radiation (positive downward) (W/m**2)
          )

       ! Get step size and time step

       nstep = get_nstep()
       dtime = get_step_size()

       ! Determine column level incoming snow and rain
       ! Assume no incident precipitation on urban wall columns (as in CanopyHydrologyMod.F90).

       do c = bounds%begc,bounds%endc
          g = col%gridcell(c)
          l = col%landunit(c)       

          if (col%itype(c) == icol_sunwall .or.  col%itype(c) == icol_shadewall) then
             forc_rain_col(c) = 0.
             forc_snow_col(c) = 0.
          else
             forc_rain_col(c) = forc_rain(c)
             forc_snow_col(c) = forc_snow(c)
          end if
       end do

       ! Water balance check

       do c = bounds%begc, bounds%endc

          ! add qflx_drain_perched and qflx_flood
          if (col%active(c)) then

             errh2o(c) = endwb(c) - begwb(c) &
                  - (forc_rain_col(c) + forc_snow_col(c)  + qflx_floodc(c) + qflx_irrig(c) &
                  - qflx_evap_tot(c) - qflx_surf(c)  - qflx_h2osfc_surf(c) &
                  - qflx_qrgwl(c) - qflx_drain(c) - qflx_drain_perched(c) - qflx_snwcp_ice(c)) * dtime
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
          g = col%gridcell(c)
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
               ' global indexc= ',GetGlobalIndex(decomp_index=indexc, clmlevel=namec), &
               ' errh2o= ',errh2o(indexc)

          if ((col%itype(indexc) == icol_roof .or. &
               col%itype(indexc) == icol_road_imperv .or. &
               col%itype(indexc) == icol_road_perv) .and. &
               abs(errh2o(indexc)) > 1.e-4_r8 .and. (nstep > 2) ) then

             write(iulog,*)'clm urban model is stopping - error is greater than 1e-4 (mm)'
             write(iulog,*)'nstep          = ',nstep
             write(iulog,*)'errh2o         = ',errh2o(indexc)
             write(iulog,*)'forc_rain      = ',forc_rain_col(indexc)
             write(iulog,*)'forc_snow      = ',forc_snow_col(indexc)
             write(iulog,*)'endwb          = ',endwb(indexc)
             write(iulog,*)'begwb          = ',begwb(indexc)
             write(iulog,*)'qflx_evap_tot  = ',qflx_evap_tot(indexc)
             write(iulog,*)'qflx_irrig     = ',qflx_irrig(indexc)
             write(iulog,*)'qflx_surf      = ',qflx_surf(indexc)
             write(iulog,*)'qflx_qrgwl     = ',qflx_qrgwl(indexc)
             write(iulog,*)'qflx_drain     = ',qflx_drain(indexc)
             write(iulog,*)'qflx_snwcp_ice = ',qflx_snwcp_ice(indexc)
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))

          else if (abs(errh2o(indexc)) > 1.e-4_r8 .and. (nstep > 2) ) then

             write(iulog,*)'clm model is stopping - error is greater than 1e-4 (mm)'
             write(iulog,*)'nstep              = ',nstep
             write(iulog,*)'errh2o             = ',errh2o(indexc)
             write(iulog,*)'forc_rain          = ',forc_rain_col(indexc)
             write(iulog,*)'forc_snow          = ',forc_snow_col(indexc)
             write(iulog,*)'endwb              = ',endwb(indexc)
             write(iulog,*)'begwb              = ',begwb(indexc)
             write(iulog,*)'qflx_evap_tot      = ',qflx_evap_tot(indexc)
             write(iulog,*)'qflx_irrig         = ',qflx_irrig(indexc)
             write(iulog,*)'qflx_surf          = ',qflx_surf(indexc)
             write(iulog,*)'qflx_h2osfc_surf   = ',qflx_h2osfc_surf(indexc)
             write(iulog,*)'qflx_qrgwl         = ',qflx_qrgwl(indexc)
             write(iulog,*)'qflx_drain         = ',qflx_drain(indexc)
             write(iulog,*)'qflx_drain_perched = ',qflx_drain_perched(indexc)
             write(iulog,*)'qflx_flood         = ',qflx_floodc(indexc)
             write(iulog,*)'qflx_snwcp_ice     = ',qflx_snwcp_ice(indexc)
             write(iulog,*)'qflx_glcice_melt   = ',qflx_glcice_melt(indexc)
             write(iulog,*)'qflx_glcice_frz    = ',qflx_glcice_frz(indexc)
             write(iulog,*)'clm model is stopping'
             call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

       ! Snow balance check

       do c = bounds%begc,bounds%endc
          if (col%active(c)) then
             g = col%gridcell(c)
             l = col%landunit(c)

             ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
             ! any given time step but only if there is at least one snow layer.  h2osno 
             ! also includes snow that is part of the soil column (an initial snow layer is 
             ! only created if h2osno > 10mm).

             if (col%snl(c) < 0) then
                snow_sources(c) = qflx_prec_grnd(c) + qflx_dew_snow(c) + qflx_dew_grnd(c)
                snow_sinks(c)  = qflx_sub_snow(c) + qflx_evap_grnd(c) + qflx_snow_melt(c) &
                     + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) + qflx_sl_top_soil(c)

                if (lun%itype(l) == istdlak) then 
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

                if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop .or. lun%itype(l) == istwet ) then
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
          if (col%active(c)) then
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
               ' global indexc= ',GetGlobalIndex(decomp_index=indexc, clmlevel=namec), &
               ' col%itype= ',col%itype(indexc), &
               ' lun%itype= ',lun%itype(col%landunit(indexc)), &
               ' errh2osno= ',errh2osno(indexc)

          if (abs(errh2osno(indexc)) > 1.e-4_r8 .and. (nstep > 2) ) then
             write(iulog,*)'clm model is stopping - error is greater than 1e-4 (mm)'
             write(iulog,*)'nstep            = ',nstep
             write(iulog,*)'errh2osno        = ',errh2osno(indexc)
             write(iulog,*)'snl              = ',col%snl(indexc)
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
          if (pft%active(p)) then
             c = pft%column(p)
             l = pft%landunit(p)
             g = pft%gridcell(p)

             ! Solar radiation energy balance
             ! Do not do this check for an urban pft since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun%urbpoi(l)) then
                errsol(p) = fsa(p) + fsr(p) &
                     - (forc_solad(g,1) + forc_solad(g,2) + forc_solai(g,1) + forc_solai(g,2))
             else
                errsol(p) = spval
             end if

             ! Longwave radiation energy balance
             ! Do not do this check for an urban pft since it will not balance on a per-column
             ! level because of interactions between columns and since a separate check is done
             ! in the urban radiation module
             if (.not. lun%urbpoi(l)) then
                errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(c)
             else
                errlon(p) = spval
             end if

             ! Surface energy balance
             ! Changed to using (eflx_lwrad_net) here instead of (forc_lwrad - eflx_lwrad_out) because
             ! there are longwave interactions between urban columns (and therefore patches). 
             ! For surfaces other than urban, (eflx_lwrad_net) equals (forc_lwrad - eflx_lwrad_out),
             ! and a separate check is done above for these terms.

             if (.not. lun%urbpoi(l)) then
                errseb(p) = sabv(p) + sabg_chk(p) + forc_lwrad(c) - eflx_lwrad_out(p) &
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
          if (pft%active(p)) then
             if ( (errsol(p) /= spval) .and. (abs(errsol(p)) > 1.e-3_r8) ) then
                found = .true.
                indexp = p
                indexg = pft%gridcell(indexp)
             end if
          end if
       end do
       if ( found  .and. (nstep > 2) ) then
          write(iulog,*)'BalanceCheck: solar radiation balance error (W/m2)'
          write(iulog,*)'nstep         = ',nstep
          write(iulog,*)'errsol        = ',errsol(indexp)
          write(iulog,*)'fsa           = ',fsa(indexp)
          write(iulog,*)'fsr           = ',fsr(indexp)
          write(iulog,*)'forc_solad(1) = ',forc_solad(indexg,1)
          write(iulog,*)'forc_solad(2) = ',forc_solad(indexg,2)
          write(iulog,*)'forc_solai(1) = ',forc_solai(indexg,1)
          write(iulog,*)'forc_solai(2) = ',forc_solai(indexg,2)
          write(iulog,*)'forc_tot      = ',forc_solad(indexg,1)+forc_solad(indexg,2) &
               +forc_solai(indexg,1)+forc_solai(indexg,2)
          write(iulog,*)'clm model is stopping'
          call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
       end if

       ! Longwave radiation energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (pft%active(p)) then
             if ( (errlon(p) /= spval) .and. (abs(errlon(p)) > 1.e-3_r8) ) then
                found = .true.
                indexp = p
             end if
          end if
       end do
       if ( found  .and. (nstep > 2) ) then
          write(iulog,*)'BalanceCheck: longwave energy balance error (W/m2)'
          write(iulog,*)'nstep        = ',nstep
          write(iulog,*)'errlon       = ',errlon(indexp)
          call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
       end if

       ! Surface energy balance check

       found = .false.
       do p = bounds%begp, bounds%endp
          if (pft%active(p)) then
             if (abs(errseb(p)) > 1.e-3_r8 ) then
                found = .true.
                indexp = p
                indexc = pft%column(indexp)
             end if
          end if
       end do
       if ( found  .and. (nstep > 2) ) then
          write(iulog,*)'BalanceCheck: surface flux energy balance error (W/m2)'
          write(iulog,*)'nstep          = ' ,nstep
          write(iulog,*)'errseb         = ' ,errseb(indexp)
          write(iulog,*)'sabv           = ' ,sabv(indexp)

          write(iulog,*)'sabg           = ' ,sabg(indexp), ((1._r8- frac_sno(indexc))*sabg_soil(indexp) + &
               frac_sno(indexc)*sabg_snow(indexp)),sabg_chk(indexp)

          write(iulog,*)'forc_tot      = '  ,forc_solad(indexg,1) + forc_solad(indexg,2) + &
               forc_solai(indexg,1) + forc_solai(indexg,2)

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
          !
          ! FIX(RF, 082914) - commented out in ED branch
          !
          ! call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
       end if

       ! Soil energy balance check

       found = .false.
       do c = bounds%begc,bounds%endc
          if (col%active(c)) then
             if (abs(errsoi_col(c)) > 1.0e-7_r8 ) then
                found = .true.
                indexc = c
             end if
          end if
       end do
       if ( found ) then
          if (abs(errsoi_col(indexc)) > 1.e-3_r8 .and. (nstep > 2) ) then
             write(iulog,*)'BalanceCheck: soil balance error (mm)'
             write(iulog,*)'nstep         = ',nstep
             write(iulog,*)'errsoi_col    = ',errsoi_col(indexc)
             write(iulog,*)'clm model is stopping'
             !
             ! FIX(RF, 082914) - commented out in ED branch
             !
             ! call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
          end if
       end if

     end associate

   end subroutine BalanceCheck

end module BalanceCheckMod
