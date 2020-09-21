module HydrologyDrainageMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates soil/snow hydrology with drainage (subsurface runoff)
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varctl        , only : iulog, use_vichydro
  use elm_varcon        , only : e_ice, denh2o, denice, rpi, spval
  use atm2lndType       , only : atm2lnd_type
  use glc2lndMod        , only : glc2lnd_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use SoilHydrologyMod  , only : WaterTable
  use TopounitDataType  , only : top_af ! atmospheric flux variables
  use LandunitType      , only : lun_pp                
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_ws, col_wf  
  use VegetationType    , only : veg_pp                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: HydrologyDrainage ! Calculates soil/snow hydrolog with drainage
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine HydrologyDrainage(bounds,    &
       num_nolakec, filter_nolakec,       &
       num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc,         &
       num_do_smb_c, filter_do_smb_c,     &
       atm2lnd_vars, glc2lnd_vars, temperature_vars,    &
       soilhydrology_vars, soilstate_vars, waterstate_vars, waterflux_vars, ep_betr)
    !
    ! !DESCRIPTION:
    ! Calculates soil/snow hydrology with drainage (subsurface runoff)
    !
    ! !USES:
    use landunit_varcon  , only : istice, istwet, istsoil, istice_mec, istcrop
    use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, icol_shadewall
    use elm_varcon       , only : denh2o, denice, secspday
    use clm_varctl       , only : glc_snow_persistence_max_days, use_vichydro, use_betr
    use domainMod        , only : ldomain
    use atm2lndType      , only : atm2lnd_type
    use elm_varpar       , only : nlevgrnd, nlevurb, nlevsoi    
    use clm_time_manager , only : get_step_size, get_nstep
    use SoilHydrologyMod , only : CLMVICMap, Drainage
    use clm_varctl       , only : use_vsfm
    use BeTRSimulationALM, only : betr_simulation_alm_type
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds               
    integer                  , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
    integer                  , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    integer                  , intent(in)    :: num_do_smb_c         ! number of bareland columns in which SMB is calculated, in column filter    
    integer                  , intent(in)    :: filter_do_smb_c(:)   ! column filter for bare land SMB columns      
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(glc2lnd_type)       , intent(in)    :: glc2lnd_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    class(betr_simulation_alm_type), intent(inout) :: ep_betr
    !
    ! !LOCAL VARIABLES:
    integer  :: g,t,l,c,j,fc               ! indices
    real(r8) :: dtime                      ! land model time step (sec)
    !-----------------------------------------------------------------------
    
    associate(                                                                  &    
         dz                     => col_pp%dz                                     , & ! Input:  [real(r8) (:,:) ]  layer thickness depth (m)                       
         ctype                  => col_pp%itype                                  , & ! Input:  [integer  (:)   ]  column type                                        

         qflx_floodg            => atm2lnd_vars%forc_flood_grc                , & ! Input:  [real(r8) (:)   ]  gridcell flux of flood water from RTM             
         forc_rain              => top_af%rain                                , & ! Input:  [real(r8) (:)   ]  rain rate (kg H2O/m**2/s, or mm liquid H2O/s)                                  
         forc_snow              => top_af%snow                                , & ! Input:  [real(r8) (:)   ]  snow rate (kg H2O/m**2/s, or mm liquid H2O/s)                                  

         glc_dyn_runoff_routing => glc2lnd_vars%glc_dyn_runoff_routing_grc    , & ! Input:  [real(r8) (:)   ]  whether we're doing runoff routing appropriate for having a dynamic icesheet

         wa                     => soilhydrology_vars%wa_col                  , & ! Input:  [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
         
         h2osoi_liq_depth_intg  => col_ws%h2osoi_liq_depth_intg , & ! Output: [real(r8) (:)   ]  grid-level depth integrated liquid soil water
         h2osoi_ice_depth_intg  => col_ws%h2osoi_ice_depth_intg , & ! Output: [real(r8) (:)   ]  grid-level depth integrated ice soil water
         h2ocan                 => col_ws%h2ocan                 , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O)                             
         h2osfc                 => col_ws%h2osfc                 , & ! Input:  [real(r8) (:)   ]  surface water (mm)                                
         h2osno                 => col_ws%h2osno                 , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                               
         begwb                  => col_ws%begwb                  , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step              
         endwb                  => col_ws%endwb                  , & ! Output: [real(r8) (:)   ]  water mass end of the time step                   
         h2osoi_ice             => col_ws%h2osoi_ice             , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
         h2osoi_liq             => col_ws%h2osoi_liq             , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osoi_vol             => col_ws%h2osoi_vol             , & ! Output: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         snow_persistence       => col_ws%snow_persistence       , & ! Output: [real(r8) (:)   ]  counter for length of time snow-covered
         total_plant_stored_h2o => col_ws%total_plant_stored_h2o , & ! Input [real(r8) (:) dynamic water stored in plants]
         qflx_evap_tot          => col_wf%qflx_evap_tot           , & ! Input:  [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg     
         qflx_irrig             => col_wf%qflx_irrig              , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)
         qflx_irr_demand        => col_wf%qflx_irr_demand         , & ! Input:  [real(r8) (:)   ]  irrigation demand sent to MOSART/WM (mm H2O /s)                       		 
         qflx_glcice_melt       => col_wf%qflx_glcice_melt        , & ! Input:  [real(r8) (:)]  ice melt (positive definite) (mm H2O/s)      
         qflx_h2osfc_surf       => col_wf%qflx_h2osfc_surf        , & ! Output: [real(r8) (:)   ]  surface water runoff (mm/s)                        
         qflx_drain_perched     => col_wf%qflx_drain_perched      , & ! Output: [real(r8) (:)   ]  sub-surface runoff from perched zwt (mm H2O /s)   
         qflx_rsub_sat          => col_wf%qflx_rsub_sat           , & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]                 
         qflx_drain             => col_wf%qflx_drain              , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf              => col_wf%qflx_surf               , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
         qflx_infl              => col_wf%qflx_infl               , & ! Output: [real(r8) (:)   ]  infiltration (mm H2O /s)                          
         qflx_qrgwl             => col_wf%qflx_qrgwl              , & ! Output: [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes            
         qflx_runoff            => col_wf%qflx_runoff             , & ! Output: [real(r8) (:)   ]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_runoff_u          => col_wf%qflx_runoff_u           , & ! Output: [real(r8) (:)   ]  Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
         qflx_runoff_r          => col_wf%qflx_runoff_r           , & ! Output: [real(r8) (:)   ]  Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_snwcp_ice         => col_wf%qflx_snwcp_ice          , & ! Output: [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]`
         qflx_glcice            => col_wf%qflx_glcice             , & ! Output: [real(r8) (:)   ]  flux of new glacier ice (mm H2O /s)               
         qflx_glcice_frz        => col_wf%qflx_glcice_frz           & ! Output: [real(r8) (:)   ]  ice growth (positive definite) (mm H2O/s)         
         )

      ! Determine time step and step size

      dtime = get_step_size()

      if (use_vichydro) then
         call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc, &
              soilhydrology_vars, waterstate_vars)
      endif

      if (use_betr) then
        call ep_betr%BeTRSetBiophysForcing(bounds, col_pp, veg_pp, 1, nlevsoi, waterstate_vars=waterstate_vars)
        call ep_betr%PreDiagSoilColWaterFlux(num_hydrologyc, filter_hydrologyc)
      endif

      if (.not. use_vsfm) then
         call Drainage(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc,&
              temperature_vars, soilhydrology_vars, soilstate_vars, &
              waterstate_vars, waterflux_vars)
      endif

      if (use_betr) then
        call ep_betr%BeTRSetBiophysForcing(bounds, col_pp, veg_pp, 1, nlevsoi, waterstate_vars=waterstate_vars, &
          waterflux_vars=waterflux_vars)
        call ep_betr%DiagDrainWaterFlux(num_hydrologyc, filter_hydrologyc)
        call ep_betr%RetrieveBiogeoFlux(bounds, 1, nlevsoi, waterflux_vars=waterflux_vars)
      endif

      do j = 1, nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof) .and. j > nlevurb) then
            else
               h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
            end if
         end do
      end do

      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         l = col_pp%landunit(c)

         if (ctype(c) == icol_roof .or. ctype(c) == icol_sunwall &
              .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_road_imperv) then
            endwb(c) = h2ocan(c) + h2osno(c)
         else
            ! add h2osfc to water balance
            endwb(c) = h2ocan(c) + h2osno(c) + h2osfc(c) + wa(c)

         end if
      end do

      do j = 1, nlevgrnd
         do fc = 1, num_nolakec
            c = filter_nolakec(fc)
            if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
                 .or. ctype(c) == icol_roof) .and. j > nlevurb) then

            else
               endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
               h2osoi_liq_depth_intg(c) = h2osoi_liq_depth_intg(c) + h2osoi_liq(c,j)
               h2osoi_ice_depth_intg(c) = h2osoi_ice_depth_intg(c) + h2osoi_ice(c,j)
            end if
         end do
      end do

      ! ---------------------------------------------------------------------------------
      ! Add stored plant water to the column water balance
      ! currently, stored plant water is only dynamic when FATES is turned on.
      ! Other orthogonal modules should not need to worry about this term,
      ! and it should be zero in all other cases and all other columns.
      ! ---------------------------------------------------------------------------------
      do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         endwb(c) = endwb(c) + total_plant_stored_h2o(c)
      end do
      
      ! Prior to summing up wetland/ice hydrology, calculate land ice contributions/sinks
      ! to this hydrology.
      ! 1) Generate SMB from capped-snow amount.  This is done over istice_mec
      !    columns, and also any other columns included in do_smb_c filter, where
      !    perennial snow has remained for at least snow_persistence_max.
      ! 2) If using glc_dyn_runoff_routing=T, zero qflx_snwcp_ice: qflx_snwcp_ice is the flux
      !    sent to ice runoff, but for glc_dyn_runoff_routing=T, we do NOT want this to be
      !    sent to ice runoff (instead it is sent to CISM).

      do c = bounds%begc,bounds%endc
         qflx_glcice_frz(c) = 0._r8
      end do
      do fc = 1,num_do_smb_c
         c = filter_do_smb_c(fc)
         l = col_pp%landunit(c)
         g = col_pp%gridcell(c)
         ! In the following, we convert glc_snow_persistence_max_days to r8 to avoid overflow
         if ( (snow_persistence(c) >= (real(glc_snow_persistence_max_days, r8) * secspday)) &
              .or. lun_pp%itype(l) == istice_mec) then
            qflx_glcice_frz(c) = qflx_snwcp_ice(c)  
            qflx_glcice(c) = qflx_glcice(c) + qflx_glcice_frz(c)
            if (glc_dyn_runoff_routing(g)) qflx_snwcp_ice(c) = 0._r8
         end if
      end do

      ! Determine wetland and land ice hydrology (must be placed here
      ! since need snow updated from CombineSnowLayers)
      do c = bounds%begc,bounds%endc
         qflx_irr_demand(c) = 0._r8
      end do

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col_pp%landunit(c)
         t = col_pp%topounit(c)
         g = col_pp%gridcell(c)

         if (lun_pp%itype(l)==istwet .or. lun_pp%itype(l)==istice      &
                                  .or. lun_pp%itype(l)==istice_mec) then

            qflx_drain(c)         = 0._r8
            qflx_drain_perched(c) = 0._r8
            qflx_h2osfc_surf(c)   = 0._r8
            qflx_surf(c)          = 0._r8
            qflx_infl(c)          = 0._r8
            qflx_qrgwl(c) = forc_rain(t) + forc_snow(t) + qflx_floodg(g) - qflx_evap_tot(c) - qflx_snwcp_ice(c) - &
                 (endwb(c)-begwb(c))/dtime

            ! With glc_dyn_runoff_routing = false (the less realistic way, typically used
            ! when NOT coupling to CISM), excess snow immediately runs off, whereas melting
            ! ice stays in place and does not run off. The reverse is true with
            ! glc_dyn_runoff_routing = true: in this case, melting ice runs off, and excess
            ! snow is sent to CISM, where it is converted to ice. These corrections are
            ! done here: 

            if (glc_dyn_runoff_routing(g) .and. lun_pp%itype(l)==istice_mec) then
               ! If glc_dyn_runoff_routing=T, add meltwater from istice_mec ice columns to the runoff.
               !    Note: The meltwater contribution is computed in PhaseChanges (part of Biogeophysics2)
               qflx_qrgwl(c) = qflx_qrgwl(c) + qflx_glcice_melt(c)
               ! Also subtract the freezing component of qflx_glcice: this ice is added to
               ! CISM's ice column rather than running off. (This is analogous to the
               ! subtraction of qflx_snwcp_ice from qflx_qrgwl above, which accounts for
               ! snow that should be put into ice runoff rather than liquid runoff. But for
               ! glc_dyn_runoff_routing=true, qflx_snwcp_ice has been zeroed out, and has
               ! been put into qflx_glcice_frz.)
               qflx_qrgwl(c) = qflx_qrgwl(c) - qflx_glcice_frz(c)
            endif

         else if (lun_pp%urbpoi(l) .and. ctype(c) /= icol_road_perv) then

            qflx_drain_perched(c) = 0._r8
            qflx_h2osfc_surf(c)   = 0._r8
            qflx_rsub_sat(c)      = spval

         end if

         qflx_runoff(c) = qflx_drain(c) + qflx_surf(c)  + qflx_h2osfc_surf(c) + qflx_qrgwl(c) + qflx_drain_perched(c)

         if ((lun_pp%itype(l)==istsoil .or. lun_pp%itype(l)==istcrop) .and. col_pp%active(c)) then
            qflx_irr_demand(c) = -1.0_r8 * ldomain%f_surf(g)*qflx_irrig(c) !surface water demand send to MOSART																																		 
         end if
         if (lun_pp%urbpoi(l)) then
            qflx_runoff_u(c) = qflx_runoff(c)
         else if (lun_pp%itype(l)==istsoil .or. lun_pp%itype(l)==istcrop) then
            qflx_runoff_r(c) = qflx_runoff(c)
         end if

      end do
																					  
    end associate

  end subroutine HydrologyDrainage

end module HydrologyDrainageMod
