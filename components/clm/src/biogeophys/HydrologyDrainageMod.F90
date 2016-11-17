module HydrologyDrainageMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates soil/snow hydrology with drainage (subsurface runoff)
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varctl        , only : iulog, use_vichydro
  use clm_varcon        , only : e_ice, denh2o, denice, rpi, spval
  use atm2lndType       , only : atm2lnd_type
  use glc2lndMod        , only : glc2lnd_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : pft

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
    use clm_varcon       , only : denh2o, denice, secspday
    use clm_varctl       , only : glc_snow_persistence_max_days, use_vichydro, use_betr
    use clm_varpar       , only : nlevgrnd, nlevurb, nlevsoi    
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
    integer  :: g,l,c,j,fc                 ! indices
    real(r8) :: dtime                      ! land model time step (sec)
    !-----------------------------------------------------------------------
    
    associate(                                                         &    
         dz                 => col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer thickness depth (m)                       
         ctype              => col%itype                             , & ! Input:  [integer  (:)   ]  column type                                        

         qflx_floodg        => atm2lnd_vars%forc_flood_grc           , & ! Input:  [real(r8) (:)   ]  gridcell flux of flood water from RTM             
         forc_rain          => atm2lnd_vars%forc_rain_downscaled_col , & ! Input:  [real(r8) (:)   ]  rain rate [mm/s]                                  
         forc_snow          => atm2lnd_vars%forc_snow_downscaled_col , & ! Input:  [real(r8) (:)   ]  snow rate [mm/s]                                  

         glc_dyn_runoff_routing => glc2lnd_vars%glc_dyn_runoff_routing_grc,& ! Input:  [real(r8) (:)   ]  whether we're doing runoff routing appropriate for having a dynamic icesheet

         wa                 => soilhydrology_vars%wa_col             , & ! Input:  [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
         
         h2ocan             => waterstate_vars%h2ocan_col            , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O)                             
         h2osfc             => waterstate_vars%h2osfc_col            , & ! Input:  [real(r8) (:)   ]  surface water (mm)                                
         h2osno             => waterstate_vars%h2osno_col            , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                               
         begwb              => waterstate_vars%begwb_col             , & ! Input:  [real(r8) (:)   ]  water mass begining of the time step              
         endwb              => waterstate_vars%endwb_col             , & ! Output: [real(r8) (:)   ]  water mass end of the time step                   
         h2osoi_ice         => waterstate_vars%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
         h2osoi_liq         => waterstate_vars%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osoi_vol         => waterstate_vars%h2osoi_vol_col        , & ! Output: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         snow_persistence   => waterstate_vars%snow_persistence_col  , & ! Output: [real(r8) (:)   ]  counter for length of time snow-covered

         qflx_evap_tot      => waterflux_vars%qflx_evap_tot_col      , & ! Input:  [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg     
         qflx_irrig         => waterflux_vars%qflx_irrig_col         , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)                       
         qflx_glcice_melt   => waterflux_vars%qflx_glcice_melt_col   , & ! Input:  [real(r8) (:)]  ice melt (positive definite) (mm H2O/s)      
         qflx_h2osfc_surf   => waterflux_vars%qflx_h2osfc_surf_col   , & ! Output: [real(r8) (:)   ]  surface water runoff (mm/s)                        
         qflx_drain_perched => waterflux_vars%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ]  sub-surface runoff from perched zwt (mm H2O /s)   
         qflx_rsub_sat      => waterflux_vars%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]                 
         qflx_drain         => waterflux_vars%qflx_drain_col         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf          => waterflux_vars%qflx_surf_col          , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
         qflx_infl          => waterflux_vars%qflx_infl_col          , & ! Output: [real(r8) (:)   ]  infiltration (mm H2O /s)                          
         qflx_qrgwl         => waterflux_vars%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ]  qflx_surf at glaciers, wetlands, lakes            
         qflx_runoff        => waterflux_vars%qflx_runoff_col        , & ! Output: [real(r8) (:)   ]  total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_runoff_u      => waterflux_vars%qflx_runoff_u_col      , & ! Output: [real(r8) (:)   ]  Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
         qflx_runoff_r      => waterflux_vars%qflx_runoff_r_col      , & ! Output: [real(r8) (:)   ]  Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
         qflx_snwcp_ice     => waterflux_vars%qflx_snwcp_ice_col     , & ! Output: [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]`
         qflx_glcice        => waterflux_vars%qflx_glcice_col        , & ! Output: [real(r8) (:)   ]  flux of new glacier ice (mm H2O /s)               
         qflx_glcice_frz    => waterflux_vars%qflx_glcice_frz_col      & ! Output: [real(r8) (:)   ]  ice growth (positive definite) (mm H2O/s)         
         )

      ! Determine time step and step size

      dtime = get_step_size()

      if (use_vichydro) then
         call CLMVICMap(bounds, num_hydrologyc, filter_hydrologyc, &
              soilhydrology_vars, waterstate_vars)
      endif

      if (use_betr) then
        call ep_betr%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, waterstate_vars=waterstate_vars)
        call ep_betr%PreDiagSoilColWaterFlux(num_hydrologyc, filter_hydrologyc)
      endif
      if (.not. use_vsfm) then
         call Drainage(bounds, num_hydrologyc, filter_hydrologyc, &
              num_urbanc, filter_urbanc,&
              temperature_vars, soilhydrology_vars, soilstate_vars, &
              waterstate_vars, waterflux_vars)
      endif

      if (use_betr) then
        call ep_betr%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, waterstate_vars=waterstate_vars, &
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
         l = col%landunit(c)

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
            end if
         end do
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
         l = col%landunit(c)
         g = col%gridcell(c)
         ! In the following, we convert glc_snow_persistence_max_days to r8 to avoid overflow
         if ( (snow_persistence(c) >= (real(glc_snow_persistence_max_days, r8) * secspday)) &
              .or. lun%itype(l) == istice_mec) then
            qflx_glcice_frz(c) = qflx_snwcp_ice(c)  
            qflx_glcice(c) = qflx_glcice(c) + qflx_glcice_frz(c)
            if (glc_dyn_runoff_routing(g)) qflx_snwcp_ice(c) = 0._r8
         end if
      end do

      ! Determine wetland and land ice hydrology (must be placed here
      ! since need snow updated from CombineSnowLayers)

      do fc = 1,num_nolakec
         c = filter_nolakec(fc)
         l = col%landunit(c)
         g = col%gridcell(c)

         if (lun%itype(l)==istwet .or. lun%itype(l)==istice      &
                                  .or. lun%itype(l)==istice_mec) then

            qflx_drain(c)         = 0._r8
            qflx_drain_perched(c) = 0._r8
            qflx_h2osfc_surf(c)   = 0._r8
            qflx_surf(c)          = 0._r8
            qflx_infl(c)          = 0._r8
            qflx_qrgwl(c) = forc_rain(c) + forc_snow(c) + qflx_floodg(g) - qflx_evap_tot(c) - qflx_snwcp_ice(c) - &
                 (endwb(c)-begwb(c))/dtime

            ! With glc_dyn_runoff_routing = false (the less realistic way, typically used
            ! when NOT coupling to CISM), excess snow immediately runs off, whereas melting
            ! ice stays in place and does not run off. The reverse is true with
            ! glc_dyn_runoff_routing = true: in this case, melting ice runs off, and excess
            ! snow is sent to CISM, where it is converted to ice. These corrections are
            ! done here: 

            if (glc_dyn_runoff_routing(g) .and. lun%itype(l)==istice_mec) then
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

         else if (lun%urbpoi(l) .and. ctype(c) /= icol_road_perv) then

            qflx_drain_perched(c) = 0._r8
            qflx_h2osfc_surf(c)   = 0._r8
            qflx_rsub_sat(c)      = spval

         end if

         qflx_runoff(c) = qflx_drain(c) + qflx_surf(c)  + qflx_h2osfc_surf(c) + qflx_qrgwl(c) + qflx_drain_perched(c)

         if ((lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) .and. col%active(c)) then
            qflx_runoff(c) = qflx_runoff(c) - qflx_irrig(c)
         end if
         if (lun%urbpoi(l)) then
            qflx_runoff_u(c) = qflx_runoff(c)
         else if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
            qflx_runoff_r(c) = qflx_runoff(c)
         end if

      end do

    end associate

  end subroutine HydrologyDrainage

end module HydrologyDrainageMod
