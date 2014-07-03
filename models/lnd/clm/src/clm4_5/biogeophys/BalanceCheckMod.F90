module BalanceCheckMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water and energy balance check.
  !
  ! !USES:
  use shr_kind_mod       , only: r8 => shr_kind_r8
  use shr_log_mod        , only: errMsg => shr_log_errMsg
  use abortutils         , only: endrun
  use clm_varctl         , only: iulog
  use decompMod          , only: bounds_type
  use GetGlobalValuesMod , only: GetGlobalIndex

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
       num_hydrologyc, filter_hydrologyc)
    !
    ! !DESCRIPTION:
    ! Initialize column-level water balance at beginning of time step
    !
    ! !USES:
    use clmtype
    use clm_varpar   , only : nlevgrnd, nlevsoi, nlevurb
    use subgridAveMod, only : p2c
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall 
    use clm_varcon   , only : icol_road_perv, icol_road_imperv
    use clm_varcon   , only : denh2o, denice
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds     ! bounds
    integer, intent(in) :: num_nolakec          ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:)    ! column filter for non-lake points
    integer, intent(in) :: num_lakec            ! number of column non-lake points in column filter
    integer, intent(in) :: filter_lakec(:)      ! column filter for non-lake points
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    !
    ! !LOCAL VARIABLES:
    integer :: c, p, f, j, fc                  ! indices
    real(r8):: h2osoi_vol
    !-----------------------------------------------------------------------

   associate(& 
   h2osfc               =>    cws%h2osfc              , & ! Input:  [real(r8) (:)]  surface water (mm)                      
   ltype                =>    lun%itype               , & ! Input:  [integer (:)]  landunit type                            
   dz                   =>    cps%dz                  , & ! Output: [real(r8) (:,:)]                                        
   h2osno               =>    cws%h2osno              , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   h2osoi_ice           =>    cws%h2osoi_ice          , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                      
   h2osoi_liq           =>    cws%h2osoi_liq          , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                  
   begwb                =>    cwbal%begwb             , & ! Output: [real(r8) (:)]  water mass begining of the time step    
   h2ocan_col           =>    pws_a%h2ocan            , & ! Output: [real(r8) (:)]  canopy water (mm H2O) (column level)    
   wa                   =>    cws%wa                  , & ! Input:  [real(r8) (:)]  water in the unconfined aquifer (mm)    
   ctype                =>    col%itype               , & ! Input:  [integer (:)]  column type                              
   zwt                  =>    cws%zwt                 , & ! Input:  [real(r8) (:)]  water table depth (m)                   
   zi                   =>    cps%zi                  , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m) 
   h2ocan_pft           =>    pws%h2ocan                & ! Input:  [real(r8) (:)]  canopy water (mm H2O) (pft-level)       
   )

    ! Determine beginning water balance for time step
    ! pft-level canopy water averaged to column
    call p2c(bounds, num_nolakec, filter_nolakec, &
         h2ocan_pft(bounds%begp:bounds%endp), &
         h2ocan_col(bounds%begc:bounds%endc))

    do f = 1, num_hydrologyc
       c = filter_hydrologyc(f)
       if(zwt(c) <= zi(c,nlevsoi)) then
          wa(c) = 5000._r8
       end if
    end do

    do f = 1, num_nolakec
       c = filter_nolakec(f)
       if (ctype(c) == icol_roof .or. ctype(c) == icol_sunwall &
          .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_road_imperv) then
         begwb(c) = h2ocan_col(c) + h2osno(c)
       else
         begwb(c) = h2ocan_col(c) + h2osno(c) + h2osfc(c) + wa(c)
       end if

    end do
    do j = 1, nlevgrnd
      do f = 1, num_nolakec
         c = filter_nolakec(f)
         if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
          .or. ctype(c) == icol_roof) .and. j > nlevurb) then
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
   subroutine BalanceCheck(bounds, num_do_smb_c, filter_do_smb_c)
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
     use clmtype
     use subgridAveMod
     use clm_atmlnd       , only : clm_a2l, a2l_downscaled_col
     use clm_time_manager , only : get_step_size, get_nstep
     use clm_varcon       , only : icol_roof, icol_sunwall, icol_shadewall
     use clm_varcon       , only : spval, icol_road_perv, icol_road_imperv, istice_mec
     use clm_varcon       , only : istdlak, istsoil,istcrop,istwet
     use clm_varctl       , only : glc_dyn_runoff_routing, create_glacier_mec_landunit
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds    ! bounds
     integer, intent(in) :: num_do_smb_c        ! number of columns in filter_do_smb_c
     integer, intent(in) :: filter_do_smb_c (:) ! column filter for points where SMB calculations are done
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c,l,g,fc                  ! indices
     real(r8) :: dtime                       ! land model time step (sec)
     integer  :: nstep                       ! time step number
     logical  :: found                       ! flag in search loop
     integer  :: indexp,indexc,indexl,indexg ! index of first found in search loop
     real(r8) :: forc_rain_col(bounds%begc:bounds%endc)      ! column level rain rate [mm/s]
     real(r8) :: forc_snow_col(bounds%begc:bounds%endc)      ! column level snow rate [mm/s]
     !-----------------------------------------------------------------------

   associate(& 
   tws                  =>    grc%tws                       , & ! Input:  [real(r8) (:)] total water storage (mm H2O)             
   area                 =>    grc%area                      , & ! Input:  [real(r8) (:)] gridcell area (km2)                      
   volr                 =>    clm_a2l%volr                  , & ! Input:  [real(r8) (:)] river water storage (m3)                 
   do_capsnow           =>    cps%do_capsnow                , & ! Input:  [logical (:)]  true => do snow capping                  
   qflx_rain_grnd_col   =>    pwf_a%qflx_rain_grnd          , & ! Input:  [real(r8) (:)]  rain on ground after interception (mm H2O/s) [+]
   qflx_snow_grnd_col   =>    pwf_a%qflx_snow_grnd          , & ! Input:  [real(r8) (:)]  snow on ground after interception (mm H2O/s) [+]
   qflx_snow_h2osfc     =>    cwf%qflx_snow_h2osfc          , & ! Input:  [real(r8) (:)]  snow falling on surface water (mm/s)    
   frac_sno_eff         =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)]  effective snow fraction                 
   qflx_h2osfc_to_ice   =>    cwf%qflx_h2osfc_to_ice        , & ! Input:  [real(r8) (:)]  conversion of h2osfc to ice             
   frac_sno             =>    cps%frac_sno                  , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   qflx_drain_perched   =>    cwf%qflx_drain_perched        , & ! Input:  [real(r8) (:)]  sub-surface runoff (mm H2O /s)          
   qflx_floodc          =>    cwf%qflx_floodc               , & ! Input:  [real(r8) (:)]  total runoff due to flooding            
   qflx_evap_soi        =>    pwf_a%qflx_evap_soi           , & ! Input:  [real(r8) (:)]  soil evaporation (mm H2O/s) (+ = to atm)
   qflx_h2osfc_surf     =>    cwf%qflx_h2osfc_surf          , & ! Input:  [real(r8) (:)] surface water runoff (mm/s)              
   qflx_snow_melt       =>    cwf%qflx_snow_melt            , & ! Input:  [real(r8) (:)]  snow melt (net)                         
   sabg_soil            =>    pef%sabg_soil                 , & ! Input:  [real(r8) (:)]  solar radiation absorbed by soil (W/m**2)
   sabg_snow            =>    pef%sabg_snow                 , & ! Input:  [real(r8) (:)]  solar radiation absorbed by snow (W/m**2)
   sabg_chk             =>    pef%sabg_chk                  , & ! Input:  [real(r8) (:)]  sum of soil/snow using current fsno, for balance check
   forc_rain            =>    a2l_downscaled_col%forc_rain  , & ! Input:  [real(r8) (:)]  rain rate [mm/s]
   forc_snow            =>    a2l_downscaled_col%forc_snow  , & ! Input:  [real(r8) (:)]  snow rate [mm/s]
   forc_lwrad           =>    a2l_downscaled_col%forc_lwrad , & ! Input:  [real(r8) (:)]  downward infrared (longwave) radiation (W/m**2)
   forc_solad           =>    clm_a2l%forc_solad            , & ! Input:  [real(r8) (:,:)]  direct beam radiation (vis=forc_sols , nir=forc_soll )
   forc_solai           =>    clm_a2l%forc_solai            , & ! Input:  [real(r8) (:,:)]  diffuse radiation     (vis=forc_solsd, nir=forc_solld)
   ltype                =>    lun%itype                     , & ! Input:  [integer (:)]  landunit type                            
   canyon_hwr           =>    lun%canyon_hwr                , & ! Input:  [real(r8) (:)]  ratio of building height to street width
   urbpoi               =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point       
   cactive              =>    col%active                    , & ! Input:  [logical (:)]  true=>do computations on this column 
   ctype                =>    col%itype                     , & ! Input:  [integer (:)]  column type                              
   endwb                =>    cwbal%endwb                   , & ! Input:  [real(r8) (:)]  water mass end of the time step         
   begwb                =>    cwbal%begwb                   , & ! Input:  [real(r8) (:)]  water mass begining of the time step    
   qflx_irrig           =>    pwf_a%qflx_irrig              , & ! Input:  [real(r8) (:)]  irrigation flux (mm H2O /s)             
   qflx_surf            =>    cwf%qflx_surf                 , & ! Input:  [real(r8) (:)]  surface runoff (mm H2O /s)              
   qflx_qrgwl           =>    cwf%qflx_qrgwl                , & ! Input:  [real(r8) (:)]  qflx_surf at glaciers, wetlands, lakes  
   qflx_drain           =>    cwf%qflx_drain                , & ! Input:  [real(r8) (:)]  sub-surface runoff (mm H2O /s)          
   qflx_runoff          =>    cwf%qflx_runoff               , & ! Input:  [real(r8) (:)]  total runoff (mm H2O /s)                
   qflx_snwcp_ice       =>    pwf_a%qflx_snwcp_ice          , & ! Input:  [real(r8) (:)]  excess snowfall due to snow capping (mm H2O /s) [+]`
   qflx_evap_tot        =>    pwf_a%qflx_evap_tot           , & ! Input:  [real(r8) (:)]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
   qflx_glcice_frz      =>    cwf%qflx_glcice_frz           , & ! Input:  [real(r8) (:)]  ice growth (mm H2O/s) [+]               
   qflx_glcice_melt     =>    cwf%qflx_glcice_melt          , & ! Input:  [real(r8) (:)]  ice melt (mm H2O/s) [             
   h2osno               =>    cws%h2osno                    , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                     
   h2osno_old           =>    cws%h2osno_old                , & ! Input:  [real(r8) (:)]  snow water (mm H2O) at previous time step
   qflx_dew_snow        =>    pwf_a%qflx_dew_snow           , & ! Input:  [real(r8) (:)]  surface dew added to snow pack (mm H2O /s) [+]
   qflx_sub_snow        =>    pwf_a%qflx_sub_snow           , & ! Input:  [real(r8) (:)]  sublimation rate from snow pack (mm H2O /s) [+]
   qflx_top_soil        =>    cwf%qflx_top_soil             , & ! Input:  [real(r8) (:)]  net water input into soil from top (mm/s)
   qflx_evap_grnd       =>    pwf_a%qflx_evap_grnd          , & ! Input:  [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]
   qflx_dew_grnd        =>    pwf_a%qflx_dew_grnd           , & ! Input:  [real(r8) (:)]  ground surface dew formation (mm H2O /s) [+]
   qflx_prec_grnd       =>    pwf_a%qflx_prec_grnd          , & ! Input:  [real(r8) (:)]  water onto ground including canopy runoff [kg/(m2 s)]
   qflx_snwcp_liq       =>    pwf_a%qflx_snwcp_liq          , & ! Input:  [real(r8) (:)]  excess liquid water due to snow capping (mm H2O /s) [+]`
   qflx_sl_top_soil     =>    cwf%qflx_sl_top_soil          , & ! Input:  [real(r8) (:)]  liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
   snl                  =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   pactive              =>    pft%active                    , & ! Input:  [logical (:)]  true=>do computations on this pft 
   fsa                  =>    pef%fsa                       , & ! Input:  [real(r8) (:)]  solar radiation absorbed (total) (W/m**2)
   fsr                  =>    pef%fsr                       , & ! Input:  [real(r8) (:)]  solar radiation reflected (W/m**2)      
   eflx_lwrad_out       =>    pef%eflx_lwrad_out            , & ! Input:  [real(r8) (:)]  emitted infrared (longwave) radiation (W/m**2)
   eflx_lwrad_net       =>    pef%eflx_lwrad_net            , & ! Input:  [real(r8) (:)]  net infrared (longwave) rad (W/m**2) [+ = to atm]
   sabv                 =>    pef%sabv                      , & ! Input:  [real(r8) (:)]  solar radiation absorbed by vegetation (W/m**2)
   sabg                 =>    pef%sabg                      , & ! Input:  [real(r8) (:)]  solar radiation absorbed by ground (W/m**2)
   eflx_sh_tot          =>    pef%eflx_sh_tot               , & ! Input:  [real(r8) (:)]  total sensible heat flux (W/m**2) [+ to atm]
   eflx_lh_tot          =>    pef%eflx_lh_tot               , & ! Input:  [real(r8) (:)]  total latent heat flux (W/m8*2)  [+ to atm]
   eflx_soil_grnd       =>    pef%eflx_soil_grnd            , & ! Input:  [real(r8) (:)]  soil heat flux (W/m**2) [+ = into soil] 
   eflx_wasteheat_pft   =>    pef%eflx_wasteheat_pft        , & ! Input:  [real(r8) (:)]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
   eflx_heat_from_ac_pft=>    pef%eflx_heat_from_ac_pft     , & ! Input:  [real(r8) (:)] sensible heat flux put back into canyon due to removal by AC (W/m**2)
   eflx_traffic_pft     =>    pef%eflx_traffic_pft          , & ! Input:  [real(r8) (:)]  traffic sensible heat flux (W/m**2)     
   qflx_runoffg         =>    gwf%qflx_runoffg              , & ! Input:  [real(r8) (:)]  total runoff at gridcell level inc land cover change flux (mm H2O /s)
   qflx_liq_dynbal      =>    gwf%qflx_liq_dynbal           , & ! Input:  [real(r8) (:)]  liq runoff due to dynamic land cover change (mm H2O /s)
   qflx_snwcp_iceg      =>    gwf%qflx_snwcp_iceg           , & ! Input:  [real(r8) (:)]  excess snowfall due to snow cap inc land cover change flux (mm H20/s)
   qflx_ice_dynbal      =>    gwf%qflx_ice_dynbal           , & ! Input:  [real(r8) (:)]  ice runoff due to dynamic land cover change (mm H2O /s)
   eflx_sh_totg         =>    gef%eflx_sh_totg              , & ! Input:  [real(r8) (:)]  total sensible heat flux at grid level (W/m**2) [+ to atm]
   eflx_dynbal          =>    gef%eflx_dynbal               , & ! Input:  [real(r8) (:)]  energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]
   snow_sources         =>    cws%snow_sources              , & ! Output: [real(r8) (:)]  snow sources (mm H2O /s)                
   snow_sinks           =>    cws%snow_sinks                , & ! Output: [real(r8) (:)]  snow sinks (mm H2O /s)                  
   errh2o               =>    cwbal%errh2o                  , & ! Output: [real(r8) (:)]  water conservation error (mm H2O)       
   errsoi_col           =>    cebal%errsoi                  , & ! Output: [real(r8) (:)]  column-level soil/lake energy conservation error (W/m**2)
   errh2osno            =>    cws%errh2osno                 , & ! Output: [real(r8) (:)]  error in h2osno (kg m-2)                
   errsol               =>    pebal%errsol                  , & ! Output: [real(r8) (:)]  solar radiation conservation error (W/m**2)
   errseb               =>    pebal%errseb                  , & ! Output: [real(r8) (:)]  surface energy conservation error (W/m**2)
   errlon               =>    pebal%errlon                  , & ! Output: [real(r8) (:)]  longwave radiation conservation error (W/m**2)
   netrad               =>    pef%netrad                      & ! Output: [real(r8) (:)]  net radiation (positive downward) (W/m**2)
   )

    ! Get step size and time step

    nstep = get_nstep()
    dtime = get_step_size()

    ! Determine column level incoming snow and rain
    ! Assume no incident precipitation on urban wall columns (as in Hydrology1Mod.F90).

    do c = bounds%begc,bounds%endc
       g = col%gridcell(c)
       l = col%landunit(c)       
        
       if (ctype(c) == icol_sunwall .or.  ctype(c) == icol_shadewall) then
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
       if (cactive(c))then
          errh2o(c) = endwb(c) - begwb(c) &
               - (forc_rain_col(c) + forc_snow_col(c)  + qflx_floodc(c) + qflx_irrig(c) &
                 - qflx_evap_tot(c) - qflx_surf(c)  - qflx_h2osfc_surf(c) &
                 - qflx_qrgwl(c) - qflx_drain(c) - qflx_drain_perched(c) - qflx_snwcp_ice(c)) * dtime

       else

          errh2o(c) = 0.0_r8

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

    if (glc_dyn_runoff_routing) then
       do fc = 1,num_do_smb_c
          c = filter_do_smb_c(fc)
          errh2o(c) = errh2o(c) + qflx_glcice_frz(c)*dtime
          errh2o(c) = errh2o(c) - qflx_glcice_melt(c)*dtime
       end do
    endif
    
    found = .false.
    do c = bounds%begc, bounds%endc
       if (abs(errh2o(c)) > 1e-7_r8) then
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

       if ((ctype(indexc) .eq. icol_roof .or. &
            ctype(indexc) .eq. icol_road_imperv .or. &
            ctype(indexc) .eq. icol_road_perv) .and. &
            abs(errh2o(indexc)) > 1.e-1 .and. (nstep > 2) ) then

          write(iulog,*)'clm urban model is stopping - error is greater than .10 (mm)'
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

       else if (abs(errh2o(indexc)) > .10_r8 .and. (nstep > 2) ) then

          write(iulog,*)'clm model is stopping - error is greater than .10 (mm)'
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
          write(iulog,*)'icemask=' , grc%icemask(col%gridcell(indexc))
          write(iulog,*)'clm model is stopping'
          call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
       end if
    end if

    ! Snow balance check

    do c = bounds%begc,bounds%endc
       if (cactive(c)) then
          g = col%gridcell(c)
          l = col%landunit(c)
          ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
          ! any given time step but only if there is at least one snow layer.  h2osno 
          ! also includes snow that is part of the soil column (an initial snow layer is 
          ! only created if h2osno > 10mm).
          if (snl(c) .lt. 0) then
             snow_sources(c) = qflx_prec_grnd(c) + qflx_dew_snow(c) + qflx_dew_grnd(c)
             snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c) + qflx_snow_melt(c) &
                  + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) + qflx_sl_top_soil(c)

             if (ltype(l) == istdlak) then 
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

                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c) ) &
                        + qflx_snow_melt(c)  + qflx_sl_top_soil(c)
                endif
             endif

             if (ltype(l) == istsoil .or. ltype(l) == istcrop .or. ltype(l) == istwet ) then
                if ( do_capsnow(c) ) then
                   snow_sources(c) = frac_sno_eff(c) * (qflx_dew_snow(c) + qflx_dew_grnd(c) ) &
                        + qflx_h2osfc_to_ice(c) + qflx_prec_grnd(c)

                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c)) &
                        + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                        + qflx_snow_melt(c) + qflx_sl_top_soil(c)
                else
                   snow_sources(c) = (qflx_snow_grnd_col(c) - qflx_snow_h2osfc(c) ) &
                        + frac_sno_eff(c) * (qflx_rain_grnd_col(c) &
                        +  qflx_dew_snow(c) + qflx_dew_grnd(c) ) + qflx_h2osfc_to_ice(c)

                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c)) &
                        + qflx_snow_melt(c) + qflx_sl_top_soil(c)
                endif
             endif
             
             if (glc_dyn_runoff_routing) then
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
       if (cactive(c)) then
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
            ' ctype= ',ctype(indexc), &
            ' ltype= ',ltype(col%landunit(indexc)), &
            ' errh2osno= ',errh2osno(indexc)
       
       if (abs(errh2osno(indexc)) > 0.1_r8 .and. (nstep > 2) ) then
          write(iulog,*)'clm model is stopping - error is greater than .10 (mm)'
          write(iulog,*)'nstep            = ',nstep
          write(iulog,*)'errh2osno        = ',errh2osno(indexc)
          write(iulog,*)'snl              = ',snl(indexc)
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
       if (pactive(p)) then
          c = pft%column(p)
          l = pft%landunit(p)
          g = pft%gridcell(p)

          ! Solar radiation energy balance
          ! Do not do this check for an urban pft since it will not balance on a per-column
          ! level because of interactions between columns and since a separate check is done
          ! in the urban radiation module
          if (.not. urbpoi(l)) then
             errsol(p) = fsa(p) + fsr(p) &
                  - (forc_solad(g,1) + forc_solad(g,2) + forc_solai(g,1) + forc_solai(g,2))
          else
             errsol(p) = spval
          end if
          
          ! Longwave radiation energy balance
          ! Do not do this check for an urban pft since it will not balance on a per-column
          ! level because of interactions between columns and since a separate check is done
          ! in the urban radiation module
          if (.not. urbpoi(l)) then
             errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(c)
          else
             errlon(p) = spval
          end if
          
          ! Surface energy balance
          ! Changed to using (eflx_lwrad_net) here instead of (forc_lwrad - eflx_lwrad_out) because
          ! there are longwave interactions between urban columns (and therefore pfts). 
          ! For surfaces other than urban, (eflx_lwrad_net) equals (forc_lwrad - eflx_lwrad_out),
          ! and a separate check is done above for these terms.
          
          if (.not. urbpoi(l)) then
             errseb(p) = sabv(p) + sabg_chk(p) + forc_lwrad(c) - eflx_lwrad_out(p) &
                         - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p)
          else
             errseb(p) = sabv(p) + sabg(p) &
                         - eflx_lwrad_net(p) &
                         - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) &
                         + eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
          end if
          netrad(p) = fsa(p) - eflx_lwrad_net(p)
       end if
    end do

    ! Solar radiation energy balance check

    found = .false.
    do p = bounds%begp, bounds%endp
       if (pactive(p)) then
          if ( (errsol(p) /= spval) .and. (abs(errsol(p)) > .10_r8) ) then
             found = .true.
             indexp = p
             indexg = pft%gridcell(indexp)
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,*)'BalanceCheck: solar radiation balance error (W/m2)'
       write(iulog,*)'nstep         = ',nstep
       write(iulog,*)'errsol        = ',errsol(p)
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
       if (pactive(p)) then
          if ( (errlon(p) /= spval) .and. (abs(errlon(p)) > .10_r8) ) then
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
       if (pactive(p)) then
          if (abs(errseb(p)) > .10_r8 ) then
             found = .true.
             indexp = p
             indexc = pft%column(indexp)
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,*)'BalanceCheck: surface flux energy balance error (W/m2)'
       write(iulog,*)'nstep          = ',nstep
       write(iulog,*)'errseb         = ',errseb(indexp)
       write(iulog,*)'sabv           = ',sabv(indexp)
       write(iulog,*)'sabg           = ',sabg(indexp), ((1._r8- frac_sno(indexc))*sabg_soil(indexp) + &
                                         frac_sno(indexc)*sabg_snow(indexp)),sabg_chk(indexp)
       write(iulog,*)'eflx_lwrad_net = ',eflx_lwrad_net(indexp)
       write(iulog,*)'eflx_sh_tot    = ',eflx_sh_tot(indexp)
       write(iulog,*)'eflx_lh_tot    = ',eflx_lh_tot(indexp)
       write(iulog,*)'eflx_soil_grnd = ',eflx_soil_grnd(indexp)
       write(iulog,*)'clm model is stopping'
       call endrun(decomp_index=indexp, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
    end if

    ! Soil energy balance check

    found = .false.
    do c = bounds%begc,bounds%endc
       if (cactive(c)) then
          if (abs(errsoi_col(c)) > 1.0e-7_r8 ) then
             found = .true.
             indexc = c
          end if
       end if
    end do
    if ( found ) then
       if (abs(errsoi_col(indexc)) > .10_r8 .and. (nstep > 2) ) then
          write(iulog,*)'BalanceCheck: soil balance error (mm)'
          write(iulog,*)'nstep         = ',nstep
          write(iulog,*)'errsoi_col    = ',errsoi_col(indexc)
          write(iulog,*)'clm model is stopping'
          call endrun(decomp_index=indexc, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
       end if
    end if

    ! Update SH and RUNOFF for dynamic land cover change energy and water fluxes
    call c2g( bounds, &
              qflx_runoff(bounds%begc:bounds%endc), qflx_runoffg(bounds%begg:bounds%endg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       qflx_runoffg(g) = qflx_runoffg(g) - qflx_liq_dynbal(g)
    enddo

    call c2g( bounds, &
              qflx_snwcp_ice(bounds%begc:bounds%endc), qflx_snwcp_iceg(bounds%begg:bounds%endg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       qflx_snwcp_iceg(g) = qflx_snwcp_iceg(g) - qflx_ice_dynbal(g)
    enddo

    call p2g( bounds, &
              eflx_sh_tot(bounds%begp:bounds%endp), eflx_sh_totg(bounds%begg:bounds%endg), &
              p2c_scale_type='unity',c2l_scale_type='urbanf',l2g_scale_type='unity')
    do g = bounds%begg, bounds%endg
       eflx_sh_totg(g) =  eflx_sh_totg(g) - eflx_dynbal(g)
    enddo

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb
    call c2g( bounds, &
         endwb(bounds%begc:bounds%endc), tws(bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    ! second add river storage as gridcell average depth
    ! 1.e-3 converts [m3/km2] to [mm]
    do g = bounds%begg, bounds%endg
       tws(g) = tws(g) + volr(g) / area(g) * 1.e-3_r8
    enddo

    end associate 
   end subroutine BalanceCheck

end module BalanceCheckMod
