module Biogeophysics2Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs the calculation of soil/snow and ground temperatures
  ! and updates surface fluxes based on the new ground temperature.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl,   only: iulog
  use decompMod   , only: bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Biogeophysics2   ! Calculate soil/snow and ground temperatures
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Biogeophysics2 (bounds, &
       num_urbanl, filter_urbanl, num_nolakec, filter_nolakec, &
       num_nolakep, filter_nolakep)
    !
    ! !DESCRIPTION:
    ! This is the main subroutine to execute the calculation of soil/snow and
    ! ground temperatures and update surface fluxes based on the new ground
    ! temperature
    !
    ! Calling sequence is:
    ! Biogeophysics2:             surface biogeophysics driver
    !    -> SoilTemperature:      soil/snow and ground temperatures
    !          -> SoilTermProp    thermal conductivities and heat capacities
    !          -> Tridiagonal     tridiagonal matrix solution
    !          -> PhaseChange     phase change of liquid/ice contents
    !
    ! (1) Snow and soil temperatures
    !     o The volumetric heat capacity is calculated as a linear combination
    !       in terms of the volumetric fraction of the constituent phases.
    !     o The thermal conductivity of soil is computed from
    !       the algorithm of Johansen (as reported by Farouki 1981), and the
    !       conductivity of snow is from the formulation used in
    !       SNTHERM (Jordan 1991).
    !     o Boundary conditions:
    !       F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
    !     o Soil / snow temperature is predicted from heat conduction
    !       in 10 soil layers and up to 5 snow layers.
    !       The thermal conductivities at the interfaces between two
    !       neighboring layers (j, j+1) are derived from an assumption that
    !       the flux across the interface is equal to that from the node j
    !       to the interface and the flux from the interface to the node j+1.
    !       The equation is solved using the Crank-Nicholson method and
    !       results in a tridiagonal system equation.
    !
    ! (2) Phase change (see PhaseChange.F90)
    !
    ! !USES:
    use clmtype
    use clm_atmlnd        , only : a2l_downscaled_col
    use clm_time_manager  , only : get_step_size
    use clm_varcon        , only : hvap, cpair, grav, vkc, tfrz, sb, icol_road_perv, &
                                   icol_roof, icol_sunwall, icol_shadewall, istsoil
    use clm_varcon        , only : istcrop
    use clm_varpar        , only : nlevsno, nlevgrnd, nlevurb, max_pft_per_col
    use SoilTemperatureMod, only : SoilTemperature
    use subgridAveMod     , only : p2c
    use perf_mod          , only : t_startf, t_stopf
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds    ! bounds
    integer, intent(in) :: num_nolakec         ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:)   ! column filter for non-lake points
    integer, intent(in) :: num_urbanl          ! number of urban landunits in clump
    integer, intent(in) :: filter_urbanl(:)    ! urban landunit filter
    integer, intent(in) :: num_nolakep         ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(:)   ! pft filter for non-lake points
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,g,j,pi,l         ! indices
    integer  :: fc,fp                ! lake filtered column and pft indices
    real(r8) :: dtime                ! land model time step (sec)
    real(r8) :: egsmax(bounds%begc:bounds%endc)      ! max. evaporation which soil can provide at one time step
    real(r8) :: egirat(bounds%begc:bounds%endc)      ! ratio of topsoil_evap_tot : egsmax
    real(r8) :: tinc(bounds%begc:bounds%endc)        ! temperature difference of two time step
    real(r8) :: xmf(bounds%begc:bounds%endc)         ! total latent heat of phase change of ground water
    real(r8) :: sumwt(bounds%begc:bounds%endc)       ! temporary
    real(r8) :: evaprat(bounds%begp:bounds%endp)     ! ratio of qflx_evap_soi/topsoil_evap_tot
    real(r8) :: save_qflx_evap_soi   ! temporary storage for qflx_evap_soi
    real(r8) :: topsoil_evap_tot(bounds%begc:bounds%endc)          ! column-level total evaporation from top soil layer
    real(r8) :: fact(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd) ! used in computing tridiagonal matrix
    real(r8) :: eflx_lwrad_del(bounds%begp:bounds%endp)            ! update due to eflx_lwrad
    real(r8) :: t_grnd0(bounds%begc:bounds%endc)    !t_grnd of previous time step
    real(r8) :: c_h2osfc(bounds%begc:bounds%endc)   !heat capacity of surface water
    real(r8) :: xmf_h2osfc(bounds%begc:bounds%endc) !latent heat of phase change of surface water
    real(r8) :: eflx_temp(bounds%begc:bounds%endc)
    real(r8) :: xmf_temp(bounds%begc:bounds%endc)
    real(r8) :: dq_temp(bounds%begc:bounds%endc)
    real(r8) :: lw_grnd
    real(r8) :: fsno_eff
    !-----------------------------------------------------------------------

   associate(& 
   forc_lwrad                =>    a2l_downscaled_col%forc_lwrad , & ! Input:  [real(r8) (:)]  downward infrared (longwave) radiation (W/m**2)
   ltype                     =>    lun%itype                     , & ! Input:  [integer (:)]  landunit type                            
   urbpoi                    =>    lun%urbpoi                    , & ! Input:  [logical (:)]  true => landunit is an urban point       
   canyon_hwr                =>   lun%canyon_hwr                 , & ! Input:  [real(r8) (:)]  ratio of building height to street width (-)
   frac_sno_eff              =>    cps%frac_sno_eff              , & ! Input:  [real(r8) (:)] eff. fraction of ground covered by snow (0 to 1)
   frac_sno                  =>    cps%frac_sno                  , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   h2osfc                    =>    cws%h2osfc                    , & ! Input:  [real(r8) (:)]  surface water (mm)                      
   frac_h2osfc               =>    cps%frac_h2osfc               , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   t_h2osfc                  =>    ces%t_h2osfc                  , & ! Input:  [real(r8) (:)]  surface water temperature               
   t_h2osfc_bef              =>    ces%t_h2osfc_bef              , & ! Input:  [real(r8) (:)]  saved surface water temperature         
   qflx_ev_snow              =>    pwf%qflx_ev_snow              , & ! Input:  [real(r8) (:)]  evaporation flux from snow (W/m**2) [+ to atm]
   qflx_ev_soil              =>    pwf%qflx_ev_soil              , & ! Input:  [real(r8) (:)]  evaporation flux from soil (W/m**2) [+ to atm]
   qflx_ev_h2osfc            =>    pwf%qflx_ev_h2osfc            , & ! Input:  [real(r8) (:)]  evaporation flux from soil (W/m**2) [+ to atm]
   sabg_soil                 =>    pef%sabg_soil                 , & ! Input:  [real(r8) (:)]  solar radiation absorbed by soil (W/m**2)
   sabg_snow                 =>    pef%sabg_snow                 , & ! Input:  [real(r8) (:)]  solar radiation absorbed by snow (W/m**2)
   ctype                     =>    col%itype                     , & ! Input:  [integer (:)]  column type                              
   npfts                     =>   col%npfts                      , & ! Input:  [integer (:)]  column's number of pfts                  
   pfti                      =>   col%pfti                       , & ! Input:  [integer (:)]  column's beginning pft index             
   snl                       =>    cps%snl                       , & ! Input:  [integer (:)]  number of snow layers                    
   do_capsnow                =>    cps%do_capsnow                , & ! Input:  [logical (:)]  true => do snow capping                  
   htvp                      =>    cps%htvp                      , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
   emg                       =>    cps%emg                       , & ! Input:  [real(r8) (:)]  ground emissivity                       
   t_grnd                    =>    ces%t_grnd                    , & ! Input:  [real(r8) (:)]  ground temperature (Kelvin)             
   dt_grnd                   =>    ces%dt_grnd                   , & ! Output: [real(r8) (:)]  change in t_grnd, last iteration (Kelvin)
   t_soisno                  =>    ces%t_soisno                  , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)             
   tssbef                    =>    ces%tssbef                    , & ! Input:  [real(r8) (:,:)]  soil/snow temperature before update   
   h2osoi_ice                =>    cws%h2osoi_ice                , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2) (new)                
   h2osoi_liq                =>    cws%h2osoi_liq                , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2) (new)            
   errsoi_col                =>    cebal%errsoi                  , & ! Output: [real(r8) (:)]  column-level soil/lake energy conservation error (W/m**2)
   eflx_building_heat        =>    cef%eflx_building_heat        , & ! Input:  [real(r8) (:)]  heat flux from urban building interior to walls, roof
   pactive                   =>    pft%active                    , & ! Input:  [logical (:)]  true=>do computations on this pft 
   pcolumn                   =>   pft%column                     , & ! Input:  [integer (:)]  pft's column index                       
   plandunit                 =>   pft%landunit                   , & ! Input:  [integer (:)]  pft's landunit index                     
   pgridcell                 =>   pft%gridcell                   , & ! Input:  [integer (:)]  pft's gridcell index                     
   frac_veg_nosno            =>    pps%frac_veg_nosno            , & ! Input:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1 now) [-]
   sabg                      =>    pef%sabg                      , & ! Input:  [real(r8) (:)]  solar radiation absorbed by ground (W/m**2)
   dlrad                     =>    pef%dlrad                     , & ! Input:  [real(r8) (:)]  downward longwave radiation below the canopy [W/m2]
   ulrad                     =>    pef%ulrad                     , & ! Input:  [real(r8) (:)]  upward longwave radiation above the canopy [W/m2]
   eflx_sh_grnd              =>    pef%eflx_sh_grnd              , & ! Input:  [real(r8) (:)]  sensible heat flux from ground (W/m**2) [+ to atm]
   eflx_sh_veg               =>    pef%eflx_sh_veg               , & ! Input:  [real(r8) (:)]  sensible heat flux from leaves (W/m**2) [+ to atm]
   qflx_evap_soi             =>    pwf%qflx_evap_soi             , & ! Input:  [real(r8) (:)]  soil evaporation (mm H2O/s) (+ = to atm)
   qflx_evap_veg             =>    pwf%qflx_evap_veg             , & ! Input:  [real(r8) (:)]  vegetation evaporation (mm H2O/s) (+ = to atm)
   qflx_tran_veg             =>    pwf%qflx_tran_veg             , & ! Input:  [real(r8) (:)]  vegetation transpiration (mm H2O/s) (+ = to atm)
   qflx_evap_can             =>    pwf%qflx_evap_can             , & ! Input:  [real(r8) (:)]  evaporation from leaves and stems (mm H2O/s) (+ = to atm)
   qflx_snwcp_liq            =>    pwf%qflx_snwcp_liq            , & ! Input:  [real(r8) (:)]  excess rainfall due to snow capping (mm H2O /s)
   qflx_snwcp_ice            =>    pwf%qflx_snwcp_ice            , & ! Input:  [real(r8) (:)]  excess snowfall due to snow capping (mm H2O /s)
   qflx_evap_tot             =>    pwf%qflx_evap_tot             , & ! Output: [real(r8) (:)]  qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
   qflx_evap_grnd            =>    pwf%qflx_evap_grnd            , & ! Output: [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]
   qflx_sub_snow             =>    pwf%qflx_sub_snow             , & ! Output: [real(r8) (:)]  sublimation rate from snow pack (mm H2O /s) [+]
   qflx_dew_snow             =>    pwf%qflx_dew_snow             , & ! Output: [real(r8) (:)]  surface dew added to snow pack (mm H2O /s) [+]
   qflx_dew_grnd             =>    pwf%qflx_dew_grnd             , & ! Output: [real(r8) (:)]  ground surface dew formation (mm H2O /s) [+]
   eflx_soil_grnd            =>    pef%eflx_soil_grnd            , & ! Output: [real(r8) (:)]  soil heat flux (W/m**2) [+ = into soil] 
   eflx_soil_grnd_u          =>    pef%eflx_soil_grnd_u          , & ! Output: [real(r8) (:)]  urban soil heat flux (W/m**2) [+ = into soil]
   eflx_soil_grnd_r          =>    pef%eflx_soil_grnd_r          , & ! Output: [real(r8) (:)]  rural soil heat flux (W/m**2) [+ = into soil]
   eflx_sh_tot               =>    pef%eflx_sh_tot               , & ! Output: [real(r8) (:)]  total sensible heat flux (W/m**2) [+ to atm]
   eflx_sh_tot_u             =>    pef%eflx_sh_tot_u             , & ! Output: [real(r8) (:)]  urban total sensible heat flux (W/m**2) [+ to atm]
   eflx_sh_tot_r             =>    pef%eflx_sh_tot_r             , & ! Output: [real(r8) (:)]  rural total sensible heat flux (W/m**2) [+ to atm]
   eflx_lh_tot               =>    pef%eflx_lh_tot               , & ! Output: [real(r8) (:)]  total latent heat flux (W/m**2)  [+ to atm]
   eflx_lh_tot_u             =>    pef%eflx_lh_tot_u             , & ! Output: [real(r8) (:)]  urban total latent heat flux (W/m**2)  [+ to atm]
   eflx_lh_tot_r             =>    pef%eflx_lh_tot_r             , & ! Output: [real(r8) (:)]  rural total latent heat flux (W/m**2)  [+ to atm]
   eflx_lwrad_out            =>    pef%eflx_lwrad_out            , & ! Output: [real(r8) (:)]  emitted infrared (longwave) radiation (W/m**2)
   eflx_lwrad_net            =>    pef%eflx_lwrad_net            , & ! Output: [real(r8) (:)]  net infrared (longwave) rad (W/m**2) [+ = to atm]
   eflx_lwrad_net_u          =>    pef%eflx_lwrad_net_u          , & ! Output: [real(r8) (:)]  urban net infrared (longwave) rad (W/m**2) [+ = to atm]
   eflx_lwrad_net_r          =>    pef%eflx_lwrad_net_r          , & ! Output: [real(r8) (:)]  rural net infrared (longwave) rad (W/m**2) [+ = to atm]
   eflx_lwrad_out_u          =>    pef%eflx_lwrad_out_u          , & ! Output: [real(r8) (:)]  urban emitted infrared (longwave) rad (W/m**2)
   eflx_lwrad_out_r          =>    pef%eflx_lwrad_out_r          , & ! Output: [real(r8) (:)]  rural emitted infrared (longwave) rad (W/m**2)
   eflx_lh_vege              =>    pef%eflx_lh_vege              , & ! Output: [real(r8) (:)]  veg evaporation heat flux (W/m**2) [+ to atm]
   eflx_lh_vegt              =>    pef%eflx_lh_vegt              , & ! Output: [real(r8) (:)]  veg transpiration heat flux (W/m**2) [+ to atm]
   eflx_lh_grnd              =>    pef%eflx_lh_grnd              , & ! Output: [real(r8) (:)]  ground evaporation heat flux (W/m**2) [+ to atm]
   cgrnds                    =>    pef%cgrnds                    , & ! Input:  [real(r8) (:)]  deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
   cgrndl                    =>    pef%cgrndl                    , & ! Input:  [real(r8) (:)]  deriv of soil latent heat flux wrt soil temp [w/m**2/k]
   errsoi_pft                =>    pebal%errsoi                  , & ! Output: [real(r8) (:)]  pft-level soil/lake energy conservation error (W/m**2)
   wtcol                     =>   pft%wtcol                      , & ! Input:  [real(r8) (:)]  pft weight relative to column           
   eflx_wasteheat_pft        =>    pef%eflx_wasteheat_pft        , & ! Input:  [real(r8) (:)]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
   eflx_heat_from_ac_pft     =>    pef%eflx_heat_from_ac_pft     , & ! Input:  [real(r8) (:)]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
   eflx_traffic_pft          =>    pef%eflx_traffic_pft            & ! Input:  [real(r8) (:)]  traffic sensible heat flux (W/m**2)     
   )

    ! Get step size

    dtime = get_step_size()

    ! Determine soil temperatures including surface soil temperature

    call t_startf('soiltemperature')
    call SoilTemperature(bounds, num_urbanl, filter_urbanl, &
                         num_nolakec, filter_nolakec, &
                         xmf(bounds%begc:bounds%endc), &
                         fact(bounds%begc:bounds%endc, :), &
                         c_h2osfc(bounds%begc:bounds%endc), &
                         xmf_h2osfc(bounds%begc:bounds%endc))
    call t_stopf('soiltemperature')

    call t_startf('bgp2_loop_1')
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       j = snl(c)+1

       ! Calculate difference in soil temperature from last time step, for
       ! flux corrections

       if (snl(c) < 0) then
          t_grnd0(c) = frac_sno_eff(c) * tssbef(c,snl(c)+1) &
               + (1 - frac_sno_eff(c) - frac_h2osfc(c)) * tssbef(c,1) &
               + frac_h2osfc(c) * t_h2osfc_bef(c)
       else
          t_grnd0(c) = (1 - frac_h2osfc(c)) * tssbef(c,1) + frac_h2osfc(c) * t_h2osfc_bef(c)
       endif

       tinc(c) = t_grnd(c) - t_grnd0(c)

       ! Determine ratio of topsoil_evap_tot

       egsmax(c) = (h2osoi_ice(c,j)+h2osoi_liq(c,j)) / dtime

       ! added to trap very small negative soil water,ice

       if (egsmax(c) < 0._r8) then
          egsmax(c) = 0._r8
       end if
    end do

    ! A preliminary pft loop to determine if corrections are required for
    ! excess evaporation from the top soil layer... Includes new logic
    ! to distribute the corrections between pfts on the basis of their
    ! evaporative demands.
    ! egirat holds the ratio of demand to availability if demand is
    ! greater than availability, or 1.0 otherwise.
    ! Correct fluxes to present soil temperature

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       eflx_sh_grnd(p) = eflx_sh_grnd(p) + tinc(c)*cgrnds(p)
       qflx_evap_soi(p) = qflx_evap_soi(p) + tinc(c)*cgrndl(p)

       ! set ev_snow, ev_soil for urban landunits here
       l = plandunit(p)
       if (urbpoi(l)) then
          qflx_ev_snow(p) = qflx_evap_soi(p)
          qflx_ev_soil(p) = 0._r8
          qflx_ev_h2osfc(p) = 0._r8
       else
          qflx_ev_snow(p) = qflx_ev_snow(p) + tinc(c)*cgrndl(p)
          qflx_ev_soil(p) = qflx_ev_soil(p) + tinc(c)*cgrndl(p)
          qflx_ev_h2osfc(p) = qflx_ev_h2osfc(p) + tinc(c)*cgrndl(p)
       endif
    end do

    ! Set the column-average qflx_evap_soi as the weighted average over all pfts
    ! but only count the pfts that are evaporating

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       topsoil_evap_tot(c) = 0._r8
       sumwt(c) = 0._r8
    end do

    do pi = 1,max_pft_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             if (pactive(p)) then
                topsoil_evap_tot(c) = topsoil_evap_tot(c) + qflx_evap_soi(p) * wtcol(p)
             end if
          end if
       end do
    end do
    call t_stopf('bgp2_loop_1')
    call t_startf('bgp2_loop_2')

    ! Calculate ratio for rescaling pft-level fluxes to meet availability

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       if (topsoil_evap_tot(c) > egsmax(c)) then
          egirat(c) = (egsmax(c)/topsoil_evap_tot(c))
       else
          egirat(c) = 1.0_r8
       end if
    end do

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)
       j = snl(c)+1

       ! Correct soil fluxes for possible evaporation in excess of top layer water
       ! excess energy is added to the sensible heat flux from soil

       if (egirat(c) < 1.0_r8) then
          save_qflx_evap_soi = qflx_evap_soi(p)
          qflx_evap_soi(p) = qflx_evap_soi(p) * egirat(c)
          eflx_sh_grnd(p) = eflx_sh_grnd(p) + (save_qflx_evap_soi - qflx_evap_soi(p))*htvp(c)
          qflx_ev_snow(p) = qflx_ev_snow(p) * egirat(c)
          qflx_ev_soil(p) = qflx_ev_soil(p) * egirat(c)
          qflx_ev_h2osfc(p) = qflx_ev_h2osfc(p) * egirat(c)
       end if

       ! Ground heat flux

       if (.not. urbpoi(l)) then
          lw_grnd=(frac_sno_eff(c)*tssbef(c,snl(c)+1)**4 &
               +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
               +frac_h2osfc(c)*t_h2osfc_bef(c)**4)

          eflx_soil_grnd(p) = ((1._r8- frac_sno_eff(c))*sabg_soil(p) + frac_sno_eff(c)*sabg_snow(p)) + dlrad(p) &
               + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
            - emg(c)*sb*lw_grnd - emg(c)*sb*t_grnd0(c)**3*(4._r8*tinc(c)) &
               - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
          
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
            eflx_soil_grnd_r(p) = eflx_soil_grnd(p)
          end if
       else
          ! For all urban columns we use the net longwave radiation (eflx_lwrad_net) since
          ! the term (emg*sb*tssbef(snl+1)**4) is not the upward longwave flux because of 
          ! interactions between urban columns.

          eflx_lwrad_del(p) = 4._r8*emg(c)*sb*t_grnd0(c)**3*tinc(c)

          ! Include transpiration term because needed for pervious road
          ! and wasteheat and traffic flux
          eflx_soil_grnd(p) = sabg(p) + dlrad(p) &
                              - eflx_lwrad_net(p) - eflx_lwrad_del(p) &
                              - (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + qflx_tran_veg(p)*hvap) &
                              + eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
          eflx_soil_grnd_u(p) = eflx_soil_grnd(p)
       end if

       ! Total fluxes (vegetation + ground)

       eflx_sh_tot(p) = eflx_sh_veg(p) + eflx_sh_grnd(p)
       qflx_evap_tot(p) = qflx_evap_veg(p) + qflx_evap_soi(p)
       eflx_lh_tot(p)= hvap*qflx_evap_veg(p) + htvp(c)*qflx_evap_soi(p)
       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         eflx_lh_tot_r(p)= eflx_lh_tot(p)
         eflx_sh_tot_r(p)= eflx_sh_tot(p)
       else if (urbpoi(l)) then
         eflx_lh_tot_u(p)= eflx_lh_tot(p)
         eflx_sh_tot_u(p)= eflx_sh_tot(p)
       end if

       ! Assign ground evaporation to sublimation from soil ice or to dew
       ! on snow or ground

       qflx_evap_grnd(p) = 0._r8
       qflx_sub_snow(p) = 0._r8
       qflx_dew_snow(p) = 0._r8
       qflx_dew_grnd(p) = 0._r8

       if (qflx_ev_snow(p) >= 0._r8) then
          ! for evaporation partitioning between liquid evap and ice sublimation, 
	  ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
	  if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0.) then
             qflx_evap_grnd(p) = max(qflx_ev_snow(p)*(h2osoi_liq(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
	  else
	     qflx_evap_grnd(p) = 0.
	  end if
          qflx_sub_snow(p) = qflx_ev_snow(p) - qflx_evap_grnd(p)
       else
          if (t_grnd(c) < tfrz) then
             qflx_dew_snow(p) = abs(qflx_ev_snow(p))
          else
             qflx_dew_grnd(p) = abs(qflx_ev_snow(p))
          end if
       end if

       ! Update the pft-level qflx_snwcp
       ! This was moved in from Hydrology2 to keep all pft-level
       ! calculations out of Hydrology2

       if (snl(c) < 0 .and. do_capsnow(c)) then
          qflx_snwcp_liq(p) = qflx_snwcp_liq(p)+frac_sno_eff(c)*qflx_dew_grnd(p)
          qflx_snwcp_ice(p) = qflx_snwcp_ice(p)+frac_sno_eff(c)*qflx_dew_snow(p)
       end if

       ! Variables needed by history tape

       qflx_evap_can(p)  = qflx_evap_veg(p) - qflx_tran_veg(p)
       eflx_lh_vege(p)   = (qflx_evap_veg(p) - qflx_tran_veg(p)) * hvap
       eflx_lh_vegt(p)   = qflx_tran_veg(p) * hvap
       eflx_lh_grnd(p)   = qflx_evap_soi(p) * htvp(c)

    end do
    call t_stopf('bgp2_loop_2')
    call t_startf('bgp2_loop_3')

    ! Soil Energy balance check

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       errsoi_pft(p) = eflx_soil_grnd(p) - xmf(c) - xmf_h2osfc(c) &
            - frac_h2osfc(c)*(t_h2osfc(c)-t_h2osfc_bef(c)) &
            *(c_h2osfc(c)/dtime)

       ! For urban sunwall, shadewall, and roof columns, the "soil" energy balance check
       ! must include the heat flux from the interior of the building.
       if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall .or. ctype(c)==icol_roof) then
          errsoi_pft(p) = errsoi_pft(p) + eflx_building_heat(c) 
       end if
    end do
    do j = -nlevsno+1,nlevgrnd
       do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          c = pcolumn(p)

          if ((ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
               .and. ctype(c) /= icol_roof) .or. ( j <= nlevurb)) then
            ! area weight heat absorbed by snow layers
             if (j >= snl(c)+1 .and. j < 1) errsoi_pft(p) = errsoi_pft(p) &
                  - frac_sno_eff(c)*(t_soisno(c,j)-tssbef(c,j))/fact(c,j)
             if (j >= 1) errsoi_pft(p) = errsoi_pft(p) &
                  - (t_soisno(c,j)-tssbef(c,j))/fact(c,j)
          end if
       end do
    end do
    call t_stopf('bgp2_loop_3')
    call t_startf('bgp2_loop_4')

    ! Outgoing long-wave radiation from vegetation + ground
    ! For conservation we put the increase of ground longwave to outgoing
    ! For urban pfts, ulrad=0 and (1-fracveg_nosno)=1, and eflx_lwrad_out and eflx_lwrad_net 
    ! are calculated in UrbanRadiation. The increase of ground longwave is added directly 
    ! to the outgoing longwave and the net longwave.

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)
       j = snl(c)+1

       if (.not. urbpoi(l)) then
          lw_grnd=(frac_sno_eff(c)*tssbef(c,snl(c)+1)**4 &
               +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
               +frac_h2osfc(c)*t_h2osfc_bef(c)**4)

          eflx_lwrad_out(p) = ulrad(p) &
               + (1-frac_veg_nosno(p))*(1.-emg(c))*forc_lwrad(c) &
               + (1-frac_veg_nosno(p))*emg(c)*sb*lw_grnd &
               + 4._r8*emg(c)*sb*t_grnd0(c)**3*tinc(c)

          eflx_lwrad_net(p) = eflx_lwrad_out(p) - forc_lwrad(c)
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
            eflx_lwrad_net_r(p) = eflx_lwrad_out(p) - forc_lwrad(c)
            eflx_lwrad_out_r(p) = eflx_lwrad_out(p)
          end if
       else
          eflx_lwrad_out(p) = eflx_lwrad_out(p) + eflx_lwrad_del(p)
          eflx_lwrad_net(p) = eflx_lwrad_net(p) + eflx_lwrad_del(p)
          eflx_lwrad_net_u(p) = eflx_lwrad_net_u(p) + eflx_lwrad_del(p)
          eflx_lwrad_out_u(p) = eflx_lwrad_out(p)
       end if
    end do

    ! lake balance for errsoi is not over pft
    ! therefore obtain column-level radiative temperature

    call p2c(bounds, num_nolakec, filter_nolakec, &
         errsoi_pft(bounds%begp:bounds%endp), &
         errsoi_col(bounds%begc:bounds%endc))
    call t_stopf('bgp2_loop_4')

    end associate 
   end subroutine Biogeophysics2

end module Biogeophysics2Mod
