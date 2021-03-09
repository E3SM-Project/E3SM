module BareGroundFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute sensible and latent fluxes and their derivatives with respect
  ! to ground temperature using ground temperatures from previous time step.
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use decompMod            , only : bounds_type
  use CH4Mod               , only : ch4_type
  use atm2lndType          , only : atm2lnd_type
  use CanopyStateType      , only : canopystate_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use SoilStateType        , only : soilstate_type
  use TemperatureType      , only : temperature_type
  use WaterfluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use TopounitDataType     , only : top_as
  use LandunitType         , only : lun_pp
  use ColumnType           , only : col_pp
  use ColumnDataType       , only : col_es, col_ef, col_ws
  use VegetationType       , only : veg_pp
  use VegetationDataType   , only : veg_es, veg_ef, veg_ws, veg_wf
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine BareGroundFluxes(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
       atm2lnd_vars, canopystate_vars, soilstate_vars, &
       frictionvel_vars, ch4_vars, energyflux_vars, temperature_vars, &
       waterflux_vars, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step.
    !
    ! !USES:
    use shr_const_mod        , only : SHR_CONST_RGAS
    use elm_varpar           , only : nlevgrnd
    use elm_varcon           , only : cpair, denice, denh2o
    use elm_varctl           , only : use_lch4
    use landunit_varcon      , only : istsoil, istcrop
    use QSatMod              , only : QSat
    use SurfaceResistanceMod , only : do_soilevap_beta
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_nolakeurbanp          ! number of pft non-lake, non-urban points in pft filter
    integer                , intent(in)    :: filter_nolakeurbanp(:)    ! patch filter for non-lake, non-urban points
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    type(ch4_type)         , intent(inout) :: ch4_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(temperature_type) , intent(inout) :: temperature_vars
    type(waterflux_type)   , intent(inout) :: waterflux_vars
    type(waterstate_type)  , intent(inout) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: dt_frac_expl = 0.5_r8 ! fraction of the time step taken explicitly
    real(r8), parameter :: dt_frac_impl = 1._r8 - dt_frac_expl ! fraction taken implicitly
    integer  :: p,c,t,g,f,j,l                    ! indices
    integer  :: filterp(bounds%endp-bounds%begp+1) ! patch filter for vegetated patches
    integer  :: fn                               ! number of values in local pft filter
    integer  :: fp                               ! lake filter pft index
    real(r8) :: dth(bounds%begp:bounds%endp)     ! diff of virtual temp. between ref. height and surface
    real(r8) :: dqh(bounds%begp:bounds%endp)     ! diff of humidity between ref. height and surface
    real(r8) :: um(bounds%begp:bounds%endp)      ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(bounds%begp:bounds%endp)   ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begp:bounds%endp) ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begp:bounds%endp)   ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begp:bounds%endp) ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(bounds%begp:bounds%endp)   ! friction velocity [m/s]
    real(r8) :: cf                               ! heat transfer coefficient from leaves [-]
    real(r8) :: ram                              ! aerodynamical resistance [s/m]
    real(r8) :: rah                              ! thermal resistance [s/m]
    real(r8) :: raw                              ! moisture resistance [s/m]
    real(r8) :: raih                             ! temporary variable [kg/m2/s]
    real(r8) :: raiw                             ! temporary variable [kg/m2/s]
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    real(r8) :: www                    ! surface soil wetness [-]
    real(r8) :: thm_adj(bounds%begp:bounds%endp) ! Adjusted thm for iteration
    real(r8) :: forc_q_adj(bounds%begp:bounds%endp) ! Adjusted forc_q for iteration
    real(r8) :: forc_u_adj(bounds%begp:bounds%endp) ! Adjusted forc_u for iteration
    real(r8) :: forc_v_adj(bounds%begp:bounds%endp) ! Adjusted forc_v for iteration
    real(r8) :: t_end(bounds%begp:bounds%endp) ! temperature output from iteration
    real(r8) :: q_end(bounds%begp:bounds%endp) ! humidity output from iteration
    real(r8) :: u_end(bounds%begp:bounds%endp) ! zonal wind output from iteration
    real(r8) :: v_end(bounds%begp:bounds%endp) ! meridional wind output from iteration
    real(r8) :: eflx_expl(bounds%begp:bounds%endp) ! sensible heat flux from explicit step
    real(r8) :: qflx_expl(bounds%begp:bounds%endp) ! vapor flux from explicit step
    real(r8) :: taux_expl(bounds%begp:bounds%endp) ! taux from explicit step
    real(r8) :: tauy_expl(bounds%begp:bounds%endp) ! tauy from explicit step
    !------------------------------------------------------------------------------

    associate(                                                          &
         snl              =>    col_pp%snl                            , & ! Input:  [integer  (:)   ]  number of snow layers
         dz               =>    col_pp%dz                             , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         forc_u           =>    top_as%ubot                           , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)
         forc_v           =>    top_as%vbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)
         u_diff           =>    top_as%u_diff                         , & ! Input:  [real(r8) (:)   ]  approximate atmosphere change to zonal wind (m/s)
         v_diff           =>    top_as%v_diff                         , & ! Input:  [real(r8) (:)   ]  approximate atmosphere change to meridional wind (m/s)
         forc_pbot        =>    top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_rho         =>    top_as%rhobot                         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
         forc_q           =>    top_as%qbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)

         frac_veg_nosno   =>    canopystate_vars%frac_veg_nosno_patch , & ! Input:  [logical  (:)   ]  true=> pft is bare ground (elai+esai = zero)

         htvp             =>    col_ef%htvp             , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg]

         watsat           =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         soilbeta         =>    soilstate_vars%soilbeta_col           , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity

         t_soisno         =>    col_es%t_soisno         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_grnd           =>    col_es%t_grnd           , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]
         thm              =>    veg_es%thm              , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
         t_h2osfc         =>    col_es%t_h2osfc         , & ! Input:  [real(r8) (:)   ]  surface water temperature

         frac_sno         =>    col_ws%frac_sno          , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         qg_snow          =>    col_ws%qg_snow           , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]
         qg_soil          =>    col_ws%qg_soil           , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]
         qg_h2osfc        =>    col_ws%qg_h2osfc         , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]
         qg               =>    col_ws%qg                , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]
         dqgdT            =>    col_ws%dqgdT             , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"
         h2osoi_ice       =>    col_ws%h2osoi_ice        , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq       =>    col_ws%h2osoi_liq        , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)

         grnd_ch4_cond    =>    ch4_vars%grnd_ch4_cond_patch          , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]

         eflx_sh_snow     =>    veg_ef%eflx_sh_snow    , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_soil     =>    veg_ef%eflx_sh_soil    , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_h2osfc   =>    veg_ef%eflx_sh_h2osfc  , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_grnd     =>    veg_ef%eflx_sh_grnd    , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_tot      =>    veg_ef%eflx_sh_tot     , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         taux             =>    veg_ef%taux            , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)
         tauy             =>    veg_ef%tauy            , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)
         dlrad            =>    veg_ef%dlrad           , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]
         ulrad            =>    veg_ef%ulrad           , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]
         cgrnds           =>    veg_ef%cgrnds          , & ! Output: [real(r8) (:)   ]  deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
         cgrndl           =>    veg_ef%cgrndl          , & ! Output: [real(r8) (:)   ]  deriv of soil latent heat flux wrt soil temp [w/m**2/k]
         cgrnd            =>    veg_ef%cgrnd           , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]

         t_ref2m          =>    veg_es%t_ref2m        , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)
         t_ref2m_r        =>    veg_es%t_ref2m_r      , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin)

         q_ref2m          =>    veg_ws%q_ref2m         , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)
         rh_ref2m_r       =>    veg_ws%rh_ref2m_r      , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)
         rh_ref2m         =>    veg_ws%rh_ref2m        , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)

         ram1             =>    frictionvel_vars%ram1_patch           , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)

         qflx_ev_snow     =>    veg_wf%qflx_ev_snow     , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]
         qflx_ev_soil     =>    veg_wf%qflx_ev_soil     , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]
         qflx_ev_h2osfc   =>    veg_wf%qflx_ev_h2osfc   , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]
         qflx_evap_soi    =>    veg_wf%qflx_evap_soi    , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_evap_tot    =>    veg_wf%qflx_evap_tot    , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
         begp             =>    bounds%begp             , &
         endp             =>    bounds%endp               &
         )

      !---------------------------------------------------
      ! Filter patches where frac_veg_nosno IS ZERO
      !---------------------------------------------------

      fn = 0
      do fp = 1,num_nolakeurbanp
         p = filter_nolakeurbanp(fp)
         if (frac_veg_nosno(p) == 0) then
            fn = fn + 1
            filterp(fn) = p
         end if
      end do

      ! Compute sensible and latent fluxes and their derivatives with respect
      ! to ground temperature using ground temperatures from previous time step

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! Initialization variables

         dlrad(p)  = 0._r8
         ulrad(p)  = 0._r8

         thm_adj(p) = thm(p)
         forc_q_adj(p) = forc_q(t)
         forc_u_adj(p) = forc_u(t) + u_diff(t)
         forc_v_adj(p) = forc_v(t) + v_diff(t)

      end do

      ! Perform stability iteration
      ! Determine friction velocity, and potential temperature and humidity
      ! profiles of the surface boundary layer

      ! We perform the integration with two half-steps, one explicit and one
      ! implicit, i.e. the trapezoidal method.

      ! Do explicit integration for the first half step, i.e. do_implicit=.false.
      call BareGroundStability(.false., dt_frac_expl, bounds, soilstate_vars, &
           fn, filterp, &
           thm_adj, forc_q_adj, forc_u_adj, forc_v_adj, t_end, q_end, &
           u_end, v_end, eflx_expl, qflx_expl, taux_expl, tauy_expl, &
           um, ustar, temp1, temp12m, temp2, temp22m, frictionvel_vars)

      ! Update atmospheric state.
      do f = 1, fn
         p = filterp(f)
         thm_adj(p) = t_end(p)
         forc_q_adj(p) = q_end(p)
         forc_u_adj(p) = u_end(p)
         forc_v_adj(p) = v_end(p)
      end do

      ! Do implicit integration.
      call BareGroundStability(.true., dt_frac_impl, bounds, soilstate_vars, &
           fn, filterp, &
           thm_adj, forc_q_adj, forc_u_adj, forc_v_adj, t_end, q_end, &
           u_end, v_end, eflx_sh_grnd(begp:endp), qflx_evap_soi(begp:endp), &
           taux(begp:endp), tauy(begp:endp), &
           um, ustar, temp1, temp12m, temp2, temp22m, frictionvel_vars)

      ! Update atmospheric state.
      do f = 1, fn
         p = filterp(f)
         thm_adj(p) = t_end(p)
         forc_q_adj(p) = q_end(p)
         forc_u_adj(p) = u_end(p)
         forc_v_adj(p) = v_end(p)
      end do

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)
         t = veg_pp%topounit(p)
         l = veg_pp%landunit(p)

         ! Determine aerodynamic resistances

         ram  = 1._r8/(ustar(p)*ustar(p)/um(p))
         rah  = 1._r8/(temp1(p)*ustar(p))
         raw  = 1._r8/(temp2(p)*ustar(p))
         raih = forc_rho(t)*cpair/rah
         if (use_lch4) then
            grnd_ch4_cond(p) = 1._r8/raw
         end if

         ! Soil evaporation resistance
         www = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
         www = min(max(www,0.0_r8),1._r8)

         dqh(p) = q_end(p) - qg(c)
         dth(p)   = t_end(p)-t_grnd(c)

         !changed by K.Sakaguchi. Soilbeta is used for evaporation
         if (dqh(p) > 0._r8) then  !dew  (beta is not applied, just like rsoil used to be)
            raiw = forc_rho(t)/(raw)
         else
            if(do_soilevap_beta())then
               ! Lee and Pielke 1992 beta is applied
               raiw    = soilbeta(c)*forc_rho(t)/(raw)
            endif
         end if

         ram1(p) = ram  !pass value to global variable

         ! Output to pft-level data structures
         ! Derivative of fluxes with respect to ground temperature
         cgrnds(p) = raih
         cgrndl(p) = raiw*dqgdT(c)
         cgrnd(p)  = cgrnds(p) + htvp(c)*cgrndl(p)

         ! Surface fluxes of momentum, sensible and latent heat
         ! using ground temperatures from previous time step
         taux(p)          = dt_frac_expl * taux_expl(p) + dt_frac_impl * taux(p)
         tauy(p)          = dt_frac_expl * taux_expl(p) + dt_frac_impl * taux(p)
         eflx_sh_grnd(p)  = dt_frac_expl * eflx_expl(p) + dt_frac_impl * eflx_sh_grnd(p)
         eflx_sh_tot(p)   = eflx_sh_grnd(p)

         ! compute sensible heat fluxes individually
         eflx_sh_snow(p)   = eflx_sh_grnd(p)*(thm_adj(p)-t_soisno(c,snl(c)+1))/dth(p)
         eflx_sh_soil(p)   = eflx_sh_grnd(p)*(thm_adj(p)-t_soisno(c,1))/dth(p)
         eflx_sh_h2osfc(p) = eflx_sh_grnd(p)*(thm_adj(p)-t_h2osfc(c))/dth(p)

         ! water fluxes from soil
         qflx_evap_soi(p)  = dt_frac_expl * qflx_expl(p) + dt_frac_impl * qflx_evap_soi(p)
         qflx_evap_tot(p)  = qflx_evap_soi(p)

         ! compute latent heat fluxes individually
         qflx_ev_snow(p)   = qflx_evap_soi(p)*(forc_q_adj(p) - qg_snow(c))/dqh(p)
         qflx_ev_soil(p)   = qflx_evap_soi(p)*(forc_q_adj(p) - qg_soil(c))/dqh(p)
         qflx_ev_h2osfc(p) = qflx_evap_soi(p)*(forc_q_adj(p) - qg_h2osfc(c))/dqh(p)

         ! 2 m height air temperature
         t_ref2m(p) = thm_adj(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

         ! 2 m height specific humidity
         q_ref2m(p) = forc_q_adj(p) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

         ! 2 m height relative humidity
         call QSat(t_ref2m(p), forc_pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)

         if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
            rh_ref2m_r(p) = rh_ref2m(p)
            t_ref2m_r(p) = t_ref2m(p)
         end if

      end do

    end associate

  end subroutine BareGroundFluxes

  subroutine BareGroundStability(do_implicit, dt_frac, bounds, soilstate_vars, &
       fn, filterp, &
       t_start, q_start, u_start, v_start, t_end, q_end, &
       u_end, v_end, eflx_sh_grnd, qflx_evap_soi, taux, tauy, &
       um, ustar, temp1, temp12m, temp2, temp22m, frictionvel_vars)
    !
    ! !DESCRIPTION:
    ! Stability iteration used to determine friction velocity and surface
    ! boundary layer profiles for a given atmospheric state.
    !
    ! !USES:
    use elm_varcon           , only : cpair, vkc, grav
    use FrictionVelocityMod  , only : FrictionVelocity, MoninObukIni
    use SurfaceResistanceMod , only : do_soilevap_beta
    !
    ! !ARGUMENTS:
    logical                , intent(in)    :: do_implicit                      ! Use backward rather than forward Euler method.
    real(r8)               , intent(in)    :: dt_frac                          ! Take this fraction of a full time step.
    type(bounds_type)      , intent(in)    :: bounds
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    integer                , intent(in)    :: fn                               ! number of values in local pft filter
    integer                , intent(in)    :: filterp(bounds%endp-bounds%begp+1) ! patch filter for vegetated patches
    real(r8)               , intent(in)    :: t_start(bounds%begp:bounds%endp) ! start-of-iteration temperature [K]
    real(r8)               , intent(in)    :: q_start(bounds%begp:bounds%endp) ! start-of-iteration specific humidity [kg/kg]
    real(r8)               , intent(in)    :: u_start(bounds%begp:bounds%endp) ! start-of-iteration zonal wind [m/s]
    real(r8)               , intent(in)    :: v_start(bounds%begp:bounds%endp) ! start-of-iteration meridional wind [m/s]
    real(r8)               , intent(out)   :: t_end(bounds%begp:bounds%endp)   ! end-of-iteration temperature [K]
    real(r8)               , intent(out)   :: q_end(bounds%begp:bounds%endp)   ! end-of-iteration specific humidity [kg/kg]
    real(r8)               , intent(out)   :: u_end(bounds%begp:bounds%endp)   ! end-of-iteration zonal wind [m/s]
    real(r8)               , intent(out)   :: v_end(bounds%begp:bounds%endp)   ! end-of-iteration meridional wind [m/s]
    real(r8)               , intent(out)   :: eflx_sh_grnd(bounds%begp:bounds%endp)  ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8)               , intent(out)   :: qflx_evap_soi(bounds%begp:bounds%endp) ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8)               , intent(out)   :: taux(bounds%begp:bounds%endp)    ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8)               , intent(out)   :: tauy(bounds%begp:bounds%endp)    ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8)               , intent(out)   :: um(bounds%begp:bounds%endp)      ! wind speed including the stablity effect [m/s]
    real(r8)               , intent(out)   :: ustar(bounds%begp:bounds%endp)   ! friction velocity [m/s]
    real(r8)               , intent(out)   :: temp1(bounds%begp:bounds%endp)   ! relation for potential temperature profile
    real(r8)               , intent(out)   :: temp12m(bounds%begp:bounds%endp) ! relation for potential temperature profile applied at 2-m
    real(r8)               , intent(out)   :: temp2(bounds%begp:bounds%endp)   ! relation for specific humidity profile
    real(r8)               , intent(out)   :: temp22m(bounds%begp:bounds%endp) ! relation for specific humidity profile applied at 2-m
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    !
    ! !LOCAL VARIABLES:
    integer, parameter  :: niters = 30           ! maximum number of iterations for surface temperature
    integer :: iter                              ! iteration index
    integer :: p,c,t,f                         ! indices
    real(r8) :: displa(bounds%begp:bounds%endp)  ! displacement height [m]
    real(r8) :: zldis(bounds%begp:bounds%endp)   ! reference height "minus" zero displacement height [m]
    real(r8) :: z0mg_patch(bounds%begp:bounds%endp)
    real(r8) :: z0hg_patch(bounds%begp:bounds%endp)
    real(r8) :: z0qg_patch(bounds%begp:bounds%endp)
    real(r8) :: ur(bounds%begp:bounds%endp)      ! wind speed at reference height [m/s]
    real(r8) :: dth(bounds%begp:bounds%endp)     ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                             ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(bounds%begp:bounds%endp)     ! diff of humidity between ref. height and surface
    real(r8) :: obu(bounds%begp:bounds%endp)     ! Monin-Obukhov length (m)
    real(r8) :: fm(bounds%begp:bounds%endp)      ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: ram                              ! aerodynamical resistance [s/m]
    real(r8) :: rah                              ! thermal resistance [s/m]
    real(r8) :: raw                              ! moisture resistance [s/m]
    real(r8) :: raih                             ! temporary variable [kg/m2/s]
    real(r8) :: raiw                             ! temporary variable [kg/m2/s]
    real(r8) :: tstar                            ! temperature scaling parameter
    real(r8) :: qstar                            ! moisture scaling parameter
    real(r8) :: thvstar                          ! virtual potential temperature scaling parameter
    real(r8) :: beta                             ! coefficient of convective velocity [-]
    real(r8) :: zeta                             ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                               ! convective velocity [m/s]

    associate(                                                      &
         zii              =>    col_pp%zii                        , & ! Input:  [real(r8) (:)   ]  convective boundary height [m]
         forc_th          =>    top_as%thbot                      , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)
         forc_rho         =>    top_as%rhobot                     , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
         tresp            =>    top_as%tresp                      , & ! Input:  [real(r8) (:)   ]  response of temperature to surface flux (K m^2/W)
         qresp            =>    top_as%qresp                      , & ! Input:  [real(r8) (:)   ]  response of spec. humidity to surface flux (m^2 s/kg)
         wsresp           =>    top_as%wsresp                     , & ! Input:  [real(r8) (:)   ]  response of wind to surface stress (m/s/Pa)
         soilbeta         =>    soilstate_vars%soilbeta_col       , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity
         t_grnd           =>    col_es%t_grnd                     , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]
         thv              =>    col_es%thv                        , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)
         qg               =>    col_ws%qg                         , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]
         forc_hgt_u_patch =>    frictionvel_vars%forc_hgt_u_patch , & ! Input:  [real(r8) (:)   ]
         z0mg_col         =>    frictionvel_vars%z0mg_col         , & ! Input:  [real(r8) (:)   ]  roughness length, momentum [m]
         z0hg_col         =>    frictionvel_vars%z0hg_col         , & ! Input:  [real(r8) (:)   ]  roughness length, sensible heat [m]
         z0qg_col         =>    frictionvel_vars%z0qg_col         , & ! Input:  [real(r8) (:)   ]  roughness length, latent heat [m]
         begp             =>    bounds%begp                       , &
         endp             =>    bounds%endp                         &
         )

      beta = 1._r8 ! previously set as a constant for all columns in CanopyTemperature()

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)

         ! Initialize variables
         displa(p) = 0._r8
         zldis(p) = forc_hgt_u_patch(p)

         t_end(p) = t_start(p)
         q_end(p) = q_start(p)
         u_end(p) = u_start(p)
         v_end(p) = v_start(p)

         dth(p)   = t_end(p)-t_grnd(c)
         dqh(p)   = q_end(p) - qg(c)
         dthv     = dth(p)*(1._r8+0.61_r8*q_end(p))+0.61_r8*forc_th(t)*dqh(p)

         ur(p)    = max(1.0_r8,sqrt(u_end(p)*u_end(p)+v_end(p)*v_end(p)))

         ! Copy column roughness to local pft-level arrays

         z0mg_patch(p) = z0mg_col(c)
         z0hg_patch(p) = z0hg_col(c)
         z0qg_patch(p) = z0qg_col(c)

         ! Initialize Monin-Obukhov length and wind speed

         call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_patch(p), um(p), obu(p))
      end do

      do iter = 1, niters

         call FrictionVelocity(begp, endp, fn, filterp, &
              displa(begp:endp), z0mg_patch(begp:endp), z0hg_patch(begp:endp), z0qg_patch(begp:endp), &
              obu(begp:endp), iter, ur(begp:endp), um(begp:endp), ustar(begp:endp), &
              temp1(begp:endp), temp2(begp:endp), temp12m(begp:endp), temp22m(begp:endp), fm(begp:endp), &
              frictionvel_vars)

         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)

            if (do_implicit) then
               ! Update winds.
               ram = 1._r8/(ustar(p)*ustar(p)/um(p))
               ! Forbid removing more than 99% of wind speed in a time step.
               ! This is mainly to avoid convergence issues since this is such a
               ! basic form of iteration in this loop...
               if ( forc_rho(t)*dt_frac*wsresp(t) / ram > 0.99 ) then
                  taux(p) = -0.99 * u_end(p) / (dt_frac * wsresp(t))
                  tauy(p) = -0.99 * v_end(p) / (dt_frac * wsresp(t))
               else
                  taux(p) = -forc_rho(t)*u_end(p)/ram
                  tauy(p) = -forc_rho(t)*v_end(p)/ram
               end if
               u_end(p) = u_start(p) + taux(p) * dt_frac * wsresp(t)
               v_end(p) = v_start(p) + tauy(p) * dt_frac * wsresp(t)
               ur(p) = max(1.0_r8,sqrt(u_end(p)*u_end(p)+v_end(p)*v_end(p)))

               ! Update air temperature.
               rah  = 1._r8/(temp1(p)*ustar(p))
               raih = forc_rho(t)*cpair/rah
               ! Limit as with wind speed.
               if ( raih*dt_frac*tresp(t) > 0.99 ) then
                  eflx_sh_grnd(p) = -0.99 * dth(p) / (dt_frac * tresp(t))
               else
                  eflx_sh_grnd(p)  = -raih*dth(p)
               end if
               t_end(p) = t_start(p) + eflx_sh_grnd(p) * dt_frac * tresp(t)
               dth(p) = t_end(p)-t_grnd(c)

               ! Update humidity.
               raw  = 1._r8/(temp2(p)*ustar(p))
               !changed by K.Sakaguchi. Soilbeta is used for evaporation
               if (dqh(p) > 0._r8) then  !dew  (beta is not applied, just like rsoil used to be)
                  raiw = forc_rho(t)/(raw)
               else
                  if(do_soilevap_beta())then
                     ! Lee and Pielke 1992 beta is applied
                     raiw    = soilbeta(c)*forc_rho(t)/(raw)
                  endif
               end if
               ! Limit as with wind speed.
               if ( raiw*dt_frac*qresp(t) > 0.99 ) then
                  qflx_evap_soi(p) = -0.99 * dqh(p) / (dt_frac * qresp(t))
               else
                  qflx_evap_soi(p) = -raiw*dqh(p)
               end if
               q_end(p) = q_start(p) + qflx_evap_soi(p) * dt_frac * qresp(t)
               dqh(p) = q_end(p) - qg(c)
            end if

            tstar = temp1(p)*dth(p)
            qstar = temp2(p)*dqh(p)
            z0hg_patch(p) = z0mg_patch(p)/exp(0.13_r8 * (ustar(p)*z0mg_patch(p)/1.5e-5_r8)**0.45_r8)
            z0qg_patch(p) = z0hg_patch(p)
            ! Adjustment to potential temperature assumes exner function is ~1,
            ! which should be true for a relatively fine atmospheric grid...
            thvstar = tstar*(1._r8+0.61_r8*q_end(p)) + 0.61_r8*(forc_th(t) + t_end(p) - t_start(p))*qstar
            zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(p) = max(ur(p),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               wc = beta*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_r8
               um(p) = sqrt(ur(p)*ur(p) + wc*wc)
            end if
            obu(p) = zldis(p)/zeta
         end do
      end do

      ! Explicit update of _end outputs.
      if (.not. do_implicit) then
         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)

            ! Update winds.
            ram = 1._r8/(ustar(p)*ustar(p)/um(p))
            ! Forbid removing more than 99% of wind speed in a time step.
            ! This is mainly to avoid convergence issues since this is such a
            ! basic form of iteration in this loop...
            if ( forc_rho(t)*dt_frac*wsresp(t) / ram > 0.99 ) then
               taux(p) = -0.99 * u_end(p) / (dt_frac * wsresp(t))
               tauy(p) = -0.99 * v_end(p) / (dt_frac * wsresp(t))
            else
               taux(p) = -forc_rho(t)*u_end(p)/ram
               tauy(p) = -forc_rho(t)*v_end(p)/ram
            end if
            u_end(p) = u_start(p) + taux(p) * dt_frac * wsresp(t)
            v_end(p) = v_start(p) + tauy(p) * dt_frac * wsresp(t)

            ! Update air temperature.
            rah  = 1._r8/(temp1(p)*ustar(p))
            raih = forc_rho(t)*cpair/rah
            ! Limit as with wind speed.
            if ( raih*dt_frac*tresp(t) > 0.99 ) then
               eflx_sh_grnd(p) = -0.99 * dth(p) / (dt_frac * tresp(t))
            else
               eflx_sh_grnd(p)  = -raih*dth(p)
            end if
            t_end(p) = t_start(p) + eflx_sh_grnd(p) * dt_frac * tresp(t)

            ! Update humidity.
            raw  = 1._r8/(temp2(p)*ustar(p))
            !changed by K.Sakaguchi. Soilbeta is used for evaporation
            if (dqh(p) > 0._r8) then  !dew  (beta is not applied, just like rsoil used to be)
               raiw = forc_rho(t)/(raw)
            else
               if(do_soilevap_beta())then
                  ! Lee and Pielke 1992 beta is applied
                  raiw    = soilbeta(c)*forc_rho(t)/(raw)
               endif
            end if
            ! Limit as with wind speed.
            if ( raiw*dt_frac*qresp(t) > 0.99 ) then
               qflx_evap_soi(p) = -0.99 * dqh(p) / (dt_frac * qresp(t))
            else
               qflx_evap_soi(p) = -raiw*dqh(p)
            end if
            q_end(p) = q_start(p) + qflx_evap_soi(p) * dt_frac * qresp(t)
         end do
      end if
    end associate

  end subroutine BareGroundStability

end module BareGroundFluxesMod
