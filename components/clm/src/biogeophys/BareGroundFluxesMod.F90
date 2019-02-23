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
  use TopounitType         , only : top_as
  use LandunitType         , only : lun_pp
  use ColumnType           , only : col_pp
  use VegetationType       , only : veg_pp
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
    use clm_varpar           , only : nlevgrnd
    use clm_varcon           , only : cpair, vkc, grav, denice, denh2o
    use clm_varctl           , only : use_lch4
    use landunit_varcon      , only : istsoil, istcrop
    use FrictionVelocityMod  , only : FrictionVelocity, MoninObukIni
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
    integer, parameter  :: niters = 3            ! maximum number of iterations for surface temperature
    integer  :: p,c,t,g,f,j,l                    ! indices
    integer  :: filterp(bounds%endp-bounds%begp+1) ! patch filter for vegetated patches
    integer  :: fn                               ! number of values in local pft filter
    integer  :: fp                               ! lake filter pft index
    integer  :: iter                             ! iteration index
    real(r8) :: zldis(bounds%begp:bounds%endp)   ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(bounds%begp:bounds%endp)  ! displacement height [m]
    real(r8) :: zeta                             ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                               ! convective velocity [m/s]
    real(r8) :: dth(bounds%begp:bounds%endp)     ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                             ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(bounds%begp:bounds%endp)     ! diff of humidity between ref. height and surface
    real(r8) :: obu(bounds%begp:bounds%endp)     ! Monin-Obukhov length (m)
    real(r8) :: ur(bounds%begp:bounds%endp)      ! wind speed at reference height [m/s]
    real(r8) :: um(bounds%begp:bounds%endp)      ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(bounds%begp:bounds%endp)   ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begp:bounds%endp) ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begp:bounds%endp)   ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begp:bounds%endp) ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(bounds%begp:bounds%endp)   ! friction velocity [m/s]
    real(r8) :: tstar                            ! temperature scaling parameter
    real(r8) :: qstar                            ! moisture scaling parameter
    real(r8) :: thvstar                          ! virtual potential temperature scaling parameter
    real(r8) :: cf                               ! heat transfer coefficient from leaves [-]
    real(r8) :: ram                              ! aerodynamical resistance [s/m]
    real(r8) :: rah                              ! thermal resistance [s/m]
    real(r8) :: raw                              ! moisture resistance [s/m]
    real(r8) :: raih                             ! temporary variable [kg/m2/s]
    real(r8) :: raiw                             ! temporary variable [kg/m2/s]
    real(r8) :: fm(bounds%begp:bounds%endp)      ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: z0mg_patch(bounds%begp:bounds%endp)
    real(r8) :: z0hg_patch(bounds%begp:bounds%endp)
    real(r8) :: z0qg_patch(bounds%begp:bounds%endp)
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    real(r8) :: www                    ! surface soil wetness [-]
    !------------------------------------------------------------------------------

    associate(                                                          &
         snl              =>    col_pp%snl                            , & ! Input:  [integer  (:)   ]  number of snow layers
         dz               =>    col_pp%dz                             , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         zii              =>    col_pp%zii                            , & ! Input:  [real(r8) (:)   ]  convective boundary height [m]
         forc_u           =>    top_as%ubot                           , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)
         forc_v           =>    top_as%vbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)
         forc_th          =>    top_as%thbot                          , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)
         forc_pbot        =>    top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_rho         =>    top_as%rhobot                         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
         forc_q           =>    top_as%qbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)

         forc_hgt_u_patch =>    frictionvel_vars%forc_hgt_u_patch     , & ! Input:

         frac_veg_nosno   =>    canopystate_vars%frac_veg_nosno_patch , & ! Input:  [logical  (:)   ]  true=> pft is bare ground (elai+esai = zero)

         htvp             =>    energyflux_vars%htvp_col              , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg]

         watsat           =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         soilbeta         =>    soilstate_vars%soilbeta_col           , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity

         t_soisno         =>    temperature_vars%t_soisno_col         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_grnd           =>    temperature_vars%t_grnd_col           , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]
         thv              =>    temperature_vars%thv_col              , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)
         thm              =>    temperature_vars%thm_patch            , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
         t_h2osfc         =>    temperature_vars%t_h2osfc_col         , & ! Input:  [real(r8) (:)   ]  surface water temperature
         beta             =>    temperature_vars%beta_col             , & ! Input:  [real(r8) (:)   ]  coefficient of conective velocity [-]

         frac_sno         =>    waterstate_vars%frac_sno_col          , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         qg_snow          =>    waterstate_vars%qg_snow_col           , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]
         qg_soil          =>    waterstate_vars%qg_soil_col           , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]
         qg_h2osfc        =>    waterstate_vars%qg_h2osfc_col         , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]
         qg               =>    waterstate_vars%qg_col                , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]
         dqgdT            =>    waterstate_vars%dqgdT_col             , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"
         h2osoi_ice       =>    waterstate_vars%h2osoi_ice_col        , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq       =>    waterstate_vars%h2osoi_liq_col        , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)

         grnd_ch4_cond    =>    ch4_vars%grnd_ch4_cond_patch          , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]

         eflx_sh_snow     =>    energyflux_vars%eflx_sh_snow_patch    , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_soil     =>    energyflux_vars%eflx_sh_soil_patch    , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_h2osfc   =>    energyflux_vars%eflx_sh_h2osfc_patch  , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_grnd     =>    energyflux_vars%eflx_sh_grnd_patch    , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_tot      =>    energyflux_vars%eflx_sh_tot_patch     , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         taux             =>    energyflux_vars%taux_patch            , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)
         tauy             =>    energyflux_vars%tauy_patch            , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)
         dlrad            =>    energyflux_vars%dlrad_patch           , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]
         ulrad            =>    energyflux_vars%ulrad_patch           , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]
         cgrnds           =>    energyflux_vars%cgrnds_patch          , & ! Output: [real(r8) (:)   ]  deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
         cgrndl           =>    energyflux_vars%cgrndl_patch          , & ! Output: [real(r8) (:)   ]  deriv of soil latent heat flux wrt soil temp [w/m**2/k]
         cgrnd            =>    energyflux_vars%cgrnd_patch           , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]

         t_ref2m          =>    temperature_vars%t_ref2m_patch        , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)
         t_ref2m_r        =>    temperature_vars%t_ref2m_r_patch      , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin)

         q_ref2m          =>    waterstate_vars%q_ref2m_patch         , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)
         rh_ref2m_r       =>    waterstate_vars%rh_ref2m_r_patch      , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)
         rh_ref2m         =>    waterstate_vars%rh_ref2m_patch        , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)

         z0mg_col         =>    frictionvel_vars%z0mg_col             , & ! Output: [real(r8) (:)   ]  roughness length, momentum [m]
         z0hg_col         =>    frictionvel_vars%z0hg_col             , & ! Output: [real(r8) (:)   ]  roughness length, sensible heat [m]
         z0qg_col         =>    frictionvel_vars%z0qg_col             , & ! Output: [real(r8) (:)   ]  roughness length, latent heat [m]
         ram1             =>    frictionvel_vars%ram1_patch           , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)

         qflx_ev_snow     =>    waterflux_vars%qflx_ev_snow_patch     , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]
         qflx_ev_soil     =>    waterflux_vars%qflx_ev_soil_patch     , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]
         qflx_ev_h2osfc   =>    waterflux_vars%qflx_ev_h2osfc_patch   , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]
         qflx_evap_soi    =>    waterflux_vars%qflx_evap_soi_patch    , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_evap_tot    =>    waterflux_vars%qflx_evap_tot_patch    , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
         begp             =>    bounds%begp                           , &
         endp             =>    bounds%endp                             &
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

         displa(p) = 0._r8
         dlrad(p)  = 0._r8
         ulrad(p)  = 0._r8

         ur(p)    = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)))
         dth(p)   = thm(p)-t_grnd(c)
         dqh(p)   = forc_q(t) - qg(c)
         dthv     = dth(p)*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*dqh(p)
         zldis(p) = forc_hgt_u_patch(p)

         ! Copy column roughness to local pft-level arrays

         z0mg_patch(p) = z0mg_col(c)
         z0hg_patch(p) = z0hg_col(c)
         z0qg_patch(p) = z0qg_col(c)

         ! Initialize Monin-Obukhov length and wind speed

         call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_patch(p), um(p), obu(p))

      end do

      ! Perform stability iteration
      ! Determine friction velocity, and potential temperature and humidity
      ! profiles of the surface boundary layer

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
            g = veg_pp%gridcell(p)

            tstar = temp1(p)*dth(p)
            qstar = temp2(p)*dqh(p)
            z0hg_patch(p) = z0mg_patch(p)/exp(0.13_r8 * (ustar(p)*z0mg_patch(p)/1.5e-5_r8)**0.45_r8)
            z0qg_patch(p) = z0hg_patch(p)
            thvstar = tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*qstar
            zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(p) = max(ur(p),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_r8
               um(p) = sqrt(ur(p)*ur(p) + wc*wc)
            end if
            obu(p) = zldis(p)/zeta
         end do

      end do ! end stability iteration

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
         taux(p)          = -forc_rho(t)*forc_u(t)/ram
         tauy(p)          = -forc_rho(t)*forc_v(t)/ram
         eflx_sh_grnd(p)  = -raih*dth(p)
         eflx_sh_tot(p)   = eflx_sh_grnd(p)

         ! compute sensible heat fluxes individually
         eflx_sh_snow(p)   = -raih*(thm(p)-t_soisno(c,snl(c)+1))
         eflx_sh_soil(p)   = -raih*(thm(p)-t_soisno(c,1))
         eflx_sh_h2osfc(p) = -raih*(thm(p)-t_h2osfc(c))

         ! water fluxes from soil
         qflx_evap_soi(p)  = -raiw*dqh(p)
         qflx_evap_tot(p)  = qflx_evap_soi(p)

         ! compute latent heat fluxes individually
         qflx_ev_snow(p)   = -raiw*(forc_q(t) - qg_snow(c))
         qflx_ev_soil(p)   = -raiw*(forc_q(t) - qg_soil(c))
         qflx_ev_h2osfc(p) = -raiw*(forc_q(t) - qg_h2osfc(c))

         ! 2 m height air temperature
         t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

         ! 2 m height specific humidity
         q_ref2m(p) = forc_q(t) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

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

end module BareGroundFluxesMod
