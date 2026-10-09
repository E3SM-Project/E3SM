module BareGroundFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute sensible and latent fluxes and their derivatives with respect
  ! to ground temperature using ground temperatures from previous time step.
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use CH4Mod               , only : ch4_type
  use CanopyStateType      , only : canopystate_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use SoilStateType        , only : soilstate_type
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
  subroutine BareGroundFluxes(bounds, num_nolu_barep, filter_nolu_barep, &
       canopystate_vars, soilstate_vars, &
       frictionvel_vars, ch4_vars)
    !
    ! !DESCRIPTION:
    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step.
    ! !USES:
    use shr_const_mod        , only : SHR_CONST_RGAS
    use shr_flux_mod         , only : shr_flux_update_stress
    use decompMod            , only : bounds_type
    use elm_varpar           , only : nlevgrnd
    use elm_varcon           , only : cpair, vkc, grav, denice, denh2o
    use elm_varctl           , only : iulog, use_lch4
    use landunit_varcon      , only : istsoil, istcrop
    use FrictionVelocityMod  , only : FrictionVelocity, MoninObukIni, &
         implicit_stress, atm_gustiness, force_land_gustiness
    use QSatMod              , only : QSat
    use SurfaceResistanceMod , only : do_soilevap_beta
    use elm_time_manager     , only : get_nstep
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_nolu_barep        ! number of pft non-lake, non-urban points in pft filter
    integer                , intent(in)    :: filter_nolu_barep(:)   ! patch filter for non-lake, non-urban bare pfts
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    type(ch4_type)         , intent(inout) :: ch4_vars
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: dtaumin = 0.01_r8     ! max limit for stress convergence [Pa]
    integer, parameter  :: itmin = 3             ! minimum number of iterations
    integer, parameter  :: itmax = 30            ! maximum number of iterations
    real(r8),PARAMETER :: beta = 1.0_r8   ! coefficient of convective velocity [-]
    integer  :: p,c,t,g,j,l                      ! indices
    integer  :: f                                ! filter pft index
    integer  :: fa                               ! active (unconverged) filter index
    integer  :: fn                               ! number of active (unconverged) patches
    integer  :: fnold                            ! previous number of active patches
    integer  :: begp, endp                       ! patch bounds
    integer  :: iter                             ! iteration index
    integer  :: iter_final                       ! number of iterations used
    integer  :: loopmax                          ! maximum number of iterations for this configuration
    real(r8) :: displa(bounds%begp:bounds%endp)     ! displacement height [m]
    real(r8) :: z0mg_patch(bounds%begp:bounds%endp) ! roughness length, momentum [m]
    real(r8) :: z0hg_patch(bounds%begp:bounds%endp) ! roughness length, sensible heat [m]
    real(r8) :: z0qg_patch(bounds%begp:bounds%endp) ! roughness length, latent heat [m]
    real(r8) :: zldis(num_nolu_barep)   ! reference height "minus" zero displacement height [m]
    real(r8) :: dth(num_nolu_barep)     ! diff of virtual temp. between ref. height and surface
    real(r8) :: dqh(num_nolu_barep)     ! diff of humidity between ref. height and surface
    real(r8) :: obu(num_nolu_barep)     ! Monin-Obukhov length (m)
    real(r8) :: ur(num_nolu_barep)      ! wind speed at reference height [m/s]
    real(r8) :: um(num_nolu_barep)      ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(num_nolu_barep)   ! relation for potential temperature profile
    real(r8) :: temp12m(num_nolu_barep) ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(num_nolu_barep)   ! relation for specific humidity profile
    real(r8) :: temp22m(num_nolu_barep) ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(num_nolu_barep)   ! friction velocity [m/s]
    real(r8) :: fm(num_nolu_barep)      ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: ugust_total(num_nolu_barep)    ! gustiness including convective velocity [m/s]
    real(r8) :: wind_speed0(num_nolu_barep)    ! Wind speed from atmosphere at start of iteration
    real(r8) :: wind_speed_adj(num_nolu_barep) ! Adjusted wind speed for iteration
    real(r8) :: tau(num_nolu_barep)            ! Stress used in iteration
    real(r8) :: tau_diff(num_nolu_barep)       ! Difference from previous iteration tau
    real(r8) :: prev_tau(num_nolu_barep)       ! Previous iteration tau
    real(r8) :: prev_tau_diff(num_nolu_barep)  ! Previous difference in iteration tau
    ! iter_filterp(1:fn) holds the patch index of patches still iterating;
    ! iter_filter_map(1:fn) holds their position in filter_nolu_barep, which is
    ! how the num_nolu_barep-sized work arrays above are addressed.
    integer  :: iter_filterp(num_nolu_barep), iter_filter_map(num_nolu_barep)
    real(r8) :: zeta   ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc     ! convective velocity [m/s]
    real(r8) :: dthv   ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: tstar   ! temperature scaling parameter
    real(r8) :: qstar   ! moisture scaling parameter
    real(r8) :: thvstar ! virtual potential temperature scaling parameter
    real(r8) :: ram     ! aerodynamical resistance [s/m]
    real(r8) :: rah     ! thermal resistance [s/m]
    real(r8) :: raw     ! moisture resistance [s/m]
    real(r8) :: raih    ! temporary variable [kg/m2/s]
    real(r8) :: raiw    ! temporary variable [kg/m2/s]
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
         wsresp           =>    top_as%wsresp                         , & ! Input:  [real(r8) (:)   ]  response of wind to surface stress (m/s/Pa)
         tau_est          =>    top_as%tau_est                        , & ! Input:  [real(r8) (:)   ]  approximate atmosphere change to zonal wind (m/s)
         ugust            =>    top_as%ugust                          , & ! Input:  [real(r8) (:)   ]  gustiness from atmosphere (m/s)
         forc_th          =>    top_as%thbot                          , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)
         forc_pbot        =>    top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_rho         =>    top_as%rhobot                         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
         forc_q           =>    top_as%qbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)

         forc_hgt_u_patch => frictionvel_vars%forc_hgt_u_patch , & ! Input:
         forc_hgt_t_patch => frictionvel_vars%forc_hgt_t_patch , & ! Input:  [real(r8) (:) ] observational height of temperature at pft level [m]
         forc_hgt_q_patch => frictionvel_vars%forc_hgt_q_patch , & ! Input:  [real(r8) (:) ] observational height of specific humidity at pft level [m]
         vds              => frictionvel_vars%vds_patch        , & ! Output: [real(r8) (:) ] dry deposition velocity term (m/s) (for SO4 NH4NO3)
         u10              => frictionvel_vars%u10_patch        , & ! Output: [real(r8) (:) ] 10-m wind (m/s) (for dust model)
         u10_elm          => frictionvel_vars%u10_elm_patch    , & ! Output: [real(r8) (:) ] 10-m wind (m/s)
         u10_with_gusts_elm=>frictionvel_vars%u10_with_gusts_elm_patch, & ! Output: [real(r8) (:) ] 10-m wind with gusts(m/s)
         va               => frictionvel_vars%va_patch         , & ! Output: [real(r8) (:) ] atmospheric wind speed plus convective velocity (m/s)
         fv               => frictionvel_vars%fv_patch         ,  & ! Output: [real(r8) (:) ] friction velocity (m/s) (for dust model)

         frac_veg_nosno   =>    canopystate_vars%frac_veg_nosno_patch , & ! Input:  [logical  (:)   ]  true=> pft is bare ground (elai+esai = zero)

         htvp             =>    col_ef%htvp             , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg]

         watsat           =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         soilbeta         =>    soilstate_vars%soilbeta_col           , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity

         t_soisno         =>    col_es%t_soisno   , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_grnd           =>    col_es%t_grnd     , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]
         thv              =>    col_es%thv        , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)
         thm              =>    veg_es%thm        , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
         t_h2osfc         =>    col_es%t_h2osfc   , & ! Input:  [real(r8) (:)   ]  surface water temperature

         frac_sno         =>    col_ws%frac_sno   , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         qg_snow          =>    col_ws%qg_snow    , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]
         qg_soil          =>    col_ws%qg_soil    , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]
         qg_h2osfc        =>    col_ws%qg_h2osfc  , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]
         qg               =>    col_ws%qg         , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]
         dqgdT            =>    col_ws%dqgdT      , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"
         h2osoi_ice       =>    col_ws%h2osoi_ice , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq       =>    col_ws%h2osoi_liq , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)

         grnd_ch4_cond    =>    ch4_vars%grnd_ch4_cond_patch , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]

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

         q_ref2m          =>    veg_ws%q_ref2m        , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)
         rh_ref2m_r       =>    veg_ws%rh_ref2m_r     , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)
         rh_ref2m         =>    veg_ws%rh_ref2m       , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)

         z0mg_col         =>    frictionvel_vars%z0mg_col   , & ! Output: [real(r8) (:)   ]  roughness length, momentum [m]
         z0hg_col         =>    frictionvel_vars%z0hg_col   , & ! Output: [real(r8) (:)   ]  roughness length, sensible heat [m]
         z0qg_col         =>    frictionvel_vars%z0qg_col   , & ! Output: [real(r8) (:)   ]  roughness length, latent heat [m]
         ram1             =>    frictionvel_vars%ram1_patch , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)

         qflx_ev_snow     =>    veg_wf%qflx_ev_snow     , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]
         qflx_ev_soil     =>    veg_wf%qflx_ev_soil     , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]
         qflx_ev_h2osfc   =>    veg_wf%qflx_ev_h2osfc   , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]
         qflx_evap_soi    =>    veg_wf%qflx_evap_soi    , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_evap_tot    =>    veg_wf%qflx_evap_tot    , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg
         num_iter         => frictionvel_vars%num_iter_patch             & ! Output: number of iterations required
         )

      !---------------------------------------------------
      ! Filter patches where frac_veg_nosno IS ZERO
      !---------------------------------------------------

      if (num_nolu_barep == 0) return

      begp = bounds%begp
      endp = bounds%endp

      if (implicit_stress) then
         loopmax = itmax
      else
         loopmax = itmin
      end if

      !$acc enter data create(displa(:), z0mg_patch(:), z0hg_patch(:), z0qg_patch(:), &
      !$acc    zldis(:), dth(:), dqh(:), obu(:), ur(:), um(:), temp1(:), temp12m(:), &
      !$acc    temp2(:), temp22m(:), ustar(:), fm(:), ugust_total(:), wind_speed0(:), &
      !$acc    wind_speed_adj(:), tau(:), tau_diff(:), prev_tau(:), prev_tau_diff(:), &
      !$acc    iter_filterp(:), iter_filter_map(:))

      ! Compute sensible and latent fluxes and their derivatives with respect
      ! to ground temperature using ground temperatures from previous time step
      !$acc parallel loop independent gang vector default(present) private(p,c,t,dthv)
      do f = 1, num_nolu_barep
         p = filter_nolu_barep(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)

         iter_filterp(f)    = p
         iter_filter_map(f) = f

         ! Initialization variables

         displa(p) = 0._r8
         dlrad(p)  = 0._r8
         ulrad(p)  = 0._r8

         ! Initialize winds for iteration.
         if (implicit_stress) then
            wind_speed0(f) = max(0.01_r8, hypot(forc_u(t), forc_v(t)))
            wind_speed_adj(f) = wind_speed0(f)
            ur(f) = max(1.0_r8, sqrt(wind_speed_adj(f)**2 + ugust(t)**2))

            prev_tau(f) = tau_est(t)
         else
            ur(f) = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)+ugust(t)*ugust(t)))
         end if
         tau_diff(f) = 1.e100_r8
         ugust_total(f) = ugust(t)

         dth(f)   = thm(p)-t_grnd(c)
         dqh(f)   = forc_q(t) - qg(c)
         dthv     = dth(f)*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*dqh(f)
         zldis(f) = forc_hgt_u_patch(p)

         ! Copy column roughness to local pft-level arrays

         z0mg_patch(p) = z0mg_col(c)
         z0hg_patch(p) = z0hg_col(c)
         z0qg_patch(p) = z0qg_col(c)

         ! Initialize Obukhov length scale and wind speed

         call MoninObukIni(ur(f), thv(c), dthv, zldis(f), z0mg_patch(p), um(f), obu(f))
         num_iter(p) = 0._r8
      end do

      ! Perform stability iteration
      ! Determine friction velocity, and potential temperature and humidity
      ! profiles of the surface boundary layer

      fn = num_nolu_barep
      iter_final = 0

      ITERATION: do iter = 1, loopmax

         call FrictionVelocity(begp, endp, fn, iter_filterp, iter_filter_map, num_nolu_barep, &
              displa(begp:endp), z0mg_patch(begp:endp), z0hg_patch(begp:endp), z0qg_patch(begp:endp), &
              obu(1:num_nolu_barep), iter, ur(1:num_nolu_barep), um(1:num_nolu_barep), &
              ugust_total(1:num_nolu_barep), ustar(1:num_nolu_barep), &
              temp1(1:num_nolu_barep), temp2(1:num_nolu_barep), temp12m(1:num_nolu_barep), &
              temp22m(1:num_nolu_barep), fm(1:num_nolu_barep), &
              frictionvel_vars)

         !$acc parallel loop independent gang vector default(present) &
         !$acc    private(p,f,c,t,ram,tstar,qstar,thvstar,zeta,wc)
         do fa = 1, fn
            p = iter_filterp(fa)
            f = iter_filter_map(fa)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)

            ! Calculate magnitude of stress and update wind speed.
            if (implicit_stress) then
               ram = 1._r8/(ustar(f)*ustar(f)/um(f))
               tau(f) = forc_rho(t)*wind_speed_adj(f)/ram
               call shr_flux_update_stress(wind_speed0(f), wsresp(t), tau_est(t), &
                    tau(f), prev_tau(f), tau_diff(f), prev_tau_diff(f), &
                    wind_speed_adj(f))
               ur(f) = max(1.0_r8, sqrt(wind_speed_adj(f)**2 + ugust(t)**2))
            end if

            tstar = temp1(f)*dth(f)
            qstar = temp2(f)*dqh(f)
            z0hg_patch(p) = z0mg_patch(p)/exp(0.13_r8 * (ustar(f)*z0mg_patch(p)/1.5e-5_r8)**0.45_r8)
            z0qg_patch(p) = z0hg_patch(p)
            thvstar = tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*qstar
            zeta = zldis(f)*vkc*grav*thvstar/(ustar(f)**2*thv(c))

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(f) = max(ur(f),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               if ((.not. atm_gustiness) .or. force_land_gustiness) then
                  wc = beta*(-grav*ustar(f)*thvstar*zii(c)/thv(c))**0.333_r8
                  ugust_total(f) = sqrt(ugust(t)**2 + wc**2)
                  um(f) = sqrt(ur(f)*ur(f) + wc*wc)
               else
                  um(f) = max(ur(f),0.1_r8)
               end if
            end if
            obu(f) = zldis(f)/zeta
         end do

         ! Test for convergence.
         ! Compact iter_filterp/iter_filter_map in place, keeping only the patches that
         ! have NOT yet converged. Sequential write-index dependency on fn.
         iter_final = iter
         if (iter >= itmin) then
            fnold = fn
            fn = 0
            !$acc parallel loop seq default(present) private(p,f) copy(fn)
            do fa = 1, fnold
               p = iter_filterp(fa)
               f = iter_filter_map(fa)
               num_iter(p) = real(iter,r8)
               if (.not. (abs(tau_diff(f)) < dtaumin)) then
                  fn = fn + 1
                  iter_filterp(fn)    = p
                  iter_filter_map(fn) = f
               end if
            end do
            if (fn == 0) exit ITERATION
         end if

      end do ITERATION ! end stability iteration

      !$acc parallel loop independent gang vector default(present) &
      !$acc    private(p,c,t,ram,rah,raw,raih,raiw,www,e_ref2m,de2mdT,qsat_ref2m,dqsat2mdT)
      do f = 1, num_nolu_barep
         p = filter_nolu_barep(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)

         ! Determine aerodynamic resistances

         ram  = 1._r8/(ustar(f)*ustar(f)/um(f))
         rah  = 1._r8/(temp1(f)*ustar(f))
         raw  = 1._r8/(temp2(f)*ustar(f))
         raih = forc_rho(t)*cpair/rah
         if (use_lch4) then
            grnd_ch4_cond(p) = 1._r8/raw
         end if

         ! Soil evaporation resistance
         www = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
         www = min(max(www,0.0_r8),1._r8)

         !changed by K.Sakaguchi. Soilbeta is used for evaporation
         if (dqh(f) > 0._r8) then  !dew  (beta is not applied, just like rsoil used to be)
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
         if (implicit_stress) then
            taux(p)          = taux(p) * (wind_speed_adj(f) / wind_speed0(f))
            tauy(p)          = tauy(p) * (wind_speed_adj(f) / wind_speed0(f))
         end if
         eflx_sh_grnd(p)  = -raih*dth(f)
         eflx_sh_tot(p)   = eflx_sh_grnd(p)

         ! compute sensible heat fluxes individually
         eflx_sh_snow(p)   = -raih*(thm(p)-t_soisno(c,snl(c)+1))
         eflx_sh_soil(p)   = -raih*(thm(p)-t_soisno(c,1))
         eflx_sh_h2osfc(p) = -raih*(thm(p)-t_h2osfc(c))

         ! water fluxes from soil
         qflx_evap_soi(p)  = -raiw*dqh(f)
         qflx_evap_tot(p)  = qflx_evap_soi(p)

         ! compute latent heat fluxes individually
         qflx_ev_snow(p)   = -raiw*(forc_q(t) - qg_snow(c))
         qflx_ev_soil(p)   = -raiw*(forc_q(t) - qg_soil(c))
         qflx_ev_h2osfc(p) = -raiw*(forc_q(t) - qg_h2osfc(c))

         ! 2 m height air temperature
         t_ref2m(p) = thm(p) + temp1(f)*dth(f)*(1._r8/temp12m(f) - 1._r8/temp1(f))

         ! 2 m height specific humidity
         q_ref2m(p) = forc_q(t) + temp2(f)*dqh(f)*(1._r8/temp22m(f) - 1._r8/temp2(f))

         ! 2 m height relative humidity
         call QSat(t_ref2m(p), forc_pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)

         if (veg_pp%is_on_soil_col(p) .or. veg_pp%is_on_crop_col(p)) then
            rh_ref2m_r(p) = rh_ref2m(p)
            t_ref2m_r(p) = t_ref2m(p)
         end if
      end do

#ifndef _OPENACC
      ! Check for convergence of stress.
      if (implicit_stress .and. get_nstep() > 0) then ! Suppress common warnings on the first time step.
         do f = 1, num_nolu_barep
            p = filter_nolu_barep(f)
            if (abs(tau_diff(f)) > dtaumin) then
               write(iulog,*)'WARNING: Stress did not converge for bare ground ',&
                    ' nstep = ',get_nstep(),' p= ',p,' prev_tau_diff= ',prev_tau_diff(f),&
                    ' tau_diff= ',tau_diff(f),' tau= ',tau(f),&
                    ' wind_speed_adj= ',wind_speed_adj(f),' iter_final= ',iter_final
            end if
         end do
      end if
#endif

      !$acc exit data delete(displa(:), z0mg_patch(:), z0hg_patch(:), z0qg_patch(:), &
      !$acc    zldis(:), dth(:), dqh(:), obu(:), ur(:), um(:), temp1(:), temp12m(:), &
      !$acc    temp2(:), temp22m(:), ustar(:), fm(:), ugust_total(:), wind_speed0(:), &
      !$acc    wind_speed_adj(:), tau(:), tau_diff(:), prev_tau(:), prev_tau_diff(:), &
      !$acc    iter_filterp(:), iter_filter_map(:))

    end associate

  end subroutine BareGroundFluxes

end module BareGroundFluxesMod
