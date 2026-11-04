module UrbanFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
  !
  ! !USES:
  use elm_varctl, only : elm_device_ready  ! compute constructs run on the host before the device copy
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use decompMod            , only : bounds_type
  use elm_varcon           , only : namel
  use elm_varctl           , only : iulog
  use abortutils           , only : endrun
  use UrbanParamsType      , only : urbanparams_type
  use UrbanParamsType      , only : urban_wasteheat_int, urban_hac_on_int, urban_hac_int,urban_hac_off_int
  use SoilStateType        , only : soilstate_type
  use FrictionVelocityType , only : frictionvel_type
  use TopounitDataType     , only : top_as
  use LandunitType         , only : lun_pp
  use LandunitDataType     , only : lun_es, lun_ef, lun_ws
  use ColumnDataType       , only : col_es, col_ef, col_ws
  use VegetationDataType   , only : veg_es, veg_ef, veg_ws, veg_wf
  use timeinfoMod  , only : nstep_mod

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanFluxes       ! Urban physics - turbulent fluxes
  !-----------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------
  subroutine UrbanFluxes (bounds,       &
       num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       urbanparams_vars, soilstate_vars,    &
      frictionvel_vars )
    !
    ! !DESCRIPTION:
    ! Turbulent and momentum fluxes from urban canyon (consisting of roof, sunwall,
    ! shadewall, pervious and impervious road).

    ! !USES:
    use shr_flux_mod         , only : shr_flux_update_stress
    use elm_varcon          , only : cpair, vkc, grav, pondmx_urban, rpi
    use elm_varcon          , only : ht_wasteheat_factor, ac_wasteheat_factor, wasteheat_limit
    use FrictionVelocityMod , only : FrictionVelocity, MoninObukIni, &
         implicit_stress, atm_gustiness, force_land_gustiness
    use QSatMod             , only : QSat
    use elm_varpar          , only : nlevurb, nlevgrnd
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_urbanl         ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:)   ! urban landunit filter
    integer                , intent(in)    :: num_urbanc         ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:)   ! urban column filter
    integer                , intent(in)    :: num_urbanp         ! number of urban patches in clump
    integer                , intent(in)    :: filter_urbanp(:)   ! urban pft filter
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fl, fa, j, k, c, l, t, p                  ! indices; fa = position in the active list
    integer  :: pu, c0, p0                                ! first patch and first column of a landunit
    integer  :: begl, endl
    integer  :: iter, iter_final, loopmax                 ! iteration index, iterations used, bound of the iteration loop
    integer  :: fnl                                       ! number of active landunits
    integer  :: n_keep, i_next, k_next                    ! compaction of the active list
    ! iter_filterl(1:fnl) holds the landunits still iterating, iter_mapl their position in filter_urbanl,
    ! which is how the work arrays below are addressed.
    integer  :: iter_filterl(1:num_urbanl), iter_mapl(1:num_urbanl)
    integer  :: next_filterl(1:num_urbanl), next_mapl(1:num_urbanl)

    ! The work arrays are only what crosses a kernel boundary. Arguments of FrictionVelocity:
    real(r8) :: obu(1:num_urbanl)                         ! Monin-Obukhov length (m)
    real(r8) :: ur(1:num_urbanl)                          ! wind speed at reference height (m/s)
    real(r8) :: um(1:num_urbanl)                          ! wind speed including the stablity effect (m/s)
    real(r8) :: ugust_total(1:num_urbanl)                 ! gustiness including convective velocity [m/s]
    real(r8) :: ustar(1:num_urbanl)                       ! friction velocity (m/s)
    real(r8) :: temp1(1:num_urbanl)                       ! relation for potential temperature profile
    real(r8) :: temp2(1:num_urbanl)                       ! relation for specific humidity profile
    real(r8) :: temp12m(1:num_urbanl)                     ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp22m(1:num_urbanl)                     ! relation for specific humidity profile applied at 2-m
    real(r8) :: fm(1:num_urbanl)                          ! needed for BGC only to diagnose 10m wind speed
    ! State of the last iteration of a landunit:
    real(r8) :: canyon_resistance(1:num_urbanl)           ! resistance to heat and moisture transfer from canyon road/walls to canyon air (s/m)
    real(r8) :: ramu(1:num_urbanl)                        ! aerodynamic resistance (s/m)
    real(r8) :: fwet_roof(1:num_urbanl)                   ! fraction of roof surface that is wet (-)
    real(r8) :: fwet_road_imperv(1:num_urbanl)            ! fraction of impervious road surface that is wet (-)
    ! Implicit stress:
    real(r8) :: wind_speed_adj(1:num_urbanl)              ! adjusted wind speed for iteration
    real(r8) :: prev_tau(1:num_urbanl)                    ! previous iteration tau
    real(r8) :: tau_diff(1:num_urbanl)                    ! difference from previous iteration tau
    real(r8) :: prev_tau_diff(1:num_urbanl)               ! previous difference in iteration tau

    ! Scalars of one landunit in the kernels
    real(r8) :: wind0, tau_l                              ! wind speed from atmosphere (m/s), stress used in iteration
    real(r8) :: rah, raw                                  ! thermal, moisture resistance (s/m)
    real(r8) :: rc, gc, gs                                ! canyon resistance, 1/rc, gc/wts_tot
    real(r8) :: thm, thv                                  ! forc_t+0.0098*forc_hgt_t, virtual potential temperature (K)
    real(r8) :: dth_l, dqh_l, dthv                        ! diff of temperature, humidity, vir. poten. temp. between ref. height and surface
    real(r8) :: zldis_l                                   ! reference height "minus" zero displacement height (m)
    real(r8) :: canyontop, canyon_u                       ! wind at canyon top, u-component of wind speed inside canyon (m/s)
    real(r8) :: wts_r, wts_p, wts_i, wts_w                ! scaled sensible heat conductance: roof, pervious road, impervious road, each wall (m/s)
    real(r8) :: wtq_r, wtq_i                              ! scaled latent heat conductance: roof, impervious road (m/s); pervious road = wts_p
    real(r8) :: wtas_l, wtaq_l                            ! sensible, latent heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wts_tot, wtq_tot                          ! sums of wtas_l, wtaq_l and the surface conductances
    real(r8) :: tstar, qstar, thvstar                     ! temperature, moisture, virtual potential temperature scaling parameter
    real(r8) :: wc, zeta                                  ! convective velocity (m/s), dimensionless height in Monin-Obukhov theory
    real(r8) :: rho, rhocp, tx, ty                        ! air density, its product with cpair, stresses of the landunit
    real(r8) :: dth_p, dq_r, dq_i, dq_p                   ! air-surface temperature difference; humidity difference: roof, roads
    real(r8) :: wh_r, wh_s, wh_d, ac_r, ac_s, ac_d        ! waste heat, heat from AC: roof, sunwall, shadewall (W/m**2)
    real(r8) :: lngth_roof                                ! length of roof (m)
    real(r8) :: e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT    ! 2 m height saturated vapor pressure, specific humidity, derivatives
    real(r8) :: rh_l                                      ! 2 m height relative humidity (%)
    real(r8) :: eflx_err_l, qflx_err_l, err_e, err_q      ! flux check residuals, the largest of the failures
    integer  :: erridx_e, erridx_q                        ! a landunit that failed the flux checks
    integer  :: spel_dc_1                                 ! spel defer-checks: elements recorded by deferred checks

    real(r8), parameter :: lapse_rate = 0.0098_r8 ! Dry adiabatic lapse rate (K/m)
    real(r8), parameter :: dtaumin = 0.01_r8      ! max limit for stress convergence [Pa]
    integer, parameter  :: itmin = 3              ! minimum number of iterations
    integer, parameter  :: itmax = 30             ! maximum number of iterations
    real(r8), parameter :: beta = 1._r8           ! coefficient of convective velocity
    real(r8), parameter :: zii  = 1000._r8        ! convective boundary layer height (m)
    ! Each urban landunit has one column (and one patch) per surface in this order: offsets from
    ! coli/pfti, checked at initialization in UrbanParamsType.
    integer, parameter  :: kroof = 0, ksunwall = 1, kshadewall = 2, kroad_imperv = 3, kroad_perv = 4
    integer, parameter  :: nsurf = 5
    !-----------------------------------------------------------------------

    associate(                                                                &
         z_0_town            =>   lun_pp%z_0_town                              , & ! Input:  [real(r8) (:)   ]  momentum roughness length of urban landunit (m)
         z_d_town            =>   lun_pp%z_d_town                              , & ! Input:  [real(r8) (:)   ]  displacement height of urban landunit (m)
         ht_roof             =>   lun_pp%ht_roof                               , & ! Input:  [real(r8) (:)   ]  height of urban roof (m)
         wtlunit_roof        =>   lun_pp%wtlunit_roof                          , & ! Input:  [real(r8) (:)   ]  weight of roof with respect to landunit
         canyon_hwr          =>   lun_pp%canyon_hwr                            , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width
         wtroad_perv         =>   lun_pp%wtroad_perv                           , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road

         forc_t              =>   top_as%tbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (K)
         forc_th             =>   top_as%thbot                              , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (K)
         forc_rho            =>   top_as%rhobot                             , & ! Input:  [real(r8) (:)   ]  air density (kg/m**3)
         forc_q              =>   top_as%qbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)
         forc_pbot           =>   top_as%pbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_u              =>   top_as%ubot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)
         forc_v              =>   top_as%vbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)
         wsresp              =>   top_as%wsresp                             , & ! Input:  [real(r8) (:)   ]  response of wind to surface stress (m/s/Pa)
         tau_est             =>   top_as%tau_est                            , & ! Input:  [real(r8) (:)   ]  approximate atmosphere change to zonal wind (m/s)
         ugust               =>   top_as%ugust                              , & ! Input:  [real(r8) (:)   ]  gustiness from atmosphere (m/s)

         wind_hgt_canyon     =>   urbanparams_vars%wind_hgt_canyon          , & ! Input:  [real(r8) (:)   ]  height above road at which wind in canyon is to be computed (m)
         eflx_traffic_factor =>   urbanparams_vars%eflx_traffic_factor      , & ! Input:  [real(r8) (:)   ]  multiplicative urban traffic factor for sensible heat flux

         rootr_road_perv     =>   soilstate_vars%rootr_road_perv_col        , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer for urban pervious road
         soilalpha_u         =>   soilstate_vars%soilalpha_u_col            , & ! Input:  [real(r8) (:)   ]  Urban factor that reduces ground saturated specific humidity (-)
         rootr               =>   soilstate_vars%rootr_patch                , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer

         t_grnd              =>   col_es%t_grnd               , & ! Input:  [real(r8) (:)   ]  ground surface temperature (K)
         t_soisno            =>   col_es%t_soisno             , & ! Input:  [real(r8) (:,:) ]  soil temperature (K)
         t_ref2m             =>   veg_es%t_ref2m            , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (K)
         t_ref2m_u           =>   veg_es%t_ref2m_u          , & ! Output: [real(r8) (:)   ]  Urban 2 m height surface air temperature (K)
         t_veg               =>   veg_es%t_veg                , & ! Output: [real(r8) (:)   ]  vegetation temperature (K)
         t_building          =>   lun_es%t_building           , & ! Output: [real(r8) (:)   ]  internal building temperature (K)
         taf                 =>   lun_es%taf                  , & ! Output: [real(r8) (:)   ]  urban canopy air temperature (K)

         frac_sno            =>   col_ws%frac_sno              , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         snow_depth          =>   col_ws%snow_depth            , & ! Input:  [real(r8) (:)   ]  snow height (m)
         dqgdT               =>   col_ws%dqgdT                 , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"
         qg                  =>   col_ws%qg                    , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface (kg/kg)
         h2osoi_ice          =>   col_ws%h2osoi_ice            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq          =>   col_ws%h2osoi_liq            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         qaf                 =>   lun_ws%qaf                   , & ! Output: [real(r8) (:)   ]  urban canopy air specific humidity (kg/kg)
         q_ref2m             =>   veg_ws%q_ref2m             , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)
         rh_ref2m            =>   veg_ws%rh_ref2m            , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)
         rh_ref2m_u          =>   veg_ws%rh_ref2m_u          , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)

         forc_hgt_u_patch    =>   frictionvel_vars%forc_hgt_u_patch         , & ! Input:  [real(r8) (:)   ]  observational height of wind at pft-level (m)
         forc_hgt_t_patch    =>   frictionvel_vars%forc_hgt_t_patch         , & ! Input:  [real(r8) (:)   ]  observational height of temperature at pft-level (m)
         ram1                =>   frictionvel_vars%ram1_patch               , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)

         htvp                =>   col_ef%htvp                  , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) (J/kg)
         eflx_urban_ac       =>   col_ef%eflx_urban_ac         , & ! Input:  [real(r8) (:)   ]  urban air conditioning flux (W/m**2)
         eflx_urban_heat     =>   col_ef%eflx_urban_heat       , & ! Input:  [real(r8) (:)   ]  urban heating flux (W/m**2)
         dlrad               =>   veg_ef%dlrad               , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy (W/m**2)
         ulrad               =>   veg_ef%ulrad               , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy (W/m**2)
         cgrnds              =>   veg_ef%cgrnds              , & ! Output: [real(r8) (:)   ]  deriv, of soil sensible heat flux wrt soil temp (W/m**2/K)
         cgrndl              =>   veg_ef%cgrndl              , & ! Output: [real(r8) (:)   ]  deriv of soil latent heat flux wrt soil temp (W/m**2/K)
         cgrnd               =>   veg_ef%cgrnd               , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp (W/m**2/K)
         eflx_sh_grnd        =>   veg_ef%eflx_sh_grnd        , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_tot         =>   veg_ef%eflx_sh_tot         , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_tot_u       =>   veg_ef%eflx_sh_tot_u       , & ! Output: [real(r8) (:)   ]  urban total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_snow        =>   veg_ef%eflx_sh_snow        , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_soil        =>   veg_ef%eflx_sh_soil        , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_h2osfc      =>   veg_ef%eflx_sh_h2osfc      , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_traffic        =>   lun_ef%eflx_traffic          , & ! Output: [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)
         eflx_wasteheat      =>   lun_ef%eflx_wasteheat        , & ! Output: [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_heat_from_ac   =>   lun_ef%eflx_heat_from_ac     , & ! Output: [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         taux                =>   veg_ef%taux                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)
         tauy                =>   veg_ef%tauy                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)

         qflx_evap_soi       =>   veg_wf%qflx_evap_soi        , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_tran_veg       =>   veg_wf%qflx_tran_veg        , & ! Output: [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         qflx_evap_veg       =>   veg_wf%qflx_evap_veg          & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)

         )
    !$acc enter data create(&
    !$acc iter_filterl(:), &
    !$acc iter_mapl(:), &
    !$acc next_filterl(:), &
    !$acc next_mapl(:), &
    !$acc obu(:), &
    !$acc ur(:), &
    !$acc um(:), &
    !$acc ugust_total(:), &
    !$acc ustar(:), &
    !$acc temp1(:), &
    !$acc temp2(:), &
    !$acc temp12m(:), &
    !$acc temp22m(:), &
    !$acc fm(:), &
    !$acc canyon_resistance(:), &
    !$acc ramu(:), &
    !$acc fwet_roof(:), &
    !$acc fwet_road_imperv(:), &
    !$acc wind_speed_adj(:), &
    !$acc prev_tau(:), &
    !$acc tau_diff(:), &
    !$acc prev_tau_diff(:) )


       begl = bounds%begl
       endl = bounds%endl

      ! Compute fluxes - Follows elm approach for bare soils (Oleson et al 2004)

      ! Initialize winds, the active list and the Monin-Obukhov length
      !$acc parallel loop independent gang vector default(present) &
      !$acc private(l, t, pu, wind0, thm, thv, dth_l, dqh_l, dthv, zldis_l) if(elm_device_ready)
      do fl = 1, num_urbanl
         l  = filter_urbanl(fl)
         t  = lun_pp%topounit(l)
         pu = lun_pp%pfti(l)

         if (implicit_stress) then
            wind0 = max(0.01_r8, hypot(forc_u(t), forc_v(t)))
            wind_speed_adj(fl) = wind0
            ur(fl) = max(1.0_r8, sqrt(wind_speed_adj(fl)**2 + ugust(t)**2))
            prev_tau(fl) = tau_est(t)
         else
            ur(fl) = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)+ugust(t)*ugust(t)))
         end if
         tau_diff(fl)     = 1.e100_r8
         ugust_total(fl)  = ugust(t)
         iter_filterl(fl) = l
         iter_mapl(fl)    = fl

         thm     = forc_t(t) + lapse_rate*forc_hgt_t_patch(pu)
         thv     = forc_th(t)*(1._r8+0.61_r8*forc_q(t))
         dth_l   = thm-taf(l)
         dqh_l   = forc_q(t)-qaf(l)
         dthv    = dth_l*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*dqh_l
         zldis_l = forc_hgt_u_patch(pu) - z_d_town(l)

         ! Initialize Monin-Obukhov length and wind speed including convective velocity
         call MoninObukIni(ur(fl), thv, dthv, zldis_l, z_0_town(l), um(fl), obu(fl))
      end do

      ! Start stability iteration

      if (implicit_stress) then
         loopmax = itmax
      else
         loopmax = itmin
      end if
      fnl = num_urbanl
      ITERATION: do iter = 1, loopmax
         ! Get friction velocity, relation for potential
         ! temperature and humidity profiles of surface boundary layer.
         call FrictionVelocity(begl, endl, &
                 fnl, iter_filterl, iter_mapl, num_urbanl, &
                 z_d_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), &
                 obu(1:num_urbanl), iter, ur(1:num_urbanl), um(1:num_urbanl), &
                 ugust_total(1:num_urbanl), ustar(1:num_urbanl), &
                 temp1(1:num_urbanl), temp2(1:num_urbanl), temp12m(1:num_urbanl), &
                 temp22m(1:num_urbanl), fm(1:num_urbanl), &
                 frictionvel_vars, landunit_index=.true.)

         ! One landunit at a time: canyon resistance, conductances of its five surfaces (columns
         ! coli(l)+k*, see UrbanParamsType's layout check), new canopy air temperature and humidity
         ! and the new stability. Sums keep the column order.
         !$acc parallel loop independent gang vector default(present) &
         !$acc private(l, fl, t, pu, c0, wind0, tau_l, rah, raw, thm, thv, dth_l, dqh_l, zldis_l, &
         !$acc   canyontop, canyon_u, rc, wts_r, wts_p, wts_i, wts_w, wtq_r, wtq_i, &
         !$acc   tstar, qstar, thvstar, zeta, wc) if(elm_device_ready)
         do fa = 1, fnl
            l  = iter_filterl(fa)
            fl = iter_mapl(fa)
            t  = lun_pp%topounit(l)
            pu = lun_pp%pfti(l)
            c0 = lun_pp%coli(l)

            ! Aerodynamic resistances to fluxes from urban canopy air to atmosphere
            ramu(fl) = 1._r8/(ustar(fl)*ustar(fl)/um(fl))
            rah = 1._r8/(temp1(fl)*ustar(fl))
            raw = 1._r8/(temp2(fl)*ustar(fl))

            ! Calculate magnitude of stress and update wind speed.
            if (implicit_stress) then
               wind0 = max(0.01_r8, hypot(forc_u(t), forc_v(t)))
               tau_l = forc_rho(t)*wind_speed_adj(fl)/ramu(fl)
               call shr_flux_update_stress(wind0, wsresp(t), tau_est(t), &
                    tau_l, prev_tau(fl), tau_diff(fl), prev_tau_diff(fl), &
                    wind_speed_adj(fl))
               ur(fl) = max(1.0_r8, sqrt(wind_speed_adj(fl)**2 + ugust(t)**2))
            end if

            ! Canyon top wind, then its u component in the canyon
            canyontop = ur(fl) * &
                 log( (ht_roof(l)-z_d_town(l)) / z_0_town(l) ) / &
                 log( (forc_hgt_u_patch(pu)-z_d_town(l)) / z_0_town(l) )

            if (canyon_hwr(l) < 0.5_r8) then  ! isolated roughness flow
               canyon_u = canyontop * exp( -0.5_r8*canyon_hwr(l)* &
                    (1._r8-(wind_hgt_canyon(l)/ht_roof(l))) )
            else if (canyon_hwr(l) < 1.0_r8) then ! wake interference flow
               canyon_u = canyontop * (1._r8+2._r8*(2._r8/rpi - 1._r8)* &
                    (ht_roof(l)/(ht_roof(l)/canyon_hwr(l)) - 0.5_r8)) * &
                    exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
            else  ! skimming flow
               canyon_u = canyontop * (2._r8/rpi) * &
                    exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
            end if

            ! The single resistance from all urban surfaces to urban canopy air, from the magnitude
            ! of the canyon wind (horizontal wind above, vertical from friction velocity; Masson 2000)
            rc = cpair * forc_rho(t) / (11.8_r8 + 4.2_r8*sqrt(canyon_u**2._r8 + ustar(fl)**2._r8))
            canyon_resistance(fl) = rc

            ! Wetness of the surfaces that have one; pervious road is always 1 and walls are always 0
            fwet_roof(fl) = surface_wetness(snow_depth(c0+kroof), &
                 h2osoi_liq(c0+kroof,1)+h2osoi_ice(c0+kroof,1), pondmx_urban, qaf(l), qg(c0+kroof))
            fwet_road_imperv(fl) = surface_wetness(snow_depth(c0+kroad_imperv), &
                 h2osoi_liq(c0+kroad_imperv,1)+h2osoi_ice(c0+kroad_imperv,1), pondmx_urban, qaf(l), qg(c0+kroad_imperv))

            ! Scaled conductances: sensible (wts_*, walls share one value), latent (wtq_*)
            wts_r = wtlunit_roof(l)/rc
            wts_p = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/rc
            wts_i = (1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/rc
            wts_w = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/rc
            wtq_r = fwet_roof(fl)*wts_r
            wtq_i = fwet_road_imperv(fl)*(1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/rc

            ! New urban canopy air temperature and specific humidity
            thm = forc_t(t) + lapse_rate*forc_hgt_t_patch(pu)
            thv = forc_th(t)*(1._r8+0.61_r8*forc_q(t))
            taf(l) = (thm/rah + t_grnd(c0+kroof)*wts_r + t_grnd(c0+ksunwall)*wts_w + t_grnd(c0+kshadewall)*wts_w &
                 + t_grnd(c0+kroad_imperv)*wts_i + t_grnd(c0+kroad_perv)*wts_p) &
                 / (1._r8/rah + wts_r + wts_w + wts_w + wts_i + wts_p)
            qaf(l) = (forc_q(t)/raw + qg(c0+kroof)*wtq_r + qg(c0+kroad_imperv)*wtq_i + qg(c0+kroad_perv)*wts_p) &
                 / (1._r8/raw + wtq_r + wtq_i + wts_p)


            ! Determine stability using new taf and qaf
            ! TODO: Some of these constants replicate what is in FrictionVelocity
            !       and BareGround fluxes should consildate. EBK
            zldis_l = forc_hgt_u_patch(pu) - z_d_town(l)
            dth_l = thm-taf(l)
            dqh_l = forc_q(t)-qaf(l)
            tstar = temp1(fl)*dth_l
            qstar = temp2(fl)*dqh_l

            thvstar = tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*qstar
            zeta = zldis_l*vkc*grav*thvstar/(ustar(fl)**2*thv)

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(fl) = max(ur(fl),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               if ((.not. atm_gustiness) .or. force_land_gustiness) then
                  wc = beta*(-grav*ustar(fl)*thvstar*zii/thv)**0.333_r8
                  ugust_total(fl) = sqrt(ugust(t)**2 + wc**2)
                  um(fl) = sqrt(ur(fl)*ur(fl) + wc*wc)
               else
                  um(fl) = max(ur(fl),0.1_r8)
               end if
            end if
            obu(fl) = zldis_l/zeta
         end do

         ! Test for convergence. Only implicit stress can converge early: with explicit stress
         ! tau_diff never changes. Keep only the landunits that have NOT converged.
         iter_final = iter
         if (implicit_stress .and. iter >= itmin) then
            ! Parallel: the order of the list doesn't matter. Most passes nobody leaves, so count first.
            n_keep = 0
            !$acc parallel loop gang vector default(present) private(fl) reduction(+:n_keep) if(elm_device_ready)
            do fa = 1, fnl
               fl = iter_mapl(fa)
               if (.not. (abs(tau_diff(fl)) < dtaumin)) n_keep = n_keep + 1
            end do
            if (n_keep < fnl) then
               i_next = 0
               !$acc parallel loop gang vector default(present) private(fl, k_next) copy(i_next) if(elm_device_ready)
               do fa = 1, fnl
                  fl = iter_mapl(fa)
                  if (.not. (abs(tau_diff(fl)) < dtaumin)) then
                     !$acc atomic capture
                     i_next = i_next + 1
                     k_next = i_next
                     !$acc end atomic
                     next_filterl(k_next) = iter_filterl(fa)
                     next_mapl(k_next) = iter_mapl(fa)
                  end if
               end do
               !$acc parallel loop gang vector default(present) if(elm_device_ready)
               do fa = 1, n_keep
                  iter_filterl(fa) = next_filterl(fa)
                  iter_mapl(fa)  = next_mapl(fa)
               end do
               fnl = n_keep
            end if
            if (n_keep == 0) exit ITERATION
         end if

      end do ITERATION ! end iteration

      ! Everything that follows from the converged state of a landunit: its five patches get their
      ! fluxes directly (patch pfti(l)+k* is on column coli(l)+k*), no column type tests.
      erridx_e = -1
      erridx_q = -1
      err_e = 0._r8
      err_q = 0._r8
      !$acc parallel loop independent gang vector default(present) &
      !$acc private(l, t, p0, c0, p, c, rc, gc, rah, raw, thm, wts_r, wts_p, wts_i, wts_w, wtq_r, wtq_i, &
      !$acc   wtas_l, wtaq_l, wts_tot, wtq_tot, gs, rho, rhocp, tx, ty, wind0, dth_p, dq_r, dq_i, dq_p, &
      !$acc   qsat_ref2m, e_ref2m, de2mdT, dqsat2mdT, rh_l, eflx_err_l, qflx_err_l, lngth_roof, &
      !$acc   wh_r, wh_s, wh_d, ac_r, ac_s, ac_d) &
      !$acc reduction(max:erridx_e, erridx_q, err_e, err_q) if(elm_device_ready)
      do fl = 1, num_urbanl
         l   = filter_urbanl(fl)
         t   = lun_pp%topounit(l)
         p0  = lun_pp%pfti(l)
         c0  = lun_pp%coli(l)
         rho = forc_rho(t)
         rhocp = rho*cpair

         ! The state of the last iteration of this landunit
         rc  = canyon_resistance(fl)
         gc  = 1._r8/rc                 ! unscaled sensible heat conductance of every surface
         rah = 1._r8/(temp1(fl)*ustar(fl))
         raw = 1._r8/(temp2(fl)*ustar(fl))
         thm = forc_t(t) + lapse_rate*forc_hgt_t_patch(p0)
         wts_r = wtlunit_roof(l)/rc
         wts_p = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/rc
         wts_i = (1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/rc
         wts_w = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/rc
         wtq_r = fwet_roof(fl)*wts_r
         wtq_i = fwet_road_imperv(fl)*(1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/rc
         wtas_l = 1._r8/rah
         wtaq_l = 1._r8/raw
         wts_tot = wtas_l + wts_r + wts_p + wts_i + wts_w + wts_w
         wtq_tot = wtaq_l + wtq_r + wts_p + wtq_i

         ! Derivative of sensible and latent heat fluxes with respect to ground temperature
         gs = gc/wts_tot
         cgrnds(p0+kroof)        = rhocp*(wtas_l + wts_p + wts_i + wts_w + wts_w)*gs
         cgrnds(p0+kroad_perv)   = rhocp*(wtas_l + wts_r + wts_i + wts_w + wts_w)*gs
         cgrnds(p0+kroad_imperv) = rhocp*(wtas_l + wts_r + wts_p + wts_w + wts_w)*gs
         cgrnds(p0+ksunwall)     = rhocp*(wtas_l + wts_r + wts_p + wts_i + wts_w)*gs
         cgrnds(p0+kshadewall)   = cgrnds(p0+ksunwall)
         cgrndl(p0+kroof)        = rho*(wtaq_l + wts_p + wtq_i)*((fwet_roof(fl)*gc)/wtq_tot)*dqgdT(c0+kroof)
         cgrndl(p0+kroad_perv)   = rho*(wtaq_l + wtq_r + wtq_i)*(gc/wtq_tot)*dqgdT(c0+kroad_perv)
         cgrndl(p0+kroad_imperv) = rho*(wtaq_l + wtq_r + wts_p)*((fwet_road_imperv(fl)*gc)/wtq_tot)*dqgdT(c0+kroad_imperv)
         cgrndl(p0+ksunwall)     = 0._r8
         cgrndl(p0+kshadewall)   = 0._r8

         ! Surface fluxes of momentum
         tx = -rho*forc_u(t)/ramu(fl)
         ty = -rho*forc_v(t)/ramu(fl)
         if (implicit_stress) then
            wind0 = max(0.01_r8, hypot(forc_u(t), forc_v(t)))
            tx = tx * (wind_speed_adj(fl) / wind0)
            ty = ty * (wind_speed_adj(fl) / wind0)
         end if

         ! Latent heat: walls have none. Pervious road evaporation is assigned to the soil term if
         ! dew or snow or if no liquid water is available in the soil column, otherwise to transpiration.
         dq_r = qaf(l) - qg(c0+kroof)
         dq_i = qaf(l) - qg(c0+kroad_imperv)
         dq_p = qaf(l) - qg(c0+kroad_perv)
         qflx_evap_soi(p0+kroof)        = -rho*(fwet_roof(fl)*gc)*dq_r
         qflx_evap_soi(p0+kroad_imperv) = -rho*(fwet_road_imperv(fl)*gc)*dq_i
         qflx_evap_soi(p0+ksunwall)     = 0._r8
         qflx_evap_soi(p0+kshadewall)   = 0._r8
         if (dq_p > 0._r8 .or. frac_sno(c0+kroad_perv) > 0._r8 .or. soilalpha_u(c0+kroad_perv) <= 0._r8) then
            qflx_evap_soi(p0+kroad_perv) = -rho*gc*dq_p
            qflx_tran_veg(p0+kroad_perv) = 0._r8
         else
            qflx_evap_soi(p0+kroad_perv) = 0._r8
            qflx_tran_veg(p0+kroad_perv) = -rho*gc*dq_p
         end if
         qflx_evap_veg(p0+kroad_perv) = qflx_tran_veg(p0+kroad_perv)

         ! 2 m height relative humidity: urban canopy air represents 2-m temperature and humidity
         call QSat(taf(l), forc_pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
         rh_l = min(100._r8, qaf(l) / qsat_ref2m * 100._r8)

         !$acc loop seq
         do k = 0, nsurf-1
            p = p0 + k
            c = c0 + k
            ram1(p)  = ramu(fl)  !pass value to global variable
            ! Upward and downward canopy longwave are zero
            ulrad(p) = 0._r8
            dlrad(p) = 0._r8
            cgrnd(p) = cgrnds(p) + cgrndl(p)*htvp(c)
            taux(p)  = tx
            tauy(p)  = ty
            ! Use new canopy air temperature
            dth_p = taf(l) - t_grnd(c)
            eflx_sh_grnd(p)   = -rhocp*gc*dth_p
            eflx_sh_snow(p)   = 0._r8
            eflx_sh_soil(p)   = 0._r8
            eflx_sh_h2osfc(p) = 0._r8
            eflx_sh_tot(p)    = eflx_sh_grnd(p)
            eflx_sh_tot_u(p)  = eflx_sh_tot(p)
            t_ref2m(p)    = taf(l)
            q_ref2m(p)    = qaf(l)
            t_ref2m_u(p)  = taf(l)
            rh_ref2m(p)   = rh_l
            rh_ref2m_u(p) = rh_l
            ! Variables needed by history tape
            t_veg(p) = forc_t(t)
         end do

         ! No roots for urban except for pervious road
         !$acc loop seq
         do j = 1, nlevgrnd
            rootr(p0+kroof,j)        = 0._r8
            rootr(p0+ksunwall,j)     = 0._r8
            rootr(p0+kshadewall,j)   = 0._r8
            rootr(p0+kroad_imperv,j) = 0._r8
            rootr(p0+kroad_perv,j)   = rootr_road_perv(c0+kroad_perv,j)
         end do

         ! Total sensible and latent heat must equal the sum of the scaled fluxes of the surfaces
         eflx_err_l = -rhocp*(wts_r*(taf(l)-t_grnd(c0+kroof)) + wts_w*(taf(l)-t_grnd(c0+ksunwall)) &
              + wts_w*(taf(l)-t_grnd(c0+kshadewall)) + wts_i*(taf(l)-t_grnd(c0+kroad_imperv)) &
              + wts_p*(taf(l)-t_grnd(c0+kroad_perv))) + (rhocp/rah)*(thm - taf(l))
         qflx_err_l = -rho*(wtq_r*dq_r + wtq_i*dq_i + wts_p*dq_p) + (rho/raw)*(forc_q(t) - qaf(l))
         if (abs(eflx_err_l) > 0.01_r8) then
            erridx_e = max(erridx_e, l)
            err_e    = max(err_e, abs(eflx_err_l))
         end if
         ! 4.e-9 kg/m**2/s = 0.01 W/m**2
         if (abs(qflx_err_l) > 4.e-9_r8) then
            erridx_q = max(erridx_q, l)
            err_q    = max(err_q, abs(qflx_err_l))
         end if

         ! Heat from air conditioning, waste heat and traffic: their inputs don't depend on the iteration
         if (urban_hac_int == urban_wasteheat_int) then
            wh_r = ac_wasteheat_factor * eflx_urban_ac(c0+kroof)      + ht_wasteheat_factor * eflx_urban_heat(c0+kroof)
            wh_s = ac_wasteheat_factor * eflx_urban_ac(c0+ksunwall)   + ht_wasteheat_factor * eflx_urban_heat(c0+ksunwall)
            wh_d = ac_wasteheat_factor * eflx_urban_ac(c0+kshadewall) + ht_wasteheat_factor * eflx_urban_heat(c0+kshadewall)
         else
            wh_r = 0._r8; wh_s = 0._r8; wh_d = 0._r8
         end if
         ! If air conditioning on, always replace heat removed with heat into canyon
         if (urban_hac_int == urban_hac_on_int .or. urban_hac_int == urban_wasteheat_int) then
            ac_r = abs(eflx_urban_ac(c0+kroof))
            ac_s = abs(eflx_urban_ac(c0+ksunwall))
            ac_d = abs(eflx_urban_ac(c0+kshadewall))
         else
            ac_r = 0._r8; ac_s = 0._r8; ac_d = 0._r8
         end if
         ! Total is the sum over walls and roof, accounting for their surface areas. Limit wasteheat
         ! to ensure that we don't get any unrealistically strong positive feedbacks due to AC in a warmer climate.
         eflx_wasteheat(l) = min(wtlunit_roof(l)*wh_r + (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(wh_s + wh_d)), wasteheat_limit)
         eflx_heat_from_ac(l) = wtlunit_roof(l)*ac_r + (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(ac_s + ac_d))
         ! Traffic heat flux only comes from impervious road
         eflx_traffic(l) = (1._r8-wtlunit_roof(l))*(1._r8-wtroad_perv(l))*eflx_traffic_factor(l)

         ! Internal building temperature
         lngth_roof = (ht_roof(l)/canyon_hwr(l))*wtlunit_roof(l)/(1._r8-wtlunit_roof(l))
         t_building(l) = (ht_roof(l)*(t_soisno(c0+kshadewall,nlevurb) + t_soisno(c0+ksunwall,nlevurb)) &
              +lngth_roof*t_soisno(c0+kroof,nlevurb))/(2._r8*ht_roof(l)+lngth_roof)
      end do

      if (erridx_e > 0) then
         write(iulog,*)'WARNING:  Total sensible heat does not equal sum of scaled heat fluxes for urban columns ',&
              ' nstep = ',nstep_mod,' landunit= ',erridx_e,' max eflx_err= ',err_e
         write(iulog,*)'elm model is stopping - error is greater than .01 W/m**2'
         call endrun(decomp_index=erridx_e, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
      end if

      if (erridx_q > 0) then
         write(iulog,*)'WARNING:  Total water vapor flux does not equal sum of scaled water vapor fluxes for urban columns ',&
              ' nstep = ',nstep_mod,' landunit= ',erridx_q,' max qflx_err= ',err_q
         write(iulog,*)'elm model is stopping - error is greater than 4.e-9 kg/m**2/s'
         call endrun(decomp_index=erridx_q, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
      end if

      ! Check for convergence of stress.
      if (implicit_stress) then
         spel_dc_1 = 0  ! spel defer-checks
         !$acc parallel loop gang vector default(present) reduction(max:spel_dc_1) if(elm_device_ready)
         do fl = 1, num_urbanl
            if (abs(tau_diff(fl)) > dtaumin) then
               if (nstep_mod > 0) spel_dc_1 = max(spel_dc_1, fl)  ! spel defer-checks
            end if
         end do
         if (spel_dc_1 > 0) then  ! spel defer-checks: the last element where the check held
            !$acc update self(prev_tau_diff, prev_tau, tau_diff, wind_speed_adj) if_present if(elm_device_ready)
            fl = spel_dc_1
            l = filter_urbanl(fl)
            ! the update of the last iteration leaves tau in prev_tau
            write(iulog,*)'WARNING: Stress did not converge for urban columns ',&
                 ' nstep = ',nstep_mod,' indexl= ',l,' prev_tau_diff= ',prev_tau_diff(fl),&
                 ' tau_diff= ',tau_diff(fl),' tau= ',prev_tau(fl),&
                 ' wind_speed_adj= ',wind_speed_adj(fl),' iter_final= ',iter_final
         end if
      end if

    !$acc exit data delete(&
    !$acc iter_filterl(:), &
    !$acc iter_mapl(:), &
    !$acc next_filterl(:), &
    !$acc next_mapl(:), &
    !$acc obu(:), &
    !$acc ur(:), &
    !$acc um(:), &
    !$acc ugust_total(:), &
    !$acc ustar(:), &
    !$acc temp1(:), &
    !$acc temp2(:), &
    !$acc temp12m(:), &
    !$acc temp22m(:), &
    !$acc fm(:), &
    !$acc canyon_resistance(:), &
    !$acc ramu(:), &
    !$acc fwet_roof(:), &
    !$acc fwet_road_imperv(:), &
    !$acc wind_speed_adj(:), &
    !$acc prev_tau(:), &
    !$acc tau_diff(:), &
    !$acc prev_tau_diff(:) )

    end associate

  end subroutine UrbanFluxes

  !-----------------------------------------------------------------------
  pure function surface_wetness(snow_depth, h2o_top, pondmx, qaf, qg) result(fwet)
    !$acc routine seq
    !
    ! !DESCRIPTION:
    ! Wet fraction of a roof or impervious road: from the snow, else the ponded water; dew makes it fully wet.
    !
    ! !ARGUMENTS:
    real(r8), value, intent(in) :: snow_depth ! snow height (m)
    real(r8), value, intent(in) :: h2o_top    ! liquid water + ice of the top layer (kg/m2)
    real(r8), value, intent(in) :: pondmx     ! ponding limit (mm)
    real(r8), value, intent(in) :: qaf        ! urban canopy air specific humidity (kg/kg)
    real(r8), value, intent(in) :: qg         ! specific humidity at the surface (kg/kg)
    real(r8) :: fwet                          ! wet fraction (-)
    !-----------------------------------------------------------------------
    if (snow_depth > 0._r8) then
       fwet = min(snow_depth/0.05_r8, 1._r8)
    else
       fwet = (max(0._r8, h2o_top)/pondmx)**0.666666666666_r8
       fwet = min(fwet,1._r8)
    end if
    if (qaf > qg) then
       fwet = 1._r8
    end if
  end function surface_wetness

end module UrbanFluxesMod
