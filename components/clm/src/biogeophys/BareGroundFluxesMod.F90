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
  use TopounitDataType     , only : top_as
  use LandunitType         , only : lun_pp
  use ColumnType           , only : col_pp
  use ColumnDataType       , only : col_es
  use ColumnDataType       , only : col_ef
  use ColumnDataType       , only : col_ws
  use VegetationType       , only : veg_pp
  use VegetationDataType   , only : veg_es
  use VegetationDataType   , only : veg_ef
  use VegetationDataType   , only : veg_wf
  use VegetationDataType   , only : veg_ws
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine BareGroundFluxes(bounds, num_nolakeurbanp, filter_nolakeurbanp, &
       atm2lnd_vars, canopystate_vars, soilstate_vars, &
       frictionvel_vars, ch4_vars &
       )
    ! !DESCRIPTION:
    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step.
    !$acc routine seq
    ! !USES:
    use shr_const_mod        , only : SHR_CONST_RGAS
    use clm_varcon           , only : cpair, vkc, grav, denice, denh2o
    use clm_varctl           , only : use_lch4
    use landunit_varcon      , only : istsoil, istcrop
    use FrictionVelocityMod  , only : FrictionVelocity, MoninObukIni
    use QSatMod              , only : QSat
    use SurfaceResistanceMod , only : do_soilevap_beta

    implicit NONE
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



    !
    ! !LOCAL VARIABLES:
    integer, parameter  :: niters = 3            ! maximum number of iterations for surface temperature
    integer  :: p,c,t,g,f,j,l,pp                 ! indices
    integer  :: filterp(bounds%endp-bounds%begp+1) ! patch filter for vegetated patches
    integer  :: fn                               ! number of values in local pft filter
    integer  :: fp                               ! lake filter pft index
    integer  :: iter                             ! iteration index
    real(r8) :: zldis(bounds%begp:bounds%endp)   ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(bounds%begp:bounds%endp)  ! displacement height [m]
    real(r8) :: zeta                             ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: beta                             ! coefficient of convective velocity [-]
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

      !---------------------------------------------------
      ! Filter patches where canopystate_vars%frac_veg_nosno_patch IS ZERO
      !---------------------------------------------------
      beta = 1._r8 ! previously set as a constant for all columns in CanopyTemperature()

      !---------------------------------------------------
      ! Filter patches where frac_veg_nosno IS ZERO
      !---------------------------------------------------
      fn = 0
      p  = 0
      fp = 0
      filterp(:) = 0
      fm(:) = 0

      do fp = 1,num_nolakeurbanp
         p = filter_nolakeurbanp(fp)

         if (canopystate_vars%frac_veg_nosno_patch(p) == 0) then
            fn = fn + 1
            filterp(fn) = p
         end if
      end do

      ! Compute sensible and latent fluxes and their derivatives with respect
      ! to ground temperature using ground temperatures from previous time step
      !$acc loop seq
      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! Initialization variables

         displa(p) = 0._r8
         veg_ef%dlrad(p)  = 0._r8
         veg_ef%ulrad(p)  = 0._r8

         ur(p)    = max(1.0_r8,sqrt(top_as%ubot(t)*top_as%ubot(t)+top_as%vbot(t)*top_as%vbot(t)))
         dth(p)   = veg_es%thm(p)-col_es%t_grnd(c)
         dqh(p)   = top_as%qbot(t) - col_ws%qg(c)
         dthv     = dth(p)*(1._r8+0.61_r8*top_as%qbot(t))+0.61_r8*top_as%thbot(t)*dqh(p)
         zldis(p) = frictionvel_vars%forc_hgt_u_patch(p)

         ! Copy column roughness to local pft-level arrays

         z0mg_patch(p) = frictionvel_vars%z0mg_col(c)
         z0hg_patch(p) = frictionvel_vars%z0hg_col(c)
         z0qg_patch(p) = frictionvel_vars%z0qg_col(c)

         ! Initialize Monin-Obukhov length and wind speed

         call MoninObukIni(ur(p), col_es%thv(c), dthv, zldis(p), z0mg_patch(p), um(p), obu(p))

      end do
      !! Perform stability iteration
      !! Determine friction velocity, and potential temperature and humidity
      !! profiles of the surface boundary layer
      iter = 0
      ITERATION : do while (iter <= niters .and. fn>0)

          call FrictionVelocity(bounds%begp, bounds%endp, fn, filterp, &
              displa, z0mg_patch,&
              z0hg_patch, z0qg_patch, &
              obu, iter+1, ur, &
              um, ustar, &
              temp1, temp2,&
              temp12m, temp22m, fm, &
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
            thvstar = tstar*(1._r8+0.61_r8*top_as%qbot(t)) + 0.61_r8*top_as%thbot(t)*qstar
            zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*col_es%thv(c))

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(p) = max(ur(p),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               wc = beta*(-grav*ustar(p)*thvstar*col_pp%zii(c)/col_es%thv(c))**0.333_r8
               um(p) = sqrt(ur(p)*ur(p) + wc*wc)
            end if
            obu(p) = zldis(p)/zeta
         end do
         iter = iter + 1

      end do ITERATION! end stability iteration

      !$acc loop seq
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
         raih = top_as%rhobot(t)*cpair/rah
         if (use_lch4) then
            ch4_vars%grnd_ch4_cond_patch(p) = 1._r8/raw
         end if

         ! Soil evaporation resistance
         www = (col_ws%h2osoi_liq(c,1)/denh2o+col_ws%h2osoi_ice(c,1)/denice)/col_pp%dz(c,1)/soilstate_vars%watsat_col(c,1)
         www = min(max(www,0.0_r8),1._r8)

         !changed by K.Sakaguchi. Soilbeta is used for evaporation
         if (dqh(p) > 0._r8) then  !dew  (beta is not applied, just like rsoil used to be)
            raiw = top_as%rhobot(t)/(raw)
         else
            if(do_soilevap_beta())then
               ! Lee and Pielke 1992 beta is applied
               raiw    = soilstate_vars%soilbeta_col(c)*top_as%rhobot(t)/(raw)
            endif
         end if

         frictionvel_vars%ram1_patch(p) = ram  !pass value to global variable

         ! Output to pft-level data structures
         ! Derivative of fluxes with respect to ground temperature
         veg_ef%cgrnds(p) = raih
         veg_ef%cgrndl(p) = raiw*col_ws%dqgdT(c)
         veg_ef%cgrnd(p)  = veg_ef%cgrnds(p) + col_ef%htvp(c)*veg_ef%cgrndl(p)

         ! Surface fluxes of momentum, sensible and latent heat
         ! using ground temperatures from previous time step
         veg_ef%taux(p)          = -top_as%rhobot(t)*top_as%ubot(t)/ram
         veg_ef%tauy(p)          = -top_as%rhobot(t)*top_as%vbot(t)/ram
         veg_ef%eflx_sh_grnd(p)  = -raih*dth(p)
         veg_ef%eflx_sh_tot(p)   = veg_ef%eflx_sh_grnd(p)

         veg_ef%eflx_sh_snow(p)   = -raih*(veg_es%thm(p)-col_es%t_soisno(c,col_pp%snl(c)+1))
         veg_ef%eflx_sh_soil(p)   = -raih*(veg_es%thm(p)-col_es%t_soisno(c,1))
         veg_ef%eflx_sh_h2osfc(p) = -raih*(veg_es%thm(p)-col_es%t_h2osfc(c))


         veg_wf%qflx_evap_soi(p)  = -raiw*dqh(p)
         veg_wf%qflx_evap_tot(p)  = veg_wf%qflx_evap_soi(p)

         veg_wf%qflx_ev_snow(p)   = -raiw*(top_as%qbot(t) - col_ws%qg_snow(c))
         veg_wf%qflx_ev_soil(p)   = -raiw*(top_as%qbot(t) - col_ws%qg_soil(c))
         veg_wf%qflx_ev_h2osfc(p) = -raiw*(top_as%qbot(t) - col_ws%qg_h2osfc(c))

         veg_es%t_ref2m(p) = veg_es%thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

         veg_ws%q_ref2m(p) = top_as%qbot(t) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

         ! 2 m height relative humidity
         call QSat(veg_es%t_ref2m(p), top_as%pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

         veg_ws%rh_ref2m(p) = min(100._r8, veg_ws%q_ref2m(p) / qsat_ref2m * 100._r8)

         if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
            veg_ws%rh_ref2m_r(p) = veg_ws%rh_ref2m(p)
            veg_es%t_ref2m_r(p) = veg_es%t_ref2m(p)
         end if

      end do

  end subroutine BareGroundFluxes

end module BareGroundFluxesMod
