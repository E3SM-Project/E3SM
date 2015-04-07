module lnd2atmMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_megan_mod        , only : shr_megan_mechcomps_n
  use clm_varpar           , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon           , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl           , only : iulog, use_c13, use_cn, use_lch4, use_voc
  use seq_drydep_mod       , only : n_drydep, drydep_method, DD_XLND
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : p2g, c2g 
  use lnd2atmType          , only : lnd2atm_type
  use atm2lndType          , only : atm2lnd_type
  use ch4Mod               , only : ch4_type
  use CNVegCarbonFluxType  , only : cnveg_carbonflux_type
  use DUSTMod              , only : dust_type
  use DryDepVelocity       , only : drydepvel_type
  use VocEmissionMod       , only : vocemis_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityMod  , only : frictionvel_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use TemperatureType      , only : temperature_type
  use WaterFluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use GridcellType         , only : grc                
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm
  public :: lnd2atm_minimal
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine lnd2atm_minimal(bounds, &
      waterstate_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)
    !
    ! !DESCRIPTION:
    ! Compute clm_l2a_inst component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !USES:
    use clm_varcon, only : sb
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    type(waterstate_type) , intent(in)    :: waterstate_inst
    type(surfalb_type)    , intent(in)    :: surfalb_inst
    type(energyflux_type) , intent(in)    :: energyflux_inst
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: g                                    ! index
    real(r8), parameter :: amC   = 12.0_r8          ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8          ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    call c2g(bounds, &
         waterstate_inst%h2osno_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%h2osno_grc    (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_inst%h2osno_grc(g) = lnd2atm_inst%h2osno_grc(g)/1000._r8
    end do

    call c2g(bounds, nlevgrnd, &
         waterstate_inst%h2osoi_vol_col (bounds%begc:bounds%endc, :), &
         lnd2atm_inst%h2osoi_vol_grc    (bounds%begg:bounds%endg, :), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_inst%albd_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_inst%albd_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_inst%albi_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_inst%albi_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_inst%eflx_lwrad_out_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%eflx_lwrad_out_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_inst%t_rad_grc(g) = sqrt(sqrt(lnd2atm_inst%eflx_lwrad_out_grc(g)/sb))
    end do

  end subroutine lnd2atm_minimal

  !------------------------------------------------------------------------
  subroutine lnd2atm(bounds, &
       atm2lnd_inst, surfalb_inst, temperature_inst, frictionvel_inst, &
       waterstate_inst, waterflux_inst, energyflux_inst,               &
       solarabs_inst, cnveg_carbonflux_inst, drydepvel_inst,  &
       vocemis_inst, dust_inst, ch4_inst, lnd2atm_inst) 
    !
    ! !DESCRIPTION:
    ! Compute lnd2atm_inst component of gridcell derived type
    !
    ! !USES:
    use ch4varcon  , only : ch4offline
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds  
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_inst
    type(surfalb_type)          , intent(in)    :: surfalb_inst
    type(temperature_type)      , intent(in)    :: temperature_inst
    type(frictionvel_type)      , intent(in)    :: frictionvel_inst
    type(waterstate_type)       , intent(inout) :: waterstate_inst
    type(waterflux_type)        , intent(in)    :: waterflux_inst
    type(energyflux_type)       , intent(in)    :: energyflux_inst
    type(solarabs_type)         , intent(in)    :: solarabs_inst
    type(cnveg_carbonflux_type) , intent(in)    :: cnveg_carbonflux_inst
    type(drydepvel_type)        , intent(in)    :: drydepvel_inst
    type(vocemis_type)          , intent(in)    :: vocemis_inst
    type(dust_type)             , intent(in)    :: dust_inst
    type(ch4_type)              , intent(in)    :: ch4_inst
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: g             ! index
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    !----------------------------------------------------
    ! lnd -> atm
    !----------------------------------------------------
    
    ! First, compute the "minimal" set of fields.
    call lnd2atm_minimal(bounds, &
         waterstate_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)

    call p2g(bounds, &
         temperature_inst%t_ref2m_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%t_ref2m_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         waterstate_inst%q_ref2m_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%q_ref2m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_inst%u10_clm_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%u_ref10m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_inst%taux_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%taux_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_inst%tauy_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%tauy_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         waterflux_inst%qflx_evap_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%qflx_evap_tot_grc     (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         solarabs_inst%fsa_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%fsa_grc    (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_inst%fv_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%fv_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_inst%ram1_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%ram1_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g( bounds, &
         energyflux_inst%eflx_sh_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%eflx_sh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity',c2l_scale_type='urbanf',l2g_scale_type='unity')
    do g = bounds%begg, bounds%endg
       lnd2atm_inst%eflx_sh_tot_grc(g) =  lnd2atm_inst%eflx_sh_tot_grc(g) - &
            energyflux_inst%eflx_dynbal_grc(g) 
    enddo

    call p2g(bounds, &
         energyflux_inst%eflx_lh_tot_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%eflx_lh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    if (use_cn) then
       call c2g(bounds, &
            cnveg_carbonflux_inst%nee_col(bounds%begc:bounds%endc), &
            lnd2atm_inst%nee_grc         (bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity')

       if (use_lch4) then
          if (.not. ch4offline) then
             ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
             do g = bounds%begg,bounds%endg
                ! nem is in g C/m2/sec
                lnd2atm_inst%nee_grc(g) = lnd2atm_inst%nee_grc(g) + lnd2atm_inst%nem_grc(g) 
             end do
          end if
       end if

       ! Convert from gC/m2/s to kgCO2/m2/s
       do g = bounds%begg,bounds%endg
          lnd2atm_inst%nee_grc(g) = lnd2atm_inst%nee_grc(g)*convertgC2kgCO2
       end do
    else
       do g = bounds%begg,bounds%endg
          lnd2atm_inst%nee_grc(g) = 0._r8
       end do
    end if

    ! drydepvel
    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
       call p2g(bounds, n_drydep, &
            drydepvel_inst%velocity_patch (bounds%begp:bounds%endp, :), &
            lnd2atm_inst%ddvel_grc        (bounds%begg:bounds%endg, :), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    endif

    ! voc emission flux
    if (use_voc .and. shr_megan_mechcomps_n>0) then
       call p2g(bounds, shr_megan_mechcomps_n, &
            vocemis_inst%vocflx_patch(bounds%begp:bounds%endp,:), &
            lnd2atm_inst%flxvoc_grc  (bounds%begg:bounds%endg,:), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    end if

    ! dust emission flux
    call p2g(bounds, ndst, &
         dust_inst%flx_mss_vrt_dst_patch(bounds%begp:bounds%endp, :), &
         lnd2atm_inst%flxdst_grc        (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')


    ! ch4 flux
    if (use_lch4) then
       call c2g( bounds,     &
            ch4_inst%ch4_surf_flux_tot_col (bounds%begc:bounds%endc), &
            lnd2atm_inst%flux_ch4_grc      (bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity' )
    end if

    !----------------------------------------------------
    ! lnd -> rof
    !----------------------------------------------------

    call c2g( bounds, &
         waterflux_inst%qflx_runoff_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%qflx_rofliq_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       lnd2atm_inst%qflx_rofliq_grc(g) = lnd2atm_inst%qflx_rofliq_grc(g) - waterflux_inst%qflx_liq_dynbal_grc(g)
    enddo

    call c2g( bounds, &
         waterflux_inst%qflx_snwcp_ice_col(bounds%begc:bounds%endc),  &
         lnd2atm_inst%qflx_rofice_grc     (bounds%begg:bounds%endg),  & 
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       lnd2atm_inst%qflx_rofice_grc(g) = lnd2atm_inst%qflx_rofice_grc(g) - waterflux_inst%qflx_ice_dynbal_grc(g)          
    enddo

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb
    ! second add river storage as gridcell average depth (1.e-3 converts [m3/km2] to [mm])
    ! TODO - this was in BalanceCheckMod - not sure where it belongs?

    call c2g( bounds, &
         waterstate_inst%endwb_col(bounds%begc:bounds%endc), &
         waterstate_inst%tws_grc  (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       waterstate_inst%tws_grc(g) = waterstate_inst%tws_grc(g) + atm2lnd_inst%volr_grc(g) / grc%area(g) * 1.e-3_r8
    enddo

  end subroutine lnd2atm

end module lnd2atmMod
