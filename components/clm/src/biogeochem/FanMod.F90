module FanMod
#ifdef _PYMOD_
  use qsatmod
#else
  use shr_const_mod
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use QSatMod                         , only : QSat
#endif
  implicit none

#ifdef _REALPR_
  integer, parameter :: r8 = 8
#endif
  
#ifdef _PYMOD_
  public

  real(r8), parameter :: SHR_CONST_BOLTZ   = 1.38065e-23
  real(r8), parameter :: SHR_CONST_AVOGAD  = 6.02214e26
  real(r8), parameter :: SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ
  real(r8), parameter :: SHR_CONST_MWDAIR  = 28.966
  real(r8), parameter :: SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR

#else
  private
  public update_org_n
  public eval_fluxes_storage
  public update_npool
  public update_3pool
  public update_4pool
  public update_urea
#endif
  
  ! Indices in flux arrays, soil:
  integer, parameter, public :: iflx_air = 1, & ! flux to air
       iflx_soild = 2, & ! diffusion to soil
       iflx_no3 = 3, &   ! nitrification
       iflx_soilq = 4, & ! percolation to soil
       iflx_roff = 5, &  ! surface runoff
       iflx_to_tan = 6   ! conversion to tan (from urea) 
  ! Number of different fluxes, the minimum size for flux vectors:
  integer, parameter, public :: num_fluxes = 6

  ! Indices in flux arrays, storage:
  integer, parameter, public :: iflx_air_barns = 1, &
       iflx_air_stores = 2, &
       iflx_appl = 3, &
       iflx_to_store = 4
  ! Indices in the organic N pools and fluxes
  integer, parameter, public :: ind_avail = 1, ind_resist = 2, ind_unavail = 3 
  
  ! nominal depth where the soil TAN concentration vanishes:
  real(r8), parameter, public :: soildepth_reservoir = 0.05_r8

  integer, parameter, public :: err_bad_theta = 1, err_negative_tan = 2, err_negative_flux = 3, &
       err_balance_tan = 4, err_balance_nitr = 5, err_nan = 6, err_bad_subst = 7

  integer, parameter, public :: subst_tan = 1, subst_urea = 2

  logical, parameter, public :: debug_fan = .false.
  
  real(r8), parameter, public :: water_relax_t = 24*3600.0_r8
  
contains

  ! accessor functions are only needed for the python interface... 
  !
  integer function ind_soild() result(ind)
    ind = iflx_soild
  end function ind_soild

  integer function ind_soilq() result(ind)
    ind = iflx_soilq
  end function ind_soilq

  integer function ind_air() result(ind)
    ind = iflx_air
  end function ind_air

  integer function ind_no3() result(ind)
    ind = iflx_no3
  end function ind_no3

  integer function ind_roff() result(ind)
    ind = iflx_roff
  end function ind_roff

  integer function ind_air_barns() result(ind)
    ind = iflx_air_barns
  end function ind_air_barns

  integer function ind_air_stores() result(ind)
    ind = iflx_air_stores
  end function ind_air_stores

  integer function ind_appl() result(ind)
    ind = iflx_appl
  end function ind_appl

  integer function ind_to_store() result(ind)
    ind = iflx_to_store
  end function ind_to_store
  
  function eval_diffusivity_liq_mq(theta, thetasat, tg) result(diff)
    ! Evaluate the aquous phase diffusivity for TAN in soil according to the Millington &
    ! Quirk model.
    implicit none
    real(r8), intent(in) :: theta, thetasat, tg
    real(r8) :: diff
    
    real(r8) :: kaq_base
    real(r8), parameter :: pw = 10.0_r8 / 3.0_r8
    
    kaq_base = 9.8e-10_r8 * 1.03_r8 ** (Tg-273.0_r8)
    diff = kaq_base * (theta**pw) / (thetasat**2)

  end function eval_diffusivity_liq_mq

  function eval_diffusivity_gas_mq(theta, thetasat, tg) result(diff)
    ! Evaluate the gas phase diffusivity for NH3 in soil according to the Millington &
    ! Quirk model.
    implicit none
    real(r8), intent(in) :: theta, thetasat, tg
    real(r8) :: diff
    
    real(r8) :: soilair, dair
    real(r8), parameter :: pw = 10.0_r8 / 3.0_r8
    real(r8), parameter :: mNH3 = 17., mair = 29, vNH3 = 14.9, vair = 20.1, press = 1.0
    
    soilair = thetasat - theta
    !dair = 1.7e-5_r8 * 1.03_r8**(Tg-293.0_r8)
    !dair = 18e-6_r8
    !dair = 1.4e-5
    
    ! Base rate from Perry's Chemical Engineer's Handbook, 8th ed.
    dair = (0.001 * tg**1.75 * sqrt(1/mNH3 + 1/mair)) / (press * (vair**(1./3) * vNH3**(1./3))**2) * 1e-4
    diff = dair * (soilair**pw) / (thetasat**2)

  end function eval_diffusivity_gas_mq

  function eval_diffusivity_gas_m03(theta, thetasat, tg) result(diff)
    ! Evaluate the gas phase diffusivity for NH3 in soil according to the method of
    ! Moldrup (2003).
    implicit none
    real(r8), intent(in) :: theta, thetasat, tg
    real(r8) :: diff
    
    real(r8) :: soilair, dair
    real(r8), parameter :: pw = 10.0_r8 / 3.0_r8
    real(r8), parameter :: bsw = 5.0_r8, m03_T = 2.0, m03_W = 3.0 / bsw
    
    soilair = thetasat - theta
    dair = 1.7e-5 * 1.03**(Tg-293.0_r8)
   
    diff = dair * soilair**m03_T * (soilair/thetasat)**m03_W
    
  end function eval_diffusivity_gas_m03

  function eval_diffusivity_liq_m03(theta, thetasat, tg) result(diff)
    ! Evaluate the aquous phase diffusivity for TAN in soil according to the method of Moldrup (2003).
    implicit none
    real(r8), intent(in) :: theta, thetasat, tg
    real(r8) :: diff
    
    real(r8) :: kaq_base
    real(r8), parameter :: pw = 10.0_r8 / 3.0_r8, bsw = 5.0_r8, m03_T = 2.0_r8, m03_W = 0.3333_r8*bsw - 1.0
    
    kaq_base = 9.8e-10 * 1.03 ** (Tg-273.0_r8)
    
    diff = kaq_base * theta**m03_T * (theta/thetasat)**m03_W

  end function eval_diffusivity_liq_m03
  

  ! The following three subroutines are meant for evaluating the net NH3 flux between soil and atmosphere as
  !
  ! F = g * (N - beta * cair)
  !
  ! where F is the net upwards flux of NH3, g is a conductance (m/s), cair is the
  ! atmospheric concentration of NH3, and beta is a unitless constant. N denotes the TAN
  ! concentration in soil in different phases as detailed below.
  
  subroutine get_volat_coefs_liq(ratm, tg, theta, thetasat, Hconc, depth, conductance, beta)
    ! 
    ! Evaluate the conductance g with N = [NH4+ (aq)] + [NH3 (aq)] in soilwater and
    ! assuming that [NH3 (g)] < N in soil.
    ! 
    ! For cair = 0, this subroutine is actually equivalent to get_volat_soil_leachn.
    ! 
    implicit none
    real(r8), intent(in) :: ratm ! resistance between the soil surface and bulk atmosphere
    real(r8), intent(in) :: tg   ! ground temperature
    real(r8), intent(in) :: theta ! volumetric water content 
    real(r8), intent(in) :: thetasat ! volumetric water content at saturation
    real(r8), intent(in) :: Hconc ! hydrogen ion concentration, -log10(pH)
    real(r8), intent(in) :: depth ! thickenss of the soil layer, m
    real(r8), intent(out) :: conductance ! as defined above  
    real(r8), intent(out) :: beta ! as defined above
    
    real(r8) :: dz, henry_eff, dsl, dsg, rsl, rsg, grad, cond, air
    
    dz = 0.5*depth 
    dsl = eval_diffusivity_liq_mq(theta, thetasat, Tg)
    dsg = eval_diffusivity_gas_mq(theta, thetasat, Tg)
    
    henry_eff = get_henry_eff(Tg, Hconc)
    beta = 1/henry_eff
    air = thetasat - theta

    rsg = dz / (dsg*air)
    rsl = dz / (dsl*theta)
    
    conductance = henry_eff*(henry_eff*rsl + rsg)/(henry_eff*ratm*rsl + henry_eff*rsg*rsl + ratm*rsg)
    
  end subroutine get_volat_coefs_liq

  subroutine get_volat_coefs_bulk(ratm, tg, theta, thetasat, Hconc, depth, conductance, beta)
    !
    ! Evaluate the conductance g with N = [NH4+ (aq)] + [NH3 (aq)] + [NH3 (g)] measured per volume of soil.
    !
    ! More accurate than get_volat_coefs_liq but in practice rarely different because when
    ! measured in mass per volume, most of the TAN is normally in soil water.
    
    implicit none
    real(r8), intent(in) :: ratm, tg, theta, thetasat, Hconc, depth
    real(r8), intent(out) :: conductance, beta
    
    real(r8) :: dz, henry_eff, dsl, dsg, rsl, rsg, cond, air
    
    dz = 0.5*depth 
    dsl = eval_diffusivity_liq_mq(theta, thetasat, Tg)
    dsg = eval_diffusivity_gas_mq(theta, thetasat, Tg)
    
    henry_eff = get_henry_eff(Tg, Hconc)
    beta = (air*henry_eff + theta)/henry_eff
    air = thetasat - theta

    rsg = dz / (dsg*air)
    rsl = dz / (dsl*theta)

    conductance = henry_eff*(henry_eff*rsl + rsg) &
         /(air*henry_eff**2*ratm*rsl + air*henry_eff**2*rsg*rsl + air*henry_eff*ratm*rsg &
         + henry_eff*ratm*rsl*theta + henry_eff*rsg*rsl*theta + ratm*rsg*theta)
    
  end subroutine get_volat_coefs_bulk

  subroutine get_volat_coefs_3p(ratm, tg, theta, thetasat, Hconc, depth, cond, beta)
    !
    ! Evaluate the conductance g including absorbed "solid" NH4+,
    !
    ! N = [NH3 (aq)] + [NHH4+ (aq)] + [NH3 (g)] + [NH4+ (s)].
    !
    ! The partitioning between absorbed and aquoeous NH4+ is evaluted with a linear isotherm
    ! such that [NH4+ (s)] / [NH4+ (aq)] = kxc, where kxc is a constant set below.
    !
    implicit none
    real(r8), intent(in) :: ratm, tg, theta, thetasat, Hconc, depth
    real(r8), intent(out) :: cond, beta
    
    real(r8) :: dz, henry_eff, dsl, dsg, rsl, rsg, air, solid
    real(r8), parameter :: kxc = 0.3_r8
    
    dz = 0.5*depth 
    dsl = eval_diffusivity_liq_mq(theta, thetasat, Tg)
    dsg = eval_diffusivity_gas_mq(theta, thetasat, Tg)
    
    henry_eff = get_henry_eff(Tg, Hconc)
    solid = 1 - thetasat
    air = thetasat - theta
    beta = (air*henry_eff + kxc*solid + theta)/henry_eff

    rsg = dz / (dsg*air)
    rsl = dz / (dsl*theta)

    cond = henry_eff*(henry_eff*rsl + rsg) &
         / (air*henry_eff**2*ratm*rsl + air*henry_eff**2*rsg*rsl + air*henry_eff*ratm*rsg &
         + henry_eff*kxc*solid*ratm*rsl + henry_eff*kxc*solid*rsg*rsl + henry_eff*ratm*rsl*theta &
         + henry_eff*rsg*rsl*theta + kxc*solid*ratm*rsg + ratm*rsg*theta)
    
  end subroutine get_volat_coefs_3p

  function get_volat_soil_leachn(ratm, tg, theta, thetasat, Hconc, depth) result(rate)
    ! Evaluate the instanteneous volatilization rate from soil as done in the LEACHN
    ! model. Includes gas and aquoues phase diffusion within soil and gas/liquid
    ! partitioning.
    real(r8), intent(in) :: ratm, tg, theta, thetasat, Hconc, depth
    real(r8) :: rate
    
    real(r8) :: dz, henry_eff, dsl, dsg, dstot, gs, gatm_eff
    
    dz = 0.5*depth 

    henry_eff = get_henry_eff(Tg, Hconc)
    gatm_eff = henry_eff / ratm 
    dsl = eval_diffusivity_liq_mq(theta, thetasat, Tg)
    dsg = eval_diffusivity_gas_mq(theta, thetasat, Tg)
    dstot = dsl*theta + dsg*henry_eff*(thetasat-theta)
    gs = dstot / dz
    rate = gs*gatm_eff / (gs + gatm_eff)
        
  end function get_volat_soil_leachn

  
  real(r8) function get_henry_eff(tg, Hconc) result(henry)
    ! Evaluate the "effective Henry constant" for ammonia, i.e the ratio
    ! H* = [NH3 (g)] / [NH3 (aq) + NH4+ (aq)]
    ! given a fixed H+ concentration in the solution.
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: Hconc ! H+ concentration, mol / l

    real(r8) :: KNH4, KH
    real(r8), parameter :: Tref = 298.15_r8
        
    KNH4 = 5.67_r8 * 1e-10_r8 * exp(-6286.0_r8 * (1.0_r8/Tg - 1.0_r8/Tref))
    KH = 4.59_r8 * Tg * exp(4092_r8 * (1.0_r8/Tg - 1.0_r8/Tref))
    henry = 1.0_r8 / (1.0_r8 + KH + KH*Hconc/KNH4)
    
  end function get_henry_eff

  real(r8) function eval_no3prod(theta, Tg, Hconc) result(kNO3)
    ! Evaluate nitrification rate as in the Riddick et al. (2016) paper.
    real(r8), intent(in) :: theta ! volumetric soil water m/m
    real(r8), intent(in) :: Tg    ! soil temperature, K
    real(r8), intent(in) :: Hconc ! hydrogen ion concentration mol/l

    real(r8) :: KNH4, KH, gas, stf, wmr, smrf, mNH4

    real(r8), parameter :: soil_dens = 1050.0_r8 ! Soil density, kg/m3
    real(r8), parameter :: water_dens = 1000.0_r8
    real(r8), parameter :: rmax = 1.16e-6_r8   ! Maximum rate of nitrification, s-1
    real(r8), parameter :: tmax = 313.0       ! Maximm temperature of microbial activity, K
    real(r8), parameter :: topt = 301.0       ! Optimal temperature of microbial acticity, K
    real(r8), parameter :: asg  = 2.4_8        ! a_sigma, empirical factor
    real(r8), parameter :: wmr_crit = 0.12_r8  ! Critical water content, g/g
    real(r8), parameter :: smrf_b = 2          ! Parameter in soil moisture response function
    real(r8), parameter :: Tref = 298.15_r8
    
    KNH4 = 5.67_r8 * 1e-10_r8 * exp(-6286.0_r8 * (1.0_r8/Tg - 1.0_r8/Tref))
    KH = 4.59_r8 * Tg * exp(4092_r8 * (1.0_r8/Tg - 1.0_r8/Tref))
    gas = 1.0_r8 / (1.0_r8 + KH + KH*Hconc/KNH4)
    mNH4 = gas * KH*Hconc/KNH4

    ! soil temperature function
    stf = (max(1e-3_r8, tmax-Tg) / (tmax-topt))**asg * exp(asg * (Tg-topt)/(tmax-topt))

    ! gravimetric soil water
    wmr = theta * water_dens / soil_dens

    ! soil moisture response function
    
    smrf = 1.0_r8 - exp(-(wmr/wmr_crit)**smrf_b)
    !if stf < 1e-9 or smrf < 1e-9:
    !    print theta
    !    1/0
    kNO3 = 2.0_r8 * rmax * mNH4 / (1.0_r8/stf + 1.0_r8/smrf)
    
  end function eval_no3prod
  
  subroutine eval_fluxes_slurry(water, mtan, Hconc, tg, ratm, theta, thetasat, perc, runoff, cnc_nh3_air, fluxes)
    ! Evaluate nitrogen fluxes for a partly infiltrated layer of slurry.
    ! The state of infiltration is detemined from the amounts water on surface and in soil.
    ! Positive flux means loss of TAN.
    implicit none
    real(r8), intent(in) :: water(2) ! water in (surface , subsurface), m
    real(r8), intent(in) :: mtan   ! TAN, mass units / m2, surface + subsurface
    real(r8), intent(out) :: fluxes(5) ! TAN fluxes, see top of the module
    real(r8), intent(in) :: Hconc    ! H+ concentration, -log10(pH)
    real(r8), intent(in) :: tg       ! soil temperature, K
    real(r8), intent(in) :: ratm     ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta    ! volumetric soil water in "clean" soil
    real(r8), intent(in) :: thetasat ! volumetric soil water at saturation
    real(r8), intent(in) :: perc     ! percolation water flux thourgh the bottom of volatilization layer, m/s
    real(r8), intent(in) :: runoff   ! surface runoff, m/s
    real(r8), intent(in) :: cnc_nh3_air ! atmospheric NH3 concentration, mass units / m3
    !real(r8), intent(in) :: dt       ! timestep

    real(r8) :: water_tot, cnc, air, depth_soilsat, diffusivity_water, diffusivity_satsoil, halfwater, insoil, r1, dz2
    real(r8) :: r2, volat_rate, kno3, henry_eff, depth_lower

    water_tot = water(1) + water(2)

    air = thetasat - theta
    ! depth of the saturated soil layer below the surface pool
    depth_soilsat = water(2) / air 

    cnc = mtan / water_tot

    fluxes(iflx_roff) = cnc * runoff
    fluxes(iflx_soilq) = cnc * perc


    diffusivity_water = 9.8e-10_r8 * 1.03_r8 ** (tg - 273.0_r8)
    diffusivity_satsoil = eval_diffusivity_liq_mq(thetasat, thetasat, tg) * thetasat

    halfwater = 0.5_r8 * water_tot

    ! Calculate the internal resistance r1 of the slurry/soil layer by integrating the
    ! diffusivity for distance that covers half of the slurry water.
    if (water(1) < halfwater) then
       ! contribution from both pool and the saturated soil.
       insoil = (halfwater - water(1)) / thetasat
       r1 = water(1) / diffusivity_water + insoil / diffusivity_satsoil
    else
       ! pool only
       r1 = halfwater / diffusivity_water
    end if
    
    depth_lower = max(soildepth_reservoir, depth_soilsat*1.5)
    ! Diffusion to deeper soil over distance dz2
    dz2 = depth_lower - 0.5*depth_soilsat
    r2 = 0.5 * depth_soilsat/diffusivity_satsoil + dz2 / eval_diffusivity_liq_mq(theta, thetasat, tg)
    !print *, 'r2', r2, diffusivity_satsoil, dz2, depth_soilsat
    fluxes(iflx_soild) = cnc / r2

    henry_eff = get_henry_eff(tg, Hconc)
    volat_rate = 1.0_r8 / (r1 + ratm / henry_eff) ! conductance from aqueous TAN in slurry to NH3 in atmosphere
    fluxes(iflx_air) = max(volat_rate*(cnc - cnc_nh3_air), 0.0_r8)
    
    ! nitrification
    kno3 = eval_no3prod(thetasat, tg, Hconc)
    fluxes(iflx_no3) = kno3 * mtan

    !fluxes(3:) = 0
    
  end subroutine eval_fluxes_slurry

  subroutine eval_fluxes_soil(mtan, water_manure, Hconc, tg, ratm, theta, thetasat, perc, &
       & runoff, cnc_nh3_air, soildepth, fluxes, substance, status)
    !
    ! Evaluate nitrogen fluxes from a soil layer. Use for all cases except the partly
    ! infiltrated slurry (above). Fluxes can be evaluated either for urea or TAN: for
    ! urea, only the aqueous phase fluxes are evaluated and nitrification is set to zero.
    ! 
    implicit none
    real(r8), intent(in) :: mtan ! TAN, mass units / m2
    real(r8), intent(in) :: water_manure ! water in the soil pool *in addition to* background soil water
    real(r8), intent(out) :: fluxes(5)   ! nitrogen fluxes, mass units / m2 / s, see top of module
    real(r8), intent(in) :: Hconc        ! Hydrogen ion concentration, mol/l
    real(r8), intent(in) :: tg           ! soil temperature, K
    real(r8), intent(in) :: ratm         ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta        ! volumetric soil water in "clean" soil, m/m
    real(r8), intent(in) :: thetasat     ! volumetric soil water at saturation
    real(r8), intent(in) :: perc         ! downwards water percolation rate at the bottom of layer, m/s, > 0
    real(r8), intent(in) :: runoff       ! runoff water flux, m / s
    real(r8), intent(in) :: cnc_nh3_air  ! NH3 concentration in air, mass units / m3
    real(r8), intent(in) :: soildepth    ! thickness of the volatlization layer
    integer, intent(in) :: substance     ! subst_tan or subst_urea.
    integer, intent(out) :: status      ! error flag
    
    
    real(r8) :: water_tot, cnc, air, henry_eff, dsl, dsg, dstot, dz2, no3_rate, volat_rate, theta_tot, beta

    water_tot = water_manure + theta*soildepth
    if (water_tot < 1e-9) then
       fluxes = 0.0
       return
    end if
    
    theta_tot = water_tot / soildepth
    if (theta_tot > thetasat) then
       status = err_bad_theta
       return
    end if

    cnc = mtan / water_tot

    air = thetasat - theta_tot
    beta = 0.0
    
    if (substance == subst_tan) then

       volat_rate = get_volat_soil_leachn(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth)
       !fluxes(iflx_air) = max((cnc-cnc_nh3_air) * volat_rate, 0.0_r8)

       !call get_volat_coefs_liq(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth, volat_rate, beta)
       !fluxes(iflx_air) = volat_rate * (cnc - beta*cnc_nh3_air)

       call get_volat_coefs_bulk(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth, volat_rate, beta)
       !call get_volat_coefs_3p(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth, volat_rate, beta)
       fluxes(iflx_air) = volat_rate * (mtan/soildepth - beta*cnc_nh3_air)

       henry_eff = get_henry_eff(tg, Hconc)
       dsg = eval_diffusivity_gas_mq(theta_tot, thetasat, tg)
       no3_rate = eval_no3prod(theta_tot, tg, Hconc)
    else if (substance == subst_urea) then
       fluxes(iflx_air) = 0.0_r8
       henry_eff = 0.0_r8
       dsg = 0.0_r8
       no3_rate = 0.0_r8
    else
       status = err_bad_subst
       return
    end if
    
    ! Downwards diffusion
    ! soil diffusivities: liquid, gas, bulk
    dsl = eval_diffusivity_liq_mq(theta_tot, thetasat, tg)
    dstot = dsl*theta_tot + dsg*henry_eff*air
    dz2 = soildepth_reservoir - 0.5 * soildepth
    !print *, 'dz2:', dz2
    fluxes(iflx_soild) = cnc * dstot / dz2
    fluxes(iflx_no3) = mtan * no3_rate
    fluxes(iflx_soilq) = cnc * perc
    fluxes(iflx_roff) = cnc * runoff

    status = 0
    !fluxes(4:) = 0
    
  end subroutine eval_fluxes_soil

  subroutine eval_fluxes_soilroff(mtan, water_manure, Hconc, tg, ratm, theta, thetasat, perc, &
       & runoff, cnc_nh3_air, soildepth, fluxes, substance, status)
    !
    ! Evaluate nitrogen fluxes from a soil layer. Use for all cases except the partly
    ! infiltrated slurry (above). Fluxes can be evaluated either for urea or TAN: for
    ! urea, only the aqueous phase fluxes are evaluated and nitrification is set to zero.
    ! 
    implicit none
    real(r8), intent(in) :: mtan ! TAN, mass units / m2
    real(r8), intent(in) :: water_manure ! water in the soil pool *in addition to* background soil water
    real(r8), intent(out) :: fluxes(5)   ! nitrogen fluxes, mass units / m2 / s, see top of module
    real(r8), intent(in) :: Hconc        ! Hydrogen ion concentration, mol/l
    real(r8), intent(in) :: tg           ! soil temperature, K
    real(r8), intent(in) :: ratm         ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta        ! volumetric soil water in "clean" soil, m/m
    real(r8), intent(in) :: thetasat     ! volumetric soil water at saturation
    real(r8), intent(in) :: perc         ! downwards water percolation rate at the bottom of layer, m/s, > 0
    real(r8), intent(in) :: runoff       ! runoff water flux, m / s
    real(r8), intent(in) :: cnc_nh3_air  ! NH3 concentration in air, mass units / m3
    real(r8), intent(in) :: soildepth    ! thickness of the volatlization layer
    integer, intent(in) :: substance     ! subst_tan or subst_urea.
    integer, intent(out) :: status      ! error flag
    
    
    real(r8) :: water_tot, cnc, air, henry_eff, dsl, dsg, dstot, dz2, no3_rate, volat_rate, theta_tot, beta
    real(r8) :: cnc_srfg, cnc_srfaq, cnc_soilaq, cnc_soilg, dz, rsl, rsg, qrf

    water_tot = water_manure + theta*soildepth
    if (water_tot < 1e-9) then
       fluxes = 0.0
       return
    end if
    
    theta_tot = water_tot / soildepth
    if (theta_tot > thetasat) then
       status = err_bad_theta
       return
    end if

    cnc = mtan / soildepth
    qrf = runoff
    air = thetasat - theta_tot
    beta = 0.0

    dz = 0.5*soildepth 
    dsl = eval_diffusivity_liq_mq(theta_tot, thetasat, Tg)
    dsg = eval_diffusivity_gas_mq(theta_tot, thetasat, Tg)
    rsg = dz / (dsg*air)
    rsl = dz / (dsl*theta_tot)

    if (substance == subst_tan) then
       henry_eff = get_henry_eff(tg, Hconc)
       no3_rate = eval_no3prod(theta_tot, tg, Hconc)
       cnc_soilg = cnc / (theta_tot/henry_eff + air)
    else if (substance == subst_urea) then
       henry_eff = 0.0
       no3_rate = 0.0_r8
       cnc_soilg = 0.0_r8
    else
       status = err_bad_subst
       return
    end if

    cnc_srfg = (henry_eff * ratm * cnc * (henry_eff*rsl + rsg)) &
         / (air*henry_eff**2*rsl*(ratm + rsg) + air*henry_eff*ratm*rsg*(qrf*rsl + 1) &
         + henry_eff*rsl*theta_tot*(ratm + rsg) + ratm*rsg*theta_tot*(qrf*rsl + 1))
    
    cnc_srfaq = ratm * cnc * (henry_eff*rsl + rsg)&
         / (air*henry_eff**2*rsl*(ratm + rsg) + air*henry_eff*ratm*rsg*(qrf*rsl + 1) &
         + henry_eff*rsl*theta_tot*(ratm + rsg) + ratm*rsg*theta_tot*(qrf*rsl + 1))
    
    fluxes(iflx_air) = cnc_srfg / ratm
    fluxes(iflx_roff) = runoff * cnc_srfaq

    dstot = dsl*theta_tot + dsg*henry_eff*air
    dz2 = soildepth_reservoir - 0.5 * soildepth
    
    cnc_soilaq = cnc / (theta_tot + air*henry_eff)
    
    fluxes(iflx_soild) = (cnc_soilaq * dsl) / dz2 + (cnc_soilg * dsg) / dz2


    fluxes(iflx_no3) = mtan * no3_rate
    fluxes(iflx_soilq) = cnc_soilaq * perc

    status = 0
    !fluxes(4:) = 0
    
  end subroutine eval_fluxes_soilroff


  subroutine eval_fluxes_soil_3p(mtan, water_manure, Hconc, tg, ratm, theta, thetasat, perc, &
       & runoff, cnc_nh3_air, soildepth, fluxes, substance, status)
    !
    ! Evaluate nitrogen fluxes from a soil layer. Use for all cases except the partly
    ! infiltrated slurry (above). Fluxes can be evaluated either for urea or TAN: for
    ! urea, only the aqueous phase fluxes are evaluated and nitrification is set to zero.
    ! 
    implicit none
    real(r8), intent(in) :: mtan ! TAN, mass units / m2
    real(r8), intent(in) :: water_manure ! water in the soil pool *in addition to* background soil water
    real(r8), intent(out) :: fluxes(5)   ! nitrogen fluxes, mass units / m2 / s, see top of module
    real(r8), intent(in) :: Hconc        ! Hydrogen ion concentration, mol/l
    real(r8), intent(in) :: tg           ! soil temperature, K
    real(r8), intent(in) :: ratm         ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta        ! volumetric soil water in "clean" soil, m/m
    real(r8), intent(in) :: thetasat     ! volumetric soil water at saturation
    real(r8), intent(in) :: perc         ! downwards water percolation rate at the bottom of layer, m/s, > 0
    real(r8), intent(in) :: runoff       ! runoff water flux, m / s
    real(r8), intent(in) :: cnc_nh3_air  ! NH3 concentration in air, mass units / m3
    real(r8), intent(in) :: soildepth    ! thickness of the volatlization layer
    integer, intent(in) :: substance     ! subst_tan or subst_urea.
    integer, intent(out) :: status      ! error flag
    
    
    real(r8) :: water_tot, cnc, air, henry_eff, dsl, dsg, dstot, dz2, no3_rate, volat_rate, theta_tot, beta
    real(r8) :: part_liq, part_gas, part_nh4, part_solid, part_avail
    
    water_tot = water_manure + theta*soildepth
    if (water_tot < 1e-9) then
       fluxes = 0.0
       return
    end if
    
    theta_tot = water_tot / soildepth
    if (theta_tot > thetasat) then
       status = err_bad_theta
       return
    end if

    cnc = mtan / water_tot

    air = thetasat - theta_tot
    beta = 0.0
    dz2 = soildepth_reservoir - 0.5 * soildepth
    
    if (substance == subst_tan) then
       call get_fluxes_3p(mtan/soildepth, ratm, theta_tot, thetasat, tg, 0.5*soildepth, dz2, Hconc, &
            fluxes(iflx_air), fluxes(iflx_soild), part_liq, part_gas, part_nh4, part_solid)
       !print *, 'part:', part_liq, part_gas, part_nh4, part_solid
       !volat_rate = get_volat_soil_leachn(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth)
       !fluxes(iflx_air) = max((cnc-cnc_nh3_air) * volat_rate, 0.0_r8)

       !call get_volat_coefs_liq(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth, volat_rate, beta)
       !fluxes(iflx_air) = volat_rate * (cnc - beta*cnc_nh3_air)

       !call get_volat_coefs_bulk(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth, volat_rate, beta)
       !call get_volat_coefs_3p(ratm, tg, water_tot/soildepth, thetasat, Hconc, soildepth, volat_rate, beta)
       !fluxes(iflx_air) = volat_rate * (mtan/soildepth - beta*cnc_nh3_air)
       part_avail = 1 - part_solid
       no3_rate = thetasat * part_avail * eval_no3prod(theta_tot, tg, Hconc) 
       no3_rate = eval_no3prod(theta_tot, tg, Hconc) 
       fluxes(iflx_no3) = no3_rate * mtan! * ( 1 - part_solid * (1-thetasat))
       cnc = mtan / soildepth * part_liq
    else if (substance == subst_urea) then
       fluxes(iflx_air) = 0.0_r8
       henry_eff = 0.0_r8
       dsg = 0.0_r8
       no3_rate = 0.0_r8
       part_avail = 1.0_r8
       dsl = eval_diffusivity_liq_mq(theta_tot, thetasat, tg)
       dstot = dsl*theta_tot + dsg*henry_eff*air
       fluxes(iflx_soild) = cnc * dstot / dz2
       fluxes(iflx_no3) = 0.0
    else
       status = err_bad_subst
       return
    end if
    
    ! Downwards diffusion
    ! soil diffusivities: liquid, gas, bulk
    
    !print *, 'dz2:', dz2
    fluxes(iflx_soilq) = cnc * perc 
    fluxes(iflx_roff) = cnc * runoff

    status = 0
    !fluxes(4:) = 0

  contains
    
    subroutine get_fluxes_3p(tan_soil, ratm, theta, thetasat, tg, dzup, dzdown, Hconc, &
         flux_atm, flux_soild, part_soil_liq, part_soil_gas, part_soil_nh4, part_soil_solid)
      real(r8), intent(in) :: tan_soil, ratm, theta, thetasat, tg, dzup, dzdown, Hconc
      real(r8), intent(out) :: flux_atm, flux_soild, part_soil_liq, part_soil_gas, part_soil_nh4, part_soil_solid
      
      real(r8) :: rsl, rsg, knh4, kh, air, solid, comp
      real(r8), parameter :: Tref = 298.15_r8
      real(r8), parameter :: kx = 1.0
      
      KNH4 = 5.67_r8 * 1e-10_r8 * exp(-6286.0_r8 * (1.0_r8/Tg - 1.0_r8/Tref))
      KH = 4.59_r8 * Tg * exp(4092_r8 * (1.0_r8/Tg - 1.0_r8/Tref))

      air = thetasat - theta
      solid = 1 - thetasat
      
      rsl = dzup / (theta * eval_diffusivity_liq_mq(theta, thetasat, tg))
      rsg = dzup / (air * eval_diffusivity_gas_mq(theta, thetasat, tg))

      comp = knh4*(cnc_nh3_air*rsg*rsl*(Hconc*kh*kx*solid + Hconc*kh*theta + &
           air*knh4 + kh*knh4*theta) + ratm*tan_soil*(Hconc*kh*rsg + kh*knh4*rsg + &
           knh4*rsl))/(Hconc**2*kh**2*kx*ratm*rsg*solid + Hconc**2*kh**2*ratm*rsg*theta + &
           Hconc*air*kh*knh4*ratm*rsg + Hconc*kh**2*knh4*kx*ratm*rsg*solid + &
           2*Hconc*kh**2*knh4*ratm*rsg*theta + Hconc*kh*knh4*kx*rsl*solid*(ratm + rsg) + &
           Hconc*kh*knh4*rsl*theta*(ratm + rsg) + air*kh*knh4**2*ratm*rsg + &
           air*knh4**2*rsl*(ratm + rsg) + kh**2*knh4**2*ratm*rsg*theta + &
           kh*knh4**2*rsl*theta*(ratm + rsg))

      flux_atm = comp / ratm

      part_soil_liq = kh*(Hconc + knh4)/(Hconc*kh*kx*solid + Hconc*kh*theta + air*knh4 + kh*knh4*theta)
      part_soil_nh4 = Hconc*kh/(Hconc*kh*kx*solid + Hconc*kh*theta + air*knh4 + kh*knh4*theta)
      part_soil_gas = knh4/(Hconc*kh*kx*solid + Hconc*kh*theta + air*knh4 + kh*knh4*theta)
      part_soil_solid = Hconc*kh*kx/(Hconc*kh*kx*solid + Hconc*kh*theta + air*knh4 + kh*knh4*theta)
      !flux_soild = &
      !     tan_soil*(kh*rsg*(Hconc + knh4) + knh4*rsl)&
      !     /(rsg*rsl*(Hconc*kh*kx*solid + H*kh*theta + epsilon*knh4 + kh*knh4*theta))

      rsl = dzdown / (theta * eval_diffusivity_liq_mq(theta, thetasat, tg))
      rsg = dzdown / (air * eval_diffusivity_gas_mq(theta, thetasat, tg))

      flux_soild = tan_soil * part_soil_liq / rsl + tan_soil * part_soil_gas / rsg
      
    end subroutine get_fluxes_3p
    
  end subroutine eval_fluxes_soil_3p

  subroutine partition_to_layer(water, theta, thetasat, soildepth, fraction_in, fraction_down, fraction_runoff)
    ! Evaluate the fraction of water volume that can be accommodated (before saturation)
    ! by soil layer with current water content theta.
    implicit none
    real(r8), intent(in) :: water ! water to be added to the layer, m
    real(r8), intent(in) :: theta, thetasat ! vol. soil water, current and saturation, m/m
    real(r8), intent(in) :: soildepth ! thickness of the layer, m
    real(r8), intent(out) :: fraction_in ! fraction that fits in
    real(r8), intent(out) :: fraction_down ! fraction excess
    real(r8), intent(out) :: fraction_runoff ! = 1 if theta is > 99% saturation otherwise 0

    real(r8) :: watvol, watvol_sat, vol_to_layer, vol_to_deeper, vol_avail
    
    watvol = theta*soildepth
    watvol_sat = thetasat*soildepth

    if (watvol < 0.99*watvol_sat) then
       vol_avail = watvol_sat - watvol
       vol_to_layer = min(water, vol_avail)
       vol_to_deeper = water - vol_to_layer
       fraction_down = vol_to_deeper / water
       fraction_runoff = 0.0
       fraction_in = 1.0_r8 - fraction_down
    else
       fraction_down = 0.0
       fraction_in = 0.0
       fraction_runoff = 1.0_r8
    end if
    
  end subroutine partition_to_layer

  subroutine age_pools_soil(ndep, dt, pools, mtan, garbage)
    implicit none
    real(r8), intent(in) :: ndep ! flux of TAN input, gN/m2/s
    real(r8), intent(inout) :: mtan(:) ! TAN pools for each age range. gN/m2
    real(r8), intent(in) :: dt ! timestep, s
    real(r8), intent(in) :: pools(:) ! age spans covered by each bin, seconds. size(mtan) = size(pools)
    real(r8), intent(out) :: garbage ! TAN removed from the oldest pool. gN/m2

    real(r8) :: flux_out(size(mtan))

    flux_out = mtan / pools

    ! new nitrogen
    mtan(1) = mtan(1) + ndep * dt
    ! transfer nitrogen from fresh to old pools
    mtan = mtan - flux_out * dt
    if (size(mtan) > 1) mtan(2:) = mtan(2:) + flux_out(:size(mtan)-1) * dt
    ! provided that the oldest pool has wide enough age range, the amount transferred out
    ! should be small.
    garbage = flux_out(size(mtan)) * dt
    
  end subroutine age_pools_soil

  subroutine age_pools_slurry(ndep, dt, water_slurry, tan_slurry, tan_soil, pools, garbage)
    implicit none
    real(r8), intent(in) :: ndep ! flux of TAN input, gN/m2/s
    real(r8), intent(in) :: dt ! timestep, s
    ! water in slurry pool, on surface (1) and below surface (2) not including the background water content (theta)
    real(r8), intent(in) :: water_slurry(2)
    real(r8), intent(inout) :: tan_slurry  ! TAN in slurry pool, gN/m2
    real(r8), intent(inout) :: tan_soil(:) ! TAN in soil pools, gN/m2
    real(r8), intent(in) :: pools(:) ! age spans covered by each pool, including the S0 surface slurry. seconds.
    real(r8), intent(out) :: garbage       ! garbage TAN, see above

    real(r8) :: fract_slurry_soil, flux_to_soilpools
    
    fract_slurry_soil = water_slurry(2) / sum(water_slurry)
    flux_to_soilpools = tan_slurry * fract_slurry_soil / pools(1)
    tan_slurry = tan_slurry + (ndep - flux_to_soilpools) * dt

    call age_pools_soil(flux_to_soilpools, dt, pools(2:), tan_soil, garbage)
    
  end subroutine age_pools_slurry

  subroutine update_3pool(tg, ratm, theta, thetasat, precip, evap, qbot, watertend, runoff, tandep, tanprod, cnc_nh3_air, &
       depth_slurry, poolranges, tanpools, fluxes, garbage, dt, status)
    !
    ! Evalute fluxes and update TAN pools for the 3-pool slurry model with a partly
    ! infiltrated, freshly infiltrated and aged TAN pools.
    !
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: ratm ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta ! volumetric soil water in soil column (unaffected by slurry)
    real(r8), intent(in) :: thetasat ! vol. soil water at saturation
    real(r8), intent(in) :: precip   ! precipitation, m/s
    real(r8), intent(in) :: evap   ! ground evaporation, m/s
    real(r8), intent(in) :: qbot   ! specific humidity (kg/kg) at lowest atmospheric model level
    real(r8), intent(in) :: watertend ! time derivative of theta
    real(r8), intent(in) :: runoff    ! surface runoff flux, m/s
    real(r8), intent(in) :: tandep    ! TAN input flux, gN/m2/s
    real(r8), intent(in) :: tanprod   ! TAN produced in the column, added to aged TAN pool
    !real(r8), intent(in) :: infiltr_slurry ! Slurry infiltration rate, m/s
    real(r8), intent(in) :: depth_slurry ! Initial slurry depth, m
    real(r8), intent(in) :: cnc_nh3_air ! NH3 concentration in air, gN/m3
    real(r8), intent(in) :: poolranges(3)  ! age ranges of TAN pools S0, S1, S2, sec. Slurry infiltration time is inferred from S0. 
    real(r8), intent(inout) :: tanpools(3) ! TAN pools gN/m2
    real(r8), intent(out) :: fluxes(5,3) ! TAN fluxes, gN/m2/s (type of flux, pool)
    real(r8), intent(out) :: garbage     ! over-aged TAN occurring during the step, gN/m.
    real(r8), intent(in) :: dt ! timestep, sec, >0
    integer, intent(out) :: status ! return status, 0 = good
    
    real(r8) :: infiltr_slurry, infiltrated, percolated, evap_slurry, water_slurry(2), perc_slurry_mean, waterloss
    real(r8) :: percolation, water_soil, age_prev, water_in_layer, tanpools_old(3)
    integer :: indpl
    
    real(r8), parameter :: dz_layer = 0.02 ! thickness of the volatilization layer, m
    ! H+ concentration in each pool 
    real(r8), parameter :: Hconc(3) = (/10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8)/)

    if (theta > thetasat) then
       status = err_bad_theta
       return
    end if

    tanpools_old = tanpools

    ! Pool S0
    !
    evap_slurry = get_evap_pool(tg, ratm, qbot)
    infiltr_slurry = max(depth_slurry / poolranges(1), precip)
    infiltrated = depth_slurry * infiltr_slurry / (infiltr_slurry + evap_slurry)
    ! Slurry water (in addition to soil water, theta) on surface and in soil. Represents
    ! mean over pool S0.
    water_slurry = (/0.5*depth_slurry, 0.5*infiltrated/) 
    ! The excess water assumed to have percolated down from the volat. layer.
    percolated = max(infiltrated - dz_layer*(thetasat-theta), 0.0)
    ! Percolation rate out of volat layer, average over the pool S0.
    perc_slurry_mean = percolated / poolranges(1)

    call eval_fluxes_slurry(water_slurry, tanpools(1), Hconc(1), tg, ratm, theta, thetasat, perc_slurry_mean, &
         runoff, cnc_nh3_air, fluxes(:,1))

    if (any(isnan(fluxes))) then
       status = err_nan * 10
    end if

    call update_pools(tanpools(1:1), fluxes(1:5,1:1), dt, 1, 5)

    if (any(tanpools < -1e-15)) then
       status = err_negative_tan
       return
    end if
    
    ! Pool aging & input
    !
    call age_pools_slurry(tandep, dt, water_slurry, tanpools(1), tanpools(2:3), poolranges, garbage)
    ! TAN produced (mineralization) goes to directly the old TAN pool. 
    tanpools(3) = tanpools(3) + tanprod*dt

    ! Soil bins S1 and S2
    !
    age_prev = 0 ! for water evaluations, consider beginning of S1 as the starting point
    water_in_layer = infiltrated - percolated ! water in layer just after slurry has infiltrated
    do indpl = 2, 3
       ! water content lost during the aging
       waterloss = water_in_layer * (waterfunction(age_prev) - waterfunction(age_prev+poolranges(indpl)))
       percolation = eval_perc(waterloss, evap, precip, watertend, poolranges(indpl))
       ! water content at the mean age of the pool
       water_soil = water_in_layer * waterfunction(age_prev + 0.5*poolranges(indpl))
       call eval_fluxes_soil(tanpools(indpl), water_soil, Hconc(indpl), tg, &
            & ratm, theta, thetasat, percolation, runoff, cnc_nh3_air, &
            & dz_layer, fluxes(:,indpl), subst_tan, status)
       if (status /= 0) return
       age_prev = age_prev + poolranges(indpl)
    end do
    
    call update_pools(tanpools(2:3), fluxes(1:5,2:3), dt, 2, 5)
    
    if (any(tanpools < -1e-15)) then
       status = err_negative_tan * 10
       return
       !end if
    end if

    if (any(isnan(fluxes))) then
       status = err_nan * 100
    end if

    if (abs(sum(tanpools - tanpools_old) - (-sum(fluxes) + tandep + tanprod)*dt + garbage) &
         > max(sum(tanpools_old)*1e-2, 1e-4)) then
       status = err_balance_tan
       return
    end if
    
    status = 0
    
  end subroutine update_3pool

  subroutine update_4pool(tg, ratm, theta, thetasat, precip, evap, qbot, watertend, runoff, tandep, tanprod, cnc_nh3_air, &
       depth_slurry, poolranges, tanpools, Hconc, fluxes, garbage, dt, status)
    !
    ! Experimental, as above but with an additional long-lived TAN pool.
    !
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: ratm ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta ! volumetric soil water in soil column (unaffected by slurry)
    real(r8), intent(in) :: thetasat ! vol. soil water at saturation
    real(r8), intent(in) :: precip   ! precipitation, m/s
    real(r8), intent(in) :: evap   ! ground evaporation, m/s
    real(r8), intent(in) :: qbot   ! specific humidity (kg/kg) at lowest atmospheric model level
    real(r8), intent(in) :: watertend ! time derivative of theta
    real(r8), intent(in) :: runoff    ! surface runoff flux, m/s
    real(r8), intent(in) :: tandep    ! TAN input flux, gN/m2/s
    real(r8), intent(in) :: tanprod   ! TAN produced in the column, added to aged TAN pool
    !real(r8), intent(in) :: infiltr_slurry ! Slurry infiltration rate, m/s
    real(r8), intent(in) :: depth_slurry ! Initial slurry depth, m
    real(r8), intent(in) :: cnc_nh3_air ! NH3 concentration in air, gN/m3
    real(r8), intent(in) :: poolranges(4)  ! age ranges of TAN pools S0, S1, S2, sec. Slurry infiltration time is inferred from S0. 
    real(r8), intent(inout) :: tanpools(4) ! TAN pools gN/m2
    real(r8), intent(out) :: fluxes(5,4) ! TAN fluxes, gN/m2/s (type of flux, pool)
    real(r8), intent(in) :: Hconc(4) ! H+ concentration
    real(r8), intent(out) :: garbage     ! over-aged TAN occurring during the step, gN/m.
    real(r8), intent(in) :: dt ! timestep, sec, >0
    integer, intent(out) :: status ! return status, 0 = good
    
    real(r8) :: infiltr_slurry, infiltrated, percolated, evap_slurry, water_slurry(2), perc_slurry_mean, waterloss
    real(r8) :: percolation, water_soil, age_prev, water_in_layer, tanpools_old(4)
    integer :: indpl
    
    real(r8), parameter :: dz_layer = 0.02 ! thickness of the volatilization layer, m
    ! H+ concentration in each pool 
    !real(r8), parameter :: Hconc(4) = (/10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-7_r8)/)

    if (theta > thetasat) then
       status = err_bad_theta
       return
    end if

    tanpools_old = tanpools

    ! Pool S0
    !
    evap_slurry = get_evap_pool(tg, ratm, qbot)
    infiltr_slurry = max(depth_slurry / poolranges(1), precip)
    infiltrated = depth_slurry * infiltr_slurry / (infiltr_slurry + evap_slurry)
    ! Slurry water (in addition to soil water, theta) on surface and in soil. Represents
    ! mean over pool S0.
    water_slurry = (/0.5*depth_slurry, 0.5*infiltrated/) 
    ! The excess water assumed to have percolated down from the volat. layer.
    percolated = max(infiltrated - dz_layer*(thetasat-theta), 0.0)
    ! Percolation rate out of volat layer, average over the pool S0.
    perc_slurry_mean = percolated / poolranges(1)

    call eval_fluxes_slurry(water_slurry, tanpools(1), Hconc(1), tg, ratm, theta, thetasat, perc_slurry_mean, &
         runoff, cnc_nh3_air, fluxes(:,1))

    if (any(isnan(fluxes))) then
       status = err_nan * 10
    end if
    
    !tanpools(1) = tanpools(1) - sum(fluxes(:,1)) * dt
    call update_pools(tanpools(1:1), fluxes(1:5,1:1), dt, 1, 5)

    if (any(tanpools < -1e-15)) then
       status = err_negative_tan
       return
    end if
    
    ! Pool aging & input
    !
    call age_pools_slurry(tandep, dt, water_slurry, tanpools(1), tanpools(2:), poolranges, garbage)
    ! TAN produced (mineralization) goes to directly the old TAN pool. 
    tanpools(4) = tanpools(4) + tanprod*dt

    ! Soil bins S1 and S2
    !
    age_prev = 0 ! for water evaluations, consider beginning of S1 as the starting point
    water_in_layer = infiltrated - percolated ! water in layer just after slurry has infiltrated
    do indpl = 2, 4
       ! water content lost during the aging
       waterloss = water_in_layer * (waterfunction(age_prev) - waterfunction(age_prev+poolranges(indpl)))
       percolation = eval_perc(waterloss, evap, precip, watertend, poolranges(indpl))
       !print *, water_in_layer*waterfunction(age_prev), water_in_layer*waterfunction(age_prev+poolranges(indpl))
       ! water content at the mean age of the pool
       water_soil = water_in_layer * waterfunction(age_prev + 0.5*poolranges(indpl))
       !print *, tanpools(indpl), water_soil, Hconc(indpl), tg, &
       !     & ratm, theta, thetasat, percolation, runoff, cnc_nh3_air, &
       !     & dz_layer
       call eval_fluxes_soil(tanpools(indpl), water_soil, Hconc(indpl), tg, &
            & ratm, theta, thetasat, percolation, runoff, cnc_nh3_air, &
            & dz_layer, fluxes(:,indpl), subst_tan, status)
       if (status /= 0) return
       !print *, fluxes(4,indpl), percolation, waterloss
       !fluxes(:,indpl) = 0
       !fluxes(5:,indpl) = 0
       !tanpools(indpl) = tanpools(indpl) - sum(fluxes(:,indpl)) * dt
       age_prev = age_prev + poolranges(indpl)
    end do
    
    call update_pools(tanpools(2:), fluxes(1:5,2:), dt, 3, 5)
    
    if (any(tanpools < -1e-15)) then
       !if (any(tanpools < -1e-3)) then
       status = err_negative_tan * 10
       return
       !end if
    end if

    if (any(isnan(fluxes))) then
       status = err_nan * 100
    end if

    if (abs(sum(tanpools - tanpools_old) - (-sum(fluxes) + tandep + tanprod)*dt + garbage) > max(sum(tanpools_old)*1e-2, 1e-4)) then
       !print *, sum(tanpools_old), sum(tanpools), sum(tanpools - tanpools_old)
       !print *, sum(fluxes), tandep*dt, tanprod*dt
       status = err_balance_tan
       return
    end if
    !print *, sum(tanpools), sum(tanpools - tanpools_old) !+ garbage
    
    status = 0
    
  end subroutine update_4pool
  
  subroutine update_npool(tg, ratm, theta, thetasat, precip, evap, qbot, watertend, runoff, tandep, tanprod, &
       water_init, cnc_nh3_air, poolranges, Hconc, dz_layer, tanpools, fluxes, garbage, dt, status, numpools)
    !
    ! Evaluate fluxes and update TAN pools for a model with arbitrary number of pools
    ! divided by age and pH.
    ! 
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: ratm ! atmospheric resistance, s/m
    real(r8), intent(in) :: theta ! volumetric soil water in soil column (unaffected by slurry)
    real(r8), intent(in) :: thetasat ! vol. soil water at saturation
    real(r8), intent(in) :: precip   ! precipitation, m/s
    real(r8), intent(in) :: evap   ! ground evaporation, m/s
    real(r8), intent(in) :: qbot   ! specific humidity (kg/kg) at lowest atmospheric model level
    real(r8), intent(in) :: watertend ! time derivative of theta*dz
    real(r8), intent(in) :: runoff    ! surface runoff flux, m/s
    real(r8), intent(in) :: tandep    ! TAN input flux, gN/m2/s
    real(r8), intent(in) :: tanprod(numpools)   ! flux of TAN produced (from urea/organic n) in the column
    real(r8), intent(in) :: water_init ! Initial water volume in the affected patch, m
    real(r8), intent(in) :: cnc_nh3_air ! NH3 concentration in air, gN/m3
    real(r8), intent(in) :: poolranges(numpools)  ! age ranges of TAN pools (npools)
    real(r8), intent(in) :: Hconc(numpools)       ! H+ concentration, mol/l (npools)
    real(r8), intent(in) :: dz_layer ! thickness of the volatilization layer, m
    real(r8), intent(inout) :: tanpools(numpools) ! TAN pools gN/m2 (npools)
    real(r8), intent(out) :: fluxes(5,numpools) ! TAN fluxes, gN/m2/s (type of flux, pool)
    real(r8), intent(out) :: garbage     ! "over-aged" TAN produced during the step, gN/m.
    real(r8), intent(in) :: dt ! timestep, sec, >0
    integer, intent(out) :: status ! 0 == OK
    integer, intent(in) :: numpools
    
    real(r8) :: fraction_layer, fraction_reservoir, fraction_runoff, waterloss, direct_runoff
    real(r8) :: percolation, water_soil, age_prev, tandep_remaining, direct_percolation, water_into_layer
    real(r8) :: tanpools_old(size(tanpools)), imbalance
    integer :: indpl
    
    real(r8), parameter :: water_relax_t = 24*3600.0_r8
    logical :: fixed

    tanpools_old = tanpools
    
    if (theta > thetasat) then
       !print *, 'bad theta update_npool 1', theta, thetasat
       status = err_bad_theta
       return
    end if
    
    ! Initial water excess goes to runoff if the surface is close to saturation, otherwise to the soil.
    !
    call partition_to_layer(water_init, theta, thetasat, dz_layer, fraction_layer, fraction_reservoir, fraction_runoff)
    direct_runoff = fraction_runoff * tandep
    direct_percolation = fraction_reservoir * tandep
    tandep_remaining = tandep - direct_runoff - direct_percolation
    water_into_layer = water_init * (1.0_r8 - fraction_reservoir - fraction_runoff)
    if (tandep_remaining < -1e-15) then
       status = err_negative_tan + 10
       return
    end if
    if (water_into_layer/dz_layer + theta > thetasat+1e-5) then
       status = err_bad_theta
       return
    end if

    if(any(isnan(tanpools))) then
       status = err_nan 
       return
    end if
    
    ! Pool aging & TAN input
    !
    call age_pools_soil(tandep_remaining, dt, poolranges, tanpools, garbage)
    ! TAN produced (mineralization) goes to directly the old TAN pool. 
    
    if (any(tanpools < 0)) then
       if (any(tanpools < -1e-15)) then
          print *, '<0 2', tanpools_old, 'new', tanpools, 'fx', &
               sum(fluxes(:,1)) * dt, sum(fluxes(:,2)) * dt, sum(fluxes(:,3)) * dt
          status = err_negative_tan + 10000
          return
       else
          where(tanpools<0) tanpools = 0.0
       end if
    end if

    imbalance = abs((sum(tanpools) - sum(tanpools_old)) - ((tandep_remaining)*dt+garbage))
    if (imbalance > max(1e-14, 0.001*sum(tanpools_old))) then
       print *, imbalance, 'old', tanpools_old, 'new', tanpools, 'g', garbage, tandep_remaining, &
            'diff', sum(tanpools_old-tanpools), &
            'depprod', tandep+sum(tanprod), garbage
       status = err_balance_tan*10
       return
    end if
    
    age_prev = 0 ! for water evaluations, consider beginning of S1 as the starting point
    do indpl = 1, size(tanpools)
       ! water content lost during the aging
       waterloss = water_into_layer * (waterfunction(age_prev) - waterfunction(age_prev+poolranges(indpl)))
       percolation = eval_perc(waterloss, evap, precip, watertend, poolranges(indpl))
       ! water content at the middle of the age range
       water_soil = water_into_layer * waterfunction(age_prev + 0.5*poolranges(indpl))
       !call eval_fluxes_soil(tanpools(indpl), water_soil, Hconc(indpl), tg, &
       !call eval_fluxes_soil_3p(tanpools(indpl), water_soil, Hconc(indpl), tg, &
       call eval_fluxes_soilroff(tanpools(indpl), water_soil, Hconc(indpl), tg, &
            & ratm, theta, thetasat, percolation, runoff, cnc_nh3_air, &
            & dz_layer, fluxes(:,indpl), subst_tan, status)
       if (status /= 0) then
          return
       end if
       age_prev = age_prev + poolranges(indpl)
    end do
    
    call update_pools(tanpools, fluxes, dt, numpools, 5, fixed)
    !print *, 'pools just after', tanpools
    tanpools = tanpools + tanprod*dt
    if(any(isnan(tanpools))) then
       status = err_nan+100
       return
    end if
    
    if (any(isnan(fluxes))) then
       status = err_nan + 1000
    end if
    
    if (any(tanpools < -1e-15)) then
       status = err_negative_tan + 1000
       print *, '<0 2', tanpools_old, tanpools, sum(fluxes(:,1)) * dt, sum(fluxes(:,2)) * dt
       
       return
    end if
    
    
    if (abs(sum(tanpools - tanpools_old) + (sum(fluxes)-tandep_remaining-sum(tanprod))*dt + garbage) &
         > max(sum(tanpools_old)*1e-2, 1d-2)) then
       print *, tanpools, tanpools_old, 'fx', fluxes*dt, 'dp', tandep_remaining*dt, tanprod*dt, 'g', garbage, &
            'ib', sum(tanpools-tanpools_old), (sum(fluxes)-tandep_remaining-sum(tanprod))*dt + garbage, 'fix', fixed
       status = err_balance_tan
       return
    end if

    ! Add the "direct" fluxes to the fluxes of the first pool 
    fluxes(iflx_roff, 1) = fluxes(iflx_roff, 1) + direct_runoff
    fluxes(iflx_soilq, 1) = fluxes(iflx_soilq, 1) + direct_percolation
    
    if (any(fluxes < -1e-6)) then
       print *, fluxes
       status = err_negative_flux
       return
    end if
    
    status = 0
    
  end subroutine update_npool

  subroutine update_pools(tanpools, fluxes, dt, np, nf, fixed)
    ! Update tan pools using the fluxes and an ad-hoc scheme agains negative TAN masses.
    implicit none
    real(r8), intent(inout) :: tanpools(np), fluxes(nf,np)
    real(r8), intent(in) :: dt
    integer, intent(in) :: np, nf
    logical, intent(out), optional :: fixed

    integer :: ip
    real(r8) :: sumflux, ff
    logical :: fixed_

    fixed_ = .false.
    do ip = 1, np
       sumflux = sum(fluxes(:,ip))*dt
       if (sumflux > tanpools(ip)) then
          if (sumflux > 1e-15) then
             fixed_ = .true.
             ff = tanpools(ip) / sumflux
             fluxes(:,ip) = fluxes(:,ip) * ff
             sumflux = tanpools(ip)
             sumflux = sum(fluxes(:,ip))*dt
          else
             sumflux = 0.0
          end if
       end if
       tanpools(ip) = tanpools(ip) - sumflux
    end do
    if (present(fixed)) fixed = fixed_
    
  end subroutine update_pools
  
  function get_evap_pool(tg, ratm, qbot) result(evap)
    ! Evaluate evaporation rate for surface water given spcific humidity at the reference
    ! height.
    implicit none
    real(r8), intent(in) :: tg, ratm, qbot
    real(r8) :: evap ! m/s

    real(r8) :: es, esdt, qs, qsdt, dens, flux
    real(r8), parameter :: press = 101300.0_r8
    
    call qsat(tg, press, es, esdt, qs, qsdt)
    if (qbot > qs) then
       evap = 0
       return
    end if

    dens = press / (SHR_CONST_RDAIR * tg)
    flux = dens * (qs - qbot) / ratm ! kg/s/m2 == mm/s
    evap = flux*1e-3
    
  end function get_evap_pool

  function waterfunction(pool_age) result(water)
    implicit none
    real(r8), intent(in) :: pool_age ! sec
    real(r8) :: water

    water = exp(-pool_age / water_relax_t)
  end function waterfunction

  function eval_perc(waterloss, evap, precip, watertend, dt) result(rate)
    !
    ! Evaluate the downwards water flux at the layer bottom given the infiltration and
    ! evaporation fluxes.
    implicit none
    real(r8), intent(in) :: waterloss ! total water loss during dt, m
    real(r8), intent(in) :: evap ! average evaporation rate, m/s
    real(r8), intent(in) :: precip ! average infiltration rate, m/s
    real(r8), intent(in) :: watertend ! background water tendency, m/s
    real(r8), intent(in) :: dt ! timespan, s

    real(r8) :: rate ! percolation rate, m/s
    real(r8) :: perc_base, perc_adj
    
    perc_base = waterloss / dt
    perc_adj = perc_base + precip - evap - watertend
    rate = max(perc_adj, 0.0)

  end function eval_perc

  subroutine eval_fluxes_storage(nitr_input, tempr_outside, windspeed, fract_direct, &
       volat_coef_barns, volat_coef_stores, &
       tan_fract_excr, fluxes_nitr, fluxes_tan, status)
    !
    ! Evaluate nitrogen fluxes in animal housings and storage. Only volatilization losses
    ! are assumed. The volatilization fluxes are assumed to depend linearly on the TAN
    ! fluxes entering the housings or storage. The base coefficients are given as
    ! arguments and adjusted according the model of Gyldenkaerne et al.
    ! 
    implicit none
    real(r8), intent(in) :: nitr_input ! total nitrogen excreted by animals in housings
    real(r8), intent(in) :: tempr_outside ! K
    real(r8), intent(in) :: windspeed ! m/s
    real(r8), intent(in) :: fract_direct
    real(r8), intent(in) :: volat_coef_barns, volat_coef_stores
    real(r8), intent(in) :: tan_fract_excr ! fraction of NH4 nitrogen in excreted N
    real(r8), intent(out) :: fluxes_nitr(4), fluxes_tan(4)
    integer, intent(out) :: status

    ! parameters for the Gyldenkaerne et al. parameterization
    real(r8), parameter :: Tfloor_barns = 4.0_8, Tfloor_stores = 1.0_8
    real(r8), parameter :: Tmin_barns = 0_8
    real(r8), parameter :: Tmax_barns = 12.5_8
    real(r8), parameter :: tempr_D = 3.0
    real(r8), parameter :: Vmin_barns = 0.2_8
    real(r8), parameter :: Vmax_barns = 0.228_8
    real(r8), parameter :: pA = 0.89_8, pB = 0.26_8  
    
    real(r8) :: flux_avail, flux_avail_tan, tempr_stores, tempr_barns, vent_barns, flux_direct, flux_direct_tan, &
         & flux_barn, flux_store, tempr_C
    
    fluxes_nitr = 0.0_r8
    fluxes_tan = 0.0_r8
    tempr_C = tempr_outside - 273
    
    tempr_barns = max(tempr_C+tempr_D, Tfloor_barns)
    if (tempr_C < Tmin_barns) then
       vent_barns = Vmin_barns
    else if (tempr_C > Tmax_barns) then
       vent_barns = Vmax_barns
    else
       vent_barns = Vmin_barns + tempr_C * (Vmax_barns-Vmin_barns) / (Tmax_barns - Tmin_barns)
    end if

    flux_avail = nitr_input
    flux_avail_tan = nitr_input * tan_fract_excr

    if (flux_avail < -1e-15 .or. flux_avail_tan < -1e-15) then
       status = err_negative_flux
       return
    end if
    
    flux_barn = flux_avail_tan * volat_coef_barns * tempr_barns**pA * vent_barns**pB

    fluxes_tan(iflx_air_barns) = flux_barn
    fluxes_nitr(iflx_air_barns) = flux_barn
    
    flux_avail = flux_avail - flux_barn
    flux_avail_tan = flux_avail_tan - flux_barn

    if (flux_avail < 0 .or. flux_avail_tan < 0) then
       status = err_negative_flux
       return
    end if
    
    flux_direct = fract_direct * flux_avail
    flux_avail = flux_avail - flux_direct
    flux_direct_tan = flux_avail_tan * fract_direct
    flux_avail_tan = flux_avail_tan - flux_direct_tan

    fluxes_tan(iflx_appl) = flux_direct_tan
    fluxes_nitr(iflx_appl) = flux_direct

    tempr_stores = max(Tfloor_stores, tempr_C)
    flux_store = flux_avail_tan &
         & * volat_coef_stores * tempr_stores**pA * windspeed**pB
    flux_store = min(flux_avail_tan, flux_store)

    fluxes_tan(iflx_air_stores) = flux_store
    fluxes_nitr(iflx_air_stores) = flux_store
    
    flux_avail = flux_avail - flux_store
    flux_avail_tan = flux_avail_tan - flux_store
    if (flux_avail < 0) then
       !print *, 'stores'
       status = err_negative_flux
       return
    end if

    fluxes_nitr(iflx_to_store) = flux_avail
    fluxes_tan(iflx_to_store) = flux_avail_tan
    
    if (abs(sum(fluxes_nitr) - nitr_input) > 1e-5*nitr_input) then
       !print *, fluxes_nitr, sum(fluxes_nitr), nitr_input
       status = err_balance_nitr
       return
    end if

    if (abs(sum(fluxes_tan) - nitr_input*tan_fract_excr) > 1e-5*nitr_input) then
       status = err_balance_tan
       return
    end if

    if (any(fluxes_nitr < 0) .or. any(fluxes_tan < 0)) then
       !print *, 'final'
       status = err_negative_flux
       return
    end if
    
    status = 0
    
  end subroutine eval_fluxes_storage

  subroutine update_org_n(flux_input, tg, pools, dt, tanprod, soilflux)
    !
    ! Evaluate the decomposition/mineralization N fluxes from the available, resistant and
    ! unavailable N fractions, and update the organic N pools. In addition, evaluate the
    ! flux of organic N into the soil pools according to a fixed time constant set below.
    implicit none
    real(r8), intent(in) :: flux_input(3) ! organic N entering the pools. gN/m2/s. For
                                          ! indices see at top of the module.
    real(r8), intent(in) :: tg          ! ground temperature, K
    real(r8), intent(inout) :: pools(3) ! organic N pools
    real(r8), intent(in) :: dt       ! timestep, sec
    real(r8), intent(out) :: tanprod(3)      ! Flux of TAN formed, both pools
    real(r8), intent(out) :: soilflux  ! Flux of organic nitrogen to soil

    real(r8) :: rate_res, rate_avail, TR
    real(r8), parameter :: ka1 = 8.94e-7_r8, ka2 = 6.38e-8 ! 1/s
    real(r8), parameter :: tr1 = 0.0106, tr2 = 0.12979 
    real(r8), parameter :: org_to_soil_time = 365*24*3600.0_r8

    real(r8) :: soilfluxes(3)
    
    TR = tr1 * exp(tr2 * (tg-273.15_r8))
    tanprod(ind_avail) = ka1 * TR * pools(ind_avail)
    tanprod(ind_resist) = ka2 * TR * pools(ind_resist)
    tanprod(ind_unavail) = 0.0
    soilfluxes = pools * 1.0_r8 / org_to_soil_time

    pools = pools + (flux_input - tanprod - soilfluxes) * dt
    soilflux = sum(soilfluxes)
    
  end subroutine update_org_n

  subroutine update_urea(tg, theta, thetasat, precip, evap, watertend, runoff, &
       ndep, pools, fluxes, garbage, ranges, dt, status, numpools)
    !
    ! Evaluate fluxes and update the urea pools. The procedure is similar to updating the
    ! soil TAN pools, but NO3 and volatilization fluxes do not occur.
    ! 
    implicit none
    real(r8), intent(in) :: tg ! soil temperature, K
    real(r8), intent(in) :: theta ! volumetric soil water in soil column (background)
    real(r8), intent(in) :: thetasat ! vol. soil water at saturation
    real(r8), intent(in) :: precip   ! precipitation, m/s
    real(r8), intent(in) :: evap   ! ground evaporation, m/s
    real(r8), intent(in) :: watertend ! time derivative of theta*dz
    real(r8), intent(in) :: runoff    ! surface runoff flux, m/s
    real(r8), intent(in) :: ndep
    real(r8), intent(inout) :: pools(numpools)
    real(r8), intent(out) :: fluxes(6, numpools) ! one extra for the to_tan flux
    real(r8), intent(in) :: ranges(numpools)
    real(r8), intent(out) :: garbage
    real(r8), intent(in) :: dt
    integer, intent(out) :: status
    integer, intent(in) :: numpools

    real(r8), parameter :: rate = 4.83e-6 ! 1/s
    real(r8), parameter :: missing = 1e36 ! for the parameters not needed for urea fluxes
    real(r8), parameter :: dz_layer = 0.02 ! thickness of the volatilization layer, m
    
    real(r8) :: age_prev, percolation, old_total, balance
    integer :: indpl

    old_total = sum(pools)
    
    call age_pools_soil(ndep, dt, ranges, pools, garbage)

    age_prev = 0 
    do indpl = 1, numpools
       percolation = eval_perc(0.0_r8, evap, precip, watertend, ranges(indpl))
       call eval_fluxes_soilroff(pools(indpl), 0.0_r8, missing, tg, &
       !call eval_fluxes_soil(pools(indpl), 0.0_r8, missing, tg, &
            & missing, theta, thetasat, percolation, runoff, missing, &
            & dz_layer, fluxes(1:5,indpl), subst_urea, status)
       if (status /= 0) then
          return
       end if
       fluxes(iflx_to_tan, indpl) = rate*pools(indpl)
       age_prev = age_prev + ranges(indpl)
    end do

    ! Here goes also flux_tan!
    call update_pools(pools, fluxes, dt, numpools, 6)

    balance = sum(pools) - old_total
    if (abs(balance - (ndep-sum(fluxes))*dt + garbage) > 1e-9) then
       print *, balance, 'f', sum(fluxes)*dt, ndep*dt, (ndep-sum(fluxes))*dt-garbage - balance, &
            'p', pools, 'g', garbage
       status = err_balance_nitr
       return
    end if

    status = 0
    
  end subroutine update_urea
  
  subroutine get_storage_fluxes_tan_ar(manure_excr, tempr_outside, windspeed, fract_direct, &
     & flux_direct, flux_direct_tan, flux_barn, flux_store, flux_resid, flux_resid_tan, &
     & volat_target_barns, volat_target_stores, volat_coef_barns, volat_coef_stores, tan_fract_excr, nn)

    real(8), intent(in), dimension(nn) :: manure_excr, tempr_outside, windspeed, fract_direct
    real(8), intent(out), dimension(nn) :: flux_barn, flux_store, flux_direct, flux_resid, &
         & flux_direct_tan, flux_resid_tan
    real(8), intent(in) :: volat_target_barns, volat_target_stores, volat_coef_barns, volat_coef_stores, tan_fract_excr
    integer, intent(in) :: nn

    integer :: ii, status
    real(r8) :: fluxes_nitr(4), fluxes_tan(4)
    
    do ii = 1, nn
       call eval_fluxes_storage(manure_excr(ii), tempr_outside(ii), windspeed(ii), fract_direct(ii), &
            & volat_coef_barns, volat_coef_stores, tan_fract_excr, &
            & fluxes_nitr, fluxes_tan, status)

       if (status /= 0) then
          print *, 'Status = ', status
          return
       end if

       flux_direct(ii) = fluxes_nitr(iflx_appl)
       flux_direct_tan(ii) = fluxes_tan(iflx_appl)
       flux_barn(ii) = fluxes_tan(iflx_air_barns)
       flux_store(ii) = fluxes_tan(iflx_air_stores)
       flux_resid(ii) = fluxes_nitr(iflx_to_store)
       flux_resid_tan(ii) = fluxes_tan(iflx_to_store)
       !print *, '1', fluxes_nitr(iflx_appl), flux_direct(ii)
       
    end do
  end subroutine get_storage_fluxes_tan_ar
  
end module FanMod
