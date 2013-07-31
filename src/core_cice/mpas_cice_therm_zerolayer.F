module cice_therm_zerolayer

  use mpas_grid_types

  implicit none

  private

  public :: therm_zerolayer, &
            init_semtner_forcing, &
            init_semtner_test_case

  logical :: ltest = .false.

  ! thermodynamic constants
  real(kind=RKIND), parameter :: &
       cal_to_joule                 = 4.184_RKIND, &
       joule_to_cal                 = 1.0_RKIND / cal_to_joule, &
       kcal_to_cal                  = 1000.0_RKIND, &
       cal_to_kcal                  = 1.0_RKIND / kcal_to_cal, &
       m_to_cm                      = 100.0_RKIND, &
       cm_to_m                      = 1.0_RKIND / m_to_cm, &
       iceConductivitySemtner       = 4.86e-3_RKIND, & ! cal cm-1 s-1 C-1
       snowConductivitySemtner      = 7.4e-4_RKIND, &  ! cal cm-1 s-1 C-1
       iceConductivity              = iceConductivitySemtner  * cal_to_joule * m_to_cm, & ! J m-1 s-1 C-1
       snowConductivity             = snowConductivitySemtner * cal_to_joule * m_to_cm, & ! J m-1 s-1 C-1
       conductivityCorrection       = 1.065_RKIND, &
       iceConductivityEffective     = iceConductivity  * conductivityCorrection, &
       snowConductivityEffective    = snowConductivity * conductivityCorrection, &
       freshWaterFreezingPoint      = 273.15_RKIND, &
       stefanBoltzmann              = 5.79484e-8_RKIND, & ! J m-2 K-4 s-1
       airSaturationHumidityCoeff1  = 11637800.0_RKIND, &
       airSaturationHumidityCoeff2  = 5897.8_RKIND, &
       iceAlbedo                    = 0.64_RKIND, &
       snowAlbedoConstant           = 0.85_RKIND, &
       internalAbsorptionCorrection = 0.4_RKIND, &
       fractionPenetratingRadiation = 0.17_RKIND, &
       liquidusSlope                = 0.054_RKIND, &
       iceEnergyDensitySemtner      = 72.0_RKIND, & ! cal cm-3
       snowEnergyDensitySemtner     = 26.2_RKIND, & ! cal cm-3
       basalEnergyDensitySemtner    = 64.0_RKIND, & ! cal cm-3
       iceEnergyDensity             = iceEnergyDensitySemtner   * cal_to_joule * m_to_cm**3, & ! J m-3
       snowEnergyDensity            = snowEnergyDensitySemtner  * cal_to_joule * m_to_cm**3, & ! J m-3
       basalEnergyDensity           = basalEnergyDensitySemtner * cal_to_joule * m_to_cm**3, & ! J m-3
       snowMeltingTemperature       = 0.0_RKIND, &
       iceMeltingTemperature        = -0.1_RKIND

  !-------------------------------------------------------------
  ! Semtner forcing
  !-------------------------------------------------------------

    real(kind=RKIND), dimension(12), parameter :: shortwaveFluxMonthly = & ! kcal cm-2
         (/ 0.0_RKIND,  0.0_RKIND,  1.9_RKIND, &
            9.9_RKIND, 17.7_RKIND, 19.2_RKIND, &
           13.6_RKIND,  9.0_RKIND,  3.7_RKIND, &
            0.4_RKIND,  0.0_RKIND,  0.0_RKIND/)

    real(kind=RKIND), dimension(12), parameter :: longwaveFluxMonthly = & ! kcal cm-2
         (/10.4_RKIND, 10.3_RKIND, 10.3_RKIND, &
           11.6_RKIND, 15.1_RKIND, 18.0_RKIND, &
           19.1_RKIND, 18.7_RKIND, 16.5_RKIND, &
           13.9_RKIND, 11.2_RKIND, 10.9_RKIND/)
    
    real(kind=RKIND), dimension(12), parameter :: sensibleFluxMonthly = & ! kcal cm-2
         (/ 1.18_RKIND,  0.76_RKIND,  0.72_RKIND, &
            0.29_RKIND, -0.45_RKIND, -0.39_RKIND, &
           -0.30_RKIND, -0.40_RKIND, -0.17_RKIND, &
            0.10_RKIND,  0.56_RKIND,  0.79_RKIND/)

    real(kind=RKIND), dimension(12), parameter :: latentFluxMonthly = & ! kcal cm-2
         (/ 0.00_RKIND, -0.02_RKIND, -0.03_RKIND, &
           -0.09_RKIND, -0.46_RKIND, -0.70_RKIND, &
           -0.64_RKIND, -0.66_RKIND, -0.39_RKIND, &
           -0.19_RKIND, -0.01_RKIND, -0.01_RKIND/)

    real(kind=RKIND), dimension(12), parameter :: snowAlbedoMonthly = &
         (/0.83_RKIND, 0.83_RKIND, 0.83_RKIND, &
           0.81_RKIND, 0.82_RKIND, 0.78_RKIND, &
           0.64_RKIND, 0.69_RKIND, 0.84_RKIND, &
           0.85_RKIND, 0.85_RKIND, 0.85_RKIND/)

    real(kind=RKIND), dimension(-1:14) :: &
         shortwaveFluxMonthlyExpand, &
         longwaveFluxMonthlyExpand, &
         sensibleFluxMonthlyExpand, &
         latentFluxMonthlyExpand, &
         snowAlbedoMonthlyExpand

    real(kind=RKIND), dimension(-1:14) :: &
         semtner_data_times

    real(kind=RKIND), dimension(0:12,4) :: &
         shortwaveFluxMonthlyIntCoeff, &
         longwaveFluxMonthlyIntCoeff, &
         sensibleFluxMonthlyIntCoeff, &
         latentFluxMonthlyIntCoeff, &
         snowAlbedoMonthlyIntCoeff

    integer, dimension(12), parameter :: days_in_month_real = &
         (/31,28,31,30,31,30,31,31,30,31,30,31/)

    integer, dimension(12), parameter :: days_in_month = &
         (/30,30,30,30,30,30,30,30,30,30,30,30/)


    integer, dimension(12), parameter :: cumulative_days_in_month = &
!         (/0,31,59,90,120,151,181,212,243,273,304,334/)
         (/0,30,60,90,120,150,180,210,240,270,300,330/)

    integer, parameter :: &
         seconds_in_day = 24.0_RKIND * 3600.0_RKIND, &
!         seconds_in_year = seconds_in_day * 365.0_RKIND
         seconds_in_year = seconds_in_day * 360.0_RKIND

    integer, parameter :: nSnowPeriods = 4

    real(kind=RKIND), dimension(nSnowPeriods,2) :: &
         snowfallTimes

    real(kind=RKIND), dimension(nSnowPeriods) :: &
         snowfallDuration, &
         snowfallAccumulationRate    

contains

  !-------------------------------------------------------------

  subroutine therm_zerolayer(mesh, icestate, atmos_forcing, ocean_forcing, dt, time)

    type(mesh_type), intent(inout) :: mesh

    type(icestate_type), pointer :: icestate
    type(atmos_forcing_type), pointer :: atmos_forcing
    type(ocean_forcing_type), pointer :: ocean_forcing

    real(kind=RKIND), intent(in) :: &
         dt, &
         time

    real(kind=RKIND) :: &
         cumulativeIceAblation, &
         cumulativeBasalAblation, &
         cumulativeBasalAccretion
         
    integer :: &
         iCell

    !call run_semtner_test_case()

    do iCell = 1, mesh % nCells

       if (icestate % iceAreaCell % array(iCell) > 0.0_RKIND) then

          if (iCell == 7134 .or. iCell == 7133) ltest = .true.

          call therm_zerolayer_column( &
               icestate % iceAreaCell % array(iCell), &
               icestate % iceVolumeCell % array(iCell), &
               icestate % snowVolumeCell % array(iCell), &
               icestate % surfaceTemperature % array(iCell), &
               atmos_forcing % shortwaveDown % array(iCell), &
               atmos_forcing % longwaveDown % array(iCell), &
               atmos_forcing % airDensity % array(iCell), &
               atmos_forcing % airSpecificHumidity % array(iCell), &
               atmos_forcing % airPotentialTemperature % array(iCell), &
               atmos_forcing % sensibleTransferCoefficient % array(iCell), &
               atmos_forcing % latentTransferCoefficient % array(iCell), &
               atmos_forcing % snowfallRate % array(iCell), &
               ocean_forcing % seaSurfaceSalinity % array(iCell), &
               cumulativeIceAblation, &
               cumulativeBasalAblation, &
               cumulativeBasalAccretion, &
               dt, &
               time)

          ltest = .false.

          !write(*,*) iCell, icestate % surfaceTemperature % array(iCell)

       else

          call add_new_ice( &
               icestate % iceAreaCell % array(iCell), &
               icestate % iceVolumeCell % array(iCell), &
               icestate % snowVolumeCell % array(iCell), &
               icestate % surfaceTemperature % array(iCell), &
               atmos_forcing % airPotentialTemperature % array(iCell))
          
       endif ! iceAreaCell

    enddo ! iCell

    !write(*,*) minloc(icestate % surfaceTemperature % array), minval(icestate % surfaceTemperature % array)

!stop
  end subroutine therm_zerolayer

  !-------------------------------------------------------------

  subroutine add_new_ice( &
       iceAreaCell, &
       iceVolumeCell, &
       snowVolumeCell, &
       surfaceTemperature, &
       airPotentialTemperature)

    real(kind=RKIND), intent(inout) :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         surfaceTemperature

    real(kind=RKIND), intent(in) :: &
         airPotentialTemperature

    if (airPotentialTemperature < -2.0_RKIND) then

       iceAreaCell = 1.0_RKIND
       iceVolumeCell = 0.1_RKIND
       snowVolumeCell = 0.0_RKIND
       surfaceTemperature = -2.0_RKIND

    endif

  end subroutine add_new_ice

  !-------------------------------------------------------------

  subroutine therm_zerolayer_column( &
       iceAreaCell, &
       iceVolumeCell, &
       snowVolumeCell, &
       surfaceTemperature, &
       shortwaveDownForcing, &
       longwaveDownForcing, &
       airDensity, &
       airSpecificHumidity, &
       airPotentialTemperature, &
       sensibleTransferCoefficient, &
       latentTransferCoefficient, &
       snowfallRateForcing, &
       seaSurfaceSalinity, &
       cumulativeIceAblation, &
       cumulativeBasalAblation, &
       cumulativeBasalAccretion, &
       dt, &
       time)

    real(kind=RKIND), intent(inout) :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         surfaceTemperature, &
         cumulativeIceAblation, &
         cumulativeBasalAblation, &
         cumulativeBasalAccretion

    real(kind=RKIND), intent(in) :: &
         shortwaveDownForcing, &
         longwaveDownForcing, &
         airDensity, &
         airSpecificHumidity, &
         airPotentialTemperature, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         snowfallRateForcing, &
         seaSurfaceSalinity, &
         dt, &
         time

    real(kind=RKIND) :: &
         iceThickness, &
         snowThickness, &
         newSurfaceTemperature, &
         shortwaveDown, &
         longwaveDown, &
         snowAlbedo, &
         shortwaveAbsorbed, &
         bottomTemperature, &
         oceanicFlux, &
         snowfallRate, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         dsensibleHeatFlux_dTsf, &
         dlatentHeatFlux_dTsf

    ! determine thicknesses
    iceThickness  = iceVolumeCell  / iceAreaCell
    snowThickness = snowVolumeCell / iceAreaCell

    ! calculate the bottom temperature
    call calculate_basal_properties( &
         bottomTemperature, &
         oceanicFlux, &
         seaSurfaceSalinity)

    shortwaveDown = shortwaveDownForcing
    longwaveDown  = longwaveDownForcing
    snowAlbedo    = snowAlbedoConstant
    snowfallRate  = snowfallRateForcing

    ! calculate absorbed shortwave flux
    call absorbed_shortwave( &
         shortwaveAbsorbed, &
         shortwaveDown, &
         snowThickness, &
         snowAlbedo)

    ! calculate the new surface temperature
    call calculate_surface_temperature( &
         newSurfaceTemperature, &
         surfaceTemperature, &
         longwaveDown, &
         shortwaveAbsorbed,  &
         bottomTemperature, &
         iceThickness, &
         snowThickness, &
         airDensity, &
         airPotentialTemperature, &
         airSpecificHumidity, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient)

    ! sensible and latent heats with new temperatures
    call sensible_latent_heat_fluxes( &
         newSurfaceTemperature, &
         airDensity, &
         airPotentialTemperature, &
         airSpecificHumidity, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         dsensibleHeatFlux_dTsf, &
         dlatentHeatFlux_dTsf)

    ! surface melting
    call surface_melting( &
         newSurfaceTemperature, &
         iceThickness, &
         snowThickness, &
         longwaveDown, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         shortwaveAbsorbed, &
         bottomTemperature, &
         cumulativeIceAblation, &
         time, &
         dt)

    surfaceTemperature = newSurfaceTemperature

    ! basal melting or accretion
    call basal_melting_accretion( &
         iceThickness, &
         snowThickness, &
         surfaceTemperature, &
         bottomTemperature, &
         oceanicFlux, &
         cumulativeBasalAblation, &
         cumulativeBasalAccretion, &
         dt)

    ! snowfall accumulation
    call snowfallAccumulation( &
         snowThickness, &
         snowfallRate, &
         surfaceTemperature, &
         dt)
    
    ! all ice gone
    if (iceThickness == 0.0_RKIND) then
       snowThickness = 0.0_RKIND
       iceAreaCell = 0.0_RKIND
       surfaceTemperature = 0.0_RKIND
    end if

    ! determine volumes
    iceVolumeCell  = iceThickness  * iceAreaCell
    snowVolumeCell = snowThickness * iceAreaCell 

  end subroutine therm_zerolayer_column

  !-------------------------------------------------------------

  subroutine snowfallAccumulation( &
       snowThickness, &
       snowfallRate, &
       surfaceTemperature, &
       dt)

    real(kind=RKIND), intent(inout) :: &
         snowThickness

    real(kind=RKIND), intent(in) :: &    
         snowfallRate, &
         surfaceTemperature, &
         dt

    real(kind=RKIND), parameter :: &
         densitySnow = 330.0_RKIND
    
    ! only accumulate if freezing surface
    if ((surfaceTemperature < snowMeltingTemperature .and. snowThickness > 0.0_RKIND) .or. &
        (surfaceTemperature < iceMeltingTemperature .and. snowThickness == 0.0_RKIND)) then

       snowThickness = snowThickness + (snowfallRate / densitySnow) * dt

    endif

  end subroutine snowfallAccumulation

  !-------------------------------------------------------------

  subroutine basal_melting_accretion( &
       iceThickness, &
       snowThickness, &
       surfaceTemperature, &
       bottomTemperature, &
       oceanicFlux, &
       cumulativeBasalAblation, &
       cumulativeBasalAccretion, &
       dt)

    real(kind=RKIND), intent(inout) :: &
         iceThickness, &
         cumulativeBasalAblation, &
         cumulativeBasalAccretion

    real(kind=RKIND), intent(in) :: &
         snowThickness, &
         surfaceTemperature, &
         bottomTemperature, &
         oceanicFlux, &
         dt

    real(kind=RKIND) :: &
         conductivityCoefficient, &
         conductiveFlux, &
         growthEnergy, &
         iceThicknessChange

    call conductivity_coefficient(conductivityCoefficient, iceThickness, snowThickness)
    
    conductiveFlux = conductivityCoefficient * (surfaceTemperature - bottomTemperature)
    
    growthEnergy = (oceanicFlux - conductiveFlux) * dt
    
    iceThicknessChange = max(growthEnergy / basalEnergyDensity, -iceThickness)

    iceThickness = iceThickness + iceThicknessChange

    cumulativeBasalAblation  = cumulativeBasalAblation  + abs(min(iceThicknessChange, 0.0_RKIND))
    cumulativeBasalAccretion = cumulativeBasalAccretion + max(iceThicknessChange, 0.0_RKIND)

  end subroutine basal_melting_accretion

  !-------------------------------------------------------------

  subroutine surface_melting( &
       newSurfaceTemperature, &
       iceThickness, &
       snowThickness, &
       longwaveDown, &
       sensibleHeatFlux, &
       latentHeatFlux, &
       shortwaveAbsorbed, &
       bottomTemperature, &
       cumulativeIceAblation, &
       time, &
       dt)

    real(kind=RKIND), intent(inout) :: &
         newSurfaceTemperature, &
         iceThickness, &
         snowThickness, &
         cumulativeIceAblation

    real(kind=RKIND), intent(in) :: &
         longwaveDown, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         shortwaveAbsorbed, &
         bottomTemperature, &
         time, &
         dt

    real(kind=RKIND) :: &
         surfaceTemperatureKelvin, &
         conductivityCoefficient, &
         conductiveFlux, &
         atmosphericFlux, &
         meltingEnergy, &
         snowThicknessChange, &
         iceThicknessChange

    logical :: &
         lmelt
    
    lmelt = .false.

    ! surface temperature
    if (newSurfaceTemperature > snowMeltingTemperature .and. snowThickness > 0.0_RKIND) then
       newSurfaceTemperature = snowMeltingTemperature
       lmelt = .true.
    else if (newSurfaceTemperature > iceMeltingTemperature .and. snowThickness == 0.0_RKIND) then
       newSurfaceTemperature = iceMeltingTemperature
       lmelt = .true.
    endif

    if (lmelt) then
       
       surfaceTemperatureKelvin = newSurfaceTemperature + freshWaterFreezingPoint
       
       call conductivity_coefficient(conductivityCoefficient, iceThickness, snowThickness)
       
       conductiveFlux = conductivityCoefficient * (newSurfaceTemperature - bottomTemperature)
    
       atmosphericFlux = longwaveDown + sensibleHeatFlux + latentHeatFlux + shortwaveAbsorbed - &
            stefanBoltzmann * surfaceTemperatureKelvin**4
    
       meltingEnergy = max((atmosphericFlux - conductiveFlux) * dt, 0.0_RKIND)

       snowThicknessChange = min(meltingEnergy / snowEnergyDensity*0.9, snowThickness)

       snowThickness = snowThickness - snowThicknessChange

       meltingEnergy = meltingEnergy - snowThicknessChange * snowEnergyDensity
       
       iceThicknessChange = min(meltingEnergy / iceEnergyDensity, iceThickness)

       iceThickness = iceThickness - iceThicknessChange

       cumulativeIceAblation = cumulativeIceAblation + iceThicknessChange

       if (snowThickness == 0.0_RKIND) newSurfaceTemperature = iceMeltingTemperature

    endif

  end subroutine surface_melting
  
  !-------------------------------------------------------------

  subroutine calculate_surface_temperature( &
       newSurfaceTemperature, &
       surfaceTemperature, &
       longwaveDown, &
       shortwaveAbsorbed, &
       bottomTemperature, &
       iceThickness, &
       snowThickness, &
       airDensity, &
       airPotentialTemperature, &
       airSpecificHumidity, &
       sensibleTransferCoefficient, &
       latentTransferCoefficient)

    real(kind=RKIND), intent(out) :: &
         newSurfaceTemperature
    
    real(kind=RKIND), intent(in) :: &
         surfaceTemperature, &
         longwaveDown, &
         shortwaveAbsorbed, &
         bottomTemperature, &
         iceThickness, &
         snowThickness, &
         airDensity, &
         airPotentialTemperature, &
         airSpecificHumidity, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient

    real(kind=RKIND) :: &
         conductivityCoefficient, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         dsensibleHeatFlux_dTsf, &
         dlatentHeatFlux_dTsf, &
         f, &
         fp

    integer :: &
         nit

    integer, parameter :: &
         nit_max = 10
 
    newSurfaceTemperature = 0.0_RKIND!surfaceTemperature

    call conductivity_coefficient(&
         conductivityCoefficient, &
         iceThickness, &
         snowThickness)

    do nit = 1, nit_max

       call sensible_latent_heat_fluxes( &
            newSurfaceTemperature, &
            airDensity, &
            airPotentialTemperature, &
            airSpecificHumidity, &
            sensibleTransferCoefficient, &
            latentTransferCoefficient, &
            sensibleHeatFlux, &
            latentHeatFlux, &
            dsensibleHeatFlux_dTsf, &
            dlatentHeatFlux_dTsf)
       
       f = surface_temperature_function( &
            newSurfaceTemperature, &
            longwaveDown, &
            shortwaveAbsorbed, &
            sensibleHeatFlux, &
            latentHeatFlux, &
            conductivityCoefficient, &
            bottomTemperature)

       fp = surface_temperature_derivative( &
            newSurfaceTemperature, &
            dsensibleHeatFlux_dTsf, &
            dlatentHeatFlux_dTsf, &
            conductivityCoefficient)

       !if (ltest) write(*,*) nit, newSurfaceTemperature - f / fp, newSurfaceTemperature, longwaveDown, shortwaveAbsorbed, sensibleHeatFlux, latentHeatFlux, &
       !     -stefanBoltzmann * (newSurfaceTemperature + freshWaterFreezingPoint)**4, -conductivityCoefficient * (bottomTemperature - newSurfaceTemperature)

       !if (ltest) write(*,*) nit, newSurfaceTemperature, &
       !     airDensity, &
       !     airPotentialTemperature, &
       !     airSpecificHumidity, &
       !     sensibleTransferCoefficient, &
       !     latentTransferCoefficient, &
       !     sensibleHeatFlux, latentHeatFlux

       newSurfaceTemperature = newSurfaceTemperature - f / fp

    enddo

    !stop

  end subroutine calculate_surface_temperature

  !-------------------------------------------------------------

  function surface_temperature_function( &
       surfaceTemperature, &
       longwaveDown, &
       shortwaveAbsorbed, &
       sensibleHeatFlux, &
       latentHeatFlux, &
       conductivityCoefficient, &
       bottomTemperature) result(f)

    real(kind=RKIND), intent(in) :: &
         surfaceTemperature, &
         longwaveDown, &
         shortwaveAbsorbed, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         conductivityCoefficient, &
         bottomTemperature

    real(kind=RKIND) :: f

    real(kind=RKIND) :: &
         surfaceTemperatureKelvin

    surfaceTemperatureKelvin = surfaceTemperature + freshWaterFreezingPoint

    f = longwaveDown + shortwaveAbsorbed + sensibleHeatFlux + latentHeatFlux - &
        stefanBoltzmann * surfaceTemperatureKelvin**4 - &
        conductivityCoefficient * (surfaceTemperature - bottomTemperature)

  end function surface_temperature_function

  !-------------------------------------------------------------
  
  function surface_temperature_derivative( &
       surfaceTemperature, &
       dsensibleHeatFlux_dTsf, &
       dlatentHeatFlux_dTsf, &
       conductivityCoefficient) result(fp)

    real(kind=RKIND), intent(in) :: &
         surfaceTemperature, &
         dsensibleHeatFlux_dTsf, &
         dlatentHeatFlux_dTsf, &
         conductivityCoefficient

    real(kind=RKIND) :: fp

    real(kind=RKIND) :: &
         surfaceTemperatureKelvin

    surfaceTemperatureKelvin = surfaceTemperature + freshWaterFreezingPoint

    fp = dsensibleHeatFlux_dTsf + dlatentHeatFlux_dTsf - &
         4.0_RKIND * stefanBoltzmann * surfaceTemperatureKelvin**3 - &
         conductivityCoefficient

  end function surface_temperature_derivative

  !-------------------------------------------------------------
  
  subroutine sensible_latent_heat_fluxes( &
       surfaceTemperature, &
       airDensity, &
       airPotentialTemperature, &
       airSpecificHumidity, &
       sensibleTransferCoefficient, &
       latentTransferCoefficient, &
       sensibleHeatFlux, &
       latentHeatFlux, &
       dsensibleHeatFlux_dTsf, &
       dlatentHeatFlux_dTsf)

    real(kind=RKIND), intent(in) :: &
         surfaceTemperature              ! ice/snow surface temperature

    real(kind=RKIND), intent(in) :: &
         airDensity                  , & ! air density (kg/m^3)
         airPotentialTemperature     , & ! air potential temperature  (K)
         airSpecificHumidity         , & ! specific humidity (kg/kg)
         sensibleTransferCoefficient , & ! transfer coefficient for sensible heat
         latentTransferCoefficient       ! transfer coefficient for latent heat

    real(kind=RKIND), intent(out) :: &
         sensibleHeatFlux            , & ! surface downward sensible heat (W m-2)
         latentHeatFlux              , & ! surface downward latent heat (W m-2)
         dsensibleHeatFlux_dTsf      , &
         dlatentHeatFlux_dTsf

    real(kind=RKIND) :: &
         surfaceTemperatureKelvin        , & ! ice/snow surface temperature (K)
         saturatedSpecificHumidity       , & ! saturated surface specific humidity (kg/kg)
         airSaturationHumidity           , & ! the saturation humidity of air (kg/m^3)
         dsaturatedSpecificHumidity_dTsf 

    ! ice surface temperature in Kelvin
    surfaceTemperatureKelvin = surfaceTemperature + freshWaterFreezingPoint

    ! saturation humidity
    airSaturationHumidity = airSaturationHumidityCoeff1 * exp(-airSaturationHumidityCoeff2 / surfaceTemperatureKelvin)
    saturatedSpecificHumidity = airSaturationHumidity / airDensity

    dsaturatedSpecificHumidity_dTsf = (saturatedSpecificHumidity * airSaturationHumidityCoeff2) / surfaceTemperatureKelvin**2

    ! downward latent and sensible heat fluxes
    sensibleHeatFlux = sensibleTransferCoefficient * (airPotentialTemperature - surfaceTemperature)
    latentHeatFlux   = latentTransferCoefficient   * (airSpecificHumidity     - saturatedSpecificHumidity)

    ! derivatives of latent and sensible heats wrt surface temperature
    dsensibleHeatFlux_dTsf = -sensibleTransferCoefficient
    dlatentHeatFlux_dTsf   = -latentTransferCoefficient * dsaturatedSpecificHumidity_dTsf
 
  end subroutine sensible_latent_heat_fluxes

  !-------------------------------------------------------------

  subroutine absorbed_shortwave( &
       shortwaveAbsorbed, &
       shortwaveDown, &
       snowThickness, &
       snowAlbedo)

    real(kind=RKIND), intent(out) :: &
         shortwaveAbsorbed

    real(kind=RKIND), intent(in) :: &
         shortwaveDown, &
         snowThickness, &
         snowAlbedo

    real(kind=RKIND) :: &
         effectiveAlbedo

    if (snowThickness > 0.0_RKIND) then

       shortwaveAbsorbed = shortwaveDown * (1.0_RKIND - snowAlbedo) 

    else

       shortwaveAbsorbed = shortwaveDown * (1.0_RKIND - iceAlbedo) * (1.0_RKIND - fractionPenetratingRadiation) + &
                           shortwaveDown * (1.0_RKIND - iceAlbedo) * fractionPenetratingRadiation * (1.0_RKIND - internalAbsorptionCorrection)

    endif
    
  end subroutine absorbed_shortwave

  !-------------------------------------------------------------

  subroutine conductivity_coefficient(conductivityCoefficient, iceThickness, snowThickness)

    real(kind=RKIND), intent(out) :: &
         conductivityCoefficient

    real(kind=RKIND), intent(in) :: &    
         iceThickness, &
         snowThickness

    if (snowThickness > 0.0_RKIND) then
       
       conductivityCoefficient = snowConductivityEffective / &
            (snowThickness + ((iceThickness * snowConductivityEffective) / iceConductivityEffective))
       
    else
       
       conductivityCoefficient = iceConductivityEffective / iceThickness
       
    endif

  end subroutine conductivity_coefficient

  !-------------------------------------------------------------

  subroutine calculate_basal_properties( &
       bottomTemperature, &
       oceanicFlux, &
       seaSurfaceSalinity)
    
    real(kind=RKIND), intent(out) :: &
         bottomTemperature, &
         oceanicFlux
    
    real(kind=RKIND), intent(in) :: &
         seaSurfaceSalinity

    real(kind=RKIND), parameter :: &
         oceanicFluxSemtner = -1.5_RKIND ! kcal cm-2 year-1

    !bottomTemperature = -liquidusSlope * seaSurfaceSalinity
    bottomTemperature = -2.0_RKIND

    oceanicFlux = (oceanicFluxSemtner * kcal_to_cal * cal_to_joule * m_to_cm**2) / seconds_in_year

  end subroutine calculate_basal_properties

  !-------------------------------------------------------------
  ! Semtner forcing test case
  !-------------------------------------------------------------

  subroutine run_semtner_test_case()

    real(kind=RKIND), parameter :: &
         dt = 3600.0_RKIND

    real(kind=RKIND) :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         surfaceTemperature

    real(kind=RKIND) :: &
         time, &
         shortwaveDownForcing, &
         longwaveDownForcing, &
         airDensity, &
         airSpecificHumidity, &
         airPotentialTemperature, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         snowfallRateForcing, &
         seaSurfaceSalinity, &
         cumulativeIceAblation, &
         cumulativeBasalAblation, &
         cumulativeBasalAccretion, &
         year

    integer, parameter :: &
         nsteps = 65 * 12 * 30 * 24

    integer :: &
         istep

    shortwaveDownForcing = 0.0_RKIND
    longwaveDownForcing = 0.0_RKIND
    airDensity = 0.0_RKIND
    airSpecificHumidity = 0.0_RKIND
    airPotentialTemperature = 0.0_RKIND
    sensibleTransferCoefficient = 0.0_RKIND
    latentTransferCoefficient = 0.0_RKIND
    snowfallRateForcing = 0.0_RKIND
    seaSurfaceSalinity = 0.0_RKIND

    iceAreaCell = 1.0_RKIND
    iceVolumeCell = 1.0_RKIND
    snowVolumeCell = 0.0_RKIND
    surfaceTemperature = -10.0_RKIND

    cumulativeIceAblation = 0.0_RKIND
    cumulativeBasalAblation = 0.0_RKIND
    cumulativeBasalAccretion = 0.0_RKIND

    open(11,file='therm_zerolayer.txt')

    open(12,file='therm_zerolayer_cum.txt')

    do istep = 1, nsteps

       time = istep * dt

       call therm_zerolayer_column( &
            iceAreaCell, &
            iceVolumeCell, &
            snowVolumeCell, &
            surfaceTemperature, &
            shortwaveDownForcing, &
            longwaveDownForcing, &
            airDensity, &
            airSpecificHumidity, &
            airPotentialTemperature, &
            sensibleTransferCoefficient, &
            latentTransferCoefficient, &
            snowfallRateForcing, &
            seaSurfaceSalinity, &
            cumulativeIceAblation, &
            cumulativeBasalAblation, &
            cumulativeBasalAccretion, &
            dt, &
            time)

       if (mod(istep,24) == 0) then
          write(11,*) time, iceVolumeCell, snowVolumeCell, surfaceTemperature
       end if

       if (mod(istep,24*360) == 0) then

          year = ceiling(real(istep,RKIND)/real(24*365,RKIND))

          write(12,*) year, cumulativeIceAblation, cumulativeBasalAblation, cumulativeBasalAccretion, &
               cumulativeIceAblation + cumulativeBasalAblation - cumulativeBasalAccretion

          cumulativeIceAblation = 0.0_RKIND
          cumulativeBasalAblation = 0.0_RKIND
          cumulativeBasalAccretion = 0.0_RKIND

       endif

    enddo ! istep

    close(11)

    close(12)

  end subroutine run_semtner_test_case

  !-------------------------------------------------------------

  subroutine init_semtner_forcing()

    real(kind=RKIND) :: &
         shortwave_forcing, &
         longwave_forcing, &
         sensible_forcing, &
         latent_forcing, &
         snowalbedo_forcing

    integer :: &
         nMonth, &
         nMonthUse, &
         i

    real(kind=RKIND) :: &
         time

    do nMonth = -1, 14

       nMonthUse = modulo(nMonth-1,12) + 1

       semtner_data_times(nMonth) = (real(cumulative_days_in_month(nMonthUse), RKIND) + &
                                     0.5_RKIND * real(days_in_month(nMonthUse), RKIND)) * &
                                     seconds_in_day

    enddo ! nMonth

    semtner_data_times(-1:0)  = semtner_data_times(-1:0)  - seconds_in_year
    semtner_data_times(13:14) = semtner_data_times(13:14) + seconds_in_year

    semtner_data_times = semtner_data_times / real(seconds_in_year)

    ! expanded input data
    shortwaveFluxMonthlyExpand(1:12)  = shortwaveFluxMonthly
    shortwaveFluxMonthlyExpand(-1:0)  = shortwaveFluxMonthly(11:12)
    shortwaveFluxMonthlyExpand(13:14) = shortwaveFluxMonthly(1:2)

    longwaveFluxMonthlyExpand(1:12)  = longwaveFluxMonthly
    longwaveFluxMonthlyExpand(-1:0)  = longwaveFluxMonthly(11:12)
    longwaveFluxMonthlyExpand(13:14) = longwaveFluxMonthly(1:2)

    sensibleFluxMonthlyExpand(1:12)  = sensibleFluxMonthly
    sensibleFluxMonthlyExpand(-1:0)  = sensibleFluxMonthly(11:12)
    sensibleFluxMonthlyExpand(13:14) = sensibleFluxMonthly(1:2)

    latentFluxMonthlyExpand(1:12)  = latentFluxMonthly
    latentFluxMonthlyExpand(-1:0)  = latentFluxMonthly(11:12)
    latentFluxMonthlyExpand(13:14) = latentFluxMonthly(1:2)

    snowAlbedoMonthlyExpand(1:12)  = snowAlbedoMonthly
    snowAlbedoMonthlyExpand(-1:0)  = snowAlbedoMonthly(11:12)
    snowAlbedoMonthlyExpand(13:14) = snowAlbedoMonthly(1:2)

    ! convert units
    shortwaveFluxMonthlyExpand = shortwaveFluxMonthlyExpand * kcal_to_cal * cal_to_joule * m_to_cm**2
    longwaveFluxMonthlyExpand  = longwaveFluxMonthlyExpand  * kcal_to_cal * cal_to_joule * m_to_cm**2
    sensibleFluxMonthlyExpand  = sensibleFluxMonthlyExpand  * kcal_to_cal * cal_to_joule * m_to_cm**2
    latentFluxMonthlyExpand    = latentFluxMonthlyExpand    * kcal_to_cal * cal_to_joule * m_to_cm**2

    ! divide by month duration to get flux
    do nMonth = -1, 14

       nMonthUse = modulo(nMonth-1,12) + 1

       shortwaveFluxMonthlyExpand(nMonth) = shortwaveFluxMonthlyExpand(nMonth) / (days_in_month_real(nMonthUse) * seconds_in_day)
       longwaveFluxMonthlyExpand(nMonth)  = longwaveFluxMonthlyExpand(nMonth)  / (days_in_month_real(nMonthUse) * seconds_in_day)
       sensibleFluxMonthlyExpand(nMonth)  = sensibleFluxMonthlyExpand(nMonth)  / (days_in_month_real(nMonthUse) * seconds_in_day)
       latentFluxMonthlyExpand(nMonth)    = latentFluxMonthlyExpand(nMonth)    / (days_in_month_real(nMonthUse) * seconds_in_day)

    enddo ! nMonth

    ! interpolation coefficients
    call calculate_interpolation_coefficients(shortwaveFluxMonthlyIntCoeff, shortwaveFluxMonthlyExpand, semtner_data_times)
    call calculate_interpolation_coefficients(longwaveFluxMonthlyIntCoeff,  longwaveFluxMonthlyExpand,  semtner_data_times)
    call calculate_interpolation_coefficients(sensibleFluxMonthlyIntCoeff,  sensibleFluxMonthlyExpand,  semtner_data_times)
    call calculate_interpolation_coefficients(latentFluxMonthlyIntCoeff,    latentFluxMonthlyExpand,    semtner_data_times)
    call calculate_interpolation_coefficients(snowAlbedoMonthlyIntCoeff,    snowAlbedoMonthlyExpand,    semtner_data_times)

    ! init snowfall rate
    call init_semtner_snowfall()

    call plot_semtner_forcing()
    
  end subroutine init_semtner_forcing

  !-------------------------------------------------------------

  subroutine plot_semtner_forcing()

    integer :: &
         nMonth, &
         iday
    
    real(kind=RKIND) :: &
         shortwave_forcing, &
         longwave_forcing, &
         sensible_forcing, &
         latent_forcing, &
         snowalbedo_forcing, &
         snowfallRate, &
         time

    open(11,file="semtner/shortwaveFluxSemtnerData.txt")
    open(12,file="semtner/longwaveFluxSemtnerData.txt")
    open(13,file="semtner/sensibleFluxSemtnerData.txt")
    open(14,file="semtner/latentFluxSemtnerData.txt")
    open(15,file="semtner/snowAlbedoSemtnerData.txt")

    do nMonth = 1, 12

       write(11,*) (semtner_data_times(nMonth) * seconds_in_year) / seconds_in_day, shortwaveFluxMonthly(nMonth)
       write(12,*) (semtner_data_times(nMonth) * seconds_in_year) / seconds_in_day, longwaveFluxMonthly(nMonth)
       write(13,*) (semtner_data_times(nMonth) * seconds_in_year) / seconds_in_day, sensibleFluxMonthly(nMonth)
       write(14,*) (semtner_data_times(nMonth) * seconds_in_year) / seconds_in_day, latentFluxMonthly(nMonth)
       write(15,*) (semtner_data_times(nMonth) * seconds_in_year) / seconds_in_day, snowAlbedoMonthly(nMonth)
       
    enddo ! nMonth

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

    open(11,file="semtner/shortwaveFluxSemtner.txt")
    open(12,file="semtner/longwaveFluxSemtner.txt")
    open(13,file="semtner/sensibleFluxSemtner.txt")
    open(14,file="semtner/latentFluxSemtner.txt")
    open(15,file="semtner/snowAlbedoSemtner.txt")

    open(16,file="semtner/snowfallRateSemtner.txt")    

    do iday = 0, 360*24 

       time = real(iday,RKIND) * 3600.0_RKIND

       call interpolate_forcing(shortwave_forcing,  time, shortwaveFluxMonthlyIntCoeff, semtner_data_times)
       call interpolate_forcing(longwave_forcing,   time, longwaveFluxMonthlyIntCoeff,  semtner_data_times)
       call interpolate_forcing(sensible_forcing,   time, sensibleFluxMonthlyIntCoeff,  semtner_data_times)
       call interpolate_forcing(latent_forcing,     time, latentFluxMonthlyIntCoeff,    semtner_data_times)
       call interpolate_forcing(snowalbedo_forcing, time, snowAlbedoMonthlyIntCoeff,    semtner_data_times)

       write(11,*) real(iday,RKIND) / 24.0_RKIND, shortwave_forcing
       write(12,*) real(iday,RKIND) / 24.0_RKIND, longwave_forcing
       write(13,*) real(iday,RKIND) / 24.0_RKIND, sensible_forcing
       write(14,*) real(iday,RKIND) / 24.0_RKIND, latent_forcing
       write(15,*) real(iday,RKIND) / 24.0_RKIND, snowalbedo_forcing

       call semtner_snowfall(snowfallRate, time)

       write(16,*) real(iday,RKIND) / 24.0_RKIND, snowfallRate

    enddo ! iday

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)

  end subroutine plot_semtner_forcing
  
  !-------------------------------------------------------------

  subroutine calculate_interpolation_coefficients(MonthlyIntCoeff, MonthlyExpand, semtner_data_times)

    use cice_dynamics_shared, only: solve_linear_basis_system

    real(kind=RKIND), dimension(0:12,4), intent(out) :: &
         MonthlyIntCoeff

    real(kind=RKIND), dimension(-1:14), intent(inout) :: &
         MonthlyExpand

    real(kind=RKIND), dimension(-1:14), intent(in) :: &
         semtner_data_times

    real(kind=RKIND), dimension(4,4) :: &
         leftMatrix

    real(kind=RKIND) :: &
         test_forcing

    integer :: &
         iPeriod

    do iPeriod = 0, 12
       
       leftMatrix(1,:) = semtner_data_times(iPeriod-1)
       leftMatrix(2,:) = semtner_data_times(iPeriod)
       leftMatrix(3,:) = semtner_data_times(iPeriod+1)
       leftMatrix(4,:) = semtner_data_times(iPeriod+2)

       leftMatrix(:,1) = leftMatrix(:,1)**3
       leftMatrix(:,2) = leftMatrix(:,2)**2
       leftMatrix(:,3) = leftMatrix(:,3)
       leftMatrix(:,4) = 1.0_RKIND
 
       call solve_linear_basis_system(leftMatrix, MonthlyExpand(iPeriod-1:iPeriod+2), MonthlyIntCoeff(iPeriod,:))

    enddo ! iPeriod

  end subroutine calculate_interpolation_coefficients

  !-------------------------------------------------------------

  subroutine semtner_forcing( &
       shortwaveDown, &
       longwaveDown, &
       sensibleHeatFlux, &
       latentHeatFlux, &
       snowAlbedo, &
       time)
    
    real(kind=RKIND), intent(out) :: &
         shortwaveDown, &
         longwaveDown, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         snowAlbedo
    
    real(kind=RKIND), intent(in) :: &
         time

    real(kind=RKIND) :: &
         time_forcing

    time_forcing = mod(time,real(seconds_in_year,RKIND))

    call interpolate_forcing(shortwaveDown,    time_forcing, shortwaveFluxMonthlyIntCoeff, semtner_data_times)
    call interpolate_forcing(longwaveDown,     time_forcing, longwaveFluxMonthlyIntCoeff,  semtner_data_times)
    call interpolate_forcing(sensibleHeatFlux, time_forcing, sensibleFluxMonthlyIntCoeff,  semtner_data_times)
    call interpolate_forcing(latentHeatFlux,   time_forcing, latentFluxMonthlyIntCoeff,    semtner_data_times)
    call interpolate_forcing(snowAlbedo,       time_forcing, snowAlbedoMonthlyIntCoeff,    semtner_data_times)

  end subroutine semtner_forcing

  !-------------------------------------------------------------

  subroutine interpolate_forcing(forcing, time, MonthlyIntCoeff, semtner_data_times)

    real(kind=RKIND), intent(out) :: &
         forcing

    real(kind=RKIND), intent(in) :: &
         time

    real(kind=RKIND), dimension(0:12,4), intent(out) :: &
         MonthlyIntCoeff

    real(kind=RKIND), dimension(-1:14), intent(in) :: &
         semtner_data_times

    real(kind=RKIND) :: &
         time_year

    integer :: &
         iMonth, &
         iPeriod

    time_year = time / real(seconds_in_year)

    do iMonth = 0, 13

       if (time_year > semtner_data_times(iMonth) .and. time_year <= semtner_data_times(iMonth+1)) then

          iPeriod = iMonth
          exit

       endif

    enddo ! iMonth

    forcing = MonthlyIntCoeff(iPeriod, 1) * time_year**3 + &
              MonthlyIntCoeff(iPeriod, 2) * time_year**2 + &
              MonthlyIntCoeff(iPeriod, 3) * time_year    + &
              MonthlyIntCoeff(iPeriod, 4)

  end subroutine interpolate_forcing

  !-------------------------------------------------------------
  ! semtner snowfall 
  !-------------------------------------------------------------

  subroutine semtner_snowfall(snowfallRate, time)

    real(kind=RKIND), intent(out) :: &
         snowfallRate

    real(kind=RKIND), intent(in) :: &    
         time

    integer :: &
         iSnowPeriod

    real(kind=RKIND) :: &
         time_forcing

    time_forcing = mod(time, real(seconds_in_year,RKIND))

    snowfallRate = 0.0_RKIND

    do iSnowPeriod = 1, nSnowPeriods

       if (time_forcing > snowfallTimes(iSnowPeriod,1) .and. time_forcing <= snowfallTimes(iSnowPeriod,2)) then

          snowfallRate = snowfallAccumulationRate(iSnowPeriod)
          exit

       endif

    enddo ! iSnowPeriod

  end subroutine semtner_snowfall

  !-------------------------------------------------------------

  subroutine init_semtner_snowfall()

    integer :: &
         iSnowPeriod

    ! set times

    ! 20 August to 30 October
    snowfallTimes(1,1) = (cumulative_days_in_month(8) + 19) * seconds_in_day
    snowfallTimes(1,2) = (cumulative_days_in_month(9) + 30) * seconds_in_day

    ! 1 November to 31 December
    snowfallTimes(2,1) = (cumulative_days_in_month(11))      * seconds_in_day
    snowfallTimes(2,2) = (cumulative_days_in_month(12) + 30) * seconds_in_day

    ! 1 January to 30 April
    snowfallTimes(3,1) = 0.0_RKIND
    snowfallTimes(3,2) = (cumulative_days_in_month(4) + 30) * seconds_in_day

    ! 1 May to 31 May
    snowfallTimes(4,1) = cumulative_days_in_month(5) * seconds_in_day
    snowfallTimes(4,2) = cumulative_days_in_month(6) * seconds_in_day

    do iSnowPeriod = 1, nSnowPeriods
       snowfallDuration(iSnowPeriod) = snowfallTimes(iSnowPeriod,2) - snowfallTimes(iSnowPeriod,1)
    enddo ! iSnowPeriod

    ! set accumulation rate
    snowfallAccumulationRate(1) = (30.0_RKIND * cm_to_m) / snowfallDuration(1)

    snowfallAccumulationRate(2) = (5.0_RKIND * cm_to_m) / (snowfallDuration(2) + snowfallDuration(3))

    snowfallAccumulationRate(3) = (5.0_RKIND * cm_to_m) / (snowfallDuration(2) + snowfallDuration(3))

    snowfallAccumulationRate(4) = (5.0_RKIND * cm_to_m) / snowfallDuration(4) 

  end subroutine init_semtner_snowfall

  !-------------------------------------------------------------
  ! init state for semtner test case
  !-------------------------------------------------------------

  subroutine init_semtner_test_case(mesh, icestate)

    type(mesh_type), intent(in) :: mesh

    type(icestate_type), pointer :: icestate

    integer :: &
         iCell

    do iCell = 1, mesh % nCells

       icestate % iceAreaCell % array(iCell) = 1.0_RKIND
       
       icestate % surfaceTemperature % array(iCell) = -1.0_RKIND
       
       icestate % iceVolumeCell % array(iCell) = 0.1_RKIND
       icestate % snowVolumeCell % array(iCell) = 0.0_RKIND

    enddo ! iCell

  end subroutine init_semtner_test_case

  !-------------------------------------------------------------

end module cice_therm_zerolayer

