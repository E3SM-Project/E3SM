#define test_1d 0

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_therm_zerolayer
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_therm_zerolayer

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping

  implicit none

  private
  save

  public :: &
       cice_zerolayer_thermodynamics, &
       cice_zero_layer_thermo_semtner_unit_test

  logical, private :: ltest = .false.

  ! thermodynamic constants
  real(kind=RKIND), parameter, private :: &
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
       iceAlbedo                    = 0.5_RKIND, &!0.64_RKIND, &
       snowAlbedoConstant           = 0.7_RKIND, &!0.85_RKIND, &
       internalAbsorptionCorrection = 0.0_RKIND, &!0.4_RKIND, &
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

!-----------------------------------------------------------------------
! Semtner forcing
!-----------------------------------------------------------------------

    real(kind=RKIND), dimension(12), parameter, private :: shortwaveFluxMonthly = & ! kcal cm-2
         (/ 0.0_RKIND,  0.0_RKIND,  1.9_RKIND, &
            9.9_RKIND, 17.7_RKIND, 19.2_RKIND, &
           13.6_RKIND,  9.0_RKIND,  3.7_RKIND, &
            0.4_RKIND,  0.0_RKIND,  0.0_RKIND/)

    real(kind=RKIND), dimension(12), parameter, private :: longwaveFluxMonthly = & ! kcal cm-2
         (/10.4_RKIND, 10.3_RKIND, 10.3_RKIND, &
           11.6_RKIND, 15.1_RKIND, 18.0_RKIND, &
           19.1_RKIND, 18.7_RKIND, 16.5_RKIND, &
           13.9_RKIND, 11.2_RKIND, 10.9_RKIND/)
    
    real(kind=RKIND), dimension(12), parameter, private :: sensibleFluxMonthly = & ! kcal cm-2
         (/ 1.18_RKIND,  0.76_RKIND,  0.72_RKIND, &
            0.29_RKIND, -0.45_RKIND, -0.39_RKIND, &
           -0.30_RKIND, -0.40_RKIND, -0.17_RKIND, &
            0.10_RKIND,  0.56_RKIND,  0.79_RKIND/)

    real(kind=RKIND), dimension(12), parameter, private :: latentFluxMonthly = & ! kcal cm-2
         (/ 0.00_RKIND, -0.02_RKIND, -0.03_RKIND, &
           -0.09_RKIND, -0.46_RKIND, -0.70_RKIND, &
           -0.64_RKIND, -0.66_RKIND, -0.39_RKIND, &
           -0.19_RKIND, -0.01_RKIND, -0.01_RKIND/)

    real(kind=RKIND), dimension(12), parameter, private :: snowAlbedoMonthly = &
         (/0.83_RKIND, 0.83_RKIND, 0.83_RKIND, &
           0.81_RKIND, 0.82_RKIND, 0.78_RKIND, &
           0.64_RKIND, 0.69_RKIND, 0.84_RKIND, &
           0.85_RKIND, 0.85_RKIND, 0.85_RKIND/)

    real(kind=RKIND), dimension(-1:14), private :: &
         shortwaveFluxMonthlyExpand, &
         longwaveFluxMonthlyExpand, &
         sensibleFluxMonthlyExpand, &
         latentFluxMonthlyExpand, &
         snowAlbedoMonthlyExpand

    real(kind=RKIND), dimension(-1:14), private :: &
         semtner_data_times

    real(kind=RKIND), dimension(0:12,4), private :: &
         shortwaveFluxMonthlyIntCoeff, &
         longwaveFluxMonthlyIntCoeff, &
         sensibleFluxMonthlyIntCoeff, &
         latentFluxMonthlyIntCoeff, &
         snowAlbedoMonthlyIntCoeff

    integer, dimension(12), parameter, private :: days_in_month_real = &
         (/31,28,31,30,31,30,31,31,30,31,30,31/)

    integer, dimension(12), parameter, private :: days_in_month = &
         (/30,30,30,30,30,30,30,30,30,30,30,30/)


    integer, dimension(12), parameter, private :: cumulative_days_in_month = &
!         (/0,31,59,90,120,151,181,212,243,273,304,334/)
         (/0,30,60,90,120,150,180,210,240,270,300,330/)

    integer, parameter, private :: &
         seconds_in_day = 24.0_RKIND * 3600.0_RKIND, &
!         seconds_in_year = seconds_in_day * 365.0_RKIND
         seconds_in_year = seconds_in_day * 360.0_RKIND

    integer, parameter, private :: nSnowPeriods = 4

    real(kind=RKIND), dimension(nSnowPeriods,2), private :: &
         snowfallTimes

    real(kind=RKIND), dimension(nSnowPeriods), private :: &
         snowfallDuration, &
         snowfallAccumulationRate   

    integer, dimension(1), private :: iCell_test = (/-1/)

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_therm_zerolayer
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_zerolayer_thermodynamics(block, clock)!{{{

    use cice_constants, only: &
         iceAreaMinimum

    use cice_debug, only: &
         cice_time_output_variable_real

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type(MPAS_Clock_type), intent(in) :: &
         clock

    type (MPAS_pool_type), pointer :: &
         mesh, &
         configs, &
         tracers, &
         atmos_forcing, &
         ocean_coupling

    real(kind=RKIND) :: &
         cumulativeIceAblation, &
         cumulativeBasalAblation, &
         cumulativeBasalAccretion

    real(kind=RKIND), pointer :: &
         dt
         
    integer :: &
         iCell, &
         iCategory

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), dimension(:), pointer :: &
         shortwaveDown, &
         longwaveDown, &
         airDensity, &
         airSpecificHumidity, &
         airPotentialTemperature, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         snowfallRate, &
         seaSurfaceSalinity

#if test_1d == 1
    ! testing
    real(kind=RKIND), dimension(:), pointer :: &
         latCell
    ! testing
#endif

    configs => block % configs
    call MPAS_pool_get_config(configs, "config_dt", dt)
    
    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(atmos_forcing, "shortwaveDown", shortwaveDown)
    call MPAS_pool_get_array(atmos_forcing, "longwaveDown", longwaveDown)
    call MPAS_pool_get_array(atmos_forcing, "airDensity", airDensity)
    call MPAS_pool_get_array(atmos_forcing, "airSpecificHumidity", airSpecificHumidity)
    call MPAS_pool_get_array(atmos_forcing, "airPotentialTemperature", airPotentialTemperature)
    call MPAS_pool_get_array(atmos_forcing, "sensibleTransferCoefficient", sensibleTransferCoefficient)
    call MPAS_pool_get_array(atmos_forcing, "latentTransferCoefficient", latentTransferCoefficient)
    call MPAS_pool_get_array(atmos_forcing, "snowfallRate", snowfallRate)
    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)

#if test_1d == 1
    call MPAS_pool_get_array(mesh, "latCell", latCell) ! testing 

    ! testing
    if (iCell_test(1) == -1) then
       iCell_test = maxloc(latCell)
    endif
    ! testing
#endif

    do iCell = 1, nCellsSolve
       
       do iCategory = 1, nCategories

#if test_1d == 1
          ! testing
          if (iCell == iCell_test(1)) then
          ! testing
#endif   

          if (iceAreaCategory(1,iCategory,iCell) > iceAreaMinimum) then

             if (iceVolumeCategory(1,iCategory,iCell) <= 0.0_RKIND) then
                write(stderrUnit,*) "Thermodynamics: volume too small"
                write(stderrUnit,*) "iCell: ", iCell, " iCategory: ", iCategory
                write(stderrUnit,*) "iceAreaCategory:   ", iceAreaCategory(1,iCategory,iCell)
                write(stderrUnit,*) "iceVolumeCategory: ", iceVolumeCategory(1,iCategory,iCell)

                iceAreaCategory(1,iCategory,iCell) = 0.0_RKIND
                iceVolumeCategory(1,iCategory,iCell) = 0.0_RKIND
                snowVolumeCategory(1,iCategory,iCell) = 0.0_RKIND
                surfaceTemperature(1,iCategory,iCell) = 0.0_RKIND

                !stop
             endif
             
             call therm_zerolayer_column( &
                  iceAreaCategory(1,iCategory,iCell), &
                  iceVolumeCategory(1,iCategory,iCell), &
                  snowVolumeCategory(1,iCategory,iCell), &
                  surfaceTemperature(1,iCategory,iCell), &
                  shortwaveDown(iCell), &
                  longwaveDown(iCell), &
                  airDensity(iCell), &
                  airSpecificHumidity(iCell), &
                  airPotentialTemperature(iCell), &
                  sensibleTransferCoefficient(iCell), &
                  latentTransferCoefficient(iCell), &
                  snowfallRate(iCell), &
                  seaSurfaceSalinity(iCell), &
                  cumulativeIceAblation, &
                  cumulativeBasalAblation, &
                  cumulativeBasalAccretion, &
                  dt)
             
          else     

             call add_new_ice( &
                  iceAreaCategory(1,iCategory,iCell), &
                  iceVolumeCategory(1,iCategory,iCell), &
                  snowVolumeCategory(1,iCategory,iCell), &
                  surfaceTemperature(1,iCategory,iCell), &
                  airPotentialTemperature(iCell))

          endif ! iceAreaCell

#if test_1d == 1
          call writeout_zero_layer_thermo_point(&
               iCell, &
               mesh, &
               iceAreaCategory(1,iCategory,iCell), &
               iceVolumeCategory(1,iCategory,iCell), &
               snowVolumeCategory(1,iCategory,iCell), &
               surfaceTemperature(1,iCategory,iCell))
          
       endif ! testing 
#endif

       enddo ! iCategory

    enddo ! iCell

    !call cice_time_output_variable_real(mesh, iceAreaCategory(1,1,:), clock, "iceAreaCategory")
    !call cice_time_output_variable_real(mesh, iceVolumeCategory(1,1,:), clock, "iceVolumeCategory")

  end subroutine cice_zerolayer_thermodynamics!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  add_new_ice
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine add_new_ice( &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       surfaceTemperature, &
       airPotentialTemperature)!{{{

    real(kind=RKIND), intent(inout) :: &
         iceAreaCategory, &    !< Input/Output: 
         iceVolumeCategory, &  !< Input/Output: 
         snowVolumeCategory, & !< Input/Output: 
         surfaceTemperature    !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         airPotentialTemperature !< Input: 

    if (airPotentialTemperature < -2.0_RKIND) then

       iceAreaCategory = 1.0_RKIND
       iceVolumeCategory = 0.01_RKIND
       snowVolumeCategory = 0.0_RKIND
       surfaceTemperature = -2.0_RKIND

    endif

  end subroutine add_new_ice!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  therm_zerolayer_column
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine therm_zerolayer_column( &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
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
       dt)!{{{

    real(kind=RKIND), intent(inout) :: &
         iceAreaCategory, &         !< Input/Output: 
         iceVolumeCategory, &       !< Input/Output: 
         snowVolumeCategory, &      !< Input/Output: 
         surfaceTemperature, &      !< Input/Output: 
         cumulativeIceAblation, &   !< Input/Output: 
         cumulativeBasalAblation, & !< Input/Output: 
         cumulativeBasalAccretion   !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         shortwaveDownForcing, &        !< Input: 
         longwaveDownForcing, &         !< Input: 
         airDensity, &                  !< Input: 
         airSpecificHumidity, &         !< Input: 
         airPotentialTemperature, &     !< Input: 
         sensibleTransferCoefficient, & !< Input: 
         latentTransferCoefficient, &   !< Input: 
         snowfallRateForcing, &         !< Input: 
         seaSurfaceSalinity, &          !< Input: 
         dt                             !< Input: 

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
    iceThickness  = iceVolumeCategory  / iceAreaCategory
    snowThickness = snowVolumeCategory / iceAreaCategory

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
    if (iceThickness < 0.0_RKIND) then
       iceThickness = 0.0_RKIND
       snowThickness = 0.0_RKIND
       iceAreaCategory = 0.0_RKIND
       surfaceTemperature = 0.0_RKIND
    end if

    ! determine volumes
    iceVolumeCategory  = iceThickness  * iceAreaCategory
    snowVolumeCategory = snowThickness * iceAreaCategory

  end subroutine therm_zerolayer_column!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  snowfallAccumulation
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine snowfallAccumulation( &
       snowThickness, &
       snowfallRate, &
       surfaceTemperature, &
       dt)!{{{

    real(kind=RKIND), intent(inout) :: &
         snowThickness !< Input/Output: 

    real(kind=RKIND), intent(in) :: &    
         snowfallRate, &       !< Input: 
         surfaceTemperature, & !< Input: 
         dt                    !< Input: 

    real(kind=RKIND), parameter :: &
         densitySnow = 330.0_RKIND
    
    if (snowfallRate < 0.0_RKIND) then
       write(*,*) "Negative snowfall rate!",snowfallRate
       stop
    endif

    ! only accumulate if freezing surface
    if ((surfaceTemperature < snowMeltingTemperature .and. snowThickness > 0.0_RKIND) .or. &
        (surfaceTemperature < iceMeltingTemperature .and. snowThickness == 0.0_RKIND)) then

       snowThickness = snowThickness + (snowfallRate / densitySnow) * dt

    endif

  end subroutine snowfallAccumulation!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  basal_melting_accretion
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine basal_melting_accretion( &
       iceThickness, &
       snowThickness, &
       surfaceTemperature, &
       bottomTemperature, &
       oceanicFlux, &
       cumulativeBasalAblation, &
       cumulativeBasalAccretion, &
       dt)!{{{

    real(kind=RKIND), intent(inout) :: &
         iceThickness, &            !< Input/Output: 
         cumulativeBasalAblation, & !< Input/Output: 
         cumulativeBasalAccretion   !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         snowThickness, &      !< Input: 
         surfaceTemperature, & !< Input: 
         bottomTemperature, &  !< Input: 
         oceanicFlux, &        !< Input: 
         dt                    !< Input: 

    real(kind=RKIND) :: &
         conductivityCoefficient, &
         conductiveFlux, &
         growthEnergy, &
         iceThicknessChange

    if (iceThickness > 0.0_RKIND) then

       call conductivity_coefficient(conductivityCoefficient, iceThickness, snowThickness)
       
       conductiveFlux = conductivityCoefficient * (surfaceTemperature - bottomTemperature)
       
       growthEnergy = (oceanicFlux - conductiveFlux) * dt
       
       iceThicknessChange = max(growthEnergy / basalEnergyDensity, -iceThickness)
       
       iceThickness = iceThickness + iceThicknessChange
       
       cumulativeBasalAblation  = cumulativeBasalAblation  + abs(min(iceThicknessChange, 0.0_RKIND))
       cumulativeBasalAccretion = cumulativeBasalAccretion + max(iceThicknessChange, 0.0_RKIND)

    endif

  end subroutine basal_melting_accretion!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_melting
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

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
       dt)!{{{

    real(kind=RKIND), intent(inout) :: &
         newSurfaceTemperature, & !< Input/Output: 
         iceThickness, &          !< Input/Output: 
         snowThickness, &         !< Input/Output: 
         cumulativeIceAblation    !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         longwaveDown, &      !< Input: 
         sensibleHeatFlux, &  !< Input: 
         latentHeatFlux, &    !< Input: 
         shortwaveAbsorbed, & !< Input: 
         bottomTemperature, & !< Input: 
         dt                   !< Input: 

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
    if (newSurfaceTemperature >= snowMeltingTemperature .and. snowThickness > 0.0_RKIND) then
       newSurfaceTemperature = snowMeltingTemperature
       lmelt = .true.
    else if (newSurfaceTemperature >= iceMeltingTemperature .and. snowThickness == 0.0_RKIND) then
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

  end subroutine surface_melting!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculate_surface_temperature
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

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
       latentTransferCoefficient)!{{{

    real(kind=RKIND), intent(out) :: &
         newSurfaceTemperature !< Output: 
    
    real(kind=RKIND), intent(in) :: &
         surfaceTemperature, &          !< Input: 
         longwaveDown, &                !< Input: 
         shortwaveAbsorbed, &           !< Input: 
         bottomTemperature, &           !< Input: 
         iceThickness, &                !< Input: 
         snowThickness, &               !< Input: 
         airDensity, &                  !< Input: 
         airPotentialTemperature, &     !< Input: 
         airSpecificHumidity, &         !< Input: 
         sensibleTransferCoefficient, & !< Input: 
         latentTransferCoefficient      !< Input: 

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

  end subroutine calculate_surface_temperature!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_temperature_function
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function surface_temperature_function( &
       surfaceTemperature, &
       longwaveDown, &
       shortwaveAbsorbed, &
       sensibleHeatFlux, &
       latentHeatFlux, &
       conductivityCoefficient, &
       bottomTemperature) result(f)!{{{

    real(kind=RKIND), intent(in) :: &
         surfaceTemperature, &      !< Input: 
         longwaveDown, &            !< Input: 
         shortwaveAbsorbed, &       !< Input: 
         sensibleHeatFlux, &        !< Input: 
         latentHeatFlux, &          !< Input: 
         conductivityCoefficient, & !< Input: 
         bottomTemperature          !< Input: 

    real(kind=RKIND) :: f

    real(kind=RKIND) :: &
         surfaceTemperatureKelvin

    surfaceTemperatureKelvin = surfaceTemperature + freshWaterFreezingPoint

    f = longwaveDown + shortwaveAbsorbed + sensibleHeatFlux + latentHeatFlux - &
        stefanBoltzmann * surfaceTemperatureKelvin**4 - &
        conductivityCoefficient * (surfaceTemperature - bottomTemperature)

  end function surface_temperature_function!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_temperature_derivative
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------
  
  function surface_temperature_derivative( &
       surfaceTemperature, &
       dsensibleHeatFlux_dTsf, &
       dlatentHeatFlux_dTsf, &
       conductivityCoefficient) result(fp)!{{{

    real(kind=RKIND), intent(in) :: &
         surfaceTemperature, &     !< Input: 
         dsensibleHeatFlux_dTsf, & !< Input: 
         dlatentHeatFlux_dTsf, &   !< Input: 
         conductivityCoefficient   !< Input: 

    real(kind=RKIND) :: fp

    real(kind=RKIND) :: &
         surfaceTemperatureKelvin

    surfaceTemperatureKelvin = surfaceTemperature + freshWaterFreezingPoint

    fp = dsensibleHeatFlux_dTsf + dlatentHeatFlux_dTsf - &
         4.0_RKIND * stefanBoltzmann * surfaceTemperatureKelvin**3 - &
         conductivityCoefficient

  end function surface_temperature_derivative!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  sensible_latent_heat_fluxes
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------
  
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
       dlatentHeatFlux_dTsf)!{{{

    real(kind=RKIND), intent(in) :: &
         surfaceTemperature              !< Input: ice/snow surface temperature

    real(kind=RKIND), intent(in) :: &
         airDensity                  , & !< Input: air density (kg/m^3)
         airPotentialTemperature     , & !< Input: air potential temperature  (K)
         airSpecificHumidity         , & !< Input: specific humidity (kg/kg)
         sensibleTransferCoefficient , & !< Input: transfer coefficient for sensible heat
         latentTransferCoefficient       !< Input: transfer coefficient for latent heat

    real(kind=RKIND), intent(out) :: &
         sensibleHeatFlux            , & !< Output: surface downward sensible heat (W m-2)
         latentHeatFlux              , & !< Output: surface downward latent heat (W m-2)
         dsensibleHeatFlux_dTsf      , & !< Output: 
         dlatentHeatFlux_dTsf            !< Output: 

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
 
  end subroutine sensible_latent_heat_fluxes!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  absorbed_shortwave
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine absorbed_shortwave( &
       shortwaveAbsorbed, &
       shortwaveDown, &
       snowThickness, &
       snowAlbedo)!{{{

    real(kind=RKIND), intent(out) :: &
         shortwaveAbsorbed !< Output: 

    real(kind=RKIND), intent(in) :: &
         shortwaveDown, & !< Input: 
         snowThickness, & !< Input: 
         snowAlbedo       !< Input: 

    real(kind=RKIND) :: &
         effectiveAlbedo

    if (snowThickness > 0.0_RKIND) then

       shortwaveAbsorbed = shortwaveDown * (1.0_RKIND - snowAlbedo) 

    else

       shortwaveAbsorbed = shortwaveDown * (1.0_RKIND - iceAlbedo) * &
            ((1.0_RKIND - fractionPenetratingRadiation) + fractionPenetratingRadiation * (1.0_RKIND - internalAbsorptionCorrection))

    endif
    
  end subroutine absorbed_shortwave!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  conductivity_coefficient
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine conductivity_coefficient(&
       conductivityCoefficient, &
       iceThickness, &
       snowThickness)!{{{

    real(kind=RKIND), intent(out) :: &
         conductivityCoefficient !< Output: 

    real(kind=RKIND), intent(in) :: &    
         iceThickness, & !< Input: 
         snowThickness   !< Input: 

    if (snowThickness > 0.0_RKIND) then
       
       conductivityCoefficient = snowConductivityEffective / &
            (snowThickness + ((iceThickness * snowConductivityEffective) / iceConductivityEffective))
       
    else
       
       conductivityCoefficient = iceConductivityEffective / iceThickness
       
    endif

  end subroutine conductivity_coefficient!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculate_basal_properties
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calculate_basal_properties( &
       bottomTemperature, &
       oceanicFlux, &
       seaSurfaceSalinity)!{{{
    
    real(kind=RKIND), intent(out) :: &
         bottomTemperature, & !< Output: 
         oceanicFlux          !< Output: 
    
    real(kind=RKIND), intent(in) :: &
         seaSurfaceSalinity !< Input: 

    real(kind=RKIND), parameter :: &
         oceanicFluxSemtner = -1.5_RKIND ! kcal cm-2 year-1

    !bottomTemperature = -liquidusSlope * seaSurfaceSalinity
    bottomTemperature = -liquidusSlope * 34.0_RKIND

    oceanicFlux = (oceanicFluxSemtner * kcal_to_cal * cal_to_joule * m_to_cm**2) / seconds_in_year

  end subroutine calculate_basal_properties!}}}

!-----------------------------------------------------------------------
! Semtner forcing test case
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_zero_layer_thermo_semtner_unit_test
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_zero_layer_thermo_semtner_unit_test()!{{{

    call init_semtner_forcing()    

    call run_semtner_test_case()

  end subroutine cice_zero_layer_thermo_semtner_unit_test!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  run_semtner_test_case
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine run_semtner_test_case()!{{{

    real(kind=RKIND), parameter :: &
         dt = 3600.0_RKIND

    real(kind=RKIND) :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
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

    iceAreaCategory = 1.0_RKIND
    iceVolumeCategory = 1.0_RKIND
    snowVolumeCategory = 0.0_RKIND
    surfaceTemperature = -10.0_RKIND

    cumulativeIceAblation = 0.0_RKIND
    cumulativeBasalAblation = 0.0_RKIND
    cumulativeBasalAccretion = 0.0_RKIND

    open(11,file='therm_zerolayer.txt')

    open(12,file='therm_zerolayer_cum.txt')

    do istep = 1, nsteps

       time = istep * dt

       call therm_zerolayer_column( &
            iceAreaCategory, &
            iceVolumeCategory, &
            snowVolumeCategory, &
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
            dt)

       if (mod(istep,24) == 0) then
          write(11,*) time, iceVolumeCategory, snowVolumeCategory, surfaceTemperature
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

  end subroutine run_semtner_test_case!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_semtner_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_semtner_forcing()!{{{

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
    
  end subroutine init_semtner_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  plot_semtner_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine plot_semtner_forcing()!{{{

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

    open(11,file="shortwaveFluxSemtnerData.txt")
    open(12,file="longwaveFluxSemtnerData.txt")
    open(13,file="sensibleFluxSemtnerData.txt")
    open(14,file="latentFluxSemtnerData.txt")
    open(15,file="snowAlbedoSemtnerData.txt")

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

    open(11,file="shortwaveFluxSemtner.txt")
    open(12,file="longwaveFluxSemtner.txt")
    open(13,file="sensibleFluxSemtner.txt")
    open(14,file="latentFluxSemtner.txt")
    open(15,file="snowAlbedoSemtner.txt")

    open(16,file="snowfallRateSemtner.txt")    

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

  end subroutine plot_semtner_forcing!}}}
  
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculate_interpolation_coefficients
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calculate_interpolation_coefficients(&
       MonthlyIntCoeff, &
       MonthlyExpand, &
       semtner_data_times)!{{{

    use cice_numerics, only: &
         cice_solve_linear_basis_system

    real(kind=RKIND), dimension(0:12,4), intent(out) :: &
         MonthlyIntCoeff !< Output: 

    real(kind=RKIND), dimension(-1:14), intent(inout) :: &
         MonthlyExpand !< Input/Output: 

    real(kind=RKIND), dimension(-1:14), intent(in) :: &
         semtner_data_times !< Input: 

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
 
       call cice_solve_linear_basis_system(leftMatrix, MonthlyExpand(iPeriod-1:iPeriod+2), MonthlyIntCoeff(iPeriod,:))

    enddo ! iPeriod

  end subroutine calculate_interpolation_coefficients!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  semtner_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine semtner_forcing( &
       shortwaveDown, &
       longwaveDown, &
       sensibleHeatFlux, &
       latentHeatFlux, &
       snowAlbedo, &
       time)!{{{
    
    real(kind=RKIND), intent(out) :: &
         shortwaveDown, &    !< Output: 
         longwaveDown, &     !< Output: 
         sensibleHeatFlux, & !< Output: 
         latentHeatFlux, &   !< Output: 
         snowAlbedo          !< Output: 
    
    real(kind=RKIND), intent(in) :: &
         time !< Input: 

    real(kind=RKIND) :: &
         time_forcing

    time_forcing = mod(time,real(seconds_in_year,RKIND))

    call interpolate_forcing(shortwaveDown,    time_forcing, shortwaveFluxMonthlyIntCoeff, semtner_data_times)
    call interpolate_forcing(longwaveDown,     time_forcing, longwaveFluxMonthlyIntCoeff,  semtner_data_times)
    call interpolate_forcing(sensibleHeatFlux, time_forcing, sensibleFluxMonthlyIntCoeff,  semtner_data_times)
    call interpolate_forcing(latentHeatFlux,   time_forcing, latentFluxMonthlyIntCoeff,    semtner_data_times)
    call interpolate_forcing(snowAlbedo,       time_forcing, snowAlbedoMonthlyIntCoeff,    semtner_data_times)

  end subroutine semtner_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  interpolate_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine interpolate_forcing(&
       forcing, &
       time, &
       MonthlyIntCoeff, &
       semtner_data_times)!{{{

    real(kind=RKIND), intent(out) :: &
         forcing !< Output: 

    real(kind=RKIND), intent(in) :: &
         time !< Input: 
!nj-bgc+
    real(kind=RKIND), dimension(0:12,4), intent(in) :: &
         MonthlyIntCoeff !< Output: 
!nj-bgc-

    real(kind=RKIND), dimension(-1:14), intent(in) :: &
         semtner_data_times !< Input: 

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

  end subroutine interpolate_forcing!}}}

!-----------------------------------------------------------------------
! semtner snowfall 
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  semtner_snowfall
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine semtner_snowfall(&
       snowfallRate, &
       time)!{{{

    real(kind=RKIND), intent(out) :: &
         snowfallRate !< Output: 

    real(kind=RKIND), intent(in) :: &    
         time !< Input: 

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

  end subroutine semtner_snowfall!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_semtner_snowfall
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_semtner_snowfall()!{{{

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

  end subroutine init_semtner_snowfall!}}}

!-----------------------------------------------------------------------
! init state for semtner test case
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_semtner_test_case
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_semtner_test_case(&
       mesh, &
       tracers)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers !< Input/Output: 

    integer :: &
         iCell, &
         iCategory

    integer, pointer :: &
         nCells, &
         nCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature
    
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    do iCategory = 1, nCategories
       do iCell = 1, nCells
          
          iceAreaCategory(iCategory,iCell) = 1.0_RKIND
          
          surfaceTemperature(iCategory,iCell) = -1.0_RKIND
          
          iceVolumeCategory(iCategory,iCell)  = 0.1_RKIND
          snowVolumeCategory(iCategory,iCell) = 0.0_RKIND
          
       enddo ! iCell
    enddo ! iCategory

  end subroutine init_semtner_test_case!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  writeout_zero_layer_thermo_point
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine writeout_zero_layer_thermo_point(&
       iCell, &
       mesh, &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       surfaceTemperature)!{{{

    integer, intent(in) :: &
         iCell !< Input: 

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 
    
    real(kind=RKIND), intent(in) :: &
         iceAreaCategory, &    !< Input: 
         iceVolumeCategory, &  !< Input: 
         snowVolumeCategory, & !< Input: 
         surfaceTemperature    !< Input: 

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    integer, dimension(1) :: &
        max_loc 

    !call MPAS_pool_get_array(mesh, "latCell", latCell)

    !max_loc = maxloc(latCell)

    !if (iCell == max_loc(1)) then
    
       open(55,file="semtner_output.txt",position="append")

       write(55,*) iceAreaCategory, iceVolumeCategory, snowVolumeCategory, surfaceTemperature
       
       close(55)

    !endif

  end subroutine writeout_zero_layer_thermo_point!}}}

!-----------------------------------------------------------------------

end module cice_therm_zerolayer

