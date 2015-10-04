!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_forcing
!
!> \brief A core forcing module example 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  An example of a forcing module that might be implemented in a core
!>  using the MPAS_forcing module
!
!-----------------------------------------------------------------------

module cice_forcing

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_forcing
  use mpas_stream_manager

  implicit none

  private
  public :: &
       cice_forcing_init, &
       cice_forcing_get, &
       cice_forcing_write_restart, &
       cice_reset_coupler_fluxes

  type (MPAS_forcing_group_type), pointer :: ciceForcingGroups 

  ! forcing parameters
  real (kind=RKIND), parameter :: &
       fracShortwaveVisibleDirect  = 0.28_RKIND, & ! fraction of incoming shortwave in visible direct band
       fracShortwaveVisibleDiffuse = 0.24_RKIND, & ! fraction of incoming shortwave in visible diffuse band
       fracShortwaveIRDirectDown   = 0.31_RKIND, & ! fraction of incoming shortwave in near IR direct band
       fracShortwaveIRDiffuseDown  = 0.17_RKIND    ! fraction of incoming shortwave in near IR diffuse band

  ! precipitation factor
  real(kind=RKIND) :: &
       precipitationFactor

contains

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_forcing_init
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_forcing_init(domain)

    type (domain_type) :: domain

    ! init the atmospheric forcing
    call init_atmospheric_forcing(domain)

    ! init the ocean forcing
    call init_oceanic_forcing(domain)

  end subroutine cice_forcing_init

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing(domain)

    type (domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_atmospheric_forcing_type

    ! set the precipitation factor
    call init_precipitation_factor(domain)

    call MPAS_pool_get_config(domain % configs, "config_atmospheric_forcing_type", config_atmospheric_forcing_type)

    select case (trim(config_atmospheric_forcing_type))
    case ("CORE")
       call init_atmospheric_forcing_CORE(domain)
    case default
       write(stderrUnit,*) "Atmospheric forcing type unknown: ", trim(config_atmospheric_forcing_type)
       call MPAS_dmpar_global_abort('Unknown config_atmospheric_forcing_type')
    end select

  end subroutine init_atmospheric_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing_CORE
!
!> \brief Initialize the forcing objects
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  This routine calls the MPAS_forcing module subroutines that initializes
!>  the forcings type
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing_CORE(domain)

    type(domain_type) :: domain

    real(kind=RKIND), pointer :: &
         config_dt

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_forcing_restart_file

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalSixHourly, &
         forcingReferenceTimeSixHourly, &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    ! get atmospheric forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_forcing_start_time", config_forcing_start_time)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_start", config_forcing_cycle_start)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_duration", config_forcing_cycle_duration)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(domain % configs, "config_forcing_restart_file", config_forcing_restart_file)

    ! create the six hourly forcing group
    call MPAS_forcing_init_group(&
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         domain, &
         config_forcing_start_time, &
         "none", &
         config_forcing_cycle_duration, &
         config_do_restart, &
         config_forcing_restart_file, &
         .false.)

    forcingIntervalSixHourly = "06:00:00"
    forcingReferenceTimeSixHourly = "2000-01-01_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "airTemperature", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "airTemperature", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "airSpecificHumidity", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "airSpecificHumidity", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "uAirVelocity", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "uAirVelocity", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "vAirVelocity", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "vAirVelocity", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly, &
         "next")

    call MPAS_forcing_init_field_data(&
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         domain % streamManager)

    ! create the monthly forcing group
    call MPAS_forcing_init_group(&
         ciceForcingGroups, &
         "cice_atmospheric_forcing_monthly", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart, &
         config_forcing_restart_file)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "2000-01-15_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_monthly", &
         "cloudFraction", &
         "LYqMonthlyForcing", &
         "atmos_forcing", &
         "cloudFraction", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_monthly", &
         "rainfallRate", &
         "LYqMonthlyForcing", &
         "atmos_coupling", &
         "rainfallRate", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         ciceForcingGroups, &
         "cice_atmospheric_forcing_monthly", &
         domain % streamManager)

  end subroutine init_atmospheric_forcing_CORE

!-----------------------------------------------------------------------
! runtime
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_forcing
!
!> \brief Retrieve forcing data during time stepping
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  This routine calls the MPAS_forcing routine that will perform the 
!>  forcing data aquisition and interpolation during timestepping
!
!-----------------------------------------------------------------------

  subroutine cice_forcing_get(&
       streamManager, &
       domain, &
       simulationClock, &
       firstTimeStep)

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    logical, intent(in) :: &
         firstTimeStep

    call atmospheric_forcing(&
         streamManager, &
         domain, &
         simulationClock)

    call oceanic_forcing(&
         streamManager, &
         domain, &
         simulationClock, &
         firstTimeStep)

  end subroutine cice_forcing_get

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  atmospheric_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine atmospheric_forcing(&
       streamManager, &
       domain, &
       simulationClock)

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    type (block_type), pointer :: block

    real(kind=RKIND), pointer :: &
         config_dt

    character(len=strKIND), pointer :: &
         config_atmospheric_forcing_type

    ! configurations
    call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)
    call mpas_pool_get_config(domain % configs, 'config_atmospheric_forcing_type', config_atmospheric_forcing_type)

    ! use the forcing layer to get data
    if (trim(config_atmospheric_forcing_type) == "CORE") then

       call MPAS_forcing_get_forcing(&
            ciceForcingGroups, &
            "cice_atmospheric_forcing_sixhrly", &
            streamManager, &
            config_dt)

       call MPAS_forcing_get_forcing(&
            ciceForcingGroups, &
            "cice_atmospheric_forcing_monthly", &
            streamManager, &
            config_dt)

    endif

    block => domain % blocklist
    do while (associated(block))

       ! convert the input forcing variables to the coupling variables       
       select case (trim(config_atmospheric_forcing_type))
       case ("CORE")
          call prepare_atmospheric_coupling_variables_CORE(block)
       end select
       
       ! perform post coupling operations
       call post_atmospheric_coupling(block)
       
       ! perform post forcing
       call post_atmospheric_forcing(block)

       block => block % next
    end do
    
  end subroutine atmospheric_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_atmospheric_coupling_variables_CORE
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine prepare_atmospheric_coupling_variables_CORE(block)

    use cice_constants, only: &
         ciceFreshWaterFreezingPoint

    use mpas_constants, only: pii

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: &
         mesh, &
         atmosCoupling, &
         atmosForcing, &
         oceanCoupling, &
         tracers_aggregate
    
    real(kind=RKIND), dimension(:), pointer :: &
         airLevelHeight, &
         airPotentialTemperature, &
         airTemperature, &
         airSpecificHumidity, &
         airDensity, &
         shortwaveDown, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         longwaveDown, &
         rainfallRate, &
         snowfallRate, &
         cloudFraction, &
         seaSurfaceTemperature, &
         surfaceTemperatureCell, &
         lonCell, &
         latCell, &
         iceAreaCell

    type (MPAS_time_type) :: &
         currentForcingTime

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)    
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "lonCell", lonCell)    
    call MPAS_pool_get_array(mesh, "latCell", latCell)   

    call MPAS_pool_get_array(atmosCoupling, "airLevelHeight", airLevelHeight)
    call MPAS_pool_get_array(atmosCoupling, "airPotentialTemperature", airPotentialTemperature)
    call MPAS_pool_get_array(atmosCoupling, "airTemperature", airTemperature)
    call MPAS_pool_get_array(atmosCoupling, "airSpecificHumidity", airSpecificHumidity)
    call MPAS_pool_get_array(atmosCoupling, "airDensity", airDensity)
    call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
    call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
    call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
    call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)
    call MPAS_pool_get_array(atmosCoupling, "longwaveDown", longwaveDown)
    call MPAS_pool_get_array(atmosCoupling, "rainfallRate", rainfallRate)
    call MPAS_pool_get_array(atmosCoupling, "snowfallRate", snowfallRate)

    call MPAS_pool_get_array(atmosForcing, "cloudFraction", cloudFraction)
    call MPAS_pool_get_array(atmosForcing, "shortwaveDown", shortwaveDown)

    call MPAS_pool_get_array(oceanCoupling, "seaSurfaceTemperature", seaSurfaceTemperature)

    call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "surfaceTemperatureCell", surfaceTemperatureCell)

    ! get the current time
    call MPAS_forcing_get_forcing_time(&
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         currentForcingTime)

    do iCell = 1, nCells

       ! limit air temperature values where ice is present
       if (iceAreaCell(iCell) > 0.1_RKIND) then
          airTemperature(iCell) = min(airTemperature(iCell), ciceFreshWaterFreezingPoint + 0.1_RKIND)
       endif
       
       ! prevent supersaturated humidity
       call limit_specific_humidity(&
            airTemperature(iCell), &
            airSpecificHumidity(iCell))

    enddo ! iCell

    ! shortwave
    call shortwave_down(&
         shortwaveDown, &
         lonCell, &
         latCell, &
         cloudFraction, &
         airSpecificHumidity, &
         currentForcingTime)

    do iCell = 1, nCells

       shortwaveVisibleDirectDown(iCell)  = shortwaveDown(iCell) * fracShortwaveVisibleDirect
       shortwaveVisibleDiffuseDown(iCell) = shortwaveDown(iCell) * fracShortwaveVisibleDiffuse
       shortwaveIRDirectDown(iCell)       = shortwaveDown(iCell) * fracShortwaveIRDirectDown
       shortwaveIRDiffuseDown(iCell)      = shortwaveDown(iCell) * fracShortwaveIRDiffuseDown

       ! ensure physically realistic values
       cloudFraction(iCell)       = max(min(cloudFraction(iCell),1.0_RKIND),0.0_RKIND)
       shortwaveDown(iCell)       = max(shortwaveDown(iCell),0.0_RKIND)
       rainfallRate(iCell)        = max(rainfallRate(iCell),0.0_RKIND)
       airSpecificHumidity(iCell) = max(airSpecificHumidity(iCell),0.0_RKIND)

       ! atmospheric level height
       airLevelHeight(iCell) = 10.0_RKIND

       ! air potential temperature
       airPotentialTemperature(iCell) = airTemperature(iCell)

       ! air density
       airDensity(iCell) = 1.3_RKIND
       
    enddo ! iCell

    ! longwave radiation
    call longwave_rosati_miyakoda(&
         longwaveDown, &
         cloudFraction, &
         iceAreaCell, &
         surfaceTemperatureCell, &
         seaSurfaceTemperature, &
         airSpecificHumidity, &
         airTemperature)

    do iCell = 1, nCells

       ! precipitation
       call precipitation(&
            rainfallRate(iCell), &
            snowfallRate(iCell), &
            airTemperature(iCell))

    enddo ! iCell

  end subroutine prepare_atmospheric_coupling_variables_CORE

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  post_atmospheric_coupling
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine post_atmospheric_coupling(block)

    use cice_mesh, only: &
         cice_latlon_vector_rotation_forward

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: &
         mesh, &
         atmosCoupling, &
         atmosForcing

    real(kind=RKIND), dimension(:), pointer :: &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         shortwaveDown, &
         latCell, &
         lonCell, &
         xCell, &
         yCell, &
         zCell

    real(kind=RKIND), pointer :: &
         sphere_radius

    logical, pointer :: &
         config_rotate_cartesian_grid

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)    
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

    call MPAS_pool_get_config(block % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)

    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "lonCell", lonCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDirectDown", shortwaveVisibleDirectDown)
    call MPAS_pool_get_array(atmosCoupling, "shortwaveVisibleDiffuseDown", shortwaveVisibleDiffuseDown)
    call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDirectDown", shortwaveIRDirectDown)
    call MPAS_pool_get_array(atmosCoupling, "shortwaveIRDiffuseDown", shortwaveIRDiffuseDown)
    call MPAS_pool_get_array(atmosCoupling, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmosCoupling, "vAirVelocity", vAirVelocity)

    call MPAS_pool_get_array(atmosForcing, "windSpeed", windSpeed)
    call MPAS_pool_get_array(atmosForcing, "shortwaveDown", shortwaveDown)

    do iCell = 1, nCells

       ! wind speed
       windSpeed(iCell) = sqrt(uAirVelocity(iCell)**2 + vAirVelocity(iCell)**2)

       ! rotate velocities from geographical to local
       call cice_latlon_vector_rotation_forward(&
            uAirVelocity(iCell), &
            vAirVelocity(iCell), &
            uAirVelocity(iCell), &
            vAirVelocity(iCell), &
            latCell(iCell), &
            lonCell(iCell), &
            xCell(iCell), &
            yCell(iCell), &
            zCell(iCell), &
            sphere_radius, &
            config_rotate_cartesian_grid)

       ! shortwave - comment out for bit for bit tests
       !shortwaveDown(iCell) = &
       !     shortwaveVisibleDirectDown(iCell) + &
       !     shortwaveVisibleDiffuseDown(iCell) + &
       !     shortwaveIRDirectDown(iCell) + &
       !     shortwaveIRDiffuseDown(iCell)
       
    enddo ! iCell

  end subroutine post_atmospheric_coupling

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  post_atmospheric_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine post_atmospheric_forcing(block)

    use cice_constants, only: &
         ciceAirSpecificHeat, &
         ciceLatentHeatSublimation

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: &
         mesh, &
         atmosCoupling, &
         atmosForcing

    real(kind=RKIND), dimension(:), pointer :: &
         airDensity, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         uAirStress, &
         vAirStress

    real(kind=RKIND) :: &
         airStressCoefficient

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)    
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(atmosCoupling, "airDensity", airDensity)
    call MPAS_pool_get_array(atmosCoupling, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmosCoupling, "vAirVelocity", vAirVelocity)

    call MPAS_pool_get_array(atmosForcing, "windSpeed", windSpeed)
    call MPAS_pool_get_array(atmosForcing, "sensibleTransferCoefficient", sensibleTransferCoefficient)
    call MPAS_pool_get_array(atmosForcing, "latentTransferCoefficient", latentTransferCoefficient)
    call MPAS_pool_get_array(atmosForcing, "uAirStress", uAirStress)
    call MPAS_pool_get_array(atmosForcing, "vAirStress", vAirStress)

    do iCell = 1, nCells

       ! transfer coefficients
       sensibleTransferCoefficient(iCell) = 1.20e-3_RKIND * ciceAirSpecificHeat       * airDensity(iCell) * windSpeed(iCell)
       latentTransferCoefficient(iCell)   = 1.50e-3_RKIND * ciceLatentHeatSublimation * airDensity(iCell) * windSpeed(iCell)
       
       ! air stresses
       airStressCoefficient = 0.0012_RKIND * airDensity(iCell) * windSpeed(iCell)
       
       uAirStress(iCell) = uAirVelocity(iCell) * airStressCoefficient
       vAirStress(iCell) = vAirVelocity(iCell) * airStressCoefficient
       
    enddo ! iCell

  end subroutine post_atmospheric_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  limit_specific_humidity
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 17th Febuary 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine limit_specific_humidity(&
       airTemperature, &
       airSpecificHumidity)

    use cice_constants, only: &
         ciceFreshWaterFreezingPoint, &
         cicePuny

    real(kind=RKIND), intent(in) :: &
         airTemperature

    real(kind=RKIND), intent(inout) :: &
         airSpecificHumidity

    real(kind=RKIND) :: &
         airSpecificHumidityMax

    ! convert air temperature from Kelvin to Celcius
    airSpecificHumidityMax = airTemperature - ciceFreshWaterFreezingPoint

    airSpecificHumidityMax = 2.0_RKIND + &
         ((0.7859_RKIND + 0.03477_RKIND * airSpecificHumidityMax) / &
          (1.0_RKIND    + 0.00412_RKIND * airSpecificHumidityMax)) + &
          0.00422_RKIND * airSpecificHumidityMax

    ! vapor pressure
    airSpecificHumidityMax = 10.0_RKIND**airSpecificHumidityMax ! saturated 
    airSpecificHumidityMax = max(airSpecificHumidityMax,cicePuny) ! prevent division by zero

    ! specific humidity
    airSpecificHumidityMax = &
         (0.622_RKIND * airSpecificHumidityMax) / &
         (1.0e5_RKIND - 0.378_RKIND * airSpecificHumidityMax)

    ! limit the specific humidity
    airSpecificHumidity = min(airSpecificHumidity, airSpecificHumidityMax)

  end subroutine limit_specific_humidity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  longwave_rosati_miyakoda
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 17th Febuary 2015
!> \details
!>  based on 
!   Rosati, A. and K. Miyakoda (1988), 
!   A general-circulation model for upper ocean simulation, 
!   J. Physical Oceanography, 18, 1601-1626, 
!   doi:10.1175/1520-0485(1988)018<1601:AGCMFU>2.0.CO;2 
!
!-----------------------------------------------------------------------

  subroutine longwave_rosati_miyakoda(&
       longwaveDown, &
       cloudFraction, &
       iceAreaCell, &
       surfaceTemperature, &
       seaSurfaceTemperature, &
       airSpecificHumidity, &
       airTemperature)

    use cice_constants, only: &
         ciceStefanBoltzmann, &
         ciceFreshWaterFreezingPoint, &
         ciceIceSnowEmissivity

    real(kind=RKIND), dimension(:), intent(out) :: &
         longwaveDown

    real(kind=RKIND), dimension(:), intent(in) :: &
         cloudFraction, &
         iceAreaCell, &
         surfaceTemperature, &
         seaSurfaceTemperature, &
         airSpecificHumidity, &
         airTemperature

    real(kind=RKIND) :: &
         clearSkyFraction, &
         combinedSurfaceTemperature, &
         vapourPressureSqrt, &
         airPotentialTemperature, &
         airSeaTemperatureDifferenceTerm

    integer :: &
         iCell

    do iCell = 1, size(longwaveDown)

       ! get a clear sky fraction
       clearSkyFraction = 1.0_RKIND - 0.8_RKIND * cloudFraction(iCell)
    
       ! combined ice and ocean temperature in Kelvin
       combinedSurfaceTemperature = &
            surfaceTemperature(iCell)    * iceAreaCell(iCell)               + &
            seaSurfaceTemperature(iCell) * (1.0_RKIND - iceAreaCell(iCell)) + &
            ciceFreshWaterFreezingPoint

       ! square root of the vapour pressure
       vapourPressureSqrt = sqrt((1000.0_RKIND * airSpecificHumidity(iCell)) / &
                                 (0.622_RKIND + 0.378_RKIND * airSpecificHumidity(iCell)))

       ! potential temperature (CICE comment: get this from stability?)
       airPotentialTemperature = airTemperature(iCell)

       ! unknown
       !airSeaTemperatureDifferenceTerm = airPotentialTemperature**3 * &
       !     (airPotentialTemperature * (0.39_RKIND - 0.05_RKIND * vapourPressureSqrt) * clearSkyFraction + &
       !     4.0_RKIND * (combinedSurfaceTemperature - airPotentialTemperature))

       ! Different version for bfb to CICE
       airSeaTemperatureDifferenceTerm = airPotentialTemperature*airPotentialTemperature*airPotentialTemperature * &
            (airPotentialTemperature * (0.39_RKIND - 0.05_RKIND * vapourPressureSqrt) * clearSkyFraction + &
            4.0_RKIND * (combinedSurfaceTemperature - airPotentialTemperature))

       ! final longwave calculation from stefan-boltzmann law
       longwaveDown(iCell) = ciceIceSnowEmissivity * ciceStefanBoltzmann * &
                             (combinedSurfaceTemperature**4 - airSeaTemperatureDifferenceTerm)

    enddo ! iCell

  end subroutine longwave_rosati_miyakoda

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine longwave_parkinson_and_washington(&
       longWavedown, &
       airTemperature, &
       cloudFraction)

    use cice_constants, only: &
         ciceStefanBoltzmann, &
         ciceFreshWaterFreezingPoint

    real(kind=RKIND), intent(out) :: &
         longWavedown

    real(kind=RKIND), intent(in) :: &    
         airTemperature, &
         cloudFraction

    ! Longwave down
    ! Parkinson, C. L. and W. M. Washington (1979),
    ! Large-scale numerical-model of sea ice,
    ! JGR, 84, 311-337, doi:10.1029/JC084iC01p00311

    longWavedown = &
         ciceStefanBoltzmann * airTemperature**4 * &
         (1.0_RKIND - 0.261_RKIND * exp(-7.77e-4_RKIND * (ciceFreshWaterFreezingPoint - airTemperature)**2)) * &
         (1.0_RKIND + 0.275_RKIND * cloudFraction)

  end subroutine longwave_parkinson_and_washington

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  convert precipitation units to kg/m^2 s
!
!-----------------------------------------------------------------------

  subroutine init_precipitation_factor(domain)

    use cice_constants, only: &
         ciceSecondsPerYear, &
         ciceSecondsPerDay

    type(domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_precipitation_units

    call MPAS_pool_get_config(domain % configs, "config_forcing_precipitation_units", config_forcing_precipitation_units)

    select case (trim(config_forcing_precipitation_units))
    case ("mm_per_month")
       precipitationFactor = 12.0_RKIND / real(ciceSecondsPerYear,RKIND)
    case ("mm_per_day")
       precipitationFactor = 1.0_RKIND / real(ciceSecondsPerDay,RKIND)
    case ("mm_per_sec","mks")
       precipitationFactor = 1.0_RKIND
    case default
       write(stderrUnit,*) "Unknown config_precipitation_units"
       call MPAS_dmpar_global_abort('Unknown config_precipitation_units')
    end select

  end subroutine init_precipitation_factor

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine precipitation(&
       rainfallRate, &
       snowfallRate, &
       airTemperature)

    use cice_constants, only: &
         ciceSecondsPerYear, &
         ciceFreshWaterFreezingPoint

    real(kind=RKIND), intent(inout) :: &
         rainfallRate

    real(kind=RKIND), intent(out) :: &    
         snowfallRate
!nj-bgc+
    real(kind=RKIND), intent(in) :: & 
!nj-bgc-
         airTemperature
   
    rainfallRate = rainfallRate * precipitationFactor

    ! divide total precipitation between rain and snow
    snowfallRate = 0.0_RKIND

    if (airTemperature < ciceFreshWaterFreezingPoint) then

       snowfallRate = rainfallRate
       rainfallRate = 0.0_RKIND

    endif

  end subroutine precipitation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  shortwave_down
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine shortwave_down(&
       shortwaveDown, &
       longitudeIn, &
       latitude, &
       cloudFraction, &
       airSpecificHumidity, &
       currentTime)

    use cice_constants, only: &
         ciceDegreesToRadians, &
         ciceSecondsPerDay

    use mpas_constants, only: &
         pii

    real(kind=RKIND), dimension(:), intent(out) :: &
         shortwaveDown

    real(kind=RKIND), dimension(:), intent(in) :: &    
         longitudeIn, &
         latitude, &
         cloudFraction, &
         airSpecificHumidity

    type(MPAS_time_type), intent(in) :: &
         currentTime

    integer :: &
         dayOfYear, H, M, S, S_d, S_n

    real(kind=RKIND) :: &
         longitude, &
         secondsToday, &
         solarTime, &
         hourAngle, &
         declination, &
         cosZ, &
         e, &
         d, &
         sw0

    integer :: &
         iCell

    call mpas_get_time(currentTime, DoY=dayOfYear, H=H, M=M, S=S, S_n=S_n, S_d=S_d)

    secondsToday = real(H,RKIND)*3600.0_RKIND + &
                   real(M,RKIND)*60.0_RKIND + &
                   real(S,RKIND) + &
                   real(S_n,RKIND)/real(S_d,RKIND)

    do iCell = 1, size(shortwaveDown)

       ! longitude needs to be [-pi,pi] not [0,2pi]
       longitude = longitudeIn(iCell)
       if (longitude > pii) longitude = longitude - 2.0_RKIND * pii

       solarTime = mod(real(secondsToday,kind=RKIND),ciceSecondsPerDay)/3600.0_RKIND + 12.0_RKIND*sin(0.5_RKIND*longitude)

       hourAngle = (12.0_RKIND - solarTime)*pii/12.0_RKIND

       ! solar declinatiom
       declination = 23.44_RKIND*cos((172.0_RKIND-real(dayOfYear,RKIND)) * 2.0_RKIND*pii/365.0_RKIND)*ciceDegreesToRadians

       ! solar zenith angle
       cosZ = sin(latitude(iCell))*sin(declination) + cos(latitude(iCell))*cos(declination)*cos(hourAngle)
       cosZ = max(cosZ,0.0_RKIND)

       e = 1.0e5_RKIND*airSpecificHumidity(iCell)/(0.622_RKIND + 0.378_RKIND*airSpecificHumidity(iCell))

       d = (cosZ + 2.7_RKIND)*e*1.0e-5_RKIND+1.085_RKIND*cosZ+0.1_RKIND

       sw0 = 1353.0_RKIND*cosZ**2/d
       sw0 = max(sw0,0.0_RKIND)

       shortwaveDown(iCell) = sw0 * (1.0_RKIND - 0.6_RKIND * cloudFraction(iCell)**3) 

    enddo ! iCell

  end subroutine shortwave_down

!-----------------------------------------------------------------------
! ocean forcing
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_oceanic_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 17th Febuary 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_oceanic_forcing(domain)

    type (domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_sst_type

    call MPAS_pool_get_config(domain % configs, "config_forcing_sst_type", config_forcing_sst_type)

    select case (trim(config_forcing_sst_type))
    case ("ncar")
       call init_oceanic_forcing_ncar(domain)
    case default
       write(stderrUnit,*) "Oceanic forcing type unknown: ", trim(config_forcing_sst_type)
       call MPAS_dmpar_global_abort('Unknown config_forcing_sst_type')
    end select

  end subroutine init_oceanic_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_oceanic_forcing_ncar
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 17th Febuary 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_oceanic_forcing_ncar(domain)

    type (domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_restart_file

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    ! get oceanic forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(domain % configs, "config_forcing_restart_file", config_forcing_restart_file)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "2000-01-15_00:00:00"

    ! create the sea surface temperature forcing group
    call MPAS_forcing_init_group(&
         ciceForcingGroups, &
         "cice_sst_forcing_monthly", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart, &
         config_forcing_restart_file)

    ! sea surface temperature
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_sst_forcing_monthly", &
         "seaSurfaceTemperature", &
         "NCARMonthlySSTForcing", &
         "ocean_coupling", &
         "seaSurfaceTemperature", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         ciceForcingGroups, &
         "cice_sst_forcing_monthly", &
         domain % streamManager)

    ! create the monthly forcing group
    call MPAS_forcing_init_group(&
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart, &
         config_forcing_restart_file)

    ! sea surface salinity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         "seaSurfaceSalinity", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceSalinity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! u ocean velocity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         "uOceanVelocity", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "uOceanVelocity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! v ocean velocity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         "vOceanVelocity", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "vOceanVelocity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! u surface tilt
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         "seaSurfaceTiltU", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceTiltU", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! v surface tilt
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         "seaSurfaceTiltV", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceTiltV", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! ocean mixed layer depth
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         "oceanMixedLayerDepth", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "oceanMixedLayerDepth", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! deep ocean heat flux convergence
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         "oceanHeatFluxConvergence", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "oceanHeatFluxConvergence", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         ciceForcingGroups, &
         "cice_oceanic_forcing_monthly", &
         domain % streamManager)

  end subroutine init_oceanic_forcing_ncar

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  oceanic_forcing
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine oceanic_forcing(&
       streamManager, &
       domain, &
       simulationClock, &
       firstTimeStep)

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    logical, intent(in) :: &
         firstTimeStep

    type(block_type), pointer :: &
         block

    real(kind=RKIND), pointer :: &
         config_dt

    character(len=strKIND), pointer :: &
         config_forcing_sst_type

    logical, pointer :: &
         config_do_restart

    ! configurations
    call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)
    call mpas_pool_get_config(domain % configs, 'config_forcing_sst_type', config_forcing_sst_type)
    call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

    ! use the forcing layer to get data
    if (trim(config_forcing_sst_type) == 'ncar') then

       call MPAS_forcing_get_forcing(&
            ciceForcingGroups, &
            "cice_oceanic_forcing_monthly", &
            streamManager, &
            config_dt)

       ! only get sst data on first timestep
       if (firstTimeStep .and. .not. config_do_restart) then

          call MPAS_forcing_get_forcing(&
               ciceForcingGroups, &
               "cice_sst_forcing_monthly", &
               streamManager, &
               config_dt)

       endif
       
    endif

    block => domain % blocklist
    do while (associated(block))

       ! convert the input forcing variables to the coupling variables       
       select case (trim(config_forcing_sst_type))
       case ("ncar")
          call prepare_oceanic_coupling_variables_ncar(block, firstTimeStep)
       end select

       block => block % next
    end do

  end subroutine oceanic_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_oceanic_coupling_variables_ncar
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine prepare_oceanic_coupling_variables_ncar(block, firstTimeStep)

    use ice_colpkg, only: &
         colpkg_sea_freezing_temperature

    use cice_mesh, only: &
         cice_latlon_vector_rotation_forward

    type (block_type), pointer :: block

    logical, intent(in) :: &
         firstTimeStep

    type (mpas_pool_type), pointer :: &
         mesh, &
         ocean_coupling

    real(kind=RKIND), dimension(:), pointer :: &
         seaSurfaceTemperature, &
         seaSurfaceSalinity, &
         oceanMixedLayerDepth, &
         seaFreezingTemperature

    real(kind=RKIND), dimension(:), pointer :: &
         uOceanVelocity, &
         vOceanVelocity, &
         seaSurfaceTiltU, &
         seaSurfaceTiltV

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         lonCell, &
         xCell, &
         yCell, &
         zCell

    character(len=strKIND), pointer :: &
         config_sea_freezing_temperature_type

    logical, pointer :: &
         config_do_restart, &
         config_rotate_cartesian_grid

    real(kind=RKIND), pointer :: &
         sphere_radius

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_config(block % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(block % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)
    call MPAS_pool_get_config(block % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh) 
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)
    call MPAS_pool_get_array(ocean_coupling, "oceanMixedLayerDepth", oceanMixedLayerDepth)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

    do iCell = 1, nCells

       ! ensure physical realism
       seaSurfaceSalinity(iCell)   = max(seaSurfaceSalinity(iCell), 0.0_RKIND)
       oceanMixedLayerDepth(iCell) = max(oceanMixedLayerDepth(iCell), 0.0_RKIND)

       ! sea freezing temperature 
       seaFreezingTemperature(iCell) = colpkg_sea_freezing_temperature(seaSurfaceSalinity(iCell))

    enddo ! iCell

    ! only update sea surface temperature on first non-restart timestep
    if (firstTimeStep .and. .not. config_do_restart) then

       do iCell = 1, nCells

          ! sea surface temperature
          seaSurfaceTemperature(iCell) = max(seaSurfaceTemperature(iCell), seaFreezingTemperature(iCell))

       enddo ! iCell

    endif

    ! rotate currents from geographical
    call MPAS_pool_get_array(ocean_coupling, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "vOceanVelocity", vOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTiltU", seaSurfaceTiltU)
    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTiltV", seaSurfaceTiltV)

    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "lonCell", lonCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    do iCell = 1, nCells

       call cice_latlon_vector_rotation_forward(&
            uOceanVelocity(iCell), &
            vOceanVelocity(iCell), &
            uOceanVelocity(iCell), &
            vOceanVelocity(iCell), &
            latCell(iCell), &
            lonCell(iCell), &
            xCell(iCell), &
            yCell(iCell), &
            zCell(iCell), &
            sphere_radius, &
            config_rotate_cartesian_grid)

       call cice_latlon_vector_rotation_forward(&
            seaSurfaceTiltU(iCell), &
            seaSurfaceTiltV(iCell), &
            seaSurfaceTiltU(iCell), &
            seaSurfaceTiltV(iCell), &
            latCell(iCell), &
            lonCell(iCell), &
            xCell(iCell), &
            yCell(iCell), &
            zCell(iCell), &
            sphere_radius, &
            config_rotate_cartesian_grid)

    enddo ! iCell

  end subroutine prepare_oceanic_coupling_variables_ncar

!-----------------------------------------------------------------------
! coupler fluxes initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_reset_coupler_fluxes
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 18th March 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_reset_coupler_fluxes(domain)

    type(domain_type) :: domain

    call reset_atmospheric_coupler_fluxes(domain)
    call reset_ocean_coupler_fluxes(domain)

  end subroutine cice_reset_coupler_fluxes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  reset_atmospheric_coupler_fluxes
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 18th March 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine reset_atmospheric_coupler_fluxes(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         atmosFluxes, &
         shortwave, &
         velocitySolver, &
         atmosCoupling, &
!nj-bgc+
         biogeochemistry
!nj-bgc-

    real(kind=RKIND), dimension(:), pointer :: &
         sensibleHeatFlux, &
         latentHeatFlux, &
         absorbedShortwaveFlux, &
         longwaveUp, &
         evaporativeHeatFlux, &
         airStressCellU, &
         airStressCellV, &
         atmosReferenceSpeed10m, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m

    real(kind=RKIND), dimension(:,:), pointer :: &
!nj-bgc+
         atmosBioFluxes, &
         atmosBlackCarbonFlux, &
         atmosDustFlux
!nj-bgc-
         

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

       call MPAS_pool_get_array(atmosFluxes, "sensibleHeatFlux", sensibleHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "latentHeatFlux", latentHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "longwaveUp", longwaveUp)
       call MPAS_pool_get_array(atmosFluxes, "evaporativeHeatFlux", evaporativeHeatFlux)

       call MPAS_pool_get_array(velocitySolver, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolver, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(shortwave, "absorbedShortwaveFlux", absorbedShortwaveFlux)

       call mpas_pool_get_array(atmosCoupling, 'atmosReferenceSpeed10m', atmosReferenceSpeed10m)
       call mpas_pool_get_array(atmosCoupling, 'atmosReferenceTemperature2m', atmosReferenceTemperature2m) 
       call mpas_pool_get_array(atmosCoupling, 'atmosReferenceHumidity2m', atmosReferenceHumidity2m)

!nj-bgc+
       call MPAS_pool_get_array(biogeochemistry, "atmosBioFluxes", atmosBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
       call MPAS_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)
!nj-bgc-
     
       sensibleHeatFlux            = 0.0_RKIND
       latentHeatFlux              = 0.0_RKIND
       longwaveUp                  = 0.0_RKIND
       evaporativeHeatFlux         = 0.0_RKIND
       
       airStressCellU              = 0.0_RKIND
       airStressCellV              = 0.0_RKIND
       
       absorbedShortwaveFlux       = 0.0_RKIND
       
       atmosReferenceSpeed10m      = 0.0_RKIND
       atmosReferenceTemperature2m = 0.0_RKIND
       atmosReferenceHumidity2m    = 0.0_RKIND
!nj-bgc+
       atmosBioFluxes              = 0.0_RKIND
       atmosBlackCarbonFlux        = 0.0_RKIND
       atmosDustFlux               = 0.0_RKIND
!nj-bgc-

       block => block % next
    end do
    
  end subroutine reset_atmospheric_coupler_fluxes
    
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  reset_ocean_coupler_fluxes
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 18th March 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine reset_ocean_coupler_fluxes(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         oceanFluxes, &
!nj-bgc+
         biogeochemistry
!nj-bgc-

    real(kind=RKIND), dimension(:), pointer :: &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
!nj-bgc+
         oceanNitrateFlux, &
         oceanSilicateFlux, &
         oceanAmmoniumFlux, &
         oceanDMSFlux, &
         oceanDMSPpFlux, &
         oceanDMSPdFlux, &
         oceanHumicsFlux, &
         oceanDustIronFlux

    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanBioFluxes, &
         oceanAlgaeFlux, &
         oceanDOCFlux, &
         oceanDICFlux, &
         oceanDONFlux, &
         oceanParticulateIronFlux, &
         oceanDissolvedIronFlux
!nj-bgc-
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)
!nj-bgc+
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
!nj-bgc-

       call MPAS_pool_get_array(oceanFluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanShortwaveFlux", oceanShortwaveFlux)

       !call MPAS_pool_get_array(aerosols, "oceanAerosolFlux", oceanAerosolFlux)
!bh-bgc+
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeFlux", oceanAlgaeFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDOCFlux", oceanDOCFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDICFlux", oceanDICFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDONFlux", oceanDONFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanNitrateFlux", oceanNitrateFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanSilicateFlux", oceanSilicateFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumFlux", oceanAmmoniumFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPpFlux", oceanDMSPpFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSPdFlux", oceanDMSPdFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDMSFlux", oceanDMSFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanHumicsFlux", oceanHumicsFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronFlux",oceanParticulateIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronFlux",oceanDissolvedIronFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanDustIronFlux",oceanDustIronFlux)
!nj-bgc-

       oceanFreshWaterFlux      = 0.0_RKIND
       oceanSaltFlux            = 0.0_RKIND
       oceanHeatFlux            = 0.0_RKIND
       oceanShortwaveFlux       = 0.0_RKIND
!nj-bgc+
       oceanBioFluxes           = 0.0_RKIND  
       oceanDOCFlux             = 0.0_RKIND
       oceanDICFlux             = 0.0_RKIND
       oceanDONFlux             = 0.0_RKIND
       oceanNitrateFlux         = 0.0_RKIND
       oceanSilicateFlux        = 0.0_RKIND
       oceanAmmoniumFlux        = 0.0_RKIND
       oceanDMSFlux             = 0.0_RKIND
       oceanDMSPpFlux           = 0.0_RKIND
       oceanDMSPdFlux           = 0.0_RKIND
       oceanHumicsFlux          = 0.0_RKIND
       oceanParticulateIronFlux = 0.0_RKIND
       oceanDissolvedIronFlux   = 0.0_RKIND
       oceanDustIronFlux        = 0.0_RKIND
!nj-bgc-

       block => block % next
    end do
    
  end subroutine reset_ocean_coupler_fluxes
  
!-----------------------------------------------------------------------
! restart
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_forcing_write_restart(domain)

    type(domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_restart_file

    call MPAS_pool_get_config(domain % configs, "config_forcing_restart_file", config_forcing_restart_file)

    call MPAS_forcing_write_restart_times(ciceForcingGroups, config_forcing_restart_file)
    
  end subroutine cice_forcing_write_restart

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine forcing_testing(block)

    use cice_debug, only: cice_debug_output_array

    type(block_type) :: block

    type(MPAS_pool_type), pointer :: pool

    real(kind=RKIND), dimension(:), pointer :: array

    character(len=200), parameter :: &
         poolName  = "ocean_coupling", &
         fieldName = "seaSurfaceTemperature"

    write(stderrUnit,*) trim(poolName), " ", trim(fieldName)

    call MPAS_pool_get_subpool(block % structs, trim(poolName), pool)
    
    call MPAS_pool_get_array(pool, trim(fieldName), array)

    call cice_debug_output_array(array)

  end subroutine forcing_testing

!-----------------------------------------------------------------------

end module cice_forcing
