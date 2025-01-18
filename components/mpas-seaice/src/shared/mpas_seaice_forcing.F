!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_forcing
!
!> \brief A core forcing module example
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  An example of a forcing module that might be implemented in a core
!>  using the MPAS_forcing module
!
!-----------------------------------------------------------------------

module seaice_forcing

  use mpas_kind_types
  use mpas_timer
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_forcing
  use mpas_stream_manager
  use mpas_log, only: mpas_log_write

  implicit none

  private
  public :: &
       seaice_forcing_init, &
       seaice_forcing_get, &
       seaice_forcing_write_restart_times, &
       seaice_reset_coupler_fluxes,  &
       post_atmospheric_coupling,  &
       post_atmospheric_forcing,   &
       post_oceanic_coupling

  type (MPAS_forcing_group_type), pointer, public :: seaiceForcingGroups

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
!  seaice_forcing_init
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_forcing_init(domain, clock)

    type (domain_type) :: domain

    type (MPAS_Clock_type) :: clock

    logical, pointer :: &
         config_use_forcing, &
         config_use_data_icebergs, &
         config_use_column_biogeochemistry, &
         config_use_zaerosols, &
         config_use_atm_dust_file, &
         config_use_iron_solubility_file

    call MPAS_pool_get_config(domain % configs, "config_use_forcing", config_use_forcing)
    call MPAS_pool_get_config(domain % configs, "config_use_data_icebergs", config_use_data_icebergs)
    call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_atm_dust_file", config_use_atm_dust_file)
    call MPAS_pool_get_config(domain % configs, "config_use_iron_solubility_file", config_use_iron_solubility_file)

    if (config_use_forcing) then

       ! init the atmospheric forcing
       call init_atmospheric_forcing(domain)

       ! init the ocean forcing
       call init_oceanic_forcing(domain)

       ! init the ocean bgc forcing
       if (config_use_column_biogeochemistry) &
          call init_oceanic_bgc_forcing(domain)

    endif

    if (config_use_zaerosols .and. (config_use_iron_solubility_file .or. &
       config_use_atm_dust_file) )  &
       call init_atm_iron_bgc_forcing(domain, clock)

    ! init the data iceberg forcing
    if (config_use_data_icebergs) call init_data_iceberg_forcing(domain, clock)

  end subroutine seaice_forcing_init

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
    case ("ISPOL")
       call init_atmospheric_forcing_ISPOL(domain)
    case ("hourly_velocities")
       call init_atmospheric_forcing_hourly_velocities(domain)
    case default
       call mpas_log_write("Atmospheric forcing type unknown: "//trim(config_atmospheric_forcing_type), MPAS_LOG_CRIT)
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
         config_forcing_cycle_duration

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

    ! create the six hourly forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         domain, &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_do_restart, &
         .false.)

    forcingIntervalSixHourly = "06:00:00"
    forcingReferenceTimeSixHourly = "2000-01-01_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
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
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
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
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
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
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         "vAirVelocity", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "vAirVelocity", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly, &
         "next")

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

    ! create the monthly forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_monthly", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "0001-01-15_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_monthly", &
         "cloudFraction", &
         "LYqMonthlyForcing", &
         "atmos_forcing", &
         "cloudFraction", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_monthly", &
         "rainfallRate", &
         "LYqMonthlyForcing", &
         "atmos_coupling", &
         "rainfallRate", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_monthly", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_atmospheric_forcing_CORE

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing_ISPOL
!
!> \brief Initialize the forcing objects
!> \author Nicole Jeffery, LANL
!> \date Nov, 2022
!> \details
!> This routine calls the MPAS_forcing module subroutines that initializes
!> a 1D forcing dataset. The data corresponds to the Weddell Sea (67.9S, 54W)
!> from June 16, 2004 to December 31, 2004, and is from the 2004
!> Ice Station Polarstern (ISPOL) field experiment.
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing_ISPOL(domain)

    type(domain_type) :: domain

    real(kind=RKIND), pointer :: &
         config_dt

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalSixHourly, &
         forcingReferenceTimeSixHourly, &
         forcingIntervalDaily, &
         forcingReferenceTimeDaily

    ! get atmospheric forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_forcing_start_time", config_forcing_start_time)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_start", config_forcing_cycle_start)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_duration", config_forcing_cycle_duration)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    ! create the six hourly forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         domain, &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_do_restart, &
         .false.)

    forcingIntervalSixHourly = "06:00:00"
    forcingReferenceTimeSixHourly = "2004-01-01_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         "shortwaveDown", &
         "ISPOLSixHourlyForcing", &
         "atmos_forcing", &
         "shortwaveDown", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         "longwaveDown", &
         "ISPOLSixHourlyForcing", &
         "atmos_coupling", &
         "longwaveDown", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly, &
         "next")

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

    ! create the daily forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_daily", &
         domain, &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_do_restart)

    forcingIntervalDaily = "00-00-01_00:00:00"
    forcingReferenceTimeDaily = "2004-01-01_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_daily", &
         "airTemperature", &
         "ISPOLDailyForcing", &
         "atmos_coupling", &
         "airTemperature", &
         "linear", &
         forcingReferenceTimeDaily, &
         forcingIntervalDaily, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_daily", &
         "airSpecificHumidity", &
         "ISPOLDailyForcing", &
         "atmos_coupling", &
         "airSpecificHumidity", &
         "linear", &
         forcingReferenceTimeDaily, &
         forcingIntervalDaily, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_daily", &
         "uAirVelocity", &
         "ISPOLDailyForcing", &
         "atmos_coupling", &
         "uAirVelocity", &
         "linear", &
         forcingReferenceTimeDaily, &
         forcingIntervalDaily, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_daily", &
         "vAirVelocity", &
         "ISPOLDailyForcing", &
         "atmos_coupling", &
         "vAirVelocity", &
         "linear", &
         forcingReferenceTimeDaily, &
         forcingIntervalDaily, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_daily", &
         "rainfallRate", &
         "ISPOLDailyForcing", &
         "atmos_coupling", &
         "rainfallRate", &
         "linear", &
         forcingReferenceTimeDaily, &
         forcingIntervalDaily, &
         "next")

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_daily", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_atmospheric_forcing_ISPOL

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atmospheric_forcing_hourly_velocities
!
!> \brief Initialize the forcing objects
!> \author Adrian K. Turner, LANL
!> \date 18th November 2021
!> \details
!>  This routine calls the MPAS_forcing module subroutines that initializes
!>  the forcings type
!
!-----------------------------------------------------------------------

  subroutine init_atmospheric_forcing_hourly_velocities(domain)

    type(domain_type) :: domain

    real(kind=RKIND), pointer :: &
         config_dt

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalHourlyVelocities, &
         forcingReferenceTimeHourlyVelocities, &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    ! get atmospheric forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_forcing_start_time", config_forcing_start_time)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_start", config_forcing_cycle_start)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_duration", config_forcing_cycle_duration)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    ! create the six hourly forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_hourly_velocities", &
         domain, &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_do_restart, &
         .false.)

    forcingIntervalHourlyVelocities = "01:00:00"
    forcingReferenceTimeHourlyVelocities = "2000-01-01_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_hourly_velocities", &
         "uAirVelocity", &
         "airHourlyVelocitiesForcing", &
         "atmos_coupling", &
         "uAirVelocity", &
         "linear", &
         forcingReferenceTimeHourlyVelocities, &
         forcingIntervalHourlyVelocities, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_hourly_velocities", &
         "vAirVelocity", &
         "airHourlyVelocitiesForcing", &
         "atmos_coupling", &
         "vAirVelocity", &
         "linear", &
         forcingReferenceTimeHourlyVelocities, &
         forcingIntervalHourlyVelocities, &
         "next")

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_hourly_velocities", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_atmospheric_forcing_hourly_velocities

!-----------------------------------------------------------------------
! runtime
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_forcing
!
!> \brief Retrieve forcing data during time stepping
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  This routine calls the MPAS_forcing routine that will perform the
!>  forcing data aquisition and interpolation during timestepping
!
!-----------------------------------------------------------------------

  subroutine seaice_forcing_get(&
       streamManager, &
       domain, &
       simulationClock, &
       firstTimeStep)

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    logical, intent(in) :: &
         firstTimeStep

    logical, pointer :: &
         config_use_forcing, &
         config_use_data_icebergs, &
         config_use_zaerosols, &
         config_use_atm_dust_file, &
         config_use_iron_solubility_file

    call MPAS_pool_get_config(domain % configs, "config_use_forcing", config_use_forcing)
    call MPAS_pool_get_config(domain % configs, "config_use_data_icebergs", config_use_data_icebergs)
    call MPAS_pool_get_config(domain % configs, "config_use_zaerosols", config_use_zaerosols)
    call MPAS_pool_get_config(domain % configs, "config_use_atm_dust_file", config_use_atm_dust_file)
    call MPAS_pool_get_config(domain % configs, "config_use_iron_solubility_file", config_use_iron_solubility_file)

    if (config_use_forcing) then

       call atmospheric_forcing(&
            streamManager, &
            domain, &
            simulationClock)

       call oceanic_forcing(&
            streamManager, &
            domain, &
            simulationClock, &
            firstTimeStep)

    endif

    if (config_use_zaerosols .and. (config_use_iron_solubility_file .or. &
      config_use_atm_dust_file) ) &
       call atmospheric_aerosols_forcing(&
            streamManager, &
            domain, &
            simulationClock, &
            firstTimeStep)

    ! data iceberg forcing
    if (config_use_data_icebergs) then

       call data_iceberg_forcing(&
            streamManager, &
            domain)

    endif

  end subroutine seaice_forcing_get

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
            seaiceForcingGroups, &
            "seaice_atmospheric_forcing_sixhrly", &
            streamManager, &
            config_dt)

       call MPAS_forcing_get_forcing(&
            seaiceForcingGroups, &
            "seaice_atmospheric_forcing_monthly", &
            streamManager, &
            config_dt)

    else if (trim(config_atmospheric_forcing_type) == "hourly_velocities") then

       call MPAS_forcing_get_forcing(&
            seaiceForcingGroups, &
            "seaice_atmospheric_forcing_hourly_velocities", &
            streamManager, &
            config_dt)

    else if (trim(config_atmospheric_forcing_type) == "ISPOL") then

       call MPAS_forcing_get_forcing(&
            seaiceForcingGroups, &
            "seaice_atmospheric_forcing_sixhrly", &
            streamManager, &
            config_dt)

       call MPAS_forcing_get_forcing(&
            seaiceForcingGroups, &
            "seaice_atmospheric_forcing_daily", &
            streamManager, &
            config_dt)

    endif

    block => domain % blocklist
    do while (associated(block))

       ! convert the input forcing variables to the coupling variables
       select case (trim(config_atmospheric_forcing_type))
       case ("CORE")
          call prepare_atmospheric_coupling_variables_CORE(block)
       case ("ISPOL")
          call prepare_atmospheric_coupling_variables_ISPOL(block)
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
!  atmospheric_aerosols_forcing
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date
!> \details: reads data file for dust fluxes and dust-iron solubility
!>
!
!-----------------------------------------------------------------------

  subroutine atmospheric_aerosols_forcing(&
       streamManager, &
       domain, &
       simulationClock, &
       firstTimeStep)

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    logical, intent(in) :: &
         firstTimeStep

    type (block_type), pointer :: block

    real(kind=RKIND), pointer :: &
         config_dt

    ! configurations
    call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

    call MPAS_forcing_get_forcing(seaiceForcingGroups, &
        'seaice_atm_bgc_forcing_monthly', streamManager, config_dt)

    block => domain % blocklist
    do while (associated(block))

       ! convert the input forcing variables to the coupling variables
       call prepare_atmospheric_dust_coupling_variables(block)

       block => block % next
    end do

  end subroutine atmospheric_aerosols_forcing

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

    use seaice_constants, only: &
         seaiceFreshWaterFreezingPoint, &
         pii

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

    logical, pointer :: &
         config_limit_air_temperatures

    type (MPAS_time_type) :: &
         currentForcingTime

    real(kind=RKIND) :: &
         secondsToday

    integer, pointer :: &
         nCellsSolve

    integer :: &
         dayOfYear, &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_config(block % configs, "config_limit_air_temperatures", config_limit_air_temperatures)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

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
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         currentForcingTime)

    ! get the number of seconds so far today
    call get_seconds_today(&
         currentForcingTime, &
         secondsToday, &
         dayOfYear)

    do iCell = 1, nCellsSolve

       ! limit air temperature values where ice is present
       if (iceAreaCell(iCell) > 0.1_RKIND .and. config_limit_air_temperatures) then
          airTemperature(iCell) = min(airTemperature(iCell), seaiceFreshWaterFreezingPoint + 0.1_RKIND)
       endif

       ! prevent supersaturated humidity
       call limit_specific_humidity(&
            airTemperature(iCell), &
            airSpecificHumidity(iCell))

       ! shortwave
       call shortwave_down(&
            shortwaveDown(iCell), &
            lonCell(iCell), &
            latCell(iCell), &
            cloudFraction(iCell), &
            airSpecificHumidity(iCell), &
            secondsToday, &
            dayOfYear)

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

       ! longwave radiation
       call longwave_rosati_miyakoda(&
            longwaveDown(iCell), &
            cloudFraction(iCell), &
            iceAreaCell(iCell), &
            surfaceTemperatureCell(iCell), &
            seaSurfaceTemperature(iCell), &
            airSpecificHumidity(iCell), &
            airTemperature(iCell))

       ! precipitation
       call precipitation(&
            rainfallRate(iCell), &
            snowfallRate(iCell), &
            airTemperature(iCell))

    enddo ! iCell

  end subroutine prepare_atmospheric_coupling_variables_CORE

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_atmospheric_coupling_variables_ISPOL
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date
!> \details
!> Defines the cell fields corresponding to the 1D forcing fields.
!
!-----------------------------------------------------------------------

  subroutine prepare_atmospheric_coupling_variables_ISPOL(block)

    use seaice_constants, only: &
         seaiceFreshWaterFreezingPoint, &
         pii

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
         rainfallRate, &
         snowfallRate, &
         cloudFraction, &
         lonCell, &
         latCell, &
         iceAreaCell

    type (MPAS_time_type) :: &
         currentForcingTime

    real(kind=RKIND) :: &
         secondsToday

    integer, pointer :: &
         nCellsSolve

    integer :: &
         dayOfYear, &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

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
    call MPAS_pool_get_array(atmosCoupling, "rainfallRate", rainfallRate)
    call MPAS_pool_get_array(atmosCoupling, "snowfallRate", snowfallRate)

    call MPAS_pool_get_array(atmosForcing, "cloudFraction", cloudFraction)
    call MPAS_pool_get_array(atmosForcing, "shortwaveDown", shortwaveDown)

    call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)

    ! get the current time
    call MPAS_forcing_get_forcing_time(&
         seaiceForcingGroups, &
         "seaice_atmospheric_forcing_sixhrly", &
         currentForcingTime)

    ! get the number of seconds so far today
    call get_seconds_today(&
         currentForcingTime, &
         secondsToday, &
         dayOfYear)

    do iCell = 1, nCellsSolve

       ! not sure if we need this
       cloudFraction(iCell)       = 0.5_RKIND

       ! limit air temperature values where ice is present
       if (iceAreaCell(iCell) > 0.1_RKIND) then
          airTemperature(iCell) = min(airTemperature(iCell), seaiceFreshWaterFreezingPoint + 0.1_RKIND)
       endif

       ! prevent supersaturated humidity
       call limit_specific_humidity(&
            airTemperature(iCell), &
            airSpecificHumidity(iCell))

       shortwaveVisibleDirectDown(iCell)  = shortwaveDown(iCell) * fracShortwaveVisibleDirect
       shortwaveVisibleDiffuseDown(iCell) = shortwaveDown(iCell) * fracShortwaveVisibleDiffuse
       shortwaveIRDirectDown(iCell)       = shortwaveDown(iCell) * fracShortwaveIRDirectDown
       shortwaveIRDiffuseDown(iCell)      = shortwaveDown(iCell) * fracShortwaveIRDiffuseDown

       ! ensure physically realistic values
       shortwaveDown(iCell)       = max(shortwaveDown(iCell),0.0_RKIND)
       rainfallRate(iCell)        = max(rainfallRate(iCell),0.0_RKIND)
       airSpecificHumidity(iCell) = max(airSpecificHumidity(iCell),0.0_RKIND)

       ! atmospheric level height
       airLevelHeight(iCell) = 10.0_RKIND ! for velocities, 2-m for air T and specific humidity

       ! air potential temperature
       airPotentialTemperature(iCell) = airTemperature(iCell)

       ! air density
       airDensity(iCell) = 1.3_RKIND

       ! precipitation

       ! divide total precipitation between rain and snow
       snowfallRate(iCell) = 0.0_RKIND

       if (airTemperature(iCell) < seaiceFreshWaterFreezingPoint) then

          snowfallRate(iCell) = rainfallRate(iCell)
          rainfallRate(iCell) = 0.0_RKIND

       endif

    enddo ! iCell

  end subroutine prepare_atmospheric_coupling_variables_ISPOL

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_atmospheric_dust_coupling_variables
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date
!> \details
!> Defines the atm flux fields for dust if reading from a file
!
!-----------------------------------------------------------------------

  subroutine prepare_atmospheric_dust_coupling_variables(block)

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: &
         mesh, &
         atmosCoupling, &
         atmosForcing, &
         biogeochemistry

    real(kind=RKIND), dimension(:,:), pointer :: &
         atmosDustFlux, &
         atmosWetDustFlux, &
         atmosDryDustFlux

    real(kind=RKIND), dimension(:), pointer :: &
         dust_FLUZ_WET, &
         dust_FLUZ_DRY

    integer, pointer :: &
         nCellsSolve, &
         maxBCType, &
         nzAerosols

    integer :: &
         iCell, &
         iBioTracers, &
         nDust

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

    ! zaerosols

    call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
    call MPAS_pool_get_dimension(mesh, "maxBCType", maxBCType)
    call MPAS_pool_get_dimension(mesh, "nzAerosols", nzAerosols)
    call MPAS_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)
    call MPAS_pool_get_array(biogeochemistry, "atmosWetDustFlux", atmosWetDustFlux)
    call MPAS_pool_get_array(biogeochemistry, "atmosDryDustFlux", atmosDryDustFlux)
    call MPAS_pool_get_array(atmosForcing, "dust_FLUZ_WET", dust_FLUZ_WET)
    call MPAS_pool_get_array(atmosForcing, "dust_FLUZ_DRY", dust_FLUZ_DRY)

    nDust = MAX(1,nzAerosols - maxBCType)
    do iCell = 1, nCellsSolve
       do iBioTracers = maxBCType + 1, nzAerosols
          atmosDustFlux(iBioTracers,iCell) = (dust_FLUZ_WET(iCell) + dust_FLUZ_DRY(iCell))/nDust
          atmosWetDustFlux(iBioTracers,iCell) = dust_FLUZ_WET(iCell)/nDust
          atmosDryDustFlux(iBiotracers,iCell) = dust_FLUZ_DRY(iCell)/nDust
      enddo
    end do

  end subroutine prepare_atmospheric_dust_coupling_variables

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

    use seaice_mesh, only: &
         seaice_latlon_vector_rotation_forward

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
         config_rotate_cartesian_grid, &
         on_a_sphere

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

    call MPAS_pool_get_config(block % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
    call MPAS_pool_get_config(mesh, "on_a_sphere", on_a_sphere)

    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

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

    do iCell = 1, nCellsSolve

       ! wind speed
       windSpeed(iCell) = sqrt(uAirVelocity(iCell)**2 + vAirVelocity(iCell)**2)

    enddo ! iCell

    if (on_a_sphere .and. config_rotate_cartesian_grid) then

       do iCell = 1, nCellsSolve

          ! rotate velocities from geographical to local
          call seaice_latlon_vector_rotation_forward(&
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

       enddo ! iCell

    endif !

    do iCell = 1, nCellsSolve

       ! shortwave - comment out for bit for bit tests
       shortwaveDown(iCell) = &
            shortwaveVisibleDirectDown(iCell) + &
            shortwaveVisibleDiffuseDown(iCell) + &
            shortwaveIRDirectDown(iCell) + &
            shortwaveIRDiffuseDown(iCell)

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

    use seaice_constants, only: &
         seaiceAirSpecificHeat, &
         seaiceLatentHeatSublimation

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
         nCellsSolve

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

    call MPAS_pool_get_array(atmosCoupling, "airDensity", airDensity)
    call MPAS_pool_get_array(atmosCoupling, "uAirVelocity", uAirVelocity)
    call MPAS_pool_get_array(atmosCoupling, "vAirVelocity", vAirVelocity)

    call MPAS_pool_get_array(atmosForcing, "windSpeed", windSpeed)
    call MPAS_pool_get_array(atmosForcing, "sensibleTransferCoefficient", sensibleTransferCoefficient)
    call MPAS_pool_get_array(atmosForcing, "latentTransferCoefficient", latentTransferCoefficient)
    call MPAS_pool_get_array(atmosForcing, "uAirStress", uAirStress)
    call MPAS_pool_get_array(atmosForcing, "vAirStress", vAirStress)

    do iCell = 1, nCellsSolve

       ! transfer coefficients
       sensibleTransferCoefficient(iCell) = 1.20e-3_RKIND * seaiceAirSpecificHeat       * airDensity(iCell) * windSpeed(iCell)
       latentTransferCoefficient(iCell)   = 1.50e-3_RKIND * seaiceLatentHeatSublimation * airDensity(iCell) * windSpeed(iCell)

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

    use seaice_constants, only: &
         seaiceFreshWaterFreezingPoint, &
         seaicePuny

    real(kind=RKIND), intent(in) :: &
         airTemperature

    real(kind=RKIND), intent(inout) :: &
         airSpecificHumidity

    real(kind=RKIND) :: &
         airSpecificHumidityMax

    ! convert air temperature from Kelvin to Celcius
    airSpecificHumidityMax = airTemperature - seaiceFreshWaterFreezingPoint

    airSpecificHumidityMax = 2.0_RKIND + &
         ((0.7859_RKIND + 0.03477_RKIND * airSpecificHumidityMax) / &
          (1.0_RKIND    + 0.00412_RKIND * airSpecificHumidityMax)) + &
          0.00422_RKIND * airSpecificHumidityMax

    ! vapor pressure
    airSpecificHumidityMax = 10.0_RKIND**airSpecificHumidityMax ! saturated
    airSpecificHumidityMax = max(airSpecificHumidityMax,seaicePuny) ! prevent division by zero

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

    use seaice_constants, only: &
         seaiceStefanBoltzmann, &
         seaiceFreshWaterFreezingPoint, &
         seaiceIceSnowEmissivity

    real(kind=RKIND), intent(out) :: &
         longwaveDown

    real(kind=RKIND), intent(in) :: &
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

    ! get a clear sky fraction
    clearSkyFraction = 1.0_RKIND - 0.8_RKIND * cloudFraction

    ! combined ice and ocean temperature in Kelvin
    combinedSurfaceTemperature = &
         surfaceTemperature    * iceAreaCell               + &
         seaSurfaceTemperature * (1.0_RKIND - iceAreaCell) + &
         seaiceFreshWaterFreezingPoint

    ! square root of the vapour pressure
    vapourPressureSqrt = sqrt((1000.0_RKIND * airSpecificHumidity) / &
                              (0.622_RKIND + 0.378_RKIND * airSpecificHumidity))

    ! potential temperature (CICE comment: get this from stability?)
    airPotentialTemperature = airTemperature

    ! unknown
    !airSeaTemperatureDifferenceTerm = airPotentialTemperature**3 * &
    !     (airPotentialTemperature * (0.39_RKIND - 0.05_RKIND * vapourPressureSqrt) * clearSkyFraction + &
    !     4.0_RKIND * (combinedSurfaceTemperature - airPotentialTemperature))

    ! Different version for bfb to CICE
    airSeaTemperatureDifferenceTerm = airPotentialTemperature*airPotentialTemperature*airPotentialTemperature * &
         (airPotentialTemperature * (0.39_RKIND - 0.05_RKIND * vapourPressureSqrt) * clearSkyFraction + &
         4.0_RKIND * (combinedSurfaceTemperature - airPotentialTemperature))

    ! final longwave calculation from stefan-boltzmann law
    longwaveDown = seaiceIceSnowEmissivity * seaiceStefanBoltzmann * &
                          (combinedSurfaceTemperature**4 - airSeaTemperatureDifferenceTerm)

  end subroutine longwave_rosati_miyakoda

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  longwave_parkinson_and_washington
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

    use seaice_constants, only: &
         seaiceStefanBoltzmann, &
         seaiceFreshWaterFreezingPoint

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
         seaiceStefanBoltzmann * airTemperature**4 * &
         (1.0_RKIND - 0.261_RKIND * exp(-7.77e-4_RKIND * (seaiceFreshWaterFreezingPoint - airTemperature)**2)) * &
         (1.0_RKIND + 0.275_RKIND * cloudFraction)

  end subroutine longwave_parkinson_and_washington

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_precipitation_factor
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>  convert precipitation units to kg/m^2 s
!
!-----------------------------------------------------------------------

  subroutine init_precipitation_factor(domain)

    use seaice_constants, only: &
         seaiceSecondsPerYear, &
         seaiceSecondsPerDay

    type(domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_precipitation_units

    call MPAS_pool_get_config(domain % configs, "config_forcing_precipitation_units", config_forcing_precipitation_units)

    select case (trim(config_forcing_precipitation_units))
    case ("mm_per_month")
       precipitationFactor = 12.0_RKIND / real(seaiceSecondsPerYear,RKIND)
    case ("mm_per_day")
       precipitationFactor = 1.0_RKIND / real(seaiceSecondsPerDay,RKIND)
    case ("mm_per_sec","mks")
       precipitationFactor = 1.0_RKIND
    case default
       call mpas_log_write("Unknown config_precipitation_units: "//trim(config_forcing_precipitation_units), MPAS_LOG_CRIT)
    end select

  end subroutine init_precipitation_factor

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  precipitation
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

    use seaice_constants, only: &
         seaiceFreshWaterFreezingPoint

    real(kind=RKIND), intent(inout) :: &
         rainfallRate

    real(kind=RKIND), intent(out) :: &
         snowfallRate

    real(kind=RKIND), intent(in) :: &
         airTemperature

    rainfallRate = rainfallRate * precipitationFactor

    ! divide total precipitation between rain and snow
    snowfallRate = 0.0_RKIND

    if (airTemperature < seaiceFreshWaterFreezingPoint) then

       snowfallRate = rainfallRate
       rainfallRate = 0.0_RKIND

    endif

  end subroutine precipitation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_seconds_today
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine get_seconds_today(&
       currentTime, &
       secondsToday, &
       dayOfYear)

    type(MPAS_time_type), intent(in) :: &
         currentTime

    real(kind=RKIND), intent(out) :: &
         secondsToday

    integer, intent(out) :: &
         dayOfYear

    integer :: &
         H, M, S, S_d, S_n

    call mpas_get_time(currentTime, DoY=dayOfYear, H=H, M=M, S=S, S_n=S_n, S_d=S_d)

    secondsToday = real(H,RKIND)*3600.0_RKIND + &
                   real(M,RKIND)*60.0_RKIND + &
                   real(S,RKIND) + &
                   real(S_n,RKIND)/real(S_d,RKIND)

  end subroutine get_seconds_today

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
       secondsToday, &
       dayOfYear)

    use seaice_constants, only: &
         seaiceDegreesToRadians, &
         seaiceSecondsPerDay, &
         pii

    real(kind=RKIND), intent(out) :: &
         shortwaveDown

    real(kind=RKIND), intent(in) :: &
         longitudeIn, &
         latitude, &
         cloudFraction, &
         airSpecificHumidity

    real(kind=RKIND), intent(in) :: &
         secondsToday

    integer, intent(in) :: &
         dayOfYear

    real(kind=RKIND) :: &
         longitude, &
         solarTime, &
         hourAngle, &
         declination, &
         cosZ, &
         e, &
         d, &
         sw0

    ! longitude needs to be [-pi,pi] not [0,2pi]
    longitude = longitudeIn
    if (longitude > pii) longitude = longitude - 2.0_RKIND * pii

    solarTime = mod(real(secondsToday,kind=RKIND),seaiceSecondsPerDay)/3600.0_RKIND + 12.0_RKIND*sin(0.5_RKIND*longitude)

    hourAngle = (12.0_RKIND - solarTime)*pii/12.0_RKIND

    ! solar declinatiom
    declination = 23.44_RKIND*cos((172.0_RKIND-real(dayOfYear,RKIND)) * 2.0_RKIND*pii/365.0_RKIND)*seaiceDegreesToRadians

    ! solar zenith angle
    cosZ = sin(latitude)*sin(declination) + cos(latitude)*cos(declination)*cos(hourAngle)
    cosZ = max(cosZ,0.0_RKIND)

    e = 1.0e5_RKIND*airSpecificHumidity/(0.622_RKIND + 0.378_RKIND*airSpecificHumidity)

    d = (cosZ + 2.7_RKIND)*e*1.0e-5_RKIND+1.085_RKIND*cosZ+0.1_RKIND

    sw0 = 1353.0_RKIND*cosZ**2/d
    sw0 = max(sw0,0.0_RKIND)

    shortwaveDown = sw0 * (1.0_RKIND - 0.6_RKIND * cloudFraction**3)

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
    case ("ISPOL")
       call init_oceanic_forcing_ISPOL(domain)
    case ("hourly_velocities")
       call init_oceanic_forcing_hourly_velocities(domain)
    case ("som_varying_thickness")
       call init_oceanic_forcing_som_varying_thickness(domain)
    case default
       call mpas_log_write("Oceanic forcing type unknown: "//trim(config_forcing_sst_type), MPAS_LOG_CRIT)
    end select

  end subroutine init_oceanic_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ocean_bgc_forcing
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 3rd March 2023
!> \details
!> Uses an ocean biogeochemistry forcing file
!  for bgc enabled sea ice only test cases.
!-----------------------------------------------------------------------

  subroutine init_oceanic_bgc_forcing(domain)

    type (domain_type) :: domain

    character(len=strKIND), pointer :: &
         config_forcing_bgc_type

    call MPAS_pool_get_config(domain % configs, "config_forcing_bgc_type", config_forcing_bgc_type)

    select case (trim(config_forcing_bgc_type))
    case ("ISPOL")
       call init_oceanic_bgc_forcing_ISPOL(domain)
    end select

  end subroutine init_oceanic_bgc_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_atm_iron_bgc_forcing
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 5th April 2024
!> \details
!> Uses monthly iron solubility maps from CAM6-MIMI (Yan Feng)
!  to convert dust to dissolved iron for bgc.
!-----------------------------------------------------------------------

  subroutine init_atm_iron_bgc_forcing(domain, clock)

    type (domain_type) :: domain

    type (MPAS_Clock_type) :: clock

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    ! get forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    ! create the dust iron solubility forcing group

    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_atm_bgc_forcing_monthly", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "0001-01-15_00:00:00"

    ! iron solubility fraction of wet dust iron
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         'seaice_atm_bgc_forcing_monthly', &
         'IRON_Zolubility_wet', &
         'DustIronMonthlyForcing', &
         'atmos_forcing', &
         'IRON_Zolubility_wet', &
         'linear', &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! iron solubility fraction of dry dust iron
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         'seaice_atm_bgc_forcing_monthly', &
         'IRON_Zolubility_dry', &
         'DustIronMonthlyForcing', &
         'atmos_forcing', &
         'IRON_Zolubility_dry', &
         'linear', &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! atmospheric wet dust flux
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         'seaice_atm_bgc_forcing_monthly', &
         'dust_FLUZ_WET', &
         'DustIronMonthlyForcing', &
         'atmos_forcing', &
         'dust_FLUZ_WET', &
         'linear', &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! atmospheric dry dust flux
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         'seaice_atm_bgc_forcing_monthly', &
         'dust_FLUZ_DRY', &
         'DustIronMonthlyForcing', &
         'atmos_forcing', &
         'dust_FLUZ_DRY', &
         'linear', &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! iron fraction in dust
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         'seaice_atm_bgc_forcing_monthly', &
         'IRON_in_duzt_fraction', &
         'DustIronMonthlyForcing', &
         'atmos_forcing', &
         'IRON_in_duzt_fraction', &
         'linear', &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         'seaice_atm_bgc_forcing_monthly', &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_atm_iron_bgc_forcing

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

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    ! get oceanic forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "0001-01-15_00:00:00"

    ! create the sea surface temperature forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart)

    ! sea surface temperature
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         "seaSurfaceTemperature", &
         "NCARMonthlySSTForcing", &
         "ocean_coupling", &
         "seaSurfaceTemperature", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

    ! create the monthly forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart)

    ! sea surface salinity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
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
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
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
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
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
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
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
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
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
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
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
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "oceanHeatFluxConvergence", &
         "NCARMonthlyForcing", &
         "ocean_coupling", &
         "oceanHeatFluxConvergence", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_oceanic_forcing_ncar

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_oceanic_forcing_ISPOL
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date October 2022
!> \details
!> The physical data corresponds to the Weddell Sea (67.9S, 54W) location from a
!> POP 1-degree (gx1v3) field derived from an ocean-ice simulation
!> forced with Ice Station Polarstern (ISPOL) 2004 atmospheric data.
!> Biogeochemical data is World Ocean Atlas (NCEP 2013) surface nutrient data.
!
!-----------------------------------------------------------------------

  subroutine init_oceanic_forcing_ISPOL(domain)

    type (domain_type) :: domain

    logical, pointer :: &
         config_do_restart

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly, &
         forcingIntervaldaily, &
         forcingReferenceTimedaily

    ! get oceanic forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    ! get ispol forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_forcing_start_time", config_forcing_start_time)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_start", config_forcing_cycle_start)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_duration", config_forcing_cycle_duration)

    ! create the sea surface temperature forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         domain, &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_do_restart)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "0001-01-15_00:00:00"

    ! sea surface temperature
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "seaSurfaceTemperature", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceTemperature", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! sea surface salinity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "seaSurfaceSalinity", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceSalinity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! u ocean velocity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "uOceanVelocity", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "uOceanVelocity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! v ocean velocity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "vOceanVelocity", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "vOceanVelocity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! ocean mixed layer depth
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "oceanMixedLayerDepth", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "oceanMixedLayerDepth", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! ocean sea surface tilt
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "seaSurfaceTiltU", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceTiltU", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "seaSurfaceTiltV", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceTiltV", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! deep ocean heat flux convergence
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "oceanHeatFluxConvergence", &
         "ISPOLMonthlyForcing", &
         "ocean_coupling", &
         "oceanHeatFluxConvergence", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_oceanic_forcing_ISPOL

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_oceanic_forcing_ISPOL
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date March 2023
!> \details
!>  surface nutrient data corresponds to the Weddell Sea (67.9S, 54W) location from
!>  the World Ocean Atlas (NCEP 2013).
!
!-----------------------------------------------------------------------

  subroutine init_oceanic_bgc_forcing_ISPOL(domain)

    type (domain_type) :: domain

    logical, pointer :: &
         config_do_restart

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly, &
         forcingIntervaldaily, &
         forcingReferenceTimedaily

    ! get oceanic forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    ! get ispol forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_forcing_start_time", config_forcing_start_time)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_start", config_forcing_cycle_start)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_duration", config_forcing_cycle_duration)


    ! ocean surface biogeochemistry fields

    ! create the sea surface biogeochemistry forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_oceanic_biogeochemistry_forcing_daily", &
         domain, &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_do_restart)

    forcingIntervalDaily = "00-00-01_00:00:00"
    forcingReferenceTimeDaily = "2004-01-01_00:00:00"

    ! sea surface nitrate
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_biogeochemistry_forcing_daily", &
         "oceanNitrateConc1D", &
         "ISPOLDailyWOAForcing", &
         "1D_ocean_biogeochemistry_forcing", &
         "oceanNitrateConc1D", &
         "linear", &
         forcingReferenceTimeDaily, &
         forcingIntervalDaily)

    ! sea surface silicate
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_biogeochemistry_forcing_daily", &
         "oceanSilicateConc1D", &
         "ISPOLDailyWOAForcing", &
         "1D_ocean_biogeochemistry_forcing", &
         "oceanSilicateConc1D", &
         "linear", &
         forcingReferenceTimeDaily, &
         forcingIntervalDaily)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_oceanic_biogeochemistry_forcing_daily", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_oceanic_bgc_forcing_ISPOL

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_oceanic_forcing_hourly_velocities
!
!> \brief Initialize the forcing objects
!> \author Adrian K. Turner, LANL
!> \date 18th November 2021
!> \details
!>  This routine calls the MPAS_forcing module subroutines that initializes
!>  the forcings type
!
!-----------------------------------------------------------------------

  subroutine init_oceanic_forcing_hourly_velocities(domain)

    type(domain_type) :: domain

    real(kind=RKIND), pointer :: &
         config_dt

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalHourlyVelocities, &
         forcingReferenceTimeHourlyVelocities, &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    ! get atmospheric forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_forcing_start_time", config_forcing_start_time)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_start", config_forcing_cycle_start)
    call MPAS_pool_get_config(domain % configs, "config_forcing_cycle_duration", config_forcing_cycle_duration)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    ! create the six hourly forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_hourly_velocities", &
         domain, &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_do_restart, &
         .false.)

    forcingIntervalHourlyVelocities = "01:00:00"
    forcingReferenceTimeHourlyVelocities = "2000-01-01_00:00:00"

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_hourly_velocities", &
         "uOceanVelocity", &
         "oceanHourlyVelocitiesForcing", &
         "ocean_coupling", &
         "uOceanVelocity", &
         "linear", &
         forcingReferenceTimeHourlyVelocities, &
         forcingIntervalHourlyVelocities, &
         "next")

    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_hourly_velocities", &
         "vOceanVelocity", &
         "oceanHourlyVelocitiesForcing", &
         "ocean_coupling", &
         "vOceanVelocity", &
         "linear", &
         forcingReferenceTimeHourlyVelocities, &
         forcingIntervalHourlyVelocities, &
         "next")

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_hourly_velocities", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_oceanic_forcing_hourly_velocities

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_oceanic_forcing_som_varying_thickness
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 14th September 2022
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_oceanic_forcing_som_varying_thickness(domain)

    type (domain_type) :: domain

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    ! get oceanic forcing configuration options
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "0001-01-15_00:00:00"

    ! create the sea surface temperature forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart)

    ! sea surface temperature
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         "seaSurfaceTemperature", &
         "SOMVaryingThicknessInitForcing", &
         "ocean_coupling", &
         "seaSurfaceTemperature", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! sea surface salinity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         "seaSurfaceSalinity", &
         "SOMVaryingThicknessInitForcing", &
         "ocean_coupling", &
         "seaSurfaceSalinity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! ocean mixed layer depth
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         "oceanMixedLayerDepth", &
         "SOMVaryingThicknessInitForcing", &
         "ocean_coupling", &
         "oceanMixedLayerDepth", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_init", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

    ! create the monthly forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         domain, &
         '0000-01-01_00:00:00', &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart)

    ! u ocean velocity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "uOceanVelocity", &
         "SOMVaryingThicknessMonthlyForcing", &
         "ocean_coupling", &
         "uOceanVelocity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! v ocean velocity
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "vOceanVelocity", &
         "SOMVaryingThicknessMonthlyForcing", &
         "ocean_coupling", &
         "vOceanVelocity", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! u surface tilt
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "seaSurfaceTiltU", &
         "SOMVaryingThicknessMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceTiltU", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! v surface tilt
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "seaSurfaceTiltV", &
         "SOMVaryingThicknessMonthlyForcing", &
         "ocean_coupling", &
         "seaSurfaceTiltV", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    ! deep ocean heat flux convergence
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         "oceanHeatFluxConvergence", &
         "SOMVaryingThicknessMonthlyForcing", &
         "ocean_coupling", &
         "oceanHeatFluxConvergence", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_oceanic_forcing_monthly", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_oceanic_forcing_som_varying_thickness

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
         config_forcing_sst_type, &
         config_forcing_bgc_type

    logical, pointer :: &
         config_do_restart, &
         config_use_column_biogeochemistry

    ! configurations
    call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)
    call mpas_pool_get_config(domain % configs, 'config_forcing_sst_type', config_forcing_sst_type)
    call mpas_pool_get_config(domain % configs, 'config_forcing_bgc_type', config_forcing_bgc_type)
    call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
    call mpas_pool_get_config(domain % configs, 'config_use_column_biogeochemistry', config_use_column_biogeochemistry)

    ! use the forcing layer to get data
    if (trim(config_forcing_sst_type) == 'ncar' .or. &
        trim(config_forcing_sst_type) == 'som_varying_thickness') then

       call MPAS_forcing_get_forcing(&
            seaiceForcingGroups, &
            "seaice_oceanic_forcing_monthly", &
            streamManager, &
            config_dt)

       ! only get sst data on first timestep
       if (firstTimeStep .and. .not. config_do_restart) then

          call MPAS_forcing_get_forcing(&
               seaiceForcingGroups, &
               "seaice_oceanic_forcing_init", &
               streamManager, &
               config_dt)

       endif

    else if (trim(config_forcing_sst_type) == 'hourly_velocities') then

       call MPAS_forcing_get_forcing(&
            seaiceForcingGroups, &
            "seaice_oceanic_forcing_hourly_velocities", &
            streamManager, &
            config_dt)

    else if (trim(config_forcing_sst_type) == 'ISPOL') then

       call MPAS_forcing_get_forcing(&
            seaiceForcingGroups, &
            "seaice_oceanic_forcing_monthly", &
            streamManager, &
            config_dt)
    endif
    if (config_use_column_biogeochemistry) then
       if (trim(config_forcing_bgc_type) == 'ISPOL') then

          call MPAS_forcing_get_forcing(&
               seaiceForcingGroups, &
               "seaice_oceanic_biogeochemistry_forcing_daily", &
               streamManager, &
               config_dt)
       endif
    endif

    block => domain % blocklist
    do while (associated(block))

       ! convert the input forcing variables to the coupling variables
       select case (trim(config_forcing_sst_type))
       case ("ncar","som_varying_thickness")
          call prepare_oceanic_coupling_variables_ncar(block, firstTimeStep)
       case ("ISPOL")
          call prepare_oceanic_coupling_variables_ISPOL(block, firstTimeStep)
       end select

       if (config_use_column_biogeochemistry) then
          ! convert the input forcing variables to the coupling variables
          select case (trim(config_forcing_bgc_type))
          case ("ISPOL")
             call prepare_oceanic_coupling_bgc_variables_ISPOL(block, firstTimeStep)
          end select
       end if

       call post_oceanic_coupling(block)

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

    use icepack_intfc, only: &
         icepack_sea_freezing_temperature
    use seaice_column, only: &
         seaice_column_sea_freezing_temperature

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

    character(len=strKIND), pointer :: &
         config_sea_freezing_temperature_type, &
         config_column_physics_type

    logical, pointer :: &
         config_do_restart

    integer, pointer :: &
         nCellsSolve, &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_config(block % configs, "config_column_physics_type", config_column_physics_type)
    call MPAS_pool_get_config(block % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(block % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)
    call MPAS_pool_get_array(ocean_coupling, "oceanMixedLayerDepth", oceanMixedLayerDepth)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

    if (trim(config_column_physics_type) == "icepack") then
       do iCell = 1, nCellsSolve

          ! ensure physical realism
          seaSurfaceSalinity(iCell)   = max(seaSurfaceSalinity(iCell), 0.0_RKIND)
          oceanMixedLayerDepth(iCell) = max(oceanMixedLayerDepth(iCell), 0.0_RKIND)

          ! sea freezing temperature
          seaFreezingTemperature(iCell) = icepack_sea_freezing_temperature(seaSurfaceSalinity(iCell))

       enddo ! iCell
    else if (trim(config_column_physics_type) == "column_package") then
       do iCell = 1, nCellsSolve

          ! ensure physical realism
          seaSurfaceSalinity(iCell)   = max(seaSurfaceSalinity(iCell), 0.0_RKIND)
          oceanMixedLayerDepth(iCell) = max(oceanMixedLayerDepth(iCell), 0.0_RKIND)

          ! sea freezing temperature
          seaFreezingTemperature(iCell) = seaice_column_sea_freezing_temperature(seaSurfaceSalinity(iCell))

       enddo ! iCell
    endif ! config_column_physics_type

    ! only update sea surface temperature on first non-restart timestep
    if (firstTimeStep .and. .not. config_do_restart) then

       do iCell = 1, nCellsSolve

          ! sea surface temperature
          seaSurfaceTemperature(iCell) = max(seaSurfaceTemperature(iCell), seaFreezingTemperature(iCell))

       enddo ! iCell

    endif

  end subroutine prepare_oceanic_coupling_variables_ncar

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_oceanic_coupling_variables_ISPOL
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date  October, 2022
!> \details
!> Prepares the ocean coupling variables from the ISPOL 1D forcing
!
!-----------------------------------------------------------------------

  subroutine prepare_oceanic_coupling_variables_ISPOL(block, firstTimeStep)

    use ice_colpkg, only: &
         colpkg_sea_freezing_temperature

    use icepack_intfc, only: &
         icepack_sea_freezing_temperature

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

    character(len=strKIND), pointer :: &
         config_sea_freezing_temperature_type, &
         config_column_physics_type

    logical, pointer :: &
         config_do_restart

    integer, pointer :: &
         nCellsSolve, &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_config(block % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(block % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)
    call MPAS_pool_get_array(ocean_coupling, "oceanMixedLayerDepth", oceanMixedLayerDepth)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

    if (trim(config_column_physics_type) == "icepack") then
       do iCell = 1, nCellsSolve

          ! ensure physical realism
          seaSurfaceSalinity(iCell)   = max(seaSurfaceSalinity(iCell), 0.0_RKIND)
          oceanMixedLayerDepth(iCell) = max(oceanMixedLayerDepth(iCell), 0.0_RKIND)

          ! sea freezing temperature
          seaFreezingTemperature(iCell) = icepack_sea_freezing_temperature(seaSurfaceSalinity(iCell))

       enddo ! iCell
    else if (trim(config_column_physics_type) == "column_package") then
       do iCell = 1, nCellsSolve

          ! ensure physical realism
          seaSurfaceSalinity(iCell)   = max(seaSurfaceSalinity(iCell), 0.0_RKIND)
          oceanMixedLayerDepth(iCell) = max(oceanMixedLayerDepth(iCell), 0.0_RKIND)

          ! sea freezing temperature
          seaFreezingTemperature(iCell) = colpkg_sea_freezing_temperature(seaSurfaceSalinity(iCell))

       enddo ! iCell
    endif ! config_column_physics_type

    ! only update sea surface temperature on first non-restart timestep
    if (firstTimeStep .and. .not. config_do_restart) then

       do iCell = 1, nCellsSolve

          ! sea surface temperature
          seaSurfaceTemperature(iCell) = max(seaSurfaceTemperature(iCell), seaFreezingTemperature(iCell))

       enddo ! iCell

    endif

  end subroutine prepare_oceanic_coupling_variables_ISPOL

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_oceanic_coupling_variables_ISPOL
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date  March, 2023
!> \details
!> Defines the ocean bgc cell fields corresponding to the 1D forcing fields.
!
!-----------------------------------------------------------------------

  subroutine prepare_oceanic_coupling_bgc_variables_ISPOL(block, firstTimeStep)

    use seaice_constants, only: &
         oceanAmmoniumISPOL, oceanDMSISPOL, oceanDMSPISPOL, oceanDiatomsISPOL, &
         oceanSmallAlgaeISPOL, oceanPhaeocystisISPOL, oceanPolysaccharidsISPOL, &
         oceanLipidsISPOL, oceanDICISPOL, oceanProteinsISPOL, oceanDissolvedIronISPOL, &
         oceanParticulateIronISPOL, oceanHumicsISPOL, oceanProteinsCarbonISPOL

    type (block_type), pointer :: block

    logical, intent(in) :: &
         firstTimeStep

    type (mpas_pool_type), pointer :: &
         mesh, &
         oceanBGCForcing, &
         biogeochemistry

    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanDONConc, &
         oceanDICConc, &
         oceanDOCConc, &
         oceanAlgaeConc, &
         oceanParticulateIronConc, &
         oceanDissolvedIronConc

    real(kind=RKIND), dimension(:), pointer :: &
         oceanNitrateConc, &
         oceanSilicateConc, &
         oceanAmmoniumConc, &
         oceanDMSConc, &
         oceanDMSPConc, &
         oceanHumicsConc

    real(kind=RKIND), pointer :: &
         oceanNitrateConc1D, &
         oceanSilicateConc1D

    logical, pointer :: &
         config_do_restart, &
         config_use_column_biogeochemistry

    integer, pointer :: &
         nCellsSolve, &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_config(block % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_subpool(block % structs, "biogeochemistry",biogeochemistry)
    call MPAS_pool_get_subpool(block % structs, "1D_ocean_biogeochemistry_forcing", &
            oceanBGCForcing)

    call MPAS_pool_get_array(oceanBGCForcing, "oceanNitrateConc1D", oceanNitrateConc1D)
    call MPAS_pool_get_array(oceanBGCForcing, "oceanSilicateConc1D", oceanSilicateConc1D)
    call MPAS_pool_get_array(biogeochemistry, "oceanSilicateConc", oceanSilicateConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanNitrateConc", oceanNitrateConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeConc", oceanAlgaeConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumConc", oceanAmmoniumConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanDMSConc", oceanDMSConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanDMSPConc", oceanDMSPConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanHumicsConc", oceanHumicsConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanDONConc", oceanDONConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanDICConc", oceanDICConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanDOCConc", oceanDOCConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronConc", oceanParticulateIronConc)
    call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronConc", oceanDissolvedIronConc)

    do iCell = 1, nCellsSolve
       ! define cell quantities
       oceanNitrateConc(iCell) = oceanNitrateConc1D
       oceanSilicateConc(iCell) = oceanSilicateConc1D
    end do

    if (firstTimeStep) then
       do iCell = 1, nCellsSolve
          oceanAmmoniumConc(iCell) = oceanAmmoniumISPOL
          oceanDMSConc(iCell) = oceanDMSISPOL
          oceanDMSPConc(iCell) = oceanDMSPISPOL
          oceanAlgaeConc(:,iCell) = oceanDiatomsISPOL
          oceanAlgaeConc(2,iCell) = oceanSmallAlgaeISPOL
          oceanAlgaeConc(3,iCell) = oceanPhaeocystisISPOL
          oceanDOCConc(:,iCell) = oceanPolysaccharidsISPOL
          oceanDOCConc(2,iCell) = oceanLipidsISPOL
          oceanDOCConc(3,iCell) = oceanProteinsCarbonISPOL
          oceanDONConc(:,iCell) = oceanProteinsISPOL
          oceanDICConc(:,iCell) = oceanDICISPOL
          oceanHumicsConc(iCell) = oceanHumicsISPOL
          oceanParticulateIronConc(:,iCell) = oceanParticulateIronISPOL
          oceanDissolvedIronConc(:,iCell) = oceanDissolvedIronISPOL
       enddo ! iCell
    endif

  end subroutine prepare_oceanic_coupling_bgc_variables_ISPOL

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  post_oceanic_coupling
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine post_oceanic_coupling(block)

    use seaice_mesh, only: &
         seaice_latlon_vector_rotation_forward

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: &
         mesh, &
         ocean_coupling

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

    logical, pointer :: &
         on_a_sphere, &
         config_rotate_cartesian_grid

    real(kind=RKIND), pointer :: &
         sphere_radius

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_config(block % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_config(mesh, "on_a_sphere", on_a_sphere)
    if (on_a_sphere .and. config_rotate_cartesian_grid) then

       call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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

       ! make rotated ocean velocity available across all cells
       do iCell = 1, nCells

          call seaice_latlon_vector_rotation_forward(&
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

          call seaice_latlon_vector_rotation_forward(&
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

    endif ! on_a_sphere

  end subroutine post_oceanic_coupling

!-----------------------------------------------------------------------
! data iceberg forcing
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_data_iceberg_forcing
!
!> \brief
!> \author Darin Comeau, LANL
!> \date 27th August 2018
!> \details This initializes the data iceberg forcing group.
!>
!-----------------------------------------------------------------------

  subroutine init_data_iceberg_forcing(domain, clock)

    type (domain_type) :: domain

    type (MPAS_Clock_type) :: clock

    logical, pointer :: &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    type (MPAS_Time_Type) :: currTime
    character(len=strKIND) :: timeStamp
    integer :: ierr

    ! get configuration options
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    forcingIntervalMonthly = "00-01-00_00:00:00"
    forcingReferenceTimeMonthly = "0001-01-15_00:00:00"

    currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
    timeStamp = '0000'//trim(timeStamp(5:))

    ! create own data iceberg forcing group
    call MPAS_forcing_init_group(&
         seaiceForcingGroups, &
         "seaice_data_iceberg_forcing_monthly", &
         domain, &
         timeStamp, &
         '0000-01-01_00:00:00', &
         '0001-00-00_00:00:00', &
         config_do_restart)

    ! iceberg freshwater fluxes
    call MPAS_forcing_init_field(&
         domain % streamManager, &
         seaiceForcingGroups, &
         "seaice_data_iceberg_forcing_monthly", &
         "bergFreshwaterFluxData", &
         "dataIcebergForcing", &
         "berg_forcing", &
         "bergFreshwaterFluxData", &
         "linear", &
         forcingReferenceTimeMonthly, &
         forcingIntervalMonthly)

    call MPAS_forcing_init_field_data(&
         seaiceForcingGroups, &
         "seaice_data_iceberg_forcing_monthly", &
         domain % streamManager, &
         config_do_restart, &
         .false.)

  end subroutine init_data_iceberg_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  data_iceberg_forcing
!
!> \brief
!> \author Darin Comeau, LANL
!> \date 27th Aug 2018
!> \details This routine is the timestep for getting and setting fluxes
!> from data icebergs.
!
!-----------------------------------------------------------------------

  subroutine data_iceberg_forcing(&
       streamManager, &
       domain)

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    ! Arguments for verbose debugging
    ! type (MPAS_time_type) :: currentForcingTime
    ! character(len=strKIND) :: currentForcingTimeStr

    real(kind=RKIND), pointer :: &
         config_dt

    call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

    ! For verbose debugging.
    ! Uncomment the lines below and the arguments above to print the forcing time
    ! to the log file, to ensure the forcing times match the simulation times.

    ! call MPAS_forcing_get_forcing_time(&
    !    seaiceForcingGroups, &
    !    "seaice_data_iceberg_forcing_monthly", &
    !    currentForcingTime)

    ! call MPAS_get_time(currentForcingTime, dateTimeString=currentForcingTimeStr)
    ! call MPAS_log_write('Get Data icebergs at: '//trim(currentForcingTimeStr))

    ! use the forcing layer to get
    call MPAS_forcing_get_forcing(&
       seaiceForcingGroups, &
       "seaice_data_iceberg_forcing_monthly", &
       streamManager, &
       config_dt)

    call get_data_iceberg_fluxes(domain)

  end subroutine data_iceberg_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_data_iceberg_fluxes
!
!> \brief   Initialize with icebergs calving mass
!> \author  Darin Comeau, LANL
!> \date    20 Aug 2018
!> \details This routine is intended to initialize the set data iceberg
!> meltwater and latent heat fluxes from a forcing file of monthly
!> climatologies.
!
!-----------------------------------------------------------------------

  subroutine get_data_iceberg_fluxes(domain)

    use seaice_constants, only: &
         seaiceLatentHeatMelting ! latent heat of melting of fresh ice (J/kg)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         berg_forcing, &
         berg_fluxes

    logical, pointer :: &
         config_scale_dib_by_removed_ice_runoff

    integer, pointer :: &
         nCellsSolve

    real(kind=RKIND), dimension(:), pointer :: &
         bergFreshwaterFluxData, & ! iceberg freshwater flux read in from file (kg/m^2/s)
         bergFreshwaterFlux, & ! iceberg freshwater flux for ocean (kg/m^2/s)
         bergLatentHeatFlux ! iceberg latent heat flux for ocean (J/m^2/s)

    real(kind=RKIND), pointer :: &
         areaIntegAnnMeanDataIcebergIceShelfFreshwaterFlux, & ! area integrated, annual mean freshwater flux from icebegs and ice shelves (kg/s)
         runningMeanRemovedIceRunoff ! the area integrated, running mean of removed ice runoff (kg/s)

    integer :: &
         iCell

    real(kind=RKIND) :: &
         scaling

    ! dc including as parameters here so as not to create new namelist options
    real(kind=RKIND), parameter :: &
         specificHeatFreshIce = 2106.0_RKIND, & ! specific heat of fresh ice J * kg^-1 * K^-1
         bergTemperature = -4.0_RKIND           ! iceberg temperature, assumed constant


    call MPAS_pool_get_config(domain % configs, "config_scale_dib_by_removed_ice_runoff", &
                              config_scale_dib_by_removed_ice_runoff)

    if (config_scale_dib_by_removed_ice_runoff) then
       block => domain % blocklist
       call MPAS_pool_get_subpool(block % structs, "berg_forcing", berg_forcing)
       call mpas_pool_get_array(berg_forcing, 'runningMeanRemovedIceRunoff', runningMeanRemovedIceRunoff)
       call mpas_pool_get_array(berg_forcing, 'areaIntegAnnMeanDataIcebergIceShelfFreshwaterFlux', &
                                areaIntegAnnMeanDataIcebergIceShelfFreshwaterFlux)
       scaling = runningMeanRemovedIceRunoff / areaIntegAnnMeanDataIcebergIceShelfFreshwaterFlux
    else
       scaling = 1.0_RKIND
    end if

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "berg_forcing", berg_forcing)
       call MPAS_pool_get_subpool(block % structs, "berg_fluxes", berg_fluxes)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_array(berg_forcing, "bergFreshwaterFluxData", bergFreshwaterFluxData)
       call MPAS_pool_get_array(berg_fluxes, "bergFreshwaterFlux", bergFreshwaterFlux)
       call MPAS_pool_get_array(berg_fluxes, "bergLatentHeatFlux", bergLatentHeatFlux)

       do iCell = 1, nCellsSolve

          bergFreshwaterFlux(iCell) = scaling * bergFreshwaterFluxData(iCell)
          bergLatentHeatFlux(iCell) = -scaling * bergFreshwaterFluxData(iCell) * &
                                     (seaiceLatentHeatMelting - specificHeatFreshIce*bergTemperature)

       enddo

       block => block % next
    enddo

  end subroutine get_data_iceberg_fluxes

!-----------------------------------------------------------------------
! coupler fluxes initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_reset_coupler_fluxes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_reset_coupler_fluxes(domain)

    type(domain_type) :: domain

    logical, pointer :: &
         config_use_column_physics

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_physics", config_use_column_physics)

    if (config_use_column_physics) then

       call reset_atmospheric_coupler_fluxes(domain)
       call reset_ocean_coupler_fluxes(domain)

    endif

  end subroutine seaice_reset_coupler_fluxes

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
         biogeochemistry

    real(kind=RKIND), dimension(:), pointer :: &
         sensibleHeatFlux, &
         latentHeatFlux, &
         absorbedShortwaveFlux, &
         longwaveUp, &
         evaporativeWaterFlux, &
         airStressCellU, &
         airStressCellV, &
         atmosReferenceSpeed10m, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m

    real(kind=RKIND), dimension(:,:), pointer :: &
         atmosBioFluxes, &
         atmosBlackCarbonFlux, &
         atmosDustFlux, &
         atmosWetDustFlux, &
         atmosDryDustFlux

    logical, pointer :: &
         config_use_zaerosols

    block => domain % blocklist
    do while (associated(block))

       ! physics
       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)

       call MPAS_pool_get_array(atmosFluxes, "sensibleHeatFlux", sensibleHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "latentHeatFlux", latentHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "longwaveUp", longwaveUp)
       call MPAS_pool_get_array(atmosFluxes, "evaporativeWaterFlux", evaporativeWaterFlux)

       call MPAS_pool_get_array(velocitySolver, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolver, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(shortwave, "absorbedShortwaveFlux", absorbedShortwaveFlux)

       call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceSpeed10m', atmosReferenceSpeed10m)
       call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceTemperature2m', atmosReferenceTemperature2m)
       call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceHumidity2m', atmosReferenceHumidity2m)

       sensibleHeatFlux            = 0.0_RKIND
       latentHeatFlux              = 0.0_RKIND
       longwaveUp                  = 0.0_RKIND
       evaporativeWaterFlux        = 0.0_RKIND

       airStressCellU              = 0.0_RKIND
       airStressCellV              = 0.0_RKIND

       absorbedShortwaveFlux       = 0.0_RKIND

       atmosReferenceSpeed10m      = 0.0_RKIND
       atmosReferenceTemperature2m = 0.0_RKIND
       atmosReferenceHumidity2m    = 0.0_RKIND

       ! zaerosols
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols", config_use_zaerosols)

       if (config_use_zaerosols) then

          call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

          call MPAS_pool_get_array(biogeochemistry, "atmosBioFluxes", atmosBioFluxes)
          call MPAS_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
          call MPAS_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)
          call MPAS_pool_get_array(biogeochemistry, "atmosWetDustFlux", atmosWetDustFlux)
          call MPAS_pool_get_array(biogeochemistry, "atmosDryDustFlux", atmosDryDustFlux)

          atmosBioFluxes              = 0.0_RKIND
          atmosBlackCarbonFlux        = 0.0_RKIND
          atmosDustFlux               = 0.0_RKIND
          atmosWetDustFlux            = 0.0_RKIND
          atmosDryDustFlux            = 0.0_RKIND

       endif

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
         icebergFluxes, &
         biogeochemistry

    real(kind=RKIND), dimension(:), pointer :: &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux

    real(kind=RKIND), dimension(:), pointer :: &
         bergFreshwaterFlux, &
         bergLatentHeatFlux

    real(kind=RKIND), dimension(:), pointer :: &
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

    logical, pointer :: &
         config_use_column_biogeochemistry, &
         config_use_data_icebergs

    block => domain % blocklist
    do while (associated(block))

       ! physics
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)

       call MPAS_pool_get_array(oceanFluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(oceanFluxes, "oceanShortwaveFlux", oceanShortwaveFlux)

       oceanFreshWaterFlux      = 0.0_RKIND
       oceanSaltFlux            = 0.0_RKIND
       oceanHeatFlux            = 0.0_RKIND
       oceanShortwaveFlux       = 0.0_RKIND

       ! biogeochemistry
       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

       if (config_use_column_biogeochemistry) then

          call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

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

       endif

       ! data icebergs
       call MPAS_pool_get_config(block % configs, "config_use_data_icebergs", config_use_data_icebergs)

       if (config_use_data_icebergs) then

          call MPAS_pool_get_subpool(block % structs, "berg_fluxes", icebergFluxes)

          call MPAS_pool_get_array(icebergFluxes, "bergFreshwaterFlux", bergFreshwaterFlux)
          call MPAS_pool_get_array(icebergFluxes, "bergLatentHeatFlux", bergLatentHeatFlux)

          bergFreshwaterFlux = 0.0_RKIND
          bergLatentHeatFlux = 0.0_RKIND

       endif

       block => block % next
    end do

  end subroutine reset_ocean_coupler_fluxes

!-----------------------------------------------------------------------
! restart
!-----------------------------------------------------------------------
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_forcing_write_restart_times
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_forcing_write_restart_times(domain)

    type(domain_type) :: domain

    logical, pointer :: &
         config_use_forcing, &
         config_use_prescribed_ice, &
         config_use_prescribed_ice_forcing, &
         config_use_data_icebergs, &
         config_use_atm_dust_file, &
         config_use_iron_solubility_file

    call MPAS_pool_get_config(domain % configs, "config_use_forcing", config_use_forcing)
    call MPAS_pool_get_config(domain % configs, "config_use_prescribed_ice", config_use_prescribed_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_prescribed_ice_forcing", config_use_prescribed_ice_forcing)
    call MPAS_pool_get_config(domain % configs, "config_use_data_icebergs", config_use_data_icebergs)
    call MPAS_pool_get_config(domain % configs, "config_use_atm_dust_file", config_use_atm_dust_file)
    call MPAS_pool_get_config(domain % configs, "config_use_iron_solubility_file", config_use_iron_solubility_file)

    if (config_use_forcing .or. &
        (config_use_prescribed_ice .and. config_use_prescribed_ice_forcing) .or. &
        config_use_data_icebergs .or. &
        config_use_atm_dust_file .or. config_use_iron_solubility_file) then

       call MPAS_forcing_write_restart_times(seaiceForcingGroups)

    endif

  end subroutine seaice_forcing_write_restart_times

!-----------------------------------------------------------------------

end module seaice_forcing
