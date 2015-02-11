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

  use mpas_grid_types
  use mpas_timekeeping
  use mpas_forcing
  use mpas_stream_manager

  implicit none

  private
  public :: &
       cice_forcing_init, &
       cice_forcing_get, &
       cice_forcing_write_restart

  type (MPAS_forcing_group_type), pointer :: ciceForcingGroups 

  ! forcing parameters
  real (kind=RKIND), parameter :: &
       fracShortwaveVisibleDirect  = 0.28_RKIND, & ! fraction of incoming shortwave in visible direct band
       fracShortwaveVisibleDiffuse = 0.24_RKIND, & ! fraction of incoming shortwave in visible diffuse band
       fracShortwaveIRDirectDown   = 0.31_RKIND, & ! fraction of incoming shortwave in near IR direct band
       fracShortwaveIRDiffuseDown  = 0.17_RKIND    ! fraction of incoming shortwave in near IR diffuse band

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

  subroutine cice_forcing_init(domain, streamManager)

    type (domain_type) :: domain

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    open(55,file="forcing_text.txt")

    call init_atmospheric_forcing(domain, streamManager)

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

  subroutine init_atmospheric_forcing(domain, streamManager)

    type (domain_type) :: domain

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    character(len=strKIND), pointer :: &
         config_atmospheric_forcing_type

    call MPAS_pool_get_config(domain % configs, "config_atmospheric_forcing_type", config_atmospheric_forcing_type)

    select case (trim(config_atmospheric_forcing_type))
    case ("CORE")
       call init_atmospheric_forcing_CORE(domain, streamManager)
    case default
       write(*,*) "CICE atmospheric forcing type unknown: ", trim(config_atmospheric_forcing_type)
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

  subroutine init_atmospheric_forcing_CORE(domain, streamManager)

    type(domain_type) :: domain

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    real(kind=RKIND), pointer :: &
         config_dt

    character(len=strKIND), pointer :: &
         config_forcing_start_time, &
         config_forcing_cycle_start, &
         config_forcing_cycle_duration, &
         config_forcing_restart_file

    character(len=strKIND) :: &
         forcingIntervalSixHourly, &
         forcingReferenceTimeSixHourly, &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    logical, pointer :: &
         config_do_restart

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
         streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "airTemperature", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "airTemperature", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly)

    call MPAS_forcing_init_field(&
         streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "airSpecificHumidity", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "airSpecificHumidity", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly)

    call MPAS_forcing_init_field(&
         streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "uAirVelocity", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "uAirVelocity", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly)

    call MPAS_forcing_init_field(&
         streamManager, &
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         "vAirVelocity", &
         "LYqSixHourlyForcing", &
         "atmos_coupling", &
         "vAirVelocity", &
         "linear", &
         forcingReferenceTimeSixHourly, &
         forcingIntervalSixHourly)

    call MPAS_forcing_init_field_data(&
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         streamManager)

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
         streamManager, &
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
         streamManager, &
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
         streamManager)

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
       simulationClock)

    type (MPAS_streamManager_type), intent(inout) :: streamManager

    type (domain_type) :: domain

    type (MPAS_clock_type) :: simulationClock

    call atmospheric_forcing(&
         streamManager, &
         domain, &
         simulationClock)

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

    use cice_debug, only: cice_debug_output_array

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
    if (trim(config_atmospheric_forcing_type) /= "none") then

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

       call forcing_testing(block)

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

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: &
         mesh, &
         atmosCoupling, &
         atmosForcing
    
    real(kind=RKIND), dimension(:), pointer :: &
         airLevelHeight, &
         airPotentialTemperature, &
         airTemperature, &
         airSpecificHumidity, &
         airDensity, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         longwaveDown, &
         rainfallRate, &
         snowfallRate, &
         cloudFraction, &
         lonCell, &
         latCell

    type (MPAS_time_type) :: &
         currentForcingTime

    real(kind=RKIND) :: &
         shortwaveDown

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)    
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

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

    ! get the current time
    call MPAS_forcing_get_forcing_time(&
         ciceForcingGroups, &
         "cice_atmospheric_forcing_sixhrly", &
         currentForcingTime)

    do iCell = 1, nCells

       ! atmospheric level height
       airLevelHeight(iCell) = 10.0_RKIND

       ! air potential temperature
       airPotentialTemperature(iCell) = &
            airTemperature(iCell) - ciceFreshWaterFreezingPoint

       ! air density
       airDensity(iCell) = 1.3_RKIND

       ! shortwave
       call shortwave_down(&
            shortwaveDown, &
            lonCell(iCell), &
            latCell(iCell), &
            cloudFraction(iCell), &
            airSpecificHumidity(iCell), &
            currentForcingTime)

       shortwaveVisibleDirectDown(iCell)  = shortwaveDown * fracShortwaveVisibleDirect
       shortwaveVisibleDiffuseDown(iCell) = shortwaveDown * fracShortwaveVisibleDiffuse
       shortwaveIRDirectDown(iCell)       = shortwaveDown * fracShortwaveIRDirectDown
       shortwaveIRDiffuseDown(iCell)      = shortwaveDown * fracShortwaveIRDiffuseDown

       ! longwave radiation
       call longwave_parkinson_and_washington(&
            longWavedown(iCell), &
            airTemperature(iCell), &
            cloudFraction(iCell))

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

       ! wind speed
       windSpeed(iCell) = sqrt(uAirVelocity(iCell)**2 + vAirVelocity(iCell)**2)

       ! shortwave
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

    real(kind=RKIND), intent(out) :: & 
         airTemperature

    rainfallRate = rainfallRate * (12.0_RKIND / ciceSecondsPerYear) ! convert mm per year

    ! divide total precipitation between rain and snow
    snowfallRate = 0.0_RKIND

    if (airTemperature < ciceFreshWaterFreezingPoint) then

       snowfallRate = rainfallRate
       rainfallRate = 0.0_RKIND

    endif

  end subroutine precipitation

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

    real(kind=RKIND), intent(out) :: &
         shortwaveDown

    real(kind=RKIND), intent(in) :: &    
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

    call mpas_get_time(currentTime, DoY=dayOfYear, H=H, M=M, S=S, S_n=S_n, S_d=S_d)

    secondsToday = real(H,RKIND)*3600.0_RKIND + &
                   real(M,RKIND)*60.0_RKIND + &
                   real(S,RKIND) + &
                   real(S_n,RKIND)/real(S_d,RKIND)

    ! longitude needs to be [-pi,pi] not [0,2pi]
    longitude = longitudeIn
    if (longitude > pii) longitude = longitude - 2.0_RKIND * pii

    solarTime = mod(real(secondsToday,kind=RKIND),ciceSecondsPerDay)/3600.0_RKIND + 12.0_RKIND*sin(0.5_RKIND*longitude)

    hourAngle = (12.0_RKIND - solarTime)*pii/12.0_RKIND

    ! solar declinatiom
    declination = 23.44_RKIND*cos((172.0_RKIND-real(dayOfYear,RKIND)) * 2.0_RKIND*pii/365.0_RKIND)*ciceDegreesToRadians
    declination = 0.0_RKIND

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
! restart
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
!  variable_interval_midmonth_forward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function variable_interval_fifteenth_forward(currentTime) result(variableInterval)

    type(MPAS_Time_type), intent(in) :: &
         currentTime !< Input: 

    type(MPAS_TimeInterval_type) :: &
         variableInterval

    call mpas_set_timeInterval(variableInterval, MM=1)

  end function variable_interval_fifteenth_forward

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  variable_interval_midmonth_forward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function variable_interval_fifteenth_backward(currentTime) result(variableInterval)

    type(MPAS_Time_type), intent(in) :: &
         currentTime !< Input: 

    type(MPAS_TimeInterval_type) :: &
         variableInterval

    call mpas_set_timeInterval(variableInterval, MM=1)

  end function variable_interval_fifteenth_backward

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  variable_interval_midmonth_forward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function variable_interval_midmonth_forward(currentTime) result(variableInterval)

    type(MPAS_Time_type), intent(in) :: &
         currentTime !< Input: 

    type(MPAS_TimeInterval_type) :: &
         variableInterval

    integer :: &
         Year, &
         Month

    type(MPAS_Time_type) :: &
         startOfMonth, &
         nextMonth, &
         nextMonthEnd

    type(MPAS_TimeInterval_type) :: &
         monthInterval

    call mpas_get_time(currentTime, Year, Month)

    call mpas_set_time(startOfMonth, Year, Month)

    call mpas_set_timeInterval(monthInterval, MM=1)

    nextMonth = currentTime + monthInterval

    call mpas_get_time(nextMonth, Year, Month)

    Month = Month + 1
    if (Month > 12) then
       Month = 1
       Year = Year + 1
    endif
    call mpas_set_time(nextMonthEnd, Year, Month)

    variableInterval = (nextMonthEnd - startOfMonth) / 2

  end function variable_interval_midmonth_forward

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  variable_interval_midmonth_forward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function variable_interval_midmonth_backward(currentTime) result(variableInterval)

    type(MPAS_Time_type), intent(in) :: &
         currentTime !< Input: 

    type(MPAS_TimeInterval_type) :: &
         variableInterval

    integer :: &
         Year, &
         Month

    type(MPAS_Time_type) :: &
         endOfMonth, &
         previousMonth, &
         previousMonthStart

    type(MPAS_TimeInterval_type) :: &
         monthInterval

    call mpas_get_time(currentTime, Year, Month)

    Month = Month + 1
    if (Month > 12) then
       Month = 1
       Year = Year + 1
    endif
    call mpas_set_time(endOfMonth, Year, Month)

    call mpas_set_timeInterval(monthInterval, MM=1)

    previousMonth = currentTime - monthInterval

    call mpas_get_time(previousMonth, Year, Month)

    call mpas_set_time(previousMonthStart, Year, Month)

    variableInterval = (endOfMonth - previousMonthStart) / 2   

  end function variable_interval_midmonth_backward

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

    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", pool)
    
    call MPAS_pool_get_array(pool, "airTemperature", array)

    call cice_debug_output_array(array)

  end subroutine forcing_testing

!-----------------------------------------------------------------------

end module cice_forcing
