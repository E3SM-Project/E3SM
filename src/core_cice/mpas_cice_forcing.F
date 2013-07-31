module cice_forcing

  use mpas_grid_types
  use mpas_timekeeping

  implicit none

  private

  public :: &
       init_forcing, &
       get_forcing
         
  integer, parameter :: &
       epoch = 1900

  type (MPAS_Time_type) :: &
       epochTime

  type (MPAS_Time_type) :: &
       startYearTime, &
       endYearTime

  integer :: currentYear
  real(kind=RKIND) :: &
       currentYearDuration, &
       currentYearEpoch

  ! individual forcing field derived type
  type forcing_field

     character(len=10) :: forcingType            ! 'interval', 'monthly'
     character(len=200) :: filenameDirectory
     character(len=200) :: filenamePrefix
     character(len=200) :: filenameSuffix
     integer :: currentYear
     integer :: pointsPerYear
     real(kind=RKIND) :: interval                ! forcing data interval
     real(kind=RKIND) :: positionInInterval      ! [0.0,1.0]
     character(len=10) :: filetype               ! 'binary', 'netcdf'

     real(kind=RKIND), dimension(2) :: dataTimes ! times of the stored data
     integer, dimension(2) :: dataRecordNumbers
     integer, dimension(2) :: dataYears
     integer, dimension(2) :: dataOrder
     integer :: ndata
     real(kind=RKIND), dimension(:,:), allocatable :: data

  end type forcing_field

  ! all forcing fields
  type forcing_fields

     integer :: nForcingFields
     character(len=200) :: forcingDirectory
     type(forcing_field), dimension(:), allocatable :: forcingFieldsArray

  end type forcing_fields

  type(forcing_fields) :: atmosForcingFields

  real(kind=RKIND) :: elapsed_time_forcing

contains

  !-------------------------------------------------------------
  ! Specific forcing
  !-------------------------------------------------------------

  subroutine init_forcing(mesh, atmos_forcing)

    type(mesh_type), intent(in) :: mesh

    type(atmos_forcing_type), pointer :: atmos_forcing

    call init_forcing_generic()

    call init_atmos_forcing(mesh)

    elapsed_time_forcing = 0.0_RKIND

    call POP_grid_basis(mesh, &
                        atmos_forcing % uPOPAlongGridVector % array, &
                        atmos_forcing % vPOPAlongGridVector % array)

  end subroutine init_forcing

  !-------------------------------------------------------------

  subroutine init_atmos_forcing(mesh)

    type(mesh_type), intent(in) :: mesh

    atmosForcingFields % nForcingFields = 8
    allocate(atmosForcingFields % forcingFieldsArray(atmosForcingFields % nForcingFields))

    atmosForcingFields % forcingDirectory = "/Users/akt/Work/MPAS-CICE/forcing/DATA/atm/gx3v2m/NCAR_bulk"

    ! air temperature
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(1), "interval", &
         "/NCEP/4XDAILY/STATES/", "t_10.", ".dat", "binaryPOP", 1460, 1.0_RKIND)

    ! u air velocity
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(2), "interval", &
         "/NCEP/4XDAILY/STATES/", "u_10.", ".dat", "binaryPOP", 1460, 1.0_RKIND)

    ! v air velocity
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(3), "interval", &
         "/NCEP/4XDAILY/STATES/", "v_10.", ".dat", "binaryPOP", 1460, 1.0_RKIND)

    ! air density
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(4), "interval", &
         "/NCEP/4XDAILY/STATES/", "dn10.", ".dat", "binaryPOP", 1460, 1.0_RKIND)

    ! air humidity
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(5), "interval", &
         "/NCEP/4XDAILY/STATES/", "q_10.", ".dat", "binaryPOP", 1460, 1.0_RKIND)

    ! shortwave down
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(6), "monthly", &
         "/ISCCPM/MONTHLY/RADFLX/", "swdn.", ".dat", "binaryPOP", 12, 0.5_RKIND)

    ! cloudFraction
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(7), "monthly", &
         "/ISCCPM/MONTHLY/RADFLX/", "cldf.", ".dat", "binaryPOP", 12, 0.5_RKIND)

    ! precipitation
    call init_forcing_field(mesh, atmosForcingFields % forcingFieldsArray(8), "monthly", &
         "/MXA/MONTHLY/PRECIP/", "prec.", ".dat", "binaryPOP", 12, 0.5_RKIND)

  end subroutine init_atmos_forcing

  !-------------------------------------------------------------

  subroutine get_forcing(mesh, clock, atmos_forcing)

    type(mesh_type), intent(in) :: mesh
    type (MPAS_Clock_type), intent(in) :: clock
    type (atmos_forcing_type), pointer :: atmos_forcing

    integer :: &
         year, &
         month, &
         dayOfMonth
    
    real(kind=RKIND) :: &
         forcingTime, &
         forcingTimeEpoch

    elapsed_time_forcing = elapsed_time_forcing + 3600.0_RKIND

    call get_forcing_times(forcingTime, forcingTimeEpoch, year, month, dayOfMonth, clock)

    call get_atmos_forcing(mesh, atmos_forcing, forcingTime, forcingTimeEpoch, year, month, dayOfMonth)

  end subroutine get_forcing

  !-------------------------------------------------------------

  subroutine get_atmos_forcing(mesh, atmos_forcing, forcingTime, forcingTimeEpoch, year, month, dayOfMonth)

    use cice_testing, only: plot_pop_grid_cells, plot_cells_spherical
    
    type(mesh_type), intent(in) :: mesh
    type (atmos_forcing_type), pointer :: atmos_forcing

    integer, intent(in) :: &
         year, &
         month, &
         dayOfMonth
    
    real(kind=RKIND), intent(in) :: &
         forcingTime, &
         forcingTimeEpoch

    integer :: iCell

    ! get the forcing fields
    ! air temperature
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(1), &
         atmos_forcing % airTemperature % array)

    ! u air velocity
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(2), &
         atmos_forcing % uAirVelocityPOP % array)

    ! v air velocity
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(3), &
         atmos_forcing % vAirVelocityPOP % array)

    ! air density
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(4), &
         atmos_forcing % airDensity % array)

    ! air humidity
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(5), &
         atmos_forcing % airSpecificHumidity % array)

    ! shortwave down
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(6), &
         atmos_forcing % shortwaveDown % array)

    ! cloudFraction
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(7), &
         atmos_forcing % cloudFraction % array)

    ! precipitation
    call get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, &
         atmosForcingFields % forcingDirectory, &
         atmosForcingFields % forcingFieldsArray(8), &
         atmos_forcing % precipitationRate % array)

    !iCell = POPiCell(mesh,69,59)
    !write(56,*) elapsed_time_forcing/3600.0_RKIND, &
    !     atmos_forcing % airTemperature % array(iCell), &
    !     atmos_forcing % uAirVelocityPOP % array(iCell), &
    !     atmos_forcing % vAirVelocityPOP % array(iCell), &
    !     atmos_forcing % airDensity % array(iCell), &
    !     atmos_forcing % airSpecificHumidity % array(iCell), &
    !     atmos_forcing % shortwaveDown % array(iCell), &         
    !     atmos_forcing % cloudFraction % array(iCell), &
    !     atmos_forcing % precipitationRate % array(iCell)

    ! plot 'fort.56' using 1:2 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:2 with lines
    ! plot 'fort.56' using 1:3 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:3 with lines
    ! plot 'fort.56' using 1:4 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:4 with lines
    ! plot 'fort.56' using 1:5 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:5 with lines
    ! plot 'fort.56' using 1:6 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:6 with lines
    ! plot 'fort.56' using 1:7 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:7 with lines
    ! plot 'fort.56' using 1:8 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:8 with lines
    ! plot 'fort.56' using 1:9 with linespoints, '/Users/akt/Work/MPAS-CICE/cice_sphere/rundir/fort.55' using 1:9 with lines

    call rotate_vector_POP_to_cartesian(atmos_forcing % uAirVelocity % array,        &   
                                        atmos_forcing % vAirVelocity % array,        &
                                        atmos_forcing % uAirVelocityPOP % array,     &  
                                        atmos_forcing % vAirVelocityPOP % array,     &
                                        atmos_forcing % uPOPAlongGridVector % array, &
                                        atmos_forcing % vPOPAlongGridVector % array)


    ! calculate other quantities
    do iCell = 1, mesh % nCells

       call calculate_derived_forcing_variables( &
            atmos_forcing % longWavedown % array(iCell), &
            atmos_forcing % snowfallRate % array(iCell), &
            atmos_forcing % airPotentialTemperature % array(iCell), &
            atmos_forcing % windSpeed % array(iCell), &
            atmos_forcing % sensibleTransferCoefficient % array(iCell), &
            atmos_forcing % latentTransferCoefficient % array(iCell), &
            atmos_forcing % uAirStress % array(iCell), &
            atmos_forcing % vAirStress % array(iCell), &
            atmos_forcing % airTemperature % array(iCell), &
            atmos_forcing % cloudFraction % array(iCell), &
            atmos_forcing % precipitationRate % array(iCell), &
            atmos_forcing % uAirVelocity % array(iCell), &
            atmos_forcing % vAirVelocity % array(iCell), &
            atmos_forcing % airDensity % array(iCell))

    enddo ! iCell

  end subroutine get_atmos_forcing

  !-------------------------------------------------------------
  ! Physical conversions
  !-------------------------------------------------------------

  subroutine calculate_derived_forcing_variables( &
       longWavedown, &
       snowfallRate, &
       airPotentialTemperature, &
       windSpeed, &
       sensibleTransferCoefficient, &
       latentTransferCoefficient, &
       uAirStress, &
       vAirStress, &
       airTemperature, &
       cloudFraction, &
       precipitationRate, &
       uAirVelocity, &
       vAirVelocity, &
       airDensity)
    
    real(kind=RKIND), intent(out) :: &
         longWavedown, &
         snowfallRate, &
         airPotentialTemperature, &
         windSpeed, &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         uAirStress, &
         vAirStress

    real(kind=RKIND), intent(inout) :: &
         precipitationRate

    real(kind=RKIND), intent(in) :: &
         airTemperature, &
         cloudFraction, &
         uAirVelocity, &
         vAirVelocity, &
         airDensity

    real(kind=RKIND), parameter :: &
         stefanBoltzmann         = 5.79484e-8_RKIND , & ! J m-2 K-4 s-1
         freshWaterFreezingPoint = 273.15_RKIND     , & ! freezing temp of fresh ice (K)
         airSpecificHeat         = 1005.0_RKIND     , & ! specific heat of air (J/kg/K)
         latentHeatSublimation   = 2.835e6_RKIND        ! latent heat, sublimation freshwater (J/kg)

    real(kind=RKIND), parameter :: &
         seconds_per_year        = 365.0_RKIND * 24.0_RKIND * 3600.0_RKIND

    real(kind=RKIND) :: &
         airStressCoefficient

    ! Longwave down
    ! Parkinson, C. L. and W. M. Washington (1979),
    ! Large-scale numerical-model of sea ice,
    ! JGR, 84, 311-337, doi:10.1029/JC084iC01p00311

    longWavedown = stefanBoltzmann * airTemperature**4 * &
                   (1.0_RKIND - 0.261_RKIND * exp(-7.77e-4_RKIND * (freshWaterFreezingPoint - airTemperature)**2)) * &
                   (1.0_RKIND + 0.275_RKIND * cloudFraction)

    ! snowfall rate
    precipitationRate = precipitationRate * (12.0_RKIND / seconds_per_year) ! convert mm per year

    snowfallRate = 0.0_RKIND

    if (airTemperature < freshWaterFreezingPoint) then

       snowfallRate = precipitationRate
       precipitationRate = 0.0_RKIND

    endif

    ! air potential temperature
    airPotentialTemperature = airTemperature - freshWaterFreezingPoint

    ! wind speed
    windSpeed = sqrt(uAirVelocity**2 + vAirVelocity**2)

    ! transfer coefficients
    sensibleTransferCoefficient = 1.20e-3_RKIND * airSpecificHeat        * airDensity * windSpeed
    latentTransferCoefficient   = 1.50e-3_RKIND * latentHeatSublimation  * airDensity * windSpeed

    ! air stresses
    airStressCoefficient = 0.0012_RKIND * airDensity * windSpeed

    uAirStress = uAirVelocity * airStressCoefficient
    vAirStress = vAirVelocity * airStressCoefficient

  end subroutine calculate_derived_forcing_variables

  !-------------------------------------------------------------
  ! Generic forcing routines
  !-------------------------------------------------------------

  subroutine init_forcing_generic()

    ! set currentYear
    currentYear = -1e6
    currentYearDuration = 0.0_RKIND

    ! set epoch time
    call mpas_set_time(epochTime, YYYY=epoch, MM=1, DD=1)

  end subroutine init_forcing_generic

  !-------------------------------------------------------------

  subroutine init_forcing_field(mesh, forcingField, forcingType, &
                                filenameDirectory, filenamePrefix, filenameSuffix, filetype, &
                                pointsPerYear, positionInInterval)

    type(mesh_type), intent(in) :: mesh
    type (forcing_field), intent(inout) :: forcingField

    character(len=*), intent(in) :: forcingType, filenameDirectory, filenamePrefix, filenameSuffix, filetype
    integer, intent(in) :: pointsPerYear
    real(kind=RKIND) :: positionInInterval

    forcingField % forcingType = trim(forcingType)

    forcingField % filenameDirectory = trim(filenameDirectory)
    forcingField % filenamePrefix    = trim(filenamePrefix)
    forcingField % filenameSuffix    = trim(filenameSuffix)

    forcingField % positionInInterval = positionInInterval
    forcingField % filetype           = trim(filetype)
    forcingField % pointsPerYear      = pointsPerYear

    forcingField % dataTimes(:) = -1.0e30_RKIND
    forcingField % dataRecordNumbers(:) = -1
    forcingField % dataYears(:) = 0    

    forcingField % dataOrder(1) = 1
    forcingField % dataOrder(2) = 2
    
    forcingField % currentYear = 0  

    forcingField % ndata = mesh % nCells

    allocate(forcingField % data(1:forcingField % ndata,1:2))

  end subroutine init_forcing_field

  !-------------------------------------------------------------
  
  subroutine get_forcing_times(forcingTime, forcingTimeEpoch, year, month, dayOfMonth, clock)

    real(kind=RKIND), intent(out) :: &
         forcingTime, &
         forcingTimeEpoch

    integer, intent(out) :: &
         year, &
         month, &
         dayOfMonth

    type (MPAS_Clock_type), intent(in) :: clock

    type (MPAS_Time_Type) :: currTime
    integer :: ierr

    type (MPAS_TimeInterval_type) :: &
         epochTimeDifference, &
         startYearTimeDifference, &
         currentYearInterval, &
         epochYearInterval

    ! get the current time
    currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime,     &
                       YYYY=year,              &
                       MM  =month,             &
                       DD  =dayOfMonth,        &
                       ierr=ierr)

    ! current year duration
    if (currentYear /= year) then

       ! start year time
       call mpas_set_time(startYearTime, YYYY=year, MM=1, DD=1)

       ! start year time
       call mpas_set_time(endYearTime, YYYY=year+1, MM=1, DD=1)

       ! year duration
       currentYearInterval = endYearTime - startYearTime
       call mpas_get_timeInterval(currentYearInterval, dt=currentYearDuration)

       ! epoch duration
       epochYearInterval = startYearTime - epochTime
       call mpas_get_timeInterval(epochYearInterval, dt=currentYearEpoch)

       currentYear = year

    endif

    ! epoch time
    epochTimeDifference = currTime - epochTime
    call mpas_get_timeInterval(epochTimeDifference, dt=forcingTimeEpoch)

    ! forcing time
    startYearTimeDifference = currTime - startYearTime
    call mpas_get_timeInterval(startYearTimeDifference, dt=forcingTime)

  end subroutine get_forcing_times

  !-------------------------------------------------------------

  function calc_interval(forcingField) result(interval)

    type (forcing_field), intent(in) :: forcingField
    real(kind=RKIND) :: interval
    
    interval = currentYearDuration / real(forcingField % pointsPerYear,RKIND)

  end function calc_interval

  !-------------------------------------------------------------

  subroutine get_individual_forcing(mesh, forcingTime, forcingTimeEpoch, year, month, dayOfMonth, forcingDirectory, forcingField, forcingArray)

    use mpas_timekeeping

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), intent(in) :: &
         forcingTime, &
         forcingTimeEpoch

    integer, intent(in) :: year, month, dayOfMonth

    character(len=*), intent(in) :: forcingDirectory

    type (forcing_field), intent(inout) :: forcingField

    real(kind=RKIND), dimension(:), intent(out) :: forcingArray

    integer, dimension(2) :: desiredYears
    integer, dimension(2) :: desiredRecordNumbers
    real(kind=RKIND), dimension(2) :: desiredTimes
    
    integer, dimension(2) :: &
         existingDataLocation

    logical, dimension(2) :: &
         existingHasData

    real(kind=RKIND), dimension(2) :: &
         interpolationCoeff

    integer :: &
         nDesired, &
         nExist, &
         iCell

    iCell = POPiCell(mesh,69,59)

    call get_desired_record_positions(forcingTime, year, month, dayOfMonth, forcingField, &
                                      desiredYears, desiredRecordNumbers, desiredTimes)
    !write(59,*) desiredYears, desiredRecordNumbers, desiredTimes

    ! get desired data
    existingDataLocation = 0
    existingHasData      = .false.

    ! fill starting data
    do nDesired = 1, 2
       do nExist = 1, 2
          
          if (forcingField % dataYears(nExist)         == desiredYears(nDesired) .and. &
              forcingField % dataRecordNumbers(nExist) == desiredRecordNumbers(nDesired)) then

             ! use the preexisting data
             existingDataLocation(nDesired) = nExist
             existingHasData(nExist)        = .true.

             forcingField % dataOrder(nDesired) = nExist 
             
             exit
             
          endif
          
       enddo ! nExist
    enddo ! nDesired

    do nDesired = 1, 2

       if (existingDataLocation(nDesired) == 0) then

          ! where do it read it into?
          if (.not. existingHasData(1)) then

             ! read into 1
             call read_data(mesh, forcingDirectory, forcingField, 1, desiredYears(nDesired), desiredRecordNumbers(nDesired), nDesired)

             forcingField % dataTimes(1)         = desiredTimes(nDesired)
             forcingField % dataRecordNumbers(1) = desiredRecordNumbers(nDesired)
             forcingField % dataYears(1)         = desiredYears(nDesired)
             forcingField % dataOrder(nDesired) = 1

             existingHasData(1) = .true.
             !write(60,*) "r:1", nDesired, POPiCell(mesh,70,60), forcingField % data(POPiCell(mesh,70,60),forcingField % dataOrder(nDesired))

          else if (.not. existingHasData(2)) then

             ! read into 2
             call read_data(mesh, forcingDirectory, forcingField, 2, desiredYears(nDesired), desiredRecordNumbers(nDesired), nDesired)

             forcingField % dataTimes(2)         = desiredTimes(nDesired)
             forcingField % dataRecordNumbers(2) = desiredRecordNumbers(nDesired)
             forcingField % dataYears(2)         = desiredYears(nDesired)
             forcingField % dataOrder(nDesired) = 2

             existingHasData(2) = .true.
             !write(60,*) "r:2", nDesired, POPiCell(mesh,70,60), forcingField % data(POPiCell(mesh,70,60),forcingField % dataOrder(nDesired))

          else

             ! major error

          endif

       endif

    enddo ! nDesired

    ! interpolate

    interpolationCoeff(1) = (desiredTimes(2) - forcingTimeEpoch) / (desiredTimes(2) - desiredTimes(1))
    interpolationCoeff(2) = (forcingTimeEpoch - desiredTimes(1)) / (desiredTimes(2) - desiredTimes(1))

    !write(61,*) interpolationCoeff(1), desiredTimes(1), desiredTimes(2), forcingTimeEpoch

    forcingArray = forcingField % data(:,forcingField % dataOrder(1)) * interpolationCoeff(1) + &
                   forcingField % data(:,forcingField % dataOrder(2)) * interpolationCoeff(2)

    !write(57,*) forcingArray(iCell), &
    !     forcingField % data(iCell,forcingField % dataOrder(1)), interpolationCoeff(1), &
    !     forcingField % data(iCell,forcingField % dataOrder(2)), interpolationCoeff(2)

  end subroutine get_individual_forcing

  !-------------------------------------------------------------

  subroutine get_desired_record_positions(secOfYear, year, month, dayOfMonth, forcingField, desiredYears, desiredRecordNumbers, desiredTimes)

    real(kind=RKIND), intent(in) :: secOfYear

    integer, intent(in) :: &
         year, &
         month, &
         dayOfMonth
   
    type (forcing_field), intent(in) :: forcingField

    integer, dimension(2), intent(out) :: &
         desiredYears, &
         desiredRecordNumbers

    real(kind=RKIND), dimension(2), intent(out) :: &
         desiredTimes

    real(kind=RKIND) :: &
         dataPosition

    integer :: &
         intervalNumber

    integer :: nDesire
         
    if (trim(forcingField % forcingType) == "interval") then

       intervalNumber = ceiling(secOfYear / calc_interval(forcingField))
       dataPosition = forcingTimeInYear(intervalNumber, forcingField % positionInInterval, calc_interval(forcingField))

    else if (trim(forcingField % forcingType) == "monthly") then

       intervalNumber = month
       dataPosition = forcingTimeInYearMonthly(year, month, forcingField % positionInInterval)
    
    endif

    desiredYears = year

    if (dataPosition > secOfYear) then

       desiredRecordNumbers(1) = intervalNumber - 1
       desiredRecordNumbers(2) = intervalNumber

    else

       desiredRecordNumbers(1) = intervalNumber
       desiredRecordNumbers(2) = intervalNumber + 1
       
    endif

    do nDesire = 1, 2
       if (desiredRecordNumbers(nDesire) < 1) then
          desiredRecordNumbers(nDesire) = desiredRecordNumbers(nDesire) + forcingField % pointsPerYear
          desiredYears(nDesire) = desiredYears(nDesire) - 1
       endif
       
       if (desiredRecordNumbers(nDesire) > forcingField % pointsPerYear) then
          desiredRecordNumbers(nDesire) = desiredRecordNumbers(nDesire) - forcingField % pointsPerYear
          desiredYears(nDesire) = desiredYears(nDesire) + 1
       endif

       if (trim(forcingField % forcingType) == "interval") then

          desiredTimes(nDesire) = forcingTimeSinceEpoch(desiredRecordNumbers(nDesire), forcingField % positionInInterval, &
                                                        calc_interval(forcingField), desiredYears(nDesire))

       else if (trim(forcingField % forcingType) == "monthly") then

          desiredTimes(nDesire) = forcingTimeSinceEpochMonthly(desiredYears(nDesire), desiredRecordNumbers(nDesire), &
                                                               forcingField % positionInInterval)

       endif

    enddo ! nDesire

  end subroutine get_desired_record_positions

  !-------------------------------------------------------------

  function forcingTimeInYear(intervalNumber, positionInInterval, interval) result(forcingTime)

    integer, intent(in) :: intervalNumber

    real(kind=RKIND), intent(in) :: &
         positionInInterval, &
         interval

    real(kind=RKIND) :: forcingTime

    forcingTime = (real(intervalNumber-1,RKIND) + positionInInterval) * interval

  end function forcingTimeInYear

  !-------------------------------------------------------------

  function forcingTimeSinceEpoch(intervalNumber, positionInInterval, interval, year) result(forcingTime)

    integer, intent(in) :: &
         intervalNumber, &
         year

    real(kind=RKIND), intent(in) :: &
         positionInInterval, &
         interval

    real(kind=RKIND) :: forcingTime

    type (MPAS_Time_type) :: &
         startOfThisYear

    type (MPAS_TimeInterval_type) :: &
         startThisYearEpoch

    call mpas_set_time(startOfThisYear, YYYY=year, MM=1, DD=1)    

    startThisYearEpoch = startOfThisYear - epochTime
    call mpas_get_timeInterval(startThisYearEpoch, dt=forcingTime)

    forcingTime = forcingTime + forcingTimeInYear(intervalNumber, positionInInterval, interval)

  end function forcingTimeSinceEpoch

  !-------------------------------------------------------------

  function forcingTimeInYearMonthly(year, month, positionInInterval) result(forcingTime)
    
    integer, intent(in) :: &
         year, &
         month

    real(kind=RKIND), intent(in) :: &
         positionInInterval

    real(kind=RKIND) :: forcingTime

    type (MPAS_Time_type) :: &
         startOfThisYear, &
         thisForcingTime

    type (MPAS_TimeInterval_type) :: &
         forcingTimeInterval

    call mpas_set_time(startOfThisYear, YYYY=year, MM=1, DD=1)  
    call mpas_set_time(thisForcingTime, YYYY=year, MM=month, DD=15)

    forcingTimeInterval = thisForcingTime - startOfThisYear

    call mpas_get_timeInterval(forcingTimeInterval, dt=forcingTime)

  end function forcingTimeInYearMonthly

  !-------------------------------------------------------------

  function forcingTimeSinceEpochMonthly(year, month, positionInInterval) result(forcingTime)

    integer, intent(in) :: &
         year, &
         month

    real(kind=RKIND), intent(in) :: &
         positionInInterval

    real(kind=RKIND) :: forcingTime

    type (MPAS_Time_type) :: &
         thisForcingTime

    type (MPAS_TimeInterval_type) :: &
         forcingTimeInterval

    call mpas_set_time(thisForcingTime, YYYY=year, MM=month, DD=15)

    forcingTimeInterval = thisForcingTime - epochTime

    call mpas_get_timeInterval(forcingTimeInterval, dt=forcingTime)

  end function forcingTimeSinceEpochMonthly

  !-------------------------------------------------------------

  subroutine read_data(mesh, forcingDirectory, forcingField, slot, desiredYear, desiredRecordNumber, position)

    type(mesh_type), intent(in) :: mesh

    character(len=*), intent(in) :: forcingDirectory

    type (forcing_field), intent(inout) :: forcingField

    integer, intent(in) :: &
         slot, &
         desiredYear, &
         desiredRecordNumber, &
         position

    if (forcingField % filetype == trim("binaryPOP")) then
    
       call read_data_binaryPOP(mesh, forcingDirectory, forcingField, slot, desiredYear, desiredRecordNumber, position)

    endif

  end subroutine read_data

  !-------------------------------------------------------------

  subroutine read_data_binaryPOP(mesh, forcingDirectory, forcingField, slot, desiredYear, desiredRecordNumber, position)

    type(mesh_type), intent(in) :: mesh

    character(len=*), intent(in) :: forcingDirectory

    type (forcing_field), intent(inout) :: forcingField

    integer, intent(in) :: &
         slot, &
         desiredYear, &
         desiredRecordNumber, &
         position

    character(len=404) :: filename

    real(kind=RKIND), dimension(:,:), allocatable :: arrayPOP

    integer :: iCell, recordLength, usedRecordNumber

    logical :: lexist

    allocate(arrayPOP(1:mesh % POP_nx, 1:mesh % POP_ny))

    recordLength = mesh % POP_nx * mesh % POP_ny * 8

    usedRecordNumber = desiredRecordNumber

    call get_filename(filename, usedRecordNumber, forcingDirectory, forcingField, desiredYear, position)

    open(55,file=trim(filename),form='unformatted',access='direct',recl=recordLength,action='read',status='old')

    read(55,rec=usedRecordNumber) arrayPOP
 
    close(55)

    forcingField % data(:,slot) = 273.15_RKIND!1.0e30_RKIND
    
    do iCell = 1, mesh % nCells

       forcingField % data(iCell,slot) = arrayPOP(mesh % POPindxi % array(iCell),mesh % POPindxj % array(iCell))

       !write(*,*) iCell, mesh % POPindxi % array(iCell),mesh % POPindxj % array(iCell)

    enddo ! iCell

    deallocate(arrayPOP)

  end subroutine read_data_binaryPOP

  !-------------------------------------------------------------

  subroutine get_filename(filename, usedRecordNumber, forcingDirectory, forcingField, desiredYear, position)

    character(len=*), intent(out) :: filename

    integer, intent(out) :: usedRecordNumber

    character(len=*), intent(in) :: forcingDirectory

    type (forcing_field), intent(in) :: forcingField

    integer, intent(in) :: desiredYear

    integer, intent(in) :: position

    logical :: lexist

    write(filename,fmt='(a,a,a,i4.4,a)') trim(forcingDirectory), trim(forcingField % filenameDirectory), trim(forcingField % filenamePrefix), desiredYear, trim(forcingField % filenameSuffix)

    inquire(file=trim(filename),exist=lexist)

    if (.not. lexist) then

       if (forcingField % forcingType == "interval") then
          
          if (position == 1) then

             write(filename,fmt='(a,a,a,i4.4,a)') trim(forcingDirectory), trim(forcingField % filenameDirectory), trim(forcingField % filenamePrefix), desiredYear+1, trim(forcingField % filenameSuffix)
             
             usedRecordNumber = 1

          else if (position == 2) then

             write(filename,fmt='(a,a,a,i4.4,a)') trim(forcingDirectory), trim(forcingField % filenameDirectory), trim(forcingField % filenamePrefix), desiredYear-1, trim(forcingField % filenameSuffix)
             
             usedRecordNumber = forcingField % pointsPerYear

          endif

       else if (forcingField % forcingType == "monthly") then

          if (position == 1) then

             write(filename,fmt='(a,a,a,i4.4,a)') trim(forcingDirectory), trim(forcingField % filenameDirectory), trim(forcingField % filenamePrefix), desiredYear+1, trim(forcingField % filenameSuffix)
             
             usedRecordNumber = forcingField % pointsPerYear

          else if (position == 2) then

             write(filename,fmt='(a,a,a,i4.4,a)') trim(forcingDirectory), trim(forcingField % filenameDirectory), trim(forcingField % filenamePrefix), desiredYear-1, trim(forcingField % filenameSuffix)
             
             usedRecordNumber = 1

          endif

       endif

    endif

  end subroutine get_filename

  !-------------------------------------------------------------

  function POPiCell(mesh, i, j) result(iCellOut)

    type(mesh_type), intent(in) :: mesh
    
    integer, intent(in) :: i, j
    integer :: iCellOut

    integer :: iCell

    iCellOut = -1

    do iCell = 1, mesh % nCells

       if (mesh % POPindxi % array(iCell) == i .and. mesh % POPindxj % array(iCell) == j) then

          iCellOut = iCell
          exit

       endif

    enddo ! iCell

  end function POPiCell

  !-------------------------------------------------------------
  
  subroutine rotate_vector_POP_to_cartesian(uArrayCartesian,     vArrayCartesian, &
                                            uArrayPOP,           vArrayPOP,       &
                                            uPOPAlongGridVector, vPOPAlongGridVector)

    real(kind=RKIND), dimension(:), intent(out) :: &
         uArrayCartesian, &
         vArrayCartesian

    real(kind=RKIND), dimension(:), intent(in) :: &
         uArrayPOP, &
         vArrayPOP

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         uPOPAlongGridVector, &
         vPOPAlongGridVector

    integer :: &
         nPoints, &
         iPoint

    nPoints = size(uArrayCartesian)

    do iPoint = 1, nPoints

       uArrayCartesian(iPoint) = uArrayPOP(iPoint) * uPOPAlongGridVector(1,iPoint) + &
                                 vArrayPOP(iPoint) * vPOPAlongGridVector(1,iPoint)
       vArrayCartesian(iPoint) = uArrayPOP(iPoint) * uPOPAlongGridVector(2,iPoint) + &
                                 vArrayPOP(iPoint) * vPOPAlongGridVector(2,iPoint)

    enddo ! iPoint

  end subroutine rotate_vector_POP_to_cartesian

  !-------------------------------------------------------------

  subroutine POP_grid_basis(mesh, uPOPAlongGridVector, vPOPAlongGridVector)

    use cice_dynamics_shared, only: grid_rotation_forward, &
                                    grid_rotation_forward_nochange, &
                                    dot_product_3space

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         uPOPAlongGridVector, &
         vPOPAlongGridVector

    real(kind=RKIND) :: &
         xCell, &
         yCell, &
         zCell, &
         xEdge1, &
         yEdge1, &
         zEdge1, &
         xEdge2, &
         yEdge2, &
         zEdge2

    real(kind=RKIND), dimension(3) :: &
         eastwardsVector, &
         uVector, &
         vVector

    integer :: &
         iCell, &
         iPOPCell, &
         jPOPCell, &
         nEdgesOnCell, &
         iEdgeOnCell, &
         iCellAdjacent, &
         iPOPAdjacent, &
         jPOPAdjacent, &
         iEdgeOnCellU2, &
         iEdgeOnCellU1, &
         iEdgeOnCellV2, &
         iEdgeOnCellV1, &
         iEdgeU2, &
         iEdgeU1, &
         iEdgeV2, &
         iEdgeV1
    
    do iCell = 1, mesh % nCells

       ! eastwards vector
       call grid_rotation_forward(xCell, yCell, zCell, &
                                  mesh % xCell % array(iCell), mesh % yCell % array(iCell), mesh % zCell % array(iCell))

       eastwardsVector(1) = -yCell
       eastwardsVector(2) =  xCell
       eastwardsVector(3) = 0.0_RKIND
       eastwardsVector = eastwardsVector / &
                         sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)

       ! cell location in POP coordinates
       iPOPCell = mesh % POPindxi % array(iCell)
       jPOPCell = mesh % POPindxj % array(iCell)       

       nEdgesOnCell = mesh % nEdgesOnCell % array(iCell)

       iEdgeOnCellU2 = -1
       iEdgeOnCellU1 = -1
       iEdgeOnCellV2 = -1
       iEdgeOnCellV1 = -1

       do iEdgeOnCell = 1, nEdgesOnCell

          iCellAdjacent = mesh % cellsOnCell % array(iEdgeOnCell,iCell)

          if (iCellAdjacent < mesh % nCells + 1) then

             iPOPAdjacent = mesh % POPindxi % array(iCellAdjacent)
             jPOPAdjacent = mesh % POPindxj % array(iCellAdjacent)

             if (iPOPAdjacent == iPOPCell+1 .and. jPOPAdjacent == jPOPCell) then
                
                iEdgeOnCellU2 = iEdgeOnCell
                iEdgeOnCellU1 = iEdgeOnCell + 2
                iEdgeOnCellV2 = iEdgeOnCell + 1
                iEdgeOnCellV1 = iEdgeOnCell + 3
                exit

             else if (iPOPAdjacent == iPOPCell .and. jPOPAdjacent == jPOPCell+1) then

                iEdgeOnCellU2 = iEdgeOnCell + 3
                iEdgeOnCellU1 = iEdgeOnCell + 1
                iEdgeOnCellV2 = iEdgeOnCell
                iEdgeOnCellV1 = iEdgeOnCell + 2
                exit

             else if (iPOPAdjacent == iPOPCell-1 .and. jPOPAdjacent == jPOPCell) then

                iEdgeOnCellU2 = iEdgeOnCell + 2
                iEdgeOnCellU1 = iEdgeOnCell
                iEdgeOnCellV2 = iEdgeOnCell + 3
                iEdgeOnCellV1 = iEdgeOnCell + 1
                exit

             else if (iPOPAdjacent == iPOPCell .and. jPOPAdjacent == jPOPCell-1) then

                iEdgeOnCellU2 = iEdgeOnCell + 1
                iEdgeOnCellU1 = iEdgeOnCell + 3
                iEdgeOnCellV2 = iEdgeOnCell + 2
                iEdgeOnCellV1 = iEdgeOnCell
                exit

             endif

          endif ! iCellAdjacent exists

       enddo ! iEdgeOnCell

       if (iEdgeOnCellU2 > nEdgesOnCell) iEdgeOnCellU2 = iEdgeOnCellU2 - nEdgesOnCell
       if (iEdgeOnCellU1 > nEdgesOnCell) iEdgeOnCellU1 = iEdgeOnCellU1 - nEdgesOnCell
       if (iEdgeOnCellV2 > nEdgesOnCell) iEdgeOnCellV2 = iEdgeOnCellV2 - nEdgesOnCell
       if (iEdgeOnCellV1 > nEdgesOnCell) iEdgeOnCellV1 = iEdgeOnCellV1 - nEdgesOnCell

       if (iEdgeOnCellU2 == -1 .or. iEdgeOnCellU1 == -1 .or. iEdgeOnCellV2 == -1 .or. iEdgeOnCellV1 == -1) then

          write(*,*) "POP_grid_basis failure!"
          write(*,*) iEdgeOnCellU2, iEdgeOnCellU1, iEdgeOnCellV2, iEdgeOnCellV1
          stop

       endif

       iEdgeU2 = mesh % edgesOnCell % array(iEdgeOnCellU2,iCell)
       iEdgeU1 = mesh % edgesOnCell % array(iEdgeOnCellU1,iCell)
       iEdgeV2 = mesh % edgesOnCell % array(iEdgeOnCellV2,iCell)
       iEdgeV1 = mesh % edgesOnCell % array(iEdgeOnCellV1,iCell)

       !write(*,*) iCell, iEdgeOnCellU2, iEdgeOnCellU1, iEdgeOnCellV2, iEdgeOnCellV1, &
       !           iEdgeU2, iEdgeU1, iEdgeV2, iEdgeV1
       
       ! U vector
       call grid_rotation_forward(xEdge2, yEdge2, zEdge2, &
                                  mesh % xEdge % array(iEdgeU2), mesh % yEdge % array(iEdgeU2), mesh % zEdge % array(iEdgeU2))
       
       call grid_rotation_forward(xEdge1, yEdge1, zEdge1, &
                                  mesh % xEdge % array(iEdgeU1), mesh % yEdge % array(iEdgeU1), mesh % zEdge % array(iEdgeU1))
       
       uVector(1) = xEdge2 - xEdge1
       uVector(2) = yEdge2 - yEdge1
       uVector(3) = zEdge2 - zEdge1
       uVector = uVector / sqrt(uVector(1)**2 + uVector(2)**2 + uVector(3)**2)

       call dot_product_3space(uPOPAlongGridVector(1,iCell), &
                               uVector(1),         uVector(2),         uVector(3), &
                               eastwardsVector(1), eastwardsVector(2), eastwardsVector(3)) 

       uPOPAlongGridVector(2,iCell) = sign(1.0_RKIND,uVector(3))&
                                    * sqrt(1.0_RKIND - min(uPOPAlongGridVector(1,iCell),1.0_RKIND)**2)

       !write(*,*) iCell, uPOPAlongGridVector(1,iCell), uPOPAlongGridVector(2,iCell)
       
       ! V vector
       call grid_rotation_forward(xEdge2, yEdge2, zEdge2, &
                                  mesh % xEdge % array(iEdgeV2), mesh % yEdge % array(iEdgeV2), mesh % zEdge % array(iEdgeV2))
       
       call grid_rotation_forward(xEdge1, yEdge1, zEdge1, &
                                  mesh % xEdge % array(iEdgeV1), mesh % yEdge % array(iEdgeV1), mesh % zEdge % array(iEdgeV1))
       
       vVector(1) = xEdge2 - xEdge1
       vVector(2) = yEdge2 - yEdge1
       vVector(3) = zEdge2 - zEdge1
       vVector = vVector / sqrt(vVector(1)**2 + vVector(2)**2 + vVector(3)**2)

       call dot_product_3space(vPOPAlongGridVector(1,iCell), &
                               vVector(1),         vVector(2),         vVector(3), &
                               eastwardsVector(1), eastwardsVector(2), eastwardsVector(3)) 

       vPOPAlongGridVector(2,iCell) = sign(1.0_RKIND,vVector(3))&
                                    * sqrt(1.0_RKIND - min(vPOPAlongGridVector(1,iCell),1.0_RKIND)**2)

    enddo ! iCell

    !stop

  end subroutine POP_grid_basis

  !-------------------------------------------------------------

end module cice_forcing
