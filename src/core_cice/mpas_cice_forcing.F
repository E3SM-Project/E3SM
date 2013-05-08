module cice_forcing

  use mpas_grid_types

  implicit none

  private

  public :: &
       init_forcing, &
       get_forcing_test

  real(kind=RKIND), parameter :: &
       seconds_in_minute = 60.0_RKIND, &
       seconds_in_hour = 60.0_RKIND * seconds_in_minute, &
       seconds_in_day = 24.0_RKIND * seconds_in_hour

  ! individual forcing field derived type
  type forcing_field

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

  type(forcing_field) :: airTemperatureForcing


contains

  !-------------------------------------------------------------

  subroutine init_forcing(mesh)

    type(mesh_type), intent(in) :: mesh

    call init_forcing_field(mesh, airTemperatureForcing, &
         "/Users/akt/Work/MPAS-CICE/forcing/DATA/atm/gx3v2m/NCAR_bulk/NCEP/4XDAILY/STATES/t_10.", &
         ".dat", "binaryPOP", 1460, 0.0_RKIND)

  end subroutine init_forcing

  !-------------------------------------------------------------

  subroutine init_forcing_field(mesh, forcingField, &
                                filenamePrefix, filenameSuffix, filetype, &
                                pointsPerYear, positionInInterval)

    type(mesh_type), intent(in) :: mesh
    type (forcing_field), intent(inout) :: forcingField

    character(len=*), intent(in) :: filenamePrefix, filenameSuffix, filetype
    integer, intent(in) :: pointsPerYear
    real(kind=RKIND) :: positionInInterval

    forcingField % filenamePrefix = trim(filenamePrefix)
    forcingField % filenameSuffix = trim(filenameSuffix)

    forcingField % positionInInterval = positionInInterval
    forcingField % filetype           = trim(filetype)
    forcingField % pointsPerYear      = pointsPerYear

    forcingField % dataTimes(:) = -1.0e30_RKIND
    forcingField % dataRecordNumbers(:) = -1
    forcingField % dataYears(:) = 0    

    forcingField % dataOrder(1) = 1
    forcingField % dataOrder(2) = 2
    
    forcingField % currentYear = 0  

    forcingField % interval = (days_in_year() * seconds_in_day) / real(forcingField % pointsPerYear,RKIND)
    
    forcingField % ndata = mesh % nCells

    allocate(forcingField % data(1:forcingField % ndata,1:2))

  end subroutine init_forcing_field

  !-------------------------------------------------------------

  subroutine get_forcing_test(mesh, clock, atmos_forcing)

    use mpas_timekeeping

    type(mesh_type), intent(in) :: mesh
    type (MPAS_Clock_type), intent(in) :: clock
    type (atmos_forcing_type), pointer :: atmos_forcing

    call get_forcing(mesh, clock, airTemperatureForcing, atmos_forcing % airTemperature % array)

  end subroutine get_forcing_test

  !-------------------------------------------------------------

  subroutine get_forcing(mesh, clock, forcingField, forcingArray)

    use mpas_timekeeping

    type(mesh_type), intent(in) :: mesh

    type (MPAS_Clock_type), intent(in) :: clock
    type (forcing_field), intent(inout) :: forcingField

    real(kind=RKIND), dimension(:), intent(out) :: forcingArray

    type (MPAS_Time_Type) :: currTime
    integer :: ierr

    integer, dimension(2) :: desiredYears
    integer, dimension(2) :: desiredRecordNumbers
    real(kind=RKIND), dimension(2) :: desiredTimes
    
    integer, dimension(2) :: &
         existingDataLocation

    logical, dimension(2) :: &
         existingHasData

    real(kind=RKIND), dimension(2) :: &
         interpolationCoeff

    integer :: YYYY, DoY, H, M, S, S_n, S_d

    real(kind=RKIND) :: &
         forcingTime, &
         forcingTime1900

    integer :: &
         nDesired, &
         nExist

    ! get the current time
    currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime, YYYY=YYYY, DoY=DoY, H=H, M=M, S=S, S_n=S_n, S_d=S_d, ierr=ierr)

    call get_forcing_time(forcingTime, DoY, H, M, S, S_n, S_d)
    call get_forcing_time_1900(forcingTime1900, YYYY, DoY, H, M, S, S_n, S_d)

    call get_desired_record_positions(forcingTime, YYYY, forcingField, desiredYears, desiredRecordNumbers, desiredTimes)
    
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

    !write(*,*) existingDataLocation, existingHasData
    !write(*,*) existingHasData

    do nDesired = 1, 2

       if (existingDataLocation(nDesired) == 0) then

          ! where do it read it into?
          if (.not. existingHasData(1)) then

             ! read into 1
             !write(*,*) "r:1"

             call read_data(mesh, forcingField, 1, desiredYears(nDesired), desiredRecordNumbers(nDesired))

             forcingField % dataTimes(1)         = desiredTimes(nDesired)
             forcingField % dataRecordNumbers(1) = desiredRecordNumbers(nDesired)
             forcingField % dataYears(1)         = desiredYears(nDesired)
             forcingField % dataOrder(nDesired) = 1

             existingHasData(1) = .true.

          else if (.not. existingHasData(2)) then

             ! read into 2
             !write(*,*) "r:2"

             call read_data(mesh, forcingField, 2, desiredYears(nDesired), desiredRecordNumbers(nDesired))

             forcingField % dataTimes(2)         = desiredTimes(nDesired)
             forcingField % dataRecordNumbers(2) = desiredRecordNumbers(nDesired)
             forcingField % dataYears(2)         = desiredYears(nDesired)
             forcingField % dataOrder(nDesired) = 2

             existingHasData(2) = .true.

          else

             ! major error

          endif

       endif

    enddo ! nDesired

    ! interpolate

    interpolationCoeff(1) = (desiredTimes(2) - forcingTime1900) / (desiredTimes(2) - desiredTimes(1))
    interpolationCoeff(2) = (forcingTime1900 - desiredTimes(1)) / (desiredTimes(2) - desiredTimes(1))

    forcingArray = forcingField % data(:,forcingField % dataOrder(1)) * interpolationCoeff(1) + &
                   forcingField % data(:,forcingField % dataOrder(2)) * interpolationCoeff(2)

    write(56,*) nint(forcingTime/3600.0_RKIND), &
         forcingField % data(mesh % cellID_ij % array(26,49),forcingField % dataOrder(1)), &
         forcingField % data(mesh % cellID_ij % array(26,49),forcingField % dataOrder(2))


  end subroutine get_forcing

  !-------------------------------------------------------------

  subroutine get_forcing_time(forcingTime, DoY, H, M, S, S_n, S_d) !!! need proper version

    real(kind=RKIND), intent(out) :: forcingTime

    integer, intent(in) :: &
         DoY, H, M, S, S_n, S_d

    forcingTime = real(DoY-1,RKIND) * seconds_in_day + &
                  real(H,RKIND) * seconds_in_hour + &
                  real(S,RKIND) * seconds_in_minute + &
                  real(S_n,RKIND) / real(S_d,RKIND)

  end subroutine get_forcing_time

  !-------------------------------------------------------------

  subroutine get_forcing_time_1900(forcingTime, YYYY, DoY, H, M, S, S_n, S_d) !!! need proper version

    real(kind=RKIND), intent(out) :: forcingTime

    integer, intent(in) :: &
         YYYY, DoY, H, M, S, S_n, S_d

    forcingTime = real(YYYY-1900,RKIND) * days_in_year() * seconds_in_day + &
                  real(DoY-1,RKIND) * seconds_in_day + &
                  real(H,RKIND) * seconds_in_hour + &
                  real(S,RKIND) * seconds_in_minute + &
                  real(S_n,RKIND) / real(S_d,RKIND)

  end subroutine get_forcing_time_1900

  !-------------------------------------------------------------
  
  function days_in_year() result(days) !!! need proper version

    real(kind=RKIND) :: days

    days = 365.0_RKIND

  end function days_in_year

  !-------------------------------------------------------------

  subroutine get_desired_record_positions(secOfYear, year, forcingField, desiredYears, desiredRecordNumbers, desiredTimes)

    real(kind=RKIND), intent(in) :: secOfYear

    integer, intent(in) :: year
   
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
         
    intervalNumber = ceiling(secOfYear / forcingField % interval)

    dataPosition = forcingTimeInYear(intervalNumber, forcingField % positionInInterval, forcingField % interval)

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

       desiredTimes(nDesire) = forcingTimeSince1900(desiredRecordNumbers(nDesire), forcingField % positionInInterval, &
                                                    forcingField % interval, desiredYears(nDesire))

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

  function forcingTimeSince1900(intervalNumber, positionInInterval, interval, year) result(forcingTime)

    integer, intent(in) :: &
         intervalNumber, &
         year

    real(kind=RKIND), intent(in) :: &
         positionInInterval, &
         interval

    real(kind=RKIND) :: forcingTime

    forcingTime = (real(intervalNumber-1,RKIND) + positionInInterval) * interval

    forcingTime = forcingTime + real(year - 1900,RKIND) * days_in_year() * seconds_in_day

  end function forcingTimeSince1900

  !-------------------------------------------------------------

  subroutine read_data(mesh, forcingField, position, desiredYear, desiredRecordNumber)

    type(mesh_type), intent(in) :: mesh

    type (forcing_field), intent(inout) :: forcingField

    integer, intent(in) :: &
         position, &
         desiredYear, &
         desiredRecordNumber

    if (forcingField % filetype == trim("binaryPOP")) then
    
       call read_data_binaryPOP(mesh, forcingField, position, desiredYear, desiredRecordNumber)

    endif

  end subroutine read_data

  !-------------------------------------------------------------

  subroutine read_data_binaryPOP(mesh, forcingField, position, desiredYear, desiredRecordNumber)

    type(mesh_type), intent(in) :: mesh

    type (forcing_field), intent(inout) :: forcingField

    integer, intent(in) :: &
         position, &
         desiredYear, &
         desiredRecordNumber

    character(len=404) :: filename

    real(kind=RKIND), dimension(:,:), allocatable :: arrayPOP

    integer :: i, j, recordLength

    allocate(arrayPOP(1:mesh % POP_nx, 1:mesh % POP_ny))

    write(filename,fmt='(a,i4.4,a)') trim(forcingField % filenamePrefix), desiredYear, trim(forcingField % filenameSuffix)

    recordLength = mesh % POP_nx * mesh % POP_ny * 8

    open(55,file=trim(filename),access='direct',form='unformatted',recl=recordLength)

    read(55,rec=desiredRecordNumber) arrayPOP
 
    close(55)

    !write(*,*) minval(arrayPOP), maxval(arrayPOP), arrayPOP(49,26), mesh % cellID_ij % array(49,26)

    forcingField % data(:,position) = 273.15_RKIND!1.0e30_RKIND
    
    do j = 1, mesh % POP_ny
    do i = 1, mesh % POP_nx

       forcingField % data(mesh % cellID_ij % array(i,j),position) = arrayPOP(i,j)

    enddo ! i
    enddo ! j

    deallocate(arrayPOP)

  end subroutine read_data_binaryPOP

  !-------------------------------------------------------------

end module cice_forcing
