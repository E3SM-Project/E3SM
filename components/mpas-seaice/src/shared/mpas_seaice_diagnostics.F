!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_diagnostics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_diagnostics

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_io_units
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       seaice_initialize_time_diagnostics, &
       seaice_set_time_diagnostics, &
       seaice_check_state, &
       seaice_set_testing_system_test_arrays, &
       seaice_load_balance_timers

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_initialize_time_diagnostics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_initialize_time_diagnostics(domain)!{{{

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         diagnosticsPool

    character(len=strKIND), pointer :: &
         xtime, &
         simulationStartTime

    logical, pointer :: &
         config_do_restart

    type(MPAS_Time_Type) :: &
         startTime

    integer :: &
         ierr

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)

       ! current time
       call MPAS_pool_get_array(diagnosticsPool, "xtime", xtime)
       startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr)
       call MPAS_get_time(startTime, dateTimeString=xtime)

       ! simulation start time
       call MPAS_pool_get_config(block % configs, "config_do_restart", config_do_restart)

       if (.not. config_do_restart) then
          call MPAS_pool_get_array(diagnosticsPool, "simulationStartTime", simulationStartTime)
          simulationStartTime = xtime
       endif

       block => block % next
    end do

  end subroutine seaice_initialize_time_diagnostics!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_initialize_time_diagnostics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_set_time_diagnostics(domain)!{{{

    use seaice_constants, only: &
         seaiceDaysPerSecond

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         diagnosticsPool

    character(len=strKIND), pointer :: &
         xtime, &
         simulationStartTime

    type(MPAS_Time_Type) :: &
         currTime, &
         xtime_timeType, &
         simulationStartTime_timeType

    real(kind=RKIND), pointer :: &
         daysSinceStartOfSim

    integer :: &
         ierr

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)

       ! set xtime
       call MPAS_pool_get_array(diagnosticsPool, "xtime", xtime)
       currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
       call mpas_get_time(curr_time=currTime, dateTimeString=xtime)

       ! compute time since start of simulation, in days
       call mpas_pool_get_array(diagnosticsPool, 'simulationStartTime', simulationStartTime)
       call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim', daysSinceStartOfSim)
       call mpas_set_time(xtime_timeType, dateTimeString=xtime)
       call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
       call mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType,dt=daysSinceStartOfSim)
       daysSinceStartOfSim = daysSinceStartOfSim*seaiceDaysPerSecond

       block => block % next
    end do

  end subroutine seaice_set_time_diagnostics!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_check_state
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 10th March 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_check_state(domain)

    use seaice_constants, only: &
         seaiceRadiansToDegrees

    type (domain_type), intent(in) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         diagnosticsPool, &
         meshPool, &
         tracersPool, &
         velocitySolverPool

    integer, dimension(:), pointer :: &
         indexToCellID, &
         indexToVertexID

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         lonCell, &
         latVertex, &
         lonVertex, &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceSalinity

    integer, pointer :: &
         nCellsSolve, &
         nVerticesSolve, &
         nCategories

    integer :: &
         iCell, &
         iVertex, &
         iCategory

    character(len=strKIND), pointer :: &
         xtime

    character(len=strKIND) :: &
         filename

    integer :: &
         errorUnit, &
         iostat

    logical :: &
         errorFlag

    logical, pointer :: &
         config_check_state

    ! limit values
    real(kind=RKIND), parameter :: &
         iceVolumeLimit        = 50.0_RKIND, &
         snowVolumeLimit       = 20.0_RKIND, &
         iceSalinityLowerLimit =  0.0_RKIND, &
         iceSalinityUpperLimit = 50.0_RKIND, &
         velocityLimit         =  5.0_RKIND

    call MPAS_pool_get_config(domain % configs, "config_check_state", config_check_state)

    if (config_check_state) then

       errorFlag = .false.

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)

          call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
          call MPAS_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)
          call MPAS_pool_get_array(tracersPool, "iceSalinity", iceSalinity, 1)

          call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

          call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
          call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
          call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)

          ! cell centre checks
          do iCell = 1, nCellsSolve

             do iCategory = 1, nCategories

                ! seaice thickness check
                if (iceVolumeCategory(1,iCategory,iCell) > iceVolumeLimit) errorFlag = .true.

                ! snow on seaice thickness check
                if (snowVolumeCategory(1,iCategory,iCell) > snowVolumeLimit) errorFlag = .true.

                ! salinity check
                if (iceSalinity(1,iCategory,iCell) < iceSalinityLowerLimit .or. &
                    iceSalinity(1,iCategory,iCell) > iceSalinityUpperLimit) errorFlag = .true.

             enddo ! iCategory

          enddo ! iCell

          ! vertex checks
          do iVertex = 1, nVerticesSolve

             ! velocity check
             if (uVelocity(iVertex) > velocityLimit .or. vVelocity(iVertex) > velocityLimit) errorFlag = .true.

          enddo ! iVertex

          block => block % next
       enddo

       if (errorFlag) then

          ! open the error file
          call mpas_new_unit(errorUnit)
          write(filename,fmt='(a,i6.6,a)') "mpas_seaice_state_test_", domain % dminfo % my_proc_id, ".log"
          open(errorUnit, file=trim(filename), position="append", iostat=iostat)
          if (iostat /= 0) then
             call mpas_log_write("seaice_check_state: problem opening state check: $i", &
                  messageType=MPAS_LOG_CRIT, intArgs=(/iostat/))
          endif

          block => domain % blocklist
          do while (associated(block))

             call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)
             call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
             call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
             call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)

             call MPAS_pool_get_array(meshPool, "indexToCellID", indexToCellID)
             call MPAS_pool_get_array(meshPool, "indexToVertexID", indexToVertexID)

             call MPAS_pool_get_array(meshPool, "latCell", latCell)
             call MPAS_pool_get_array(meshPool, "lonCell", lonCell)
             call MPAS_pool_get_array(meshPool, "latVertex", latVertex)
             call MPAS_pool_get_array(meshPool, "lonVertex", lonVertex)

             call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
             call MPAS_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)
             call MPAS_pool_get_array(tracersPool, "iceSalinity", iceSalinity, 1)

             call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
             call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

             call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
             call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
             call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)

             call MPAS_pool_get_array(diagnosticsPool, "xtime", xtime)

             ! cell centre checks
             do iCell = 1, nCellsSolve

                do iCategory = 1, nCategories

                   ! seaice thickness check
                   if (iceVolumeCategory(1,iCategory,iCell) > iceVolumeLimit) then

                      write(errorUnit,fmt=10) "Time=", trim(xtime), ", iCell=", indexToCellID(iCell), ", iCategory=", iCategory, &
                           ", lat=", latCell(iCell)*seaiceRadiansToDegrees, ", lon=", lonCell(iCell)*seaiceRadiansToDegrees, &
                           ", iceVolumeCategory=", iceVolumeCategory(1,iCategory,iCell)

                   endif ! seaice thickness check

                   ! snow on seaice thickness check
                   if (snowVolumeCategory(1,iCategory,iCell) > snowVolumeLimit) then

                      write(errorUnit,fmt=10) "Time=", trim(xtime), ", iCell=", indexToCellID(iCell), ", iCategory=", iCategory, &
                           ", lat=", latCell(iCell)*seaiceRadiansToDegrees, ", lon=", lonCell(iCell)*seaiceRadiansToDegrees, &
                           ", snowVolumeCategory=", snowVolumeCategory(1,iCategory,iCell)

                   endif ! snow on seaice thickness check

                   ! salinity check
                   if (iceSalinity(1,iCategory,iCell) < iceSalinityLowerLimit .or. &
                       iceSalinity(1,iCategory,iCell) > iceSalinityUpperLimit) then

                      write(errorUnit,fmt=10) "Time=", trim(xtime), ", iCell=", indexToCellID(iCell), ", iCategory=", iCategory, &
                           ", lat=", latCell(iCell)*seaiceRadiansToDegrees, ", lon=", lonCell(iCell)*seaiceRadiansToDegrees, &
                           ", iceSalinity=", iceSalinity(1,iCategory,iCell)

                   endif ! salinity check

                enddo ! iCategory

             enddo ! iCell

             ! vertex checks
             do iVertex = 1, nVerticesSolve

                ! velocity check
                if (uVelocity(iVertex) >= velocityLimit .or. vVelocity(iVertex) >= velocityLimit) then

                   write(errorUnit,fmt=20) "Time=", trim(xtime), ", iVertex=", indexToVertexID(iVertex), &
                        ", lat=", latVertex(iVertex)*seaiceRadiansToDegrees, ", lon=", lonVertex(iVertex)*seaiceRadiansToDegrees, &
                        ", uVelocity=", uVelocity(iVertex), ", vVelocity=", vVelocity(iVertex)

                endif ! velocity check

             enddo ! iVertex

             block => block % next
          enddo

          close(errorUnit)
          call mpas_release_unit(errorUnit)

       endif ! errorFlag

    endif ! config_check_state

    ! output formats
10  format(a,a,a,i9,a,i5,a,e12.4,a,e12.4,a,e12.4)
20  format(a,a,a,i9,a,e12.4,a,e12.4,a,e12.4,a,e12.4)

  end subroutine seaice_check_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_set_testing_system_test_arrays
!
!> \brief Set up arrays for testing the testing system
!> \author Adrian K. Turner, LANL
!> \date 11th May 2017
!> \details
!> This routine sets up the testing system test arrays that will cause
!> failures in each test.
!
!-----------------------------------------------------------------------

  subroutine seaice_set_testing_system_test_arrays(domain)

    type(domain_type), intent(in) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         testingSystemTestPool

    real(kind=RKIND), dimension(:), pointer :: &
         testArrayRegression, &
         testArrayParallelism, &
         testArrayRestartability

    logical, pointer :: &
         config_testing_system_test, &
         config_do_restart

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    ! random number variables
    integer :: &
         seedSize, &
         clock, &
         i

    integer, dimension(:), allocatable :: &
         seed

    call MPAS_pool_get_config(domain % configs, "config_testing_system_test", config_testing_system_test)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (config_testing_system_test) then

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

          call MPAS_pool_get_subpool(block % structs, "testing_system_test", testingSystemTestPool)

          ! regression
          call random_seed(size=seedSize)
          allocate(seed(seedSize))

          call system_clock(count=clock)

          seed = clock + 37 * (/ (i - 1, i = 1, seedSize) /)
          call random_seed(PUT=seed)

          deallocate(seed)

          call MPAS_pool_get_array(testingSystemTestPool, "testArrayRegression", testArrayRegression)

          do iCell = 1, nCellsSolve
             call random_number(testArrayRegression(iCell))
          enddo ! iCell

          ! parallelism
          call MPAS_pool_get_array(testingSystemTestPool, "testArrayParallelism", testArrayParallelism)
          testArrayParallelism = real(domain % dminfo % my_proc_id, RKIND)

          ! restartability
          call MPAS_pool_get_array(testingSystemTestPool, "testArrayRestartability", testArrayRestartability)
          if (config_do_restart) then
             testArrayRestartability = 1.0_RKIND
          endif

          block => block % next
       enddo

    endif

  end subroutine seaice_set_testing_system_test_arrays

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_load_balance_timers
!
!> \brief Time load in-balance
!> \author Adrian K. Turner, LANL
!> \date 12th December 2017
!> \details
!  Add a call to MPI_Barrier to allow timing of load imbalance
!
!-----------------------------------------------------------------------

  subroutine seaice_load_balance_timers(domain, label)

#ifdef _MPI
#ifndef NOMPIMOD
   use mpi
#endif
#endif
   use mpas_timer

    type(domain_type), intent(in) :: &
         domain

    character(len=*), intent(in) :: &
         label

#ifdef _MPI

    logical, pointer :: &
         config_load_balance_timers

    integer :: &
         mpi_ierr

    call MPAS_pool_get_config(domain % configs, "config_load_balance_timers", config_load_balance_timers)
    if (config_load_balance_timers) then
       call mpas_timer_start("halo barrier "//trim(label))
       call MPI_Barrier(domain % dminfo % comm, mpi_ierr)
       call mpas_timer_stop("halo barrier "//trim(label))
    endif

#endif

  end subroutine seaice_load_balance_timers

  !-----------------------------------------------------------------------

end module seaice_diagnostics
