!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_diagnostics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_diagnostics

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping

  implicit none

  private
  save

  public :: &
       cice_initialize_time_diagnostics, &
       cice_set_time_diagnostics, &
       cice_check_state, &
       cice_warning_message

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_initialize_time_diagnostics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_initialize_time_diagnostics(domain)!{{{

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

  end subroutine cice_initialize_time_diagnostics!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_initialize_time_diagnostics
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_set_time_diagnostics(domain)!{{{

    use cice_constants, only: &
         ciceDaysPerSecond

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
       daysSinceStartOfSim = daysSinceStartOfSim*ciceDaysPerSecond

       block => block % next
    end do

  end subroutine cice_set_time_diagnostics!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_check_state
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 10th March 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_check_state(domain)

    use cice_constants, only: &
         ciceRadiansToDegrees

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
             write(stderrUnit,*) "cice_check_state: problem opening state check: ", iostat
             call MPAS_dmpar_global_abort("MPAS-seaice: problem opening state check")
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
                           ", lat=", latCell(iCell)*ciceRadiansToDegrees, ", lon=", lonCell(iCell)*ciceRadiansToDegrees, &
                           ", iceVolumeCategory=", iceVolumeCategory(1,iCategory,iCell)

                   endif ! seaice thickness check

                   ! snow on seaice thickness check
                   if (snowVolumeCategory(1,iCategory,iCell) > snowVolumeLimit) then

                      write(errorUnit,fmt=10) "Time=", trim(xtime), ", iCell=", indexToCellID(iCell), ", iCategory=", iCategory, &
                           ", lat=", latCell(iCell)*ciceRadiansToDegrees, ", lon=", lonCell(iCell)*ciceRadiansToDegrees, &
                           ", snowVolumeCategory=", snowVolumeCategory(1,iCategory,iCell)

                   endif ! snow on seaice thickness check

                   ! salinity check
                   if (iceSalinity(1,iCategory,iCell) < iceSalinityLowerLimit .or. &
                       iceSalinity(1,iCategory,iCell) > iceSalinityUpperLimit) then

                      write(errorUnit,fmt=10) "Time=", trim(xtime), ", iCell=", indexToCellID(iCell), ", iCategory=", iCategory, &
                           ", lat=", latCell(iCell)*ciceRadiansToDegrees, ", lon=", lonCell(iCell)*ciceRadiansToDegrees, &
                           ", iceSalinity=", iceSalinity(1,iCategory,iCell)

                   endif ! salinity check

                enddo ! iCategory

             enddo ! iCell

             ! vertex checks
             do iVertex = 1, nVerticesSolve

                ! velocity check
                if (uVelocity(iVertex) >= velocityLimit .or. vVelocity(iVertex) >= velocityLimit) then

                   write(errorUnit,fmt=20) "Time=", trim(xtime), ", iVertex=", indexToVertexID(iVertex), &
                        ", lat=", latVertex(iVertex)*ciceRadiansToDegrees, ", lon=", lonVertex(iVertex)*ciceRadiansToDegrees, &
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

  end subroutine cice_check_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_warning_message
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 13th May 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_warning_message(domain, message)

    type(domain_type), intent(in) :: &
         domain

    character(len=*), intent(in) :: &
         message

    integer :: &
         errorUnit, &
         iostat

    character(len=strKIND) :: &
         filename

    type(MPAS_pool_type), pointer :: &
         diagnosticsPool

    character(len=strKIND), pointer :: &
         xtime

    call mpas_new_unit(errorUnit)
    write(filename,fmt='(a,i6.6,a)') "mpas_seaice_state_test_", domain % dminfo % my_proc_id, ".log"
    open(errorUnit, file=trim(filename), position="append", iostat=iostat)
    if (iostat /= 0) then
       write(stderrUnit,*) "cice_warning_message: Problem opening state check file", iostat
       call MPAS_dmpar_global_abort("MPAS-seaice: Problem opening state check file")
    endif

    call MPAS_pool_get_subpool(domain % blocklist % structs, "diagnostics", diagnosticsPool)
    call MPAS_pool_get_array(diagnosticsPool, "xtime", xtime)

    write(errorUnit,fmt='(a,a,a,a)') "Time=", trim(xtime), ", ", trim(message)

    close(errorUnit)
    call mpas_release_unit(errorUnit)

  end subroutine cice_warning_message

!-----------------------------------------------------------------------

end module cice_diagnostics
