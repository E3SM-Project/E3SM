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
       cice_check_state
       
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

    type (domain_type), intent(in) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         diagnosticsPool, &
         meshPool, &
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

    integer, pointer :: &
         nCellsSolve, &
         nVerticesSolve

    integer :: &
         iCell, &
         iVertex

    character(len=strKIND), pointer :: &
         xtime

    character(len=strKIND) :: &
         filename

    integer :: &
         errorUnit

    logical :: &
         errorFlag

    logical, pointer :: &
         config_check_state

    ! limit values
    real(kind=RKIND), parameter :: &
         velocityLimit = 0.0_RKIND

    call MPAS_pool_get_config(domain % configs, "config_check_state", config_check_state)

    if (config_check_state) then

       errorFlag = .false.

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

          call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

          call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
          call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

          ! cell centre checks
          !do iCell = 1, nCellsSolve



          !enddo ! iCell

          ! vertex checks
          do iVertex = 1, nVerticesSolve

             ! velocity check
             if (uVelocity(iVertex) >= velocityLimit .or. vVelocity(iVertex) >= velocityLimit) errorFlag = .true.

          enddo ! iVertex

          block => block % next
       enddo

       if (errorFlag) then

          ! open the error file
          call mpas_new_unit(errorUnit)
          write(filename,fmt='(a,i6.6,a)') "mpas_seaice_state_test_", domain % dminfo % my_proc_id, ".log"
          open(errorUnit, file=trim(filename))

          block => domain % blocklist
          do while (associated(block))

             call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)
             call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
             call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

             call MPAS_pool_get_array(meshPool, "indexToCellID", indexToCellID)
             call MPAS_pool_get_array(meshPool, "indexToVertexID", indexToVertexID)

             call MPAS_pool_get_array(meshPool, "latCell", latCell)
             call MPAS_pool_get_array(meshPool, "lonCell", lonCell)
             call MPAS_pool_get_array(meshPool, "latVertex", latVertex)
             call MPAS_pool_get_array(meshPool, "lonVertex", lonVertex)

             call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
             call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

             call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
             call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

             call MPAS_pool_get_array(diagnosticsPool, "xtime", xtime)

             ! cell centre checks
             !do iCell = 1, nCellsSolve



             !enddo ! iCell

             ! vertex checks
             do iVertex = 1, nVerticesSolve

                ! velocity check
                if (uVelocity(iVertex) >= velocityLimit .or. vVelocity(iVertex) >= velocityLimit) then

                   write(errorUnit,fmt=10) "Time=", trim(xtime), ", iVertex=", indexToVertexID(iVertex), &
                        ", lat=", latVertex(iVertex), ", lon=", lonVertex(iVertex), &
                        ", uVelocity=", uVelocity(iVertex), ", vVelocity=", vVelocity(iVertex)

                   10 format(a,a,a,i9,a,e12.4,a,e12.4,a,e12.4,a,e12.4)

                endif ! velocity check

             enddo ! iVertex

             block => block % next
          enddo

          close(errorUnit)
          call mpas_release_unit(errorUnit)

       endif ! errorFlag

    endif ! config_check_state

  end subroutine cice_check_state

!-----------------------------------------------------------------------

end module cice_diagnostics
