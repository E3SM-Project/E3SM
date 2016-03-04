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
       cice_set_time_diagnostics
       
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

!-----------------------------------------------------------------------

end module cice_diagnostics
