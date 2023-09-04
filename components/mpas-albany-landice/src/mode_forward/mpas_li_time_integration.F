! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_time_integration
!
!> \brief MPAS land ice time integration driver
!> \author Matt Hoffman
!> \date   17 April 2011
!> \details
!>  This module contains the main driver routine for calling
!>  time integration schemes
!
!-----------------------------------------------------------------------

module li_time_integration

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_log

   use li_time_integration_fe_rk
   use li_setup
   use li_constants

   implicit none
   private

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: li_timestep

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************
   contains


!***********************************************************************
!
!  routine li_timestep
!
!> \brief   Advance model state forward in time by the specified time step
!> \author  Matt Hoffman
!> \date    20 April 2012
!> \details
!>  This routine advances model state forward in time by the specified time step.
!> Input: domain - current model state in time level 1 plus mesh data
!> Output: domain - upon exit, time level 2 contains
!>                  model state advanced forward in time by dt seconds
!-----------------------------------------------------------------------
   subroutine li_timestep(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      ! Pools pointers
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      character (len=StrKIND), pointer :: xtime
      real (kind=RKIND), pointer :: daysSinceStart
      character (len=StrKIND), pointer :: simulationStartTime
      character (len=StrKIND), pointer :: config_time_integration
      logical, pointer :: config_adaptive_timestep
      real (kind=RKIND), pointer :: deltat !< the current time step in seconds, on each block

      ! Other local variables
      type (MPAS_TimeInterval_type) :: timeStepInterval !< the current time step as an interval
      real (kind=RKIND) :: dtSeconds !< the current time step in seconds, local variable
      type (MPAS_Time_Type) :: currTime  !< current time as time type
      character(len=StrKIND) :: timeStamp !< current time as a string
      type (MPAS_Time_type) :: simulationStartTime_timeType

      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_time_integration', config_time_integration)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, err_tmp)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)

      ! ===
      ! === Non-adaptive timestep: Get dt in seconds
      ! ===
      if (.not. config_adaptive_timestep) then
         ! Get the interval at this point in time - will be fixed for nonadaptive timestep, but need to get it out of the clock
         timeStepInterval = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
         err = ior(err,err_tmp)
         ! Convert the clock's time interval into a dt in seconds to be used by the time stepper,
         ! using the currTime as the start time for this interval.
         ! (We want to do this conversion before doint the timestep and advancing the clock because the dt
         !  in seconds may change as the base time changes, and we want the old time as the base time.
         !  For example, the number of seconds in a year will be longer in a leap year.
         !  That is why nonadaptive timesteps have to be handled before the timestep, while
         !  adaptive timesteps are handled after the timestep.
         !  It may be possible to have them handled in the same place within li_tendency.F if we want to embed it that deeply.)
         call mpas_get_timeInterval(timeStepInterval, StartTimeIn=currTime, dt=dtSeconds, ierr=err_tmp)
         err = ior(err,err_tmp)
      else
         ! initialize the dt to 0 when using the adaptive timestepper
         !   it will be set by the time stepper but needs to avoid triggering an error in tend_layerThickness_fo_upwind
         dtSeconds = 0.0_RKIND
      endif

      ! Assign current value to variabels in all blocks
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         deltat = dtSeconds

         block => block % next
      end do


      ! ===
      ! === Perform timestep
      ! ===
      !call mpas_log_write('Using ' // trim(config_time_integration) // ' time integration.')
      select case (config_time_integration)
      case ('forward_euler')
         call li_time_integrator_forwardeuler_rungekutta(domain, err_tmp)
      case ('runge_kutta')
         call li_time_integrator_forwardeuler_rungekutta(domain, err_tmp)
      case default
         call mpas_log_write(trim(config_time_integration) // ' is not a valid land ice time integration option.', MPAS_LOG_ERR)
         err_tmp = 1
      end select
      err = ior(err,err_tmp)


      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_timestep.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
   end subroutine li_timestep


end module li_time_integration
