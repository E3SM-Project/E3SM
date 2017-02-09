! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration
!
!> \brief MPAS ocean time integration driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the time integration scheme
!
!-----------------------------------------------------------------------

module ocn_time_integration

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_io_units

   use ocn_constants
   use ocn_time_integration_rk4
   use ocn_time_integration_split

   implicit none
   private
   save

   public :: ocn_timestep, &
             ocn_timestep_init

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

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

    logical :: rk4On, splitOn

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_timestep
!
!> \brief MPAS ocean time integration driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine handles a single timestep for the ocean. It determines
!>  the time integrator that will be used for the run, and calls the
!>  appropriate one.
!
!-----------------------------------------------------------------------

   subroutine ocn_timestep(domain, dt, timeStamp)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:))
   !                 plus mesh meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (dm_info) :: dminfo
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: diagnosticsPool, statePool, meshPool

      character (len=StrKIND), pointer :: xtime
      real (kind=RKIND), pointer :: daysSinceStartOfSim
      character (len=StrKIND), pointer :: simulationStartTime
      type (MPAS_Time_type) :: xtime_timeType, simulationStartTime_timeType


      if (rk4On) then
         call ocn_time_integrator_rk4(domain, dt)
      elseif (splitOn) then
         call ocn_time_integrator_split(domain, dt)
     endif

     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

        call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)

        xtime = timeStamp

        ! compute time since start of simulation, in days
        call mpas_pool_get_array(diagnosticsPool, 'simulationStartTime', simulationStartTime)
        call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim',daysSinceStartOfSim)
        call mpas_set_time(xtime_timeType, dateTimeString=xtime)
        call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
        call mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType,dt=daysSinceStartOfSim)

        !$omp single
        daysSinceStartOfSim = daysSinceStartOfSim*days_per_second
        !$omp end single

        block => block % next
     end do

   end subroutine ocn_timestep!}}}

   subroutine ocn_timestep_init(err)!{{{

      integer, intent(out) :: err

      character (len=StrKIND), pointer :: config_time_integrator

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_time_integrator', config_time_integrator)

      rk4On = .false.
      splitOn = .false.

      if (trim(config_time_integrator) == 'RK4') then
          rk4On = .true.
      elseif (trim(config_time_integrator) == 'split_explicit' &
          .or.trim(config_time_integrator) == 'unsplit_explicit') then
          splitOn = .true.
      else
          err = 1
          write (stderrUnit,*) 'Incorrect choice for config_time_integrator:', trim(config_time_integrator)
          write (stderrUnit,*) '   choices are: RK4, split_explicit, unsplit_explicit'
      endif


   end subroutine ocn_timestep_init!}}}

end module ocn_time_integration

! vim: foldmethod=marker
