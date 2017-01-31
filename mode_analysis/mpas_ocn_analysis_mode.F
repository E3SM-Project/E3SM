! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_analysis_mode
!
!> \brief Main driver for MPAS ocean analysis core in post-processing mode
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This module contains the drivers for the MPAS ocean analysis core in
!>  post-processing mode.  None of these routines are used in run-time mode
!>  by the ocean core.
!
!-----------------------------------------------------------------------

module ocn_analysis_mode

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_stream_manager
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_io_units
   use mpas_decomp
   use mpas_abort, only : mpas_dmpar_global_abort

   use ocn_analysis_driver
   use ocn_init_routines
   use ocn_diagnostics
   use ocn_equation_of_state
   use ocn_constants

   private

   public :: ocn_analysis_mode_init, ocn_analysis_mode_run, ocn_analysis_mode_finalize
   public :: ocn_analysis_mode_setup_clock

   contains

!***********************************************************************
!
!  function ocn_analysis_mode_init
!
!> \brief   Initialize MPAS-Ocean analysis post-processing core
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This function calls all initializations required to begin a
!>  simulation with MPAS-Ocean in analysis mode.
!
!-----------------------------------------------------------------------

   function ocn_analysis_mode_init(domain, startTimeStamp) result(ierr)!{{{

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
      integer :: ierr

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: diagnosticsPool

      integer :: err_tmp

      type (MPAS_timeInterval_type) :: timeStep

      ! remove dt later
      real (kind=RKIND) :: dt
      character (len=StrKIND), pointer :: xtime, simulationStartTime
      real (kind=RKIND), pointer :: daysSinceStartOfSim
      type (MPAS_Time_type) :: xtime_timeType, simulationStartTime_timeType
      type (MPAS_Time_Type) :: startTime

      ierr = 0

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, err_tmp)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp)
      ierr = ior(ierr, err_tmp)

      ! Setup ocean config pool
      call ocn_constants_init(domain % configs, domain % packages)

      !
      ! Read input data for model
      !
      call mpas_timer_start('io_read')
      call MPAS_stream_mgr_read(domain % streamManager, streamID='mesh', whence=MPAS_STREAM_NEAREST, ierr=err_tmp)
      call MPAS_stream_mgr_read(domain % streamManager, streamID='input', ierr=err_tmp)
      call mpas_timer_stop('io_read')
      call mpas_timer_start('io_reset_alarms')
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='input', ierr=err_tmp)
      call mpas_timer_stop('io_reset_alarms')

      ! Initialize submodules before initializing blocks.
      call ocn_equation_of_state_init(err_tmp)
      ierr = ior(ierr, err_tmp)

      call ocn_analysis_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)

      call ocn_init_routines_vert_coord(domain)

      call ocn_init_routines_compute_max_level(domain)

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)

         call ocn_init_routines_block(block, dt, ierr)
         if(ierr.eq.1) then
             call mpas_dmpar_global_abort('MPAS-ocean: ERROR: An error was encountered in ocn_init_routines_block')
         endif

         xtime = startTimeStamp

         ! Set simulationStartTime only if that variable is not read from the restart file.
         call mpas_pool_get_array(diagnosticsPool, 'simulationStartTime', simulationStartTime)
         if (trim(simulationStartTime)=="no_date_available") then
            simulationStartTime = startTimeStamp
         end if

         ! compute time since start of simulation, in days
         call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim',daysSinceStartOfSim)
         call mpas_set_time(xtime_timeType, dateTimeString=xtime)
         call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
         call mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType,dt=daysSinceStartOfSim)
         daysSinceStartOfSim = daysSinceStartOfSim*days_per_second

         block => block % next
      end do

   end function ocn_analysis_mode_init!}}}

!***********************************************************************
!
!  routine ocn_analysis_mode_setup_clock
!
!> \brief   Initialize timer variables
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine initializes all timer variables for analysis mode
!
!-----------------------------------------------------------------------

   function ocn_analysis_mode_setup_clock(core_clock, configs) result(ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character(len=StrKIND) :: restartTimeStamp
      integer :: err_tmp
      character (len=StrKIND), pointer :: config_start_time, config_stop_time
      character (len=StrKIND), pointer :: config_run_duration
      character (len=StrKIND), pointer :: config_dt, config_restart_timestamp_name

      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      if ( trim(config_start_time) == "file" ) then
         open(22,file=config_restart_timestamp_name,form='formatted',status='old')
         read(22,*) restartTimeStamp
         close(22)
         call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp, ierr=err_tmp)
      else
         call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=err_tmp)
      end if

      call mpas_set_timeInterval(timeStep, timeString="0000_00:00:00", ierr=err_tmp)
      call mpas_set_timeInterval(runDuration, timeString="0000_00:00:00", ierr=err_tmp)
      call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=err_tmp)

   end function ocn_analysis_mode_setup_clock!}}}

!***********************************************************************
!
!  function ocn_analysis_mode_run
!
!> \brief   Main driver for MPAS-Ocean analysis mode
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This function includes the calls to perform analysis of the input file.
!
!-----------------------------------------------------------------------

   function ocn_analysis_mode_run(domain) result(ierr)!{{{

      type (domain_type), intent(inout) :: domain

      integer :: itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: err, ierr

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: scratchPool

      type (MPAS_timeInterval_type) :: timeStep
      character (len=StrKIND), pointer :: config_dt
      logical, pointer :: config_write_output_on_startup

      ierr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_dt', config_dt)
      call mpas_pool_get_config(ocnConfigs, 'config_write_output_on_startup', config_write_output_on_startup)

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(stderrUnit,*) 'Initial time ', trim(timeStamp)

      ! fill in diagnostics variables
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 1)
         block_ptr => block_ptr % next
      end do

      if (config_write_output_on_startup) then
         call ocn_analysis_compute_startup(domain, err)
         call mpas_timer_start('io_write')
         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_timer_stop('io_write')
      endif

   end function ocn_analysis_mode_run!}}}

!***********************************************************************
!
!  function ocn_analysis_mode_finalize
!
!> \brief   Finalize MPAS-Ocean analysis mode
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This function finalizes the MPAS-Ocean core that was run with analysis mode.
!
!-----------------------------------------------------------------------

   function ocn_analysis_mode_finalize(domain) result(iErr)!{{{

      type (domain_type), intent(inout) :: domain
      integer :: ierr

      iErr = 0

      call mpas_destroy_clock(domain % clock, ierr)

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

   end function ocn_analysis_mode_finalize!}}}

end module ocn_analysis_mode

! vim: foldmethod=marker
