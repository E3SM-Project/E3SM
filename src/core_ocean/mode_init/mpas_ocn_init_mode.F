! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_mode
!
!> \brief Main driver for MPAS ocean core
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains initialization and timestep drivers for
!>  the MPAS ocean core.
!
!-----------------------------------------------------------------------

module ocn_init_mode

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_stream_manager
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_io_units
   use mpas_constants
   use mpas_decomp

   use ocn_init_routines

   use ocn_equation_of_state

   use ocn_constants

   use ocn_init_spherical_utils

   !use ocn_init_TEMPLATE
   use ocn_init_baroclinic_channel
   use ocn_init_lock_exchange
   use ocn_init_internal_waves
   use ocn_init_overflow
   use ocn_init_global_ocean
   use ocn_init_cvmix_WSwSBF
   use ocn_init_iso
   use ocn_init_soma
   use ocn_init_ziso
   use ocn_init_sub_ice_shelf_2D
   use ocn_init_periodic_planar
   use ocn_init_ecosys_column
   use ocn_init_sea_mount
   use ocn_init_isomip
   use ocn_init_isomip_plus

   implicit none
   private

   public :: ocn_init_mode_init, ocn_init_mode_run, ocn_init_mode_finalize
   public :: ocn_init_mode_setup_clock, ocn_init_mode_validate_configuration

   contains

!***********************************************************************
!
!  function ocn_init_mode_init
!
!> \brief   Initialize MPAS-Ocean core in init mode
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!>  This function calls all initializations required to start MPAS-Ocean in
!>  init mode.
!
!-----------------------------------------------------------------------

   function ocn_init_mode_init(domain, startTimeStamp) result(ierr)!{{{

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
      integer :: ierr

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block

      integer :: err_tmp
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: maxDensity, maxDensity_global
      real (kind=RKIND), dimension(:), pointer :: meshDensity
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime
      type (MPAS_TimeInterval_type) :: timeStep

      logical, pointer :: config_do_restart, config_filter_btr_mode, config_conduct_tests
      logical, pointer :: config_write_stats_on_startup
      character (len=StrKIND), pointer :: config_vert_coord_movement, config_pressure_gradient_type
      real (kind=RKIND), pointer :: config_maxMeshDensity

      ierr = 0

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, err_tmp)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp)
      ierr = ior(ierr, err_tmp)

      ! Setup ocean config pool
      call ocn_constants_init(domain % configs, domain % packages)

      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort("MPAS-ocean: ERROR: Failed validation...")
      end if

      !
      ! Read input data for model
      !
      call mpas_timer_start('io_read', .false.)
      call MPAS_stream_mgr_read(domain % streamManager, streamID='input_init', ierr=err_tmp)
      call mpas_timer_stop('io_read')

      call mpas_timer_start('reset_io_alarms', .false.)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='input_init', ierr=err_tmp)
   !  call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
      call mpas_timer_stop('reset_io_alarms')

      ! Read the remaining input streams
      call mpas_timer_start('io_read', .false.)
      call mpas_stream_mgr_read(domain % streamManager, ierr=err_tmp)
      ierr = ior(ierr, err_tmp)
      call mpas_timer_stop('io_read')
      call mpas_timer_start('reset_io_alarms', .false.)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=err_tmp)
      ierr = ior(ierr, err_tmp)
      call mpas_timer_stop('reset_io_alarms')

      ! Initialize submodules before initializing blocks.
      call ocn_equation_of_state_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      if(ierr.eq.1) then
          call mpas_dmpar_global_abort('MPAS-ocean: ERROR: An error was encountered while initializing the MPAS-Ocean init mode')
      endif

      !
      ! Initialize core
      !
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         xtime = startTimeStamp
         block => block % next
      end do

      ! Expand sphere if it needs to be expanded
      call ocn_init_expand_sphere(domain, domain % streamManager, a, ierr)

   end function ocn_init_mode_init!}}}

!***********************************************************************
!
!  function ocn_init_mode_setup_clock
!
!> \brief   Setup MPAS-Ocean clock
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!>  This function initializes the MPAS-Ocean clock for the init mode.
!
!-----------------------------------------------------------------------
   function ocn_init_mode_setup_clock(core_clock, configs) result(ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character(len=StrKIND) :: restartTimeStamp
      character(len=StrKIND), pointer :: config_start_time, config_stop_time, config_run_duration
      character(len=StrKIND), pointer :: config_dt, config_restart_timestamp_name
      integer :: err_tmp

      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      call mpas_set_time(startTime, dateTimeString=config_start_time, ierr=err_tmp)
      call mpas_set_timeInterval(timeStep, timeString=config_dt, ierr=err_tmp)
      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=err_tmp)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
            if(startTime + runduration /= stopTime) then
               write(stderrUnit,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=err_tmp)
      else
         write(stderrUnit, *) ' Warning: config_run_duration and config_start_time were "none", setting run duration to 1 second.'
         call mpas_set_timeInterval(runDuration, timeString="0000_00:00:01", ierr=err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=err_tmp)
      end if

   end function ocn_init_mode_setup_clock!}}}

!***********************************************************************
!
!  function ocn_init_mode_run
!
!> \brief   MPAS-Ocean init mode run step
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!>  This function sets up the initial configuration using the MPAS-Ocean init
!>  mode.
!
!-----------------------------------------------------------------------

   function ocn_init_mode_run(domain) result(iErr)!{{{

      type (domain_type), intent(inout) :: domain
      integer :: iErr

      integer :: itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp

      type (mpas_pool_type), pointer :: averagePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: forcingPool

      type (MPAS_timeInterval_type) :: timeStep

      character (len=StrKIND), pointer :: config_init_configuration

      ierr = 0

      ! Eventually, dt should be domain specific
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)
      write(stderrUnit, *) ' Generating configuration: ' // trim(config_init_configuration)

      call ocn_init_setup_baroclinic_channel(domain, ierr)
      call ocn_init_setup_lock_exchange(domain, ierr)
      call ocn_init_setup_internal_waves(domain, ierr)
      call ocn_init_setup_overflow(domain, ierr)
      call ocn_init_setup_global_ocean(domain, ierr)
      call ocn_init_setup_cvmix_WSwSBF(domain, ierr)
      call ocn_init_setup_iso(domain, ierr)
      call ocn_init_setup_soma(domain, ierr)
      call ocn_init_setup_ziso(domain, ierr)
      call ocn_init_setup_sub_ice_shelf_2D(domain, ierr)
      call ocn_init_setup_periodic_planar(domain, ierr)
      call ocn_init_setup_ecosys_column(domain, ierr)
      call ocn_init_setup_sea_mount(domain, ierr)
      call ocn_init_setup_isomip(domain, ierr)
      call ocn_init_setup_isomip_plus(domain, ierr)
      !call ocn_init_setup_TEMPLATE(domain, ierr)

      write(stderrUnit, *) ' Completed setup of: ' // trim(config_init_configuration)
      call mpas_timer_start('io_write', .false.)
      call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
      call mpas_timer_stop('io_write')
      call mpas_timer_start('reset_io_alarms', .false.)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
      call mpas_timer_stop('reset_io_alarms')

   end function ocn_init_mode_run!}}}

!***********************************************************************
!
!  function ocn_init_mode_finalize
!
!> \brief   MPAS-Ocean init mode run step
!> \author  Doug Jacobsen
!> \date    06/15/2015
!> \details
!>  This function sets up the initial configuration using the MPAS-Ocean init
!>  mode.
!
!-----------------------------------------------------------------------

   function ocn_init_mode_finalize(domain) result(iErr)!{{{

      type (domain_type), intent(inout) :: domain
      integer :: ierr

      iErr = 0

      call mpas_destroy_clock(domain % clock, ierr)

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

   end function ocn_init_mode_finalize!}}}

!***********************************************************************
!
!  routine ocn_init_validate_configuration
!
!> \brief   Configuration validation routine
!> \author  Doug Jacobsen
!> \date    03/20/2015
!> \details
!>  This routine is used to validate the namelist options against the
!>  configuration definition. It will call the validate routines for each of the
!>  configurations to ensure namelist options are set in a valid way.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_mode_validate_configuration(configPool, packagePool, iocontext, iErr)!{{{

      type (mpas_pool_type), intent(inout) :: configPool !< Input: Pool with namelist options
      type (mpas_pool_type), intent(inout) :: packagePool !< Input: Pool with packages
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer, intent(out) :: iErr !< Output: Error core

      logical, pointer :: cullCellsActive

      logical, pointer :: config_write_cull_cell_mask

      integer :: err_tmp

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_write_cull_cell_mask', config_write_cull_cell_mask)
      call mpas_pool_get_package(packagePool, 'cullCellsActive', cullCellsActive)

      if ( config_write_cull_cell_mask ) then
         cullCellsActive = .true.
      end if

      call ocn_init_validate_baroclinic_channel(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_lock_exchange(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_internal_waves(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_overflow(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_global_ocean(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_cvmix_WSwSBF(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_iso(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_soma(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_ziso(configPool, packagePool, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_sub_ice_shelf_2D(configPool, packagePool,iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_periodic_planar(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_ecosys_column(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_sea_mount(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_isomip(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      call ocn_init_validate_isomip_plus(configPool, packagePool, iocontext, iErr=err_tmp)
      iErr = ior(iErr, err_tmp)
      ! call ocn_init_validate_TEMPLATE(configPool, packagePool, iocontext, iErr=err_tmp)
      ! iErr = ior(iErr, err_tmp)
   end subroutine ocn_init_mode_validate_configuration!}}}

end module ocn_init_mode

! vim: foldmethod=marker
