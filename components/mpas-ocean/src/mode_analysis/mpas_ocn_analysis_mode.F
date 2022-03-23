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
!> \brief Driver for MPAS ocean analysis core in post-processing mode
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This module contains the drivers for the MPAS ocean analysis core
!>  in post-processing mode.  None of these routines are used in
!>  run-time mode by the ocean core.
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
   use mpas_log
   use mpas_decomp

   use ocn_analysis_driver
   use ocn_init_routines
   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_equation_of_state
   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_surface_land_ice_fluxes

   implicit none
   private

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------

   public :: ocn_analysis_mode_init, &
             ocn_analysis_mode_run, &
             ocn_analysis_mode_finalize, &
             ocn_analysis_mode_setup_clock

!***********************************************************************

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

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain    !< [inout] All ocean model state and variables

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      character(len=*), intent(out) :: &
         startTimeStamp   !< [out] start time in char string

      integer :: ierr     !< [out] returned error code

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: &
         block           ! structure containing data for subdomain block

      integer :: err_tmp ! internal error flag for routine calls

      real (kind=RKIND) :: &
         dt              ! model time step in seconds

      type (MPAS_timeInterval_type) :: &
         timeStep        ! model time step in time interval form

      type (MPAS_Time_type) :: &
         xtime_timeType, &! temporary for computing time since start
         startTime,      &! start time of simulation clock
         simulationStartTime_timeType ! start time again

      ! End preamble
      !-------------
      ! Begin code

      ierr = 0 ! initialize return error to success

      ! Set startTimeStamp based on start time of the simulation clock
      startTime = mpas_get_clock_time(domain%clock, MPAS_START_TIME, &
                                      err_tmp)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp)
      ierr = ior(ierr, err_tmp)
      if (err_tmp /= 0) call mpas_log_write( &
                        'Error setting start time ', MPAS_LOG_ERR)

      ! Set up ocean constants
      call ocn_constants_init(domain%configs, domain%packages)

      ! Setup ocean config pool
      call ocn_config_init(domain%configs)

      ! Read input data for model
      call mpas_timer_start('io_read')
      call MPAS_stream_mgr_read(domain%streamManager, streamID='mesh', &
                             whence=MPAS_STREAM_NEAREST, ierr=err_tmp)
      call mpas_timer_stop('io_read')
      if (err_tmp /= 0) call mpas_log_write( &
                        'Error reading mesh file ', MPAS_LOG_CRIT)

      call mpas_timer_start('io_read')
      call MPAS_stream_mgr_read(domain%streamManager, streamID='input',&
                                ierr=err_tmp)
      call mpas_timer_stop('io_read')
      if (err_tmp /= 0) call mpas_log_write( &
                        'Error reading input file ', MPAS_LOG_CRIT)

      call mpas_timer_start('io_reset_alarms')
      call MPAS_stream_mgr_reset_alarms(domain%streamManager, &
                                        streamID='input', ierr=err_tmp)
      call mpas_timer_stop('io_reset_alarms')
      if (err_tmp /= 0) call mpas_log_write( &
                        'Error resetting io alarms ', MPAS_LOG_ERR)

      ! Initialize all mesh data before other inits
      call ocn_meshCreate(domain)

      ! Initialize submodules before initializing blocks.
      call ocn_equation_of_state_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)

      call ocn_analysis_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      call ocn_init_routines_block(block, dt, err_tmp)
      ierr = ior(ierr, err_tmp)
      if (err_tmp /= 0) then
         call mpas_log_write( &
            'An error was encountered in ocn_init_routines_block', &
            MPAS_LOG_CRIT)
      endif

      xtime = startTimeStamp
      ! Set simulationStartTime only if that variable is not read from
      ! the restart file.
      if (trim(simulationStartTime)=="no_date_available") then
         simulationStartTime = startTimeStamp
      end if

      !$omp master
      ! compute time since start of simulation, in days
      call mpas_set_time(xtime_timeType, dateTimeString=xtime)
      call mpas_set_time(simulationStartTime_timeType, &
                         dateTimeString=simulationStartTime)
      call mpas_get_timeInterval( &
           xtime_timeType - simulationStartTime_timeType, &
           dt=daysSinceStartOfSim)
      daysSinceStartOfSim = daysSinceStartOfSim*days_per_second

      !$omp end master
      !$omp barrier

   !--------------------------------------------------------------------

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
      type (mpas_pool_type), pointer :: scratchPool

      type (MPAS_timeInterval_type) :: timeStep
      character (len=StrKIND), pointer :: config_dt
      logical, pointer :: config_write_output_on_startup

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft

      real (kind=RKIND), dimension(:), pointer :: &
        layerThickness, normalVelocity

      real (kind=RKIND), dimension(:), pointer :: ssh

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      ierr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_dt', config_dt)
      call mpas_pool_get_config(ocnConfigs, 'config_write_output_on_startup', config_write_output_on_startup)

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      call mpas_log_write('Initial time ' // trim(timeStamp))

      ! fill in diagnostics variables
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

         call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
         call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
         call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

         if (landIcePressureOn) then 
            call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
            call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
         endif

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
         call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThickness, normalVelocity)
         !$acc update device (normalTransportVelocity, &
         !$acc                normalGMBolusVelocity)
         !$acc enter data copyin(atmosphericPressure, seaIcePressure)
         !$acc enter data copyin(ssh)
         !$acc enter data copyin(activeTracers)
         !$acc update device(tracersSurfaceValue)
         if ( associated(frazilSurfacePressure) ) then
            !$acc enter data copyin(frazilSurfacePressure)
         endif
         if (landIcePressureOn) then
            !$acc enter data copyin(landIcePressure)
            !$acc enter data copyin(landIceDraft)
         endif
#endif
         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, scratchPool, tracersPool, 1)
#ifdef MPAS_OPENACC
         !$acc update host(layerThickEdge)
         !$acc update host(relativeVorticity, circulation)
         !$acc update host(vertTransportVelocityTop, &
         !$acc             vertGMBolusVelocityTop, &
         !$acc             relativeVorticityCell, &
         !$acc             divergence, &
         !$acc             kineticEnergyCell, &
         !$acc             tangentialVelocity, &
         !$acc             vertVelocityTop)
         !$acc update host(normRelVortEdge, normPlanetVortEdge, &
         !$acc             normalizedRelativeVorticityCell)
         !$acc update host (surfacePressure)
         !$acc update host(zMid, zTop)
         !$acc exit data copyout(ssh)
         !$acc update host(tracersSurfaceValue)
         !$acc update host(normalVelocitySurfaceLayer)
         !$acc exit data delete (atmosphericPressure, seaIcePressure)
         if ( associated(frazilSurfacePressure) ) then
            !$acc exit data delete(frazilSurfacePressure)
         endif
         if (landIcePressureOn) then
            !$acc exit data delete(landIcePressure)
            !$acc exit data delete(landIceDraft)
         endif
         !$acc exit data delete (activeTracers)
         !$acc exit data delete (layerThickness, normalVelocity)
         !$acc update host(density, potentialDensity, displacedDensity)
         !$acc update host(thermExpCoeff,  &
         !$acc&            salineContractCoeff)
         !$acc update host(montgomeryPotential, pressure)
         !$acc update host(RiTopOfCell, &
         !$acc             BruntVaisalaFreqTop)
         !$acc update host(tracersSurfaceLayerValue, &
         !$acc             indexSurfaceLayerDepth, &
         !$acc             normalVelocitySurfaceLayer, &
         !$acc             sfcFlxAttCoeff, &
         !$acc             surfaceFluxAttenuationCoefficientRunoff)
#endif
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
