! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_analysis_driver
!
!> \brief Driver for MPAS ocean analysis core
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This is the driver for the MPAS ocean core.
!
!-----------------------------------------------------------------------

module ocn_analysis_driver

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use mpas_stream_manager

   use ocn_constants
   use ocn_global_stats
   use ocn_surface_area_weighted_averages
   use ocn_layer_volume_weighted_averages
   use ocn_zonal_mean
   use ocn_okubo_weiss
   use ocn_water_mass_census
   use ocn_meridional_heat_transport
   use ocn_test_compute_interval
   use ocn_high_frequency_output
   use ocn_time_filters
   use ocn_lagrangian_particle_tracking
   use ocn_eliassen_palm
   use ocn_mixed_layer_depths
   use ocn_time_series_stats
   use ocn_pointwise_stats
   use ocn_debug_diagnostics
   use ocn_regional_stats
   use ocn_rpn_calculator
   use ocn_transect_transport
   use ocn_eddy_product_variables
!   use ocn_TEM_PLATE

   implicit none
   private
   save

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

   public :: ocn_analysis_setup_packages, &
             ocn_analysis_bootstrap, &
             ocn_analysis_init, &
             ocn_analysis_compute_startup, &
             ocn_analysis_compute, &
             ocn_analysis_write, &
             ocn_analysis_restart, &
             ocn_analysis_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


   character (len=*), parameter :: initReadTimerPrefix = 'init_read_'
   character (len=*), parameter :: initTimerPrefix = 'init_'
   character (len=*), parameter :: computeTimerPrefix = 'compute_'
   character (len=*), parameter :: computeStartupTimerPrefix = 'compute_startup_'
   character (len=*), parameter :: writeTimerPrefix = 'write_'
   character (len=*), parameter :: alarmTimerPrefix = 'reset_alarm_'
   character (len=*), parameter :: restartTimerPrefix = 'restart_'
   character (len=*), parameter :: finalizeTimerPrefix = 'finalize_'
   character (len=*), parameter :: computeAlarmSuffix = 'CMPALRM'

   character (len=StrKIND), parameter :: timeSeriesDailyTAG = 'Daily'
   character (len=StrKIND), parameter :: timeSeriesMonthlyTAG = 'Monthly'
   character (len=StrKIND), parameter :: timeSeriesClimatologyTAG = 'Climatology'
   character (len=StrKIND), parameter :: timeSeriesCustomTAG = 'Custom'
   character (len=StrKIND), parameter :: regionalStatsDailyTAG = 'Daily'
   character (len=StrKIND), parameter :: regionalStatsMonthlyTAG = 'Monthly'
   character (len=StrKIND), parameter :: regionalStatsWeeklyTAG = 'Weekly'
   character (len=StrKIND), parameter :: regionalStatsCustomTAG = 'Custom'

   type (mpas_pool_type), pointer :: analysisMemberList

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_analysis_setup_packages
!
!> \brief   Setup packages for MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine is intended to configure the packages for all
!>   ocean analysis members.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_setup_packages(configPool, packagePool, iocontext, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      integer :: err_tmp

      character (len=StrKIND) :: configName, packageName
      logical, pointer :: config_AM_enable
      logical, pointer :: AMPackageActive
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      err = 0

      call mpas_pool_create_pool(analysisMemberList)
      call mpas_pool_add_config(analysisMemberList, 'globalStats', 1)
      call mpas_pool_add_config(analysisMemberList, 'testComputeInterval', 1)
      call mpas_pool_add_config(analysisMemberList, 'layerVolumeWeightedAverage', 1)
      call mpas_pool_add_config(analysisMemberList, 'meridionalHeatTransport', 1)
      call mpas_pool_add_config(analysisMemberList, 'okuboWeiss', 1)
      call mpas_pool_add_config(analysisMemberList, 'surfaceAreaWeightedAverages', 1)
      call mpas_pool_add_config(analysisMemberList, 'waterMassCensus', 1)
      call mpas_pool_add_config(analysisMemberList, 'zonalMean', 1)
      call mpas_pool_add_config(analysisMemberList, 'highFrequencyOutput', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeFilters', 1)
      call mpas_pool_add_config(analysisMemberList, 'lagrPartTrack', 1)
      call mpas_pool_add_config(analysisMemberList, 'eliassenPalm', 1)
      call mpas_pool_add_config(analysisMemberList, 'mixedLayerDepths', 1)
      call mpas_pool_add_config(analysisMemberList, 'rpnCalculator', 1)
      call mpas_pool_add_config(analysisMemberList, 'regionalStats', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeSeriesStatsDaily', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeSeriesStatsMonthly', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeSeriesStatsClimatology', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeSeriesStatsCustom', 1)
      call mpas_pool_add_config(analysisMemberList, 'regionalStatsDaily', 1)
      call mpas_pool_add_config(analysisMemberList, 'regionalStatsWeekly', 1)
      call mpas_pool_add_config(analysisMemberList, 'regionalStatsMonthly', 1)
      call mpas_pool_add_config(analysisMemberList, 'regionalStatsCustom', 1)
      call mpas_pool_add_config(analysisMemberList, 'pointwiseStats', 1)
      call mpas_pool_add_config(analysisMemberList, 'debugDiagnostics', 1)
      call mpas_pool_add_config(analysisMemberList, 'transectTransport', 1)
      call mpas_pool_add_config(analysisMemberList, 'eddyProductVariables', 1)
!     call mpas_pool_add_config(analysisMemberList, 'temPlate', 1)

      ! DON'T EDIT BELOW HERE

      ! Iterate over all analysis members to setup packages
      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(configPool, configName, config_AM_enable)

         if ( config_AM_enable ) then
            packageName = poolItr % memberName(1:nameLength) // 'AMPKGActive'
            call mpas_pool_get_package(packagePool, packageName, AMPackageActive)
            AMPackageActive = .true.
         end if
      end do

   end subroutine ocn_analysis_setup_packages!}}}

!***********************************************************************
!
!  routine ocn_analysis_bootstrap
!
!> \brief   Bootstrap analysis members (pre-init configuration)
!> \author  Doug Jacobsen
!> \date    10/08/2015
!> \details
!>  This routine will read either a restart or an input stream for each analysis member.
!>  The stream names that will be read are controlled via the analysis member's
!>     - config_AM_${AM}_restart_stream
!>     - config_AM_${AM}_input_stream
!>  namelist options.
!>
!>  If the AM doesn't specify either of these, it will be ignored. If the AM
!>  specifies only the restart stream, it will only be read if the config_do_restart flag
!>  for the model is set to true. If the AM specifies both, the restart_stream will be read if
!>  config_do_restart is true, and the input_stream will be read if config_do_restart is false.
!>
!>  After this call, alarms on both streams are reset.
!>
!>  Additionally, if a bootstrap subroutine has been defined properly for the
!>  analysis member, it will be called here.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_bootstrap(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      integer :: err_tmp

      character (len=StrKIND) :: configName, alarmName, restartStreamName, inputStreamName, timerName
      logical, pointer :: config_AM_enable, config_do_restart
      character (len=StrKIND), pointer :: config_AM_restart_stream, config_AM_input_stream
      integer :: nameLength
      type (mpas_pool_iterator_type) :: poolItr

      logical :: streamFound
      character  (len=StrKIND) :: referenceTimeString, outputIntervalString
      type (MPAS_Time_Type) :: referenceTime
      type (MPAS_TimeInterval_type) :: alarmTimeStep

      integer :: poolErrorLevel

      err = 0

      poolErrorLevel = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)

      call mpas_timer_start('analysis_bootstrap')

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'

         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            timerName = trim(initReadTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)

            write(stderrUnit,*) '      Bootstrapping AM ' // poolItr % memberName(1:nameLength)
            call ocn_bootstrap_analysis_members(domain, poolItr % memberName(1:nameLength), ierr=err)

            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_restart_stream'
            nullify(config_AM_restart_stream)
            call mpas_pool_get_config(domain % configs, configName, config_AM_restart_stream)

            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_input_stream'
            nullify(config_AM_input_stream)
            call mpas_pool_get_config(domain % configs, configName, config_AM_input_stream)

            ! Verify the restart stream exists
            if ( associated(config_AM_restart_stream) ) then
               if ( trim(config_AM_restart_stream) == 'none' ) then
                  ! If the stream is set to 'none' nullify the config, so it doesn't get read in
                  nullify(config_AM_restart_stream)
               else if ( .not. mpas_stream_mgr_stream_exists(domain % streamManager, config_AM_restart_stream) ) then
                  call mpas_dmpar_global_abort('MPAS-ocean: ERROR: Stream named ''' // trim(config_AM_restart_stream) // &
                                               ''' does not exist in config for analysis member ''' // &
                                               trim(poolItr % memberName(1:nameLength)) // '''')
               end if
            end if

            ! Verify the input stream exists
            if ( associated(config_AM_input_stream) ) then
               if ( trim(config_AM_input_stream) == 'none' ) then
                  ! If the stream is set to 'none' nullify the config, so it doesn't get read in
                  nullify(config_AM_input_stream)
               else if ( .not. mpas_stream_mgr_stream_exists(domain % streamManager, config_AM_input_stream) ) then
                  call mpas_dmpar_global_abort('MPAS-ocean: ERROR: Stream named ''' // trim(config_AM_input_stream) // &
                                               ''' does not exist in config for analysis member ''' // &
                                               trim(poolItr % memberName(1:nameLength)) // '''')
               end if
            end if

            ! Handle reading of streams that exist.
            if ( associated(config_AM_restart_stream) .and. associated(config_AM_input_stream) ) then

               if ( config_do_restart ) then
                  call mpas_stream_mgr_read(domain % streamManager, streamID=config_AM_restart_stream, ierr=err)
               else
                  call mpas_stream_mgr_read(domain % streamManager, streamID=config_AM_input_stream, ierr=err)
               end if
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=config_AM_restart_stream, &
                                                 direction=MPAS_STREAM_INPUT, ierr=err)
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=config_AM_input_stream, &
                                                 direction=MPAS_STREAM_INPUT, ierr=err)
            else if ( associated(config_AM_restart_stream) ) then
               if ( config_do_restart ) then
                  call mpas_stream_mgr_read(domain % streamManager, streamID=config_AM_restart_stream, ierr=err)
               end if
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=config_AM_restart_stream, &
                                                 direction=MPAS_STREAM_INPUT, ierr=err)
            else if ( associated(config_AM_input_stream) ) then
               if ( .not. config_do_restart ) then
                  call mpas_stream_mgr_read(domain % streamManager, streamID=config_AM_input_stream, ierr=err)
               end if
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=config_AM_input_stream, &
                                                 direction=MPAS_STREAM_INPUT, ierr=err)
            end if
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_bootstrap')

      call mpas_pool_set_error_level(poolErrorLevel)

   end subroutine ocn_analysis_bootstrap!}}}

!***********************************************************************
!
!  routine ocn_analysis_init
!
!> \brief   Initialize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all initializations required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_init(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      integer :: err_tmp

      character (len=StrKIND) :: configName, alarmName, streamName, timerName
      logical, pointer :: config_AM_enable
      character (len=StrKIND), pointer :: config_AM_compute_interval, config_AM_output_stream, config_start_time
      integer :: nameLength
      type (mpas_pool_iterator_type) :: poolItr

      logical :: streamFound
      character  (len=StrKIND) :: referenceTimeString, outputIntervalString
      type (MPAS_Time_Type) :: referenceTime
      type (MPAS_TimeInterval_type) :: alarmTimeStep

      err = 0

      call mpas_timer_start('analysis_init')

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            write(stderrUnit,*) '      Initializing AM ' // poolItr % memberName(1:nameLength)
            timerName = trim(initTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call ocn_init_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)

            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_compute_interval'
            call mpas_pool_get_config(domain % configs, configName, config_AM_compute_interval)

            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_output_stream'
            call mpas_pool_get_config(domain % configs, configName, config_AM_output_stream)

            if ( config_AM_compute_interval == 'dt' ) then
               alarmTimeStep = mpas_get_clock_timestep(domain % clock, err_tmp)
               call mpas_get_timeInterval(alarmTimeStep, timeString=config_AM_compute_interval, ierr=err_tmp)
            end if

            ! Verify stream exists before trying to use output_interval
            if ( config_AM_output_stream /= 'none' ) then
               streamFound = .false.

               call mpas_stream_mgr_begin_iteration(domain % streamManager)
               do while ( mpas_stream_mgr_get_next_stream(domain % streamManager, streamName) )
                  if ( trim(streamName) == trim(config_AM_output_stream) ) then
                     streamFound = .true.
                  end if
               end do

               if ( .not. streamFound ) then
                  call mpas_dmpar_global_abort('MPAS-ocean: ERROR: Stream ' // trim(config_AM_output_stream) // &
                                               ' does not exist. Exiting...')
               end if
            end if

            if ( config_AM_compute_interval == 'output_interval' .and. config_AM_output_stream == 'none') then
               call mpas_dmpar_global_abort('MPAS-ocean: ERROR: Analysis member has compute_interval of ''output_interval'' ' // &
                                            'without an output stream.')
            end if

            if ( config_AM_compute_interval /= 'output_interval' ) then
               alarmName = poolItr % memberName(1:nameLength) // computeAlarmSuffix
               call mpas_set_timeInterval(alarmTimeStep, timeString=config_AM_compute_interval, ierr=err_tmp)
               if ( config_AM_output_stream /= 'none' ) then
                  call MPAS_stream_mgr_get_property(domain % streamManager, config_AM_output_stream, &
                                                    MPAS_STREAM_PROPERTY_REF_TIME, referenceTimeString, err_tmp)
                  call mpas_set_time(referenceTime, dateTimeString=referenceTimeString, ierr=err_tmp)
               else
                  call mpas_pool_get_config(domain % configs, 'config_start_time', config_start_time)

                  ! TODO FIXME I'm not sure what it's supposed to be
                  !            but 'file' is causing the code to fail
                  if (trim(config_start_time) == 'file') then
                    ! FIXME big kludge
                    call mpas_set_time(referenceTime, &
                      dateTimeString='0000-01-01_00:00:00', ierr=err_tmp)
                  else
                    ! FIXME this is what it was without the if-else
                    !       I suppose it's supposed to actually read it from
                    !       the file first
                    call mpas_set_time(referenceTime, dateTimeString=config_start_time, ierr=err_tmp)
                  end if

               end if
               call mpas_add_clock_alarm(domain % clock, alarmName, referenceTime, alarmTimeStep, ierr=err_tmp)
               call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_tmp)
            end if

            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_init')

   end subroutine ocn_analysis_init!}}}

!***********************************************************************
!
!  routine ocn_analysis_compute_startup
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_compute_startup(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      integer :: timeLevel, err_tmp

      character (len=StrKIND) :: configName, timerName
      character (len=StrKIND), pointer :: config_AM_output_stream
      logical, pointer :: config_AM_enable, config_AM_write_on_startup, config_AM_compute_on_startup
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      err = 0

      call mpas_timer_start('analysis_compute_startup')

      timeLevel=1

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_compute_on_startup'
            call mpas_pool_get_config(domain % configs, configName, config_AM_compute_on_startup)
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_write_on_startup'
            call mpas_pool_get_config(domain % configs, configName, config_AM_write_on_startup)

            if ( config_AM_compute_on_startup ) then
               timerName = trim(computeStartupTimerPrefix) // poolItr % memberName(1:nameLength)
               write(stderrUnit,*) '      Computing AM ' // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName)
               call ocn_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
               err = ior(err, err_tmp)
            end if

            if ( config_AM_write_on_startup ) then
               configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_output_stream'
               call mpas_pool_get_config(domain % configs, configName, config_AM_output_stream)
               if ( config_AM_output_stream /= 'none' ) then
                  write(stderrUnit,*) '      Writing AM ' // poolItr % memberName(1:nameLength)
                  call mpas_stream_mgr_write(domain % streamManager, streamID=config_AM_output_stream, &
                                             forceWriteNow=.true., ierr=err_tmp)
               end if
               if (.not. config_AM_compute_on_startup) then
                  write(stderrUnit, *) ' *** WARNING: write_on_startup called without compute_on_startup for analysis member: ' &
                    // poolItr % memberName(1:nameLength) // '.'
               end if
            end if
         end if
      end do

      call mpas_timer_stop('analysis_compute_startup')

   end subroutine ocn_analysis_compute_startup!}}}

!***********************************************************************
!
!  routine ocn_analysis_compute
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_compute(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      integer :: timeLevel, err_tmp

      character (len=StrKIND) :: configName, alarmName, timerName
      character (len=StrKIND), pointer :: config_AM_output_stream, config_AM_compute_interval
      logical, pointer :: config_AM_enable
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      err = 0

      call mpas_timer_start('analysis_compute')

      timeLevel=1

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_compute_interval'
            call mpas_pool_get_config(domain % configs, configName, config_AM_compute_interval)
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_output_stream'
            call mpas_pool_get_config(domain % configs, configName, config_AM_output_stream)

            ! Build name of alarm for analysis member
            alarmName = poolItr % memberName(1:nameLength) // computeAlarmSuffix
            timerName = trim(computeTimerPrefix) // poolItr % memberName(1:nameLength)

            ! Compute analysis member just before output
            if ( config_AM_compute_interval == 'output_interval' .and. config_AM_output_stream /= 'none') then
               if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID=config_AM_output_stream, &
                    direction=MPAS_STREAM_OUTPUT, ierr=err_tmp) ) then
                  write(stderrUnit,*) '      Computing AM ' // poolItr % memberName(1:nameLength)
                  call mpas_timer_start(timerName)
                  call ocn_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
                  call mpas_timer_stop(timerName)
               end if
            else if ( mpas_is_alarm_ringing(domain % clock, alarmName, ierr=err_tmp) ) then
               call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_tmp)
               write(stderrUnit,*) '      Computing AM ' // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName)
               call ocn_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
            end if
         end if
      end do

      call mpas_timer_stop('analysis_compute')

   end subroutine ocn_analysis_compute!}}}

!***********************************************************************
!
!  routine ocn_analysis_restart
!
!> \brief   Save restart for MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all subroutines required to prepare to save
!>  the restart state for the MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_restart(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      integer :: err_tmp

      character (len=StrKIND) :: configName, timerName
      type (mpas_pool_iterator_type) :: poolItr
      logical, pointer :: config_AM_enable
      integer :: nameLength

      err = 0

      call mpas_timer_start('analysis_restart')

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            write(stderrUnit,*) '      Preparing AM ' // poolItr % memberName(1:nameLength) // ' for restart write'
            timerName = trim(restartTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call ocn_restart_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_restart')

   end subroutine ocn_analysis_restart!}}}

!***********************************************************************
!
!  routine ocn_analysis_write
!
!> \brief   Driver for MPAS-Ocean analysis output
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all output writing subroutines required for the
!>  MPAS-Ocean analysis driver.
!>  At this time this is just a stub, and all analysis output is written
!>  to the output file specified by config_output_name.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_write(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

      integer :: err_tmp

      character (len=StrKIND) :: configName, timerName, outputTimeString
      character (len=StrKIND), pointer :: config_AM_output_stream
      character (len=StrKIND), pointer :: config_AM_backwardOffset, config_AM_forwardOffset
      logical, pointer :: config_AM_enable
      type (mpas_pool_iterator_type) :: poolItr
      type (mpas_time_type) :: outputTime, nowTime
      type (mpas_timeinterval_type) :: offsetInt
      integer :: nameLength

      integer :: poolErrorLevel

      err = 0

      call mpas_timer_start('analysis_write')
      nowTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr=err_tmp)

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then

            poolErrorLevel = mpas_pool_get_error_level()
            call mpas_pool_set_error_level(MPAS_POOL_SILENT)

            nullify(config_AM_backwardOffset)
            nullify(config_AM_forwardOffset)
            write(stderrUnit,*) '      Writing AM ' // poolItr % memberName(1:nameLength)
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_output_stream'
            call mpas_pool_get_config(domain % configs, configName, config_AM_output_stream)
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_backward_output_offset'
            call mpas_pool_get_config(domain % configs, configName, config_AM_backwardOffset)
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_forward_output_offset'
            call mpas_pool_get_config(domain % configs, configName, config_AM_forwardOffset)

            call mpas_pool_set_error_level(poolErrorLevel)
            if ( config_AM_output_stream /= 'none' ) then
               timerName = trim(writeTimerPrefix) // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName)

               if ( associated(config_AM_backwardOffset) ) then
                  if ( associated(config_AM_forwardOffset) ) then
                     write(stderrUnit, *) 'WARNING: Both backward and forward offsets are set for AM ' // poolItr % memberName(1:nameLength)
                     write(stderrUnit, *) '         will only use backward offset. Forward offset will be ignored.'
                  end if
                  call mpas_set_timeinterval(offsetInt, timeString=config_AM_backwardOffset, ierr=err_tmp)
                  outputTime = nowTime - offsetInt
                  call mpas_get_time(outputTime, dateTimeString=outputTimeString, ierr=err_tmp)
               else if ( associated(config_AM_forwardOffset) ) then
                  call mpas_set_timeinterval(offsetInt, timeString=config_AM_backwardOffset, ierr=err_tmp)
                  outputTime = nowTime + offsetInt
                  call mpas_stream_mgr_write(domain % streamManager, streamID=config_AM_output_stream, ierr=err_tmp)
               else
                  outputTime = nowTime
               end if

               call mpas_get_time(outputTime, dateTimeString=outputTimeString, ierr=err_tmp)
               call mpas_stream_mgr_write(domain % streamManager, streamID=config_AM_output_stream, writeTime=outputTimeString, ierr=err_tmp)
               call mpas_timer_stop(timerName)
               timerName = trim(alarmTimerPrefix) // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName)
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=config_AM_output_stream, ierr=err_tmp)
               call mpas_timer_stop(timerName)
            end if
         end if
      end do

      call mpas_timer_stop('analysis_write')

   end subroutine ocn_analysis_write!}}}

!***********************************************************************
!
!  routine ocn_analysis_finalize
!
!> \brief   Finalize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all finalize routines required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_finalize(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      integer :: err_tmp

      character (len=StrKIND) :: configName, timerName
      logical, pointer :: config_AM_enable
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      err = 0

      call mpas_timer_start('analysis_finalize')

      call mpas_pool_begin_iteration(analysisMemberList)

      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            write(stderrUnit,*) '      Finalizing AM ' // poolItr % memberName(1:nameLength)
            timerName = trim(finalizeTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call ocn_finalize_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_finalize')

   end subroutine ocn_analysis_finalize!}}}

!***********************************************************************
!
!  routine ocn_bootstrap_analysis_members
!
!> \brief Analysis member initialization driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct init routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_bootstrap_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0
      err_tmp = 0

      nameLength = len_trim(analysisMemberName)

      !if ( analysisMemberName(1:nameLength) == 'testComputeInterval' ) then
      !   call ocn_bootstrap_test_compute_interval(domain, err_tmp)
      if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call ocn_bootstrap_time_series_stats(domain, timeSeriesDailyTAG, &
           err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call ocn_bootstrap_time_series_stats(domain, timeSeriesMonthlyTAG, &
           err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call ocn_bootstrap_time_series_stats(domain, timeSeriesClimatologyTAG, &
           err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call ocn_bootstrap_time_series_stats(domain, timeSeriesCustomTAG, &
           err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call ocn_bootstrap_pointwise_stats(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'transectTransport' ) then
        call ocn_init_transect_transport(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'temPlate' ) then
!        call ocn_init_TEM_PLATE(domain, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine ocn_bootstrap_analysis_members!}}}

!***********************************************************************
!
!  routine ocn_init_analysis_members
!
!> \brief Analysis member initialization driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct init routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'globalStats' ) then
         call ocn_init_global_stats(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'testComputeInterval' ) then
         call ocn_init_test_compute_interval(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'layerVolumeWeightedAverage' ) then
         call ocn_init_layer_volume_weighted_averages(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'meridionalHeatTransport' ) then
         call ocn_init_meridional_heat_transport(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'okuboWeiss' ) then
         call ocn_init_okubo_weiss(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
         call ocn_init_surface_area_weighted_averages(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'waterMassCensus' ) then
         call ocn_init_water_mass_census(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'zonalMean' ) then
         call ocn_init_zonal_mean(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call ocn_init_high_frequency_output(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeFilters' ) then
         call ocn_init_time_filters(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'lagrPartTrack' ) then
         call ocn_init_lagrangian_particle_tracking(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'eliassenPalm' ) then
         call ocn_init_eliassen_palm(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'mixedLayerDepths' ) then
         call ocn_init_mixed_layer_depths(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call ocn_init_pointwise_stats(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'debugDiagnostics' ) then
         call ocn_init_debug_diagnostics(domain, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'transectTransport' ) then
        call ocn_init_transect_transport(domain, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'eddyProductVariables' ) then
        call ocn_init_eddy_product_variables(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'temPlate' ) then
!        call ocn_init_TEM_PLATE(domain, err_tmp)
      ! rpn is third to last
      else if ( analysisMemberName(1:nameLength) == 'rpnCalculator' ) then
         call ocn_init_rpn_calculator(domain, err_tmp)
      ! regional is second to last
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsDaily' ) then
         call ocn_init_regional_stats(domain, regionalStatsDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsMonthly' ) then
         call ocn_init_regional_stats(domain, regionalStatsMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsWeekly' ) then
         call ocn_init_regional_stats(domain, regionalStatsWeeklyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsCustom' ) then
         call ocn_init_regional_stats(domain, regionalStatsCustomTAG, err_tmp)
      ! time is last
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call ocn_init_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call ocn_init_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call ocn_init_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call ocn_init_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine ocn_init_analysis_members!}}}

!***********************************************************************
!
!  routine ocn_compute_analysis_members
!
!> \brief Analysis member compute driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct compute routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_compute_analysis_members(domain, timeLevel, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      integer, intent(in) :: timeLevel !< Input: Time level to compute with in analysis member
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'globalStats' ) then
         call ocn_compute_global_stats(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'testComputeInterval' ) then
         call ocn_compute_test_compute_interval(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'layerVolumeWeightedAverage' ) then
         call ocn_compute_layer_volume_weighted_averages(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'meridionalHeatTransport' ) then
         call ocn_compute_meridional_heat_transport(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'okuboWeiss' ) then
         call ocn_compute_okubo_weiss(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
         call ocn_compute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'waterMassCensus' ) then
         call ocn_compute_water_mass_census(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'zonalMean' ) then
         call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call ocn_compute_high_frequency_output(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeFilters' ) then
         call ocn_compute_time_filters(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'lagrPartTrack' ) then
         call ocn_compute_lagrangian_particle_tracking(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'eliassenPalm' ) then
         call ocn_compute_eliassen_palm(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'mixedLayerDepths' ) then
         call ocn_compute_mixed_layer_depths(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call ocn_compute_pointwise_stats(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'debugDiagnostics' ) then
         call ocn_compute_debug_diagnostics(domain, timeLevel, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'transectTransport' ) then
        call ocn_compute_transect_transport(domain, timeLevel, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'eddyProductVariables' ) then
        call ocn_compute_eddy_product_variables(domain, timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'temPlate' ) then
!        call ocn_compute_TEM_PLATE(domain, timeLevel, err_tmp)
      ! rpn is third to last
      else if ( analysisMemberName(1:nameLength) == 'rpnCalculator' ) then
         call ocn_compute_rpn_calculator(domain, timeLevel, err_tmp)
      ! regional is second to last
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsDaily' ) then
         call ocn_compute_regional_stats(domain, timeLevel, regionalStatsDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsMonthly' ) then
         call ocn_compute_regional_stats(domain, timeLevel, regionalStatsMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsWeekly' ) then
         call ocn_compute_regional_stats(domain, timeLevel, regionalStatsWeeklyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsCustom' ) then
         call ocn_compute_regional_stats(domain, timeLevel, regionalStatsCustomTAG, err_tmp)
      ! time is last
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call ocn_compute_time_series_stats(domain, timeLevel, &
           timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call ocn_compute_time_series_stats(domain, timeLevel, &
           timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call ocn_compute_time_series_stats(domain, timeLevel, &
           timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call ocn_compute_time_series_stats(domain, timeLevel, &
           timeSeriesCustomTAG, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine ocn_compute_analysis_members!}}}

!***********************************************************************
!
!  routine ocn_restart_analysis_members
!
!> \brief Analysis member restart driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct restart routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_restart_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'globalStats' ) then
         call ocn_restart_global_stats(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'testComputeInterval' ) then
         call ocn_restart_test_compute_interval(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'layerVolumeWeightedAverage' ) then
         call ocn_restart_layer_volume_weighted_averages(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'meridionalHeatTransport' ) then
         call ocn_restart_meridional_heat_transport(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'okuboWeiss' ) then
         call ocn_restart_okubo_weiss(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
         call ocn_restart_surface_area_weighted_averages(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'waterMassCensus' ) then
         call ocn_restart_water_mass_census(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'zonalMean' ) then
         call ocn_restart_zonal_mean(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call ocn_restart_high_frequency_output(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeFilters' ) then
         call ocn_restart_time_filters(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'lagrPartTrack' ) then
         call ocn_restart_lagrangian_particle_tracking(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'eliassenPalm' ) then
         call ocn_restart_eliassen_palm(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'mixedLayerDepths' ) then
         call ocn_restart_mixed_layer_depths(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call ocn_restart_pointwise_stats(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'debugDiagnostics' ) then
         call ocn_restart_debug_diagnostics(domain, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'transectTransport' ) then
        call ocn_restart_transect_transport(domain, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'eddyProductVariables' ) then
        call ocn_restart_eddy_product_variables(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'temPlate' ) then
!        call ocn_restart_TEM_PLATE(domain, err_tmp)
      ! rpn is third to last
      else if ( analysisMemberName(1:nameLength) == 'rpnCalculator' ) then
         call ocn_restart_rpn_calculator(domain, err_tmp)
      ! regional is second to last
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsDaily' ) then
         call ocn_restart_regional_stats(domain, regionalStatsDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsMonthly' ) then
         call ocn_restart_regional_stats(domain, regionalStatsMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsWeekly' ) then
         call ocn_restart_regional_stats(domain, regionalStatsWeeklyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsCustom' ) then
         call ocn_restart_regional_stats(domain, regionalStatsCustomTAG, err_tmp)
      ! time is last
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call ocn_restart_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call ocn_restart_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call ocn_restart_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call ocn_restart_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine ocn_restart_analysis_members!}}}

!***********************************************************************
!
!  routine ocn_finalize_analysis_members
!
!> \brief Analysis member finalize driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct finalize routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine ocn_finalize_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'globalStats' ) then
         call ocn_finalize_global_stats(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'testComputeInterval' ) then
         call ocn_finalize_test_compute_interval(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'layerVolumeWeightedAverage' ) then
         call ocn_finalize_layer_volume_weighted_averages(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'meridionalHeatTransport' ) then
         call ocn_finalize_meridional_heat_transport(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'okuboWeiss' ) then
         call ocn_finalize_okubo_weiss(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
         call ocn_finalize_surface_area_weighted_averages(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'waterMassCensus' ) then
         call ocn_finalize_water_mass_census(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'zonalMean' ) then
         call ocn_finalize_zonal_mean(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call ocn_finalize_high_frequency_output(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeFilters' ) then
         call ocn_finalize_time_filters(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'lagrPartTrack' ) then
         call ocn_finalize_lagrangian_particle_tracking(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'eliassenPalm' ) then
         call ocn_finalize_eliassen_palm(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'mixedLayerDepths' ) then
         call ocn_finalize_mixed_layer_depths(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call ocn_finalize_pointwise_stats(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'debugDiagnostics' ) then
         call ocn_finalize_debug_diagnostics(domain, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'transectTransport' ) then
        call ocn_finalize_transect_transport(domain, err_tmp)
     else if ( analysisMemberName(1:nameLength) == 'eddyProductVariables' ) then
        call ocn_finalize_eddy_product_variables(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'temPlate' ) then
!        call ocn_finalize_TEM_PLATE(domain, err_tmp)
      ! rpn is third to last
      else if ( analysisMemberName(1:nameLength) == 'rpnCalculator' ) then
         call ocn_finalize_rpn_calculator(domain, err_tmp)
      ! regional is second to last
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsDaily' ) then
         call ocn_finalize_regional_stats(domain, regionalStatsDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsMonthly' ) then
         call ocn_finalize_regional_stats(domain, regionalStatsMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsWeekly' ) then
         call ocn_finalize_regional_stats(domain, regionalStatsWeeklyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatsCustom' ) then
         call ocn_finalize_regional_stats(domain, regionalStatsCustomTAG, err_tmp)
      ! time is last
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call ocn_finalize_time_series_stats(domain, timeSeriesDailyTAG, &
           err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call ocn_finalize_time_series_stats(domain, timeSeriesMonthlyTAG, &
           err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call ocn_finalize_time_series_stats(domain, timeSeriesClimatologyTAG, &
           err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call ocn_finalize_time_series_stats(domain, timeSeriesCustomTAG, &
           err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine ocn_finalize_analysis_members!}}}

end module ocn_analysis_driver

! vim: foldmethod=marker
