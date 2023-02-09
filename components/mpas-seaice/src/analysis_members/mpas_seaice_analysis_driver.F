!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_analysis_driver
!
!> \brief Driver for MPAS ocean analysis core
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This is the driver for the MPAS ocean core.
!
!-----------------------------------------------------------------------

module seaice_analysis_driver

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use mpas_stream_manager
   use mpas_log, only: mpas_log_write

   use seaice_high_frequency_output
   use seaice_temperatures
   use seaice_thicknesses
   use seaice_regional_statistics
   use seaice_ridging_diagnostics
   use seaice_conservation_check
   use seaice_geographical_vectors
   use seaice_load_balance
   use seaice_maximum_ice_presence
   use seaice_miscellaneous
   use seaice_area_variables
   use seaice_pond_diagnostics
   use seaice_ice_present
   use seaice_time_series_stats
   use seaice_pointwise_stats
   use seaice_unit_conversion
   use seaice_ice_shelves
!   use seaice_TEM_PLATE

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

   public :: seaice_analysis_setup_packages, &
             seaice_analysis_bootstrap, &
             seaice_analysis_init, &
             seaice_analysis_compute_startup, &
             seaice_analysis_precompute, &
             seaice_analysis_compute, &
             seaice_analysis_write, &
             seaice_analysis_restart, &
             seaice_analysis_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


   character (len=*), parameter :: initReadTimerPrefix = 'init_read_'
   character (len=*), parameter :: initTimerPrefix = 'init_'
   character (len=*), parameter :: precomputeTimerPrefix = 'precompute_'
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

   type (mpas_pool_type), pointer :: analysisMemberList

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_analysis_setup_packages
!
!> \brief   Setup packages for MPAS-Seaice analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine is intended to configure the packages for all
!>   ocean analysis members.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_setup_packages(configPool, packagePool, iocontext, err)!{{{

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
      call mpas_pool_add_config(analysisMemberList, 'highFrequencyOutput', 1)
      call mpas_pool_add_config(analysisMemberList, 'temperatures', 1)
      call mpas_pool_add_config(analysisMemberList, 'thicknesses', 1)
      call mpas_pool_add_config(analysisMemberList, 'regionalStatistics', 1)
      call mpas_pool_add_config(analysisMemberList, 'ridgingDiagnostics', 1)
      call mpas_pool_add_config(analysisMemberList, 'conservationCheck', 1)
      call mpas_pool_add_config(analysisMemberList, 'geographicalVectors', 1)
      call mpas_pool_add_config(analysisMemberList, 'loadBalance', 1)
      call mpas_pool_add_config(analysisMemberList, 'maximumIcePresence', 1)
      call mpas_pool_add_config(analysisMemberList, 'miscellaneous', 1)
      call mpas_pool_add_config(analysisMemberList, 'areaVariables', 1)
      call mpas_pool_add_config(analysisMemberList, 'pondDiagnostics', 1)
      call mpas_pool_add_config(analysisMemberList, 'unitConversion', 1)
      call mpas_pool_add_config(analysisMemberList, 'pointwiseStats', 1)
      call mpas_pool_add_config(analysisMemberList, 'iceShelves', 1)
      call mpas_pool_add_config(analysisMemberList, 'icePresent', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeSeriesStatsDaily', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeSeriesStatsMonthly', 1)
      call mpas_pool_add_config(analysisMemberList, 'timeSeriesStatsClimatology', 1)
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

   end subroutine seaice_analysis_setup_packages!}}}

!***********************************************************************
!
!  routine seaice_analysis_bootstrap
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

   subroutine seaice_analysis_bootstrap(domain, err)!{{{

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

            !call mpas_log_write('      Bootstrapping AM ' // poolItr % memberName(1:nameLength))
            call seaice_bootstrap_analysis_members(domain, poolItr % memberName(1:nameLength), ierr=err)

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
                  call mpas_log_write('MPAS-ocean: ERROR: Stream named ''' // trim(config_AM_restart_stream) // &
                                      ''' does not exist in config for analysis member ''' // &
                                      trim(poolItr % memberName(1:nameLength)) // '''', MPAS_LOG_CRIT)
               end if
            end if

            ! Verify the input stream exists
            if ( associated(config_AM_input_stream) ) then
               if ( trim(config_AM_input_stream) == 'none' ) then
                  ! If the stream is set to 'none' nullify the config, so it doesn't get read in
                  nullify(config_AM_input_stream)
               else if ( .not. mpas_stream_mgr_stream_exists(domain % streamManager, config_AM_input_stream) ) then
                  call mpas_log_write('MPAS-ocean: ERROR: Stream named ''' // trim(config_AM_input_stream) // &
                                      ''' does not exist in config for analysis member ''' // &
                                      trim(poolItr % memberName(1:nameLength)) // '''', MPAS_LOG_CRIT)
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

   end subroutine seaice_analysis_bootstrap!}}}

!***********************************************************************
!
!  routine seaice_analysis_init
!
!> \brief   Initialize MPAS-Seaice analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all initializations required for the
!>  MPAS-Seaice analysis driver.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_init(domain, err)!{{{

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
            !call mpas_log_write('      Initializing AM ' // poolItr % memberName(1:nameLength))
            timerName = trim(initTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call seaice_init_analysis_members(domain, poolItr % memberName, err_tmp)
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
                  call mpas_log_write('MPAS-ocean: ERROR: Stream ' // trim(config_AM_output_stream) // &
                                      ' does not exist. Exiting...', MPAS_LOG_CRIT)
               end if
            end if

            if ( config_AM_compute_interval == 'output_interval' .and. config_AM_output_stream == 'none') then
               call mpas_log_write('MPAS-ocean: ERROR: Analysis member has compute_interval of ''output_interval'' ' // &
                                   'without an output stream.', MPAS_LOG_CRIT)
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

   end subroutine seaice_analysis_init!}}}

!***********************************************************************
!
!  routine seaice_analysis_compute_startup
!
!> \brief   Driver for MPAS-Seaice analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-Seaice analysis driver.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_compute_startup(domain, err)!{{{

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
               !call mpas_log_write('      Computing AM ' // poolItr % memberName(1:nameLength))
               call mpas_timer_start(timerName)
               call seaice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
               err = ior(err, err_tmp)
            end if

            if ( config_AM_write_on_startup ) then
               configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_output_stream'
               call mpas_pool_get_config(domain % configs, configName, config_AM_output_stream)
               if ( config_AM_output_stream /= 'none' ) then
                  !call mpas_log_write('      Writing AM ' // poolItr % memberName(1:nameLength))
                  call mpas_stream_mgr_write(domain % streamManager, streamID=config_AM_output_stream, &
                                             forceWriteNow=.true., ierr=err_tmp)
               end if
               if (.not. config_AM_compute_on_startup) then
                  call mpas_log_write(' *** WARNING: write_on_startup called without compute_on_startup for analysis member: ' &
                    // poolItr % memberName(1:nameLength) // '.')
               end if
            end if
         end if
      end do

      call mpas_timer_stop('analysis_compute_startup')

   end subroutine seaice_analysis_compute_startup!}}}

!***********************************************************************
!
!  routine seaice_analysis_precompute
!
!> \brief   Driver for MPAS-Seaice analysis computations
!> \author  MPAS-Seaice/Seaice development team
!> \date    November 2013
!> \details
!>  This routine calls all pre timestep computation subroutines
!>  required for the MPAS-Seaice analysis driver.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_precompute(domain, err)!{{{

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

      call mpas_timer_start('analysis_precompute', .false.)

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
            timerName = trim(precomputeTimerPrefix) // poolItr % memberName(1:nameLength)

            ! Compute analysis member just before output
            if ( config_AM_compute_interval == 'output_interval' .and. config_AM_output_stream /= 'none') then
               if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, &
                    streamID=config_AM_output_stream, direction=MPAS_STREAM_OUTPUT, ierr=err_tmp) ) then
                  call mpas_timer_start(timerName, .false.)
                  call seaice_precompute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
                  call mpas_timer_stop(timerName)
               end if
            else if ( mpas_is_alarm_ringing(domain % clock, alarmName, ierr=err_tmp) ) then
               !call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_tmp) ! not needed in precompute
               call mpas_timer_start(timerName, .false.)
               call seaice_precompute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
            end if
         end if
      end do

      call mpas_timer_stop('analysis_precompute')

   end subroutine seaice_analysis_precompute!}}}

!***********************************************************************
!
!  routine seaice_analysis_compute
!
!> \brief   Driver for MPAS-Seaice analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-Seaice analysis driver.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_compute(domain, err)!{{{

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
                  !call mpas_log_write('      Computing AM ' // poolItr % memberName(1:nameLength))
                  call mpas_timer_start(timerName)
                  call seaice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
                  call mpas_timer_stop(timerName)
               end if
            else if ( mpas_is_alarm_ringing(domain % clock, alarmName, ierr=err_tmp) ) then
               call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_tmp)
               !call mpas_log_write('      Computing AM ' // poolItr % memberName(1:nameLength))
               call mpas_timer_start(timerName)
               call seaice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
            end if
         end if
      end do

      call mpas_timer_stop('analysis_compute')

   end subroutine seaice_analysis_compute!}}}

!***********************************************************************
!
!  routine seaice_analysis_restart
!
!> \brief   Save restart for MPAS-Seaice analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all subroutines required to prepare to save
!>  the restart state for the MPAS-Seaice analysis driver.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_restart(domain, err)!{{{

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
            !call mpas_log_write('      Preparing AM ' // poolItr % memberName(1:nameLength) // ' for restart write')
            timerName = trim(restartTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call seaice_restart_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_restart')

   end subroutine seaice_analysis_restart!}}}

!***********************************************************************
!
!  routine seaice_analysis_write
!
!> \brief   Driver for MPAS-Seaice analysis output
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all output writing subroutines required for the
!>  MPAS-Seaice analysis driver.
!>  At this time this is just a stub, and all analysis output is written
!>  to the output file specified by config_output_name.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_write(domain, err)!{{{

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
            !call mpas_log_write('      Writing AM ' // poolItr % memberName(1:nameLength))
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
                     call mpas_log_write( &
                          'Both backward and forward offsets are set for AM ' // poolItr % memberName(1:nameLength), &
                          MPAS_LOG_WARN)
                     call mpas_log_write('         will only use backward offset. Forward offset will be ignored.', MPAS_LOG_WARN)
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
               call mpas_stream_mgr_write(domain % streamManager, &
                    streamID=config_AM_output_stream, writeTime=outputTimeString, ierr=err_tmp)
               call mpas_timer_stop(timerName)
               timerName = trim(alarmTimerPrefix) // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName)
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=config_AM_output_stream, ierr=err_tmp)
               call mpas_timer_stop(timerName)
            end if
         end if
      end do

      call mpas_timer_stop('analysis_write')

   end subroutine seaice_analysis_write!}}}

!***********************************************************************
!
!  routine seaice_analysis_finalize
!
!> \brief   Finalize MPAS-Seaice analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all finalize routines required for the
!>  MPAS-Seaice analysis driver.
!
!-----------------------------------------------------------------------

   subroutine seaice_analysis_finalize(domain, err)!{{{

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
            !call mpas_log_write('      Finalizing AM ' // poolItr % memberName(1:nameLength))
            timerName = trim(finalizeTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call seaice_finalize_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_finalize')

   end subroutine seaice_analysis_finalize!}}}

!***********************************************************************
!
!  routine seaice_bootstrap_analysis_members
!
!> \brief Analysis member initialization driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct init routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine seaice_bootstrap_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0
      err_tmp = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call seaice_bootstrap_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call seaice_bootstrap_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'thicknesses' ) then
         call seaice_bootstrap_thicknesses(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call seaice_bootstrap_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call seaice_bootstrap_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call seaice_bootstrap_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call seaice_bootstrap_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call seaice_bootstrap_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call seaice_bootstrap_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call seaice_bootstrap_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call seaice_bootstrap_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call seaice_bootstrap_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call seaice_bootstrap_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call seaice_bootstrap_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call seaice_bootstrap_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call seaice_bootstrap_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call seaice_bootstrap_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call seaice_bootstrap_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call seaice_bootstrap_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call seaice_bootstrap_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call seaice_bootstrap_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call seaice_bootstrap_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine seaice_bootstrap_analysis_members!}}}

!***********************************************************************
!
!  routine seaice_init_analysis_members
!
!> \brief Analysis member initialization driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct init routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine seaice_init_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call seaice_init_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call seaice_init_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'thicknesses' ) then
         call seaice_init_thicknesses(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call seaice_init_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call seaice_init_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call seaice_init_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call seaice_init_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call seaice_init_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call seaice_init_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call seaice_init_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call seaice_init_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call seaice_init_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call seaice_init_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call seaice_init_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call seaice_init_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call seaice_init_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call seaice_init_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call seaice_init_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call seaice_init_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call seaice_init_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call seaice_init_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call seaice_init_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine seaice_init_analysis_members!}}}

!***********************************************************************
!
!  routine seaice_precompute_analysis_members
!
!> \brief Analysis member compute driver
!> \author Adrian K. Turner
!> \date 09/09/2015
!> \details
!>  This private routine calls the correct precompute routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine seaice_precompute_analysis_members(domain, timeLevel, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      integer, intent(in) :: timeLevel !< Input: Time level to compute with in analysis member
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      err_tmp = 0 ! needed since time series stats is commented out

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call seaice_precompute_high_frequency_output(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call seaice_precompute_temperatures(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'thicknesses' ) then
         call seaice_precompute_thicknesses(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call seaice_precompute_regional_statistics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call seaice_precompute_ridging_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call seaice_precompute_conservation_check(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call seaice_precompute_geographical_vectors(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call seaice_precompute_load_balance(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call seaice_precompute_maximum_ice_presence(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call seaice_precompute_miscellaneous(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call seaice_precompute_area_variables(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call seaice_precompute_pond_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call seaice_precompute_pointwise_stats(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call seaice_precompute_unit_conversion(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call seaice_precompute_ice_shelves(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call seaice_precompute_ice_present(domain, '', timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
      !   call seaice_precompute_time_series_stats(domain, timeSeriesDailyTAG, timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
      !   call seaice_precompute_time_series_stats(domain, timeSeriesMonthlyTAG, timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
      !   call seaice_precompute_time_series_stats(domain, timeSeriesClimatologyTAG, timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
      !   call seaice_precompute_time_series_stats(domain, timeSeriesCustomTAG, timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call seaice_precompute_surface_area_weighted_averages(domain, '', timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call seaice_precompute_TEMP_LATE(domain, '', timeLevel, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine seaice_precompute_analysis_members!}}}

!***********************************************************************
!
!  routine seaice_compute_analysis_members
!
!> \brief Analysis member compute driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct compute routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine seaice_compute_analysis_members(domain, timeLevel, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      integer, intent(in) :: timeLevel !< Input: Time level to compute with in analysis member
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call seaice_compute_high_frequency_output(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call seaice_compute_temperatures(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'thicknesses' ) then
         call seaice_compute_thicknesses(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call seaice_compute_regional_statistics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call seaice_compute_ridging_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call seaice_compute_conservation_check(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call seaice_compute_geographical_vectors(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call seaice_compute_load_balance(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call seaice_compute_maximum_ice_presence(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call seaice_compute_miscellaneous(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call seaice_compute_area_variables(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call seaice_compute_pond_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call seaice_compute_pointwise_stats(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call seaice_compute_unit_conversion(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call seaice_compute_ice_shelves(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call seaice_compute_ice_present(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call seaice_compute_time_series_stats(domain, timeLevel, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call seaice_compute_time_series_stats(domain, timeLevel, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call seaice_compute_time_series_stats(domain, timeLevel, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call seaice_compute_time_series_stats(domain, timeLevel, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call seaice_compute_surface_area_weighted_averages(domain, '', timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call seaice_compute_TEMP_LATE(domain, '', timeLevel, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine seaice_compute_analysis_members!}}}

!***********************************************************************
!
!  routine seaice_restart_analysis_members
!
!> \brief Analysis member restart driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct restart routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine seaice_restart_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call seaice_restart_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call seaice_restart_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'thicknesses' ) then
         call seaice_restart_thicknesses(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call seaice_restart_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call seaice_restart_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call seaice_restart_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call seaice_restart_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call seaice_restart_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call seaice_restart_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call seaice_restart_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call seaice_restart_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call seaice_restart_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call seaice_restart_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call seaice_restart_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call seaice_restart_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call seaice_restart_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call seaice_restart_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call seaice_restart_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call seaice_restart_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call seaice_restart_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call seaice_restart_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call seaice_restart_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine seaice_restart_analysis_members!}}}

!***********************************************************************
!
!  routine seaice_finalize_analysis_members
!
!> \brief Analysis member finalize driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct finalize routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine seaice_finalize_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call seaice_finalize_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call seaice_finalize_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'thicknesses' ) then
         call seaice_finalize_thicknesses(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call seaice_finalize_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call seaice_finalize_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call seaice_finalize_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call seaice_finalize_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call seaice_finalize_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call seaice_finalize_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call seaice_finalize_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call seaice_finalize_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call seaice_finalize_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call seaice_finalize_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call seaice_finalize_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call seaice_finalize_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call seaice_finalize_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call seaice_finalize_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call seaice_finalize_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call seaice_finalize_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call seaice_finalize_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call seaice_finalize_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call seaice_finalize_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine seaice_finalize_analysis_members!}}}

end module seaice_analysis_driver

! vim: foldmethod=marker
