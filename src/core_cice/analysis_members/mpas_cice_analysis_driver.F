! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_analysis_driver
!
!> \brief Driver for MPAS ocean analysis core
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This is the driver for the MPAS ocean core.
!
!-----------------------------------------------------------------------

module cice_analysis_driver

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use mpas_stream_manager

   use cice_high_frequency_output
   use cice_temperatures
   use cice_regional_statistics
   use cice_ridging_diagnostics
   use cice_conservation_check
   use cice_geographical_vectors
   use cice_load_balance
   use cice_maximum_ice_presence
   use cice_miscellaneous
   use cice_area_variables
   use cice_pond_diagnostics
   use cice_ice_present
   use cice_time_series_stats
   use cice_pointwise_stats
   use cice_unit_conversion
   use cice_ice_shelves
!   use cice_TEM_PLATE

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

   public :: cice_analysis_setup_packages, &
             cice_analysis_bootstrap, &
             cice_analysis_init, &
             cice_analysis_compute_startup, &
             cice_analysis_precompute, &
             cice_analysis_compute, &
             cice_analysis_write, &
             cice_analysis_restart, &
             cice_analysis_finalize

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
!  routine cice_analysis_setup_packages
!
!> \brief   Setup packages for MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine is intended to configure the packages for all
!>   ocean analysis members.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_setup_packages(configPool, packagePool, iocontext, err)!{{{

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

   end subroutine cice_analysis_setup_packages!}}}

!***********************************************************************
!
!  routine cice_analysis_bootstrap
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

   subroutine cice_analysis_bootstrap(domain, err)!{{{

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

            !write(stderrUnit,*) '      Bootstrapping AM ' // poolItr % memberName(1:nameLength)
            call cice_bootstrap_analysis_members(domain, poolItr % memberName(1:nameLength), ierr=err)

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

   end subroutine cice_analysis_bootstrap!}}}

!***********************************************************************
!
!  routine cice_analysis_init
!
!> \brief   Initialize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all initializations required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_init(domain, err)!{{{

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
            !write(stderrUnit,*) '      Initializing AM ' // poolItr % memberName(1:nameLength)
            timerName = trim(initTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call cice_init_analysis_members(domain, poolItr % memberName, err_tmp)
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

   end subroutine cice_analysis_init!}}}

!***********************************************************************
!
!  routine cice_analysis_compute_startup
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_compute_startup(domain, err)!{{{

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
               !write(stderrUnit,*) '      Computing AM ' // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName)
               call cice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
               err = ior(err, err_tmp)
            end if

            if ( config_AM_write_on_startup ) then
               configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_output_stream'
               call mpas_pool_get_config(domain % configs, configName, config_AM_output_stream)
               if ( config_AM_output_stream /= 'none' ) then
                  !write(stderrUnit,*) '      Writing AM ' // poolItr % memberName(1:nameLength)
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

   end subroutine cice_analysis_compute_startup!}}}

!***********************************************************************
!
!  routine cice_analysis_precompute
!
!> \brief   Driver for MPAS-CICE analysis computations
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine calls all pre timestep computation subroutines
!>  required for the MPAS-CICE analysis driver.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_precompute(domain, err)!{{{

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
                  call cice_precompute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
                  call mpas_timer_stop(timerName)
               end if
            else if ( mpas_is_alarm_ringing(domain % clock, alarmName, ierr=err_tmp) ) then
               !call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_tmp) ! not needed in precompute
               call mpas_timer_start(timerName, .false.)
               call cice_precompute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
            end if
         end if
      end do

      call mpas_timer_stop('analysis_precompute')

   end subroutine cice_analysis_precompute!}}}

!***********************************************************************
!
!  routine cice_analysis_compute
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_compute(domain, err)!{{{

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
                  !write(stderrUnit,*) '      Computing AM ' // poolItr % memberName(1:nameLength)
                  call mpas_timer_start(timerName)
                  call cice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
                  call mpas_timer_stop(timerName)
               end if
            else if ( mpas_is_alarm_ringing(domain % clock, alarmName, ierr=err_tmp) ) then
               call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_tmp)
               !write(stderrUnit,*) '      Computing AM ' // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName)
               call cice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
            end if
         end if
      end do

      call mpas_timer_stop('analysis_compute')

   end subroutine cice_analysis_compute!}}}

!***********************************************************************
!
!  routine cice_analysis_restart
!
!> \brief   Save restart for MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all subroutines required to prepare to save
!>  the restart state for the MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_restart(domain, err)!{{{

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
            !write(stderrUnit,*) '      Preparing AM ' // poolItr % memberName(1:nameLength) // ' for restart write'
            timerName = trim(restartTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call cice_restart_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_restart')

   end subroutine cice_analysis_restart!}}}

!***********************************************************************
!
!  routine cice_analysis_write
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

   subroutine cice_analysis_write(domain, err)!{{{

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
            !write(stderrUnit,*) '      Writing AM ' // poolItr % memberName(1:nameLength)
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
                     write(stderrUnit, *) &
                          'WARNING: Both backward and forward offsets are set for AM ' // poolItr % memberName(1:nameLength)
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

   end subroutine cice_analysis_write!}}}

!***********************************************************************
!
!  routine cice_analysis_finalize
!
!> \brief   Finalize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all finalize routines required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_finalize(domain, err)!{{{

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
            !write(stderrUnit,*) '      Finalizing AM ' // poolItr % memberName(1:nameLength)
            timerName = trim(finalizeTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName)
            call cice_finalize_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_finalize')

   end subroutine cice_analysis_finalize!}}}

!***********************************************************************
!
!  routine cice_bootstrap_analysis_members
!
!> \brief Analysis member initialization driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct init routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine cice_bootstrap_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0
      err_tmp = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call cice_bootstrap_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_bootstrap_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_bootstrap_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_bootstrap_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_bootstrap_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_bootstrap_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call cice_bootstrap_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call cice_bootstrap_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call cice_bootstrap_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call cice_bootstrap_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call cice_bootstrap_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call cice_bootstrap_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call cice_bootstrap_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call cice_bootstrap_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call cice_bootstrap_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call cice_bootstrap_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call cice_bootstrap_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call cice_bootstrap_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call cice_bootstrap_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_bootstrap_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_bootstrap_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine cice_bootstrap_analysis_members!}}}

!***********************************************************************
!
!  routine cice_init_analysis_members
!
!> \brief Analysis member initialization driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct init routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine cice_init_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call cice_init_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_init_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_init_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_init_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_init_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_init_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call cice_init_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call cice_init_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call cice_init_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call cice_init_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call cice_init_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call cice_init_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call cice_init_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call cice_init_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call cice_init_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call cice_init_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call cice_init_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call cice_init_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call cice_init_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_init_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_init_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine cice_init_analysis_members!}}}

!***********************************************************************
!
!  routine cice_precompute_analysis_members
!
!> \brief Analysis member compute driver
!> \author Adrian K. Turner
!> \date 09/09/2015
!> \details
!>  This private routine calls the correct precompute routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine cice_precompute_analysis_members(domain, timeLevel, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      integer, intent(in) :: timeLevel !< Input: Time level to compute with in analysis member
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call cice_precompute_high_frequency_output(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_precompute_temperatures(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_precompute_regional_statistics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_precompute_ridging_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_precompute_conservation_check(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_precompute_geographical_vectors(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call cice_precompute_load_balance(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call cice_precompute_maximum_ice_presence(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call cice_precompute_miscellaneous(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call cice_precompute_area_variables(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call cice_precompute_pond_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call cice_precompute_pointwise_stats(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call cice_precompute_unit_conversion(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call cice_precompute_ice_shelves(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call cice_precompute_ice_present(domain, '', timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
      !   call cice_precompute_time_series_stats(domain, timeSeriesDailyTAG, timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
      !   call cice_precompute_time_series_stats(domain, timeSeriesMonthlyTAG, timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
      !   call cice_precompute_time_series_stats(domain, timeSeriesClimatologyTAG, timeLevel, err_tmp)
      !else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
      !   call cice_precompute_time_series_stats(domain, timeSeriesCustomTAG, timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_precompute_surface_area_weighted_averages(domain, '', timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_precompute_TEMP_LATE(domain, '', timeLevel, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine cice_precompute_analysis_members!}}}

!***********************************************************************
!
!  routine cice_compute_analysis_members
!
!> \brief Analysis member compute driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct compute routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine cice_compute_analysis_members(domain, timeLevel, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      integer, intent(in) :: timeLevel !< Input: Time level to compute with in analysis member
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call cice_compute_high_frequency_output(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_compute_temperatures(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_compute_regional_statistics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_compute_ridging_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_compute_conservation_check(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_compute_geographical_vectors(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call cice_compute_load_balance(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call cice_compute_maximum_ice_presence(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call cice_compute_miscellaneous(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call cice_compute_area_variables(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call cice_compute_pond_diagnostics(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call cice_compute_pointwise_stats(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call cice_compute_unit_conversion(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call cice_compute_ice_shelves(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call cice_compute_ice_present(domain, '', timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call cice_compute_time_series_stats(domain, timeLevel, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call cice_compute_time_series_stats(domain, timeLevel, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call cice_compute_time_series_stats(domain, timeLevel, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call cice_compute_time_series_stats(domain, timeLevel, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_compute_surface_area_weighted_averages(domain, '', timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_compute_TEMP_LATE(domain, '', timeLevel, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine cice_compute_analysis_members!}}}

!***********************************************************************
!
!  routine cice_restart_analysis_members
!
!> \brief Analysis member restart driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct restart routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine cice_restart_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call cice_restart_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_restart_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_restart_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_restart_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_restart_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_restart_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call cice_restart_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call cice_restart_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call cice_restart_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call cice_restart_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call cice_restart_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call cice_restart_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call cice_restart_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call cice_restart_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call cice_restart_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call cice_restart_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call cice_restart_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call cice_restart_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call cice_restart_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_restart_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_restart_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine cice_restart_analysis_members!}}}

!***********************************************************************
!
!  routine cice_finalize_analysis_members
!
!> \brief Analysis member finalize driver
!> \author Doug Jacobsen
!> \date 07/01/2015
!> \details
!>  This private routine calls the correct finalize routine for each analysis member.
!
!-----------------------------------------------------------------------
   subroutine cice_finalize_analysis_members(domain, analysisMemberName, iErr)!{{{
      type (domain_type), intent(inout) :: domain !< Input: Domain information
      character (len=*), intent(in) :: analysisMemberName !< Input: Name of analysis member
      integer, intent(out) :: iErr !< Output: Error code

      integer :: nameLength, err_tmp

      iErr = 0

      nameLength = len_trim(analysisMemberName)

      if ( analysisMemberName(1:nameLength) == 'highFrequencyOutput' ) then
         call cice_finalize_high_frequency_output(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_finalize_temperatures(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_finalize_regional_statistics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_finalize_ridging_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_finalize_conservation_check(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_finalize_geographical_vectors(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'loadBalance' ) then
         call cice_finalize_load_balance(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'maximumIcePresence' ) then
         call cice_finalize_maximum_ice_presence(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'miscellaneous' ) then
         call cice_finalize_miscellaneous(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'areaVariables' ) then
         call cice_finalize_area_variables(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pondDiagnostics' ) then
         call cice_finalize_pond_diagnostics(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'pointwiseStats' ) then
         call cice_finalize_pointwise_stats(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'unitConversion' ) then
         call cice_finalize_unit_conversion(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'iceShelves' ) then
         call cice_finalize_ice_shelves(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'icePresent' ) then
         call cice_finalize_ice_present(domain, '', err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsDaily' ) then
         call cice_finalize_time_series_stats(domain, timeSeriesDailyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsMonthly' ) then
         call cice_finalize_time_series_stats(domain, timeSeriesMonthlyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsClimatology' ) then
         call cice_finalize_time_series_stats(domain, timeSeriesClimatologyTAG, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'timeSeriesStatsCustom' ) then
         call cice_finalize_time_series_stats(domain, timeSeriesCustomTAG, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_finalize_surface_area_weighted_averages(domain, '', err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_finalize_TEMP_LATE(domain, '', err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine cice_finalize_analysis_members!}}}

end module cice_analysis_driver

! vim: foldmethod=marker
