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
!> \brief Driver for MPAS sea ice analysis core
!> \author MPAS-Ocean/CICE development team
!> \date   November 2013
!> \details
!>  This is the driver for the MPAS sea-ice core.
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
!   use cice_TEMP_LATE

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


   character (len=*), parameter :: initTimerPrefix = 'init_'
   character (len=*), parameter :: computeTimerPrefix = 'compute_'
   character (len=*), parameter :: precomputeTimerPrefix = 'precompute_'
   character (len=*), parameter :: writeTimerPrefix = 'write_'
   character (len=*), parameter :: alarmTimerPrefix = 'reset_alarm_'
   character (len=*), parameter :: restartTimerPrefix = 'restart_'
   character (len=*), parameter :: finalizeTimerPrefix = 'finalize_'
   character (len=*), parameter :: computeAlarmSuffix = 'CMPALRM'
   type (mpas_pool_type), pointer :: analysisMemberList

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_analysis_setup_packages
!
!> \brief   Setup packages for MPAS-CICE analysis driver
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine is intended to configure the packages for all
!>   sea-ice analysis members.
!
!-----------------------------------------------------------------------

   subroutine cice_analysis_setup_packages(configPool, packagePool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

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
   !   write(6,*) ' debug:  mpas_pool_add_config'
   !  call mpas_pool_add_config(analysisMemberList, 'surfaceAreaWeightedAverages', 1)
   !  call mpas_pool_add_config(analysisMemberList, 'TEMPLATE', 1)

      ! DON'T EDIT BELOW HERE

      ! Iterate over all analysis members to setup packages
      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         !write(6,*) ' debug : ', configName
         call mpas_pool_get_config(configPool, configName, config_AM_enable)

         if ( config_AM_enable ) then
            packageName = poolItr % memberName(1:nameLength) // 'AMPKGActive'
            !write(6,*) ' debug : ', packageName
            call mpas_pool_get_package(packagePool, packageName, AMPackageActive)
            AMPackageActive = .true.
         end if
      end do

   end subroutine cice_analysis_setup_packages!}}}

!***********************************************************************
!
!  routine cice_analysis_init
!
!> \brief   Initialize MPAS-CICE analysis driver
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine calls all initializations required for the
!>  MPAS-CICE analysis driver.
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
      character (len=StrKIND), pointer :: config_AM_compute_interval, config_AM_stream_name
      integer :: nameLength
      type (mpas_pool_iterator_type) :: poolItr

      logical :: streamFound
      character  (len=StrKIND) :: referenceTimeString, outputIntervalString
      type (MPAS_Time_Type) :: referenceTime 
      type (MPAS_TimeInterval_type) :: alarmTimeStep

      err = 0

      call mpas_timer_start('analysis_init', .false.)

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            timerName = trim(initTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName, .false.)
            call cice_init_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)

            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_compute_interval'
            call mpas_pool_get_config(domain % configs, configName, config_AM_compute_interval)

            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_stream_name'
            call mpas_pool_get_config(domain % configs, configName, config_AM_stream_name)

            if ( config_AM_compute_interval == 'dt' ) then
               alarmTimeStep = mpas_get_clock_timestep(domain % clock, err_tmp)
               call mpas_get_timeInterval(alarmTimeStep, timeString=config_AM_compute_interval, ierr=err_tmp)
            end if

            ! Verify stream exists before trying to use output_interval
            if ( config_AM_stream_name /= 'none' ) then
               streamFound = .false.
   
               call mpas_stream_mgr_begin_iteration(domain % streamManager)
               do while ( mpas_stream_mgr_get_next_stream(domain % streamManager, streamName) )
                  if ( trim(streamName) == trim(config_AM_stream_name) ) then
                     streamFound = .true.
                  end if
               end do
   
               if ( .not. streamFound ) then
                  call mpas_dmpar_global_abort('ERROR: Stream ' // trim(config_AM_stream_name) // ' does not exist. Exiting...')
               end if
            end if
            

            if ( config_AM_compute_interval /= 'output_interval' .and. config_AM_stream_name /= 'none') then
               alarmName = poolItr % memberName(1:nameLength) // computeAlarmSuffix
               call mpas_set_timeInterval(alarmTimeStep, timeString=config_AM_compute_interval, ierr=err_tmp)
               call MPAS_stream_mgr_get_property(domain % streamManager, config_AM_stream_name, MPAS_STREAM_PROPERTY_REF_TIME, referenceTimeString, err_tmp)
               call mpas_set_time(referenceTime, dateTimeString=referenceTimeString, ierr=err_tmp)
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
!> \brief   Driver for MPAS-CICE analysis computations
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-CICE analysis driver.
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
      character (len=StrKIND), pointer :: config_AM_stream_name
      logical, pointer :: config_AM_enable, config_AM_write_on_startup, config_AM_compute_on_startup
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      err = 0

      call mpas_timer_start('analysis_compute', .false.)

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
               timerName = trim(computeTimerPrefix) // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName, .false.)
               call cice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
               call mpas_timer_stop(timerName)
               err = ior(err, err_tmp)

               if ( config_AM_write_on_startup ) then
                  configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_stream_name'
                  call mpas_pool_get_config(domain % configs, configName, config_AM_stream_name)
                  if ( config_AM_stream_name /= 'none' ) then
                     call mpas_stream_mgr_write(domain % streamManager, streamID=config_AM_stream_name, forceWriteNow=.true., ierr=err_tmp)
                  end if
               end if
            else
               if ( config_AM_write_on_startup ) then
                  write(stderrUnit, *) ' *** WARNING: write_on_startup called without compute_on_startup for analysis member: ' &
                                        // poolItr % memberName(1:nameLength) // '. Skipping output...'
               end if
            end if
         end if
      end do

      call mpas_timer_stop('analysis_compute')

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
      character (len=StrKIND), pointer :: config_AM_stream_name, config_AM_compute_interval
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
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_stream_name'
            call mpas_pool_get_config(domain % configs, configName, config_AM_stream_name)

            ! Build name of alarm for analysis member
            alarmName = poolItr % memberName(1:nameLength) // computeAlarmSuffix
            timerName = trim(precomputeTimerPrefix) // poolItr % memberName(1:nameLength)

            ! Compute analysis member just before output
            if ( config_AM_compute_interval == 'output_interval' .and. config_AM_stream_name /= 'none') then
               if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID=config_AM_stream_name, direction=MPAS_STREAM_OUTPUT, ierr=err_tmp) ) then
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
!> \brief   Driver for MPAS-CICE analysis computations
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-CICE analysis driver.
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
      character (len=StrKIND), pointer :: config_AM_stream_name, config_AM_compute_interval
      logical, pointer :: config_AM_enable
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      err = 0

      call mpas_timer_start('analysis_compute', .false.)

      timeLevel=1

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_compute_interval'
            call mpas_pool_get_config(domain % configs, configName, config_AM_compute_interval)
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_stream_name'
            call mpas_pool_get_config(domain % configs, configName, config_AM_stream_name)

            ! Build name of alarm for analysis member
            alarmName = poolItr % memberName(1:nameLength) // computeAlarmSuffix
            timerName = trim(computeTimerPrefix) // poolItr % memberName(1:nameLength)

            ! Compute analysis member just before output
            if ( config_AM_compute_interval == 'output_interval' .and. config_AM_stream_name /= 'none') then
               if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID=config_AM_stream_name, direction=MPAS_STREAM_OUTPUT, ierr=err_tmp) ) then
                  call mpas_timer_start(timerName, .false.)
                  call cice_compute_analysis_members(domain, timeLevel, poolItr % memberName, err_tmp)
                  call mpas_timer_stop(timerName)
               end if
            else if ( mpas_is_alarm_ringing(domain % clock, alarmName, ierr=err_tmp) ) then

               call mpas_reset_clock_alarm(domain % clock, alarmName, ierr=err_tmp)
               call mpas_timer_start(timerName, .false.)
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
!> \brief   Save restart for MPAS-CICE analysis driver
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine calls all subroutines required to prepare to save
!>  the restart state for the MPAS-CICE analysis driver.
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

      call mpas_timer_start('analysis_restart', .false.)

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            timerName = trim(restartTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName, .false.)
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
!> \brief   Driver for MPAS-CICE analysis output
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine calls all output writing subroutines required for the
!>  MPAS-CICE analysis driver.
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

      type (domain_type), intent(inout) :: domain

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

      character (len=StrKIND) :: configName, timerName
      character (len=StrKIND), pointer :: config_AM_stream_name
      logical, pointer :: config_AM_enable
      type (mpas_pool_iterator_type) :: poolItr
      integer :: nameLength

      err = 0

      call mpas_timer_start('analysis_write', .false.)

      call mpas_pool_begin_iteration(analysisMemberList)
      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_stream_name'
            call mpas_pool_get_config(domain % configs, configName, config_AM_stream_name)
            if ( config_AM_stream_name /= 'none' ) then
               timerName = trim(writeTimerPrefix) // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName, .false.)
               call mpas_stream_mgr_write(domain % streamManager, streamID=config_AM_stream_name, ierr=err_tmp)
               call mpas_timer_stop(timerName)
               timerName = trim(alarmTimerPrefix) // poolItr % memberName(1:nameLength)
               call mpas_timer_start(timerName, .false.)
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=config_AM_stream_name, ierr=err_tmp)
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
!> \brief   Finalize MPAS-CICE analysis driver
!> \author  MPAS-Ocean/CICE development team
!> \date    November 2013
!> \details
!>  This routine calls all finalize routines required for the
!>  MPAS-CICE analysis driver.
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

      call mpas_timer_start('analysis_finalize', .false.)

      call mpas_pool_begin_iteration(analysisMemberList)

      do while ( mpas_pool_get_next_member(analysisMemberList, poolItr) )
         nameLength = len_trim(poolItr % memberName)
         configName = 'config_AM_' // poolItr % memberName(1:nameLength) // '_enable'
         call mpas_pool_get_config(domain % configs, configName, config_AM_enable)

         if ( config_AM_enable ) then
            timerName = trim(finalizeTimerPrefix) // poolItr % memberName(1:nameLength)
            call mpas_timer_start(timerName, .false.)
            call cice_finalize_analysis_members(domain, poolItr % memberName, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop(timerName)
         end if
      end do

      call mpas_timer_stop('analysis_finalize')

   end subroutine cice_analysis_finalize!}}}

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
         call cice_init_high_frequency_output(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_init_temperatures(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_init_regional_statistics(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_init_ridging_diagnostics(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_init_conservation_check(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_init_geographical_vectors(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_init_surface_area_weighted_averages(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_init_TEMP_LATE(domain, err_tmp)
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
         call cice_precompute_high_frequency_output(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_precompute_temperatures(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_precompute_regional_statistics(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_precompute_ridging_diagnostics(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_precompute_conservation_check(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_precompute_geographical_vectors(domain, timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_precompute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_precompute_TEMP_LATE(domain, timeLevel, err_tmp)
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
         call cice_compute_high_frequency_output(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_compute_temperatures(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_compute_regional_statistics(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_compute_ridging_diagnostics(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_compute_conservation_check(domain, timeLevel, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_compute_geographical_vectors(domain, timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_compute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_compute_TEMP_LATE(domain, timeLevel, err_tmp)
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
         call cice_restart_high_frequency_output(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_restart_temperatures(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_restart_regional_statistics(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_restart_ridging_diagnostics(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_restart_conservation_check(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_restart_geographical_vectors(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_restart_surface_area_weighted_averages(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_restart_TEMP_LATE(domain, err_tmp)
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
         call cice_finalize_high_frequency_output(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'temperatures' ) then
         call cice_finalize_temperatures(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'regionalStatistics' ) then
         call cice_finalize_regional_statistics(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'ridgingDiagnostics' ) then
         call cice_finalize_ridging_diagnostics(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'conservationCheck' ) then
         call cice_finalize_conservation_check(domain, err_tmp)
      else if ( analysisMemberName(1:nameLength) == 'geographicalVectors' ) then
         call cice_finalize_geographical_vectors(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'surfaceAreaWeightedAverages' ) then
!        call cice_finalize_surface_area_weighted_averages(domain, err_tmp)
!     else if ( analysisMemberName(1:nameLength) == 'TEMPLATE' ) then
!        call cice_finalize_TEMP_LATE(domain, err_tmp)
      end if

      iErr = ior(iErr, err_tmp)

   end subroutine cice_finalize_analysis_members!}}}

!-----------------------------------------------------------------------

end module cice_analysis_driver

! vim: foldmethod=marker
