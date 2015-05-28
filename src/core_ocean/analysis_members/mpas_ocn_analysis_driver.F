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

   use mpas_grid_types
   use mpas_timekeeping
   use mpas_timer
   use mpas_stream_manager

   use ocn_constants
   use ocn_global_stats
   use ocn_surface_area_weighted_averages
   use ocn_zonal_mean
   use ocn_okubo_weiss
   use ocn_test_compute_interval
!   use ocn_TEMPLATE

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

   character (len=*), parameter :: test_compute_intervalAlarmID = 'test_compute_intervalAlarm'

   ! To use a compute interval using config_TEMPLATE_compute_interval, copy the following.
   ! Otherwise, compute is only called just before write, as managed by the stream.
   !character (len=*), parameter :: TEMPLATEAlarmID = 'TEMPLATEAlarm'

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

   subroutine ocn_analysis_setup_packages(configPool, packagePool, err)!{{{

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
      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_zonal_mean 
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_test_compute_interval
!      logical, pointer :: config_use_TEMPLATE

      err = 0


      call mpas_pool_get_config(configPool, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_setup_packages_global_stats(configPool, packagePool, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(configPool, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
          call ocn_setup_packages_surface_area_weighted_averages(configPool, packagePool, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(configPool, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_setup_packages_zonal_mean(configPool, packagePool, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(configPool, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         call ocn_setup_packages_okubo_weiss(configPool, packagePool, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(configPool, 'config_use_test_compute_interval', config_use_test_compute_interval)
      if (config_use_test_compute_interval) then
         call ocn_setup_packages_test_compute_interval(configPool, packagePool, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(configPool, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         call ocn_setup_packages_TEMPLATE(configPool, packagePool, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_analysis_setup_packages!}}}

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

   subroutine ocn_analysis_init(domain, stream_manager, err)!{{{

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
      type (mpas_streamManager_type), intent(inout) :: stream_manager

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
      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_test_compute_interval
!      logical, pointer :: config_use_TEMPLATE

      type (MPAS_Clock_type), pointer :: clock
      character  (len=StrKIND) :: referenceTimeString
      type (MPAS_Time_Type) :: referenceTime 
      type (MPAS_TimeInterval_type) :: alarmTimeStep

      character (len=StrKIND), pointer :: config_test_compute_interval_compute_interval
   ! To use a compute interval using config_TEMPLATE_compute_interval, copy the following.
   ! Otherwise, compute is only called just before write, as managed by the stream.
!      character (len=StrKIND), pointer :: config_TEMPLATE_compute_interval

      err = 0

      clock => domain % clock

      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_init_global_stats(domain, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
          call ocn_init_surface_area_weighted_averages(domain, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_init_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         call ocn_init_okubo_weiss(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_test_compute_interval', config_use_test_compute_interval)
      if (config_use_test_compute_interval) then
         call ocn_init_test_compute_interval(domain, err_tmp)
         err = ior(err, err_tmp)

         call mpas_pool_get_config(domain % configs, 'config_test_compute_interval_compute_interval',config_test_compute_interval_compute_interval)
         if (config_test_compute_interval_compute_interval=='dt') config_test_compute_interval_compute_interval='00-00-00_00:00:01'
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_test_compute_interval_compute_interval, ierr=err_tmp)

         call MPAS_stream_mgr_get_property(stream_manager, 'testCompIntOutput', MPAS_STREAM_PROPERTY_REF_TIME, referenceTimeString, err_tmp)
         call mpas_set_time(referenceTime, dateTimeString=referenceTimeString, ierr=err_tmp)
         call mpas_add_clock_alarm(clock, test_compute_intervalAlarmID, referenceTime, alarmTimeStep, ierr=err_tmp)
         call mpas_stream_mgr_reset_alarms(stream_manager, streamID='testCompIntOutput', ierr=err_tmp)

      endif

      ! call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
      ! if (config_use_TEMPLATE) then
      !    call ocn_init_TEMPLATE(domain, err_tmp)
      !    err = ior(err, err_tmp)

      !    ! To use a compute interval using config_TEMPLATE_compute_interval, uncomment the following.
      !    ! Otherwise, compute is only called just before write, as managed by the stream.
      !    call mpas_pool_get_config(domain % configs, 'config_TEMPLATE_compute_interval',config_TEMPLATE_compute_interval)
      !    if (config_TEMPLATE_compute_interval=='dt') config_TEMPLATE_compute_interval='00-00-00_00:00:01'
      !    call mpas_set_timeInterval(alarmTimeStep, timeString=config_TEMPLATE_compute_interval, ierr=err_tmp)

      !    call MPAS_stream_mgr_get_property(stream_manager, 'testCompIntOutput', MPAS_STREAM_PROPERTY_REF_TIME, referenceTimeString, err_tmp)
      !    call mpas_set_time(referenceTime, dateTimeString=referenceTimeString, ierr=err_tmp)
      !    call mpas_add_clock_alarm(clock, TEMPLATEAlarmID, referenceTime, alarmTimeStep, ierr=err_tmp)
      !    call mpas_stream_mgr_reset_alarms(stream_manager, streamID='TEMPLATEOutput', ierr=err_tmp)

      ! endif

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

   subroutine ocn_analysis_compute_startup(domain, stream_manager, err)!{{{

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
      type (mpas_streamManager_type), intent(inout) :: stream_manager

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

      logical, pointer :: config_use_global_stats, config_global_stats_compute_startup
      logical, pointer :: config_use_AM_sfc_area_weighted_avg, config_AM_sfc_area_weighted_avg_compute_startup
      logical, pointer :: config_use_zonal_mean, config_zonal_mean_compute_startup
      logical, pointer :: config_use_okubo_weiss, config_okubo_weiss_compute_startup
      logical, pointer :: config_use_test_compute_interval, config_test_compute_interval_compute_startup
!      logical, pointer :: config_use_TEMPLATE, config_TEMPLATE_compute_startup

      err = 0

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      call mpas_pool_get_config(domain % configs, 'config_global_stats_compute_startup', config_global_stats_compute_startup)
      if (config_use_global_stats.and.config_global_stats_compute_startup) then
         call ocn_compute_global_stats(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='globalStatsOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif 

      call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      call mpas_pool_get_config(domain % configs, 'config_AM_sfc_area_weighted_avg_compute_startup', config_AM_sfc_area_weighted_avg_compute_startup)
      if (config_use_AM_sfc_area_weighted_avg.and.config_AM_sfc_area_weighted_avg_compute_startup) then
         call ocn_compute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='sfcAreaOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      call mpas_pool_get_config(domain % configs, 'config_zonal_mean_compute_startup', config_zonal_mean_compute_startup)
      if (config_use_zonal_mean.and.config_zonal_mean_compute_startup) then
         call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='zonalMeanOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      call mpas_pool_get_config(domain % configs, 'config_okubo_weiss_compute_startup', config_okubo_weiss_compute_startup)
      if (config_use_okubo_weiss.and.config_okubo_weiss_compute_startup) then
         call ocn_compute_okubo_weiss(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='okuboWeissOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_test_compute_interval', config_use_test_compute_interval)
      call mpas_pool_get_config(domain % configs, 'config_test_compute_interval_compute_startup', config_test_compute_interval_compute_startup)
      if (config_use_test_compute_interval.and.config_test_compute_interval_compute_startup) then
         call ocn_compute_test_compute_interval(domain, timeLevel, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      call mpas_pool_get_config(domain % configs, 'config_TEMPLATE_compute_startup', config_TEMPLATE_compute_startup)
!      if (config_use_TEMPLATE.and.config_TEMPLATE_compute_startup) then

!         ! This line is to compute on start-up:
!         call ocn_compute_TEMPLATE(domain, timeLevel, err_tmp)

!         ! This line is to write on start-up:
!         call mpas_timer_start('io_write', .false.)
!         call mpas_stream_mgr_write(stream_manager, streamID='TEMPLATEOutput', forceWriteNow=.true., ierr=err_tmp)
!         call mpas_timer_stop('io_write')
!         err = ior(err, err_tmp)
!      endif

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

   subroutine ocn_analysis_compute(domain, stream_manager, err)!{{{

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

      type (MPAS_streamManager_type), intent(inout) :: stream_manager
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
      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_test_compute_interval
!      logical, pointer :: config_use_TEMPLATE
      type (MPAS_Clock_type), pointer :: clock

      err = 0
      clock => domain % clock

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
         ! Compute analysis member just before output:
         if (mpas_stream_mgr_ringing_alarms(stream_manager, streamID='globalStatsOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_global_stats(domain, timeLevel, err_tmp)
         endif 
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
         ! Compute analysis member just before output:
         if (mpas_stream_mgr_ringing_alarms(stream_manager, streamID='sfcAreaOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
         endif
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         ! Compute analysis member just before output:
         if (mpas_stream_mgr_ringing_alarms(stream_manager, streamID='zonalMeanOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
         endif
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         ! Compute analysis member just before output:
         if (mpas_stream_mgr_ringing_alarms(stream_manager, streamID='okuboWeissOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_okubo_weiss(domain, timeLevel, err_tmp)
         endif
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_test_compute_interval', config_use_test_compute_interval)
      if (config_use_test_compute_interval) then
         ! Compute analysis member when called for compute:
         if (mpas_is_alarm_ringing(clock, test_compute_intervalAlarmID, ierr=err_tmp)) then
            call mpas_reset_clock_alarm(clock, test_compute_intervalAlarmID, ierr=err_tmp)
            call ocn_compute_test_compute_interval(domain, timeLevel, err_tmp)
         endif

         ! Compute analysis member just before output:
         !if (mpas_stream_mgr_ringing_alarms(stream_manager, streamID='test_compute_intervalOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
         !   call ocn_compute_test_compute_interval(domain, timeLevel, err_tmp)
         !endif

      endif

      ! call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
      ! if (config_use_TEMPLATE) then
           ! When adding a new analysis member, you may choose when to compute: 
           ! 1. using config_TEMPLATE_compute_interval
           ! 2. just before stream output
           ! or both.  Uncomment these sections accordingly.

      !    ! Compute analysis member when called for compute, using config_TEMPLATE_compute_interval
      !    if (mpas_is_alarm_ringing(clock, TEMPLATEAlarmID, ierr=err_tmp)) then
      !       call mpas_reset_clock_alarm(clock, TEMPLATEAlarmID, ierr=err_tmp)
      !       call ocn_compute_TEMPLATE(domain, timeLevel, err_tmp)
      !    endif

      !    ! Compute analysis member just before stream output:
      !    !if (mpas_stream_mgr_ringing_alarms(stream_manager, streamID='TEMPLATEOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
      !    !   call ocn_compute_TEMPLATE(domain, timeLevel, err_tmp)
      !    !endif

      ! endif

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
      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_test_compute_interval
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_restart_global_stats(domain, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
          call ocn_restart_surface_area_weighted_averages(domain, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_restart_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         call ocn_restart_okubo_weiss(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_test_compute_interval', config_use_test_compute_interval)
      if (config_use_test_compute_interval) then
         call ocn_restart_test_compute_interval(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         call ocn_restart_TEMPLATE(domain, err_tmp)
!         err = ior(err, err_tmp)
!      endif

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

   subroutine ocn_analysis_write(stream_manager, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (MPAS_streamManager_type), intent(inout) :: stream_manager

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

      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg 
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_test_compute_interval
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(stream_manager, streamID='globalStatsOutput', ierr=err_tmp)
          call mpas_timer_stop('io_write')
          call mpas_timer_start('io_reset_alarms', .false.)
          call mpas_stream_mgr_reset_alarms(stream_manager, streamID='globalStatsOutput', ierr=err_tmp)
          call mpas_timer_stop('io_reset_alarms')
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(stream_manager, streamID='sfcAreaOutput', ierr=err_tmp)
          call mpas_timer_stop('io_write')
          call mpas_timer_start('io_reset_alarms', .false.)
          call mpas_stream_mgr_reset_alarms(stream_manager, streamID='sfcAreaOutput', ierr=err_tmp)
          call mpas_timer_stop('io_reset_alarms')
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='zonalMeanOutput', ierr=err_tmp)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('io_reset_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(stream_manager, streamID='zonalMeanOutput', ierr=err_tmp)
         call mpas_timer_stop('io_reset_alarms')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='okuboWeissOutput', ierr=err_tmp)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('io_reset_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(stream_manager, streamID='okuboWeissOutput', ierr=err_tmp)
         call mpas_timer_stop('io_reset_alarms')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_test_compute_interval', config_use_test_compute_interval)
      if (config_use_test_compute_interval) then
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(stream_manager, streamID='test_compute_intervalOutput', ierr=err_tmp)
          call mpas_timer_stop('io_write')
          call mpas_timer_start('io_reset_alarms', .false.)
          call mpas_stream_mgr_reset_alarms(stream_manager, streamID='test_compute_intervalOutput', ierr=err_tmp)
          call mpas_timer_stop('io_reset_alarms')
         err = ior(err, err_tmp)
      endif

      ! call mpas_pool_get_config(ocnConfigs, 'config_use_TEMPLATE', config_use_TEMPLATE)
      ! if (config_use_TEMPLATE) then
      !     call mpas_timer_start('io_write', .false.)
      !     call mpas_stream_mgr_write(stream_manager, streamID='TEMPLATEOutput', ierr=err_tmp)
      !     call mpas_timer_stop('io_write')
      !     call mpas_timer_start('io_reset_alarms', .false.)
      !     call mpas_stream_mgr_reset_alarms(stream_manager, streamID='TEMPLATEOutput', ierr=err_tmp)
      !     call mpas_timer_stop('io_reset_alarms')
      !    err = ior(err, err_tmp)
      ! endif

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
      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_test_compute_interval
!      logical, pointer :: config_use_TEMPLATE

      err = 0

     call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
     if (config_use_global_stats) then
        call ocn_finalize_global_stats(domain, err_tmp)
        err = ior(err, err_tmp)
     endif

     call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
     if (config_use_AM_sfc_area_weighted_avg) then
        call ocn_finalize_surface_area_weighted_averages(domain, err_tmp)
        err = ior(err, err_tmp)
     endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_finalize_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         call ocn_finalize_okubo_weiss(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_test_compute_interval', config_use_test_compute_interval)
      if (config_use_test_compute_interval) then
         call ocn_finalize_test_compute_interval(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         call ocn_finalize_TEMPLATE(domain, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_analysis_finalize!}}}

end module ocn_analysis_driver

! vim: foldmethod=marker
