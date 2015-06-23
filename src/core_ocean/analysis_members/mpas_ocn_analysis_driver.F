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
             ocn_analysis_compute_w_alarms, &
             ocn_analysis_write, &
             ocn_analysis_restart, &
             ocn_analysis_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

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
      logical, pointer :: config_use_AM_water_mass_census
      logical, pointer :: config_use_AM_layer_volume_weighted_avg
      logical, pointer :: config_use_meridional_heat_transport
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

      call mpas_pool_get_config(configPool, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
      if (config_use_AM_layer_volume_weighted_avg) then
          call ocn_setup_packages_layer_volume_weighted_averages(configPool, packagePool, err_tmp)
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

       call mpas_pool_get_config(configPool, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       if (config_use_AM_water_mass_census) then
          call ocn_setup_packages_water_mass_census(configPool, packagePool, err_tmp)
          err = ior(err, err_tmp)
       endif

      call mpas_pool_get_config(configPool, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
      if (config_use_meridional_heat_transport) then
         call ocn_setup_packages_meridional_heat_transport(configPool, packagePool, err_tmp)
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
      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_AM_layer_volume_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_AM_water_mass_census
      logical, pointer :: config_use_meridional_heat_transport
!      logical, pointer :: config_use_TEMPLATE

      err = 0

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

      call mpas_pool_get_config(domain % configs, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
      if (config_use_AM_layer_volume_weighted_avg) then
          call ocn_init_layer_volume_weighted_averages(domain, err_tmp)
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

       call mpas_pool_get_config(domain % configs, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       if (config_use_AM_water_mass_census) then
          call ocn_init_water_mass_census(domain, err_tmp)
          err = ior(err, err_tmp)
       endif

      call mpas_pool_get_config(domain % configs, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
      if (config_use_meridional_heat_transport) then
         call ocn_init_meridional_heat_transport(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         call ocn_init_TEMPLATE(domain, err_tmp)
!         err = ior(err, err_tmp)
!      endif

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

      logical, pointer :: config_use_global_stats, config_global_stats_compute_startup
      logical, pointer :: config_use_AM_sfc_area_weighted_avg, config_AM_sfc_area_weighted_avg_compute_startup
      logical, pointer :: config_use_AM_layer_volume_weighted_avg, config_AM_layer_volume_weighted_avg_compute_startup
      logical, pointer :: config_use_zonal_mean, config_zonal_mean_compute_startup
      logical, pointer :: config_use_okubo_weiss, config_okubo_weiss_compute_startup
      logical, pointer :: config_use_AM_water_mass_census, config_AM_water_mass_census_compute_startup
      logical, pointer :: config_use_meridional_heat_transport, config_meridional_heat_transport_compute_startup
!      logical, pointer :: config_use_TEMPLATE, config_TEMPLATE_compute_startup

      err = 0

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      call mpas_pool_get_config(domain % configs, 'config_global_stats_compute_startup', config_global_stats_compute_startup)
      if (config_use_global_stats.and.config_global_stats_compute_startup) then
         call ocn_compute_global_stats(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(domain % streamManager, streamID='globalStatsOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      call mpas_pool_get_config(domain % configs, 'config_AM_sfc_area_weighted_avg_compute_startup', config_AM_sfc_area_weighted_avg_compute_startup)
      if (config_use_AM_sfc_area_weighted_avg.and.config_AM_sfc_area_weighted_avg_compute_startup) then
         call ocn_compute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(domain % streamManager, streamID='sfcAreaOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
      call mpas_pool_get_config(domain % configs, 'config_AM_layer_volume_weighted_avg_compute_startup', config_AM_layer_volume_weighted_avg_compute_startup)
      if (config_use_AM_layer_volume_weighted_avg.and.config_AM_layer_volume_weighted_avg_compute_startup) then
         call ocn_compute_layer_volume_weighted_averages(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(domain % streamManager, streamID='layerVolOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      call mpas_pool_get_config(domain % configs, 'config_zonal_mean_compute_startup', config_zonal_mean_compute_startup)
      if (config_use_zonal_mean.and.config_zonal_mean_compute_startup) then
         call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(domain % streamManager, streamID='zonalMeanOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      call mpas_pool_get_config(domain % configs, 'config_okubo_weiss_compute_startup', config_okubo_weiss_compute_startup)
      if (config_use_okubo_weiss.and.config_okubo_weiss_compute_startup) then
         call ocn_compute_okubo_weiss(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(domain % streamManager, streamID='okuboWeissOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

       call mpas_pool_get_config(domain % configs, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       call mpas_pool_get_config(domain % configs, 'config_AM_water_mass_census_compute_startup', config_AM_water_mass_census_compute_startup)
       if (config_use_AM_water_mass_census.and.config_AM_water_mass_census_compute_startup) then
          call ocn_compute_water_mass_census(domain, timeLevel, err_tmp)
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(domain % streamManager, streamID='waterMassOutput', forceWriteNow=.true., ierr=err_tmp)
          call mpas_timer_stop('io_write')
          err = ior(err, err_tmp)
       endif

      call mpas_pool_get_config(domain % configs, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
      call mpas_pool_get_config(domain % configs, 'config_meridional_heat_transport_compute_startup', config_meridional_heat_transport_compute_startup)
      if (config_use_meridional_heat_transport.and.config_meridional_heat_transport_compute_startup) then
         call ocn_compute_meridional_heat_transport(domain, timeLevel, err_tmp)
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(domain % streamManager, streamID='MerHeatTransOutput', forceWriteNow=.true., ierr=err_tmp)
         call mpas_timer_stop('io_write')
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      call mpas_pool_get_config(domain % configs, 'config_TEMPLATE_compute_startup', config_TEMPLATE_compute_startup)
!      if (config_use_TEMPLATE.and.config_TEMPLATE_compute_startup) then
!         call mpas_timer_start('io_write', .false.)
!         call ocn_compute_TEMPLATE(domain, timeLevel, err_tmp)
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

      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_AM_layer_volume_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_AM_water_mass_census
      logical, pointer :: config_use_meridional_heat_transport
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
         call ocn_compute_global_stats(domain, timeLevel, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
         call ocn_compute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
      if (config_use_AM_layer_volume_weighted_avg) then
         call ocn_compute_layer_volume_weighted_averages(domain, timeLevel, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         call ocn_compute_okubo_weiss(domain, timeLevel, err_tmp)
      endif

       call mpas_pool_get_config(domain % configs, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       if (config_use_AM_water_mass_census) then
          call ocn_compute_water_mass_census(domain, timeLevel, err_tmp)
       endif

      call mpas_pool_get_config(domain % configs, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
      if (config_use_meridional_heat_transport) then
         call ocn_compute_meridional_heat_transport(domain, timeLevel, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         call ocn_compute_TEMPLATE(domain, timeLevel, err_tmp)
!      endif

   end subroutine ocn_analysis_compute!}}}

!***********************************************************************
!
!  routine ocn_analysis_compute_w_alarms
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details
!>  This routine calls all computation subroutines required for the
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_compute_w_alarms(domain, err)!{{{

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
      logical, pointer :: config_use_global_stats
      logical, pointer :: config_use_AM_sfc_area_weighted_avg
      logical, pointer :: config_use_AM_layer_volume_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_AM_water_mass_census
      logical, pointer :: config_use_meridional_heat_transport
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
         if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='globalStatsOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_global_stats(domain, timeLevel, err_tmp)
         endif
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
         if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='sfcAreaOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_surface_area_weighted_averages(domain, timeLevel, err_tmp)
         endif
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
      if (config_use_AM_layer_volume_weighted_avg) then
         if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='layerVolOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_layer_volume_weighted_averages(domain, timeLevel, err_tmp)
         endif
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='zonalMeanOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
         endif
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='okuboWeissOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_okubo_weiss(domain, timeLevel, err_tmp)
         endif
      endif

       call mpas_pool_get_config(domain % configs, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       if (config_use_AM_water_mass_census) then
         if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='waterMassOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
             call ocn_compute_water_mass_census(domain, timeLevel, err_tmp)
          endif
       endif

      call mpas_pool_get_config(domain % configs, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
      if (config_use_meridional_heat_transport) then
        if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='MerHeatTransOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)) then
            call ocn_compute_meridional_heat_transport(domain, timeLevel, err_tmp)
         endif
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!        if (mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='TEMPLATEOutput', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp) then
!            call ocn_compute_TEMPLATE(domain, timeLevel, err_tmp)
!         endif
!      endif

   end subroutine ocn_analysis_compute_w_alarms!}}}

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
      logical, pointer :: config_use_AM_layer_volume_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_AM_water_mass_census
      logical, pointer :: config_use_meridional_heat_transport
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

      call mpas_pool_get_config(domain % configs, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
      if (config_use_AM_layer_volume_weighted_avg) then
          call ocn_restart_layer_volume_weighted_averages(domain, err_tmp)
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

       call mpas_pool_get_config(domain % configs, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       if (config_use_AM_water_mass_census) then
          call ocn_restart_water_mass_census(domain, err_tmp)
          err = ior(err, err_tmp)
       endif

      call mpas_pool_get_config(domain % configs, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
      if (config_use_meridional_heat_transport) then
         call ocn_restart_meridional_heat_transport(domain, err_tmp)
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

   subroutine ocn_analysis_write(streamManager, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (MPAS_streamManager_type), intent(inout) :: streamManager

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
      logical, pointer :: config_use_AM_layer_volume_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_AM_water_mass_census
      logical, pointer :: config_use_meridional_heat_transport
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(streamManager, streamID='globalStatsOutput', ierr=err_tmp)
          call mpas_timer_stop('io_write')
          call mpas_timer_start('io_reset_alarms', .false.)
          call mpas_stream_mgr_reset_alarms(streamManager, streamID='globalStatsOutput', ierr=err_tmp)
          call mpas_timer_stop('io_reset_alarms')
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_AM_sfc_area_weighted_avg', config_use_AM_sfc_area_weighted_avg)
      if (config_use_AM_sfc_area_weighted_avg) then
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(streamManager, streamID='sfcAreaOutput', ierr=err_tmp)
          call mpas_timer_stop('io_write')
          call mpas_timer_start('io_reset_alarms', .false.)
          call mpas_stream_mgr_reset_alarms(streamManager, streamID='sfcAreaOutput', ierr=err_tmp)
          call mpas_timer_stop('io_reset_alarms')
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
      if (config_use_AM_layer_volume_weighted_avg) then
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(streamManager, streamID='layerVolOutput', ierr=err_tmp)
          call mpas_timer_stop('io_write')
          call mpas_timer_start('io_reset_alarms', .false.)
          call mpas_stream_mgr_reset_alarms(streamManager, streamID='layerVolOutput', ierr=err_tmp)
          call mpas_timer_stop('io_reset_alarms')
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(streamManager, streamID='zonalMeanOutput', ierr=err_tmp)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('io_reset_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(streamManager, streamID='zonalMeanOutput', ierr=err_tmp)
         call mpas_timer_stop('io_reset_alarms')
         err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_okubo_weiss', config_use_okubo_weiss)
      if (config_use_okubo_weiss) then
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(streamManager, streamID='okuboWeissOutput', ierr=err_tmp)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('io_reset_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(streamManager, streamID='okuboWeissOutput', ierr=err_tmp)
         call mpas_timer_stop('io_reset_alarms')
         err = ior(err, err_tmp)
      endif

       call mpas_pool_get_config(ocnConfigs, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       if (config_use_AM_water_mass_census) then
           call mpas_timer_start('io_write', .false.)
           call mpas_stream_mgr_write(streamManager, streamID='waterMassOutput', ierr=err_tmp)
           call mpas_timer_stop('io_write')
           call mpas_timer_start('io_reset_alarms', .false.)
           call mpas_stream_mgr_reset_alarms(streamManager, streamID='waterMassOutput', ierr=err_tmp)
           call mpas_timer_stop('io_reset_alarms')
          err = ior(err, err_tmp)
       endif

     call mpas_pool_get_config(ocnConfigs, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
     if (config_use_meridional_heat_transport) then
         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(streamManager, streamID='MerHeatTransOutput', ierr=err_tmp)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('io_reset_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(streamManager, streamID='MerHeatTransOutput', ierr=err_tmp)
         call mpas_timer_stop('io_reset_alarms')
        err = ior(err, err_tmp)
     endif

!      call mpas_pool_get_config(ocnConfigs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!          call mpas_timer_start('io_write', .false.)
!          call mpas_stream_mgr_write(streamManager, streamID='TEMPLATEOutput', ierr=err_tmp)
!          call mpas_timer_stop('io_write')
!          call mpas_timer_start('io_reset_alarms', .false.)
!          call mpas_stream_mgr_reset_alarms(streamManager, streamID='TEMPLATEOutput', ierr=err_tmp)
!          call mpas_timer_stop('io_reset_alarms')
!         err = ior(err, err_tmp)
!      endif

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
      logical, pointer :: config_use_AM_layer_volume_weighted_avg
      logical, pointer :: config_use_zonal_mean
      logical, pointer :: config_use_okubo_weiss
      logical, pointer :: config_use_AM_water_mass_census
      logical, pointer :: config_use_meridional_heat_transport
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

     call mpas_pool_get_config(domain % configs, 'config_use_AM_layer_volume_weighted_avg', config_use_AM_layer_volume_weighted_avg)
     if (config_use_AM_layer_volume_weighted_avg) then
        call ocn_finalize_layer_volume_weighted_averages(domain, err_tmp)
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

       call mpas_pool_get_config(domain % configs, 'config_use_AM_water_mass_census', config_use_AM_water_mass_census)
       if (config_use_AM_water_mass_census) then
          call ocn_finalize_water_mass_census(domain, err_tmp)
          err = ior(err, err_tmp)
       endif

      call mpas_pool_get_config(domain % configs, 'config_use_meridional_heat_transport', config_use_meridional_heat_transport)
      if (config_use_meridional_heat_transport) then
         call ocn_finalize_meridional_heat_transport(domain, err_tmp)
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
