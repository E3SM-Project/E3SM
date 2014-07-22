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

   use ocn_constants
   use ocn_global_stats
   use ocn_zonal_mean
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
             ocn_analysis_init_alarms, &
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

   ! start alarm IDs at 11 because lower values are assigned in mpas_ocn_mpas_core.F
   character (len=*), parameter :: compute_global_stats_AlarmID = 'globalStats'
   character (len=*), parameter :: compute_zonal_mean_AlarmID = 'zonalMean'
!  character (len=*), parameter :: compute_TEMPLATE_AlarmID = 'templateName'

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
      logical, pointer :: config_use_global_stats, config_use_zonal_mean 
!      logical, pointer :: config_use_TEMPLATE

      err = 0


      call mpas_pool_get_config(configPool, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_setup_packages_global_stats(configPool, packagePool, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(configPool, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_setup_packages_zonal_mean(configPool, packagePool, err_tmp)
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
      logical, pointer :: config_use_global_stats, config_use_zonal_mean
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_init_global_stats(domain, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_init_zonal_mean(domain, err_tmp)
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
!  routine ocn_analysis_init_alarms
!
!> \brief   Initialize alarms MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    March 2014
!> \details 
!>  This routine calls all alarm initializations required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_analysis_init_alarms(clock, startTime, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (MPAS_Time_Type), intent(in) :: startTime

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (MPAS_Clock_type), intent(inout) :: clock

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

      logical, pointer :: config_use_global_stats, config_use_zonal_mean
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_init_alarms_global_stats(clock, startTime, compute_global_stats_AlarmID, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_init_alarms_zonal_mean(clock, startTime, compute_zonal_mean_AlarmID, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(ocnConfigs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         call ocn_init_alarms_TEMPLATE(clock, startTime, compute_TEMPLATE_AlarmID, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_analysis_init_alarms!}}}

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
      logical, pointer :: config_use_zonal_mean, config_zonal_mean_compute_startup
!      logical, pointer :: config_use_TEMPLATE, config_TEMPLATE_compute_startup

      err = 0

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      call mpas_pool_get_config(domain % configs, 'config_global_stats_compute_startup', config_global_stats_compute_startup)
      if (config_use_global_stats.and.config_global_stats_compute_startup) then
         call ocn_compute_global_stats(domain, timeLevel, err_tmp)
         err = ior(err, err_tmp)
      endif 

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      call mpas_pool_get_config(domain % configs, 'config_zonal_mean_compute_startup', config_zonal_mean_compute_startup)
      if (config_use_zonal_mean.and.config_zonal_mean_compute_startup) then
         call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      call mpas_pool_get_config(domain % configs, 'config_TEMPLATE_compute_startup', config_TEMPLATE_compute_startup)
!      if (config_use_TEMPLATE.and.config_TEMPLATE_compute_startup) then
!         call ocn_compute_TEMPLATE(domain, timeLevel, err_tmp)
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

      logical, pointer :: config_use_global_stats, config_use_zonal_mean
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
         call ocn_compute_global_stats(domain, timeLevel, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
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

   subroutine ocn_analysis_compute_w_alarms(clock, domain, err)!{{{

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

      type (MPAS_Clock_type), intent(inout) :: clock
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
      logical, pointer :: config_use_global_stats, config_use_zonal_mean
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      timeLevel=1
      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
         if (mpas_is_alarm_ringing(clock, compute_global_stats_AlarmID, ierr=err_tmp)) then
            call mpas_reset_clock_alarm(clock, compute_global_stats_AlarmID, ierr=err_tmp)
            call ocn_compute_global_stats(domain, timeLevel, err_tmp)
         endif 
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         if (mpas_is_alarm_ringing(clock, compute_zonal_mean_AlarmID, ierr=err_tmp)) then
            call mpas_reset_clock_alarm(clock, compute_zonal_mean_AlarmID, ierr=err_tmp)
            call ocn_compute_zonal_mean(domain, timeLevel, err_tmp)
         endif
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         if (mpas_is_alarm_ringing(clock, compute_TEMPLATE_AlarmID, ierr=err_tmp)) then
!            call mpas_reset_clock_alarm(clock, compute_TEMPLATE_AlarmID, ierr=err_tmp)
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
      logical, pointer :: config_use_global_stats, config_use_zonal_mean
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_restart_global_stats(domain, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_restart_zonal_mean(domain, err_tmp)
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

      logical, pointer :: config_use_global_stats, config_use_zonal_mean
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
      if (config_use_global_stats) then
          call ocn_write_global_stats(domain, err_tmp)
          err = ior(err, err_tmp)
      endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_write_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      call mpas_pool_get_config(domain % configs, 'config_use_TEMPLATE', config_use_TEMPLATE)
!      if (config_use_TEMPLATE) then
!         call ocn_write_TEMPLATE(domain, err_tmp)
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
      logical, pointer :: config_use_global_stats, config_use_zonal_mean
!      logical, pointer :: config_use_TEMPLATE

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_use_global_stats', config_use_global_stats)
     if (config_use_global_stats) then
        call ocn_finalize_global_stats(domain, err_tmp)
        err = ior(err, err_tmp)
     endif

      call mpas_pool_get_config(domain % configs, 'config_use_zonal_mean', config_use_zonal_mean)
      if (config_use_zonal_mean) then
         call ocn_finalize_zonal_mean(domain, err_tmp)
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
