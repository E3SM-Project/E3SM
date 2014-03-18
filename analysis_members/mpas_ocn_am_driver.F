! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_am_driver
!
!> \brief Driver for MPAS ocean analysis core
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This is the driver for the MPAS ocean core.
!
!-----------------------------------------------------------------------

module ocn_am_driver

   use mpas_grid_types
   use mpas_configure
   use mpas_timekeeping

   use ocn_am_global_stats
   use ocn_am_zonal_mean
!   use ocn_am_TEMPLATE

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

   public :: ocn_am_driver_setup_packages, &
             ocn_am_driver_init, &
             ocn_am_driver_init_alarms, &
             ocn_am_driver_compute_startup, &
             ocn_am_driver_compute, &
             ocn_am_driver_compute_w_alarms, &
             ocn_am_driver_write, &
             ocn_am_driver_restart, &
             ocn_am_driver_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   ! start alarm IDs at 11 because lower values are assigned in mpas_ocn_mpas_core.F
   integer, parameter :: compute_global_stats_AlarmID = 11
   integer, parameter :: compute_zonal_mean_AlarmID = 12
!   integer, parameter :: compute_TEMPLATE_AlarmID = (increment previous)

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_am_driver_setup_packages
!
!> \brief   Setup packages for MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine is intended to configure the packages for all
!>   ocean analysis members.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_setup_packages(err)!{{{

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

      err = 0

     if (config_ocn_am_global_stats) then
         call ocn_am_setup_packages_global_stats(err_tmp)
         err = ior(err, err_tmp)
     endif

      if (config_ocn_am_zonal_mean) then
         call ocn_am_setup_packages_zonal_mean(err_tmp)
         err = ior(err, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE) then
!         call ocn_am_setup_packages_TEMPLATE(err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_am_driver_setup_packages!}}}

!***********************************************************************
!
!  routine ocn_am_driver_init
!
!> \brief   Initialize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all initializations required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_init(domain, err)!{{{

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

      err = 0

     if (config_ocn_am_global_stats) then
         call ocn_am_init_global_stats(domain, err_tmp)
         err = ior(err, err_tmp)
     endif

      if (config_ocn_am_zonal_mean) then
         call ocn_am_init_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE) then
!         call ocn_am_init_TEMPLATE(domain, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_am_driver_init!}}}

!***********************************************************************
!
!  routine ocn_am_driver_init_alarms
!
!> \brief   Initialize alarms MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    March 2014
!> \details 
!>  This routine calls all alarm initializations required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_init_alarms(clock, startTime, err)!{{{

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

      err = 0

     if (config_ocn_am_global_stats) then
         call ocn_am_init_alarms_global_stats(clock, startTime, compute_global_stats_AlarmID, err_tmp)
         err = ior(err, err_tmp)
     endif

      if (config_ocn_am_zonal_mean) then
         call ocn_am_init_alarms_zonal_mean(clock, startTime, compute_zonal_mean_AlarmID, err_tmp)
         err = ior(err, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE) then
      !         call ocn_am_init_alarms_TEMPLATE(clock, startTime, compute_TEMPLATE_AlarmID, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_am_driver_init_alarms!}}}

!***********************************************************************
!
!  routine ocn_am_driver_compute_startup
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all computation subroutines required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_compute_startup(domain, err)!{{{

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
      type (state_type), pointer :: state

      err = 0

      timeLevel=1
      if (config_ocn_am_global_stats.and.config_ocn_am_global_stats_compute_startup) then
         call ocn_am_compute_global_stats(domain, timeLevel, err_tmp)
         err = ior(err, err_tmp)
      endif 

      if (config_ocn_am_zonal_mean.and.config_ocn_am_zonal_mean_compute_startup) then
         call ocn_am_compute_zonal_mean(domain, timeLevel, err_tmp)
         err = ior(err, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE.and.config_ocn_am_TEMPLATE_compute_startup) then
!         call ocn_am_compute_TEMPLATE(domain, timeLevel, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_am_driver_compute_startup!}}}

!***********************************************************************
!
!  routine ocn_am_driver_compute
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all computation subroutines required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_compute(domain, err)!{{{

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
      type (state_type), pointer :: state

      err = 0

      timeLevel=1
      if (config_ocn_am_global_stats) then
         call ocn_am_compute_global_stats(domain, timeLevel, err_tmp)
      endif

      if (config_ocn_am_zonal_mean) then
         call ocn_am_compute_zonal_mean(domain, timeLevel, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE) then
!         call ocn_am_compute_TEMPLATE(domain, timeLevel, err_tmp)
!      endif

   end subroutine ocn_am_driver_compute!}}}

!***********************************************************************
!
!  routine ocn_am_driver_compute_w_alarms
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all computation subroutines required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_compute_w_alarms(clock, domain, err)!{{{

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
      type (state_type), pointer :: state

      err = 0

      timeLevel=1
      if (config_ocn_am_global_stats) then
         if (mpas_is_alarm_ringing(clock, compute_global_stats_AlarmID, ierr=err_tmp)) then
            call mpas_reset_clock_alarm(clock, compute_global_stats_AlarmID, ierr=err_tmp)
            call ocn_am_compute_global_stats(domain, timeLevel, err_tmp)
         endif 
      endif

      if (config_ocn_am_zonal_mean) then
         if (mpas_is_alarm_ringing(clock, compute_zonal_mean_AlarmID, ierr=err_tmp)) then
            call mpas_reset_clock_alarm(clock, compute_zonal_mean_AlarmID, ierr=err_tmp)
            call ocn_am_compute_zonal_mean(domain, timeLevel, err_tmp)
         endif
      endif

!      if (config_ocn_am_TEMPLATE) then
!         if (mpas_is_alarm_ringing(clock, compute_TEMPLATE_AlarmID, ierr=err_tmp)) then
!            call mpas_reset_clock_alarm(clock, compute_TEMPLATE_AlarmID, ierr=err_tmp)
!            call ocn_am_compute_TEMPLATE(domain, timeLevel, err_tmp)
!         endif
!      endif

   end subroutine ocn_am_driver_compute_w_alarms!}}}

!***********************************************************************
!
!  routine ocn_am_driver_restart
!
!> \brief   Save restart for MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all subroutines required to prepare to save
!>  the restart state for the MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_restart(domain, err)!{{{

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

      err = 0

     if (config_ocn_am_global_stats) then
         call ocn_am_restart_global_stats(domain, err_tmp)
         err = ior(err, err_tmp)
     endif

      if (config_ocn_am_zonal_mean) then
         call ocn_am_restart_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE) then
!         call ocn_am_restart_TEMPLATE(domain, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_am_driver_restart!}}}

!***********************************************************************
!
!  routine ocn_am_driver_write
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

   subroutine ocn_am_driver_write(domain, err)!{{{

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

      err = 0

     if (config_ocn_am_global_stats) then
         call ocn_am_write_global_stats(domain, err_tmp)
         err = ior(err, err_tmp)
     endif

      if (config_ocn_am_zonal_mean) then
         call ocn_am_write_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE) then
!         call ocn_am_write_TEMPLATE(domain, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_am_driver_write!}}}

!***********************************************************************
!
!  routine ocn_am_driver_finalize
!
!> \brief   Finalize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all finalize routines required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine ocn_am_driver_finalize(domain, err)!{{{

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

      err = 0

     if (config_ocn_am_global_stats) then
        call ocn_am_finalize_global_stats(domain, err_tmp)
        err = ior(err, err_tmp)
     endif

      if (config_ocn_am_zonal_mean) then
         call ocn_am_finalize_zonal_mean(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

!      if (config_ocn_am_TEMPLATE) then
!         call ocn_am_finalize_TEMPLATE(domain, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine ocn_am_driver_finalize!}}}

end module ocn_am_driver

! vim: foldmethod=marker
