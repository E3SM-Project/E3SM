! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  oac_driver
!
!> \brief Driver for MPAS ocean analysis core
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This is the driver for the MPAS ocean core.
!
!-----------------------------------------------------------------------

module oac_driver

   use mpas_grid_types
   use mpas_configure

   use oac_global_stats
!   use oac_TEMPLATE

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

   public :: oac_driver_init, &
             oac_driver_compute, &
             oac_driver_write, &
             oac_driver_restart, &
             oac_driver_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine oac_driver_init
!
!> \brief   Initialize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all initializations required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine oac_driver_init(err)!{{{

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

     if (config_oac_global_stats) then
         call oac_init_global_stats(err_tmp)
         err = ior(err, err_tmp)
     endif

!      if (config_oac_TEMPLATE) then
!         call oac_init_TEMPLATE(err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine oac_driver_init!}}}

!***********************************************************************
!
!  routine oac_driver_compute
!
!> \brief   Driver for MPAS-Ocean analysis computations
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all computation subroutines required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine oac_driver_compute(domain, err)!{{{

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
      state => domain % blocklist % state % time_levs(timeLevel) % state

     if (config_oac_global_stats) then
         call oac_compute_global_stats(domain % blocklist % amGlobalStats, state, domain % dminfo, timeLevel, err_tmp)
         err = ior(err, err_tmp)
     endif

!      if (config_oac_TEMPLATE) then
!         call oac_compute_TEMPLATE(domain % blocklist % amTemplate, state, domain % dminfo, timeLevel, err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine oac_driver_compute!}}}

!***********************************************************************
!
!  routine oac_driver_restart
!
!> \brief   Save restart for MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all subroutines required to prepare to save
!>  the restart state for the MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine oac_driver_restart(err)!{{{

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

     if (config_oac_global_stats) then
         call oac_restart_global_stats(err_tmp)
         err = ior(err, err_tmp)
     endif

!      if (config_oac_TEMPLATE) then
!         call oac_restart_TEMPLATE(err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine oac_driver_restart!}}}

!***********************************************************************
!
!  routine oac_driver_write
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

   subroutine oac_driver_write(err)!{{{

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

     if (config_oac_global_stats) then
         call oac_write_global_stats(err_tmp)
         err = ior(err, err_tmp)
     endif

!      if (config_oac_TEMPLATE) then
!         call oac_write_TEMPLATE(err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine oac_driver_write!}}}

!***********************************************************************
!
!  routine oac_driver_finalize
!
!> \brief   Finalize MPAS-Ocean analysis driver
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all finalize routines required for the 
!>  MPAS-Ocean analysis driver.
!
!-----------------------------------------------------------------------

   subroutine oac_driver_finalize(err)!{{{

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

     if (config_oac_global_stats) then
        call oac_finalize_global_stats(err_tmp)
         err = ior(err, err_tmp)
     endif

!      if (config_oac_TEMPLATE) then
!         call oac_finalize_TEMPLATE(err_tmp)
!         err = ior(err, err_tmp)
!      endif

   end subroutine oac_driver_finalize!}}}

end module oac_driver

! vim: foldmethod=marker
