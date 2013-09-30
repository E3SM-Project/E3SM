! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_time_integration
!
!> \brief MPAS land ice time integration driver
!> \author Matt Hoffman
!> \date   17 April 2011
!> \version SVN:$Id:$
!> \details
!>  This module contains the main driver routine for calling
!>  time integration schemes
!
!-----------------------------------------------------------------------

module li_time_integration

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use li_time_integration_fe

   implicit none
   private

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

   public :: li_timestep

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************
   contains


!***********************************************************************
!
!  routine li_timestep
!
!> \brief   Advance model state forward in time by the specified time step
!> \author  Matt Hoffman
!> \date    20 April 2012
!> \version SVN:$Id$
!> \details
!>  This routine advances model state forward in time by the specified time step.
!> Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
!>                 plus grid meta-data
!> Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
!>                  model state advanced forward in time by dt seconds
!-----------------------------------------------------------------------
   subroutine li_timestep(domain, dt, timeStamp, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: dt    !< Input:  time step
      character(len=*), intent(in) :: timeStamp  !< Input: current time stamp

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

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
      type (block_type), pointer :: block
      type (state_type), pointer :: state
      integer :: err_tmp

      err = 0
      err_tmp = 0

      !write(*,*) 'Using ', trim(config_time_integration), ' time integration.'
      select case (config_time_integration)
      case ('forward_euler')
         call li_time_integrator_forwardeuler(domain, dt, err_tmp)
      case ('rk4')
         write(0,*) trim(config_time_integration), ' is not currently supported.'
         call mpas_dmpar_abort(domain % dminfo)
         err_tmp = 1
      case default
         write(0,*) trim(config_time_integration), ' is not a valid land ice time integration option.'
         err_tmp = 1
      end select
      err = ior(err,err_tmp)

      block => domain % blocklist
      do while (associated(block))
         ! Assign the time stamp for this time step
         block % state % time_levs(2) % state % xtime % scalar = timeStamp

!         ! Abort the simulation if NaNs occur in the velocity field
!         if (isNaN(sum(block % state % time_levs(2) % state % u % array))) then
!            write(0,*) 'Abort: NaN detected'
!            call mpas_dmpar_abort(dminfo)
!         endif

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in li_timestep."
      endif

   !--------------------------------------------------------------------
   end subroutine li_timestep


end module li_time_integration
