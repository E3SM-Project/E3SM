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
!> \details
!>  This module contains the main driver routine for calling
!>  time integration schemes
!
!-----------------------------------------------------------------------

module li_time_integration

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_log
   use li_time_integration_fe
   use li_setup

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
!> \details
!>  This routine advances model state forward in time by the specified time step.
!> Input: domain - current model state in time level 1 plus mesh data
!> Output: domain - upon exit, time level 2 contains
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
      type (mpas_pool_type), pointer :: statePool
      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: config_time_integration
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_time_integration', config_time_integration)

      !write(*,*) 'Using ', trim(config_time_integration), ' time integration.'
      select case (config_time_integration)
      case ('forward_euler')
         call li_time_integrator_forwardeuler(domain, dt, err_tmp)
      case ('rk4')
         call mpas_log_write(trim(config_time_integration) // ' is not currently supported.', MPAS_LOG_ERR)
         call mpas_dmpar_abort(domain % dminfo)
         err_tmp = 1
      case default
         call mpas_log_write(trim(config_time_integration) // ' is not a valid land ice time integration option.', MPAS_LOG_ERR)
         err_tmp = 1
      end select
      err = ior(err,err_tmp)

      block => domain % blocklist
      do while (associated(block))
         ! Assign the time stamp for this time step
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_array(statePool, 'xtime', xtime, timeLevel=1)  ! xtime only has one time level, but stating is explicitly here to avoid confusion later.
         xtime = timeStamp

!         ! Abort the simulation if NaNs occur in the velocity field
!         if (isNaN(sum(block % state % time_levs(2) % state % u % array))) then
!            write(0,*) 'Abort: NaN detected'
!            call mpas_dmpar_abort(dminfo)
!         endif

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_timestep.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
   end subroutine li_timestep


end module li_time_integration
