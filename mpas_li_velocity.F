! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_velocity
!
!> \MPAS land-ice velocity driver
!> \author William Lipscomb
!> \date   10 January 2012
!> \details
!>  This module contains the routines for calling dycores
!>  for calculating ice velocity. 
!>
!
!-----------------------------------------------------------------------


module li_velocity

   use mpas_grid_types
   use mpas_configure
   use li_velocity_external
   use li_sia
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

   public :: li_velocity_init, &
             li_velocity_finalize, &
             li_velocity_block_init, &
             li_velocity_solve

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine li_velocity_init
!
!> \brief   Initializes velocity solver
!> \author  Xylar Asay-Davis
!> \date    18 January 2012
!> \details
!>  This routine initializes the ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_init(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

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
      character (len=StrKIND), pointer :: config_velocity_solver
      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      write(*,*) 'Using ', trim(config_velocity_solver), ' dynamical core.'
      select case (config_velocity_solver)
      case ('sia')
          call li_sia_init(domain, err)
      case ('L1L2', 'FO', 'Stokes')
          call li_velocity_external_init(domain, err)
      case default
          write(0,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
      end select

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in li_velocity_init."
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_init



!***********************************************************************
!
!  routine li_velocity_block_init
!
!> \brief   Initializes velocity solver
!> \author  William Lipscomb
!> \date    10 January 2012
!> \details 
!>  This routine initializes each block of the ice velocity solver.
!
!-----------------------------------------------------------------------
   subroutine li_velocity_block_init(block, err)

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

      type (block_type), intent(inout) :: &
         block          !< Input/Output: block object 

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
      character (len=StrKIND), pointer :: config_velocity_solver

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      select case (config_velocity_solver)
      case ('sia')
          call li_sia_block_init(block, err)
      case ('L1L2', 'FO', 'Stokes')
          call li_velocity_external_block_init(block, err)
      case default
          write(*,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
          return
      end select

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in li_velocity_block_init."
      endif

   !--------------------------------------------------------------------

   end subroutine li_velocity_block_init



!***********************************************************************
!
!  routine li_velocity_solve
!
!> \brief   Solver for calculating ice velocity 
!> \author  William Lipscomb
!> \date    10 January 2012
!> \details 
!>  This routine calls velocity solvers.
!
!-----------------------------------------------------------------------
   subroutine li_velocity_solve(meshPool, statePool, timeLevel, err)

      use li_mask

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in) :: timeLevel  !< Input: Time level on which to calculate diagnostic variables

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool          !< Input: state information 

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
      ! pointers to get from pools
      character (len=StrKIND), pointer :: config_velocity_solver
      integer, pointer :: nEdges
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      integer, dimension(:), pointer :: edgeMask
      ! truly local variables
      integer :: iEdge

      err = 0

      ! Get variables from pools
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'edgeMask', edgeMask, timeLevel=timeLevel)


      select case (config_velocity_solver)
      case ('sia')
          call li_sia_solve(meshPool, statePool, timeLevel, err)
      case ('L1L2', 'FO', 'Stokes')
          call li_velocity_external_solve(meshPool, statePool, timeLevel, err)
      case default
          write(0,*) 'Error: ', trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
          return
      end select

      ! Check if the velocity solver has returned a velocity on any non-dynamic edges
      do iEdge = 1, nEdges
         if ( li_mask_is_ice(edgeMask(iEdge)) .and.                    &
              (.not.(li_mask_is_dynamic_ice(edgeMask(iEdge)))) .and.   &
              (maxval(abs(normalVelocity(:,iEdge))) /= 0.0_RKIND)      &
            ) then
            err = 1
            !!!normalVelocity(:,iEdge) = 0.0_RKIND  ! this is a hack because the rest of the code requires this, but this condition should really cause a fatal error.
         endif
      enddo
      if (err == 1) then
         write(0,*) 'Error: Velocity has been calculated on non-dynamic edges.  There is a problem with the velocity solver.' !!!  Velocity on those edges have been set to 0, but this should be a fatal error.'
         err = 1  ! a hack to let the code continue until this can be fixed in the velocity solver
      end if

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in li_velocity_solve."
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_solve



!***********************************************************************
!
!  routine li_velocity_finalize
!
!> \brief   Finalizes velocity solver
!> \author  Xylar Asay-Davis
!> \date    18 January 2012
!> \details
!>  This routine finalizes the ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_finalize(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain   !< Input/Output: domain object

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
      character (len=StrKIND), pointer :: config_velocity_solver

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      select case (config_velocity_solver)
      case ('sia')
          call li_sia_finalize(domain, err)
      case ('L1L2', 'FO', 'Stokes')
          call li_velocity_external_finalize(err)
      case default
          write(*,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
          return
      end select

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in li_velocity_finalize."
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_finalize



!***********************************************************************

end module li_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
