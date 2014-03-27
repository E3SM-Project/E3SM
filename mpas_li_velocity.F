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
!!!   use li_lifev
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
!!!      case ('L1L2')
!!!          call li_lifev_init(domain, err)
!!!      case ('FO')
!!!          call li_lifev_init(domain, err)
!!!      case ('Stokes')
!!!          call li_lifev_init(domain, err)
!!!          call li_phg_init(domain, err)
      case default
          write(*,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
          return
      end select

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in li_timestep."
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
!!!      case ('L1L2')
!!!          call li_lifev_block_init(block, err)
!!!      case ('FO')
!!!          call li_lifev_block_init(block, err)
!!!      case ('Stokes')
!!!          call li_lifev_block_init(block, err)
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

      use li_sia

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
      character (len=StrKIND), pointer :: config_velocity_solver
!      integer :: iEdge, nEdges
!      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
!      integer, dimension(:), pointer :: edgeMask

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

!      nEdges = mesh % nEdges
!      normalVelocity => state % normalVelocity % array
!      edgeMask => state % edgeMask % array

      select case (config_velocity_solver)
      case ('sia')
          call li_sia_solve(meshPool, statePool, timeLevel, err)
!!!      case ('L1L2')
!!!          call li_lifev_solve(mesh, state, timeLevel, err)
!!!      case ('FO')
!!!          call li_lifev_solve(mesh, state, timeLevel, err)
!!!      case ('Stokes')
!!!          call li_lifev_solve(mesh, state, timeLevel, err)
      case default
          write(*,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
          return
      end select


!!!      do iEdge = 1, nEdges
!!!         if ( MASK_IS_THIN_ICE(edgeMask(iEdge)) .and. (maxval(abs(normalVelocity(:,iEdge))) /= 0.0_RKIND) ) then
!!!            err = 1
!!!            normalVelocity(:,iEdge) = 0.0_RKIND  ! this is a hack because the rest of the code requires this, but this condition should really cause a fatal error.
!!!         endif
!!!      enddo
!!!      if (err == 1) then
!!!         write(0,*) 'Velocity has been calculated on non-dynamic edges.  There is a problem with the velocity solver.  Velocity on those edges have been set to 0, but this should be a fatal error.'
!!!         err = 0  ! a hack to let the code continue until this can be fixed in the velocity solver
!!!      end if

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
!!!      case ('L1L2')
!!!          call li_lifev_finalize(domain, err)
!!!      case ('FO')
!!!          call li_lifev_finalize(domain, err)
!!!      case ('Stokes')
!!!          call li_lifev_finalize(domain, err)
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
