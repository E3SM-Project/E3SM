! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  land_ice_velocity
!
!> \MPAS land-ice velocity driver
!> \author William Lipscomb
!> \date   10 January 2012
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for calling dycores
!>  for calculating ice velocity. 
!>
!
!-----------------------------------------------------------------------


module li_velocity

   use mpas_grid_types
   use mpas_configure
!!!   use land_ice_lifev
   use li_sia

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
!> \version SVN:$Id$
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

      err = 0
      
      write(*,*) 'Using ', trim(config_dycore), ' dynamical core.'
      select case (config_dycore)
      case ('sia')
          call land_ice_sia_init(domain, err)
!!!      case ('L1L2')
!!!          call land_ice_lifev_init(domain, err)
!!!      case ('FO')
!!!          call land_ice_lifev_init(domain, err)
!!!      case ('Stokes')
!!!          call land_ice_lifev_init(domain, err)
!!!          call land_ice_phg_init(domain, err)
      case default
          write(*,*) trim(config_dycore), ' is not a valid land ice dycore option.'
          err = 1
          return
      end select

   !--------------------------------------------------------------------
   end subroutine land_ice_velocity_init



!***********************************************************************
!
!  routine li_velocity_block_init
!
!> \brief   Initializes velocity solver
!> \author  William Lipscomb
!> \date    10 January 2012
!> \version SVN:$Id$
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

      err = 0

      select case (config_dycore)
      case ('sia')
          call li_sia_block_init(block, err)
!!!      case ('L1L2')
!!!          call land_ice_lifev_block_init(block, err)
!!!      case ('FO')
!!!          call land_ice_lifev_block_init(block, err)
!!!      case ('Stokes')
!!!          call land_ice_lifev_block_init(block, err)
      case default
          write(*,*) trim(config_dycore), ' is not a valid land ice dycore option.'
          err = 1
          return
      end select


   !--------------------------------------------------------------------

   end subroutine land_ice_velocity_block_init



!***********************************************************************
!
!  routine li_velocity_solve
!
!> \brief   Solver for calculating ice velocity 
!> \author  William Lipscomb
!> \date    10 January 2012
!> \version SVN:$Id$
!> \details 
!>  This routine calls velocity solvers.
!
!-----------------------------------------------------------------------
   subroutine land_ice_velocity_solve(mesh, state, err)

      use li_sia

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         state          !< Input: state information 

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
      integer :: iEdge, nEdges
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      integer, dimension(:), pointer :: edgeMask

      err = 0

      nEdges = mesh % nEdges
      normalVelocity => state % normalVelocity % array
      edgeMask => state % edgeMask % array

      select case (config_dycore)
      case ('SIA')
          call li_sia_solve(mesh, state, err)
!!!      case ('L1L2')
!!!          call land_ice_lifev_solve(mesh, state, err)
!!!      case ('FO')
!!!          call land_ice_lifev_solve(mesh, state, err)
!!!      case ('Stokes')
!!!          call land_ice_lifev_solve(mesh, state, err)
      case default
          write(*,*) trim(config_dycore), ' is not a valid land ice dycore option.'
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

   !--------------------------------------------------------------------
   end subroutine land_ice_velocity_solve



!***********************************************************************
!
!  routine li_velocity_finalize
!
!> \brief   Finalizes velocity solver
!> \author  Xylar Asay-Davis
!> \date    18 January 2012
!> \version SVN:$Id$
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

      err = 0

      select case (config_dycore)
      case ('SIA')
          call land_ice_sia_finalize(domain, err)
!!!      case ('L1L2')
!!!          call land_ice_lifev_finalize(domain, err)
!!!      case ('FO')
!!!          call land_ice_lifev_finalize(domain, err)
!!!      case ('Stokes')
!!!          call land_ice_lifev_finalize(domain, err)
      case default
          write(*,*) trim(config_dycore), ' is not a valid land ice dycore option.'
          err = 1
          return
      end select

   !--------------------------------------------------------------------
   end subroutine land_ice_velocity_finalize



!***********************************************************************

end module li_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
