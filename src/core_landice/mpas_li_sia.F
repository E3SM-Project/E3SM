! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_sia
!
!> \MPAS land-ice SIA velocity driver
!> \author Matt Hoffman
!> \date   16 March 2012
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for calculating velocity using the shallow ice approximation.
!>
!
!-----------------------------------------------------------------------

module li_sia

   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar
   use li_mask

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
   public :: li_sia_init, &
             li_sia_finalize, &
             li_sia_block_init, &
             li_sia_solve

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------



!***********************************************************************

contains

!***********************************************************************
!
!  routine li_sia_init
!
!> \brief   Initializes SIA velocity solver
!> \author  Matt Hoffman/Xylar Asay-Davis
!> \date    16 March 2012
!> \version SVN:$Id$
!> \details
!>  This routine initializes the SIA ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_sia_init(domain, err)

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


      ! No init is needed.
      err = 0

   !--------------------------------------------------------------------

   end subroutine li_sia_init



!***********************************************************************
!
!  routine li_sia_block_init
!
!> \brief   Initializes blocks for SIA velocity solver
!> \author  Matt Hoffman/Xylar Asay-Davis
!> \date    16 March 2012
!> \version SVN:$Id$
!> \details 
!>  This routine initializes each block of the SIA ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_sia_block_init(block, err)

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

      ! No block init needed.
      err = 0

   !--------------------------------------------------------------------
   end subroutine li_sia_block_init



!***********************************************************************
!
!  subroutine li_sia_solve
!
!> \brief   Computes velocity using Shallow Ice Appoximation
!> \author  Matt Hoffman
!> \date    21 May 2012
!> \version SVN:$Id$
!> \details 
!>  This routine computes the normal velocity on edges for each layer 
!>  using the Shallow Ice Approximation.  It calculates ice thickness on 
!>  on an edge using the average of the two neighboring cells (2nd order).
!
!-----------------------------------------------------------------------
   subroutine li_sia_solve(mesh, state, err)
      use mpas_constants, only: gravity

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

      real (kind=RKIND), dimension(:), pointer :: thickness, layerCenterSigma, dcEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:), pointer :: edgeMask
      integer :: nVertLevels, nEdges, iLevel, iEdge, cell1, cell2
      real (kind=RKIND) :: basalVelocity, slopeOnEdge, &
               layerCenterHeightOnEdge, thicknessEdge
      real (kind=RKIND) :: rhoi         ! ice density
      real (kind=RKIND) :: ratefactor   ! flow law parameter, A
      real (kind=RKIND) :: n            ! flow law exponent, n

      err = 0

      ! Set needed variables and pointers
      nEdges = mesh % nEdges
      nVertLevels = mesh % nVertLevels

      dcEdge => mesh % dcEdge % array
      cellsOnEdge => mesh % cellsOnEdge % array
      layerCenterSigma => mesh % layerCenterSigma % array

      normalVelocity => state % normalVelocity % array
      thickness => state % thickness % array
      edgeMask => state % edgeMask % array

      ! Get parameters specified in the namelist
      rhoi = config_ice_density
      n = config_flowLawExponent
      ! Calculate ratefactor (A) at edge - TODO This should be calculated external to this subroutine and as a function of temperature
      ratefactor = config_default_flowParamA  ! units of s^{-1} Pa^{-n}

      basalVelocity = 0.0_RKIND  ! Assume no sliding


      ! Loop over edges
      do iEdge = 1, nEdges
         ! Only calculate velocity for edges that are part of the dynamic ice sheet.(thick ice)
         ! Also, the velocity calculation should be valid for non-ice edges (i.e. returns 0).
         if ( li_mask_is_dynamic_ice(edgeMask(iEdge)) ) then
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)
             ! Calculate slope at edge
             ! This could/should be calculated externally to this subroutine
             slopeOnEdge = (thickness(cell1) - thickness(cell2) ) / dcEdge(iEdge) 
             ! Calculate thickness on edge - 2nd order
             thicknessEdge = (thickness(cell1) + thickness(cell2) ) * 0.5_RKIND
             ! Loop over layers
             do iLevel = 1, nVertLevels
                ! Determine the height of each layer above the bed
                layerCenterHeightOnEdge = thicknessEdge * (1.0_RKIND - layerCenterSigma(iLevel) )
                ! Calculate SIA velocity
                normalVelocity(iLevel,iEdge) = basalVelocity + &
                    0.5_RKIND * ratefactor * (rhoi * gravity)**n * slopeOnEdge**n * &
                    (thicknessEdge**(n+1) - (thicknessEdge - layerCenterHeightOnEdge)**(n+1))
             end do  ! Levels
         endif
      end do  ! edges  
      


     ! === error check
     if (err > 0) then
         write (0,*) "An error has occurred in li_sia_solve."
     endif

   !--------------------------------------------------------------------

   end subroutine li_sia_solve




!***********************************************************************
!
!  routine li_sia_finalize
!
!> \brief   finalizes SIA velocity solver
!> \author  Matt Hoffman/Xylar Asay-Davis
!> \date    16 March 2012
!> \version SVN:$Id$
!> \details
!>  This routine initializes the SIA ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_sia_finalize(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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


   !--------------------------------------------------------------------

   end subroutine li_sia_finalize



   ! private subroutines

 


!***********************************************************************

end module li_sia

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
