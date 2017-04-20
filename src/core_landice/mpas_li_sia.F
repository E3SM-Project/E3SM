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
!> \details
!>  This module contains the routines for calculating velocity using the shallow ice approximation.
!>
!
!-----------------------------------------------------------------------

module li_sia

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use li_mask
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
!> \details 
!>  This routine computes the normal velocity on edges for each layer 
!>  using the Shallow Ice Approximation.  It calculates ice thickness on 
!>  on an edge using the average of the two neighboring cells (2nd order).
!
!-----------------------------------------------------------------------
   subroutine li_sia_solve(meshPool, statePool, timeLevel, err)
      use mpas_constants, only: gravity

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in) :: &
         timeLevel     !< Input: time level from which to calculate velocity

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

      real (kind=RKIND), dimension(:), pointer :: thickness, layerCenterSigma, dcEdge, dvEdge, upperSurface, upperSurfaceVertex
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge
      integer, dimension(:), pointer :: edgeMask
      integer, pointer :: nVertLevels, nEdges, nVertices, vertexDegree
      integer :: iLevel, iEdge, iCell, iVertex, cell1, cell2, cellIndex
      real (kind=RKIND) :: basalVelocity, slopeOnEdge, &
               normalSlopeOnEdge, tangentSlopeOnEdge,  &
               layerCenterHeightOnEdge, thicknessEdge, hVertexAccum
      real (kind=RKIND), pointer :: rhoi         ! ice density
      real (kind=RKIND), pointer :: ratefactor   ! flow law parameter, A
      real (kind=RKIND), pointer :: n            ! flow law exponent, n

      err = 0

      ! Set needed variables and pointers
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'thickness', thickness, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'edgeMask', edgeMask, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'upperSurface', upperSurface, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'upperSurfaceVertex', upperSurfaceVertex, timeLevel=timeLevel)


      ! Get parameters specified in the namelist
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', n)
      call mpas_pool_get_config(liConfigs, 'config_default_flowParamA', ratefactor)  ! units of s^{-1} Pa^{-n}
      ! Calculate ratefactor (A) at edge - TODO This should be calculated external to this subroutine and as a function of temperature

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
             normalSlopeOnEdge = (upperSurface(cell1) - upperSurface(cell2) ) / dcEdge(iEdge) 
             tangentSlopeOnEdge = ( upperSurfaceVertex(verticesOnEdge(1,iEdge)) - upperSurfaceVertex(verticesOnEdge(2,iEdge)) ) / dvEdge(iEdge)
             slopeOnEdge = (normalSlopeOnEdge**2 + tangentSlopeOnEdge**2)**0.5
             ! Calculate thickness on edge - 2nd order
             thicknessEdge = (thickness(cell1) + thickness(cell2) ) * 0.5_RKIND
             ! Loop over layers
             do iLevel = 1, nVertLevels
                ! Determine the height of each layer above the bed
                layerCenterHeightOnEdge = thicknessEdge * (1.0_RKIND - layerCenterSigma(iLevel) )
                ! Calculate SIA velocity
                normalVelocity(iLevel,iEdge) = basalVelocity + &
                    0.5_RKIND * ratefactor * (rhoi * gravity)**n * slopeOnEdge**(n-1) * normalSlopeOnEdge * &
                    (thicknessEdge**(n+1) - (thicknessEdge - layerCenterHeightOnEdge)**(n+1))
             end do  ! Levels
         else
             normalVelocity(:,iEdge) = 0.0_RKIND
         endif
      end do  ! edges  

     ! === error check
     if (err > 0) then
         call mpas_log_write("An error has occurred in li_sia_solve.", MPAS_LOG_ERR)
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
