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

   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar
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

      use mpas_geometry_utils, only: mpas_calculate_barycentric_weights_for_points

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
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      integer :: iCell, iLevel, i, iVertex, err_tmp
      integer, pointer :: nVertices
      character (len=StrKIND), pointer :: config_sia_tangent_slope_calculation
      integer, dimension(:,:), pointer :: baryCellsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: baryWeightsOnVertex
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      type (field1dInteger), pointer :: vertexIndicesField

      ! No block init needed.
      err = 0
      err_tmp = 0

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_array(meshPool, 'baryCellsOnVertex', baryCellsOnVertex)
      call mpas_pool_get_array(meshPool, 'baryWeightsOnVertex', baryWeightsOnVertex)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_config(liConfigs, 'config_sia_tangent_slope_calculation', config_sia_tangent_slope_calculation)
      call mpas_pool_get_field(scratchPool, 'vertexIndices', vertexIndicesField)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)

      ! The SIA solver may need to setup these weights for calculating upperSurfaceVertex
      if (trim(config_sia_tangent_slope_calculation) == 'from_vertex_barycentric') then
         call mpas_allocate_scratch_field(vertexIndicesField, .true.)
         do iVertex = 1, nVertices
            vertexIndicesField % array(iVertex) = iVertex
         enddo
         call mpas_calculate_barycentric_weights_for_points(meshPool, &
                xVertex(1:nVertices), yVertex(1:nVertices), zVertex(1:nVertices), &
                vertexIndicesField % array(1:nVertices), &
                baryCellsOnVertex(:, 1:nVertices), baryWeightsOnVertex(:, 1:nVertices), err_tmp)
         err = ior(err, err_tmp)
         call mpas_deallocate_scratch_field(vertexIndicesField, .true.)
      endif


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

      real (kind=RKIND), dimension(:), pointer :: thickness, layerInterfaceSigma
      real (kind=RKIND), dimension(:), pointer :: slopeEdge, normalSlopeEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, flowParamA
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:), pointer :: edgeMask, cellMask
      integer, pointer :: nVertInterfaces, nEdges
      integer :: iLevel, iEdge
      integer :: cell1, cell2
      real (kind=RKIND) :: thicknessEdge, flwaLevelEdge
      real (kind=RKIND) :: positionIndependentFactor ! The portion of the velocity calculation that is completely independent of position
      real (kind=RKIND) :: levelIndependentFactor ! The portion of the velocity calculation that depends on horizontal location but not vertical position
      real (kind=RKIND), pointer :: rhoi         ! ice density
      real (kind=RKIND), pointer :: n            ! flow law exponent, n
      integer :: cell1_is_dynamic, cell2_is_dynamic

      err = 0

      ! Set needed variables and pointers
      call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'flowParamA', flowParamA, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'thickness', thickness, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'edgeMask', edgeMask, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'cellMask', cellMask, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'slopeEdge', slopeEdge, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'normalSlopeEdge', normalSlopeEdge, timeLevel=timeLevel)


      ! Get parameters specified in the namelist
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', n)

      positionIndependentFactor = -0.5_RKIND * (rhoi * gravity)**n  ! could be calculated once on init

      ! Loop over edges
      do iEdge = 1, nEdges

         ! Only calculate velocity for edges that are part of the dynamic ice sheet.(thick ice)
         ! Also, the velocity calculation should be valid for non-ice edges (i.e. returns 0).
         if ( li_mask_is_dynamic_ice(edgeMask(iEdge)) ) then
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)
             cell1_is_dynamic = li_mask_is_dynamic_ice_int(cellMask(cell1))
             cell2_is_dynamic = li_mask_is_dynamic_ice_int(cellMask(cell2))

             ! Calculate thickness on edge - 2nd order
             thicknessEdge = (thickness(cell1) + thickness(cell2) ) * 0.5_RKIND
!             thicknessEdge = (thickness(cell1) * cell1_is_dynamic + thickness(cell2) * cell2_is_dynamic) / (cell1_is_dynamic + cell2_is_dynamic)  ! this version does an upwind thickness on margin edges only.  Most Halfar error stats are higher by about 10-20%
             ! Also tried upwind everywhere [for dome can be hacked with: thicknessEdge = max(thickness(cell1), thickness(cell2) ]  This results in Halfar errors that are about 5x larger than centered difference

             levelIndependentFactor = slopeEdge(iEdge)**(n-1) * normalSlopeEdge(iEdge) * thicknessEdge**(n+1)

             normalVelocity(nVertInterfaces, iEdge) = 0.0_RKIND  ! Assume no sliding

             do iLevel = nVertInterfaces-1, 1, -1  ! Loop upwards from second lowest level to surface
                ! Calculate flwa on edge for this level - 2nd order, except can't do centered difference into areas where flwa may not be valid, so excluding the downwind flwa value in non-dynamic cells
!                flwaLevelEdge = (flowParamA(iLevel, cell1) + flowParamA(iLevel, cell2) ) * 0.5_RKIND
                flwaLevelEdge = (flowParamA(iLevel, cell1) * cell1_is_dynamic +  &
                   flowParamA(iLevel, cell2) * cell2_is_dynamic) / &
                   (cell1_is_dynamic + cell2_is_dynamic)

                ! Calculate SIA velocity for this layer interface by adding on incremental velocity for the layer below
                ! (This requires that flwa be constant over that layer, which it is.)
                normalVelocity(iLevel, iEdge) = normalVelocity(iLevel+1, iEdge) +  &
                   positionIndependentFactor * levelIndependentFactor * flwaLevelEdge *  &
                   ( (layerInterfaceSigma(iLevel))**(n+1) - (layerInterfaceSigma(iLevel+1))**(n+1) )
             end do
         else
             normalVelocity(:,iEdge) = 0.0_RKIND  ! zero velocity on non-dynamic edges
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
