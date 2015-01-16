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

      real (kind=RKIND), dimension(:), pointer :: thickness, layerInterfaceSigma, dcEdge, dvEdge
      real (kind=RKIND), dimension(:), pointer :: slopeEdge, normalSlopeEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge
      integer, dimension(:), pointer :: edgeMask
      integer, pointer :: nVertInterfaces, nEdges, nVertices, vertexDegree
      integer :: iLevel, iEdge, iCell, iVertex, cell1, cell2, cellIndex
      real (kind=RKIND) :: basalVelocity, &
               layerInterfaceHeightOnEdge, thicknessEdge, hVertexAccum
      real (kind=RKIND), pointer :: rhoi         ! ice density
      real (kind=RKIND), pointer :: ratefactor   ! flow law parameter, A
      real (kind=RKIND), pointer :: n            ! flow law exponent, n

      err = 0

      ! Set needed variables and pointers
      call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'thickness', thickness, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'edgeMask', edgeMask, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'slopeEdge', slopeEdge, timeLevel=timeLevel)
      call mpas_pool_get_array(statePool, 'normalSlopeEdge', normalSlopeEdge, timeLevel=timeLevel)


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
             ! Calculate thickness on edge - 2nd order
             thicknessEdge = (thickness(cell1) + thickness(cell2) ) * 0.5_RKIND
             ! Loop over layers
             do iLevel = 1, nVertInterfaces
                ! Determine the height of each layer above the bed
                layerInterfaceHeightOnEdge = thicknessEdge * (1.0_RKIND - layerInterfaceSigma(iLevel) )
                ! Calculate SIA velocity
                normalVelocity(iLevel,iEdge) = basalVelocity + &
                    0.5_RKIND * ratefactor * (rhoi * gravity)**n * slopeEdge(iEdge)**(n-1) * normalSlopeEdge(iEdge) * &
                    (thicknessEdge**(n+1) - (thicknessEdge - layerInterfaceHeightOnEdge)**(n+1))
             end do  ! Levels
         else
             normalVelocity(:,iEdge) = 0.0_RKIND
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
