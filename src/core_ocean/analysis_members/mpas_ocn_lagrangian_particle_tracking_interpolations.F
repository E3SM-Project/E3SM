! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  ocn_lagrangian_particle_tracking_interpolations
!
!> \brief   LIGHT Vector reconstruction and filtering module
!> \author  Phillip J. Wolfram
!> \date    07/21/2015
!> \details
!> This module provides routines for performing vector interpolations
!> and spatial filtering.
!
!-----------------------------------------------------------------------
module ocn_lagrangian_particle_tracking_interpolations

  use mpas_derived_types
  use mpas_constants
  use mpas_rbf_interpolation
  use mpas_geometry_utils
  use mpas_vector_reconstruction
  use mpas_dmpar

  implicit none

  contains

!***********************************************************************
!
!  routine ocn_vertex_reconstruction
!
!> \brief   Reconstruct vertex velocity driver / interface
!> \author  Phillip Wolfram
!> \date    03/27/2014
!> \details
!>  Purpose: reconstruct vector field at vertex locations based on
!>           particular choice of reconstruction method
!>  Input: mesh meta data and vector component data residing at cell edges
!>         initialize_weights logical is to determine if weights should be initialized
!>  Output: reconstructed vector field (measured in X,Y,Z) located at vertices
!-----------------------------------------------------------------------
  subroutine ocn_vertex_reconstruction(filterNum, meshPool, scratchPool, particleCellPool, layerThickness, u, & !{{{
                                       uvReconstructX, uvReconstructY, uvReconstructZ )

    implicit none

    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    type (mpas_pool_type), pointer, intent(in) :: scratchPool !< Input: Scratch variables
    type (mpas_pool_type), pointer, intent(in) :: particleCellPool !< Input: particlefield variables
    integer, intent(in) :: filterNum  ! filtering strength employed
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness !< Input: layerThickness on cells
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: u !< Input: Velocity field on edges (normalVelocity)
    type (field2DReal), pointer, intent(inout) :: uvReconstructX !< Output: X Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(inout) :: uvReconstructY !< Output: Y Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(inout) :: uvReconstructZ !< Output: Z Component of velocity reconstructed to vertices

    ! could add additional reconstruction techniques here with switch if desired

    ! assumption is made that mpas_init_reconstruct was previously called
    call ocn_RBFvertex(meshPool, filterNum, layerThickness, u, uvReconstructX, uvReconstructY, uvReconstructZ, .false., &
                       scratchPool, particleCellPool)

  end subroutine ocn_vertex_reconstruction!}}}

!***********************************************************************
!
!  routine ocn_RBFvertex
!
!> \brief   Reconstruct vertex velocity using linear interpolation of
!>          RBFs reconstruction at cell centers
!> \author  Phillip Wolfram, Todd Ringler
!> \date    03/26/2014
!> \details
!>  Purpose: reconstruct vector field at vertex locations based on radial basis functions
!>  Input: mesh meta data and vector component data residing at cell edges
!>         initialize_weights logical is to determine if weights should be initialized
!>  Output: reconstructed vector field (measured in X,Y,Z) located at vertices
!-----------------------------------------------------------------------
  subroutine ocn_RBFvertex(meshPool, filterNum, layerThickness, u, uvReconstructX, uvReconstructY, uvReconstructZ, & !{{{
                           initialize_weights, scratchPool, particleCellPool)

    implicit none

    ! inputs
    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    type (mpas_pool_type), pointer, intent(in) :: scratchPool
    type (mpas_pool_type), pointer, intent(in) :: particleCellPool !< Input: particlefield variables
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: u !< Input: Velocity field on edges
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness !< Input: layerThickness on cells
    integer, intent(in) :: filterNum  !< number of times to filter
    logical, intent(in) :: initialize_weights !< Input: Determine if weights for RBF should be pre-computed

    ! outputs
    type (field2DReal), pointer, intent(inout) :: uvReconstructX !< Output: X Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(inout) :: uvReconstructY !< Output: Y Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(inout) :: uvReconstructZ !< Output: Z Component of velocity reconstructed to vertices

    ! local / temporary arrays needed in the compute procedure
    type (field2DReal), pointer :: &
      ucReconstructX, ucReconstructY, ucReconstructZ, ucReconstructZonal, ucReconstructMeridional ! cell center values
    type (field2DReal), pointer :: ucStore, vcStore, wcStore
    type (field2DInteger), pointer :: boundaryVertex, boundaryCell, boundaryCellGlobal, boundaryVertexGlobal

    ! get pointers
    call mpas_pool_get_field(scratchPool, 'ucReconstructX', ucReconstructX)
    call mpas_pool_get_field(scratchPool, 'ucReconstructY', ucReconstructY)
    call mpas_pool_get_field(scratchPool, 'ucReconstructZ', ucReconstructZ)
    call mpas_pool_get_field(scratchPool, 'ucReconstructZonal', ucReconstructZonal)
    call mpas_pool_get_field(scratchPool, 'ucReconstructMeridional', ucReconstructMeridional)
    call mpas_pool_get_field(scratchPool, 'boundaryVertexGlobal', boundaryVertexGlobal)

    ! allocate memory
    call mpas_allocate_scratch_field(ucReconstructX, .True.)
    call mpas_allocate_scratch_field(ucReconstructY, .True.)
    call mpas_allocate_scratch_field(ucReconstructZ, .True.)
    call mpas_allocate_scratch_field(ucReconstructZonal, .True.)
    call mpas_allocate_scratch_field(ucReconstructMeridional, .True.)
    call mpas_allocate_scratch_field(boundaryVertexGlobal, .True.)

    ucReconstructX % array = 0.0_RKIND
    ucReconstructY % array = 0.0_RKIND
    ucReconstructZ % array = 0.0_RKIND

    ! initialize weights (should be pre-initialized)
    if (initialize_weights) then
      call mpas_init_reconstruct(meshPool)
    end if

    ! get cell center reconstructed RBF values
    call mpas_reconstruct(meshPool, u, ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array, &
      ucReconstructZonal % array, ucReconstructMeridional % array)

    ! need to do exchange for uc components (we don't use Zonal / Meridional for this calculation)
    call mpas_dmpar_exch_halo_field(ucReconstructX)
    call mpas_dmpar_exch_halo_field(ucReconstructY)
    call mpas_dmpar_exch_halo_field(ucReconstructZ)

    ! get boundaries
    call mpas_pool_get_field(meshPool,'boundaryVertex', boundaryVertex)
    call mpas_pool_get_field(meshPool,'boundaryCell', boundaryCell)


    if (filternum > 0) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! filter the cell velocity field
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call ocn_second_order_shapiro_filter_ops(filterNum, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! store filter data &
      ! write data to file for output
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_pool_get_field(particleCellPool, 'filteredVelocityU', ucStore)
      call mpas_pool_get_field(particleCellPool, 'filteredVelocityV', vcStore)
      call mpas_pool_get_field(particleCellPool, 'filteredVelocityW', wcStore)
      ucStore % array = ucReconstructX % array
      vcStore % array = ucReconstructY % array
      wcStore % array = ucReconstructZ % array

    end if

    ! interpolate to vertex locations for use in Wachspress
    call ocn_vector_cell_center_to_vertex(meshPool, boundaryVertex % array, boundaryCell % array, &
      ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array, &
      uvReconstructX % array, uvReconstructY % array, uvReconstructZ % array)

    ! handle boundary vertices (should be zero).  Can potentially remove if mpas_init_block initializes to 0 vs -1e34
    boundaryVertexGlobal % array = boundaryVertex % array
    call mpas_dmpar_exch_halo_field(boundaryVertexGlobal)
    ! definite change between these fields!
    uvReconstructX % array = uvReconstructX % array * (1.0_RKIND - boundaryVertexGlobal % array)
    uvReconstructY % array = uvReconstructY % array * (1.0_RKIND - boundaryVertexGlobal % array)
    uvReconstructZ % array = uvReconstructZ % array * (1.0_RKIND - boundaryVertexGlobal % array)

    ! do halo exchanges
    call mpas_dmpar_exch_halo_field(uvReconstructX)
    call mpas_dmpar_exch_halo_field(uvReconstructY)
    call mpas_dmpar_exch_halo_field(uvReconstructZ)

    ! deallocate memory
    call mpas_deallocate_scratch_field(ucReconstructX, .True.)
    call mpas_deallocate_scratch_field(ucReconstructY, .True.)
    call mpas_deallocate_scratch_field(ucReconstructZ, .True.)
    call mpas_deallocate_scratch_field(ucReconstructZonal, .True.)
    call mpas_deallocate_scratch_field(ucReconstructMeridional, .True.)
    call mpas_deallocate_scratch_field(boundaryVertexGlobal, .True.)

  end subroutine ocn_RBFvertex!}}}

!***********************************************************************
!
!  routine ocn_vector_cell_center_to_vertex
!
!> \brief   Interpolate cell center values to vertex values
!> \author  Phillip Wolfram
!> \date    05/27/2014
!> \details
!>  Purpose: interpolate vector field at vertex locations from cell center values
!>           using Barycentric (via Wachspress) interpolation
!>  Input: cell center data and mesh information
!>  Output: interpolated vertex values
!-----------------------------------------------------------------------
  subroutine ocn_vector_cell_center_to_vertex(meshPool, boundaryVertex, boundaryCell, & !{{{
      ucReconstructX, ucReconstructY, ucReconstructZ, &
      uvReconstructX, uvReconstructY, uvReconstructZ)

    implicit none

    ! input variables
    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: ucReconstructX, & !< Input: X Cell center values
                                                              ucReconstructY, & !< Input: Y Cell center values
                                                              ucReconstructZ    !< Input: z Cell center values
    integer, dimension(:,:), pointer, intent(in) :: boundaryVertex, boundaryCell !< Input: Boundary flags

    ! output variables
    real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: uvReconstructX !< Output: Vertex Reconstructed X Velocity Component
    real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: uvReconstructY !< Output: Vertex Reconstructed Y Velocity Component
    real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: uvReconstructZ !< Output: Vertex Reconstructed Z Velocity Component

    ! local variables
    integer, pointer :: nVerticesSolve, nCells, vertexDegree, nVertLevels
    integer :: aVertex, aCell, aLevel
    integer, dimension(:,:), pointer :: cellsOnVertex
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
    real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
    ! temporary arrays needed in the (to be constructed) init procedure
    ! note that lambda is going to be constant for this and could be cached
    real (kind=RKIND), dimension(:), allocatable :: lambda
    real (kind=RKIND), dimension(:,:), allocatable :: pointVertex
    real (kind=RKIND), dimension(3) :: pointInterp
    real (kind=RKIND) :: xp,yp,zp , sumArea, kiteArea
    logical, pointer :: is_periodic
    real(kind=RKIND), pointer :: x_period, y_period

    uvReconstructX = 0.0_RKIND
    uvReconstructY = 0.0_RKIND
    uvReconstructZ = 0.0_RKIND

    call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

    allocate(lambda(vertexDegree), pointVertex(3,vertexDegree))

    ! setup pointers
    call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
    call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'zCell', zCell)

    call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
    call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
    call mpas_pool_get_array(meshPool, 'zVertex', zVertex)

    call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)

    call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
    call mpas_pool_get_config(meshPool, 'x_period', x_period)
    call mpas_pool_get_config(meshPool, 'y_period', y_period)

    ! loop over all vertices
    do aVertex = 1, nVerticesSolve
      ! could precompute the list as an optimization
      ! really, condition is any boundaryVertex in column greater than 0
      if(any(boundaryVertex(:,aVertex) < 1)) then
        ! get vertex location and cell center locations
        do aCell = 1, vertexDegree
          ! logical could be moved outside of code block as an optimization
          ! (then essentially would have two nearly identical code blocks...)
          if (is_periodic) then
            ! fix periodicity with respect to pointInterp (xVertex)
            pointVertex(1,aCell) = mpas_fix_periodicity(xCell(cellsOnVertex(aCell, aVertex)), xVertex(aVertex), x_period)
            pointVertex(2,aCell) = mpas_fix_periodicity(yCell(cellsOnVertex(aCell, aVertex)), yVertex(aVertex), y_period)
            pointVertex(3,aCell) = zCell(cellsOnVertex(aCell, aVertex))
          else
            pointVertex(1,aCell) = xCell(cellsOnVertex(aCell, aVertex))
            pointVertex(2,aCell) = yCell(cellsOnVertex(aCell, aVertex))
            pointVertex(3,aCell) = zCell(cellsOnVertex(aCell, aVertex))
          end if
        end do
        ! vertex point for reconstruction
        pointInterp(1) = xVertex(aVertex)
        pointInterp(2) = yVertex(aVertex)
        pointInterp(3) = zVertex(aVertex)
        ! get interpolation constants (could be cached / optimized with areaBin)
        lambda = mpas_wachspress_coordinates(vertexDegree, pointVertex , pointInterp, meshPool)
      else
        lambda = 0.0_RKIND
      end if

      do aLevel = 1, nVertLevels
        if(boundaryVertex(aLevel,aVertex) < 1) then
          ! perform interpolation
          uvReconstructX(aLevel,aVertex) = sum(ucReconstructX(aLevel,cellsOnVertex(:,aVertex)) * lambda)
          uvReconstructY(aLevel,aVertex) = sum(ucReconstructY(aLevel,cellsOnVertex(:,aVertex)) * lambda)
          uvReconstructZ(aLevel,aVertex) = sum(ucReconstructZ(aLevel,cellsOnVertex(:,aVertex)) * lambda)
        end if
      end do

      ! need to specify boundary conditions for the vertexes (outside this subroutine)

    end do

    deallocate(lambda, pointVertex)

  end subroutine ocn_vector_cell_center_to_vertex!}}}

!***********************************************************************
!
!  routine ocn_vector_vertex_to_cell_center
!
!> \brief   Interpolate vertex values to cell center
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: interpolate vector field at cell center locations from vertex values
!>           using Wachspress interpolation
!>  Input: vertex vector data and mesh information
!>  Output: interpolated cell values
!-----------------------------------------------------------------------
  subroutine ocn_vector_vertex_to_cell_center(meshPool, & !{{{
      uvReconstructX, uvReconstructY, uvReconstructZ, &
      ucReconstructX, ucReconstructY, ucReconstructZ)

    implicit none

    ! input variables
    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: uvReconstructX, & !< Input: Vertex x values
                                                              uvReconstructY, & !< Input: Vertex y values
                                                              uvReconstructZ    !< Input: Vertex z values

    ! output variables
    real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: ucReconstructX !< Output: X Component of velocity
                                                                              !<         reconstructed to cells
    real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: ucReconstructY !< Output: Y Component of velocity
                                                                              !<         reconstructed to cells
    real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: ucReconstructZ !< Output: Z Component of velocity
                                                                              !<         reconstructed to cells

    ! local variables
    integer, pointer :: nCellsSolve, nVertLevels
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: aVertex, aCell, aLevel, nLocalVertices
    integer, dimension(:,:), pointer :: verticesOnCell
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
    ! temporary arrays needed in the (to be constructed) init procedure
    ! note that lambda is going to be constant for this and could be cached
    real (kind=RKIND), dimension(:), allocatable :: lambda
    real (kind=RKIND), dimension(3) :: pointInterp
    real (kind=RKIND), dimension(:,:), allocatable :: pointVertex
    real (kind=RKIND) :: xp,yp,zp

    ucReconstructX = 0.0_RKIND
    ucReconstructY = 0.0_RKIND
    ucReconstructZ = 0.0_RKIND

    ! setup pointers
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'zCell', zCell)

    call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
    call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
    call mpas_pool_get_array(meshPool, 'zVertex', zVertex)


    ! loop over all vertices
    do aCell = 1, nCellsSolve
      ! could precompute the list as an optimization to
      ! remove the following lines !{{{
      nLocalVertices = nEdgesOnCell(aCell)
      ! really, condition is any boundaryVertex in column greater than 0
      allocate(lambda(nLocalVertices), pointVertex(3,nLocalVertices))
      ! get vertex location and cell center locations
      do aVertex = 1, nLocalVertices
        pointVertex(1,aVertex) = xVertex(verticesOnCell(aVertex, aCell))
        pointVertex(2,aVertex) = yVertex(verticesOnCell(aVertex, aCell))
        pointVertex(3,aVertex) = zVertex(verticesOnCell(aVertex, aCell))
      end do
      ! vertex point for reconstruction
      pointInterp(1) = xCell(aCell)
      pointInterp(2) = yCell(aCell)
      pointInterp(3) = zCell(aCell)
      ! get interpolation constants (should be cached as an optimization!)
      lambda = mpas_wachspress_coordinates(nLocalVertices, pointVertex , pointInterp, meshPool)
      !}}}

      do aLevel = 1, nVertLevels
          ! perform interpolation
          ucReconstructX(aLevel,aCell) = sum(uvReconstructX(aLevel,verticesOnCell(1:nLocalVertices,aCell)) * lambda)
          ucReconstructY(aLevel,aCell) = sum(uvReconstructY(aLevel,verticesOnCell(1:nLocalVertices,aCell)) * lambda)
          ucReconstructZ(aLevel,aCell) = sum(uvReconstructZ(aLevel,verticesOnCell(1:nLocalVertices,aCell)) * lambda)
      end do

      deallocate(lambda, pointVertex)
    end do

  end subroutine ocn_vector_vertex_to_cell_center !}}}

!***********************************************************************
!
!  routine ocn_second_order_shapiro_filter_ops
!
!> \brief   Do Ntimes simple shapiro filtering operations, but make
!>          higher order
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: multiple applications of digital shapiro filter (discrete Laplacian)
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
  subroutine ocn_second_order_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      type (field2DInteger), pointer, intent(in) :: boundaryVertex, boundaryCell
      type (field2DReal), pointer, intent(inout) :: ucReconstructX, ucReconstructY, ucReconstructZ ! cell center values
      integer, intent(in) :: Ntimes  ! number of filter applications
      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness

      type (field2DReal), pointer :: ucStore, vcStore, wcStore

      call mpas_pool_get_field(scratchPool,'ucX',ucStore)
      call mpas_pool_get_field(scratchPool,'ucY',vcStore)
      call mpas_pool_get_field(scratchPool,'ucZ',wcStore)
      call mpas_allocate_scratch_field(ucStore,.True.)
      call mpas_allocate_scratch_field(vcStore,.True.)
      call mpas_allocate_scratch_field(wcStore,.True.)


      call ocn_multiple_vector_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ)
      ucStore % array = 2.0_RKIND*ucReconstructX % array
      vcStore % array = 2.0_RKIND*ucReconstructY % array
      wcStore % array = 2.0_RKIND*ucReconstructZ % array
      call ocn_multiple_vector_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ)
      ucStore % array = ucStore % array - ucReconstructX % array
      vcStore % array = vcStore % array - ucReconstructY % array
      wcStore % array = wcStore % array - ucReconstructZ % array

      ! move temporary storage into final storage
      ucReconstructX % array = ucStore % array
      ucReconstructY % array = vcStore % array
      ucReconstructZ % array = wcStore % array

      ! deallocate temporary memory
      call mpas_deallocate_scratch_field(ucStore,.True.)
      call mpas_deallocate_scratch_field(vcStore,.True.)
      call mpas_deallocate_scratch_field(wcStore,.True.)

    end subroutine ocn_second_order_shapiro_filter_ops !}}}

!***********************************************************************
!
!  routine ocn_multiple_vector_shapiro_filter_ops
!
!> \brief   Do Ntimes simple shapiro filtering operations
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: multiple applications of digital shapiro filter (discrete Laplacian)
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
  subroutine ocn_multiple_vector_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      type (field2DInteger), pointer, intent(in) :: boundaryVertex, boundaryCell
      type (field2DReal), pointer, intent(inout) :: ucReconstructX, ucReconstructY, ucReconstructZ ! cell center values
      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness
      integer, intent(in) :: Ntimes  ! number of filter applications

      ! local variables
      integer atime

      do atime = 1,Ntimes
        !call ocn_simple_vector_shapiro_filter(meshPool, scratchPool, boundaryVertex, boundaryCell, &
        !  ucReconstructX, ucReconstructY, ucReconstructZ)
        call ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell % array, layerThickness, &
                                                ucReconstructX % array)
        call ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell % array, layerThickness, &
                                                ucReconstructY % array)
        call ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell % array, layerThickness, &
                                                ucReconstructZ % array)

      end do

    end subroutine ocn_multiple_vector_shapiro_filter_ops !}}}

!***********************************************************************
!
!  routine ocn_simple_vector_laplacian_filter
!
!> \brief   Do 1 pass of simple laplacian filter
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: one pass of digital shapiro filter (discrete Laplacian)
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
    subroutine ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell, layerThickness, ucReconstruct) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      integer, dimension(:,:), pointer, intent(in) :: boundaryCell
      real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: ucReconstruct
      real (kind=RKIND), dimension(:,:), pointer, intent(in) ::  layerThickness

      ! local variables
      type (field2DReal), pointer :: ucTemp
      integer :: aCell, aNeigh, aLevel
      integer, pointer :: nCellsSolve, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND) :: volSum, cellVol

      ! allocate scratch memory
      call mpas_pool_get_field(scratchPool, 'ucTemp', ucTemp)
      call mpas_allocate_scratch_field(ucTemp,.True.)

      ! get values from pools
      call mpas_pool_get_dimension(meshPool,'nCellsSolve',nCellsSolve)
      call mpas_pool_get_dimension(meshPool,'nVertLevels',nVertLevels)
      call mpas_pool_get_array(meshPool,'nEdgesOnCell',nEdgesOnCell)
      call mpas_pool_get_array(meshPool,'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool,'areaCell',areaCell)

      ucTemp % array = 0.0_RKIND

      ! perform laplacian filtering
      do aCell = 1, nCellsSolve
        do aLevel = 1, nVertLevels
          volSum = nEdgesOnCell(aCell) * layerThickness(aLevel,aCell) * areaCell(aCell) * (1-boundaryCell(aLevel,aCell))
          ucTemp % array(aLevel, aCell) = ucReconstruct(aLevel, aCell) * volSum
          if (volSum /= 0 ) then
            ! loop over all neighbors
            do aNeigh = 1, nEdgesOnCell(aCell)
              cellVol = layerThickness(aLevel,cellsOnCell(aNeigh,aCell)) * areaCell(cellsOnCell(aNeigh,aCell)) &
                * (1-boundaryCell(aLevel, cellsOnCell(aNeigh,aCell)))
              volSum = volSum + cellVol
              ucTemp % array(aLevel, aCell) = ucTemp % array(aLevel, aCell) + ucReconstruct(aLevel,cellsOnCell(aNeigh,aCell)) &
                                            * cellVol
            end do
            ucTemp % array(aLevel, aCell) = ucTemp % array(aLevel, aCell) / volSum
          end if
        end do
      end do

      ! exchange halo values
      call mpas_dmpar_exch_halo_field(ucTemp)

      ! replace input values with filtered values
      ucReconstruct = ucTemp % array

      ! deallocate scratch memory
      call mpas_deallocate_scratch_field(ucTemp , .True.)

    end subroutine ocn_simple_vector_laplacian_filter !}}}

!***********************************************************************
!
!  routine ocn_simple_vector_shapiro_filter
!
!> \brief   Do 1 pass of simple shapiro filter
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: one pass of digital shapiro filter to vertexes, back to cells
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
    subroutine ocn_simple_vector_shapiro_filter(meshPool, scratchPool, boundaryVertex, boundaryCell, &
        ucReconstructX, ucReconstructY, ucReconstructZ) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      type (field2DInteger), pointer, intent(in) :: boundaryVertex, boundaryCell
      type (field2DReal), pointer, intent(inout) :: ucReconstructX, ucReconstructY, ucReconstructZ ! cell center values

      ! local variables
      type (field2DReal), pointer :: uvX , uvY, uvZ ! cell center values

      ! allocate scratch memory
      call mpas_pool_get_field(scratchPool, 'uvX', uvX)
      call mpas_pool_get_field(scratchPool, 'uvY', uvY)
      call mpas_pool_get_field(scratchPool, 'uvZ', uvZ)
      call mpas_allocate_scratch_field(uvX,.True.)
      call mpas_allocate_scratch_field(uvY,.True.)
      call mpas_allocate_scratch_field(uvZ,.True.)

      uvX % array = 0.0_RKIND
      uvY % array = 0.0_RKIND
      uvZ % array = 0.0_RKIND

      ! perform filtering

      ! CC -> vertices
      call ocn_vector_cell_center_to_vertex(meshPool, boundaryVertex % array, boundaryCell % array, &
        ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array, &
        uvX % array, uvY % array, uvZ % array)
      ! do halo exchanges
      call mpas_dmpar_exch_halo_field(uvX)
      call mpas_dmpar_exch_halo_field(uvY)
      call mpas_dmpar_exch_halo_field(uvZ)
      ! vertices -> CC
      call ocn_vector_vertex_to_cell_center(meshPool, &
        uvX % array, uvY % array, uvZ % array, &
        ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array)
      ! do halo exchanges
      call mpas_dmpar_exch_halo_field(ucReconstructX)
      call mpas_dmpar_exch_halo_field(ucReconstructY)
      call mpas_dmpar_exch_halo_field(ucReconstructZ)

      ! N.B., effect of forgetting halo exchange may be subtle for a single pass

      ! deallocate scratch memory
      call mpas_deallocate_scratch_field(uvX , .True.)
      call mpas_deallocate_scratch_field(uvY , .True.)
      call mpas_deallocate_scratch_field(uvZ , .True.)

    end subroutine ocn_simple_vector_shapiro_filter !}}}

end module ocn_lagrangian_particle_tracking_interpolations

