! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_cice_advection_incremental_remap
!
!> \brief  incremental remapping transport scheme for cice
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  This module contains routines for transport of sea-ice mass and tracers
!>  using incremental remapping
!
!-----------------------------------------------------------------------

module cice_advection_incremental_remap

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timer

   use cice_advection_incremental_remap_tracers, only: &
        tracer_type, &
        tracersHead, &
        cice_add_tracers_to_linked_list, &
        cice_set_tracer_array_pointers, &
        cice_update_tracer_halo

   implicit none
   private
   save

   public :: &
        cice_init_advection_incremental_remap, &
        cice_run_advection_incremental_remap

   ! variables private to this module

   type :: geometric_avg_cell_type
      real(kind=RKIND), dimension(:), pointer :: x    => null()   ! grid cell average of x
      real(kind=RKIND), dimension(:), pointer :: y    => null()   ! grid cell average of y
      real(kind=RKIND), dimension(:), pointer :: xx   => null()   ! grid cell average of x^2
      real(kind=RKIND), dimension(:), pointer :: xy   => null()   ! grid cell average of x*y
      real(kind=RKIND), dimension(:), pointer :: yy   => null()   ! grid cell average of y^2
      real(kind=RKIND), dimension(:), pointer :: xxx  => null()   ! grid cell average of x^3
      real(kind=RKIND), dimension(:), pointer :: xxy  => null()   ! grid cell average of x^2*y
      real(kind=RKIND), dimension(:), pointer :: xyy  => null()   ! grid cell average of x*y^2
      real(kind=RKIND), dimension(:), pointer :: yyy  => null()   ! grid cell average of y^3
      real(kind=RKIND), dimension(:), pointer :: xxxx => null()   ! grid cell average of x^4
      real(kind=RKIND), dimension(:), pointer :: xxxy => null()   ! grid cell average of x^3*y
      real(kind=RKIND), dimension(:), pointer :: xxyy => null()   ! grid cell average of x^2*y^2
      real(kind=RKIND), dimension(:), pointer :: xyyy => null()   ! grid cell average of x*y^3
      real(kind=RKIND), dimension(:), pointer :: yyyy => null()   ! grid cell average of y^4
   end type geometric_avg_cell_type

   ! parameters private to this module
   real(kind=RKIND), parameter :: &
        eps11 = 1.0e-11_RKIND           ! small number (= puny in CICE)

   real(kind=RKIND), parameter :: &
        w1TriangleQP = 1.09951743655321885e-01_RKIND, &  ! weighting factors for triangle quadrature
        w2TriangleQP = 2.23381589678011389e-01_RKIND     ! (3rd and 4th degree polynomials)

   real(kind=RKIND), parameter :: &                      ! geometric factors for triangle quadrature
        q1TriangleQP = 9.15762135097710761e-02_RKIND, &  ! (3rd and 4th degree polynomials)
        q2TriangleQP = 8.16847572980458514e-01_RKIND, &
        q3TriangleQP = 1.08103018168070275e-01_RKIND, &
        q4TriangleQP = 4.45948490915965612e-01_RKIND

   integer, pointer :: &
        nQuadPoints   ! number of triangle quadrature points
                      ! = 6 for exact integration of polynomials up to degree 4
                      ! = 3 for exact integration of polynomials up to degree 2

   real(kind=RKIND), dimension(:), allocatable ::   &
        weightQuadPoint                               ! weighting factor for each quadrature point

   ! local diagnostic indices and IDs
   ! These are included for writing detailed diagnostic output for individual cells, edges and vertices
   integer :: ctest, etest, vtest   ! local cell, edge and vertex indices for diagnostics
   logical :: ctestOnProc,  etestOnProc,  vtestOnProc     ! true if ctestGlobal, etc. are on this processor node
                                                          ! Note: Can be true for more than one node (because of halo cells)
   integer :: ctestBlockID, etestBlockID, vtestBlockID    ! ID for blocks that own ctest, etest and vtest

   ! global indices for diagnostic cells, edges and vertices
   ! Uncomment to choose appropriate values for specific meshes

   ! default values that should work on any mesh
   integer, parameter :: ctestGlobal = 1
   integer, parameter :: etestGlobal = 1
   integer, parameter :: vtestGlobal = 1

   ! planar quad
!   integer, parameter :: ctestGlobal = 3244    ! test cell for planar quad ('square') test case
!   integer, parameter :: etestGlobal = 6529    ! test edge for planar quad test case
!   integer, parameter :: vtestGlobal = 3366    ! test vertex for planar quad test case

   ! planar hex
!   integer, parameter :: ctestGlobal = 3240     ! test cell for planar hex test case
!   integer, parameter :: etestGlobal = 9801     ! test edge for planer hex test case
!   integer, parameter :: vtestGlobal = 6720     ! test vertex for planar hex test case

   ! spherical quad (gx3)
!   integer, parameter :: ctestGlobal = 1        ! test cell for spherical quad test case (gx3, rotated, near SP)
!   integer, parameter :: etestGlobal = 13898    ! test edge for spherical quad test case
!   integer, parameter :: vtestGlobal = 1        ! test vertex for spherical quad test case
!!   integer, parameter :: ctestGlobal = 5411    ! test cell for spherical quad test case (gx3_noland, near equator)
!!   integer, parameter :: etestGlobal = 5706    ! test edge for spherical quad test case
!!   integer, parameter :: vtestGlobal = 5511    ! test vertex for spherical quad test case

   ! spherical quad (gx1)
!   integer, parameter :: ctestGlobal = 2         ! test cell for spherical quad test case (gx1, rotated, near SP)
!   integer, parameter :: etestGlobal = 124953    ! test edge for spherical quad test case
!   integer, parameter :: vtestGlobal = 2         ! test vertex for spherical quad test case
!!   integer, parameter :: ctestGlobal = 51157    ! test cell for spherical quad test case (gx1, not rotated)
!!   integer, parameter :: etestGlobal = 126267   ! test edge for spherical quad test case
!!   integer, parameter :: vtestGlobal = 52059    ! test vertex for spherical quad test case

   ! spherical hex (2562 mesh)
!   integer, parameter :: ctestGlobal = 3          ! test cell for spherical hex test case (North Pole on 2562 grid)
!   integer, parameter :: etestGlobal = 473        ! test edge for spherical hex test case
!   integer, parameter :: vtestGlobal = 12         ! test vertex for spherical hex test case
!!   integer, parameter :: ctestGlobal = 2354       ! test cell for spherical hex test case (rotated)
!!   integer, parameter :: etestGlobal = 3473       ! test edge for spherical hex test case
!!   integer, parameter :: vtestGlobal = 2585       ! test vertex for spherical hex test case

   integer, parameter :: iCatTest = 1
   integer, parameter :: iLayerTest = 1

   ! verbose options; set to true for extensive diagnostic output
   ! WHL: I often debug with verboseInit/Run/Construct/Geometry/Fluxes/Global set to true.
!   logical, parameter :: verboseInit = .true.
!   logical, parameter :: verboseRun = .true.
!   logical, parameter :: verboseConstruct = .true.
!   logical, parameter :: verboseGeometry = .true.
!   logical, parameter :: verboseFluxes = .true.
!   logical, parameter :: verboseGlobal = .true.
   logical, parameter :: verboseInit = .false.
   logical, parameter :: verboseRun = .false.
   logical, parameter :: verboseConstruct = .false.
   logical, parameter :: verboseGeometry = .false.
   logical, parameter :: verboseFluxes = .false.
   logical, parameter :: verboseGlobal = .false.
   logical, parameter :: verboseGeomAvg = .false.

   logical :: first_call = .true.  ! set to false after the first call

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cice_init_advection_incremental_remap
!
!> \brief MPAS-CICE incremental remapping initialization
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  Given input mesh quantities, this routine computes some geometric
!>  quantities required for incremental remapping.
!
!-----------------------------------------------------------------------

  subroutine cice_init_advection_incremental_remap(domain)

    use mpas_rbf_interpolation, only: mpas_rbf_interp_initialize
    use mpas_vector_reconstruction, only: mpas_init_reconstruct

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    ! local variables

    type(block_type), pointer :: block

    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: rotatedMeshPool
    type (mpas_pool_type), pointer :: incrementalRemapPool
    type (mpas_pool_type), pointer :: configPool

    real(kind=RKIND), dimension(:,:), pointer ::  &
         xVertexOnCell, yVertexOnCell,  & ! local x/y coordinates of vertices relative to cell center
         xVertexOnEdge, yVertexOnEdge     ! local x/y coordinates of vertices relative to edge

    type(geometric_avg_cell_type), pointer :: &
         geomAvgCell          ! derived type holding geometric averages

    integer, pointer ::   &
         nCells,          &   ! number of cells
         nVertices,       &   ! number of vertices
         nEdges,          &   ! number of edges
         nCellsSolve,     &   ! number of locally owned cells
         nVerticesSolve,  &   ! number of locally owned vertices
         nEdgesSolve,     &   ! number of locally owned edges
         maxEdges,        &   ! max number of edges per cell
         vertexDegree         ! number of edges that meet at each vertex

    integer, dimension(:), pointer :: &
         nEdgesOnCell,      & ! number of edges for each cell
         remapEdge            ! = 1 if IR fluxes need to be computed across an edge, else = 0
                              ! = 1 for all edges of locally owned cells, provided a cell exists on each side

    integer, dimension(:,:), pointer ::  &
         cellsOnCell,       & ! index for each neighbor cell of a given cell
         edgesOnCell,       & ! index for each edge of a given cell
         verticesOnCell,    & ! index for each vertex of a given cell
         verticesOnEdge,    & ! index for each vertex of a given edge
         cellsOnEdge,       & ! index for each cell on a given edge
         cellsOnVertex,     & ! index for each cell on a given vertex
         edgesOnVertex,     & ! index for each edge neighbor of a given vertex
         cellsOnEdgeRemap,  & ! cell index for each cell neighbor of a given edge,
                              !  including cells that share a vertex with the edge
         edgesOnEdgeRemap     ! edge index for nearest edge neighbors of a given edge
                              ! (i.e., those sharing a vertex with the edge)

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell,      &  ! cell area
         dcEdge,        &  ! distance between the 2 cell centers on each side of an edge
         dvEdge            ! distance between the 2 vertices of an edge

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, yCell, zCell,         & ! global coordinates of cell centers
         xVertex, yVertex, zVertex,   & ! global coordinates of vertices
         xEdge, yEdge, zEdge            ! global coordinates of edge midpoints

    real(kind=RKIND), dimension(:), pointer :: &
         xCellRotate, yCellRotate, zCellRotate,         & ! rotated global coordinates of cell centers
         xVertexRotate, yVertexRotate, zVertexRotate,   & ! rotated global coordinates of vertices
         xEdgeRotate, yEdgeRotate, zEdgeRotate            ! rotated global coordinates of edge midpoints

    real(kind=RKIND), dimension(:), pointer :: &
         minLengthEdgesOnVertex         ! minimum length of the edges on each vertex
                                        ! used for CFL diagnostics

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         transGlobalToCell,    &  ! 3x3 matrix for transforming vectors from global to cell-based coordinates
         transGlobalToVertex,  &  ! 3x3 matrix for transforming vectors from global to vertex-based coordinates
         transGlobalToEdge        ! 3x3 matrix for transforming vectors from global to edge-based coordinates

    !Note: The matrices transCellToGlobal, transVertexToGlobal and transEdgeToGlobal are not currently used,
    !      but they are computed just in case.
    real(kind=RKIND), dimension(:,:,:), pointer :: &
         transCellToGlobal,    &  ! 3x3 matrix for transforming vectors from cell-based to global coordinates
         transVertexToGlobal,  &  ! 3x3 matrix for transforming vectors from vertex-based to global coordinates
         transEdgeToGlobal        ! 3x3 matrix for transforming vectors from edge-based to global coordinates

    real(kind=RKIND), dimension(:,:,:), pointer ::  &
         coeffsReconstruct        ! coefficients for reconstructing edge-based fields at cell centers

    logical, pointer :: &
         on_a_sphere,                & ! if true, then the mesh lives on the surface of a sphere; else in a plane
         config_rotate_cartesian_grid  ! if true, then rotate the true North and South Poles to the equator

    ! IDs for diagnostics
    integer, dimension(:), pointer ::  &
         indexToCellID,   & ! global index for each cell on a block
         indexToVertexID, & ! global index for each vertex on a block
         indexToEdgeID      ! global index for each edge on a block

    integer :: i, n, iCell, iEdge, iVertex, iCellNeighbor, iEdgeNeighbor
    integer :: iVertexOnCell, iEdgeOnCell, iCellOnEdge, iCellOnCell, iVertexOnEdge, iEdgeOnEdge, iEdgeOnVertex
    integer :: nCellsOnEdgeRemap, nEdgesOnEdgeRemap, nVerticesOnEdgeRemap

    character(len=strKIND) :: &
         errorMessage ! error message for global abort

    ! check the number of halos against the vertex degree
    call check_halo_layer_number(domain)

    if (verboseInit) then
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'In cice_init_advection_incremental_remap'
       write(stderrUnit,*) 'Add tracers to linked list'
    endif

    ! construct a linked list of tracers
    call cice_add_tracers_to_linked_list(domain)

    ! assign pointers
    block => domain % blocklist

    ! loop over the blocks in the domain

    do while (associated(block))

       if (verboseInit) write(stderrUnit,*) 'Get arrays from pools'

       ! get pools
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'rotated_mesh', rotatedMeshPool)
       call mpas_pool_get_subpool(block % structs, 'incremental_remap', incrementalRemapPool)

       ! get config options
       configPool => block % configs
       call mpas_pool_get_config(configPool, 'config_rotate_cartesian_grid', config_rotate_cartesian_grid)
       call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

       if (verboseInit) write(stderrUnit,*) 'rotate =', config_rotate_cartesian_grid

       ! get mesh dimensions
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
       call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
       call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
       call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

       ! get some mesh arrays
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
       call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
       call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
       call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
       call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'zCell', zCell)
       call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
       call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
       call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
       call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
       call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
       call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)
       call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)

       ! If on a rotated spherical grid, get arrays for rotated cell, vertex and edge vectors.
       ! These vectors will be rotated to put the North and South Poles at the equator.
       ! After rotation, we will assign the standard pointers (without the 'Rotate' suffix)
       !  to these vectors so that subroutine arguments are the same with or without rotation.

       if (config_rotate_cartesian_grid .and. on_a_sphere) then
          call mpas_pool_get_array(rotatedMeshPool, 'xCellRotate', xCellRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'yCellRotate', yCellRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'zCellRotate', zCellRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'xVertexRotate', xVertexRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'yVertexRotate', yVertexRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'zVertexRotate', zVertexRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'xEdgeRotate', xEdgeRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'yEdgeRotate', yEdgeRotate)
          call mpas_pool_get_array(rotatedMeshPool, 'zEdgeRotate', zEdgeRotate)
       endif

       ! get IR arrays
       allocate(geomAvgCell)
       call get_incremental_remap_geometry_pointers(geomAvgCell, block)

       call mpas_pool_get_array(incrementalRemapPool, 'transCellToGlobal', transCellToGlobal)
       call mpas_pool_get_array(incrementalRemapPool, 'transGlobalToCell', transGlobalToCell)
       call mpas_pool_get_array(incrementalRemapPool, 'transVertexToGlobal', transVertexToGlobal)
       call mpas_pool_get_array(incrementalRemapPool, 'transGlobalToVertex', transGlobalToVertex)
       call mpas_pool_get_array(incrementalRemapPool, 'transEdgeToGlobal', transEdgeToGlobal)
       call mpas_pool_get_array(incrementalRemapPool, 'transGlobalToEdge', transGlobalToEdge)

       call mpas_pool_get_array(incrementalRemapPool, 'xVertexOnCell', xVertexOnCell)
       call mpas_pool_get_array(incrementalRemapPool, 'yVertexOnCell', yVertexOnCell)
       call mpas_pool_get_array(incrementalRemapPool, 'xVertexOnEdge', xVertexOnEdge)
       call mpas_pool_get_array(incrementalRemapPool, 'yVertexOnEdge', yVertexOnEdge)

       call mpas_pool_get_array(incrementalRemapPool, 'remapEdge', remapEdge)
       call mpas_pool_get_array(incrementalRemapPool, 'cellsOnEdgeRemap', cellsOnEdgeRemap)
       call mpas_pool_get_array(incrementalRemapPool, 'edgesOnEdgeRemap', edgesOnEdgeRemap)
       call mpas_pool_get_array(incrementalRemapPool, 'minLengthEdgesOnVertex', minLengthEdgesOnVertex)

       ! Find the local index, block ID and proc ID associated with global diagnostic cells, edges and vertices.
       ! Note: If ctest/etest/vtest lie on a processor other than the head node, then diagnostic output will be written
       !       only if the code is built in debug mode.

       if (ctestGlobal > 0) then

          call get_local_index(&
               block,          &
               nCells,         &
               indexToCellID,  &
               ctestGlobal,    &
               ctest,          &
               ctestBlockID,   &
               ctestOnProc)

          if (verboseInit) then
             if (ctestOnProc) then
                write(stderrUnit,*) 'ctestGlobal, local ctest, block:', ctestGlobal, ctest, ctestBlockID
             else
                write(stderrUnit,*) 'Not on this proc: ctestGlobal =', ctestGlobal
             endif
          endif

       endif   ! ctestGlobal > 0

       if (etestGlobal > 0) then

          call get_local_index(&
               block,          &
               nEdges,         &
               indexToEdgeID,  &
               etestGlobal,    &
               etest,          &
               etestBlockID,   &
               etestOnProc)

          if (verboseInit) then
             if (etestOnProc) then
                write(stderrUnit,*) 'etestGlobal, local etest, block:', etestGlobal, etest, etestBlockID
             else
                write(stderrUnit,*) 'Not on this proc: etestGlobal =', etestGlobal
             endif
          endif

       endif   ! etestGlobal > 0

       if (vtestGlobal > 0) then

          call get_local_index(&
               block,          &
               nVertices,      &
               indexToVertexID,&
               vtestGlobal,    &
               vtest,          &
               vtestBlockID,   &
               vtestOnProc)

          if (verboseInit) then
             if (vtestOnProc) then
                write(stderrUnit,*) 'vtestGlobal, local vtest, block:', vtestGlobal, vtest, vtestBlockID
             else
                write(stderrUnit,*) 'Not on this proc: vtestGlobal =', vtestGlobal
             endif
          endif

       endif   ! vtestGlobal > 0

       ! If on a rotated spherical grid, then rotate the global vectors pointing to cells, vertices and edges

       if (config_rotate_cartesian_grid .and. on_a_sphere) then

          call rotate_global_vectors(&
               xCell,       yCell,       zCell, &
               xCellRotate, yCellRotate, zCellRotate, &
               nCells)

          call rotate_global_vectors(&
               xVertex,       yVertex,       zVertex, &
               xVertexRotate, yVertexRotate, zVertexRotate, &
               nVertices)

          call rotate_global_vectors(&
               xEdge,       yEdge,       zEdge, &
               xEdgeRotate, yEdgeRotate, zEdgeRotate, &
               nEdges)

          if (verboseInit) then
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'Rotated the global Cell, Vertex and Edge vectors'
             if (ctestOnProc .and. block % localBlockID == ctestBlockID) then
                iCell = ctest
                write(stderrUnit,*) 'iCell, xCell, yCell, zCell:', &
                                     iCell, xCell(iCell), yCell(iCell), zCell(iCell)
                write(stderrUnit,*) 'Rotated iCell, xCell, yCell, zCell:', &
                                     iCell, xCellRotate(iCell), yCellRotate(iCell), zCellRotate(iCell)
             endif
          endif

          ! We do not need the unrotated vectors again in this subroutine, so reassign the local pointers
          !  to point to the rotated vectors.

          xCell => xCellRotate
          yCell => yCellRotate
          zCell => zCellRotate
          xEdge => xEdgeRotate
          yEdge => yEdgeRotate
          zEdge => zEdgeRotate
          xVertex => xVertexRotate
          yVertex => yVertexRotate
          zVertex => zVertexRotate

       endif  ! config_rotate_cartesian_grid and on_a_sphere

       ! set up some matrices for transformation between global coordinates and local coordinates

       if (on_a_sphere) then

          call define_local_to_global_transformations(&
               xCell, yCell, zCell,  &
               nCells,               &
               transCellToGlobal,    &
               transGlobalToCell)

          call define_local_to_global_transformations(&
               xVertex, yVertex, zVertex,  &
               nVertices,                  &
               transVertexToGlobal,        &
               transGlobalToVertex)

          call define_local_to_global_transformations(&
               xEdge, yEdge, zEdge,  &
               nEdges,               &
               transEdgeToGlobal,    &
               transGlobalToEdge)

          if (verboseInit .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
             iCell = ctest
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'Defined local to global transformations'
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'iCell, transCellToGlobal'
             do n = 1, 3
                write(stderrUnit,*) transCellToGlobal(n,:,iCell)
             enddo
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'iCell, transGlobalToCell'
             do n = 1, 3
                write(stderrUnit,*) transGlobalToCell(n,:,iCell)
             enddo
          endif

       endif   ! on a sphere

       ! compute local vertex coordinates relative to cell center
       ! Note: This calculation could be done in subroutine get_geometry_incremental_remap,
       !       but I kept it separate in case we want to make it available to other parts of the code.
       call get_vertex_on_cell_coordinates(&
            nCells,                       &
            nEdgesOnCell,                 &
            verticesOnCell,               &
            xVertex, yVertex, zVertex,    &
            xCell,   yCell,   zCell,      &
            on_a_sphere,                  &
            transGlobalToCell,            &
            xVertexOnCell, yVertexOnCell)

       ! compute other geometric data structures used for IR
       call get_geometry_incremental_remap(&
            nEdges,                    &
            nVertices,                 &
            nCells,                    &
            nCellsSolve,               &
            nEdgesOnCell,              &
            maxEdges,                  &
            xCell,   yCell,   zCell,   &
            xVertex, yVertex, zVertex, &
            xEdge,   yEdge,   zEdge,   &
            indexToEdgeID,             &
            vertexDegree,              &
            cellsOnEdge,               &
            verticesOnEdge,            &
            verticesOnCell,            &
            edgesOnCell,               &
            edgesOnVertex,             &
            cellsOnVertex,             &
            xVertexOnCell, yVertexOnCell, &
            on_a_sphere,               &
            transGlobalToEdge,         &
            transGlobalToVertex,       &
            remapEdge,                 &
            cellsOnEdgeRemap,          &
            edgesOnEdgeRemap,          &
            xVertexOnEdge, yVertexOnEdge, &
            minLengthEdgesOnVertex)

       ! Write diagnostics for the test cell
       if (verboseInit .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then

          iCell = ctest
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'iCell (local ID, global ID) =', iCell, indexToCellID(iCell)

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'cellsOnCell:'
          do iCellOnCell = 1, nEdgesOnCell(iCell)
             write(stderrUnit,*) iCellOnCell, cellsOnCell(iCellOnCell, iCell), &
                  indexToCellID(cellsOnCell(iCellOnCell,iCell))
          enddo

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'edgesOnCell:'
          do iEdgeOnCell = 1, nEdgesOnCell(iCell)
             write(stderrUnit,*) iEdgeOnCell, edgesOnCell(iEdgeOnCell, iCell), &
                  indexToEdgeID(edgesOnCell(iEdgeOnCell,iCell))
          enddo

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'verticesOnCell:'
          do iVertexOnCell = 1, nEdgesOnCell(iCell)
             write(stderrUnit,*) iVertexOnCell, verticesOnCell(iVertexOnCell,iCell), &
                  indexToVertexID(verticesOnCell(iVertexOnCell,iCell))
          enddo

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'xVertexOnCell, yVertexOnCell:'
          do iVertexOnCell = 1, nEdgesOnCell(iCell)
             write(stderrUnit,*) iVertexOnCell, xVertexOnCell(iVertexOnCell, iCell), &
                  yVertexOnCell(iVertexOnCell, iCell)
          enddo

       endif   ! ctestOnProc & ctestBlockID

       ! Write diagnostics for the test edge
       if (verboseInit .and. etestOnProc .and. block % localBlockID == etestBlockID) then

          iEdge = etest
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'iEdge (local ID, global ID) =', iEdge, indexToEdgeID(iEdge)

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'cellsOnEdge:'
          do iCellOnEdge = 1, 2
             write(stderrUnit,*) iCellOnEdge, cellsOnEdge(iCellOnEdge, iEdge), &
                  indexToCellID(cellsOnEdge(iCellOnEdge, iEdge))
          enddo

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'verticesOnEdge:'
          do iVertexOnEdge = 1, 2
             write(stderrUnit,*) iVertexOnEdge, verticesOnEdge(iVertexOnEdge, iEdge), &
                  indexToVertexID(verticesOnEdge(iVertexOnEdge, iEdge))
          enddo

          if (vertexDegree == 3) then   ! hex mesh
             nCellsOnEdgeRemap = 4
             nEdgesOnEdgeRemap = 4
             nVerticesOnEdgeRemap = 6
          elseif (vertexDegree == 4) then   ! quad mesh
             nCellsOnEdgeRemap = 6
             nEdgesOnEdgeRemap = 6
             nVerticesOnEdgeRemap = 8
          endif

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'cellsOnEdgeRemap:'
          do iCellOnEdge = 1, nCellsOnEdgeRemap
             iCellNeighbor = cellsOnEdgeRemap(iCellOnEdge,iEdge)
             if (iCellNeighbor >= 1) then
                write(stderrUnit,*) iCellOnEdge, iCellNeighbor, indexToCellID(iCellNeighbor)
             else  ! no global ID if local ID = 0
                write(stderrUnit,*) iCellOnEdge, iCellNeighbor
             endif
          enddo

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'edgesOnEdgeRemap:'
          do iEdgeOnEdge = 1, nEdgesOnEdgeRemap
             iEdgeNeighbor = edgesOnEdgeRemap(iEdgeOnEdge,iEdge)
             if (iEdgeNeighbor >= 1) then
                write(stderrUnit,*) iEdgeOnEdge, iEdgeNeighbor, indexToEdgeID(iEdgeNeighbor)
             else  ! no global ID if local ID = 0
                write(stderrUnit,*) iEdgeOnEdge, iEdgeNeighbor
             endif
          enddo

          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'xVertexOnEdge, yVertexOnEdge:', indexToEdgeID(iEdge)
          write(stderrUnit,*) 'vertex, x, y:'
          do iVertexOnEdge = 1, nVerticesOnEdgeRemap
             write(stderrUnit,*) iVertexOnEdge, xVertexOnEdge(iVertexOnEdge,iEdge), &
                  yVertexOnEdge(iVertexOnEdge,iEdge)
          enddo

          write(stderrUnit,*) ' '

       endif  ! etestOnProc & etestBlockID

       ! Write diagnostics for the test vertex
       if (verboseInit .and. vtestOnProc .and. block % localBlockID == vtestBlockID) then

          iVertex = vtest
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'iVertex (local ID, global ID) =', iVertex, indexToVertexID(iVertex)
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'edgesOnVertex:'
          do iEdgeOnVertex = 1, vertexDegree
             write(stderrUnit,*) iEdgeOnVertex, edgesOnVertex(iEdgeOnVertex,iVertex), &
                  indexToEdgeID(edgesOnVertex(iEdgeOnVertex,iVertex))
          enddo

       endif

       ! compute geometric means of x, y, xx, xy, yy, etc.

       call compute_geometric_cell_averages(&
            xVertexOnCell, yVertexOnCell,  &
            nCells,                        &
            maxEdges,                      &
            edgesOnCell,                   &
            nEdgesOnCell,                  &
            dcEdge,                        &
            dvEdge,                        &
            geomAvgCell)

       if (verboseInit .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
          iCell = ctest
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'Grid cell averages, iCell =', indexToCellID(iCell)
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'x', geomAvgCell % x(iCell)
          write(stderrUnit,*) 'y', geomAvgCell % y(iCell)
          write(stderrUnit,*) 'xx', geomAvgCell % xx(iCell)
          write(stderrUnit,*) 'xy', geomAvgCell % xy(iCell)
          write(stderrUnit,*) 'yy', geomAvgCell % yy(iCell)
          write(stderrUnit,*) 'xxx', geomAvgCell % xxx(iCell)
          write(stderrUnit,*) 'xxy', geomAvgCell % xxy(iCell)
          write(stderrUnit,*) 'xyy', geomAvgCell % xyy(iCell)
          write(stderrUnit,*) 'yyy', geomAvgCell % yyy(iCell)
          write(stderrUnit,*) 'xxxx', geomAvgCell % xxxx(iCell)
          write(stderrUnit,*) 'xxxy', geomAvgCell % xxxy(iCell)
          write(stderrUnit,*) 'xxyy', geomAvgCell % xxyy(iCell)
          write(stderrUnit,*) 'xyyy', geomAvgCell % xyyy(iCell)
          write(stderrUnit,*) 'yyyy', geomAvgCell % yyyy(iCell)
       endif   ! ctestOnProc & ctestBlockID

       ! compute coefficients needed for reconstructing tracer gradients at cell centers
       ! Note: This field gets a halo update outside the block loop.
       ! Note: The 3D gradient vector that is computed using these coefficients lives in
       !       Cartesian x-y-z space. If working on a rotated sea ice mesh, a rotation is needed
       !       before transforming to the local tangent plane.  This is done in the compute_gradient
       !       subroutines below.
       ! Note: The operator subroutines assume this field is called 'coeffs_reconstruct'.
       !       It would be better to call it 'coeffsReconstruct', but this would require a change
       !       in the operators.

       call mpas_rbf_interp_initialize(meshPool)
       call mpas_init_reconstruct(meshPool)
       call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffsReconstruct)

       if (verboseInit .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
          iCell = ctest
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'Reconstruction coefficients, iCell =', indexToCellID(iCell)
          write(stderrUnit,*) ' '
          do iEdgeOnCell = 1, nEdgesOnCell(iCell)
             write(stderrUnit,*) iEdgeOnCell, coeffsReconstruct(:,iEdgeOnCell,iCell)
          enddo
       endif   ! ctestOnProc & ctestBlockID

       ! clean up
       deallocate(geomAvgCell)

       block => block % next

    enddo  ! associated(block)

    ! Do a halo update for the reconstruction coefficients.
    ! Otherwise the coefficients are zero in halo cells.
    call MPAS_dmpar_field_halo_exch(domain, 'coeffs_reconstruct')

    ! get the number of points for triangle quadrature
    ! = 6 for exact integration of polynomials up to degree 4
    ! = 3 for exact integration of polynomials up to degree 2
    call mpas_pool_get_dimension(meshPool, 'nQuadPoints', nQuadPoints)

    ! set weights for triangle quadrature

    allocate(weightQuadPoint(nQuadPoints))
    if (nQuadPoints == 3) then
       weightQuadPoint(:) = 1.0_RKIND / 3.0_RKIND
    elseif (nQuadPoints == 6) then
       weightQuadPoint(1:3) = w1TriangleQP    ! 0.010995...
       weightQuadPoint(4:6) = w2TriangleQP    ! 0.022338...
    else
       write(stderrUnit,*) 'ERROR: must have nQuadPoints = 3 or 6 for IR transport'
       call mpas_dmpar_global_abort('MPAS-seaice: IR must have nQuadPoints = 3 or 6 for IR transport')
    endif

    ! check for non-positive cell areas
    ! Note: As of Sept. 2015, the mesh converter tool had a bug that can give negative areas
    !       for some cells on a spherical quad mesh.
    do iCell = 1, nCellsSolve
       if (areaCell(iCell) <= 0.0_RKIND) then
          write(stderrUnit,*) ' '
          write(errorMessage,*) 'Bad cell area, iCell, area:', iCell, areaCell(iCell)
          write(stderrUnit,*) trim(errorMessage)
!!          write(stderrUnit,*) 'cellsOnCell:', cellsOnCell(:,iCell)
          call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
       endif
    enddo

    ! uncomment to get a list of cells and edges with local and global IDs
    if (verboseInit) then
       write(stderrUnit,*) 'nCells, nEdges, nVertices:', nCells, nEdges, nVertices
       write(stderrUnit,*) 'nCellsSolve, nEdgesSolve, nVerticesSolve:', nCellsSolve, nEdgesSolve, nVerticesSolve
       write(stderrUnit,*) ' '

       !  write(stderrUnit,*) 'Cells (local ID, global ID)'
       do iCell = 1, nCells
          !  write(stderrUnit,*) iCell, indexToCellID(iCell)
       enddo

       !  write(stderrUnit,*) 'Edges (local ID, global ID)'
       do iEdge = 1, nEdges
          !  write(stderrUnit,*) iEdge, indexToEdgeID(iEdge)
       enddo
    endif

  end subroutine cice_init_advection_incremental_remap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine check_halo_layer_number
!
!> \brief Checks the number of halo layers against vertex degree
!> \author Adrian K. Turner
!> \date   5th November 2015
!> \details
!>
!-----------------------------------------------------------------------

  subroutine check_halo_layer_number(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    integer, pointer :: &
         config_num_halos, &
         vertexDegree

    integer :: &
         numHalosOptimal

    ! get the number of halo layers
    call MPAS_pool_get_config(domain % blocklist % configs, "config_num_halos", config_num_halos)

    ! get the vertex degree
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "vertexDegree", vertexDegree)

    ! determine the optimal number of halo layers
    if (vertexDegree == 3) then
       numHalosOptimal = 2 ! SCVT grids
    else if (vertexDegree == 4) then
       numHalosOptimal = 3 ! Quad grids
    endif

    ! check the number of halo layers against the vertex degree
    if (config_num_halos < numHalosOptimal) then

       ! too few layers
       write(stderrUnit,*) "Error: too few halo layers for incremental remapping"
       write(stderrUnit,*) "-- vertexDegree:     ", vertexDegree
       write(stderrUnit,*) "-- config_num_halos: ", config_num_halos
       write(stderrUnit,*) "-- numHalosOptimal:  ", numHalosOptimal
       call MPAS_dmpar_global_abort("MPAS-seaice: IR too few halo layers for incremental remapping")

    else if (config_num_halos > numHalosOptimal) then

       ! too many layers
       write(stderrUnit,*) "Warning: too many halo layers for incremental remapping"
       write(stderrUnit,*) "-- vertexDegree:     ", vertexDegree
       write(stderrUnit,*) "-- config_num_halos: ", config_num_halos
       write(stderrUnit,*) "-- numHalosOptimal:  ", numHalosOptimal

    endif

  end subroutine check_halo_layer_number

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine get_local_index
!
!> \brief  find local index, blockID and procID given a global index
!> \author William Lipscomb
!> \date   June 2015
!> \details
!>  This routine finds the local index, block ID and proc ID for an element
!>  (edge, vertex or cell), given the global index.
!-----------------------------------------------------------------------

  subroutine get_local_index(&
       block,              &
       nElements,          &
       indexToElementID,   &
       testElementGlobal,  &
       testElementLocal,   &
       testBlockID,        &
       testOnProc)

    type(block_type), intent(in) :: block

    integer, intent(in) :: &
         nElements             !< Input:  number of elements (cells, vertices or edges)

    integer, dimension(:), intent(in) ::  &
         indexToElementID      !< Input:  global index corresponding to each element on a block

    integer, intent(in) ::  &
         testElementGlobal     !< Input:  global index for test element

    integer, intent(out) ::  &
         testElementLocal, &   !< Output: local index for test element
         testBlockId           !< Output: block ID for test element

    logical, intent(out) ::  &
         testOnProc            !< Output: true if testElementGlobal is on this processor node

    integer :: iElement

    ! initialize output
    testElementLocal = 0
    testBlockID = -999  ! not 0, because 0 is an actual block ID
    testOnProc = .false.

    ! loop over elements and check whether the test element is on this block and proc

    do iElement = 1, nElements
       if (indexToElementID(iElement) == testElementGlobal) then
          testElementLocal = iElement
          testBlockID = block % localBlockID
          testOnProc = .true.
       endif
    enddo

  end subroutine get_local_index

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine rotate_global_vectors
!
!> \brief  rotate global vectors through 90 degrees
!> \author William Lipscomb
!> \date   June 2015
!> \details
!>  This routine rotates global vectors through 90 degrees such that the
!>  true poles have z = 0 and the true equator has x = 0.
!>  It is similar to cice_grid_rotation_forward, except that it operates
!>  on a global array instead of a single vector.
!-----------------------------------------------------------------------

  subroutine rotate_global_vectors(&
       xElement,       yElement,       zElement, &
       xElementRotate, yElementRotate, zElementRotate, &
       nElements)

    real(kind=RKIND), dimension(:), intent(in) :: &
         xElement, yElement, zElement   !< Input: components of a global vector pointing from the center of a sphere
                                        !         to an element (cell, vertex or edge) on its surface

    real(kind=RKIND), dimension(:), intent(out) :: &
         xElementRotate, yElementRotate, zElementRotate   !< Output: global vector rotated by 90 degrees about the y-axis

    integer, intent(in) :: nElements    ! number of elements

    integer :: iElement

    ! Rotate the vectors by 90 degrees, such that the true North Pole goes to (-1,0,0)
    ! and the true South Pole to (1,0,0).  The true equator will have x = 0.

    do iElement = 1, nElements

       xElementRotate(iElement) = -zElement(iElement)
       yElementRotate(iElement) =  yElement(iElement)
       zElementRotate(iElement) =  xElement(iElement)

    enddo

  end subroutine rotate_global_vectors

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine define_local_to_global_transformations
!
!> \brief  compute matrices for transforming between local and global coordinates
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine computes 3x3 matrices for transforming between global xyz coordinates
!>  and local coordinates based at a cell center, vertex or edge.
!
!-----------------------------------------------------------------------

  subroutine define_local_to_global_transformations(&
       xGlobal, yGlobal, zGlobal,  &
       nElements,                  &
       transLocalToGlobal,         &
       transGlobalToLocal)

    real(kind=RKIND), dimension(nElements), intent(in) ::  &
         xGlobal, yGlobal, zGlobal  !< Input: global coordinates of points on the surface of the sphere

    integer, intent(in) :: &
         nElements                  !< Input: number of points (typically nCells, nVertices or nEdges)

    real(kind=RKIND), dimension(3, 3, nElements), intent(out) ::  &
         transLocalToGlobal         !< Output: 3 x 3 matrix for transforming from local to global coordinates

    real(kind=RKIND), dimension(3, 3, nElements), intent(out) ::  &
         transGlobalToLocal         !< Output: 3 x 3 matrix for transforming from global to local coordinates

    ! local variables

    integer :: iElement

    real(kind=RKIND), dimension(3,3) :: &
         unitVectorGlobal,       &  ! unit vectors for global coordinate system
         unitVectorLocal            ! unit vectors for local coordinate system
                                    ! (expressed in global coordinates)

    transLocalToGlobal(:,:,:) = 0.0_RKIND
    transGlobalToLocal(:,:,:) = 0.0_RKIND

    ! set unit vectors for the global coordinate system

    ! global x
    unitVectorGlobal(1,1) = 1.0_RKIND
    unitVectorGlobal(2,1) = 0.0_RKIND
    unitVectorGlobal(3,1) = 0.0_RKIND

    ! global y
    unitVectorGlobal(1,2) = 0.0_RKIND
    unitVectorGlobal(2,2) = 1.0_RKIND
    unitVectorGlobal(3,2) = 0.0_RKIND

    ! global z
    unitVectorGlobal(1,3) = 0.0_RKIND
    unitVectorGlobal(2,3) = 0.0_RKIND
    unitVectorGlobal(3,3) = 1.0_RKIND

    do iElement = 1, nElements

       ! compute unit vectors for the local coordinate system
       ! Note: If using a rotating Cartesian grid, then the input global vectors must already have been rotated.

       ! local vertical
       ! this is the same as the vector from the center of the sphere to the point on the surface
       unitVectorLocal(1,3) = xGlobal(iElement)
       unitVectorLocal(2,3) = yGlobal(iElement)
       unitVectorLocal(3,3) = zGlobal(iElement)
       call unit_vector_3d(unitVectorLocal(:,3))

       if (abs(xGlobal(iElement)**2 + yGlobal(iElement)**2) > eps11) then  ! not at a pole

          ! local east
          ! this is obtained by a 90-degree rotation of (xGlobal, yGlobal) in a plane of constant latitude

          unitVectorLocal(1,1) = -yGlobal(iElement)
          unitVectorLocal(2,1) =  xGlobal(iElement)
          unitVectorLocal(3,1) =  0.0_RKIND
          call unit_vector_3d(unitVectorLocal(:,1))

          ! local north
          ! this is found by taking the cross product of the local vertical and local east unit vectors
          call cross_product_3d(unitVectorLocal(:,3), unitVectorLocal(:,1), unitVectorLocal(:,2))

       else  ! xGlobal = yGlobal = 0; at a pole

          ! take local east and local north to lie in a plane parallel to the equator, as in unitVectorGlobal
          if (zGlobal(iElement) > 0.0_RKIND) then  ! North Pole
             unitVectorLocal(:,1) = unitVectorGlobal(:,1)
             unitVectorLocal(:,2) = unitVectorGlobal(:,2)
          else   ! South Pole (reverse the order of east and north to obtain right-handed coordinates)
             unitVectorLocal(:,1) = unitVectorGlobal(:,2)
             unitVectorLocal(:,2) = unitVectorGlobal(:,1)
          endif

       endif   ! at the pole or not

       ! compute matrix for transforming from local to global coordinates
       transLocalToGlobal(:,:,iElement) = matmul(transpose(unitVectorGlobal), unitVectorLocal)

       ! compute matrix for transforming from local to global coordinates
       transGlobalToLocal(:,:,iElement) = matmul(transpose(unitVectorLocal), unitVectorGlobal)

    enddo   ! iElement

  end subroutine define_local_to_global_transformations

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine get_geometry_incremental_remap
!
!> \brief MPAS-CICE find various mesh-related quantities for incremental remapping
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  This routine finds several mesh quantities needed for incremental remapping:
!>  (1) cellsOnEdgeRemap, which is like cellsOnEdge but includes additional cells
!>  (2) edgesOnEdgeRemap, which is like edgesOnEdge but includes fewer edges
!>  (3) xVertexOnEdge, yVertexOnEdge, the local x/y coordinates of the vertices
!>      of each edge included in edgesOnEdgeRemap
!>  (4) minLengthEdgesOnVertex, which is used to detect potential CFL violations
!>  It also checks that cell 1 is in the left half-plane (and cell 2 in the right
!>   half-plane) for each edge.
!>  See the diagrams below for labeling conventions.
!-----------------------------------------------------------------------

  subroutine get_geometry_incremental_remap(&
       nEdges,                    &
       nVertices,                 &
       nCells,                    &
       nCellsSolve,               &
       nEdgesOnCell,              &
       maxEdges,                  &
       xCell,   yCell,   zCell,   &
       xVertex, yVertex, zVertex, &
       xEdge,   yEdge,   zEdge,   &
       indexToEdgeID,             &
       vertexDegree,              &
       cellsOnEdge,               &
       verticesOnEdge,            &
       verticesOnCell,            &
       edgesOnCell,               &
       edgesOnVertex,             &
       cellsOnVertex,             &
       xVertexOnCell, yVertexOnCell, &
       on_a_sphere,               &
       transGlobalToEdge,         &
       transGlobalToVertex,       &
       remapEdge,                 &
       cellsOnEdgeRemap,          &
       edgesOnEdgeRemap,          &
       xVertexOnEdge, yVertexOnEdge, &
       minLengthEdgesOnVertex)

       !--------------------------------------------------------------------------------------------------------!
       !                                                                                                        !
       ! The IR edge-numbering and cell-numbering conventions are as follows for a hex mesh (vertexDegree = 3): !
       ! The main edge (E0) has vertices V1 and V2 and is shared by cells C1 and C2.                            !
       ! C1 lies in the left half-plane of the vector pointing from C1 to C2, and C2 in the right half-plane.   !
       ! Edges E1 and E2 lie on C1, and edges E3 and E4 lie on C2.                                              !
       ! C3 is the cell sharing V1 with the main edge, and C4 is the cell sharing V2 with the main edge.        !
       ! Vertices 3-6 lie on the distant ends of edges 1-4, respectively.                                       !
       !                                                                                                        !
       !                       V3                       V4                                                      !
       !                        \                       /                                                       !
       !                         \                     /                                                        !
       !                          \E1      C1       E2/                                                         !
       !                           \                 /                                                          !
       !                    C3      \V1___________V2/      C4                                                   !
       !                            /      E0       \                                                           !
       !                           /                 \                                                          !
       !                          /E3      C2       E4\                                                         !
       !                         /                     \                                                        !
       !                        /                       \                                                       !
       !                       V5                       V6                                                      !
       !                                                                                                        !
       ! On a quad mesh (vertexDegree = 4), we add the 2 edges that are colinear with the main edge             !
       ! (E5 containing V1/V7 and E6 containing V2/V8).                                                         !
       ! We also add two cells that share a vertex with the main edge (C5 sharing V1 and C6 sharing V2).        !
       !                                                                                                        !
       !                             V3            V4                                                           !
       !                             |             |                                                            !
       !                             |             |                                                            !
       !                    C3       |E1    C1     |E2     C4                                                   !
       !                             |             |                                                            !
       !                V7 _________ |V1___________|V2_________V8                                               !
       !                       E5    |      E0     |     E6                                                     !
       !                             |             |                                                            !
       !                    C5       |E3    C2     |E4     C6                                                   !
       !                             |             |                                                            !
       !                             |             |                                                            !
       !                             V5            V6                                                           !
       !                                                                                                        !
       !--------------------------------------------------------------------------------------------------------!

    integer, intent(in) :: &
         nEdges,         & !< Input: number of edges
         nVertices,      & !< Input: number of vertices
         nCells,         & !< Input: number of locally owned cells
         nCellsSolve,    & !< Input: number of cells
         maxEdges,       & !< Input: max number of edges (and vertices) per cell
         vertexDegree      !< Input: number of edges that meet at each vertex

    integer, dimension(:), intent(in) ::  &
         nEdgesOnCell      !< Input: number of edges of each cell

    real(kind=RKIND), dimension(:), intent(in) ::  &
         xCell,   yCell,   zCell,    & !< Input: global coordinates of cells
         xVertex, yVertex, zVertex,  & !< Input: global coordinates of vertices
         xEdge,   yEdge,   zEdge       !< Input: global coordinates of edges

    integer, dimension(:), intent(in) ::  &
         indexToEdgeID     !< Input: global index for each edge on a block (diagnostic only)

    integer, dimension(:,:), intent(in) :: &
         cellsOnEdge,    & !< Input: index for each cell sharing a given edge
         verticesOnEdge, & !< Input: index for each vertex of a given edge
         verticesOnCell, & !< Input: index for each vertex of a given cell
         edgesOnVertex,  & !< Input: index for each edge neighbor of a given vertex
         edgesOnCell,    & !< Input: index for each edge of a given cell
         cellsOnVertex     !< Input: cell index for each cell of a given vertex

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnCell,  & !< Output: x (east) coordinate of vertex relative to cell center in local tangent plane
         yVertexOnCell     !< Output: y (north) coordinate of vertex relative to cell center in local tangent plane

    logical, intent(in) ::  &
         on_a_sphere       !< Input: if true, then the mesh lives on the surface of a sphere; else in a plane

    real(kind=RKIND), dimension(3, 3, nEdges), intent(in) ::  &
         transGlobalToEdge    !< Input:  3 x 3 matrix for transforming from global to edge-based coordinates

    real(kind=RKIND), dimension(3, 3, nEdges), intent(in) ::  &
         transGlobalToVertex  !< Input:  3 x 3 matrix for transforming from global to vertex-based coordinates

    !Note: remapEdge could be a logical array, but logical fields are not supported (as of June 2015)
    integer, dimension(:), intent(out) :: &
         remapEdge          !< Output: = 1 if IR fluxes need to be computed across an edge, else = 0
                            !         (i.e., = 1 for all edges of locally owned cells)

    integer, dimension(:,:), intent(out) ::  &
         cellsOnEdgeRemap, &!< Output: index for each cell neighbor of a given edge,
                            !          including cells that share a vertex with the edge
         edgesOnEdgeRemap   !< Output: index for nearest edge neighbors of a given edge

    real(kind=RKIND), dimension(:,:), intent(out) ::  &
         xVertexOnEdge,  &  !< Output: x (east) coordinate of vertex relative to edge midpoint in local tangent plane
         yVertexOnEdge      !< Output: y (north) coordinate of vertex relative to edge midpoint in local tangent plane

    real(kind=RKIND), dimension(:), intent(out) :: &
         minLengthEdgesOnVertex     !< Output: minimum length of the edges on each vertex

    ! local variables

    integer :: iEdge, iVertex, iCell, iMainVertex, iSideVertex, iVertex1, iVertex2
    integer :: iVertexOnEdge, iCellOnVertex, iVertexOnCell, iEdgeOnEdge, iEdgeOnCell, iCellOnEdge, iEdgeOnVertex
    integer :: iMainEdgeOnCell
    integer :: iCellNeighbor, iEdgeNeighbor

    real(kind=RKIND), dimension(2) :: cellCenter   ! x/y coordinates of neighboring cell center
    real(kind=RKIND), dimension(3) :: workVector   ! x/y/z coordinates of a global vector

    real(kind=RKIND) :: xVector, yVector, edgeLength

    real(kind=RKIND), dimension(2,2) :: &
         edgeVertex    ! x/y coordinates of each of 2 vertices on an edge, relative to the edge midpoint

    integer :: count, m, m1, m2, n, n1

    integer :: nEdgesOnEdgeRemap

    logical :: newRemapEdge  ! true if this edge is to be added to edgesOnEdgeRemap

    character(len=strKIND) :: &
         errorMessage ! error message for abort

    ! Identify edges of locally owned cells.  These are the edges that require IR fluxes.

    remapEdge(:) = 0
    do iCell = 1, nCellsSolve
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(iEdgeOnCell,iCell)
          remapEdge(iEdge) = 1
       enddo
    enddo   ! iCell

    ! Exclude edges that do not have cell neighbors on each side.
    ! This amounts to a no-flux condition at any global boundary
    ! (unless halo cells exist beyond the boundary).

    do iEdge = 1, nEdges
       if (remapEdge(iEdge) == 1) then
          do iCellOnEdge = 1, 2
             iCell = cellsOnEdge(iCellOnEdge,iEdge)
             if (iCell < 1 .or. iCell > nCells) then
                remapEdge(iEdge) = 0
!!                write(stderrUnit,*) 'Excluding edge, iedge (local, global):', iEdge, indexToEdgeID(iEdge)
!!                write(stderrUnit,*) 'Cell neighbors:', cellsOnEdge(1,iEdge), cellsOnEdge(2,iEdge)
                exit
             endif   ! cell neighbor does not exist
          enddo   ! iCellOnEdge
       endif
    enddo

    if (verboseInit) then
!!       write(stderrUnit,*) 'Edges not used for IR:'
       do iEdge = 1, nEdges
          if (remapEdge(iEdge) == 0) then
!!             write(stderrUnit,*) iEdge
          endif
       enddo
    endif

    ! Make sure cell 1 is in the left half-plane (implying that cell 2 is in the right half-plane)
    !  for each edge.  This should be true for MPAS-conforming meshes.
    ! If this check fails, then the user will need to apply the mpas mesh converter to the input grid file.

    if (on_a_sphere) then

       ! loop over edges
       do iEdge = 1, nEdges
          if (remapEdge(iEdge) == 1) then

             do iVertexOnEdge = 1, 2

                ! Compute coordinates of this vertex (relative to the edge) in global coordinates
                iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

                workVector(1) = xVertex(iVertex) - xEdge(iEdge)
                workVector(2) = yVertex(iVertex) - yEdge(iEdge)
                workVector(3) = zVertex(iVertex) - zEdge(iEdge)

                ! transform to coordinates of iEdge
                edgeVertex(1,iVertexOnEdge) = transGlobalToEdge(1,1,iEdge) * workVector(1)  &
                                            + transGlobalToEdge(1,2,iEdge) * workVector(2)  &
                                            + transGlobalToEdge(1,3,iEdge) * workVector(3)

                edgeVertex(2,iVertexOnEdge) = transGlobalToEdge(2,1,iEdge) * workVector(1)  &
                                            + transGlobalToEdge(2,2,iEdge) * workVector(2)  &
                                            + transGlobalToEdge(2,3,iEdge) * workVector(3)

             enddo   ! iVertexOnEdge

             ! compute global coordinates of the center of neighboring cell 1

             iCell = cellsOnEdge(1,iEdge)

             if (iCell >= 1 .and. iCell <= nCells) then  ! neighbor cell exists for this edge

                workVector(1) = xCell(iCell) - xEdge(iEdge)
                workVector(2) = yCell(iCell) - yEdge(iEdge)
                workVector(3) = zCell(iCell) - zEdge(iEdge)

                ! transform to coordinates of iEdge
                cellCenter(1) = transGlobalToEdge(1,1,iEdge) * workVector(1)  &
                              + transGlobalToEdge(1,2,iEdge) * workVector(2)  &
                              + transGlobalToEdge(1,3,iEdge) * workVector(3)

                cellCenter(2) = transGlobalToEdge(2,1,iEdge) * workVector(1)  &
                              + transGlobalToEdge(2,2,iEdge) * workVector(2)  &
                              + transGlobalToEdge(2,3,iEdge) * workVector(3)

                ! find whether this cell center is in the left half-plane of the edge
                if (point_in_half_plane(edgeVertex(:,1), edgeVertex(:,2), cellCenter(:))) then
                   ! cell 1 is in the left half-plane, as required
                else
                   write(errorMessage,*) 'cell 1 is not in the left half-plane: iEdge (local, global ID) =', &
                        iEdge, indexToEdgeID(iEdge)
                   write(stderrUnit,*) 'ERROR: '//trim(errorMessage)
                   write(stderrUnit,*) 'edgeVertex(:,1):', edgeVertex(:,1)
                   write(stderrUnit,*) 'edgeVertex(:,2):', edgeVertex(:,2)
                   write(stderrUnit,*) 'cellCenter(:):', cellCenter(:)
                   call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
                endif

             endif  ! neighbor cell exists

          endif  ! remapEdge = 1
       enddo  ! iEdge

    else   ! on a plane

       do iEdge = 1, nEdges
          if (remapEdge(iEdge) == 1) then

             do iVertexOnEdge = 1, 2

                iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

                edgeVertex(1,iVertexOnEdge) = xVertex(iVertex) - xEdge(iEdge)
                edgeVertex(2,iVertexOnEdge) = yVertex(iVertex) - yEdge(iEdge)

             enddo   ! iVertexOnEdge

             iCell = cellsOnEdge(1,iEdge)

             if (iCell >= 1 .and. iCell <= nCells) then  ! neighbor cell exists for this edge

                cellCenter(1) = xCell(iCell) - xEdge(iEdge)
                cellCenter(2) = yCell(iCell) - yEdge(iEdge)

                ! find whether this cell center is in the left half plane of the edge
                if (point_in_half_plane(edgeVertex(:,1), edgeVertex(:,2), cellCenter(:))) then
                   ! cell 1 is in the left half-plane, as required
                else
                   write(errorMessage,*) 'cell 1 is not in the left half-plane: iEdge (local, global ID) =', &
                        iEdge, indexToEdgeID(iEdge)
                   write(stderrUnit,*) 'ERROR: '//trim(errorMessage)
                   write(stderrUnit,*) 'edgeVertex(:,1):', edgeVertex(:,1)
                   write(stderrUnit,*) 'edgeVertex(:,2):', edgeVertex(:,2)
                   write(stderrUnit,*) 'cellCenter(:):', cellCenter(:)
                   call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
                endif

             endif  ! neighbor cell exists

          endif   ! remapEdge = 1
       enddo   ! iEdge

    endif   ! on a sphere

    ! Compute indices of each edge neighbor of each edge: i.e., each edge that shares a vertex
    !  with a given edge (4 edge neighbors for hexes, 6 for quads).
    ! Also compute indices of each cell neighbor of each edge: i.e., each cell that includes a vertex
    !  of a given edge (4 cell neighbors for hexes, 6 for quads).
    ! Note: This can get tricky at a global boundary. If a cell or edge at a given orientation
    !        to the main edge does not exist, it receives index 0, nCells+1 or nEdges+1.
    !       I have tried to write this in a way that we find all the cells and edges that exist,
    !        even if some neighbor cells and edges do not exist.

    edgesOnEdgeRemap(:,:) = 0
    cellsOnEdgeRemap(:,:) = 0

    ! loop over edges
    do iEdge = 1, nEdges
       if (remapEdge(iEdge) == 1) then  ! this edge is a boundary of a locally owned cell
                                        ! (but might be on the global boundary, hence missing some neighbors)

          ! identify the 2 cells that share the edge
          do iCellOnEdge = 1, 2
             iCell = cellsOnEdge(iCellOnEdge,iEdge)
             cellsOnEdgeRemap(iCellOnEdge,iEdge) = iCell   ! C1, then C2
          enddo

          ! Identify the side edges (4 for a hex mesh, 6 for a quad mesh) that share a single vertex with iEdge
          ! Note: I first coded a method looping over edgesOnVertex for vertices 1 and 2.
          !       However, that method fails if there are non-existent edges at global boundaries that are
          !        pushed to the end of the array, as can happen using the MPAS mesh converter tool (as of Aug. 2015).
          !       The following method (using edgesOnCell) is less elegant, but it finds all the remapping edges
          !        in the required order as shown in the diagrams above.

          ! loop over the 2 cells on the edge
          do iCellOnEdge = 1, 2

             iCell = cellsOnEdge(iCellOnEdge,iEdge)

             if (iCell >= 1 .and. iCell <= nCells) then  ! this cell exists

                ! Of the several edges of cell 1, identify the one that is iEdge

                do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                   if (edgesOnCell(iEdgeOnCell,iCell) == iEdge) then
                      iMainEdgeOnCell = iEdgeOnCell
                      exit
                   endif
                enddo

                ! Having found the index of iEdgeOnCell corresponding to iEdge, pick up the
                !  neighboring edges on iCell that share a vertex with iEdge (E1 and E2 for C1,
                !  (E3 and E4 for C2).

                if (iCellOnEdge == 1) then

                   ! find E1
                   iEdgeOnCell = iMainEdgeOnCell - 1
                   if (iEdgeOnCell < 1) iEdgeOnCell = iEdgeOnCell + nEdgesOnCell(iCell)
                   edgesOnEdgeRemap(1,iEdge) = edgesOnCell(iEdgeOnCell,iCell)

                   ! find E2
                   iEdgeOnCell = iMainEdgeOnCell + 1
                   if (iEdgeOnCell > nEdgesOnCell(iCell)) iEdgeOnCell = iEdgeOnCell - nEdgesOnCell(iCell)
                   edgesOnEdgeRemap(2,iEdge) = edgesOnCell(iEdgeOnCell,iCell)

                elseif (iCellOnEdge == 2) then

                   ! find E3
                   iEdgeOnCell = iMainEdgeOnCell + 1
                   if (iEdgeOnCell > nEdgesOnCell(iCell)) iEdgeOnCell = iEdgeOnCell - nEdgesOnCell(iCell)
                   edgesOnEdgeRemap(3,iEdge) = edgesOnCell(iEdgeOnCell,iCell)

                   ! find E4
                   iEdgeOnCell = iMainEdgeOnCell - 1
                   if (iEdgeOnCell < 1) iEdgeOnCell = iEdgeOnCell + nEdgesOnCell(iCell)
                   edgesOnEdgeRemap(4,iEdge) = edgesOnCell(iEdgeOnCell,iCell)

                endif  ! iCellOnEdge

             endif  ! this cell exists

          enddo  ! iCellOnEdge

          ! For quad meshes, we still need to pick up E5 and E6, which share a single vertex with iEdge and cells 1 and 2

          if (vertexDegree == 4) then  ! quad mesh

             do iVertexOnEdge = 1, 2

                iVertex = verticesOnEdge(iVertexOnEdge,iEdge)  ! V1 or V2

                ! Of the 4 edges meeting at this vertex, find the one that has not yet been picked up
                do iEdgeOnVertex = 1, vertexDegree  ! loop through candidate edges for E5/E6
                   iEdgeNeighbor = edgesOnVertex(iEdgeOnVertex,iVertex)
                   if (iEdgeNeighbor >= 1 .and. iEdgeNeighbor <= nEdges) then  ! this is a real edge
                      newRemapEdge = .true.  ! true until found to be false
                      do iEdgeOnEdge = 1, 4  ! loop over the 4 remap edges found so far
                         if (iEdgeNeighbor == edgesOnEdgeRemap(iEdgeOnEdge,iEdge) .or. &
                              iEdgeNeighbor == iEdge) then  ! already picked up this edge
                            newRemapEdge = .false.
                            exit   ! exit loop over iEdgeOnEdgeRemap
                         endif
                      enddo   ! iEdgeOnEdgeRemap
                      ! If we get this far and newRemapEdge = true, then we have found the missing edge (E5 for V1, E6 for V2)
                      if (newRemapEdge) then
                         edgesOnEdgeRemap(iVertexOnEdge+4,iEdge) = iEdgeNeighbor
                         exit   ! exit loop over iEdgeOnVertex
                      endif
                   endif   ! this is a real edge
                enddo

             enddo   ! iVertexOnEdge

          endif   ! vertexDegree = 4 (quad mesh)

          ! At this point we have C1 and C2, along with all the side edges
          ! Finally, find C3, C4 and (on quad mesh) C5 and C6
          ! V1 lies on C3 (and C5); V2 lies on C4 (and C6)

          if (vertexDegree == 3) then  ! hex mesh

             ! initialize cellsOnEdgeRemap to nCells+1 (in case no cell is found below)
             cellsOnEdgeRemap(3:4,iEdge) = nCells + 1

             do iVertexOnEdge = 1, 2

                ! if E1 (or E2) exists, find the other cell (apart from C1) that lies on this edge
                ! else if E3 (or E4) exists, find the other cell (apart from C2) that lies on this edge

                iEdgeNeighbor = edgesOnEdgeRemap(iVertexOnEdge,iEdge)  ! E1 or E2
                if (iEdgeNeighbor < 1 .or. iEdgeNeighbor > nEdges) then  ! edges does not exist; choose E3 or E4 instead
                   iEdgeNeighbor = edgesOnEdgeRemap(iVertexOnEdge+2,iEdge)  ! E3 or E4
                endif

                ! look for a cell that lies on this edge and is not C1 or C2
                do iCellOnEdge = 1, 2
                   iCellNeighbor = cellsOnEdge(iCellOnEdge,iEdgeNeighbor)
                   if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then   ! this cell exists
                      if (iCellNeighbor /= cellsOnEdgeRemap(1,iEdge) .and. &
                          iCellNeighbor /= cellsOnEdgeRemap(2,iEdge)) then ! not C1 or C2
                         cellsOnEdgeRemap(iVertexOnEdge+2,iEdge) = iCellNeighbor
                      endif
                   endif
                enddo

             enddo  ! iVertexOnEdge

          elseif (vertexDegree == 4) then   ! quad mesh

             ! initialize cellsOnEdgeRemap to nCells+1 (in case no cell is found below)
             cellsOnEdgeRemap(3:6,iEdge) = nCells + 1

             do iVertexOnEdge = 1, 2

                ! if E1 (or E2) exists, find the other cell (apart from C1) that lies on this edge
                iEdgeNeighbor = edgesOnEdgeRemap(iVertexOnEdge,iEdge)  ! E1 or E2
                if (iEdgeNeighbor >= 1 .and. iEdgeNeighbor <= nEdges) then  ! E1 (or E2) exists
                   ! look for a cell that lies on this edge and is not C1
                   do iCellOnEdge = 1, 2
                      iCellNeighbor = cellsOnEdge(iCellOnEdge,iEdgeNeighbor)
                      if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then ! this cell exists
                         if (iCellNeighbor /= cellsOnEdgeRemap(1,iEdge)) then    ! not C1
                            cellsOnEdgeRemap(iVertexOnEdge+2,iEdge) = iCellNeighbor  ! C3 or C4
                         endif
                      endif
                   enddo  ! iCellOnEdge
                endif  ! E1 (or E2) exists

                ! if E3 (or E4) exists, find the other cell (apart from C2) that lies on this edge
                iEdgeNeighbor = edgesOnEdgeRemap(iVertexOnEdge+2,iEdge)  ! E3 or E4
                if (iEdgeNeighbor >= 1 .and. iEdgeNeighbor <= nEdges) then  ! E3 (or E4) exists
                   ! look for a cell that lies on this edge and is not C2
                   do iCellOnEdge = 1, 2
                      iCellNeighbor = cellsOnEdge(iCellOnEdge,iEdgeNeighbor)
                      if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then ! this cell exists
                         if (iCellNeighbor /= cellsOnEdgeRemap(2,iEdge)) then    ! not C2
                            cellsOnEdgeRemap(iVertexOnEdge+4,iEdge) = iCellNeighbor  ! C5 or C6
                         endif
                      endif

                   enddo  ! iCellOnEdge

                endif  ! E1 (or E2) exists

             enddo  ! iVertexOnEdge

          endif  ! vertexDegree

       endif  ! remapEdge = 1
    enddo  ! iEdge

    ! compute vertex coordinates relative to each edge

    if (vertexDegree == 3) then   ! hex mesh
       nEdgesOnEdgeRemap = 4
    elseif (vertexDegree == 4) then   ! quad mesh
       nEdgesOnEdgeRemap = 6
    endif

    if (on_a_sphere) then

       ! first loop: compute the local x/y coordinates of the 2 vertices of each edge
       ! Note: This loop is not limited to edges with remapEdge = 1, because
       !       some values in cells without remapEdge = 1 are needed in the second loop.
       ! Note: Previously, the global vectors (x/y/zVertex - x/y/zEdge) were transformed
       !       to the plane, one vertex at a time. This resulted in edge midpoints that
       !       did not lie exactly halfway between the 2 vertices. With the calculation
       !       below, the edge midpoints (by construction) exactly bisect the segment
       !       joining the 2 vertices in the plane.

       do iEdge = 1, nEdges

          ! Compute a global vector pointing from vertex 1 to vertex 2 of the edge

          iVertex1 = verticesOnEdge(1,iEdge)
          iVertex2 = verticesOnEdge(2,iEdge)

          workVector(1) = xVertex(iVertex2) - xVertex(iVertex1)
          workVector(2) = yVertex(iVertex2) - yVertex(iVertex1)
          workVector(3) = zVertex(iVertex2) - zVertex(iVertex1)

          ! Transform this vector from global to edge-based coordinates
          ! (Ignore the small component in the z direction)

          xVector = transGlobalToEdge(1,1,iEdge) * workVector(1)  &
                  + transGlobalToEdge(1,2,iEdge) * workVector(2)  &
                  + transGlobalToEdge(1,3,iEdge) * workVector(3)

          yVector = transGlobalToEdge(2,1,iEdge) * workVector(1)  &
                  + transGlobalToEdge(2,2,iEdge) * workVector(2)  &
                  + transGlobalToEdge(2,3,iEdge) * workVector(3)

          ! Compute xVertexToEdge and yVertexToEdge such that the edge midpoint lies at (0,0)
          !  and is exactly halfway between the 2 vertices.

          xVertexOnEdge(1,iEdge) = -0.5_RKIND * xVector
          yVertexOnEdge(1,iEdge) = -0.5_RKIND * yVector

          xVertexOnEdge(2,iEdge) =  0.5_RKIND * xVector
          yVertexOnEdge(2,iEdge) =  0.5_RKIND * yVector

       enddo   ! iEdge

       ! Second loop: For each edge, compute the coordinates of the more distant vertices of the neighboring edges
       ! (i.e., those edges that share a single vertex with iEdge).
       ! This is done in a separate loop because we make use of the nearest-vertex coordinates computed above.
       ! That is, the distance from an edge midpoint to a more distant vertex (V3, V4, ... V8)
       !  is computed as the sum of (1) the distance to the nearby vertex (V1 or V2) and
       !  (2) the distance from the nearby vertex to the distant vertex, where distance (2)
       !  is taken in the reference frame of the neighbor edge.
       ! Note: There are other ways to compute these coordinates, but this particular method ensures
       !       consistency in the angles between departure trajectories and cell edges, when viewed from
       !       the point of view of different edge midpoints. Without this consistency, it is possible
       !       for a departure triangle to have a greater area upon leaving a cell than it has upon entry,
       !       which can result in negative masses.

       do iEdge = 1, nEdges
          if (remapEdge(iEdge) == 1) then

             do iEdgeOnEdge = 1, nEdgesOnEdgeRemap   ! 4 neighbor edges for hexes, 6 for quads

                iEdgeNeighbor = edgesOnEdgeRemap(iEdgeOnEdge,iEdge)

                if (iEdgeNeighbor >= 1 .and. iEdgeNeighbor <= nEdges) then  ! neighbor edge exists

                   ! identify which vertex of the neighbor edge is shared with iEdge; the other is a new vertex
                   do n = 1, 2      ! vertices of iEdge
                      do m = 1, 2   ! vertices of iEdgeNeighbor

                         if (verticesOnEdge(m,iEdgeNeighbor) == verticesOnEdge(n,iEdge)) then  ! shared vertex

                            ! save the index of the shared vertex, relative to iEdge
                            n1 = n

                            ! save the index of the shared and distant vertex, relative to the neighbor edge
                            if (m == 1) then
                               m1 = 1
                               m2 = 2
                            else ! m = 2
                               m1 = 2
                               m2 = 1
                            endif

                            iMainVertex = verticesOnEdge(n1,iEdge)
                            iSideVertex = verticesOnEdge(m2,iEdgeNeighbor)
                            exit

                         endif   ! shared vertex
                      enddo   ! vertices of iEdgeNeighbor
                   enddo   ! vertices of iEdge

                   ! Compute x/yVertexOnEdge for the more distant vertex

                   iVertexOnEdge = 2 + iEdgeOnEdge  ! vertex index = edge index + 2

                   xVertexOnEdge(iVertexOnEdge,iEdge) = xVertexOnEdge(n1,iEdge)  &
                                                     + (xVertexOnEdge(m2,iEdgeNeighbor) - xVertexOnEdge(m1,iEdgeNeighbor))
                   yVertexOnEdge(iVertexOnEdge,iEdge) = yVertexOnEdge(n1,iEdge)  &
                                                     + (yVertexOnEdge(m2,iEdgeNeighbor) - yVertexOnEdge(m1,iEdgeNeighbor))

                endif  ! neighbor edge exists

             enddo  ! iEdgeOnEdge

          endif   ! remapEdge
       enddo   ! iEdge

    else  ! on a plane

       do iEdge = 1, nEdges
          if (remapEdge(iEdge) == 1) then

             ! first, compute the local x/y coordinates of the 2 vertices of iEdge

             do iVertexOnEdge = 1, 2

                iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

                ! Compute a vector pointing from the edge to the vertex

                xVertexOnEdge(iVertexOnEdge,iEdge) = xVertex(iVertex) - xEdge(iEdge)
                yVertexOnEdge(iVertexOnEdge,iEdge) = yVertex(iVertex) - yEdge(iEdge)

             enddo   ! iVertexOnEdge

             ! next, compute the coordinates of the vertices of the neighboring edges
             ! (i.e., those that share a vertex with iEdge)

             count = 2   ! 2 vertices already computed above

             do iEdgeOnEdge = 1, nEdgesOnEdgeRemap   ! 4 neighbor edges for hexes, 6 for quads

                iEdgeNeighbor = edgesOnEdgeRemap(iEdgeOnEdge,iEdge)

                if (iEdgeNeighbor >= 1 .and. iEdgeNeighbor <= nEdges) then  ! neighbor edge exists

                   ! identify the side vertex of the neighbor edge
                   do n = 1, 2    ! vertices of iEdge
                      do m = 1, 2   ! vertices of iEdgeNeighbor

                         if (verticesOnEdge(m,iEdgeNeighbor) == verticesOnEdge(n,iEdge)) then  ! shared vertex

                            ! identify the other vertex of this neighbor edge
                            if (m == 1) then
                               m2 = 2  ! the other vertex of iEdgeNeighbor
                            else ! m = 2
                               m2 = 1
                            endif
                            iVertex = verticesOnEdge(m2,iEdgeNeighbor)
                            count = count + 1
                            exit

                         endif   ! shared vertex

                      enddo   ! vertices of iEdgeNeighbor
                   enddo   ! vertices of iEdge

                   ! Compute a global vector pointing from the midpoint of iEdge to this vertex
                   xVertexOnEdge(count,iEdge) = xVertex(iVertex) - xEdge(iEdge)
                   yVertexOnEdge(count,iEdge) = yVertex(iVertex) - yEdge(iEdge)

                endif   ! neighbor edge exists

             enddo  ! iEdgeOnEdge

          endif   ! remapEdge
       enddo   ! iEdge

    endif  ! on a sphere

    ! Find the minimum length of the edges that meet at each vertex.
    ! This length is used to detect potential CFL violations.

    ! initialize
    minLengthEdgesOnVertex(:) = huge(0.0_RKIND)

    ! loop over vertices
    do iVertex = 1, nVertices

       ! loop over the edges on this vertex (3 for hexes, 4 for quads)
       do iEdgeOnVertex = 1, vertexDegree

          iEdge = edgesOnVertex(iEdgeOnVertex,iVertex)

          if (iEdge >= 1 .and. iEdge <= nEdges) then  ! this edge exists

             ! compute the length of the edge
             ! Note: On a sphere, this length will not exactly equal the edge length
             !       projected onto a local tangent plane, but should be pretty close.

             iVertex1 = verticesOnEdge(1,iEdge)
             iVertex2 = verticesOnEdge(2,iEdge)

             workVector(1) = xVertex(iVertex2) - xVertex(iVertex1)
             workVector(2) = yVertex(iVertex2) - yVertex(iVertex1)
             workVector(3) = zVertex(iVertex2) - zVertex(iVertex1)

             edgeLength = sqrt(workVector(1)**2 + workVector(2)**2 + workVector(3)**2)

             ! adjust the minimum as appropriate
             if (edgeLength < minLengthEdgesOnVertex(iVertex)) &
                  minLengthEdgesOnVertex(iVertex) = edgeLength

!!             if (verboseInit .and. vtestOnProc .and. block % localBlockID == vtestBlockID) then
!!                write(stderrUnit,*) 'iVertex, iEdge, edgeLength, minLength:', &
!!                     iVertex, iEdge, edgeLength, minLengthEdgesOnVertex(iVertex)
!!             endif

          endif  ! edge exists

       enddo  ! iEdgeOnVertex

    enddo  ! iVertex

  end subroutine get_geometry_incremental_remap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine get_vertex_on_cell_coordinates
!
!> \brief MPAS-CICE get vertex coordinates relative to cell centers
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine computes the coordinates of each vertex of each cell,
!>  in the local tangent plane of the cell center.
!
!-----------------------------------------------------------------------

  subroutine get_vertex_on_cell_coordinates(&
       nCells,                       &
       nEdgesOnCell,                 &
       verticesOnCell,               &
       xVertex, yVertex, zVertex,    &
       xCell,   yCell,   zCell,      &
       on_a_sphere,                  &
       transGlobalToCell,            &
       xVertexOnCell, yVertexOnCell)

    integer, intent(in) ::  &
         nCells            !< Input: number of cells

    integer, dimension(:), intent(in) ::  &
         nEdgesOnCell      !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         verticesOnCell    !< Input: vertex index for each vertex of a given cell

    real(kind=RKIND), dimension(:), intent(in) ::  &
         xVertex, yVertex, zVertex,  & !< Input: global coordinates of vertices
         xCell,   yCell,   zCell       !< Input: global coordinates of edges

    logical, intent(in) :: &
         on_a_sphere        !< Input: if true, then the mesh lives on the surface of a sphere; else in a plane

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         transGlobalToCell  !< Input: 3x3 matrix for transforming vectors from global to cell-based coordinates

    real(kind=RKIND), dimension(:,:), intent(out) ::  &
         xVertexOnCell,  &  !< Output: x (east) coordinate of vertex relative to cell center in local tangent plane
         yVertexOnCell      !< Output: y (north) coordinate of vertex relative to cell center in local tangent plane

    ! local variables

    integer :: iCell, iVertex, iVertexOnCell, iVertexP1, iVertexP2, iVertexCell

    real(kind=RKIND), dimension(3) :: workVector

    real(kind=RKIND), dimension(2,2) :: &
         edgeVertex         ! x/y coordinates of each vertex on an edge

    real(kind=RKIND), dimension(2) :: &
         nextVertex         ! x/y coordinates of a vertex

    real(kind=RKIND), dimension(3) :: workVector1, workVector2, workVector3
    real(kind=RKIND) :: dotProduct, crossProduct

    ! compute vertex coordinates relative to cell center

    if (on_a_sphere) then

       do iCell = 1, nCells

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             iVertex = verticesOnCell(iVertexOnCell, iCell)

             ! Compute global vector pointing from cell center to vertex

             workVector(1) = xVertex(iVertex) - xCell(iCell)
             workVector(2) = yVertex(iVertex) - yCell(iCell)
             workVector(3) = zVertex(iVertex) - zCell(iCell)

             ! Transform this vector from global to local cell-based coordinates

             xVertexOnCell(iVertexOnCell,iCell) = transGlobalToCell(1,1,iCell) * workVector(1)  &
                                                + transGlobalToCell(1,2,iCell) * workVector(2)  &
                                                + transGlobalToCell(1,3,iCell) * workVector(3)

             yVertexOnCell(iVertexOnCell,iCell) = transGlobalToCell(2,1,iCell) * workVector(1)  &
                                                + transGlobalToCell(2,2,iCell) * workVector(2)  &
                                                + transGlobalToCell(2,3,iCell) * workVector(3)

             ! The z component would be computed as follows, but is thrown away.
             ! It is small compared to the x and y components.
!             workVector(3) = transGlobalToCell(3,1,iCell) * workVector(1)  &
!                           + transGlobalToCell(3,2,iCell) * workVector(2)  &
!                           + transGlobalToCell(3,3,iCell) * workVector(3)

!             if (iCell == ctest) then
!                write(stderrUnit,*) 'Local:', &
!                     iVertexOnCell, xVertexOnCell(iVertexOnCell,iCell), yVertexOnCell(iVertexOnCell,iCell), workVector(3)
!             endif

          enddo  ! VertexOnCell

       enddo     ! iCell

    else  ! on a plane

       do iCell = 1, nCells

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             iVertex = verticesOnCell(iVertexOnCell, iCell)

             xVertexOnCell(iVertexOnCell,iCell) = xVertex(iVertex) - xCell(iCell)
             yVertexOnCell(iVertexOnCell,iCell) = yVertex(iVertex) - yCell(iCell)

          enddo ! iVertexOnCell

       enddo ! iCell

    endif

    ! Bug check - Make sure each cell is convex with vertices ordered counterclockwise.
    ! This should be true for conforming MPAS meshes.

    ! First check for CCW ordering in 3D spherical coordinates

    workVector1(:) = 0.0_RKIND
    workVector2(:) = 0.0_RKIND
    workVector3(:) = 0.0_RKIND

    if (on_a_sphere) then

       do iCell = 1, nCells

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             ! set the next vertex
             iVertexP1 = iVertexOnCell + 1
             if (iVertexP1 > nEdgesOnCell(iCell)) iVertexP1 = 1

             ! vector from cell center to vertex
             iVertex = verticesOnCell(iVertexOnCell,iCell)
             workVector1(1) = xVertex(iVertex) - xCell(iCell)
             workVector1(2) = yVertex(iVertex) - yCell(iCell)
             workVector1(3) = zVertex(iVertex) - zCell(iCell)

             ! vector from cell center to the next vertex
             iVertex = verticesOnCell(iVertexP1,iCell)
             workVector2(1) = xVertex(iVertex) - xCell(iCell)
             workVector2(2) = yVertex(iVertex) - yCell(iCell)
             workVector2(3) = zVertex(iVertex) - zCell(iCell)

             ! cross product of vector1 and vector 2
             call cross_product_3d(workVector1, workVector2, workVector)

             ! vector from center of Earth to cell center
             workVector3(1) = xCell(iCell)
             workVector3(2) = yCell(iCell)
             workVector3(3) = zCell(iCell)

             ! dot product of the result with vector3

             call dot_product_3d(workVector, workVector3, dotProduct)

             if (dotProduct >= 0.0_RKIND) then
                continue
             else
!                write(stderrUnit,*) 'Violation: iCell, zCell, dotProduct:', iCell, zCell(iCell), dotProduct
                write(stderrUnit,*) 'ERROR: IR requires convex cells with vertices arranged CCW around the cell'
                write(stderrUnit,*) 'Violation for cell:', iCell
                write(stderrUnit,*) 'iVertexOnCell, xVertex, yVertex, zVertex:'
                do iVertexCell = 1, nEdgesOnCell(iCell)
                   iVertex = verticesOnCell(iVertexCell,iCell)
                   write(stderrUnit,*) iVertexCell, xVertex(iVertex), yVertex(iVertex), zVertex(iVertex)
                enddo
                call mpas_dmpar_global_abort('MPAS-seaice: IR requires convex cells with vertices arranged CCW around the cell')
             endif

          enddo   ! iVertexOnCell

       enddo   ! iCell

       if (verboseInit) write(stderrUnit,*) 'CCW ordering of cell vertices on a sphere is OK'

    endif  ! on a sphere

    ! Bug check: Verify CCW ordering on the plane. This makes sure the projection from the sphere to the plane is correct.

    do iCell = 1, nCells

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          ! set the next vertex
          iVertexP1 = iVertexOnCell + 1
          if (iVertexP1 > nEdgesOnCell(iCell)) iVertexP1 = 1

          ! vector from cell center to vertex
          workVector1(1) = xVertexOnCell(iVertexOnCell,iCell)
          workVector1(2) = yVertexOnCell(iVertexOnCell,iCell)

          ! vector from cell center to the next vertex
          workVector2(1) = xVertexOnCell(iVertexP1,iCell)
          workVector2(2) = yVertexOnCell(iVertexP1,iCell)

          ! cross product of vector1 and vector 2
          call cross_product_2d(workVector1(1:2), workVector2(1:2), crossProduct)

          if (crossProduct >= 0.0_RKIND) then
             continue
          else
!             write(stderrUnit,*) 'Violation: iCell, zCell, crossProduct:', iCell, zCell(iCell), crossProduct
             write(stderrUnit,*) 'ERROR: IR requires convex cells with vertices arranged CCW around the cell'
             write(stderrUnit,*) 'Violation for cell:', iCell
             write(stderrUnit,*) 'iVertexOnCell, xVertexOnCell, yVertexOnCell:'
             do iVertexCell = 1, nEdgesOnCell(iCell)
                write(stderrUnit,*) iVertexCell, xVertexOnCell(iVertexCell,iCell), yVertexOnCell(iVertexCell,iCell)
             enddo
             call mpas_dmpar_global_abort(&
                  'MPAS-seaice: IR requires convex cells with vertices arranged CCW around the cell: check')
          endif

       enddo   ! iVertexOnCell

    enddo   ! iCell

    if (verboseInit) write(stderrUnit,*) 'CCW ordering of cell vertices on local tangent plane is OK'

  end subroutine get_vertex_on_cell_coordinates

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine get_incremental_remap_geometry_pointers
!
!> \brief MPAS-CICE get pointers to geometric quantities used for IR
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine assigns pointers to various geometric quantities used
!>  for incremental remapping.
!
!-----------------------------------------------------------------------

  subroutine get_incremental_remap_geometry_pointers(geomAvgCell, block)

    type(geometric_avg_cell_type), pointer :: &
         geomAvgCell       ! derived type holding geometric averages

    type(block_type), intent(in) :: block

    ! local variables

    type (mpas_pool_type), pointer :: incrementalRemapPool

    ! get IR pool
    call MPAS_pool_get_subpool(block % structs, 'incremental_remap', incrementalRemapPool)

    ! get arrays and assign pointers

    call MPAS_pool_get_array(incrementalRemapPool, 'xAvgCell',     geomAvgCell % x)
    call MPAS_pool_get_array(incrementalRemapPool, 'yAvgCell',     geomAvgCell % y)
    call MPAS_pool_get_array(incrementalRemapPool, 'xxAvgCell',    geomAvgCell % xx)
    call MPAS_pool_get_array(incrementalRemapPool, 'xyAvgCell',    geomAvgCell % xy)
    call MPAS_pool_get_array(incrementalRemapPool, 'yyAvgCell',    geomAvgCell % yy)
    call MPAS_pool_get_array(incrementalRemapPool, 'xxxAvgCell',   geomAvgCell % xxx)
    call MPAS_pool_get_array(incrementalRemapPool, 'xxyAvgCell',   geomAvgCell % xxy)
    call MPAS_pool_get_array(incrementalRemapPool, 'xyyAvgCell',   geomAvgCell % xyy)
    call MPAS_pool_get_array(incrementalRemapPool, 'yyyAvgCell',   geomAvgCell % yyy)
    call MPAS_pool_get_array(incrementalRemapPool, 'xxxxAvgCell',  geomAvgCell % xxxx)
    call MPAS_pool_get_array(incrementalRemapPool, 'xxxyAvgCell',  geomAvgCell % xxxy)
    call MPAS_pool_get_array(incrementalRemapPool, 'xxyyAvgCell',  geomAvgCell % xxyy)
    call MPAS_pool_get_array(incrementalRemapPool, 'xyyyAvgCell',  geomAvgCell % xyyy)
    call MPAS_pool_get_array(incrementalRemapPool, 'yyyyAvgCell',  geomAvgCell % yyyy)

  end subroutine get_incremental_remap_geometry_pointers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine compute_geometric_cell_averages
!
!> \brief MPAS-CICE compute some geometric averages over each cell for IR
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine computes the average values of x, y, x^2, xy, y^2, etc.
!>  for each grid cell, as required for incremental remapping
!
!-----------------------------------------------------------------------

  subroutine compute_geometric_cell_averages(&
       xVertexOnCell, yVertexOnCell,  &
       nCells,                        &
       maxEdges,                      &
       edgesOnCell,                   &
       nEdgesOnCell,                  &
       dcEdge,                        &
       dvEdge,                        &
       geomAvgCell)

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnCell,  & !< Input: x (east) coordinate of vertex relative to cell center in local tangent plane
         yVertexOnCell     !< Input: y (north) coordinate of vertex relative to cell center in local tangent plane

    integer, intent(in) ::  &
         nCells,        &  !< Input: number of cells
         maxEdges          !< Input: maximum number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         edgesOnCell       !< Input: edge index for each edge of a given cell

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell      !< Input: number of edges for each cell

    real(kind=RKIND), dimension(:), intent(in) :: &
         dcEdge,        &  !< Input: distance between the 2 cell centers on each side of an edge
         dvEdge            !< Input: distance between the 2 vertices of an edge

    type(geometric_avg_cell_type), pointer :: &
         geomAvgCell       !< Output: derived type holding geometric averages

    ! local variables

    real(kind=RKIND), dimension(maxEdges) :: &
         fracEdgeArea      ! fractional area of the triangle joining the cell center to the 2 vertices of each edge
                           ! (relative to the total cell area)
    real(kind=RKIND) :: &
         sumArea           ! sum of edge areas for a cell

    real(kind=RKIND) :: &
         x1Vertex, y1Vertex, x2Vertex, y2Vertex, x3Vertex, y3Vertex  ! coordinates of triangle vertices

    real(kind=RKIND), dimension(6) ::  &
         xQP, yQP,      &  ! x and y coordinates of 6 quadrature points per triangle
         wQP               ! weight assigned to each quadrature point

    real(kind=RKIND) ::    &
         xavg,    yavg,    &    ! geometric averages for a given triangle
         xxavg,   xyavg,   yyavg,    &
         xxxavg,  xxyavg,  xyyavg,  yyyavg,  &
         xxxxavg, xxxyavg, xxyyavg, xyyyavg, yyyyavg

    integer :: n, iCell, iEdge, iEdgeOnCell, iVertexOnCell

    ! loop over cells
    do iCell = 1, nCells

       fracEdgeArea(:) = 0.0_RKIND
       sumArea = 0.0_RKIND

       ! compute the area of the triangle joining the cell center to the 2 vertices of each edge
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(iEdgeOnCell, iCell)
          fracEdgeArea(iEdgeOnCell) = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)
          sumArea = sumArea + fracEdgeArea(iEdgeOnCell)

       enddo

       ! convert this area to a fraction of the total cell area
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          fracEdgeArea(iEdgeOnCell) = fracEdgeArea(iEdgeOnCell) / sumArea

       enddo

       ! initialize geometric averages for this cell
       geomAvgCell % x(iCell)    = 0.0_RKIND
       geomAvgCell % y(iCell)    = 0.0_RKIND
       geomAvgCell % xx(iCell)   = 0.0_RKIND
       geomAvgCell % xy(iCell)   = 0.0_RKIND
       geomAvgCell % yy(iCell)   = 0.0_RKIND
       geomAvgCell % xxx(iCell)  = 0.0_RKIND
       geomAvgCell % xxy(iCell)  = 0.0_RKIND
       geomAvgCell % xyy(iCell)  = 0.0_RKIND
       geomAvgCell % yyy(iCell)  = 0.0_RKIND
       geomAvgCell % xxxx(iCell) = 0.0_RKIND
       geomAvgCell % xxxy(iCell) = 0.0_RKIND
       geomAvgCell % xxyy(iCell) = 0.0_RKIND
       geomAvgCell % xyyy(iCell) = 0.0_RKIND
       geomAvgCell % yyyy(iCell) = 0.0_RKIND

       ! loop over edges to compute the cell-average values of x, y, xx, xy, yy, etc.

       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          ! find the coordinates of the triangle vertices associated with this edge
          ! Note: Vertex indices lead edge indices, so edge n has vertices n and n+1.
          x1Vertex = 0.0_RKIND   ! cell center
          y1Vertex = 0.0_RKIND

          iVertexOnCell = iEdgeOnCell   ! one vertex of this edge
          x2Vertex = xVertexOnCell(iVertexOnCell, iCell)
          y2Vertex = yVertexOnCell(iVertexOnCell, iCell)

          iVertexOnCell = iEdgeOnCell+1  ! the other vertex of this edge
          if (iVertexOnCell > nEdgesOnCell(iCell)) iVertexOnCell = 1
          x3Vertex = xVertexOnCell(iVertexOnCell, iCell)
          y3Vertex = yVertexOnCell(iVertexOnCell, iCell)

          ! find the coordinates of 6 quadrature points
          ! Note: 6 points are sufficient for exact integration of 4th order polynomials

          xQP(1) = q1TriangleQP * x1Vertex + q1TriangleQP * x2Vertex + q2TriangleQP * x3Vertex
          yQP(1) = q1TriangleQP * y1Vertex + q1TriangleQP * y2Vertex + q2TriangleQP * y3Vertex

          xQP(2) = q1TriangleQP * x1Vertex + q2TriangleQP * x2Vertex + q1TriangleQP * x3Vertex
          yQP(2) = q1TriangleQP * y1Vertex + q2TriangleQP * y2Vertex + q1TriangleQP * y3Vertex

          xQP(3) = q2TriangleQP * x1Vertex + q1TriangleQP * x2Vertex + q1TriangleQP * x3Vertex
          yQP(3) = q2TriangleQP * y1Vertex + q1TriangleQP * y2Vertex + q1TriangleQP * y3Vertex

          xQP(4) = q3TriangleQP * x1Vertex + q4TriangleQP * x2Vertex + q4TriangleQP * x3Vertex
          yQP(4) = q3TriangleQP * y1Vertex + q4TriangleQP * y2Vertex + q4TriangleQP * y3Vertex

          xQP(5) = q4TriangleQP * x1Vertex + q3TriangleQP * x2Vertex + q4TriangleQP * x3Vertex
          yQP(5) = q4TriangleQP * y1Vertex + q3TriangleQP * y2Vertex + q4TriangleQP * y3Vertex

          xQP(6) = q4TriangleQP * x1Vertex + q4TriangleQP * x2Vertex + q3TriangleQP * x3Vertex
          yQP(6) = q4TriangleQP * y1Vertex + q4TriangleQP * y2Vertex + q3TriangleQP * y3Vertex

          ! set weight of each quadrature point
          wQP(1:3) = w1TriangleQP
          wQP(4:6) = w2TriangleQP

          if (verboseGeomAvg .and. ctestOnProc) then
             if (iCell == ctest) then
                write(stderrUnit,*) ' '
                write(stderrUnit,*) 'edge =', iEdgeOnCell
                write(stderrUnit,*) 'x1, y1:', x1Vertex, y1Vertex
                write(stderrUnit,*) 'x2, y2:', x2Vertex, y2Vertex
                write(stderrUnit,*) 'x3, y3:', x3Vertex, y3Vertex
                write(stderrUnit,*) 'quad point x, y, w:'
                do n = 1, 6
                   write(stderrUnit,*) n, xQP(n), yQP(n), wQP(n)
                enddo
             endif
          endif

          ! compute geometric averages for this triangle

          xavg = 0.0_RKIND
          yavg = 0.0_RKIND
          xxavg = 0.0_RKIND
          xyavg = 0.0_RKIND
          yyavg = 0.0_RKIND
          xxxavg = 0.0_RKIND
          xxyavg = 0.0_RKIND
          xyyavg = 0.0_RKIND
          yyyavg = 0.0_RKIND
          xxxxavg = 0.0_RKIND
          xxxyavg = 0.0_RKIND
          xxyyavg = 0.0_RKIND
          xyyyavg = 0.0_RKIND
          yyyyavg = 0.0_RKIND

          ! sum over quadrature points
          do n = 1, 6
             xavg    = xavg    + wQP(n) * xqp(n)
             yavg    = yavg    + wQP(n) * yqp(n)
             xxavg   = xxavg   + wQP(n) * xqp(n)**2
             xyavg   = xyavg   + wQP(n) * xqp(n) * yqp(n)
             yyavg   = yyavg   + wQP(n) * yqp(n)**2
             xxxavg  = xxxavg  + wQP(n) * xqp(n)**3
             xxyavg  = xxyavg  + wQP(n) * xqp(n)**2 * yqp(n)
             xyyavg  = xyyavg  + wQP(n) * xqp(n) * yqp(n)**2
             yyyavg  = yyyavg  + wQP(n) * yqp(n)**3
             xxxxavg = xxxxavg + wQP(n) * xqp(n)**4
             xxxyavg = xxxyavg + wQP(n) * xqp(n)**3 * yqp(n)
             xxyyavg = xxyyavg + wQP(n) * xqp(n)**2 * yqp(n)**2
             xyyyavg = xyyyavg + wQP(n) * xqp(n) * yqp(n)**3
             yyyyavg = yyyyavg + wQP(n) * yqp(n)**4
          enddo

          if (verboseGeomAvg .and. ctestOnProc) then
             if (iCell == ctest) then
                write(stderrUnit,*) ' '
                write(stderrUnit,*) 'Averages for edge ', iEdgeOnCell
                write(stderrUnit,*) 'x:', xavg
                write(stderrUnit,*) 'y:', yavg
                write(stderrUnit,*) 'xx:', xxavg
                write(stderrUnit,*) 'xy:', xyavg
                write(stderrUnit,*) 'yy:', yyavg
                write(stderrUnit,*) 'xxx:', xxxavg
                write(stderrUnit,*) 'xxy:', xxyavg
                write(stderrUnit,*) 'xyy:', xyyavg
                write(stderrUnit,*) 'yyy:', yyyavg
                write(stderrUnit,*) 'xxxx:', xxxxavg
                write(stderrUnit,*) 'xxxy:', xxxyavg
                write(stderrUnit,*) 'xxyy:', xxyyavg
                write(stderrUnit,*) 'xyyy:', xyyyavg
                write(stderrUnit,*) 'yyyy:', yyyyavg
             endif
          endif

          ! add this triangle contribution to the cumulative sum for the cell,
          ! weighted by the fractional area of the triangle
          geomAvgCell % x   (iCell) = geomAvgCell % x   (iCell) + fracEdgeArea(iEdgeOnCell) * xavg
          geomAvgCell % y   (iCell) = geomAvgCell % y   (iCell) + fracEdgeArea(iEdgeOnCell) * yavg
          geomAvgCell % xx  (iCell) = geomAvgCell % xx  (iCell) + fracEdgeArea(iEdgeOnCell) * xxavg
          geomAvgCell % xy  (iCell) = geomAvgCell % xy  (iCell) + fracEdgeArea(iEdgeOnCell) * xyavg
          geomAvgCell % yy  (iCell) = geomAvgCell % yy  (iCell) + fracEdgeArea(iEdgeOnCell) * yyavg
          geomAvgCell % xxx (iCell) = geomAvgCell % xxx (iCell) + fracEdgeArea(iEdgeOnCell) * xxxavg
          geomAvgCell % xxy (iCell) = geomAvgCell % xxy (iCell) + fracEdgeArea(iEdgeOnCell) * xxyavg
          geomAvgCell % xyy (iCell) = geomAvgCell % xyy (iCell) + fracEdgeArea(iEdgeOnCell) * xyyavg
          geomAvgCell % yyy (iCell) = geomAvgCell % yyy (iCell) + fracEdgeArea(iEdgeOnCell) * yyyavg
          geomAvgCell % xxxx(iCell) = geomAvgCell % xxxx(iCell) + fracEdgeArea(iEdgeOnCell) * xxxxavg
          geomAvgCell % xxxy(iCell) = geomAvgCell % xxxy(iCell) + fracEdgeArea(iEdgeOnCell) * xxxyavg
          geomAvgCell % xxyy(iCell) = geomAvgCell % xxyy(iCell) + fracEdgeArea(iEdgeOnCell) * xxyyavg
          geomAvgCell % xyyy(iCell) = geomAvgCell % xyyy(iCell) + fracEdgeArea(iEdgeOnCell) * xyyyavg
          geomAvgCell % yyyy(iCell) = geomAvgCell % yyyy(iCell) + fracEdgeArea(iEdgeOnCell) * yyyyavg

       enddo  ! iEdgeOnCell

    enddo     ! iCell

  end subroutine compute_geometric_cell_averages

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cice_run_advection_incremental_remap
!
!> \brief MPAS-CICE incremental remapping driver
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine loops over blocks, setting up tracer data structures
!>  and calling the main incremental remap subroutine for each block
!
!-----------------------------------------------------------------------

  subroutine cice_run_advection_incremental_remap(&
       domain,      &
       clock,       &
       timeLevelIn, &
       updateHaloInitIn,  &
       updateHaloFinalIn)

    use mpas_dmpar

    type (domain_type), intent(inout) :: &
         domain    !< Input/Output: domain info

    type (MPAS_Clock_type), intent(in) :: &
         clock     !< Input: clock

    integer, intent(in), optional :: &
         timeLevelIn    !< Input: time level of input fields (1 or 2)

    logical, intent(in), optional :: &
         updateHaloInitIn,   & !< Input: if true, do initial halo updates for velocity, mass and tracer fields
         updateHaloFinalIn     !< Input: if true, do final halo updates for mass and tracer fields

    ! Local variables

    integer :: timeLevel    ! time level of input fields (1 or 2)

    logical :: &
         updateHaloInit,   & ! if true, do initial halo updates for velocity, mass and tracer fields
         updateHaloFinal     ! if true, do final halo updates for mass and tracer fields

    type (mpas_pool_type), pointer :: velocitySolverPool
    type (mpas_pool_type), pointer :: incrementalRemapPool
    type (mpas_pool_type), pointer :: tracerPool

    type(block_type), pointer :: block

    type (dm_info), pointer :: dminfo

    type(tracer_type), pointer :: thisTracer

    real(kind=RKIND), dimension(:,:), pointer ::  &
         workCategoryCell        ! work array with dimension(nCategories,nCells)

    logical, pointer :: &
         configConservationCheck, & ! namelist configuration whether perform conservation check
         configMonotonicityCheck    ! namelist configuration whether perform monotonicity check

    ! dynamics time step
    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    ! assign pointers
    dminfo => domain % dminfo

    ! Set time level
    if (present(timeLevelIn)) then
       timeLevel = timeLevelIn
    else
       timeLevel = 1        ! old time by default
    endif

    ! Set logical variables for halo updates

    if (present(updateHaloInitIn)) then
       updateHaloInit = updateHaloInitIn
    else
       updateHaloInit = .true.    ! do initial updates to be safe
    endif

    if (present(updateHaloFinalIn)) then
       updateHaloFinal = updateHaloFinalIn
    else
       updateHaloFinal = .false.  ! another module can take care of halo updates later
    endif

    if (verboseRun) then
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'In cice_incremental remap: timeLevel =', timeLevel
    endif

    ! halo updates
    ! NOTE: In order to use incremental_remap_block to update mass and tracers
    !       correctly for all locally owned cells in a block, the input mass
    !       and tracer fields must be up to date in a halo of 2 or more layers.
    !TODO - Make sure there are at least two layers of halo cells.
    !       (Is this part of the field info?)

    if (updateHaloInit) then

       if (verboseRun) write(stderrUnit,*) 'Doing initial halo updates for velocity, mass and tracers'

       call mpas_timer_start("incr remap init halo")
       ! velocity

       call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
       call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')

       ! mass and tracers
       call cice_update_tracer_halo(tracersHead, domain, timeLevel)
       call mpas_timer_stop("incr remap init halo")

    endif  ! halo update

    ! Convert the ice and snow volume to thickness.
    ! The IR code assumes that all tracers except the first (fractional ice area,
    !  which is a mass-like field) are true tracers, and not the product of area*tracer.
    ! Note: The volume_to_thickness calls could go inside subroutine incremental_remap_block,
    !       saving some subpool calls. But the thickness_to_volume calls have to go after
    !       the monotonicity check, which is outside incremental_remap_block.
    !       For symmetry, I've put both the volume_to_thickness and thickness_to_volume calls
    !       outside incremental_remap_block.
    ! TODO: Switch iceVolume and snowVolume to iceThickness and snowThickness?
    !       This would require some code rewrites outside the IR.

    ! loop over blocks

    block => domain % blocklist

    do while (associated(block))

       if (verboseRun) write(stderrUnit,*) 'Convert volume to thickness'

       ! get tracer pool and work array
       call mpas_pool_get_subpool(block % structs, 'tracers', tracerPool)
       call mpas_pool_get_subpool(block % structs, 'incremental_remap', incrementalRemapPool)
       call mpas_pool_get_array(incrementalRemapPool, 'workCategoryCell', workCategoryCell)

       ! set pointers to tracer arrays
       call cice_set_tracer_array_pointers(tracersHead, block, timeLevel)

       ! loop over tracers
       thisTracer => tracersHead

       do while (associated(thisTracer))

          if (trim(thisTracer % tracerName) == 'iceAreaCategory') then

             ! copy the ice area (assumed to be 2D) to a scratch array
             ! Note: iceArea always precedes iceVolume and snowVolume in the tracer list

             workCategoryCell(:,:) = thisTracer % array2D(:,:)

          elseif (trim(thisTracer % tracerName) == 'iceVolumeCategory' .or. &
                  trim(thisTracer % tracerName) == 'snowVolumeCategory') then

             ! convert volume to thickness
             call volume_to_thickness(&
                  workCategoryCell(:,:),  &
                  thisTracer % array2D(:,:))

          endif

          thisTracer => thisTracer % next

       enddo   ! associated(thisTracer)

       block => block % next

    enddo  ! associated(block)

    ! Call the main IR algorithm for each block

    block => domain % blocklist

    call mpas_timer_start("incr remap blocks")
    do while (associated(block))

       if (verboseRun) write(stderrUnit,*) 'Call incremental_remap_block: block =', block % localBlockID

       ! get dynamics time step
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_array(velocitySolverPool, "dynamicsTimeStep", dynamicsTimeStep)

       ! do IR for this block
       call incremental_remap_block(&
            domain, &
            block, &
            dynamicsTimeStep, &
            tracersHead)

       block  => block % next

    enddo  ! associated(block)
    call mpas_timer_stop("incr remap blocks")

    ! Optional check for conservation of mass and mass*tracer
    ! Note: This check must be done outside the block loop because it requires global sums

    call MPAS_pool_get_config(domain % configs, "config_conservation_check", configConservationCheck)
    if (configConservationCheck) then

       if (verboseRun) write(stderrUnit,*) 'Check conservation'

       call mpas_timer_start("incr remap tracer cons check")
       call check_tracer_conservation(dminfo, tracersHead)
       call mpas_timer_stop("incr remap tracer cons check")

    endif

    ! Optional check for tracer monotonicity
    ! Note: This check should be done outside the block loop because it requires a halo update.
    !       (However, monotonicity should hold in most cases without the halo update.)

    call MPAS_pool_get_config(domain % configs, "config_monotonicity_check", configMonotonicityCheck)
    if (configMonotonicityCheck) then

       if (verboseRun) write(stderrUnit,*) 'Check monotonicity'

       call mpas_timer_start("incr remap tracer mono check")
       call check_tracer_monotonicity(domain, tracersHead)
       call mpas_timer_stop("incr remap tracer mono check")

    endif

    ! Convert the ice and snow thickness back to volume.
    ! Note: This conversion must be done after the monotonicity check, which assumes
    !       that iceVolumeCategory and snowVolumeCategory are actually thicknesses.

    ! loop over blocks

    block => domain % blocklist

    do while (associated(block))

       if (verboseRun) write(stderrUnit,*) 'Convert thickness to volume'

       ! get tracer pool and work array
       call mpas_pool_get_subpool(block % structs, 'tracers', tracerPool)
       call mpas_pool_get_subpool(block % structs, 'incremental_remap', incrementalRemapPool)
       call mpas_pool_get_array(incrementalRemapPool, 'workCategoryCell', workCategoryCell)

       ! set pointers to tracer arrays
       call cice_set_tracer_array_pointers(tracersHead, block, timeLevel)

       ! loop over tracers
       thisTracer => tracersHead

       do while (associated(thisTracer))

          if (trim(thisTracer % tracerName) == 'iceAreaCategory') then

             ! copy the ice area (assumed to be 2D) to a scratch array
             ! Note: iceArea always precedes iceVolume and snowVolume in the tracer list

             workCategoryCell(:,:) = thisTracer % array2D(:,:)

          elseif (trim(thisTracer % tracerName) == 'iceVolumeCategory' .or. &
                  trim(thisTracer % tracerName) == 'snowVolumeCategory') then

             ! convert thickness to volume
             call thickness_to_volume(&
                  workCategoryCell(:,:),  &
                  thisTracer % array2D(:,:))

          endif

          thisTracer => thisTracer % next

       enddo   ! associated(thisTracer)

       block  => block % next

    enddo  ! associated(block)

    ! final halo updates

    if (updateHaloFinal) then

       if (verboseRun) write(stderrUnit,*) 'Doing final halo updates for mass and tracers'

       ! mass and tracers
       call mpas_timer_start("incr remap final halo")
       call cice_update_tracer_halo(tracersHead, domain, timeLevel)
       call mpas_timer_stop("incr remap final halo")

    endif  ! halo update

    if (verboseRun) write(stderrUnit,*) 'Done in IR'

  end subroutine cice_run_advection_incremental_remap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine incremental_remap_block
!
!> \brief MPAS-CICE incremental remapping driver
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine carries out the incremental remapping algorithm for a block,
!>  returning updated mass and tracer fields
!
!-----------------------------------------------------------------------
! Some notes on incremental remapping:
!
! References:
!
! Dukowicz, J. K., and J. R. Baumgardner, 2000: Incremental
!  remapping as a transport/advection algorithm, J. Comput. Phys.,
!  160, 318-335.
!
! Lipscomb, W. H., and E. C. Hunke, 2004: Modeling sea ice
!  transport using incremental remapping, Mon. Wea. Rev., 132,
!  1341-1354.
!
! The remapping routine is designed to transport a generic mass-like
! field (in CICE, the ice fractional area) along with an arbitrary number
! of tracers in two dimensions.  The velocity components are assumed
! to lie at cell vertices and the transported scalars at cell centers.
! Incremental remapping has the following desirable properties:
!
! (1) Tracer monotonicity is preserved.  That is, no new local
!     extrema are produced in fields like ice thickness or internal
!     energy.
! (2) The reconstucted mass and tracer fields vary linearly in x and y.
!     This means that remapping is second-order accurate in space,
!     except where horizontal gradients are limited to preserve
!     monotonicity.
! (3) There are economies of scale.  Transporting a single field
!     is fairly expensive, but additional fields have a relatively
!     low marginal cost.
!
! The following generic conservation equations may be solved:
!
!            dm/dt = del*(u*m)             (0)
!       d(m*T1)/dt = del*(u*m*T1)          (1)
!    d(m*T1*T2)/dt = del*(u*m*T1*T2)       (2)
! d(m*T1*T2*T3)/dt = del*(u*m*T1*T2*T3)    (3)
!
! where d is a partial derivative, del is the 2D divergence operator,
! u is the horizontal velocity, m is the mass density field, and
! T1, T2, and T3 are tracers.
!
! We say that T1 has one parent (the mass field m), T2 has two parents,
! and T3 has three parents.
!
! In CICE, m corresponds to ice fractional concentration.
! The paradigm for T1 is ice or snow thickness, and the paradigm for T2
! (assuming constant density) is ice or snow enthalpy.
! If density is variable, then density is T2 and enthalpy is T3.
!
! Requirements for input fields::
! (1) The input mass and tracer fields must be up to date in two or more layers
!     of halo cells. Else the output mass and tracers will not be correct
!     for all locally owned cells in the block.
! (2) The velocity must be up to date for all vertices of locally owned cells,
!     plus the vertices in the first halo layer.
! (3) In the tracer linked list, all T1 tracers must appear before their T2 children,
!     and all T2 tracers before their T3 children.
!
!-----------------------------------------------------------------------

  subroutine incremental_remap_block(&
       domain,      &
       block,       &
       dt,          &
       tracersHead, &
       timeLevelIn)

    ! in/out arguments
    type(domain_type), intent(in) :: &
         domain       !< Input: domain

    type(block_type), intent(inout) :: &
         block        !< Input: block info

    real(kind=RKIND), intent(in) :: &
         dt           !< Input: time step

    type(tracer_type), pointer :: &
         tracersHead  !< Input/output: pointer to first element of linked list of tracers
                      ! The pointer stays attached to the first tracer, but all tracers are updated

    integer, intent(in), optional :: timeLevelIn

    ! local arguments

    ! pools
    type (mpas_pool_type), pointer :: tracerPool
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: velocityPool
    type (mpas_pool_type), pointer :: incrementalRemapPool
    type (mpas_pool_type), pointer :: configPool

    ! state variables
    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity,   &     ! velocity components at vertices
         vVelocity

    ! mesh quantities
    integer, pointer :: &
         nCells,          & ! number of cells on block
         nVertices,       & ! number of vertices on block
         nEdges,          & ! number of edges on block
         nCellsSolve,     & ! number of locally owned cells on block
         maxEdges,        & ! max number of edges per cell
         vertexDegree,    & ! number of edges that meet at each vertex
         nTriPerEdgeRemap   ! max number of departure triangles per edge

    logical, pointer ::  &
         on_a_sphere,                & ! if true, then the mesh lives on the surface of a sphere; else in a plane
         config_rotate_cartesian_grid  ! if true, then the North and South Poles are rotated to the equator

    ! mesh arrays

    integer, dimension(:), pointer :: &
         nEdgesOnCell       ! number of edges for each cell

    integer, dimension(:,:), pointer ::  &
         edgesOnCell,  &    ! edge index for each edge of a given cell
         cellsOnEdge,  &    ! cell index for each of two cells sharing a given edge
         cellsOnCell,  &    ! cell index for each edge neighbor of a given cell
         verticesOnCell,  & ! vertex index for each vertex of a given cell
         verticesOnEdge     ! vertex index for each vertex of a given edge

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell,                  &  ! cell area
         dcEdge,                    &  ! distance between the 2 cell centers on each side of an edge
         dvEdge,                    &  ! distance between the 2 vertices of an edge
         minLengthEdgesOnVertex        ! minimum length of the edges on each vertex

    real(kind=RKIND), dimension(:), pointer :: &
         xCell,   yCell,   zCell       ! global coordinates of cells; diagnostic only

    !Note: The indexTo*ID arrays are for diagnostics only
    integer, dimension(:), pointer ::  &
         indexToCellID,   & ! global index for each cell on a block
         indexToVertexID, & ! global index for each vertex on a block
         indexToEdgeID      ! global index for each edge on a block

    ! IR-specific arrays

    type(geometric_avg_cell_type), pointer :: &
         geomAvgCell        ! derived type holding geometric averages

    integer, dimension(:), pointer :: &
         remapEdge          ! = 1 if IR fluxes need to be computed across an edge, else = 0
                            ! (i.e., = 1 for all edges of locally owned cells)

    integer, dimension(:,:), pointer ::  &
         cellsOnEdgeRemap,& ! index for each cell neighbor of a given edge,
                            !  including cells that share a vertex with the edge
         edgesOnEdgeRemap   ! edge index for nearest edge neighbors of a given edge
                            ! (i.e., those sharing a vertex with the edge)

    integer, dimension(:,:), pointer ::  &
         iCellTriangle      ! index of cell containing a departure triangle

    real(kind=RKIND), dimension(:,:), pointer ::  &
         xVertexOnCell, yVertexOnCell,    & ! local x/y coordinates of vertices relative to cell center
                                            ! 1st dimension = maxEdges, 2nd dimension = nCells
         xVertexOnEdge, yVertexOnEdge       ! local x/y coordinates of vertices relative to edge midpoint
                                            ! 1st dimension = maxVerticesPerEdgeRemap, 2nd dimension = nEdges

    real(kind=RKIND), dimension(:,:), pointer ::  &
         departurePoint,                & ! x/y coordinates of departure points relative to vertices
                                          ! 1st dimension = 2, 2nd dimension = nVertices
         triangleArea                     ! area of departure triangle
                                          ! 1st dimension = nTriPerEdgeRemap, 2nd dimension = nEdges

    real(kind=RKIND), dimension(:,:,:), pointer ::  &
         coeffsReconstruct,    & ! coefficients for reconstructing the gradient at a cell center, given normal components on edges
         xTriangle, yTriangle    ! x and y coordinates of vertices and/or quadrature points of departure triangles for each edge
                                 ! 1st dimension = nQuadPoints, 2nd dimension = nTriPerEdgeRemap, 3rd dimension = nEdges

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         transGlobalToCell       ! 3x3 matrix for transforming vectors from global to cell-based coordinates

    integer, dimension(:), pointer ::  &
         maskEdge,             & ! mask array with dimension(nEdges)
                                 ! = 1 where departure area is nonzero, else = 0
         maskCell                ! mask array with dimension(nCells)
                                 ! = 1 where ice is present, else = 0

    integer, dimension(:,:), pointer ::  &
         maskCategoryCell        ! mask array with dimension(nCategories,nCells)

    real(kind=RKIND), dimension(:,:), pointer ::  &
         workCategoryCell        ! work array with dimension(nCategories,nCells)

    integer :: timeLevel

    logical, parameter :: zapSmallMass = .true.  ! if true, remove mass values (i.e., fractional ice area in CICE)
                                                 ! Note: The default threshold is 10^(-22), or eps11^2
                                                 !TODO - Turn off if this is handled elsewhere (e.g., in column package)

    type(tracer_type), pointer :: thisTracer
    integer :: n, m, iEdge, iCat, iLayer, iCell
    integer :: nCategories, nLayers, count

    logical, pointer :: &
         configConservationCheck, & ! namelist configuration whether perform conservation check
         configMonotonicityCheck, & ! namelist configuration whether perform monotonicity check
         configRecoverTracerMeansCheck       ! namelist configuration whether perform recover tracer means check

    if (present(timeLevelIn)) then
       timeLevel = timeLevelIn
    else
       timeLevel = 1
    endif

    ! get pools
    call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
    call mpas_pool_get_subpool(block % structs, 'velocity_solver', velocityPool)
    call mpas_pool_get_subpool(block % structs, 'tracers', tracerPool)
    call mpas_pool_get_subpool(block % structs, 'incremental_remap', incrementalRemapPool)

    ! get config options
    configPool => block % configs
    call mpas_pool_get_config(configPool, 'config_rotate_cartesian_grid', config_rotate_cartesian_grid)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! get velocity
    call mpas_pool_get_array(velocityPool, 'uVelocity', uVelocity)
    call mpas_pool_get_array(velocityPool, 'vVelocity', vVelocity)

    ! set pointers to tracer arrays
    call cice_set_tracer_array_pointers(tracersHead, block, timeLevel)

    ! get mesh dimensions
    call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
    call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
    call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
    call mpas_pool_get_dimension(meshPool, 'nTriPerEdgeRemap', nTriPerEdgeRemap)

    ! get some mesh arrays
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
    call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
    call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
    call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
    call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
    call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
    call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
    call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
    call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffsReconstruct)
    call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
    call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)
    call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)  ! x/y/zCell are diagnostic only
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'zCell', zCell)

    ! get IR geometry arrays
    allocate(geomAvgCell)
    call get_incremental_remap_geometry_pointers(geomAvgCell, block)

    call mpas_pool_get_array(incrementalRemapPool, 'xVertexOnCell', xVertexOnCell)
    call mpas_pool_get_array(incrementalRemapPool, 'yVertexOnCell', yVertexOnCell)
    call mpas_pool_get_array(incrementalRemapPool, 'xVertexOnEdge', xVertexOnEdge)
    call mpas_pool_get_array(incrementalRemapPool, 'yVertexOnEdge', yVertexOnEdge)
    call mpas_pool_get_array(incrementalRemapPool, 'departurePoint', departurePoint)
    call mpas_pool_get_array(incrementalRemapPool, 'remapEdge', remapEdge)
    call mpas_pool_get_array(incrementalRemapPool, 'cellsOnEdgeRemap', cellsOnEdgeRemap)
    call mpas_pool_get_array(incrementalRemapPool, 'edgesOnEdgeRemap', edgesOnEdgeRemap)
    call mpas_pool_get_array(incrementalRemapPool, 'xTriangle', xTriangle)
    call mpas_pool_get_array(incrementalRemapPool, 'yTriangle', yTriangle)
    call mpas_pool_get_array(incrementalRemapPool, 'iCellTriangle', iCellTriangle)
    call mpas_pool_get_array(incrementalRemapPool, 'triangleArea', triangleArea)
    call mpas_pool_get_array(incrementalRemapPool, 'transGlobalToCell',   transGlobalToCell)
    call mpas_pool_get_array(incrementalRemapPool, 'minLengthEdgesOnVertex', minLengthEdgesOnVertex)
    call mpas_pool_get_array(incrementalRemapPool, 'maskEdge', maskEdge)
    call mpas_pool_get_array(incrementalRemapPool, 'maskCell', maskCell)
    call mpas_pool_get_array(incrementalRemapPool, 'maskCategoryCell', maskCategoryCell)
    call mpas_pool_get_array(incrementalRemapPool, 'workCategoryCell', workCategoryCell)

    ! Optional initial diagnostics: Write a list of cells with ice
    ! If verboseGlobal = true, then this info is also written at the end of each time step.

    if (first_call) then
       first_call = .false.
       thisTracer => tracersHead  ! point to the mass-like field, iceAreaCategory
       if (verboseGlobal .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'Cells with initial ice:'
          write(stderrUnit,*) 'iCell (local), iCell(global), ice area, xCell, yCell, zCell:'
          if (thisTracer % ndims == 2) then  ! the pointer to iceAreaCategory should always be 2D, but check just in case
             do iCell = 1, nCellsSolve
                if (thisTracer % array2D(iCatTest,iCell) > 0.0_RKIND) then
                   write(stderrUnit,*) iCell, indexToCellID(iCell), thisTracer % array2D(iCatTest,iCell), &
                        xCell(iCell), yCell(iCell), zCell(iCell)
                endif
             enddo
          endif  ! ndims
       endif   ! verbose and iceAreaCategory
    endif   ! first_call

    !-------------------------------------------------------------------
    ! Compute masks.
    ! The tracer mask identifies cells whose values are physically meaningful.
    ! For example, the ice thickness is meaningful if and only if the area is nonzero,
    !  and the ice enthalpy is meaningful if and only if the area and thickness are nonzero.
    ! Note: There are two calls to the mask routine:
    !   (1) The first mask is used in computing local maxes and mins for an optional monotonicity check.
    !       It counts any positive ice area or thickness, however small, as physically meaningful.
    !       Once the local maxes and mins are computed, this mask is no longer needed.
    !   (2) The second mask is used below in gradient calculations. It masks out area and thickness
    !       values that are very small but nonzero, to ensure that junk tracer values are
    !       not used in computations.
    !-------------------------------------------------------------------

    if (verboseRun) write(stderrUnit,*) 'Make masks'

    call MPAS_pool_get_config(block % configs, "config_monotonicity_check", configMonotonicityCheck)
    if (configMonotonicityCheck) then

       if (verboseRun) write(stderrUnit,*) 'Compute local mins and maxes'

       call mpas_timer_start("incr remap masks")
       call make_masks(&
            nCells,      &
            maskCell,    &
            tracersHead, &
            threshold_in = 0.0_RKIND)
       call mpas_timer_stop("incr remap masks")

       call mpas_timer_start("incr remap monotonicity check")
       call tracer_local_min_max(&
            tracersHead,   &
            nCells,        &
            nCellsSolve,   &
            nEdgesOnCell,  &
            cellsOnCell)
       call mpas_timer_stop("incr remap monotonicity check")

    endif

    call mpas_timer_start("incr remap masks")
    call make_masks(&
         nCells,      &
         maskCell,    &
         tracersHead, &
         threshold_in = eps11)
    call mpas_timer_stop("incr remap masks")

    !-------------------------------------------------------------------
    ! Construct each tracer field as a linear function of x and y in each cell
    !-------------------------------------------------------------------

    if (verboseRun) write(stderrUnit,*) 'Construct fields'

    call mpas_timer_start("incr remap construct fields")
    call construct_linear_tracer_fields(&
         nCells,                        &
         maxEdges,                      &
         nEdgesOnCell,                  &
         edgesOnCell,                   &
         cellsOnCell,                   &
         cellsOnEdge,                   &
         xVertexOnCell, yVertexOnCell,  &
         dcEdge,                        &
         on_a_sphere,                   &
         config_rotate_cartesian_grid,  &
         transGlobalToCell,             &
         coeffsReconstruct,             &
         maskCell,                      &
         geomAvgCell,                   &
         tracersHead,                   &
         indexToCellID,                 &
         indexToEdgeID,                 &
         block)
    call mpas_timer_stop("incr remap construct fields")

    !-------------------------------------------------------------------
    ! Find the locations of departure points for each edge.
    ! Both the velocity components and the departure points are given in local east/north coordinates at vertices.
    !-------------------------------------------------------------------

    if (verboseRun) then
       write(stderrUnit,*) 'Find departure points'
       write(stderrUnit,*) 'Max (uvel, vvel) =', maxval(uVelocity), maxval(vVelocity)
       write(stderrUnit,*) 'Max (uvel, vvel)*dt =', maxval(uVelocity)*dt, maxval(vVelocity)*dt
    endif

    call mpas_timer_start("incr remap departure points")
    call find_departure_points(&
         domain,                              &
         nVertices,                           &
         dt,                                  &
         uVelocity,    vVelocity,             &
         minLengthEdgesOnVertex,              &
         indexToVertexID,                     &
         departurePoint,                      &
         on_a_sphere)
    call mpas_timer_stop("incr remap departure points")

    !-------------------------------------------------------------------
    ! For each edge, divide the departure region into triangles, with one
    ! or two triangles for each cell contributing a flux across the edge.
    ! Note: This is done for nEdges, which (unlike nEdgesSolve) includes
    !       all edges of locally owned cells.
    !-------------------------------------------------------------------

    if (verboseRun) write(stderrUnit,*) 'Locate departure triangles'

    call mpas_timer_start("incr remap departure triangles")
    call find_departure_triangles(&
         nEdges,                                    &
         nCells,                                    &
         maxEdges,                                  &
         nTriPerEdgeRemap,                          &
         nEdgesOnCell,                              &
         xVertexOnEdge,    yVertexOnEdge,           &
         xVertexOnCell,    yVertexOnCell,           &
         verticesOnCell,                            &
         verticesOnEdge,                            &
         edgesOnCell,                               &
         cellsOnEdge,                               &
         remapEdge,                                 &
         cellsOnEdgeRemap,                          &
         edgesOnEdgeRemap,                          &
         vertexDegree,                              &
         departurePoint,                            &
         maskEdge,                                  &
         xTriangle,        yTriangle,               &
         iCellTriangle,                             &
         triangleArea,                              &
         on_a_sphere,                               &
         transGlobalToCell,                         &
         indexToCellID,                             &
         indexToEdgeID,                             &
         indexToVertexID,                           &
         block)
    call mpas_timer_stop("incr remap departure triangles")

    if (verboseGeometry .and. etestOnProc .and. block % localBlockID == etestBlockID) then
       iEdge = etest
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'Triangles, iEdge =', indexToEdgeID(iEdge)
       do n = 1, nTriPerEdgeRemap
          iCell = iCellTriangle(n,iEdge)
          if (iCell >= 1 .and. iCell <= nCells) then
             write(stderrUnit,*) n, 'cell, area =', indexToCellID(iCell), triangleArea(n,iEdge)
             do m = 1, 3
                write(stderrUnit,*) m, xTriangle(m,n,iEdge), yTriangle(m,n,iEdge)
             enddo
          endif
       enddo
    endif

    !-------------------------------------------------------------------
    ! Compute quadrature points for each departure triangle.
    !-------------------------------------------------------------------

    call mpas_timer_start("incr remap quadrature points")
    call get_triangle_quadrature_points (&
         nEdges,             &
         nTriPerEdgeRemap,   &
         maskEdge,          &
         xTriangle, yTriangle)
    call mpas_timer_stop("incr remap quadrature points")

    if (verboseGeometry .and. etestOnProc .and. block % localBlockID == etestBlockID) then
       iEdge = etest
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'Quad points, iEdge =', indexToEdgeID(iEdge)
       do n = 1, nTriPerEdgeRemap
          if (maxval(xTriangle(:,n,iEdge)) /= 0.0_RKIND) then
             write(stderrUnit,*) 'Triangle', n
             do m = 1, nQuadPoints
                write(stderrUnit,*) m, xTriangle(m,n,iEdge), yTriangle(m,n,iEdge)
             enddo
          endif
       enddo
    endif

    !-------------------------------------------------------------------
    ! Integrate mass and tracer fluxes over each departure triangle.
    !-------------------------------------------------------------------

    call mpas_timer_start("incr remap integrate fluxes")
    call integrate_fluxes_over_triangles(&
         tracersHead,          &
         nCells,               &
         nEdges,               &
         maskEdge,             &
         xTriangle, yTriangle, &
         triangleArea,         &
         cellsOnEdge,          &
         iCellTriangle,        &
         indexToCellID,        &
         indexToEdgeID,        &
         block)
    call mpas_timer_stop("incr remap integrate fluxes")

    if (verboseFluxes .and. etestOnProc .and. block % localBlockID == etestBlockID) then
       iEdge = etest
       thisTracer => tracersHead   ! point to the mass-like field, iceAreaCategory

       do while(associated(thisTracer))
          if (trim(thisTracer % tracerName) == 'iceAreaCategory') then
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'Fluxes, iEdge, tracer =', indexToEdgeID(iEdge), trim(thisTracer % tracerName)
             if (thisTracer % ndims == 2) then
                nCategories = size(thisTracer % edgeFlux2D, 1)
                write(stderrUnit,*)'category, flux:'
                do iCat = 1, nCategories
                   write(stderrUnit,*) iCat, thisTracer % edgeFlux2D(iCat,iEdge)
                enddo
             elseif (thisTracer % ndims == 3) then
                nLayers = size(thisTracer % edgeFlux3D, 1)
                nCategories = size(thisTracer % edgeFlux3D, 2)
                write(stderrUnit,*)'layer, category, flux:'
                do iLayer = 1, nLayers
                   do iCat = 1, nCategories
                      write(stderrUnit,*) iLayer, iCat, thisTracer % edgeFlux3D(iLayer,iCat,iEdge)
                   enddo
                enddo
             endif
          endif

          thisTracer => thisTracer % next
       enddo     ! associated(thisTracer)
    endif     ! verboseFluxes & etestProcId & etestBlockID

    !-------------------------------------------------------------------
    ! Compute mass*tracer products
    ! Note: This subroutine must be called before calling sum_tracers or update_mass_and_tracers.
    !       (So it is needed even if conservationCheck = .false.)
    !-------------------------------------------------------------------

    call mpas_timer_start("incr remap mass-tracer products")
    call compute_mass_tracer_products(&
         nCells,           &
         tracersHead)
    call mpas_timer_stop("incr remap mass-tracer products")

    !-------------------------------------------------------------------
    ! Optionally, compare these mass*tracer products to the values computed
    !  by analytically integrating the reconstructed tracer fields over each cell
    !  (using the center/xGrad/yGrad values found in construct_linear_tracer_fields).
    !-------------------------------------------------------------------

    call MPAS_pool_get_config(block % configs, "config_recover_tracer_means_check", configRecoverTracerMeansCheck)
    if (configRecoverTracerMeansCheck) then

       if (verboseRun) then
          write(stderrUnit,*) 'Recover tracer means'
          write(stderrUnit,*) 'WARNING: This subroutine will considerably increase the computation time for IR'
       endif

       call mpas_timer_start("incr remap recover tracer means")
       call recover_tracer_means(&
            nCells,                        &
            maxEdges,                      &
            nEdgesOnCell,                  &
            xVertexOnCell, yVertexOnCell,  &
            edgesOnCell,                   &
            dcEdge,                        &
            dvEdge, &
            tracersHead)
       call mpas_timer_stop("incr remap recover tracer means")

    endif

    !-------------------------------------------------------------------
    ! Compute initial sums of mass*tracer over the locally owned cells
    !-------------------------------------------------------------------

    call MPAS_pool_get_config(block % configs, "config_conservation_check", configConservationCheck)
    if (configConservationCheck) then

       ! compute initial sums of mass*tracer
       call mpas_timer_start("incr remap conservation check")
       call sum_tracers(&
            tracersHead, &
            nCellsSolve, &
            areaCell,    &
            init=.true.)
       call mpas_timer_stop("incr remap conservation check")

    endif

    !-------------------------------------------------------------------
    ! Compute the new values of mass and tracers in each locally owned cell
    !-------------------------------------------------------------------

    call mpas_timer_start("incr remap update tracers")
    call update_mass_and_tracers(&
         nCellsSolve,      &
         nEdgesOnCell,     &
         edgesOnCell,      &
         cellsOnEdge,      &
         areaCell,         &
         tracersHead,      &
         indexToCellID,    &
         indexToEdgeID,    &
         block)
    call mpas_timer_stop("incr remap update tracers")

    !-------------------------------------------------------------------
    ! Compute final sums of mass*tracer over the locally owned cells
    !-------------------------------------------------------------------

    call MPAS_pool_get_config(block % configs, "config_conservation_check", configConservationCheck)
    if (configConservationCheck) then

       ! compute new mass*tracer products
       call mpas_timer_start("incr remap conservation check")

       call compute_mass_tracer_products(&
            nCells,        &
            tracersHead)

       ! compute final sums of mass*tracer
       call sum_tracers(&
            tracersHead, &
            nCellsSolve, &
            areaCell,    &
            init=.false.)

       call mpas_timer_stop("incr remap conservation check")

    endif

    ! Optionally, zero out the ice area and tracers in cells with very small ice areas
    ! TODO - Turn off this call if not needed (for example, if very small areas are zapped
    !        in the column-physics calculation).

    if (zapSmallMass) then

       if (verboseRun) write(stderrUnit,*) 'Zap small mass'

       call zap_small_mass(&
            tracersHead,   &
            nCellsSolve,   &
            maskCategoryCell)

    endif

    ! clean up
    deallocate(geomAvgCell)

    ! More optional diagnostics
    if (verboseRun .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then

       count = 0
       thisTracer => tracersHead

       do while(count < 3)  ! output for first 3 tracers
          count = count + 1

          iCell = ctest
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'New values, iCell, tracer =', indexToCellID(iCell), trim(thisTracer % tracerName)
          if (thisTracer % ndims == 2) then
             write(stderrUnit,*)'category, value:'
             nCategories = size(thisTracer % array2D, 1)
             do iCat = 1, nCategories
                write(stderrUnit,*) iCat, thisTracer % array2D(iCat,iCell)
             enddo
          elseif (thisTracer % ndims == 3) then
             write(stderrUnit,*)'layer, category, value:'
             nLayers = size(thisTracer % array3D, 1)
             nCategories = size(thisTracer % array3D, 2)
             do iCat = 1, nCategories
                do iLayer = 1, nLayers
                   write(stderrUnit,*) iLayer, iCat, thisTracer % array3D(iLayer,iCat,iCell)
                enddo
             enddo
          endif  ! ndims

          thisTracer => thisTracer % next
       enddo   ! count < 3

    endif   ! verboseRun & ctestOnProc & ctestBlockID

    ! Optional diagnostics: Write a list of cells with ice

    thisTracer => tracersHead  ! point to the mass-like field, iceAreaCategory
    if (verboseGlobal .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'Cells with ice:'
       write(stderrUnit,*) 'iCell (local), iCell(global), ice area, xCell, yCell, zCell:'
       if (thisTracer % ndims == 2) then  ! the pointer to iceAreaCategory should always be 2D, but check just in case
          do iCell = 1, nCellsSolve
             if (thisTracer % array2D(iCatTest,iCell) > 0.0_RKIND) then
                write(stderrUnit,*) iCell, indexToCellID(iCell), thisTracer % array2D(iCatTest,iCell), &
                     xCell(iCell), yCell(iCell), zCell(iCell)
             endif
          enddo
       endif  ! ndims
    endif   ! verbose and iceAreaCategory

  end subroutine incremental_remap_block

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine make_masks
!
!> \brief  integer masks for incremental remapping tracers
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine sets integer masks to indicate where tracer values
!>  are physically meaningful.
!
!-----------------------------------------------------------------------

  subroutine make_masks(&
       nCells,          &
       maskCell,        &
       tracersHead,     &
       threshold_in)

    integer, intent(in) ::  &
         nCells       ! number of cells on block

    integer, intent(out), dimension(:) :: &
         maskCell     ! = 1 for cells with initial nonzero mass, else = 0

    type(tracer_type), pointer :: &
         tracersHead  !< Input/output: pointer to first element of linked list of tracers
                      ! The 'arrayMask' variable of the tracer derived type is set here

    real(kind=RKIND), intent(in), optional :: &
         threshold_in ! threshold for deciding whether a small value is physically meaningful

    type(tracer_type), pointer :: &
         thisTracer   ! pointer that loops through linked list of tracers

    type(tracer_type), pointer :: &
         parentTracer ! pointer to parent of this tracer

    real(kind=RKIND) :: massSumCell

    real(kind=RKIND) :: &
         threshold    ! threshold for deciding whether a small value is physically meaningful

    integer :: nCategories, nLayers
    integer :: iCat, iCell, iLayer

    if (present(threshold_in)) then
       threshold = threshold_in
    else
       threshold = 0.0_RKIND
    endif

    ! Set tracer masks.
    ! The mask is set to 1 where tracer values are physically meaningful, else the mask = 0.
    ! For the mass-like field (fractional ice concentration for CICE), the value is deemed physically meaningful everywhere.
    ! For tracers with parents, the value is deemed physically meaningful wherever the parent tracer is > 0
    !  or some small threshold value.
    ! The reason for the optional threshold is as follows:
    ! - For some applications, such as monotonicity checks, we want the mask to be inclusive, and (for instance)
    !   to count very small ice area or ice/snow thickness as real. So the best threshold is 0.
    ! - For other applications, such as gradient computations, we want the mask to be exclusive, and
    !   to definitely ignore tracer values if the area or thickness is within roundoff level of zero
    !   (in which case the tracer values could be junk). So a good threshold is a small finite number like 1.e-11.
    !
    ! Note: This is different from the convention in standard CICE, where tracer values are deemed meaningful
    !       wherever the *parent* tracer mask = 1.
    !
    ! Note: An earlier version of this subroutine used 'where' statements instead of do loops.
    !       However, this led to errors when the pondArea mask was set to 1 for iCell = nCells+1 (i.e., non-existent cells),
    !        as a result of its parent, levelIceArea, being initialized to 1.0 everywhere (including non-existent cells).
    !       The masks are now initialized to 0 everywhere, and then set to 1 only in the range 1:nCells.

    thisTracer => tracersHead

    do while (associated(thisTracer))

       !  initialize
       if (associated(thisTracer % array2D)) then
          thisTracer % arrayMask2D(:,:) = 0
          nCategories = size(thisTracer % array2D, 1)
       elseif (associated(thisTracer % array3D)) then
          thisTracer % arrayMask3D(:,:,:) = 0
          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)
       endif

       if (thisTracer % nParents == 0) then  ! mass-like field

          if (associated(thisTracer % array2D)) then

             do iCell = 1, nCells

                ! set mask = 1 to indicate the value is meaningful everywhere
                thisTracer % arrayMask2D(:,iCell) = 1

                ! Set maskCell = 1 if ice exists in any category or layer; else set maskCell = 0
                massSumCell = sum(thisTracer%array2D(:,iCell))
                if (massSumCell > 0.0_RKIND) then
                   maskCell(iCell) = 1
                else
                   maskCell(iCell) = 0
                endif

             enddo

          elseif (associated(thisTracer % array3D)) then

             do iCell = 1, nCells

                ! set mask = 1 to indicate the value is meaningful everywhere
                thisTracer % arrayMask3D(:,:,iCell) = 1

                ! Set maskCell = 1 if ice exists in any category or layer; else set maskCell = 0
                massSumCell = sum(thisTracer%array3D(:,:,iCell))
                if (massSumCell > 0.0_RKIND) then
                   maskCell(iCell) = 1
                else
                   maskCell(iCell) = 0
                endif

             enddo

          endif

       else   ! nParents = 1, 2 or 3; value is meaningful where parent value > 0
              !TODO - Confirm that there are no physically meaningful parents with nonzero negative values
              !TODO - Change eps11 threshold to 0.0?

          parentTracer => thisTracer % parent

          if (associated(thisTracer % array2D)) then  ! parent must also be 2D

             do iCell = 1, nCells
                do iCat = 1, nCategories
                   if (parentTracer % array2D(iCat,iCell) > threshold) then
                      thisTracer % arrayMask2D(iCat,iCell) = 1
                   endif
                enddo
             enddo

          elseif (associated(thisTracer % array3D)) then  ! parent can be 2D or 3D

             if (parentTracer % ndims == 2) then

                do iCell = 1, nCells
                   do iCat = 1, nCategories
                      if (parentTracer % array2D(iCat,iCell) > threshold) then
                         thisTracer % arrayMask3D(:,iCat,iCell) = 1
                      endif
                   enddo
                enddo

             else   ! parent tracer is 3D

                do iCell = 1, nCells
                   do iCat = 1, nCategories
                      do iLayer = 1, nLayers
                         if (parentTracer % array3D(iLayer,iCat,iCell) > threshold) then
                            thisTracer % arrayMask3D(iLayer,iCat,iCell) = 1
                         endif
                      enddo
                   enddo
                enddo

             endif  ! parentTracer % ndims

          endif     ! associated(thisTracer % array2D)

       endif        ! thisTracer % nParents

       thisTracer => thisTracer % next

    enddo   ! associated(thisTracer)

  end subroutine make_masks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine construct_linear_tracer_fields
!
!> \brief  construct each tracer field as a linear function of x and y
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine constructs tracer fields (center value, plus x and y gradient
!>  components), such that the area-integrated value recovers the grid cell mean.
!
!-----------------------------------------------------------------------

  subroutine construct_linear_tracer_fields(&
       nCells,                       &
       maxEdges,                     &
       nEdgesOnCell,                 &
       edgesOnCell,                  &
       cellsOnCell,                  &
       cellsOnEdge,                  &
       xVertexOnCell, yVertexOnCell, &
       dcEdge,                       &
       on_a_sphere,                  &
       config_rotate_cartesian_grid, &
       transGlobalToCell,            &
       coeffsReconstruct,            &
       maskCell,                     &
       geomAvgCell,                  &
       tracersHead,                  &
       indexToCellID,                &
       indexToEdgeID,                &
       block)

    integer, intent(in) :: &
         nCells,         & !< Input: number of cells
                           !   Note: Need to construct fields in a layer of halo cells outside nCellsSolve
         maxEdges          !< Input: max number of edges per cell

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell      !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         edgesOnCell,    & !< Input: index for each edge of a given cell
         cellsOnCell,    & !< Input: index for each cell neighbor of a given cell
         cellsOnEdge       !< Input: index for each cell neighbor of a given edge

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnCell,  &  !< Input: x (east) coordinate of vertex relative to cell center in local tangent plane
         yVertexOnCell      !< Input: y (north) coordinate of vertex relative to cell center in local tangent plane

    real(kind=RKIND), dimension(:), intent(in) :: &
         dcEdge             !< Input: distance between the 2 cell centers on each side of an edge

    logical, intent(in) ::   &
         on_a_sphere,               &  !< Input: T if flow is on a sphere, F if on a plane
         config_rotate_cartesian_grid  !< Input: if true, then the North and South Poles are rotated to the equator

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         transGlobalToCell  !< Input: 3x3 matrix for transforming vectors from global to cell-based coordinates

    real(kind=RKIND), dimension(:,:,:), intent(in) ::  &
         coeffsReconstruct  !< Input: coefficients for reconstructing the gradient at a cell center,
                            !         given normal components on edges

    integer, dimension(:), intent(in) :: &
         maskCell           !< Input: = 1 for cells with ice, else = 0

    type(geometric_avg_cell_type), pointer :: &
         geomAvgCell        !< Input: derived type holding geometric averages

    type(tracer_type), pointer :: &
         tracersHead        !< Input/output: pointer to first element of linked list of tracers
                            ! unlimited gradients on input; limited gradients on output

    integer, dimension(:), intent(in) ::  &
         indexToCellID,   & !< Input: global index for each cell on a block (diagnostic only)
         indexToEdgeID      !< Input: global index for each edge on a block (diagnostic only)

    type(block_type), intent(in) :: &
         block             !< Input: local block (diagnostic only)

    ! local variables

    type(tracer_type), pointer :: &
         thisTracer,        &   ! pointer that loops through linked list of tracers
         parentTracer,      &   ! parent of thisTracer
         grandparentTracer, &   ! parent of parentTracer
         dummyTracer            ! dummy tracer with barycenter at cell centroid

    integer :: iCell, iCat, iLayer

    real(kind=RKIND) :: &
         c0, cx, cy, cxx, cxy, cyy,  &   ! coefficients in equations for barycentric coordinates
         cxxx, cxxy, cxyy, cyyy

    real(kind=RKIND) :: reciprocal

    integer :: nCategories   ! number of categories
    integer :: nLayers       ! number of layers

    ! loop through linked list of tracers

    thisTracer => tracersHead
    do while (associated(thisTracer))

       if (verboseConstruct) then
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'tracer: ', trim(thisTracer % tracerName)
       endif

       ! initialize quantities needed for tracer reconstruction
       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)

          thisTracer % center2D(:,:) = 0.0_RKIND
          thisTracer % xGrad2D (:,:) = 0.0_RKIND
          thisTracer % yGrad2D (:,:) = 0.0_RKIND
          thisTracer % xBarycenter2D(:,:) = 0.0_RKIND
          thisTracer % yBarycenter2D(:,:) = 0.0_RKIND

       elseif (thisTracer % ndims == 3) then

          nLayers     = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)

          thisTracer % center3D(:,:,:) = 0.0_RKIND
          thisTracer % xGrad3D (:,:,:) = 0.0_RKIND
          thisTracer % yGrad3D (:,:,:) = 0.0_RKIND
          thisTracer % xBarycenter3D(:,:,:) = 0.0_RKIND
          thisTracer % yBarycenter3D(:,:,:) = 0.0_RKIND

       endif   ! ndims

       ! set parent tracer
       if (thisTracer % nParents == 0) then   ! mass-like field (fractional ice concentration for CICE)

          ! This tracer has no parent, but create a dummy parent tracer with (xBarycenter, yBarycenter) at cell centroid.
          ! This ensures correct values for xBarycenter and yBarycenter below.
          allocate(dummyTracer)

          if (thisTracer % ndims == 2) then

             allocate(dummyTracer % xBarycenter2D(nCategories, nCells))
             allocate(dummyTracer % yBarycenter2D(nCategories, nCells))
             do iCell = 1, nCells
                dummyTracer % xBarycenter2D(:,iCell) = geomAvgCell % x(iCell)
                dummyTracer % yBarycenter2D(:,iCell) = geomAvgCell % y(iCell)
             enddo
             dummyTracer % ndims = 2

          elseif (thisTracer % ndims == 3) then

             allocate(dummyTracer % xBarycenter3D(nLayers, nCategories, nCells))
             allocate(dummyTracer % yBarycenter3D(nLayers, nCategories, nCells))
             do iCell = 1, nCells
                dummyTracer % xBarycenter3D(:,:,iCell) = geomAvgCell % x(iCell)
                dummyTracer % yBarycenter3D(:,:,iCell) = geomAvgCell % y(iCell)
             enddo
             dummyTracer % ndims = 3

          endif

          parentTracer => dummyTracer

       elseif (thisTracer % nParents == 1) then

          parentTracer => thisTracer % parent

       else   ! nParents = 2 or 3

          parentTracer => thisTracer % parent
          grandparentTracer => parentTracer % parent

       endif  ! nParents

       ! compute quantities needed for linear reconstruction: cell center value plus x and y gradient components

       if (thisTracer % ndims == 2) then

          ! compute unlimited gradient
          call compute_gradient_2d(&
               nCells,                       &
               maxEdges,                     &
               nEdgesOnCell,                 &
               edgesOnCell,                  &
               cellsOnCell,                  &
               cellsOnEdge,                  &
               dcEdge,                       &
               coeffsReconstruct,            &
               on_a_sphere,                  &
               config_rotate_cartesian_grid, &
               transGlobalToCell,            &
               maskCell,                     &
               thisTracer % array2D,         &
               thisTracer % arrayMask2D,     &
               thisTracer % xGrad2D,         &
               thisTracer % yGrad2D,         &
               indexToCellID,                &
               indexToEdgeID,                &
               block)

          ! limit the gradient
          call limit_tracer_gradient_2d(&
               nCells,                        &
               nEdgesOnCell,                  &
               cellsOnCell,                   &
               xVertexOnCell, yVertexOnCell,  &
               maskCell,                      &
               thisTracer % array2D,          &
               thisTracer % arrayMask2D,      &
               parentTracer % xBarycenter2D,  &
               parentTracer % yBarycenter2D,  &
               thisTracer % xGrad2D,          &
               thisTracer % yGrad2D)

          ! compute the value of the field at the geometric center of the cell

          do iCell = 1, nCells
             do iCat = 1, nCategories
                thisTracer % center2D(iCat,iCell) = thisTracer % array2D(iCat,iCell)  &
                                                  - thisTracer % xGrad2D(iCat,iCell) * parentTracer % xBarycenter2D(iCat,iCell)  &
                                                  - thisTracer % yGrad2D(iCat,iCell) * parentTracer % yBarycenter2D(iCat,iCell)
             enddo
          enddo

          if (verboseConstruct .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
             iCell = ctest
             write(stderrUnit,*) '  Limited center gradient:', thisTracer % xGrad2D(iCatTest,iCell), &
                                                               thisTracer % yGrad2D(iCatTest,iCell)
             write(stderrUnit,*) 'Center value:', thisTracer % center2D(iCatTest,iCell)
          endif

          ! if this tracer has children, then compute its barycenter, as required for reconstruction of the child tracer

          if (thisTracer % hasChild) then

             if (thisTracer % nParents == 0) then   ! this tracer is the mass-like field (fractional ice concentration for CICE)

                do iCell = 1, nCells
                   if (maskCell(iCell) == 1) then
                      do iCat = 1, nCategories

                         call compute_barycenter_coordinates(&
                              geomAvgCell,                             &
                              iCell,                                   &
                              thisTracer % xBarycenter2D(iCat,iCell),  &
                              thisTracer % yBarycenter2D(iCat,iCell),  &
                              thisTracer % array2D      (iCat,iCell),  &
                              thisTracer % center2D     (iCat,iCell),  &
                              thisTracer % xGrad2D      (iCat,iCell),  &
                              thisTracer % yGrad2D      (iCat,iCell))

                      enddo  ! iCat
                   endif    ! maskCell
                enddo     ! iCell

             elseif (thisTracer % nParents == 1) then   ! the mass-like field is the parent

                do iCell = 1, nCells
                   if (maskCell(iCell) == 1) then
                      do iCat = 1, nCategories

                         call compute_barycenter_coordinates(&
                              geomAvgCell,                             &
                              iCell,                                   &
                              thisTracer % xBarycenter2D(iCat,iCell),  &
                              thisTracer % yBarycenter2D(iCat,iCell),  &
                              parentTracer % array2D    (iCat,iCell),  &
                              parentTracer % center2D   (iCat,iCell),  &
                              parentTracer % xGrad2D    (iCat,iCell),  &
                              parentTracer % yGrad2D    (iCat,iCell),  &
                              thisTracer % array2D      (iCat,iCell),  &
                              thisTracer % center2D     (iCat,iCell),  &
                              thisTracer % xGrad2D      (iCat,iCell),  &
                              thisTracer % yGrad2D      (iCat,iCell))

                      enddo  ! iCat
                   endif    ! maskCell
                enddo     ! iCell

             elseif (thisTracer % nParents == 2) then   ! the mass-like field is the grandparent

                do iCell = 1, nCells
                   if (maskCell(iCell) == 1) then
                      do iCat = 1, nCategories

                         call compute_barycenter_coordinates(&
                              geomAvgCell,                               &
                              iCell,                                     &
                              thisTracer % xBarycenter2D  (iCat,iCell),  &
                              thisTracer % yBarycenter2D  (iCat,iCell),  &
                              grandparentTracer % array2D (iCat,iCell),  &
                              grandparentTracer % center2D(iCat,iCell),  &
                              grandparentTracer % xGrad2D (iCat,iCell),  &
                              grandparentTracer % yGrad2D (iCat,iCell),  &
                              parentTracer % array2D      (iCat,iCell),  &
                              parentTracer % center2D     (iCat,iCell),  &
                              parentTracer % xGrad2D      (iCat,iCell),  &
                              parentTracer % yGrad2D      (iCat,iCell),  &
                              thisTracer % array2D        (iCat,iCell),  &
                              thisTracer % center2D       (iCat,iCell),  &
                              thisTracer % xGrad2D        (iCat,iCell),  &
                              thisTracer % yGrad2D        (iCat,iCell))

                      enddo  ! iCat
                   endif    ! maskCell
                enddo     ! iCell

             elseif (thisTracer % nParents == 3) then

                write(stderrUnit,*) 'IR, construct_linear_tracer_fields: Tracers with 3 parents should not have children'
                call mpas_dmpar_global_abort(&
                     'MPAS-seaice: IR, construct_linear_tracer_fields: Tracers with 3 parents should not have children')

             endif   ! nParents

          endif      ! hasChild

       elseif (thisTracer % ndims == 3) then

          ! compute unlimited gradient
          call compute_gradient_3d(&
               nCells,                       &
               maxEdges,                     &
               nEdgesOnCell,                 &
               edgesOnCell,                  &
               cellsOnCell,                  &
               cellsOnEdge,                  &
               dcEdge,                       &
               coeffsReconstruct,            &
               on_a_sphere,                  &
               config_rotate_cartesian_grid, &
               transGlobalToCell,            &
               maskCell,                     &
               thisTracer % array3D,         &
               thisTracer % arrayMask3D,     &
               thisTracer % xGrad3D,         &
               thisTracer % yGrad3D,         &
               indexToCellID,                &
               indexToEdgeID,                &
               block)

          ! Note: The following argument lists are different depending on whether the parent tracer is 2D or 3D

          if (parentTracer % ndims == 2) then

             ! limit the gradient

             call limit_tracer_gradient_3d(&
                  nCells,                        &
                  nEdgesOnCell,                  &
                  cellsOnCell,                   &
                  xVertexOnCell, yVertexOnCell,  &
                  maskCell,                      &
                  thisTracer % array3D,          &
                  thisTracer % arrayMask3D,      &
                  thisTracer % xGrad3D,          &
                  thisTracer % yGrad3D,          &
                  xBarycenter2D = parentTracer % xBarycenter2D,  &
                  yBarycenter2D = parentTracer % yBarycenter2D)

             ! compute the value of the field at the geometric center of the cell

             do iCell = 1, nCells
                do iCat = 1, nCategories
                   do iLayer = 1, nLayers
                      thisTracer % center3D(iLayer,iCat,iCell) = &
                           thisTracer % array3D(iLayer,iCat,iCell)  &
                         - thisTracer % xGrad3D(iLayer,iCat,iCell) * parentTracer % xBarycenter2D(iCat,iCell)  &
                         - thisTracer % yGrad3D(iLayer,iCat,iCell) * parentTracer % yBarycenter2D(iCat,iCell)
                   enddo
                enddo
             enddo

          elseif (parentTracer % ndims == 3) then

             ! limit the gradient

             call limit_tracer_gradient_3d(&
                  nCells,                        &
                  nEdgesOnCell,                  &
                  cellsOnCell,                   &
                  xVertexOnCell, yVertexOnCell,  &
                  maskCell,                      &
                  thisTracer % array3D,          &
                  thisTracer % arrayMask3D,      &
                  thisTracer % xGrad3D,          &
                  thisTracer % yGrad3D,          &
                  xBarycenter3D = parentTracer % xBarycenter3D,  &
                  yBarycenter3D = parentTracer % yBarycenter3D)

             ! compute the value of the field at the geometric center of the cell

             do iCell = 1, nCells
                do iCat = 1, nCategories
                   do iLayer = 1, nLayers
                      thisTracer % center3D(iLayer,iCat,iCell) = &
                           thisTracer % array3D(iLayer,iCat,iCell)  &
                         - thisTracer % xGrad3D(iLayer,iCat,iCell) * parentTracer % xBarycenter3D(iLayer,iCat,iCell)  &
                         - thisTracer % yGrad3D(iLayer,iCat,iCell) * parentTracer % yBarycenter3D(iLayer,iCat,iCell)
                   enddo
                enddo
             enddo

          endif  ! parentTracer % ndims

          if (verboseConstruct .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
             iCell = ctest
             write(stderrUnit,*) '  Limited center gradient:', thisTracer % xGrad3D(iLayerTest,iCatTest,iCell), &
                                                               thisTracer % yGrad3D(iLayerTest,iCatTest,iCell)
             write(stderrUnit,*) 'Center value:', thisTracer % center3D(iLayer,iCat,iCell)
          endif

          ! if this tracer has children, then compute its barycenter, as required for reconstruction of the child tracer

          if (thisTracer % hasChild) then

             if (thisTracer % nParents == 0) then   ! this tracer is the mass-like field (fractional ice concentration for CICE)

                do iCell = 1, nCells
                   if (maskCell(iCell) == 1) then
                      do iCat = 1, nCategories
                         do iLayer = 1, nLayers

                            call compute_barycenter_coordinates(&
                                 geomAvgCell,                                    &
                                 iCell,                                          &
                                 thisTracer % xBarycenter3D(iLayer,iCat,iCell),  &
                                 thisTracer % yBarycenter3D(iLayer,iCat,iCell),  &
                                 thisTracer % array3D      (iLayer,iCat,iCell),  &
                                 thisTracer % center3D     (iLayer,iCat,iCell),  &
                                 thisTracer % xGrad3D      (iLayer,iCat,iCell),  &
                                 thisTracer % yGrad3D      (iLayer,iCat,iCell))

                         enddo  ! iLayer
                      enddo     ! iCat
                   endif       ! maskCell
                enddo        ! iCell

             elseif (thisTracer % nParents == 1) then   ! the parent is the mass-like field

                if (parentTracer % ndims == 2) then

                   do iCell = 1, nCells
                      if (maskCell(iCell) == 1) then
                         do iCat = 1, nCategories
                            do iLayer = 1, nLayers

                               call compute_barycenter_coordinates(&
                                    geomAvgCell,                                    &
                                    iCell,                                          &
                                    thisTracer % xBarycenter3D(iLayer,iCat,iCell),  &
                                    thisTracer % yBarycenter3D(iLayer,iCat,iCell),  &
                                    parentTracer % array2D           (iCat,iCell),  &
                                    parentTracer % center2D          (iCat,iCell),  &
                                    parentTracer % xGrad2D           (iCat,iCell),  &
                                    parentTracer % yGrad2D           (iCat,iCell),  &
                                    thisTracer % array3D      (iLayer,iCat,iCell),  &
                                    thisTracer % center3D     (iLayer,iCat,iCell),  &
                                    thisTracer % xGrad3D      (iLayer,iCat,iCell),  &
                                    thisTracer % yGrad3D      (iLayer,iCat,iCell))

                            enddo  ! iLayer
                         enddo     ! iCat
                      endif       ! maskCell
                   enddo        ! iCell

                elseif (parentTracer % ndims == 3) then

                   do iCell = 1, nCells
                      if (maskCell(iCell) == 1) then
                         do iCat = 1, nCategories
                            do iLayer = 1, nLayers

                               call compute_barycenter_coordinates(&
                                    geomAvgCell,                                    &
                                    iCell,                                          &
                                    thisTracer % xBarycenter3D(iLayer,iCat,iCell),  &
                                    thisTracer % yBarycenter3D(iLayer,iCat,iCell),  &
                                    parentTracer % array3D    (iLayer,iCat,iCell),  &
                                    parentTracer % center3D   (iLayer,iCat,iCell),  &
                                    parentTracer % xGrad3D    (iLayer,iCat,iCell),  &
                                    parentTracer % yGrad3D    (iLayer,iCat,iCell),  &
                                    thisTracer % array3D      (iLayer,iCat,iCell),  &
                                    thisTracer % center3D     (iLayer,iCat,iCell),  &
                                    thisTracer % xGrad3D      (iLayer,iCat,iCell),  &
                                    thisTracer % yGrad3D      (iLayer,iCat,iCell))

                            enddo  ! iLayer
                         enddo     ! iCat
                      endif       ! maskCell
                   enddo        ! iCell

                endif   ! parentTracer % ndims

             elseif (thisTracer % nParents == 2) then   ! the grandparent is the mass-like field

                if (parentTracer % ndims == 2) then   ! grandparent also must have ndims = 2

                   do iCell = 1, nCells
                      if (maskCell(iCell) == 1) then
                         do iCat = 1, nCategories
                            do iLayer = 1, nLayers

                               call compute_barycenter_coordinates(&
                                    geomAvgCell,                               &
                                    iCell,                                     &
                                    thisTracer % xBarycenter3D  (iLayer,iCat,iCell),  &
                                    thisTracer % yBarycenter3D  (iLayer,iCat,iCell),  &
                                    grandparentTracer % array2D        (iCat,iCell),  &
                                    grandparentTracer % center2D       (iCat,iCell),  &
                                    grandparentTracer % xGrad2D        (iCat,iCell),  &
                                    grandparentTracer % yGrad2D        (iCat,iCell),  &
                                    parentTracer % array2D             (iCat,iCell),  &
                                    parentTracer % center2D            (iCat,iCell),  &
                                    parentTracer % xGrad2D             (iCat,iCell),  &
                                    parentTracer % yGrad2D             (iCat,iCell),  &
                                    thisTracer % array3D        (iLayer,iCat,iCell),  &
                                    thisTracer % center3D       (iLayer,iCat,iCell),  &
                                    thisTracer % xGrad3D        (iLayer,iCat,iCell),  &
                                    thisTracer % yGrad3D        (iLayer,iCat,iCell))

                            enddo  ! iLayer
                         enddo     ! iCat
                      endif       ! maskCell
                   enddo        ! iCell

                elseif (parentTracer % ndims == 3) then

                   if (grandparentTracer % ndims == 2) then

                      do iCell = 1, nCells
                         if (maskCell(iCell) == 1) then
                            do iCat = 1, nCategories
                               do iLayer = 1, nLayers

                                  call compute_barycenter_coordinates(&
                                       geomAvgCell,                               &
                                       iCell,                                     &
                                       thisTracer % xBarycenter3D  (iLayer,iCat,iCell),  &
                                       thisTracer % yBarycenter3D  (iLayer,iCat,iCell),  &
                                       grandparentTracer % array2D        (iCat,iCell),  &
                                       grandparentTracer % center2D       (iCat,iCell),  &
                                       grandparentTracer % xGrad2D        (iCat,iCell),  &
                                       grandparentTracer % yGrad2D        (iCat,iCell),  &
                                       parentTracer % array3D      (iLayer,iCat,iCell),  &
                                       parentTracer % center3D     (iLayer,iCat,iCell),  &
                                       parentTracer % xGrad3D      (iLayer,iCat,iCell),  &
                                       parentTracer % yGrad3D      (iLayer,iCat,iCell),  &
                                       thisTracer % array3D        (iLayer,iCat,iCell),  &
                                       thisTracer % center3D       (iLayer,iCat,iCell),  &
                                       thisTracer % xGrad3D        (iLayer,iCat,iCell),  &
                                       thisTracer % yGrad3D        (iLayer,iCat,iCell))

                               enddo  ! iLayer
                            enddo     ! iCat
                         endif       ! maskCell
                      enddo        ! iCell

                   elseif (grandparentTracer % ndims == 3) then

                      do iCell = 1, nCells
                         if (maskCell(iCell) == 1) then
                            do iCat = 1, nCategories
                               do iLayer = 1, nLayers

                                  call compute_barycenter_coordinates(&
                                       geomAvgCell,                                      &
                                       iCell,                                            &
                                       thisTracer % xBarycenter3D  (iLayer,iCat,iCell),  &
                                       thisTracer % yBarycenter3D  (iLayer,iCat,iCell),  &
                                       grandparentTracer % array3D (iLayer,iCat,iCell),  &
                                       grandparentTracer % center3D(iLayer,iCat,iCell),  &
                                       grandparentTracer % xGrad3D (iLayer,iCat,iCell),  &
                                       grandparentTracer % yGrad3D (iLayer,iCat,iCell),  &
                                       parentTracer % array3D      (iLayer,iCat,iCell),  &
                                       parentTracer % center3D     (iLayer,iCat,iCell),  &
                                       parentTracer % xGrad3D      (iLayer,iCat,iCell),  &
                                       parentTracer % yGrad3D      (iLayer,iCat,iCell),  &
                                       thisTracer % array3D        (iLayer,iCat,iCell),  &
                                       thisTracer % center3D       (iLayer,iCat,iCell),  &
                                       thisTracer % xGrad3D        (iLayer,iCat,iCell),  &
                                       thisTracer % yGrad3D        (iLayer,iCat,iCell))

                               enddo  ! iLayer
                            enddo     ! iCat
                         endif       ! maskCell
                      enddo        ! iCell

                   endif   ! grandparentTracer % ndims

                elseif (thisTracer % nParents == 3) then

                   write(stderrUnit,*) 'ERROR: IR, construct_linear_tracer_fields: Tracers with 3 parents should not have children'
                   call mpas_dmpar_global_abort(&
                        'MPAS-seaice: IR, construct_linear_tracer_fields: Tracers with 3 parents should not have children')

                endif   ! parentTracer % ndims

             endif      ! thisTracer % nParents

          endif         ! thisTracer % hasChild

       endif            ! thisTracer % ndims

       ! clean up
       if (associated(dummyTracer)) then

          if (dummyTracer % ndims == 2) then
             if (associated(dummyTracer%xBarycenter2D)) deallocate(dummyTracer % xBarycenter2D)
             if (associated(dummyTracer%yBarycenter2D)) deallocate(dummyTracer % yBarycenter2D)
          elseif (dummyTracer % ndims == 3) then
             if (associated(dummyTracer%xBarycenter3D)) deallocate(dummyTracer % xBarycenter3D)
             if (associated(dummyTracer%yBarycenter3D)) deallocate(dummyTracer % yBarycenter3D)
          endif
          deallocate(dummyTracer)

       endif

       thisTracer => thisTracer % next

    enddo   ! while(associated)

  end subroutine construct_linear_tracer_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine compute_gradient_2d
!
!> \brief  compute the gradient of a 2D scalar field
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine computes the x and y coordinates of the gradient of a
!>  2D scalar field. The gradient is located at cell centers.
!
!-----------------------------------------------------------------------

  subroutine compute_gradient_2d(&
       nCells,                       &
       maxEdges,                     &
       nEdgesOnCell,                 &
       edgesOnCell,                  &
       cellsOnCell,                  &
       cellsOnEdge,                  &
       dcEdge,                       &
       coeffsReconstruct,            &
       on_a_sphere,                  &
       config_rotate_cartesian_grid, &
       transGlobalToCell,            &
       maskCell,                     &
       field,                        &
       mask,                         &
       xGrad,  yGrad,                &
       indexToCellID,                &
       indexToEdgeID,                &
       block)

    integer, intent(in) :: &
         nCells,           & !< Input: number of cells
         maxEdges            !< Input: max number of edges per cell

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell        !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         edgesOnCell,      & !< Input: index for each edge of a given cell
         cellsOnCell,      & !< Input: index for each cell neighbor of a given cell
         cellsOnEdge         !< Input: index for each cell neighbor of a given edge

    real(kind=RKIND), dimension(:), intent(in) :: &
         dcEdge              !< Input: distance between the 2 cell centers on each side of an edge

    real(kind=RKIND), dimension(:,:,:), intent(in) ::  &
         coeffsReconstruct   !< Input: coefficients for reconstructing the gradient at a cell center,
                             !         given normal components on edges

    logical, intent(in) ::   &
         on_a_sphere,               &  !< Input: T if flow is on a sphere, F if on a plane
         config_rotate_cartesian_grid  !< Input: if true, then the North and South Poles are rotated to the equator

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         transGlobalToCell   !< Input: 3x3 matrix for transforming vectors from global to cell-based coordinates

    integer, dimension(:), intent(in) :: &
         maskCell            !< Input: = 1 for cells with ice, else = 0

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         field               !< Input: 2D scalar field defined at cell centers
                             !         1st index = nLayers; 2nd index = nCategories, 3rd index = nCells

    integer, dimension(:,:), intent(in) :: &
         mask                !< Input: integer mask for parent tracer;
                             !         = 1 where field values for this tracer are physically meaningful, else = 0

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xGrad, yGrad        !< Output: x and y components of the gradient

    integer, dimension(:), intent(in) ::  &
         indexToCellID,    & !< Input: global index for each cell on a block (diagnostic only)
         indexToEdgeID       !< Input: global index for each edge on a block (diagnostic only)

    type(block_type), intent(in) :: &
         block               !< Input: local block (diagnostic only)

    ! local variables

    integer :: iCat, nCategories

    integer :: iCell, iEdge, iEdgeOnCell, iCellNeighbor

    real(kind=RKIND) :: &
         signGradient      ! = 1 or -1, depending on which direction is taken as positive at a given edge

    real(kind=RKIND), dimension(:,:), allocatable ::   &
         normalGrad,     & ! normal components of the gradient, defined on cell edges
         globalGrad,     & ! gradient at cell center, in global x/y/z coordinates
         globalGradRotate  ! rotated globalGrad vector

    real(kind=RKIND), dimension(:), allocatable ::   &
         zGrad             ! diagnostic only; should be much smaller than xGrad and yGrad

    ! find dimensions and allocate arrays
    nCategories = size(field,1)
    allocate(zGrad(nCategories))
    allocate(normalGrad(nCategories,maxEdges))
    allocate(globalGrad(nCategories,3))
    if (config_rotate_cartesian_grid) allocate(globalGradRotate(nCategories,3))

    ! initialize the gradient
    xGrad(:,:) = 0.0_RKIND
    yGrad(:,:) = 0.0_RKIND
    zGrad(:) = 0.0_RKIND  ! diagnostic only

    ! loop over cells
    do iCell = 1, nCells

       if (maskCell(iCell) == 1) then  ! ice is present in the cell

          ! initialize gradient components
          normalGrad(:,:) = 0.0_RKIND
          globalGrad(:,:) = 0.0_RKIND

          ! loop over edges of this cell
          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             iCellNeighbor = cellsOnCell(iEdgeOnCell, iCell)

             ! compute the normal component of the gradient on this edge
             if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then ! there is a cell neighbor on this edge

                iEdge = edgesOnCell(iEdgeOnCell,iCell)

                do iCat = 1, nCategories

                   ! field values in cell and its neighbor have physical meaning
                   if (mask(iCat,iCell) == 1 .and. mask(iCat,iCellNeighbor) == 1) then

                      if (iCell == cellsOnEdge(1,iEdge)) then
                         signGradient =  1.0_RKIND
                      else
                         signGradient = -1.0_RKIND
                      end if

                      normalGrad(iCat,iEdgeOnCell) = signGradient * (field(iCat,iCellNeighbor) - field(iCat,iCell)) / dcEdge(iEdge)

                   else   ! either or both field values do not have physical meaning; set gradient component = 0

                      normalGrad(iCat,iEdgeOnCell) = 0.0_RKIND

                   endif

                enddo   ! iCat

             else   ! there is no cell neighbor on this edge

                ! set gradient component = 0
                normalGrad(:,iEdgeOnCell) = 0.0_RKIND

             endif

             ! add the contribution of this normal component to the reconstructed
             ! gradient at the cell center (in global x/y/z coordinates)
             globalGrad(:,1) = globalGrad(:,1) + coeffsReconstruct(1,iEdgeOnCell,iCell) * normalGrad(:,iEdgeOnCell)
             globalGrad(:,2) = globalGrad(:,2) + coeffsReconstruct(2,iEdgeOnCell,iCell) * normalGrad(:,iEdgeOnCell)
             globalGrad(:,3) = globalGrad(:,3) + coeffsReconstruct(3,iEdgeOnCell,iCell) * normalGrad(:,iEdgeOnCell)

          enddo  ! iEdgeOnCell

          ! The global gradient vector lives in Cartesian x/y/z space.
          ! If on a rotated grid, then rotate the vector by 90 degrees, such that
          !  the true North Pole goes to (-1,0,0) and the true South Pole to (1,0,0).

          !TODO - Rotate the gradient vector if on a plane?

          if (config_rotate_cartesian_grid .and. on_a_sphere) then
             globalGradRotate(:,1) = -globalGrad(:,3)   ! xR = -z
             globalGradRotate(:,2) =  globalGrad(:,2)   ! yR = y
             globalGradRotate(:,3) =  globalGrad(:,1)   ! zR = x
             globalGrad(:,:) = globalGradRotate(:,:)
          endif

          ! transform from global x/y/z coordinates to local east/west coordinates

          if (on_a_sphere) then

             xGrad(:,iCell) = transGlobalToCell(1,1,iCell) * globalGrad(:,1)  &
                            + transGlobalToCell(1,2,iCell) * globalGrad(:,2)  &
                            + transGlobalToCell(1,3,iCell) * globalGrad(:,3)

             yGrad(:,iCell) = transGlobalToCell(2,1,iCell) * globalGrad(:,1)  &
                            + transGlobalToCell(2,2,iCell) * globalGrad(:,2)  &
                            + transGlobalToCell(2,3,iCell) * globalGrad(:,3)

             ! Note: The zGrad component is never used; it is simply computed as a diagnostic
             !TODO - Comment out zGrad computation?
             zGrad(:) = transGlobalToCell(3,1,iCell) * globalGrad(:,1)  &
                      + transGlobalToCell(3,2,iCell) * globalGrad(:,2)  &
                      + transGlobalToCell(3,3,iCell) * globalGrad(:,3)

          else  ! on a plane; do a simple copy

             xGrad(:,iCell) = globalGrad(:,1)
             yGrad(:,iCell) = globalGrad(:,2)

          endif

       endif   ! maskCell = 1

       if (verboseConstruct .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
          if (iCell == ctest) then
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'iCat, iCell =', iCatTest, indexToCellID(iCell)
             write(stderrUnit,*) 'field value =', field(iCatTest,iCell)
             write(stderrUnit,*) 'neighbor values and edge gradients:'
             do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                iCellNeighbor = cellsOnCell(iEdgeOnCell,iCell)
                iEdge = edgesOnCell(iEdgeOnCell,iCell)
                write(stderrUnit,*) iEdgeOnCell, indexToEdgeID(iEdge), &
                     field(iCatTest,iCellNeighbor), normalGrad(iCatTest,iEdgeOnCell)
             enddo
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'Unlimited center gradient:', xGrad(iCatTest,iCell), yGrad(iCatTest,iCell), zGrad(iCatTest)
          endif
       endif

    enddo  ! iCell

    ! cleanup
    deallocate(zGrad)
    deallocate(normalGrad)
    deallocate(globalGrad)
    if (config_rotate_cartesian_grid) deallocate(globalGradRotate)

  end subroutine compute_gradient_2d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine compute_gradient_3d
!
!> \brief  compute the gradient of a 3D scalar field
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine computes the x and y coordinates of the gradient of a
!>  3D scalar field. The gradient is located at cell centers.
!
!-----------------------------------------------------------------------

  subroutine compute_gradient_3d(&
       nCells,                       &
       maxEdges,                     &
       nEdgesOnCell,                 &
       edgesOnCell,                  &
       cellsOnCell,                  &
       cellsOnEdge,                  &
       dcEdge,                       &
       coeffsReconstruct,            &
       on_a_sphere,                  &
       config_rotate_cartesian_grid, &
       transGlobalToCell,            &
       maskCell,                     &
       field,                        &
       mask,                         &
       xGrad,  yGrad,                &
       indexToCellID,                &
       indexToEdgeID,                &
       block)

    integer, intent(in) :: &
         nCells,           & !< Input: number of cells
         maxEdges            !< Input: max number of edges per cell

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell        !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         edgesOnCell,      & !< Input: index for each edge of a given cell
         cellsOnCell,      & !< Input: index for each cell neighbor of a given cell
         cellsOnEdge         !< Input: index for each cell neighbor of a given edge

    real(kind=RKIND), dimension(:), intent(in) :: &
         dcEdge              !< Input: distance between the 2 cell centers on each side of an edge

    real(kind=RKIND), dimension(:,:,:), intent(in) ::  &
         coeffsReconstruct   !< Input: coefficients for reconstructing the gradient at a cell center,
                             ! given normal components on edges

    logical, intent(in) ::   &
         on_a_sphere,               &  !< Input: T if flow is on a sphere, F if on a plane
         config_rotate_cartesian_grid  !< Input: if true, then the North and South Poles are rotated to the equator

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         transGlobalToCell   !< Input: 3x3 matrix for transforming vectors from global to cell-based coordinates

    integer, dimension(:), intent(in) :: &
         maskCell            !< Input: = 1 for cells with ice, else = 0

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         field               !< Input: 3D scalar field defined at cell centers
                             !         1st index = nCategories, 2nd index = nCells

    integer, dimension(:,:,:), intent(in) :: &
         mask                !< Input: integer mask for parent tracer;
                             !         = 1 where field values for this tracer are physically meaningful, else = 0

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         xGrad, yGrad        !< Output: x and y components of the gradient

    integer, dimension(:), intent(in) ::  &
         indexToCellID,    & !< Input: global index for each cell on a block (diagnostic only)
         indexToEdgeID       !< Input: global index for each edge on a block (diagnostic only)

    type(block_type), intent(in) :: &
         block               !< Input: local block (diagnostic only)

    ! local variables

    integer :: iCat, iLayer, nCategories, nLayers

    integer :: iCell, iEdge, iEdgeOnCell, iCellNeighbor

    real(kind=RKIND), dimension(:,:,:), allocatable ::   &
         normalGrad,     & ! normal components of the gradient, defined on cell edges
         globalGrad,     & ! gradient at cell center, in global x/y/z coordinates
         globalGradRotate  ! rotated globalGrad vector

    real(kind=RKIND) :: &
         signGradient      ! = 1 or -1, depending on which direction is taken as positive at a given edge

    real(kind=RKIND), dimension(:,:), allocatable :: &
         zGrad             ! diagnostic only; should be much smaller than xGrad and yGrad

    ! find dimensions and allocate arrays
    nLayers = size(field,1)
    nCategories = size(field,2)
    allocate(zGrad(nLayers,nCategories))
    allocate(normalGrad(nLayers,nCategories,maxEdges))
    allocate(globalGrad(nLayers,nCategories,3))
    if (config_rotate_cartesian_grid) allocate(globalGradRotate(nLayers,nCategories,3))

    ! initialize the gradient
    xGrad(:,:,:) = 0.0_RKIND
    yGrad(:,:,:) = 0.0_RKIND
    zGrad(:,:) = 0.0_RKIND     ! diagnostic only

    ! loop over cells
    do iCell = 1, nCells

       if (maskCell(iCell) == 1) then  ! ice is present in the cell

          ! initialize gradient components
          normalGrad(:,:,:) = 0.0_RKIND
          globalGrad(:,:,:) = 0.0_RKIND

          ! loop over edges of this cell
          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             iCellNeighbor = cellsOnCell(iEdgeOnCell, iCell)

             ! compute the normal component of the gradient on this edge

             if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then ! there is a cell neighbor on this edge

                iEdge = edgesOnCell(iEdgeOnCell,iCell)

                do iCat = 1, nCategories
                   do iLayer = 1, nLayers

                      if (mask(iLayer,iCat,iCell) == 1 .and. mask(iLayer,iCat,iCellNeighbor) == 1) then

                         ! field values in cell and its neighbor have physical meaning; set the normal gradient
                         if (iCell == cellsOnEdge(1,iEdge)) then
                            signGradient =  1.0_RKIND
                         else
                            signGradient = -1.0_RKIND
                         end if

                         normalGrad(iLayer,iCat,iEdgeOnCell) = &
                              signGradient * (field(iLayer,iCat,iCellNeighbor) - field(iLayer,iCat,iCell)) / dcEdge(iEdge)

                      else   ! either or both field values do not have physical meaning; set gradient component = 0

                         normalGrad(iLayer,iCat,iEdgeOnCell) = 0.0_RKIND

                      endif

                   enddo   ! iCat
                enddo      ! iLayer

             else   ! there is no cell neighbor on this edge

                ! set gradient component = 0
                normalGrad(:,:,iEdgeOnCell) = 0.0_RKIND

             endif

             ! add the contribution of this normal component to the reconstructed gradient
             ! at the cell center (in global x/y/z coordinates)
             globalGrad(:,:,1) = globalGrad(:,:,1) + coeffsReconstruct(1,iEdgeOnCell,iCell) * normalGrad(:,:,iEdgeOnCell)
             globalGrad(:,:,2) = globalGrad(:,:,2) + coeffsReconstruct(2,iEdgeOnCell,iCell) * normalGrad(:,:,iEdgeOnCell)
             globalGrad(:,:,3) = globalGrad(:,:,3) + coeffsReconstruct(3,iEdgeOnCell,iCell) * normalGrad(:,:,iEdgeOnCell)

          enddo  ! iEdgeOnCell

          ! The global gradient vector lives in Cartesian x/y/z space.
          ! If on a rotated grid, then rotate the vector by 90 degrees, such that
          !  the true North Pole goes to (-1,0,0) and the true South Pole to (1,0,0).

          !TODO - Rotate the gradient vector if on a plane?

          if (config_rotate_cartesian_grid .and. on_a_sphere) then
             globalGradRotate(:,:,1) = -globalGrad(:,:,3)   ! xR = -z
             globalGradRotate(:,:,2) =  globalGrad(:,:,2)   ! yR = y
             globalGradRotate(:,:,3) =  globalGrad(:,:,1)   ! zR = x
             globalGrad(:,:,:) = globalGradRotate(:,:,:)
          endif

          ! transform from global x/y/z coordinates to local east/west coordinates

          if (on_a_sphere) then

             xGrad(:,:,iCell) = transGlobalToCell(1,1,iCell) * globalGrad(:,:,1)  &
                              + transGlobalToCell(1,2,iCell) * globalGrad(:,:,2)  &
                              + transGlobalToCell(1,3,iCell) * globalGrad(:,:,3)

             yGrad(:,:,iCell) = transGlobalToCell(2,1,iCell) * globalGrad(:,:,1)  &
                              + transGlobalToCell(2,2,iCell) * globalGrad(:,:,2)  &
                              + transGlobalToCell(2,3,iCell) * globalGrad(:,:,3)

             ! Note: The zGrad component is never used; it is simply computed as a diagnostic
             !TODO - Comment out zGrad computation?
             zGrad(:,:) = transGlobalToCell(3,1,iCell) * globalGrad(:,:,1)  &
                        + transGlobalToCell(3,2,iCell) * globalGrad(:,:,2)  &
                        + transGlobalToCell(3,3,iCell) * globalGrad(:,:,3)

          else  ! on a plane; do a simple copy

             xGrad(:,:,iCell) = globalGrad(:,:,1)
             yGrad(:,:,iCell) = globalGrad(:,:,2)

          endif

       endif   ! maskCell = 1

       if (verboseConstruct .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
          if (iCell == ctest) then
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'iLayer, iCat, iCell =', iLayerTest, iCatTest, indexToCellID(iCell)
             write(stderrUnit,*) 'field value =', field(iLayerTest,iCatTest,iCell)
             write(stderrUnit,*) 'neighbor values and edge gradients:'
             do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                iCellNeighbor = cellsOnCell(iEdgeOnCell,iCell)
                iEdge = edgesOnCell(iEdgeOnCell,iCell)
                write(stderrUnit,*) iEdgeOnCell, indexToEdgeID(iEdge), &
                     field(iLayerTest,iCatTest,iCellNeighbor), normalGrad(iLayerTest,iCatTest,iEdgeOnCell)
             enddo
             write(stderrUnit,*) ' '
             write(stderrUnit,*) 'Unlimited center gradient:', &
                  xGrad(iLayerTest,iCatTest,iCell), yGrad(iLayerTest,iCatTest,iCell), zGrad(iLayerTest,iCatTest)
          endif
       endif

    enddo  ! iCell

    ! cleanup
    deallocate(zGrad)
    deallocate(normalGrad)
    deallocate(globalGrad)
    if (config_rotate_cartesian_grid) deallocate(globalGradRotate)

  end subroutine compute_gradient_3d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine compute_barycenter_coordinates
!
!> \brief  compute the coordinates of the center of mass and related quantites
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine computes the x and y coordinates of a cell's center of mass,
!>  center of mass*tracer1, or center of mass*tracer1*tracer2.
!
!-----------------------------------------------------------------------

  subroutine compute_barycenter_coordinates(&
       geomAvgCell,                     &
       iCell,                           &
       xBarycenter,    yBarycenter,     &
       mean0, center0, xGrad0, yGrad0,  &
       mean1, center1, xGrad1, yGrad1,  &
       mean2, center2, xGrad2, yGrad2)

    type(geometric_avg_cell_type), pointer :: &
         geomAvgCell                       !< Input: derived type holding geometric averages

    integer, intent(in) :: iCell           !< Input: cell for which barycenter quantities are computed

    real(kind=RKIND), intent(out) ::  &
         xBarycenter, yBarycenter          !< Output: coordinates of barycenter relative to cell center

    real(kind=RKIND), intent(in), optional ::  &
         mean0, center0, xGrad0, yGrad0    !< Input:  mean value, center value, x/y gradients for mass-like quantity

    real(kind=RKIND), intent(in), optional ::  &
         mean1, center1, xGrad1, yGrad1    !< Input:  mean value, center value, x/y gradients for tracer 1

    real(kind=RKIND), intent(in), optional ::  &
         mean2, center2, xGrad2, yGrad2    !< Input:  mean value, center value, x/y gradients for tracer 2

    ! local variables

    real(kind=RKIND) ::  &
         c0,  &            ! coefficients in barycenter formulas
         cx,    cy,  &
         cxx,   cxy,   cyy,  &
         cxxx,  cxxy,  cxyy,  cyyy,   &
         cxxxx, cxxxy, cxxyy, cxyyy, cyyyy

    real(kind=RKIND) :: reciprocal, massTracerProd

    ! compute the barycenter coordinates, depending on how many tracer types were passed in
    ! Note: These formulas become more complex as the tracer hierarchy deepens.  See the model documentation for details.

    if (.not.present(mean0) .and. .not.present(mean1) .and. .not.present(mean2)) then

       ! barycenter = centroid
       xBarycenter = geomAvgCell % x(iCell)
       yBarycenter = geomAvgCell % y(iCell)

    elseif (.not.present(mean1) .and. .not.present(mean2)) then

       ! barycenter = center of mass
       c0 = center0
       cx = xGrad0
       cy = yGrad0
       if (abs(mean0) > 0.0_RKIND) then
          reciprocal = 1.0_RKIND / mean0
       else
          reciprocal = 0.0_RKIND
       endif

       xBarycenter = (c0 * geomAvgCell % x (iCell)  &
                    + cx * geomAvgCell % xx(iCell)  + cy * geomAvgCell % xy(iCell))  &
                    * reciprocal
       yBarycenter = (c0 * geomAvgCell % y (iCell)  &
                    + cx * geomAvgCell % xy(iCell)  + cy * geomAvgCell % yy(iCell))  &
                    * reciprocal

    elseif (.not.present(mean2)) then

       ! barycenter = center of mass*tracer1

       c0  = center0 * center1
       cx  = center0 * xGrad1  + xGrad0  * center1
       cy  = center0 * yGrad1  + yGrad0  * center1
       cxx = xGrad0  * xGrad1
       cxy = xGrad0  * yGrad1  + yGrad0  * xGrad1
       cyy = yGrad0  * yGrad1

       massTracerProd = mean0 * mean1
       if (abs(massTracerProd) > 0.0_RKIND) then
          reciprocal = 1.0_RKIND / massTracerProd
       else
          reciprocal = 0.0_RKIND
       endif

       xBarycenter = (c0  * geomAvgCell % x  (iCell) &
                    + cx  * geomAvgCell % xx (iCell) + cy  * geomAvgCell % xy (iCell)  &
                    + cxx * geomAvgCell % xxx(iCell) + cxy * geomAvgCell % xxy(iCell) + cyy * geomAvgCell % xyy(iCell)) &
                    * reciprocal
       yBarycenter = (c0  * geomAvgCell % y  (iCell) &
                    + cx  * geomAvgCell % xy (iCell) + cy  * geomAvgCell % yy (iCell)  &
                    + cxx * geomAvgCell % xxy(iCell) + cxy * geomAvgCell % xyy(iCell) + cyy * geomAvgCell % yyy(iCell)) &
                    * reciprocal

    else  ! mass field, tracer1 and tracer2 are all present

       ! barycenter = center of mass*tracer1*tracer2

       c0    = center0 * center1 * center2
       cx    = center0 * center1 * xGrad2   + center0 * xGrad1  * center2 + xGrad0  * center1 * center2
       cy    = center0 * center1 * yGrad2   + center0 * yGrad1  * center2 + yGrad0  * center1 * center2
       cxx   = center0 * xGrad1  * xGrad2   + xGrad0  * center1 * xGrad2  + xGrad0  * xGrad1  * center2
       cxy   = center0 * xGrad1  * yGrad2   + xGrad0  * yGrad1  * center2 + yGrad0  * center1 * xGrad2  &
             + center0 * yGrad1  * xGrad2   + xGrad0  * center1 * yGrad2  + yGrad0  * xGrad1  * center2
       cyy   = center0 * yGrad1  * yGrad2   + yGrad0  * center1 * yGrad2  + yGrad0  * yGrad1  * center2
       cxxx  = xGrad0  * xGrad1  * xGrad2
       cxxy  = xGrad0  * xGrad1  * yGrad2   + xGrad0  * yGrad1  * xGrad2  + yGrad0  * xGrad1  * xGrad2
       cxyy  = yGrad0  * yGrad1  * xGrad2   + yGrad0  * xGrad1  * yGrad2  + xGrad0  * yGrad1  * yGrad2
       cyyy  = yGrad0  * yGrad1  * yGrad2

       massTracerProd = mean0 * mean1 * mean2
       if (abs(massTracerProd) > 0.0_RKIND) then
          reciprocal = 1.0_RKIND / massTracerProd
       else
          reciprocal = 0.0_RKIND
       endif

       xBarycenter = (c0   * geomAvgCell % x   (iCell) &
                    + cx   * geomAvgCell % xx  (iCell) + cy   * geomAvgCell % xy  (iCell) &
                    + cxx  * geomAvgCell % xxx (iCell) + cxy  * geomAvgCell % xxy (iCell) + cyy  * geomAvgCell % xyy (iCell) &
                    + cxxx * geomAvgCell % xxxx(iCell) + cxxy * geomAvgCell % xxxy(iCell) + cxyy * geomAvgCell % xxyy(iCell) &
                                                                                          + cyyy * geomAvgCell % xyyy(iCell)) &
                    * reciprocal
       yBarycenter = (c0   * geomAvgCell % y   (iCell) &
                    + cx   * geomAvgCell % xy  (iCell) + cy   * geomAvgCell % yy  (iCell) &
                    + cxx  * geomAvgCell % xxy (iCell) + cxy  * geomAvgCell % xyy (iCell) + cyy  * geomAvgCell % yyy (iCell) &
                    + cxxx * geomAvgCell % xxxy(iCell) + cxxy * geomAvgCell % xxyy(iCell) + cxyy * geomAvgCell % xyyy(iCell) &
                                                                                          + cyyy * geomAvgCell % yyyy(iCell)) &
                    * reciprocal

    endif

  end subroutine compute_barycenter_coordinates

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine limit_tracer_gradient_2d
!
!> \brief  limit the gradient of a 2d tracer field
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine limits the gradient of a 2d tracer field so that the reconstructed
!>  value within a cell falls within the range of the values in the neighbor cells
!
!-----------------------------------------------------------------------

  subroutine limit_tracer_gradient_2d(&
       nCells,                       &
       nEdgesOnCell,                 &
       cellsOnCell,                  &
       xVertexOnCell, yVertexOnCell, &
       maskCell,                     &
       field,                        &
       fieldMask,                    &
       xBarycenter,   yBarycenter,   &
       xGrad,         yGrad)

    integer, intent(in) :: &
         nCells            !< Input: number of cells

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell      !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         cellsOnCell       !< Input: cell index for each edge neighbor of a given cell

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         field             !< Input: 2d field for which we are limiting the gradient

    integer, dimension(:,:), intent(in) ::  &
         fieldMask         !< Input: mask = 1 where field value is physically meaningful, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnCell,  & !< Input: x (east) coordinate of vertex relative to cell center in local tangent plane
         yVertexOnCell     !< Input: y (north) coordinate of vertex relative to cell center in local tangent plane

    integer, dimension(:), intent(in) :: &
         maskCell          !< Input: = 1 for cells with ice, else = 0

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xBarycenter,   &  !< Input: x (east) coordinate of barycenter of this tracer's parent
         yBarycenter       !< Input: y (north) coordinate of barycenter of this tracer's parent

    real(kind=RKIND), dimension(:,:), intent(inout) ::  &
         xGrad,         &  !< Input/output: x (east) component of gradient vector
         yGrad             !< Input/output: y (north) coordinate of barycenter of this tracer's parent
                           !  Gradient components are unlimited on input, limited on output

    ! local variables

    integer :: nCategories   ! number of categories

    integer :: iCell, iVertex, iCat, iEdgeOnCell  ! counting indices
    integer :: iCellNeighbor   ! index of neighbor cell

    real(kind=RKIND), dimension(:), allocatable ::  &
         maxNeighbor,    &    ! max value of the field in a cell and its nearest neighbors
         minNeighbor          ! min value of the field in a cell and its nearest neighbors

    real(kind=RKIND), dimension(:), allocatable ::  &
         maxLocal,    &       ! max reconstructed value of the field within a cell
         minLocal             ! min reconstructed value of the field within a cell

    real(kind=RKIND) :: &
         deviationAtVertex          ! tracer deviation (i.e., difference from mean value) at vertex, based on unlimited gradient

    real(kind=RKIND) ::  &
         gradFactor,     &          ! limiting factor for gradient, 0 < gradFactor < 1
         gradFactor1, gradFactor2   ! the two candidates for gradFactor, based on cell maxes and mins
                                    ! (the smaller value is chosen)

    nCategories = size(field, 1)

    ! allocate temporary arrays
    allocate(maxNeighbor(nCategories))
    allocate(minNeighbor(nCategories))
    allocate(maxLocal(nCategories))
    allocate(minLocal(nCategories))

    ! loop over cells
    do iCell = 1, nCells

       if (maskCell(iCell) == 1) then   ! ice is present

          ! compute max and min values of this tracer in the cell and its nearest neighbors

          ! initialize neighbor max and min to the value in iCell
          maxNeighbor(:) = field(:,iCell)
          minNeighbor(:) = field(:,iCell)

          ! loop over edges of the cell
          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             ! find the index of the neighbor cell
             iCellNeighbor = cellsOnCell(iEdgeOnCell,iCell)

             ! loop over categories
             do iCat = 1, nCategories

                ! modify the max and min, as appropriate
                if (fieldMask(iCat,iCellNeighbor) == 1) then
                   maxNeighbor(iCat) = max(maxNeighbor(iCat), field(iCat,iCellNeighbor))
                   minNeighbor(iCat) = min(minNeighbor(iCat), field(iCat,iCellNeighbor))
                endif

             enddo   ! iCat

          enddo      ! iEdgeOnCell

          ! convert the max and min to differences
          maxNeighbor(:) = maxNeighbor(:) - field(:,iCell)
          minNeighbor(:) = minNeighbor(:) - field(:,iCell)

          ! compute the max and min deviation of the reconstructed tracer within the cell

          ! initialize local max and min for each category
          !TODO - huge and -huge?
          maxLocal(:) = 0.0_RKIND
          minLocal(:) = 0.0_RKIND

          ! loop over vertices of this cell
          do iVertex = 1, nEdgesOnCell(iCell)

             ! loop over categories
             do iCat = 1, nCategories

                ! for each vertex, compute the deviation from the mean value
                !  Note: The barycenter (i.e., the location were the tracer has the value contained in the 'field' array)
                !         coincides with the geometric center only for mass-type fields (nParents = 0).
                !        This is why the x and y terms consist of the difference between
                !        the vertex coordinate and the barycenter coordinate.
                !        For nParents = 1, the barycenter is the center of mass.
                !        For nParents = 2, the barycenter is the center of mass*tracer1, etc.

                deviationAtVertex = xGrad(iCat,iCell) * (xVertexOnCell(iVertex,iCell) - xBarycenter(iCat,iCell)) &
                                  + yGrad(iCat,iCell) * (yVertexOnCell(iVertex,iCell) - yBarycenter(iCat,iCell))

                ! modify the max and min, as appropriate
                maxLocal(iCat) = max(maxLocal(iCat), deviationAtVertex)
                minLocal(iCat) = min(minLocal(iCat), deviationAtVertex)

             enddo  ! iCat

          enddo     ! iVertex

          ! compute the gradient limiting factor

          do iCat = 1, nCategories

             if (abs(maxLocal(iCat)) > abs(maxNeighbor(iCat))) then
                gradFactor1 = max (0.0_RKIND, maxNeighbor(iCat)/maxLocal(iCat))
             else
                gradFactor1 = 1.0_RKIND
             endif

             if (abs(minLocal(iCat)) > abs(minNeighbor(iCat))) then
                gradFactor2 = max (0.0_RKIND, minNeighbor(iCat)/minLocal(iCat))
             else
                gradFactor2 = 1.0_RKIND
             endif

             gradFactor = min(gradFactor1, gradFactor2)

             ! Reduce this factor slightly to avoid going out of bounds due to roundoff errors.
             ! This can happen with the current (Aug. 2015) geometry calculations when a departure trajectory
             !  and an edge are very nearly parallel. In this case their intersection is sometimes ignored, resulting
             !  in a small triangle lying very slightly outside the cell that the code thinks is its source cell.
             gradFactor = max(0.0_RKIND, gradFactor - eps11)

             ! limit the gradient components

             xGrad(iCat,iCell) = xGrad(iCat,iCell) * gradFactor
             yGrad(iCat,iCell) = yGrad(iCat,iCell) * gradFactor

          enddo   ! iCat

       endif     ! maskCell = 1

    enddo      ! iCell

    ! cleanup
    deallocate(maxNeighbor)
    deallocate(minNeighbor)
    deallocate(maxLocal)
    deallocate(minLocal)

  end subroutine limit_tracer_gradient_2d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine limit_tracer_gradient_3d
!
!> \brief  limit the gradient of a 3d tracer field
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine limits the gradient of a 3d tracer field so that the reconstructed
!>  value within a cell falls within the range of the values in the neighbor cells
!
!-----------------------------------------------------------------------

  subroutine limit_tracer_gradient_3d(&
       nCells,                       &
       nEdgesOnCell,                 &
       cellsOnCell,                  &
       xVertexOnCell, yVertexOnCell, &
       maskCell,                     &
       field,                        &
       fieldMask,                    &
       xGrad,         yGrad,         &
       xBarycenter2D, yBarycenter2D, &
       xBarycenter3D, yBarycenter3D)

    integer, intent(in) :: &
         nCells            !< Input: number of cells

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell      !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         cellsOnCell       !< Input: cell index for each edge neighbor of a given cell

    real(kind=RKIND), dimension(:,:,:), intent(in) ::  &
         field             !< Input: 3d field for which we are limiting the gradient

    integer, dimension(:,:,:), intent(in) ::  &
         fieldMask         !< Input: mask = 1 where field value is physically meaningful, = 0 elsewhere

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnCell,  & !< Input: x (east) coordinate of vertex relative to cell center in local tangent plane
         yVertexOnCell     !< Input: y (north) coordinate of vertex relative to cell center in local tangent plane

    integer, dimension(:), intent(in) :: &
         maskCell            !< Input: = 1 for cells with ice, else = 0

    real(kind=RKIND), dimension(:,:,:), intent(inout) ::  &
         xGrad,         &  !< Input/output: x (east) component of gradient vector
         yGrad             !< Input/output: y (north) coordinate of barycenter of this tracer's parent
                           !  Gradient components are unlimited on input, limited on output

    real(kind=RKIND), dimension(:,:), intent(in), optional ::  &
         xBarycenter2D, &  !< Input: x (east) coordinate of barycenter of this tracer's 2Dparent
         yBarycenter2D     !< Input: y (north) coordinate of barycenter of this tracer's 2D parent

    real(kind=RKIND), dimension(:,:,:), intent(in), optional ::  &
         xBarycenter3D, &  !< Input: x (east) coordinate of barycenter of this tracer's 3D parent
         yBarycenter3D     !< Input: y (north) coordinate of barycenter of this tracer's 3D parent

    ! local variables

    integer :: nLayers       ! number of layers
    integer :: nCategories   ! number of categories

    integer :: iCell, iVertex, iCat, iLayer, iEdgeOnCell  ! counting indices
    integer :: iCellNeighbor   ! index of neighbor cell

    real(kind=RKIND), dimension(:,:), allocatable ::  &
         maxNeighbor,    &    ! max value of the field in a cell and its nearest neighbors
         minNeighbor          ! min value of the field in a cell and its nearest neighbors

    real(kind=RKIND), dimension(:,:), allocatable ::  &
         maxLocal,    &       ! max reconstructed value of the field within a cell
         minLocal             ! min reconstructed value of the field within a cell

    real(kind=RKIND) :: &
         deviationAtVertex          ! tracer deviation (i.e., difference from mean value) at vertex, based on unlimited gradient

    real(kind=RKIND) ::  &
         gradFactor,     &          ! limiting factor for gradient, 0 < gradFactor < 1
         gradFactor1, gradFactor2   ! the two candidates for gradFactor, based on cell maxes and mins
                                    ! (the smaller value is chosen)

    ! Make sure one set of barycentric coordinates was passed in.
    ! Note: The parent tracer can be either 2D or 3D. Rather than write a separate subroutine for each case,
    !       I made the barycenter coordinates optional arguments. But these arguments are not strictly optional;
    !       one set or the other must be present.

    if (.not. (present(xBarycenter2D) .and. present(yBarycenter2D)) .and.  &
        .not. (present(xBarycenter3D) .and. present(yBarycenter3D))) then
       write(stderrUnit,*) 'ERROR: Subroutine limit_tracer_gradient_3d requires input barycenter coordinates'
       call mpas_dmpar_global_abort('MPAS-seaice: IR subroutine limit_tracer_gradient_3d requires input barycenter coordinates')
    endif

    nLayers     = size(field, 1)
    nCategories = size(field, 2)

    ! allocate temporary arrays
    allocate(maxNeighbor(nLayers,nCategories))
    allocate(minNeighbor(nLayers,nCategories))
    allocate(maxLocal(nLayers,nCategories))
    allocate(minLocal(nLayers,nCategories))

    ! loop over cells
    do iCell = 1, nCells

       if (maskCell(iCell) == 1) then  ! ice is present

          ! compute max and min values of this tracer in the cell and its nearest neighbors

          ! initialize neighbor max and min for each category
          maxNeighbor(:,:) = field(:,:,iCell)
          minNeighbor(:,:) = field(:,:,iCell)

          ! loop over edges of the cell
          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             ! find the index of the neighbor cell
             iCellNeighbor = cellsOnCell(iEdgeOnCell,iCell)

             ! loop over layers and categories
             do iCat = 1, nCategories
                do iLayer = 1, nLayers

                   ! modify the max and min, as appropriate
                   if (fieldMask(iLayer,iCat,iCellNeighbor) == 1) then
                      maxNeighbor(iLayer,iCat) = max(maxNeighbor(iLayer,iCat), field(iLayer,iCat,iCellNeighbor))
                      minNeighbor(iLayer,iCat) = min(minNeighbor(iLayer,iCat), field(iLayer,iCat,iCellNeighbor))
                   endif

                enddo   ! iCat
             enddo      ! iLayer

          enddo      ! iEdgeOnCell

          ! convert the max and min to differences
          maxNeighbor(:,:) = maxNeighbor(:,:) - field(:,:,iCell)
          minNeighbor(:,:) = minNeighbor(:,:) - field(:,:,iCell)

          ! compute the max and min deviation of the reconstructed tracer within the cell

          ! initialize local max and min for each category
          !TODO - huge and -huge?
          maxLocal(:,:) = 0.0_RKIND
          minLocal(:,:) = 0.0_RKIND

          if (present(xBarycenter2D) .and. present(yBarycenter2D)) then  ! parent tracer is 2D

             ! loop over vertices of this cell
             do iVertex = 1, nEdgesOnCell(iCell)

                ! loop over layers and categories
                do iCat = 1, nCategories
                   do iLayer = 1, nLayers

                      ! for each vertex, compute the deviation from the mean value
                      !  Note: The barycenter (i.e., the location were the tracer has the value contained in the 'field' array)
                      !         coincides with the geometric center only for mass-type fields (nParents = 0).
                      !        This is why the x and y terms consist of the difference
                      !        between the vertex coordinate and the barycenter coordinate.
                      !        For nParents = 1, the barycenter is the center of mass.
                      !        For nParents = 2, the barycenter is the center of mass*tracer1, etc.

                      deviationAtVertex = &
                           xGrad(iLayer,iCat,iCell) * (xVertexOnCell(iVertex,iCell) - xBarycenter2D(iCat,iCell)) &
                         + yGrad(iLayer,iCat,iCell) * (yVertexOnCell(iVertex,iCell) - yBarycenter2D(iCat,iCell))

                      ! modify the max and min, as appropriate
                      maxLocal(iLayer,iCat) = max(maxLocal(iLayer,iCat), deviationAtVertex)
                      minLocal(iLayer,iCat) = min(minLocal(iLayer,iCat), deviationAtVertex)

                   enddo  ! iCat
                enddo     ! iLayer

             enddo     ! iVertex

          elseif (present(xBarycenter3D) .and. present(yBarycenter3D)) then  ! parent tracer is 3D

             ! loop over vertices of this cell
             do iVertex = 1, nEdgesOnCell(iCell)

                ! loop over layers and categories
                do iCat = 1, nCategories
                   do iLayer = 1, nLayers

                      ! for each vertex, compute the deviation from the mean value
                      !  Note: The barycenter (i.e., the location were the tracer has the value contained in the 'field' array)
                      !         coincides with the geometric center only for mass-type fields (nParents = 0).
                      !        This is why the x and y terms consist of the difference
                      !        between the vertex coordinate and the barycenter coordinate.
                      !        For nParents = 1, the barycenter is the center of mass.
                      !        For nParents = 2, the barycenter is the center of mass*tracer1, etc.

                      deviationAtVertex = &
                           xGrad(iLayer,iCat,iCell) * (xVertexOnCell(iVertex,iCell) - xBarycenter3D(iLayer,iCat,iCell)) &
                         + yGrad(iLayer,iCat,iCell) * (yVertexOnCell(iVertex,iCell) - yBarycenter3D(iLayer,iCat,iCell))

                      ! modify the max and min, as appropriate
                      maxLocal(iLayer,iCat) = max(maxLocal(iLayer,iCat), deviationAtVertex)
                      minLocal(iLayer,iCat) = min(minLocal(iLayer,iCat), deviationAtVertex)

                   enddo  ! iCat
                enddo     ! iLayer

             enddo     ! iVertex

          endif   ! xBarycenter and yBarycenter are 2d or 3d

          ! compute the gradient limiting factor

          do iCat = 1, nCategories
             do iLayer = 1, nLayers

                if (abs(maxLocal(iLayer,iCat)) > abs(maxNeighbor(iLayer,iCat))) then
                   gradFactor1 = max (0.0_RKIND, maxNeighbor(iLayer,iCat)/maxLocal(iLayer,iCat))
                else
                   gradFactor1 = 1.0_RKIND
                endif

                if (abs(minLocal(iLayer,iCat)) > abs(minNeighbor(iLayer,iCat))) then
                   gradFactor2 = max (0.0_RKIND, minNeighbor(iLayer,iCat)/minLocal(iLayer,iCat))
                else
                   gradFactor2 = 1.0_RKIND
                endif

                gradFactor = min(gradFactor1, gradFactor2)

                ! Reduce this factor slightly to avoid going out of bounds due to roundoff errors.
                ! This can happen with the current (Aug. 2015) geometry calculations when a departure trajectory
                !  and an edge are very nearly parallel. In this case their intersection is sometimes ignored, resulting
                !  in a small triangle lying very slightly outside the cell that the code thinks is its source cell.

                gradFactor = max(0.0_RKIND, gradFactor - eps11)

                ! limit the gradient components

                xGrad(iLayer,iCat,iCell) = xGrad(iLayer,iCat,iCell) * gradFactor
                yGrad(iLayer,iCat,iCell) = yGrad(iLayer,iCat,iCell) * gradFactor

             enddo   ! iCat
          enddo      ! iLayer

       endif        ! maskCell = 1

    enddo         ! iCell

    ! cleanup
    deallocate(maxNeighbor)
    deallocate(minNeighbor)
    deallocate(maxLocal)
    deallocate(minLocal)

  end subroutine limit_tracer_gradient_3d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine find_departure_points
!
!> \brief  find departure points associated with each vertex
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  This routine traces backward trajectories to locate the departure
!>  points associated with each vertex.
!
!-----------------------------------------------------------------------

  subroutine find_departure_points(&
       domain,                              &
       nVertices,                           &
       dt,                                  &
       uVelocity,    vVelocity,             &
       minLengthEdgesOnVertex,              &
       indexToVertexID,                     &
       departurePoint,                      &
       on_a_sphere)

    use cice_diagnostics, only: &
         cice_warning_message

    type(domain_type), intent(in) :: &
         domain             !< Input: domain

    integer, intent(in) :: &
         nVertices          !< Input: number of vertices

    real(kind=RKIND) ::  &
         dt                 !< Input: time step

    real(kind=RKIND), dimension(:) ::   &
         uVelocity, vVelocity   !< Input: eastward and northward components of velocity field

    real(kind=RKIND), dimension(:) ::  &
         minLengthEdgesOnVertex !< Input: minimum length of the edges on each vertex

    integer, dimension(:), pointer ::  &
         indexToVertexID        !< Input: global index for each vertex

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         departurePoint         !< Output: x/y coordinates of departure points
                                !          1st index = 2, 2nd index = nVertices

    logical, intent(in) :: &
         on_a_sphere       !< Input: T if flow is on a sphere, F if on a plane

    ! local variables

    integer :: iVertex

    real(kind=RKIND), dimension(3) ::  &
         velocityGlobal   ! 3D velocity in global coordinates`

    real(kind=RKIND) :: trajectoryLength

    character(len=strKIND) :: &
         errorMessage ! error message for CFL warning

    ! initialize
    departurePoint(:,:) = 0.0_RKIND

    ! compute 3D velocity in global coordinates
    ! Note: The logic on a sphere is the same as on a plane.
    !TODO: Add an option to compute DPs more accurately using midpoint approximation

    ! loop over vertices
    do iVertex = 1, nVertices

       ! Note: Departure points are computed in vertex-based coordinates.
       !       I tried translating to global coordinates, from which they were
       !       further translated to edge-based coordinates.  But because of small
       !       projection differences, this led to negative areas.
       !       See the comments in subroutine get_geometry_incremental_remap.

       departurePoint(1,iVertex) = -uVelocity(iVertex) * dt
       departurePoint(2,iVertex) = -vVelocity(iVertex) * dt

    enddo   ! iVertex

    !TODO - The midpoint option could go here

    ! Check for potential CFL violations
    ! IR requires that the departure trajectory does not extend beyond the neighbor cells of a given edge
    ! (4 cellneighbors for hexes, 6 for quads).  Here, for each vertex we check whether the trajectory
    ! is longer than the shortest of the edges that meet at the vertex.  If so, a non-fatal warning is written.
    ! The code might run successfully, but there are no guarantees.

    do iVertex = 1, nVertices

       ! compute the length of the departure trajectory
       trajectoryLength = sqrt(departurePoint(1,iVertex)**2 + departurePoint(2,iVertex)**2)

       ! compare to the length of the shortest edge that meets at this vertex (computed at initialization)

       if (trajectoryLength > minLengthEdgesOnVertex(iVertex)) then
          write(errorMessage,*) 'Potential CFL violation in IR advection, global vertex ID =', indexToVertexID(iVertex)
          write(stderrUnit,*) 'WARNING: '//trim(errorMessage)
          call cice_warning_message(domain, errorMessage)
          write(errorMessage,*) 'Speed at vertex =', trajectoryLength/dt
          write(stderrUnit,*) trim(errorMessage)
          call cice_warning_message(domain, errorMessage)
          write(errorMessage,*) 'Maximum safe speed =', minLengthEdgesOnVertex(iVertex)/dt
          write(stderrUnit,*) trim(errorMessage)
          call cice_warning_message(domain, errorMessage)
       endif

    enddo

  end subroutine find_departure_points

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine find_departure_triangles
!
!> \brief  find the departure triangles associated with each edge
!> \author William Lipscomb
!> \date   June 2015
!> \details
!>  This routine computes the vertices of the triangles in the departure
!>  regions associated with each edge, along with the triangle areas
!>  and the cell where each triangle lies.
!
!-----------------------------------------------------------------------

  subroutine find_departure_triangles(&
       nEdges,                                    &
       nCells,                                    &
       maxEdges,                                  &
       nTriPerEdgeRemap,                          &
       nEdgesOnCell,                              &
       xVertexOnEdge,    yVertexOnEdge,           &
       xVertexOnCell,    yVertexOnCell,           &
       verticesOnCell,                            &
       verticesOnEdge,                            &
       edgesOnCell,                               &
       cellsOnEdge,                               &
       remapEdge,                                 &
       cellsOnEdgeRemap,                          &
       edgesOnEdgeRemap,                          &
       vertexDegree,                              &
       departurePointIn,                          &
       maskEdge,                                  &
       xTriangle,        yTriangle,               &
       iCellTriangle,                             &
       triangleArea,                              &
       on_a_sphere,                               &
       transGlobalToCell,                         &
       indexToCellID,                             &
       indexToEdgeID,                             &
       indexToVertexID,                           &
       block)

    !----------------------------------------------------------------------------------------------------------------!
    ! Here is a diagram of the IR geometry on a hex mesh (vertexDegree = 3).                                         !
    ! C denotes cells, V denotes vertices, E denotes edges, and D denotes departure points.                          !
    !                                                                                                                !
    ! Note the following:                                                                                            !
    ! (1) IR computes fluxes across edges. The main edge has vertices V1 and V2.                                     !
    ! (2) Cells C1 and C2 correspond to iCellOnEdge = 1 and 2, respectively. These cells share the main edge.        !
    !     C1 lies in the left half-plane of the vector pointing from vertex V1 to V2.                                !
    ! (3) There are two adjacent side cells: C3 and C4.                                                              !
    !     C3 includes V1, and C4 includes V2.                                                                        !
    ! (4) There are four adjacent edges: E1, E2, E3 and E4.                                                          !
    !     Edges E1 and E3 include V1, with E1 in the left half-plane and E3 in the right half-plane.                 !
    !     Similarly, edges E2 and E4 contain V2, with E2 in the left half-plane and E4 in the right half-plane.      !
    ! (5) Each adjacent edge includes one vertex (V1 or V2) on the main edge, and one vertex not on the main edge.   !
    !     Vertices V3 to V6 lie on edges E1 to E4. That is, the vertex index is equal to the edge index + 2.         !
    ! (6) Departure point D1 connects to V1, and departure point D2 connects to V2.                                  !
    !     D1 and D2 may lie in either half-plane. (In the diagram, both lie in the left half-plane.)                 !
    !     The line segment joining D1 or D2 may intersect one or more edges.                                         !
    ! (7) The quadrilateral defined by points V1-V2-D2-D1 can lie in up to 4 source cells                            !
    !     and can consist of up to 4 subtriangles (1 triangle in each side cell and up to 2 triangles                !
    !     in C1 and/or C2.                                                                                           !
    !                                                                                                                !
    !                      V3                         V4                                                             !
    !                        \                       /                                                               !
    !                         \   D1-----------D2   /                                                                !
    !                        E1\  |     C1     |   /E2                                                               !
    !                           \ |            |  /                                                                  !
    !                    C3      \V1___________V2/     C4                                                            !
    !                            /               \                                                                   !
    !                           /                 \                                                                  !
    !                        E3/        C2         \E4                                                               !
    !                         /                     \                                                                !
    !                        /                       \                                                               !
    !                      V5                         V6                                                             !
    !                                                                                                                !
    ! Suppose we have a quad mesh (vertexDegree = 4) instead of a hex mesh.                                          !
    ! Then the labeling is as follows.                                                                               !
    !                                                                                                                !
    !                             V3            V4                                                                   !
    !                             |             |                                                                    !
    !                         D1--|----------D2 |                                                                    !
    !                    C3    \  |E1    C1   \ |E2     C4                                                           !
    !                           \ |            \|                                                                    !
    !                 V7_________\|V1___________|V2_________V8                                                       !
    !                       E5    |             |     E6                                                             !
    !                             |             |                                                                    !
    !                    C5       |E3    C2     |E4     C6                                                           !
    !                             |             |                                                                    !
    !                             |             |                                                                    !
    !                             V5            V6                                                                   !
    !                                                                                                                !
    ! Note:                                                                                                          !
    ! (1) We have two additional cells, C3 and C5.                                                                   !
    !     C3 and C5 include V1, with C3 in the left half-plane of the edge and C5 in the right half-plane.           !
    !     Similarly, C4 and C6 include V2, with C4 in the left half-plane of the edge and C6 in the right half-plane.!
    ! (2) We have two additional edges, E5 and E6, that are (at least approximately) colinear with the main edge.    !
    !     E5 includes V1, and E6 includes V2.                                                                        !
    ! (3) We have two additional vertices, V7 and V8. These lie on edges E5 and E6, respectively.                    !
    ! (4) The maximum number of departure triangles is 6 instead of 4, if there is 1 triangle                        !
    !     in each of the 4 side cells in addition to the 2 triangles in C1 and/or C2.                                !
    !     (This can happen if E5 and E6 bend up or down away from the main edge.)                                    !
    !                                                                                                                !
    ! The values of nTriPerEdgeRemap, maxCellsPerEdgeRemap, maxEdgesPerEdgeRemap and maxVerticesPerEdgeRemap         !
    ! in the Registry (6, 6, 6 and 8 respectively) are based on the quad mesh.  The values for a hex mesh            !
    ! are 4, 4, 4 and 6.                                                                                             !
    ! TODO - Modify the setup so that these remap dimensions are based on vertexDegree?                              !
    !                                                                                                                !
    ! Most of the logic below is agnostic with respect to vertexDegree.                                              !
    ! The main difference is that on the hex mesh, side triangles lie in either C3 or C4, whereas                    !
    !  for the quad mesh these triangles can also lie in C5 or C6.                                                   !
    !                                                                                                                !
    !----------------------------------------------------------------------------------------------------------------!

    integer, intent(in) :: &
         nEdges,         &  !< Input: number of edges
         nCells,         &  !< Input: number of cells
         maxEdges,       &  !< Input: max number of edges per cell
         nTriPerEdgeRemap, &!< Input: max number of departure triangles per edge
         vertexDegree       !< Input: number of edges that meet at each vertex

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell       !< Input: number of edges per cell

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnEdge, yVertexOnEdge,  & !< Input: x/y coordinates of vertices in local edge-based coordinates
         xVertexOnCell, yVertexOnCell     !< Input: x/y coordinates of vertices in local cell-based coordinates

    integer, dimension(:), intent(in) ::  &
         remapEdge            !< Input: = 1 if IR fluxes need to be computed across an edge, else = 0
                              !        (i.e., = 1 for all edges of locally owned cells)

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell,      & !< Input: vertex index for each vertex of a given cell
         verticesOnEdge,      & !< Input: vertex index for each vertex of a given edge
         edgesOnCell,         & !< Input: edge index for each edge of a given cell
         cellsOnEdge,         & !< Input: index for each cell neighboring an edge
         cellsOnEdgeRemap,    & !< Input: index for each cell neighboring an edge
                                !  Note:  Unlike cellsOnEdge, this includes cells that share a single vertex with the edge
         edgesOnEdgeRemap       !< Input: index for each edges neighboring an edge
                                !  Note:  Unlike edgessOnEdge, this includes only edges that share a vertex with the edge

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         departurePointIn       !< Input: x/y coordinates of departure points relative to vertices
                                !         1st index = 2, 2nd index = nVertices

    integer, dimension(:), intent(out) :: &
         maskEdge               !< Output: = 1 for edges with potential nonzero fluxes, else = 0

    real(kind=RKIND), dimension(:,:,:), intent(out) ::  &
         xTriangle, yTriangle   !< Output: x and y coordinates of vertices of departure triangles,
                                !          relative to the center of the cell where they are located
                                !          1st index = nQuadPoints, 2nd index = nTriPerEdgeRemap, 3rd index = nEdges

    real(kind=RKIND), dimension(:,:), intent(out) ::  &
         triangleArea           !< Output: area of departure triangles
                                !          positive if the triangle moves across the edge from cell 1 to cell 2; else negative

    integer, dimension(:,:), intent(out) :: &
         iCellTriangle          !< Output: index of the cell containing each departure triangle

    logical, intent(in) :: &
         on_a_sphere            !< Input: T if flow is on a sphere, F if on a plane

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         transGlobalToCell      !< Input: 3x3 matrix for transforming vectors from global to cell-based coordinates

    integer, dimension(:), intent(in) ::  &
         indexToCellID,   & !< Input: global index for each cell on a block (diagnostic only)
         indexToVertexID, & !< Input: global index for each vertex on a block (diagnostic only)
         indexToEdgeID      !< Input: global index for each edge on a block (diagnostic only)

    type(block_type), intent(in) :: &
         block              !< Input: local block (diagnostic only)

    ! local variables

    integer :: &
         triangleCount          ! running total of the number of departure triangles for an edge

    logical :: &
         edgeIntersect,        &! true if the line segment joining the departure points intersects a side edge
         edgeIntersectMain      ! true if the line segment joining the departure points intersects the main edge

    real(kind=RKIND), dimension(2) ::  &
         intersectionPoint,   & ! x and y coordinates of the intersection of two lines
                                ! one line is a side edge; the other line joins the two departure points on the edge
         intersectionPointMain  ! x and y coordinates of the intersection of two lines
                                ! one line is the main (central) edge or its extension;
                                ! the other line joins the two departure points on the edge

    integer, dimension(nTriPerEdgeRemap) ::  &
         triangleVertexOnEdge, &! vertexOnEdge index for a given triangle, relative to iEdge
         triangleVertexOnCell, &! vertexOnCell index for a given triangle, relative to iCell
         fluxSign               ! sign of the flux across a given edge
                                ! (positive if going from cell 1 to cell 2, negative if in the other direction)

    real(kind=RKIND), dimension(2,2) ::  &
         edgeVertex,         & ! x/y coordinates of 2 vertices of edge
         edgeNeighborVertex    ! x/y coordinates of 2 vertices of neighboring edge

    real(kind=RKIND), dimension(2,nTriPerEdgeRemap) ::  &
         edgeVector1,        & ! 2D vectors pointing along adjacent edges; used for triangle coordinate transformations
         edgeVector2

    real(kind=RKIND), dimension(2,2) ::  &
         departurePoint        ! x/y coordinates of 2 departure points relative to an edge midpoint
                               ! Note: The input vector, departurePointIn, points from a vertex to the DP.

    real(kind=RKIND) :: quadArea, lenSquared

    logical, dimension(2) ::  &
         dpInHalfPlane         ! true if a departure point lies in the left half-plane of the edge

    integer :: iEdge, iCell, iVertex, iVertexOnCell, iEdgeOnCell, iVertexOnEdge, iTri, iEdgeNeighbor
    integer :: iCellOnEdgeRemap, iEdgeOnEdgeRemap, iVertexOnEdgeRemap
    integer :: iVertex1, iVertex2, iSideIndex, iOtherEdge, iOtherVertex
    integer :: n

    if (verboseGeometry) then
       write(stderrUnit,*) 'In find_departure_triangles'
       write(stderrUnit,*) 'nCells, maxEdges =', nCells, maxEdges

       if (etestOnProc .and. block % localBlockID == etestBlockID) then
          iEdge = etest
          write(stderrUnit,*) 'cellsOnEdgeRemap, test edge =', indexToEdgeID(iEdge)
          do n = 1, 6
             iCell = cellsOnEdgeRemap(n,iEdge)
             if (iCell >= 1) then
                write(stderrUnit,*) indexToCellID(iCell)
             endif
          enddo
       endif
    endif   ! verbose

    ! initialize
    xTriangle(:,:,:) = 0.0_RKIND
    yTriangle(:,:,:) = 0.0_RKIND
    triangleArea(:,:) = 0.0_RKIND
    iCellTriangle(:,:) = 0
    fluxSign(:) = 0
    edgeVector1(:,:) = 0.0_RKIND
    edgeVector2(:,:) = 0.0_RKIND
    triangleVertexOnEdge(:) = 0
    triangleVertexOnCell(:) = 0

    ! Compute a mask of edges that potentially have nonzero fluxes.
    ! Must have (1) remapEdge = 1 and (2) departurePointIn has nonzero length for at least one vertex.

    maskEdge(:) = 0
    do iEdge = 1, nEdges
       if (remapEdge(iEdge) == 1) then
          do iVertexOnEdge = 1, 2
             iVertex = verticesOnEdge(iVertexOnEdge,iEdge)
             lenSquared = departurePointIn(1,iVertex)**2 + departurePointIn(2,iVertex)**2
             if (lenSquared > 0.0_RKIND) maskEdge(iEdge) = 1
          enddo  ! iVertexOnEdge
       endif     ! remapEdge = 1
    enddo        ! iEdge

    ! loop over edges
    ! Note: The loop is over nEdges, which (unlike nEdgesSolve) includes all edges of locally owned cells.

    do iEdge = 1, nEdges
       if (maskEdge(iEdge) == 1) then   !TODO - Modify indentation from here to the end of the loop?

          ! Compute the coordinates of the departure points in local edge-based coordinates
          ! Note: The logic is the same on a sphere as on a plane, because in either case the departurePoint
          !        vector uses vertex-based coordinates.
          !       I tried a different method where the departurePoint vector was translated to global coordinates
          !        and then edge-based coordinates, but this led to negative areas in some cases.
          !        See comments in subroutine get_geometry_incremental_remap.

          if (verboseGeometry .and. etestOnProc .and. block % localBlockID == etestBlockID) then
             if (iEdge == etest) then
                write(stderrUnit,*) ' '
                write(stderrUnit,*) 'Departure points, iEdge =', indexToEdgeID(iEdge)
             endif
          endif

          do iVertexOnEdge = 1, 2

             iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

             ! Add the departurePoint vector computed above (in vertex-based coordinates) to x/yVertexOnEdge.
             ! The resulting vector points from the edge midpoint to the DP.
             departurePoint(1,iVertexOnEdge) = xVertexOnEdge(iVertexOnEdge,iEdge) + departurePointIn(1,iVertex)
             departurePoint(2,iVertexOnEdge) = yVertexOnEdge(iVertexOnEdge,iEdge) + departurePointIn(2,iVertex)

             if (verboseGeometry .and. etestOnProc .and. block % localBlockID == etestBlockID) then
                if (iEdge == etest) then
                   write(stderrUnit,*) indexToVertexID(iVertex), departurePoint(1,iVertexOnEdge), departurePoint(2,iVertexOnEdge)
                endif
             endif

          enddo   ! iVertexOnEdge

          ! Set the coordinates of the 2 vertices on iEdge

          iVertex1 = verticesOnEdge(1,iEdge)
          iVertex2 = verticesOnEdge(2,iEdge)

          edgeVertex(1,1) = xVertexOnEdge(1,iEdge)
          edgeVertex(2,1) = yVertexOnEdge(1,iEdge)

          edgeVertex(1,2) = xVertexOnEdge(2,iEdge)
          edgeVertex(2,2) = yVertexOnEdge(2,iEdge)

          ! Determine whether each departure point lies in the left half-plane of the edge.
          !
          !  The possibilities are:
          ! (1) Both departure points lie in the left half-plane.
          !     The departure region is a convex quadrilateral.
          !     The flux across the edge is from cell 1 to cell 2 and is defined to be positive.
          ! (2) Neither departure point lies in the left half-plane.
          !     The departure region is a convex quadrilateral.
          !     The flux across the edge is from cell 2 to cell 1 and is defined to be negative.
          ! (3) One departure point lies in the left half-plane, and the other does not.
          !     The departure region consists of two triangles, one in the half-plane but not the other.
          !     (a) The line segment connecting D1 and D2 intersects the edge. The triangle within the left half-plane
          !         contributes a positive flux; the other triangle contributes a negative flux.
          !     (b) The line segment connecting D1 and D2 does not intersect the edge. Both triangles contribute
          !         a flux of the same sign.  The flux is positive if the midpoint between D1 and D2
          !         lies in the left half-plane, else the flux is negative.
          !     Note: In case 3b the departure region is a quadrilateral, but is not necessarily convex.

          do iVertexOnEdge = 1,2
             dpInHalfPlane(iVertexOnEdge) = point_in_half_plane(edgeVertex(:,1), edgeVertex(:,2), departurePoint(:,iVertexOnEdge))
          enddo

          if (verboseGeometry .and. etestOnProc .and. block % localBlockID == etestBlockID) then
             if (iEdge == etest) then
                write(stderrUnit,*) ' '
                write(stderrUnit,*) 'edge vertex 1:', edgeVertex(:,1)
                write(stderrUnit,*) 'edge vertex 2:', edgeVertex(:,2)
                write(stderrUnit,*) 'Distance 1:', &
                     sqrt( (edgeVertex(1,1) - departurePoint(1,1))**2 + (edgeVertex(2,1) - departurePoint(2,1))**2 )
                write(stderrUnit,*) 'Distance 2:', &
                     sqrt( (edgeVertex(1,2) - departurePoint(1,2))**2 + (edgeVertex(2,2) - departurePoint(2,2))**2 )
                write(stderrUnit,*) 'Length of edge:', &
                     sqrt( (edgeVertex(1,2) - edgeVertex(1,1))**2 + (edgeVertex(2,2) - edgeVertex(2,1))**2 )
                write(stderrUnit,*) 'Length of D12:', &
                     sqrt( (departurePoint(1,2) - departurePoint(1,1))**2 + (departurePoint(2,2) - departurePoint(2,1))**2 )
                write(stderrUnit,*) 'dp1 in left half plane = ', dpInHalfPlane(1)
                write(stderrUnit,*) 'dp2 in left half plane = ', dpInHalfPlane(2)
             endif
          endif

          ! Compute the vertices and area of the departure triangles across this edge

          ! initialize the number of flux triangles for this edge
          triangleCount = 0

          ! Loop over the 4 side edges (E1, E2, E3 and E4).
          ! For each edge, determine whether the line segment joining D1 and D2 intersects the edge.
          ! If so, then add a side triangle.

          do iVertexOnEdge = 1, 2

             edgeNeighborVertex(:,1) = edgeVertex(:,iVertexOnEdge)   ! one vertex shared with main edge

             do iSideIndex = 0, 1   ! 0 = left half-plane, 1 = right half-plane

                ! choose one of the 4 side edges
                iEdgeOnEdgeRemap = iVertexOnEdge + 2*iSideIndex   ! E1 and E3 for V1; E2 and E4 for V2

                ! save the vertex of this side edge for later reference
                iVertexOnEdgeRemap = iEdgeOnEdgeRemap + 2  ! V3 and V5 for edges E1 and E3; V4 and V6 for edges V2 and V4

                ! if the side edge exists, then compute its intersection with the line segment joining D1 and D2

                iEdgeNeighbor = edgesOnEdgeRemap(iEdgeOnEdgeRemap,iEdge)

                if (iEdgeNeighbor >= 1 .and. iEdgeNeighbor <= nEdges) then  ! the side edge exists

                   edgeNeighborVertex(1,2) = xVertexOnEdge(iVertexOnEdgeRemap,iEdge)
                   edgeNeighborVertex(2,2) = yVertexOnEdge(iVertexOnEdgeRemap,iEdge)

                   ! find whether the line segment joining D1 and D2 intersects the edge
                   call find_line_intersection(&
                        departurePoint(:,1),     departurePoint(:,2),      &
                        edgeNeighborVertex(:,1), edgeNeighborVertex(:,2),  &
                        edgeIntersect,                                     &
                        intersectionPoint)

                else   ! the side edge does not exist

                   edgeIntersect = .false.

                endif

                if (edgeIntersect) then

                   ! add a triangle with vertices V1/D1/IP or V2/D2/IP
                   triangleCount = triangleCount + 1
                   xTriangle(1,triangleCount,iEdge) = edgeVertex(1,iVertexOnEdge)      ! V1 or V2
                   yTriangle(1,triangleCount,iEdge) = edgeVertex(2,iVertexOnEdge)
                   xTriangle(2,triangleCount,iEdge) = departurePoint(1,iVertexOnEdge)  ! D1 or D2
                   yTriangle(2,triangleCount,iEdge) = departurePoint(2,iVertexOnEdge)
                   xTriangle(3,triangleCount,iEdge) = intersectionPoint(1)
                   yTriangle(3,triangleCount,iEdge) = intersectionPoint(2)

                   ! identify the vertexOnEdge value (1 or 2) for this triangle
                   triangleVertexOnEdge(triangleCount) = iVertexOnEdge

                   ! identify the source cell where the triangle is located
                   iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(iVertexOnEdge+2,iEdge)  ! C3 for V1; C4 for V2

                   ! identify the vertexOnCell value (1 to nEdgesOnCell) for this triangle
                   iCell = iCellTriangle(triangleCount,iEdge)
                   do iVertexOnCell = 1, nEdgesOnCell(iCell)
                      if (verticesOnEdge(iVertexOnEdge,iEdge) == verticesOnCell(iVertexOnCell,iCell)) then
                         ! this is the vertex shared by iCell and iEdge
                         triangleVertexOnCell(triangleCount) = iVertexOnCell
                      endif
                   enddo

                   ! set the vectors that will be used as a basis for triangle vertices (spherical mesh only)

                   if (on_a_sphere) then

                      ! first the edge on which the intersection point lies
                      edgeVector1(:,triangleCount) = edgeNeighborVertex(:,2) - edgeNeighborVertex(:,1)

                      ! now the other side edge at this vertex
                      ! For vertexDegree = 3, E3 is paired with E1, and E4 is paired with E2.
                      ! For vertexDegree = 4, either E1 or E3 is paired with E5, and either E2 or E4 is paired with E6.

                      if (vertexDegree == 3) then  ! hex mesh
                         iOtherEdge = iEdgeOnEdgeRemap + 2
                         if (iOtherEdge > 4) iOtherEdge = iOtherEdge - 4 ! 4 side edges in total: E1 and E3 on V1, E2 and E4 on V2
                      elseif (vertexDegree == 4) then  ! quad mesh
                         iOtherEdge = iVertexOnEdge + 4  ! E5 for V1; E6 for V2
                      endif  ! vertexDegree

                      iOtherVertex = iOtherEdge + 2 ! V3 and V5 for edges E1 and E3;
                                                    ! V4 and V6 for E2 and E4; V7 and V8 for E5 and E6
                      edgeNeighborVertex(1,2) = xVertexOnEdge(iOtherVertex,iEdge)
                      edgeNeighborVertex(2,2) = yVertexOnEdge(iOtherVertex,iEdge)
                      edgeVector2(:,triangleCount) = edgeNeighborVertex(:,2) - edgeNeighborVertex(:,1)

                   endif  ! on a sphere

                   ! set the sign of the flux
                   if (iSideIndex == 0) then  ! IP in left half-plane
                      fluxSign(triangleCount) = 1
                   else   ! IP in right half-plane
                      fluxSign(triangleCount) = -1
                   endif

                   ! If vertexDegree = 4 (quad mesh), then check whether the line joining D1 and D2
                   ! intersects E5 or E6.  If so, break the triangle above into 2 triangles, one each in C3 and C5.

                   if (vertexDegree == 4) then

                      ! if E5 or E6 exists, then compute its intersection with the line segment joining D1 and D2

                      iEdgeOnEdgeRemap = iVertexOnEdge + 4  ! E5 or E6
                      iEdgeNeighbor = edgesOnEdgeRemap(iEdgeOnEdgeRemap,iEdge)
                      if (iEdgeNeighbor >= 1 .and. iEdgeNeighbor <= nEdges) then  ! E5 or E6 exists

                         ! set edgeNeighborVertex(:,2) to V7 or V8; edgeNeighborVertex(:,1) is unchanged
                         edgeNeighborVertex(1,2) = xVertexOnEdge(iVertexOnEdge+6,iEdge)
                         edgeNeighborVertex(2,2) = yVertexOnEdge(iVertexOnEdge+6,iEdge)

                         ! find whether the line segment joining D1 and D2 intersects E5 or E6
                         ! (not actually the main edge, but colinear with the main edge)
                         call find_line_intersection(&
                              departurePoint(:,1),     departurePoint(:,2),      &
                              edgeNeighborVertex(:,1), edgeNeighborVertex(:,2),  &
                              edgeIntersectMain,                                 &
                              intersectionPointMain)

                      else  ! E5 or E6 does not exist

                         edgeIntersectMain = .false.

                      endif

                      if (edgeIntersectMain) then

                         ! change one vertex of the triangle just computed; replace IP with IP0
                         xTriangle(3,triangleCount,iEdge) = intersectionPointMain(1)
                         yTriangle(3,triangleCount,iEdge) = intersectionPointMain(2)

                         ! change the source cell if necessary
                         if (iSideIndex == 0) then  ! IP in left half-plane, so this triangle is in right half-plane
                            iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(iVertexOnEdge+4,iEdge)  ! C5 or C6

                            ! change triangleVertexOnCell
                            iCell = iCellTriangle(triangleCount,iEdge)
                            do iVertexOnCell = 1, nEdgesOnCell(iCell)
                               if (verticesOnEdge(iVertexOnEdge,iEdge) == verticesOnCell(iVertexOnCell,iCell)) then
                                  ! this is the vertex shared by iCell and iEdge
                                  triangleVertexOnCell(triangleCount) = iVertexOnCell
                               endif
                            enddo

                            if (on_a_sphere) then

                               ! change edgeVector1 (spherical mesh only)
                               ! Currently, iVertexOnEdgeRemap = V3 or V4 (since IP is in left half-plane);
                               ! switch to V5 or V6 in right half-plane
                               ! Since edgeVector2 already points along E5 or E6, it need not be changed

                               iOtherVertex = iVertexOnEdgeRemap + 2    ! V3 -> V5; V4 -> V6
                               edgeNeighborVertex(1,2) = xVertexOnEdge(iOtherVertex,iEdge)
                               edgeNeighborVertex(2,2) = yVertexOnEdge(iOtherVertex,iEdge)
                               edgeVector1(:,triangleCount) = edgeNeighborVertex(:,2) - edgeNeighborVertex(:,1)

                            endif

                         else    ! iSideIndex = 1; IP in right half-plane, so this triangle is in left half-plane
                                 ! iCellTriangle is already correct (C3 or C4)
                                 ! Thus triangleVertexOnCell is also correct

                            if (on_a_sphere) then

                               ! change edgeVector1
                               ! Currently, iVertexOnEdgeRemap = V5 or V6 (since IP is in right half-plane);
                               ! switch to V3 or V4 in right half-plane
                               ! Since edgeVector2 already points along E5 or E6, it need not be changed

                               iOtherVertex = iVertexOnEdgeRemap - 2    ! V5 -> V3; V6 -> V4
                               edgeNeighborVertex(1,2) = xVertexOnEdge(iOtherVertex,iEdge)
                               edgeNeighborVertex(2,2) = yVertexOnEdge(iOtherVertex,iEdge)
                               edgeVector1(:,triangleCount) = edgeNeighborVertex(:,2) - edgeNeighborVertex(:,1)

                            endif  ! on a sphere

                         endif  ! iSideIndex

                         ! add a triangle with vertices V1/V2, IP0 and IP
                         triangleCount = triangleCount + 1
                         xTriangle(1,triangleCount,iEdge) = edgeVertex(1,iVertexOnEdge)    ! V1 or V2
                         yTriangle(1,triangleCount,iEdge) = edgeVertex(2,iVertexOnEdge)
                         xTriangle(2,triangleCount,iEdge) = intersectionPointMain(1)
                         yTriangle(2,triangleCount,iEdge) = intersectionPointMain(2)
                         xTriangle(3,triangleCount,iEdge) = intersectionPoint(1)
                         yTriangle(3,triangleCount,iEdge) = intersectionPoint(2)

                         ! identify the source cell
                         if (iSideIndex == 0) then  ! IP in left half-plane, where this triangle is located
                            iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(iVertexOnEdge+2,iEdge)  ! C3 or C4
                         else    ! iSideIndex = 1; IP in right half-plane, where this triangle is located
                            iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(iVertexOnEdge+4,iEdge)  ! C5 or C6
                         endif

                         ! set vertexOnEdge (same as for the previous triangle)
                         triangleVertexOnEdge(triangleCount) = triangleVertexOnEdge(triangleCount-1)

                         ! set vertexOnCell (different from previous triangle, since source cell is different)
                         iCell = iCellTriangle(triangleCount,iEdge)
                         do iVertexOnCell = 1, nEdgesOnCell(iCell)
                            if (verticesOnEdge(iVertexOnEdge,iEdge) == verticesOnCell(iVertexOnCell,iCell)) then
                               ! this is the vertex shared by iCell and iEdge
                               triangleVertexOnCell(triangleCount) = iVertexOnCell
                            endif
                         enddo

                         if (on_a_sphere) then

                            ! set the vectors that will be used as a basis for triangle vertices
                            ! edgeVector1 should point toward iVertexOnEdgeRemap, which was set above (V3, V4, V5 or V6)
                            edgeNeighborVertex(1,2) = xVertexOnEdge(iVertexOnEdgeRemap,iEdge)
                            edgeNeighborVertex(2,2) = yVertexOnEdge(iVertexOnEdgeRemap,iEdge)
                            edgeVector1(:,triangleCount) = edgeNeighborVertex(:,2) - edgeNeighborVertex(:,1)

                            ! edgeVector2 is either E5 or E6, as for the previous triangle
                            edgeVector2(:,triangleCount) = edgeVector2(:,triangleCount-1)

                         endif

                         ! set the sign of the flux (same as for the previous triangle)
                         fluxSign(triangleCount) = fluxSign(triangleCount-1)

                      else    ! no intersection with E5 or E6, but may have to change the source cell to C5 or C6

                         if (iSideIndex == 0) then  ! IP in left half-plane, where this triangle is located

                            ! iCellTriangle is already correct (C3 or C4); do nothing

                         else   ! IP in right half-plane, where this triangle is located

                            ! change the source cell
                            iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(iVertexOnEdge+4,iEdge)  ! C5 for V1; C6 for V2

                            ! change triangleVertexOnCell
                            iCell = iCellTriangle(triangleCount,iEdge)
                            do iVertexOnCell = 1, nEdgesOnCell(iCell)
                               if (verticesOnEdge(iVertexOnEdge,iEdge) == verticesOnCell(iVertexOnCell,iCell)) then
                                  ! this is the vertex shared by iCell and iEdge
                                  triangleVertexOnCell(triangleCount) = iVertexOnCell
                               endif
                            enddo

                            ! Note: edgeVectors do not need to be changed
                            !       edgeVector1 already points toward iVertexOnEdgeRemap (V3, V4, V5 or V6)
                            !       edgeVector2 points toward V7 or V8

                         endif   ! iSideIndex

                      endif   ! edgeIntersectMain (E5 or E6)

                   endif   ! vertexDegree = 4 (quad mesh)

                   ! Switch the departure point to the side intersection point
                   ! (for the purpose of computing triangle vertices in C1 and C2 below)
                   departurePoint(:,iVertexOnEdge) = intersectionPoint(:)

                endif   ! edgeIntersect (E1, E2, E3 or E4)

             enddo   ! iSideIndex

          enddo   ! iVertexOnEdge

          ! Next compute vertices of the departure triangles in C1 and/or C2
          ! Either there is a departure quad on one side of the edge, which is split into two triangles,
          !  or there is a departure triangle on each side of the edge.

          ! find whether the line segment joining D1 and D2 intersects the main edge
          call find_line_intersection(&
               departurePoint(:,1), departurePoint(:,2),  &
               edgeVertex(:,1),     edgeVertex(:,2),      &
               edgeIntersectMain,                         &
               intersectionPointMain)

          if (edgeIntersectMain) then   !  one departure triangle on each side of the main edge

             do iVertexOnEdge = 1, 2

                ! add triangle associated with this departure point; vertices are V1/V2, D1/D2 and IP0
                triangleCount = triangleCount + 1
                xTriangle(1,triangleCount,iEdge) = edgeVertex(1,iVertexOnEdge)      ! V1 or V2
                yTriangle(1,triangleCount,iEdge) = edgeVertex(2,iVertexOnEdge)
                xTriangle(2,triangleCount,iEdge) = departurePoint(1,iVertexOnEdge)  ! D1 or D2
                yTriangle(2,triangleCount,iEdge) = departurePoint(2,iVertexOnEdge)
                xTriangle(3,triangleCount,iEdge) = intersectionPointMain(1)
                yTriangle(3,triangleCount,iEdge) = intersectionPointMain(2)

                ! identify the vertexOnEdge value (1 or 2) for this triangle
                ! Although the intersection point may lie much closer to one edge than the other,
                !  the vertexOnEdge value is chosen based on iVertexOnEdge.
                triangleVertexOnEdge(triangleCount) = iVertexOnEdge

                ! choose the cell location and set the sign of the flux
                dpInHalfPlane(iVertexOnEdge) = &
                     point_in_half_plane(edgeVertex(:,1), edgeVertex(:,2), departurePoint(:,iVertexOnEdge))
                if (dpInHalfPlane(iVertexOnEdge)) then  ! departure point in left half-plane of iEdge
                   iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(1,iEdge)  ! C1
                   fluxSign(triangleCount) = 1
                else   ! departure point in right half-plane
                   iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(2,iEdge)  ! C2
                   fluxSign(triangleCount) = -1
                endif

                ! identify the vertexOnCell value (1 to nEdgesOnCell) for this triangle
                iCell = iCellTriangle(triangleCount,iEdge)
                do iVertexOnCell = 1, nEdgesOnCell(iCell)
                   if (verticesOnEdge(iVertexOnEdge,iEdge) == verticesOnCell(iVertexOnCell,iCell)) then
                      ! this is the vertex shared by iCell and iEdge
                      triangleVertexOnCell(triangleCount) = iVertexOnCell
                   endif
                enddo

                if (on_a_sphere) then  ! set the vectors that will be used as a basis for triangle vertices

                   ! first the main edge that on which the intersection point lies
                   if (iVertexOnEdge == 1) then
                      iOtherVertex = 2
                   else
                      iOtherVertex = 1
                   endif
                   edgeVector1(:,triangleCount) = edgeVertex(:,iOtherVertex) - edgeVertex(:,iVertexOnEdge)

                   ! now the side edge at this vertex
                   ! Choose E1 or E3 if iVertexOnEdge = 1, else choose E2 or E4

                   if (dpInHalfPlane(iVertexOnEdge)) then  ! departure point in left half-plane of iEdge; choose E1 or E2
                      iOtherEdge = iVertexOnEdge
                   else  ! departure point in right half-plane; choose E3 or E4
                      iOtherEdge = iVertexOnEdge + 2
                   endif
                   edgeNeighborVertex(1,2) = &
                        xVertexOnEdge(iOtherEdge+2,iEdge)   ! V3 and V5 for edges E1 and E3; V4 and V6 for edges E2 and E4
                   edgeNeighborVertex(2,2) = yVertexOnEdge(iOtherEdge+2,iEdge)
                   edgeVector2(:,triangleCount) = edgeNeighborVertex(:,2) - edgeVertex(:,iVertexOnEdge)

                endif  ! on a sphere

             enddo

          else  ! no intersection with main edge; departure quad on one side of the edge

             ! There is a possible degenerate case where both departure trajectories lie along the main edge,
             !  resulting in a departure quad with zero area.  Check for this degenerate case.
             ! Note: You might think such a triangle would be harmless. But it is not, because its area can
             !       become slightly greater than zero on transformation into cell-based coordinates.
             !       Then if either or both DPs lies outside the central cell, it can be given a value beyond
             !       the range of the gradient limiter, leading to negative mass.

             call quadrilateral_area(edgeVertex(:,1), edgeVertex(:,2), departurePoint(:,2), departurePoint(:,1), quadArea)

             if (quadArea > 0.0_RKIND) then  !TODO - change minimum quad area to eps11 or eps11**2?

                ! add triangle with vertices V1, V2 and D1
                triangleCount = triangleCount + 1
                xTriangle(1,triangleCount,iEdge) = edgeVertex(1,1)      ! V1
                yTriangle(1,triangleCount,iEdge) = edgeVertex(2,1)
                xTriangle(2,triangleCount,iEdge) = edgeVertex(1,2)      ! V2
                yTriangle(2,triangleCount,iEdge) = edgeVertex(2,2)
                xTriangle(3,triangleCount,iEdge) = departurePoint(1,1)  ! D1
                yTriangle(3,triangleCount,iEdge) = departurePoint(2,1)

                ! identify the vertexOnEdge value (1 or 2) for this triangle
                ! Somewhat arbitrarily, the triangle with vertices V1/V2/D1 is assigned to vertex 1.
                ! The other triangle (with vertices V2/D1/D2) is assigned to vertex 2.
                ! TODO: Think about handling this differently by constructing a departure quad and
                !        assigning V1/D1 to vertex 1, and V2/D2 to vertex 2.  The division into triangles
                !        would then be done later, after V1/D1 and V2/D2 are transformed to cell-based coordinates.
                !       This would be conceptually cleaner, and might come closer to conserving triangle area.

                triangleVertexOnEdge(triangleCount) = 1

                ! choose the source cell and set the sign of the flux
                dpInHalfPlane(1) = point_in_half_plane(edgeVertex(:,1), edgeVertex(:,2), departurePoint(:,1))
                if (dpInHalfPlane(1)) then  ! departure quad in C1
                   iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(1,iEdge)
                   fluxSign(triangleCount) = 1
                else
                   iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(2,iEdge)
                   fluxSign(triangleCount) = -1
                endif

                ! identify the vertexOnCell value (1 to nEdgesOnCell) for this triangle
                iCell = iCellTriangle(triangleCount,iEdge)
                do iVertexOnCell = 1, nEdgesOnCell(iCell)
                   if (verticesOnEdge(1,iEdge) == verticesOnCell(iVertexOnCell,iCell)) then
                      ! this is the vertex of iCell corresponding to vertex 1 of iEdge
                      triangleVertexOnCell(triangleCount) = iVertexOnCell
                   endif
                enddo

                if (on_a_sphere) then  ! set the vectors that will be used as a basis for triangle vertices

                   ! first the main edge that on which the intersection point lies
                   ! (recalling that this triangle is assigned to vertex 1)
                   iVertexOnEdge = 1
                   iOtherVertex = 2
                   edgeVector1(:,triangleCount) = edgeVertex(:,iOtherVertex) - edgeVertex(:,iVertexOnEdge)

                   ! now a side edge at this vertex
                   ! Choose E1 or E3, since iVertexOnEdge = 1

                   if (dpInHalfPlane(iVertexOnEdge)) then  ! departure point in left half-plane of iEdge; choose E1
                      iOtherEdge = iVertexOnEdge
                   else  ! departure point in right half-plane; choose E3
                      iOtherEdge = iVertexOnEdge + 2
                   endif
                   edgeNeighborVertex(1,2) = xVertexOnEdge(iOtherEdge+2,iEdge)   ! V3 for edge E1; V5 for edge E3
                   edgeNeighborVertex(2,2) = yVertexOnEdge(iOtherEdge+2,iEdge)
                   edgeVector2(:,triangleCount) = edgeNeighborVertex(:,2) - edgeVertex(:,iVertexOnEdge)

                endif  ! on a sphere

                ! add triangle with vertices V2, D1 and D2
                triangleCount = triangleCount + 1
                xTriangle(1,triangleCount,iEdge) = edgeVertex(1,2)      ! V2
                yTriangle(1,triangleCount,iEdge) = edgeVertex(2,2)
                xTriangle(2,triangleCount,iEdge) = departurePoint(1,1)  ! D1
                yTriangle(2,triangleCount,iEdge) = departurePoint(2,1)
                xTriangle(3,triangleCount,iEdge) = departurePoint(1,2)  ! D2
                yTriangle(3,triangleCount,iEdge) = departurePoint(2,2)

                ! identify the vertexOnEdge value (1 or 2) for this triangle
                ! Somewhat arbitrarily, this triangle is assigned to vertex 2;
                !  see comment above for the other central triangle
                triangleVertexOnEdge(triangleCount) = 2

                ! choose the source cell and set the sign of the flux
                ! Note: Since there is no intersection with the main edge, both D1 and D2 lie on the same side of the edge,
                !        so you might think they should always have the same source cell.
                !       However, suppose D1 lies directly on the edge and D2 is in the right half-plane.
                !       By definition, the edge itself is in the LHP. Thus triangle 1 (with zero area) is in the LHP
                !        and triangle 2 (with nonzero area) is in the RHP. So it is not safe to assume that triangle 2
                !        lies in the same half-plane (with the same source cell) as triangle 1.

                dpInHalfPlane(2) = point_in_half_plane(edgeVertex(:,1), edgeVertex(:,2), departurePoint(:,2))
                if (dpInHalfPlane(2)) then  ! departure quad in C1
                   iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(1,iEdge)
                   fluxSign(triangleCount) = 1
                else
                   iCellTriangle(triangleCount,iEdge) = cellsOnEdgeRemap(2,iEdge)
                   fluxSign(triangleCount) = -1
                endif

                ! identify the vertexOnCell value (1 to nEdgesOnCell) for this triangle
                iCell = iCellTriangle(triangleCount,iEdge)
                do iVertexOnCell = 1, nEdgesOnCell(iCell)
                   if (verticesOnEdge(2,iEdge) == verticesOnCell(iVertexOnCell,iCell)) then
                      ! this is the vertex of iCell corresponding to vertex 2 of iEdge
                      triangleVertexOnCell(triangleCount) = iVertexOnCell
                   endif
                enddo

                if (on_a_sphere) then  ! set the vectors that will be used as a basis for triangle vertices

                   ! first the main edge on which the intersection point lies
                   ! (recalling that this triangle is assigned to vertex 2)
                   iVertexOnEdge = 2
                   iOtherVertex = 1
                   edgeVector1(:,triangleCount) = edgeVertex(:,iOtherVertex) - edgeVertex(:,iVertexOnEdge)

                   ! now a side edge at this vertex
                   ! Choose E2 or E4, since iVertexOnEdge = 2

                   if (dpInHalfPlane(iVertexOnEdge)) then  ! departure point in left half-plane of iEdge; choose E2
                      iOtherEdge = iVertexOnEdge
                   else  ! departure point in right half-plane; choose E4
                      iOtherEdge = iVertexOnEdge + 2
                   endif
                   edgeNeighborVertex(1,2) = xVertexOnEdge(iOtherEdge+2,iEdge)   ! V4 for edge E2; V6 for edge E4
                   edgeNeighborVertex(2,2) = yVertexOnEdge(iOtherEdge+2,iEdge)
                   edgeVector2(:,triangleCount) = edgeNeighborVertex(:,2) - edgeVertex(:,iVertexOnEdge)

                endif  ! on a sphere

             endif  ! quadrilateral area > 0

          endif   ! edgeIntersectMain

          ! Above, we have found the vertices of each departure triangle relative to the edge.
          ! Now find the vertices of each triangle relative to the cell where the triangle is located.
          ! Also find the triangle area (positive if the flux is from cell 1 to cell 2, else negative).

          do iTri = 1, triangleCount

             iCell = iCellTriangle(iTri,iEdge)

             if (iCell >= 1 .and. iCell <= nCells) then  ! source cell exists

                if (verboseGeometry .and. etestOnProc .and. block % localBlockID == etestBlockID) then
                   if (iEdge == etest) then
                      write(stderrUnit,*) ' '
                      write(stderrUnit,*) 'vertices in edge-based coordinates, iTri =', iTri
                      do n = 1, 3
                         write(stderrUnit,*) n, xTriangle(n,iTri,iEdge), yTriangle(n,iTri,iEdge)
                      enddo
                      write(stderrUnit,*) 'edge-based triangle area =',  &
                           abs (0.5_RKIND * ( (xTriangle(2,iTri,iEdge) - xTriangle(1,iTri,iEdge)) * &
                                              (yTriangle(3,iTri,iEdge) - yTriangle(1,iTri,iEdge))  &
                                            - (yTriangle(2,iTri,iEdge) - yTriangle(1,iTri,iEdge)) * &
                                              (xTriangle(3,iTri,iEdge) - xTriangle(1,iTri,iEdge)) ) )
                      write(stderrUnit,*) 'source cell =', indexToCellID(iCell)
                   endif
                endif

                !TODO: Compute triangle area in edge-based coordinates?
                !Then would use the edge-based triangle area instead of the cell-based area.

                call shift_vertices_of_departure_triangle(&
                     iEdge,                         &
                     iCell,                         &
                     nEdgesOnCell,                  &
                     on_a_sphere,                   &
                     xVertexOnEdge, yVertexOnEdge,  &
                     xVertexOnCell, yVertexOnCell,  &
                     edgeVector1(:,iTri),           &
                     edgeVector2(:,iTri),           &
                     xTriangle(1:3,iTri,iEdge),     &
                     yTriangle(1:3,iTri,iEdge),     &
                     triangleVertexOnEdge(iTri),    &
                     triangleVertexOnCell(iTri),    &
                     triangleArea(iTri,iEdge))

                ! choose the sign of the flux
                !TODO - Eliminate triangles with very small areas?
                triangleArea(iTri,iEdge) = triangleArea(iTri,iEdge) * fluxSign(iTri)

                if (verboseGeometry .and. etestOnProc .and. block % localBlockID == etestBlockID) then
                   if (iEdge == etest) then
                      write(stderrUnit,*) 'vertices in cell-based coordinates, iTri =', iTri
                      do n = 1, 3
                         write(stderrUnit,*) n, xTriangle(n,iTri,iEdge), yTriangle(n,iTri,iEdge)
                      enddo
                      write(stderrUnit,*) 'cell-based triangle area =', abs(triangleArea(iTri,iEdge))
                   endif
                endif

             else   ! source cell does not exist

                ! do nothing, since we already have triangleArea = 0

             endif  ! source cell exists

          enddo   ! iTri

       endif      ! maskEdge = 1
    enddo         ! iEdge

  end subroutine find_departure_triangles

!-----------------------------------------------------------------------
!  routine shift_vertices_of_departure_triangle
!
!> \brief  find vertices of a departure triangle relative to the cell center
!> \author William Lipscomb
!> \date   August 2015
!> \details
!>  This routine translates the vertex coordinates of a departure triangle
!>  from edge-based to cell-based coordinates. It does this by locating
!>  triangle vertices in a basis defined by adjacent cell edges, and then
!>  transforming to a cell-centered basis.
!
!-----------------------------------------------------------------------

  subroutine shift_vertices_of_departure_triangle(&
       iEdge,                         &
       iCell,                         &
       nEdgesOnCell,                  &
       on_a_sphere,                   &
       xVertexOnEdge, yVertexOnEdge,  &
       xVertexOnCell, yVertexOnCell,  &
       edgeVector1,   edgeVector2,    &
       xTriangle,     yTriangle,      &
       triangleVertexOnEdge,          &
       triangleVertexOnCell,          &
       triangleArea)

    integer, intent(in) :: &
         iEdge,          &   !< Input: edge across which the triangle is transported
         iCell               !< Input: source cell where the triangle is located

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell        !< Input: number of edges per cell

    logical, intent(in) :: &
         on_a_sphere         !< Input: T if flow is on a sphere, F if on a plane

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnEdge,    & !< Input: x (east) coordinate of edge vertex relative to edge midpoint in local tangent plane
         yVertexOnEdge       !< Input: y (north) coordinate of edge vertex relative to edge midpoint in local tangent plane

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnCell,    & !< Input: x (east) coordinate of cell vertex relative to cell center in local tangent plane
         yVertexOnCell       !< Input: y (north) coordinate of cell vertex relative to cell center in local tangent plane

    ! Note: Together, edgeVector1 and edgeVector2 define a basis for locating triangle vertices.
    !       Both edges are among the edges included in the diagrams at the top of subroutine
    !       find_departure_triangles (7 for quads, 5 for hexes).
    !       The edges are assumed to be adjacent and to share a vertex; each vector points away from the common vertex.
    !       Both edges must belong to the source cell, with the second edge following the first edge
    !       when proceeding CCW around the cell. This implies that the cross product (edgeVector1 x edgeVector2) < 0.

    real(kind=RKIND), dimension(2), intent(inout) ::  &
         edgeVector1,      & !< Input: x and y coordinates of a vector pointing along an edge
         edgeVector2         !< Input: x and y coordinates for a vector along an adjacent edge
                             ! These are intent(inout) in case they need to be switched with each other below

    real(kind=RKIND), dimension(3), intent(inout) ::  &
         xTriangle, yTriangle  !< Input/output: x and y coordinates of vertices of departure triangle
                               !                On input: relative to the one of the vertices on the edge
                               !                On output: relative to the center of the cell where the triangle lies
                               ! Note: It is assumed that the first of the 3 points is located at an edge vertex.
                               !       The other two points will be departure points or intersection points
                               !       (or possibly another edge vertex).

    integer, intent(in) ::  &
         triangleVertexOnEdge, &!< Input: vertexOnEdge index for a given triangle vertex, relative to iEdge
         triangleVertexOnCell   !< Input: vertexOnCell index for a given triangle vertex, relative to iCell

    real(kind=RKIND), intent(out) ::  &
         triangleArea           !< Output: area of the triangle

    ! local variables

    integer :: iTrivert

    integer :: iVertexOnCell, iVertexOnCell_p1, iVertexOnCell_m1

    real(kind=RKIND) :: coeff_a, coeff_b, denom, crossProduct

    real(kind=RKIND), dimension(2) :: &
         edgeVectorTemp,   &  ! temporary vector
         edgeVector1OnCell, edgeVector2OnCell  ! input edgeVectors in cell-based coordinates

    character(len=strKIND) :: &
         errorMessage ! error message for abort

    ! This subroutine has two methods for transforming from edge-based to cell-based coordinates:
    !  a simple method for planar meshes and a more complex method for spherical meshes.
    !
    ! A picture may be useful in understanding the more complex method, which uses basis vectors.
    ! Consider a side triangle that shares one vertex with iEdge and lies in iCell.  Its vertices are V1, D1 and I.
    ! (This is for a quad mesh, but the same method applies to a hex mesh.)
    !
    !                V3            V4
    !                |             |
    !         D1-----|I            |
    !            \   |             |
    !             \  |E1           |E2
    !              \ |             |
    !    V7________ \V1____________V2__________V8
    !          E5    |             |     E6
    !                |             |
    !                |E3           |E4
    !                |             |
    !                |             |
    !                V5            V6
    !
    ! The 2D basis vectors are u1 = V3 - V1  (the vector pointing along E1)
    !                          u2 = V7 - V1  (the vector pointing along E5)
    !
    ! The steps for each triangle vertex (V1, D1 or I) are as follows:
    ! (1) Given the (x,y) coordinates relative to the midpoint of the main edge,
    !     subtract x/yVertexOnEdge to get coordinates relative to V1.
    !     Thus V1 will have coordinates (0,0).
    ! (2) Next, express the triangle vertex coordinates relative to V1 as a
    !     linear combination of the basis vectors u1 and u2:
    !           (x,y) = a*u1 + b*u2
    !     The formulas are as follows:
    !           a = (x*u2(2) - y*u2(1)) / (u1(1)*u2(2) - u2(1)*u1(2))
    !           b = (y*u1(1) - x*u1(2)) / (u1(1)*u2(2) - u2(1)*u1(2))
    !     Note: The denominator is just the magnitude of the cross-product u1 x u2.
    !     Since the two vectors are never parallel, it must be nonzero.
    ! (3) Given the coefficients a and b, transform to the cell-centered basis:
    !           (x_c, y_c) = a*u1_c + b*u2_c
    !     where u1_c and u2_c are the same edges as u1 and u2, but given in
    !     cell-based coordinates.
    ! One property of this method is that any point lying on an edge in edge-based
    ! coordinates (such that either a = 0 or b = 0) will also lie on that edge
    ! in cell-based coordinates.
    !
    ! You might ask: Why go to all this trouble?  For instance, why not just subtract
    !  x/yVertexOnEdge and add x/yVertexOnCell, as for a plane?
    !
    ! The answer is that on a sphere, a point lying on an edge in edge-based coordinates
    !  will not necessarily lie on that edge in cell-based coordinates. It could lie
    !  slightly outside the source cell and therefore outside the range of values
    !  determined by the gradient limiter.  This can result in small negative masses
    !  under transport.
    !
    ! For future reference: A drawback of this method--at least aesthetically--is that
    ! when we have two central triangles that are part of a departure quadrilateral,
    ! one triangle is (somewhat arbitrarily) referenced to V1 and the other to V2.
    ! This means that when the two triangles are separately transformed to the source cell,
    ! they can overlap.
    ! I am not sure whether this is a bad thing.  I don't think it will destroy the nice properties of IR
    ! (e.g., tracer monotonicity).  If we decide later that it is a bad thing,
    ! an alternative would be to write a separate procedure to deal with departure quadrilaterals.

    if (on_a_sphere) then   ! use the more complex method with an expansion in terms of basis vectors

       ! This method requires edgeVector1 x edgeVector2 < 0.
       ! If this is not already the case, then interchange the two vectors.

       ! Also check that the cross product is nonzero; else the vectors are parallel and the algorithm will not work.
       ! This can happen on a sphere for cells near the poles.
       ! TODO - Mask out cells and edges near the poles on the spherical hex mesh?

       call cross_product_2d(edgeVector1, edgeVector2, &
            crossProduct)

       if (abs(crossProduct) < eps11) then  ! the two vectors are (nearly) parallel; something is wrong
          write(errorMessage,*) 'basis vectors for IR coordinate transformations must not be parallel: iEdge, iCell =', &
               iEdge, iCell
          write(stderrUnit,*) 'ERROR: '//trim(errorMessage)
          call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
       endif

       if (crossProduct > 0.0_RKIND) then
          edgeVectorTemp(:) = edgeVector1(:)
          edgeVector1(:)    = edgeVector2(:)
          edgeVector2(:)    = edgeVectorTemp(:)
       endif

       ! Transform the two edge vectors to cell-based coordinates.
       ! Note: It is assumed that edgeVector2 follows edgeVector1 when proceeding CCW around the cell.
       !       Thus edgeVector1 is the difference between iVertexOnCell - 1 and iVertexOnCell,
       !        whereas edgeVector2 is the difference between iVertexOnCell + 1 and iVertexOnCell.
       iVertexOnCell = triangleVertexOnCell
       iVertexOnCell_m1 = iVertexOnCell - 1
       if (iVertexOnCell_m1 < 1) &
            iVertexOnCell_m1 = iVertexOnCell_m1 + nEdgesOnCell(iCell)
       iVertexOnCell_p1 = iVertexOnCell + 1
       if (iVertexOnCell_p1 > nEdgesOnCell(iCell)) &
            iVertexOnCell_p1 = iVertexOnCell_p1 - nEdgesOnCell(iCell)

       edgeVector1OnCell(1) = xVertexOnCell(iVertexOnCell_m1,iCell) - xVertexOnCell(iVertexOnCell,iCell)
       edgeVector1OnCell(2) = yVertexOnCell(iVertexOnCell_m1,iCell) - yVertexOnCell(iVertexOnCell,iCell)

       edgeVector2OnCell(1) = xVertexOnCell(iVertexOnCell_p1,iCell) - xVertexOnCell(iVertexOnCell,iCell)
       edgeVector2OnCell(2) = yVertexOnCell(iVertexOnCell_p1,iCell) - yVertexOnCell(iVertexOnCell,iCell)

       denom = edgeVector1(1)*edgeVector2(2) - edgeVector2(1)*edgeVector1(2)  ! z component of cross product

       ! loop over the 3 triangle vertices
       do iTriVert = 1, 3

          ! Given the x/y coordinates relative to iEdge, compute x/y coordinates
          ! relative to the first triangle vertex (which itself is a vertex of iEdge)

          xTriangle(iTriVert) = xTriangle(iTriVert)  &
                              - xVertexOnEdge(triangleVertexOnEdge,iEdge)

          yTriangle(iTriVert) = yTriangle(iTriVert)  &
                              - yVertexOnEdge(triangleVertexOnEdge,iEdge)

          ! Compute the coefficients a and b in the basis defined by the input edge vectors

          coeff_a = (xTriangle(iTriVert) * edgeVector2(2) - yTriangle(iTriVert) * edgeVector2(1)) / denom
          coeff_b = (yTriangle(iTriVert) * edgeVector1(1) - xTriangle(iTriVert) * edgeVector1(2)) / denom

          ! Express the triangle vertex in cell-based coordinates by taking a linear combination
          ! of the two edge vectors

          xTriangle(iTriVert) = xVertexOnCell(triangleVertexOnCell,iCell)  &
                              + coeff_a * edgeVector1OnCell(1) + coeff_b * edgeVector2OnCell(1)
          yTriangle(iTriVert) = yVertexOnCell(triangleVertexOnCell,iCell)  &
                              + coeff_a * edgeVector1OnCell(2) + coeff_b * edgeVector2OnCell(2)

       enddo  ! iTriVert

    else   ! on a plane; do a simpler calculation that gives the same result

       ! loop over triangle vertices
       do iTriVert = 1, 3

          xTriangle(iTriVert) = xTriangle(iTriVert)  &
                              - xVertexOnEdge(triangleVertexOnEdge,iEdge) &
                              + xVertexOnCell(triangleVertexOnCell,iCell)

          yTriangle(iTriVert) = yTriangle(iTriVert)  &
                              - yVertexOnEdge(triangleVertexOnEdge,iEdge) &
                              + yVertexOnCell(triangleVertexOnCell,iCell)

       enddo   ! iTriVert

    endif  ! on a sphere

    ! compute the area of the triangle

    triangleArea = abs (0.5_RKIND * ( (xTriangle(2) - xTriangle(1)) * (yTriangle(3) - yTriangle(1))  &
                                    - (yTriangle(2) - yTriangle(1)) * (xTriangle(3) - xTriangle(1)) ) )

  end subroutine shift_vertices_of_departure_triangle

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine get_triangle_quadrature_points
!
!> \brief  quadrature points for integrating over triangles
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine computes quadrature points for integrating over triangles;
!>  up to 6 points as required to integrate 4th degree polynomials exactly.
!
!-----------------------------------------------------------------------
!
! For each triangle, find the coordinates of the quadrature points needed
!  to compute integrals of polynomials up to 4th degree.
! Two options here: Either compute 3 quadrature points as required for
!  linear/quadratic polynomials, or 6 points as required for cubic/quartic.
!
! I2 = A * (f(p1) + f(p2) + f(p3))
! where A is the triangle area, and the three quadrature points are located
! halfway between the midpoint and the three vertices of the triangle.
!
! I4 = A * { w1 * [f(p1) + f(p2) + f(p3)]
!          + w2 * [f(p4) + f(p5) + f(p6)]}
!
! where w1 = 0.109951743655321885
!       w2 = 0.223381589678011389
!
!       p1 = q1*(x1,y1) +  q1*(x2,y2) + q2*(x3,y3)
!       p2 = q1*(x1,y1) +  q2*(x2,y2) + q1*(x3,y3)
!       p3 = q2*(x1,y1) +  q1*(x2,y2) + q1*(x3,y3)
!       p4 = q3*(x1,y1) +  q4*(x2,y2) + q4*(x3,y3)
!       p5 = q4*(x1,y1) +  q3*(x2,y2) + q4*(x3,y3)
!       p6 = q4*(x1,y1) +  q4*(x2,y2) + q3*(x3,y3)
!
!       (x1,y1), (x2,y2) and (x3,y3) are the triangle vertices
!
!       q1 = 0.0915762135097704655
!       q2 = 0.816847572980458514
!       q3 = 0.108103018168070275
!       q4 = 0.445948490915964113
!
!-----------------------------------------------------------------------

  subroutine get_triangle_quadrature_points (&
       nEdges,            &
       nTriPerEdgeRemap,  &
       maskEdge,          &
       xTriangle, yTriangle)

    integer, intent(in) :: &
         nEdges,           & !< Input: number of edges
         nTriPerEdgeRemap    !< Input: number of triangles per edge

    integer, dimension(:), intent(in) ::  &
         maskEdge           !< Input: = 1 if IR fluxes need to be computed across an edge, else = 0

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &  ! nQuadPoints, nTriPerEdgeRemap, nEdges
         xTriangle, yTriangle   !< Input/output: on input, coordinates of triangle vertices (indices 1-3)
                                !<               on output, coordinates of quadrature points (indices 1-3 or 1-6)

    ! local variables

    integer ::   &
         iEdge, iTri            ! counting indices

    real(kind=RKIND) ::  &
         xMidpoint, yMidpoint   ! coordinates of triangle midpoint

    real(kind=RKIND) ::  &
         x1Vertex, y1Vertex,  &     ! coordinates of triangle vertices
         x2Vertex, y2Vertex,  &
         x3Vertex, y3Vertex

    if (nQuadPoints == 3) then

       ! compute 3 quadrature points that are accurate for polynomials up to degree 2
       ! Not used in CICE (which uses 6 quadrature points), but likely to be used for land ice

       do iEdge = 1, nEdges
          if (maskEdge(iEdge) == 1) then

             do iTri = 1, nTriPerEdgeRemap

                ! coordinates of triangle midpoint
                xMidpoint = (xTriangle(1,iTri,iEdge) + xTriangle(2,iTri,iEdge) + xTriangle(3,iTri,iEdge)) / 3.0_RKIND
                yMidpoint = (xTriangle(1,iTri,iEdge) + xTriangle(2,iTri,iEdge) + xTriangle(3,iTri,iEdge)) / 3.0_RKIND

                ! coordinates of 3 quadrature points
                xTriangle(1,iTri,iEdge) = 0.5_RKIND * (xTriangle(1,iTri,iEdge) + xMidpoint)
                yTriangle(1,iTri,iEdge) = 0.5_RKIND * (yTriangle(1,iTri,iEdge) + yMidpoint)

                xTriangle(2,iTri,iEdge) = 0.5_RKIND * (xTriangle(2,iTri,iEdge) + xMidpoint)
                yTriangle(2,iTri,iEdge) = 0.5_RKIND * (yTriangle(2,iTri,iEdge) + yMidpoint)

                xTriangle(3,iTri,iEdge) = 0.5_RKIND * (xTriangle(3,iTri,iEdge) + xMidpoint)
                yTriangle(3,iTri,iEdge) = 0.5_RKIND * (yTriangle(3,iTri,iEdge) + yMidpoint)

             enddo     ! nTriPerEdgeRemap

          endif        ! maskEdge = 1
       enddo           ! iEdge

    elseif (nQuadPoints == 6) then

       ! compute coordinates of 6 quadrature points that are accurate for polynomials up to degree 4
       do iEdge = 1, nEdges
          if (maskEdge(iEdge) == 1) then

             do iTri = 1, nTriPerEdgeRemap

                ! save coordinates of the 3 vertices
                x1Vertex = xTriangle(1,iTri,iEdge)
                y1Vertex = yTriangle(1,iTri,iEdge)
                x2Vertex = xTriangle(2,iTri,iEdge)
                y2Vertex = yTriangle(2,iTri,iEdge)
                x3Vertex = xTriangle(3,iTri,iEdge)
                y3Vertex = yTriangle(3,iTri,iEdge)

                ! coordinates of 6 quadrature points

                xTriangle(1,iTri,iEdge) = q1TriangleQP * x1Vertex + q1TriangleQP * x2Vertex + q2TriangleQP * x3Vertex
                yTriangle(1,iTri,iEdge) = q1TriangleQP * y1Vertex + q1TriangleQP * y2Vertex + q2TriangleQP * y3Vertex

                xTriangle(2,iTri,iEdge) = q1TriangleQP * x1Vertex + q2TriangleQP * x2Vertex + q1TriangleQP * x3Vertex
                yTriangle(2,iTri,iEdge) = q1TriangleQP * y1Vertex + q2TriangleQP * y2Vertex + q1TriangleQP * y3Vertex

                xTriangle(3,iTri,iEdge) = q2TriangleQP * x1Vertex + q1TriangleQP * x2Vertex + q1TriangleQP * x3Vertex
                yTriangle(3,iTri,iEdge) = q2TriangleQP * y1Vertex + q1TriangleQP * y2Vertex + q1TriangleQP * y3Vertex

                xTriangle(4,iTri,iEdge) = q3TriangleQP * x1Vertex + q4TriangleQP * x2Vertex + q4TriangleQP * x3Vertex
                yTriangle(4,iTri,iEdge) = q3TriangleQP * y1Vertex + q4TriangleQP * y2Vertex + q4TriangleQP * y3Vertex

                xTriangle(5,iTri,iEdge) = q4TriangleQP * x1Vertex + q3TriangleQP * x2Vertex + q4TriangleQP * x3Vertex
                yTriangle(5,iTri,iEdge) = q4TriangleQP * y1Vertex + q3TriangleQP * y2Vertex + q4TriangleQP * y3Vertex

                xTriangle(6,iTri,iEdge) = q4TriangleQP * x1Vertex + q4TriangleQP * x2Vertex + q3TriangleQP * x3Vertex
                yTriangle(6,iTri,iEdge) = q4TriangleQP * y1Vertex + q4TriangleQP * y2Vertex + q3TriangleQP * y3Vertex

             enddo   ! nTriPerEdgeRemap

          endif      ! maskEdge = 1
       enddo         ! nEdges

    else   ! nQuadPoint /= 3 or 6

       write(stderrUnit,*) 'ERROR: must have nQuadPoints = 3 or 6'
       call mpas_dmpar_global_abort('MPAS-seaice: IR must have nQuadPoints = 3 or 6')

    endif            ! nQuadPoints

  end subroutine get_triangle_quadrature_points

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine integrate_fluxes_over_triangles
!
!> \brief  integration of mass and tracer fluxes over departure triangles
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine sums over quadrature points in each departure triangle to
!>  obtain the area-integrated fluxes of masses and tracers across cell edges.
!
!-----------------------------------------------------------------------

  subroutine integrate_fluxes_over_triangles(&
       tracersHead,          &
       nCells,               &
       nEdges,               &
       maskEdge,             &
       xTriangle, yTriangle, &
       triangleArea,         &
       cellsOnEdge,          &
       iCellTriangle,        &
       indexToCellID,        &
       indexToEdgeID,        &
       block)

    type(tracer_type), pointer :: &
         tracersHead      !< Input/output: pointer to first element of linked list of tracers
                          !                This subroutine outputs thisTracer % edgeFlux2d/3d

    integer, intent(in) :: &
         nCells, &        !< Input: number of cells
         nEdges           !< Input: number of edges

    integer, dimension(:), intent(in) ::  &
         maskEdge         !< Input: = 1 if IR fluxes need to be computed across an edge, else = 0

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &  ! nQuadPoints, nTriPerEdgeRemap, nEdges
         xTriangle, yTriangle     !< Input: coordinates of triangle quadrature points (indices 1-3 or 1-6)

    real(kind=RKIND), dimension(:,:), intent(in) :: &  ! nTriPerEdgeRemap, nEdges
         triangleArea             !< Input: triangle area
                                  !   Note: Area is positive if flux across the edge is from cell 1 to cell 2, else is negative

    integer, dimension(:,:), intent(in) :: &  ! 2, nEdges
         cellsOnEdge              !< Input: index for each cell neighboring an edge

    integer, dimension(:,:), intent(in) :: &  ! nTriPerEdgeRemap, nEdges
         iCellTriangle            !< Input: index of cell where each departure triangle is located

    integer, dimension(:), intent(in) ::  &
         indexToCellID,  & !< Input: global index for each cell on a block (diagnostic only)
         indexToEdgeID     !< Input: global index for each edge on a block (diagnostic only)

    type(block_type), intent(in) :: &
         block             !< Input: local block (diagnostic only)

    ! local variables

    type(tracer_type), pointer :: &
         thisTracer,   &   ! pointer that loops through linked list of tracers
         parentTracer, &   ! parent of thisTracer
         dummyTracer       ! dummy tracer with value of 1 everywhere

    integer ::   &
         iEdge, iCell, iCat, iLayer, iTri, iqp

    integer :: &
         nTriPerEdgeRemap,          & ! number of triangles per edge
         nCategories,               & ! number of ice thickness categories
         nLayers                      ! number of layers

    real(kind=RKIND), dimension(:), allocatable ::   &
         tracerIntegral2D         ! integral over a triangle of mass, mass*tracer, etc.
                                  ! for each category

    real(kind=RKIND), dimension(:,:), allocatable ::   &
         tracerIntegral3D         ! integral over a triangle of mass, mass*tracer, etc.
                                  ! for each category and layer

    nTriPerEdgeRemap = size(triangleArea,1)

    ! loop through linked list of tracers

    thisTracer => tracersHead
    do while (associated(thisTracer))

       if (verboseFluxes) then
          write(stderrUnit,*) ' '
          write(stderrUnit,*) 'Tracer: ', trim(thisTracer % tracerName)
       endif

       ! check that nParents <= 3
       if (thisTracer % nParents > 3) then
          write(stderrUnit,*) 'ERROR: Tracer has > 3 parents:', trim(thisTracer % tracerName)
          call mpas_dmpar_global_abort('MPAS-seaice: IR: Tracer has > 3 parents: '//trim(thisTracer % tracerName))
       endif

       ! initialize edge flux and allocate tracerIntegral array
       if (thisTracer % ndims == 2) then

          thisTracer % edgeFlux2D(:,:) = 0.0_RKIND
          nCategories = size(thisTracer % array2D, 1)
          allocate(tracerIntegral2D(nCategories))

       elseif (thisTracer % ndims == 3) then

          thisTracer % edgeFlux3D(:,:,:) = 0.0_RKIND
          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)
          allocate(tracerIntegral3D(nLayers,nCategories))

       endif

       ! set parent tracer
       if (thisTracer % nParents == 0) then   ! mass-like field (fractional ice concentration for CICE)

          ! This field has no parent, but set a parentTracer pointer to have value of 1 at each QP.
          ! This ensures correct values for triangleValue2D/3D below.

          allocate(dummyTracer)
          if (thisTracer % ndims == 2) then

             allocate(dummyTracer % triangleValue2D(nCategories,nQuadPoints,nTriPerEdgeRemap,nEdges))
             dummyTracer % triangleValue2D(:,:,:,:) = 1.0_RKIND
             dummyTracer % ndims = 2

          elseif (thisTracer % ndims == 3) then

             allocate(dummyTracer % triangleValue3D(nLayers,nCategories,nQuadPoints,nTriPerEdgeRemap,nEdges))
             dummyTracer % triangleValue3D(:,:,:,:,:) = 1.0_RKIND
             dummyTracer % ndims = 3

          endif
          parentTracer => dummyTracer

       else  ! nParents = 1, 2 or 3

          parentTracer => thisTracer % parent

       endif

       if (verboseFluxes) then
          if (thisTracer % nParents > 0) then
             write(stderrUnit,*) 'Parent: ', trim(parentTracer % tracerName)
          endif
       endif

       ! Integrate the fluxes of this tracer (for each category and layer) over each triangle of each edge.
       ! The integral formulas are described above in get_triangle_quadrature_points.

       if (thisTracer % ndims == 2) then

          do iEdge = 1, nEdges
             if (maskEdge(iEdge) == 1) then

                do iTri = 1, nTriPerEdgeRemap

                   if (triangleArea(iTri,iEdge) /= 0.0_RKIND) then

                      ! In the following arrays:
                      ! 1st index of triangleValue2D = category; 2nd index = QP
                      ! 1st index of center/xGrad/yGrad = category
                      ! 1st index of xTriangle/yTriangle = QP
                      ! 1st index of tracerIntegral2D = category
                      ! 1st index of edgeFlux2D = category

                      ! identify the cell where the triangle is located
                      iCell = iCellTriangle(iTri,iEdge)

                      ! evaluate the tracer at each quadrature point
                      do iqp = 1, nQuadPoints
                         thisTracer % triangleValue2D(:,iqp,iTri,iEdge) = &
                              thisTracer % center2D(:,iCell)                             &
                            + thisTracer % xGrad2D(:,iCell) * xTriangle(iqp,iTri,iEdge)  &
                            + thisTracer % yGrad2D(:,iCell) * yTriangle(iqp,iTri,iEdge)
                      enddo

                      ! evaluate the product mass*tracer at each quadrature point (using parent tracer info computed already)
                      ! In nParents = 0, this is mass
                      ! If nParents = 1, this is mass*tracer1
                      ! If nParents = 2, this is mass*tracer1*tracer2
                      ! If nParents = 3, this is mass*tracer1*tracer2*tracer3
                      ! Note: Parent tracer must have ndims = 2
                      thisTracer % triangleValue2D(:,:,iTri,iEdge) = parentTracer % triangleValue2D(:,:,iTri,iEdge)  &
                                                                     * thisTracer % triangleValue2D(:,:,iTri,iEdge)

                      ! integrate over the triangle by summing over quadrature points
                      tracerIntegral2D(:) = 0.0_RKIND
                      do iqp = 1, nQuadPoints
                         tracerIntegral2D(:) = tracerIntegral2D(:) &
                                             + weightQuadPoint(iqp) * thisTracer % triangleValue2D(:,iqp,iTri,iEdge)
                      enddo

                      ! increment the area-weighted flux across the edge
                      thisTracer % edgeFlux2D(:,iEdge) = thisTracer % edgeFlux2D(:,iEdge) &
                                                       + triangleArea(iTri,iEdge) * tracerIntegral2D(:)

                   endif ! triangleArea /= 0

                enddo    ! nTriPerEdgeRemap

             endif       ! maskEdge = 1
          enddo          ! nEdges

          ! Check for negative reconstructed ice area
          if (trim(thisTracer % tracerName) == 'iceAreaCategory') then
             do iEdge = 1, nEdges
                do iTri = 1, nTriPerEdgeRemap
                   if (triangleArea(iTri,iEdge) /= 0.0_RKIND) then
                      iCell = iCellTriangle(iTri,iEdge)
                      do iqp = 1, nQuadPoints
                         do iCat = 1, nCategories
                            if (thisTracer % triangleValue2D(iCat,iqp,iTri,iEdge) < 0.0_RKIND) then
                               write(stderrUnit,*)'ERROR: Negative reconstructed ice area'
                               write(stderrUnit,*)'nCells =', nCells
                               write(stderrUnit,*)'iCat, iCell, global iCell:', iCat, iCell, indexToCellID(iCell)
                               write(stderrUnit,*)'iEdge, global iEdge, iTri, iqp:', iEdge, indexToEdgeID(iEdge), iTri, iqp
                               write(stderrUnit,*)'triangle area:', triangleArea(iTri,iEdge)
                               write(stderrUnit,*)'tracer val:', thisTracer % triangleValue2D(iCat,iqp,iTri,iEdge)
                               write(stderrUnit,*)'center val:', thisTracer % center2D(iCat,iCell)
                               write(stderrUnit,*)'x gradient:', thisTracer % xGrad2D(iCat,iCell)
                               write(stderrUnit,*)'y gradient:', thisTracer % yGrad2D(iCat,iCell)
                               call mpas_dmpar_global_abort('MPAS-seaice: IR negative reconstructed ice area (nDims == 2)')
                            endif  ! negative area
                         enddo   ! iCat
                      enddo   ! iqp
                   endif   ! triangleArea > 0
                enddo   ! iTri
             enddo   ! iEdge
          endif   ! iceAreaCategory

          deallocate(tracerIntegral2D)

       elseif (thisTracer % ndims == 3) then

          do iEdge = 1, nEdges
             if (maskEdge(iEdge) == 1) then

                do iTri = 1, nTriPerEdgeRemap

                   if (triangleArea(iTri,iEdge) /= 0.0_RKIND) then

                      ! In the following arrays:
                      ! 1st index of triangleValue3D = layer; 2nd index = category; 3rd index = QP
                      ! 1st index of center/xGrad/yGrad = layer; 2nd index = category
                      ! 1st index of xTriangle/yTriangle = QP
                      ! 1st index of tracerIntegral = layer; 2nd index = category
                      ! 1st index of edgeFlux = layer; 2nd index = category

                      ! identify the cell where the triangle is located
                      iCell = iCellTriangle(iTri,iEdge)

                      ! evaluate the tracer at each quadrature point
                      do iqp = 1, nQuadPoints
                         thisTracer % triangleValue3D(:,:,iqp,iTri,iEdge) = &
                              thisTracer % center3D(:,:,iCell)   &
                            + thisTracer % xGrad3D(:,:,iCell) * xTriangle(iqp,iTri,iEdge)  &
                            + thisTracer % yGrad3D(:,:,iCell) * yTriangle(iqp,iTri,iEdge)
                      enddo

                      ! evaluate the product mass*tracer at each quadrature point (using parent tracer info computed already)
                      ! In nParents = 0, this is mass
                      ! If nParents = 1, this is mass*tracer1
                      ! If nParents = 2, this is mass*tracer1*tracer2
                      ! If nParents = 3, this is mass*tracer1*tracer2*tracer3
                      ! Note: Parent tracer can have ndims = 2 or 3
                      if (parentTracer % ndims == 2) then
                         do iLayer = 1, nLayers
                            thisTracer % triangleValue3D(iLayer,:,:,iTri,iEdge) = &
                                 parentTracer % triangleValue2D(:,:,iTri,iEdge) * &
                                 thisTracer % triangleValue3D(iLayer,:,:,iTri,iEdge)
                         enddo
                      else  ! parents has ndims = 3
                         thisTracer % triangleValue3D(:,:,:,iTri,iEdge) = parentTracer % triangleValue3D(:,:,:,iTri,iEdge) &
                                                                          * thisTracer % triangleValue3D(:,:,:,iTri,iEdge)
                      endif

                      ! integrate over the triangle by summing over quadrature points
                      tracerIntegral3D(:,:) = 0.0_RKIND
                      do iqp = 1, nQuadPoints
                         tracerIntegral3D(:,:) = tracerIntegral3D(:,:) &
                                               + weightQuadPoint(iqp) * thisTracer % triangleValue3D(:,:,iqp,iTri,iEdge)
                      enddo

                      ! increment the area-weighted flux across the edge
                      thisTracer % edgeFlux3D(:,:,iEdge) = thisTracer % edgeFlux3D(:,:,iEdge) &
                                                         + triangleArea(iTri,iEdge) * tracerIntegral3D(:,:)

                   endif ! triangleArea /= 0

                enddo    ! nTriPerEdgeRemap

             endif       ! maskEdge = 1
          enddo          ! nEdges

          ! Check for negative reconstructed ice area
          if (trim(thisTracer % tracerName) == 'iceAreaCategory') then
             do iEdge = 1, nEdges
                do iTri = 1, nTriPerEdgeRemap
                   if (triangleArea(iTri,iEdge) /= 0.0_RKIND) then
                      iCell = iCellTriangle(iTri,iEdge)
                      do iqp = 1, nQuadPoints
                         do iCat = 1, nCategories
                            do iLayer = 1, nLayers
                               if (thisTracer % triangleValue3D(iLayer,iCat,iqp,iTri,iEdge) < 0.0_RKIND) then
                                  write(stderrUnit,*)'ERROR: Negative reconstructed ice area'
                                  write(stderrUnit,*)'nCells =', nCells
                                  write(stderrUnit,*)'iLayer, iCat, iCell, global iCell:', &
                                       iLayer, iCat, iCell, indexToCellID(iCell)
                                  write(stderrUnit,*)'iEdge, global iEdge, iTri, iqp:', iEdge, indexToEdgeID(iEdge), iTri, iqp
                                  write(stderrUnit,*)'triangle area:', triangleArea(iTri,iEdge)
                                  write(stderrUnit,*)'tracer val:', thisTracer % triangleValue3D(iLayer,iCat,iqp,iTri,iEdge)
                                  write(stderrUnit,*)'center val:', thisTracer % center3D(iLayer,iCat,iCell)
                                  write(stderrUnit,*)'x gradient:', thisTracer % xGrad3D(iLayer,iCat,iCell)
                                  write(stderrUnit,*)'y gradient:', thisTracer % yGrad3D(iLayer,iCat,iCell)
                                  call mpas_dmpar_global_abort('MPAS-seaice: IR negative reconstructed ice area (nDims == 3)')
                               endif  ! negative area
                            enddo   ! iLayer
                         enddo   ! iCat
                      enddo   ! iqp
                   endif   ! triangleArea > 0
                enddo   ! iTri
             enddo   ! iEdge
          endif   ! iceAreaCategory

          deallocate(tracerIntegral3D)

       endif   ! ndims

       ! clean up

       if (associated(dummyTracer)) then
          if (dummyTracer % ndims == 2) then
             if (associated(dummyTracer % triangleValue2D)) deallocate(dummyTracer % triangleValue2D)
          elseif (dummyTracer % ndims == 3) then
             if (associated(dummyTracer % triangleValue3D)) deallocate(dummyTracer % triangleValue3D)
          endif
          deallocate(dummyTracer)
       endif

       thisTracer => thisTracer % next

    enddo   ! associated(thisTracer)

  end subroutine integrate_fluxes_over_triangles

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine compute_mass_tracer_products
!
!> \brief  compute mass*tracer products in each cell
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  This routine computes values of mass, mass*tracer, mass*tracer*tracer, etc.
!>  in each grid cell.
!
!-----------------------------------------------------------------------

  subroutine compute_mass_tracer_products(&
       nCells,        &
       tracersHead)

    integer, intent(in) :: &
         nCells            !< Input: number of cells to be updated

    type(tracer_type), pointer :: &
         tracersHead       !< Input/output: pointer to first element of linked list of tracers
                           ! old values on input; new values on output

    ! local variables

    type(tracer_type), pointer :: &
         thisTracer,   &   ! pointer that loops through linked list of tracers
         parentTracer, &   ! pointer to parent of a given tracer
         dummyTracer       ! dummy tracer with value of 1 everywhere

    integer :: &
         nCategories, &    ! number of categories
         nLayers           ! number of layers

    integer :: iCell, iCat, iLayer

    ! loop through linked list of tracers

    thisTracer => tracersHead
    do while (associated(thisTracer))

       ! get array sizes

       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)

       elseif (thisTracer % ndims == 3) then

          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)

       endif

       ! set parent tracer
       if (thisTracer % nParents == 0) then   ! mass-like field (fractional ice concentration for CICE)

          ! Set a parentTracer pointer to have a mass-tracer product of 1 in each cell.
          ! This ensures correct values for massTracerProduct below.
          allocate(dummyTracer)
          if (thisTracer % ndims == 2) then

             allocate(dummyTracer % massTracerProduct2D(nCategories,nCells))
             dummyTracer % massTracerProduct2D(:,:) = 1.0_RKIND
             dummyTracer % ndims = 2

          elseif (thisTracer % ndims == 3) then

             allocate(dummyTracer % massTracerProduct3D(nLayers,nCategories,nCells))
             dummyTracer % massTracerProduct3D(:,:,:) = 1.0_RKIND
             dummyTracer % ndims = 3

          endif  ! ndims

          parentTracer => dummyTracer

       else  ! nParents = 1, 2 or 3

          parentTracer => thisTracer % parent

       endif ! nParents

       ! evaluate the product mass*tracer
       ! If nParents = 0, this is mass
       ! If nParents = 1, this is mass*tracer1
       ! If nParents = 2, this is mass*tracer1*tracer2
       ! If nParents = 3, this is mass*tracer1*tracer2*tracer3

       if (thisTracer % ndims == 2) then

          do iCell = 1, nCells
             do iCat = 1, nCategories
                thisTracer % massTracerProduct2D(iCat,iCell) = parentTracer % massTracerProduct2D(iCat,iCell)  &
                                                             * thisTracer % array2D(iCat,iCell)
             enddo
          enddo

       elseif (thisTracer % ndims == 3) then   ! Note: Parent tracer can have ndims = 2 or 3

          if (parentTracer % ndims == 2) then

             do iCell = 1, nCells
                do iCat = 1, nCategories
                   do iLayer = 1, nLayers
                      thisTracer % massTracerProduct3D(iLayer,iCat,iCell) = parentTracer % massTracerProduct2D(iCat,iCell) &
                                                                          * thisTracer % array3D(iLayer,iCat,iCell)
                   enddo
                enddo
             enddo

          elseif (parentTracer % ndims == 3) then

             do iCell = 1, nCells
                do iCat = 1, nCategories
                   do iLayer = 1, nLayers
                      thisTracer % massTracerProduct3D(iLayer,iCat,iCell) = parentTracer % massTracerProduct3D(iLayer,iCat,iCell) &
                                                                          * thisTracer % array3D(iLayer,iCat,iCell)
                   enddo
                enddo
             enddo

          endif

       endif  ! ndims

       ! clean up
       if (associated(dummyTracer)) then
          if (dummyTracer % ndims == 2) then
             if (associated(dummyTracer % massTracerProduct2D)) deallocate(dummyTracer % massTracerProduct2D)
          elseif (dummyTracer % ndims == 3) then
             if (associated(dummyTracer % massTracerProduct3D)) deallocate(dummyTracer % massTracerProduct3D)
          endif
          deallocate(dummyTracer)
       endif

       thisTracer => thisTracer % next

    enddo  ! while(associated)

  end subroutine compute_mass_tracer_products

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine update_mass_and_tracers
!
!> \brief  new values of mass and tracers in each cell
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine updates mass and tracer values in each grid cell by
!>  summing over the fluxes across each cell edge.
!>  NOTE: This routine should be called after compute_mass_tracer_products.
!
!-----------------------------------------------------------------------

  subroutine update_mass_and_tracers(&
       nCellsSolve,      &
       nEdgesOnCell,     &
       edgesOnCell,      &
       cellsOnEdge,      &
       areaCell,         &
       tracersHead,      &
       indexToCellID,    &
       indexToEdgeID,    &
       block)

    integer, intent(in) :: &
         nCellsSolve       !< Input: number of locally owned cells to be updated

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell      !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         edgesOnCell,  &   !< Input: edge index for each edge of a given cell
         cellsOnEdge       !< Input: cell index for each of two cells sharing a given edge

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaCell          !< Input: cell area

    type(tracer_type), pointer :: &
         tracersHead       !< Input/output: pointer to first element of linked list of tracers
                           ! old tracer values on input; new values on output

    integer, dimension(:), intent(in) ::  &
         indexToCellID,  & !< Input: global index for each edge on a block (diagnostic only)
         indexToEdgeID     !< Input: global index for each edge on a block (diagnostic only)

    type(block_type), intent(in) :: &
         block             !< Input: local block (diagnostic only)

    ! local variables

    type(tracer_type), pointer :: &
         thisTracer,   &   ! pointer that loops through linked list of tracers
         parentTracer, &   ! pointer to parent of a given tracer
         dummyTracer       ! dummy tracer with value of 1 everywhere

    integer :: iCell, iEdge, iEdgeOnCell, iLayer, iCat

    integer :: &
         edgeSignOnCell    ! = 1 or -1, depending on which side of the edge a triangle is located

    integer :: &
         nCategories, &    ! number of categories
         nLayers           ! number of layers

    real(kind=RKIND), dimension(:), allocatable :: &
         fluxFromCell2D      ! net flux*dt out of a grid cell ('2D' suffix refers to tracer array dimension)

    real(kind=RKIND), dimension(:,:), allocatable :: &
         fluxFromCell3D      ! net flux*dt out of a grid cell ('3D' suffix refers to tracer array dimension)

    character(len=strKIND) :: &
         errorMessage ! error message for abort

    ! initial loop through linked list of tracers

    thisTracer => tracersHead
    do while (associated(thisTracer))

       ! get dimension sizes

       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)

       elseif (thisTracer % ndims == 3) then

          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)

       endif

       ! set parent tracer
       if (thisTracer % nParents == 0) then   ! mass-like field (fractional ice concentration for CICE)

          ! Set a parentTracer pointer to have a mass-tracer product of 1 in each cell.
          ! This ensures correct values for massTracerProduct below.
          allocate(dummyTracer)
          if (thisTracer % ndims == 2) then

             allocate(dummyTracer % massTracerProduct2D(nCategories,nCellsSolve))
             dummyTracer % massTracerProduct2D(:,:) = 1.0_RKIND
             dummyTracer % ndims = 2

          elseif (thisTracer % ndims == 3) then

             allocate(dummyTracer % massTracerProduct3D(nLayers,nCategories,nCellsSolve))
             dummyTracer % massTracerProduct3D(:,:,:) = 1.0_RKIND
             dummyTracer % ndims = 3

          endif  ! ndims

          parentTracer => dummyTracer

       else  ! nParents = 1, 2 or 3

          parentTracer => thisTracer % parent

       endif ! nParents

       thisTracer => thisTracer % next

    enddo  ! while(associated)

    ! update loop through linked list of tracers
    ! Note: Must begin with up-to-date values of massTracerProduct from previous loop

    thisTracer => tracersHead
    do while (associated(thisTracer))

       ! set parent tracer
       if (thisTracer % nParents == 0) then   ! mass-like field (fractional ice concentration for CICE)

          parentTracer => dummyTracer

       else  ! nParents = 1, 2 or 3

          parentTracer => thisTracer % parent

       endif ! nParents

       if (thisTracer % ndims == 2) then

          ! get size of 1st dimension
          nCategories = size(thisTracer % array2D, 1)

          allocate (fluxFromCell2D(nCategories))

          do iCell = 1, nCellsSolve

             fluxFromCell2D = 0.0_RKIND

             if (ctestOnProc .and. block % localBlockID == ctestBlockID) then
                if (verboseFluxes .and. iCell == ctest .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
                   write(stderrUnit,*) ' '
                   write(stderrUnit,*) 'Update: iCell, tracer =', indexToCellID(iCell), trim(thisTracer % tracerName)
                   write(stderrUnit,*) 'Edge, flux from cell:'
                endif
             endif

             do iEdgeOnCell = 1, nEdgesOnCell(iCell)

                ! find the index for this edge
                iEdge = edgesOnCell(iEdgeOnCell,iCell)

                ! choose the correct sign
                ! If this is cell 1, then a positive flux across an edge removes mass;
                !  else a positive flux adds mass
                if (iCell == cellsOnEdge(1,iEdge)) then
                   edgeSignOnCell = 1
                else
                   edgeSignOnCell = -1
                end if

                ! acccumulate the edge flux

                fluxFromCell2D(:) = fluxFromCell2D(:) + thisTracer % edgeFlux2D(:,iEdge) * edgeSignOnCell

                if (verboseFluxes .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
                   if (iCell == ctest .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
                      write(stderrUnit,*) indexToEdgeID(iEdge), thisTracer % edgeFlux2D(iCatTest,iEdge)*edgeSignOnCell
                   endif
                endif

             enddo   ! nEdgesOnCell

             ! update the tracer
             ! Note: This expression uses the old value of massTracerProduct for this tracer,
             ! but the new value for the parent tracer
             do iCat = 1, nCategories
                if (parentTracer % massTracerProduct2D(iCat,iCell) > 0.0_RKIND) then
                   thisTracer % array2D(iCat,iCell) = &
                        (thisTracer % massTracerProduct2D(iCat,iCell) - (fluxFromCell2D(iCat) / areaCell(iCell))) / &
                        parentTracer % massTracerProduct2D(iCat,iCell)
                else
                   thisTracer % array2D(iCat,iCell) = 0.0_RKIND
                endif
             enddo

             if (verboseFluxes .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
                if (iCell == ctest .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
                   write(stderrUnit,*) 'cell area:', areaCell(iCell)
                   write(stderrUnit,*) 'iceAreaCategory, massTracerProduct:', thisTracer % massTracerProduct2D(iCatTest,iCell)
                   write(stderrUnit,*) 'flux/area:', fluxFromCell2D(iCatTest) / areaCell(iCell)
                   write(stderrUnit,*) 'difference:', &
                        (thisTracer % massTracerProduct2D(iCatTest,iCell) - (fluxFromCell2D(iCatTest) / areaCell(iCell)))
                   write(stderrUnit,*) 'parent massTracerProduct:', parentTracer % massTracerProduct2D(iCatTest,iCell)
                   write(stderrUnit,*) 'new value:', thisTracer % array2D(iCatTest,iCell)
                endif
             endif

             ! update the mass-tracer product, as needed for child tracers
             !TODO - Do this only if hasChild = T?

             thisTracer % massTracerProduct2D(:,iCell) = parentTracer % massTracerProduct2D(:,iCell)  &
                                                         * thisTracer % array2D(:,iCell)

          enddo  ! nCells

          deallocate (fluxFromCell2D)

       elseif (thisTracer % ndims == 3) then

          ! get size of 1st and 2nd dimensions
          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)

          allocate (fluxFromCell3D(nLayers,nCategories))

          do iCell = 1, nCellsSolve

             fluxFromCell3D = 0.0_RKIND

             if (ctestOnProc .and. block % localBlockID == ctestBlockID) then
                if (verboseFluxes .and. iCell == ctest .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
                   write(stderrUnit,*) ' '
                   write(stderrUnit,*) 'Update: iCell, tracer =', indexToCellID(iCell), trim(thisTracer % tracerName)
                   write(stderrUnit,*) 'Edge, flux from cell:'
                endif
             endif

             do iEdgeOnCell = 1, nEdgesOnCell(iCell)

                ! find the index for this edge
                iEdge = edgesOnCell(iEdgeOnCell, iCell)

                ! choose the correct sign
                ! If this is cell 1, then a positive flux across an edge removes mass;
                !  else a positive flux adds mass
                if (iCell == cellsOnEdge(1,iEdge)) then
                   edgeSignOnCell = 1
                else
                   edgeSignOnCell = -1
                end if

                ! acccumulate the edge flux

                fluxFromCell3D(:,:) = fluxFromCell3D(:,:) + thisTracer % edgeFlux3D(:,:,iEdge) * edgeSignOnCell

                if (verboseFluxes .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
                   if (iCell == ctest .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
                      write(stderrUnit,*) indexToEdgeID(iEdge), thisTracer % edgeFlux3D(iLayerTest,iCatTest,iEdge)*edgeSignOnCell
                   endif
                endif

             enddo  ! iEdgeOnCell

             ! update the tracer
             ! Note: This expression uses the old value of massTracerProduct for this tracer,
             ! but the new value for the parent tracer

             if (parentTracer % ndims == 2) then

                do iCat = 1, nCategories
                   do iLayer = 1, nLayers
                      if (parentTracer % massTracerProduct2D(iCat,iCell) > 0.0_RKIND) then
                         thisTracer % array3D(iLayer,iCat,iCell) = (thisTracer % massTracerProduct3D(iLayer,iCat,iCell) &
                                                                       - (fluxFromCell3D(iLayer,iCat) / areaCell(iCell))) &
                                                                  / parentTracer % massTracerProduct2D(iCat,iCell)
                      else
                         thisTracer % array3D(iLayer,iCat,iCell) = 0.0_RKIND
                      endif
                   enddo
                enddo

                if (verboseFluxes .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
                   if (iCell == ctest .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
                      write(stderrUnit,*) 'cell area:', areaCell(iCell)
                      write(stderrUnit,*) 'massTracerProduct:', thisTracer % massTracerProduct3D(iLayerTest,iCatTest,iCell)
                      write(stderrUnit,*) 'parent massTracerProduct:', parentTracer % massTracerProduct2D(iCatTest,iCell)
                      write(stderrUnit,*) 'new value:', thisTracer % array3D(iLayerTest,iCatTest,iCell)
                   endif
                endif

                ! update the mass-tracer product, as needed for child tracers
                !TODO - Do this only if hasChild = T?

                do iCat = 1, nCategories

                   thisTracer % massTracerProduct3D(:,iCat,iCell) = parentTracer % massTracerProduct2D(iCat,iCell)  &
                                                                    * thisTracer % array3D(:,iCat,iCell)

                enddo

             elseif (parentTracer % ndims == 3) then

                do iCat = 1, nCategories
                   do iLayer = 1, nLayers
                      if (parentTracer % massTracerProduct3D(iLayer,iCat,iCell) > 0.0_RKIND) then
                         thisTracer % array3D(iLayer,iCat,iCell) = (thisTracer % massTracerProduct3D(iLayer,iCat,iCell) &
                                                                       - (fluxFromCell3D(iLayer,iCat) / areaCell(iCell))) &
                                                                  / parentTracer % massTracerProduct3D(iLayer,iCat,iCell)
                      else
                         thisTracer % array3D(iLayer,iCat,iCell) = 0.0_RKIND
                      endif
                   enddo
                enddo

                if (verboseFluxes .and. ctestOnProc .and. block % localBlockID == ctestBlockID) then
                   if (iCell == ctest .and. trim(thisTracer % tracerName) == 'iceAreaCategory') then
                      write(stderrUnit,*) 'cell area:', areaCell(iCell)
                      write(stderrUnit,*) 'massTracerProduct:', thisTracer % massTracerProduct3D(iLayerTest,iCatTest,iCell)
                      write(stderrUnit,*) 'parent massTracerProduct:', &
                           parentTracer % massTracerProduct3D(iLayerTest,iCatTest,iCell)
                      write(stderrUnit,*) 'new value:', thisTracer % array3D(iLayerTest,iCatTest,iCell)
                   endif
                endif

                ! update the mass-tracer product, as needed for child tracers
                !TODO - Do this only if hasChild = T?

                thisTracer % massTracerProduct3D(:,:,iCell) = parentTracer % massTracerProduct3D(:,:,iCell)  &
                                                              * thisTracer % array3D(:,:,iCell)

             endif   ! parentTracer % ndims

          enddo  ! nCells

          deallocate (fluxFromCell3D)

       endif

       ! Check for negative values

       if (thisTracer % nParents == 0) then   ! mass-like field

          if (thisTracer % ndims == 2) then

             do iCell = 1, nCellsSolve
                do iCat = 1, nCategories
                   if (thisTracer % array2D(iCat,iCell) < 0.0_RKIND) then
                      write(stderrUnit,*) ' '
                      write(errorMessage,*) 'Negative mass in IR: iCat, iCell, value:', &
                           iCat, iCell, thisTracer % array2D(iCat,iCell)
                      write(stderrUnit,*) 'ERROR: '//trim(errorMessage)
                      call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
                   endif
                enddo
             enddo

          elseif (thisTracer % ndims == 3) then

             do iCell = 1, nCellsSolve
                do iCat = 1, nCategories
                   do iLayer = 1, nLayers
                      if (thisTracer % array3D(iLayer,iCat,iCell) < 0.0_RKIND) then
                         write(stderrUnit,*) ' '
                         write(errorMessage,*) 'Negative mass in IR: iLayer, iCat, iCell, value:',  &
                              iLayer, iCat, iCell, thisTracer % array3D(iLayer,iCat,iCell)
                         write(stderrUnit,*) 'ERROR: '//trim(errorMessage)
                         call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
                      endif
                   enddo
                enddo
             enddo

          endif   ! ndims

       endif      ! nParents

       ! clean up
       if (associated(dummyTracer)) then
          if (dummyTracer % ndims == 2) then
             if (associated(dummyTracer % massTracerProduct2D)) deallocate(dummyTracer % massTracerProduct2D)
          elseif (dummyTracer % ndims == 3) then
             if (associated(dummyTracer % massTracerProduct3D)) deallocate(dummyTracer % massTracerProduct3D)
          endif
          deallocate(dummyTracer)
       endif

       thisTracer => thisTracer % next

    enddo   ! while(associated)

  end subroutine update_mass_and_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine recover_tracer_means
!
!> \brief MPAS-CICE check whether reconstructed tracer fields recover the means
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  This routine integrates mass and mass*tracer for each cell to verify
!>  that the integrals recover the mean values.
!>  It is just a check; it is not required for production runs.
!>  Note: It is coded with embedded tracer loops, such that is has a small memory footprint
!>        but runs very slowly.
!
!-----------------------------------------------------------------------

  subroutine recover_tracer_means(&
       nCells,                        &
       maxEdges,                      &
       nEdgesOnCell,                  &
       xVertexOnCell, yVertexOnCell,  &
       edgesOnCell,                   &
       dcEdge,                        &
       dvEdge, &
       tracersHead)

    integer, intent(in) ::  &
         nCells,        &  !< Input: number of cells
         maxEdges          !< Input: maximum number of edges per cell

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell      !< Input: number of edges for each cell

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         xVertexOnCell,  & !< Input: x (east) coordinate of vertex relative to cell center in local tangent plane
         yVertexOnCell     !< Input: y (north) coordinate of vertex relative to cell center in local tangent plane

    integer, dimension(:,:), intent(in) ::  &
         edgesOnCell       !< Input: edge index for each edge of a given cell

    real(kind=RKIND), dimension(:), intent(in) :: &
         dcEdge,        &  !< Input: distance between the 2 cell centers on each side of an edge
         dvEdge            !< Input: distance between the 2 vertices of an edge

    type(tracer_type), pointer :: &
         tracersHead       !< Input/output: pointer to first element of linked list of tracers
                           ! The pointer stays attached to the first tracer, but all tracers are updated
    ! local variables

    type(tracer_type), pointer :: &
         thisTracer,        &   ! pointer that loops through linked list of tracers
         parentTracer,      &   ! parent of thisTracer
         dummyTracer            ! dummy tracer with barycenter at cell centroid

    real(kind=RKIND), dimension(maxEdges) :: &
         fracEdgeArea      ! fractional area of the triangle joining the cell center to the 2 vertices of each edge
                           ! (relative to the total cell area)
    real(kind=RKIND) :: &
         sumArea           ! sum of edge areas for a cell

    real(kind=RKIND) :: &
         x1Vertex, y1Vertex, x2Vertex, y2Vertex, x3Vertex, y3Vertex  ! coordinates of triangle vertices

    real(kind=RKIND) :: difference  ! difference between two values

    real(kind=RKIND), dimension(nQuadPoints) ::  &
         xQP, yQP,      &  ! x and y coordinates of quadrature points in each triangle
         wQP               ! weight assigned to each quadrature point

    integer :: iCell, iEdge, iCat, iLayer, iEdgeOnCell, iVertexOnCell, iqp, nLayers, nCategories

    ! loop over tracers, allocating some local arrays

    thisTracer => tracersHead

    do while (associated(thisTracer))

       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)
          allocate(thisTracer % massTracerProductCell2D(nCategories))
          allocate(thisTracer % massTracerProductTriangle2D(nCategories))
          allocate(thisTracer % valueQP2D(nCategories,6))

       elseif (thisTracer % ndims == 3) then

          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)
          allocate(thisTracer % massTracerProductCell3D(nLayers,nCategories))
          allocate(thisTracer % massTracerProductTriangle3D(nLayers,nCategories))
          allocate(thisTracer % valueQP3D(nLayers,nCategories,6))

       endif   ! ndims

       if (thisTracer % nParents == 0) then   ! mass-like field; set up dummy tracer with value of 1 at each QP

          allocate(dummyTracer)
          dummyTracer % ndims = 2
          allocate(dummyTracer % valueQP2D(nCategories,6))   ! hardwire nQuadPoints = 6
          dummyTracer % valueQP2D(:,:) = 1.0_RKIND
          thisTracer % parent => dummyTracer

       else   ! nParents = 1, 2 or 3

          ! do nothing

       endif  ! nParents

       thisTracer => thisTracer % next

    enddo   ! associated(thisTracer)

    ! Note: The following loop over cells contains embedded loops over tracers. This is not good for performance.
    !       However, since this subroutine is just for bug-checking and can be turned off for production,
    !        I am trading speed for simplicity (and reduced memory footprint).
    !       It could be sped up by defining some additional tracer fields and moving the tracer loops to the outside.

    ! loop over cells
    do iCell = 1, nCells

       if (verboseConstruct .and. ctestOnProc) then
          if (iCell == ctest) then
             write(stderrUnit,*) 'Reconstructed mass*tracer products for cell ', iCell
          endif
       endif

       ! loop over tracers, initializing the mass*tracer products for this cell

       thisTracer => tracersHead

       do while(associated(thisTracer))

          if (thisTracer % ndims == 2) then
             thisTracer % massTracerProductCell2D(:) = 0.0_RKIND
          elseif (thisTracer % ndims == 3) then
             thisTracer % massTracerProductCell3D(:,:) = 0.0_RKIND
          endif

          thisTracer => thisTracer % next

       enddo  ! associate(thisTracer)

       ! loop over edges of the cell
       ! For each edge, compute the area of the triangle joining the cell center to the 2 vertices of each edge
       ! Accumulate the sum

       fracEdgeArea(:) = 0.0_RKIND
       sumArea = 0.0_RKIND

       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(iEdgeOnCell, iCell)
          fracEdgeArea(iEdgeOnCell) = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)
          sumArea = sumArea + fracEdgeArea(iEdgeOnCell)

       enddo

       ! loop over edges again
       ! Each edge is associated with a triangle joining the cell center to the 2 vertices of the edge;
       !  do various calculations on each triangle

       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          ! convert the area of this triangle to a fraction of the total cell area
          fracEdgeArea(iEdgeOnCell) = fracEdgeArea(iEdgeOnCell) / sumArea

          if (verboseConstruct .and. ctestOnProc) then
             if (iCell == ctest) then
                write(stderrUnit,*) ' '
                write(stderrUnit,*) 'Edge triangle, fracEdgeArea:', iEdgeOnCell, fracEdgeArea(iEdgeOnCell)
             endif
          endif

          ! find the coordinates of the triangle vertices associated with this edge
          ! Note: Vertex indices lead edge indices, so edge n has vertices n and n+1.
          x1Vertex = 0.0_RKIND   ! cell center
          y1Vertex = 0.0_RKIND

          iVertexOnCell = iEdgeOnCell   ! one vertex of this edge
          x2Vertex = xVertexOnCell(iVertexOnCell, iCell)
          y2Vertex = yVertexOnCell(iVertexOnCell, iCell)

          iVertexOnCell = iEdgeOnCell+1  ! the other vertex of this edge
          if (iVertexOnCell > nEdgesOnCell(iCell)) iVertexOnCell = 1
          x3Vertex = xVertexOnCell(iVertexOnCell, iCell)
          y3Vertex = yVertexOnCell(iVertexOnCell, iCell)

          if (verboseConstruct .and. ctestOnProc) then
             if (iCell == ctest) then
                write(stderrUnit,*) 'Vertices:', x1Vertex, y1Vertex, x2Vertex, y2Vertex, x3Vertex, y3Vertex
             endif
          endif

          ! find the coordinates of 6 quadrature points
          ! Note: 6 points are sufficient for exact integration of polynomials up to 4th order

          xQP(1) = q1TriangleQP * x1Vertex + q1TriangleQP * x2Vertex + q2TriangleQP * x3Vertex
          yQP(1) = q1TriangleQP * y1Vertex + q1TriangleQP * y2Vertex + q2TriangleQP * y3Vertex

          xQP(2) = q1TriangleQP * x1Vertex + q2TriangleQP * x2Vertex + q1TriangleQP * x3Vertex
          yQP(2) = q1TriangleQP * y1Vertex + q2TriangleQP * y2Vertex + q1TriangleQP * y3Vertex

          xQP(3) = q2TriangleQP * x1Vertex + q1TriangleQP * x2Vertex + q1TriangleQP * x3Vertex
          yQP(3) = q2TriangleQP * y1Vertex + q1TriangleQP * y2Vertex + q1TriangleQP * y3Vertex

          xQP(4) = q3TriangleQP * x1Vertex + q4TriangleQP * x2Vertex + q4TriangleQP * x3Vertex
          yQP(4) = q3TriangleQP * y1Vertex + q4TriangleQP * y2Vertex + q4TriangleQP * y3Vertex

          xQP(5) = q4TriangleQP * x1Vertex + q3TriangleQP * x2Vertex + q4TriangleQP * x3Vertex
          yQP(5) = q4TriangleQP * y1Vertex + q3TriangleQP * y2Vertex + q4TriangleQP * y3Vertex

          xQP(6) = q4TriangleQP * x1Vertex + q4TriangleQP * x2Vertex + q3TriangleQP * x3Vertex
          yQP(6) = q4TriangleQP * y1Vertex + q4TriangleQP * y2Vertex + q3TriangleQP * y3Vertex

          ! set weight of each quadrature point
          wQP(1:3) = w1TriangleQP
          wQP(4:6) = w2TriangleQP

          ! loop over tracers

          thisTracer => tracersHead

          do while (associated(thisTracer))

             if (verboseConstruct .and. ctestOnProc) then
                if (iCell == ctest) then

                endif
             endif

             parentTracer => thisTracer % parent

             ! Integrate mass*tracer over this triangle by computing a sum of mass*tracer over quadrature points

             if (thisTracer % ndims == 2) then

                ! initialize
                thisTracer % massTracerProductTriangle2D(:) = 0.0_RKIND

                ! sum over quadrature points
                ! For each QP, evaluate this tracer at the quadrature point, then multiply by the value of the parent tracer.
                ! For nParents = 0, this is the value of the mass-like field at the QP.
                ! For nParents = 1, it is mass*tracer1
                ! For nParents = 2, it is mass*tracer1*tracer2
                ! For nParents = 3, it is mass*tracer1*tracer2*tracer3

                do iqp = 1, 6

                   thisTracer % valueQP2D(:,iqp) = thisTracer % center2D(:,iCell)         &
                                                 + thisTracer % xGrad2D(:,iCell)*xQP(iqp) &
                                                 + thisTracer % yGrad2D(:,iCell)*yQP(iqp)

                   ! reset the current tracer value to the product of the tracer and its parent(s)
                   ! Note: The parent of a 2D tracer always has ndims = 2
                   thisTracer % valueQP2D(:,iqp) = thisTracer % valueQP2D(:,iqp) * parentTracer % valueQP2D(:,iqp)

                   ! increment the sum over QPs
                   thisTracer % massTracerProductTriangle2D(:) = thisTracer % massTracerProductTriangle2D(:)    &
                                                               + wQP(iqp) * thisTracer % valueQP2D(:,iqp)

                enddo  ! iqp

                if (verboseConstruct .and. ctestOnProc) then
                   if (iCell == ctest) then
                      write(stderrUnit,*) '   Tracer, triangle value: ', &
                           trim(thisTracer % tracerName), thisTracer % massTracerProductTriangle2D(iCatTest)
                   endif
                endif

                ! add this triangle contribution to the cumulative sum for the cell,
                ! weighted by the fractional area of the triangle

                thisTracer % massTracerProductCell2D(:) = thisTracer % massTracerProductCell2D(:)   &
                                                        + fracEdgeArea(iEdgeOnCell) * thisTracer % massTracerProductTriangle2D(:)

             elseif (thisTracer % ndims == 3) then

                ! initialize
                thisTracer % massTracerProductTriangle3D(:,:) = 0.0_RKIND

                ! sum over quadrature points
                ! For each QP, evaluate this tracer at the quadrature point, then multiply by the value of the parent tracer.
                ! For nParents = 0, this is the value of the mass-like field at the QP.
                ! For nParents = 1, it is mass*tracer1
                ! For nParents = 2, it is mass*tracer1*tracer2
                ! For nParents = 3, it is mass*tracer1*tracer2*tracer3

                do iqp = 1, 6

                   thisTracer % valueQP3D(:,:,iqp) = thisTracer % center3D(:,:,iCell)         &
                                                   + thisTracer % xGrad3D(:,:,iCell)*xQP(iqp) &
                                                   + thisTracer % yGrad3D(:,:,iCell)*yQP(iqp)

                   ! Note: parent tracer can have ndims = 2 or 3

                   if (parentTracer % ndims == 2) then

                      do iLayer = 1, nLayers

                         ! reset the current tracer value to the product of the tracer and its parent(s)
                         thisTracer % valueQP3D(iLayer,:,iqp) = &
                              thisTracer % valueQP3D(iLayer,:,iqp) * parentTracer % valueQP2D(:,iqp)

                         ! increment the sum over QPs
                         thisTracer % massTracerProductTriangle3D(iLayer,:) = &
                              thisTracer % massTracerProductTriangle3D(iLayer,:) &
                            + wQP(iqp) * thisTracer % valueQP3D(iLayer,:,iqp)

                      enddo   ! nLayers

                   elseif (parentTracer % ndims == 3) then

                      ! reset the current tracer value to the product of the tracer and its parent(s)
                      thisTracer % valueQP3D(:,:,iqp) = thisTracer % valueQP3D(:,:,iqp) * parentTracer % valueQP3D(:,:,iqp)

                      ! increment the sum over QPs
                      thisTracer % massTracerProductTriangle3D(:,:) = thisTracer % massTracerProductTriangle3D(:,:) &
                                                                    + wQP(iqp) * thisTracer % valueQP3D(:,:,iqp)

                   endif  ! parentTracer % ndims

                enddo  ! iqp

                if (verboseConstruct .and. ctestOnProc) then
                   if (iCell == ctest) then
                      write(stderrUnit,*) '   Tracer, triangle value:', &
                           trim(thisTracer % tracerName), thisTracer % massTracerProductTriangle3D(iLayerTest,iCatTest)
                   endif
                endif

                ! add this triangle contribution to the cumulative sum for the cell,
                ! weighted by the fractional area of the triangle

                thisTracer % massTracerProductCell3D(:,:) = &
                     thisTracer % massTracerProductCell3D(:,:) &
                   + fracEdgeArea(iEdgeOnCell) * thisTracer % massTracerProductTriangle3D(:,:)

             endif   ! thisTracer % ndims

             thisTracer => thisTracer % next

          enddo   ! associated(thisTracer)

       enddo  ! iEdgeOnCell

       ! We now have the integral of mass*tracer over the cell for each tracer
       ! Loop over tracers again, and compare to the known mass*tracer product computed above.

       thisTracer => tracersHead   ! start by pointing to the mass-like field

       do while (associated(thisTracer))

          if (thisTracer % ndims == 2) then

             do iCat = 1, nCategories

                !TODO - Figure out a criterion based on nParents?
                ! Relative errors are largest for tracers with 3 parents (e.g., pond depth)
                difference = thisTracer % massTracerProductCell2D(iCat) - thisTracer % massTracerProduct2D(iCat,iCell)
                if (abs(difference) > eps11 * abs(thisTracer % massTracerProduct2D(iCat,iCell)) .and. abs(difference) > eps11) then
                   write(stderrUnit,*) 'WARNING: Reconstructed tracer does not recover the mean, iCat, iCell, tracer =', &
                        iCat, iCell, trim(thisTracer % tracerName)
                   write(stderrUnit,*) 'Product of means:     ', thisTracer % massTracerProduct2D(iCat,iCell)
                   write(stderrUnit,*) 'Reconstructed product:', thisTracer % massTracerProductCell2D(iCat)
                   write(stderrUnit,*) 'Difference:           ', difference
                endif

             enddo

             if (verboseConstruct .and. ctestOnProc) then
                if (iCell == ctest) then
                   write(stderrUnit,*) '   Tracer, reconstructed and desired value of mass*tracer: ', &
                        trim(thisTracer % tracerName),  &
                        thisTracer % massTracerProductCell2D(iCatTest), &
                        thisTracer % massTracerProduct2D(iCatTest,iCell)
                endif
             endif

          elseif (thisTracer % ndims == 3) then

             do iCat = 1, nCategories
                do iLayer = 1, nLayers

                   difference = &
                        thisTracer % massTracerProductCell3D(iLayer,iCat) - thisTracer % massTracerProduct3D(iLayer,iCat,iCell)
                   if (abs(difference) > eps11 * abs(thisTracer % massTracerProduct3D(iLayer,iCat,iCell)) .and. &
                        abs(difference) > eps11) then
                      write(stderrUnit,*) 'WARNING: Reconstructed tracer does not recover the mean, iCat, iCell, tracer =', &
                           iCat, iCell, trim(thisTracer % tracerName)
                      write(stderrUnit,*) 'Product of means:     ', thisTracer % massTracerProduct3D(iLayer,iCat,iCell)
                      write(stderrUnit,*) 'Reconstructed product:', thisTracer % massTracerProductCell3D(iLayer,iCat)
                      write(stderrUnit,*) 'Difference:           ', difference
                   endif

                enddo   ! iLayer
             enddo   ! iCat

             if (verboseConstruct .and. ctestOnProc) then
                if (iCell == ctest) then
                   write(stderrUnit,*) '   Tracer, reconstructed and desired value of  mass*tracer: ', &
                        trim(thisTracer % tracerName),  &
                        thisTracer % massTracerProductCell3D(iLayerTest,iCatTest), &
                        thisTracer % massTracerProduct3D(iLayerTest,iCatTest,iCell)
                endif
             endif

          endif  ! ndims

          thisTracer => thisTracer % next

       enddo   ! associated(thisTracer)

    enddo     ! iCell

    ! clean up

    thisTracer => tracersHead

    do while (associated(thisTracer))

       if (thisTracer % ndims == 2) then

          deallocate(thisTracer % massTracerProductTriangle2D)
          deallocate(thisTracer % massTracerProductCell2D)
          deallocate(thisTracer % valueQP2D)

       elseif (thisTracer % ndims == 3) then

          deallocate(thisTracer % massTracerProductTriangle3D)
          deallocate(thisTracer % massTracerProductCell3D)
          deallocate(thisTracer % valueQP3D)

       endif  ! ndims

       if (associated(dummyTracer)) then
          if (dummyTracer % ndims == 2) then
             if (associated(dummyTracer % valueQP2D)) deallocate(dummyTracer % valueQP2D)
          elseif (dummyTracer % ndims == 3) then
             if (associated(dummyTracer % valueQP3D)) deallocate(dummyTracer % valueQP3D)
          endif
          deallocate(dummyTracer)
       endif

       thisTracer => thisTracer % next

    enddo   ! associated(thisTracer)

  end subroutine recover_tracer_means

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine sum_tracers
!
!> \brief  sum mass*tracer products over all sums in block
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  This routine computes sums of mass, mass*tracer, mass*tracer*tracer, etc.
!>  over all grid cells in the block.
!>  NOTE: This routine should be called after compute_mass_tracer_products.
!
!-----------------------------------------------------------------------

  subroutine sum_tracers(&
       tracersHead,   &
       nCellsSolve,   &
       areaCell,      &
       init)

    type(tracer_type), pointer :: &
         tracersHead  !< Input/output: pointer to first element of linked list of tracers
                      ! The pointer stays attached to the first tracer, but all tracer sums are updated

    integer, intent(in) ::  &
         nCellsSolve  !< Input: number of locally owned cells

    real(kind=RKIND), dimension(:), intent(in) ::   &
         areaCell     !< Input: grid cell area

    logical, intent(in) :: &
         init         !< Input: If true, then compute initial tracer sums
                      !         Else compute final tracer sums

    ! local variables
    type(tracer_type), pointer :: thisTracer

    integer :: nCategories, nLayers

    integer :: iCell

    integer :: iCat, iLayer

    real(kind=RKIND), dimension(:), allocatable ::  &
         tracerSum2D

    real(kind=RKIND), dimension(:,:), allocatable ::  &
         tracerSum3D

    thisTracer => tracersHead

    do while (associated(thisTracer))

       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)
          allocate(tracerSum2D(nCategories))

          tracerSum2D(:) = 0.0_RKIND
          do iCell = 1, nCellsSolve
             do iCat = 1, nCategories
                tracerSum2D(iCat) = tracerSum2D(iCat) + areaCell(iCell) * thisTracer % massTracerProduct2D(iCat,iCell)
             enddo
          enddo

          if (init) then
             thisTracer % globalSumInit2D(:) = tracerSum2D(:)
          else
             thisTracer % globalSumFinal2D(:) = tracerSum2D(:)
          endif

          deallocate(tracerSum2D)

          if (verboseRun) then
             if (trim(thisTracer % tracerName) == 'iceAreaCategory' .or. trim(thisTracer % tracerName) == 'iceVolumeCategory') then
                write(stderrUnit,*) ' '
                write(stderrUnit,*) 'Check conservation: ', trim(thisTracer % tracerName)
                if (init) then
                   write(stderrUnit,*) 'Initial sum on block:', thisTracer % globalSumInit2D(1)
                else
                   write(stderrUnit,*) 'Final sum on block:', thisTracer % globalSumFinal2D(1)
                endif
             endif
          endif   ! verbose

       elseif (thisTracer % ndims == 3) then

          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)
          allocate(tracerSum3D(nLayers, nCategories))
          tracerSum3D(:,:) = 0.0_RKIND

          do iCell = 1, nCellsSolve
             do iCat = 1, nCategories
                do iLayer = 1, nLayers
                   tracerSum3D(iLayer,iCat) = &
                        tracerSum3D(iLayer,iCat) + areaCell(iCell) * thisTracer % massTracerProduct3D(iLayer,iCat,iCell)
                enddo
             enddo
          enddo

          if (init) then
             thisTracer % globalSumInit3D(:,:) = tracerSum3D(:,:)
          else
             thisTracer % globalSumFinal3D(:,:) = tracerSum3D(:,:)
          endif

          deallocate(tracerSum3D)

          if (verboseRun) then
             if (trim(thisTracer % tracerName) == 'iceAreaCategory' .or. trim(thisTracer % tracerName) == 'iceVolumeCategory') then
                write(stderrUnit,*) ' '
                write(stderrUnit,*) 'Check conservation: ', trim(thisTracer % tracerName)
                if (init) then
                   write(stderrUnit,*) 'Initial sum on block:', thisTracer % globalSumInit3D(1,1)
                else
                   write(stderrUnit,*) 'Final sum on block:', thisTracer % globalSumFinal3D(1,1)
                endif
             endif
          endif   ! verbose

       endif   ! ndims

       thisTracer => thisTracer % next

    enddo   ! associated(thisTracer)

  end subroutine sum_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine check_tracer_conservation
!
!> \brief  check for conservation of global mass*tracer
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  This routine checks that mass*tracer values are globally conserved.
!   NOTE: Includes global sums, so cannot be called inside a block loop.
!
!-----------------------------------------------------------------------

  subroutine check_tracer_conservation(dminfo, tracersHead)

    type (dm_info), intent(in) :: dminfo   !< Input: domain info

    type(tracer_type), pointer :: &
         tracersHead    !< Input/output: pointer to first element of linked list of tracers
                        ! The pointer stays attached to the first tracer, but tracer sums are updated
    ! local variables
    type(tracer_type), pointer :: thisTracer

    real(kind=RKIND) :: difference, ratio, localSum

    integer :: nCategories, nLayers

    integer :: iCat, iLayer

    character(len=strKIND) :: &
         errorMessage ! error message for abort

    thisTracer => tracersHead

    do while (associated(thisTracer))

       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)

          do iCat = 1, nCategories

             ! global sums
             call mpas_timer_start("tracer cons check 2D glbsum")
             localSum = thisTracer % globalSumInit2D(iCat)
             call mpas_dmpar_sum_real (dminfo, localSum, thisTracer % globalSumInit2D(iCat))
             localSum = thisTracer % globalSumFinal2D(iCat)
             call mpas_dmpar_sum_real (dminfo, localSum, thisTracer % globalSumFinal2D(iCat))
             call mpas_timer_stop("tracer cons check 2D glbsum")

             if (verboseRun .and. iCat == 1) then
                if (trim(thisTracer % tracerName) == 'iceAreaCategory' .or. &
                    trim(thisTracer % tracerName) == 'iceVolumeCategory') then
                   write(stderrUnit,*) ' '
                   write(stderrUnit,*) 'Tracer conservation, init global sum, final global sum: ', trim(thisTracer % tracerName), &
                        thisTracer % globalSumInit2D(iCat), thisTracer % globalSumFinal2D(iCat)
                endif
             endif

             ! check for equality of initial and final sums

             if (abs(thisTracer % globalSumInit2D(iCat)) > eps11) then

                difference = thisTracer % globalSumFinal2D(iCat) - thisTracer % globalSumInit2D(iCat)
                ratio = difference / thisTracer % globalSumInit2D(iCat)

                !TODO - Verify that eps11 will work as a threshold
                if (abs(ratio) > eps11) then
                   write(errorMessage,*) 'Tracer conservation error: tracer, category =', trim(thisTracer % tracerName), iCat
                   write(stderrUnit,*) trim(errorMessage)
                   write(stderrUnit,*) 'Init value, final value:', &
                        thisTracer % globalSumInit2D(iCat), thisTracer % globalSumFinal2D(iCat)
                   write(stderrUnit,*) 'Difference, ratio:', difference, ratio
                   call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
                endif

             endif

          enddo    ! iCat

       elseif (thisTracer % ndims == 3) then

          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)

          do iCat = 1, nCategories
             do iLayer = 1, nLayers

                ! global sums
                call mpas_timer_start("tracer cons check 3D glbsum")
                localSum = thisTracer % globalSumInit3D(iLayer,iCat)
                call mpas_dmpar_sum_real (dminfo, localSum, thisTracer % globalSumInit3D(iLayer,iCat))
                localSum = thisTracer % globalSumFinal3D(iLayer,iCat)
                call mpas_dmpar_sum_real (dminfo, localSum, thisTracer % globalSumFinal3D(iLayer,iCat))
                call mpas_timer_stop("tracer cons check 3D glbsum")

                if (verboseRun .and. iCat == 1 .and. iLayer == 1) then
                   if (trim(thisTracer % tracerName) == 'iceAreaCategory' .or. &
                       trim(thisTracer % tracerName) == 'iceVolumeCategory') then
                      write(stderrUnit,*) ' '
                      write(stderrUnit,*) 'Tracer conservation, init global sum, final global sum: ', &
                           trim(thisTracer % tracerName), &
                           thisTracer % globalSumInit3D(iLayer,iCat), thisTracer % globalSumFinal3D(iLayer,iCat)
                   endif
                endif

                if (abs(thisTracer % globalSumInit3D(iLayer,iCat)) > eps11) then

                   difference = thisTracer % globalSumFinal3D(iLayer,iCat) - thisTracer % globalSumInit3D(iLayer,iCat)
                   ratio = difference/ thisTracer % globalSumInit3D(iLayer,iCat)

                   !TODO - Verify that eps11 will work as a threshold
                   if (abs(ratio) > eps11) then
                      write(errorMessage,*) 'Tracer conservation error: tracer, category, layer =', &
                           trim(thisTracer % tracerName), iCat, iLayer
                      write(stderrUnit,*) trim(errorMessage)
                      write(stderrUnit,*) 'Init value, final value:', &
                           thisTracer % globalSumInit3D(iLayer,iCat), thisTracer % globalSumFinal3D(iLayer,iCat)
                      write(stderrUnit,*) 'Difference, ratio:', difference, ratio
                      call mpas_dmpar_global_abort('MPAS-seaice: IR: '//trim(errorMessage))
                   endif

                endif

             enddo   ! iLayer
          enddo   ! iCat

       endif  ! ndims

       thisTracer => thisTracer % next

    enddo     ! while(associated)

  end subroutine check_tracer_conservation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine tracer_local_min_max
!
!> \brief  find the local min and max values for tracers
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  This routine computes for each tracer the maximum and minimum values
!>  in a grid cell and its nearest neighbors.
!
!-----------------------------------------------------------------------

  subroutine tracer_local_min_max(&
       tracersHead,   &
       nCells,        &
       nCellsSolve,   &
       nEdgesOnCell,  &
       cellsOnCell)

    type(tracer_type), pointer :: &
         tracersHead   !< Input/output: pointer to first element of linked list of tracers
                       ! The pointer stays attached to the first tracer, but all tracer sums are updated

    integer, intent(in) ::  &
         nCells,     & !< Input: number of cells
         nCellsSolve   !< Input: number of cells to be updated

    integer, dimension(:), intent(in) ::  &
         nEdgesOnCell  !< Input: number of edges per cell

    integer, dimension(:,:), intent(in) ::  &
         cellsOnCell   !< Input: cell index for each edge neighbor of a given cell

    ! local variables
    type(tracer_type), pointer :: thisTracer

    integer :: nCategories, nLayers

    integer :: iCell, iCellOnCell, iCellNeighbor, iCat, iLayer

    ! This routine is called to support a monotonicity check following advection. If the code is working,
    !  the tracer values in a given cell at the new time should always lie within the range of the tracer
    !  values in that cell's neighborhood at the old time. This routines finds the max and min within
    !  the local neighborhood (a cell plus its nearest neighbors).
    ! If a tracer has a parent tracer with a value of zero (e.g., enthalpy in a category with zero
    !  ice/snow thickness), the tracer value is masked as unmeaningful and is not included in min/max.

    thisTracer => tracersHead

    do while (associated(thisTracer))

       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)

          !TODO - huge and -huge?
          thisTracer % localMin2D(:,:) = 0.0_RKIND
          thisTracer % localMax2D(:,:) = 0.0_RKIND

          do iCell = 1, nCellsSolve
             do iCat = 1, nCategories

                ! initialize to the values in the cell itself
                if (thisTracer % arrayMask2D(iCat,iCell) == 1) then   ! the tracer value in this cell has physical meaning
                   thisTracer % localMin2D(iCat,iCell) = thisTracer % array2D(iCat,iCell)
                   thisTracer % localMax2D(iCat,iCell) = thisTracer % array2D(iCat,iCell)
                endif

                do iCellOnCell = 1, nEdgesOnCell(iCell)

                   iCellNeighbor = cellsOnCell(iCellOnCell, iCell)

                   if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then       ! the neighbor cell exists
                      ! the tracer value in the neighbor cell has physical meaning
                      if (thisTracer % arrayMask2D(iCat,iCellNeighbor) == 1) then
                         if (thisTracer % array2D(iCat,iCellNeighbor) < thisTracer % localMin2D(iCat,iCell)) then
                            thisTracer % localMin2D(iCat,iCell) = thisTracer % array2D(iCat,iCellNeighbor)
                         endif
                         if (thisTracer % array2D(iCat,iCellNeighbor) > thisTracer % localMax2D(iCat,iCell)) then
                            thisTracer % localMax2D(iCat,iCell) = thisTracer % array2D(iCat,iCellNeighbor)
                         endif
                      endif   ! mask = 1
                   endif      ! cell exists

                enddo   ! iCellOnCell

             enddo   ! iCat
          enddo   ! iCell

       elseif (thisTracer % ndims == 3) then

          nLayers = size(thisTracer % array3D, 1)
          nCategories = size(thisTracer % array3D, 2)

          !TODO - huge and -huge?
          thisTracer % localMin3D(:,:,:) = 0.0_RKIND
          thisTracer % localMax3D(:,:,:) = 0.0_RKIND

          do iCell = 1, nCellsSolve
             do iCat = 1, nCategories
                do iLayer = 1, nLayers

                   ! initialize to the values in the cell itself
                   if (thisTracer % arrayMask3D(iLayer,iCat,iCell) == 1) then ! the tracer value in this cell has physical meaning
                      thisTracer % localMin3D(iLayer,iCat,iCell) = thisTracer % array3D(iLayer,iCat,iCell)
                      thisTracer % localMax3D(iLayer,iCat,iCell) = thisTracer % array3D(iLayer,iCat,iCell)
                   endif

                   do iCellOnCell = 1, nEdgesOnCell(iCell)

                      iCellNeighbor = cellsOnCell(iCellOnCell, iCell)

                      if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then       ! the neighbor cell exists
                         ! the tracer value in the neighbor cell has physical meaning
                         if (thisTracer % arrayMask3D(iLayer,iCat,iCellNeighbor) == 1) then
                            if (thisTracer % array3D(iLayer,iCat,iCellNeighbor) < thisTracer % localMin3D(iLayer,iCat,iCell)) then
                               thisTracer % localMin3D(iLayer,iCat,iCell) = thisTracer % array3D(iLayer,iCat,iCellNeighbor)
                            endif
                            if (thisTracer % array3D(iLayer,iCat,iCellNeighbor) > thisTracer % localMax3D(iLayer,iCat,iCell)) then
                               thisTracer % localMax3D(iLayer,iCat,iCell) = thisTracer % array3D(iLayer,iCat,iCellNeighbor)
                            endif
                         endif   ! mask = 1
                      endif      ! cell exists

                   enddo   ! iCellOnCell

                enddo   ! iLayer
             enddo   ! iCat
          enddo   ! iCell

       endif   ! ndims

       thisTracer => thisTracer % next

    enddo   ! while (associated)

  end subroutine tracer_local_min_max

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine check_tracer_monotonicity
!
!> \brief  check tracer monotonicity
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  Given the initial range of tracer values in a cell and its nearest neighbors,
!>  verify that the new values are within this range.
!>  Note: In rare cases, IR gives tracer values outside the range of the local max/min,
!>        but values should always be within the range of the extended max/min
!>        computed here.
!>  Note: This routine includes halo updates and should not be called from
!>        inside a block loop.
!
!-----------------------------------------------------------------------

  !TODO - Recompute the tracer mask before calling this subroutine?
  !       Might not work to do so without the extended min/max, because local min/max
  !        will not be computed correctly for initially ice-free cells

  subroutine check_tracer_monotonicity(domain, tracersHead)

    type (domain_type), intent(inout) :: &
         domain        !< Input: domain info

    type(tracer_type), pointer :: &
         tracersHead   !< Input/output: pointer to first element of linked list of tracers
                       ! The pointer stays attached to the first tracer, but all tracer sums are updated

    ! local variables

    type (mpas_pool_type), pointer :: meshPool
    type(MPAS_pool_type), pointer :: tracerMonotonicityPool

    integer, pointer :: &
         nCells,          & ! number of cells on block
         nCellsSolve        ! number of locally owned cells on block

    integer, dimension(:), pointer :: &
         nEdgesOnCell       ! number of edges for each cell

    integer, dimension(:,:), pointer ::  &
         cellsOnCell        ! cell index for each edge neighbor of a given cell

    type(tracer_type), pointer :: &
         thisTracer

    type(dm_info), pointer :: &
         dminfo

    type(block_type), pointer :: &
         block

    real(kind=RKIND) :: toleranceMin  ! error tolerance below local min
    real(kind=RKIND) :: toleranceMax  ! error tolerance above local max

    integer :: nCategories, nLayers

    integer :: iCell, iCellOnCell, iCellNeighbor, iCat, iLayer

    logical, parameter :: &
         extendedMinMax = .true.  ! if true, then extend the local min/max by one layer of cells
                                  ! may be necessary for complex fields with large CFL

    ! get pools

    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_subpool(domain % blocklist % structs, 'tracer_monotonicity', tracerMonotonicityPool)

    ! get grid quantities

    call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

    ! loop over tracers

    thisTracer => tracersHead

    do while (associated(thisTracer))

       if (thisTracer % nParents > 0) then   ! Note: monotonicity holds for tracers but not for the mass-like field

          if (thisTracer % ndims == 2) then

             nCategories = size(thisTracer % array2D, 1)

             if (extendedMinMax) then      ! extend local max/min by 1 layer of neighbor cells

                ! halo update for local max/min

                call mpas_timer_start("tracer mono check 2D halo")

                call MPAS_dmpar_field_halo_exch(domain, trim(thisTracer % tracerName) // 'LocalMin')
                call MPAS_dmpar_field_halo_exch(domain, trim(thisTracer % tracerName) // 'LocalMax')

                call mpas_timer_stop("tracer mono check 2D halo")

                ! loop over blocks

                block => domain % blocklist

                do while(associated(block))

                   do iCell = 1, nCellsSolve
                      do iCat = 1, nCategories

                         if (thisTracer % arrayMask2D(iCat,iCell) == 1) then   ! the tracer value in this cell has physical meaning

                            do iCellOnCell = 1, nEdgesOnCell(iCell)

                               iCellNeighbor = cellsOnCell(iCellOnCell, iCell)

                               if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then       ! the neighbor cell exists
                                  ! the tracer value in the neighbor cell has physical meaning
                                  if (thisTracer % arrayMask2D(iCat,iCellNeighbor) == 1) then

                                     if (thisTracer % localMin2D(iCat,iCellNeighbor) < thisTracer % localMin2D(iCat,iCell)) then
                                        thisTracer % localMin2D(iCat,iCell) = thisTracer % localMin2D(iCat,iCellNeighbor)
                                     endif

                                     if (thisTracer % localMax2D(iCat,iCellNeighbor) > thisTracer % localMax2D(iCat,iCell)) then
                                        thisTracer % localMax2D(iCat,iCell) = thisTracer % localMax2D(iCat,iCellNeighbor)
                                     endif

                                  endif  ! tracer value in neighbor cell has physical meaning
                               endif     ! neighbor cell exists

                            enddo   ! iCellOnCell

                         endif      ! tracer value has physical meaning

                      enddo   ! iCat
                   enddo   ! iCell

                   block  => block % next

                end do  ! associated(block)

             endif    ! extendedMinMax

             ! loop over blocks

             block => domain % blocklist

             do while(associated(block))

                ! check for monotonicity;
                ! compare the new value in each cell with the max/min of the old values in the cell neighborhood

                do iCell = 1, nCellsSolve
                   do iCat = 1, nCategories

                      if (thisTracer % arrayMask2D(iCat,iCell) == 1) then   ! the tracer value in this cell has physical meaning

                         toleranceMin = eps11 * max(1.0_RKIND, abs(thisTracer % localMin2D(iCat,iCell)))
                         toleranceMax = eps11 * max(1.0_RKIND, abs(thisTracer % localMax2D(iCat,iCell)))

                         if (thisTracer % array2D(iCat,iCell) < thisTracer % localMin2D(iCat,iCell) - toleranceMin) then

                            write(stderrUnit,*) ' '
                            write(stderrUnit,*) 'IR advection, monotonicity violation!'
                            write(stderrUnit,*) 'Tracer, iCat, iCell =', trim(thisTracer % tracerName), iCat, iCell
                            write(stderrUnit,*) 'Old min value in neighborhood =', thisTracer % localMin2D(iCat,iCell)
                            write(stderrUnit,*) 'New value in cell =', thisTracer % array2D(iCat,iCell)
                            write(stderrUnit,*) 'Tolerance, difference:', toleranceMin,  &
                                 thisTracer % localMin2D(iCat,iCell) - thisTracer % array2D(iCat,iCell)
                            call mpas_dmpar_global_abort('MPAS-seaice: IR advection, monotonicity violation (min, nDims == 2)')

                         elseif (thisTracer % array2D(iCat,iCell) > thisTracer % localMax2D(iCat,iCell) + toleranceMax) then

                            write(stderrUnit,*) ' '
                            write(stderrUnit,*) 'IR advection, monotonicity violation!'
                            write(stderrUnit,*) 'Tracer, iCat, iCell =', trim(thisTracer % tracerName), iCat, iCell
                            write(stderrUnit,*) 'Old max value in neighborhood =', thisTracer % localMax2D(iCat,iCell)
                            write(stderrUnit,*) 'New value in cell =', thisTracer % array2D(iCat,iCell)
                            write(stderrUnit,*) 'Tolerance, difference:', toleranceMax,  &
                                 thisTracer % array2D(iCat,iCell) - thisTracer % localMax2D(iCat,iCell)
                            call mpas_dmpar_global_abort('MPAS-seaice: IR advection, monotonicity violation (max, nDims == 2)')

                         endif

                      endif   ! tracer value has physical meaning

                   enddo   ! iCat
                enddo   ! iCell

                block => block % next

             enddo   ! associated(block)

          elseif (thisTracer % ndims == 3) then

             nLayers = size(thisTracer % array3D, 1)
             nCategories = size(thisTracer % array3D, 2)

             if (extendedMinMax) then      ! extend local max/min by 1 layer of neighbor cells

                ! halo update for local max/min

                call mpas_timer_start("tracer mono check 3D halo")

                call MPAS_dmpar_field_halo_exch(domain, trim(thisTracer % tracerName) // 'LocalMin')
                call MPAS_dmpar_field_halo_exch(domain, trim(thisTracer % tracerName) // 'LocalMax')

                call mpas_timer_stop("tracer mono check 3D halo")

                ! loop over blocks

                block => domain % blocklist

                do while(associated(block))

                   do iCell = 1, nCellsSolve
                      do iCat = 1, nCategories
                         do iLayer = 1, nLayers

                            ! the tracer value in this cell has physical meaning
                            if (thisTracer % arrayMask3D(iLayer,iCat,iCell) == 1) then

                               do iCellOnCell = 1, nEdgesOnCell(iCell)

                                  iCellNeighbor = cellsOnCell(iCellOnCell, iCell)

                                  if (iCellNeighbor >= 1 .and. iCellNeighbor <= nCells) then       ! the neighbor cell exists
                                     ! the tracer value in the neighbor cell has physical meaning
                                     if (thisTracer % arrayMask3D(iLayer,iCat,iCellNeighbor) == 1) then

                                        if (thisTracer % localMin3D(iLayer,iCat,iCellNeighbor) < &
                                            thisTracer % localMin3D(iLayer,iCat,iCell)) then
                                           thisTracer % localMin3D(iLayer,iCat,iCell) = &
                                                thisTracer % localMin3D(iLayer,iCat,iCellNeighbor)
                                        endif

                                        if (thisTracer % localMax3D(iLayer,iCat,iCellNeighbor) > &
                                            thisTracer % localMax3D(iLayer,iCat,iCell)) then
                                           thisTracer % localMax3D(iLayer,iCat,iCell) = &
                                                thisTracer % localMax3D(iLayer,iCat,iCellNeighbor)
                                        endif

                                     endif  ! tracer value in neighbor cell has physical meaning
                                  endif     ! neighbor cell exists

                               enddo   ! iCellOnCell

                            endif   ! tracer value has physical meaning

                         enddo   ! iLayer
                      enddo   ! iCat
                   enddo   ! iCell

                   block => block % next

                enddo    ! associated(block)

             endif    ! extendedMinMax

             ! loop over blocks

             block => domain % blocklist

             do while(associated(block))

                ! check for monotonicity;
                ! compare the new value in each cell with the max/min of the old values in the cell neighborhood

                do iCell = 1, nCellsSolve
                   do iCat = 1, nCategories
                      do iLayer = 1, nLayers

                         ! the tracer value in this cell has physical meaning
                         if (thisTracer % arrayMask3D(iLayer,iCat,iCell) == 1) then

                            toleranceMin = eps11 * max(1.0_RKIND, abs(thisTracer % localMin3D(iLayer,iCat,iCell)))
                            toleranceMax = eps11 * max(1.0_RKIND, abs(thisTracer % localMax3D(iLayer,iCat,iCell)))

                            if (thisTracer % array3D(iLayer,iCat,iCell) < &
                                thisTracer % localMin3D(iLayer,iCat,iCell) - toleranceMin) then

                               write(stderrUnit,*) ' '
                               write(stderrUnit,*) 'IR advection, monotonicity violation!'
                               write(stderrUnit,*) 'Tracer, iLayer,iCat, iCell =', &
                                    trim(thisTracer % tracerName), iLayer,iCat, iCell
                               write(stderrUnit,*) 'Old min value in neighborhood =', thisTracer % localMin3D(iLayer,iCat,iCell)
                               write(stderrUnit,*) 'New value in cell =', thisTracer % array3D(iLayer,iCat,iCell)
                               write(stderrUnit,*) 'Tolerance, difference:', toleranceMin,  &
                                    thisTracer % localMin3D(iLayer,iCat,iCell) - thisTracer % array3D(iLayer,iCat,iCell)
                               call mpas_dmpar_global_abort('MPAS-seaice: IR advection, monotonicity violation (min, nDims == 3)')

                            elseif (thisTracer % array3D(iLayer,iCat,iCell) > &
                                    thisTracer % localMax3D(iLayer,iCat,iCell) + toleranceMax) then

                               write(stderrUnit,*) ' '
                               write(stderrUnit,*) 'IR advection, monotonicity violation!'
                               write(stderrUnit,*) 'Tracer, iLayer,iCat, iCell =', &
                                    trim(thisTracer % tracerName), iLayer,iCat, iCell
                               write(stderrUnit,*) 'Old max value in neighborhood =', thisTracer % localMax3D(iLayer,iCat,iCell)
                               write(stderrUnit,*) 'New value in cell =', thisTracer % array3D(iLayer,iCat,iCell)
                               write(stderrUnit,*) 'Tolerance, difference:', toleranceMax,  &
                                    thisTracer % array3D(iLayer,iCat,iCell) - thisTracer % localMax3D(iLayer,iCat,iCell)
                               call mpas_dmpar_global_abort('MPAS-seaice: IR advection, monotonicity violation (max, nDims == 3)')

                            endif

                         endif     ! tracer value has physical meaning

                      enddo   ! iLayer
                   enddo   ! iCell
                enddo   ! iCat

                block => block % next

             enddo     ! associated(block)

          endif        ! ndims
       endif           ! nparents

       thisTracer => thisTracer % next

    enddo      ! while(associated)

  end subroutine check_tracer_monotonicity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine zap_small_mass
!
!> \brief  set mass and tracers to zero in cells with very small mass
!> \author William Lipscomb
!> \date   August 2015
!> \details
!>  This routine zeroes out the mass and tracer values in cells with a
!>  very small mass.

!-----------------------------------------------------------------------

  subroutine zap_small_mass(&
       tracersHead,   &
       nCellsSolve,   &
       maskCategoryCell)

    type(tracer_type), pointer :: &
         tracersHead     !< Input/output: pointer to first element of linked list of tracers
                         ! The pointer stays attached to the first tracer

    integer, intent(in) ::  &
         nCellsSolve     !< Input: number of locally owned cells

    integer, dimension(:,:), intent(out) ::  &
         maskCategoryCell   !< Output: work array with dimensions (nCategories,nCells)
                            !  Here, it is a mask indicating which cells and categories are zeroed out

    ! local variables

    type(tracer_type), pointer :: thisTracer

    integer :: nCategories

    integer :: iCell, iCat

    ! Make the threshold very small (<< eps11) so as not to sweep potential problems under the rug

    real(kind=RKIND), parameter :: &
!!         smallMassThreshold = eps11
         smallMassThreshold = 1.0e-22_RKIND

    maskCategoryCell(:,:) = 0

    ! Starting with the mass-like field, loop over each cell and category.
    ! Zap those whose mass is below a small threshold.

    thisTracer => tracersHead   ! point to mass-like field

    if (thisTracer % ndims == 2) then

       nCategories = size(thisTracer % array2D, 1)

       do iCell = 1, nCellsSolve
          do iCat = 1, nCategories

             if (thisTracer % array2D(iCat,iCell) > 0.0_RKIND .and. &
                 thisTracer % array2D(iCat,iCell) < smallMassThreshold) then

!!                write(stderrUnit,*) 'Zap mass: iCell, iCat =', iCell, iCat

                ! zero out the mass and mask for this category and cell
                thisTracer % array2D(iCat,iCell) = 0.0_RKIND
                thisTracer % arrayMask2D(iCat,iCell) = 0

                ! for future reference, set a mask = 1 for this category and cell
                maskCategoryCell(iCat,iCell) = 1

             endif

          enddo  ! iCat
       enddo   ! iCell

    elseif (thisTracer % ndims == 3) then

       ! mass-like field must be 2D
       write(stderrUnit,*) 'ERROR: mass-like field must be 2D'
       call mpas_dmpar_global_abort('MPAS-seaice: IR: mass-like field must be 2D')

    endif   ! ndims

    ! Now loop through the tracers, setting tracer = 0 wherever iceCellCategoryMask = 1

    thisTracer => thisTracer % next

    do while (associated(thisTracer))

       if (thisTracer % ndims == 2) then

          nCategories = size(thisTracer % array2D, 1)

          do iCell = 1, nCellsSolve
             do iCat = 1, nCategories

                if (maskCategoryCell(iCat,iCell) == 1) then

                   ! zero out the tracer and mask
                   thisTracer % array2D(iCat,iCell) = 0.0_RKIND
                   thisTracer % arrayMask2D(iCat,iCell) = 0

!!                   write(stderrUnit,*) 'Zap tracer: name, iCell, iCat = ', trim(thisTracer % tracerName), iCell, iCat

                endif

             enddo   ! iCat
          enddo  ! iCell

       elseif (thisTracer % ndims == 3) then

          nCategories = size(thisTracer % array3D, 2)

          do iCell = 1, nCellsSolve
             do iCat = 1, nCategories

                if (maskCategoryCell(iCat,iCell) == 1) then

                   ! zero the tracer value and mask in each layer
                   thisTracer % array3D(:,iCat,iCell) = 0.0_RKIND
                   thisTracer % arrayMask3D(:,iCat,iCell) = 0

!!                   write(stderrUnit,*) 'Zap tracer: name, iCell, iCat = ', trim(thisTracer % tracerName), iCell, iCat

                endif

             enddo
          enddo

       endif   ! ndims

       thisTracer => thisTracer % next

    enddo  ! associated(thisTracer)

  end subroutine zap_small_mass

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine unit_vector_3d
!
!> \brief  normalize a 3D vector to have unit length
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  Given a 3D vector, return a vector with unit length
!
!-----------------------------------------------------------------------

  subroutine unit_vector_3d(vector)

    real(kind=RKIND), dimension(3), intent(inout) ::  &
         vector        !< Input/output: on input, any 3D Cartesian vector
                       !                on output, a vector of unit length

    ! local variables

    real(kind=RKIND) :: magnitude

    magnitude = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)

    if (magnitude > 0.0_RKIND) then
       vector(:) = vector(:) / magnitude
    else
       vector(:) = 0.0_RKIND
    endif

  end subroutine unit_vector_3d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine find_line_intersection
!
!> \brief  find the intersection point of two line segments in a plane
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  Given the Cartesian coordinates of the endpoints of 2 line segments in a plane,
!>  find the point where the two line segments intersect (if they intersect).
!>  If the lines extending these segments intersect (even if the segment themselves
!>  do not), this routine returns the vertices of their intersection point.
!-----------------------------------------------------------------------

  subroutine find_line_intersection(&
       p1, p2,             &
       p3, p4,             &
       lineIntersect,      &
       intersectionPoint,  &
       scalar)

    real(kind=RKIND), dimension(2), intent(in) :: &
         p1, p2,     & !< Input: x and y coordinates of the endpoints of the first line segment
         p3, p4        !< Input: x and y coordinates of the endpoints of the second line segment

    logical, intent(out) :: &
         lineIntersect !< Output: true if the lines intersect, else false
                       !    Note: If the lines share one or both endpoints, they do not intersect

    real(kind=RKIND), dimension(2), intent(out) :: &
         intersectionPoint !< Output:  x and y coordinates of the intersection point
                           !           returns (huge,huge) if the lines are parallel

    real(kind=RKIND), intent(out), optional :: &
         scalar        !< Output: a number describing where the intersection point lies on the second line
                       !          0 <= scalar <= 1 if the IP lies on the segment itself
                       !          scalar < 0 if the IP lies nearer p3 than p4
                       !          scalar > 1 if the IP lies nearer p4 than p3

    real(kind=RKIND) :: x1, x2, x3, x4, y1, y2, y3, y4
    real(kind=RKIND) :: rx, ry, sx, sy, t1, t2, rsCross, rsCrossMin

    ! initialize
    lineIntersect = .false.

    x1 = p1(1); y1 = p1(2)
    x2 = p2(1); y2 = p2(2)
    x3 = p3(1); y3 = p3(2)
    x4 = p4(1); y4 = p4(2)

    ! vector r = (x2,y2) - (x1,y1)
    rx = x2 - x1
    ry = y2 - y1

    ! vector s = (x4,y4) - (x3,y3)
    sx = x4 - x3
    sy = y4 - y3

    ! z component of cross product, r x s
    rsCross = rx*sy - ry*sx

    !TODO - Adjust this tolerance for single precision?
    ! Min value for cross product (if smaller than this, the lines are deemed to be parallel)
    rsCrossMin = eps11 * sqrt( (rx*rx + ry*ry) * (sx*sx + sy*sy) )

    ! Compute t1 and t2 such that (x1,y1) + t1*(rx,ry) = (x3,y3) + t2*(sx,sy)
    ! If both t1 and t2 are in the range [0,1], then the line segments intersect

    if (abs(rsCross) > rsCrossMin) then

       t1 = (sy*(x3-x1) - sx*(y3-y1)) / rsCross
       t2 = (ry*(x3-x1) - rx*(y3-y1)) / rsCross

       intersectionPoint(1) = x1 + t1*rx
       intersectionPoint(2) = y1 + t1*ry

       ! Note: The algorithm uses  > and < rather than >= and <=
       !       This means that if a vertex of one segment lies exactly on the other segment,
       !        the segments are deemed not to intersect.
       !       In this way we avoid creating triangles with degenerate vertices and zero area.
!!       if (t1 >= 0.0_RKIND .and. t1 <= 1.0_RKIND .and. t2 >= 0.0_RKIND .and. t2 <= 1.0_RKIND) then
       if (t1 > 0.0_RKIND .and. t1 < 1.0_RKIND .and. t2 > 0.0_RKIND .and. t2 < 1.0_RKIND) then
          lineIntersect = .true.
       endif

       if (present(scalar)) scalar = t2

    else   ! r x s <= rsCrossMin; the two lines are deemed to be parallel

       intersectionPoint(1) = huge(0.0_RKIND)
       intersectionPoint(2) = huge(0.0_RKIND)

       if (present(scalar)) scalar = huge(0.0_RKIND)

    endif

  end subroutine find_line_intersection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine dot_product_2d
!
!> \brief  compute the dot product of two 2D vectors
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  Given two 2D vectors, return the dot product
!
!-----------------------------------------------------------------------

  subroutine dot_product_2d(&
       vector1, vector2,    &
       dotProduct)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         vector1, vector2        !< Input: x/y/z coordinates of two 3D Cartesian vectorx

    real(kind=RKIND), intent(out) ::  &
         dotProduct              !< Output: dot product

    dotProduct = vector1(1)*vector2(1) + vector1(2)*vector2(2)

  end subroutine dot_product_2d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine dot_product_3d
!
!> \brief  compute the dot product of two 3D vectors
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  Given two 3D vectors, return the dot product
!
!-----------------------------------------------------------------------

  subroutine dot_product_3d(&
       vector1, vector2,      &
       dotProduct)

    real(kind=RKIND), dimension(3), intent(in) ::  &
         vector1, vector2        !< Input: x/y/z coordinates of two 3D Cartesian vectorx

    real(kind=RKIND), intent(out) ::  &
         dotProduct              !< Output: dot product

    dotProduct = vector1(1)*vector2(1) + vector1(2)*vector2(2) + vector1(3)*vector2(3)

  end subroutine dot_product_3d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cross_product_2d
!
!> \brief  compute the scalar value of the z component of the cross product
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  Given the components of two 2D vectors in the x-y plane, return the
!>  scalar equal to the z component of the cross product
!
!-----------------------------------------------------------------------

  subroutine cross_product_2d(&
       vector1, vector2,      &
       scalarOut)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         vector1, vector2        !< Input: x/y coordinates of two 2D Cartesian vectors

    real(kind=RKIND), intent(out) ::  &
         scalarOut               !< Output: z component of the cross product

    scalarOut = vector1(1)*vector2(2) - vector1(2)*vector2(1)

  end subroutine cross_product_2d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cross_product_3d
!
!> \brief  compute the cross product of two 3D vectors
!> \author William Lipscomb
!> \date   April 2015
!> \details
!>  Given two 3D vectors, return the cross product
!
!-----------------------------------------------------------------------

  subroutine cross_product_3d(&
       vector1, vector2,      &
       vectorOut)

    real(kind=RKIND), dimension(3), intent(in) ::  &
         vector1, vector2        !< Input: x/y/z coordinates of two 3D Cartesian vectorx

    real(kind=RKIND), dimension(3), intent(out) ::  &
         vectorOut               !< Output: x/y/z coordinates of the cross product

    vectorOut(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
    vectorOut(2) = vector1(3)*vector2(1) - vector1(1)*vector2(3)
    vectorOut(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)

  end subroutine cross_product_3d

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine quadrilateral_area
!
!> \brief  compute the area of a quadrilateral in a plane
!> \author William Lipscomb
!> \date   August 2015
!> \details
!>  Given four vertices of a quadrilateral, return its area.
!
!-----------------------------------------------------------------------

  subroutine quadrilateral_area(&
       vertex1, vertex2, vertex3, vertex4, &
       quadrilateralArea)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         vertex1, vertex2, vertex3, vertex4    !< Input: x/y coordinates of four quadrilateral vertices
                                               ! Note: These should be in order going around the perimeter
                                               !       (either clockwise or counterclockwise)
    real(kind=RKIND), intent(out) ::  &
         quadrilateralArea        !< Output: area of a quadrilateral

    real(kind=RKIND), dimension(2) :: triangleArea

    ! Compute quadrilateral area by splitting into 2 triangles

    call triangle_area(vertex1, vertex2, vertex3, triangleArea(1))

    call triangle_area(vertex1, vertex3, vertex4, triangleArea(2))

    quadrilateralArea = triangleArea(1) + triangleArea(2)

  end subroutine quadrilateral_area

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine triangle_area
!
!> \brief  compute the area of a triangle in a plane
!> \author William Lipscomb
!> \date   August 2015
!> \details
!>  Given three vertices of a triangle, return its area.
!
!-----------------------------------------------------------------------

  subroutine triangle_area(&
       vertex1, vertex2, vertex3, &
       triangleArea)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         vertex1, vertex2, vertex3    !< Input: x/y coordinates of three triangle vertices

    real(kind=RKIND), intent(out) ::  &
         triangleArea                 !< Output: area of a triangle

    triangleArea = abs (0.5_RKIND * ( (vertex2(1) - vertex1(1)) * (vertex3(2) - vertex1(2))  &
                                    - (vertex2(2) - vertex1(2)) * (vertex3(1) - vertex1(1)) ) )

  end subroutine triangle_area

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine point_in_half_plane
!
!> \brief  determine whether a given point is in the half-plane of a given line
!> \author William Lipscomb
!> \date   May 2015
!> \details
!>  Given the 2D Cartesian coordinates of two points on a line, along with a
!>  third point, determine whether the third point lies in the left half-plane
!>  of the line.
!
!-----------------------------------------------------------------------

  logical function point_in_half_plane(point, lineStart, lineEnd)

    real(kind=RKIND), dimension(2), intent(in) ::  &
         point          !< Input: x/y coordinates of an input point

    real(kind=RKIND), dimension(2), intent(in) ::  &
         lineStart, &   !< Input: x/y coordinates of starting point on a line segment
         lineEnd        !< Input: x/y coordinates of ending point on a line segment

    ! local variables

    real(kind=RKIND), dimension(2) :: &
         vector1, vector2   ! 2D work vectors in the x/y plane

    real(kind=RKIND) ::  &
         crossProduct       ! z component of the cross product

    ! vector from start point to end point of the input line
    vector1(:) = lineEnd(:) - lineStart(:)

    ! vector from start point of line to the input point
    vector2(:) = point(:) - lineStart(:)

    ! cross product
    call cross_product_2d(vector1, vector2, crossProduct)

    ! determine whether the point is in the half plane of the line
    if (crossProduct >= 0.0_RKIND) then
       point_in_half_plane = .true.
    else
       point_in_half_plane = .false.
    endif

  end function point_in_half_plane

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine volume_to_thickness
!
!> \brief  convert ice or snow volume to thickness
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  Given the area and volume of ice or snow, compute the thickness.
!
!-----------------------------------------------------------------------

  subroutine volume_to_thickness(&
       area,    &
       volumeToThickness)

    !Note: For the current version of CICE (as of Sept. 2015), ice and snow volume are in the tracer list.
    !      However, volume is not a tracer, but the product of area and thickness, where thickness is a tracer.
    !      This subroutine converts volume to thickness for purposes of IR.

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         area               !< Input: fractional area; first index is nCategories, second index is nCells

    real(kind=RKIND), dimension(:,:), intent(inout) ::  &
         volumeToThickness  !< Input/output: volume = area*thickness on input; thickness on output

    integer :: nCategories, nCells

    integer :: iCat, iCell

    nCategories = size(area,1)
    nCells = size(area,2)

    do iCell = 1, nCells
       do iCat = 1, nCategories

          if (area(iCat,iCell) > 0.0_RKIND) then
             volumeToThickness(iCat,iCell) = volumeToThickness(iCat,iCell) / area(iCat,iCell)
          else
             volumeToThickness(iCat,iCell) = 0.0_RKIND
          endif

       enddo
    enddo

  end subroutine volume_to_thickness

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine thickness_to_volume
!
!> \brief  convert ice or snow thickness to volume
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  Given the area and thickness of ice or snow, compute the volume.
!
!-----------------------------------------------------------------------

  subroutine thickness_to_volume(&
       area,       &
       thicknessToVolume)

    !Note: For the current version of CICE (as of Sept. 2015), ice and snow volume are in the tracer list.
    !      However, volume is not a tracer, but the product of area and thickness, where thickness is a tracer.
    !      This subroutine converts thickness back to volume for purposes of IR.

    real(kind=RKIND), dimension(:,:), intent(in) ::  &
         area             !< Input: fractional area; first index is nCategories, second index is nCells

    real(kind=RKIND), dimension(:,:), intent(inout) ::  &
         thicknessToVolume   !< Input/output: thickness on input; volume = area*thickness on output

    integer :: nCategories, nCells

    integer :: iCat, iCell

    nCategories = size(area,1)
    nCells = size(area,2)

    do iCell = 1, nCells
       do iCat = 1, nCategories

          thicknessToVolume(iCat,iCell) = area(iCat,iCell) * thicknessToVolume(iCat,iCell)

       enddo
    enddo

  end subroutine thickness_to_volume

end module cice_advection_incremental_remap
