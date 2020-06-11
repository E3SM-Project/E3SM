!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! \file mpas_ocn_mesh.F
!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_mesh
!
!>  \brief MPAS ocean mesh structure with GPU support
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> This module creates and maintains a primary ocean mesh structure
!> and ensures all mesh variables are copied to an accelerator device
!> if needed. Currently it consists of pointers to the existing MPAS mesh pool
!> variables, but is intended to eventually replace the mesh pool later.
!
!-------------------------------------------------------------------------------

module ocn_mesh

   ! module dependencies
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   implicit none
   private

   !----------------------------------------------------------------------------
   !
   ! Public parameters
   !
   !----------------------------------------------------------------------------
   !{{{

   integer, public :: &! mesh, array sizes
      nCellsAll, &! total number of local (owned+halo) cells in primary
      nEdgesAll, &! total number of local edge midpoints
      nVerticesAll, &! total number of local cells in dual (cell vertices)
      nCellsOwned, &! number of cells    owned by the local domain
      nEdgesOwned, &! number of edges    owned by the local domain
      nVerticesOwned, &! number of vertices owned by the local domain
      maxEdges, &! largest number of edges any polygon has
      maxEdges2, &! 2x the largest number of edges any polygon has
      vertexDegree, &! number of cells or edges touching each vertex
      nVertLevels, &! number of vertical levels
      nVertLevelsP1, &! number of vertical interfaces (levels plus one)
      nAdvectionCells  ! largest number of advection cells for any edge

   integer, public, dimension(:), allocatable :: &
      nCellsHalo, &! number of owned+halo(n) cells in local domain
      nEdgesHalo, &! number of owned+halo(n) edges in local domain
      nVerticesHalo      ! number of owned+halo(n) vertices in local domain

   integer, public, dimension(:), pointer :: &
      nEdgesOnEdge, &! number of edges connected to each edge point
      nEdgesOnCell, &! number of edges associated with each cell center
      maxLevelCell, &! max ocean level at bottom of cell
      maxLevelEdgeTop, &! max ocean level at top    of edge column
      maxLevelEdgeBot, &! max ocean level at bottom of edge column
      maxLevelVertexTop, &! max ocean level at top    of each vertex
      maxLevelVertexBot, &! max ocean level at bottom of each vertex
      nAdvCellsForEdge, &! number of cells contrib to advection at edge
      indexToCellID, &! global ID of each local cell
      indexToEdgeID, &! global ID of each local edge
      indexToVertexID    ! global ID of each local vertex

   integer, public, dimension(:, :), pointer :: &
      edgesOnEdge, &! index of edges connected to each edge
      cellsOnEdge, &! index of cells connected to each edge
      verticesOnEdge, &! index of vertices connected to each edge
      cellsOnCell, &! index of cells connected to each cell
      edgesOnCell, &! index of edges connected to each cell
      verticesOnCell, &! index of vertices connected to each cell
      cellsOnVertex, &! index of cells connected to each vertex
      edgesOnVertex, &! index of edges connected to each vertex
      kiteIndexOnCell, &! index of kite associated with each cell
      advCellsForEdge  ! index of cells contrib to advective flux at edge

   real(kind=RKIND), public, dimension(:), pointer :: &
      latCell, &! latitude  of cell centers
      lonCell, &! longitude of cell centers
      xCell, &! Cartesian x coord of cell center
      yCell, &! Cartesian y coord of cell center
      zCell, &! Cartesian z coord of cell center
      latEdge, &! latitude  of edge
      lonEdge, &! longitude of edge
      xEdge, &! Cartesian x coord of edge
      yEdge, &! Cartesian y coord of edge
      zEdge, &! Cartesian z coord of edge
      latVertex, &! latitude  of vertex
      lonVertex, &! longitude of vertex
      xVertex, &! Cartesian coord of vertex
      yVertex, &! Cartesian y coord of vertex
      zVertex, &! Cartesian z coord of vertex
      fEdge, &! Coriolus parameter at edge
      fVertex, &! Coriolus parameter at vertex
      fCell, &! Coriolus parameter at cell center
      dcEdge, &! length of edge = dist between cells across edge
      dvEdge, &! length of edge = dist between vertices along edge
      areaCell, &! area of each cell
      areaTriangle, &! area of each cell on dual grid
      bottomDepth, &! ocean bottom depth at each cell center
      refBottomDepth, &! ocean depth at bottom of cell for reference profile
      refBottomDepthTopOfCell, &! depth at top of cell for reference profile
      vertCoordMovementWeights, &! weights for distributing height perturb
      meshScalingDel2, &! mesh scaling factor for use in del2 diffusion
      meshScalingDel4, &! mesh scaling factor for use in del4 diffusion
      meshDensity, &! density of mesh
      angleEdge        ! angle the edge normal makes with local east

   ! Multiplicative masks and vectors for various conditions
   real(kind=RKIND), public, dimension(:, :), allocatable :: &
      edgeMask, &! mask to denote active edges    with depth
      cellMask, &! mask to denote active cells    with depth
      vertexMask, &! mask to denote active vertices with depth
      boundaryEdge, &! mask for boundary edges    at each level
      boundaryCell, &! mask for boundary cells    at each level
      boundaryVertex, &! mask for boundary vertices at each level
      edgeSignOnCell, &! sign of edge contributions to a cell
      edgeSignOnVertex, &! sign of edge contributions to a vertex
      highOrderAdvectionMask ! mask for high order advection contributions

   real(kind=RKIND), public, dimension(:, :), pointer :: &
      weightsOnEdge, &! weights on each edge
      kiteAreasOnVertex, &! real (vertexDegree nVertices)
      edgeTangentVectors, &! tangent unit vector at edge
      edgeNormalVectors, &! normal  unit vector at edge
      localVerticalUnitVectors, &! local unit vector iin vertical
      advCoefs, &! mesh-based advection coefficients
      advCoefs3rd         ! mesh-based advection coeffs for high order

   real(kind=RKIND), public, dimension(:, :, :), pointer :: &
      derivTwo, &! 2nd derivative of edge reconstruction polynomial
      cellTangentPlane, &! two vectors defining tangent plane at cell center
      coeffs_reconstruct  ! coeffs for reconstructing vectors at cell centers
   !}}}

   !----------------------------------------------------------------------------
   !
   ! Public member functions
   !
   !----------------------------------------------------------------------------
   !{{{

   public :: &
      ocn_meshCreate, &
      ocn_meshUpdateFields, &
      ocn_meshDestroy
   !}}}

!*******************************************************************************

contains

!*******************************************************************************
!
!  ocn_meshCreate
!
!> \brief Creates the primary ocean mesh data structure on both host and device
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> This module creates and maintains a primary ocean mesh structure
!> and ensures all mesh variables are copied to an accelerator device
!> if needed.
!
!-------------------------------------------------------------------------------

   subroutine ocn_meshCreate(domain) !{{{

      ! Input arguments

      type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

      ! Local variables

      integer :: &
         blockCount               ! counter for number of blocks

      type(block_type), pointer :: &
         block                    ! variables in current subblock

      type(mpas_pool_type), pointer :: &
         meshPool                 ! mesh variables in MPAS pool structure

      ! scalar pointers for retrieval, but convert to actual scalars in struct
      integer, pointer :: &! mesh dimensions
         nCellsTmp, &!
         nEdgesTmp, &!
         nVerticesTmp, &!
         maxEdgesTmp, &!
         maxEdges2Tmp, &!
         vertexDegreeTmp, &!
         nVertLevelsTmp, &!
         nVertLevelsP1Tmp, &!
         nAdvectionCellsTmp  !

      ! temporary pointers for converting index arrays
      integer, dimension(:), pointer :: &
         nCellsArrayTmp, &
         nEdgesArrayTmp, &
         nVerticesArrayTmp

      ! temporary pointers for converting masks
      integer i, k, n          ! loop indices
      integer, dimension(:, :), pointer :: &
         edgeMaskTmp, &
         vertexMaskTmp, &
         cellMaskTmp, &
         highOrderAdvectionMaskTmp, &
         edgeSignOnCellTmp, &
         edgeSignOnVertexTmp, &
         boundaryEdgeTmp, &
         boundaryVertexTmp, &
         boundaryCellTmp

      !***
      !*** end of preamble, begin code
      !***

      !$omp master

      blockCount = 0
      block => domain%blocklist
      do while (associated(block))

         !*** we now only support one block per task
         blockCount = blockCount + 1
         if (blockCount > 1) then
            call mpas_log_write( &
               'ocn_meshCreate: more than one block no longer supported', &
               MPAS_LOG_CRIT)
         endif

         ! retrieve the mpas mesh pool
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)

         ! set all mesh dimensions
         call mpas_pool_get_dimension(meshPool, 'nCells', &
                                      nCellsTmp)
         call mpas_pool_get_dimension(meshPool, 'nEdges', &
                                      nEdgesTmp)
         call mpas_pool_get_dimension(meshPool, 'nVertices', &
                                      nVerticesTmp)
         call mpas_pool_get_dimension(meshPool, 'maxEdges', &
                                      maxEdgesTmp)
         call mpas_pool_get_dimension(meshPool, 'maxEdges2', &
                                      maxEdges2Tmp)
         call mpas_pool_get_dimension(meshPool, 'vertexDegree', &
                                      vertexDegreeTmp)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', &
                                      nVertLevelsTmp)
         call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', &
                                      nVertLevelsP1Tmp)
         call mpas_pool_get_dimension(meshPool, 'nAdvectionCells', &
                                      nAdvectionCellsTmp)
         call mpas_pool_get_dimension(meshPool, 'nCellsArray', &
                                      nCellsArrayTmp)
         call mpas_pool_get_dimension(meshPool, 'nEdgesArray', &
                                      nEdgesArrayTmp)
         call mpas_pool_get_dimension(meshPool, 'nVerticesArray', &
                                      nVerticesArrayTmp)

         ! translate scalar pointers to scalars in new mesh structure
         maxEdges = maxEdgesTmp
         maxEdges2 = maxEdges2Tmp
         vertexDegree = vertexDegreeTmp
         nVertLevels = nVertLevelsTmp
         nVertLevelsP1 = nVertLevelsP1Tmp
         nAdvectionCells = nAdvectionCellsTmp

         ! convert previous index limits into new arrangement
         nCellsAll = nCellsTmp
         nEdgesAll = nEdgesTmp
         nVerticesAll = nVerticesTmp

         n = size(nCellsArrayTmp)
         allocate (nCellsHalo(n - 1))
         nCellsOwned = nCellsArrayTmp(1)
         do i = 2, n
            nCellsHalo(i - 1) = nCellsArrayTmp(i)
         end do

         n = size(nEdgesArrayTmp)
         allocate (nEdgesHalo(n - 1))
         nEdgesOwned = nEdgesArrayTmp(1)
         do i = 2, n
            nEdgesHalo(i - 1) = nEdgesArrayTmp(i)
         end do

         n = size(nVerticesArrayTmp)
         allocate (nVerticesHalo(n - 1))
         nVerticesOwned = nVerticesArrayTmp(1)
         do i = 2, n
            nVerticesHalo(i - 1) = nVerticesArrayTmp(i)
         end do

         ! set a lot of connectivity info
         call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', &
                                  nEdgesOnEdge)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', &
                                  nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', &
                                  maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', &
                                  maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', &
                                  maxLevelEdgeBot)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', &
                                  maxLevelVertexTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', &
                                  maxLevelVertexBot)
         call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge', &
                                  nAdvCellsForEdge)
         call mpas_pool_get_array(meshPool, 'indexToCellID', &
                                  indexToCellID)
         call mpas_pool_get_array(meshPool, 'indexToEdgeID', &
                                  indexToEdgeID)
         call mpas_pool_get_array(meshPool, 'indexToVertexID', &
                                  indexToVertexID)
         call mpas_pool_get_array(meshPool, 'edgesOnEdge', &
                                  edgesOnEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', &
                                  cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'verticesOnEdge', &
                                  verticesOnEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', &
                                  cellsOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', &
                                  edgesOnCell)
         call mpas_pool_get_array(meshPool, 'verticesOnCell', &
                                  verticesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnVertex', &
                                  cellsOnVertex)
         call mpas_pool_get_array(meshPool, 'edgesOnVertex', &
                                  edgesOnVertex)
         call mpas_pool_get_array(meshPool, 'kiteIndexOnCell', &
                                  kiteIndexOnCell)
         call mpas_pool_get_array(meshPool, 'advCellsForEdge', &
                                  advCellsForEdge)

         ! now set a number of physics and numerical properties of mesh
         call mpas_pool_get_array(meshPool, 'latCell', &
                                  latCell)
         call mpas_pool_get_array(meshPool, 'lonCell', &
                                  lonCell)
         call mpas_pool_get_array(meshPool, 'xCell', &
                                  xCell)
         call mpas_pool_get_array(meshPool, 'yCell', &
                                  yCell)
         call mpas_pool_get_array(meshPool, 'zCell', &
                                  zCell)
         call mpas_pool_get_array(meshPool, 'latEdge', &
                                  latEdge)
         call mpas_pool_get_array(meshPool, 'lonEdge', &
                                  lonEdge)
         call mpas_pool_get_array(meshPool, 'xEdge', &
                                  xEdge)
         call mpas_pool_get_array(meshPool, 'yEdge', &
                                  yEdge)
         call mpas_pool_get_array(meshPool, 'zEdge', &
                                  zEdge)
         call mpas_pool_get_array(meshPool, 'latVertex', &
                                  latVertex)
         call mpas_pool_get_array(meshPool, 'lonVertex', &
                                  lonVertex)
         call mpas_pool_get_array(meshPool, 'xVertex', &
                                  xVertex)
         call mpas_pool_get_array(meshPool, 'yVertex', &
                                  yVertex)
         call mpas_pool_get_array(meshPool, 'zVertex', &
                                  zVertex)
         call mpas_pool_get_array(meshPool, 'fEdge', &
                                  fEdge)
         call mpas_pool_get_array(meshPool, 'fVertex', &
                                  fVertex)
         call mpas_pool_get_array(meshPool, 'fCell', &
                                  fCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', &
                                  dcEdge)
         call mpas_pool_get_array(meshPool, 'dvEdge', &
                                  dvEdge)
         call mpas_pool_get_array(meshPool, 'areaCell', &
                                  areaCell)
         call mpas_pool_get_array(meshPool, 'areaTriangle', &
                                  areaTriangle)
         call mpas_pool_get_array(meshPool, 'weightsOnEdge', &
                                  weightsOnEdge)
         call mpas_pool_get_array(meshPool, 'bottomDepth', &
                                  bottomDepth)
         call mpas_pool_get_array(meshPool, 'refBottomDepth', &
                                  refBottomDepth)
         call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', &
                                  refBottomDepthTopOfCell)
         call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', &
                                  vertCoordMovementWeights)
         call mpas_pool_get_array(meshPool, 'meshScalingDel2', &
                                  meshScalingDel2)
         call mpas_pool_get_array(meshPool, 'meshScalingDel4', &
                                  meshScalingDel4)
         call mpas_pool_get_array(meshPool, 'meshDensity', &
                                  meshDensity)
         call mpas_pool_get_array(meshPool, 'angleEdge', &
                                  angleEdge)

         call mpas_pool_get_array(meshPool, 'weightsOnEdge', &
                                  weightsOnEdge)
         call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', &
                                  kiteAreasOnVertex)
         call mpas_pool_get_array(meshPool, 'edgeTangentVectors', &
                                  edgeTangentVectors)
         call mpas_pool_get_array(meshPool, 'edgeNormalVectors', &
                                  edgeNormalVectors)
         call mpas_pool_get_array(meshPool, 'localVerticalUnitVectors', &
                                  localVerticalUnitVectors)
         call mpas_pool_get_array(meshPool, 'advCoefs', &
                                  advCoefs)
         call mpas_pool_get_array(meshPool, 'advCoefs3rd', &
                                  advCoefs3rd)
         call mpas_pool_get_array(meshPool, 'derivTwo', &
                                  derivTwo)
         call mpas_pool_get_array(meshPool, 'cellTangentPlane', &
                                  cellTangentPlane)
         call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', &
                                  coeffs_reconstruct)

         ! For masks, we wish to convert to real multiplicative masks
         ! so retrieve integer version and set real version in struct.
         ! Once these are converted in Registry, we can eliminate conversion.
         call mpas_pool_get_array(meshPool, 'edgeMask', &
                                  edgeMaskTmp)
         call mpas_pool_get_array(meshPool, 'vertexMask', &
                                  vertexMaskTmp)
         call mpas_pool_get_array(meshPool, 'cellMask', &
                                  cellMaskTmp)
         call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', &
                                  highOrderAdvectionMaskTmp)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell', &
                                  edgeSignOnCellTmp)
         call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', &
                                  edgeSignOnVertexTmp)
         call mpas_pool_get_array(meshPool, 'boundaryEdge', &
                                  boundaryEdgeTmp)
         call mpas_pool_get_array(meshPool, 'boundaryVertex', &
                                  boundaryVertexTmp)
         call mpas_pool_get_array(meshPool, 'boundaryCell', &
                                  boundaryCellTmp)

         allocate ( &
            edgeMask(nVertLevels, nEdgesAll+1), &
            cellMask(nVertLevels, nCellsAll+1), &
            vertexMask(nVertLevels, nVerticesAll+1), &
            boundaryEdge(nVertLevels, nEdgesAll+1), &
            boundaryCell(nVertLevels, nCellsAll+1), &
            boundaryVertex(nVertLevels, nVerticesAll+1), &
            edgeSignOnCell(maxEdges, nCellsAll+1), &
            edgeSignOnVertex(maxEdges, nVerticesAll+1), &
            highOrderAdvectionMask(nVertLevels, nEdgesAll+1))

         do n = 1, nCellsAll+1
         do k = 1, nVertLevels
            cellMask(k, n) = real(cellMaskTmp(k, n), RKIND)
            boundaryCell(k, n) = real(boundaryCellTmp(k, n), RKIND)
         end do
         end do

         do n = 1, nCellsAll+1
         do k = 1, maxEdges
            edgeSignOnCell(k, n) = real(edgeSignOnCellTmp(k, n), RKIND)
         end do
         end do

         do n = 1, nEdgesAll+1
         do k = 1, nVertLevels
            edgeMask(k, n) = real(edgeMaskTmp(k, n), RKIND)
            boundaryEdge(k, n) = real(boundaryEdgeTmp(k, n), RKIND)
            highOrderAdvectionMask(k, n) = &
               real(highOrderAdvectionMaskTmp(k, n), RKIND)
         end do
         end do

         do n = 1, nVerticesAll+1
         do k = 1, nVertLevels
            vertexMask(k, n) = real(vertexMaskTmp(k, n), RKIND)
            boundaryVertex(k, n) = real(boundaryVertexTmp(k, n), RKIND)
         end do
         end do

         do n = 1, nVerticesAll+1
         do k = 1, vertexDegree
            edgeSignOnVertex(k, n) = real(edgeSignOnVertexTmp(k, n), RKIND)
         end do
         end do

         ! Copy structure and components to accelerator
         ! Must manually do a deep copy with structure itself copied first.
         !$acc enter data copyin( &
         !$acc                   nCellsAll,                &
         !$acc                   nEdgesAll,                &
         !$acc                   nVerticesAll,             &
         !$acc                   nCellsOwned,              &
         !$acc                   nEdgesOwned,              &
         !$acc                   nVerticesOwned,           &
         !$acc                   maxEdges,                 &
         !$acc                   maxEdges2,                &
         !$acc                   vertexDegree,             &
         !$acc                   nVertLevels,              &
         !$acc                   nVertLevelsP1,            &
         !$acc                   nAdvectionCells,          &
         !$acc                   nEdgesHalo,               &
         !$acc                   nCellsHalo,               &
         !$acc                   nVerticesHalo,            &
         !$acc                   nEdgesOnEdge,             &
         !$acc                   nEdgesOnCell,             &
         !$acc                   maxLevelCell,             &
         !$acc                   maxLevelEdgeTop,          &
         !$acc                   maxLevelEdgeBot,          &
         !$acc                   maxLevelVertexTop,        &
         !$acc                   maxLevelVertexBot,        &
         !$acc                   nAdvCellsForEdge,         &
         !$acc                   indexToCellID,            &
         !$acc                   indexToEdgeID,            &
         !$acc                   indexToVertexID,          &
         !$acc                   edgesOnEdge,              &
         !$acc                   cellsOnEdge,              &
         !$acc                   verticesOnEdge,           &
         !$acc                   cellsOnCell,              &
         !$acc                   edgesOnCell,              &
         !$acc                   verticesOnCell,           &
         !$acc                   cellsOnVertex,            &
         !$acc                   edgesOnVertex,            &
         !$acc                   kiteIndexOnCell,          &
         !$acc                   advCellsForEdge,          &
         !$acc                   latCell,                  &
         !$acc                   lonCell,                  &
         !$acc                   xCell,                    &
         !$acc                   yCell,                    &
         !$acc                   zCell,                    &
         !$acc                   latEdge,                  &
         !$acc                   lonEdge,                  &
         !$acc                   xEdge,                    &
         !$acc                   yEdge,                    &
         !$acc                   zEdge,                    &
         !$acc                   latVertex,                &
         !$acc                   lonVertex,                &
         !$acc                   xVertex,                  &
         !$acc                   yVertex,                  &
         !$acc                   zVertex,                  &
         !$acc                   fEdge,                    &
         !$acc                   fVertex,                  &
         !$acc                   fCell,                    &
         !$acc                   dcEdge,                   &
         !$acc                   dvEdge,                   &
         !$acc                   areaCell,                 &
         !$acc                   areaTriangle,             &
         !$acc                   bottomDepth,              &
         !$acc                   refBottomDepth,           &
         !$acc                   refBottomDepthTopOfCell,  &
         !$acc                   vertCoordMovementWeights, &
         !$acc                   meshScalingDel2,          &
         !$acc                   meshScalingDel4,          &
         !$acc                   meshDensity,              &
         !$acc                   angleEdge,                &
         !$acc                   edgeMask,                 &
         !$acc                   cellMask,                 &
         !$acc                   vertexMask,               &
         !$acc                   boundaryEdge,             &
         !$acc                   boundaryCell,             &
         !$acc                   boundaryVertex,           &
         !$acc                   edgeSignOnCell,           &
         !$acc                   edgeSignOnVertex,         &
         !$acc                   highOrderAdvectionMask,   &
         !$acc                   weightsOnEdge,            &
         !$acc                   kiteAreasOnVertex,        &
         !$acc                   edgeTangentVectors,       &
         !$acc                   edgeNormalVectors,        &
         !$acc                   localVerticalUnitVectors, &
         !$acc                   advCoefs,                 &
         !$acc                   advCoefs3rd,              &
         !$acc                   derivTwo,                 &
         !$acc                   cellTangentPlane,         &
         !$acc                   coeffs_reconstruct)

         block => block%next
      end do

      !$omp end master
      !$omp barrier

!-------------------------------------------------------------------------------

   end subroutine ocn_meshCreate !}}}

!*******************************************************************************
!
!  ocn_meshDestroy
!
!> \brief Destroy mesh structure and removes from device
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> This module removes the mesh variables from the device and invalidates
!> all pointers in the mesh structure.
!
!-------------------------------------------------------------------------------

   subroutine ocn_meshDestroy(err) !{{{

      ! Input variables

      ! Since the ocnMesh is currently a public module variable, no inputs
      ! here, but eventually may want to treat ocnMesh as a specific
      ! instantiation instead and pass via args everywhere. If so, need an
      ! input mesh here

      ! Output variables

      integer, intent(out) :: &
         err                   ! returned error flag

      ! Local variables

      !***
      !*** end of preamble, begin code
      !***

      err = 0

      !$omp barrier
      !$omp master

      ! First remove data from the device. Must remove components first,
      ! then remove the mesh type itself.

      !$acc exit data delete(nCellsAll,         &
      !$acc                  nEdgesAll,         &
      !$acc                  nVerticesAll,      &
      !$acc                  nCellsOwned,       &
      !$acc                  nEdgesOwned,       &
      !$acc                  nVerticesOwned,    &
      !$acc                  maxEdges,          &
      !$acc                  maxEdges2,         &
      !$acc                  vertexDegree,      &
      !$acc                  nVertLevels,       &
      !$acc                  nVertLevelsP1,     &
      !$acc                  nAdvectionCells,   &
      !$acc                  nEdgesHalo,        &
      !$acc                  nCellsHalo,        &
      !$acc                  nVerticesHalo,     &
      !$acc                  nEdgesOnEdge,      &
      !$acc                  nEdgesOnCell,      &
      !$acc                  maxLevelCell,      &
      !$acc                  maxLevelEdgeTop,   &
      !$acc                  maxLevelEdgeBot,   &
      !$acc                  maxLevelVertexTop, &
      !$acc                  maxLevelVertexBot, &
      !$acc                  nAdvCellsForEdge,  &
      !$acc                  indexToCellID,     &
      !$acc                  indexToEdgeID,     &
      !$acc                  indexToVertexID,   &
      !$acc                  edgesOnEdge,       &
      !$acc                  cellsOnEdge,       &
      !$acc                  verticesOnEdge,    &
      !$acc                  cellsOnCell,       &
      !$acc                  edgesOnCell,       &
      !$acc                  verticesOnCell,    &
      !$acc                  cellsOnVertex,     &
      !$acc                  edgesOnVertex,     &
      !$acc                  kiteIndexOnCell,   &
      !$acc                  advCellsForEdge,   &
      !$acc                  latCell,           &
      !$acc                  lonCell,           &
      !$acc                  xCell,             &
      !$acc                  yCell,             &
      !$acc                  zCell,             &
      !$acc                  latEdge,           &
      !$acc                  lonEdge,           &
      !$acc                  xEdge,             &
      !$acc                  yEdge,             &
      !$acc                  zEdge,             &
      !$acc                  latVertex,         &
      !$acc                  lonVertex,         &
      !$acc                  xVertex,           &
      !$acc                  yVertex,           &
      !$acc                  zVertex,           &
      !$acc                  fEdge,             &
      !$acc                  fVertex,           &
      !$acc                  fCell,             &
      !$acc                  dcEdge,            &
      !$acc                  dvEdge,            &
      !$acc                  areaCell,          &
      !$acc                  areaTriangle,      &
      !$acc                  bottomDepth,       &
      !$acc                  refBottomDepth,    &
      !$acc                  refBottomDepthTopOfCell, &
      !$acc                  vertCoordMovementWeights, &
      !$acc                  meshScalingDel2,   &
      !$acc                  meshScalingDel4,   &
      !$acc                  meshDensity,       &
      !$acc                  angleEdge,         &
      !$acc                  edgeMask,          &
      !$acc                  cellMask,          &
      !$acc                  vertexMask,        &
      !$acc                  boundaryEdge,      &
      !$acc                  boundaryCell,      &
      !$acc                  boundaryVertex,    &
      !$acc                  edgeSignOnCell,    &
      !$acc                  edgeSignOnVertex,  &
      !$acc                  highOrderAdvectionMask,&
      !$acc                  weightsOnEdge,     &
      !$acc                  kiteAreasOnVertex, &
      !$acc                  edgeTangentVectors,&
      !$acc                  edgeNormalVectors, &
      !$acc                  localVerticalUnitVectors, &
      !$acc                  advCoefs,          &
      !$acc                  advCoefs3rd,       &
      !$acc                  derivTwo,          &
      !$acc                  cellTangentPlane,  &
      !$acc                  coeffs_reconstruct)

      ! Reset all scalars to zero
      nCellsAll = 0
      nEdgesAll = 0
      nVerticesAll = 0
      nCellsOwned = 0
      nEdgesOwned = 0
      nVerticesOwned = 0
      maxEdges = 0
      maxEdges2 = 0
      vertexDegree = 0
      nVertLevels = 0
      nVertLevelsP1 = 0
      nAdvectionCells = 0

      ! Now nullify all pointers to invalidate fields
      ! If this becomes the only mesh structure and mesh pool is eliminated,
      !  then we will want to deallocate here instead of nullify.

      nullify (nEdgesOnEdge, &
               nEdgesOnCell, &
               maxLevelCell, &
               maxLevelEdgeTop, &
               maxLevelEdgeBot, &
               maxLevelVertexTop, &
               maxLevelVertexBot, &
               nAdvCellsForEdge, &
               indexToCellID, &
               indexToEdgeID, &
               indexToVertexID, &
               edgesOnEdge, &
               cellsOnEdge, &
               verticesOnEdge, &
               cellsOnCell, &
               edgesOnCell, &
               verticesOnCell, &
               cellsOnVertex, &
               edgesOnVertex, &
               kiteIndexOnCell, &
               advCellsForEdge, &
               latCell, &
               lonCell, &
               xCell, &
               yCell, &
               zCell, &
               latEdge, &
               lonEdge, &
               xEdge, &
               yEdge, &
               zEdge, &
               latVertex, &
               lonVertex, &
               xVertex, &
               yVertex, &
               zVertex, &
               fEdge, &
               fVertex, &
               fCell, &
               dcEdge, &
               dvEdge, &
               areaCell, &
               areaTriangle, &
               bottomDepth, &
               refBottomDepth, &
               refBottomDepthTopOfCell, &
               vertCoordMovementWeights, &
               meshScalingDel2, &
               meshScalingDel4, &
               meshDensity, &
               angleEdge, &
               weightsOnEdge, &
               kiteAreasOnVertex, &
               edgeTangentVectors, &
               edgeNormalVectors, &
               localVerticalUnitVectors, &
               advCoefs, &
               advCoefs3rd, &
               derivTwo, &
               cellTangentPlane, &
               coeffs_reconstruct)

      deallocate (nEdgesHalo, &
                  nCellsHalo, &
                  nVerticesHalo, &
                  edgeMask, &
                  cellMask, &
                  vertexMask, &
                  boundaryEdge, &
                  boundaryCell, &
                  boundaryVertex, &
                  edgeSignOnCell, &
                  edgeSignOnVertex, &
                  highOrderAdvectionMask)

      !$omp end master

!-------------------------------------------------------------------------------

   end subroutine ocn_meshDestroy !}}}

!*******************************************************************************
!
!  ocn_meshUpdateFields
!
!> \brief Updates fields on an accelerator device
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> Many mesh fields are computed or input later in the initialization
!> phase after the meshCreate call. This routine updates these fields
!> on the device. The routine is not needed if the original mesh pool is
!> eliminated and this mesh becomes the only mesh structure and can be
!> updated directly.
!
!-------------------------------------------------------------------------------

   subroutine ocn_meshUpdateFields(domain) !{{{

      ! Input arguments

      type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

      ! Local variables

      type(block_type), pointer :: &
         block                    ! variables in current subblock

      type(mpas_pool_type), pointer :: &
         meshPool                 ! mesh variables in MPAS pool structure

      ! temporary pointers for converting masks
      integer i, k, n          ! loop indices
      integer, dimension(:, :), pointer :: &
         edgeMaskTmp, &
         vertexMaskTmp, &
         cellMaskTmp, &
         highOrderAdvectionMaskTmp, &
         edgeSignOnCellTmp, &
         edgeSignOnVertexTmp, &
         boundaryEdgeTmp, &
         boundaryVertexTmp, &
         boundaryCellTmp

      !***
      !*** end of preamble, begin code
      !***

      !$omp master

      ! already check during create that there should only be one block
      block => domain%blocklist

      ! retrieve the mpas mesh pool
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)

      ! Mesh dimensions should not have changed so no updates

      ! Because these fields are pointers into meshPool, they do not need
      ! to be updated - they capture the updates automatically.
      !call mpas_pool_get_array(meshPool, 'nEdgesOnEdge',      &
      !                                    nEdgesOnEdge)
      !call mpas_pool_get_array(meshPool, 'nEdgesOnCell',      &
      !                                    nEdgesOnCell)
      !call mpas_pool_get_array(meshPool, 'maxLevelCell',      &
      !                                    maxLevelCell)
      !call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',   &
      !                                    maxLevelEdgeTop)
      !call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot',   &
      !                                    maxLevelEdgeBot)
      !call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', &
      !                                    maxLevelVertexTop)
      !call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', &
      !                                    maxLevelVertexBot)
      !call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge',  &
      !                                    nAdvCellsForEdge)
      !call mpas_pool_get_array(meshPool, 'indexToCellID',     &
      !                                    indexToCellID)
      !call mpas_pool_get_array(meshPool, 'indexToEdgeID',     &
      !                                    indexToEdgeID)
      !call mpas_pool_get_array(meshPool, 'indexToVertexID',   &
      !                                    indexToVertexID)
      !call mpas_pool_get_array(meshPool, 'edgesOnEdge',       &
      !                                    edgesOnEdge)
      !call mpas_pool_get_array(meshPool, 'cellsOnEdge',       &
      !                                    cellsOnEdge)
      !call mpas_pool_get_array(meshPool, 'verticesOnEdge',    &
      !                                    verticesOnEdge)
      !call mpas_pool_get_array(meshPool, 'cellsOnCell',       &
      !                                    cellsOnCell)
      !call mpas_pool_get_array(meshPool, 'edgesOnCell',       &
      !                                    edgesOnCell)
      !call mpas_pool_get_array(meshPool, 'verticesOnCell',    &
      !                                    verticesOnCell)
      !call mpas_pool_get_array(meshPool, 'cellsOnVertex',     &
      !                                    cellsOnVertex)
      !call mpas_pool_get_array(meshPool, 'edgesOnVertex',     &
      !                                    edgesOnVertex)
      !call mpas_pool_get_array(meshPool, 'kiteIndexOnCell',   &
      !                                    kiteIndexOnCell)
      !call mpas_pool_get_array(meshPool, 'advCellsForEdge',   &
      !                                    advCellsForEdge)

      ! these are also pointers that do not require updating
      ! now set a number of physics and numerical properties of mesh
      !call mpas_pool_get_array(meshPool, 'latCell',           &
      !                                    latCell)
      !call mpas_pool_get_array(meshPool, 'lonCell',           &
      !                                    lonCell)
      !call mpas_pool_get_array(meshPool, 'xCell',             &
      !                                    xCell)
      !call mpas_pool_get_array(meshPool, 'yCell',             &
      !                                    yCell)
      !call mpas_pool_get_array(meshPool, 'zCell',             &
      !                                    zCell)
      !call mpas_pool_get_array(meshPool, 'latEdge',           &
      !                                    latEdge)
      !call mpas_pool_get_array(meshPool, 'lonEdge',           &
      !                                    lonEdge)
      !call mpas_pool_get_array(meshPool, 'xEdge',             &
      !                                    xEdge)
      !call mpas_pool_get_array(meshPool, 'yEdge',             &
      !                                    yEdge)
      !call mpas_pool_get_array(meshPool, 'zEdge',             &
      !                                    zEdge)
      !call mpas_pool_get_array(meshPool, 'latVertex',         &
      !                                    latVertex)
      !call mpas_pool_get_array(meshPool, 'lonVertex',         &
      !                                    lonVertex)
      !call mpas_pool_get_array(meshPool, 'xVertex',           &
      !                                    xVertex)
      !call mpas_pool_get_array(meshPool, 'yVertex',           &
      !                                    yVertex)
      !call mpas_pool_get_array(meshPool, 'zVertex',           &
      !                                    zVertex)
      !call mpas_pool_get_array(meshPool, 'fEdge',             &
      !                                    fEdge)
      !call mpas_pool_get_array(meshPool, 'fVertex',           &
      !                                    fVertex)
      !call mpas_pool_get_array(meshPool, 'fCell',             &
      !                                    fCell)
      !call mpas_pool_get_array(meshPool, 'dcEdge',            &
      !                                    dcEdge)
      !call mpas_pool_get_array(meshPool, 'dvEdge',            &
      !                                    dvEdge)
      !call mpas_pool_get_array(meshPool, 'areaCell',          &
      !                                    areaCell)
      !call mpas_pool_get_array(meshPool, 'areaTriangle',      &
      !                                    areaTriangle)
      !call mpas_pool_get_array(meshPool, 'weightsOnEdge',     &
      !                                    weightsOnEdge)
      !call mpas_pool_get_array(meshPool, 'bottomDepth',       &
      !                                    bottomDepth)
      !call mpas_pool_get_array(meshPool, 'refBottomDepth',    &
      !                                    refBottomDepth)
      !call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell',   &
      !                                    refBottomDepthTopOfCell)
      !call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights',  &
      !                                    vertCoordMovementWeights)
      !call mpas_pool_get_array(meshPool, 'meshScalingDel2',   &
      !                                    meshScalingDel2)
      !call mpas_pool_get_array(meshPool, 'meshScalingDel4',   &
      !                                    meshScalingDel4)
      !call mpas_pool_get_array(meshPool, 'meshDensity',       &
      !                                    meshDensity)
      !call mpas_pool_get_array(meshPool, 'angleEdge',         &
      !                                    angleEdge)
      !
      !call mpas_pool_get_array(meshPool, 'weightsOnEdge',             &
      !                                    weightsOnEdge)
      !call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex',         &
      !                                    kiteAreasOnVertex)
      !call mpas_pool_get_array(meshPool, 'edgeTangentVectors',        &
      !                                    edgeTangentVectors)
      !call mpas_pool_get_array(meshPool, 'edgeNormalVectors',         &
      !                                    edgeNormalVectors)
      !call mpas_pool_get_array(meshPool, 'localVerticalUnitVectors',  &
      !                                    localVerticalUnitVectors)
      !call mpas_pool_get_array(meshPool, 'advCoefs',                  &
      !                                    advCoefs)
      !call mpas_pool_get_array(meshPool, 'advCoefs3rd',               &
      !                                    advCoefs3rd)
      !call mpas_pool_get_array(meshPool, 'derivTwo',                  &
      !                                    derivTwo)
      !call mpas_pool_get_array(meshPool, 'cellTangentPlane',          &
      !                                    cellTangentPlane)
      !call mpas_pool_get_array(meshPool, 'coeffs_reconstruct',        &
      !                                    coeffs_reconstruct)

      ! For masks, we converted to real masks, so need to recompute for
      ! updated values. Probably only the edgeSign masks need updating, but
      ! do them all to be sure.
      call mpas_pool_get_array(meshPool, 'edgeMask', &
                               edgeMaskTmp)
      call mpas_pool_get_array(meshPool, 'vertexMask', &
                               vertexMaskTmp)
      call mpas_pool_get_array(meshPool, 'cellMask', &
                               cellMaskTmp)
      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', &
                               highOrderAdvectionMaskTmp)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', &
                               edgeSignOnCellTmp)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', &
                               edgeSignOnVertexTmp)
      call mpas_pool_get_array(meshPool, 'boundaryEdge', &
                               boundaryEdgeTmp)
      call mpas_pool_get_array(meshPool, 'boundaryVertex', &
                               boundaryVertexTmp)
      call mpas_pool_get_array(meshPool, 'boundaryCell', &
                               boundaryCellTmp)

      do n = 1, nCellsAll+1
      do k = 1, nVertLevels
         cellMask(k, n) = real(cellMaskTmp(k, n), RKIND)
         boundaryCell(k, n) = real(boundaryCellTmp(k, n), RKIND)
      end do
      end do

      do n = 1, nCellsAll+1
      do k = 1, maxEdges
         edgeSignOnCell(k, n) = real(edgeSignOnCellTmp(k, n), RKIND)
      end do
      end do

      do n = 1, nEdgesAll+1
      do k = 1, nVertLevels
         edgeMask(k, n) = real(edgeMaskTmp(k, n), RKIND)
         boundaryEdge(k, n) = real(boundaryEdgeTmp(k, n), RKIND)
         highOrderAdvectionMask(k, n) = &
            real(highOrderAdvectionMaskTmp(k, n), RKIND)
      end do
      end do

      do n = 1, nVerticesAll+1
      do k = 1, nVertLevels
         vertexMask(k, n) = real(vertexMaskTmp(k, n), RKIND)
         boundaryVertex(k, n) = real(boundaryVertexTmp(k, n), RKIND)
      end do
      end do

      do n = 1, nVerticesAll+1
      do k = 1, vertexDegree
         edgeSignOnVertex(k, n) = &
            real(edgeSignOnVertexTmp(k, n), RKIND)
      end do
      end do

      ! Go ahead and update all fields on device to be safe.
      ! NOTE: if we end up computing some fields on the device during
      !       init, the update must go the opposite direction (update host)

      !$acc update device(nCellsAll,                &
      !$acc               nEdgesAll,                &
      !$acc               nVerticesAll,             &
      !$acc               nCellsOwned,              &
      !$acc               nEdgesOwned,              &
      !$acc               nVerticesOwned,           &
      !$acc               maxEdges,                 &
      !$acc               maxEdges2,                &
      !$acc               vertexDegree,             &
      !$acc               nVertLevels,              &
      !$acc               nVertLevelsP1,            &
      !$acc               nAdvectionCells,          &
      !$acc               nEdgesHalo,               &
      !$acc               nCellsHalo,               &
      !$acc               nVerticesHalo,            &
      !$acc               nEdgesOnEdge,             &
      !$acc               nEdgesOnCell,             &
      !$acc               maxLevelCell,             &
      !$acc               maxLevelEdgeTop,          &
      !$acc               maxLevelEdgeBot,          &
      !$acc               maxLevelVertexTop,        &
      !$acc               maxLevelVertexBot,        &
      !$acc               nAdvCellsForEdge,         &
      !$acc               indexToCellID,            &
      !$acc               indexToEdgeID,            &
      !$acc               indexToVertexID,          &
      !$acc               edgesOnEdge,              &
      !$acc               cellsOnEdge,              &
      !$acc               verticesOnEdge,           &
      !$acc               cellsOnCell,              &
      !$acc               edgesOnCell,              &
      !$acc               verticesOnCell,           &
      !$acc               cellsOnVertex,            &
      !$acc               edgesOnVertex,            &
      !$acc               kiteIndexOnCell,          &
      !$acc               advCellsForEdge,          &
      !$acc               latCell,                  &
      !$acc               lonCell,                  &
      !$acc               xCell,                    &
      !$acc               yCell,                    &
      !$acc               zCell,                    &
      !$acc               latEdge,                  &
      !$acc               lonEdge,                  &
      !$acc               xEdge,                    &
      !$acc               yEdge,                    &
      !$acc               zEdge,                    &
      !$acc               latVertex,                &
      !$acc               lonVertex,                &
      !$acc               xVertex,                  &
      !$acc               yVertex,                  &
      !$acc               zVertex,                  &
      !$acc               fEdge,                    &
      !$acc               fVertex,                  &
      !$acc               fCell,                    &
      !$acc               dcEdge,                   &
      !$acc               dvEdge,                   &
      !$acc               areaCell,                 &
      !$acc               areaTriangle,             &
      !$acc               bottomDepth,              &
      !$acc               refBottomDepth,           &
      !$acc               refBottomDepthTopOfCell,  &
      !$acc               vertCoordMovementWeights, &
      !$acc               meshScalingDel2,          &
      !$acc               meshScalingDel4,          &
      !$acc               meshDensity,              &
      !$acc               angleEdge,                &
      !$acc               edgeMask,                 &
      !$acc               cellMask,                 &
      !$acc               vertexMask,               &
      !$acc               boundaryEdge,             &
      !$acc               boundaryCell,             &
      !$acc               boundaryVertex,           &
      !$acc               edgeSignOnCell,           &
      !$acc               edgeSignOnVertex,         &
      !$acc               highOrderAdvectionMask,   &
      !$acc               weightsOnEdge,            &
      !$acc               kiteAreasOnVertex,        &
      !$acc               edgeTangentVectors,       &
      !$acc               edgeNormalVectors,        &
      !$acc               localVerticalUnitVectors, &
      !$acc               advCoefs,                 &
      !$acc               advCoefs3rd,              &
      !$acc               derivTwo,                 &
      !$acc               cellTangentPlane,         &
      !$acc               coeffs_reconstruct)

      !$omp end master
      !$omp barrier

!-------------------------------------------------------------------------------

   end subroutine ocn_meshUpdateFields !}}}

!*******************************************************************************

end module ocn_mesh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
