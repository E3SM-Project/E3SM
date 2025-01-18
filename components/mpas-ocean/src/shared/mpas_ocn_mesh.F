!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! \file mpas_ocn_mesh.F
!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use mpas_dmpar
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_log

   use ocn_config
   use ocn_constants

   implicit none
   private

   !----------------------------------------------------------------------------
   !
   ! Public parameters
   !
   !----------------------------------------------------------------------------
   !{{{

   logical, public :: &
      onSphere        ! this mesh is on the sphere

   real (kind=RKIND), public :: &
      sphereRadius    ! radius of sphere for spherical meshes

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
      nVertLevelsP1 ! number of vertical interfaces (levels plus one)

   integer, public, dimension(:), allocatable :: &
      nCellsHalo, &! number of owned+halo(n) cells in local domain
      nEdgesHalo, &! number of owned+halo(n) edges in local domain
      nVerticesHalo      ! number of owned+halo(n) vertices in local domain

   integer, public, dimension(:), pointer :: &
      nEdgesOnEdge, &! number of edges connected to each edge point
      nEdgesOnCell, &! number of edges associated with each cell center
      minLevelCell, &! max ocean level at bottom of cell
      minLevelEdgeTop, &! max ocean level at top    of edge column
      minLevelEdgeBot, &! max ocean level at bottom of edge column
      minLevelVertexTop, &! max ocean level at top    of each vertex
      minLevelVertexBot, &! max ocean level at bottom of each vertex
      maxLevelCell, &! max ocean level at bottom of cell
      maxLevelEdgeTop, &! max ocean level at top    of edge column
      maxLevelEdgeBot, &! max ocean level at bottom of edge column
      maxLevelVertexTop, &! max ocean level at top    of each vertex
      maxLevelVertexBot, &! max ocean level at bottom of each vertex
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
      kiteIndexOnCell ! index of kite associated with each cell

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
      angleEdge, &! angle the edge normal makes with local east
      distanceToCoast ! distance to nearest coast cell 

   real(kind=RKIND), public, dimension(:), allocatable :: &
      invAreaCell   ! inverse of area of each cell

   ! Multiplicative masks and vectors for various conditions
   real(kind=RKIND), public, dimension(:, :), allocatable :: &
      edgeMask, &! mask to denote active edges    with depth
      cellMask, &! mask to denote active cells    with depth
      vertexMask, &! mask to denote active vertices with depth
      boundaryEdge, &! mask for boundary edges    at each level
      boundaryCell, &! mask for boundary cells    at each level
      boundaryVertex, &! mask for boundary vertices at each level
      edgeSignOnCell, &! sign of edge contributions to a cell
      edgeSignOnVertex ! sign of edge contributions to a vertex

   real(kind=RKIND), public, dimension(:, :), pointer :: &
      weightsOnEdge, &! weights on each edge
      kiteAreasOnVertex, &! real (vertexDegree nVertices)
      edgeTangentVectors, &! tangent unit vector at edge
      edgeNormalVectors, &! normal  unit vector at edge
      localVerticalUnitVectors ! local unit vector iin vertical

   real(kind=RKIND), public, dimension(:, :, :), pointer :: &
      cellTangentPlane, &! two vectors defining tangent plane at cell center
      coeffs_reconstruct  ! coeffs for reconstructing vectors at cell centers

   ! Derived mesh values needed
   real(kind=RKIND), public, dimension(:, :), allocatable :: &
      edgeAreaFractionOfCell

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

!***********************************************************************

contains

!***********************************************************************
!
!  ocn_meshCreate
!
!> \brief Creates the ocean mesh data structure on both host and device
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> This module creates and maintains public ocean mesh data
!> and ensures all mesh variables are copied to an accelerator device
!> if needed.
!
!-----------------------------------------------------------------------

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

      real (kind=RKIND) :: &
         maxDensityLocal, maxDensityGlobal ! temps for mesh density

      ! scalar pointers for retrieval, but convert to actual scalars in struct
      logical, pointer :: &
         on_a_sphere

      real (kind=RKIND), pointer :: &
         sphere_radius

      integer, pointer :: &! mesh dimensions
         nCellsTmp, &!
         nEdgesTmp, &!
         nVerticesTmp, &!
         maxEdgesTmp, &!
         maxEdges2Tmp, &!
         vertexDegreeTmp, &!
         nVertLevelsTmp, &!
         nVertLevelsP1Tmp !

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
         edgeSignOnCellTmp, &
         edgeSignOnVertexTmp, &
         boundaryEdgeTmp, &
         boundaryVertexTmp, &
         boundaryCellTmp

      !***
      !*** end of preamble, begin code
      !***

      ! We only support one block so test for condition here
      blockCount = 0
      block => domain%blocklist
      do while (associated(block))
         blockCount = blockCount + 1
         if (blockCount > 1) then
            call mpas_log_write( &
               'ocn_meshCreate: more than one block no longer supported', &
               MPAS_LOG_CRIT)
         endif
         block => block%next
      end do

      ! Reset to original block
      block => domain%blocklist

      ! retrieve the mpas mesh pool
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)

      !-----------------------------------------------------------------
      ! first set pointers/values for all mesh variables
      ! many variables already initialized based on read of mesh file
      !-----------------------------------------------------------------

      ! set all mesh properties
      call mpas_pool_get_config(meshPool, 'on_a_sphere', &
                                           on_a_sphere)
      call mpas_pool_get_config(meshPool, 'sphere_radius', &
                                           sphere_radius)

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
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', &
                                   nCellsArrayTmp)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', &
                                   nEdgesArrayTmp)
      call mpas_pool_get_dimension(meshPool, 'nVerticesArray', &
                                   nVerticesArrayTmp)

      ! translate scalar pointers to scalars in new mesh structure
      onSphere = on_a_sphere
      sphereRadius = sphere_radius
      maxEdges = maxEdgesTmp
      maxEdges2 = maxEdges2Tmp
      vertexDegree = vertexDegreeTmp
      nVertLevels = nVertLevelsTmp
      nVertLevelsP1 = nVertLevelsP1Tmp

      ! convert previous index limits into new halo definitions
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

      ! set pointers for a lot of connectivity info
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', &
                               nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', &
                               nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', &
                               minLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', &
                               minLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', &
                               minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'minLevelVertexTop', &
                               minLevelVertexTop)
      call mpas_pool_get_array(meshPool, 'minLevelVertexBot', &
                               minLevelVertexBot)
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
      call mpas_pool_get_array(meshPool, 'distanceToCoast', &
                               distanceToCoast)

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
      call mpas_pool_get_array(meshPool, 'cellTangentPlane', &
                               cellTangentPlane)
      call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', &
                               coeffs_reconstruct)

      ! For masks, we wish to convert to real multiplicative masks
      ! so retrieve integer version pointers and allocate real masks
      ! Once these are converted in Registry, we can eliminate this.
      call mpas_pool_get_array(meshPool, 'edgeMask', &
                               edgeMaskTmp)
      call mpas_pool_get_array(meshPool, 'vertexMask', &
                               vertexMaskTmp)
      call mpas_pool_get_array(meshPool, 'cellMask', &
                               cellMaskTmp)
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
            invAreaCell(nCellsAll+1))

      !-----------------------------------------------------------------
      ! Now that all pointers are set and mesh variables allocated
      ! we initialize other mesh quantities
      !-----------------------------------------------------------------

      ! Start by initializing vertical mesh, min/max cells and
      ! sign/index fields
      call ocn_meshVertCoordInit(domain)
      call ocn_meshMinMaxLevel()
      call ocn_meshSignIndexFields()

      ! Compute max mesh density and mesh scaling
      if (config_maxMeshDensity < 0.0_RKIND) then
         maxDensityLocal = maxval(meshDensity)
         call mpas_dmpar_max_real(domain%dminfo, maxDensityLocal, &
                                                 maxDensityGlobal)
         config_maxMeshDensity = maxDensityGlobal 
      endif
      call ocn_meshScaling()
   
      ! Routine calls above initialized the real mask arrays
      ! so convert the meshPool integer pointers here
      do n = 1, nCellsAll+1
      do k = 1, nVertLevels
         cellMaskTmp(k,n)     = nint(cellMask(k,n))
         boundaryCellTmp(k,n) = nint(boundaryCell(k,n))
      end do
      end do

      do n = 1, nCellsAll+1
      do k = 1, maxEdges
         edgeSignOnCellTmp(k, n) = nint(edgeSignOnCell(k,n))
      end do
      end do

      do n = 1, nEdgesAll+1
      do k = 1, nVertLevels
         edgeMaskTmp(k, n) = nint(edgeMask(k,n))
         boundaryEdgeTmp(k, n) = nint(boundaryEdge(k,n))
      end do
      end do

      do n = 1, nVerticesAll+1
      do k = 1, nVertLevels
         vertexMaskTmp(k, n) = nint(vertexMask(k,n))
         boundaryVertexTmp(k, n) = nint(boundaryVertex(k,n))
      end do
      end do

      do n = 1, nVerticesAll+1
      do k = 1, vertexDegree
         edgeSignOnVertexTmp(k,n) = nint(edgeSignOnVertex(k,n))
      end do
      end do

      ! Compute area weighting for some interpolation
      call ocn_meshAreaWeights()
      areaCell(nCellsAll+1) = -1.0e34_RKIND

      ! Compute the inverse of areaCell
      do n = 1, nCellsAll
         invAreaCell(n) = 1.0_RKIND / areaCell(n)
      end do
      invAreaCell(nCellsAll+1) = 0.0_RKIND

      !-----------------------------------------------------------------
      ! Copy mesh data to accelerator if needed.
      !-----------------------------------------------------------------

         !$acc enter data copyin( &
         !$acc                   onSphere,                 &
         !$acc                   sphereRadius,             &
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
         !$acc                   nEdgesHalo,               &
         !$acc                   nCellsHalo,               &
         !$acc                   nVerticesHalo,            &
         !$acc                   nEdgesOnEdge,             &
         !$acc                   nEdgesOnCell,             &
         !$acc                   minLevelCell,             &
         !$acc                   minLevelEdgeTop,          &
         !$acc                   minLevelEdgeBot,          &
         !$acc                   minLevelVertexTop,        &
         !$acc                   minLevelVertexBot,        &
         !$acc                   maxLevelCell,             &
         !$acc                   maxLevelEdgeTop,          &
         !$acc                   maxLevelEdgeBot,          &
         !$acc                   maxLevelVertexTop,        &
         !$acc                   maxLevelVertexBot,        &
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
         !$acc                   invAreaCell,              &
         !$acc                   areaTriangle,             &
         !$acc                   bottomDepth,              &
         !$acc                   refBottomDepth,           &
         !$acc                   refBottomDepthTopOfCell,  &
         !$acc                   vertCoordMovementWeights, &
         !$acc                   meshScalingDel2,          &
         !$acc                   meshScalingDel4,          &
         !$acc                   meshDensity,              &
         !$acc                   angleEdge,                &
         !$acc                   distanceToCoast,          &
         !$acc                   edgeMask,                 &
         !$acc                   cellMask,                 &
         !$acc                   vertexMask,               &
         !$acc                   boundaryEdge,             &
         !$acc                   boundaryCell,             &
         !$acc                   boundaryVertex,           &
         !$acc                   edgeSignOnCell,           &
         !$acc                   edgeSignOnVertex,         &
         !$acc                   weightsOnEdge,            &
         !$acc                   kiteAreasOnVertex,        &
         !$acc                   edgeTangentVectors,       &
         !$acc                   edgeNormalVectors,        &
         !$acc                   localVerticalUnitVectors, &
         !$acc                   cellTangentPlane,         &
         !$acc                   coeffs_reconstruct)

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

      ! First remove data from the device. Must remove components first,
      ! then remove the mesh type itself.

      !$acc exit data delete(onSphere,          &
      !$acc                  sphereRadius,      &
      !$acc                  nCellsAll,         &
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
      !$acc                  nEdgesHalo,        &
      !$acc                  nCellsHalo,        &
      !$acc                  nVerticesHalo,     &
      !$acc                  nEdgesOnEdge,      &
      !$acc                  nEdgesOnCell,      &
      !$acc                  minLevelCell,      &
      !$acc                  minLevelEdgeTop,   &
      !$acc                  minLevelEdgeBot,   &
      !$acc                  minLevelVertexTop, &
      !$acc                  minLevelVertexBot, &
      !$acc                  maxLevelCell,      &
      !$acc                  maxLevelEdgeTop,   &
      !$acc                  maxLevelEdgeBot,   &
      !$acc                  maxLevelVertexTop, &
      !$acc                  maxLevelVertexBot, &
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
      !$acc                  invAreaCell,       &
      !$acc                  areaTriangle,      &
      !$acc                  bottomDepth,       &
      !$acc                  refBottomDepth,    &
      !$acc                  refBottomDepthTopOfCell, &
      !$acc                  vertCoordMovementWeights, &
      !$acc                  meshScalingDel2,   &
      !$acc                  meshScalingDel4,   &
      !$acc                  meshDensity,       &
      !$acc                  angleEdge,         &
      !$acc                  distanceToCoast,   &
      !$acc                  edgeMask,          &
      !$acc                  cellMask,          &
      !$acc                  vertexMask,        &
      !$acc                  boundaryEdge,      &
      !$acc                  boundaryCell,      &
      !$acc                  boundaryVertex,    &
      !$acc                  edgeSignOnCell,    &
      !$acc                  edgeSignOnVertex,  &
      !$acc                  weightsOnEdge,     &
      !$acc                  kiteAreasOnVertex, &
      !$acc                  edgeTangentVectors,&
      !$acc                  edgeNormalVectors, &
      !$acc                  localVerticalUnitVectors, &
      !$acc                  cellTangentPlane,  &
      !$acc                  coeffs_reconstruct)

      ! Reset all scalars to zero
      onSphere  = .false.
      sphereRadius = 0.0_RKIND
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

      ! Now nullify all pointers to invalidate fields
      ! If this becomes the only mesh structure and mesh pool is eliminated,
      !  then we will want to deallocate here instead of nullify.

      nullify (nEdgesOnEdge, &
               nEdgesOnCell, &
               minLevelCell, &
               minLevelEdgeTop, &
               minLevelEdgeBot, &
               minLevelVertexTop, &
               minLevelVertexBot, &
               maxLevelCell, &
               maxLevelEdgeTop, &
               maxLevelEdgeBot, &
               maxLevelVertexTop, &
               maxLevelVertexBot, &
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
               distanceToCoast, &
               weightsOnEdge, &
               kiteAreasOnVertex, &
               edgeTangentVectors, &
               edgeNormalVectors, &
               localVerticalUnitVectors, &
               cellTangentPlane, &
               coeffs_reconstruct)

#ifdef MPAS_OPENACC
      !$acc exit data delete(edgeAreaFractionOfCell)
#endif
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
                  edgeAreaFractionOfCell, &
                  invAreaCell)

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
      integer k, n          ! loop indices
      integer, dimension(:, :), pointer :: &
         edgeMaskTmp, &
         vertexMaskTmp, &
         cellMaskTmp, &
         edgeSignOnCellTmp, &
         edgeSignOnVertexTmp, &
         boundaryEdgeTmp, &
         boundaryVertexTmp, &
         boundaryCellTmp

      !***
      !*** end of preamble, begin code
      !***

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
      !call mpas_pool_get_array(meshPool, 'distanceToCoast',   &
      !                                    distanceToCoast)
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

      ! Allocatable array was updated in ocn_init_routines_setup_sign_and_index_fields.  Now the
      ! pointer needs to be udpated
      do n = 1, nCellsAll+1
      do k = 1, maxEdges
         edgeSignOnCellTmp(k, n) = int(edgeSignOnCell(k, n))
      end do
      end do

      do n = 1, nEdgesAll+1
      do k = 1, nVertLevels
         edgeMask(k, n) = real(edgeMaskTmp(k, n), RKIND)
         boundaryEdge(k, n) = real(boundaryEdgeTmp(k, n), RKIND)
      end do
      end do

      do n = 1, nVerticesAll+1
      do k = 1, nVertLevels
         vertexMask(k, n) = real(vertexMaskTmp(k, n), RKIND)
         boundaryVertex(k, n) = real(boundaryVertexTmp(k, n), RKIND)
      end do
      end do

      ! Allocatable array was updated in ocn_init_routines_setup_sign_and_index_fields.  Now the
      ! pointer needs to be udpated
      do n = 1, nVerticesAll+1
      do k = 1, vertexDegree
         edgeSignOnVertexTmp(k, n) = &
            int(edgeSignOnVertex(k, n))
      end do
      end do

      ! Go ahead and update all fields on device to be safe.
      ! NOTE: if we end up computing some fields on the device during
      !       init, the update must go the opposite direction (update host)

      !$acc update device(onSphere,                 &
      !$acc               sphereRadius,             &
      !$acc               nCellsAll,                &
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
      !$acc               nEdgesHalo,               &
      !$acc               nCellsHalo,               &
      !$acc               nVerticesHalo,            &
      !$acc               nEdgesOnEdge,             &
      !$acc               nEdgesOnCell,             &
      !$acc               minLevelCell,             &
      !$acc               minLevelEdgeTop,          &
      !$acc               minLevelEdgeBot,          &
      !$acc               minLevelVertexTop,        &
      !$acc               minLevelVertexBot,        &
      !$acc               maxLevelCell,             &
      !$acc               maxLevelEdgeTop,          &
      !$acc               maxLevelEdgeBot,          &
      !$acc               maxLevelVertexTop,        &
      !$acc               maxLevelVertexBot,        &
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
      !$acc               invAreaCell,              &
      !$acc               areaTriangle,             &
      !$acc               bottomDepth,              &
      !$acc               refBottomDepth,           &
      !$acc               refBottomDepthTopOfCell,  &
      !$acc               vertCoordMovementWeights, &
      !$acc               meshScalingDel2,          &
      !$acc               meshScalingDel4,          &
      !$acc               meshDensity,              &
      !$acc               angleEdge,                &
      !$acc               distanceToCoast,          &
      !$acc               edgeMask,                 &
      !$acc               cellMask,                 &
      !$acc               vertexMask,               &
      !$acc               boundaryEdge,             &
      !$acc               boundaryCell,             &
      !$acc               boundaryVertex,           &
      !$acc               edgeSignOnCell,           &
      !$acc               edgeSignOnVertex,         &
      !$acc               weightsOnEdge,            &
      !$acc               kiteAreasOnVertex,        &
      !$acc               edgeTangentVectors,       &
      !$acc               edgeNormalVectors,        &
      !$acc               localVerticalUnitVectors, &
      !$acc               cellTangentPlane,         &
      !$acc               coeffs_reconstruct)

!-------------------------------------------------------------------------------

   end subroutine ocn_meshUpdateFields !}}}

!***********************************************************************
!
!  routine ocn_MeshMinMaxLevel
!
!> \brief  initialize max level and boundary mask variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes max level and boundary mask variables
!
!-----------------------------------------------------------------------

   subroutine ocn_meshMinMaxLevel()!{{{

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, iEdge, iVertex, &! loop indices for cells,edges,vertices
         i, k                    ! loop indices for neighbors, vertical

      ! End preamble
      !-------------
      ! Begin code

      ! minLevelEdgeTop is the minimum (shallowest) of surrounding cells
      ! if the water column is dry, set minLevelCell outside of bounds
      ! to prevent dry cells from determining minimum
      do iCell = 1, nCellsAll+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell)=nVertLevels+1 
      end do
      do iEdge = 1, nEdgesAll
         minLevelEdgeTop(iEdge) = &
            min( minLevelCell(cellsOnEdge(1,iEdge)), &
                 minLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      minLevelEdgeTop(nEdgesAll+1) = 0

      ! minLevelEdgeBot is the maximum (deepest) of surrounding cells
      ! prevent dry cells from determining maximum
      do iCell = 1, nCellsAll+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell) = 1 
      end do
      do iEdge = 1, nEdgesAll
         minLevelEdgeBot(iEdge) = &
            max( minLevelCell(cellsOnEdge(1,iEdge)), &
                 minLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      minLevelEdgeBot(nEdgesAll+1) = 0

      ! minLevelVertexBot is the maximum (deepest) of surrounding cells
      do iVertex = 1,nVerticesAll
         minLevelVertexBot(iVertex) = &
              minLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            minLevelVertexBot(iVertex) = &
               max( minLevelVertexBot(iVertex), &
                    minLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      minLevelVertexBot(nVerticesAll+1) = 0

      ! minLevelVertexTop is the minimum (shallowest) of surrounding cells
      ! if the water column is dry, set minLevelCell outside of bounds
      ! to prevent dry cells from determining minimum
      do iCell = 1, nCellsAll+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell) = nVertLevels+1 
      end do
      do iVertex = 1,nVerticesAll
         minLevelVertexTop(iVertex) = minLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            minLevelVertexTop(iVertex) = &
               min( minLevelVertexTop(iVertex), &
                    minLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      minLevelVertexTop(nVerticesAll+1) = 0

      ! if the water column is dry, set minLevelCell = 1
      do iCell = 1, nCellsAll+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell) = 1
      end do
      
      ! maxLevelEdgeTop is the minimum (shallowest) of surrounding cells
      do iEdge = 1, nEdgesAll
         maxLevelEdgeTop(iEdge) = &
            min( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      maxLevelEdgeTop(nEdgesAll+1) = 0

      ! maxLevelEdgeBot is the maximum (deepest) of surrounding cells
      do iEdge = 1, nEdgesAll
         maxLevelEdgeBot(iEdge) = &
            max( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      maxLevelEdgeBot(nEdgesAll+1) = 0

      ! maxLevelVertexBot is the maximum (deepest) of surrounding cells
      do iVertex = 1,nVerticesAll
         maxLevelVertexBot(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            maxLevelVertexBot(iVertex) = &
               max( maxLevelVertexBot(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      maxLevelVertexBot(nVerticesAll+1) = 0

      ! maxLevelVertexTop is the minimum (shallowest) of surrounding cells
      do iVertex = 1,nVerticesAll
         maxLevelVertexTop(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            maxLevelVertexTop(iVertex) = &
               min( maxLevelVertexTop(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      maxLevelVertexTop(nVerticesAll+1) = 0

      ! set boundary edge
      boundaryEdge(:,1:nEdgesAll+1)=1.0_RKIND
      edgeMask(:,1:nEdgesAll+1)=0.0_RKIND

      do iEdge = 1, nEdgesAll
      do k= minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
         boundaryEdge(k,iEdge)=0.0_RKIND
         edgeMask(k,iEdge)=1.0_RKIND
      end do
      end do

      ! Find cells and vertices that have an edge on the boundary
      boundaryCell(:,1:nCellsAll+1) = 0.0_RKIND
      cellMask(:,1:nCellsAll+1) = 0.0_RKIND
      boundaryVertex(:,1:nVerticesAll+1) = 0.0_RKIND
      vertexMask(:,1:nVerticesAll+1) = 0.0_RKIND

      do iEdge = 1, nEdgesAll
      do k = 1, nVertLevels
         if (boundaryEdge(k,iEdge) == 1) then
            boundaryCell(k,cellsOnEdge(1,iEdge)) = 1.0_RKIND
            boundaryCell(k,cellsOnEdge(2,iEdge)) = 1.0_RKIND
            boundaryVertex(k,verticesOnEdge(1,iEdge)) = 1.0_RKIND
            boundaryVertex(k,verticesOnEdge(2,iEdge)) = 1.0_RKIND
         endif
      end do
      end do

      do iCell = 1, nCellsAll
      do k = 1, nVertLevels
         if ( minLevelCell(iCell) <= k .and. &
              maxLevelCell(iCell) >= k ) then
            cellMask(k, iCell) = 1.0_RKIND
         end if
      end do
      end do

      do iVertex = 1, nVerticesAll
      do k = 1, nVertLevels
         if ( minLevelVertexTop(iVertex) <= k .and. &
              maxLevelVertexBot(iVertex) >= k ) then
            vertexMask(k, iVertex) = 1.0_RKIND
         end if
      end do
      end do

      ! Note: We do not update halos on maxLevel* variables.  We want
      ! the outside edge of a halo to be zero on each processor.

   !--------------------------------------------------------------------

   end subroutine ocn_meshMinMaxLevel !}}}

!***********************************************************************
!
!  routine ocn_meshSignIndexFields
!
!> \brief   set up sign and index fields
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes edgeSignOnCell, edgeSignOnVertex, and
!>  kiteIndexOnCell.
!
!-----------------------------------------------------------------------

   subroutine ocn_meshSignIndexFields()!{{{

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      integer :: iCell, iEdge, iVertex, i, j

      ! End preamble
      !-------------
      ! Begin code

      ! Initialize to zero
      edgeSignOnCell   = 0.0_RKIND
      edgeSignOnVertex = 0.0_RKIND
      kiteIndexOnCell  = 0

      do iCell = 1, nCellsAll
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            iVertex = verticesOnCell(i, iCell)

            ! Vector points from cell 1 to cell 2
            if (iCell == cellsOnEdge(1, iEdge)) then
               edgeSignOnCell(i, iCell) = -1.0_RKIND
            else
               edgeSignOnCell(i, iCell) =  1.0_RKIND
            end if

            do j = 1, vertexDegree
               if (cellsOnVertex(j,iVertex) == iCell) then
                  kiteIndexOnCell(i,iCell) = j
               end if
            end do
         end do
      end do

      do iVertex = 1, nVerticesAll
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)

            ! Vector points from vertex 1 to vertex 2
            if (iVertex == verticesOnEdge(1,iEdge)) then
               edgeSignOnVertex(i,iVertex) = -1.0_RKIND
            else
               edgeSignOnVertex(i,iVertex) =  1.0_RKIND
            end if
         end do
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_meshSignIndexFields !}}}

!***********************************************************************
!
!  routine ocn_meshVertCoordInit
!
!> \brief  initialize vertical coordinate variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes vertical coordinate variables, including
!>  zlevel-type variables and adjusted initial conditions for partial
!>  bottom cells.
!
!-----------------------------------------------------------------------

   subroutine ocn_meshVertCoordInit(domain)!{{{

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain !< [inout] ocean state

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block ! all ocn data in a block

      type (mpas_pool_type), pointer :: &
         verticalMeshPool, &! vertical mesh data
         statePool,        &! ocean state (layer thickness)
         forcingPool        ! forcing data (land ice mask)

      integer :: &
         iCell, k,         &! loop indices for cell, vertical loops
         kmin, kmax         ! min, max active levels in column

      logical :: &
         consistentSSH,    &! flag if SSH consistent with depth
         landIceFlag,      &! flag if land ice mask exists
         checkSSH           ! flag for checking SSH in cell column

      real (kind=RKIND) &
         thickDepth         ! depth based on sum of layer thicknesses

      integer, dimension(:), pointer :: &
         landIceMask        ! land ice mask

      real (kind=RKIND), dimension(:), pointer :: &
         refZMid,                  &! reference depth at mid-layer
         refLayerThickness          ! reference profile layer thickness

      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThickness     ! current layer thickness

      ! End preamble
      !-------------
      ! Begin code

      ! Extract pools from domain
      block => domain % blocklist
      call mpas_pool_get_subpool(block%structs, 'state', statePool)
      call mpas_pool_get_subpool(block%structs, 'verticalMesh', &
                                                 verticalMeshPool)
      call mpas_pool_get_subpool(block%structs, 'forcing', forcingPool)

      ! Extract needed variables from pools
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThickness, 1)

      call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
      call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', &
                                                  refLayerThickness)

      call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
      if (associated(landIceMask)) then
         landIceFlag = .true.  ! ice shelves
      else
         landIceFlag = .false. ! ice shelves
      endif

      ! Initialize reference z coordinates based on refBottomDepth
      refBottomDepthTopOfCell(1) = 0.0_RKIND
      do k = 1, nVertLevels
         refBottomDepthTopOfCell(k+1) = refBottomDepth(k)
         refLayerThickness(k) = refBottomDepth(k) - &
                                refBottomDepthTopOfCell(k)
         refZMid(k) = - refBottomDepthTopOfCell(k) - &
                        refLayerThickness(k)/2.0_RKIND
      end do

      ! Initialization of vertCoordMovementWeights. 
      ! This determines how SSH perturbations are distributed throughout
      ! the column. Check for invalid choice of vert coord movement.

      call mpas_log_write(' Vertical coordinate movement is: ' // &
                          trim(config_vert_coord_movement))


      select case (trim(config_vert_coord_movement))
      case ('fixed')

         vertCoordMovementWeights = 0.0_RKIND
         vertCoordMovementWeights(1) = 1.0_RKIND

      case ('uniform_stretching')

         vertCoordMovementWeights = 1.0_RKIND

      case ('tapered')

         ! Set weight tapering:
         !  1.0 shallower than config_vert_taper_weight_depth_1
         !  linear in between
         !  0.0 deeper than config_vert_taper_weight_depth_2
         do k = 1, nVertLevels
            vertCoordMovementWeights(k) = 1.0_RKIND + &
                 (refZMid(k) + config_vert_taper_weight_depth_1) / &
                 (config_vert_taper_weight_depth_2 - &
                  config_vert_taper_weight_depth_1)
            ! limit range to (0,1)
            vertCoordMovementWeights(k) = max( 0.0_RKIND, &
                      min( 1.0_RKIND, vertCoordMovementWeights(k) ) )
         end do

      case ('impermeable_interfaces')

         ! weights are never used, no init needed

      case ('user_specified')

         ! weights should be input? probably should check

      case default

         call mpas_log_write( &
            'Incorrect choice of config_vert_coord_movement.', &
             MPAS_LOG_CRIT)

      end select

      !-----------------------------------------------------------------
      ! Check consistency and validity of options
      !-----------------------------------------------------------------

      ! SSH consistency
      if (config_check_ssh_consistency) then

         ! Check if abs(ssh)>20m.  If so, print warning.
         consistentSSH = .true.
         do iCell = 1,nCellsAll

            ! Only check for open ocean when land ice is being used
            checkSSH = .true.
            if (landIceFlag) then
               if (landIceMask(iCell) /= 0) checkSSH = .false.
            endif

            ! Check whether thickness is consistent with bottom depth
            if (checkSSH) then
               kmin = minLevelCell(iCell)
               kmax = maxLevelCell(iCell)

               thickDepth = sum(layerThickness(kmin:kmax,iCell))

               if (abs(thickDepth - bottomDepth(iCell)) > 20.0_RKIND) then
                  consistentSSH = .false.

                  call mpas_log_write( &
                     ' Warning: Sea surface height outside of acceptable range (20m)', &
                     MPAS_LOG_ERR)
                  call mpas_log_write( &
                     ' iCell: $i, minLevelCell(iCell), maxLevelCell(iCell): $i, $i, bottomDepth(iCell): $r, sum(h): $r', &
                     intArgs=(/iCell, minLevelCell(iCell), maxLevelCell(iCell) /), &
                     realArgs=(/ bottomDepth(iCell), thickDepth /) )
               endif ! depth diff check
            endif ! checkSSH
         end do ! cell loop

         if (.not. consistentSSH) call mpas_log_write( &
            'Warning: SSH not consistent - ' // &
            'initial layerThickness does not match bottomDepth.')

      endif ! config_check_ssh_consistency

      ! Z-level consistency
      if (config_check_zlevel_consistency) then
         do iCell = 1,nCellsAll
            ! Check that bottomDepth and maxLevelCell match.  Some
            ! older meshs do not have the bottomDepth variable.
            kmax = maxLevelCell(iCell)
            if (bottomDepth(iCell) > refBottomDepth(kmax) .or. &
                bottomDepth(iCell) < refBottomDepthTopOfCell(kmax)) then
               call mpas_log_write( &
                 ' fatal error: bottomDepth and maxLevelCell do not match:',  &
                 MPAS_LOG_ERR)
               call mpas_log_write( &
                 ' iCell: $i, maxLevelCell(iCell): $i, bottomDepth(iCell): $r', &
                 intArgs=(/iCell, kmax /), &
                 realArgs=(/ bottomDepth(iCell) /) )
            endif
         enddo ! cell loop
      endif ! check_zlevel_consistency

      ! pressure gradient compatibility
      if (config_vert_coord_movement /= 'impermeable_interfaces' .and. &
          config_pressure_gradient_type == 'MontgomeryPotential') then
         call mpas_log_write( 'Incompatible choice of ' // &
            'impermeable vertical interfaces and Montgomery Potential',&
            MPAS_LOG_CRIT)
      end if

      ! barotropic filter mode compatibility
      if (config_filter_btr_mode .and. &
          config_vert_coord_movement /= 'fixed')then
         call mpas_log_write( &
          'filter_btr_mode has only been tested with ' // &
          'config_vert_coord_movement=fixed.', &
          MPAS_LOG_CRIT)
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_meshVertCoordInit !}}}

!***********************************************************************
!
!  routine ocn_meshAreaWeights
!
!> \brief  set up area weighting
!> \author Mark Petersen
!> \date   June 2020
!> \details
!>  This routine initializes edgeAreaFractionOfCell which is defined as
!>  the fractional area of cell iCell encompassed by the triangle with
!>  edge iEdge connected to the cell center of cell iCell. On a perfect
!>  planar hex mesh this is always 1/6. This weighting is used to
!>  interpolate scalars from edges to cell centers.
!>  Use 2*edgeAreaFractionOfCell to interpolate normal vectors at
!>  edges to vector norms at cell centers.
!
!-----------------------------------------------------------------------

   subroutine ocn_meshAreaWeights()!{{{

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, iEdge, &! cell, edge addresses
         i,            &! edge on cell index
         ierr           ! internal error flag

      ! End preamble
      !-------------
      ! Begin code

      ! Allocate space
      allocate (edgeAreaFractionOfCell(maxEdges,nCellsAll), STAT=ierr)
#ifdef MPAS_OPENACC
      !$acc enter data create(edgeAreaFractionOfCell)
#endif
      if (ierr /= 0) call mpas_log_write( &
         'Error allocating edgeAreaFractionOfCell in ocn_mesh', &
          MPAS_LOG_CRIT)

      do iCell = 1, nCellsAll
      do i = 1, nEdgesOnCell(iCell)
         iEdge = edgesOnCell(i, iCell)
         edgeAreaFractionOfCell(i, iCell) = 0.25_RKIND* &
               dcEdge(iEdge)*dvEdge(iEdge)/areaCell(iCell)
      end do
      end do
#ifdef MPAS_OPENACC
      !$acc update device(edgeAreaFractionOfCell)
#endif

   !--------------------------------------------------------------------

   end subroutine ocn_meshAreaWeights !}}}

!***********************************************************************
!
!  routine ocn_meshScaling
!
!> \brief  Computes mesh scaling variables for mixing
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes meshScalingDel2 and meshScalingDel4
!>  that are used to scale coefficients with mesh size.
!
!-----------------------------------------------------------------------

   subroutine ocn_meshScaling() !{{{

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iEdge,       &! loop index/address for edges
         cell1, cell2  ! cell addresses for neighbors across edge

      real (kind=RKIND) :: &
         cellWidth,   &! reference cell width for scaling
         avgDensity    ! avg mesh density across edge

      ! End preamble
      !-------------
      ! Begin code

      if (config_hmix_scaleWithMesh) then
         if (config_hmix_use_ref_cell_width) then
            ! Mesh scaling is set by areaCell and and an input
            ! reference widht config_hmix_ref_cell_width
            ! See description of config_hmix_ref_cell_width in
            ! Registry.xml for more detail.
            do iEdge = 1, nEdgesAll
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! Effective cell width at edge iEdge, assuming
               ! neighboring cells are circles for this calculation.
               cellWidth = 2.0_RKIND* &
                           sqrt((areaCell(cell1) + areaCell(cell2))/ &
                                2.0_RKIND/pi)
               meshScalingDel2(iEdge) =  cellWidth/ &
                                         config_hmix_ref_cell_width
               meshScalingDel4(iEdge) = (cellWidth/ &
                                         config_hmix_ref_cell_width)**3
            end do

         else ! not using ref cell width
            ! Mesh scaling is set by meshDensity. This is both confusing
            ! and inconvenient, as the flags like config_mom_del2 need
            ! to be reset for every resolution. It is kept for backwards
            ! compatibility, but should become defunct.
            ! del2 scales as dc**1, del4 scales as dc**3
            do iEdge = 1, nEdgesAll
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               avgDensity = (meshDensity(cell1) + meshDensity(cell2))/ &
                            2.0_RKIND
               meshScalingDel2(iEdge) = 1.0_RKIND/ &
                  (avgDensity/config_maxMeshDensity)**(1.0_RKIND/4.0_RKIND)
               meshScalingDel4(iEdge) = 1.0_RKIND/ &
                  (avgDensity/config_maxMeshDensity)**(3.0_RKIND/4.0_RKIND)
            end do
         end if

      else ! no scaling with mesh

         ! If config_hmix_scaleWithMesh is false, hmix coefficients
         ! do not vary with cell size but remain constant across domain.
         do iEdge = 1, nEdgesAll
            meshScalingDel2(iEdge) = 1.0_RKIND
            meshScalingDel4(iEdge) = 1.0_RKIND
         end do
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_meshScaling !}}}

!***********************************************************************

end module ocn_mesh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
