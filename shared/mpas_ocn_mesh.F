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
!!  \author Rob Aulwes and Phil Jones
!!  \date   14 Jan 2020
!!  \details
!!  This module creates and maintains a primary ocean mesh structure
!!  and ensures all mesh variables are copied to an accelerator device
!!  if needed. Currently it consists of pointers to the existing MPAS mesh pool
!!  variables, but is intended to eventually replace the mesh pool later.
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

   type, public :: meshType !< mesh data structure with all needed mesh fields

      integer ::         &! mesh, array sizes
         nCellsAll,      &! total number of local (owned+halo) cells in primary
         nEdgesAll,      &! total number of local edge midpoints
         nVerticesAll,   &! total number of local cells in dual (cell vertices)
         nCellsOwned,    &! number of cells    owned by the local domain
         nEdgesOwned,    &! number of edges    owned by the local domain
         nVerticesOwned, &! number of vertices owned by the local domain
         maxEdges,       &! largest number of edges any polygon has
         maxEdges2,      &! 2x the largest number of edges any polygon has
         vertexDegree,   &! number of cells or edges touching each vertex
         nVertLevels,    &! number of vertical levels
         nVertLevelsP1,  &! number of vertical interfaces (levels plus one)
         nAdvectionCells  ! largest number of advection cells for any edge

      integer, dimension(:), allocatable :: & 
         nCellsHalo,       &! number of owned+halo(n) cells in local domain
         nEdgesHalo,       &! number of owned+halo(n) edges in local domain
         nVerticesHalo      ! number of owned+halo(n) vertices in local domain

      integer, dimension(:), pointer :: & 
         nEdgesOnEdge,     &! number of edges connected to each edge point
         nEdgesOnCell,     &! number of edges associated with each cell center
         maxLevelCell,     &! max ocean level at bottom of cell
         maxLevelEdgeTop,  &! max ocean level at top    of edge column
         maxLevelEdgeBot,  &! max ocean level at bottom of edge column
         maxLevelVertexTop,&! max ocean level at top    of each vertex
         maxLevelVertexBot,&! max ocean level at bottom of each vertex
         nAdvCellsForEdge, &! number of cells contrib to advection at edge
         indexToCellID,    &! global ID of each local cell
         indexToEdgeID,    &! global ID of each local edge
         indexToVertexID    ! global ID of each local vertex

      integer, dimension(:,:), pointer :: &
         edgesOnEdge,    &! index of edges connected to each edge
         cellsOnEdge,    &! index of cells connected to each edge
         verticesOnEdge, &! index of vertices connected to each edge
         cellsOnCell,    &! index of cells connected to each cell
         edgesOnCell,    &! index of edges connected to each cell
         verticesOnCell, &! index of vertices connected to each cell
         cellsOnVertex,  &! index of cells connected to each vertex
         edgesOnVertex,  &! index of edges connected to each vertex
         kiteIndexOnCell,&! index of kite associated with each cell
         advCellsForEdge  ! index of cells contrib to advective flux at edge

      real (kind=RKIND), dimension(:), pointer :: &
         latCell,        &! latitude  of cell centers
         lonCell,        &! longitude of cell centers
         xCell,          &! Cartesian x coord of cell center
         yCell,          &! Cartesian y coord of cell center
         zCell,          &! Cartesian z coord of cell center
         latEdge,        &! latitude  of edge
         lonEdge,        &! longitude of edge
         xEdge,          &! Cartesian x coord of edge
         yEdge,          &! Cartesian y coord of edge
         zEdge,          &! Cartesian z coord of edge
         latVertex,      &! latitude  of vertex
         lonVertex,      &! longitude of vertex
         xVertex,        &! Cartesian coord of vertex
         yVertex,        &! Cartesian y coord of vertex
         zVertex,        &! Cartesian z coord of vertex
         fEdge,          &! Coriolus parameter at edge
         fVertex,        &! Coriolus parameter at vertex
         fCell,          &! Coriolus parameter at cell center
         dcEdge,         &! length of edge = dist between cells across edge
         dvEdge,         &! length of edge = dist between vertices along edge
         areaCell,       &! area of each cell
         areaTriangle,   &! area of each cell on dual grid
         bottomDepth,    &! ocean bottom depth at each cell center
         refBottomDepth, &! ocean depth at bottom of cell for reference profile
         refBottomDepthTopOfCell, &! depth at top of cell for reference profile
         vertCoordMovementWeights, &! weights for distributing height perturb
         meshScalingDel2,&! mesh scaling factor for use in del2 diffusion
         meshScalingDel4,&! mesh scaling factor for use in del4 diffusion
         meshScaling,    &! general scaling of mesh size
         meshDensity,    &! density of mesh
         angleEdge        ! angle the edge normal makes with local east

      ! Multiplicative masks and vectors for various conditions
      real (kind=RKIND), dimension(:,:), allocatable :: &
         edgeMask,       &! mask to denote active edges    with depth
         cellMask,       &! mask to denote active cells    with depth
         vertexMask,     &! mask to denote active vertices with depth
         boundaryEdge,   &! mask for boundary edges    at each level
         boundaryCell,   &! mask for boundary cells    at each level
         boundaryVertex, &! mask for boundary vertices at each level
         edgeSignOnCell, &! sign of edge contributions to a cell
         edgeSignOnVertex,&! sign of edge contributions to a vertex
         highOrderAdvectionMask ! mask for high order advection contributions

      real (kind=RKIND), dimension(:,:), pointer :: &
         weightsOnEdge,     &! weights on each edge
         kiteAreasOnVertex, &! real (vertexDegree nVertices)
         edgeTangentVectors,&! tangent unit vector at edge
         edgeNormalVectors, &! normal  unit vector at edge
         localVerticalUnitVectors, &! local unit vector iin vertical
         advCoefs,          &! mesh-based advection coefficients
         advCoefs3rd         ! mesh-based advection coeffs for high order

      real (kind=RKIND), dimension(:,:,:), pointer :: &
         derivTwo,          &! 2nd derivative of edge reconstruction polynomial
         cellTangentPlane,  &! two vectors defining tangent plane at cell center
         coeffs_reconstruct  ! coeffs for reconstructing vectors at cell centers

   end type meshType

   ! default mesh must be public for now to provide identical namespace
   ! on CPU, accelerator
   type (meshType), public :: &
      ocnMesh              ! default mesh

   !----------------------------------------------------------------------------
   !
   ! Public member functions
   !
   !----------------------------------------------------------------------------

   public ::                     &
      ocn_meshCreate,            &
      ocn_meshUpdateFields,      &
      ocn_meshDestroy

!*******************************************************************************

contains

!*******************************************************************************
!
!  ocn_meshCreate
!
!> \brief Creates the primary ocean mesh data structure on both host and device
!!  \author Rob Aulwes and Phil Jones
!!  \date   14 Jan 2020
!!  \details
!!  This module creates and maintains a primary ocean mesh structure
!!  and ensures all mesh variables are copied to an accelerator device
!!  if needed.
!
!-------------------------------------------------------------------------------

   subroutine ocn_meshCreate(domain)

      ! Input arguments

      type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

      !{{{
      ! Local variables

      integer :: &
         blockCount               ! counter for number of blocks

      type(block_type), pointer :: &
         block                    ! variables in current subblock

      type (mpas_pool_type), pointer :: &
         meshPool                 ! mesh variables in MPAS pool structure

      ! scalar pointers for retrieval, but convert to actual scalars in struct
      integer, pointer ::   &! mesh dimensions
         nCellsTmp,         &!
         nEdgesTmp,         &!
         nVerticesTmp,      &!
         maxEdgesTmp,       &!
         maxEdges2Tmp,      &!
         vertexDegreeTmp,   &!
         nVertLevelsTmp,    &!
         nVertLevelsP1Tmp,  &!
         nAdvectionCellsTmp  !

      ! temporary pointers for converting index arrays
      integer, dimension(:), pointer :: &
         nCellsArrayTmp,                &
         nEdgesArrayTmp,                &
         nVerticesArrayTmp

      ! temporary pointers for converting masks
      integer i,k,n          ! loop indices
      integer, dimension(:,:), pointer :: &
         edgeMaskTmp,                     &
         vertexMaskTmp,                   &
         cellMaskTmp,                     &
         highOrderAdvectionMaskTmp,       &
         edgeSignOnCellTmp,               &
         edgeSignOnVertexTmp,             &
         boundaryEdgeTmp,                 &
         boundaryVertexTmp,               &
         boundaryCellTmp

      !***
      !*** end of preamble, begin code
      !***

      !$omp master

      blockCount = 0
      block => domain % blocklist
      do while ( associated(block) )

         !*** we now only support one block per task
         blockCount = blockCount + 1
         if (blockCount > 1) then
            call mpas_log_write(&
               'ocn_meshCreate: more than one block no longer supported', &
                MPAS_LOG_CRIT)
         endif

         ! retrieve the mpas mesh pool
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

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
         ocnMesh%maxEdges        = maxEdgesTmp
         ocnMesh%maxEdges2       = maxEdges2Tmp
         ocnMesh%vertexDegree    = vertexDegreeTmp
         ocnMesh%nVertLevels     = nVertLevelsTmp
         ocnMesh%nVertLevelsP1   = nVertLevelsP1Tmp
         ocnMesh%nAdvectionCells = nAdvectionCellsTmp

         ! convert previous index limits into new arrangement
         ocnMesh%nCellsAll       = nCellsTmp
         ocnMesh%nEdgesAll       = nEdgesTmp
         ocnMesh%nVerticesAll    = nVerticesTmp

         n = size(nCellsArrayTmp)
         allocate(ocnMesh%nCellsHalo(n-1))
         ocnMesh%nCellsOwned = nCellsArrayTmp(1)
         do i=2,n
            ocnMesh%nCellsHalo(i-1) = nCellsArrayTmp(i) 
         end do

         n = size(nEdgesArrayTmp)
         allocate(ocnMesh%nEdgesHalo(n-1))
         ocnMesh%nEdgesOwned = nEdgesArrayTmp(1)
         do i=2,n
            ocnMesh%nEdgesHalo(i-1) = nEdgesArrayTmp(i) 
         end do

         n = size(nVerticesArrayTmp)
         allocate(ocnMesh%nVerticesHalo(n-1))
         ocnMesh%nVerticesOwned = nVerticesArrayTmp(1)
         do i=2,n
            ocnMesh%nVerticesHalo(i-1) = nVerticesArrayTmp(i) 
         end do

         ! set a lot of connectivity info
         call mpas_pool_get_array(meshPool, 'nEdgesOnEdge',      &
                                     ocnMesh%nEdgesOnEdge)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell',      &
                                     ocnMesh%nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell',      &
                                     ocnMesh%maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',   &
                                     ocnMesh%maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot',   &
                                     ocnMesh%maxLevelEdgeBot)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', &
                                     ocnMesh%maxLevelVertexTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', &
                                     ocnMesh%maxLevelVertexBot)
         call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge',  &
                                     ocnMesh%nAdvCellsForEdge)
         call mpas_pool_get_array(meshPool, 'indexToCellID',     &
                                     ocnMesh%indexToCellID)
         call mpas_pool_get_array(meshPool, 'indexToEdgeID',     &
                                     ocnMesh%indexToEdgeID)
         call mpas_pool_get_array(meshPool, 'indexToVertexID',   &
                                     ocnMesh%indexToVertexID)
         call mpas_pool_get_array(meshPool, 'edgesOnEdge',       &
                                     ocnMesh%edgesOnEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge',       &
                                     ocnMesh%cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'verticesOnEdge',    &
                                     ocnMesh%verticesOnEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnCell',       &
                                     ocnMesh%cellsOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell',       &
                                     ocnMesh%edgesOnCell)
         call mpas_pool_get_array(meshPool, 'verticesOnCell',    &
                                     ocnMesh%verticesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnVertex',     &
                                     ocnMesh%cellsOnVertex)
         call mpas_pool_get_array(meshPool, 'edgesOnVertex',     &
                                     ocnMesh%edgesOnVertex)
         call mpas_pool_get_array(meshPool, 'kiteIndexOnCell',   &
                                     ocnMesh%kiteIndexOnCell)
         call mpas_pool_get_array(meshPool, 'advCellsForEdge',   &
                                     ocnMesh%advCellsForEdge)

         ! now set a number of physics and numerical properties of mesh
         call mpas_pool_get_array(meshPool, 'latCell',           &
                                     ocnMesh%latCell)
         call mpas_pool_get_array(meshPool, 'lonCell',           &
                                     ocnMesh%lonCell)
         call mpas_pool_get_array(meshPool, 'xCell',             &
                                     ocnMesh%xCell)
         call mpas_pool_get_array(meshPool, 'yCell',             &
                                     ocnMesh%yCell)
         call mpas_pool_get_array(meshPool, 'zCell',             &
                                     ocnMesh%zCell)
         call mpas_pool_get_array(meshPool, 'latEdge',           &
                                     ocnMesh%latEdge)
         call mpas_pool_get_array(meshPool, 'lonEdge',           &
                                     ocnMesh%lonEdge)
         call mpas_pool_get_array(meshPool, 'xEdge',             &
                                     ocnMesh%xEdge)
         call mpas_pool_get_array(meshPool, 'yEdge',             &
                                     ocnMesh%yEdge)
         call mpas_pool_get_array(meshPool, 'zEdge',             &
                                     ocnMesh%zEdge)
         call mpas_pool_get_array(meshPool, 'latVertex',         &
                                     ocnMesh%latVertex)
         call mpas_pool_get_array(meshPool, 'lonVertex',         &
                                     ocnMesh%lonVertex)
         call mpas_pool_get_array(meshPool, 'xVertex',           &
                                     ocnMesh%xVertex)
         call mpas_pool_get_array(meshPool, 'yVertex',           &
                                     ocnMesh%yVertex)
         call mpas_pool_get_array(meshPool, 'zVertex',           &
                                     ocnMesh%zVertex)
         call mpas_pool_get_array(meshPool, 'fEdge',             &
                                     ocnMesh%fEdge)
         call mpas_pool_get_array(meshPool, 'fVertex',           &
                                     ocnMesh%fVertex)
         call mpas_pool_get_array(meshPool, 'fCell',             &
                                     ocnMesh%fCell)
         call mpas_pool_get_array(meshPool, 'dcEdge',            &
                                     ocnMesh%dcEdge)
         call mpas_pool_get_array(meshPool, 'dvEdge',            &
                                     ocnMesh%dvEdge)
         call mpas_pool_get_array(meshPool, 'areaCell',          &
                                     ocnMesh%areaCell)
         call mpas_pool_get_array(meshPool, 'areaTriangle',      &
                                     ocnMesh%areaTriangle)
         call mpas_pool_get_array(meshPool, 'weightsOnEdge',     &
                                     ocnMesh%weightsOnEdge)
         call mpas_pool_get_array(meshPool, 'bottomDepth',       &
                                     ocnMesh%bottomDepth)
         call mpas_pool_get_array(meshPool, 'refBottomDepth',    &
                                     ocnMesh%refBottomDepth)
         call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell',   &
                                     ocnMesh%refBottomDepthTopOfCell)
         call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights',  &
                                     ocnMesh%vertCoordMovementWeights)
         call mpas_pool_get_array(meshPool, 'meshScalingDel2',   &
                                     ocnMesh%meshScalingDel2)
         call mpas_pool_get_array(meshPool, 'meshScalingDel4',   &
                                     ocnMesh%meshScalingDel4)
         call mpas_pool_get_array(meshPool, 'meshScaling',       &
                                     ocnMesh%meshScaling)
         call mpas_pool_get_array(meshPool, 'meshDensity',       &
                                     ocnMesh%meshDensity)
         call mpas_pool_get_array(meshPool, 'angleEdge',         &
                                     ocnMesh%angleEdge)
           
         call mpas_pool_get_array(meshPool, 'weightsOnEdge',             &
                                     ocnMesh%weightsOnEdge)
         call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex',         &
                                     ocnMesh%kiteAreasOnVertex)
         call mpas_pool_get_array(meshPool, 'edgeTangentVectors',        &
                                     ocnMesh%edgeTangentVectors)
         call mpas_pool_get_array(meshPool, 'edgeNormalVectors',         &
                                     ocnMesh%edgeNormalVectors)
         call mpas_pool_get_array(meshPool, 'localVerticalUnitVectors',  &
                                     ocnMesh%localVerticalUnitVectors)
         call mpas_pool_get_array(meshPool, 'advCoefs',                  &
                                     ocnMesh%advCoefs)
         call mpas_pool_get_array(meshPool, 'advCoefs3rd',               &
                                     ocnMesh%advCoefs3rd)
         call mpas_pool_get_array(meshPool, 'derivTwo',                  &
                                     ocnMesh%derivTwo)
         call mpas_pool_get_array(meshPool, 'cellTangentPlane',          &
                                     ocnMesh%cellTangentPlane)
         call mpas_pool_get_array(meshPool, 'coeffs_reconstruct',        &
                                     ocnMesh%coeffs_reconstruct)

         ! For masks, we wish to convert to real multiplicative masks
         ! so retrieve integer version and set real version in struct.
         ! Once these are converted in Registry, we can eliminate conversion.
         call mpas_pool_get_array(meshPool, 'edgeMask',        &
                                             edgeMaskTmp)
         call mpas_pool_get_array(meshPool, 'vertexMask',        &
                                             vertexMaskTmp)
         call mpas_pool_get_array(meshPool, 'cellMask',        &
                                             cellMaskTmp)
         call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask',        &
                                             highOrderAdvectionMaskTmp)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell',  &
                                             edgeSignOnCellTmp)
         call mpas_pool_get_array(meshPool, 'edgeSignOnVertex',  &
                                             edgeSignOnVertexTmp)
         call mpas_pool_get_array(meshPool, 'boundaryEdge',  &
                                             boundaryEdgeTmp)
         call mpas_pool_get_array(meshPool, 'boundaryVertex',  &
                                             boundaryVertexTmp)
         call mpas_pool_get_array(meshPool, 'boundaryCell',  &
                                             boundaryCellTmp)

         allocate( &
            ocnMesh%edgeMask      (ocnMesh%nVertLevels, ocnMesh%nEdgesAll),    &
            ocnMesh%cellMask      (ocnMesh%nVertLevels, ocnMesh%nCellsAll),    &
            ocnMesh%vertexMask    (ocnMesh%nVertLevels, ocnMesh%nVerticesAll), &
            ocnMesh%boundaryEdge  (ocnMesh%nVertLevels, ocnMesh%nEdgesAll),    &
            ocnMesh%boundaryCell  (ocnMesh%nVertLevels, ocnMesh%nCellsAll),    &
            ocnMesh%boundaryVertex(ocnMesh%nVertLevels, ocnMesh%nVerticesAll), &
            ocnMesh%edgeSignOnCell(ocnMesh%maxEdges,    ocnMesh%nCellsAll),    &
            ocnMesh%edgeSignOnVertex(ocnMesh%maxEdges,  ocnMesh%nVerticesAll), &
            ocnMesh%highOrderAdvectionMask(ocnMesh%nVertLevels, &
                                           ocnMesh%nEdgesAll))

         do n=1,ocnMesh%nCellsAll
         do k=1,ocnMesh%nVertLevels
            ocnMesh%cellMask(k,n)       = real(cellMaskTmp(k,n),       RKIND)
            ocnMesh%boundaryCell(k,n)   = real(boundaryCellTmp(k,n),   RKIND)
         end do
         end do

         do n=1,ocnMesh%nCellsAll
         do k=1,ocnMesh%maxEdges
            ocnMesh%edgeSignOnCell(k,n) = real(edgeSignOnCellTmp(k,n), RKIND)
         end do
         end do

         do n=1,ocnMesh%nEdgesAll
         do k=1,ocnMesh%nVertLevels
            ocnMesh%edgeMask(k,n)     = real(edgeMaskTmp(k,n),     RKIND)
            ocnMesh%boundaryEdge(k,n) = real(boundaryEdgeTmp(k,n), RKIND)
            ocnMesh%highOrderAdvectionMask(k,n) = & 
                              real(highOrderAdvectionMaskTmp(k,n), RKIND)
         end do
         end do

         do n=1,ocnMesh%nVerticesAll
         do k=1,ocnMesh%nVertLevels
            ocnMesh%vertexMask(k,n)       = real(vertexMaskTmp(k,n),     RKIND)
            ocnMesh%boundaryVertex(k,n)   = real(boundaryVertexTmp(k,n), RKIND)
         end do
         end do

         do n=1,ocnMesh%nVerticesAll
         do k=1,ocnMesh%maxEdges
            ocnMesh%edgeSignOnVertex(k,n) = &
                                          real(edgeSignOnVertexTmp(k,n), RKIND)
         end do
         end do

         ! Copy structure and components to accelerator
         ! Must manually do a deep copy with structure itself copied first.
         !$acc enter data copyin(ocnMesh,                          &
         !$acc                   ocnMesh%nCellsAll,                &
         !$acc                   ocnMesh%nEdgesAll,                &
         !$acc                   ocnMesh%nVerticesAll,             &
         !$acc                   ocnMesh%nCellsOwned,              &
         !$acc                   ocnMesh%nEdgesOwned,              &
         !$acc                   ocnMesh%nVerticesOwned,           &
         !$acc                   ocnMesh%maxEdges,                 &
         !$acc                   ocnMesh%maxEdges2,                &
         !$acc                   ocnMesh%vertexDegree,             &
         !$acc                   ocnMesh%nVertLevels,              &
         !$acc                   ocnMesh%nVertLevelsP1,            &
         !$acc                   ocnMesh%nAdvectionCells,          &
         !$acc                   ocnMesh%nEdgesHalo,               &
         !$acc                   ocnMesh%nCellsHalo,               &
         !$acc                   ocnMesh%nVerticesHalo,            &
         !$acc                   ocnMesh%nEdgesOnEdge,             &
         !$acc                   ocnMesh%nEdgesOnCell,             &
         !$acc                   ocnMesh%maxLevelCell,             &
         !$acc                   ocnMesh%maxLevelEdgeTop,          &
         !$acc                   ocnMesh%maxLevelEdgeBot,          &
         !$acc                   ocnMesh%maxLevelVertexTop,        &
         !$acc                   ocnMesh%maxLevelVertexBot,        &
         !$acc                   ocnMesh%nAdvCellsForEdge,         &
         !$acc                   ocnMesh%indexToCellID,            &
         !$acc                   ocnMesh%indexToEdgeID,            &
         !$acc                   ocnMesh%indexToVertexID,          &
         !$acc                   ocnMesh%edgesOnEdge,              &
         !$acc                   ocnMesh%cellsOnEdge,              &
         !$acc                   ocnMesh%verticesOnEdge,           &
         !$acc                   ocnMesh%cellsOnCell,              &
         !$acc                   ocnMesh%edgesOnCell,              &
         !$acc                   ocnMesh%verticesOnCell,           &
         !$acc                   ocnMesh%cellsOnVertex,            &
         !$acc                   ocnMesh%edgesOnVertex,            &
         !$acc                   ocnMesh%kiteIndexOnCell,          &
         !$acc                   ocnMesh%advCellsForEdge,          &
         !$acc                   ocnMesh%latCell,                  &
         !$acc                   ocnMesh%lonCell,                  &
         !$acc                   ocnMesh%xCell,                    &
         !$acc                   ocnMesh%yCell,                    &
         !$acc                   ocnMesh%zCell,                    &
         !$acc                   ocnMesh%latEdge,                  &
         !$acc                   ocnMesh%lonEdge,                  &
         !$acc                   ocnMesh%xEdge,                    &
         !$acc                   ocnMesh%yEdge,                    &
         !$acc                   ocnMesh%zEdge,                    &
         !$acc                   ocnMesh%latVertex,                &
         !$acc                   ocnMesh%lonVertex,                &
         !$acc                   ocnMesh%xVertex,                  &
         !$acc                   ocnMesh%yVertex,                  &
         !$acc                   ocnMesh%zVertex,                  &
         !$acc                   ocnMesh%fEdge,                    &
         !$acc                   ocnMesh%fVertex,                  &
         !$acc                   ocnMesh%fCell,                    &
         !$acc                   ocnMesh%dcEdge,                   &
         !$acc                   ocnMesh%dvEdge,                   &
         !$acc                   ocnMesh%areaCell,                 &
         !$acc                   ocnMesh%areaTriangle,             &
         !$acc                   ocnMesh%bottomDepth,              &
         !$acc                   ocnMesh%refBottomDepth,           &
         !$acc                   ocnMesh%refBottomDepthTopOfCell,  &
         !$acc                   ocnMesh%vertCoordMovementWeights, &
         !$acc                   ocnMesh%meshScalingDel2,          &
         !$acc                   ocnMesh%meshScalingDel4,          &
         !$acc                   ocnMesh%meshScaling,              &
         !$acc                   ocnMesh%meshDensity,              &
         !$acc                   ocnMesh%angleEdge,                &
         !$acc                   ocnMesh%edgeMask,                 &
         !$acc                   ocnMesh%cellMask,                 &
         !$acc                   ocnMesh%vertexMask,               &
         !$acc                   ocnMesh%boundaryEdge,             &
         !$acc                   ocnMesh%boundaryCell,             &
         !$acc                   ocnMesh%boundaryVertex,           &
         !$acc                   ocnMesh%edgeSignOnCell,           &
         !$acc                   ocnMesh%edgeSignOnVertex,         &
         !$acc                   ocnMesh%highOrderAdvectionMask,   &
         !$acc                   ocnMesh%weightsOnEdge,            &
         !$acc                   ocnMesh%kiteAreasOnVertex,        &
         !$acc                   ocnMesh%edgeTangentVectors,       &
         !$acc                   ocnMesh%edgeNormalVectors,        &
         !$acc                   ocnMesh%localVerticalUnitVectors, &
         !$acc                   ocnMesh%advCoefs,                 &
         !$acc                   ocnMesh%advCoefs3rd,              &
         !$acc                   ocnMesh%derivTwo,                 &
         !$acc                   ocnMesh%cellTangentPlane,         &
         !$acc                   ocnMesh%coeffs_reconstruct)

         block => block % next
      end do

      !$omp end master
      !$omp barrier

!-------------------------------------------------------------------------------

   end subroutine ocn_meshCreate

!*******************************************************************************
!
!  ocn_meshDestroy
!
!> \brief Destroy mesh structure and removes from device
!!  \author Rob Aulwes and Phil Jones
!!  \date   14 Jan 2020
!!  \details
!!  This module removes the mesh variables from the device and invalidates
!!  all pointers in the mesh structure.
!
!-------------------------------------------------------------------------------

   subroutine ocn_meshDestroy(err)

      ! Input variables

      ! Since the ocnMesh is currently a public module variable, no inputs 
      ! here, but eventually may want to treat ocnMesh as a specific 
      ! instantiation instead and pass via args everywhere. If so, need an
      ! input mesh here

      ! Output variables

      integer, intent(out) :: &
         err                   ! returned error flag

      !{{{

      ! Local variables

      !***
      !*** end of preamble, begin code
      !***

      err = 0 

      !$omp barrier
      !$omp master

      ! First remove data from the device. Must remove components first,
      ! then remove the mesh type itself.

      !$acc exit data delete(ocnMesh%nCellsAll,         &
      !$acc                  ocnMesh%nEdgesAll,         &
      !$acc                  ocnMesh%nVerticesAll,      &
      !$acc                  ocnMesh%nCellsOwned,       &
      !$acc                  ocnMesh%nEdgesOwned,       &
      !$acc                  ocnMesh%nVerticesOwned,    &
      !$acc                  ocnMesh%maxEdges,          &
      !$acc                  ocnMesh%maxEdges2,         &
      !$acc                  ocnMesh%vertexDegree,      &
      !$acc                  ocnMesh%nVertLevels,       &
      !$acc                  ocnMesh%nVertLevelsP1,     &
      !$acc                  ocnMesh%nAdvectionCells,   &
      !$acc                  ocnMesh%nEdgesHalo,        &
      !$acc                  ocnMesh%nCellsHalo,        &
      !$acc                  ocnMesh%nVerticesHalo,     &
      !$acc                  ocnMesh%nEdgesOnEdge,      &
      !$acc                  ocnMesh%nEdgesOnCell,      &
      !$acc                  ocnMesh%maxLevelCell,      &
      !$acc                  ocnMesh%maxLevelEdgeTop,   &
      !$acc                  ocnMesh%maxLevelEdgeBot,   &
      !$acc                  ocnMesh%maxLevelVertexTop, &
      !$acc                  ocnMesh%maxLevelVertexBot, &
      !$acc                  ocnMesh%nAdvCellsForEdge,  &
      !$acc                  ocnMesh%indexToCellID,     &
      !$acc                  ocnMesh%indexToEdgeID,     &
      !$acc                  ocnMesh%indexToVertexID,   &
      !$acc                  ocnMesh%edgesOnEdge,       &
      !$acc                  ocnMesh%cellsOnEdge,       &
      !$acc                  ocnMesh%verticesOnEdge,    &
      !$acc                  ocnMesh%cellsOnCell,       &
      !$acc                  ocnMesh%edgesOnCell,       &
      !$acc                  ocnMesh%verticesOnCell,    &
      !$acc                  ocnMesh%cellsOnVertex,     &
      !$acc                  ocnMesh%edgesOnVertex,     &
      !$acc                  ocnMesh%kiteIndexOnCell,   &
      !$acc                  ocnMesh%advCellsForEdge,   &
      !$acc                  ocnMesh%latCell,           &
      !$acc                  ocnMesh%lonCell,           &
      !$acc                  ocnMesh%xCell,             &
      !$acc                  ocnMesh%yCell,             &
      !$acc                  ocnMesh%zCell,             &
      !$acc                  ocnMesh%latEdge,           &
      !$acc                  ocnMesh%lonEdge,           &
      !$acc                  ocnMesh%xEdge,             &
      !$acc                  ocnMesh%yEdge,             &
      !$acc                  ocnMesh%zEdge,             &
      !$acc                  ocnMesh%latVertex,         &
      !$acc                  ocnMesh%lonVertex,         &
      !$acc                  ocnMesh%xVertex,           &
      !$acc                  ocnMesh%yVertex,           &
      !$acc                  ocnMesh%zVertex,           &
      !$acc                  ocnMesh%fEdge,             &
      !$acc                  ocnMesh%fVertex,           &
      !$acc                  ocnMesh%fCell,             &
      !$acc                  ocnMesh%dcEdge,            &
      !$acc                  ocnMesh%dvEdge,            &
      !$acc                  ocnMesh%areaCell,          &
      !$acc                  ocnMesh%areaTriangle,      &
      !$acc                  ocnMesh%bottomDepth,       &
      !$acc                  ocnMesh%refBottomDepth,    &
      !$acc                  ocnMesh%refBottomDepthTopOfCell, &
      !$acc                  ocnMesh%vertCoordMovementWeights, &
      !$acc                  ocnMesh%meshScalingDel2,   &
      !$acc                  ocnMesh%meshScalingDel4,   &
      !$acc                  ocnMesh%meshScaling,       &
      !$acc                  ocnMesh%meshDensity,       &
      !$acc                  ocnMesh%angleEdge,         &
      !$acc                  ocnMesh%edgeMask,          &
      !$acc                  ocnMesh%cellMask,          &
      !$acc                  ocnMesh%vertexMask,        &
      !$acc                  ocnMesh%boundaryEdge,      &
      !$acc                  ocnMesh%boundaryCell,      &
      !$acc                  ocnMesh%boundaryVertex,    &
      !$acc                  ocnMesh%edgeSignOnCell,    &
      !$acc                  ocnMesh%edgeSignOnVertex,  &
      !$acc                  ocnMesh%highOrderAdvectionMask,&
      !$acc                  ocnMesh%weightsOnEdge,     &
      !$acc                  ocnMesh%kiteAreasOnVertex, &
      !$acc                  ocnMesh%edgeTangentVectors,&
      !$acc                  ocnMesh%edgeNormalVectors, &
      !$acc                  ocnMesh%localVerticalUnitVectors, &
      !$acc                  ocnMesh%advCoefs,          &
      !$acc                  ocnMesh%advCoefs3rd,       &
      !$acc                  ocnMesh%derivTwo,          &
      !$acc                  ocnMesh%cellTangentPlane,  &
      !$acc                  ocnMesh%coeffs_reconstruct,&
      !$acc                  ocnMesh)

      ! Reset all scalars to zero
      ocnMesh%nCellsAll       = 0
      ocnMesh%nEdgesAll       = 0
      ocnMesh%nVerticesAll    = 0
      ocnMesh%nCellsOwned     = 0
      ocnMesh%nEdgesOwned     = 0
      ocnMesh%nVerticesOwned  = 0
      ocnMesh%maxEdges        = 0
      ocnMesh%maxEdges2       = 0
      ocnMesh%vertexDegree    = 0
      ocnMesh%nVertLevels     = 0
      ocnMesh%nVertLevelsP1   = 0
      ocnMesh%nAdvectionCells = 0

      ! Now nullify all pointers to invalidate fields
      ! If this becomes the only mesh structure and mesh pool is eliminated,
      !  then we will want to deallocate here instead of nullify.

      nullify(ocnMesh%nEdgesOnEdge,             &
              ocnMesh%nEdgesOnCell,             &
              ocnMesh%maxLevelCell,             &
              ocnMesh%maxLevelEdgeTop,          &
              ocnMesh%maxLevelEdgeBot,          &
              ocnMesh%maxLevelVertexTop,        &
              ocnMesh%maxLevelVertexBot,        &
              ocnMesh%nAdvCellsForEdge,         &
              ocnMesh%indexToCellID,            &
              ocnMesh%indexToEdgeID,            &
              ocnMesh%indexToVertexID,          &
              ocnMesh%edgesOnEdge,              &
              ocnMesh%cellsOnEdge,              &
              ocnMesh%verticesOnEdge,           &
              ocnMesh%cellsOnCell,              &
              ocnMesh%edgesOnCell,              &
              ocnMesh%verticesOnCell,           &
              ocnMesh%cellsOnVertex,            &
              ocnMesh%edgesOnVertex,            &
              ocnMesh%kiteIndexOnCell,          &
              ocnMesh%advCellsForEdge,          &
              ocnMesh%latCell,                  &
              ocnMesh%lonCell,                  &
              ocnMesh%xCell,                    &
              ocnMesh%yCell,                    &
              ocnMesh%zCell,                    &
              ocnMesh%latEdge,                  &
              ocnMesh%lonEdge,                  &
              ocnMesh%xEdge,                    &
              ocnMesh%yEdge,                    &
              ocnMesh%zEdge,                    &
              ocnMesh%latVertex,                &
              ocnMesh%lonVertex,                &
              ocnMesh%xVertex,                  &
              ocnMesh%yVertex,                  &
              ocnMesh%zVertex,                  &
              ocnMesh%fEdge,                    &
              ocnMesh%fVertex,                  &
              ocnMesh%fCell,                    &
              ocnMesh%dcEdge,                   &
              ocnMesh%dvEdge,                   &
              ocnMesh%areaCell,                 &
              ocnMesh%areaTriangle,             &
              ocnMesh%bottomDepth,              &
              ocnMesh%refBottomDepth,           &
              ocnMesh%refBottomDepthTopOfCell,  &
              ocnMesh%vertCoordMovementWeights, &
              ocnMesh%meshScalingDel2,          &
              ocnMesh%meshScalingDel4,          &
              ocnMesh%meshScaling,              &
              ocnMesh%meshDensity,              &
              ocnMesh%angleEdge,                &
              ocnMesh%weightsOnEdge,            &
              ocnMesh%kiteAreasOnVertex,        &
              ocnMesh%edgeTangentVectors,       &
              ocnMesh%edgeNormalVectors,        &
              ocnMesh%localVerticalUnitVectors, &
              ocnMesh%advCoefs,                 &
              ocnMesh%advCoefs3rd,              &
              ocnMesh%derivTwo,                 &
              ocnMesh%cellTangentPlane,         &
              ocnMesh%coeffs_reconstruct)
      
      deallocate(ocnMesh%nEdgesHalo,       &
                 ocnMesh%nCellsHalo,       &
                 ocnMesh%nVerticesHalo,    &
                 ocnMesh%edgeMask,         &
                 ocnMesh%cellMask,         &
                 ocnMesh%vertexMask,       &
                 ocnMesh%boundaryEdge,     &
                 ocnMesh%boundaryCell,     &
                 ocnMesh%boundaryVertex,   &
                 ocnMesh%edgeSignOnCell,   &
                 ocnMesh%edgeSignOnVertex, &
                 ocnMesh%highOrderAdvectionMask)

      !$omp end master

!-------------------------------------------------------------------------------

   end subroutine ocn_meshDestroy

!*******************************************************************************
!
!  ocn_meshUpdateFields
!
!> \brief Updates fields on an accelerator device
!!  \author Rob Aulwes and Phil Jones
!!  \date   14 Jan 2020
!!  \details
!!  Many mesh fields are computed or input later in the initialization 
!!  phase after the meshCreate call. This routine updates these fields
!!  on the device. The routine is not needed if the original mesh pool is
!!  eliminated and this mesh becomes the only mesh structure and can be
!!  updated directly.
!
!-------------------------------------------------------------------------------

   subroutine ocn_meshUpdateFields(domain)

      ! Input arguments

      type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

      !{{{
      ! Local variables

      type(block_type), pointer :: &
         block                    ! variables in current subblock

      type (mpas_pool_type), pointer :: &
         meshPool                 ! mesh variables in MPAS pool structure

      ! scalar pointers for retrieval, but convert to actual scalars in struct
      integer, pointer ::   &! mesh dimensions
         nCellsTmp,         &!
         nEdgesTmp,         &!
         nVerticesTmp,      &!
         maxEdgesTmp,       &!
         maxEdges2Tmp,      &!
         vertexDegreeTmp,   &!
         nVertLevelsTmp,    &!
         nVertLevelsP1Tmp,  &!
         nAdvectionCellsTmp  !

      ! temporary pointers for converting index arrays
      integer, dimension(:), pointer :: &
         nCellsArrayTmp,                &
         nEdgesArrayTmp,                &
         nVerticesArrayTmp

      ! temporary pointers for converting masks
      integer i,k,n          ! loop indices
      integer, dimension(:,:), pointer :: &
         edgeMaskTmp,                     &
         vertexMaskTmp,                   &
         cellMaskTmp,                     &
         highOrderAdvectionMaskTmp,       &
         edgeSignOnCellTmp,               &
         edgeSignOnVertexTmp,             &
         boundaryEdgeTmp,                 &
         boundaryVertexTmp,               &
         boundaryCellTmp

      !***
      !*** end of preamble, begin code
      !***

      !$omp master

      ! already check during create that there should only be one block
      block => domain % blocklist

      ! retrieve the mpas mesh pool
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      ! reset all mesh dimensions, just in case
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
      ocnMesh%maxEdges        = maxEdgesTmp
      ocnMesh%maxEdges2       = maxEdges2Tmp
      ocnMesh%vertexDegree    = vertexDegreeTmp
      ocnMesh%nVertLevels     = nVertLevelsTmp
      ocnMesh%nVertLevelsP1   = nVertLevelsP1Tmp
      ocnMesh%nAdvectionCells = nAdvectionCellsTmp

      ! convert previous index limits into new arrangement
      ocnMesh%nCellsAll       = nCellsTmp
      ocnMesh%nEdgesAll       = nEdgesTmp
      ocnMesh%nVerticesAll    = nVerticesTmp

      n = size(nCellsArrayTmp)
      ocnMesh%nCellsOwned = nCellsArrayTmp(1)
      do i=2,n
         ocnMesh%nCellsHalo(i-1) = nCellsArrayTmp(i) 
      end do

      n = size(nEdgesArrayTmp)
      ocnMesh%nEdgesOwned = nEdgesArrayTmp(1)
      do i=2,n
         ocnMesh%nEdgesHalo(i-1) = nEdgesArrayTmp(i) 
      end do

      n = size(nVerticesArrayTmp)
      ocnMesh%nVerticesOwned = nVerticesArrayTmp(1)
      do i=2,n
         ocnMesh%nVerticesHalo(i-1) = nVerticesArrayTmp(i) 
      end do

      ! Because these are pointers into meshPool, they do not need
      ! to be updated - they capture the updates automatically. 
      !call mpas_pool_get_array(meshPool, 'nEdgesOnEdge',      &
      !                            ocnMesh%nEdgesOnEdge)
      !call mpas_pool_get_array(meshPool, 'nEdgesOnCell',      &
      !                            ocnMesh%nEdgesOnCell)
      !call mpas_pool_get_array(meshPool, 'maxLevelCell',      &
      !                            ocnMesh%maxLevelCell)
      !call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop',   &
      !                            ocnMesh%maxLevelEdgeTop)
      !call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot',   &
      !                            ocnMesh%maxLevelEdgeBot)
      !call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', &
      !                            ocnMesh%maxLevelVertexTop)
      !call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', &
      !                            ocnMesh%maxLevelVertexBot)
      !call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge',  &
      !                            ocnMesh%nAdvCellsForEdge)
      !call mpas_pool_get_array(meshPool, 'indexToCellID',     &
      !                            ocnMesh%indexToCellID)
      !call mpas_pool_get_array(meshPool, 'indexToEdgeID',     &
      !                            ocnMesh%indexToEdgeID)
      !call mpas_pool_get_array(meshPool, 'indexToVertexID',   &
      !                            ocnMesh%indexToVertexID)
      !call mpas_pool_get_array(meshPool, 'edgesOnEdge',       &
      !                            ocnMesh%edgesOnEdge)
      !call mpas_pool_get_array(meshPool, 'cellsOnEdge',       &
      !                            ocnMesh%cellsOnEdge)
      !call mpas_pool_get_array(meshPool, 'verticesOnEdge',    &
      !                            ocnMesh%verticesOnEdge)
      !call mpas_pool_get_array(meshPool, 'cellsOnCell',       &
      !                            ocnMesh%cellsOnCell)
      !call mpas_pool_get_array(meshPool, 'edgesOnCell',       &
      !                            ocnMesh%edgesOnCell)
      !call mpas_pool_get_array(meshPool, 'verticesOnCell',    &
      !                            ocnMesh%verticesOnCell)
      !call mpas_pool_get_array(meshPool, 'cellsOnVertex',     &
      !                            ocnMesh%cellsOnVertex)
      !call mpas_pool_get_array(meshPool, 'edgesOnVertex',     &
      !                            ocnMesh%edgesOnVertex)
      !call mpas_pool_get_array(meshPool, 'kiteIndexOnCell',   &
      !                            ocnMesh%kiteIndexOnCell)
      !call mpas_pool_get_array(meshPool, 'advCellsForEdge',   &
      !                            ocnMesh%advCellsForEdge)

      ! these are also pointers that do not require updating
      ! now set a number of physics and numerical properties of mesh
      !call mpas_pool_get_array(meshPool, 'latCell',           &
      !                            ocnMesh%latCell)
      !call mpas_pool_get_array(meshPool, 'lonCell',           &
      !                            ocnMesh%lonCell)
      !call mpas_pool_get_array(meshPool, 'xCell',             &
      !                            ocnMesh%xCell)
      !call mpas_pool_get_array(meshPool, 'yCell',             &
      !                            ocnMesh%yCell)
      !call mpas_pool_get_array(meshPool, 'zCell',             &
      !                            ocnMesh%zCell)
      !call mpas_pool_get_array(meshPool, 'latEdge',           &
      !                            ocnMesh%latEdge)
      !call mpas_pool_get_array(meshPool, 'lonEdge',           &
      !                            ocnMesh%lonEdge)
      !call mpas_pool_get_array(meshPool, 'xEdge',             &
      !                            ocnMesh%xEdge)
      !call mpas_pool_get_array(meshPool, 'yEdge',             &
      !                            ocnMesh%yEdge)
      !call mpas_pool_get_array(meshPool, 'zEdge',             &
      !                            ocnMesh%zEdge)
      !call mpas_pool_get_array(meshPool, 'latVertex',         &
      !                            ocnMesh%latVertex)
      !call mpas_pool_get_array(meshPool, 'lonVertex',         &
      !                            ocnMesh%lonVertex)
      !call mpas_pool_get_array(meshPool, 'xVertex',           &
      !                            ocnMesh%xVertex)
      !call mpas_pool_get_array(meshPool, 'yVertex',           &
      !                            ocnMesh%yVertex)
      !call mpas_pool_get_array(meshPool, 'zVertex',           &
      !                            ocnMesh%zVertex)
      !call mpas_pool_get_array(meshPool, 'fEdge',             &
      !                            ocnMesh%fEdge)
      !call mpas_pool_get_array(meshPool, 'fVertex',           &
      !                            ocnMesh%fVertex)
      !call mpas_pool_get_array(meshPool, 'fCell',             &
      !                            ocnMesh%fCell)
      !call mpas_pool_get_array(meshPool, 'dcEdge',            &
      !                            ocnMesh%dcEdge)
      !call mpas_pool_get_array(meshPool, 'dvEdge',            &
      !                            ocnMesh%dvEdge)
      !call mpas_pool_get_array(meshPool, 'areaCell',          &
      !                            ocnMesh%areaCell)
      !call mpas_pool_get_array(meshPool, 'areaTriangle',      &
      !                            ocnMesh%areaTriangle)
      !call mpas_pool_get_array(meshPool, 'weightsOnEdge',     &
      !                            ocnMesh%weightsOnEdge)
      !call mpas_pool_get_array(meshPool, 'bottomDepth',       &
      !                            ocnMesh%bottomDepth)
      !call mpas_pool_get_array(meshPool, 'refBottomDepth',    &
      !                            ocnMesh%refBottomDepth)
      !call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell',   &
      !                            ocnMesh%refBottomDepthTopOfCell)
      !call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights',  &
      !                            ocnMesh%vertCoordMovementWeights)
      !call mpas_pool_get_array(meshPool, 'meshScalingDel2',   &
      !                            ocnMesh%meshScalingDel2)
      !call mpas_pool_get_array(meshPool, 'meshScalingDel4',   &
      !                            ocnMesh%meshScalingDel4)
      !call mpas_pool_get_array(meshPool, 'meshScaling',       &
      !                            ocnMesh%meshScaling)
      !call mpas_pool_get_array(meshPool, 'meshDensity',       &
      !                            ocnMesh%meshDensity)
      !call mpas_pool_get_array(meshPool, 'angleEdge',         &
      !                            ocnMesh%angleEdge)
      !  
      !call mpas_pool_get_array(meshPool, 'weightsOnEdge',             &
      !                            ocnMesh%weightsOnEdge)
      !call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex',         &
      !                            ocnMesh%kiteAreasOnVertex)
      !call mpas_pool_get_array(meshPool, 'edgeTangentVectors',        &
      !                            ocnMesh%edgeTangentVectors)
      !call mpas_pool_get_array(meshPool, 'edgeNormalVectors',         &
      !                            ocnMesh%edgeNormalVectors)
      !call mpas_pool_get_array(meshPool, 'localVerticalUnitVectors',  &
      !                            ocnMesh%localVerticalUnitVectors)
      !call mpas_pool_get_array(meshPool, 'advCoefs',                  &
      !                            ocnMesh%advCoefs)
      !call mpas_pool_get_array(meshPool, 'advCoefs3rd',               &
      !                            ocnMesh%advCoefs3rd)
      !call mpas_pool_get_array(meshPool, 'derivTwo',                  &
      !                            ocnMesh%derivTwo)
      !call mpas_pool_get_array(meshPool, 'cellTangentPlane',          &
      !                            ocnMesh%cellTangentPlane)
      !call mpas_pool_get_array(meshPool, 'coeffs_reconstruct',        &
      !                            ocnMesh%coeffs_reconstruct)

      ! For masks, we converted to real masks, so need to recompute for
      ! updated values.
      call mpas_pool_get_array(meshPool, 'edgeMask',        &
                                          edgeMaskTmp)
      call mpas_pool_get_array(meshPool, 'vertexMask',        &
                                          vertexMaskTmp)
      call mpas_pool_get_array(meshPool, 'cellMask',        &
                                          cellMaskTmp)
      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask',        &
                                          highOrderAdvectionMaskTmp)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell',  &
                                          edgeSignOnCellTmp)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex',  &
                                          edgeSignOnVertexTmp)
      call mpas_pool_get_array(meshPool, 'boundaryEdge',  &
                                          boundaryEdgeTmp)
      call mpas_pool_get_array(meshPool, 'boundaryVertex',  &
                                          boundaryVertexTmp)
      call mpas_pool_get_array(meshPool, 'boundaryCell',  &
                                          boundaryCellTmp)

      do n=1,ocnMesh%nCellsAll
      do k=1,ocnMesh%nVertLevels
         ocnMesh%cellMask(k,n)       = real(cellMaskTmp(k,n),       RKIND)
         ocnMesh%boundaryCell(k,n)   = real(boundaryCellTmp(k,n),   RKIND)
      end do
      end do

      do n=1,ocnMesh%nCellsAll
      do k=1,ocnMesh%maxEdges
         ocnMesh%edgeSignOnCell(k,n) = real(edgeSignOnCellTmp(k,n), RKIND)
      end do
      end do

      do n=1,ocnMesh%nEdgesAll
      do k=1,ocnMesh%nVertLevels
         ocnMesh%edgeMask(k,n)     = real(edgeMaskTmp(k,n),     RKIND)
         ocnMesh%boundaryEdge(k,n) = real(boundaryEdgeTmp(k,n), RKIND)
         ocnMesh%highOrderAdvectionMask(k,n) = & 
                           real(highOrderAdvectionMaskTmp(k,n), RKIND)
      end do
      end do

      do n=1,ocnMesh%nVerticesAll
      do k=1,ocnMesh%nVertLevels
         ocnMesh%vertexMask(k,n)       = real(vertexMaskTmp(k,n),     RKIND)
         ocnMesh%boundaryVertex(k,n)   = real(boundaryVertexTmp(k,n), RKIND)
      end do
      end do

      do n=1,ocnMesh%nVerticesAll
      do k=1,ocnMesh%maxEdges
         ocnMesh%edgeSignOnVertex(k,n) = &
                                       real(edgeSignOnVertexTmp(k,n), RKIND)
      end do
      end do

      ! Go ahead and update all fields on device to be safe.
      ! NOTE: if we end up computing some fields on the device during
      !       init, the update must go the opposite direction (update host)

      !$acc update device(ocnMesh%nCellsAll,                &
      !$acc               ocnMesh%nEdgesAll,                &
      !$acc               ocnMesh%nVerticesAll,             &
      !$acc               ocnMesh%nCellsOwned,              &
      !$acc               ocnMesh%nEdgesOwned,              &
      !$acc               ocnMesh%nVerticesOwned,           &
      !$acc               ocnMesh%maxEdges,                 &
      !$acc               ocnMesh%maxEdges2,                &
      !$acc               ocnMesh%vertexDegree,             &
      !$acc               ocnMesh%nVertLevels,              &
      !$acc               ocnMesh%nVertLevelsP1,            &
      !$acc               ocnMesh%nAdvectionCells,          &
      !$acc               ocnMesh%nEdgesHalo,               &
      !$acc               ocnMesh%nCellsHalo,               &
      !$acc               ocnMesh%nVerticesHalo,            &
      !$acc               ocnMesh%nEdgesOnEdge,             &
      !$acc               ocnMesh%nEdgesOnCell,             &
      !$acc               ocnMesh%maxLevelCell,             &
      !$acc               ocnMesh%maxLevelEdgeTop,          &
      !$acc               ocnMesh%maxLevelEdgeBot,          &
      !$acc               ocnMesh%maxLevelVertexTop,        &
      !$acc               ocnMesh%maxLevelVertexBot,        &
      !$acc               ocnMesh%nAdvCellsForEdge,         &
      !$acc               ocnMesh%indexToCellID,            &
      !$acc               ocnMesh%indexToEdgeID,            &
      !$acc               ocnMesh%indexToVertexID,          &
      !$acc               ocnMesh%edgesOnEdge,              &
      !$acc               ocnMesh%cellsOnEdge,              &
      !$acc               ocnMesh%verticesOnEdge,           &
      !$acc               ocnMesh%cellsOnCell,              &
      !$acc               ocnMesh%edgesOnCell,              &
      !$acc               ocnMesh%verticesOnCell,           &
      !$acc               ocnMesh%cellsOnVertex,            &
      !$acc               ocnMesh%edgesOnVertex,            &
      !$acc               ocnMesh%kiteIndexOnCell,          &
      !$acc               ocnMesh%advCellsForEdge,          &
      !$acc               ocnMesh%latCell,                  &
      !$acc               ocnMesh%lonCell,                  &
      !$acc               ocnMesh%xCell,                    &
      !$acc               ocnMesh%yCell,                    &
      !$acc               ocnMesh%zCell,                    &
      !$acc               ocnMesh%latEdge,                  &
      !$acc               ocnMesh%lonEdge,                  &
      !$acc               ocnMesh%xEdge,                    &
      !$acc               ocnMesh%yEdge,                    &
      !$acc               ocnMesh%zEdge,                    &
      !$acc               ocnMesh%latVertex,                &
      !$acc               ocnMesh%lonVertex,                &
      !$acc               ocnMesh%xVertex,                  &
      !$acc               ocnMesh%yVertex,                  &
      !$acc               ocnMesh%zVertex,                  &
      !$acc               ocnMesh%fEdge,                    &
      !$acc               ocnMesh%fVertex,                  &
      !$acc               ocnMesh%fCell,                    &
      !$acc               ocnMesh%dcEdge,                   &
      !$acc               ocnMesh%dvEdge,                   &
      !$acc               ocnMesh%areaCell,                 &
      !$acc               ocnMesh%areaTriangle,             &
      !$acc               ocnMesh%bottomDepth,              &
      !$acc               ocnMesh%refBottomDepth,           &
      !$acc               ocnMesh%refBottomDepthTopOfCell,  &
      !$acc               ocnMesh%vertCoordMovementWeights, &
      !$acc               ocnMesh%meshScalingDel2,          &
      !$acc               ocnMesh%meshScalingDel4,          &
      !$acc               ocnMesh%meshScaling,              &
      !$acc               ocnMesh%meshDensity,              &
      !$acc               ocnMesh%angleEdge,                &
      !$acc               ocnMesh%edgeMask,                 &
      !$acc               ocnMesh%cellMask,                 &
      !$acc               ocnMesh%vertexMask,               &
      !$acc               ocnMesh%boundaryEdge,             &
      !$acc               ocnMesh%boundaryCell,             &
      !$acc               ocnMesh%boundaryVertex,           &
      !$acc               ocnMesh%edgeSignOnCell,           &
      !$acc               ocnMesh%edgeSignOnVertex,         &
      !$acc               ocnMesh%highOrderAdvectionMask,   &
      !$acc               ocnMesh%weightsOnEdge,            &
      !$acc               ocnMesh%kiteAreasOnVertex,        &
      !$acc               ocnMesh%edgeTangentVectors,       &
      !$acc               ocnMesh%edgeNormalVectors,        &
      !$acc               ocnMesh%localVerticalUnitVectors, &
      !$acc               ocnMesh%advCoefs,                 &
      !$acc               ocnMesh%advCoefs3rd,              &
      !$acc               ocnMesh%derivTwo,                 &
      !$acc               ocnMesh%cellTangentPlane,         &
      !$acc               ocnMesh%coeffs_reconstruct)

      !$omp end master
      !$omp barrier

!-------------------------------------------------------------------------------

   end subroutine ocn_meshUpdateFields

!*******************************************************************************
   
end module ocn_mesh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
