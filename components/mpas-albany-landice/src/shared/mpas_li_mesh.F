!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! \file mpas_li_mesh.F
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
!  li_mesh
!
!>  \brief MALI mesh structure with GPU support
!> \author Rob Aulwes and Phil Jones, modified for MALI
!>         by Trevor Hillebrand and Matthew Hoffman
!> \date   14 Jan 2020, added to MALI 2023
!> \details
!> This module creates and maintains a primary land ice mesh structure
!> and ensures all mesh variables are copied to an accelerator device
!> if needed. Currently it consists of pointers to the existing MPAS mesh pool
!> variables, but is intended to eventually replace the mesh pool later.
!
!-------------------------------------------------------------------------------

module li_mesh

   ! module dependencies
   use mpas_dmpar
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_log

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
      nCells, &! total number of local (owned+halo) cells in primary
      nEdges, &! total number of local edge midpoints
      nVertices, &! total number of local cells in dual (cell vertices)
      nCellsSolve, &! number of cells    owned by the local domain
      nEdgesSolve, &! number of edges    owned by the local domain
      nVerticesSolve, &! number of vertices owned by the local domain
      maxEdges, &! largest number of edges any polygon has
      maxEdges2, &! 2x the largest number of edges any polygon has
      vertexDegree, &! number of cells or edges touching each vertex
      nVertLevels ! number of vertical levels

   integer, public, dimension(:), allocatable :: &
      nCellsHalo, &! number of owned+halo(n) cells in local domain
      nEdgesHalo, &! number of owned+halo(n) edges in local domain
      nVerticesHalo      ! number of owned+halo(n) vertices in local domain

   integer, public, dimension(:), pointer :: &
      nEdgesOnEdge, &! number of edges connected to each edge point
      nEdgesOnCell, & ! number of edges associated with each cell center
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
      edgesOnVertex, & ! index of edges connected to each vertex
      edgeSignOnCell, &! sign of edge contributions to a cell
      edgeSignOnVertex ! sign of edge contributions to a vertex

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
      dcEdge, &! length of edge = dist between cells across edge
      dvEdge, &! length of edge = dist between vertices along edge
      areaCell, &! area of each cell
      areaTriangle, &! area of each cell on dual grid
      meshDensity, &! density of mesh
      angleEdge ! angle the edge normal makes with local east

   real(kind=RKIND), public, dimension(:), allocatable :: &
      invAreaCell   ! inverse of area of each cell

   ! Multiplicative masks and vectors for various conditions
   integer, public, dimension(:), allocatable :: &
      boundaryCell ! mask for boundary cells    at each level

   real(kind=RKIND), public, dimension(:, :), pointer :: &
      weightsOnEdge, &! weights on each edge
      kiteAreasOnVertex, &! real (vertexDegree nVertices)
      edgeNormalVectors, &! normal  unit vector at edge
      localVerticalUnitVectors ! local unit vector iin vertical

   real(kind=RKIND), public, dimension(:, :, :), pointer :: &
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
      li_meshCreate, &
      li_meshDestroy
   !}}}

!***********************************************************************

contains

!***********************************************************************
!
!  li_meshCreate
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

   subroutine li_meshCreate(domain) !{{{

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
         nVertLevelsTmp !
         !nVertLevelsP1Tmp !

      ! temporary pointers for converting index arrays
      integer, dimension(:), pointer :: &
         nCellsArrayTmp, &
         nEdgesArrayTmp, &
         nVerticesArrayTmp

      ! temporary pointers for converting masks
      integer i, k, n          ! loop indices
      integer, dimension(:, :), pointer :: &
         edgeSignOnCellTmp, &
         edgeSignOnVertexTmp

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
               'li_meshCreate: more than one block no longer supported', &
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

      ! convert previous index limits into new halo definitions
      nCells = nCellsTmp
      nEdges = nEdgesTmp
      nVertices = nVerticesTmp

      n = size(nCellsArrayTmp)
      allocate (nCellsHalo(n - 1))
      nCellsSolve = nCellsArrayTmp(1)
      do i = 2, n
         nCellsHalo(i - 1) = nCellsArrayTmp(i)
      end do

      n = size(nEdgesArrayTmp)
      allocate (nEdgesHalo(n - 1))
      nEdgesSolve = nEdgesArrayTmp(1)
      do i = 2, n
         nEdgesHalo(i - 1) = nEdgesArrayTmp(i)
      end do

      n = size(nVerticesArrayTmp)
      allocate (nVerticesHalo(n - 1))
      nVerticesSolve = nVerticesArrayTmp(1)
      do i = 2, n
         nVerticesHalo(i - 1) = nVerticesArrayTmp(i)
      end do

      ! set pointers for a lot of connectivity info
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', &
                               nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', &
                               nEdgesOnCell)
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
      call mpas_pool_get_array(meshPool, 'meshDensity', &
                               meshDensity)
      call mpas_pool_get_array(meshPool, 'angleEdge', &
                               angleEdge)
      call mpas_pool_get_array(meshPool, 'weightsOnEdge', &
                               weightsOnEdge)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', &
                               kiteAreasOnVertex)
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
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', &
                               edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', &
                               edgeSignOnVertex)

      allocate ( &
            boundaryCell(nCells+1), &
            invAreaCell(nCells+1))

      !-----------------------------------------------------------------
      ! Now that all pointers are set and mesh variables allocated
      ! we initialize other mesh quantities
      !-----------------------------------------------------------------

      ! Start by initializing vertical mesh, min/max cells and
      ! sign/index fields
      call meshSignIndexFields()

      areaCell(nCells+1) = -1.0e34_RKIND

      ! Compute the inverse of areaCell
      do n = 1, nCells
         invAreaCell(n) = 1.0_RKIND / areaCell(n)
      end do
      invAreaCell(nCells+1) = 0.0_RKIND

!-------------------------------------------------------------------------------

   end subroutine li_meshCreate !}}}

!*******************************************************************************
!
!  li_meshDestroy
!
!> \brief Destroy mesh structure and removes from device
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> This module removes the mesh variables from the device and invalidates
!> all pointers in the mesh structure.
!
!-------------------------------------------------------------------------------

   subroutine li_meshDestroy(err) !{{{

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

      ! Reset all scalars to zero
      onSphere  = .false.
      sphereRadius = 0.0_RKIND
      nCells = 0
      nEdges = 0
      nVertices = 0
      nCellsSolve = 0
      nEdgesSolve = 0
      nVerticesSolve = 0
      maxEdges = 0
      maxEdges2 = 0
      vertexDegree = 0
      nVertLevels = 0

      ! Now nullify all pointers to invalidate fields
      ! If this becomes the only mesh structure and mesh pool is eliminated,
      !  then we will want to deallocate here instead of nullify.

      nullify (nEdgesOnEdge, &
               nEdgesOnCell, &
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
               dcEdge, &
               dvEdge, &
               areaCell, &
               areaTriangle, &
               meshDensity, &
               angleEdge, &
               weightsOnEdge, &
               kiteAreasOnVertex, &
               edgeNormalVectors, &
               localVerticalUnitVectors, &
               cellTangentPlane, &
               edgeSignOnCell, &
               edgeSignOnVertex, &
               coeffs_reconstruct)

      deallocate (nEdgesHalo, &
                  nCellsHalo, &
                  nVerticesHalo, &
                  boundaryCell, &
                  invAreaCell)

!-------------------------------------------------------------------------------

   end subroutine li_meshDestroy !}}}

!***********************************************************************
!
!  routine li_meshSignIndexFields
!
!> \brief   set up sign and index fields
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes edgeSignOnCell, edgeSignOnVertex, and
!>  kiteIndexOnCell.
!
!-----------------------------------------------------------------------

   subroutine meshSignIndexFields()!{{{

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

      do iCell = 1, nCells
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            iVertex = verticesOnCell(i, iCell)

            ! Vector points from cell 1 to cell 2
            if (iCell == cellsOnEdge(1, iEdge)) then
               edgeSignOnCell(i, iCell) = -1.0_RKIND
            else
               edgeSignOnCell(i, iCell) =  1.0_RKIND
            end if

         end do
      end do

      do iVertex = 1, nVertices
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

   end subroutine meshSignIndexFields !}}}
!***********************************************************************

end module li_mesh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
