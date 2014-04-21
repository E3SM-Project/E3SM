module cice_mesh

  use mpas_grid_types

  implicit none

  private
  public :: &
       init_boundary, &
       calc_cell_vertices_at_vertex, &
       normal_vectors, &
       dot_product_3space, &
       grid_rotation_forward, &
       grid_rotation_forward_nochange, &
       latlon_vector_rotation_forward, &
       latlon_vector_rotation_backward, &
       latlon_from_xyz, &
       interpolate_cell_to_vertex

contains

  !-------------------------------------------------------------
  
  subroutine init_boundary(block)

    type (block_type), intent(inout) :: block

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary

    integer, dimension(:), pointer :: &
         interiorVertex

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call interior_vertices(mesh, interiorVertex)

  end subroutine init_boundary

  !-------------------------------------------------------------

  subroutine interior_vertices(mesh, interiorVertex)

    type (MPAS_pool_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         interiorVertex

    integer :: &
         nInteriorAdjacentCells

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    integer, pointer :: &
         nVertices, &
         vertexDegree, &
         nCells

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)

    ! boundary vertices
    do iVertex = 1, nVertices

       interiorVertex(iVertex) = 0
       
       nInteriorAdjacentCells = 0

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree, iVertex)

          if (iCell >= 1 .and. iCell <= nCells) then
             
             nInteriorAdjacentCells = nInteriorAdjacentCells + 1
             
          endif
          
       enddo ! iVertexDegree

       ! vertex points we directly calculate on
       if (nInteriorAdjacentCells == vertexDegree) then

          interiorVertex(iVertex) = 1
          
       endif
       
    enddo ! iVertex

  end subroutine interior_vertices

  !-------------------------------------------------------------

  subroutine boundary_cells(mesh, boundaryCell, boundaryCell2, interiorCell) 

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         boundaryCell, &
         boundaryCell2, &
         interiorCell

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnCell
         
    integer :: &
         iCell, &
         iCellOnCell, &
         iAdjacentCell, &
         nAdjacentCells

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnCell", cellsOnCell)

    do iCell = 1, nCells

       boundaryCell(iCell) = 0
       interiorCell(iCell) = 0

       nAdjacentCells = 0

       do iCellOnCell = 1, nEdgesOnCell(iCell)

          iAdjacentCell = cellsOnCell(iCellOnCell,iCell)

          if (iAdjacentCell <= nCells) then

             nAdjacentCells = nAdjacentCells + 1

          endif

          if (nAdjacentCells == nEdgesOnCell(iCell)) interiorCell(iCell) = 1
          if (nAdjacentCells <  nEdgesOnCell(iCell)) boundaryCell(iCell) = 1

       enddo ! iCellOnCell

    enddo ! iCell

    do iCell = 1, nCells

       boundaryCell2(iCell) = 0

       if (boundaryCell(iCell) == 0) then

          do iCellOnCell = 1, nEdgesOnCell(iCell)

             iAdjacentCell = cellsOnCell(iCellOnCell,iCell)

             if (boundaryCell(iAdjacentCell) == 1) then

                boundaryCell2(iCell) = 1

             endif

          enddo ! iCellOnCell

       endif

    enddo ! iCell

  end subroutine boundary_cells

  !-------------------------------------------------------------

  subroutine boundary_vertices(mesh, boundaryVertex, interiorVertex, interiorVertex2, boundaryCell)

    type(MPAS_pool_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         boundaryVertex, & ! boundary vertices of boundary vertex triangles
         interiorVertex, & ! all normal shaped interior vertices
         interiorVertex2

    integer, dimension(:), intent(in) :: &
         boundaryCell

    integer :: &
         nInteriorAdjacentCells

    integer, pointer :: &
         nVertices, &
         nCells, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    logical :: &
         linterior

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)

    ! boundary vertices
    do iVertex = 1, nVertices

       boundaryVertex(iVertex) = 0
       interiorVertex(iVertex) = 0
       
       nInteriorAdjacentCells = 0

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree, iVertex)
          
          if (iCell >= 1 .and. iCell <= nCells) then
             
             nInteriorAdjacentCells = nInteriorAdjacentCells + 1
             
          endif
          
       enddo ! iVertexDegree

       ! boundary vertices that have triangles
       if (nInteriorAdjacentCells == 2) then

          boundaryVertex(iVertex) = 1
          
       endif

       ! vertex points we directly calculate on
       if (nInteriorAdjacentCells == 3) then

          interiorVertex(iVertex) = 1
          
       endif
       
    enddo ! iVertex


    ! second interior vertices
    do iVertex = 1, nVertices

       interiorVertex2(iVertex) = 0

       linterior = .true.

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree, iVertex)
          
          if (boundaryCell(iCell) == 1) then
             
             linterior = .false.
             
          endif
          
       enddo ! iVertexDegree

       if (linterior) interiorVertex2(iVertex) = 1

    enddo ! iVertex

  end subroutine boundary_vertices

  !-------------------------------------------------------------

  subroutine boundary_edges(mesh, boundaryEdge)

    type(MPAS_pool_type), intent(in) :: mesh

    integer, dimension(:), intent(out) :: &
         boundaryEdge

    integer :: &
         nInteriorAdjacentCells

    integer, pointer :: &
         nEdges, &
         nCells

    integer, dimension(:,:), pointer :: &
         cellsOnEdge

    integer :: &
         iCell, &
         iEdge, &
         iCellsOnEdge
    
    ! init variables
    call MPAS_pool_get_array(mesh, "nEdges", nEdges)
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)

    ! boundary edges
    do iEdge = 1, nEdges
       
       do iCellsOnEdge = 1, 2
          
          iCell = cellsOnEdge(iCellsOnEdge,iEdge)
          
          if (iCell > nCells) then
             boundaryEdge(iEdge) = 1
          else
             boundaryEdge(iEdge) = 0
          endif
          
       enddo ! iCellsOnEdge

    enddo ! iEdge
    
  end subroutine boundary_edges

  !-------------------------------------------------------------
  ! mesh searches
  !-------------------------------------------------------------

   subroutine calc_cell_vertices_at_vertex(mesh, cellVerticesAtVertex)

    type(MPAS_pool_type), intent(in) :: mesh

    integer, dimension(:,:), intent(out) :: &
         cellVerticesAtVertex

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex, &
         verticesOnCell

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         iVertexOnCell, &
         jVertex
    
    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)

    do iVertex = 1, nVertices

       do iVertexDegree = 1, vertexDegree

          cellVerticesAtVertex(iVertexDegree,iVertex) = 0

          iCell = cellsOnVertex(iVertexDegree, iVertex)

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             jVertex = verticesOnCell(iVertexOnCell,iCell)

             if (iVertex == jVertex) then

                cellVerticesAtVertex(iVertexDegree,iVertex) = iVertexOnCell

             endif

          enddo ! iVertexOnCell

       enddo ! iVertexDegree

    enddo ! iVertex

   end subroutine calc_cell_vertices_at_vertex


  !-------------------------------------------------------------
  ! normal vectors
  !-------------------------------------------------------------

  subroutine normal_vectors(mesh, &
                            normalVectorPolygon, normalVectorTriangle, &
                            latCellRotated,      latVertexRotated,     &
                            interiorVertex,      gridCurvature)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon, &
         normalVectorTriangle

    real(kind=RKIND), dimension(:), intent(out) :: &
         latCellRotated, &
         latVertexRotated

    integer, dimension(:), intent(in) :: &
         interiorVertex

    character(len=*), intent(in) :: &
         gridCurvature

    if (trim(gridCurvature) == "planar") then

       call normal_vectors_planar_polygon(mesh, normalVectorPolygon)
       call normal_vectors_planar_triangle(mesh, normalVectorTriangle, interiorVertex)

    else if (trim(gridCurvature) == "spherical") then

       call normal_vectors_spherical_polygon_metric(mesh, normalVectorPolygon, latCellRotated)
       call normal_vectors_spherical_triangle_metric(mesh, normalVectorTriangle, latVertexRotated, interiorVertex)

    endif

  end subroutine normal_vectors

  !-------------------------------------------------------------

  subroutine normal_vectors_planar_polygon(mesh, normalVectorPolygon)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xEdge, &
         yEdge, &
         xCell, &
         yCell

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge

    real(kind=RKIND) :: &
         dx, dy

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    do iCell = 1, nCells

       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          
          iEdge = edgesOnCell(iEdgeOnCell,iCell)

          dx = xEdge(iEdge) - xCell(iCell)
          dy = yEdge(iEdge) - yCell(iCell)

          normalVectorPolygon(1,iEdgeOnCell,iCell) = dx / sqrt(dx**2 + dy**2)
          normalVectorPolygon(2,iEdgeOnCell,iCell) = dy / sqrt(dx**2 + dy**2)

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine normal_vectors_planar_polygon

  !-------------------------------------------------------------

  subroutine normal_vectors_planar_triangle(mesh, normalVectorTriangle, interiorVertex)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         interiorVertex

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         edgesOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         xEdge, &
         yEdge, &
         xVertex, &
         yVertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge

    real(kind=RKIND) :: &
         dx, dy

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)

    do iVertex = 1, nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          do iVertexDegree = 1, vertexDegree
             
             iEdge = edgesOnVertex(iVertexDegree,iVertex)
             
             dx = xEdge(iEdge) - xVertex(iVertex)
             dy = yEdge(iEdge) - yVertex(iVertex)
             
             normalVectorTriangle(1,iVertexDegree,iVertex) = dx / sqrt(dx**2 + dy**2)
             normalVectorTriangle(2,iVertexDegree,iVertex) = dy / sqrt(dx**2 + dy**2)
             
          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_planar_triangle

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon_metric(mesh, normalVectorPolygon, latCellRotated)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon

    real(kind=RKIND), dimension(:), intent(out) :: &
         latCellRotated
    
    real(kind=RKIND), dimension(3) :: &
         cellCentreRotated, &
         cellCentreEquator, &
         edgeRotated, &
         vertexRotated1, &
         vertexRotated2, &
         edgeEquator, &
         vertexEquator1, &
         vertexEquator2, &
         vertexVector, &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         lonCellRotated, &
         normalGreatCircleNorm

    real(kind=RKIND), dimension(3,3) :: &
         yRotationMatrix, &
         zRotationMatrix

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell, &
         verticesOnEdge, &
         cellsOnEdge

    real(kind=RKIND), pointer :: &
         sphere_radius
    
    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell, &
         xEdge, &
         yEdge, &
         zEdge, &
         xVertex, &
         yVertex, &
         zVertex

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertex1, &
         iVertex2

    real(kind=RKIND), parameter :: &
         rad_to_deg = 180.0_RKIND / 3.14159265359_RKIND

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "zEdge", zEdge)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    yRotationMatrix = 0.0_RKIND
    zRotationMatrix = 0.0_RKIND

    yRotationMatrix(2,2) = 1.0_RKIND
    zRotationMatrix(3,3) = 1.0_RKIND

    do iCell = 1, nCells

       ! rotate the cell centre to rotated geographical grid
       call grid_rotation_forward(cellCentreRotated(1), cellCentreRotated(2), cellCentreRotated(3), &
                                  xCell(iCell),         yCell(iCell),         zCell(iCell))

       ! calculate lon and lat of the cell centre in the rotated coordinate system
       lonCellRotated        = atan2(cellCentreRotated(2), cellCentreRotated(1))
       latCellRotated(iCell) = asin(cellCentreRotated(3) / sphere_radius)

       ! rotate the cell centre onto the rotated equator
       !yRotationMatrix(1,1) =  cos(latCellRotated(iCell))
       !yRotationMatrix(1,3) =  sin(latCellRotated(iCell))
       !yRotationMatrix(3,1) = -sin(latCellRotated(iCell))
       !yRotationMatrix(3,3) =  cos(latCellRotated(iCell))

       !zRotationMatrix(1,1) =  cos(-lonCellRotated)
       !zRotationMatrix(1,2) = -sin(-lonCellRotated)
       !zRotationMatrix(2,1) =  sin(-lonCellRotated)
       !zRotationMatrix(2,2) =  cos(-lonCellRotated)

       yRotationMatrix(1,1) = 1.0_RKIND
       yRotationMatrix(1,3) = 0.0_RKIND
       yRotationMatrix(3,1) = 0.0_RKIND
       yRotationMatrix(3,3) = 1.0_RKIND

       zRotationMatrix(1,1) = 1.0_RKIND
       zRotationMatrix(1,2) = 0.0_RKIND
       zRotationMatrix(2,1) = 0.0_RKIND
       zRotationMatrix(2,2) = 1.0_RKIND

       ! this should be (r, 0, 0)
       !cellCentreEquator = matmul(yRotationMatrix, matmul(zRotationMatrix, cellCentreRotated))
       
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          
          iEdge = edgesOnCell(iEdgeOnCell,iCell)

          ! vector in great circle plane
          iVertex1 = verticesOnEdge(1,iEdge)
          iVertex2 = verticesOnEdge(2,iEdge)

          ! perform grid rotations
          call grid_rotation_forward(edgeRotated(1), edgeRotated(2), edgeRotated(3), &
                                     xEdge(iEdge),   yEdge(iEdge),   zEdge(iEdge))
          
          call grid_rotation_forward(vertexRotated1(1), vertexRotated1(2), vertexRotated1(3), &
                                     xVertex(iVertex1), yVertex(iVertex1), zVertex(iVertex1))

          call grid_rotation_forward(vertexRotated2(1), vertexRotated2(2), vertexRotated2(3), &
                                     xVertex(iVertex2), yVertex(iVertex2), zVertex(iVertex2))

          ! rotate to equator in new coords
          edgeEquator    = matmul(yRotationMatrix, matmul(zRotationMatrix, edgeRotated))
          vertexEquator1 = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated1))
          vertexEquator2 = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated2))

          ! vector joining vertices
          vertexVector = vertexEquator2 - vertexEquator1

          ! form great circle plane normal vector
          call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                    vertexVector(1),      vertexVector(2),      vertexVector(3), &
                                    edgeEquator(1),       edgeEquator(2),       edgeEquator(3))

          if (iCell == cellsOnEdge(2,iEdge)) then
             normalGreatCircle = -1.0_RKIND * normalGreatCircle
          endif

          ! normalize normal vector
          normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
          normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
          normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
          normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm

          ! eastwards vector at edge
          eastwardsVector(1) = -edgeEquator(2)
          eastwardsVector(2) =  edgeEquator(1)
          eastwardsVector(3) = 0.0_RKIND
          eastwardsVector = eastwardsVector / &
                            sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)

          call dot_product_3space(normalVectorPolygon(1,iEdgeOnCell,iCell), &
                                  normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                  eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
          normalVectorPolygon(2,iEdgeOnCell,iCell) = sign(1.0_RKIND,normalGreatCircle(3))&
                                                   * sqrt(1.0_RKIND - max(min(normalVectorPolygon(1,iEdgeOnCell,iCell),1.0_RKIND),-1.0_RKIND)**2)

       enddo ! iEdgeOnCell
       
    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon_metric

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon(mesh, normalVectorPolygon)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon

    real(kind=RKIND), dimension(3) :: &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         xEdge0,   yEdge0,   zEdge0,   &
         xVertex1, yVertex1, zVertex1, &
         xVertex2, yVertex2, zVertex2

    real(kind=RKIND) :: &
         normalGreatCircleNorm

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell, &
         verticesOnEdge, &
         cellsOnEdge

    real(kind=RKIND), dimension(:), pointer :: &
         xEdge, &
         yEdge, &
         zEdge, &
         xVertex, &
         yVertex, &
         zVertex

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertex1, &
         iVertex2

    real(kind=RKIND), parameter :: &
         rad_to_deg = 180.0_RKIND / 3.14159265359_RKIND

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "zEdge", zEdge)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    do iCell = 1, nCells

       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          
          iEdge = edgesOnCell(iEdgeOnCell,iCell)

          ! vector in great circle plane
          iVertex1 = verticesOnEdge(1,iEdge)
          iVertex2 = verticesOnEdge(2,iEdge)

          ! perform grid rotations
          call grid_rotation_forward_nochange(xEdge0, yEdge0, zEdge0, &
                                     xEdge(iEdge), yEdge(iEdge), zEdge(iEdge))
          
          call grid_rotation_forward_nochange(xVertex1, yVertex1, zVertex1, &
                                     xVertex(iVertex1), yVertex(iVertex1), zVertex(iVertex1))

          call grid_rotation_forward_nochange(xVertex2, yVertex2, zVertex2, &
                                     xVertex(iVertex2), yVertex(iVertex2), zVertex(iVertex2))

          ! form great circle plane normal vector
          call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                    xVertex2 - xVertex1,  yVertex2 - yVertex1,  zVertex2 - zVertex1, &
                                    xEdge0,               yEdge0,               zEdge0)

          if (iCell == cellsOnEdge(2,iEdge)) then
             normalGreatCircle = -1.0_RKIND * normalGreatCircle
          endif

          ! normalize normal vector
          normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
          normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
          normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
          normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm

          ! eastwards vector at edge
          eastwardsVector(1) = -yEdge0
          eastwardsVector(2) =  xEdge0
          eastwardsVector(3) = 0.0_RKIND
          eastwardsVector = eastwardsVector / &
                            sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)

          call dot_product_3space(normalVectorPolygon(1,iEdgeOnCell,iCell), &
                                  normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                  eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
          normalVectorPolygon(2,iEdgeOnCell,iCell) = sign(1.0_RKIND,normalGreatCircle(3))&
                                                   * sqrt(1.0_RKIND - max(min(normalVectorPolygon(1,iEdgeOnCell,iCell),1.0_RKIND),-1.0_RKIND)**2)

       enddo ! iEdgeOnCell
       
    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_triangle_metric(mesh, normalVectorTriangle, latVertexRotated, interiorVertex)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle

    real(kind=RKIND), dimension(:), intent(out) :: &
         latVertexRotated

    integer, dimension(:), intent(in) :: &
         interiorVertex

    real(kind=RKIND), dimension(3) :: &
         vertexRotated, &
         vertexEquator, &
         edgeRotated, &
         cellRotated1, &
         cellRotated2, &
         edgeEquator, &
         cellEquator1, &
         cellEquator2, &
         cellVector, &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         lonVertexRotated, &
         normalGreatCircleNorm

    real(kind=RKIND), dimension(3,3) :: &
         yRotationMatrix, &
         zRotationMatrix

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         edgesOnVertex, &
         cellsOnEdge, &
         verticesOnEdge

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex, &
         xEdge, &
         yEdge, &
         zEdge, &
         xCell, &
         yCell, &
         zCell

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCell1, &
         iCell2

    real(kind=RKIND), parameter :: &
         pi = 3.14159265359_RKIND, &
         deg_to_rad = pi / 180.0_RKIND, &
         rad_to_deg = 180.0_RKIND / pi

    yRotationMatrix = 0.0_RKIND
    zRotationMatrix = 0.0_RKIND

    yRotationMatrix(2,2) = 1.0_RKIND
    zRotationMatrix(3,3) = 1.0_RKIND

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "zEdge", zEdge)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    do iVertex = 1, nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          ! rotate the cell centre to rotated geographical grid
          call grid_rotation_forward(vertexRotated(1), vertexRotated(2), vertexRotated(3), &
                                     xVertex(iVertex), yVertex(iVertex), zVertex(iVertex))
 
          ! calculate lon and lat of the cell centre in the rotated coordinate system
          lonVertexRotated          = atan2(vertexRotated(2), vertexRotated(1))
          latVertexRotated(iVertex) = asin(vertexRotated(3) / sphere_radius)
          
          ! rotate the cell centre onto the rotated equator
          yRotationMatrix(1,1) =  cos(latVertexRotated(iVertex))
          yRotationMatrix(1,3) =  sin(latVertexRotated(iVertex))
          yRotationMatrix(3,1) = -sin(latVertexRotated(iVertex))
          yRotationMatrix(3,3) =  cos(latVertexRotated(iVertex))
          
          zRotationMatrix(1,1) =  cos(-lonVertexRotated)
          zRotationMatrix(1,2) = -sin(-lonVertexRotated)
          zRotationMatrix(2,1) =  sin(-lonVertexRotated)
          zRotationMatrix(2,2) =  cos(-lonVertexRotated)

          ! this should be (r, 0, 0)
          !vertexEquator = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated))
          
          do iVertexDegree = 1, vertexDegree
             
             iEdge = edgesOnVertex(iVertexDegree,iVertex)

             ! vector in great circle plane
             iCell1 = cellsOnEdge(1,iEdge)
             iCell2 = cellsOnEdge(2,iEdge)

             ! perform grid rotations
             call grid_rotation_forward(edgeRotated(1), edgeRotated(2), edgeRotated(3), &
                                        xEdge(iEdge),   yEdge(iEdge),   zEdge(iEdge))

             call grid_rotation_forward(cellRotated1(1), cellRotated1(2), cellRotated1(3), &
                                        xCell(iCell1),   yCell(iCell1),   zCell(iCell1))

             call grid_rotation_forward(cellRotated2(1), cellRotated2(2), cellRotated2(3), &
                                        xCell(iCell2),   yCell(iCell2),   zCell(iCell2))

             ! rotate to equator in new coords
             edgeEquator  = matmul(yRotationMatrix, matmul(zRotationMatrix, edgeRotated))
             cellEquator1 = matmul(yRotationMatrix, matmul(zRotationMatrix, cellRotated1))
             cellEquator2 = matmul(yRotationMatrix, matmul(zRotationMatrix, cellRotated2))

             ! vector joining vertices
             cellVector = cellEquator2 - cellEquator1

             ! form great circle plane normal vector
             call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       cellVector(1),        cellVector(2),        cellVector(3), &
                                       edgeEquator(1),       edgeEquator(2),       edgeEquator(3))

             if (iVertex == verticesOnEdge(1,iEdge)) then
                normalGreatCircle = -1.0_RKIND * normalGreatCircle
             endif

             ! normalize normal vector
             normalGreatCircleNorm = sqrt(normalGreatCircle(1)**2 + normalGreatCircle(2)**2 + normalGreatCircle(3)**2)
             normalGreatCircle(1) = normalGreatCircle(1) / normalGreatCircleNorm
             normalGreatCircle(2) = normalGreatCircle(2) / normalGreatCircleNorm
             normalGreatCircle(3) = normalGreatCircle(3) / normalGreatCircleNorm
             
             ! eastwards vector at edge
             eastwardsVector(1) = -edgeEquator(2)
             eastwardsVector(2) =  edgeEquator(1)
             eastwardsVector(3) = 0.0_RKIND
             eastwardsVector = eastwardsVector / &
                               sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)
             
             call dot_product_3space(normalVectorTriangle(1,iVertexDegree,iVertex), &
                                     normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                     eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangle(2,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3)) &
                                                           * sqrt(1.0_RKIND - max(min(normalVectorTriangle(1,iVertexDegree,iVertex),1.0_RKIND),-1.0_RKIND)**2)

          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle_metric

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_triangle(mesh, normalVectorTriangle, interiorVertex)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle

    integer, dimension(:), intent(in) :: &
         interiorVertex

    real(kind=RKIND), dimension(3) :: &
         normalGreatCircle, &
         eastwardsVector

    real(kind=RKIND) :: &
         xEdge0, yEdge0, zEdge0,  &
         xCell1, yCell1, zCell1, &
         xCell2, yCell2, zCell2

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         edgesOnVertex, &
         cellsOnEdge, &
         verticesOnEdge

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell, &
         xEdge, &
         yEdge, &
         zEdge

    integer :: &
         iVertex, &
         iVertexDegree, &
         iEdge, &
         iCell1, &
         iCell2

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "zEdge", zEdge)

    do iVertex = 1, nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          do iVertexDegree = 1, vertexDegree
             
             iEdge = edgesOnVertex(iVertexDegree,iVertex)

             ! vector in great circle plane
             iCell1 = cellsOnEdge(1,iEdge)
             iCell2 = cellsOnEdge(2,iEdge)

             ! perform grid rotations
             call grid_rotation_forward(xEdge0,       yEdge0,       zEdge0, &
                                        xEdge(iEdge), yEdge(iEdge), zEdge(iEdge))

             call grid_rotation_forward(xCell1,        yCell1,        zCell1, &
                                        xCell(iCell1), yCell(iCell1), zCell(iCell1))

             call grid_rotation_forward(xCell2,        yCell2,        zCell2, &
                                        xCell(iCell2), yCell(iCell2), zCell(iCell2))

             ! form great circle plane normal vector
             call cross_product_3space(normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       xCell2 - xCell1,      yCell2 - yCell1,      zCell2 - zCell1,      &
                                       xEdge0,               yEdge0,               zEdge0)

             if (iVertex == verticesOnEdge(2,iEdge)) then
                normalGreatCircle = -1.0_RKIND * normalGreatCircle
             endif
             
             ! eastwards vector at edge
             eastwardsVector(1) = -yEdge0
             eastwardsVector(2) =  xEdge0
             eastwardsVector(3) = 0.0_RKIND
             eastwardsVector = eastwardsVector / &
                               sqrt(eastwardsVector(1)**2 + eastwardsVector(2)**2)
             
             call dot_product_3space(normalVectorTriangle(1,iVertexDegree,iVertex), &
                                     normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                     eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangle(2,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3))&
                                                           * sqrt(1.0_RKIND - min(normalVectorTriangle(1,iVertexDegree,iVertex),1.0_RKIND)**2)
             
          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle

  !-------------------------------------------------------------

  subroutine dot_product_3space(dot, x1, y1, z1, x2, y2, z2)

    real(kind=RKIND), intent(out) :: &
         dot

    real(kind=RKIND), intent(in) :: &
         x1, y1, z1, &
         x2, y2, z2

    dot = x1 * x2 + y1 * y2 + z1 * z2

  end subroutine dot_product_3space

  !-------------------------------------------------------------

  subroutine cross_product_3space(cp1, cp2, cp3, x1, y1, z1, x2, y2, z2)

    real(kind=RKIND), intent(out) :: &
         cp1, cp2, cp3

    real(kind=RKIND), intent(in) :: &
         x1, y1, z1, &
         x2, y2, z2

    cp1 = y1 * z2 - z1 * y2
    cp2 = z1 * x2 - x1 * z2
    cp3 = x1 * y2 - y1 * x2

  end subroutine cross_product_3space

  !-------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon_2(mesh, normalVectorPolygon, latCellRotated)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon ! TWO maxEdges nCells

    real(kind=RKIND), dimension(:), intent(out) :: &
         latCellRotated
    
    integer :: &
         iCell, &
         iEdgeOnCell, &
         iCell2, &
         iEdge

    real(kind=RKIND), dimension(3) :: &
         normalVector3D, &
         unitVectorEast, &
         unitVectorNorth

    real(kind=RKIND) :: &
         xCell0, yCell0, zCell0, &
         xCell2, yCell2, zCell2, &
         xEdge,  yEdge,  zEdge, &
         vectorMagnitude

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    real(kind=RKIND), pointer :: &
         sphere_radius
    
    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    ! loop over all cells
    do iCell = 1, nCells

       ! get rotated three space positions
       call grid_rotation_forward(&
            xCell0,       yCell0,       zCell0, &
            xCell(iCell), yCell(iCell), zCell(iCell))

       ! get the rotated latitude
       latCellRotated(iCell) = asin(zCell0 / sphere_radius)

       ! loop over edges of cell
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          ! get the polygon side normal unit vector
          call get_polygon_side_normal_3D_vector(mesh, iCell, iEdgeOnCell, xCell0, yCell0, zCell0, normalVector3D)

          ! project the 3D vector onto the local 2D plane
          call project_3D_vector_onto_local_2D(normalVectorPolygon(:,iEdgeOnCell,iCell), normalVector3D, xCell0, yCell0, zCell0)

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon_2

  !-------------------------------------------------------------
  
  subroutine normal_vectors_spherical_triangle_2(mesh, normalVectorTriangle, latVertexRotated, interiorVertex)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle ! TWO vertexDegree nVertices

    real(kind=RKIND), dimension(:), intent(out) :: &
         latVertexRotated

    integer, dimension(:), intent(in) :: &
         interiorVertex

    integer :: &
         iVertex, &
         iVertexDegree

    real(kind=RKIND), dimension(3) :: &
         normalVector3D

    real(kind=RKIND) :: &
         xVertex0, yVertex0, zVertex0

    integer, pointer :: &
         nVertices, &
         vertexDegree

    real(kind=RKIND), pointer :: &
         sphere_radius
    
    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "sphere_radius", sphere_radius)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    ! loop over interior vertices
    do iVertex = 1, nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          ! get rotated three space positions
          call grid_rotation_forward(&
               xVertex0,         yVertex0,         zVertex0, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex))

          ! get the rotated latitude
          latVertexRotated(iVertex) = asin(zVertex0 / sphere_radius)
          
          ! loop over edges of triangle
          do iVertexDegree = 1, vertexDegree

             ! get the triangle side normal unit vector
             call get_triangle_side_normal_3D_vector(mesh, iVertex, iVertexDegree, xVertex0, yVertex0, zVertex0, normalVector3D)

             ! project the 3D vector onto the local 2D plane
             call project_3D_vector_onto_local_2D(normalVectorTriangle(:,iVertexDegree,iVertex), normalVector3D, xVertex0, yVertex0, zVertex0)

          enddo ! iVertexDegree

       end if ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle_2

  !-------------------------------------------------------------

  subroutine project_3D_vector_onto_local_2D(normalVector2D, normalVector3D, xPoint, yPoint, zPoint)

    real(kind=RKIND), dimension(2), intent(out) :: &
         normalVector2D

    real(kind=RKIND), dimension(3), intent(in) :: &
         normalVector3D

    real(kind=RKIND), intent(in) :: &
         xPoint, yPoint, zPoint

    real(kind=RKIND), dimension(3) :: &
         unitVectorEast, &
         unitVectorNorth

    real(kind=RKIND) :: &
         vectorMagnitude

    ! calculate the local eastern and northern unit vectors at the point
    call local_eastern_and_northern_unit_vectors(&
         xPoint, yPoint, zPoint, &
         unitVectorEast, &
         unitVectorNorth)

    ! project normal 3D vector onto U 2D direction - eastwards
    call dot_product_3space(&
         normalVector2D(1), &
         normalVector3D(1), normalVector3D(2), normalVector3D(3), &
         unitVectorEast(1), unitVectorEast(2), unitVectorEast(3))

    ! project normal 3D vector onto V 2D direction - northwards
    call dot_product_3space(&
         normalVector2D(2), &
         normalVector3D(1),  normalVector3D(2),  normalVector3D(3), &
         unitVectorNorth(1), unitVectorNorth(2), unitVectorNorth(3))

    ! magnitude of projected vector
    vectorMagnitude = sqrt(normalVector2D(1)**2 + normalVector2D(2)**2)

    ! normalize projected vector
    normalVector2D(1) = normalVector2D(1) / vectorMagnitude
    normalVector2D(2) = normalVector2D(2) / vectorMagnitude

  end subroutine project_3D_vector_onto_local_2D

  !-------------------------------------------------------------

  subroutine get_polygon_side_normal_3D_vector(mesh, iCell, iEdgeOnCell, xCell0, yCell0, zCell0, normalVector3D)

    type(MPAS_pool_type), intent(in) :: mesh

    integer, intent(in) :: &
         iCell, &
         iEdgeOnCell

    real(kind=RKIND), intent(in) :: &
         xCell0, yCell0, zCell0

    real(kind=RKIND), dimension(3), intent(out) :: &
         normalVector3D

    integer :: &
         iCell2, &
         iEdge

    real(kind=RKIND) :: &
         xCell2, yCell2, zCell2, &
         xEdge0, yEdge0, zEdge0

    integer, pointer :: &
         nCells

    integer, dimension(:,:), pointer :: &
         cellsOnCell, &
         edgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell, &
         xEdge, &
         yEdge, &
         zEdge

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "cellsOnCell", cellsOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "zEdge", zEdge)
    
    ! get cell other side of edge
    iCell2 = cellsOnCell(iEdgeOnCell, iCell)

    ! check if opposite cell exists
    if (iCell2 <= nCells) then

       ! interior

       ! get rotated three space positions
       call grid_rotation_forward(&
            xCell2,        yCell2,        zCell2, &
            xCell(iCell2), yCell(iCell2), zCell(iCell2))

       normalVector3D(1) = xCell2 - xCell0
       normalVector3D(2) = yCell2 - yCell0
       normalVector3D(3) = zCell2 - zCell0

    else

       ! at domain edge

       iEdge = edgesOnCell(iEdgeOnCell, iCell)

       ! get rotated three space positions
       call grid_rotation_forward(&
            xEdge0,       yEdge0,       zEdge0, &
            xEdge(iEdge), yEdge(iEdge), zEdge(iEdge))

       normalVector3D(1) = xEdge0 - xCell0
       normalVector3D(2) = yEdge0 - yCell0
       normalVector3D(3) = zEdge0 - zCell0

    endif

  end subroutine get_polygon_side_normal_3D_vector

  !-------------------------------------------------------------

  subroutine get_triangle_side_normal_3D_vector(mesh, iVertex, iVertexDegree, xVertex0, yVertex0, zVertex0, normalVector3D)

    type(MPAS_pool_type), intent(in) :: mesh

    integer, intent(in) :: &
         iVertex, &
         iVertexDegree

    real(kind=RKIND), intent(in) :: &
         xVertex0, yVertex0, zVertex0

    real(kind=RKIND), dimension(3), intent(out) :: &
         normalVector3D

    integer :: &
         iEdge, &
         iVertexOnEdge, &
         iVertex2
         
    real(kind=RKIND) :: &
         xVertex2, yVertex2, zVertex2

    integer, dimension(:,:), pointer :: &
         edgesOnVertex, &
         verticesOnEdge

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    ! get vertex other side of edge
    iEdge = edgesOnVertex(iVertexDegree, iVertex)

    do iVertexOnEdge = 1, 2
       iVertex2 = verticesOnEdge(iVertexOnEdge, iEdge)
       if (iVertex2 /= iVertex) exit
    enddo ! iVertexOnEdge
    
    ! get rotated three space positions
    call grid_rotation_forward(&
         xVertex2,          yVertex2,          zVertex2, &
         xVertex(iVertex2), yVertex(iVertex2), zVertex(iVertex2))
    
    normalVector3D(1) = xVertex2 - xVertex0
    normalVector3D(2) = yVertex2 - yVertex0
    normalVector3D(3) = zVertex2 - zVertex0

  end subroutine get_triangle_side_normal_3D_vector
  
  !-------------------------------------------------------------

  subroutine local_eastern_and_northern_unit_vectors(&
       xPoint, yPoint, zPoint, &
       unitVectorEast, &
       unitVectorNorth)

    real(kind=RKIND), intent(in) :: &
         xPoint, yPoint, zPoint

    real(kind=RKIND), dimension(3), intent(out) :: &
         unitVectorEast, &
         unitVectorNorth

    real(kind=RKIND) :: &
         vectorMagnitude

    ! determine unit vector at cell centre pointing in U directions - local eastwards
    unitVectorEast(1) = -yPoint
    unitVectorEast(2) =  xPoint
    unitVectorEast(3) =  0.0_RKIND

    vectorMagnitude = sqrt(unitVectorEast(1)**2 + unitVectorEast(2)**2 + unitVectorEast(3)**2)

    unitVectorEast(1) = unitVectorEast(1) / vectorMagnitude
    unitVectorEast(2) = unitVectorEast(2) / vectorMagnitude
    unitVectorEast(3) = unitVectorEast(3) / vectorMagnitude

    ! determine unit vector at cell centre pointing in V directions - local northwards
    if (zPoint /= 0.0_RKIND) then

       unitVectorNorth(1) = -xPoint
       unitVectorNorth(2) = -yPoint
       unitVectorNorth(3) = (xPoint**2 + yPoint**2) / zPoint

       vectorMagnitude = sqrt(unitVectorNorth(1)**2 + unitVectorNorth(2)**2 + unitVectorNorth(3)**2)

       unitVectorNorth(1) = unitVectorNorth(1) / vectorMagnitude
       unitVectorNorth(2) = unitVectorNorth(2) / vectorMagnitude
       unitVectorNorth(3) = unitVectorNorth(3) / vectorMagnitude

    else   

       unitVectorNorth(1) = 0.0_RKIND
       unitVectorNorth(2) = 0.0_RKIND
       unitVectorNorth(3) = 1.0_RKIND

    endif

  end subroutine local_eastern_and_northern_unit_vectors

  !-------------------------------------------------------------
  ! rotated lat-lon grid
  !-------------------------------------------------------------

  subroutine grid_rotation_forward(xp, yp, zp, x, y, z)

    ! rotate xyz coordinates from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    real(kind=RKIND), intent(in) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    !call grid_rotation_forward_nochange(xp, yp, zp, x, y, z)

    xp = -z
    yp = y
    zp = x

  end subroutine grid_rotation_forward

  !-------------------------------------------------------------

  subroutine grid_rotation_forward_nochange(xp, yp, zp, x, y, z)

    ! leave the grid unchanged - for testing

    real(kind=RKIND), intent(out) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    real(kind=RKIND), intent(in) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    xp = x
    yp = y
    zp = z

  end subroutine grid_rotation_forward_nochange

  !-------------------------------------------------------------

  subroutine grid_rotation_backward(x, y, z, xp, yp, zp)

    ! rotate xyz coordinates from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    real(kind=RKIND), intent(in) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    !call grid_rotation_backward_nochange(x, y, z, xp, yp, zp)

    x = zp
    y = yp
    z = -xp

  end subroutine grid_rotation_backward

  !-------------------------------------------------------------

  subroutine grid_rotation_backward_nochange(x, y, z, xp, yp, zp)

    ! rotate xyz coordinates from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         x,  & ! x position of point on geographical grid
         y,  & ! y position of point on geographical grid
         z     ! z position of point on geographical grid

    real(kind=RKIND), intent(in) :: &
         xp, & ! x position of point on rotated grid
         yp, & ! y position of point on rotated grid
         zp    ! z position of point on rotated grid

    x = xp
    y = yp
    z = zp

  end subroutine grid_rotation_backward_nochange

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_forward(vxp, vyp, vzp, vx, vy, vz)

    ! rotate a xyz vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    real(kind=RKIND), intent(in) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid
    
    !call xyz_vector_rotation_forward_nochange(vxp, vyp, vzp, vx, vy, vz) 

    vxp = -vz
    vyp =  vy
    vzp =  vx

  end subroutine xyz_vector_rotation_forward

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_forward_nochange(vxp, vyp, vzp, vx, vy, vz)

    ! rotate a xyz vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    real(kind=RKIND), intent(in) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid
    
    vxp = vx
    vyp = vy
    vzp = vz

  end subroutine xyz_vector_rotation_forward_nochange

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_backward(vx, vy, vz, vxp, vyp, vzp)

    ! rotate xyz vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid

    real(kind=RKIND), intent(in) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    !call xyz_vector_rotation_backward_nochange(vx, vy, vz, vxp, vyp, vzp)

    vx =  vzp
    vy =  vyp
    vz = -vxp

  end subroutine xyz_vector_rotation_backward

  !-------------------------------------------------------------

  subroutine xyz_vector_rotation_backward_nochange(vx, vy, vz, vxp, vyp, vzp)

    ! rotate xyz vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         vx,  & ! x component of velocity vector on geographical grid
         vy,  & ! y component of velocity vector on geographical grid
         vz     ! z component of velocity vector on geographical grid

    real(kind=RKIND), intent(in) :: &
         vxp, & ! x component of velocity vector on rotated grid
         vyp, & ! y component of velocity vector on rotated grid
         vzp    ! z component of velocity vector on rotated grid

    vx = vxp
    vy = vyp
    vz = vzp

  end subroutine xyz_vector_rotation_backward_nochange

  !-------------------------------------------------------------

  subroutine latlon_vector_rotation_forward(up, vp, u, v, lat, lon, x, y, z, r)

    ! rotate a latlon vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         up,   & ! u component of velocity in rotated grid
         vp      ! v component of velocity in rotated grid

    real(kind=RKIND), intent(in) :: &
         u,    & ! u component of velocity on geographical grid
         v,    & ! v component of velocity on geographical grid
         lat,  & ! latitude of point on geographical grid
         lon,  & ! longitude of point on geographical grid
         x,    & ! x position of point on geographical grid
         y,    & ! y position of point on geographical grid
         z,    & ! z position of point on geographical grid
         r       ! radius of the earth 

    real(kind=RKIND) :: &
         xp,   & ! x position of point on rotated grid
         yp,   & ! y position of point on rotated grid
         zp,   & ! z position of point on rotated grid
         latp, & ! latitude of point on rotated grid
         lonp, & ! longitude of point on rotated grid
         vx,   & ! x component of velocity vector on geographical grid
         vy,   & ! y component of velocity vector on geographical grid
         vz,   & ! z component of velocity vector on geographical grid
         vxp,  & ! x component of velocity vector on rotated grid
         vyp,  & ! y component of velocity vector on rotated grid
         vzp     ! z component of velocity vector on rotated grid

    ! perform rotation of the point from geographical grid to rotated grid
    call grid_rotation_forward(xp, yp, zp, x, y, z)

    ! calculate latitude and longitude of the point in the rotated grid
    call latlon_from_xyz(latp, lonp, xp, yp, zp, r)

    ! convert lat lon vector to xyz vector on geographical grid
    call latlon_vector_to_xyz_vector(vx, vy, vz, u, v, lat, lon)

    ! perform rotation of geographical xyz vector to rotated grid
    call xyz_vector_rotation_forward(vxp, vyp, vzp, vx, vy, vz)

    ! convert xyz vector to lat lon vector on rotated grid
    call xyz_vector_to_latlon_vector(up, vp, vxp, vyp, vzp, latp, lonp)

  end subroutine latlon_vector_rotation_forward

  !-------------------------------------------------------------

  subroutine latlon_vector_rotation_backward(u, v, up, vp, lat, lon, x, y, z, r)

    ! rotate latlon vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         u,    & ! u component of velocity on geographical grid
         v       ! v component of velocity on geographical grid

    real(kind=RKIND), intent(in) :: &
         up,   & ! u component of velocity in rotated grid
         vp,   & ! v component of velocity in rotated grid
         lat,  & ! latitude of point on geographical grid
         lon,  & ! longitude of point on geographical grid
         x,    & ! x position of point on geographical grid
         y,    & ! y position of point on geographical grid
         z,    & ! z position of point on geographical grid
         r       ! radius of the earth 

    real(kind=RKIND) :: &
         xp,   & ! x position of point on rotated grid
         yp,   & ! y position of point on rotated grid
         zp,   & ! z position of point on rotated grid
         latp, & ! latitude of point on rotated grid
         lonp, & ! longitude of point on rotated grid
         vx,   & ! x component of velocity vector on geographical grid
         vy,   & ! y component of velocity vector on geographical grid
         vz,   & ! z component of velocity vector on geographical grid
         vxp,  & ! x component of velocity vector on rotated grid
         vyp,  & ! y component of velocity vector on rotated grid
         vzp     ! z component of velocity vector on rotated grid

    ! perform rotation of the point from geographical grid to rotated grid
    call grid_rotation_forward(xp, yp, zp, x, y, z)

    ! calculate latitude and longitude of the point in the rotated grid
    call latlon_from_xyz(latp, lonp, xp, yp, zp, r)

    ! convert lat lon vector to xyz vector on rotated grid
    call latlon_vector_to_xyz_vector(vxp, vyp, vzp, up, vp, latp, lonp)

    ! perform rotation of rotated xyz vector to geographical grid
    call xyz_vector_rotation_backward(vx, vy, vz, vxp, vyp, vzp)

    ! convert xyz vector to lat lon vector on geographical grid
    call xyz_vector_to_latlon_vector(u, v, vx, vy, vz, lat, lon)
     
  end subroutine latlon_vector_rotation_backward

  !-------------------------------------------------------------

  subroutine latlon_vector_to_xyz_vector(vx, vy, vz, u, v, lat, lon)

    ! convert a latlon vector to a xyz vector

    real(kind=RKIND), intent(out) :: &
         vx,  & ! x component of velocity vector
         vy,  & ! y component of velocity vector
         vz     ! z component of velocity vector

    real(kind=RKIND), intent(in) :: &
         u,   & ! u component of velocity
         v,   & ! v component of velocity
         lat, & ! latitude of point 
         lon    ! longitude of point 

    vx = (-u) * sin(lon) - v * sin(lat) * cos(lon)
    vy =   u  * cos(lon) - v * sin(lat) * sin(lon)
    vz =                   v * cos(lat)

  end subroutine latlon_vector_to_xyz_vector

  !-------------------------------------------------------------

  subroutine xyz_vector_to_latlon_vector(u, v, vx, vy, vz, lat, lon)

    ! convert a xyz vector vector to a latlon vector

    real(kind=RKIND), intent(out) :: &
         u,   & ! u component of velocity
         v      ! v component of velocity

    real(kind=RKIND), intent(in) :: &
         vx,  & ! x component of velocity vector
         vy,  & ! y component of velocity vector
         vz,  & ! z component of velocity vector
         lat, & ! latitude of point 
         lon    ! longitude of point 

    u = (-sin(lon)) * vx + &
          cos(lon)  * vy
    
    !write(*,*) 1, (-sin(lat)) * cos(lon), vx, (-sin(lat)) * cos(lon) * vx, (-sin(lat)) * cos(lon),  vy, (-sin(lat)) * cos(lon) * vy, cos(lat),  vz, cos(lat)  * vz

    v = (-sin(lat)) * cos(lon) * vx + &
        (-sin(lat)) * sin(lon) * vy + &
          cos(lat)  * vz

  end subroutine xyz_vector_to_latlon_vector

  !-------------------------------------------------------------

  subroutine latlon_from_xyz(lat, lon, x, y, z, r)

    ! given xyz coordinates determine the latitude and longitude

    real(kind=RKIND), intent(out) :: &
         lat, & ! latitude of point 
         lon    ! longitude of point 

    real(kind=RKIND), intent(in) :: &
         x,   & ! x position of point
         y,   & ! y position of point
         z,   & ! z position of point
         r      ! radius of earth

    lon = atan2(y, x)
    lat = asin(z/r)

  end subroutine latlon_from_xyz

  !-------------------------------------------------------------
  ! interpolaton
  !-------------------------------------------------------------

  subroutine interpolate_cell_to_vertex(mesh, variableVertex, variableCell)

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         variableVertex

    real(kind=RKIND), dimension(:), intent(in) :: &
         variableCell

    real(kind=RKIND) :: &
         totalArea

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell, &
         areaTriangle

    real(kind=RKIND), dimension(:,:), pointer :: &
         kiteAreasOnVertex

    ! init variables
    call MPAS_pool_get_array(mesh, "nVertices", nVertices)
    call MPAS_pool_get_array(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)
    call MPAS_pool_get_array(mesh, "kiteAreasOnVertex", kiteAreasOnVertex)

    do iVertex = 1, nVertices
    
       variableVertex(iVertex) = 0.0_RKIND
       totalArea = 0.0_RKIND

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex) 

          !variableVertex(iVertex) = variableVertex(iVertex) + kiteAreasOnVertex(iVertexDegree, iVertex) * variableCell(iCell)

          variableVertex(iVertex) = variableVertex(iVertex) + areaCell(iCell) * variableCell(iCell)
          totalArea = totalArea + areaCell(iCell)

       enddo ! iVertexDegree

       variableVertex(iVertex) = variableVertex(iVertex) / totalArea

       !variableVertex(iVertex) = variableVertex(iVertex) / areaTriangle(iVertex)

    enddo ! iVertex    

  end subroutine interpolate_cell_to_vertex

  !-------------------------------------------------------------

end module cice_mesh
