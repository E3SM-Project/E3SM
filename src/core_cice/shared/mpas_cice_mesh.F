!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_mesh
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_mesh

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar

  implicit none

  private
  save

  public :: &
       cice_init_mesh, &
       cice_cell_vertices_at_vertex, &
       cice_normal_vectors, &
       cice_normal_vectors_polygon, &
       cice_dot_product_3space, &
       cice_grid_rotation_forward, &
       cice_latlon_vector_rotation_forward, &
       cice_latlon_vector_rotation_backward, &
       cice_latlon_from_xyz, &
       cice_interpolate_cell_to_vertex, &
       cice_project_3D_vector_onto_local_2D

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_mesh
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 7th August 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_mesh(&
       domain)!{{{

    type (domain_type), intent(inout) :: &
         domain !< Input/Output: 

    ! rescale mesh quantities
    call rescale_mesh(domain)

    ! calculate coriolois parameter
    call coriolis_parameter(domain)

    ! boundary variables
    call init_boundary(domain)

  end subroutine cice_init_mesh!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  rescale_mesh
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 7th August 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine rescale_mesh(&
       domain)!{{{

    use mpas_constants, only: &
         a ! This is the excellently named Earth Radius - couldn't you tell?

    type (domain_type), intent(inout) :: &
         domain !< Input/Output: 

    logical, pointer :: &
         on_a_sphere

    real(kind=RKIND), pointer :: &
         oldRadius

    real(kind=RKIND) :: &
         newRadius

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool

    integer, pointer :: &
         nCells, &
         nEdges, &
         nVertices, &
         vertexDegree

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell, &
         areaCell, &
         xEdge, &
         yEdge, &
         zEdge, &
         dvEdge, &
         dcEdge, &
         xVertex, &
         yVertex, &
         zVertex, &
         areaTriangle

    real(kind=RKIND), dimension(:,:), pointer :: &
         kiteAreasOnVertex

    integer, dimension(:,:), pointer :: &
         cellsOnVertex
  
    real(kind=RKIND) :: &
         oldX, &
         oldY, &
         oldZ, &
         norm

    integer :: &
         iCell, &
         iVertex, &
         iEdge, &
         iVertexDegree
    
    block => domain % blocklist
    do while(associated(block))

       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

       call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

       if (on_a_sphere) then

          ! get the old radius
          call mpas_pool_get_config(meshPool, 'sphere_radius', oldRadius)

          ! get the new radius
          newRadius = a

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
          call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
          call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

          call mpas_pool_get_array(meshPool, 'xCell', xCell)
          call mpas_pool_get_array(meshPool, 'yCell', yCell)
          call mpas_pool_get_array(meshPool, 'zCell', zCell)
          call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

          call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
          call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
          call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
          call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
          call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

          call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
          call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
          call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
          call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
          call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
          call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

          ! rescale cell quantities
          do iCell = 1, nCells

             oldX = xCell(iCell)
             oldY = yCell(iCell)
             oldZ = zCell(iCell)

             norm = sqrt(oldX**2 + oldY**2 + oldZ**2)

             xCell(iCell) = (oldX / norm) * newRadius
             yCell(iCell) = (oldY / norm) * newRadius
             zCell(iCell) = (oldZ / norm) * newRadius
             areaCell(iCell) = (areaCell(iCell) / (oldRadius**2) )* newRadius**2

          end do

          ! rescale vertex quantities
          do iVertex = 1, nVertices

             oldX = xVertex(iVertex)
             oldY = yVertex(iVertex)
             oldZ = zVertex(iVertex)

             norm = sqrt(oldX**2 + oldY**2 + oldZ**2)

             xVertex(iVertex) = (oldX / norm) * newRadius
             yVertex(iVertex) = (oldY / norm) * newRadius
             zVertex(iVertex) = (oldZ / norm) * newRadius
             areaTriangle(iVertex) = 0.0_RKIND

             do iVertexDegree = 1, vertexDegree
                if (cellsOnVertex(iVertexDegree, iVertex) < nCells+1) then
                   kiteAreasOnVertex(iVertexDegree, iVertex) = ( kiteAreasOnVertex(iVertexDegree, iVertex) / oldRadius**2) * newRadius**2
                else
                   kiteAreasOnVertex(iVertexDegree, iVertex) = 0.0_RKIND
                end if
                areaTriangle(iVertex) = areaTriangle(iVertex) + kiteAreasOnVertex(iVertexDegree, iVertex)
             end do

          end do

          ! Expand edge quantities
          do iEdge = 1, nEdges

             oldX = xEdge(iEdge)
             oldY = yEdge(iEdge)
             oldZ = zEdge(iEdge)

             norm = sqrt(oldX**2 + oldY**2 + oldZ**2)

             xEdge(iEdge) = (oldX / norm) * newRadius
             yEdge(iEdge) = (oldY / norm) * newRadius
             zEdge(iEdge) = (oldZ / norm) * newRadius
             dvEdge(iEdge) = (dvEdge(iEdge) / oldRadius) * newRadius
             dcEdge(iEdge) = (dcEdge(iEdge) / oldRadius) * newRadius

          end do

       endif

       block => block % next
    end do

  end subroutine rescale_mesh!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  coriolis_parameter
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 7th August 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine coriolis_parameter(&
       domain)!{{{

    use mpas_constants, only: &
         omega

    type (domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         latVertex, &
         fVertex

    integer, pointer :: &
         nVertices

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_dimension(meshPool, "nVertices", nVertices)

       call MPAS_pool_get_array(meshPool, "latVertex", latVertex)
       call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

       do iVertex = 1, nVertices
       
          fVertex(iVertex) = 2.0_RKIND * omega * sin(latVertex(iVertex))

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine coriolis_parameter!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_boundary
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------
  
  subroutine init_boundary(&
       domain)!{{{

    type (domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary

    integer, dimension(:), pointer :: &
         blockIDout

    type (MPAS_pool_type), pointer :: &
         boundaryField

    type(field1DInteger), pointer :: &
         interiorVertexField, &
         interiorCellField, &
         interiorEdgeField

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
       
       call MPAS_pool_get_array(boundary, "blockIDout", blockIDout)
       blockIDout = block % blockID
       
       call interior_vertices(mesh, boundary)
       call interior_cells(mesh, boundary)
       call interior_edges(mesh, boundary)

       block => block % next
    end do

    ! exchange fields
    call MPAS_pool_get_subpool(domain % blocklist % structs, "boundary", boundaryField)

    call MPAS_pool_get_field(boundaryField, "interiorVertex", interiorVertexField)
    call MPAS_pool_get_field(boundaryField, "interiorCell", interiorCellField)
    call MPAS_pool_get_field(boundaryField, "interiorEdge", interiorEdgeField)

    call MPAS_dmpar_exch_halo_field(interiorVertexField)
    call MPAS_dmpar_exch_halo_field(interiorCellField)
    call MPAS_dmpar_exch_halo_field(interiorEdgeField)

  end subroutine init_boundary!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  interior_vertices
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine interior_vertices(&
       mesh, &
       boundary)!{{{

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh

    type (MPAS_pool_type), pointer :: &
         boundary !< Input: 

    integer :: &
         nInteriorAdjacentCells

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    integer, pointer :: &
         nVerticesSolve, &
         vertexDegree, &
         nCells

    integer, dimension(:), pointer :: &
         interiorVertex

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)

    ! boundary vertices
    do iVertex = 1, nVerticesSolve

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

  end subroutine interior_vertices!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  interior_cells
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine interior_cells(&
       mesh, &
       boundary)!{{{

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type (MPAS_pool_type), pointer, intent(inout) :: &
         boundary !< Input: 

    integer :: &
         iCell, &
         iEdgeOnCell

    integer, pointer :: &
         nCellsSolve, &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnCell

    integer, dimension(:), pointer :: &
         interiorCell

    call MPAS_pool_get_array(boundary, "interiorCell", interiorCell)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnCell", cellsOnCell)

    do iCell = 1, nCellsSolve

       interiorCell(iCell) = 1

       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          if (cellsOnCell(iEdgeOnCell,iCell) > nCells) then

             interiorCell(iCell) = 0

          endif

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine interior_cells!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  interior_cells
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine interior_edges(&
       mesh, &
       boundary)!{{{

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type (MPAS_pool_type), pointer, intent(inout) :: &
         boundary !< Input: 

    integer, pointer :: &
         nCells, &
         nEdgesSolve

    integer, dimension(:,:), pointer :: &
         cellsOnEdge

    integer, dimension(:), pointer :: &
         interiorEdge

    integer :: &
         iEdge, &
         iCellOnEdge
    
    call MPAS_pool_get_array(boundary, "interiorEdge", interiorEdge)

    call MPAS_pool_get_dimension(mesh, "nEdgesSolve", nEdgesSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)

    do iEdge = 1, nEdgesSolve

       interiorEdge(iEdge) = 1

       do iCellOnEdge = 1, 2

          if (cellsOnEdge(iCellOnEdge,iEdge) > nCells) then

             interiorEdge(iEdge) = 0

          endif

       enddo ! iCellOnEdge

    end do ! iEdge

  end subroutine interior_edges

!-----------------------------------------------------------------------
! mesh searches
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_cell_vertices_at_vertex
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_cell_vertices_at_vertex(&
       mesh, &
       cellVerticesAtVertex)!{{{
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    integer, dimension(:,:), intent(out) :: &
         cellVerticesAtVertex !< Output: 

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
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

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

  end subroutine cice_cell_vertices_at_vertex!}}}

!-----------------------------------------------------------------------
! normal vectors
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_normal_vectors
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_normal_vectors(&
       mesh, &
       normalVectorPolygon, &
       normalVectorTriangle, &
       interiorVertex, &
       rotateCartesianGrid, &
       removeMetricTerms, &
       latCellRotated, &
       latVertexRotated)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon, & !< Output: 
         normalVectorTriangle   !< Output: 

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         removeMetricTerms !< Input: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         latCellRotated, & !< Output: 
         latVertexRotated  !< Output: 

    call cice_normal_vectors_polygon(&
         mesh, &
         normalVectorPolygon, &
         rotateCartesianGrid, &
         removeMetricTerms, &
         latCellRotated)
    
    call cice_normal_vectors_triangle(&
         mesh, &
         normalVectorTriangle, &
         interiorVertex, &
         rotateCartesianGrid, &
         removeMetricTerms, &
         latVertexRotated)

  end subroutine cice_normal_vectors!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_normal_vectors_polygon
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_normal_vectors_polygon(&
       mesh, &
       normalVectorPolygon, &
       rotateCartesianGrid, &
       removeMetricTerms, &
       latCellRotated)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon !< Output: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         removeMetricTerms !< Input: 

    real(kind=RKIND), dimension(:), optional, intent(out) :: &
         latCellRotated !< Output: 

    logical, pointer :: &
         on_a_sphere

    call MPAS_pool_get_config(mesh, "on_a_sphere", on_a_sphere)

    if (on_a_sphere) then
       call normal_vectors_spherical_polygon_metric(mesh, rotateCartesianGrid, removeMetricTerms, normalVectorPolygon, latCellRotated)
    else 
       call normal_vectors_planar_polygon(mesh, normalVectorPolygon)
       if (present(latCellRotated)) latCellRotated = 0.0_RKIND
    endif

  end subroutine cice_normal_vectors_polygon!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_normal_vectors_triangle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_normal_vectors_triangle(&
       mesh, &
       normalVectorTriangle, &
       interiorVertex, &
       rotateCartesianGrid, &
       removeMetricTerms, &
       latVertexRotated)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle !< Output: 

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         removeMetricTerms !< Input: 

    real(kind=RKIND), dimension(:), optional, intent(out) :: &
         latVertexRotated !< Output: 

    logical, pointer :: &
         on_a_sphere

    call MPAS_pool_get_config(mesh, "on_a_sphere", on_a_sphere)

    if (on_a_sphere) then
       call normal_vectors_spherical_triangle_metric(mesh, rotateCartesianGrid, removeMetricTerms, normalVectorTriangle, interiorVertex, latVertexRotated)
    else 
       call normal_vectors_planar_triangle(mesh, normalVectorTriangle, interiorVertex)
       if (present(latVertexRotated)) latVertexRotated = 0.0_RKIND
    endif

  end subroutine cice_normal_vectors_triangle!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_planar_polygon
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine normal_vectors_planar_polygon(&
       mesh, &
       normalVectorPolygon)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon !< Output: 

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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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

  end subroutine normal_vectors_planar_polygon!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_planar_triangle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine normal_vectors_planar_triangle(&
       mesh, &
       normalVectorTriangle, &
       interiorVertex)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle !< Output: 

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input: 

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
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

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

  end subroutine normal_vectors_planar_triangle!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_spherical_polygon_metric
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon_metric(&
       mesh, &
       rotateCartesianGrid, &
       removeMetricTerms, &
       normalVectorPolygon, &
       latCellRotatedOut)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         removeMetricTerms !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon !< Output: 

    real(kind=RKIND), dimension(:), optional, intent(out) :: &
         latCellRotatedOut !< Output: 
    
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
         latCellRotated, &
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

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
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
       call cice_grid_rotation_forward(&
            cellCentreRotated(1), cellCentreRotated(2), cellCentreRotated(3), &
            xCell(iCell),         yCell(iCell),         zCell(iCell), &
            rotateCartesianGrid)

       ! calculate lon and lat of the cell centre in the rotated coordinate system
       lonCellRotated = atan2(cellCentreRotated(2), cellCentreRotated(1))
       latCellRotated = asin(cellCentreRotated(3) / sphere_radius)

       ! rotate the cell centre onto the rotated equator
       if (removeMetricTerms) then
          yRotationMatrix(1,1) =  cos(latCellRotated)
          yRotationMatrix(1,3) =  sin(latCellRotated)
          yRotationMatrix(3,1) = -sin(latCellRotated)
          yRotationMatrix(3,3) =  cos(latCellRotated)
          
          zRotationMatrix(1,1) =  cos(-lonCellRotated)
          zRotationMatrix(1,2) = -sin(-lonCellRotated)
          zRotationMatrix(2,1) =  sin(-lonCellRotated)
          zRotationMatrix(2,2) =  cos(-lonCellRotated)
       else
          yRotationMatrix(1,1) = 1.0_RKIND
          yRotationMatrix(1,3) = 0.0_RKIND
          yRotationMatrix(3,1) = 0.0_RKIND
          yRotationMatrix(3,3) = 1.0_RKIND
          
          zRotationMatrix(1,1) = 1.0_RKIND
          zRotationMatrix(1,2) = 0.0_RKIND
          zRotationMatrix(2,1) = 0.0_RKIND
          zRotationMatrix(2,2) = 1.0_RKIND
       endif

       ! this should be (r, 0, 0)
       !cellCentreEquator = matmul(yRotationMatrix, matmul(zRotationMatrix, cellCentreRotated))
       
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          
          iEdge = edgesOnCell(iEdgeOnCell,iCell)

          ! vector in great circle plane
          iVertex1 = verticesOnEdge(1,iEdge)
          iVertex2 = verticesOnEdge(2,iEdge)

          ! perform grid rotations
          call cice_grid_rotation_forward(&
               edgeRotated(1), edgeRotated(2), edgeRotated(3), &
               xEdge(iEdge),   yEdge(iEdge),   zEdge(iEdge), &
               rotateCartesianGrid)
          
          call cice_grid_rotation_forward(&
               vertexRotated1(1), vertexRotated1(2), vertexRotated1(3), &
               xVertex(iVertex1), yVertex(iVertex1), zVertex(iVertex1), &
               rotateCartesianGrid)

          call cice_grid_rotation_forward(&
               vertexRotated2(1), vertexRotated2(2), vertexRotated2(3), &
               xVertex(iVertex2), yVertex(iVertex2), zVertex(iVertex2), &
               rotateCartesianGrid)

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

          call cice_dot_product_3space(normalVectorPolygon(1,iEdgeOnCell,iCell), &
                                       normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
          normalVectorPolygon(2,iEdgeOnCell,iCell) = sign(1.0_RKIND,normalGreatCircle(3))&
               * sqrt(1.0_RKIND - max(min(normalVectorPolygon(1,iEdgeOnCell,iCell),1.0_RKIND),-1.0_RKIND)**2)

       enddo ! iEdgeOnCell

       if (present(latCellRotatedOut)) latCellRotatedOut(iCell) = latCellRotated
       
    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon_metric!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_spherical_polygon
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon(&
       mesh, &
       rotateCartesianGrid, &
       normalVectorPolygon)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon !< Output: 

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

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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
          call cice_grid_rotation_forward(&
               xEdge0,       yEdge0,       zEdge0, &
               xEdge(iEdge), yEdge(iEdge), zEdge(iEdge), &
               rotateCartesianGrid)
          
          call cice_grid_rotation_forward(&
               xVertex1,          yVertex1,          zVertex1, &
               xVertex(iVertex1), yVertex(iVertex1), zVertex(iVertex1), &
               rotateCartesianGrid)

          call cice_grid_rotation_forward(&
               xVertex2,          yVertex2,          zVertex2, &
               xVertex(iVertex2), yVertex(iVertex2), zVertex(iVertex2), &
               rotateCartesianGrid)

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

          call cice_dot_product_3space(normalVectorPolygon(1,iEdgeOnCell,iCell), &
                                       normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                       eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
          normalVectorPolygon(2,iEdgeOnCell,iCell) = sign(1.0_RKIND,normalGreatCircle(3))&
               * sqrt(1.0_RKIND - max(min(normalVectorPolygon(1,iEdgeOnCell,iCell),1.0_RKIND),-1.0_RKIND)**2)

       enddo ! iEdgeOnCell
       
    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_spherical_triangle_metric
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine normal_vectors_spherical_triangle_metric(&
       mesh, &
       rotateCartesianGrid, &
       removeMetricTerms, &
       normalVectorTriangle, &
       interiorVertex, &
       latVertexRotatedOut)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input: 
         removeMetricTerms !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle !< Output: 

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input: 

    real(kind=RKIND), dimension(:), optional, intent(out) :: &
         latVertexRotatedOut !< Output: 

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
         latVertexRotated, &
         lonVertexRotated, &
         normalGreatCircleNorm

    real(kind=RKIND), dimension(3,3) :: &
         yRotationMatrix, &
         zRotationMatrix

    integer, pointer :: &
         nVerticesSolve, &
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

    yRotationMatrix = 0.0_RKIND
    zRotationMatrix = 0.0_RKIND

    yRotationMatrix(2,2) = 1.0_RKIND
    zRotationMatrix(3,3) = 1.0_RKIND

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "zEdge", zEdge)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    do iVertex = 1, nVerticesSolve

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (present(latVertexRotatedOut)) latVertexRotatedOut(iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          ! rotate the cell centre to rotated geographical grid
          call cice_grid_rotation_forward(&
               vertexRotated(1), vertexRotated(2), vertexRotated(3), &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               rotateCartesianGrid)
 
          ! calculate lon and lat of the cell centre in the rotated coordinate system
          lonVertexRotated = atan2(vertexRotated(2), vertexRotated(1))
          latVertexRotated = asin(vertexRotated(3) / sphere_radius)
          
          ! rotate the cell centre onto the rotated equator
          if (removeMetricTerms) then
             yRotationMatrix(1,1) =  cos(latVertexRotated)
             yRotationMatrix(1,3) =  sin(latVertexRotated)
             yRotationMatrix(3,1) = -sin(latVertexRotated)
             yRotationMatrix(3,3) =  cos(latVertexRotated)
             
             zRotationMatrix(1,1) =  cos(-lonVertexRotated)
             zRotationMatrix(1,2) = -sin(-lonVertexRotated)
             zRotationMatrix(2,1) =  sin(-lonVertexRotated)
             zRotationMatrix(2,2) =  cos(-lonVertexRotated)
          else
             yRotationMatrix(1,1) = 1.0_RKIND
             yRotationMatrix(1,3) = 0.0_RKIND
             yRotationMatrix(3,1) = 0.0_RKIND
             yRotationMatrix(3,3) = 1.0_RKIND
             
             zRotationMatrix(1,1) = 1.0_RKIND
             zRotationMatrix(1,2) = 0.0_RKIND
             zRotationMatrix(2,1) = 0.0_RKIND
             zRotationMatrix(2,2) = 1.0_RKIND
          endif

          ! this should be (r, 0, 0)
          !vertexEquator = matmul(yRotationMatrix, matmul(zRotationMatrix, vertexRotated))
          
          do iVertexDegree = 1, vertexDegree
             
             iEdge = edgesOnVertex(iVertexDegree,iVertex)

             ! vector in great circle plane
             iCell1 = cellsOnEdge(1,iEdge)
             iCell2 = cellsOnEdge(2,iEdge)

             ! perform grid rotations
             call cice_grid_rotation_forward(&
                  edgeRotated(1), edgeRotated(2), edgeRotated(3), &
                  xEdge(iEdge),   yEdge(iEdge),   zEdge(iEdge), &
                  rotateCartesianGrid)

             call cice_grid_rotation_forward(&
                  cellRotated1(1), cellRotated1(2), cellRotated1(3), &
                  xCell(iCell1),   yCell(iCell1),   zCell(iCell1), &
                  rotateCartesianGrid)

             call cice_grid_rotation_forward(&
                  cellRotated2(1), cellRotated2(2), cellRotated2(3), &
                  xCell(iCell2),   yCell(iCell2),   zCell(iCell2), &
                  rotateCartesianGrid)

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
             
             call cice_dot_product_3space(normalVectorTriangle(1,iVertexDegree,iVertex), &
                                          normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                          eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangle(2,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3)) &
                  * sqrt(1.0_RKIND - max(min(normalVectorTriangle(1,iVertexDegree,iVertex),1.0_RKIND),-1.0_RKIND)**2)

          enddo ! iVertexDegree

          if (present(latVertexRotatedOut)) latVertexRotatedOut(iVertex) = latVertexRotated

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle_metric!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_spherical_triangle
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine normal_vectors_spherical_triangle(&
       mesh, &
       rotateCartesianGrid, &
       normalVectorTriangle, &
       interiorVertex)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle !< Output: 

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input: 

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
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

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
             call cice_grid_rotation_forward(&
                  xEdge0,       yEdge0,       zEdge0, &
                  xEdge(iEdge), yEdge(iEdge), zEdge(iEdge), &
                  rotateCartesianGrid)

             call cice_grid_rotation_forward(&
                  xCell1,        yCell1,        zCell1, &
                  xCell(iCell1), yCell(iCell1), zCell(iCell1), &
                  rotateCartesianGrid)

             call cice_grid_rotation_forward(&
                  xCell2,        yCell2,        zCell2, &
                  xCell(iCell2), yCell(iCell2), zCell(iCell2), &
                  rotateCartesianGrid)

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
             
             call cice_dot_product_3space(normalVectorTriangle(1,iVertexDegree,iVertex), &
                                          normalGreatCircle(1), normalGreatCircle(2), normalGreatCircle(3), &
                                          eastwardsVector(1),   eastwardsVector(2),   eastwardsVector(3))
          
             normalVectorTriangle(2,iVertexDegree,iVertex) = sign(1.0_RKIND,normalGreatCircle(3))&
                  * sqrt(1.0_RKIND - min(normalVectorTriangle(1,iVertexDegree,iVertex),1.0_RKIND)**2)
             
          enddo ! iVertexDegree

       endif ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_dot_product_3space
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_dot_product_3space(dot, x1, y1, z1, x2, y2, z2)!{{{

    real(kind=RKIND), intent(out) :: &
         dot !< Output: 

    real(kind=RKIND), intent(in) :: &
         x1, y1, z1, & !< Input: 
         x2, y2, z2    !< Input: 

    dot = x1 * x2 + y1 * y2 + z1 * z2

  end subroutine cice_dot_product_3space!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cross_product_3space
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cross_product_3space(cp1, cp2, cp3, x1, y1, z1, x2, y2, z2)!{{{

    real(kind=RKIND), intent(out) :: &
         cp1, cp2, cp3 !< Output: 

    real(kind=RKIND), intent(in) :: &
         x1, y1, z1, & !< Input: 
         x2, y2, z2    !< Input: 

    cp1 = y1 * z2 - z1 * y2
    cp2 = z1 * x2 - x1 * z2
    cp3 = x1 * y2 - y1 * x2

  end subroutine cross_product_3space!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_spherical_polygon_2
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine normal_vectors_spherical_polygon_2(&
       mesh, &
       rotateCartesianGrid, &
       normalVectorPolygon, &
       latCellRotated)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorPolygon !< Output: TWO maxEdges nCells

    real(kind=RKIND), dimension(:), intent(out) :: &
         latCellRotated !< Output: 
    
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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    ! loop over all cells
    do iCell = 1, nCells

       ! get rotated three space positions
       call cice_grid_rotation_forward(&
            xCell0,       yCell0,       zCell0, &
            xCell(iCell), yCell(iCell), zCell(iCell), &
            rotateCartesianGrid)

       ! get the rotated latitude
       latCellRotated(iCell) = asin(zCell0 / sphere_radius)

       ! loop over edges of cell
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          ! get the polygon side normal unit vector
          call get_polygon_side_normal_3D_vector(mesh, iCell, iEdgeOnCell, xCell0, yCell0, zCell0, rotateCartesianGrid, normalVector3D)

          ! project the 3D vector onto the local 2D plane
          call cice_project_3D_vector_onto_local_2D(normalVectorPolygon(:,iEdgeOnCell,iCell), normalVector3D, xCell0, yCell0, zCell0)

          ! normalize the projected vector
          vectorMagnitude = sqrt(normalVectorPolygon(1,iEdgeOnCell,iCell)**2 + normalVectorPolygon(2,iEdgeOnCell,iCell)**2)
          normalVectorPolygon(1,iEdgeOnCell,iCell) = normalVectorPolygon(1,iEdgeOnCell,iCell) / vectorMagnitude
          normalVectorPolygon(2,iEdgeOnCell,iCell) = normalVectorPolygon(2,iEdgeOnCell,iCell) / vectorMagnitude

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine normal_vectors_spherical_polygon_2!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  normal_vectors_spherical_triangle_2
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------
  
  subroutine normal_vectors_spherical_triangle_2(&
       mesh, &
       rotateCartesianGrid, &
       normalVectorTriangle, &
       latVertexRotated, &
       interiorVertex)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         normalVectorTriangle !< Output: TWO vertexDegree nVertices

    real(kind=RKIND), dimension(:), intent(out) :: &
         latVertexRotated !< Output: 

    integer, dimension(:), intent(in) :: &
         interiorVertex !< Input: 

    integer :: &
         iVertex, &
         iVertexDegree

    real(kind=RKIND), dimension(3) :: &
         normalVector3D

    real(kind=RKIND) :: &
         xVertex0, yVertex0, zVertex0, &
         vectorMagnitude

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
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    ! loop over interior vertices
    do iVertex = 1, nVertices

       normalVectorTriangle(:,:,iVertex) = 0.0_RKIND

       if (interiorVertex(iVertex) == 1) then

          ! get rotated three space positions
          call cice_grid_rotation_forward(&
               xVertex0,         yVertex0,         zVertex0, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               rotateCartesianGrid)

          ! get the rotated latitude
          latVertexRotated(iVertex) = asin(zVertex0 / sphere_radius)
          
          ! loop over edges of triangle
          do iVertexDegree = 1, vertexDegree

             ! get the triangle side normal unit vector
             call get_triangle_side_normal_3D_vector(mesh, iVertex, iVertexDegree, xVertex0, yVertex0, zVertex0, rotateCartesianGrid, normalVector3D)

             ! project the 3D vector onto the local 2D plane
             call cice_project_3D_vector_onto_local_2D(normalVectorTriangle(:,iVertexDegree,iVertex), normalVector3D, xVertex0, yVertex0, zVertex0)

             ! normalize the projected vector
             vectorMagnitude = sqrt(normalVectorTriangle(1,iVertexDegree,iVertex)**2 + normalVectorTriangle(2,iVertexDegree,iVertex)**2)
             normalVectorTriangle(1,iVertexDegree,iVertex) = normalVectorTriangle(1,iVertexDegree,iVertex) / vectorMagnitude
             normalVectorTriangle(2,iVertexDegree,iVertex) = normalVectorTriangle(2,iVertexDegree,iVertex) / vectorMagnitude

          enddo ! iVertexDegree

       end if ! interiorVertex

    enddo ! iVertex

  end subroutine normal_vectors_spherical_triangle_2!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_project_3D_vector_onto_local_2D
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_project_3D_vector_onto_local_2D(&
       normalVector2D, &
       normalVector3D, &
       xPoint, &
       yPoint, &
       zPoint)!{{{

    real(kind=RKIND), dimension(2), intent(out) :: &
         normalVector2D !< Output: 

    real(kind=RKIND), dimension(3), intent(in) :: &
         normalVector3D !< Input: 

    real(kind=RKIND), intent(in) :: &
         xPoint, & !< Input: 
         yPoint, & !< Input: 
         zPoint    !< Input: 

    real(kind=RKIND), dimension(3) :: &
         unitVectorEast, &
         unitVectorNorth

    ! calculate the local eastern and northern unit vectors at the point
    call local_eastern_and_northern_unit_vectors(&
         xPoint, yPoint, zPoint, &
         unitVectorEast, &
         unitVectorNorth)

    ! project normal 3D vector onto U 2D direction - eastwards
    call cice_dot_product_3space(&
         normalVector2D(1), &
         normalVector3D(1), normalVector3D(2), normalVector3D(3), &
         unitVectorEast(1), unitVectorEast(2), unitVectorEast(3))

    ! project normal 3D vector onto V 2D direction - northwards
    call cice_dot_product_3space(&
         normalVector2D(2), &
         normalVector3D(1),  normalVector3D(2),  normalVector3D(3), &
         unitVectorNorth(1), unitVectorNorth(2), unitVectorNorth(3))

  end subroutine cice_project_3D_vector_onto_local_2D!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_polygon_side_normal_3D_vector
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine get_polygon_side_normal_3D_vector(&
       mesh, &
       iCell, &
       iEdgeOnCell, &
       xCell0, &
       yCell0, &
       zCell0, &
       rotateCartesianGrid, &
       normalVector3D)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    integer, intent(in) :: &
         iCell, &    !< Input: 
         iEdgeOnCell !< Input: 

    real(kind=RKIND), intent(in) :: &
         xCell0, & !< Input: 
         yCell0, & !< Input: 
         zCell0    !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

    real(kind=RKIND), dimension(3), intent(out) :: &
         normalVector3D !< Output: 

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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

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
       call cice_grid_rotation_forward(&
            xCell2,        yCell2,        zCell2, &
            xCell(iCell2), yCell(iCell2), zCell(iCell2), &
            rotateCartesianGrid)

       normalVector3D(1) = xCell2 - xCell0
       normalVector3D(2) = yCell2 - yCell0
       normalVector3D(3) = zCell2 - zCell0

    else

       ! at domain edge

       iEdge = edgesOnCell(iEdgeOnCell, iCell)

       ! get rotated three space positions
       call cice_grid_rotation_forward(&
            xEdge0,       yEdge0,       zEdge0, &
            xEdge(iEdge), yEdge(iEdge), zEdge(iEdge), &
            rotateCartesianGrid)

       normalVector3D(1) = xEdge0 - xCell0
       normalVector3D(2) = yEdge0 - yCell0
       normalVector3D(3) = zEdge0 - zCell0

    endif

  end subroutine get_polygon_side_normal_3D_vector!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_triangle_side_normal_3D_vector
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine get_triangle_side_normal_3D_vector(&
       mesh, &
       iVertex, &
       iVertexDegree, &
       xVertex0, &
       yVertex0, &
       zVertex0, &
       rotateCartesianGrid, &
       normalVector3D)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    integer, intent(in) :: &
         iVertex, &    !< Input: 
         iVertexDegree !< Input: 

    real(kind=RKIND), intent(in) :: &
         xVertex0, & !< Input: 
         yVertex0, & !< Input: 
         zVertex0    !< Input: 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

    real(kind=RKIND), dimension(3), intent(out) :: &
         normalVector3D !< Output: 

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
    call cice_grid_rotation_forward(&
         xVertex2,          yVertex2,          zVertex2, &
         xVertex(iVertex2), yVertex(iVertex2), zVertex(iVertex2), &
         rotateCartesianGrid)
    
    normalVector3D(1) = xVertex2 - xVertex0
    normalVector3D(2) = yVertex2 - yVertex0
    normalVector3D(3) = zVertex2 - zVertex0

  end subroutine get_triangle_side_normal_3D_vector!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  local_eastern_and_northern_unit_vectors
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine local_eastern_and_northern_unit_vectors(&
       xPoint, &
       yPoint, &
       zPoint, &
       unitVectorEast, &
       unitVectorNorth)!{{{

    real(kind=RKIND), intent(in) :: &
         xPoint, & !< Input: 
         yPoint, & !< Input: 
         zPoint    !< Input: 

    real(kind=RKIND), dimension(3), intent(out) :: &
         unitVectorEast, & !< Output: 
         unitVectorNorth   !< Output: 

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

       if (zPoint < 0.0_RKIND) then

          ! have calculated due south vector if in southern hemisphere
          unitVectorNorth(1) = -unitVectorNorth(1)
          unitVectorNorth(2) = -unitVectorNorth(2)
          unitVectorNorth(3) = -unitVectorNorth(3)
       
       endif

    else   

       unitVectorNorth(1) = 0.0_RKIND
       unitVectorNorth(2) = 0.0_RKIND
       unitVectorNorth(3) = 1.0_RKIND

    endif

  end subroutine local_eastern_and_northern_unit_vectors!}}}

!-----------------------------------------------------------------------
! rotated lat-lon grid
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_grid_rotation_forward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_grid_rotation_forward(xp, yp, zp, x, y, z, rotateCartesianGrid)!{{{

    ! rotate xyz coordinates from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         xp, & !< Output: x position of point on rotated grid
         yp, & !< Output: y position of point on rotated grid
         zp    !< Output: z position of point on rotated grid

    real(kind=RKIND), intent(in) :: &
         x,  & !< Input: x position of point on geographical grid
         y,  & !< Input: y position of point on geographical grid
         z     !< Input: z position of point on geographical grid

    logical, intent(in) :: &
         rotateCartesianGrid

    if (rotateCartesianGrid) then

       xp = -z
       yp = y
       zp = x

    else

       xp = x
       yp = y
       zp = z

    endif

  end subroutine cice_grid_rotation_forward!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  grid_rotation_backward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine grid_rotation_backward(x, y, z, xp, yp, zp, rotateCartesianGrid)!{{{

    ! rotate xyz coordinates from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         x,  & !< Output: x position of point on geographical grid
         y,  & !< Output: y position of point on geographical grid
         z     !< Output: z position of point on geographical grid

    real(kind=RKIND), intent(in) :: &
         xp, & !< Input: x position of point on rotated grid
         yp, & !< Input: y position of point on rotated grid
         zp    !< Input: z position of point on rotated grid

    logical, intent(in) :: &
         rotateCartesianGrid

    if (rotateCartesianGrid) then

       x = zp
       y = yp
       z = -xp

    else

       x = xp
       y = yp
       z = zp

    endif

  end subroutine grid_rotation_backward!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  xyz_vector_rotation_forward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine xyz_vector_rotation_forward(vxp, vyp, vzp, vx, vy, vz, rotateCartesianGrid)!{{{

    ! rotate a xyz vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         vxp, & !< Output: x component of velocity vector on rotated grid
         vyp, & !< Output: y component of velocity vector on rotated grid
         vzp    !< Output: z component of velocity vector on rotated grid

    real(kind=RKIND), intent(in) :: &
         vx,  & !< Input: x component of velocity vector on geographical grid
         vy,  & !< Input: y component of velocity vector on geographical grid
         vz     !< Input: z component of velocity vector on geographical grid
    
    logical, intent(in) :: &
         rotateCartesianGrid

    if (rotateCartesianGrid) then
       
       vxp = -vz
       vyp =  vy
       vzp =  vx
       
    else

       vxp = vx
       vyp = vy
       vzp = vz

    endif

  end subroutine xyz_vector_rotation_forward!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  xyz_vector_rotation_backward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine xyz_vector_rotation_backward(vx, vy, vz, vxp, vyp, vzp, rotateCartesianGrid)!{{{

    ! rotate xyz vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         vx,  & !< Output: x component of velocity vector on geographical grid
         vy,  & !< Output: y component of velocity vector on geographical grid
         vz     !< Output: z component of velocity vector on geographical grid

    real(kind=RKIND), intent(in) :: &
         vxp, & !< Input: x component of velocity vector on rotated grid
         vyp, & !< Input: y component of velocity vector on rotated grid
         vzp    !< Input: z component of velocity vector on rotated grid

    logical, intent(in) :: &
         rotateCartesianGrid

    if (rotateCartesianGrid) then

       vx =  vzp
       vy =  vyp
       vz = -vxp

    else

       vx = vxp
       vy = vyp
       vz = vzp

    endif

  end subroutine xyz_vector_rotation_backward!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_latlon_vector_rotation_forward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_latlon_vector_rotation_forward(up, vp, u, v, lat, lon, x, y, z, r, rotateCartesianGrid)!{{{

    ! rotate a latlon vector from geographical grid to rotated grid with poles on real equator

    real(kind=RKIND), intent(out) :: &
         up,   & !< Output: u component of velocity in rotated grid
         vp      !< Output: v component of velocity in rotated grid

    real(kind=RKIND), intent(in) :: &
         u,    & !< Input: u component of velocity on geographical grid
         v,    & !< Input: v component of velocity on geographical grid
         lat,  & !< Input: latitude of point on geographical grid
         lon,  & !< Input: longitude of point on geographical grid
         x,    & !< Input: x position of point on geographical grid
         y,    & !< Input: y position of point on geographical grid
         z,    & !< Input: z position of point on geographical grid
         r       !< Input: radius of the earth 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

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
    call cice_grid_rotation_forward(xp, yp, zp, x, y, z, rotateCartesianGrid)

    ! calculate latitude and longitude of the point in the rotated grid
    call cice_latlon_from_xyz(latp, lonp, xp, yp, zp, r)

    ! convert lat lon vector to xyz vector on geographical grid
    call latlon_vector_to_xyz_vector(vx, vy, vz, u, v, lat, lon)

    ! perform rotation of geographical xyz vector to rotated grid
    call xyz_vector_rotation_forward(vxp, vyp, vzp, vx, vy, vz, rotateCartesianGrid)

    ! convert xyz vector to lat lon vector on rotated grid
    call xyz_vector_to_latlon_vector(up, vp, vxp, vyp, vzp, latp, lonp)

  end subroutine cice_latlon_vector_rotation_forward!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_latlon_vector_rotation_backward
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_latlon_vector_rotation_backward(u, v, up, vp, lat, lon, x, y, z, r, rotateCartesianGrid)!{{{

    ! rotate latlon vector from rotated grid with poles on real equator to geographical grid

    real(kind=RKIND), intent(out) :: &
         u,    & !< Output: u component of velocity on geographical grid
         v       !< Output: v component of velocity on geographical grid

    real(kind=RKIND), intent(in) :: &
         up,   & !< Input: u component of velocity in rotated grid
         vp,   & !< Input: v component of velocity in rotated grid
         lat,  & !< Input: latitude of point on geographical grid
         lon,  & !< Input: longitude of point on geographical grid
         x,    & !< Input: x position of point on geographical grid
         y,    & !< Input: y position of point on geographical grid
         z,    & !< Input: z position of point on geographical grid
         r       !< Input: radius of the earth 

    logical, intent(in) :: &
         rotateCartesianGrid !< Input: 

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
    call cice_grid_rotation_forward(xp, yp, zp, x, y, z, rotateCartesianGrid)

    ! calculate latitude and longitude of the point in the rotated grid
    call cice_latlon_from_xyz(latp, lonp, xp, yp, zp, r)

    ! convert lat lon vector to xyz vector on rotated grid
    call latlon_vector_to_xyz_vector(vxp, vyp, vzp, up, vp, latp, lonp)

    ! perform rotation of rotated xyz vector to geographical grid
    call xyz_vector_rotation_backward(vx, vy, vz, vxp, vyp, vzp, rotateCartesianGrid)

    ! convert xyz vector to lat lon vector on geographical grid
    call xyz_vector_to_latlon_vector(u, v, vx, vy, vz, lat, lon)
     
  end subroutine cice_latlon_vector_rotation_backward!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  latlon_vector_to_xyz_vector
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine latlon_vector_to_xyz_vector(vx, vy, vz, u, v, lat, lon)!{{{

    ! convert a latlon vector to a xyz vector

    real(kind=RKIND), intent(out) :: &
         vx,  & !< Output: x component of velocity vector
         vy,  & !< Output: y component of velocity vector
         vz     !< Output: z component of velocity vector

    real(kind=RKIND), intent(in) :: &
         u,   & !< Input: u component of velocity
         v,   & !< Input: v component of velocity
         lat, & !< Input: latitude of point 
         lon    !< Input: longitude of point 

    vx = (-u) * sin(lon) - v * sin(lat) * cos(lon)
    vy =   u  * cos(lon) - v * sin(lat) * sin(lon)
    vz =                   v * cos(lat)

  end subroutine latlon_vector_to_xyz_vector!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  xyz_vector_to_latlon_vector
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine xyz_vector_to_latlon_vector(u, v, vx, vy, vz, lat, lon)!{{{

    ! convert a xyz vector vector to a latlon vector

    real(kind=RKIND), intent(out) :: &
         u,   & !< Output: u component of velocity
         v      !< Output: v component of velocity

    real(kind=RKIND), intent(in) :: &
         vx,  & !< Input: x component of velocity vector
         vy,  & !< Input: y component of velocity vector
         vz,  & !< Input: z component of velocity vector
         lat, & !< Input: latitude of point 
         lon    !< Input: longitude of point 

    u = (-sin(lon)) * vx + &
          cos(lon)  * vy

    v = (-sin(lat)) * cos(lon) * vx + &
        (-sin(lat)) * sin(lon) * vy + &
          cos(lat)  * vz

  end subroutine xyz_vector_to_latlon_vector!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_latlon_from_xyz
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_latlon_from_xyz(lat, lon, x, y, z, r)!{{{

    ! given xyz coordinates determine the latitude and longitude

    real(kind=RKIND), intent(out) :: &
         lat, & !< Output: latitude of point 
         lon    !< Output: longitude of point 

    real(kind=RKIND), intent(in) :: &
         x,   & !< Input: x position of point
         y,   & !< Input: y position of point
         z,   & !< Input: z position of point
         r      !< Input: radius of earth

    lon = atan2(y, x)
    lat = asin(z/r)

  end subroutine cice_latlon_from_xyz!}}}

!-----------------------------------------------------------------------
! interpolaton
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_interpolate_cell_to_vertex
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_interpolate_cell_to_vertex(&
       mesh, &
       variableVertex, &
       variableCell)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         variableVertex !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         variableCell !< Input: 

    real(kind=RKIND) :: &
         totalArea

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell

    integer, pointer :: &
         nVerticesSolve, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell, &
         areaTriangle

    real(kind=RKIND), dimension(:,:), pointer :: &
         kiteAreasOnVertex

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVerticesSolve", nVerticesSolve)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)
    call MPAS_pool_get_array(mesh, "kiteAreasOnVertex", kiteAreasOnVertex)

    do iVertex = 1, nVerticesSolve
    
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

  end subroutine cice_interpolate_cell_to_vertex!}}}

!-----------------------------------------------------------------------

end module cice_mesh
