!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_variational_shared
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 24 October 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_variational_shared

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       cice_calc_local_coords, &
       cice_calc_variational_metric_terms, &
       cice_wrapped_index

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_calc_local_coords
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_calc_local_coords(&
       mesh, &
       xLocal, &
       yLocal, &
       rotateCartesianGrid)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    logical, intent(in) :: &
         rotateCartesianGrid !< Input:

    logical, pointer :: &
         on_a_sphere

    call MPAS_pool_get_config(mesh, "on_a_sphere", on_a_sphere)

    if (on_a_sphere) then
       call calc_local_coords_spherical(&
            mesh, &
            xLocal, &
            yLocal, &
            rotateCartesianGrid)
    else
       call calc_local_coords_planar(&
            mesh, &
            xLocal, &
            yLocal)
    endif

  end subroutine cice_calc_local_coords!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calc_local_coords_planar
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calc_local_coords_planar(&
       mesh, &
       xLocal, &
       yLocal)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         xCell, &
         yCell

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)

    do iCell = 1, nCells

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          xLocal(iVertexOnCell,iCell) = xVertex(iVertex) - xCell(iCell)
          yLocal(iVertexOnCell,iCell) = yVertex(iVertex) - yCell(iCell)

       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine calc_local_coords_planar!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calc_local_coords_spherical
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine calc_local_coords_spherical(&
       mesh, &
       xLocal, &
       yLocal, &
       rotateCartesianGrid)!{{{

    use cice_mesh, only: &
         cice_project_3D_vector_onto_local_2D, &
         cice_grid_rotation_forward

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    logical, intent(in) :: &
         rotateCartesianGrid !< Input:

    real(kind=RKIND), dimension(3) :: &
         normalVector3D

    real(kind=RKIND), dimension(2) :: &
         normalVector2D

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex, &
         xCell, &
         yCell, &
         zCell

    real(kind=RKIND) :: &
         xCellRotated, &
         yCellRotated, &
         zCellRotated

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    do iCell = 1, nCells

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          call cice_grid_rotation_forward(&
               normalVector3D(1), normalVector3D(2), normalVector3D(3), &
               xVertex(iVertex),  yVertex(iVertex),  zVertex(iVertex), &
               rotateCartesianGrid)

          call cice_grid_rotation_forward(&
               xCellRotated, yCellRotated, zCellRotated, &
               xCell(iCell), yCell(iCell), zCell(iCell), &
               rotateCartesianGrid)

          call cice_project_3D_vector_onto_local_2D(&
               normalVector2D, &
               normalVector3D, &
               xCellRotated, &
               yCellRotated, &
               zCellRotated)

          xLocal(iVertexOnCell,iCell) = normalVector2D(1)
          yLocal(iVertexOnCell,iCell) = normalVector2D(2)

       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine calc_local_coords_spherical!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_calc_variational_metric_terms
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_calc_variational_metric_terms(&
       mesh, &
       tanLatVertexRotatedOverRadius, &
       rotateCartesianGrid, &
       includeMetricTerms)

    use cice_mesh, only: &
         cice_grid_rotation_forward

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         tanLatVertexRotatedOverRadius !< Output:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         includeMetricTerms     !< Input:

    integer, pointer :: &
         nVertices

    integer :: &
         iVertex

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex

    real(kind=RKIND), pointer :: &
         sphere_radius

    real(kind=RKIND) :: &
         xVertexRotated, &
         yVertexRotated, &
         zVertexRotated, &
         latVertexRotated

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)

    if (includeMetricTerms) then

       do iVertex = 1, nVertices

          call cice_grid_rotation_forward(&
               xVertexRotated,   yVertexRotated,   zVertexRotated, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               rotateCartesianGrid)

          latVertexRotated = asin(zVertexRotated / sphere_radius)

          tanLatVertexRotatedOverRadius(iVertex) = tan(latVertexRotated) / sphere_radius

       enddo ! iVertex

    else

       do iVertex = 1, nVertices

          tanLatVertexRotatedOverRadius(iVertex) = 0.0_RKIND

       enddo ! iVertex

    endif

  end subroutine cice_calc_variational_metric_terms

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_wrapped_index
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  function cice_wrapped_index(&
       input, &
       nelements) &
       result(output)!{{{

    integer, intent(in) :: &
         input, &  !< Input:
         nelements !< Input:

    integer :: output

    output = modulo(input - 1, nelements) + 1

  end function cice_wrapped_index!}}}

!-----------------------------------------------------------------------

end module cice_velocity_solver_variational_shared
