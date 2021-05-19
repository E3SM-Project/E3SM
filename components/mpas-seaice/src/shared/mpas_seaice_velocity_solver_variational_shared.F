!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_variational_shared
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 24 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_variational_shared

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       seaice_calc_local_coords, &
       seaice_calc_variational_metric_terms, &
       seaice_wrapped_index

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_calc_local_coords
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_calc_local_coords(&
       xLocal, &
       yLocal, &
       nCells, &
       nEdgesOnCell, &
       verticesOnCell, &
       xVertex, &
       yVertex, &
       zVertex, &
       xCell, &
       yCell, &
       zCell, &
       rotateCartesianGrid, &
       onASphere)!{{{

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex, & !< Input:
         xCell, &   !< Input:
         yCell, &   !< Input:
         zCell      !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         onASphere              !< Input:

    if (onASphere) then
       call calc_local_coords_spherical(&
            xLocal, &
            yLocal, &
            nCells, &
            nEdgesOnCell, &
            verticesOnCell, &
            xVertex, &
            yVertex, &
            zVertex, &
            xCell, &
            yCell, &
            zCell, &
            rotateCartesianGrid)
    else
       call calc_local_coords_planar(&
            xLocal, &
            yLocal, &
            nCells, &
            nEdgesOnCell, &
            verticesOnCell, &
            xVertex, &
            yVertex, &
            zVertex, &
            xCell, &
            yCell, &
            zCell)
    endif

  end subroutine seaice_calc_local_coords!}}}

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
       xLocal, &
       yLocal, &
       nCells, &
       nEdgesOnCell, &
       verticesOnCell, &
       xVertex, &
       yVertex, &
       zVertex, &
       xCell, &
       yCell, &
       zCell)!{{{

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex, & !< Input:
         xCell, &   !< Input:
         yCell, &   !< Input:
         zCell      !< Input:

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell

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
       xLocal, &
       yLocal, &
       nCells, &
       nEdgesOnCell, &
       verticesOnCell, &
       xVertex, &
       yVertex, &
       zVertex, &
       xCell, &
       yCell, &
       zCell, &
       rotateCartesianGrid)!{{{

    use seaice_mesh, only: &
         seaice_project_3D_vector_onto_local_2D, &
         seaice_grid_rotation_forward

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex, & !< Input:
         xCell, &   !< Input:
         yCell, &   !< Input:
         zCell      !< Input:

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

    real(kind=RKIND) :: &
         xCellRotated, &
         yCellRotated, &
         zCellRotated

    do iCell = 1, nCells

       call seaice_grid_rotation_forward(&
            xCellRotated, yCellRotated, zCellRotated, &
            xCell(iCell), yCell(iCell), zCell(iCell), &
            rotateCartesianGrid)

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          call seaice_grid_rotation_forward(&
               normalVector3D(1), normalVector3D(2), normalVector3D(3), &
               xVertex(iVertex),  yVertex(iVertex),  zVertex(iVertex), &
               rotateCartesianGrid)

          call seaice_project_3D_vector_onto_local_2D(&
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
!  seaice_calc_variational_metric_terms
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_calc_variational_metric_terms(&
       tanLatVertexRotatedOverRadius, &
       nVertices, &
       xVertex, &
       yVertex, &
       zVertex, &
       sphereRadius, &
       rotateCartesianGrid, &
       includeMetricTerms)

    use seaice_mesh, only: &
         seaice_grid_rotation_forward

    real(kind=RKIND), dimension(:), intent(out) :: &
         tanLatVertexRotatedOverRadius !< Output:

    integer, intent(in) :: &
         nVertices !< Input:

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex    !< Input:

    real(kind=RKIND), pointer :: &
         sphereRadius !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         includeMetricTerms     !< Input:

    integer :: &
         iVertex

    real(kind=RKIND) :: &
         xVertexRotated, &
         yVertexRotated, &
         zVertexRotated, &
         latVertexRotated

    if (includeMetricTerms) then

       do iVertex = 1, nVertices

          call seaice_grid_rotation_forward(&
               xVertexRotated,   yVertexRotated,   zVertexRotated, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               rotateCartesianGrid)

          latVertexRotated = asin(zVertexRotated / sphereRadius)

          tanLatVertexRotatedOverRadius(iVertex) = tan(latVertexRotated) / sphereRadius

       enddo ! iVertex

    else

       do iVertex = 1, nVertices

          tanLatVertexRotatedOverRadius(iVertex) = 0.0_RKIND

       enddo ! iVertex

    endif

  end subroutine seaice_calc_variational_metric_terms

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_wrapped_index
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  function seaice_wrapped_index(&
       input, &
       nelements) &
       result(output)!{{{

    integer, intent(in) :: &
         input, &  !< Input:
         nelements !< Input:

    integer :: output

    output = modulo(input - 1, nelements) + 1

  end function seaice_wrapped_index!}}}

!-----------------------------------------------------------------------

end module seaice_velocity_solver_variational_shared
