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
       seaice_calc_variational_metric_terms, &
       seaice_cell_vertices_at_vertex

contains

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

  subroutine seaice_calc_variational_metric_terms(domain)

    use seaice_mesh, only: &
         seaice_grid_rotation_forward

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocityVariationalPool

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius

    integer, pointer :: &
         nVertices

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex

    real(kind=RKIND), pointer :: &
         sphereRadius

    logical, pointer :: &
         config_rotate_cartesian_grid, &
         config_include_metric_terms

    integer :: &
         iVertex

    real(kind=RKIND) :: &
         xVertexRotated, &
         yVertexRotated, &
         zVertexRotated, &
         latVertexRotated

    call MPAS_pool_get_config(domain % configs, "config_include_metric_terms", config_include_metric_terms)

    if (config_include_metric_terms) then

       call MPAS_pool_get_config(domain % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

          call MPAS_pool_get_dimension(meshPool, "nVertices", nVertices)

          call MPAS_pool_get_config(meshPool, "sphere_radius", sphereRadius)

          call MPAS_pool_get_array(meshPool, "xVertex", xVertex)
          call MPAS_pool_get_array(meshPool, "yVertex", yVertex)
          call MPAS_pool_get_array(meshPool, "zVertex", zVertex)

          call MPAS_pool_get_array(velocityVariationalPool, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)

          do iVertex = 1, nVertices

             call seaice_grid_rotation_forward(&
                  xVertexRotated,   yVertexRotated,   zVertexRotated, &
                  xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
                  config_rotate_cartesian_grid)

             latVertexRotated = asin(zVertexRotated / sphereRadius)

             tanLatVertexRotatedOverRadius(iVertex) = tan(latVertexRotated) / sphereRadius

          enddo ! iVertex

          blockPtr => blockPtr % next
       enddo

    else

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

          call MPAS_pool_get_dimension(meshPool, "nVertices", nVertices)

          call MPAS_pool_get_array(velocityVariationalPool, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)

          do iVertex = 1, nVertices

             tanLatVertexRotatedOverRadius(iVertex) = 0.0_RKIND

          enddo ! iVertex

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine seaice_calc_variational_metric_terms

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_cell_vertices_at_vertex
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_cell_vertices_at_vertex(domain)!{{{

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocityVariationalPool

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnVertex, &
         verticesOnCell

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         iVertexOnCell, &
         jVertex

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_dimension(meshPool, "nVertices", nVertices)
       call MPAS_pool_get_dimension(meshPool, "vertexDegree", vertexDegree)

       call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
       call MPAS_pool_get_array(meshPool, "verticesOnCell", verticesOnCell)
       call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)

       call MPAS_pool_get_array(velocityVariationalPool, "cellVerticesAtVertex", cellVerticesAtVertex)

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

       blockPtr => blockPtr % next
    enddo

  end subroutine seaice_cell_vertices_at_vertex!}}}

!-----------------------------------------------------------------------

end module seaice_velocity_solver_variational_shared
