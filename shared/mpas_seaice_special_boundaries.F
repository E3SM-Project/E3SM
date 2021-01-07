!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_special_boundaries
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 23rd December 2020
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_special_boundaries

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       seaice_init_special_boundaries, &
       seaice_periodic_velocity

  logical, pointer :: &
       usePeriodicVelocities

  integer, parameter :: &
       VERTEX_BOUNDARY_NONE = 0, &
       VERTEX_BOUNDARY_PERIODIC = 1

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_special_boundaries
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 23rd December 2020
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_special_boundaries(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    call init_periodic_boundaries(domain)

  end subroutine seaice_init_special_boundaries

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_periodic_boundaries
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 23rd December 2020
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_periodic_boundaries(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(mpas_pool_type), pointer :: &
         specialBoundariesPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex

    integer, dimension(:), pointer :: &
         vertexBoundaryType, &
         vertexBoundarySource, &
         vertexBoundarySourceLocal, &
         indexToVertexID

    integer, dimension(:), allocatable :: &
         globalToLocalID

    integer, pointer :: &
         nVertices

    integer :: &
         iVertex

    call MPAS_pool_get_config(domain % configs, "config_use_periodic_velocities", usePeriodicVelocities)

    if (usePeriodicVelocities) then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_dimension(blockPtr % dimensions, "nVertices", nVertices)

          call MPAS_pool_get_subpool(blockPtr % structs, "special_boundaries", specialBoundariesPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)

          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundaryType", vertexBoundaryType)
          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundarySource", vertexBoundarySource)
          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundarySourceLocal", vertexBoundarySourceLocal)

          call MPAS_pool_get_array(meshPool, "xVertex", xVertex)
          call MPAS_pool_get_array(meshPool, "yVertex", yVertex)
          call MPAS_pool_get_array(meshPool, "indexToVertexID", indexToVertexID)

          allocate(globalToLocalID(nVertices))

          do iVertex = 1, nVertices

             globalToLocalID(indexToVertexID(iVertex)) = iVertex

          enddo ! iVertex

          do iVertex = 1, nVertices

             vertexBoundarySourceLocal(iVertex) = globalToLocalID(vertexBoundarySource(iVertex))

          enddo ! iVertex

          deallocate(globalToLocalID)

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine init_periodic_boundaries

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_periodic_velocity
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 23rd December 2020
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_periodic_velocity(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(mpas_pool_type), pointer :: &
         velocitySolverPool, &
         specialBoundariesPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         xVertex, &
         yVertex

    integer, dimension(:), pointer :: &
         vertexBoundaryType, &
         vertexBoundarySource, &
         indexToVertexID

    integer, pointer :: &
         nVertices

    integer :: &
         iVertex, &
         iVertexSource

    if (usePeriodicVelocities) then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_dimension(blockPtr % dimensions, "nVertices", nVertices)

          call MPAS_pool_get_subpool(blockPtr % structs, "special_boundaries", specialBoundariesPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)

          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundaryType", vertexBoundaryType)
          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundarySourceLocal", vertexBoundarySource)

          call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

          call MPAS_pool_get_array(meshPool, "xVertex", xVertex)
          call MPAS_pool_get_array(meshPool, "yVertex", yVertex)
          call MPAS_pool_get_array(meshPool, "indexToVertexID", indexToVertexID)

          do iVertex = 1, nVertices

             if (vertexBoundaryType(iVertex) == VERTEX_BOUNDARY_PERIODIC) then

                iVertexSource = vertexBoundarySource(iVertex)

                uVelocity(iVertex) = uVelocity(iVertexSource)
                vVelocity(iVertex) = vVelocity(iVertexSource)

             endif

          enddo ! iVertex

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine seaice_periodic_velocity

  !-----------------------------------------------------------------------

end module seaice_special_boundaries
