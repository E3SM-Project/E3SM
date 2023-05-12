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
       seaice_set_special_boundaries_velocity, &
       seaice_set_special_boundaries_velocity_masks, &
       seaice_set_special_boundaries_tracers

  logical, pointer :: &
       useSpecialBoundariesVelocity, &
       useSpecialBoundariesVelocityMasks, &
       useSpecialBoundariesTracers

  ! vertex flags
  integer, parameter :: &
       VELOCITY_BOUNDARY_NONE     = 0, &
       VELOCITY_BOUNDARY_PERIODIC = 1, &
       VELOCITY_BOUNDARY_REVERSE  = 2, &
       VELOCITY_BOUNDARY_ZERO     = 3

  ! cell flags
  integer, parameter :: &
       TRACER_BOUNDARY_NONE = 0, &
       TRACER_BOUNDARY_ZERO = 1, &
       TRACER_BOUNDARY_SET  = 2

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

    call init_special_boundaries_velocity(domain)

    call init_special_boundaries_tracers(domain)

  end subroutine seaice_init_special_boundaries

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_special_boundaries_velocity
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 23rd December 2020
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_special_boundaries_velocity(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(mpas_pool_type), pointer :: &
         specialBoundariesPool, &
         meshPool

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

    call MPAS_pool_get_config(domain % configs, "config_use_special_boundaries_velocity", useSpecialBoundariesVelocity)
    call MPAS_pool_get_config(domain % configs, "config_use_special_boundaries_velocity_masks", useSpecialBoundariesVelocityMasks)

    if (useSpecialBoundariesVelocity) then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_dimension(blockPtr % dimensions, "nVertices", nVertices)

          call MPAS_pool_get_subpool(blockPtr % structs, "special_boundaries", specialBoundariesPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)

          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundaryType", vertexBoundaryType)
          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundarySource", vertexBoundarySource)
          call MPAS_pool_get_array(specialBoundariesPool, "vertexBoundarySourceLocal", vertexBoundarySourceLocal)

          call MPAS_pool_get_array(meshPool, "indexToVertexID", indexToVertexID)

          allocate(globalToLocalID(nVertices))

          do iVertex = 1, nVertices

             globalToLocalID(indexToVertexID(iVertex)) = iVertex

          enddo ! iVertex

          do iVertex = 1, nVertices

             if (vertexBoundaryType(iVertex) /= VELOCITY_BOUNDARY_NONE) then

                vertexBoundarySourceLocal(iVertex) = globalToLocalID(vertexBoundarySource(iVertex))

             end if

          enddo ! iVertex

          deallocate(globalToLocalID)

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine init_special_boundaries_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_special_boundaries_tracers
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 16th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_special_boundaries_tracers(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(mpas_pool_type), pointer :: &
         specialBoundariesPool, &
         meshPool

    integer, dimension(:), pointer :: &
         tracerBoundaryType, &
         tracerBoundarySource, &
         tracerBoundarySourceLocal, &
         indexToCellID

    integer, dimension(:), allocatable :: &
         globalToLocalID

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    call MPAS_pool_get_config(domain % configs, "config_use_special_boundaries_tracers", useSpecialBoundariesTracers)

    if (useSpecialBoundariesTracers) then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_dimension(blockPtr % dimensions, "nCells", nCells)

          call MPAS_pool_get_subpool(blockPtr % structs, "special_boundaries", specialBoundariesPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)

          call MPAS_pool_get_array(specialBoundariesPool, "tracerBoundaryType", tracerBoundaryType)
          call MPAS_pool_get_array(specialBoundariesPool, "tracerBoundarySource", tracerBoundarySource)
          call MPAS_pool_get_array(specialBoundariesPool, "tracerBoundarySourceLocal", tracerBoundarySourceLocal)

          call MPAS_pool_get_array(meshPool, "indexToCellID", indexToCellID)

          allocate(globalToLocalID(nCells))

          do iCell = 1, nCells

             globalToLocalID(indexToCellID(iCell)) = iCell

          enddo ! iCell

          do iCell = 1, nCells

             if (tracerBoundaryType(iCell) /= TRACER_BOUNDARY_NONE .and. &
                 tracerBoundaryType(iCell) /= TRACER_BOUNDARY_ZERO) then

                tracerBoundarySourceLocal(iCell) = globalToLocalID(tracerBoundarySource(iCell))

             endif

          enddo ! iCell

          deallocate(globalToLocalID)

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine init_special_boundaries_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_set_special_boundaries_velocity
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 23rd December 2020
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_set_special_boundaries_velocity(domain)

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

    if (useSpecialBoundariesVelocity) then

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

          do iVertex = 1, nVertices

             if (vertexBoundaryType(iVertex) == VELOCITY_BOUNDARY_PERIODIC) then

                iVertexSource = vertexBoundarySource(iVertex)

                uVelocity(iVertex) = uVelocity(iVertexSource)
                vVelocity(iVertex) = vVelocity(iVertexSource)

             else if (vertexBoundaryType(iVertex) == VELOCITY_BOUNDARY_REVERSE) then

                iVertexSource = vertexBoundarySource(iVertex)

                uVelocity(iVertex) = -uVelocity(iVertexSource)
                vVelocity(iVertex) = -vVelocity(iVertexSource)

             else if (vertexBoundaryType(iVertex) == VELOCITY_BOUNDARY_ZERO) then

                uVelocity(iVertex) = 0.0
                vVelocity(iVertex) = 0.0

             endif

          enddo ! iVertex

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine seaice_set_special_boundaries_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_set_special_boundaries_velocity_masks
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 23rd December 2020
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_set_special_boundaries_velocity_masks(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(mpas_pool_type), pointer :: &
         velocitySolverPool, &
         specialBoundariesPool

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveStress, &
         solveVelocitySpecialBoundaries, &
         solveStressSpecialBoundaries

    integer, pointer :: &
         nCells, &
         nVertices

    integer :: &
         iCell, &
         iVertex

    if (useSpecialBoundariesVelocityMasks) then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_dimension(blockPtr % dimensions, "nCells", nCells)
          call MPAS_pool_get_dimension(blockPtr % dimensions, "nVertices", nVertices)

          call MPAS_pool_get_subpool(blockPtr % structs, "special_boundaries", specialBoundariesPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)

          call MPAS_pool_get_array(specialBoundariesPool, "solveVelocitySpecialBoundaries", solveVelocitySpecialBoundaries)
          call MPAS_pool_get_array(specialBoundariesPool, "solveStressSpecialBoundaries", solveStressSpecialBoundaries)

          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

          do iVertex = 1, nVertices
             solveVelocity(iVertex) = solveVelocitySpecialBoundaries(iVertex)
          enddo ! iVertex

          do iCell = 1, nCells
             solveStress(iCell) = solveStressSpecialBoundaries(iCell)
          enddo ! iCell

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine seaice_set_special_boundaries_velocity_masks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_set_special_boundaries_tracers
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 16th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_set_special_boundaries_tracers(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(mpas_pool_type), pointer :: &
         specialBoundariesPool, &
         tracersPool

    integer, dimension(:), pointer :: &
         tracerBoundaryType, &
         tracerBoundarySource

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, pointer :: &
         nCells

    integer :: &
         iCell, &
         iCellSource

    if (useSpecialBoundariesTracers) then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_dimension(blockPtr % dimensions, "nCells", nCells)

          call MPAS_pool_get_subpool(blockPtr % structs, "special_boundaries", specialBoundariesPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "tracers", tracersPool)

          call MPAS_pool_get_array(specialBoundariesPool, "tracerBoundaryType", tracerBoundaryType)
          call MPAS_pool_get_array(specialBoundariesPool, "tracerBoundarySourceLocal", tracerBoundarySource)

          call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
          call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
          call MPAS_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)

          do iCell = 1, nCells

             if (tracerBoundaryType(iCell) == TRACER_BOUNDARY_ZERO) then

                iceAreaCategory(:,:,iCell)    = 0.0_RKIND
                iceVolumeCategory(:,:,iCell)  = 0.0_RKIND
                snowVolumeCategory(:,:,iCell) = 0.0_RKIND

             else if (tracerBoundaryType(iCell) == TRACER_BOUNDARY_SET) then

                iCellSource = tracerBoundarySource(iCell)

                iceAreaCategory(:,:,iCell)    = iceAreaCategory(:,:,iCellSource)
                iceVolumeCategory(:,:,iCell)  = iceVolumeCategory(:,:,iCellSource)
                snowVolumeCategory(:,:,iCell) = snowVolumeCategory(:,:,iCellSource)

             endif

          enddo ! iCell

          blockPtr => blockPtr % next
       enddo

    endif

  end subroutine seaice_set_special_boundaries_tracers

  !-----------------------------------------------------------------------

end module seaice_special_boundaries
