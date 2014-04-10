module cice_advection

  use mpas_grid_types

  use cice_testing

  implicit none

  private
  public :: init_advection, &
            run_advection

contains

  !--------------------------------------------------------------------------
  
  subroutine init_advection(mesh)

    use mpas_tracer_advection_helpers, only: &
         mpas_initialize_deriv_two, &
         mpas_tracer_advection_coefficients

    use mpas_configure, only: &
         config_horiz_tracer_adv_order
    
    type (mesh_type), intent(inout) :: mesh

    integer :: &
         err, &
         err1

    integer :: &
         i1, i2, i3

    !do i1 = 1, mesh % maxEdges2
    !   do i2 = 1, 2
    !      do i3 = 1, mesh % nEdges
    !         write(*,*) i1, i2, i3, mesh % derivTwo % array(i1,i2,i3)
    !      enddo
    !   enddo
    !enddo

    !call mpas_initialize_deriv_two(&
    !     mesh, &
    !     mesh % derivTwo % array, &
    !     err)

    !FIFTEEN TWO nEdges
    !write(*,*) &
    !     size(mesh % derivTwo % array,1), &
    !     size(mesh % derivTwo % array,2), &
    !     size(mesh % derivTwo % array,3), &
    !     mesh % maxEdges2, &
    !     mesh % nEdges

    !do i1 = 1, mesh % maxEdges2
    !   do i2 = 1, 2
    !      do i3 = 1, mesh % nEdges
    !         write(*,*) i1, i2, i3, mesh % derivTwo % array(i1,i2,i3)
    !      enddo
    !   enddo
    !enddo

    call mpas_tracer_advection_coefficients(&
         mesh, &
         config_horiz_tracer_adv_order, &
         mesh % derivTwo % array, &
         mesh % advCoefs % array, &
         mesh % advCoefs3rd % array, &
         mesh % nAdvCellsForEdge % array, &
         mesh % advCellsForEdge % array, &
         err1, &
         mesh % maxLevelCell % array, &
         mesh % highOrderAdvectionMask % array, &
         mesh % boundaryCell % array)

    !stop

  end subroutine init_advection

  !--------------------------------------------------------------------------

  subroutine run_advection(block, dt, nstep)

    use cice_diagnostics, only: plot_pop_grid_cells, output_meridonal_strip

    type (block_type), intent(inout) :: block
    real (kind=RKIND), intent(in) :: dt
    integer, intent(in) :: nstep

    type (mesh_type), pointer :: mesh
    type (hexdyn_type), pointer :: hexdyn
    type (normal_type), pointer :: normal

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         tracers, &              ! (nTracers, nVertLevels, nCells)
         tracersOnArea

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         tendency                ! (nTracers, nVertLevels, nCells)

    real(kind=RKIND), dimension(:,:), allocatable :: &
         lowerLevelTracer,     & ! (nVertLevels, nCells)
         lowerLevelTracerPrev, & ! (nVertLevels, nCells)
         lowerLevelTracerNext, & ! (nVertLevels, nCells)
         lowerLevelTendency      ! (nVertLevels, nCells)

    real(kind=RKIND), dimension(:,:), allocatable :: &
         verticalVelocity,     & ! (nVertLevels+1, nEdges)
         verticalCellSize        ! (nVertLevels, nCells)

    integer :: &
         nTracers, &
         nVertLevels

    integer :: &
         iCell, &
         iVertLevel, &
         iTracer

    ! perform pointer assignments
    mesh    => block % mesh
    hexdyn  => block % hexdyn
    normal  => block % normal

    ! init non changing things
    call init_run_advection(&
         mesh, &
         verticalVelocity, &
         verticalCellSize)

    nVertLevels = 1

    allocate(lowerLevelTracer     (1,mesh % nCells))
    allocate(lowerLevelTracerPrev (1,mesh % nCells))
    allocate(lowerLevelTracerNext (1,mesh % nCells))
    allocate(lowerLevelTendency   (1,mesh % nCells))

    ! calculate edge velocity
    call edge_from_vertex_velocity(&
       mesh, &
       hexdyn % uVelocity % array, &
       hexdyn % vVelocity % array, &
       normal % normalVectorPolygon % array, &
       hexdyn % edgeVelocity % array)

    ! perform first level
    tracers => block % icestate % tracers % array
    nTracers = size(tracers, 1)
    allocate(tendency(nTracers,1,mesh % nCells))

    do iCell = 1, mesh % nCells
       do iVertLevel = 1, nVertLevels

          lowerLevelTracer     (iVertLevel, iCell) = 1.0_RKIND
          lowerLevelTracerPrev (iVertLevel, iCell) = 1.0_RKIND
          lowerLevelTendency   (iVertLevel, iCell) = 0.0_RKIND

          do iTracer = 1, nTracers
             
             tendency(iTracer, iVertLevel, iCell) = 0.0_RKIND

          enddo ! iTracer

       enddo ! iVertLevel
    enddo ! iCell

    ! perform the first level of advection
    call perform_advection(&
       mesh, &
       dt, &
       tracers, &
       tendency, &
       hexdyn % edgeVelocity % array, &
       lowerLevelTracer, &
       lowerLevelTracerPrev, &
       lowerLevelTendency, &
       verticalVelocity, &
       verticalCellSize)

    ! perform second level




    ! cleanup

    deallocate(lowerLevelTracer)
    deallocate(lowerLevelTracerPrev)
    deallocate(lowerLevelTendency)

    call cleanup_run_advection(&
         verticalVelocity, &
         verticalCellSize)

  end subroutine run_advection

  !--------------------------------------------------------------------------

  subroutine perform_advection(&
       mesh, &
       dt, &
       tracers, &
       tendency, &
       edgeVelocity, &
       lowerLevelTracer, &
       lowerLevelTracerPrev, &
       lowerLevelTendency, &
       verticalVelocity, &
       verticalCellSize)

    use mpas_tracer_advection_mono, only: &
         mpas_tracer_advection_mono_tend

    type (mesh_type), intent(in) :: mesh

    real(kind=RKIND), intent(in) :: &
         dt

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers                 ! (nTracers, nVertLevels, nCells)

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         tendency                ! (nTracers, nVertLevels, nCells)

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity            ! (nEdges)

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         lowerLevelTracer,     & ! (nVertLevels, nCells)
         lowerLevelTracerPrev, & ! (nVertLevels, nCells)
         lowerLevelTendency      ! (nVertLevels, nCells)

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity,     & ! (nVertLevels+1, nEdges)
         verticalCellSize        ! (nVertLevels, nCells)

    real(kind=RKIND), dimension(:,:), allocatable :: &
         edgeFlux                ! (nVertLevels, nEdges)

    integer :: &
         iEdge, &
         iCell, &
         iTracer, &
         iVertLevel, &
         iCell1, &
         iCell2, &
         iAdvectionCell

    integer :: &
         nTracers, &
         nVertLevels, &
         nAdvectionCells

    nTracers    = size(tracers, 1)
    nVertLevels = size(tracers, 2)

    ! calculate edgeVelocity * lowerLevelTracerAtEdge
    allocate(edgeFlux(nVertLevels, mesh % nEdges))

    do iEdge = 1, mesh % nEdges

       iCell1 = mesh % cellsOnEdge % array(1,iEdge)
       iCell2 = mesh % cellsOnEdge % array(2,iEdge)

       if (iCell1 > mesh % nCells .or. iCell2 > mesh % nCells) then

          ! edge of domain
          do iVertLevel = 1, nVertLevels
             edgeFlux(iVertLevel,iEdge) = 0.0_RKIND
          enddo ! iVertLevel

       else

          do iVertLevel = 1, nVertLevels
             edgeFlux(iVertLevel,iEdge) = edgeVelocity(iEdge) * &
                  0.5_RKIND * (lowerLevelTracer(iVertLevel,iCell1) + lowerLevelTracer(iVertLevel,iCell2))
          enddo ! iVertLevel

       endif

    enddo ! iEdge

    ! run the advection routine
    call mpas_tracer_advection_mono_tend(&
         tracers,                         & ! (nTracers, nVertLevels, nCells)
         mesh % advCoefs % array,         & ! (nAdvectionCells nEdges)
         mesh % advCoefs3rd % array,      & ! (nAdvectionCells nEdges)
         mesh % nAdvCellsForEdge % array, & ! (nEdges)
         mesh % advCellsForEdge % array,  & ! (nAdvectionCells nEdges)
         edgeFlux,                        & ! (nVertLevels, nEdges)
         verticalVelocity,                & ! (nVertLevels+1, nEdges)
         lowerLevelTracer,                & ! (nVertLevels, nCells)
         verticalCellSize,                & ! (nVertLevels, nCells)
         dt,                              & ! 
         mesh,                            & ! 
         lowerLevelTendency,              & ! (nVertLevels, nCells)
         tendency,                        & ! (nTracers, nVertLevels, nCells)
         highOrderAdvectionMask_in =      & ! (nVertLevels nEdges)
         mesh % highOrderAdvectionMask % array, &
         lvertical_in=.false.)

    ! apply the advective changes
    do iCell = 1, mesh % nCells
       do iVertLevel = 1, nVertLevels
          do iTracer = 1, nTracers

             tracers(iTracer,iVertLevel,iCell) = (tracers(iTracer,iVertLevel,iCell) * lowerLevelTracerPrev(iVertLevel,iCell) &
                                               + tendency(iTracer,iVertLevel,iCell) * dt) &
                                               / lowerLevelTracer(iVertLevel,iCell)

          enddo ! iTracer
       enddo ! iVertLevel
    enddo ! iCell

    deallocate(edgeFlux)

  end subroutine perform_advection

  !--------------------------------------------------------------------------

  subroutine edge_from_vertex_velocity(&
       mesh, &
       uVelocity, &
       vVelocity, &
       normalVectorPolygon, &
       edgeVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon

    real(kind=RKIND), dimension(:), intent(out) :: &
         edgeVelocity

    real(kind=RKIND) :: &
         uVelocityEdge, &
         vVelocityEdge

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge, &
         iVertexOnEdge, &
         iVertex

    ! loop over cells
    do iCell = 1, mesh % nCells

       ! loop over cell edges
       do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)

          iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

          ! determine if velocity points outwards
          if (mesh % cellsOnEdge % array(1,iEdge) == iCell) then

             ! find u,v velocity at edge
             uVelocityEdge = 0.0_RKIND
             vVelocityEdge = 0.0_RKIND

             do iVertexOnEdge = 1, 2

                iVertex = mesh % verticesOnEdge % array(iVertexOnEdge, iEdge)

                uVelocityEdge = uVelocityEdge + uVelocity(iVertex)
                vVelocityEdge = vVelocityEdge + vVelocity(iVertex)

             enddo ! iVertexOnEdge

             uVelocityEdge = uVelocityEdge / 2.0_RKIND
             vVelocityEdge = vVelocityEdge / 2.0_RKIND

             ! rotate u,v velocity to normal
             edgeVelocity(iEdge) = &
                  uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell)

          endif

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine edge_from_vertex_velocity

  !--------------------------------------------------------------------------

  subroutine init_run_advection(mesh, verticalVelocity, verticalCellSize)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), allocatable, intent(inout) :: &
         verticalVelocity, &
         verticalCellSize

    integer :: &
         iEdge, &
         iCell

    allocate(verticalVelocity(2, mesh % nEdges))
    allocate(verticalCellSize(1, mesh % nCells))

    do iEdge = 1, mesh % nEdges

       verticalVelocity(:,iEdge) = 0.0_RKIND       

    enddo ! iEdge

    do iCell = 1, mesh % nCells

       verticalCellSize(:,iCell) = 1.0_RKIND
    
    enddo ! iCell

  end subroutine init_run_advection

  !--------------------------------------------------------------------------

  subroutine cleanup_run_advection(verticalVelocity, verticalCellSize)

    real(kind=RKIND), dimension(:,:), allocatable, intent(inout) :: &
         verticalVelocity, &
         verticalCellSize

    deallocate(verticalVelocity)
    deallocate(verticalCellSize)

  end subroutine cleanup_run_advection

  !--------------------------------------------------------------------------


end module cice_advection
