module cice_advection

  use mpas_grid_types

  implicit none

  private
  public :: init_advection, &
            run_advection

contains

  !--------------------------------------------------------------------------
  
  subroutine init_advection(block, mesh)

    use mpas_tracer_advection_helpers, only: &
         mpas_initialize_deriv_two, &
         mpas_tracer_advection_coefficients

    use mpas_configure, only: &
         config_horiz_tracer_adv_order
    
    type (block_type), intent(inout) :: block
    type (mesh_type), intent(inout) :: mesh

    integer :: &
         err, &
         err1

    call mpas_initialize_deriv_two(&
         mesh, &
         mesh % derivTwo % array, &
         err)

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

  end subroutine init_advection

  !--------------------------------------------------------------------------

  subroutine run_advection(block, dt, nstep)

    type (block_type), intent(inout) :: block
    real (kind=RKIND), intent(in) :: dt
    integer, intent(in) :: nstep

    real(kind=RKIND), dimension(:,:), allocatable :: &
         normalThicknessFlux, &
         verticalVelocity, &
         layerThickness, &
         verticalCellSize, &
         tend_layerThickness

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         tendency

    integer :: &
         iCell, &
         iCellOnEdge, &
         iEdge, &
         iTracer, &
         nTracers

    nTracers = size(block % icestate % tracers % array, 1)

    allocate(normalThicknessFlux(1, block % mesh % nEdges))
    allocate(verticalVelocity(2, block % mesh % nCells))
    allocate(layerThickness(1, block % mesh % nCells))
    allocate(verticalCellSize(1, block % mesh % nCells))
    allocate(tend_layerThickness(1, block % mesh % nCells))

    allocate(tendency(nTracers, 1, block % mesh % nCells))

    ! edge velocity
    call edge_from_vertex_velocity(&
         block % mesh, &
         block % hexdyn % uVelocity % array, &
         block % hexdyn % vVelocity % array, &
         block % normal % normalVectorPolygon % array, &
         normalThicknessFlux(1,:))

    ! other quantities
    do iCell = 1, block % mesh % nCells

       verticalVelocity(:,iCell)    = 0.0_RKIND
       layerThickness(:,iCell)      = 1.0_RKIND
       verticalCellSize(:,iCell)    = 1.0_RKIND
       tend_layerThickness(:,iCell) = 0.0_RKIND
    
    enddo ! iCell

    ! run the advection routine
    call mpas_tracer_advection_std_tend(&
         block % icestate % tracers % array, &
         block % mesh % advCoefs % array, &
         block % mesh % advCoefs3rd % array, &
         block % mesh % nAdvCellsForEdge % array, &
         block % mesh % advCellsForEdge % array, &
         normalThicknessFlux, &
         verticalVelocity, &
         layerThickness, &
         verticalCellSize, &
         dt, &
         block % mesh, &
         tend_layerThickness, &
         tendency, &
         highOrderAdvectionMask_in=block % mesh % highOrderAdvectionMask % array)

    ! apply the advective changes
    do iCell = 1, block % mesh % nCells
       do iTracer = 1, nTracers

          block % icestate % tracers % array(iTracer,1,iCell) = block % icestate % tracers % array(iTracer,1,iCell) + &
               tendency(iTracer,1,iCell) * dt

       enddo ! iTracers
    enddo ! iCell

    ! cleanup
    deallocate(normalThicknessFlux)
    deallocate(verticalVelocity)
    deallocate(layerThickness)
    deallocate(verticalCellSize)
    deallocate(tend_layerThickness)

    deallocate(tendency)

  end subroutine run_advection

  !--------------------------------------------------------------------------

  subroutine edge_from_vertex_velocity(mesh, uVelocity, vVelocity, normalVectorPolygon, edgeVelocity)

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

end module cice_advection
