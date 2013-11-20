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

    integer :: &
         nVertLevels = , &
         iVertLevel, &
         iCell, &
         iEdge

    allocate(normalThicknessFlux(nVertLevels, mesh % nEdges))
    allocate(verticalVelocity(nVertLevels+1,  mesh % nCells))
    allocate(layerThickness(nVertLevels, mesh % nCells))
    allocate(verticalCellSize(nVertLevels, mesh % nCells))
    allocate(tend_layerThickness(nVertLevels, mesh % nCells))

    ! edge velocity
    do iEdge = 1, mesh % nEdges
       do iVertLevel = 1, nVertLevels

          normalThicknessFlux(iVertLevel,iEdge) = ??

       enddo ! iVertLevel
    enddo ! iEdge

    ! other quantities
    do iCell = 1, mesh % nCells
       do iVertLevel = 1, nVertLevels

          verticalVelocity(iVertLevel,iCell)    = 0.0_RKIND
          layerThickness(iVertLevel,iCell)      = 1.0_RKIND
          verticalCellSize(iVertLevel,iCell)    = 1.0_RKIND
          tend_layerThickness(iVertLevel,iCell) = 0.0_RKIND
    
       enddo ! iVertLevel

       verticalVelocity(nVertLevels+1,iCell) = 0.0_RKIND

    enddo ! iCell

    ! run the advection routine
    call mpas_tracer_advection_std_tend(&
         tracers, &
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
         tend, &
         highOrderAdvectionMask_in=block % mesh % highOrderAdvectionMask % array)

    ! apply the advective changes
    do iCell = 1, mesh % nCells
       do iVertLevel = 1, nVertLevels
          do iTracers = 1, nTracers

             tracers

          enddo ! iTracers
       enddo ! iVertLevel
    enddo ! iCell

    ! cleanup
    deallocate(normalThicknessFlux)
    deallocate(verticalVelocity)
    deallocate(layerThickness)
    deallocate(verticalCellSize)
    deallocate(tend_layerThickness)

  end subroutine run_advection

  !--------------------------------------------------------------------------

end module cice_advection
