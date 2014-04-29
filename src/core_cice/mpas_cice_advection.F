module cice_advection

  use mpas_grid_types

  implicit none

  private
  public :: &
       init_advection, &
       advection
  
  integer, parameter :: &
       nTracerVariables = 7

  ! tracerConnectivity(iTracer, 1) ! child tracer name in code
  ! tracerConnectivity(iTracer, 2) ! parent tracer name in code
  character(len=200), dimension(nTracerVariables, 2) :: &
       tracerConnectivity

contains

  !--------------------------------------------------------------------------
  ! Initialization
  !--------------------------------------------------------------------------
  
  subroutine init_advection(block)

    use mpas_tracer_advection_helpers, only: &
         mpas_tracer_advection_coefficients

    type(block_type), intent(inout) :: &
         block

    type (MPAS_pool_type), pointer :: &
         mesh, &
         configs

    integer, pointer :: &
         config_horiz_tracer_adv_order

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         derivTwo

    real(kind=RKIND), dimension(:,:), pointer :: &
         advCoefs, &
         advCoefs3rd

    integer, dimension(:), pointer :: &
         nAdvCellsForEdge, &
         maxLevelCell

    integer, dimension(:,:), pointer :: &
         advCellsForEdge, &
         highOrderAdvectionMask, &
         boundaryCell

    integer :: &
         err

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    configs => block % configs
    call MPAS_pool_get_config(configs, "config_horiz_tracer_adv_order", config_horiz_tracer_adv_order)

    call MPAS_pool_get_array(mesh, "derivTwo", derivTwo)
    call MPAS_pool_get_array(mesh, "advCoefs", advCoefs)
    call MPAS_pool_get_array(mesh, "advCoefs3rd", advCoefs3rd)
    call MPAS_pool_get_array(mesh, "nAdvCellsForEdge", nAdvCellsForEdge)
    call MPAS_pool_get_array(mesh, "advCellsForEdge", advCellsForEdge)
    call MPAS_pool_get_array(mesh, "maxLevelCell", maxLevelCell)
    call MPAS_pool_get_array(mesh, "highOrderAdvectionMask", highOrderAdvectionMask)
    call MPAS_pool_get_array(mesh, "boundaryCell", boundaryCell)

    err = 0

    call mpas_tracer_advection_coefficients(&
         mesh, &
         config_horiz_tracer_adv_order, &
         derivTwo, &
         advCoefs, &
         advCoefs3rd, &
         nAdvCellsForEdge, &
         advCellsForEdge, &
         err, &
         maxLevelCell, &
         highOrderAdvectionMask, &
         boundaryCell)

    call define_tracer_connectivity()

  end subroutine init_advection

  !--------------------------------------------------------------------------
  
  subroutine define_tracer_connectivity()

    ! child tracer                                 ! parent tracer to child
    tracerConnectivity(1,1) = 'iceAreaCell'        ; tracerConnectivity(1,2) = 'none'
    tracerConnectivity(2,1) = 'surfaceTemperature' ; tracerConnectivity(2,2) = 'iceAreaCell'
    tracerConnectivity(3,1) = 'iceVolumeCell'      ; tracerConnectivity(3,2) = 'iceAreaCell'
    tracerConnectivity(4,1) = 'snowVolumeCell'     ; tracerConnectivity(4,2) = 'iceAreaCell'
    tracerConnectivity(5,1) = 'iceEnthalpy'        ; tracerConnectivity(5,2) = 'iceVolumeCell'
    tracerConnectivity(6,1) = 'iceSalinity'        ; tracerConnectivity(6,2) = 'iceVolumeCell'
    tracerConnectivity(7,1) = 'snowEnthalpy'       ; tracerConnectivity(7,2) = 'snowVolumeCell'

  end subroutine define_tracer_connectivity

  !--------------------------------------------------------------------------
  ! time stepping
  !--------------------------------------------------------------------------

  subroutine advection(block, dt)

    type (block_type), intent(inout) :: block

    real(kind=RKIND), intent(in) :: &
         dt

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracer_tendencies, &
         hexdyn, &
         normal

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon

    integer :: &
         iTracerVariable

    real(kind=RKIND), dimension(:), allocatable :: &
         edgeVelocity

    real(kind=RKIND), dimension(:,:), allocatable :: &    
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracer_tendencies", tracer_tendencies)
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "normal", normal)

    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)
    call MPAS_pool_get_array(normal, "normalVectorPolygon", normalVectorPolygon)

    call prepare_advection(&
         mesh, &
         tracers, &
         uVelocity, &
         vVelocity, &
         normalVectorPolygon, &
         edgeVelocity, &
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize)

    do iTracerVariable = 1, nTracerVariables

       call run_advection_variable(&
            mesh, &
            tracers, &
            tracer_tendencies, &
            tracerConnectivity(iTracerVariable,1), &
            tracerConnectivity(iTracerVariable,2), &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            dt)

    enddo ! iTracerVariable

    call finalize_advection(&
         mesh, &
         tracers, &
         edgeVelocity, &
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize)

  end subroutine advection

  !--------------------------------------------------------------------------
  
  subroutine advection_variable(&
       mesh, &
       tracers, &
       tracer_tendencies, &
       childTracerName, &
       parentTracerName, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       dt)

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh
  
    type (MPAS_pool_type), pointer :: &
         tracers, &
         tracer_tendencies
    
    character(len=*), intent(in) :: &
         childTracerName, &
         parentTracerName

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, &
         verticalCellSize
    
    real(kind=RKIND), intent(in) :: &
         dt

    type (mpas_pool_field_info_type) :: childFieldInfo

    call mpas_pool_get_field_info(tracers, trim(childTracerName), childFieldInfo)

    if (childFieldInfo % nDims == 3) then

       call advection_variable_3D(&
            mesh, &
            tracers, &
            tracer_tendencies, &
            childTracerName, &
            parentTracerName, &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            dt)
       
    else if (childFieldInfo % nDims == 4) then

       call advection_variable_4D(&
            mesh, &
            tracers, &
            tracer_tendencies, &
            childTracerName, &
            parentTracerName, &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            dt)

    endif

  end subroutine advection_variable

  !--------------------------------------------------------------------------

  subroutine advection_variable_3D(&
       mesh, &
       tracers, &
       tracer_tendencies, &
       childTracerName, &
       parentTracerName, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       dt)

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh
    
    type (MPAS_pool_type), pointer :: &
         tracers, &
         tracer_tendencies
    
    character(len=*), intent(in) :: &
         childTracerName, &
         parentTracerName

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, &
         verticalCellSize
    
    real(kind=RKIND), intent(in) :: &
         dt

    real (kind=RKIND), dimension(:,:,:), pointer :: &
         childTracerOld, &
         childTracerNew, &
         childTendency

    real (kind=RKIND), dimension(:,:,:), pointer :: &
         parentTracerOld, &
         parentTracerNew, &
         parentTendency    

    ! child tracer
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(childTracerName), childTendency)

    ! parent tracer
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerOld, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(parentTracerName), parentTendency)

    call advection_subvariable(&
         mesh, &
         childTracerOld(:,:,:), &
         childTracerNew(:,:,:), &
         childTendency(:,:,:), &
         parentTracerOld(:,:,:), &
         parentTracerNew(:,:,:), &
         parentTendency(:,:,:), &
         edgeVelocity, &
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize, &
         dt)

  end subroutine advection_variable_3D

  !--------------------------------------------------------------------------

  subroutine advection_variable_4D(&
       mesh, &
       tracers, &
       tracer_tendencies, &
       childTracerName, &
       parentTracerName, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       dt)
    
    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh

    type (MPAS_pool_type), pointer :: &
         tracers, &
         tracer_tendencies
    
    character(len=*), intent(in) :: &
         childTracerName, &
         parentTracerName

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, &
         verticalCellSize
    
    real(kind=RKIND), intent(in) :: &
         dt

    real (kind=RKIND), dimension(:,:,:,:), pointer :: &
         childTracerOld, &
         childTracerNew, &
         childTendency

    real (kind=RKIND), dimension(:,:,:), pointer :: &
         parentTracerOld, &
         parentTracerNew, &
         parentTendency  

    integer :: &
         iTracerDimension

    ! child tracer
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(childTracerName), childTendency)

    ! parent tracer
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(parentTracerName), parentTendency)
    
    do iTracerDimension = 1, size(childTracerNew,1)

       call advection_subvariable(&
            mesh, &
            childTracerOld(iTracerDimension,:,:,:), &
            childTracerNew(iTracerDimension,:,:,:), &
            childTendency(iTracerDimension,:,:,:), &
            parentTracerOld(:,:,:), &
            parentTracerNew(:,:,:), &
            parentTendency(:,:,:), &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            dt)
       
    enddo ! iTracerDimension
    
  end subroutine advection_variable_4D

  !--------------------------------------------------------------------------

  subroutine advection_subvariable(&
       mesh, &
       childTracerOld, &
       childTracerNew, &
       childTendency, &
       parentTracerOld, &
       parentTracerNew, &
       parentTendency, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       dt)

    use mpas_tracer_advection_mono, only: &
         mpas_tracer_advection_mono_tend

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
         childTracerOld, &
         childTracerNew, &
         childTendency

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         parentTracerOld, &
         parentTracerNew, &
         parentTendency

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, &
         verticalCellSize
    
    real(kind=RKIND), intent(in) :: &
         dt

    integer, pointer :: &
         nEdges, &
         nCells, &
         nCategories

    integer, dimension(:), pointer :: &
         nAdvCellsForEdge

    integer, dimension(:,:), pointer :: &
         cellsOnEdge, &
         advCellsForEdge, &
         highOrderAdvectionMask

    real(kind=RKIND), dimension(:,:), pointer :: &
         advCoefs, &
         advCoefs3rd

    integer :: &
         iCell, &
         iEdge, &
         iCategory, &
         iTracerDimension

    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)   
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)    
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "nAdvCellsForEdge", nAdvCellsForEdge)  
    call MPAS_pool_get_array(mesh, "advCellsForEdge", advCellsForEdge)  
    call MPAS_pool_get_array(mesh, "highOrderAdvectionMask", highOrderAdvectionMask)  
    call MPAS_pool_get_array(mesh, "advCoefs", advCoefs)  
    call MPAS_pool_get_array(mesh, "advCoefs3rd", advCoefs3rd)  

    ! edge flux of parent tracer
    do iEdge = 1, nEdges

       if (cellsOnEdge(1,iEdge) > nCells .or. cellsOnEdge(2,iEdge) > nCells) then

          ! edge of domain
          do iCategory = 1, nCategories
             edgeFlux(iCategory,iEdge) = 0.0_RKIND
          enddo ! iCategory

       else

          do iCategory = 1, nCategories
             edgeFlux(iCategory,iEdge) = edgeVelocity(iEdge) * &
                  0.5_RKIND * (parentTracerNew(1,iCategory,cellsOnEdge(1,iEdge)) + parentTracerNew(1,iCategory,cellsOnEdge(2,iEdge)))
          enddo ! iCategory

       endif

    enddo ! iEdge

    ! run the advection routine
    call mpas_tracer_advection_mono_tend(&
         childTracerOld,                  & ! (nTracers, nCategories, nCells)
         advCoefs,                        & ! (nAdvectionCells nEdges)
         advCoefs3rd,                     & ! (nAdvectionCells nEdges)
         nAdvCellsForEdge,                & ! (nEdges)
         advCellsForEdge,                 & ! (nAdvectionCells nEdges)
         edgeFlux,                        & ! (nCategories, nEdges)
         verticalVelocity,                & ! (nCategories+1, nEdges)
         parentTracerNew(1,:,:),          & ! (nCategories, nCells)
         verticalCellSize,                & ! (nCategories, nCells)
         dt,                              & ! 
         mesh,                            & ! 
         parentTendency(1,:,:),           & ! (nCategories, nCells)
         childTendency,                   & ! (nTracers, nCategories, nCells)
         highOrderAdvectionMask_in =      & ! (nCategories nEdges)
         highOrderAdvectionMask)!,          &
         !lvertical_in=.false.)

    ! apply the advective changes
    do iCell = 1, nCells
       do iCategory = 1, nCategories
          do iTracerDimension = 1, size(childTracerNew,1)

             childTracerNew(iTracerDimension,iCategory,iCell) = &
                  (childTracerOld(iTracerDimension,iCategory,iCell) * parentTracerOld(1,iCategory,iCell) &
                  + childTendency(iTracerDimension,iCategory,iCell) * dt) &
                  / parentTracerNew(1,iCategory,iCell)

          enddo ! iTracer
       enddo ! iCategory
    enddo ! iCell

  end subroutine advection_subvariable

  !--------------------------------------------------------------------------

  subroutine prepare_advection(&
       mesh, &
       tracers, &
       uVelocity, &
       vVelocity, &
       normalVectorPolygon, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize)
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorPolygon

    real(kind=RKIND), dimension(:), allocatable :: &
         edgeVelocity

    real(kind=RKIND), dimension(:,:), allocatable :: &    
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize

    integer, pointer :: &
         nEdges, &
         nCategories

    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    allocate(edgeVelocity(nEdges))
    allocate(edgeFlux(nCategories,nEdges))
    allocate(verticalVelocity(nCategories+1,nEdges))
    allocate(verticalCellSize(nCategories,nEdges))

    ! prepare vertical quantities
    verticalVelocity = 0.0_RKIND
    verticalCellSize = 1.0_RKIND
    
    call edge_from_vertex_velocity(&
         mesh, &
         uVelocity, &
         vVelocity, &
         normalVectorPolygon, &
         edgeVelocity)

    call prepare_tracers(&
         mesh, &
         tracers)

  end subroutine prepare_advection

  !--------------------------------------------------------------------------

  subroutine finalize_advection(&
       mesh, &
       tracers, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:), allocatable :: &
         edgeVelocity

    real(kind=RKIND), dimension(:,:), allocatable :: &    
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize

    call finalize_tracers(&
         mesh, &
         tracers)

    deallocate(edgeVelocity)
    deallocate(edgeFlux)
    deallocate(verticalVelocity)
    deallocate(verticalCellSize)

  end subroutine finalize_advection

  !--------------------------------------------------------------------------

  subroutine edge_from_vertex_velocity(&
       mesh, &
       uVelocity, &
       vVelocity, &
       normalVectorPolygon, &
       edgeVelocity)

    type(MPAS_pool_type), intent(in) :: mesh

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

    integer, pointer :: &
         nCells
    
    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell, &
         cellsOnEdge, &
         verticesOnEdge

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)

    ! loop over cells
    do iCell = 1, nCells

       ! loop over cell edges
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)

          iEdge = edgesOnCell(iEdgeOnCell,iCell)

          ! determine if velocity points outwards
          if (cellsOnEdge(1,iEdge) == iCell) then

             ! find u,v velocity at edge
             uVelocityEdge = 0.0_RKIND
             vVelocityEdge = 0.0_RKIND

             do iVertexOnEdge = 1, 2

                iVertex = verticesOnEdge(iVertexOnEdge, iEdge)

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

  subroutine prepare_tracers(&
       mesh, &
       tracers)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, pointer :: &
         nCells, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    do iCell = 1, nCells
       do iCategory = 1, nCategories

          if (iceAreaCategory(1,iCategory,iCell) > 0.0_RKIND) then

             iceVolumeCategory(1,iCategory,iCell) = iceVolumeCategory(1,iCategory,iCell) / &
                                                    iceAreaCategory(1,iCategory,iCell)

             snowVolumeCategory(1,iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell) / &
                                                     iceAreaCategory(1,iCategory,iCell)

          endif ! area > 0

       enddo ! iCategory
    enddo ! iCell

  end subroutine prepare_tracers

  !--------------------------------------------------------------------------

  subroutine finalize_tracers(&
       mesh, &
       tracers)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    type(MPAS_pool_type), pointer :: &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, pointer :: &
         nCells, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_shift_time_levels(tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    do iCell = 1, nCells
       do iCategory = 1, nCategories

          if (iceAreaCategory(1,iCategory,iCell) > 0.0_RKIND) then

             iceVolumeCategory(1,iCategory,iCell) = iceVolumeCategory(1,iCategory,iCell) * &
                                                    iceAreaCategory(1,iCategory,iCell)

             snowVolumeCategory(1,iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell) * &
                                                     iceAreaCategory(1,iCategory,iCell)

          endif ! area > 0

       enddo ! iCategory
    enddo ! iCell

  end subroutine finalize_tracers

  !--------------------------------------------------------------------------

end module cice_advection
