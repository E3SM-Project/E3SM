!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_advection
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_advection

  use mpas_grid_types

  implicit none

  private
  save 

  public :: &
       cice_init_advection, &
       cice_run_advection

  ! minimum sea ice area
  real(kind=RKIND), parameter, private :: &
       iceAreaMinimum = 1.0e-11_RKIND
  
  ! maximum number of tracers
  integer, parameter, private :: &
       nTracerVariables = 7

  ! single tracer connectivity. Defines the relationship between the child and its parent
  type, private :: tracerConnectivity
    
     character(len=200) :: &
          childTracerName  = "", & ! name of the child tracer
          parentTracerName = ""    ! name of the parent tracer to the child (or 'none' if none exists)

     integer :: &
          defined = 0 ! flag whether the tracer is defined for this entry
 
  end type tracerConnectivity

  ! array of all the tracer connectivities
  type(tracerConnectivity), dimension(nTracerVariables), private :: &
       tracerConnectivities

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------
  
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_advection
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_advection(block)!{{{

    use mpas_tracer_advection_helpers, only: &
         mpas_tracer_advection_coefficients
    
    use mpas_tracer_advection_mono, only: &
         mpas_tracer_advection_mono_init

    use cice_mesh, only: &
         cice_normal_vectors_polygon

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type (MPAS_pool_type), pointer :: &
         mesh, &
         configs, &
         velocity_solver

    integer, pointer :: &
         config_horiz_tracer_adv_order, &
         config_vert_tracer_adv_order

    real(kind=RKIND), pointer :: &
         config_coef_3rd_order

    logical, pointer :: &
         config_rotate_cartesian_grid

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         derivTwo, &
         normalVectorEdge

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
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)

    configs => block % configs
    call MPAS_pool_get_config(configs, "config_horiz_tracer_adv_order", config_horiz_tracer_adv_order)
    call MPAS_pool_get_config(configs, "config_vert_tracer_adv_order", config_vert_tracer_adv_order)
    call MPAS_pool_get_config(configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
    call MPAS_pool_get_config(configs, "config_coef_3rd_order", config_coef_3rd_order)

    call MPAS_pool_get_array(mesh, "derivTwo", derivTwo)
    call MPAS_pool_get_array(mesh, "advCoefs", advCoefs)
    call MPAS_pool_get_array(mesh, "advCoefs3rd", advCoefs3rd)
    call MPAS_pool_get_array(mesh, "nAdvCellsForEdge", nAdvCellsForEdge)
    call MPAS_pool_get_array(mesh, "advCellsForEdge", advCellsForEdge)
    call MPAS_pool_get_array(mesh, "maxLevelCell", maxLevelCell)
    call MPAS_pool_get_array(mesh, "highOrderAdvectionMask", highOrderAdvectionMask)
    call MPAS_pool_get_array(mesh, "boundaryCell", boundaryCell)

    call MPAS_pool_get_array(velocity_solver, "normalVectorEdge", normalVectorEdge)

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

    call mpas_tracer_advection_mono_init(&
         3, &
         config_horiz_tracer_adv_order, &
         config_vert_tracer_adv_order, &
         config_coef_3rd_order, &
         .false., .false., err)

    call define_tracer_connectivities()
    
    call cice_normal_vectors_polygon(&
         mesh, &
         normalVectorEdge, &
         config_rotate_cartesian_grid, &
         .false.)

  end subroutine cice_init_advection!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------
  
  subroutine define_tracer_connectivities()!{{{

    call init_tracer_connectivity(tracerConnectivities)

    ! This subroutine defines the tracer connectivities i.e. the child/parent relationship
    ! More senior parent tracers must be defined first
    ! The first tracer must have a parent of 'none'

    call add_tracer_connectivity(tracerConnectivities, "iceAreaCategory",    "none")
    call add_tracer_connectivity(tracerConnectivities, "surfaceTemperature", "iceAreaCategory")
    call add_tracer_connectivity(tracerConnectivities, "iceVolumeCategory",  "surfaceTemperature")
    call add_tracer_connectivity(tracerConnectivities, "snowVolumeCategory", "iceVolumeCategory")

    !call add_tracer_connectivity(tracerConnectivities, "iceAreaCategory",    "none")
    !call add_tracer_connectivity(tracerConnectivities, "surfaceTemperature", "iceAreaCategory")
    !call add_tracer_connectivity(tracerConnectivities, "iceVolumeCategory",  "iceAreaCategory")
    !call add_tracer_connectivity(tracerConnectivities, "snowVolumeCategory", "iceAreaCategory")
    !call add_tracer_connectivity(tracerConnectivities, "iceEnthalpy",        "iceVolumeCategory")
    !call add_tracer_connectivity(tracerConnectivities, "iceSalinity",        "iceVolumeCategory")
    !call add_tracer_connectivity(tracerConnectivities, "snowEnthalpy",       "iceVolumeCategory")

  end subroutine define_tracer_connectivities!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine add_tracer_connectivity(&
       tracerConnectivities, &
       childTracerName, &
       parentTracerName)!{{{

    type(tracerConnectivity), dimension(:), intent(inout) :: &
         tracerConnectivities !< Input/Output: 

    character(len=*), intent(in) :: &
         childTracerName, & !< Input: 
         parentTracerName   !< Input: 

    integer :: &
         iTracerVariable

    do iTracerVariable = 1, nTracerVariables

       if (tracerConnectivities(iTracerVariable) % defined == 0) then

          tracerConnectivities(iTracerVariable) % defined = 1
          tracerConnectivities(iTracerVariable) % childTracerName = trim(childTracerName)
          tracerConnectivities(iTracerVariable) % parentTracerName = trim(parentTracerName)

          exit

       endif

    enddo ! iTracerVariable

  end subroutine add_tracer_connectivity!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_tracer_connectivity(&
       tracerConnectivities)!{{{

    type(tracerConnectivity), dimension(:), intent(inout) :: &
         tracerConnectivities !< Input/Output: 

    integer :: iTracerVariable

    do iTracerVariable = 1, nTracerVariables    

            tracerConnectivities(iTracerVariable) % childTracerName = ""
            tracerConnectivities(iTracerVariable) % parentTracerName = ""
            tracerConnectivities(iTracerVariable) % defined = 0

    enddo ! iTracerVariable

  end subroutine init_tracer_connectivity!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine write_tracer_connectivities(unit)!{{{

    integer, intent(in) :: &
          unit !< Input: 

    integer :: &
         iTracerVariable

    write(unit,*) "-------------------------------------"

    do iTracerVariable = 1, nTracerVariables

       write(unit,*) iTracerVariable, &
            trim(tracerConnectivities(iTracerVariable) % childTracerName), &
            trim(tracerConnectivities(iTracerVariable) % parentTracerName), &
            tracerConnectivities(iTracerVariable) % defined

    enddo ! iTracerVariable

    write(unit,*) "-------------------------------------"    

  end subroutine write_tracer_connectivities!}}}

!-----------------------------------------------------------------------
! time stepping
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_advection
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_run_advection(block, dt)!{{{

    type (block_type), intent(inout) :: &
         block !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    type (MPAS_pool_type), pointer :: &
         configs, &
         mesh, &
         tracers, &
         tracer_tendencies, &
         tracer_edge_fluxes, &
         icestate, &
         velocity_solver

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorEdge

    integer :: &
         iTracerVariable

    real(kind=RKIND), dimension(:), allocatable :: &
         edgeVelocity

    real(kind=RKIND), dimension(:,:), allocatable :: &    
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize

    configs => block % configs

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracer_tendencies", tracer_tendencies)
    call MPAS_pool_get_subpool(block % structs, "tracer_edge_fluxes", tracer_edge_fluxes)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
    call MPAS_pool_get_array(velocity_solver, "normalVectorEdge", normalVectorEdge)

    call prepare_advection(&
         configs, &
         mesh, &
         tracers, &
         icestate, &
         uVelocity, &
         vVelocity, &
         normalVectorEdge, &
         edgeVelocity, &
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize)

    do iTracerVariable = 1, nTracerVariables

       if (tracerConnectivities(iTracerVariable) % defined == 1) then
       
          !write(*,*) iTracerVariable, trim(tracerConnectivities(iTracerVariable) % childTracerName)

          call run_advection_variable(&
               mesh, &
               tracers, &
               tracer_tendencies, &
               tracer_edge_fluxes, &
               icestate, &
               trim(tracerConnectivities(iTracerVariable) % childTracerName), &
               trim(tracerConnectivities(iTracerVariable) % parentTracerName), &
               edgeVelocity, &
               edgeFlux, &
               verticalVelocity, &
               verticalCellSize, &
               dt)

       endif

    enddo ! iTracerVariable

    call finalize_advection(&
         configs, &
         mesh, &
         tracers, &
         icestate, &
         edgeVelocity, &
         edgeFlux, &
         verticalVelocity, &
         verticalCellSize)

  end subroutine cice_run_advection!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  run_advection_variable
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------
  
  subroutine run_advection_variable(&
       mesh, &
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       icestate, &
       childTracerName, &
       parentTracerName, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       dt)!{{{

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 
  
    type (MPAS_pool_type), pointer :: &
         tracers, &            !< Input/Output: 
         tracer_tendencies, &  !< Input/Output: 
         tracer_edge_fluxes, & !< Input/Output: 
         icestate              !< Input/Output: 
    
    character(len=*), intent(in) :: &
         childTracerName, & !< Input: 
         parentTracerName   !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input: 
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, & !< Input: 
         verticalCellSize    !< Input: 
    
    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    type (mpas_pool_field_info_type) :: childFieldInfo

    call mpas_pool_get_field_info(tracers, trim(childTracerName), childFieldInfo)

    if (childFieldInfo % nDims == 3) then

       call run_advection_variable_3D(&
            mesh, &
            tracers, &
            tracer_tendencies, &
            tracer_edge_fluxes, &
            icestate, &
            childTracerName, &
            parentTracerName, &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            dt)
       
    else if (childFieldInfo % nDims == 4) then

       call run_advection_variable_4D(&
            mesh, &
            tracers, &
            tracer_tendencies, &
            tracer_edge_fluxes, &
            icestate, &
            childTracerName, &
            parentTracerName, &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            dt)

    endif

  end subroutine run_advection_variable!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  run_advection_variable_3D
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine run_advection_variable_3D(&
       mesh, &
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       icestate, &
       childTracerName, &
       parentTracerName, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       dt)!{{{

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 
    
    type (MPAS_pool_type), pointer :: &
         tracers, &            !< Input/Output: 
         tracer_tendencies, &  !< Input/Output: 
         tracer_edge_fluxes, & !< Input/Output: 
         icestate              !< Input/Output: 
    
    character(len=*), intent(in) :: &
         childTracerName, & !< Input: 
         parentTracerName   !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input: 
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, & !< Input: 
         verticalCellSize    !< Input: 
    
    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    real (kind=RKIND), dimension(:,:,:), pointer :: &
         childTracerOld, &
         childTracerNew, &
         childTendency, &
         childEdgeFlux

    real (kind=RKIND), dimension(:,:,:), pointer :: &
         parentTracerOld, &
         parentTracerNew, &
         parentTendency, &
         parentEdgeFlux

    real (kind=RKIND), dimension(:,:,:), allocatable :: &
         parentTracerOldNone, &
         parentTracerNewNone, &
         parentTendencyNone, &
         parentEdgeFluxNone

    integer, pointer :: &
         ONE, &
         nCategories, &
         nCells, &
         nEdges

    integer :: &
         iCategory, &
         iCell, &
         iEdgeOnCell, &
         iEdge

    if (trim(parentTracerName) == 'none') then

       ! child tracer
       call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerOld, 1)
       call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 2)
       call MPAS_pool_get_array(tracer_tendencies, trim(childTracerName), childTendency)
       call MPAS_pool_get_array(tracer_edge_fluxes, trim(childTracerName), childEdgeFlux)

       ! parent tracer
       call MPAS_pool_get_dimension(mesh, "ONE", ONE)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nCells", nCells)
       call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)

       allocate(parentTracerOldNone(ONE,nCategories,nCells))
       allocate(parentTracerNewNone(ONE,nCategories,nCells))
       allocate(parentTendencyNone(ONE,nCategories,nCells))
       allocate(parentEdgeFluxNone(ONE,nCategories,nEdges))

       call prepare_none_parent_tracer(&
            mesh, &
            icestate, &
            childTracerOld, &
            edgeVelocity, &
            parentTracerOldNone, &
            parentTracerNewNone, &
            parentTendencyNone, &
            parentEdgeFluxNone)

       call run_advection_subvariable(&
            mesh, &
            icestate, &
            childTracerOld(:,:,:), &
            childTracerNew(:,:,:), &
            childTendency(:,:,:), &
            childEdgeFlux(:,:,:), &
            parentTracerOldNone(:,:,:), &
            parentTracerNewNone(:,:,:), &
            parentTendencyNone(:,:,:), &
            parentEdgeFluxNone(:,:,:), &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            .true., &
            dt)

       deallocate(parentTracerOldNone)
       deallocate(parentTracerNewNone)
       deallocate(parentTendencyNone)
       deallocate(parentEdgeFluxNone)

    else

       ! child tracer
       call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerOld, 1)
       call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 2)
       call MPAS_pool_get_array(tracer_tendencies, trim(childTracerName), childTendency)
       call MPAS_pool_get_array(tracer_edge_fluxes, trim(childTracerName), childEdgeFlux)

       ! parent tracer
       call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerOld, 1)
       call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 2)
       call MPAS_pool_get_array(tracer_tendencies, trim(parentTracerName), parentTendency)
       call MPAS_pool_get_array(tracer_edge_fluxes, trim(parentTracerName), parentEdgeFlux)

       call run_advection_subvariable(&
            mesh, &
            icestate, &
            childTracerOld(:,:,:), &
            childTracerNew(:,:,:), &
            childTendency(:,:,:), &
            childEdgeFlux(:,:,:), &
            parentTracerOld(:,:,:), &
            parentTracerNew(:,:,:), &
            parentTendency(:,:,:), &
            parentEdgeFlux(:,:,:), &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            .false., &
            dt)

    endif

  end subroutine run_advection_variable_3D!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_none_parent_tracer
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine prepare_none_parent_tracer(&
       mesh, &
       icestate, &
       childTracerOld, &
       edgeVelocity, &
       parentTracerOldNone, &
       parentTracerNewNone, &
       parentTendencyNone, &
       parentEdgeFluxNone)!{{{

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh, &  !< Input: 
         icestate !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         childTracerOld !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &    
         edgeVelocity !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
         parentTracerOldNone, & !< Input/Output: 
         parentTracerNewNone, & !< Input/Output: 
         parentTendencyNone, &  !< Input/Output: 
         parentEdgeFluxNone     !< Input/Output: 

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnCell, &
         cellsOnEdge, &
         iceAreaCategoryMask

    integer, pointer :: &
         nCells, &
         nEdges, &
         nCategories

    integer :: &
         iCell, &
         iCategory, &
         iEdgeOnCell, &
         iEdge

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnCell", cellsOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)

    call MPAS_pool_get_array(icestate, "iceAreaCategoryMask", iceAreaCategoryMask)

    do iCell = 1, nCells
       do iCategory = 1, nCategories

          iceAreaCategoryMask(iCategory,iCell) = 0
          parentTracerOldNone(1,iCategory,iCell) = 0.0_RKIND ! need this?
          parentTracerNewNone(1,iCategory,iCell) = 0.0_RKIND ! need this?

          if (childTracerOld(1,iCategory,iCell) > iceAreaMinimum) then
             iceAreaCategoryMask(iCategory,iCell) = 1
             parentTracerOldNone(1,iCategory,iCell) = 1.0_RKIND
             parentTracerNewNone(1,iCategory,iCell) = 1.0_RKIND
             parentTendencyNone(1,iCategory,iCell) = 0.0_RKIND 
          endif
          
          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             if (childTracerOld(1,iCategory,cellsOnCell(iEdgeOnCell,iCell)) > iceAreaMinimum) then
                iceAreaCategoryMask(iCategory,iCell) = 1
                parentTracerOldNone(1,iCategory,iCell) = 1.0_RKIND
                parentTracerNewNone(1,iCategory,iCell) = 1.0_RKIND
                parentTendencyNone(1,iCategory,iCell) = 0.0_RKIND 
                exit
             endif

          enddo ! iEdgeOnCell

       enddo ! iCategory
    enddo ! iCell

    do iEdge = 1, nEdges
       do iCategory = 1, nCategories
          if (iceAreaCategoryMask(iCategory,cellsOnEdge(1,iEdge)) == 1 .or. &
              iceAreaCategoryMask(iCategory,cellsOnEdge(2,iEdge)) == 1) then
             parentEdgeFluxNone(1,iCategory,iEdge) = edgeVelocity(iEdge)
          endif
       enddo ! iCategory
    enddo ! iEdge

  end subroutine prepare_none_parent_tracer!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  run_advection_variable_4D
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine run_advection_variable_4D(&
       mesh, &
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       icestate, &
       childTracerName, &
       parentTracerName, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       dt)!{{{
    
    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type (MPAS_pool_type), pointer :: &
         tracers, &            !< Input/Output: 
         tracer_tendencies, &  !< Input/Output: 
         tracer_edge_fluxes, & !< Input/Output: 
         icestate              !< Input/Output: 
    
    character(len=*), intent(in) :: &
         childTracerName, & !< Input: 
         parentTracerName   !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input: 
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, & !< Input: 
         verticalCellSize    !< Input: 
    
    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    real (kind=RKIND), dimension(:,:,:,:), pointer :: &
         childTracerOld, &
         childTracerNew, &
         childTendency, &
         childEdgeFlux

    real (kind=RKIND), dimension(:,:,:), pointer :: &
         parentTracerOld, &
         parentTracerNew, &
         parentTendency, &
         parentEdgeFlux

    integer :: &
         iTracerDimension

    ! child tracer
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(childTracerName), childTendency)
    call MPAS_pool_get_array(tracer_edge_fluxes, trim(childTracerName), childEdgeFlux)

    ! parent tracer
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(parentTracerName), parentTendency)
    call MPAS_pool_get_array(tracer_edge_fluxes, trim(parentTracerName), parentEdgeFlux)
    
    do iTracerDimension = 1, size(childTracerNew,1)

       call run_advection_subvariable(&
            mesh, &
            icestate, &
            childTracerOld(iTracerDimension,:,:,:), & ! Inefficient!
            childTracerNew(iTracerDimension,:,:,:), & ! Inefficient!
            childTendency(iTracerDimension,:,:,:), &  ! Inefficient!
            childEdgeFlux(iTracerDimension,:,:,:), &  ! Inefficient!
            parentTracerOld(:,:,:), &
            parentTracerNew(:,:,:), &
            parentTendency(:,:,:), &
            parentEdgeFlux(:,:,:), &
            edgeVelocity, &
            edgeFlux, &
            verticalVelocity, &
            verticalCellSize, &
            .false., &
            dt)
       
    enddo ! iTracerDimension
    
  end subroutine run_advection_variable_4D!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  run_advection_subvariable
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine run_advection_subvariable(&
       mesh, &
       icestate, &
       childTracerOld, &
       childTracerNew, &
       childTendency, &
       childEdgeFlux, &
       parentTracerOld, &
       parentTracerNew, &
       parentTendency, &
       parentEdgeFlux, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize, &
       isIceArea, &
       dt)!{{{

    use mpas_tracer_advection_mono, only: &
         mpas_tracer_advection_mono_tend_seaice

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type (MPAS_pool_type), pointer :: &
         icestate !< Input/Output: 

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
         childTracerOld, & !< Input/Output: 
         childTracerNew, & !< Input/Output: 
         childTendency, &  !< Input/Output: 
         childEdgeFlux     !< Input/Output: 

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
         parentTracerOld, & !< Input/Output: 
         parentTracerNew, & !< Input/Output: 
         parentTendency, &  !< Input/Output: 
         parentEdgeFlux     !< Input/Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input: 
         
    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         verticalVelocity, & !< Input: 
         verticalCellSize    !< Input: 

    logical, intent(in) :: &
         isIceArea !< Input: 
    
    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

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

    integer, dimension(:,:), pointer :: &
         iceAreaCategoryMask

    integer :: &
         iCell, &
         iEdge, &
         iCategory, &
         iTracerDimension

    integer, parameter :: &
         iCellPoint = 81*40

    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)   
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)    
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "nAdvCellsForEdge", nAdvCellsForEdge)  
    call MPAS_pool_get_array(mesh, "advCellsForEdge", advCellsForEdge)  
    call MPAS_pool_get_array(mesh, "highOrderAdvectionMask", highOrderAdvectionMask)  
    call MPAS_pool_get_array(mesh, "advCoefs", advCoefs)  
    call MPAS_pool_get_array(mesh, "advCoefs3rd", advCoefs3rd)  
    
    call MPAS_pool_get_array(icestate, "iceAreaCategoryMask", iceAreaCategoryMask)
    
    ! zero out child tendancy
    childTendency(:,:,:) = 0.0_RKIND

    ! run the advection routine
    call mpas_tracer_advection_mono_tend_seaice(&
         childTracerOld,                  & ! (nTracers, nCategories, nCells)
         advCoefs,                        & ! (nAdvectionCells nEdges)
         advCoefs3rd,                     & ! (nAdvectionCells nEdges)
         nAdvCellsForEdge,                & ! (nEdges)
         advCellsForEdge,                 & ! (nAdvectionCells nEdges)
         parentEdgeFlux(1,:,:),           & ! (nCategories, nEdges)
         verticalVelocity,                & ! (nCategories+1, nEdges)
         parentTracerOld(1,:,:),          & ! (nCategories, nCells)
         verticalCellSize,                & ! (nCategories, nCells)
         dt,                              & ! 
         mesh,                            & ! 
         parentTendency(1,:,:),           & ! (nCategories, nCells)
         childTendency                    & ! (nTracers, nCategories, nCells)
         !highOrderAdvectionMask_in =      & ! (nCategories nEdges)
         !highOrderAdvectionMask,          &
         )

    ! apply the advective changes
    do iCell = 1, nCells
       do iCategory = 1, nCategories

          if (iceAreaCategoryMask(iCategory,iCell) == 1) then
             
             do iTracerDimension = 1, size(childTracerNew,1)

                ! calculate the new tracer
                childTracerNew(iTracerDimension,iCategory,iCell) = &
                     childTracerOld(iTracerDimension,iCategory,iCell) * parentTracerOld(iTracerDimension,iCategory,iCell) + childTendency(iTracerDimension,iCategory,iCell) * dt
                
                ! store the old child tracer multiplied by the old parent. This will be the parent tracer of the next child
                childTracerOld(iTracerDimension,iCategory,iCell) = childTracerOld(iTracerDimension,iCategory,iCell) * parentTracerOld(iTracerDimension,iCategory,iCell)

             enddo ! iTracer

          endif

       enddo ! iCategory
    enddo ! iCell

    ! if th etracer is ice area update the ice area mask
    if (isIceArea) then
       do iCell = 1, nCells
          do iCategory = 1, nCategories

             if (childTracerNew(1,iCategory,iCell) > iceAreaMinimum) then
                iceAreaCategoryMask(iCategory,iCell) = 1
             else
                iceAreaCategoryMask(iCategory,iCell) = 0
             endif

          enddo ! iCell
       enddo ! iCell
    endif

  end subroutine run_advection_subvariable!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_advection
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine prepare_advection(&
       configs, &
       mesh, &
       tracers, &
       icestate, &
       uVelocity, &
       vVelocity, &
       normalVectorEdge, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize)!{{{
    
    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, &  !< Input: 
         mesh        !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers, & !< Input/Output: 
         icestate   !< Input/Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input: 
         vVelocity    !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorEdge !< Input: 

    real(kind=RKIND), dimension(:), allocatable, intent(inout) :: &
         edgeVelocity !< Input/Output: 

    real(kind=RKIND), dimension(:,:), allocatable, intent(inout) :: &    
         edgeFlux, &         !< Input/Output: 
         verticalVelocity, & !< Input/Output: 
         verticalCellSize    !< Input/Output: 

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
         normalVectorEdge, &
         edgeVelocity)

    call prepare_tracers(&
         configs, &
         mesh, &
         tracers, &
         icestate)

  end subroutine prepare_advection!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  finalize_advection
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine finalize_advection(&
       configs, &
       mesh, &
       tracers, &
       icestate, &
       edgeVelocity, &
       edgeFlux, &
       verticalVelocity, &
       verticalCellSize)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers, & !< Input/Output: 
         icestate   !< Input/Output: 

    real(kind=RKIND), dimension(:), allocatable, intent(inout) :: &
         edgeVelocity !< Input/Output: 

    real(kind=RKIND), dimension(:,:), allocatable, intent(inout) :: &    
         edgeFlux, &         !< Input/Output: 
         verticalVelocity, & !< Input/Output: 
         verticalCellSize    !< Input/Output: 

    call finalize_tracers(&
         configs, &
         mesh, &
         tracers, &
         icestate)

    deallocate(edgeVelocity)
    deallocate(edgeFlux)
    deallocate(verticalVelocity)
    deallocate(verticalCellSize)

  end subroutine finalize_advection!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  edge_from_vertex_velocity
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine edge_from_vertex_velocity(&
       mesh, &
       uVelocity, &
       vVelocity, &
       normalVectorEdge, &
       edgeVelocity)!{{{

    type(MPAS_pool_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input: 
         vVelocity    !< Input: 

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorEdge !< Input: 

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
                  uVelocityEdge * normalVectorEdge(1,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorEdge(2,iEdgeOnCell,iCell)

          endif

       enddo ! iEdgeOnCell

    enddo ! iCell

  end subroutine edge_from_vertex_velocity!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  prepare_tracers
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine prepare_tracers(&
       configs, &
       mesh, &
       tracers, &
       icestate)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, & !< Input: 
         mesh       !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers, & !< Input/Output: 
         icestate   !< Input/Output: 

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    integer, dimension(:,:), pointer :: &
         iceAreaCategoryMask

    integer, pointer :: &
         nCells, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    logical, pointer :: &
         config_convert_volume_to_thickness

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(icestate, "iceAreaCategoryMask", iceAreaCategoryMask)

    call MPAS_pool_get_config(configs, "config_convert_volume_to_thickness", config_convert_volume_to_thickness)
    
    if (config_convert_volume_to_thickness) then

       do iCell = 1, nCells
          do iCategory = 1, nCategories

             if (iceAreaCategoryMask(iCategory,iCell) == 1) then

                iceVolumeCategory(1,iCategory,iCell) = iceVolumeCategory(1,iCategory,iCell) / &
                     iceAreaCategory(1,iCategory,iCell)

                snowVolumeCategory(1,iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell) / &
                     iceAreaCategory(1,iCategory,iCell)

             endif ! area > 0

          enddo ! iCategory
       enddo ! iCell

    endif ! config_convert_volume_to_thickness

  end subroutine prepare_tracers!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  finalize_tracers
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine finalize_tracers(&
       configs, &
       mesh, &
       tracers, &
       icestate)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, & !< Input: 
         mesh       !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers, & !< Input/Output: 
         icestate   !< Input/Output: 

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         surfaceTemperature, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, dimension(:,:), pointer :: &
         iceAreaCategoryMask

    integer, pointer :: &
         nCells, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    logical, pointer :: &
         config_convert_volume_to_thickness

    call MPAS_pool_shift_time_levels(tracers)

    call scale_tracers_back(mesh, tracers, icestate)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(icestate, "iceAreaCategoryMask", iceAreaCategoryMask)

    call MPAS_pool_get_config(configs, "config_convert_volume_to_thickness", config_convert_volume_to_thickness)

    if (config_convert_volume_to_thickness) then

       do iCell = 1, nCells
          do iCategory = 1, nCategories

             if (iceAreaCategoryMask(iCategory,iCell) == 1) then

                iceVolumeCategory(1,iCategory,iCell) = iceVolumeCategory(1,iCategory,iCell) * &
                     iceAreaCategory(1,iCategory,iCell)

                snowVolumeCategory(1,iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell) * &
                     iceAreaCategory(1,iCategory,iCell)

             endif ! area > 0

          enddo ! iCategory
       enddo ! iCell

    endif ! config_convert_volume_to_thickness

  end subroutine finalize_tracers!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  scale_tracers_back
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine scale_tracers_back(&
       mesh, &
       tracers, &
       icestate)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers, & !< Input/Output: 
         icestate   !< Input/Output: 

    integer :: &
         iTracerVariable

    type (mpas_pool_field_info_type) :: childFieldInfo

    do iTracerVariable = nTracerVariables, 1, -1

       if (tracerConnectivities(iTracerVariable) % defined == 1 .and. &
           trim(tracerConnectivities(iTracerVariable) % parentTracerName) /= "none") then

          call mpas_pool_get_field_info(tracers, trim(tracerConnectivities(iTracerVariable) % childTracerName), childFieldInfo)
          
          if (childFieldInfo % nDims == 3) then
             
             call scale_tracers_back_3D(&
                  mesh, &
                  tracers, &
                  icestate, &
                  trim(tracerConnectivities(iTracerVariable) % childTracerName), &
                  trim(tracerConnectivities(iTracerVariable) % parentTracerName))
             
          else if (childFieldInfo % nDims == 4) then
             
             call scale_tracers_back_4D(&
                  mesh, &
                  tracers, &
                  icestate, &
                  trim(tracerConnectivities(iTracerVariable) % childTracerName), &
                  trim(tracerConnectivities(iTracerVariable) % parentTracerName))
             
          endif
          
       endif

    enddo ! iTracer

  end subroutine scale_tracers_back!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  scale_tracers_back_3D
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine scale_tracers_back_3D(&
       mesh, &
       tracers, &
       icestate, &
       childTracerName, &
       parentTracerName)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers, & !< Input/Output: 
         icestate   !< Input/Output: 

    character(len=*), intent(in) :: &
         childTracerName, & !< Input: 
         parentTracerName   !< Input: 

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         childTracerNew, &
         parentTracerNew

    integer, dimension(:,:), pointer :: &
         iceAreaCategoryMask

    integer, pointer :: &
         nCells, &
         nCategories

    integer :: &
         iCell, &
         iCategory, &
         iTracer

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
         
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 1)
    call MPAS_pool_get_array(icestate, "iceAreaCategoryMask", iceAreaCategoryMask)

    do iCell = 1, nCells
       do iCategory = 1, nCategories
          do iTracer = 1, size(childTracerNew,1)

             if (iceAreaCategoryMask(iCategory,iCell) == 1) then

                childTracerNew(iTracer,iCategory,iCell) = &
                     childTracerNew(iTracer,iCategory,iCell) / parentTracerNew(1,iCategory,iCell)

             else

                childTracerNew(iTracer,iCategory,iCell) = 0.0_RKIND

             endif

          enddo ! iTracer
       enddo ! iCategory
    enddo ! iCell

  end subroutine scale_tracers_back_3D!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  scale_tracers_back_4D
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine scale_tracers_back_4D(&
       mesh, &
       tracers, &
       icestate, &
       childTracerName, &
       parentTracerName)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    type(MPAS_pool_type), pointer :: &
         tracers, & !< Input/Output: 
         icestate   !< Input/Output: 

    character(len=*), intent(in) :: &
         childTracerName, & !< Input: 
         parentTracerName   !< Input: 

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         parentTracerNew

    real(kind=RKIND), dimension(:,:,:,:), pointer :: &
         childTracerNew

    integer, dimension(:,:), pointer :: &
         iceAreaCategoryMask

    integer, pointer :: &
         nCells, &
         nCategories

    integer :: &
         iCell, &
         iCategory, &
         iTracer, &
         iTracerDimension

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
         
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 1)
    call MPAS_pool_get_array(icestate, "iceAreaCategoryMask", iceAreaCategoryMask)


    do iCell = 1, nCells
       do iCategory = 1, nCategories
          do iTracer = 1, size(childTracerNew,2)
             do iTracerDimension = 1, size(childTracerNew,1)

                if (iceAreaCategoryMask(iCategory,iCell) == 1) then

                   childTracerNew(iTracerDimension,iTracer,iCategory,iCell) = &
                        childTracerNew(iTracerDimension,iTracer,iCategory,iCell) / parentTracerNew(1,iCategory,iCell)

                else

                   childTracerNew(iTracerDimension,iTracer,iCategory,iCell) = 0.0_RKIND                   

                endif

             enddo ! iTracerDimension
          enddo ! iTracer
       enddo ! iCategory
    enddo ! iCell

  end subroutine scale_tracers_back_4D!}}}

!-----------------------------------------------------------------------

end module cice_advection
