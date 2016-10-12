!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_advection_upwind
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_advection_upwind

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_dmpar

  implicit none

  private
  save

  public :: &
       cice_init_advection_upwind, &
       cice_run_advection_upwind

  ! maximum number of tracers
  integer, parameter, private :: &
       nTracerVariables = 7

  ! single tracer connectivity. Defines the relationship between the child and its parent
  type, private :: tracerConnectivity

     character(len=200) :: &
          childTracerName  = "", & ! name of the child tracer
          parentTracerName = ""    ! name of the parent tracer to the child (or 'none' if none exists)

     logical :: &
          updateHalo ! logical flag whether the tracer undergoes a halo update before being advected

     integer :: &
          defined = 0 ! flag whether the tracer is defined for this entry

     real(kind=RKIND) :: &
          childTracerMinimum, &
          parentTracerMinimum

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
!  cice_init_advection_upwind
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_advection_upwind(domain)!{{{

    !use mpas_tracer_advection_helpers, only: &
    !     mpas_tracer_advection_coefficients

    !use mpas_tracer_advection_mono, only: &
    !     mpas_tracer_advection_mono_init

    use cice_mesh, only: &
         cice_normal_vectors_polygon

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

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
         config_rotate_cartesian_grid, &
         config_use_velocity_solver

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         derivTwo

    real(kind=RKIND), dimension(:,:,:), pointer :: &
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

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)

       configs => block % configs
       call MPAS_pool_get_config(configs, "config_use_velocity_solver", config_use_velocity_solver)
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

       !call mpas_tracer_advection_coefficients(&
       !     mesh, &
       !     config_horiz_tracer_adv_order, &
       !     derivTwo, &
       !     advCoefs, &
       !     advCoefs3rd, &
       !     nAdvCellsForEdge, &
       !     advCellsForEdge, &
       !     err, &
       !     maxLevelCell, &
       !     highOrderAdvectionMask, &
       !     boundaryCell)

       !call mpas_tracer_advection_mono_init(&
       !     3, &
       !     config_horiz_tracer_adv_order, &
       !     config_vert_tracer_adv_order, &
       !     config_coef_3rd_order, &
       !     .false., .true., err)

       call define_tracer_connectivities(&
            config_use_velocity_solver)

       call cice_normal_vectors_polygon(&
            mesh, &
            normalVectorEdge, &
            config_rotate_cartesian_grid, &
            .false.)

       block => block % next
    enddo

  end subroutine cice_init_advection_upwind!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  define_tracer_connectivities
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine define_tracer_connectivities(config_use_velocity_solver)!{{{

    use cice_constants, only: &
         iceAreaMinimum, &
         iceThicknessMinimum, &
         snowThicknessMinimum

    logical, intent(in) :: &
         config_use_velocity_solver

    call init_tracer_connectivity(tracerConnectivities)

    ! This subroutine defines the tracer connectivities i.e. the child/parent relationship
    ! More senior parent tracers must be defined first
    ! The first tracer must have a parent of 'none'

    call add_tracer_connectivity(tracerConnectivities, "iceAreaCategory",    "none")
    call add_tracer_connectivity(tracerConnectivities, "surfaceTemperature", "iceAreaCategory")
    call add_tracer_connectivity(tracerConnectivities, "iceVolumeCategory",  "surfaceTemperature")
    call add_tracer_connectivity(tracerConnectivities, "snowVolumeCategory", "iceVolumeCategory")

    !if (config_use_velocity_solver) then
    !   call add_tracer_connectivity(tracerConnectivities, "iceAreaCategory",    "none",            updateHalo = .false., childTracerMinimum=iceAreaMinimum)
    !   call add_tracer_connectivity(tracerConnectivities, "iceVolumeCategory",  "iceAreaCategory", updateHalo = .false., childTracerMinimum=iceThicknessMinimum)
    !   call add_tracer_connectivity(tracerConnectivities, "snowVolumeCategory", "iceAreaCategory", updateHalo = .false., childTracerMinimum=snowThicknessMinimum)
    !else
    !   call add_tracer_connectivity(tracerConnectivities, "iceAreaCategory",    "none",            childTracerMinimum=iceAreaMinimum)
    !   call add_tracer_connectivity(tracerConnectivities, "iceVolumeCategory",  "iceAreaCategory", childTracerMinimum=iceThicknessMinimum)
    !   call add_tracer_connectivity(tracerConnectivities, "snowVolumeCategory", "iceAreaCategory", childTracerMinimum=snowThicknessMinimum)
    !endif

    !call add_tracer_connectivity(tracerConnectivities, "surfaceTemperature", "iceAreaCategory")
    !call add_tracer_connectivity(tracerConnectivities, "iceEnthalpy",        "iceVolumeCategory")
    !call add_tracer_connectivity(tracerConnectivities, "iceSalinity",        "iceVolumeCategory")
    !call add_tracer_connectivity(tracerConnectivities, "snowEnthalpy",       "iceVolumeCategory")

    call add_parent_tracer_minimums(tracerConnectivities)

  end subroutine define_tracer_connectivities!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  add_tracer_connectivity
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
       parentTracerName, &
       updateHalo, &
       childTracerMinimum)!{{{

    type(tracerConnectivity), dimension(:), intent(inout) :: &
         tracerConnectivities !< Input/Output:

    character(len=*), intent(in) :: &
         childTracerName, & !< Input:
         parentTracerName   !< Input:

    logical, optional, intent(in) :: &
         updateHalo !< Input:

    real(kind=RKIND), optional, intent(in) :: &
         childTracerMinimum !< Input:

    integer :: &
         iTracerVariable

    do iTracerVariable = 1, nTracerVariables

       if (tracerConnectivities(iTracerVariable) % defined == 0) then

          tracerConnectivities(iTracerVariable) % defined = 1
          tracerConnectivities(iTracerVariable) % childTracerName = trim(childTracerName)
          tracerConnectivities(iTracerVariable) % parentTracerName = trim(parentTracerName)

          if (present(updateHalo)) then
             tracerConnectivities(iTracerVariable) % updateHalo = updateHalo
          else
             tracerConnectivities(iTracerVariable) % updateHalo = .true.
          endif

          if (present(childTracerMinimum)) then
             tracerConnectivities(iTracerVariable) % childTracerMinimum = childTracerMinimum
          else
             tracerConnectivities(iTracerVariable) % childTracerMinimum = 0.0_RKIND
          endif

          exit

       endif

    enddo ! iTracerVariable

  end subroutine add_tracer_connectivity!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_tracer_connectivity
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
       tracerConnectivities(iTracerVariable) % childTracerMinimum = 0.0_RKIND
       tracerConnectivities(iTracerVariable) % parentTracerMinimum = 0.0_RKIND

    enddo ! iTracerVariable

  end subroutine init_tracer_connectivity!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  write_tracer_connectivities
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

    write(unit,*) "--------------------------------------------------------------------"//&
                  "--------------------------------------------------------------------"

    write(unit,fmt='(a3,a3,a40,a3,a40,a3,a7,a3,a15,a3,a15,a3)') "#", " | ", "Child name", " | ", "Parent name", " | ", &
         "Defined", " | ", "Child minimum", " | ", "Parent Minimum" , " | "
    write(unit,fmt='(a3,a3,a40,a3,a40,a3,a7,a3,a15,a3,a15,a3)') "-", " | ", "----------", " | ", "-----------", " | ", &
         "-------", " | ", "-------------", " | ", "--------------" , " | "

    do iTracerVariable = 1, nTracerVariables

       write(unit,fmt='(i3,a3,a40,a3,a40,a3,i7,a3,e15.2,a3,e15.2,a3)') iTracerVariable, " | ", &
            trim(tracerConnectivities(iTracerVariable) % childTracerName), " | ", &
            trim(tracerConnectivities(iTracerVariable) % parentTracerName), " | ", &
            tracerConnectivities(iTracerVariable) % defined, " | ", &
            tracerConnectivities(iTracerVariable) % childTracerMinimum, " | ", &
            tracerConnectivities(iTracerVariable) % parentTracerMinimum, " | "

    enddo ! iTracerVariable

    write(unit,*) "--------------------------------------------------------------------"//&
                  "--------------------------------------------------------------------"

  end subroutine write_tracer_connectivities!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  add_parent_tracer_minimums
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date September 16th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine add_parent_tracer_minimums(tracerConnectivities)!{{{

    type(tracerConnectivity), dimension(:), intent(inout) :: &
         tracerConnectivities !< Input/Output:

    integer :: &
         iTracerVariable, &
         iTracerVariable2

    do iTracerVariable = 1, nTracerVariables

       do iTracerVariable2 = 1, nTracerVariables

          if (trim(tracerConnectivities(iTracerVariable) % parentTracerName) == &
              trim(tracerConnectivities(iTracerVariable2) % childTracerName)) &
               tracerConnectivities(iTracerVariable) % parentTracerMinimum = &
               tracerConnectivities(iTracerVariable2) % childTracerMinimum

          if (trim(tracerConnectivities(iTracerVariable) % parentTracerName) == "none") &
               tracerConnectivities(iTracerVariable) % parentTracerMinimum = 0.0_RKIND

       enddo ! iTracerVariable2

    enddo ! iTracerVariable

  end subroutine add_parent_tracer_minimums!}}}

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

  subroutine cice_run_advection_upwind(domain, clock)!{{{

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         tracers, &
         tracer_tendencies, &
         tracer_edge_fluxes, &
         tracer_conservation, &
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
         edgeFlux

    integer, dimension(:), pointer :: &
         interiorEdge

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracer_tendencies", tracer_tendencies)
       call MPAS_pool_get_subpool(block % structs, "tracer_edge_fluxes", tracer_edge_fluxes)
       call MPAS_pool_get_subpool(block % structs, "tracer_conservation", tracer_conservation)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)

       call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocity_solver, "normalVectorEdge", normalVectorEdge)

       call prepare_advection(&
            clock, &
            block % configs, &
            mesh, &
            tracers, &
            tracer_tendencies, &
            tracer_edge_fluxes, &
            tracer_conservation, &
            uVelocity, &
            vVelocity, &
            normalVectorEdge, &
            edgeVelocity, &
            edgeFlux)

       block => block % next
    end do

    ! advection halo exchange
    call halo_exchange_advection(domain)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "tracer_tendencies", tracer_tendencies)
       call MPAS_pool_get_subpool(block % structs, "tracer_edge_fluxes", tracer_edge_fluxes)
       call MPAS_pool_get_subpool(block % structs, "tracer_conservation", tracer_conservation)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)

       call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocity_solver, "normalVectorEdge", normalVectorEdge)
       call MPAS_pool_get_array(velocity_solver, "dynamicsTimeStep", dynamicsTimeStep)

       call MPAS_pool_get_array(boundary, "interiorEdge", interiorEdge)

       do iTracerVariable = 1, nTracerVariables

          !write(*,*) trim(tracerConnectivities(iTracerVariable) % childTracerName), trim(tracerConnectivities(iTracerVariable) % parentTracerName)

          if (tracerConnectivities(iTracerVariable) % defined == 1) then

             call run_advection_variable(&
                  clock, &
                  mesh, &
                  tracers, &
                  tracer_tendencies, &
                  tracer_edge_fluxes, &
                  trim(tracerConnectivities(iTracerVariable) % childTracerName), &
                  trim(tracerConnectivities(iTracerVariable) % parentTracerName), &
                  tracerConnectivities(iTracerVariable) % childTracerMinimum, &
                  tracerConnectivities(iTracerVariable) % parentTracerMinimum, &
                  tracerConnectivities(iTracerVariable) % updateHalo, &
                  edgeVelocity, &
                  edgeFlux, &
                  interiorEdge, &
                  dynamicsTimeStep)

          endif

       enddo ! iTracerVariable

       call finalize_advection(&
            clock, &
            block % configs, &
            mesh, &
            tracers, &
            tracer_conservation, &
            edgeVelocity, &
            edgeFlux)

       block => block % next
    end do

  end subroutine cice_run_advection_upwind!}}}

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
       clock, &
       mesh, &
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       childTracerName, &
       parentTracerName, &
       childTracerMinimum, &
       parentTracerMinimum, &
       updateHalo, &
       edgeVelocity, &
       edgeFlux, &
       InteriorEdge, &
       dt)!{{{

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type (MPAS_pool_type), pointer :: &
         tracers, &           !< Input/Output:
         tracer_tendencies, & !< Input/Output:
         tracer_edge_fluxes   !< Input/Output:

    character(len=*), intent(in) :: &
         childTracerName, & !< Input:
         parentTracerName   !< Input:

    real(kind=RKIND), intent(in) :: &
         childTracerMinimum, & !< Input:
         parentTracerMinimum   !< Input:

    logical, intent(in) :: &
         updateHalo !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input:

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output:

    integer, dimension(:), intent(in) :: &
         InteriorEdge !< Input:

    real(kind=RKIND), intent(in) :: &
         dt !< Input:

    type (mpas_pool_field_info_type) :: childFieldInfo

    call mpas_pool_get_field_info(tracers, trim(childTracerName), childFieldInfo)

    if (childFieldInfo % nDims == 3) then

       call run_advection_variable_3D(&
            clock, &
            mesh, &
            tracers, &
            tracer_tendencies, &
            tracer_edge_fluxes, &
            childTracerName, &
            parentTracerName, &
            childTracerMinimum, &
            parentTracerMinimum, &
            updateHalo, &
            edgeVelocity, &
            edgeFlux, &
            InteriorEdge, &
            dt)

    else if (childFieldInfo % nDims == 4) then

       call run_advection_variable_4D(&
            clock, &
            mesh, &
            tracers, &
            tracer_tendencies, &
            tracer_edge_fluxes, &
            childTracerName, &
            parentTracerName, &
            childTracerMinimum, &
            parentTracerMinimum, &
            updateHalo, &
            edgeVelocity, &
            edgeFlux, &
            InteriorEdge, &
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
       clock, &
       mesh, &
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       childTracerName, &
       parentTracerName, &
       childTracerMinimum, &
       parentTracerMinimum, &
       updateHalo, &
       edgeVelocity, &
       edgeFlux, &
       InteriorEdge, &
       dt)!{{{

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type (MPAS_pool_type), pointer :: &
         tracers, &           !< Input/Output:
         tracer_tendencies, & !< Input/Output:
         tracer_edge_fluxes   !< Input/Output:

    character(len=*), intent(in) :: &
         childTracerName, & !< Input:
         parentTracerName   !< Input:

    real(kind=RKIND), intent(in) :: &
         childTracerMinimum, & !< Input:
         parentTracerMinimum   !< Input:

    logical, intent(in) :: &
         updateHalo !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input:

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output:

    integer, dimension(:), intent(in) :: &
         InteriorEdge !< Input:

    real(kind=RKIND), intent(in) :: &
         dt !< Input:

    type(field3DReal), pointer :: &
         childTracerOld

    real (kind=RKIND), dimension(:,:,:), pointer :: &
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

    ! child tracer
    call MPAS_pool_get_field(tracers, trim(childTracerName), childTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(childTracerName) // "Tend", childTendency)
    call MPAS_pool_get_array(tracer_edge_fluxes, trim(childTracerName) // "EdgeFlux", childEdgeFlux)

    ! halo update
    !if (updateHalo) call mpas_dmpar_exch_halo_field(childTracerOld)

    if (trim(parentTracerName) == 'none') then

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
            clock, &
            mesh, &
            childTracerOld % array, &
            edgeVelocity, &
            parentTracerOldNone, &
            parentTracerNewNone, &
            parentTendencyNone, &
            parentEdgeFluxNone)

       call run_advection_subvariable(&
            clock, &
            mesh, &
            trim(childTracerName), &
            childTracerOld % array(:,:,:), &
            childTracerNew(:,:,:), &
            childTendency(:,:,:), &
            childEdgeFlux(:,:,:), &
            parentTracerOldNone(:,:,:), &
            parentTracerNewNone(:,:,:), &
            parentTendencyNone(:,:,:), &
            parentEdgeFluxNone(:,:,:), &
            childTracerMinimum, &
            parentTracerMinimum, &
            edgeVelocity, &
            edgeFlux, &
            InteriorEdge, &
            dt)

       deallocate(parentTracerOldNone)
       deallocate(parentTracerNewNone)
       deallocate(parentTendencyNone)
       deallocate(parentEdgeFluxNone)

    else

       ! parent tracer
       call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerOld, 1)
       call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 2)
       call MPAS_pool_get_array(tracer_tendencies, trim(parentTracerName) // "Tend", parentTendency)
       call MPAS_pool_get_array(tracer_edge_fluxes, trim(parentTracerName) // "EdgeFlux", parentEdgeFlux)

       call run_advection_subvariable(&
            clock, &
            mesh, &
            trim(childTracerName), &
            childTracerOld % array(:,:,:), &
            childTracerNew(:,:,:), &
            childTendency(:,:,:), &
            childEdgeFlux(:,:,:), &
            parentTracerOld(:,:,:), &
            parentTracerNew(:,:,:), &
            parentTendency(:,:,:), &
            parentEdgeFlux(:,:,:), &
            childTracerMinimum, &
            parentTracerMinimum, &
            edgeVelocity, &
            edgeFlux, &
            InteriorEdge, &
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
       clock, &
       mesh, &
       childTracerOld, &
       edgeVelocity, &
       parentTracerOldNone, &
       parentTracerNewNone, &
       parentTendencyNone, &
       parentEdgeFluxNone)!{{{

    use cice_constants, only: &
         iceAreaMinimum

    use cice_debug, only: &
         cice_time_output_variable_int

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

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
         cellsOnEdge

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

    do iCell = 1, nCells
       do iCategory = 1, nCategories

          parentTracerOldNone(1,iCategory,iCell) = 0.0_RKIND ! need this?
          parentTracerNewNone(1,iCategory,iCell) = 0.0_RKIND ! need this?

          parentTendencyNone(1,iCategory,iCell) = 0.0_RKIND

          if (childTracerOld(1,iCategory,iCell) > iceAreaMinimum) then
             parentTracerOldNone(1,iCategory,iCell) = 1.0_RKIND
             parentTracerNewNone(1,iCategory,iCell) = 1.0_RKIND
             !parentTendencyNone(1,iCategory,iCell) = 0.0_RKIND
          endif

          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             if (childTracerOld(1,iCategory,cellsOnCell(iEdgeOnCell,iCell)) > iceAreaMinimum) then
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
          if (childTracerOld(1,iCategory,cellsOnEdge(1,iEdge)) > iceAreaMinimum .or. &
              childTracerOld(1,iCategory,cellsOnEdge(2,iEdge)) > iceAreaMinimum) then
             parentEdgeFluxNone(1,iCategory,iEdge) = edgeVelocity(iEdge)
          else
             parentEdgeFluxNone(1,iCategory,iEdge) = 0.0_RKIND
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
       clock, &
       mesh, &
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       childTracerName, &
       parentTracerName, &
       childTracerMinimum, &
       parentTracerMinimum, &
       updateHalo, &
       edgeVelocity, &
       edgeFlux, &
       InteriorEdge, &
       dt)!{{{

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type (MPAS_pool_type), pointer :: &
         tracers, &           !< Input/Output:
         tracer_tendencies, & !< Input/Output:
         tracer_edge_fluxes   !< Input/Output:

    character(len=*), intent(in) :: &
         childTracerName, & !< Input:
         parentTracerName   !< Input:

    real(kind=RKIND), intent(in) :: &
         childTracerMinimum, &
         parentTracerMinimum

    logical, intent(in) :: &
         updateHalo !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input:

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output:

    integer, dimension(:), intent(in) :: &
         InteriorEdge

    real(kind=RKIND), intent(in) :: &
         dt !< Input:

    type(field4DReal), pointer :: &
         childTracerOld

    real (kind=RKIND), dimension(:,:,:,:), pointer :: &
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
    call MPAS_pool_get_field(tracers, trim(childTracerName), childTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(childTracerName) // "Tend", childTendency)
    call MPAS_pool_get_array(tracer_edge_fluxes, trim(childTracerName) // "EdgeFlux", childEdgeFlux)

    ! halo update
    !if (updateHalo) call mpas_dmpar_exch_halo_field(childTracerOld)

    ! parent tracer
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerOld, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 2)
    call MPAS_pool_get_array(tracer_tendencies, trim(parentTracerName) // "Tend", parentTendency)
    call MPAS_pool_get_array(tracer_edge_fluxes, trim(parentTracerName) // "EdgeFlux", parentEdgeFlux)

    do iTracerDimension = 1, size(childTracerNew,1)

       call run_advection_subvariable(&
            clock, &
            mesh, &
            trim(childTracerName), &
            childTracerOld % array(iTracerDimension,:,:,:), & ! Inefficient!
            childTracerNew(iTracerDimension,:,:,:), & ! Inefficient!
            childTendency(iTracerDimension,:,:,:), &  ! Inefficient!
            childEdgeFlux(iTracerDimension,:,:,:), &  ! Inefficient!
            parentTracerOld(:,:,:), &
            parentTracerNew(:,:,:), &
            parentTendency(:,:,:), &
            parentEdgeFlux(:,:,:), &
            childTracerMinimum, &
            parentTracerMinimum, &
            edgeVelocity, &
            edgeFlux, &
            InteriorEdge, &
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
       clock, &
       mesh, &
       childTracerName, &
       childTracerOld, &
       childTracerNew, &
       childTendency, &
       childEdgeFlux, &
       parentTracerOld, &
       parentTracerNew, &
       parentTendency, &
       parentEdgeFlux, &
       childTracerMinimum, &
       parentTracerMinimum, &
       edgeVelocity, &
       edgeFlux, &
       InteriorEdge, &
       dt)!{{{

    !use mpas_tracer_advection_mono, only: &
    !     mpas_tracer_advection_mono_tend_seaice

    use cice_debug, only: &
         cice_time_output_variable_int, &
         cice_time_output_variable_int_2D, &
         cice_time_output_variable_real, &
         cice_time_output_variable_real_2D

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type (MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    character(len=*), intent(in) :: &
         childTracerName

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

    real(kind=RKIND), intent(in) :: &
         childTracerMinimum, & !< Input
         parentTracerMinimum   !< Input

    real(kind=RKIND), dimension(:), intent(in) :: &
         edgeVelocity !< Input:

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         edgeFlux !< Input/Output:

    integer, dimension(:), intent(in) :: &
         InteriorEdge !< Input:

    real(kind=RKIND), intent(in) :: &
         dt !< Input:

    integer, pointer :: &
         nEdges, &
         nCellsSolve, &
         nCategories

    !integer, dimension(:), pointer :: &
    !     nAdvCellsForEdge

    integer, dimension(:,:), pointer :: &
         cellsOnEdge!, &
         !advCellsForEdge, &
         !highOrderAdvectionMask

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnCell, &
         edgesOnCell

    !real(kind=RKIND), dimension(:,:), pointer :: &
    !     advCoefs, &
    !     advCoefs3rd

    integer :: &
         iCell, &
         iEdge, &
         iCategory, &
         iTracerDimension

    integer, parameter :: &
         iCellTest = 7062

    integer :: &
         iCellOnCell, &
         iEdgeOnCell

    logical :: &
         monotonicityCheck

    monotonicityCheck = .true.

    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)
    call MPAS_pool_get_dimension(mesh, "nCells", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    !call MPAS_pool_get_array(mesh, "nAdvCellsForEdge", nAdvCellsForEdge)
    !call MPAS_pool_get_array(mesh, "advCellsForEdge", advCellsForEdge)
    !call MPAS_pool_get_array(mesh, "highOrderAdvectionMask", highOrderAdvectionMask)
    !call MPAS_pool_get_array(mesh, "advCoefs", advCoefs)
    !call MPAS_pool_get_array(mesh, "advCoefs3rd", advCoefs3rd)

    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "cellsOnCell", cellsOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)

    ! zero out child tendancy
    childTendency(:,:,:) = 0.0_RKIND
    childTracerNew(:,:,:) = 0.0_RKIND

    ! run the advection routine
    call upwind_tendencies(&
         mesh, &
         childTracerOld, &         ! (nTracers, nCategories, nCells)
         parentEdgeFlux(1,:,:), &  ! (nCategories, nEdges)
         parentTracerOld(1,:,:), & ! (nCategories, nCells)
         parentTracerMinimum, &
         childTendency, &          ! (nTracers, nCategories, nCells)
         childEdgeFlux, &          ! (nTracers, nCategories, nEdges)
         InteriorEdge, &
         childTracerName)

    !call mpas_tracer_advection_mono_tend_seaice(&
    !     childTracerOld,                  & ! (nTracers, nCategories, nCells)
    !     advCoefs,                        & ! (nAdvectionCells nEdges)
    !     advCoefs3rd,                     & ! (nAdvectionCells nEdges)
    !     nAdvCellsForEdge,                & ! (nEdges)
    !     advCellsForEdge,                 & ! (nAdvectionCells nEdges)
    !     parentEdgeFlux(1,:,:),           & ! (nCategories, nEdges)
    !     parentTracerOld(1,:,:),          & ! (nCategories, nCells)
    !     dt,                              & !
    !     mesh,                            & !
    !     parentTendency(1,:,:),           & ! (nCategories, nCells)
    !     parentTracerMaskOld(1,:,:),      & ! (nCategories, nCells)
    !     parentTracerMaskNew(1,:,:),      & ! (nCategories, nCells)
    !     childTendency,                   & ! (nTracers, nCategories, nCells)
    !     childEdgeFlux,                   & ! (nTracers, nCategories, nEdges)
    !     monotonicityCheck, &
    !     childTracerName)

    ! calculate the new tracer
    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories

          if (parentTracerNew(1,iCategory,iCell) > parentTracerMinimum) then

             do iTracerDimension = 1, size(childTracerNew,1)

                ! calculate the new tracer
                childTracerNew(iTracerDimension,iCategory,iCell) = &
                     childTracerOld(iTracerDimension,iCategory,iCell) * parentTracerOld(1,iCategory,iCell) + childTendency(iTracerDimension,iCategory,iCell) * dt

                ! store the old child tracer multiplied by the old parent. This will be the parent tracer of the next child
                childTracerOld(iTracerDimension,iCategory,iCell) = childTracerOld(iTracerDimension,iCategory,iCell) * parentTracerOld(iTracerDimension,iCategory,iCell)

             enddo ! iTracer

          endif ! parentTracerMaskNew

       enddo ! iCategory
    enddo ! iCell

  end subroutine run_advection_subvariable!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date September 16th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine upwind_tendencies(&
       mesh, &
       childTracerOld, &      ! (nTracers, nCategories, nCells)
       parentEdgeFlux, &      ! (nCategories, nEdges)
       parentTracerOld, &     ! (nCategories, nCells)
       parentTracerMinimum, &
       childTendency, &       ! (nTracers, nCategories, nCells)
       childEdgeFlux, &       ! (nTracers, nCategories, nEdges)
       interiorEdge, &
       childTracerName)

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         childTracerOld

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         parentEdgeFlux, &
         parentTracerOld

    real(kind=RKIND), intent(in) :: &
         parentTracerMinimum

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         childTendency, &
         childEdgeFlux

    integer, dimension(:), intent(in) :: &
         interiorEdge

    character(len=*), intent(in) :: &
         childTracerName

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iCategory, &
         iEdge, &
         cell1, &
         cell2, &
         iTracer

    integer, pointer :: &
         nCellsSolve, &
         nEdges, &
         nCategories, &
         maxEdges

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell, &
         cellsOnEdge

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell, &
         dvEdge

    real(kind=RKIND) :: &
         invAreaCell1, &
         flux_upwind

    integer, dimension(:,:), allocatable :: &
         edgeSignOnCell

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)

    allocate(edgeSignOnCell(maxEdges, nCellsSolve))
    do iCell = 1, nCellsSolve
       do iEdgeOnCell = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(iEdgeOnCell, iCell)
          if(iCell == cellsOnEdge(1, iEdge)) then
             edgeSignOnCell(iEdgeOnCell, iCell) = -1
          else
             edgeSignOnCell(iEdgeOnCell, iCell) =  1
          end if
       end do
    end do

    do iTracer = 1, size(childTracerOld,dim=1)

       do iCell = 1, nCellsSolve
          do iCategory = 1, nCategories
             childTendency(iTracer,iCategory,iCell) = 0.0_RKIND
          enddo ! iCategory
       enddo ! iCell

       do iEdge = 1, nEdges
          do iCategory = 1, nCategories
             childEdgeFlux(iTracer,iCategory,iEdge) = 0.0_RKIND
          enddo ! iCategory
       enddo ! iCell

       do iCell = 1, nCellsSolve

          invAreaCell1 = 1.0_RKIND / areaCell(iCell)

          do iEdgeOnCell = 1, nEdgesOnCell(iCell)

             iEdge = edgesOnCell(iEdgeOnCell, iCell)
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)

             if (interiorEdge(iEdge) == 1) then

                do iCategory = 1, nCategories

                   if (parentTracerOld(iCategory,cell1) > parentTracerMinimum .or. &
                       parentTracerOld(iCategory,cell2) > parentTracerMinimum) then

                      flux_upwind = dvEdge(iEdge) * (max(0.0_RKIND,parentEdgeFlux(iCategory,iEdge))*childTracerOld(iTracer,iCategory,cell1) + &
                                                     min(0.0_RKIND,parentEdgeFlux(iCategory,iEdge))*childTracerOld(iTracer,iCategory,cell2))

                      childTendency(iTracer,iCategory,iCell) = childTendency(iTracer,iCategory,iCell) + &
                           edgeSignOncell(iEdgeOnCell,iCell) * flux_upwind * invAreaCell1

                      childEdgeFlux(iTracer,iCategory,iEdge) = flux_upwind / dvEdge(iEdge)

                      !childTendency(iTracer,iCategory,iCell) = 0.0_RKIND
                      !childEdgeFlux(iTracer,iCategory,iEdge) = 0.0_RKIND

                   endif

                enddo ! iCategory

             endif ! interiorEdge

          enddo ! iEdgeOnCell

       enddo ! iCell

    enddo ! iTracer

    deallocate(edgeSignOnCell)

  end subroutine upwind_tendencies

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
         nCells, &
         nCellsSolve

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell, &
         cellsOnEdge, &
         verticesOnEdge

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)

    edgeVelocity = 0.0_RKIND

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

    ! loop over cells
    !do iEdge = 1, nEdges

       ! get edges surrounding owned cells
    !   if (cellsOnEdge(1,iEdge) <= nCellsSolve .or. &
    !       cellsOnEdge(2,iEdge) <= nCellsSolve) then

    !      iCell = cellsOnEdge(1,iEdge)
    !      do iEdgeOnCell1 = 1, nEdgesOnCell(iCell)
    !         if (cellsOnCell(iEdgeOnCell,iCell) == cellsOnEdge(2,iEdge)) then
    !            iEdgeOnCell = iEdgeOnCell1
    !            exit
    !         endif
    !      enddo ! iEdgeOnCell

          ! find u,v velocity at edge
    !      uVelocityEdge = 0.0_RKIND
    !      vVelocityEdge = 0.0_RKIND

    !      do iVertexOnEdge = 1, 2

    !         iVertex = verticesOnEdge(iVertexOnEdge, iEdge)

    !         uVelocityEdge = uVelocityEdge + uVelocity(iVertex)
    !         vVelocityEdge = vVelocityEdge + vVelocity(iVertex)

    !      enddo ! iVertexOnEdge

    !      uVelocityEdge = uVelocityEdge / 2.0_RKIND
    !      vVelocityEdge = vVelocityEdge / 2.0_RKIND

          ! rotate u,v velocity to normal
    !      edgeVelocity(iEdge) = &
    !           uVelocityEdge * normalVectorEdge(1,iEdgeOnCell,iCell) + &
    !           vVelocityEdge * normalVectorEdge(2,iEdgeOnCell,iCell)

    !   endif ! interior block edge

    !enddo ! iEdge



  end subroutine edge_from_vertex_velocity!}}}

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
       clock, &
       configs, &
       mesh, &
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       tracer_conservation, &
       uVelocity, &
       vVelocity, &
       normalVectorEdge, &
       edgeVelocity, &
       edgeFlux)!{{{

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, &  !< Input:
         mesh        !< Input:

    type(MPAS_pool_type), pointer :: &
         tracers, &            !< Input/Output:
         tracer_tendencies, &  !< Input/Output:
         tracer_edge_fluxes, & !< Input/Output:
         tracer_conservation   !< Input/Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input:
         vVelocity    !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         normalVectorEdge !< Input:

    real(kind=RKIND), dimension(:), allocatable, intent(inout) :: &
         edgeVelocity !< Input/Output:

    real(kind=RKIND), dimension(:,:), allocatable, intent(inout) :: &
         edgeFlux !< Input/Output:

    integer, pointer :: &
         nEdges, &
         nCategories

    call initialize_timelevel_variables(&
         tracers, &
         tracer_tendencies, &
         tracer_edge_fluxes, &
         tracer_conservation)

    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    allocate(edgeVelocity(nEdges))
    allocate(edgeFlux(nCategories,nEdges))

    call edge_from_vertex_velocity(&
         mesh, &
         uVelocity, &
         vVelocity, &
         normalVectorEdge, &
         edgeVelocity)

    call prepare_tracers(&
         clock, &
         configs, &
         mesh, &
         tracers, &
         tracer_conservation)

  end subroutine prepare_advection!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  initialize_timelevel_variables
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date August 18th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine initialize_timelevel_variables(&
       tracers, &
       tracer_tendencies, &
       tracer_edge_fluxes, &
       tracer_conservation)

    type(MPAS_pool_type), intent(inout) :: &
         tracers, &
         tracer_tendencies, &
         tracer_edge_fluxes, &
         tracer_conservation

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory2, &
         iceVolumeCategory2, &
         snowVolumeCategory2, &
         surfaceTemperature2, &
         iceEnthalpy2, &
         iceSalinity2, &
         snowEnthalpy2

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryCons1, &
         iceVolumeCategoryCons1, &
         snowVolumeCategoryCons1, &
         surfaceTemperatureCons1, &
         iceEnthalpyCons1, &
         iceSalinityCons1, &
         snowEnthalpyCons1

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryCons2, &
         iceVolumeCategoryCons2, &
         snowVolumeCategoryCons2, &
         surfaceTemperatureCons2, &
         iceEnthalpyCons2, &
         iceSalinityCons2, &
         snowEnthalpyCons2

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategoryTend1, &
         iceVolumeCategoryTend1, &
         snowVolumeCategoryTend1, &
         surfaceTemperatureTend1, &
         iceEnthalpyTend1, &
         iceSalinityTend1, &
         snowEnthalpyTend1

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategoryEdgeFlux1, &
         iceVolumeCategoryEdgeFlux1, &
         snowVolumeCategoryEdgeFlux1, &
         surfaceTemperatureEdgeFlux1, &
         iceEnthalpyEdgeFlux1, &
         iceSalinityEdgeFlux1, &
         snowEnthalpyEdgeFlux1

    ! second time level
    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory2, 2)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory2, 2)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory2, 2)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature2, 2)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy2, 2)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity2, 2)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy2, 2)

    call MPAS_pool_get_array(tracer_conservation, "iceAreaCategoryCons", iceAreaCategoryCons2, 2)
    call MPAS_pool_get_array(tracer_conservation, "iceVolumeCategoryCons", iceVolumeCategoryCons2, 2)
    call MPAS_pool_get_array(tracer_conservation, "snowVolumeCategoryCons", snowVolumeCategoryCons2, 2)
    call MPAS_pool_get_array(tracer_conservation, "surfaceTemperatureCons", surfaceTemperatureCons2, 2)
    call MPAS_pool_get_array(tracer_conservation, "iceEnthalpyCons", iceEnthalpyCons2, 2)
    call MPAS_pool_get_array(tracer_conservation, "iceSalinityCons", iceSalinityCons2, 2)
    call MPAS_pool_get_array(tracer_conservation, "snowEnthalpyCons", snowEnthalpyCons2, 2)

    iceAreaCategory2 = 0.0_RKIND
    iceVolumeCategory2 = 0.0_RKIND
    snowVolumeCategory2 = 0.0_RKIND
    surfaceTemperature2 = 0.0_RKIND
    iceEnthalpy2 = 0.0_RKIND
    iceSalinity2 = 0.0_RKIND
    snowEnthalpy2 = 0.0_RKIND

    iceAreaCategoryCons2 = 0.0_RKIND
    iceVolumeCategoryCons2 = 0.0_RKIND
    snowVolumeCategoryCons2 = 0.0_RKIND
    surfaceTemperatureCons2 = 0.0_RKIND
    iceEnthalpyCons2 = 0.0_RKIND
    iceSalinityCons2 = 0.0_RKIND
    snowEnthalpyCons2 = 0.0_RKIND

    ! first time level
    call MPAS_pool_get_array(tracer_conservation, "iceAreaCategoryCons", iceAreaCategoryCons1, 1)
    call MPAS_pool_get_array(tracer_conservation, "iceVolumeCategoryCons", iceVolumeCategoryCons1, 1)
    call MPAS_pool_get_array(tracer_conservation, "snowVolumeCategoryCons", snowVolumeCategoryCons1, 1)
    call MPAS_pool_get_array(tracer_conservation, "surfaceTemperatureCons", surfaceTemperatureCons1, 1)
    call MPAS_pool_get_array(tracer_conservation, "iceEnthalpyCons", iceEnthalpyCons1, 1)
    call MPAS_pool_get_array(tracer_conservation, "iceSalinityCons", iceSalinityCons1, 1)
    call MPAS_pool_get_array(tracer_conservation, "snowEnthalpyCons", snowEnthalpyCons1, 1)

    call MPAS_pool_get_array(tracer_tendencies, "iceAreaCategoryTend", iceAreaCategoryTend1, 1)
    call MPAS_pool_get_array(tracer_tendencies, "iceVolumeCategoryTend", iceVolumeCategoryTend1, 1)
    call MPAS_pool_get_array(tracer_tendencies, "snowVolumeCategoryTend", snowVolumeCategoryTend1, 1)
    call MPAS_pool_get_array(tracer_tendencies, "surfaceTemperatureTend", surfaceTemperatureTend1, 1)
    call MPAS_pool_get_array(tracer_tendencies, "iceEnthalpyTend", iceEnthalpyTend1, 1)
    call MPAS_pool_get_array(tracer_tendencies, "iceSalinityTend", iceSalinityTend1, 1)
    call MPAS_pool_get_array(tracer_tendencies, "snowEnthalpyTend", snowEnthalpyTend1, 1)

    call MPAS_pool_get_array(tracer_edge_fluxes, "iceAreaCategoryEdgeFlux", iceAreaCategoryEdgeFlux1, 1)
    call MPAS_pool_get_array(tracer_edge_fluxes, "iceVolumeCategoryEdgeFlux", iceVolumeCategoryEdgeFlux1, 1)
    call MPAS_pool_get_array(tracer_edge_fluxes, "snowVolumeCategoryEdgeFlux", snowVolumeCategoryEdgeFlux1, 1)
    call MPAS_pool_get_array(tracer_edge_fluxes, "surfaceTemperatureEdgeFlux", surfaceTemperatureEdgeFlux1, 1)
    call MPAS_pool_get_array(tracer_edge_fluxes, "iceEnthalpyEdgeFlux", iceEnthalpyEdgeFlux1, 1)
    call MPAS_pool_get_array(tracer_edge_fluxes, "iceSalinityEdgeFlux", iceSalinityEdgeFlux1, 1)
    call MPAS_pool_get_array(tracer_edge_fluxes, "snowEnthalpyEdgeFlux", snowEnthalpyEdgeFlux1, 1)

    iceAreaCategoryCons1 = 0.0_RKIND
    iceVolumeCategoryCons1 = 0.0_RKIND
    snowVolumeCategoryCons1 = 0.0_RKIND
    surfaceTemperatureCons1 = 0.0_RKIND
    iceEnthalpyCons1 = 0.0_RKIND
    iceSalinityCons1 = 0.0_RKIND
    snowEnthalpyCons1 = 0.0_RKIND

    iceAreaCategoryTend1 = 0.0_RKIND
    iceVolumeCategoryTend1 = 0.0_RKIND
    snowVolumeCategoryTend1 = 0.0_RKIND
    surfaceTemperatureTend1 = 0.0_RKIND
    iceEnthalpyTend1 = 0.0_RKIND
    iceSalinityTend1 = 0.0_RKIND
    snowEnthalpyTend1 = 0.0_RKIND

    iceAreaCategoryEdgeFlux1 = 0.0_RKIND
    iceVolumeCategoryEdgeFlux1 = 0.0_RKIND
    snowVolumeCategoryEdgeFlux1 = 0.0_RKIND
    surfaceTemperatureEdgeFlux1 = 0.0_RKIND
    iceEnthalpyEdgeFlux1 = 0.0_RKIND
    iceSalinityEdgeFlux1 = 0.0_RKIND
    snowEnthalpyEdgeFlux1 = 0.0_RKIND

  end subroutine initialize_timelevel_variables

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  halo_exchange_advection
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date January 13th 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine halo_exchange_advection(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(MPAS_pool_type), pointer :: &
         tracers_fields

    type (MPAS_pool_field_info_type) :: &
         childFieldInfo

    type(field3DReal), pointer :: tracer3DReal
    type(field4DReal), pointer :: tracer4DReal

    integer :: &
         iTracerVariable

    call MPAS_pool_get_subpool(domain % blocklist % structs, "tracers", tracers_fields)

    do iTracerVariable = 1, nTracerVariables

       if (tracerConnectivities(iTracerVariable) % defined == 1 .and. &
           tracerConnectivities(iTracerVariable) % updateHalo) then

          call MPAS_pool_get_field_info(tracers_fields, trim(tracerConnectivities(iTracerVariable) % childTracerName), childFieldInfo)

          select case (childFieldInfo % nDims)
          case (3)
             call MPAS_pool_get_field(tracers_fields, trim(tracerConnectivities(iTracerVariable) % childTracerName), tracer3DReal)
             call MPAS_dmpar_exch_halo_field(tracer3DReal)
          case (4)
             call MPAS_pool_get_field(tracers_fields, trim(tracerConnectivities(iTracerVariable) % childTracerName), tracer4DReal)
             call MPAS_dmpar_exch_halo_field(tracer4DReal)
          end select

       endif

    enddo ! iTracerVariable

  end subroutine halo_exchange_advection

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
       clock, &
       configs, &
       mesh, &
       tracers, &
       tracer_conservation, &
       edgeVelocity, &
       edgeFlux)!{{{

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, &
         mesh !< Input:

    type(MPAS_pool_type), pointer :: &
         tracers, &          !< Input/Output:
         tracer_conservation !< Input/Output:

    real(kind=RKIND), dimension(:), allocatable, intent(inout) :: &
         edgeVelocity !< Input/Output:

    real(kind=RKIND), dimension(:,:), allocatable, intent(inout) :: &
         edgeFlux !< Input/Output:

    call finalize_tracers(&
         clock, &
         configs, &
         mesh, &
         tracers, &
         tracer_conservation)

    deallocate(edgeVelocity)
    deallocate(edgeFlux)

  end subroutine finalize_advection!}}}

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
       clock, &
       configs, &
       mesh, &
       tracers, &
       tracer_conservation)!{{{

    use cice_constants, only: &
         iceAreaMinimum

    use cice_debug, only: cice_time_output_variable_real

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, & !< Input:
         mesh       !< Input:

    type(MPAS_pool_type), pointer :: &
         tracers, &          !< Input/Output:
         tracer_conservation !< Input/Output:

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature, &
         iceVolumeCategoryOther

    integer, pointer :: &
         nCellsSolve, &
         nCells, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    logical, pointer :: &
         config_convert_volume_to_thickness, &
         config_conservation_check

    integer, dimension(:), pointer :: &
         itimestep

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_config(configs, "config_convert_volume_to_thickness", config_convert_volume_to_thickness)
    call MPAS_pool_get_config(configs, "config_conservation_check", config_conservation_check)

    if (config_convert_volume_to_thickness) then

       do iCell = 1, nCells
          do iCategory = 1, nCategories

             if (iceAreaCategory(1,iCategory,iCell) > iceAreaMinimum) then

                iceVolumeCategory(1,iCategory,iCell) = iceVolumeCategory(1,iCategory,iCell) / &
                     iceAreaCategory(1,iCategory,iCell)

                snowVolumeCategory(1,iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell) / &
                     iceAreaCategory(1,iCategory,iCell)

             endif ! area > 0

          enddo ! iCategory
       enddo ! iCell

    endif ! config_convert_volume_to_thickness

    if (config_conservation_check) &
         call initial_conservation(&
            mesh, &
            tracers, &
            tracer_conservation)

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
       clock, &
       configs, &
       mesh, &
       tracers, &
       tracer_conservation)!{{{

    use cice_constants, only: &
         iceAreaMinimum

    use cice_debug, only: cice_time_output_variable_real

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs, & !< Input:
         mesh       !< Input:

    type(MPAS_pool_type), pointer :: &
         tracers, &          !< Input/Output:
         tracer_conservation !< Input/Output:

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         surfaceTemperature, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceVolumeCategoryOther

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         iCategory

    logical, pointer :: &
         config_convert_volume_to_thickness, &
         config_limit_ice_concentration, &
         config_conservation_check, &
         config_clean_tracers

    integer, dimension(:), pointer :: &
         itimestep

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    call MPAS_pool_shift_time_levels(tracers)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_config(configs, "config_convert_volume_to_thickness", config_convert_volume_to_thickness)
    call MPAS_pool_get_config(configs, "config_limit_ice_concentration", config_limit_ice_concentration)
    call MPAS_pool_get_config(configs, "config_conservation_check", config_conservation_check)
    call MPAS_pool_get_config(configs, "config_clean_tracers", config_clean_tracers)

    if (config_conservation_check) &
         call final_conservation(&
            mesh, &
            tracers, &
            tracer_conservation)

    if (config_limit_ice_concentration) &
         call limit_ice_concentration(mesh, tracers)

    call scale_tracers_back(mesh, tracers)

    if (config_convert_volume_to_thickness) then

       do iCell = 1, nCellsSolve
          do iCategory = 1, nCategories

             if (iceAreaCategory(1,iCategory,iCell) > iceAreaMinimum) then

                iceVolumeCategory(1,iCategory,iCell) = iceVolumeCategory(1,iCategory,iCell) * &
                     iceAreaCategory(1,iCategory,iCell)

                snowVolumeCategory(1,iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell) * &
                     iceAreaCategory(1,iCategory,iCell)

             endif ! area > 0

          enddo ! iCategory
       enddo ! iCell

    endif ! config_convert_volume_to_thickness

    if (config_clean_tracers) &
         call clean_tracers(clock, mesh, tracers)

  end subroutine finalize_tracers!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  limit_ice_concentration
!
!> \brief Limit the ice concentraion to lie in interval [0,1]
!> \author Adrian K. Turner, LANL
!> \date July 24th 2014
!> \details Instead of mechanical redistribution of areas greater than 1.0
!>  limit ice concentration before other tracers calculated after advection
!
!-----------------------------------------------------------------------

  subroutine limit_ice_concentration(mesh, tracers)

    type(MPAS_pool_type), intent(in) :: &
         mesh

    type(MPAS_pool_type) :: &
         tracers

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
        iceAreaCategory

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)

    do iCell = 1, nCellsSolve

       iceAreaCategory(:,:,iCell) = min(max(iceAreaCategory(:,:,iCell),0.0_RKIND),1.0_RKIND)

    enddo ! iCell

  end subroutine limit_ice_concentration

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
       tracers)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type(MPAS_pool_type), pointer :: &
         tracers !< Input/Output:

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
                  trim(tracerConnectivities(iTracerVariable) % childTracerName), &
                  trim(tracerConnectivities(iTracerVariable) % parentTracerName))

          else if (childFieldInfo % nDims == 4) then

             call scale_tracers_back_4D(&
                  mesh, &
                  tracers, &
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
       childTracerName, &
       parentTracerName)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type(MPAS_pool_type), pointer :: &
         tracers !< Input/Output:

    character(len=*), intent(in) :: &
         childTracerName, & !< Input:
         parentTracerName   !< Input:

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         childTracerNew, &
         parentTracerNew

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         iCategory, &
         iTracer

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 1)

    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories
          do iTracer = 1, size(childTracerNew,1)

             if (parentTracerNew(1,iCategory,iCell) > 0.0_RKIND) then

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
       childTracerName, &
       parentTracerName)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type(MPAS_pool_type), pointer :: &
         tracers !< Input/Output:

    character(len=*), intent(in) :: &
         childTracerName, & !< Input:
         parentTracerName   !< Input:

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         parentTracerNew

    real(kind=RKIND), dimension(:,:,:,:), pointer :: &
         childTracerNew

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         iCategory, &
         iTracer, &
         iTracerDimension

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(tracers, trim(childTracerName), childTracerNew, 1)
    call MPAS_pool_get_array(tracers, trim(parentTracerName), parentTracerNew, 1)

    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories
          do iTracer = 1, size(childTracerNew,2)
             do iTracerDimension = 1, size(childTracerNew,1)

                if (parentTracerNew(1,iCategory,iCell) > 0.0_RKIND) then

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

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  clean_tracers
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date August 4th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine clean_tracers(clock, mesh, tracers)

    use cice_constants, only: &
         iceAreaMinimum

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(MPAS_pool_type), pointer :: &
         mesh

    type(MPAS_pool_type), intent(inout), pointer :: &
         tracers

    integer :: &
         iCell, &
         iCategory

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    type(MPAS_Time_type) :: &
         curr_time

    character(len=strKIND) :: &
         dateTimeString

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories

          ! fail for negative ice area
          if (iceAreaCategory(1,iCategory,iCell) < 0.0_RKIND) then
             curr_time = mpas_get_clock_time(clock, MPAS_NOW)
             call mpas_get_time(curr_time, dateTimeString=dateTimeString)
             write(*,*) "Negative ice area! iCell, iCategory, time: ", iCell, iCategory, dateTimeString
             write(*,*) "iceAreaCategory:    ", iceAreaCategory(1,iCategory,iCell)
             write(*,*) "iceVolumeCategory:  ", iceVolumeCategory(1,iCategory,iCell)
             write(*,*) "snowVolumeCategory: ", snowVolumeCategory(1,iCategory,iCell)
             write(*,*) "surfaceTemperature: ", surfaceTemperature(1,iCategory,iCell)
             stop
          endif

          ! fail for negative ice volume
          if (iceVolumeCategory(1,iCategory,iCell) < 0.0_RKIND) then
             curr_time = mpas_get_clock_time(clock, MPAS_NOW)
             call mpas_get_time(curr_time, dateTimeString=dateTimeString)
             write(*,*) "Negative ice volume! iCell, iCategory, time: ", iCell, iCategory, dateTimeString
             write(*,*) "iceAreaCategory:    ", iceAreaCategory(1,iCategory,iCell)
             write(*,*) "iceVolumeCategory:  ", iceVolumeCategory(1,iCategory,iCell)
             write(*,*) "snowVolumeCategory: ", snowVolumeCategory(1,iCategory,iCell)
             write(*,*) "surfaceTemperature: ", surfaceTemperature(1,iCategory,iCell)
             stop
          endif

          ! fail for negative snow volume
          if (snowVolumeCategory(1,iCategory,iCell) < 0.0_RKIND) then
             curr_time = mpas_get_clock_time(clock, MPAS_NOW)
             call mpas_get_time(curr_time, dateTimeString=dateTimeString)
             write(*,*) "Negative snow volume! iCell, iCategory, time: ", iCell, iCategory, dateTimeString
             write(*,*) "iceAreaCategory:    ", iceAreaCategory(1,iCategory,iCell)
             write(*,*) "iceVolumeCategory:  ", iceVolumeCategory(1,iCategory,iCell)
             write(*,*) "snowVolumeCategory: ", snowVolumeCategory(1,iCategory,iCell)
             write(*,*) "surfaceTemperature: ", surfaceTemperature(1,iCategory,iCell)
             stop
          endif

          ! remove small areas
          if (iceAreaCategory(1,iCategory,iCell) < iceAreaMinimum) then
             iceAreaCategory(1,iCategory,iCell) = 0.0_RKIND
             iceVolumeCategory(1,iCategory,iCell) = 0.0_RKIND
             snowVolumeCategory(1,iCategory,iCell) = 0.0_RKIND
             surfaceTemperature(1,iCategory,iCell) = 0.0_RKIND
          end if

          ! remove small volume
          if (iceVolumeCategory(1,iCategory,iCell) <= 0.0_RKIND) then
             iceAreaCategory(1,iCategory,iCell) = 0.0_RKIND
             iceVolumeCategory(1,iCategory,iCell) = 0.0_RKIND
             snowVolumeCategory(1,iCategory,iCell) = 0.0_RKIND
             surfaceTemperature(1,iCategory,iCell) = 0.0_RKIND
          endif

       enddo ! iCategory
    enddo ! iCell

  end subroutine clean_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  initial_conservation
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date July 29th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine initial_conservation(&
       mesh, &
       tracers, &
       tracer_conservation)

    type(MPAS_pool_type), pointer :: &
         mesh

    type(MPAS_pool_type), intent(inout), pointer :: &
         tracers, &
         tracer_conservation

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         childTracer, &
         parentTracer, &
         childParentTracer

    integer :: &
         iTracerVariable, &
         iCell, &
         iTracer, &
         iCategory

    integer, pointer :: &
         nCategories, &
         nCellsSolve

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         tracerAggregate

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    do iTracerVariable = 1, nTracerVariables

       if (tracerConnectivities(iTracerVariable) % defined == 1) then

          call MPAS_pool_get_array(tracers, trim(tracerConnectivities(iTracerVariable) % childTracerName), childTracer, 1)
          call MPAS_pool_get_array(tracers, trim(tracerConnectivities(iTracerVariable) % childTracerName), childParentTracer, 2)

          call MPAS_pool_get_array(tracer_conservation, trim(tracerConnectivities(iTracerVariable) % childTracerName) // "Cons", tracerAggregate, 1)

          tracerAggregate = 0.0_RKIND

          if (trim(tracerConnectivities(iTracerVariable) % parentTracerName) == "none") then

             do iCell = 1, nCellsSolve
                do iCategory = 1, nCategories
                   do iTracer = 1, size(tracerAggregate,1)

                      tracerAggregate(iTracer,iCategory) = tracerAggregate(iTracer,iCategory) + &
                           childTracer(iTracer,iCategory,iCell) * areaCell(iCell)

                      childParentTracer(iTracer,iCategory,iCell) = &
                           childTracer(iTracer,iCategory,iCell)

                   enddo ! iTracer
                enddo ! iCategory
             enddo ! iCell

          else

             call MPAS_pool_get_array(tracers, trim(tracerConnectivities(iTracerVariable) % parentTracerName), parentTracer, 2)

             do iCell = 1, nCellsSolve
                do iCategory = 1, nCategories
                   do iTracer = 1, size(tracerAggregate,1)

                      tracerAggregate(iTracer,iCategory) = tracerAggregate(iTracer,iCategory) + &
                           childTracer(iTracer,iCategory,iCell) * parentTracer(iTracer,iCategory,iCell) * areaCell(iCell)

                      childParentTracer(iTracer,iCategory,iCell) = &
                           parentTracer(iTracer,iCategory,iCell) * childTracer(iTracer,iCategory,iCell)

                   enddo ! iTracer
                enddo ! iCategory
             enddo ! iCell

          end if ! tracer poarent none

       end if ! tracer defined

    enddo ! iTracerVariable

    ! erase

    do iTracerVariable = 1, nTracerVariables

       if (tracerConnectivities(iTracerVariable) % defined == 1) then

          call MPAS_pool_get_array(tracers, trim(tracerConnectivities(iTracerVariable) % childTracerName), childParentTracer, 2)
          childParentTracer = 0.0_RKIND

       endif ! tracer defined

    enddo ! iTracerVariable

  end subroutine initial_conservation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  final_conservation
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date July 29th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine final_conservation(&
       mesh, &
       tracers, &
       tracer_conservation)

    type(MPAS_pool_type), pointer :: &
         mesh

    type(MPAS_pool_type), intent(inout), pointer :: &
         tracers, &
         tracer_conservation

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         tracer

    real(kind=RKIND), dimension(:,:), pointer :: &
         tracerAggregate, &
         initialTracerAggregate

    integer, pointer :: &
         nCategories, &
         nCellsSolve

    integer :: &
         iTracerVariable, &
         iCell, &
         iCategory, &
         iTracer

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    character(len=30) :: &
         strTracerName, &
         strNameHeader

    character(len=10) :: &
         strCategoryHeader, &
         strTracerHeader

    character(len=18) :: &
         strInitHeader, &
         strFinalHeader, &
         strDiffHeader, &
         strRatioHeader

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    strNameHeader     = "Tracer Name"
    strCategoryHeader = "iCategory"
    strTracerHeader   = "iTracer"
    strInitHeader     = "Initial"
    strFinalHeader    = "Final"
    strDiffHeader     = "Difference"
    strRatioHeader    = "Ratio"

    write(65,fmt='(a30,3x,a10,3x,a10,3x,a18,3x,a18,3x,a18,3x,a18)') &
         adjustl(strNameHeader), adjustl(strCategoryHeader), adjustl(strTracerHeader), &
         adjustl(strInitHeader), adjustl(strFinalHeader), adjustl(strDiffHeader), adjustl(strRatioHeader)

    do iTracerVariable = 1, nTracerVariables

       if (tracerConnectivities(iTracerVariable) % defined == 1) then

          call MPAS_pool_get_array(tracers, trim(tracerConnectivities(iTracerVariable) % childTracerName), tracer, 1)

          call MPAS_pool_get_array(tracer_conservation, trim(tracerConnectivities(iTracerVariable) % childTracerName) // "Cons", tracerAggregate, 2)
          call MPAS_pool_get_array(tracer_conservation, trim(tracerConnectivities(iTracerVariable) % childTracerName) // "Cons", initialTracerAggregate, 1)

          tracerAggregate = 0.0_RKIND

          do iCell = 1, nCellsSolve
             do iCategory = 1, nCategories
                do iTracer = 1, size(tracerAggregate,1)

                   tracerAggregate(iTracer,iCategory) = tracerAggregate(iTracer,iCategory) + &
                        tracer(iTracer,iCategory,iCell) * areaCell(iCell)

                enddo ! iTracer
             enddo ! iCategory
          enddo ! iCell

          do iCategory = 1, nCategories
             do iTracer = 1, size(tracerAggregate,1)
                write(strTracerName,fmt='(a30)') trim(tracerConnectivities(iTracerVariable) % childTracerName)
                if (initialTracerAggregate(iTracer,iCategory) /= 0.0_RKIND) then
                   write(65,fmt='(a30,3x,i10,3x,i10,3x,e18.10,3x,e18.10,3x,e18.10,3x,e18.10)') &
                        adjustl(strTracerName), iCategory, iTracer, &
                        initialTracerAggregate(iTracer,iCategory), tracerAggregate(iTracer,iCategory), &
                        tracerAggregate(iTracer,iCategory) - initialTracerAggregate(iTracer,iCategory), &
                        (tracerAggregate(iTracer,iCategory) - initialTracerAggregate(iTracer,iCategory)) / initialTracerAggregate(iTracer,iCategory)
                else
                   write(65,fmt='(a30,3x,i10,3x,i10,3x,e18.10,3x,e18.10,3x,e18.10,3x,a)') &
                        adjustl(strTracerName), iCategory, iTracer, &
                        initialTracerAggregate(iTracer,iCategory), tracerAggregate(iTracer,iCategory), &
                        tracerAggregate(iTracer,iCategory) - initialTracerAggregate(iTracer,iCategory), &
                        "----"
                endif
             enddo ! iTracer
          enddo ! iCategory

       endif ! tracer defined

    enddo ! iTracerVariable

  end subroutine final_conservation

!-----------------------------------------------------------------------

end module cice_advection_upwind
