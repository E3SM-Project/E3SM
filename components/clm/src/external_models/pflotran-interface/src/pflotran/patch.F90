module Patch_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Grid_module
  use Coupler_module
  use Observation_module
  use Integral_Flux_module
  use Strata_module
  use Region_module
  use Reaction_Aux_module
  use Dataset_Base_class
  use Material_module
  use Field_module
  use Characteristic_Curves_module
  use Surface_Field_module
  use Surface_Material_module
  use Surface_Auxiliary_module

  use Auxiliary_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: patch_type

    PetscInt :: id

    ! These arrays will be used by all modes, mode-specific arrays should
    ! go in the auxiliary data stucture for that mode
    PetscInt, pointer :: imat(:)
    PetscInt, pointer :: imat_internal_to_external(:)
    PetscInt, pointer :: sat_func_id(:)

    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)
    PetscReal, pointer :: internal_tran_coefs(:,:,:)
    PetscReal, pointer :: boundary_tran_coefs(:,:,:)
    PetscReal, pointer :: internal_flow_fluxes(:,:)
    PetscReal, pointer :: boundary_flow_fluxes(:,:)
    ! fluid fluxes in moles/sec
    PetscReal, pointer :: ss_flow_fluxes(:,:)
    ! volumetric flux (m^3/sec) for liquid phase needed for transport
    PetscReal, pointer :: ss_flow_vol_fluxes(:,:)
    PetscReal, pointer :: internal_tran_fluxes(:,:)
    PetscReal, pointer :: boundary_tran_fluxes(:,:)
    PetscReal, pointer :: ss_tran_fluxes(:,:)

    ! for TH surface/subsurface
    PetscReal, pointer :: boundary_energy_flux(:,:)

    type(grid_type), pointer :: grid

    type(region_list_type), pointer :: region_list

    type(coupler_list_type), pointer :: boundary_condition_list
    type(coupler_list_type), pointer :: initial_condition_list
    type(coupler_list_type), pointer :: source_sink_list

    type(material_property_type), pointer :: material_properties
    type(material_property_ptr_type), pointer :: material_property_array(:)
    class(characteristic_curves_type), pointer :: characteristic_curves
    type(characteristic_curves_ptr_type), pointer :: characteristic_curves_array(:)

    type(strata_list_type), pointer :: strata_list
    type(observation_list_type), pointer :: observation_list
    type(integral_flux_list_type), pointer :: integral_flux_list

    ! Pointers to objects in mother realization object
    type(field_type), pointer :: field
    type(reaction_type), pointer :: reaction
    class(dataset_base_type), pointer :: datasets

    type(auxiliary_type) :: aux

    type(patch_type), pointer :: next

    PetscInt :: surf_or_subsurf_flag  ! Flag to identify if the current patch
                                      ! is a surface or subsurface (default)
    type(surface_material_property_type), pointer :: surf_material_properties
    type(surface_material_property_ptr_type), pointer :: surf_material_property_array(:)
    type(surface_field_type), pointer :: surf_field
    type(surface_auxiliary_type) :: surf_aux

  end type patch_type

  ! pointer data structure required for making an array of patch pointers in F90
  type, public :: patch_ptr_type
    type(patch_type), pointer :: ptr           ! pointer to the patch_type
  end type patch_ptr_type

  type, public :: patch_list_type
    PetscInt :: num_patch_objects
    type(patch_type), pointer :: first
    type(patch_type), pointer :: last
    type(patch_ptr_type), pointer :: array(:)
  end type patch_list_type

  PetscInt, parameter, public :: INT_VAR = 0
  PetscInt, parameter, public :: REAL_VAR = 1

  interface PatchGetVariable
    module procedure PatchGetVariable1
    module procedure PatchGetVariable2
  end interface

  public :: PatchCreate, PatchDestroy, PatchCreateList, PatchDestroyList, &
            PatchAddToList, PatchConvertListToArray, PatchProcessCouplers, &
            PatchUpdateAllCouplerAuxVars, PatchInitAllCouplerAuxVars, &
            PatchLocalizeRegions, PatchUpdateUniformVelocity, &
            PatchGetVariable, PatchGetVariableValueAtCell, &
            PatchSetVariable, PatchCouplerInputRecord, &
            PatchInitConstraints, &
            PatchCountCells, PatchGetIvarsFromKeyword, &
            PatchGetVarNameFromKeyword, &
            PatchCalculateCFL1Timestep, &
            PatchGetCellCenteredVelocities, &
            PatchGetCompMassInRegion, &
            PatchGetWaterMassInRegion, &
            PatchGetCompMassInRegionAssign

contains

! ************************************************************************** !

function PatchCreate()
  !
  ! Allocates and initializes a new Patch object
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_type), pointer :: PatchCreate

  type(patch_type), pointer :: patch

  allocate(patch)

  patch%id = 0
  patch%surf_or_subsurf_flag = SUBSURFACE
  nullify(patch%imat)
  nullify(patch%imat_internal_to_external)
  nullify(patch%sat_func_id)
  nullify(patch%internal_velocities)
  nullify(patch%boundary_velocities)
  nullify(patch%internal_tran_coefs)
  nullify(patch%boundary_tran_coefs)
  nullify(patch%internal_flow_fluxes)
  nullify(patch%boundary_flow_fluxes)
  nullify(patch%internal_tran_fluxes)
  nullify(patch%boundary_tran_fluxes)
  nullify(patch%ss_flow_fluxes)
  nullify(patch%ss_tran_fluxes)
  nullify(patch%ss_flow_vol_fluxes)
  nullify(patch%boundary_energy_flux)

  nullify(patch%grid)

  allocate(patch%region_list)
  call RegionInitList(patch%region_list)

  allocate(patch%boundary_condition_list)
  call CouplerInitList(patch%boundary_condition_list)
  allocate(patch%initial_condition_list)
  call CouplerInitList(patch%initial_condition_list)
  allocate(patch%source_sink_list)
  call CouplerInitList(patch%source_sink_list)

  nullify(patch%material_properties)
  nullify(patch%material_property_array)

  nullify(patch%characteristic_curves)
  nullify(patch%characteristic_curves_array)

  allocate(patch%observation_list)
  call ObservationInitList(patch%observation_list)
  allocate(patch%integral_flux_list)
  call IntegralFluxInitList(patch%integral_flux_list)
  allocate(patch%strata_list)
  call StrataInitList(patch%strata_list)

  call AuxInit(patch%aux)

  nullify(patch%field)
  nullify(patch%reaction)
  nullify(patch%datasets)

  nullify(patch%next)

  nullify(patch%surf_material_properties)
  nullify(patch%surf_material_property_array)
  nullify(patch%surf_field)
  call SurfaceAuxInit(patch%surf_aux)

  PatchCreate => patch

end function PatchCreate

! ************************************************************************** !

function PatchCreateList()
  !
  ! PatchListCreate: Creates a patch list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_list_type), pointer :: PatchCreateList

  type(patch_list_type), pointer :: patch_list

  allocate(patch_list)
  nullify(patch_list%first)
  nullify(patch_list%last)
  nullify(patch_list%array)
  patch_list%num_patch_objects = 0

  PatchCreateList => patch_list

end function PatchCreateList

! ************************************************************************** !

subroutine PatchAddToList(new_patch,patch_list)
  !
  ! Adds a new patch to list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_type), pointer :: new_patch
  type(patch_list_type) :: patch_list

  if (associated(new_patch)) then
     patch_list%num_patch_objects = patch_list%num_patch_objects + 1
     new_patch%id = patch_list%num_patch_objects
     if (.not.associated(patch_list%first)) patch_list%first => new_patch
     if (associated(patch_list%last)) patch_list%last%next => new_patch
     patch_list%last => new_patch
  end if
end subroutine PatchAddToList

! ************************************************************************** !

subroutine PatchConvertListToArray(patch_list)
  !
  ! Creates an array of pointers to the
  ! patchs in the patch list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_list_type) :: patch_list

  PetscInt :: count
  type(patch_type), pointer :: cur_patch


  allocate(patch_list%array(patch_list%num_patch_objects))

  cur_patch => patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    patch_list%array(cur_patch%id)%ptr => cur_patch
    cur_patch => cur_patch%next
  enddo

end subroutine PatchConvertListToArray

! ************************************************************************** !

subroutine PatchLocalizeRegions(patch,regions,option)
  !
  ! Localizes regions within each patch
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Output_Aux_module
  use Region_module

  implicit none

  type(patch_type) :: patch
  type(region_list_type) :: regions
  type(option_type) :: option

  type(region_type), pointer :: cur_region
  type(region_type), pointer :: patch_region

  cur_region => regions%first
  do
    if (.not.associated(cur_region)) exit
    patch_region => RegionCreate(cur_region)
    call RegionAddToList(patch_region,patch%region_list)
    cur_region => cur_region%next
  enddo

  !geh: All grids must be localized through GridLocalizeRegions.  Patch
  !     should not differentiate between structured/unstructured, etc.
  call GridLocalizeRegions(patch%grid,patch%region_list,option)

end subroutine PatchLocalizeRegions

! ************************************************************************** !
subroutine PatchProcessCouplers(patch,flow_conditions,transport_conditions, &
                                option)

  !
  ! Assigns conditions and regions to couplers
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Material_module
  use Condition_module
  use Transport_Constraint_module
  use Connection_module

  implicit none

  type(patch_type) :: patch
  type(condition_list_type) :: flow_conditions
  type(tran_condition_list_type) :: transport_conditions
  type(option_type) :: option

  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation, next_observation
  type(integral_flux_type), pointer :: integral_flux

  PetscInt :: temp_int, isub
  PetscInt :: nphase
  PetscErrorCode :: ierr

  ! boundary conditions
  coupler => patch%boundary_condition_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in boundary condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call printErrMsg(option)
    endif
    if (associated(patch%grid%structured_grid)) then
      if (coupler%region%num_cells > 0 .and. &
          (coupler%region%iface == 0 .and. &
           .not.associated(coupler%region%faces))) then
        option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '", which is tied to a boundary condition, has not &
                 &been assigned a face in the structured grid. '
        call printErrMsg(option)
      endif
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in &
                           &BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name, &
                                      transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
           option%io_buffer = 'Transport condition "' // &
                   trim(coupler%tran_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in &
                           &BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => patch%initial_condition_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in initial condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call printErrMsg(option)
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in initial condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in ' // &
                           'INITIAL_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name, &
                                      transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
          option%io_buffer = 'Transport condition "' // &
                   trim(coupler%tran_condition_name) // &
                   '" in initial condition "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in &
                           &INITIAL_CONDITION: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in source/sink "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call printErrMsg(option)
    endif
   
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call printErrMsg(option)
        endif
        ! check to ensure that a rate subcondition exists
        if (.not.associated(coupler%flow_condition%rate)) then
          temp_int = 0
          if (temp_int == 0) then
            option%io_buffer = 'FLOW_CONDITIONs associated with &
              &SOURCE_SINKs must have a RATE or WELL expression within them.'
            call printErrMsg(option)
          endif
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in &
                           &SOURCE_SINK: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name, &
                                      transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
          option%io_buffer = 'Transport condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call printErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in &
                           &SOURCE_SINK: ' // trim(coupler%name) // '.'
        call printErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo

!----------------------------
! AUX

  ! strata
  ! connect pointers from strata to regions
  strata => patch%strata_list%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    if (len_trim(strata%region_name) > 0) then
      strata%region => RegionGetPtrFromList(strata%region_name, &
                                                  patch%region_list)
      if (.not.associated(strata%region)) then
        option%io_buffer = 'Region "' // trim(strata%region_name) // &
                 '" in strata not found in region list'
        call printErrMsg(option)
      endif
      if (strata%active) then
        ! pointer to material
        ! gb: Depending on a surface/subsurface patch, use corresponding
        !     material properties
        if (patch%surf_or_subsurf_flag == SUBSURFACE) then
          strata%material_property => &
            MaterialPropGetPtrFromArray(strata%material_property_name, &
                                        patch%material_property_array)
          if (.not.associated(strata%material_property)) then
            option%io_buffer = 'Material "' // &
                              trim(strata%material_property_name) // &
                              '" not found in material list'
            call printErrMsg(option)
          endif
        endif

        if (patch%surf_or_subsurf_flag == SURFACE) then
          strata%surf_material_property => &
            SurfaceMaterialPropGetPtrFromArray(strata%material_property_name, &
                                            patch%surf_material_property_array)
          if (.not.associated(strata%surf_material_property)) then
            option%io_buffer = 'Material "' // &
                              trim(strata%material_property_name) // &
                              '" not found in material list'
            call printErrMsg(option)
          endif
        endif

      endif
    else
      nullify(strata%region)
      nullify(strata%material_property)
    endif
    strata => strata%next
  enddo

  ! connectivity between initial conditions, boundary conditions,
  ! srcs/sinks, etc and grid
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%initial_condition_list)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%boundary_condition_list)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%source_sink_list)

  ! linkage of observation to regions and couplers must take place after
  ! connection list have been created.
  ! observation
  observation => patch%observation_list%first
  do
    if (.not.associated(observation)) exit
    next_observation => observation%next
    select case(observation%itype)
      case(OBSERVATION_SCALAR)
        ! pointer to region
        observation%region => RegionGetPtrFromList(observation%linkage_name, &
                                                    patch%region_list)
        if (.not.associated(observation%region)) then
          option%io_buffer = 'Region "' // &
                   trim(observation%linkage_name) // &
                 '" in observation point "' // &
                 trim(observation%name) // &
                 '" not found in region list'
          call printErrMsg(option)
        endif
        call MPI_Allreduce(observation%region%num_cells,temp_int, &
                           ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                           option%mycomm,ierr)
        if (temp_int == 0) then
          option%io_buffer = 'Region "' // trim(observation%region%name) // &
            '" is used in an observation point but lies outside the &
            &model domain.'
          call printErrMsg(option)
        endif
        if (observation%region%num_cells == 0) then
          ! remove the observation object
          call ObservationRemoveFromList(observation,patch%observation_list)
        endif
      case(OBSERVATION_FLUX)
        coupler => CouplerGetPtrFromList(observation%linkage_name, &
                                         patch%boundary_condition_list)
        if (associated(coupler)) then
          observation%connection_set => coupler%connection_set
        else
          option%io_buffer = 'Boundary Condition "' // &
                   trim(observation%linkage_name) // &
                   '" not found in Boundary Condition list'
          call printErrMsg(option)
        endif
        if (observation%connection_set%num_connections == 0) then
          ! cannot remove from list, since there must be a global reduction
          ! across all procs
          ! therefore, just nullify connection set
          nullify(observation%connection_set)
        endif
    end select
    observation => next_observation
  enddo

  ! linkage of observation to regions and couplers must take place after
  ! connection list have been created.
  ! observation
  integral_flux => patch%integral_flux_list%first
  do
    if (.not.associated(integral_flux)) exit
    integral_flux%connections => &
      PatchGetConnectionsFromCoords(patch,integral_flux%coordinates, &
                                    integral_flux%name,option)
    call IntegralFluxSizeStorage(integral_flux,option)
    integral_flux => integral_flux%next
    option%flow%store_fluxes = PETSC_TRUE
    option%transport%store_fluxes = PETSC_TRUE
  enddo

  temp_int = ConnectionGetNumberInList(patch%grid%internal_connection_set_list)
  temp_int = max(temp_int,1)
  nphase = max(option%nphase,option%transport%nphase)

  ! all simulations
  allocate(patch%internal_velocities(nphase,temp_int))
  patch%internal_velocities = 0.d0

  ! flow
  if (option%nflowdof > 0) then
    if (option%flow%store_fluxes .or. &
        (patch%surf_or_subsurf_flag == SURFACE)) then
      allocate(patch%internal_flow_fluxes(option%nflowdof,temp_int))
      patch%internal_flow_fluxes = 0.d0
    endif
  endif

  ! transport
  if (option%ntrandof > 0) then
    allocate(patch%internal_tran_coefs(option%ntrandof,nphase,temp_int))
    patch%internal_tran_coefs = 0.d0
    if (option%transport%store_fluxes) then
      allocate(patch%internal_tran_fluxes(option%ntrandof,temp_int))
      patch%internal_tran_fluxes = 0.d0
    endif
  endif

  temp_int = CouplerGetNumConnectionsInList(patch%boundary_condition_list)

  if (temp_int > 0) then
    ! all simulations
    allocate(patch%boundary_velocities(nphase,temp_int))
    patch%boundary_velocities = 0.d0
    ! flow
    if (option%nflowdof > 0) then
      if (option%flow%store_fluxes .or. &
          (patch%surf_or_subsurf_flag == SURFACE)) then
        allocate(patch%boundary_flow_fluxes(option%nflowdof,temp_int))
        patch%boundary_flow_fluxes = 0.d0
      endif
      ! surface/subsurface storage
      if (option%iflowmode == TH_MODE) then
        allocate(patch%boundary_energy_flux(2,temp_int))
        patch%boundary_energy_flux = 0.d0
      endif
    endif
    ! transport
    if (option%ntrandof > 0) then
      allocate(patch%boundary_tran_coefs(option%ntrandof,nphase, &
                                         temp_int))
      patch%boundary_tran_coefs = 0.d0
      if (option%transport%store_fluxes) then
        allocate(patch%boundary_tran_fluxes(option%ntrandof,temp_int))
        patch%boundary_tran_fluxes = 0.d0
      endif
    endif
  endif

  temp_int = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (temp_int > 0) then
    ! flow
    if (option%nflowdof > 0) then
      allocate(patch%ss_flow_fluxes(option%nflowdof,temp_int))
      patch%ss_flow_fluxes = 0.d0
      if (option%nwells > 0) then
        ! needed by wells
        allocate(patch%ss_flow_vol_fluxes(nphase,temp_int))
        patch%ss_flow_vol_fluxes = 0.d0
      end if
    endif
    ! transport
    if (option%ntrandof > 0) then
      allocate(patch%ss_tran_fluxes(option%ntrandof,temp_int))
      patch%ss_tran_fluxes = 0.d0
      ! only needed by transport
      allocate(patch%ss_flow_vol_fluxes(nphase,temp_int))
      patch%ss_flow_vol_fluxes = 0.d0
    endif
  endif

#ifdef WELL_CLASS
  !create well communicators - if no wells exit without any operation
  call PatchCreateWellComms(patch,option)
#endif

end subroutine PatchProcessCouplers

! ************************************************************************** !
#ifdef WELL_CLASS
subroutine PatchCreateWellComms(patch,option)
  !
  ! Create well groups and communicators
  ! Wells are defined as source_sink couplers
  !
  ! Author: Paolo Orsini - OpenGoSim
  ! Date: 6/10/2015
  !

  use Option_module
  !use Connection_module ! do I need this??
  !use Well_Base_class

  implicit none


  type(patch_type) :: patch
  type(option_type) :: option

  PetscInt :: num_couplers,num_wells, comm_size
  PetscInt :: i_well, well_idx, i_rank
  PetscInt, pointer :: rnk_idx(:)
  PetscInt :: temp_int
  PetscInt :: recvbuf
  PetscErrorCode :: ierr
  PetscInt, pointer :: snd_well_ranks(:), rcv_well_ranks(:)
  PetscInt, pointer :: clp_to_well(:)
  PetscInt :: nw
  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list

  type wells_loc_type
    PetscInt :: cpl_id
    PetscInt :: num_ranks
    PetscMPIInt :: comm
    PetscMPIInt :: group
    PetscInt, pointer :: ranks(:)
  end type wells_loc_type

  type(wells_loc_type), pointer :: wells_loc(:)

  ! to be redesigned so that:
  ! number of well spec can be repeated to define more wells:
  ! WELL_SPEC read from the input deck can
  ! be used to create more coupler_wells (new type of coupler or extended type)

  num_couplers = patch%source_sink_list%num_couplers
  allocate(clp_to_well(num_couplers))
  clp_to_well = -1

  num_wells = 0
  ! count the number of wells and create couplers to wells map
  coupler => patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    if ( associated(coupler%well) ) then
      num_wells = num_wells + 1
      clp_to_well(coupler%id) = num_wells
    end if
    coupler => coupler%next
  enddo

  if (num_wells == 0) return

  allocate(snd_well_ranks(num_wells))
  ! flag with -1 processors/ranks that do not share a well
  snd_well_ranks = -1
  comm_size=option%mycommsize
  allocate(rcv_well_ranks(num_wells*comm_size))

  allocate(wells_loc(num_wells))
  do i_well=1,num_wells
    wells_loc(i_well)%cpl_id=-9
    wells_loc(i_well)%num_ranks = 0
    wells_loc(i_well)%comm = 0
    wells_loc(i_well)%group = 0
    nullify(wells_loc(i_well)%ranks)
  end do

  ! detecting and saving ranks containing each well
  coupler => patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    if ( associated(coupler%well) ) then
      ! empty well coupler (num_connections=0) can exist because of
      ! global regions that do not exist in the current processor/sub-domain
      ! only non-empty well are added to the well group/communicator
      if (coupler%connection_set%num_connections > 0) then
        nw = clp_to_well(coupler%id)
        snd_well_ranks(nw) = option%myrank
        wells_loc(nw)%cpl_id = coupler%id ! well to coupler map
      end if
    end if
    coupler => coupler%next
  enddo

#ifdef WELL_DEBUG
  print *, "rank= ", option%myrank
  print *, "before MPI_Allgather"
  do i_well=1,num_wells
    print *, "i_well",i_well,"snd=",snd_well_ranks(i_well),"  rank= ",option%myrank
  end do
  !print *, "snd1= ", snd_well_ranks(1),"  snd2= ",snd_well_ranks(2),"  rank= ",option%myrank
#endif
  ! if in the current processor/rank there are no wells snd_well_ranks = -1
  call MPI_Allgather(snd_well_ranks, num_wells, MPI_INTEGER, rcv_well_ranks, &
                     num_wells, MPI_INTEGER, option%mycomm, ierr);
#ifdef WELL_DEBUG
  well_idx = 1
  do i_rank=1,comm_size
    print *, "rcv rank segment = ", i_rank,"  rank= ",option%myrank
    do i_well=1,num_wells
      print *, "rcv_idx= ", well_idx, "rcv= ", rcv_well_ranks(well_idx),&
               "  rank= ",option%myrank
      well_idx = well_idx + 1
    end do
  end do
 ! print *, "rcv1= ", rcv_well_ranks(1),"  rcv2= ",rcv_well_ranks(2),"  rank= ",option%myrank
 ! print *, "rcv3= ", rcv_well_ranks(3),"  rcv4= ",rcv_well_ranks(4),"  rank= ",option%myrank
 ! print *, "rcv5= ", rcv_well_ranks(5),"  rcv6= ",rcv_well_ranks(6),"  rank= ",option%myrank
 ! print *, "rcv7= ", rcv_well_ranks(7),"  rcv8= ",rcv_well_ranks(8),"  rank= ",option%myrank
  print *, "after MPI_Allgather"
#endif

  ! First, it defines how many processors/ranks share each well,
  ! then, it loops again to allocate and load the ranks.
  ! This could be done with a list to avoid the double loop,
  ! but this is not an expensive operation ( max_num_wells= 100-200)
  well_idx=1
  do i_rank=1,comm_size
    do i_well=1,num_wells
      if ( rcv_well_ranks(well_idx) /= -1) then
        wells_loc(i_well)%num_ranks = wells_loc(i_well)%num_ranks + 1
      end if
      well_idx = well_idx + 1
    end do
  end do

  do i_well=1,num_wells
    if (wells_loc(i_well)%num_ranks > 0) then
      allocate(wells_loc(i_well)%ranks(wells_loc(i_well)%num_ranks))
      wells_loc(i_well)%ranks=-11
    end if
  end do

  allocate(rnk_idx(num_wells))
  rnk_idx=1
  do i_well=1,num_wells
    do i_rank=1,comm_size
      well_idx = i_well + (i_rank-1)*num_wells
      if (rcv_well_ranks(well_idx) /= -1 .and. rcv_well_ranks(well_idx)>=0) then
        wells_loc(i_well)%ranks(rnk_idx(i_well)) = rcv_well_ranks(well_idx)
        rnk_idx(i_well) = rnk_idx(i_well) + 1
      else if (rcv_well_ranks(well_idx)/=-1.and.rcv_well_ranks(well_idx)<0) then
        option%io_buffer = 'WELL processing: well ranks must be greater than 0'
        call printErrMsg(option)
      end if
    end do
  end do
  deallocate(rnk_idx)
  nullify(rnk_idx)

#ifdef WELL_DEBUG
  print *,"num_wells=", num_wells, "rank= ",option%myrank
  do i_well=1,num_wells
    print *,"well id= ",i_well,"num rank=",wells_loc(i_well)%num_ranks, &
             "rank= ",option%myrank
    do i_rank=1,wells_loc(i_well)%num_ranks
      print *,"well id= ",i_well,"i_rank=", wells_loc(i_well)%ranks(i_rank), &
            "rank= ",option%myrank
    end do
   ! print *,"well id= ",i_well,"ranks=", wells_loc(i_well)%ranks(1), &
   !                      wells_loc(i_well)%ranks(2), &
   !                     "rank= ",option%myrank
  end do
#endif

  ! create well groups and communicators for in each process
  ! each well must have at least one rank associated
  ! each well must have at least one group and one comm
  ! all processes call MPI_GROUP_INCL and MPI_COMM_CREATE
  do i_well=1,num_wells
    call MPI_GROUP_INCL(option%mygroup,wells_loc(i_well)%num_ranks,&
                        wells_loc(i_well)%ranks,wells_loc(i_well)%group,ierr)
    call MPI_COMM_CREATE(option%mycomm, wells_loc(i_well)%group, &
                         wells_loc(i_well)%comm, ierr)
  end do

  ! assign well group and communicators to well-couplers
  coupler => patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    if ( associated(coupler%well) ) then
#ifdef WELL_DEBUG
      print *, "well= ", coupler%id, "  num_conn= ", &
                coupler%connection_set%num_connections,"rank= ", option%myrank
      print *, coupler%id, wells_loc(coupler%id)%group
      print *, coupler%id, wells_loc(coupler%id)%comm
      print *, coupler%id, coupler%well%group
      print *, coupler%id, coupler%well%comm
#endif
      nw = clp_to_well(coupler%id)
      coupler%well%group = wells_loc(nw)%group
      coupler%well%comm = wells_loc(nw)%comm
    end if
    coupler => coupler%next
  end do

#ifdef WELL_DEBUG
  ! testing
  coupler => patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    !if (associated(coupler%well)) then
    ! only the non-empty wells in the rank can use the well communicators
    ! only ranks with a non-empty portion of the well are included in the
    ! well communicators. Only ranks belonging to the well can use the well comm
    if ( associated(coupler%well) ) then
      if (coupler%connection_set%num_connections > 0) then
        if (coupler%id==1) then
          print *, "before MPI_ALLREDUCE, well_id= ",coupler%id, "rank= ", option%myrank
          call MPI_ALLREDUCE(1, recvbuf, 1, MPI_INTEGER,MPI_SUM, &
                           coupler%well%comm, ierr)
          print *, "result reduce well 1 =", recvbuf, "rank= ", option%myrank
        end if
        if (coupler%id==2) then
          print *, "before MPI_ALLREDUCE, well_id= ",coupler%id, "rank= ", option%myrank
          call MPI_ALLREDUCE(2, recvbuf, 1, MPI_INTEGER,MPI_SUM, &
                             coupler%well%comm, ierr)
          print *, "result reduce well 2 =", recvbuf, "rank= ", option%myrank
        end if
      end if
    end if
    coupler => coupler%next
  end do
#endif

 ! deallocate local variables
  deallocate(clp_to_well)
  nullify(clp_to_well)
  deallocate(snd_well_ranks)
  nullify(snd_well_ranks)
  deallocate(rcv_well_ranks)
  nullify(rcv_well_ranks)

  do i_well=1,size(wells_loc(:))
    if(associated(wells_loc(i_well)%ranks)) then
      deallocate(wells_loc(i_well)%ranks)
      nullify(wells_loc(i_well)%ranks)
    end if
  end do
  deallocate(wells_loc)
  nullify(wells_loc)

  nullify(coupler)

end subroutine PatchCreateWellComms
#endif
! ************************************************************************** !

subroutine PatchInitAllCouplerAuxVars(patch,option)
  !
  ! Initializes coupler auxillary variables
  ! within list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Reaction_Aux_module

  implicit none

  type(patch_type), pointer :: patch
  type(option_type) :: option

  PetscBool :: force_update_flag = PETSC_TRUE

  call PatchInitCouplerAuxVars(patch%initial_condition_list,patch, &
                               option)
  call PatchInitCouplerAuxVars(patch%boundary_condition_list,patch, &
                               option)
  call PatchInitCouplerAuxVars(patch%source_sink_list,patch, &
                               option)

  !geh: This should not be included in PatchUpdateAllCouplerAuxVars
  ! as it will result in excessive updates to initial conditions
  ! that are not necessary after the simulation has started time stepping.
  call PatchUpdateCouplerAuxVars(patch,patch%initial_condition_list, &
                                 force_update_flag,option)
  call PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

end subroutine PatchInitAllCouplerAuxVars

! ************************************************************************** !

subroutine PatchInitCouplerAuxVars(coupler_list,patch,option)
  !
  ! Initializes coupler auxillary variables
  ! within list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Connection_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Condition_module
  use Transport_Constraint_module


  implicit none

  type(coupler_list_type), pointer :: coupler_list
  type(patch_type), pointer :: patch
  type(option_type) :: option

  PetscInt :: num_connections
  PetscBool :: force_update_flag

  type(coupler_type), pointer :: coupler
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  PetscInt :: idof
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: temp_int

  if (.not.associated(coupler_list)) return

  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit

    if (associated(coupler%connection_set)) then
      num_connections = coupler%connection_set%num_connections

      ! FLOW
      if (associated(coupler%flow_condition)) then
        ! determine whether flow_condition is transient
        coupler%flow_condition%is_transient = &
          FlowConditionIsTransient(coupler%flow_condition)
        if (coupler%itype == INITIAL_COUPLER_TYPE .or. &
            coupler%itype == BOUNDARY_COUPLER_TYPE) then

          if (associated(coupler%flow_condition%pressure) .or. &
              associated(coupler%flow_condition%concentration) .or. &
              associated(coupler%flow_condition%saturation) .or. &
              associated(coupler%flow_condition%temperature) &
              ) then

            ! allocate arrays that match the number of connections
            select case(option%iflowmode)

              case(TH_MODE)
                temp_int = 2
                select case(coupler%flow_condition%pressure%itype)
                  case(CONDUCTANCE_BC,HET_CONDUCTANCE_BC)
                    temp_int = temp_int + 1
                end select
                allocate(coupler%flow_aux_real_var(temp_int,num_connections))
                allocate(coupler%flow_aux_int_var(1,num_connections))
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case default
            end select

          else if (associated(coupler%flow_condition%rate)) then
            option%io_buffer = 'Flow condition "' // &
              trim(coupler%flow_condition%name) // '" can only be used in a &
              &SOURCE_SINK since a rate is prescribed.'
            call printErrMsg(option)
          endif ! associated(coupler%flow_condition%pressure)

        else if (coupler%itype == SRC_SINK_COUPLER_TYPE) then

          if (associated(coupler%flow_condition%rate)) then

            select case(coupler%flow_condition%rate%itype)
              case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS, &
                   VOLUMETRIC_RATE_SS,MASS_RATE_SS, &
                   HET_VOL_RATE_SS,HET_MASS_RATE_SS)
                select case(option%iflowmode)
                  case(TH_MODE)
                    allocate(coupler%flow_aux_real_var(option%nflowdof,num_connections))
                    coupler%flow_aux_real_var = 0.d0
                  case default
                    string = GetSubConditionName(coupler%flow_condition%rate%itype)
                    option%io_buffer='Source/Sink of rate%itype = "' // &
                      trim(adjustl(string)) // '", not implemented in this mode.'
                    call printErrMsg(option)
                end select
              case default
                string = GetSubConditionName(coupler%flow_condition%rate%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition,'rate', &
                                            string)
                call printErrMsg(option)
            end select
          endif ! associated(coupler%flow_condition%rate)
        endif ! coupler%itype == SRC_SINK_COUPLER_TYPE
      endif ! associated(coupler%flow_condition)
    endif ! associated(coupler%connection_set)

    ! TRANSPORT
    if (associated(coupler%tran_condition)) then
      cur_constraint_coupler => &
        coupler%tran_condition%constraint_coupler_list
      do
        if (.not.associated(cur_constraint_coupler)) exit
        ! Setting option%iflag = 0 ensures that the "mass_balance" array
        ! is not allocated.
        option%iflag = 0
        ! Only allocate the XXX_auxvar objects if they have not been allocated.
        ! Since coupler%tran_condition is a pointer to a separate list of
        ! tran conditions, the XXX_auxvar object may already be allocated.
        if (.not.associated(cur_constraint_coupler%global_auxvar)) then
          allocate(cur_constraint_coupler%global_auxvar)
          call GlobalAuxVarInit(cur_constraint_coupler%global_auxvar,option)
        endif
        if (.not.associated(cur_constraint_coupler%rt_auxvar)) then
          allocate(cur_constraint_coupler%rt_auxvar)
          call RTAuxVarInit(cur_constraint_coupler%rt_auxvar,patch%reaction, &
                            option)
        endif
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif

    coupler => coupler%next
  enddo

end subroutine PatchInitCouplerAuxVars

! ************************************************************************** !

subroutine PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)
  !
  ! Updates auxiliary variables associated
  ! with couplers in list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module

  implicit none

  type(patch_type) :: patch
  PetscBool :: force_update_flag
  type(option_type) :: option

  PetscInt :: iconn

  !geh: no need to update initial conditions as they only need updating
  !     once as performed in PatchInitCouplerAuxVars()
  call PatchUpdateCouplerAuxVars(patch,patch%boundary_condition_list, &
                                 force_update_flag,option)
  call PatchUpdateCouplerAuxVars(patch,patch%source_sink_list, &
                                 force_update_flag,option)

!  stop
end subroutine PatchUpdateAllCouplerAuxVars

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVars(patch,coupler_list,force_update_flag, &
                                     option)
  !
  ! Updates auxiliary variables associated
  ! with couplers in list
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !
  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module


  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class

  implicit none

  type(patch_type) :: patch
  type(coupler_list_type), pointer :: coupler_list
  PetscBool :: force_update_flag
  type(option_type) :: option

  type(coupler_type), pointer :: coupler
  type(flow_condition_type), pointer :: flow_condition

  if (.not.associated(coupler_list)) return

  coupler => coupler_list%first

  do
    if (.not.associated(coupler)) exit

    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then

      flow_condition => coupler%flow_condition
      if (force_update_flag .or. flow_condition%is_transient) then
        select case(option%iflowmode)
          case(TH_MODE)
            call PatchUpdateCouplerAuxVarsTH(patch,coupler,option)
        end select
      endif
    endif

    ! TRANSPORT
    ! nothing for transport at this point in time
    coupler => coupler%next
  enddo

end subroutine PatchUpdateCouplerAuxVars

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsTH(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for TH_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Utility_module, only : DeallocateArray

  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  PetscBool :: apply_temp_cond
  PetscInt :: rate_scale_type

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: iphase

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    !geh: this is a fix for an Intel compiler bug. Not sure why Intel cannot
    !     access flow_condition%iphase directly....
    iphase = flow_condition%iphase
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = iphase
!    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
!                                                        flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC,SPILLOVER_BC)
        select type(selector =>flow_condition%pressure%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(TH_PRESSURE_DOF,1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerFromDataset(coupler,option, &
                                               patch%grid,selector, &
                                               TH_PRESSURE_DOF)
          class default
            option%io_buffer = 'Unknown dataset class (TH%' // &
              'pressure%itype,DIRICHLET_BC)'
            call printErrMsg(option)
        end select
      case(HYDROSTATIC_BC,SEEPAGE_BC,CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
      case(HET_DIRICHLET_BC,HET_SEEPAGE_BC,HET_CONDUCTANCE_BC)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%pressure%dataset,TH_PRESSURE_DOF,option)
        if (flow_condition%pressure%itype == HET_CONDUCTANCE_BC) then
          coupler%flow_aux_real_var(TH_CONDUCTANCE_DOF,1:num_connections) = &
            flow_condition%pressure%aux_real(1)
        endif
      case(HET_SURF_SEEPAGE_BC)
        ! Do nothing, since this BC type is only used for coupling of
        ! surface-subsurface model
      case default
        string = &
          GetSubConditionName(flow_condition%pressure%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH pressure',string)
        call printErrMsg(option)
    end select
    if (associated(flow_condition%temperature)) then
      select case(flow_condition%temperature%itype)
        case(DIRICHLET_BC,ZERO_GRADIENT_BC)
          select type(selector =>flow_condition%temperature%dataset)
            class is(dataset_ascii_type)
              if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
                 (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
                 flow_condition%temperature%itype /= DIRICHLET_BC)) then
                coupler%flow_aux_real_var(TH_TEMPERATURE_DOF, &
                                          1:num_connections) = &
                  selector%rarray(1)
              endif
            class is(dataset_gridded_hdf5_type)
              call PatchUpdateCouplerFromDataset(coupler,option, &
                                                 patch%grid,selector, &
                                                 TH_TEMPERATURE_DOF)
            class default
              option%io_buffer = 'Unknown dataset class (TH%' // &
                'temperature%itype,DIRICHLET_BC)'
              call printErrMsg(option)
          end select
        case (HET_DIRICHLET_BC)
          call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                  flow_condition%temperature%dataset, &
                  TH_TEMPERATURE_DOF,option)
        case default
          string = &
            GetSubConditionName(flow_condition%temperature%itype)
          option%io_buffer = &
            FlowConditionUnknownItype(flow_condition,'TH temperature',string)
          call printErrMsg(option)
      end select
    endif
    if (associated(flow_condition%energy_flux)) then
      coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,1:num_connections) = &
        flow_condition%energy_flux%dataset%rarray(1)
    endif
  endif

  apply_temp_cond = PETSC_FALSE
  if (associated(flow_condition%temperature) .and. associated(flow_condition%pressure)) then
    if (flow_condition%pressure%itype /= HYDROSTATIC_BC) then
      apply_temp_cond = PETSC_TRUE
    else
      if (flow_condition%temperature%itype /= DIRICHLET_BC) then
        apply_temp_cond = PETSC_TRUE
      endif
    endif
  else
    apply_temp_cond = PETSC_TRUE
  endif

  if (associated(flow_condition%temperature) .and. apply_temp_cond) then
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        select type(selector =>flow_condition%temperature%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(TH_TEMPERATURE_DOF, &
                                      1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerFromDataset(coupler,option, &
                                               patch%grid,selector, &
                                               TH_TEMPERATURE_DOF)
          class default
            option%io_buffer = 'Unknown dataset class (TH%' // &
              'pressure%itype,DIRICHLET_BC)'
            call printErrMsg(option)
        end select
      case (HET_DIRICHLET_BC)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%temperature%dataset, &
                TH_TEMPERATURE_DOF,option)
      case default
        string = &
          GetSubConditionName(flow_condition%temperature%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH temperature',string)
        call printErrMsg(option)
    end select
  endif

  if (associated(flow_condition%energy_flux)) then
    select case(flow_condition%energy_flux%itype)
      case(NEUMANN_BC)
        select type(selector =>flow_condition%energy_flux%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(TH_TEMPERATURE_DOF, &
                                      1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerFromDataset(coupler,option, &
                                               patch%grid,selector, &
                                               TH_TEMPERATURE_DOF)
          class default
            option%io_buffer = 'Unknown dataset class (TH%' // &
              'pressure%itype,NEUMANN_BC)'
            call printErrMsg(option)
        end select
      case default
        string = &
          GetSubConditionName(flow_condition%energy_flux%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH energy flux',string)
        call printErrMsg(option)
    end select
  endif

  !geh: we set this flag to ensure that we are not scaling mass and energy
  !     differently
  rate_scale_type = 0
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case (HET_MASS_RATE_SS,HET_VOL_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                                            flow_condition%rate%dataset, &
                                            TH_PRESSURE_DOF,option)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,flow_condition%rate%isubtype, &
                                  option)
        rate_scale_type = flow_condition%rate%isubtype
      case(MASS_RATE_SS,VOLUMETRIC_RATE_SS)
      ! do nothing here
      case default
        string = &
          GetSubConditionName(flow_condition%rate%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH rate',string)
        call printErrMsg(option)
    end select
  endif
  if (associated(flow_condition%energy_rate)) then
    select case (flow_condition%energy_rate%itype)
      case (ENERGY_RATE_SS)
        !geh: this is pointless as %dataset%rarray(1) is reference in TH,
        !     not the flow_aux_real_var!
        coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%energy_rate%dataset%rarray(1)
      case (SCALED_ENERGY_RATE_SS)
        if (rate_scale_type == 0) then
          call PatchScaleSourceSink(patch,coupler, &
                                    flow_condition%energy_rate%isubtype,option)
        else if (rate_scale_type == flow_condition%energy_rate%isubtype) then
          !geh: do nothing as it is taken care of later.
        else
          option%io_buffer = 'MASS and ENERGY scaling mismatch in ' // &
            'FLOW_CONDITION "' // trim(flow_condition%name) // '".'
          call printErrMsg(option)
        endif
        !geh: do nothing as the
      case (HET_ENERGY_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%energy_rate%dataset, &
                TH_TEMPERATURE_DOF,option)
      case default
        string = &
          GetSubConditionName(flow_condition%energy_rate%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH energy rate',string)
        call printErrMsg(option)
    end select
  endif

end subroutine PatchUpdateCouplerAuxVarsTH

! ************************************************************************** !

! ************************************************************************** !

subroutine PatchUpdateCouplerFromDataset(coupler,option,grid,dataset,dof)
  !
  ! Updates auxiliary variables from dataset.
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Grid_module
  use Coupler_module
  use Dataset_Gridded_HDF5_class

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  class(dataset_gridded_hdf5_type) :: dataset
  PetscInt :: dof

  PetscReal :: temp_real
  PetscInt :: iconn
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: x
  PetscReal :: y
  PetscReal :: z
  PetscReal :: dist(-1:3)

  do iconn = 1, coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    x = grid%x(ghosted_id)
    y = grid%y(ghosted_id)
    z = grid%z(ghosted_id)
    if (associated(coupler%connection_set%dist)) then
      dist = coupler%connection_set%dist(:,iconn)
      x = x-dist(0)*dist(1)
      y = y-dist(0)*dist(2)
      z = z-dist(0)*dist(3)
    endif
    call DatasetGriddedHDF5InterpolateReal(dataset,x,y,z,temp_real,option)
    coupler%flow_aux_real_var(dof,iconn) = temp_real
  enddo

end subroutine PatchUpdateCouplerFromDataset

! ************************************************************************** !

subroutine PatchScaleSourceSink(patch,source_sink,iscale_type,option)
  !
  ! Scales select source/sinks based on perms*volume
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  !

#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Material_Aux_class
  use Variables_module, only : PERMEABILITY_X

  implicit none

  type(patch_type) :: patch
  type(coupler_type) :: source_sink
  PetscInt :: iscale_type
  type(option_type) :: option

  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(connection_set_type), pointer :: cur_connection_set
  type(field_type), pointer :: field

  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount, x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  PetscInt :: ghosted_neighbors(27)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  field => patch%field
  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  grid => patch%grid

  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

  cur_connection_set => source_sink%connection_set

  select case(iscale_type)
    case(SCALE_BY_VOLUME)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        vec_ptr(local_id) = vec_ptr(local_id) + &
          material_auxvars(ghosted_id)%volume
      enddo
    case(SCALE_BY_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        vec_ptr(local_id) = vec_ptr(local_id) + &
          ! this function protects from error in gfortran compiler when indexing
          ! the permeability array
          MaterialAuxVarGetValue(material_auxvars(ghosted_id), &
                                 PERMEABILITY_X) * &
          material_auxvars(ghosted_id)%volume
      enddo
    case(SCALE_BY_NEIGHBOR_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        !geh: kludge for 64-bit integers.
        call GridGetGhostedNeighbors(grid,ghosted_id,DMDA_STENCIL_STAR, &
                                    x_width,y_width,z_width, &
                                    x_count,y_count,z_count, &
                                    ghosted_neighbors,option)
        ! ghosted neighbors is ordered first in x, then, y, then z
        icount = 0
        sum = 0.d0
        ! x-direction
        do while (icount < x_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_X) * &
                grid%structured_grid%dy(neighbor_ghosted_id)* &
                grid%structured_grid%dz(neighbor_ghosted_id)
        enddo
        ! y-direction
        do while (icount < x_count + y_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_X) * &
                grid%structured_grid%dx(neighbor_ghosted_id)* &
                grid%structured_grid%dz(neighbor_ghosted_id)
        enddo
        ! z-direction
        do while (icount < x_count + y_count + z_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_X) * &
                grid%structured_grid%dx(neighbor_ghosted_id)* &
                grid%structured_grid%dy(neighbor_ghosted_id)
        enddo
        vec_ptr(local_id) = vec_ptr(local_id) + sum
      enddo
    case(0)
      option%io_buffer = 'Unknown scaling type in PatchScaleSourceSink ' // &
        'for FLOW_CONDITION "' // trim(source_sink%flow_condition%name) // '".'
      call printErrMsg(option)
  end select

  call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
  call VecNorm(field%work,NORM_1,scale,ierr);CHKERRQ(ierr)
  if (scale < 1.d-40) then
    option%io_buffer = 'Zero infinity norm in PatchScaleSourceSink for ' // &
      'FLOW_CONDITION "' // trim(source_sink%flow_condition%name) // '".'
    call printErrMsg(option)
  endif
  scale = 1.d0/scale
  call VecScale(field%work,scale,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%work,vec_ptr, ierr);CHKERRQ(ierr)
  do iconn = 1, cur_connection_set%num_connections
    local_id = cur_connection_set%id_dn(iconn)
    select case(option%iflowmode)
      case(TH_MODE)
        source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
          vec_ptr(local_id)
    end select 
  enddo
  call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PatchScaleSourceSink

! ************************************************************************** !

subroutine PatchUpdateHetroCouplerAuxVars(patch,coupler,dataset_base, &
                                          isub_condition,option)
  !
  ! This subroutine updates aux vars for distributed copuler_type
  !
  ! Author: Gautam Bisht, LBL
  ! Date: 10/03/2012
  !

#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Dataset_module
  use Dataset_Map_HDF5_class
  use Dataset_Base_class
  use Dataset_Ascii_class

  implicit none

  type(patch_type) :: patch
  type(coupler_type) :: coupler
  class(dataset_base_type), pointer :: dataset_base
  PetscInt :: isub_condition
  type(option_type) :: option

  type(connection_set_type), pointer :: cur_connection_set
  type(grid_type),pointer :: grid
  PetscErrorCode :: ierr
  PetscInt :: iconn
  PetscInt :: ghosted_id,local_id
  PetscInt,pointer ::cell_ids_nat(:)
  type(flow_sub_condition_type) :: flow_sub_condition

  class(dataset_map_hdf5_type), pointer :: dataset_map_hdf5
  class(dataset_ascii_type), pointer :: dataset_ascii

  grid => patch%grid

  if (isub_condition>option%nflowdof*option%nphase) then
    option%io_buffer='ERROR: PatchUpdateHetroCouplerAuxVars  '// &
      'isub_condition > option%nflowdof*option%nphase.'
    call printErrMsg(option)
  endif

  if (option%iflowmode/=TH_MODE) then
    option%io_buffer='PatchUpdateHetroCouplerAuxVars only implemented '// &
      ' for TH mode.'
    call printErrMsg(option)
  endif

  cur_connection_set => coupler%connection_set

  select type(selector=>dataset_base)
    class is(dataset_map_hdf5_type)
      dataset_map_hdf5 => selector

      ! If called for the first time, create the map
      if (dataset_map_hdf5%first_time) then
        allocate(cell_ids_nat(cur_connection_set%num_connections))
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          cell_ids_nat(iconn)=grid%nG2A(ghosted_id)
        enddo

        call PatchCreateFlowConditionDatasetMap(patch%grid,dataset_map_hdf5,&
                cell_ids_nat,cur_connection_set%num_connections,option)

        dataset_map_hdf5%first_time = PETSC_FALSE
        deallocate(cell_ids_nat)

      endif

      ! Save the data in the array
      do iconn=1,cur_connection_set%num_connections
        coupler%flow_aux_real_var(isub_condition,iconn) = &
          dataset_map_hdf5%rarray(dataset_map_hdf5%datatocell_ids(iconn))
      enddo

    class is(dataset_ascii_type)
      dataset_ascii => selector

      do iconn=1,cur_connection_set%num_connections
        coupler%flow_aux_real_var(isub_condition,iconn) = &
          dataset_ascii%rarray(1)
      enddo

    class default
      option%io_buffer = 'Incorrect dataset class (' // &
        trim(DatasetGetClass(dataset_base)) // &
        ') for coupler "' // trim(coupler%name) // &
        '" in PatchUpdateHetroCouplerAuxVars.'
      call printErrMsg(option)
  end select

end subroutine PatchUpdateHetroCouplerAuxVars

! ************************************************************************** !

subroutine PatchCreateFlowConditionDatasetMap(grid,dataset_map_hdf5,cell_ids,ncells,option)
  !
  ! This routine creates dataset-map for flow condition
  !
  ! Author: Gautam Bisht, LBL
  ! Date: 10/26/12
  !
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Dataset_Map_HDF5_class
  use Option_module

  implicit none

  type(grid_type) :: grid
  class(dataset_map_hdf5_type) :: dataset_map_hdf5
  type(option_type) :: option
  PetscInt,pointer :: cell_ids(:)
  PetscInt :: ncells

  PetscInt, allocatable :: int_array(:)
  PetscInt :: ghosted_id,local_id
  PetscInt :: ii,count
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  PetscInt :: nloc,nglo
  PetscInt :: istart

  IS :: is_from, is_to
  Vec :: map_ids_1, map_ids_2,map_ids_3
  VecScatter ::vec_scatter
  PetscViewer :: viewer

  ! Step-1: Rearrange map dataset
  nloc = maxval(dataset_map_hdf5%mapping(2,:))
  call MPI_Allreduce(nloc,nglo,ONE_INTEGER,MPIU_INTEGER,MPI_Max, &
                     option%mycomm,ierr)
  call VecCreateMPI(option%mycomm,dataset_map_hdf5%map_dims_local(2),&
                    PETSC_DETERMINE,map_ids_1,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE,nglo,map_ids_2, &
                    ierr);CHKERRQ(ierr)
  call VecSet(map_ids_2,0.d0,ierr);CHKERRQ(ierr)

  istart = 0
  call MPI_Exscan(dataset_map_hdf5%map_dims_local(2), istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  allocate(int_array(dataset_map_hdf5%map_dims_local(2)))
  do ii=1,dataset_map_hdf5%map_dims_local(2)
    int_array(ii)=ii+istart
  enddo
  int_array=int_array-1

  call ISCreateBlock(option%mycomm,1,dataset_map_hdf5%map_dims_local(2), &
                     int_array,PETSC_COPY_VALUES,is_from,ierr);CHKERRQ(ierr)
  deallocate(int_array)

  allocate(int_array(dataset_map_hdf5%map_dims_local(2)))
  do ii=1,dataset_map_hdf5%map_dims_local(2)
    int_array(ii)=dataset_map_hdf5%mapping(2,ii)
  enddo
  int_array=int_array-1

  call ISCreateBlock(option%mycomm,1,dataset_map_hdf5%map_dims_local(2), &
                     int_array,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)
  deallocate(int_array)

  !call VecCreateSeq(PETSC_COMM_SELF,dataset_map%map_dims_global(2),map_ids_1,ierr)
  !call VecCreateSeq(PETSC_COMM_SELF,maxval(dataset_map%map(2,:)),map_ids_2,ierr)
  !call VecSet(map_ids_2,0,ierr)

  call VecScatterCreate(map_ids_1,is_from,map_ids_2,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(map_ids_1,vec_ptr,ierr);CHKERRQ(ierr)
  do ii=1,dataset_map_hdf5%map_dims_local(2)
    vec_ptr(ii)=dataset_map_hdf5%mapping(1,ii)
  enddo
  call VecRestoreArrayF90(map_ids_1,vec_ptr,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,map_ids_1,map_ids_2, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,map_ids_1,map_ids_2, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

  ! Step-2: Get ids in map dataset for cells
  allocate(int_array(ncells))
  allocate(dataset_map_hdf5%cell_ids_local(ncells))
  int_array=cell_ids-1

  call ISCreateBlock(option%mycomm,1,ncells,int_array,PETSC_COPY_VALUES,is_from, &
                     ierr);CHKERRQ(ierr)

  istart = 0
  call MPI_Exscan(ncells, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  do local_id=1,ncells
    int_array(local_id)=local_id+istart
  enddo
  int_array=int_array-1

  call ISCreateBlock(option%mycomm,1,ncells,int_array,PETSC_COPY_VALUES,is_to, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)

  !call VecCreateSeq(PETSC_COMM_SELF,ncells,map_ids_3,ierr)
  call VecCreateMPI(option%mycomm,ncells,PETSC_DETERMINE,map_ids_3, &
                    ierr);CHKERRQ(ierr)

  call VecScatterCreate(map_ids_2,is_from,map_ids_3,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,map_ids_2,map_ids_3, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,map_ids_2,map_ids_3, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

  ! Step-3: Save the datatocell_ids
  allocate(dataset_map_hdf5%datatocell_ids(ncells))
  call VecGetArrayF90(map_ids_3,vec_ptr,ierr);CHKERRQ(ierr)
  do local_id=1,ncells
    dataset_map_hdf5%datatocell_ids(local_id) = int(vec_ptr(local_id))
  enddo
  call VecRestoreArrayF90(map_ids_3,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(map_ids_1,ierr);CHKERRQ(ierr)
  call VecDestroy(map_ids_2,ierr);CHKERRQ(ierr)
  call VecDestroy(map_ids_3,ierr);CHKERRQ(ierr)

end subroutine PatchCreateFlowConditionDatasetMap

! ************************************************************************** !

subroutine PatchInitConstraints(patch,reaction,option)
  !
  ! Initializes constraint concentrations
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !

  use Reaction_Aux_module

  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  type(reaction_type), pointer :: reaction

  call PatchInitCouplerConstraints(patch%initial_condition_list, &
                                   reaction,option)

  call PatchInitCouplerConstraints(patch%boundary_condition_list, &
                                   reaction,option)

  call PatchInitCouplerConstraints(patch%source_sink_list, &
                                   reaction,option)

end subroutine PatchInitConstraints

! ************************************************************************** !

subroutine PatchInitCouplerConstraints(coupler_list,reaction,option)
  !
  ! Initializes constraint concentrations
  ! for a given coupler
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !

  use Reaction_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Transport_Constraint_module

  use EOS_Water_module

  implicit none

  type(coupler_list_type), pointer :: coupler_list
  type(option_type) :: option
  type(reaction_type), pointer :: reaction

  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  class(material_auxvar_type), allocatable :: material_auxvar
  type(coupler_type), pointer :: cur_coupler
  type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
  PetscReal :: dum1
  PetscErrorCode :: ierr

  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  material_auxvar%porosity = option%reference_porosity

  cur_coupler => coupler_list%first
  do
    if (.not.associated(cur_coupler)) exit

    if (.not.associated(cur_coupler%tran_condition)) then
      option%io_buffer = 'Null transport condition found in coupler'
      if (len_trim(cur_coupler%name) > 1) then
        option%io_buffer = trim(option%io_buffer) // &
                           ' "' // trim(cur_coupler%name) // '"'
      endif
      call printErrMsg(option)
    endif

    cur_constraint_coupler => &
      cur_coupler%tran_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      global_auxvar => cur_constraint_coupler%global_auxvar
      rt_auxvar => cur_constraint_coupler%rt_auxvar
      if (associated(cur_coupler%flow_condition)) then
        if (associated(cur_coupler%flow_condition%pressure)) then
          if (associated(cur_coupler%flow_condition%pressure%dataset)) then
            global_auxvar%pres = &
              cur_coupler%flow_condition%pressure%dataset%rarray(1)
          else
            global_auxvar%pres = option%reference_pressure
          endif
        else
          global_auxvar%pres = option%reference_pressure
        endif
        if (associated(cur_coupler%flow_condition%temperature)) then
          if (associated(cur_coupler%flow_condition%temperature%dataset)) then
            global_auxvar%temp = &
              cur_coupler%flow_condition%temperature%dataset%rarray(1)
          else
            global_auxvar%temp = option%reference_temperature
          endif
        else
          global_auxvar%temp = option%reference_temperature
        endif

#ifndef CLM_PFLOTRAN
        call EOSWaterDensity(global_auxvar%temp, &
                             global_auxvar%pres(1), &
                             global_auxvar%den_kg(1), &
                             dum1,ierr)
#else
        call EOSWaterDensity(global_auxvar%temp, &
                             max(global_auxvar%pres(1), 0.01d0), &
                             global_auxvar%den_kg(1), &
                             dum1,ierr)
#endif

      else
        global_auxvar%pres = option%reference_pressure
        global_auxvar%temp = option%reference_temperature
        global_auxvar%den_kg = option%reference_density
      endif
      global_auxvar%sat = option%reference_saturation

      call ReactionEquilibrateConstraint(rt_auxvar,global_auxvar, &
                            material_auxvar, &
                            reaction,cur_constraint_coupler%constraint_name, &
                            cur_constraint_coupler%aqueous_species, &
                            cur_constraint_coupler%free_ion_guess, &
                            cur_constraint_coupler%minerals, &
                            cur_constraint_coupler%colloids, &
                            cur_constraint_coupler%immobile_species, &
                            cur_constraint_coupler%num_iterations, &
                            PETSC_FALSE,option)
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
    cur_coupler => cur_coupler%next
  enddo

  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)

end subroutine PatchInitCouplerConstraints

! ************************************************************************** !

subroutine PatchUpdateUniformVelocity(patch,velocity,option)
  !
  ! Assigns uniform velocity in connection list
  ! darcy velocities
  !
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  !

  use Option_module
  use Coupler_module
  use Condition_module
  use Connection_module

  implicit none

  type(patch_type), pointer :: patch
  PetscReal :: velocity(3)
  type(option_type), pointer :: option

  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscReal :: vdarcy

  grid => patch%grid

  ! Internal Flux Terms -----------------------------------
  cur_connection_set => grid%internal_connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      vdarcy = dot_product(velocity, &
                           cur_connection_set%dist(1:3,iconn))
      patch%internal_velocities(1,sum_connection) = vdarcy
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      vdarcy = dot_product(velocity, &
                           cur_connection_set%dist(1:3,iconn))
      patch%boundary_velocities(1,sum_connection) = vdarcy
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine PatchUpdateUniformVelocity

! ************************************************************************** !

subroutine PatchGetVariable1(patch,field,reaction,option,output_option,vec, &
                             ivar,isubvar,isubvar2)
  !
  ! PatchGetVariable: Extracts variables indexed by ivar and isubvar from a patch
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Option_module
  use Field_module
  
  use TH_Aux_module
  use Reaction_Mineral_module
  use Reaction_module
  use Reactive_Transport_Aux_module

  use Output_Aux_module
  use Variables_module
  use Material_Aux_class

  implicit none

  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: isubvar2
  PetscInt :: iphase

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: xmass, lnQKgas, ehfac, eh0, pe0, ph0, tk
  PetscReal :: tempreal
  PetscInt :: tempint, tempint2
  PetscInt :: irate, istate, irxn, ifo2, jcomp, comp_id
  PetscInt :: ivar_temp
  PetscErrorCode :: ierr


  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE,AIR_PRESSURE, &
         LIQUID_SATURATION,GAS_SATURATION,ICE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,CAPILLARY_PRESSURE,LIQUID_DENSITY_MOL, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,STATE,ICE_DENSITY, &
         EFFECTIVE_POROSITY,LIQUID_HEAD,VAPOR_PRESSURE,SATURATION_PRESSURE, &
         MAXIMUM_PRESSURE,LIQUID_MASS_FRACTION,GAS_MASS_FRACTION, &
         OIL_PRESSURE,OIL_SATURATION,OIL_DENSITY,OIL_DENSITY_MOL,OIL_ENERGY, &
         OIL_MOBILITY,OIL_VISCOSITY)

      if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%temp
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(1)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY,GAS_VISCOSITY) ! still needs implementation
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by TH')
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              if (option%use_th_freezing) then
                vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_gas
              else
                vec_ptr(local_id) = &
                  1.d0 - patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(1)
              endif
            enddo
          case(ICE_SATURATION)
            if (option%use_th_freezing) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_ice
              enddo
            else
              call printErrMsg(option,'ICE_SATURATION not supported by without freezing option TH')
            endif
          case(ICE_DENSITY)
            if (option%use_th_freezing) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%den_ice*FMWH2O
              enddo
            else
              call printErrMsg(option,'ICE_DENSITY not supported without freezing option in TH')
            endif
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%vis
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%kvr
            enddo
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by TH')
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%u
            enddo
          case(EFFECTIVE_POROSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%transient_por
            enddo
        end select

      end if

    case(PH,PE,EH,O2,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY, &
         SECONDARY_MOLARITY,TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_RATE,MINERAL_VOLUME_FRACTION,MINERAL_SATURATION_INDEX, &
         SURFACE_CMPLX,SURFACE_CMPLX_FREE,SURFACE_SITE_DENSITY, &
         KIN_SURFACE_CMPLX,KIN_SURFACE_CMPLX_FREE, PRIMARY_ACTIVITY_COEF, &
         SECONDARY_ACTIVITY_COEF,PRIMARY_KD,TOTAL_SORBED,TOTAL_SORBED_MOBILE, &
         COLLOID_MOBILE,COLLOID_IMMOBILE,AGE,TOTAL_BULK,IMMOBILE_SPECIES, &
         GAS_CONCENTRATION,REACTION_AUXILIARY)

      select case(ivar)
        case(PH)
          if (isubvar > 0) then
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                -log10(patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(isubvar)* &
                       patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar))
            enddo
          else
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
               -log10(patch%aux%RT%auxvars(ghosted_id)%sec_act_coef(-isubvar)* &
                      patch%aux%RT%auxvars(ghosted_id)%sec_molal(-isubvar))
            enddo
          endif
        case(EH)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then
              !geh: all the below should be calculated somewhere else, not in
              !     patch.F90.  most likely reactive_transport.F90
              ph0 = &
                -log10(patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(isubvar)* &
                       patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar))
              ifo2 = reaction%species_idx%o2_gas_id
              ! compute gas partial pressure
              lnQKgas = -reaction%gas%paseqlogK(ifo2)*LOG_TO_LN
              ! activity of water
              if (reaction%gas%paseqh2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%gas%paseqh2ostoich(ifo2) * &
                    patch%aux%RT%auxvars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%gas%paseqspecid(0,ifo2)
                comp_id = reaction%gas%paseqspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%gas%paseqstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%auxvars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(comp_id))
              enddo
              tk = patch%aux%Global%auxvars(ghosted_id)%temp + &
                   273.15d0
              ehfac = IDEAL_GAS_CONSTANT*tk*LOG_TO_LN/FARADAY
              eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0
              pe0 = eh0/ehfac
              vec_ptr(local_id) = eh0
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo
        case(PE)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then
              !geh: all the below should be calculated somewhere else, not in
              !     patch.F90.  most likely reactive_transport.F90
              ph0 = &
                -log10(patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(isubvar)* &
                       patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar))
              ifo2 = reaction%species_idx%o2_gas_id
              ! compute gas partial pressure
              lnQKgas = -reaction%gas%paseqlogK(ifo2)*LOG_TO_LN
              ! activity of water
              if (reaction%gas%paseqh2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%gas%paseqh2ostoich(ifo2) * &
                    patch%aux%RT%auxvars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%gas%paseqspecid(0,ifo2)
                comp_id = reaction%gas%paseqspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%gas%paseqstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%auxvars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(comp_id))
              enddo
              tk = patch%aux%Global%auxvars(ghosted_id)%temp + &
                   273.15d0
              ehfac = IDEAL_GAS_CONSTANT*tk*LOG_TO_LN/FARADAY
              eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0
              pe0 = eh0/ehfac
              vec_ptr(local_id) = pe0
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo

        case(O2)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then
              !geh: all the below should be calculated somewhere else, not in
              !     patch.F90.  most likely reactive_transport.F90
              ifo2 = reaction%species_idx%o2_gas_id
              ! compute gas partial pressure
              lnQKgas = -reaction%gas%paseqlogK(ifo2)*LOG_TO_LN
              ! activity of water
              if (reaction%gas%paseqh2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%gas%paseqh2ostoich(ifo2) * &
                    patch%aux%RT%auxvars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%gas%paseqspecid(0,ifo2)
                comp_id = reaction%gas%paseqspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%gas%paseqstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%auxvars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(comp_id))
              enddo
              vec_ptr(local_id) = lnQKgas * LN_TO_LOG
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo
        case(PRIMARY_MOLALITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%pri_molal(isubvar)
          enddo
        case(PRIMARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) * xmass * &
              (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(SECONDARY_MOLALITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar)
          enddo
        case(SECONDARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar) * xmass * &
              (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(TOTAL_MOLALITY)
          do local_id=1,grid%nlmax
            ghosted_id =grid%nL2G(local_id)
            if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            if (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) > 0.d0) then
              vec_ptr(local_id) = &
                patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) / &
                xmass / &
                (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo
        case(TOTAL_MOLARITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%total(isubvar,iphase)
          enddo
        case(TOTAL_BULK) ! mol/m^3 bulk
          ! add in total molarity and convert to mol/m^3 bulk
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) * &
              patch%aux%Material%auxvars(ghosted_id)%porosity * &
                                                             ! mol/L -> mol/m^3
              patch%aux%Global%auxvars(ghosted_id)%sat(iphase) * 1.d-3
          enddo
        case(GAS_CONCENTRATION)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%gas_pp(isubvar)
          enddo
        case(MINERAL_VOLUME_FRACTION)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%mnrl_volfrac(isubvar)
          enddo
        case(MINERAL_RATE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%mnrl_rate(isubvar)
          enddo
        case(MINERAL_SATURATION_INDEX)
          do local_id = 1, grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              RMineralSaturationIndex(isubvar, &
                                      patch%aux%RT%auxvars(ghosted_id), &
                                      patch%aux%Global%auxvars(ghosted_id), &
                                      reaction,option)
          enddo
        case(IMMOBILE_SPECIES)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%immobile(isubvar)
          enddo
        case(SURFACE_CMPLX)
          if (associated(patch%aux%RT%auxvars(1)%eqsrfcplx_conc)) then
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                                    eqsrfcplx_conc(isubvar)
            enddo
          else
            vec_ptr = UNINITIALIZED_DOUBLE
          endif
        case(SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
              srfcplxrxn_free_site_conc(isubvar)
          enddo
        case(KIN_SURFACE_CMPLX)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
              kinsrfcplx_conc(isubvar,1)
          enddo
        case(KIN_SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
              kinsrfcplx_free_site_conc(isubvar)
          enddo
        case(PRIMARY_ACTIVITY_COEF)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(isubvar)
          enddo
        case(SECONDARY_ACTIVITY_COEF)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%sec_act_coef(isubvar)
          enddo
        case(PRIMARY_KD)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call ReactionComputeKd(isubvar,vec_ptr(local_id), &
                                   patch%aux%RT%auxvars(ghosted_id), &
                                   patch%aux%Global%auxvars(ghosted_id), &
                                   patch%aux%Material%auxvars(ghosted_id), &
                                   patch%reaction,option)
          enddo
        case(TOTAL_SORBED)
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%neqsorb > 0) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%RT%auxvars(ghosted_id)%total_sorb_eq(isubvar)
              enddo
            endif
          endif
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%nsorb > 0 .and. patch%reaction%ncollcomp > 0) then
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = patch%aux%RT%auxvars(ghosted_id)%colloid% &
                total_eq_mob(isubvar)
            enddo
          endif
        case(COLLOID_MOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            do local_id=1,grid%nlmax
              ghosted_id =grid%nL2G(local_id)
              if (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) > &
                  0.d0) then
                vec_ptr(local_id) = &
                  patch%aux%RT%auxvars(ghosted_id)% &
                    colloid%conc_mob(isubvar) / &
                  (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
              else
                vec_ptr(local_id) = 0.d0
              endif
            enddo
          else
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                                    colloid%conc_mob(isubvar)
            enddo
          endif
        case(COLLOID_IMMOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            do local_id=1,grid%nlmax
              ghosted_id =grid%nL2G(local_id)
              if (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) > &
                  0.d0) then
                vec_ptr(local_id) = &
                  patch%aux%RT%auxvars(ghosted_id)% &
                    colloid%conc_imb(isubvar) / &
                  (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
              else
                vec_ptr(local_id) = 0.d0
              endif
            enddo
          else
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                                    colloid%conc_imb(isubvar)
            enddo
          endif
        case(AGE)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then
              vec_ptr(local_id) = &
                patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) / &
                patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar2) / &
                output_option%tconv
            endif
          enddo
        case(REACTION_AUXILIARY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%auxiliary_data(isubvar)
          enddo
      end select
    case(POROSITY,MINERAL_POROSITY,VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY, &
         SOIL_REFERENCE_PRESSURE)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          MaterialAuxVarGetValue(material_auxvars(grid%nL2G(local_id)),ivar)
      enddo
    case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y, PERMEABILITY_Z, &
         PERMEABILITY_XY,PERMEABILITY_XZ,PERMEABILITY_YZ, &
         GAS_PERMEABILITY,GAS_PERMEABILITY_X,GAS_PERMEABILITY_Y, &
         GAS_PERMEABILITY_Z)
      ivar_temp = ivar
      ! only liquid permeabilities in x, y, z are stored.
      select case(ivar)
        case(PERMEABILITY,GAS_PERMEABILITY,GAS_PERMEABILITY_X)
          ivar_temp = PERMEABILITY_X
        case(GAS_PERMEABILITY_Y)
          ivar_temp = PERMEABILITY_Y
        case(GAS_PERMEABILITY_Z)
          ivar_temp = PERMEABILITY_Z
      end select
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          MaterialAuxVarGetValue(material_auxvars(grid%nL2G(local_id)), &
                                 ivar_temp)
      enddo
    case(PHASE)
      call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = vec_ptr2(grid%nL2G(local_id))
      enddo
      call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%imat_internal_to_external(abs(patch%imat(grid%nL2G(local_id))))
      enddo
    case(FRACTURE)
      vec_ptr = 0.d0
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (associated(material_auxvars(ghosted_id)%fracture)) then
          if (material_auxvars(ghosted_id)%fracture%fracture_is_on) then
            vec_ptr(local_id) = 1.d0
          endif
        endif
      enddo
    case(PROCESS_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = option%myrank
      enddo
    case(NATURAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = grid%nG2A(grid%nL2G(local_id))
      enddo
    case(RESIDUAL)
      call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
      call VecStrideGather(field%flow_r,isubvar-1,vec,INSERT_VALUES,ierr)
      call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchGetVariable'')') ivar
      call printErrMsg(option)
  end select

  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PatchGetVariable1

! ************************************************************************** !

function PatchGetVariableValueAtCell(patch,field,reaction,option, &
                                     output_option,ghosted_id, &
                                     ivar,isubvar,isubvar2)
  !
  ! Returns variables indexed by ivar,
  ! isubvar, local id from Reactive Transport type
  !
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Option_module
  use Field_module

  use TH_Aux_module

  use Reactive_Transport_Aux_module  
  use Reaction_Mineral_module
  use Reaction_module
  use Reaction_Mineral_Aux_module
  use Output_Aux_module
  use Variables_module
  use Material_Aux_class

  implicit none

  PetscReal :: PatchGetVariableValueAtCell
  type(option_type), pointer :: option
  type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: tempint, tempint2
  PetscInt :: isubvar2
  PetscInt :: iphase
  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: ivar_temp

  PetscReal :: value, xmass, lnQKgas, tk, ehfac, eh0, pe0, ph0
  PetscInt :: irate, istate, irxn, ifo2, jcomp, comp_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  value = UNINITIALIZED_DOUBLE

  ! inactive grid cell
  if (patch%imat(ghosted_id) <= 0) then
    PatchGetVariableValueAtCell = 0.d0
    return
  endif

  iphase = 1
  xmass = 1.d0
  if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) &
    xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)

  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE, &
         LIQUID_SATURATION,GAS_SATURATION,ICE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,AIR_PRESSURE,CAPILLARY_PRESSURE, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,STATE,ICE_DENSITY, &
         SECONDARY_TEMPERATURE,LIQUID_DENSITY_MOL,EFFECTIVE_POROSITY, &
         LIQUID_HEAD,VAPOR_PRESSURE,SATURATION_PRESSURE,MAXIMUM_PRESSURE, &
         LIQUID_MASS_FRACTION,GAS_MASS_FRACTION, &
         OIL_PRESSURE,OIL_SATURATION,OIL_DENSITY,OIL_DENSITY_MOL,OIL_ENERGY, &
         OIL_MOBILITY,OIL_VISCOSITY)

     if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%auxvars(ghosted_id)%temp
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%auxvars(ghosted_id)%pres(1)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%auxvars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%auxvars(ghosted_id)%den_kg(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%TH%auxvars(ghosted_id)%vis
          case(LIQUID_MOBILITY)
            value = patch%aux%TH%auxvars(ghosted_id)%kvr
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY) ! still need implementation
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by TH')
          case(GAS_SATURATION)
            if (option%use_th_freezing) then
              value = patch%aux%TH%auxvars(ghosted_id)%ice%sat_gas
            else
              value = 0.d0
            endif
          case(ICE_SATURATION)
            if (option%use_th_freezing) then
              value = patch%aux%TH%auxvars(ghosted_id)%ice%sat_ice
            endif
          case(ICE_DENSITY)
            if (option%use_th_freezing) then
              value = patch%aux%TH%auxvars(ghosted_id)%ice%den_ice*FMWH2O
            endif
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by TH')
          case(LIQUID_ENERGY)
            value = patch%aux%TH%auxvars(ghosted_id)%u
          case(SECONDARY_TEMPERATURE)
            local_id = grid%nG2L(ghosted_id)
            value = patch%aux%SC_heat%sec_heat_vars(local_id)%sec_temp(isubvar)
          case(EFFECTIVE_POROSITY)
            value = patch%aux%TH%auxvars(ghosted_id)%transient_por
        end select
      endif

    case(PH,PE,EH,O2,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY, &
         SECONDARY_MOLARITY, TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_VOLUME_FRACTION,MINERAL_RATE,MINERAL_SATURATION_INDEX, &
         SURFACE_CMPLX,SURFACE_CMPLX_FREE,SURFACE_SITE_DENSITY, &
         KIN_SURFACE_CMPLX,KIN_SURFACE_CMPLX_FREE, PRIMARY_ACTIVITY_COEF, &
         SECONDARY_ACTIVITY_COEF,PRIMARY_KD, TOTAL_SORBED, &
         TOTAL_SORBED_MOBILE,COLLOID_MOBILE,COLLOID_IMMOBILE,AGE,TOTAL_BULK, &
         IMMOBILE_SPECIES,GAS_CONCENTRATION,REACTION_AUXILIARY)

      select case(ivar)
        case(PH)
          if (isubvar > 0) then
            value = -log10(patch%aux%RT%auxvars(ghosted_id)% &
                           pri_act_coef(isubvar)* &
                           patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar))
          else
            value = -log10(patch%aux%RT%auxvars(ghosted_id)% &
                           sec_act_coef(-isubvar)* &
                           patch%aux%RT%auxvars(ghosted_id)%sec_molal(-isubvar))
          endif
        case(EH)
          ph0 = -log10(patch%aux%RT%auxvars(ghosted_id)% &
                         pri_act_coef(isubvar)* &
                         patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar))

          ifo2 = reaction%species_idx%o2_gas_id

      ! compute gas partial pressure
          lnQKgas = -reaction%gas%paseqlogK(ifo2)*LOG_TO_LN

      ! activity of water
          if (reaction%gas%paseqh2oid(ifo2) > 0) then
            lnQKgas = lnQKgas + reaction%gas%paseqh2ostoich(ifo2) * &
                    patch%aux%RT%auxvars(ghosted_id)%ln_act_h2o
          endif
          do jcomp = 1, reaction%gas%paseqspecid(0,ifo2)
            comp_id = reaction%gas%paseqspecid(jcomp,ifo2)
            lnQKgas = lnQKgas + reaction%gas%paseqstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%auxvars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(comp_id))
          enddo

          tk = patch%aux%Global%auxvars(ghosted_id)%temp+273.15d0
          ehfac = IDEAL_GAS_CONSTANT*tk*LOG_TO_LN/FARADAY
          eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0

          value = eh0

        case(PE)
          ph0 = -log10(patch%aux%RT%auxvars(ghosted_id)% &
                         pri_act_coef(isubvar)* &
                         patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar))

          ifo2 = reaction%species_idx%o2_gas_id

      ! compute gas partial pressure
          lnQKgas = -reaction%gas%paseqlogK(ifo2)*LOG_TO_LN

      ! activity of water
          if (reaction%gas%paseqh2oid(ifo2) > 0) then
            lnQKgas = lnQKgas + reaction%gas%paseqh2ostoich(ifo2) * &
                    patch%aux%RT%auxvars(ghosted_id)%ln_act_h2o
          endif
          do jcomp = 1, reaction%gas%paseqspecid(0,ifo2)
            comp_id = reaction%gas%paseqspecid(jcomp,ifo2)
            lnQKgas = lnQKgas + reaction%gas%paseqstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%auxvars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(comp_id))
          enddo

          tk = patch%aux%Global%auxvars(ghosted_id)%temp+273.15d0
          ehfac = IDEAL_GAS_CONSTANT*tk*LOG_TO_LN/FARADAY
          eh0 = ehfac*(-4.d0*ph0+lnQKgas*LN_TO_LOG+logKeh(tk))/4.d0
          pe0 = eh0/ehfac
          value = pe0

        case(O2)

      ! compute gas partial pressure
              ifo2 = reaction%species_idx%o2_gas_id
              lnQKgas = -reaction%gas%paseqlogK(ifo2)*LOG_TO_LN

      ! activity of water
              if (reaction%gas%paseqh2oid(ifo2) > 0) then
                lnQKgas = lnQKgas + reaction%gas%paseqh2ostoich(ifo2) * &
                    patch%aux%RT%auxvars(ghosted_id)%ln_act_h2o
              endif
              do jcomp = 1, reaction%gas%paseqspecid(0,ifo2)
                comp_id = reaction%gas%paseqspecid(jcomp,ifo2)
                lnQKgas = lnQKgas + reaction%gas%paseqstoich(jcomp,ifo2)* &
                      log(patch%aux%RT%auxvars(ghosted_id)%pri_molal(comp_id)* &
                        patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(comp_id))
              enddo
           value = lnQKgas * LN_TO_LOG
        case(PRIMARY_MOLALITY)
          value = patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar)
        case(PRIMARY_MOLARITY)
          value = patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar)*xmass * &
                  patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) / 1000.d0
        case(SECONDARY_MOLALITY)
          value = patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar)
        case(SECONDARY_MOLARITY)
          value = patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar)*xmass * &
                  patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) / 1000.d0
        case(TOTAL_MOLALITY)
          value = patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) / &
                  xmass / &
                  patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)*1000.d0
        case(TOTAL_MOLARITY)
          value = patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase)
        case(TOTAL_BULK) ! mol/m^3 bulk
          ! add in total molarity and convert to mol/m^3 bulk
          value = &
              patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) * &
              patch%aux%Material%auxvars(ghosted_id)%porosity * &
                                                              ! mol/L -> mol/m^3
              patch%aux%Global%auxvars(ghosted_id)%sat(iphase) * 1.d-3
        case(GAS_CONCENTRATION)
          value = patch%aux%RT%auxvars(ghosted_id)%gas_pp(isubvar)
        case(MINERAL_VOLUME_FRACTION)
          value = patch%aux%RT%auxvars(ghosted_id)%mnrl_volfrac(isubvar)
        case(MINERAL_RATE)
          value = patch%aux%RT%auxvars(ghosted_id)%mnrl_rate(isubvar)
        case(MINERAL_SATURATION_INDEX)
          value = RMineralSaturationIndex(isubvar, &
                                         patch%aux%RT%auxvars(ghosted_id), &
                                         patch%aux%Global%auxvars(ghosted_id), &
                                         reaction,option)
        case(IMMOBILE_SPECIES)
          value = patch%aux%RT%auxvars(ghosted_id)%immobile(isubvar)
        case(PRIMARY_ACTIVITY_COEF)
          value = patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(isubvar)
        case(SECONDARY_ACTIVITY_COEF)
          value = patch%aux%RT%auxvars(ghosted_id)%sec_act_coef(isubvar)
        case(PRIMARY_KD)
          call ReactionComputeKd(isubvar,value, &
                                 patch%aux%RT%auxvars(ghosted_id), &
                                 patch%aux%Global%auxvars(ghosted_id), &
                                 material_auxvars(ghosted_id), &
                                 patch%reaction,option)
        case(TOTAL_SORBED)
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%neqsorb > 0) then
              value = patch%aux%RT%auxvars(ghosted_id)%total_sorb_eq(isubvar)
            endif
          endif
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%nsorb > 0 .and. patch%reaction%ncollcomp > 0) then
            value = &
              patch%aux%RT%auxvars(ghosted_id)%colloid%total_eq_mob(isubvar)
          endif
        case(COLLOID_MOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            value = patch%aux%RT%auxvars(ghosted_id)% &
                      colloid%conc_mob(isubvar) / &
                    patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)*1000.d0
          else
            value = patch%aux%RT%auxvars(ghosted_id)%colloid%conc_mob(isubvar)
          endif
        case(COLLOID_IMMOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            value = patch%aux%RT%auxvars(ghosted_id)% &
                      colloid%conc_imb(isubvar) / &
                    patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)*1000.d0
          else
            value = patch%aux%RT%auxvars(ghosted_id)%colloid%conc_imb(isubvar)
          endif
        case(AGE)
          if (patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) > &
              0.d0) then
            value = patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) / &
            patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar2) / &
            output_option%tconv
          endif
        case(REACTION_AUXILIARY)
          value = patch%aux%RT%auxvars(ghosted_id)%auxiliary_data(isubvar)
      end select
    case(POROSITY,MINERAL_POROSITY,VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY, &
         SOIL_REFERENCE_PRESSURE)
      value = MaterialAuxVarGetValue(material_auxvars(ghosted_id),ivar)
    case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y, PERMEABILITY_Z, &
         PERMEABILITY_XY,PERMEABILITY_XZ,PERMEABILITY_YZ, &
         GAS_PERMEABILITY,GAS_PERMEABILITY_X,GAS_PERMEABILITY_Y, &
         GAS_PERMEABILITY_Z)
      ivar_temp = ivar
      ! only liquid permeabilities in x, y, z are stored.
      select case(ivar)
        case(PERMEABILITY,GAS_PERMEABILITY,GAS_PERMEABILITY_X)
          ivar_temp = PERMEABILITY_X
        case(GAS_PERMEABILITY_Y)
          ivar_temp = PERMEABILITY_Y
        case(GAS_PERMEABILITY_Z)
          ivar_temp = PERMEABILITY_Z
      end select
      value = MaterialAuxVarGetValue(material_auxvars(ghosted_id),ivar_temp)
    case(PHASE)
      call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
      value = vec_ptr2(ghosted_id)
      call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
    case(MATERIAL_ID)
      value = patch%imat_internal_to_external(abs(patch%imat(ghosted_id)))
    case(FRACTURE)
      value = 0.d0
      if (associated(material_auxvars(ghosted_id)%fracture)) then
        if (material_auxvars(ghosted_id)%fracture%fracture_is_on) then
          value = 1.d0
        endif
      endif
    case(PROCESS_ID)
      value = grid%nG2A(ghosted_id)
    case(NATURAL_ID)
      value = option%myrank
    ! Need to fix the below two cases (they assume only one component) -- SK 02/06/13
    case(SECONDARY_CONCENTRATION)
      ! Note that the units are in mol/kg
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
              sec_rt_auxvar(isubvar)%pri_molal(isubvar2)
    case(SEC_MIN_VOLFRAC)
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
              sec_rt_auxvar(isubvar)%mnrl_volfrac(isubvar2)
    case(SEC_MIN_RATE)
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
              sec_rt_auxvar(isubvar)%mnrl_rate(isubvar2)
    case(SEC_MIN_SI)
      local_id = grid%nG2L(ghosted_id)
      value = RMineralSaturationIndex(isubvar2,&
                                      patch%aux%SC_RT% &
                                      sec_transport_vars(local_id)% &
                                      sec_rt_auxvar(isubvar), &
                                      patch%aux%Global%auxvars(ghosted_id),&
                                      reaction,option)
    case(RESIDUAL)
      local_id = grid%nG2L(ghosted_id)
      call VecGetArrayF90(field%flow_r,vec_ptr2,ierr);CHKERRQ(ierr)
      value = vec_ptr2((local_id-1)*option%nflowdof+isubvar)
      call VecRestoreArrayF90(field%flow_r,vec_ptr2,ierr);CHKERRQ(ierr)
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchGetVariableValueAtCell'')') &
            ivar
      call printErrMsg(option)
  end select

  PatchGetVariableValueAtCell = value

end function PatchGetVariableValueAtCell

! ************************************************************************** !

subroutine PatchSetVariable(patch,field,option,vec,vec_format,ivar,isubvar)
  !
  ! Sets variables indexed by ivar and isubvar within a patch
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Option_module
  use Field_module
  use Variables_module
  use Material_Aux_class

  implicit none

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase, istate

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

  if (vec_format == NATURAL) then
    call printErrMsg(option,&
                     'NATURAL vector format not supported by PatchSetVariable')
  endif

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE,LIQUID_SATURATION, &
         GAS_SATURATION,AIR_PRESSURE,CAPILLARY_PRESSURE, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY, &
         LIQUID_MOBILITY,GAS_MOBILITY,STATE)

      if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%temp = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%pres = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%sat = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%den_kg = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY)
            call printErrMsg(option,'GAS_MOLE_FRACTION not supported by TH')
          case(GAS_SATURATION)
            if (option%use_th_freezing) then
              if (vec_format == GLOBAL) then
                do local_id=1,grid%nlmax
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_gas = &
                    vec_ptr(local_id)
                enddo
              else if (vec_format == LOCAL) then
                do ghosted_id=1,grid%ngmax
                  patch%aux%TH%auxvars(ghosted_id)%ice%sat_gas = vec_ptr(ghosted_id)
                enddo
              endif
            endif
          case(ICE_SATURATION)
            if (option%use_th_freezing) then
              if (vec_format == GLOBAL) then
                do local_id=1,grid%nlmax
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_ice = &
                    vec_ptr(local_id)
                enddo
              else if (vec_format == LOCAL) then
                do ghosted_id=1,grid%ngmax
                  patch%aux%TH%auxvars(ghosted_id)%ice%sat_ice = vec_ptr(ghosted_id)
                enddo
              endif
            endif
          case(ICE_DENSITY)
            if (option%use_th_freezing) then
              if (vec_format == GLOBAL) then
                do local_id=1,grid%nlmax
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%den_ice = &
                    vec_ptr(local_id)
                enddo
              else if (vec_format == LOCAL) then
                do ghosted_id=1,grid%ngmax
                  patch%aux%TH%auxvars(ghosted_id)%ice%den_ice = vec_ptr(ghosted_id)
                enddo
              endif
            endif
          case(LIQUID_VISCOSITY)
          case(GAS_VISCOSITY)
          case(LIQUID_MOLE_FRACTION)
            call printErrMsg(option,'LIQUID_MOLE_FRACTION not supported by TH')
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%TH%auxvars(grid%nL2G(local_id))%u = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%TH%auxvars(ghosted_id)%u = vec_ptr(ghosted_id)
              enddo
            endif
        end select

      endif

    case(PRIMARY_MOLALITY,TOTAL_MOLARITY,MINERAL_VOLUME_FRACTION, &
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF,IMMOBILE_SPECIES, &
         GAS_CONCENTRATION,REACTION_AUXILIARY)
      select case(ivar)
        case(PRIMARY_MOLALITY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))%pri_molal(isubvar) = &
                vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) = &
                vec_ptr(ghosted_id)
            enddo
          endif
        case(TOTAL_MOLARITY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                total(isubvar,iphase) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                total(isubvar,iphase) = vec_ptr(ghosted_id)
            enddo
          endif
        case(GAS_CONCENTRATION)
          option%io_buffer = 'Active gas concentrations cannot be set in &
            &PatchSetVariable.'
          call printErrMsg(option)
        case(MINERAL_VOLUME_FRACTION)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                mnrl_volfrac(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                mnrl_volfrac(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(IMMOBILE_SPECIES)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                immobile(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                immobile(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(PRIMARY_ACTIVITY_COEF)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                pri_act_coef(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                pri_act_coef(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(SECONDARY_ACTIVITY_COEF)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                sec_act_coef(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                sec_act_coef(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(REACTION_AUXILIARY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                auxiliary_data(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                auxiliary_data(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
      end select
    case(PRIMARY_MOLARITY,SECONDARY_MOLALITY,SECONDARY_MOLARITY,TOTAL_MOLALITY, &
         COLLOID_MOBILE,COLLOID_IMMOBILE)
      select case(ivar)
        case(PRIMARY_MOLARITY)
          call printErrMsg(option,'Setting of primary molarity at grid cell not supported.')
        case(SECONDARY_MOLALITY)
          call printErrMsg(option,'Setting of secondary molality at grid cell not supported.')
        case(SECONDARY_MOLARITY)
          call printErrMsg(option,'Setting of secondary molarity at grid cell not supported.')
        case(TOTAL_MOLALITY)
          call printErrMsg(option,'Setting of total molality at grid cell not supported.')
        case(COLLOID_MOBILE)
          call printErrMsg(option,'Setting of mobile colloid concentration at grid cell not supported.')
        case(COLLOID_IMMOBILE)
          call printErrMsg(option,'Setting of immobile colloid concentration at grid cell not supported.')
      end select
    case(POROSITY,MINERAL_POROSITY)
      if (vec_format == GLOBAL) then
        do local_id=1,grid%nlmax
          call MaterialAuxVarSetValue(material_auxvars(grid%nL2G(local_id)), &
                                      ivar,vec_ptr(local_id))
        enddo
      else if (vec_format == LOCAL) then
        do ghosted_id=1,grid%ngmax
          call MaterialAuxVarSetValue(material_auxvars(ghosted_id), &
                                      ivar,vec_ptr(ghosted_id))
        enddo
      endif
    case(VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY,SOIL_REFERENCE_PRESSURE)
      option%io_buffer = 'Setting of volume, tortuosity, ' // &
        'soil compressibility or soil reference pressure in ' // &
        '"PatchSetVariable" not supported.'
      call printErrMsg(option)
    case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y,PERMEABILITY_Z, &
         GAS_PERMEABILITY,GAS_PERMEABILITY_X,GAS_PERMEABILITY_Y, &
         GAS_PERMEABILITY_Z)
      option%io_buffer = 'Setting of permeability in "PatchSetVariable"' // &
        ' not supported.'
      call printErrMsg(option)
    case(PHASE)
      if (vec_format == GLOBAL) then
        call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
        do local_id=1,grid%nlmax
          vec_ptr2(grid%nL2G(local_id)) = vec_ptr(local_id)
        enddo
        call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
      else if (vec_format == LOCAL) then
        call VecGetArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
        vec_ptr2(1:grid%ngmax) = vec_ptr(1:grid%ngmax)
        call VecRestoreArrayF90(field%iphas_loc,vec_ptr2,ierr);CHKERRQ(ierr)
      endif
    case(MATERIAL_ID)
      !geh: this would require the creation of a permanent mapping between
      !     external and internal material ids, which we want to avoid.
      call printErrMsg(option, &
                       'Cannot set MATERIAL_ID through PatchSetVariable()')
      if (vec_format == GLOBAL) then
        do local_id=1,grid%nlmax
          patch%imat(grid%nL2G(local_id)) = int(vec_ptr(local_id))
        enddo
      else if (vec_format == LOCAL) then
        patch%imat(1:grid%ngmax) = int(vec_ptr(1:grid%ngmax))
      endif
    case(PROCESS_ID)
      call printErrMsg(option, &
                       'Cannot set PROCESS_ID through PatchSetVariable()')
    case(NATURAL_ID)
      call printErrMsg(option, &
                       'Cannot set NATURAL_ID through PatchSetVariable()')
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchSetVariable'')') ivar
      call printErrMsg(option)
  end select

  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PatchSetVariable

! ************************************************************************** !

subroutine PatchCountCells(patch,total_count,active_count)
  !
  ! Counts # of active and inactive grid cells
  !
  ! Author: Glenn Hammond
  ! Date: 06/01/10
  !

  use Option_module

  implicit none

  type(patch_type) :: patch
  PetscInt :: total_count
  PetscInt :: active_count

  type(grid_type), pointer :: grid
  PetscInt :: local_id

  grid => patch%grid

  total_count = grid%nlmax

  active_count = 0
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    active_count = active_count + 1
  enddo

end subroutine PatchCountCells

! ************************************************************************** !

subroutine PatchCalculateCFL1Timestep(patch,option,max_dt_cfl_1)
  !
  ! Calculates largest time step to preserves a
  ! CFL # of 1 in a patch
  !
  ! Author: Glenn Hammond
  ! Date: 10/06/11
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  PetscReal :: max_dt_cfl_1

  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: por_sat_ave, por_sat_min, v_darcy, v_pore_ave, v_pore_max
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase

  PetscReal :: dt_cfl_1
  PetscErrorCode :: ierr

  field => patch%field
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid

  max_dt_cfl_1 = 1.d20

  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping
      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle
      distance = cur_connection_set%dist(0,iconn)
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      do iphase = 1, option%nphase
        ! if the phase is not present in either cell, skip the connection
        if (.not.(global_auxvars(ghosted_id_up)%sat(iphase) > 0.d0 .and. &
                  global_auxvars(ghosted_id_dn)%sat(iphase) > 0.d0)) cycle
        por_sat_min = min(material_auxvars(ghosted_id_up)%porosity* &
                          global_auxvars(ghosted_id_up)%sat(iphase), &
                          material_auxvars(ghosted_id_dn)%porosity* &
                          global_auxvars(ghosted_id_dn)%sat(iphase))
        por_sat_ave = (fraction_upwind* &
                       material_auxvars(ghosted_id_up)%porosity* &
                       global_auxvars(ghosted_id_up)%sat(iphase) + &
                      (1.d0-fraction_upwind)* &
                      material_auxvars(ghosted_id_dn)%porosity* &
                      global_auxvars(ghosted_id_dn)%sat(iphase))
        v_darcy = patch%internal_velocities(iphase,sum_connection)
        v_pore_max = v_darcy / por_sat_min
        v_pore_ave = v_darcy / por_sat_ave
        !geh: I use v_por_max to ensure that we limit the cfl based on the
        !     highest velocity through the face.  If porosity*saturation
        !     varies, the pore water velocity will be highest on the side
        !     of the face with the smalled value of porosity*saturation.
        dt_cfl_1 = distance / dabs(v_pore_max)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
      enddo
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
      if (patch%imat(ghosted_id_dn) <= 0) cycle
      !geh: since on boundary, dist must be scaled by 2.d0
      distance = 2.d0*cur_connection_set%dist(0,iconn)
      do iphase = 1, option%nphase
        por_sat_ave = material_auxvars(ghosted_id_dn)%porosity* &
                      global_auxvars(ghosted_id_dn)%sat(iphase)
        v_darcy = patch%boundary_velocities(iphase,sum_connection)
        v_pore_ave = v_darcy / por_sat_ave
        dt_cfl_1 = distance / dabs(v_pore_ave)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine PatchCalculateCFL1Timestep

! ************************************************************************** !

function PatchGetVarNameFromKeyword(keyword,option)
  !
  ! Returns the name of variable defined by keyword
  !
  ! Author: Glenn Hammond
  ! Date: 07/28/11
  !

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: keyword
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: PatchGetVarNameFromKeyword
  character(len=MAXSTRINGLENGTH) :: var_name

  select case(keyword)
    case('PROCESS_ID')
      var_name = 'Processor ID'
    case('NATURAL_ID')
      var_name = 'Natural ID'
    case default
      option%io_buffer = 'Keyword "' // trim(keyword) // '" not ' // &
                         'recognized in PatchGetIvarsFromKeyword()'
      call printErrMsg(option)
  end select

  PatchGetVarNameFromKeyword = var_name

end function PatchGetVarNameFromKeyword

! ************************************************************************** !

subroutine PatchGetIvarsFromKeyword(keyword,ivar,isubvar,var_type,option)
  !
  ! Returns the ivar and isubvars for extracting
  ! datasets using PatchGet/PatchSet routines
  !
  ! Author: Glenn Hammond
  ! Date: 07/28/11
  !

  use Option_module
  use Variables_module

  implicit none

  character(len=MAXWORDLENGTH) :: keyword
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: var_type
  type(option_type) :: option

  select case(keyword)
    case('PROCESS_ID')
      ivar = PROCESS_ID
      isubvar = ZERO_INTEGER
      var_type = INT_VAR
    case('NATURAL_ID')
      ivar = NATURAL_ID
      isubvar = ZERO_INTEGER
      var_type = INT_VAR
    case default
      option%io_buffer = 'Keyword "' // trim(keyword) // '" not ' // &
                         'recognized in PatchGetIvarsFromKeyword()'
      call printErrMsg(option)
  end select

end subroutine

! ************************************************************************** !

subroutine PatchGetVariable2(patch,surf_field,option,output_option,vec, &
                             ivar,isubvar,isubvar2)
  !
  ! PatchGetVariable: Extracts variables indexed by ivar and isubvar from a patch
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  !

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Option_module
  use Output_Aux_module
  use Surface_Field_module
  use Variables_module

  implicit none

  type(option_type), pointer :: option
  !type(reaction_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(surface_field_type), pointer :: surf_field
  type(patch_type), pointer :: patch
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: isubvar2
  PetscInt :: iphase

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: xmass
  PetscReal :: tempreal
  PetscInt :: tempint
  PetscInt :: irate, istate, irxn
  PetscErrorCode :: ierr

  grid => patch%grid

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

  iphase = 1

  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%surf_aux%SurfaceGlobal%auxvars(grid%nL2G(local_id))%head(1)
      enddo
    case(SURFACE_LIQUID_TEMPERATURE)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = patch%surf_aux%SurfaceGlobal%auxvars(grid%nL2G(local_id))%temp
      enddo
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%imat_internal_to_external(abs(patch%imat(grid%nL2G(local_id))))
      enddo
    case(PROCESS_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = option%myrank
      enddo
    case(NATURAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = grid%nG2A(grid%nL2G(local_id))
      enddo
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchGetVariable'')') ivar
      call printErrMsg(option)
  end select

end subroutine PatchGetVariable2

! ************************************************************************** !

subroutine PatchGetCellCenteredVelocities(patch,iphase,velocities)
  !
  ! Calculates the Darcy velocity at the center of all cells in a patch
  !
  ! Author: Glenn Hammond
  ! Date: 01/31/14
  !
  use Connection_module
  use Coupler_module

  implicit none

  type(patch_type), pointer :: patch
  PetscInt :: iphase
  PetscReal, intent(out) :: velocities(:,:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn, num_connections
  PetscReal, allocatable :: sum_area(:,:), sum_velocity(:,:)
  PetscReal :: area(3), velocity(3)
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: local_id
  PetscInt :: i

  grid => patch%grid

  allocate(sum_velocity(3,grid%nlmax))
  allocate(sum_area(3,grid%nlmax))
  sum_velocity(:,:) = 0.d0
  sum_area(:,:) = 0.d0

  ! interior velocities
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(1:3,iconn)
      velocity = patch%internal_velocities(iphase,sum_connection)*area
      if (local_id_up > 0) then
        sum_velocity(:,local_id_up) = sum_velocity(:,local_id_up) + velocity
        sum_area(:,local_id_up) = sum_area(:,local_id_up) + dabs(area)
      endif
      if (local_id_dn > 0) then
        sum_velocity(:,local_id_dn) = sum_velocity(:,local_id_dn) + velocity
        sum_area(:,local_id_dn) = sum_area(:,local_id_dn) + dabs(area)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! boundary velocities
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(1:3,iconn)
      velocity = patch%boundary_velocities(iphase,sum_connection)*area
      sum_velocity(:,local_id) = sum_velocity(:,local_id) + velocity
      sum_area(:,local_id) = sum_area(:,local_id) + dabs(area)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! divide by total area
  do local_id=1,grid%nlmax
    do i=1,3
      if (sum_area(i,local_id) > 0.d0) then
        velocities(i,local_id) = sum_velocity(i,local_id) / &
                                 sum_area(i,local_id)
      else
        velocities(i,local_id) = 0.d0
      endif
    enddo
  enddo

  deallocate(sum_velocity)
  deallocate(sum_area)

end subroutine PatchGetCellCenteredVelocities

! ************************************************************************** !

function PatchGetConnectionsFromCoords(patch,coordinates,integral_flux_name, &
                                       option)
  !
  !
  ! Returns a list of internal and boundary connection ids for cell
  ! interfaces within a polygon.
  !
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  !
  use Option_module
  use Geometry_module
  use Utility_module
  use Connection_module
  use Coupler_module

  implicit none

  type(patch_type) :: patch
  type(point3d_type) :: coordinates(:)
  character(len=MAXWORDLENGTH) :: integral_flux_name
  type(option_type) :: option

  PetscInt, pointer :: PatchGetConnectionsFromCoords(:)

  PetscInt, pointer :: connections(:)
  type(grid_type), pointer :: grid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition

  PetscInt :: idir
  PetscInt :: icount
  PetscInt :: array_size
  PetscInt :: sum_connection
  PetscInt :: iconn
  PetscInt :: i
  PetscInt :: local_id
  PetscInt :: local_id_up
  PetscInt :: ghosted_id
  PetscInt :: ghosted_id_up
  PetscInt :: ghosted_id_dn
  PetscReal :: fraction_upwind
  PetscReal :: magnitude
  PetscReal :: v1(3), v2(3)
  PetscReal :: x, y, z
  PetscReal :: value1, value2
  PetscReal, parameter :: relative_tolerance = 1.d-6
  PetscBool :: within_tolerance
  PetscErrorCode :: ierr

  grid => patch%grid

  ! determine orientation of polygon
  if (size(coordinates) > 2) then
    v1(1) = coordinates(2)%x - coordinates(1)%x
    v1(2) = coordinates(2)%y - coordinates(1)%y
    v1(3) = coordinates(2)%z - coordinates(1)%z
    v2(1) = coordinates(2)%x - coordinates(3)%x
    v2(2) = coordinates(2)%y - coordinates(3)%y
    v2(3) = coordinates(2)%z - coordinates(3)%z
    v1 = CrossProduct(v1,v2)
    icount = 0
    idir = 0
    do i = X_DIRECTION, Z_DIRECTION
      if (v1(i) > 1.d-10) then
        icount = icount + 1
        idir = i
      endif
    enddo
  else
    v1(1) = coordinates(1)%x
    v1(2) = coordinates(1)%y
    v1(3) = coordinates(1)%z
    v2(1) = coordinates(2)%x
    v2(2) = coordinates(2)%y
    v2(3) = coordinates(2)%z
    icount = 0
    do i = X_DIRECTION, Z_DIRECTION
      if (Equal(v2(i),v1(i))) then
        idir = i
        icount = icount + 1
      endif
    enddo
    if (icount == 0) icount = 3
  endif

  if (icount > 1) then
    option%io_buffer = 'Rectangle defined in integral flux "' // &
      trim(adjustl(integral_flux_name)) // &
      '" must be aligned with structured grid coordinates axes.'
    call printErrMsg(option)
  endif

  array_size = 100
  allocate(connections(array_size))
  icount = 0
  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      ! if one of the cells is ghosted, the process stores the flux only
      ! when the upwind cell is non-ghosted.
      if (local_id_up <= 0) cycle

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      magnitude = cur_connection_set%dist(0,iconn)
      x = grid%x(ghosted_id_up) + fraction_upwind * magnitude * &
          cur_connection_set%dist(X_DIRECTION,iconn)
      y = grid%y(ghosted_id_up) + fraction_upwind * magnitude * &
          cur_connection_set%dist(Y_DIRECTION,iconn)
      z = grid%z(ghosted_id_up) + fraction_upwind * magnitude * &
          cur_connection_set%dist(Z_DIRECTION,iconn)
      select case(idir)
        case(X_DIRECTION)
          value1 = x
          value2 = coordinates(1)%x
        case(Y_DIRECTION)
          value1 = y
          value2 = coordinates(1)%y
        case(Z_DIRECTION)
          value1 = z
          value2 = coordinates(1)%z
      end select
      within_tolerance = PETSC_FALSE
      if (Equal(value1,0.d0)) then
        within_tolerance = Equal(value1,value2)
      else
        within_tolerance = dabs((value1-value2)/value1) < relative_tolerance
      endif
      if (within_tolerance .and. &
          GeometryPointInPolygon(x,y,z,idir,coordinates)) then
        icount = icount + 1
        if (icount > size(connections)) then
          call reallocateIntArray(connections,array_size)
        endif
        connections(icount) = sum_connection
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      fraction_upwind = 1.d0
      magnitude = cur_connection_set%dist(0,iconn)
      x = grid%x(ghosted_id) - fraction_upwind * magnitude * &
                               cur_connection_set%dist(X_DIRECTION,iconn)
      y = grid%y(ghosted_id) - fraction_upwind * magnitude * &
                               cur_connection_set%dist(Y_DIRECTION,iconn)
      z = grid%z(ghosted_id) - fraction_upwind * magnitude * &
                               cur_connection_set%dist(Z_DIRECTION,iconn)
      select case(idir)
        case(X_DIRECTION)
          value1 = x
          value2 = coordinates(1)%x
        case(Y_DIRECTION)
          value1 = y
          value2 = coordinates(1)%y
        case(Z_DIRECTION)
          value1 = z
          value2 = coordinates(1)%z
      end select
      within_tolerance = PETSC_FALSE
      if (Equal(value1,0.d0)) then
        within_tolerance = Equal(value1,value2)
      else
        within_tolerance = dabs((value1-value2)/value1) < relative_tolerance
      endif
      if (within_tolerance .and. &
          GeometryPointInPolygon(x,y,z,idir,coordinates)) then
        icount = icount + 1
        if (icount > size(connections)) then
          call reallocateIntArray(connections,array_size)
        endif
        connections(icount) = -1 * sum_connection
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  nullify(PatchGetConnectionsFromCoords)
  if (icount > 0) then
    allocate(PatchGetConnectionsFromCoords(icount))
    PatchGetConnectionsFromCoords = connections(1:icount)
  endif
  deallocate(connections)
  nullify(connections)

  call MPI_Allreduce(icount,i,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                     option%mycomm,ierr)
  if (i == 0) then
    option%io_buffer = 'Zero connections found for INTEGRAL_FLUX "' // &
      trim(adjustl(integral_flux_name)) // &
      '".  Please ensure that the coordinates coincide with a cell boundary.'
    call printErrMsg(option)
  endif

end function PatchGetConnectionsFromCoords

! **************************************************************************** !

subroutine PatchCouplerInputRecord(patch)
  !
  ! Prints ingested coupler information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 04/18/2016
  !
  use Coupler_module

  implicit none

  type(patch_type), pointer :: patch

  type(coupler_type), pointer :: cur_coupler
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT

  k = 0

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'INITIAL CONDITIONS'

  ! Initial conditions
  cur_coupler => patch%initial_condition_list%first
  do
    if (.not.associated(cur_coupler)) exit
    k = k + 1
    write(id,'(a29)',advance='no') 'initial condition listed: '
    write(word1,*) k
    write(id,'(a)') '#' // adjustl(trim(word1))
    write(id,'(a29)',advance='no') 'applies to region: '
    write(id,'(a)') adjustl(trim(cur_coupler%region_name))
    if (len_trim(cur_coupler%flow_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'flow condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%flow_condition_name))
    endif
    if (len_trim(cur_coupler%tran_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'transport condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%tran_condition_name))
    endif
    write(id,'(a29)') '---------------------------: '
    cur_coupler => cur_coupler%next
  enddo

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'BOUNDARY CONDITIONS'

  ! Boundary conditions
  cur_coupler => patch%boundary_condition_list%first
  do
    if (.not.associated(cur_coupler)) exit
    write(id,'(a29)',advance='no') 'boundary condition name: '
    write(id,'(a)') adjustl(trim(cur_coupler%name))
    write(id,'(a29)',advance='no') 'applies to region: '
    write(id,'(a)') adjustl(trim(cur_coupler%region_name))
    if (len_trim(cur_coupler%flow_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'flow condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%flow_condition_name))
    endif
    if (len_trim(cur_coupler%tran_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'transport condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%tran_condition_name))
    endif
    write(id,'(a29)') '---------------------------: '
    cur_coupler => cur_coupler%next
  enddo

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'SOURCE-SINKS'

  ! Source-Sink conditions
  cur_coupler => patch%source_sink_list%first
  do
    if (.not.associated(cur_coupler)) exit
    write(id,'(a29)',advance='no') 'source-sink name: '
    write(id,'(a)') adjustl(trim(cur_coupler%name))
    write(id,'(a29)',advance='no') 'applies to region: '
    write(id,'(a)') adjustl(trim(cur_coupler%region_name))
    if (len_trim(cur_coupler%flow_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'flow condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%flow_condition_name))
    endif
    if (len_trim(cur_coupler%tran_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'transport condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%tran_condition_name))
    endif
    write(id,'(a29)') '---------------------------: '
    cur_coupler => cur_coupler%next
  enddo

end subroutine PatchCouplerInputRecord

! **************************************************************************** !

subroutine PatchGetCompMassInRegion(cell_ids,num_cells,patch,option, &
                                    global_total_mass)
  !
  ! Calculates the total mass (aqueous, sorbed, and precipitated) in a region
  ! in units of mol.
  !
  ! Author: Jenn Frederick
  ! Date: 04/25/2016
  !
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Aux_module
  use Grid_module
  use Option_module
  use Reactive_Transport_Aux_module

  implicit none

  PetscInt, pointer :: cell_ids(:)
  PetscInt :: num_cells
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscReal :: global_total_mass  ! [mol]

  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(reaction_type), pointer :: reaction
  PetscReal :: aq_species_mass    ! [mol]
  PetscReal :: sorb_species_mass  ! [mol]
  PetscReal :: ppt_species_mass   ! [mol]
  PetscReal :: m3_water           ! [m^3-water]
  PetscReal :: m3_bulk            ! [m^3-bulk]
  PetscInt :: k, j, m
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  PetscReal :: local_total_mass

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  rt_auxvars => patch%aux%RT%auxvars
  reaction => patch%reaction
  local_total_mass = 0.d0
  global_total_mass = 0.d0

  ! Loop through all cells in the region:
  do k = 1,num_cells
    local_id = cell_ids(k)
    ghosted_id = patch%grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    m3_water = material_auxvars(ghosted_id)%porosity * &         ! [-]
               global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &  ! [water]
               material_auxvars(ghosted_id)%volume               ! [m^3-bulk]
    m3_bulk = material_auxvars(ghosted_id)%volume                ! [m^3-bulk]
    ! Loop through aqueous and sorbed species:
    do j = 1,reaction%ncomp
      aq_species_mass = 0.d0
      sorb_species_mass = 0.d0
      ! aqueous species; units [mol/L-water]*[m^3-water]*[1000L/m^3-water]=[mol]
      aq_species_mass = rt_auxvars(ghosted_id)%total(j,LIQUID_PHASE) * &
                        m3_water * 1.0d3
      if (associated(rt_auxvars(ghosted_id)%total_sorb_eq)) then
        ! sorbed species; units [mol/m^3-bulk]*[m^3-bulk]=[mol]
        sorb_species_mass = rt_auxvars(ghosted_id)%total_sorb_eq(j) * m3_bulk
      else
        sorb_species_mass = 0.d0
      endif
      local_total_mass = local_total_mass + aq_species_mass + &
                                            sorb_species_mass
    enddo
    ! Loop through precipitated species:
    do m = 1,reaction%mineral%nkinmnrl
      ppt_species_mass = 0.d0
      ! precip. species; units [m^3-mnrl/m^3-bulk]*[m^3-bulk]/[m^3-mnrl/mol-mnrl]=[mol]
      ppt_species_mass = rt_auxvars(ghosted_id)%mnrl_volfrac(m) * m3_bulk / &
                         reaction%mineral%kinmnrl_molar_vol(m)
      local_total_mass = local_total_mass + ppt_species_mass
    enddo
  enddo ! Cell loop

  ! Sum the local_total_mass across all processes that own the region:
  call MPI_Allreduce(local_total_mass,global_total_mass,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

end subroutine PatchGetCompMassInRegion

! **************************************************************************** !

subroutine PatchGetWaterMassInRegion(cell_ids,num_cells,patch,option, &
                                     global_water_mass)
  !
  ! Calculates the water mass in a region in kg
  !
  ! Author: Satish Karra
  ! Date: 09/20/2016
  !
  use Global_Aux_module
  use Material_Aux_class
  use Grid_module
  use Option_module

  implicit none

  PetscInt, pointer :: cell_ids(:)
  PetscInt :: num_cells
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscReal :: global_water_mass

  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal :: m3_water, kg_water
  PetscInt :: k, j, m
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  PetscReal :: local_water_mass

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  local_water_mass = 0.d0
  global_water_mass = 0.d0

  ! Loop through all cells in the region:
  do k = 1,num_cells
    local_id = cell_ids(k)
    ghosted_id = patch%grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    m3_water = material_auxvars(ghosted_id)%porosity * &         ! [-]
               global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &  ! [water]
               material_auxvars(ghosted_id)%volume               ! [m^3-bulk]
    kg_water = m3_water*global_auxvars(ghosted_id)% &            ! [m^3-water]
               den(LIQUID_PHASE)                                 ! [kg/m^3-water]
    local_water_mass = local_water_mass + kg_water
  enddo ! Cell loop

  ! Sum the local_water_mass across all processes that own the region:
  call MPI_Allreduce(local_water_mass,global_water_mass,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

end subroutine PatchGetWaterMassInRegion

! **************************************************************************** !

subroutine PatchGetCompMassInRegionAssign(region_list, &
           mass_balance_region_list,option)
  !
  ! Assigns patch%region information to the mass balance region object
  !
  ! Author: Jenn Frederick
  ! Date: 04/26/2016
  !
  use Output_Aux_module
  use Region_module
  use String_module

  implicit none

  type(region_list_type), pointer :: region_list
  type(mass_balance_region_type), pointer :: mass_balance_region_list
  type(option_type), pointer :: option

  type(region_type), pointer :: cur_region
  type(mass_balance_region_type), pointer :: cur_mbr
  PetscBool :: success

  cur_mbr => mass_balance_region_list
  do
    if (.not.associated(cur_mbr)) exit
    ! Loop through patch%region_list to find wanted region:
    cur_region => region_list%first
    do
      if (.not.associated(cur_region)) exit
      success = PETSC_TRUE
      if (StringCompareIgnoreCase(cur_region%name,cur_mbr%region_name)) exit
      success = PETSC_FALSE
      cur_region => cur_region%next
    enddo
    ! If the wanted region was not found, throw an error msg:
    if (.not.success) then
      option%io_buffer = 'Region ' // trim(cur_mbr%region_name) // ' not &
                          &found among listed regions.'
      call printErrMsg(option)
    endif
    ! Assign the mass balance region the wanted region's info:
    cur_mbr%num_cells = cur_region%num_cells
    cur_mbr%region_cell_ids => cur_region%cell_ids
    ! Go to next mass balance region
    cur_mbr => cur_mbr%next
  enddo

end subroutine PatchGetCompMassInRegionAssign

! ************************************************************************** !

subroutine PatchVerifyDatasetGriddedForFlux(dataset,coupler,option)
  !
  ! Verifies that a dataset being used to define fluxes adheres to
  ! several rules that attempt to minimize mass balance error.
  !
  ! Author: Jennifer Frederick
  ! Date: 11/04/2016
  !

  use Option_module
  use Coupler_module
  use Dataset_Gridded_HDF5_class

  implicit none

  class(dataset_gridded_hdf5_type) :: dataset
  type(coupler_type) :: coupler
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscInt :: i, dataset_size

  ! check if dataset is cell-centered:
  if (.not.dataset%is_cell_centered) then
    option%io_buffer = 'Dataset ' // trim(dataset%hdf5_dataset_name) // &
      " must be cell-centered for fluxes. You must set attribute: &
      &h5grp.attrs['Cell Centered'] = True."
    call printErrMsg(option)
  endif
  ! check if the dimensions match:
  dataset_size = 1
  do i=1,size(dataset%dims)
    dataset_size = dataset%dims(i)*dataset_size
  enddo
  if (coupler%connection_set%num_connections /= dataset_size) then
    write(string,*) dataset%dims
    write(string2,*) coupler%connection_set%num_connections
    option%io_buffer = 'Dataset ' // trim(dataset%hdf5_dataset_name) // &
      " must have a value for each cell on the boundary defined by &
      &REGION " // trim(coupler%region%name) // '. The dataset dimension &
      &is ' // adjustl(trim(string)) // ' but the number of boundary &
      &connections is ' // adjustl(trim(string2)) // '.'
    call printErrMsg(option)
  endif
  ! check if the interpolation method is STEP:
  if (.not.dataset%interpolation_method == INTERPOLATION_STEP) then
    option%io_buffer = 'Dataset ' // trim(dataset%hdf5_dataset_name) // &
      " must be assigned the STEP interpolation method for fluxes. You &
      &must set attribute: h5grp.attrs['Interpolation Method'] = &
      &np.string_('STEP')."
    call printErrMsg(option)
  endif

end subroutine PatchVerifyDatasetGriddedForFlux

! ************************************************************************** !

subroutine PatchDestroyList(patch_list)
  !
  ! Deallocates a patch list and array of patches
  !
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  !

  implicit none

  type(patch_list_type), pointer :: patch_list

  type(patch_type), pointer :: cur_patch, prev_patch

  if (.not.associated(patch_list)) return

  if (associated(patch_list%array)) deallocate(patch_list%array)
  nullify(patch_list%array)

  cur_patch => patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    prev_patch => cur_patch
    cur_patch => cur_patch%next
    call PatchDestroy(prev_patch)
  enddo

  nullify(patch_list%first)
  nullify(patch_list%last)
  patch_list%num_patch_objects = 0

  deallocate(patch_list)
  nullify(patch_list)

end subroutine PatchDestroyList

! ************************************************************************** !

subroutine PatchDestroy(patch)
  !
  ! Deallocates a patch object
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Utility_module, only : DeallocateArray

  implicit none

  type(patch_type), pointer :: patch

  call DeallocateArray(patch%imat)
  call DeallocateArray(patch%imat_internal_to_external)
  call DeallocateArray(patch%sat_func_id)
  call DeallocateArray(patch%internal_velocities)
  call DeallocateArray(patch%boundary_velocities)
  call DeallocateArray(patch%internal_tran_coefs)
  call DeallocateArray(patch%boundary_tran_coefs)
  call DeallocateArray(patch%internal_flow_fluxes)
  call DeallocateArray(patch%boundary_flow_fluxes)
  call DeallocateArray(patch%ss_flow_fluxes)
  call DeallocateArray(patch%internal_tran_fluxes)
  call DeallocateArray(patch%boundary_tran_fluxes)
  call DeallocateArray(patch%ss_tran_fluxes)
  call DeallocateArray(patch%ss_flow_vol_fluxes)

  call DeallocateArray(patch%boundary_energy_flux)


  if (associated(patch%material_property_array)) &
    deallocate(patch%material_property_array)
  nullify(patch%material_property_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%material_properties)

  if (associated(patch%characteristic_curves_array)) &
    deallocate(patch%characteristic_curves_array)
  nullify(patch%characteristic_curves_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%characteristic_curves)

  nullify(patch%surf_field)
  if (associated(patch%surf_material_property_array)) &
    deallocate(patch%surf_material_property_array)
  nullify(patch%surf_material_property_array)
  nullify(patch%surf_material_properties)

  ! solely nullify grid since destroyed in discretization
  nullify(patch%grid)
  call RegionDestroyList(patch%region_list)
  call CouplerDestroyList(patch%boundary_condition_list)
  call CouplerDestroyList(patch%initial_condition_list)
  call CouplerDestroyList(patch%source_sink_list)

  call ObservationDestroyList(patch%observation_list)
  call IntegralFluxDestroyList(patch%integral_flux_list)
  call StrataDestroyList(patch%strata_list)

  call AuxDestroy(patch%aux)
  call SurfaceAuxDestroy(patch%surf_aux)

  ! these are solely pointers, must not destroy.
  nullify(patch%reaction)
  nullify(patch%datasets)
  nullify(patch%field)

  deallocate(patch)
  nullify(patch)

end subroutine PatchDestroy

end module Patch_module
