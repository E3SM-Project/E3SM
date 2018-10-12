module Realization_Surface_class

  use Realization_Base_class
  
  use Condition_module
  use Debug_module
  use Discretization_module
  use Input_Aux_module
  use Option_module
  use Patch_module
  use Region_module
  use Surface_Field_module
  use Surface_Material_module
  use Dataset_Base_class
  use Reaction_Aux_module
  use Output_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
private


#include "petsc/finclude/petscsys.h"

  PetscReal, parameter :: eps       = 1.D-8

  type, public, extends(realization_base_type) :: realization_surface_type

    type(surface_field_type), pointer :: surf_field
    type(region_list_type), pointer :: surf_regions
    type(condition_list_type),pointer :: surf_flow_conditions
    type(tran_condition_list_type),pointer :: surf_transport_conditions
    type(surface_material_property_type), pointer :: surf_material_properties
    type(surface_material_property_ptr_type), pointer :: surf_material_property_array(:)
    type(reaction_type), pointer :: surf_reaction
    character(len=MAXSTRINGLENGTH) :: surf_filename
    character(len=MAXSTRINGLENGTH) :: subsurf_filename

    class(dataset_base_type), pointer :: datasets
    
    PetscReal :: dt_max
    PetscReal :: dt_init
    PetscReal :: dt_coupling
    
    PetscInt :: iter_count
    PetscBool :: first_time

  end type realization_surface_type

  !         123456789+123456789+123456789+1
  public :: RealizSurfCreate, &
            RealizSurfDestroy, &
            RealizSurfStrip, &
            RealizSurfAddWaypointsToList, &
            RealizSurfCreateDiscretization, &
            RealizSurfAddCoupler, &
            RealizSurfAddStrata, &
            RealizSurfLocalizeRegions, &
            RealizSurfPassFieldPtrToPatches, &
            RealizSurfProcessCouplers, &
            RealizSurfLocalToLocalWithArray, &
            RealizSurfProcessConditions, &
            RealizSurfProcessFlowConditions, &
            RealizSurfMapSurfSubsurfGrids, &
            RealizSurfInitAllCouplerAuxVars, &
            RealizSurfAllCouplerAuxVars, &
            RealizSurfProcessMatProp, &
            RealizSurfUpdate, &
!            RealizSurfCreateSurfSubsurfVec, &
!            RealizSurfUpdateSubsurfBC, &
!            RealizSurfUpdateSurfBC, &
!            RealizSurfSurf2SubsurfFlux, &
            RealizSurfGetVariable

  !TODO(intel)
!  public :: SurfaceRealizationGetVariable     
  
!  interface SurfaceRealizationGetVariable
!    module procedure :: RealizationGetVariable ! from Realization_Base_class
!  end interface
  
contains

! ************************************************************************** !

function RealizSurfCreate(option)
  ! 
  ! This routine allocates and initializes a new SurfaceRealization object
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/16/12
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none

  type(option_type), pointer :: option
  class(realization_surface_type),pointer :: RealizSurfCreate
  class(realization_surface_type),pointer :: surf_realization
  
  allocate(surf_realization)
  call RealizationBaseInit(surf_realization,option)
  surf_realization%option => option

  surf_realization%surf_field => SurfaceFieldCreate()
  !geh: debug, output_option, patch_list already allocated in 
  !     RealizationBaseInit()
  !geh: surf_realization%debug => DebugCreate()
  !geh: surf_realization%output_option => OutputOptionCreate()
  !geh: surf_realization%patch_list => PatchCreateList()

  nullify(surf_realization%surf_material_properties)
  nullify(surf_realization%surf_material_property_array)

  allocate(surf_realization%surf_regions)
  call RegionInitList(surf_realization%surf_regions)
  
  allocate(surf_realization%surf_flow_conditions)
  call FlowConditionInitList(surf_realization%surf_flow_conditions)
  allocate(surf_realization%surf_transport_conditions)
  call TranConditionInitList(surf_realization%surf_transport_conditions)
  
  nullify(surf_realization%surf_reaction)
  nullify(surf_realization%datasets)

  surf_realization%iter_count = 0
  surf_realization%dt_init = 1.d0
  surf_realization%dt_max = 1.d0
  surf_realization%dt_coupling = 0.d0
  
  surf_realization%first_time = PETSC_TRUE
  RealizSurfCreate => surf_realization

end function RealizSurfCreate

! ************************************************************************** !

subroutine RealizSurfAddCoupler(surf_realization,coupler)
  ! 
  ! This routine adds a copy of a coupler to a list
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/10/12
  ! 

  use Coupler_module

  implicit none
  
  class(realization_surface_type) :: surf_realization
  type(coupler_type), pointer :: coupler
  
  type(patch_type), pointer :: cur_patch
  type(coupler_type), pointer :: new_coupler
  
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    ! only add to flow list for now, since they will be split out later
    new_coupler => CouplerCreate(coupler)
    select case(coupler%itype)
      case(BOUNDARY_COUPLER_TYPE)
        call CouplerAddToList(new_coupler,cur_patch%boundary_condition_list)
      case(INITIAL_COUPLER_TYPE)
        call CouplerAddToList(new_coupler,cur_patch%initial_condition_list)
      case(SRC_SINK_COUPLER_TYPE)
        call CouplerAddToList(new_coupler,cur_patch%source_sink_list)
    end select
    nullify(new_coupler)
    cur_patch => cur_patch%next
  enddo

  call CouplerDestroy(coupler)
 
end subroutine RealizSurfAddCoupler

! ************************************************************************** !

subroutine RealizSurfProcessCouplers(surf_realization)
  ! 
  ! This routine sets connectivity and pointers for couplers related to
  ! surface flow.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/17/12
  ! 

  use Option_module

  implicit none
  
  class(realization_surface_type) :: surf_realization
  type(patch_type), pointer :: cur_patch
  
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    call PatchProcessCouplers(cur_patch,surf_realization%surf_flow_conditions, &
                              surf_realization%surf_transport_conditions, &
                              surf_realization%option)
    cur_patch => cur_patch%next
  enddo

end subroutine RealizSurfProcessCouplers

! ************************************************************************** !

subroutine RealizSurfProcessMatProp(surf_realization)
  ! 
  ! This routine sets up linkeage between surface material properties
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/17/12
  ! 

  use String_module
  
  implicit none
  
  class(realization_surface_type) :: surf_realization
  
  PetscBool :: found
  PetscInt :: i
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string

  type(patch_type), pointer :: cur_patch
  
  option => surf_realization%option
  
  ! organize lists
  call SurfaceMaterialPropConvertListToArray( &
                                surf_realization%surf_material_properties, &
                                surf_realization%surf_material_property_array, &
                                option)

  ! set up mirrored pointer arrays within patches to saturation functions
  ! and material properties
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    cur_patch%surf_material_properties => surf_realization%surf_material_properties
    call SurfaceMaterialPropConvertListToArray( &
                                    cur_patch%surf_material_properties, &
                                    cur_patch%surf_material_property_array, &
                                    option)
    ! create mapping of internal to external material id
    call SurfaceMaterialCreateIntToExtMapping(cur_patch%surf_material_property_array, &
                                              cur_patch%imat_internal_to_external)

    cur_patch => cur_patch%next
  enddo
  
end subroutine RealizSurfProcessMatProp

! ************************************************************************** !

subroutine RealizSurfLocalizeRegions(surf_realization)
  ! 
  ! This routine localizes surface regions within each patch
  ! (similar to RealizationLocalizeRegions)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/17/12
  ! 

  use Option_module
  use String_module
  use Grid_module

  implicit none
  
  class(realization_surface_type) :: surf_realization
  
  type(patch_type), pointer :: cur_patch
  type (region_type), pointer :: cur_region
  type(option_type), pointer :: option
  type(region_type), pointer :: patch_region

  option => surf_realization%option

  ! localize the regions on each patch
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    call PatchLocalizeRegions(cur_patch,surf_realization%surf_regions, &
                              surf_realization%option)
    cur_patch => cur_patch%next
  enddo
 
end subroutine RealizSurfLocalizeRegions

! ************************************************************************** !

subroutine RealizSurfAddStrata(surf_realization,strata)
  ! 
  ! This routine adds a copy of a strata to a list
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/17/12
  ! 

  use Strata_module

  implicit none
  
  class(realization_surface_type) :: surf_realization
  type(strata_type), pointer :: strata
  
  type(patch_type), pointer :: cur_patch
  type(strata_type), pointer :: new_strata
  
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    new_strata => StrataCreate(strata)
    call StrataAddToList(new_strata,cur_patch%strata_list)
    nullify(new_strata)
    cur_patch => cur_patch%next
  enddo
  
  call StrataDestroy(strata)
 
end subroutine RealizSurfAddStrata

! ************************************************************************** !

subroutine RealizSurfCreateDiscretization(surf_realization)
  ! 
  ! This routine creates grid
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/17/12
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Grid_module
  use Grid_Unstructured_Aux_module, only : UGridMapIndices
  use Grid_Unstructured_module, only     : UGridEnsureRightHandRule
  use Coupler_module
  use Discretization_module
  use Grid_Unstructured_Cell_module
  use DM_Kludge_module
  
  implicit none

  class(realization_surface_type) :: surf_realization
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  type(dm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr

  option => surf_realization%option
  surf_field => surf_realization%surf_field
  discretization => surf_realization%discretization

  call DiscretizationCreateDMs(discretization, option%nsurfflowdof, &
                               ZERO_INTEGER, ZERO_INTEGER, &
                               ZERO_INTEGER, ZERO_INTEGER, &
                               option)

  ! n degree of freedom, global
  call DiscretizationCreateVector(discretization,NFLOWDOF,surf_field%flow_xx, &
                                  GLOBAL,option)
  call VecSet(surf_field%flow_xx,0.d0,ierr);CHKERRQ(ierr)

  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_yy)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_dxx)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_r)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%flow_accum)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx, &
                                     surf_field%work)

  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,surf_field%mannings0, &
                                  GLOBAL,option)
  call VecSet(surf_field%mannings0,0.d0,ierr);CHKERRQ(ierr)

   call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%area)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%press_subsurf)
  call DiscretizationDuplicateVector(discretization,surf_field%mannings0, &
                                     surf_field%temp_subsurf)
  ! n degrees of freedom, local
  call DiscretizationCreateVector(discretization,NFLOWDOF,surf_field%flow_xx_loc, &
                                  LOCAL,option)
  call VecSet(surf_field%flow_xx_loc,0.d0,ierr);CHKERRQ(ierr)
  call DiscretizationDuplicateVector(discretization,surf_field%flow_xx_loc, &
                                     surf_field%work_loc)

  ! 1-dof degrees of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,surf_field%mannings_loc, &
                                  LOCAL,option)
  call VecSet(surf_field%mannings_loc,0.d0,ierr);CHKERRQ(ierr)

  grid => discretization%grid

  ! set up nG2L, NL2G, etc.
  call UGridMapIndices(grid%unstructured_grid,discretization%dm_1dof%ugdm, &
                        grid%nG2L,grid%nL2G,grid%nG2A,option)
  call GridComputeCoordinates(grid,discretization%origin_global,option, &
                              discretization%dm_1dof%ugdm) 
  call UGridEnsureRightHandRule(grid%unstructured_grid,grid%x, &
                                grid%y,grid%z,grid%nG2A,grid%nL2G,option)

  ! set up internal connectivity, distance, etc.
  call GridComputeInternalConnect(grid,option,discretization%dm_1dof%ugdm) 
  call GridComputeAreas(grid,surf_field%area,option)
  call GridPrintExtents(grid,option)

  ! Allocate vectors to hold flowrate quantities
  if (surf_realization%output_option%print_hdf5_mass_flowrate.or. &
     surf_realization%output_option%print_hdf5_energy_flowrate.or. &
     surf_realization%output_option%print_hdf5_aveg_mass_flowrate.or. &
     surf_realization%output_option%print_hdf5_aveg_energy_flowrate) then

    call VecCreateMPI(option%mycomm, &
         (option%nflowdof*MAX_FACE_PER_CELL_SURF+1)*surf_realization%patch%grid%nlmax, &
          PETSC_DETERMINE,surf_field%flowrate_inst,ierr);CHKERRQ(ierr)
    call VecSet(surf_field%flowrate_inst,0.d0,ierr);CHKERRQ(ierr)

    ! If average flowrate has to be saved, create a vector for it
    if (surf_realization%output_option%print_hdf5_aveg_mass_flowrate.or. &
       surf_realization%output_option%print_hdf5_aveg_energy_flowrate) then
      call VecCreateMPI(option%mycomm, &
          (option%nflowdof*MAX_FACE_PER_CELL_SURF+1)*surf_realization%patch%grid%nlmax, &
          PETSC_DETERMINE,surf_field%flowrate_aveg,ierr);CHKERRQ(ierr)
    call VecSet(surf_field%flowrate_aveg,0.d0,ierr);CHKERRQ(ierr)
    endif
  endif

end subroutine RealizSurfCreateDiscretization

! ************************************************************************** !

subroutine RealizSurfPassFieldPtrToPatches(surf_realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/19/12
  ! 

  use Option_module

  implicit none
  
  class(realization_surface_type) :: surf_realization

  type(patch_type), pointer :: cur_patch

  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    cur_patch%surf_field => surf_realization%surf_field
    cur_patch => cur_patch%next
  enddo
  
end subroutine RealizSurfPassFieldPtrToPatches

! ************************************************************************** !

subroutine RealizSurfProcessConditions(surf_realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/19/12
  ! 

  implicit none
  
  class(realization_surface_type) :: surf_realization
  
  if (surf_realization%option%nflowdof > 0) then
    call RealizSurfProcessFlowConditions(surf_realization)
  endif
  if (surf_realization%option%ntrandof > 0) then
    !call SurfaceRealProcessTranConditions(surf_realization)
  endif
 
end subroutine RealizSurfProcessConditions

! ************************************************************************** !

subroutine RealizSurfLocalToLocalWithArray(surf_realization,array_id)
  ! 
  ! This routine takes an F90 array that is ghosted and updates the ghosted
  ! values (similar to RealLocalToLocalWithArray)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/13/12
  ! 

  use Grid_module
  use Surface_Field_module

  implicit none

  class(realization_surface_type) :: surf_realization
  PetscInt :: array_id
  
  type(patch_type), pointer :: cur_patch
  type(grid_type), pointer :: grid
  type(surface_field_type), pointer :: surf_field

  surf_field => surf_realization%surf_field

  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    grid => cur_patch%grid
    select case(array_id)
      case(MATERIAL_ID_ARRAY)
        call GridCopyIntegerArrayToVec(grid, cur_patch%imat,surf_field%work_loc, &
                                        grid%ngmax)
    end select
    cur_patch => cur_patch%next
  enddo
  call DiscretizationLocalToLocal(surf_realization%discretization, &
                                  surf_field%work_loc, &
                                  surf_field%work_loc,ONEDOF)
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    grid => cur_patch%grid

    select case(array_id)
      case(MATERIAL_ID_ARRAY)
        call GridCopyVecToIntegerArray(grid, cur_patch%imat,surf_field%work_loc, &
                                        grid%ngmax)
    end select
    cur_patch => cur_patch%next
  enddo

end subroutine RealizSurfLocalToLocalWithArray

! ************************************************************************** !

subroutine RealizSurfProcessFlowConditions(surf_realization)
  ! 
  ! This routine sets linkage of flow conditions to dataset
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/20/12
  ! 

  use Dataset_Base_class
  use Dataset_module

  implicit none

  class(realization_surface_type) :: surf_realization
  
  type(flow_condition_type), pointer :: cur_surf_flow_condition
  type(flow_sub_condition_type), pointer :: cur_surf_flow_sub_condition
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: dataset_name
  class(dataset_base_type), pointer :: dataset
  PetscInt :: i
  
  option => surf_realization%option
  
  ! loop over flow conditions looking for linkage to datasets
  cur_surf_flow_condition => surf_realization%surf_flow_conditions%first
  do
    if (.not.associated(cur_surf_flow_condition)) exit
    string = 'flow_condition ' // trim(cur_surf_flow_condition%name)
    ! find datum dataset
    call DatasetFindInList(surf_realization%datasets, &
                           cur_surf_flow_condition%datum, &
                           cur_surf_flow_condition%default_time_storage, &
                           string,option)
    select case(option%iflowmode)
      case(TH_MODE)
        do i = 1, size(cur_surf_flow_condition%sub_condition_ptr)
           ! find dataset
          call DatasetFindInList(surf_realization%datasets, &
                 cur_surf_flow_condition%sub_condition_ptr(i)%ptr%dataset, &
                 cur_surf_flow_condition%default_time_storage, &
                 string,option)
          ! find gradient dataset
          call DatasetFindInList(surf_realization%datasets, &
                 cur_surf_flow_condition%sub_condition_ptr(i)%ptr%gradient, &
                 cur_surf_flow_condition%default_time_storage, &
                 string,option)
        enddo
      case default
        option%io_buffer='RealizSurfProcessFlowConditions not implemented in this mode'
        call printErrMsg(option)
    end select
    cur_surf_flow_condition => cur_surf_flow_condition%next
  enddo

end subroutine RealizSurfProcessFlowConditions

! ************************************************************************** !

subroutine RealizSurfInitAllCouplerAuxVars(surf_realization)
  ! 
  ! This routine initializez coupler auxillary variables within list
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/21/12
  ! 

  use Option_module

  implicit none
  
  class(realization_surface_type) :: surf_realization
  
  type(patch_type), pointer :: cur_patch

  call FlowConditionUpdate(surf_realization%surf_flow_conditions, &
                           surf_realization%option)

  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    call PatchInitAllCouplerAuxVars(cur_patch, &
                                    surf_realization%option)
    cur_patch => cur_patch%next
  enddo

end subroutine RealizSurfInitAllCouplerAuxVars

! ************************************************************************** !

subroutine RealizSurfAllCouplerAuxVars(surf_realization,force_update_flag)
  ! 
  ! This routine updates auxiliary variables associated with couplers in the
  ! list.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/18/13
  ! 

  use Option_module

  implicit none

  class(realization_surface_type) :: surf_realization
  PetscBool :: force_update_flag

  call PatchUpdateAllCouplerAuxVars(surf_realization%patch,force_update_flag, &
                                    surf_realization%option)

end subroutine RealizSurfAllCouplerAuxVars

! ************************************************************************** !

subroutine RealizSurfMapSurfSubsurfGrids(realization,surf_realization)
  ! 
  ! This routine creates vector scatter contexts between surface and subsurface
  ! grids.
  ! Algorithm:
  ! - It uses a similar logic of Matrix-Vector multiplication used in
  ! UGridMapSideSet() subroutine. The algorithm here is extended to use
  ! Matrix-Matrix mulitplication
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 01/17/12
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Grid_module
  use String_module
  use Grid_Unstructured_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Region_module
  use DM_Kludge_module, only : dm_ptr_type

  implicit none

  class(realization_subsurface_type), pointer :: realization
  class(realization_surface_type), pointer :: surf_realization

  type(option_type), pointer           :: option
  type(grid_unstructured_type),pointer :: subsurf_grid
  type(grid_unstructured_type),pointer :: surf_grid
  type(patch_type), pointer            :: cur_patch 
  type(region_type), pointer           :: cur_region, top_region
  type(region_type), pointer           :: patch_region
  type(dm_ptr_type), pointer :: dm_ptr

  Mat :: Mat_vert_to_face_subsurf
  Mat :: Mat_vert_to_face_subsurf_transp
  Mat :: Mat_vert_to_face_surf
  Mat :: Mat_vert_to_face_surf_transp
  Mat :: prod
  Vec :: subsurf_petsc_ids, surf_petsc_ids
  Vec :: subsurf_nat_ids, surf_nat_ids
  Vec :: corr_subsurf_nat_ids, corr_surf_nat_ids

  PetscViewer :: viewer


  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer ::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4)
  PetscInt :: nvertices
  PetscInt :: iface
  PetscInt :: local_id, ii, jj
  PetscInt :: cell_type
  PetscInt :: ivertex, vertex_id_local
  PetscReal :: real_array4(4)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)

  PetscErrorCode :: ierr
  PetscBool :: found
  
  found = PETSC_FALSE
  
  if (.not.associated(realization)) return
  
  option => realization%option
  subsurf_grid => realization%discretization%grid%unstructured_grid
  surf_grid    => surf_realization%discretization%grid%unstructured_grid

  ! localize the regions on each patch
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    cur_region => cur_patch%region_list%first
      do
        if (.not.associated(cur_region)) exit
        if (StringCompare(cur_region%name,'top')) then
          found = PETSC_TRUE
          top_region => cur_region
          exit
        endif
        cur_region => cur_region%next
      enddo
    cur_patch => cur_patch%next
  enddo

  if (found.eqv.PETSC_FALSE) then
    option%io_buffer = 'When running with -DSURFACE_FLOW need to specify ' // &
      ' in the inputfile explicitly region: top '
    call printErrMsg(option)
  endif

  call MatCreateAIJ(option%mycomm, &
                       top_region%num_cells, &
                       PETSC_DECIDE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_subsurf, &
                       ierr);CHKERRQ(ierr)

   call MatCreateAIJ(option%mycomm, &
                     PETSC_DECIDE, &
                     top_region%num_cells, &
                     subsurf_grid%num_vertices_global, &
                     PETSC_DETERMINE, &
                     12, &
                     PETSC_NULL_INTEGER, &
                     12, &
                     PETSC_NULL_INTEGER, &
                     Mat_vert_to_face_subsurf_transp, &
                     ierr);CHKERRQ(ierr)

  call VecCreateMPI(option%mycomm,top_region%num_cells,PETSC_DETERMINE, &
                    subsurf_petsc_ids,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,top_region%num_cells,PETSC_DETERMINE, &
                    subsurf_nat_ids,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,top_region%num_cells,PETSC_DETERMINE, &
                    corr_surf_nat_ids,ierr);CHKERRQ(ierr)
  call MatZeroEntries(Mat_vert_to_face_subsurf,ierr);CHKERRQ(ierr)
  real_array4 = 1.d0

  call VecGetArrayF90(subsurf_petsc_ids,vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(subsurf_nat_ids,vec_ptr2,ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(top_region%num_cells,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  do ii = 1, top_region%num_cells
    local_id = top_region%cell_ids(ii)
    vec_ptr(ii) = subsurf_grid%cell_ids_petsc(local_id)
    iface    = top_region%faces(ii)
    cell_type = subsurf_grid%cell_type(local_id)
    vec_ptr2(ii) = subsurf_grid%cell_ids_natural(local_id)
    !nfaces = UCellGetNFaces(cell_type,option)

    call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                    int_array4)
    ! For this matrix:
    !   irow = local face id
    !   icol = natural (global) vertex id
    do ivertex = 1, nvertices
      vertex_id_local = &
        subsurf_grid%cell_vertices(int_array4(ivertex),local_id)
      int_array4_0(ivertex) = &
        subsurf_grid%vertex_ids_natural(vertex_id_local)-1
    enddo
    call MatSetValues(Mat_vert_to_face_subsurf, &
                      1,ii-1+offset, &
                      nvertices,int_array4_0, &
                      real_array4, &
                      INSERT_VALUES,ierr);CHKERRQ(ierr)
    call MatSetValues(Mat_vert_to_face_subsurf_transp, &
                      nvertices,int_array4_0, &
                      1,ii-1+offset, &
                      real_array4, &
                      INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo

  call VecRestoreArrayF90(subsurf_petsc_ids,vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(subsurf_nat_ids,vec_ptr2,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_face_subsurf,MAT_FINAL_ASSEMBLY, &
                      ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(Mat_vert_to_face_subsurf_transp,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_face_subsurf_transp,MAT_FINAL_ASSEMBLY, &
                      ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_subsurf.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_vert_to_face_subsurf,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  !call MatTranspose(Mat_vert_to_face_subsurf,MAT_INITIAL_MATRIX, &
  !                  Mat_vert_to_face_subsurf_transp,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_subsurf_transp.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_vert_to_face_subsurf_transp,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  string = 'subsurf_petsc_ids.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecView(subsurf_petsc_ids,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call MatCreateAIJ(option%mycomm, &
                       surf_grid%nlmax, &
                       PETSC_DECIDE, &
                       PETSC_DETERMINE, &
                       subsurf_grid%num_vertices_global, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       4, &
                       PETSC_NULL_INTEGER, &
                       Mat_vert_to_face_surf, &
                       ierr);CHKERRQ(ierr)

  call MatCreateAIJ(option%mycomm, &
                    PETSC_DECIDE, &
                    surf_grid%nlmax, &
                    subsurf_grid%num_vertices_global, &
                    PETSC_DETERMINE, &
                    12, &
                    PETSC_NULL_INTEGER, &
                    12, &
                    PETSC_NULL_INTEGER, &
                    Mat_vert_to_face_surf_transp, &
                    ierr);CHKERRQ(ierr)

  call VecCreateMPI(option%mycomm,surf_grid%nlmax,PETSC_DETERMINE, &
                    surf_petsc_ids,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,surf_grid%nlmax,PETSC_DETERMINE, &
                    surf_nat_ids,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,surf_grid%nlmax,PETSC_DETERMINE, &
                    corr_subsurf_nat_ids,ierr);CHKERRQ(ierr)
  offset=0
  call MPI_Exscan(surf_grid%nlmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  CHKERRQ(ierr)

  call VecGetArrayF90(surf_petsc_ids,vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(surf_nat_ids,vec_ptr2,ierr);CHKERRQ(ierr)

  do local_id = 1, surf_grid%nlmax
    cell_type = surf_grid%cell_type(local_id)
    vec_ptr(local_id) = surf_grid%cell_ids_petsc(local_id)
    vec_ptr2(local_id)= surf_grid%cell_ids_natural(local_id)
    
    int_array4_0 = 0
    nvertices = surf_grid%cell_vertices(0,local_id)
    do ivertex = 1, nvertices
      vertex_id_local = surf_grid%cell_vertices(ivertex,local_id)
      int_array4_0(ivertex) = &
        surf_grid%vertex_ids_natural(vertex_id_local)-1
    enddo    
   call MatSetValues(Mat_vert_to_face_surf,1,local_id-1+offset, &
                     nvertices,int_array4_0,real_array4, &
                     INSERT_VALUES,ierr);CHKERRQ(ierr)
   call MatSetValues(Mat_vert_to_face_surf_transp, &
                     nvertices,int_array4_0, &
                     1,local_id-1+offset, &
                     real_array4, &
                     INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo

  call VecRestoreArrayF90(surf_petsc_ids,vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(surf_nat_ids,vec_ptr2,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_face_surf,MAT_FINAL_ASSEMBLY, &
                      ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(Mat_vert_to_face_surf_transp,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Mat_vert_to_face_surf_transp,MAT_FINAL_ASSEMBLY, &
                      ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_surf.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_vert_to_face_surf,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  string = 'surf_petsc_ids.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecView(surf_petsc_ids,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  !call MatTranspose(Mat_vert_to_face_surf,MAT_INITIAL_MATRIX, &
  !                  Mat_vert_to_face_surf_transp,ierr)

#if UGRID_DEBUG
  string = 'Mat_vert_to_face_surf_transp.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(Mat_vert_to_face_surf_transp,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  
  
  call MatMatMult(Mat_vert_to_face_subsurf,Mat_vert_to_face_surf_transp, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,prod, &
                  ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  string = 'prod.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(prod,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  

  call RealizSurfMapSurfSubsurfGrid(realization, surf_realization, prod, TWO_DIM_GRID, &
                                     surf_petsc_ids)
  call MatDestroy(prod,ierr);CHKERRQ(ierr)

  call MatMatMult(Mat_vert_to_face_surf,Mat_vert_to_face_subsurf_transp, &
                  MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,prod, &
                  ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  string = 'prod_2.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call MatView(prod,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif  
  call RealizSurfMapSurfSubsurfGrid(realization, surf_realization, prod, THREE_DIM_GRID, &
                                        subsurf_petsc_ids)

  ! For each control volume in surface mesh, get the corresponding natural ids of
  ! subsurface control volume
  dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
  call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                       subsurf_nat_ids, corr_subsurf_nat_ids, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                       subsurf_nat_ids, corr_subsurf_nat_ids, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,ONEDOF)
  call VecScatterBegin(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                       surf_nat_ids, corr_surf_nat_ids, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(dm_ptr%ugdm%scatter_bet_grids_1dof, &
                       surf_nat_ids, corr_surf_nat_ids, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr)

  ! Save the natural ids
  allocate(surf_grid%nat_ids_of_other_grid(surf_grid%nlmax))
  allocate(subsurf_grid%nat_ids_of_other_grid(top_region%num_cells))

  call VecGetArrayF90(corr_subsurf_nat_ids,vec_ptr,ierr)
  call VecGetArrayF90(corr_surf_nat_ids,vec_ptr2,ierr)
  surf_grid%nat_ids_of_other_grid = int(vec_ptr)
  subsurf_grid%nat_ids_of_other_grid = int(vec_ptr2)
  call VecRestoreArrayF90(corr_subsurf_nat_ids,vec_ptr,ierr)
  call VecRestoreArrayF90(corr_surf_nat_ids,vec_ptr2,ierr)

  ! Free up the memory
  call MatDestroy(prod,ierr);CHKERRQ(ierr)

  call MatDestroy(Mat_vert_to_face_subsurf,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_vert_to_face_subsurf_transp,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_vert_to_face_surf,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_vert_to_face_surf_transp,ierr);CHKERRQ(ierr)

  call VecDestroy(subsurf_petsc_ids,ierr);CHKERRQ(ierr)
  call VecDestroy(surf_petsc_ids,ierr);CHKERRQ(ierr)
  
end subroutine RealizSurfMapSurfSubsurfGrids

! ************************************************************************** !

subroutine RealizSurfMapSurfSubsurfGrid( &
              realization, &       !<
              surf_realization, &  !<
              prod_mat, &          !< Mat-Mat-Mult matrix
              source_grid_flag, &  !< To identify a surface or subsurface grid
              source_petsc_ids &   !< MPI-Vector containing cell ids in PETSc order
              )
  ! 
  ! This subroutine creates a single vector scatter context
  ! Algorithm:
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 01/18/12
  ! 

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Grid_module
  use String_module
  use Grid_Unstructured_module
  use Grid_Unstructured_Cell_module
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Surface_Field_module
  use Grid_Unstructured_module
  use Discretization_module
  use Grid_Unstructured_Aux_module
  use DM_Kludge_module

  implicit none

  class(realization_subsurface_type), pointer :: realization
  class(realization_surface_type), pointer :: surf_realization
  Mat :: prod_mat
  PetscInt :: source_grid_flag
  Vec :: source_petsc_ids

  Mat :: prod_loc_mat
  Vec :: source_loc_vec
  Vec :: corr_dest_ids_vec
  Vec :: corr_dest_ids_vec_ndof
  Vec :: source_petsc_ids_ndof
  IS :: is_tmp1, is_tmp2
  IS :: is_tmp3, is_tmp4
  PetscInt,pointer :: corr_v2_ids(:)
  VecScatter :: scatter
  VecScatter :: scatter_ndof

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(surface_field_type), pointer :: surf_field

  type(dm_ptr_type), pointer :: dm_ptr
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer ::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4)
  PetscReal :: real_array4(4)
  PetscInt :: ii, jj
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: ivertex, cell_id, vertex_id_local
  PetscReal :: max_value

  PetscInt, pointer :: ia_p(:), ja_p(:)
  PetscInt :: nrow,rstart,rend,icol(1)
  PetscInt :: index
  PetscInt :: vertex_id
  PetscOffset :: iia,jja,iicol
  PetscBool :: done
  PetscScalar, pointer :: aa_v(:)
  PetscInt :: row, col

  PetscErrorCode :: ierr
  PetscBool :: found
  PetscInt :: nlocal

  option     => realization%option
  field      => realization%field
  surf_field => surf_realization%surf_field
  
  if (option%mycommsize > 1) then
    ! From the MPI-Matrix get the local-matrix
    call MatMPIAIJGetLocalMat(prod_mat,MAT_INITIAL_MATRIX,prod_loc_mat, &
                              ierr);CHKERRQ(ierr)
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_loc_mat, ONE_INTEGER, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArrayF90(prod_loc_mat,aa_v,ierr);CHKERRQ(ierr)
  else
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_mat, ONE_INTEGER, PETSC_FALSE, PETSC_FALSE, &
                        nrow, ia_p, ja_p, done, ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArrayF90(prod_mat,aa_v,ierr);CHKERRQ(ierr)
  endif

  ! For each row of the local-matrix, find the column with the largest value
  allocate(corr_v2_ids(nrow))
  row = 1
  col = 0
  do ii = 1, nrow
    max_value = 0.d0
    do jj = ia_p(ii), ia_p(ii + 1) - 1
      if (aa_v(jj) > max_value) then
        corr_v2_ids(ii) = ja_p(jj)
        max_value = aa_v(jj)
      endif
    enddo
    if (max_value<3) then
      option%io_buffer = 'Atleast three vertices need to form a face'
      call printErrMsg(option)
    endif
  enddo

  offset = 0
  call MPI_Exscan(nrow,offset,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  allocate(int_array(nrow))
  do ii = 1, nrow
    int_array(ii) = (ii-1)+offset
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,option%nflowdof,nrow, &
                     int_array,PETSC_COPY_VALUES,is_tmp3,ierr);CHKERRQ(ierr)

  do ii = 1, nrow
    int_array(ii) = corr_v2_ids(ii)-1
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp2,ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,option%nflowdof,nrow, &
                     int_array,PETSC_COPY_VALUES,is_tmp4,ierr);CHKERRQ(ierr)
  deallocate(int_array)
  
  call VecCreateMPI(option%mycomm,nrow,PETSC_DETERMINE, &
                    corr_dest_ids_vec,ierr);CHKERRQ(ierr)
  call VecScatterCreate(source_petsc_ids,is_tmp2,corr_dest_ids_vec,is_tmp1, &
                        scatter,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp1,ierr);CHKERRQ(ierr)
  call ISDestroy(is_tmp2,ierr);CHKERRQ(ierr)
  
  call VecScatterBegin(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter,source_petsc_ids,corr_dest_ids_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  select case(source_grid_flag)
    case(TWO_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,ONEDOF)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids, &
                          ierr);CHKERRQ(ierr)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids_1dof, &
                          ierr);CHKERRQ(ierr)
    case(THREE_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,ONEDOF)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids, &
                          ierr);CHKERRQ(ierr)
      call VecScatterCopy(scatter,dm_ptr%ugdm%scatter_bet_grids_1dof, &
                          ierr);CHKERRQ(ierr)
  end select
  call VecScatterDestroy(scatter,ierr);CHKERRQ(ierr)

#if UGRID_DEBUG
  if (source_grid_flag==TWO_DIM_GRID) write(string,*) 'surf'
  if (source_grid_flag==THREE_DIM_GRID) write(string,*) 'subsurf'
  string = adjustl(string)
  string = 'corr_dest_ids_vec_' // trim(string) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecView(corr_dest_ids_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecDestroy(corr_dest_ids_vec,ierr);CHKERRQ(ierr)
  if (option%mycommsize>1) then
    call MatSeqAIJRestoreArrayF90(prod_loc_mat,aa_v,ierr);CHKERRQ(ierr)
    call MatDestroy(prod_loc_mat,ierr);CHKERRQ(ierr)
  else
    call MatSeqAIJRestoreArrayF90(prod_mat,aa_v,ierr);CHKERRQ(ierr)
  endif

#if UGRID_DEBUG
  if (source_grid_flag==TWO_DIM_GRID) write(string,*) 'surf'
  if (source_grid_flag==THREE_DIM_GRID) write(string,*) 'subsurf'
  string = adjustl(string)
  string = 'scatter_bet_grids_' // trim(string) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,string,viewer,ierr);CHKERRQ(ierr)
  call VecScatterView(dm_ptr%ugdm%scatter_bet_grids,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  ! Create stridded vectors
  call VecCreate(option%mycomm,corr_dest_ids_vec_ndof,ierr);CHKERRQ(ierr)
  call VecSetSizes(corr_dest_ids_vec_ndof,nrow*option%nflowdof, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(corr_dest_ids_vec_ndof,option%nflowdof, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(corr_dest_ids_vec_ndof,ierr);CHKERRQ(ierr)

  call VecGetLocalSize(source_petsc_ids,nlocal,ierr);CHKERRQ(ierr)
  call VecCreate(option%mycomm,source_petsc_ids_ndof,ierr);CHKERRQ(ierr)
  call VecSetSizes(source_petsc_ids_ndof,nlocal*option%nflowdof, &
                  PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(source_petsc_ids_ndof,option%nflowdof, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(source_petsc_ids_ndof,ierr);CHKERRQ(ierr)

  ! Create stridded vectors-scatter context
  call VecScatterCreate(source_petsc_ids_ndof,is_tmp4, &
                        corr_dest_ids_vec_ndof,is_tmp3, &
                        scatter_ndof,ierr);CHKERRQ(ierr)

  ! Save the stridded vectors-scatter context
  select case(source_grid_flag)
    case(TWO_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(surf_realization%discretization,NFLOWDOF)
    case(THREE_DIM_GRID)
      dm_ptr => DiscretizationGetDMPtrFromIndex(realization%discretization,NFLOWDOF)
  end select
  call VecScatterCopy(scatter_ndof,dm_ptr%ugdm%scatter_bet_grids_ndof, &
                      ierr);CHKERRQ(ierr)

  ! Cleanup
  call VecScatterDestroy(scatter_ndof,ierr);CHKERRQ(ierr)
  call VecDestroy(source_petsc_ids_ndof,ierr);CHKERRQ(ierr)
  call VecDestroy(corr_dest_ids_vec_ndof,ierr);CHKERRQ(ierr)

end subroutine RealizSurfMapSurfSubsurfGrid

! ************************************************************************** !

subroutine RealizSurfDestroy(surf_realization)
  ! 
  ! This routine destroys RealizSurf object
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/16/12
  ! 

  implicit none
  
  class(realization_surface_type), pointer :: surf_realization
  
  if (.not.associated(surf_realization)) return
  
  !geh: deallocate everything in base
  call RealizationBaseStrip(surf_realization)
  
  call SurfaceFieldDestroy(surf_realization%surf_field)

  call OutputOptionDestroy(surf_realization%output_option)
  
  call RegionDestroyList(surf_realization%surf_regions)
  
  call FlowConditionDestroyList(surf_realization%surf_flow_conditions)
  call TranConditionDestroyList(surf_realization%surf_transport_conditions)
  
  call PatchDestroyList(surf_realization%patch_list)
  
  if (associated(surf_realization%debug)) deallocate(surf_realization%debug)
  nullify(surf_realization%debug)
  
  if (associated(surf_realization%surf_material_property_array)) &
    deallocate(surf_realization%surf_material_property_array)
  nullify(surf_realization%surf_material_property_array)
  call SurfaceMaterialPropertyDestroy(surf_realization%surf_material_properties)
  
  call DiscretizationDestroy(surf_realization%discretization)

  if (associated(surf_realization)) deallocate(surf_realization)
  nullify(surf_realization)

end subroutine RealizSurfDestroy


! ************************************************************************** !

subroutine RealizSurfStrip(surf_realization)
  ! 
  ! This routine destroys RealizSurf object
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/16/12
  ! 

  implicit none
  
  class(realization_surface_type), pointer :: surf_realization
  
  if (.not.associated(surf_realization)) return
  
  !geh: deallocate everything in base
  call RealizationBaseStrip(surf_realization)
  
  call SurfaceFieldDestroy(surf_realization%surf_field)

  call OutputOptionDestroy(surf_realization%output_option)
  
  call RegionDestroyList(surf_realization%surf_regions)
  
  call FlowConditionDestroyList(surf_realization%surf_flow_conditions)
  call TranConditionDestroyList(surf_realization%surf_transport_conditions)
  
  call PatchDestroyList(surf_realization%patch_list)
  
  if (associated(surf_realization%debug)) deallocate(surf_realization%debug)
  nullify(surf_realization%debug)
  
  if (associated(surf_realization%surf_material_property_array)) &
    deallocate(surf_realization%surf_material_property_array)
  nullify(surf_realization%surf_material_property_array)
  call SurfaceMaterialPropertyDestroy(surf_realization%surf_material_properties)
  
  call DiscretizationDestroy(surf_realization%discretization)

  call ReactionDestroy(surf_realization%reaction,surf_realization%option)

end subroutine RealizSurfStrip

! ************************************************************************** !

subroutine RealizSurfUpdate(surf_realization)
  ! 
  ! This routine updates parameters in realization (eg. conditions, bcs, srcs)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/22/12
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  
  class(realization_surface_type) :: surf_realization

  PetscBool :: force_update_flag = PETSC_FALSE

  ! must update conditions first
  call FlowConditionUpdate(surf_realization%surf_flow_conditions, &
                           surf_realization%option)

  call RealizSurfAllCouplerAuxVars(surf_realization,force_update_flag)

end subroutine RealizSurfUpdate

! ************************************************************************** !

subroutine RealizSurfGetVariable(surf_realization,vec,ivar,isubvar,isubvar1)
  ! 
  ! This routine extracts variables indexed by ivar and isubvar from surface
  ! realization.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/22/12
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Surface_Field_module

  implicit none

  class(realization_surface_type) :: surf_realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1

  call PatchGetVariable(surf_realization%patch, &
                       surf_realization%surf_field, &
                       !surf_realization%reaction, &
                       surf_realization%option, &
                       surf_realization%output_option, &
                       vec,ivar,isubvar,isubvar1)

end subroutine RealizSurfGetVariable

! ************************************************************************** !

subroutine RealizSurfAddWaypointsToList(surf_realization,waypoint_list)
  ! 
  ! This routine creates waypoints assocated with source/sink, boundary
  ! condition, etc. and adds to a list
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/15/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Waypoint_module
  use Time_Storage_module  

  implicit none
  
  class(realization_surface_type) :: surf_realization
  type(waypoint_list_type) :: waypoint_list

  type(flow_condition_type), pointer :: cur_flow_condition
  type(flow_sub_condition_type), pointer :: sub_condition
  type(waypoint_type), pointer :: waypoint, cur_waypoint
  type(option_type), pointer :: option
  PetscInt :: itime, isub_condition
  PetscReal :: temp_real, final_time
  PetscReal, pointer :: times(:)

  option => surf_realization%option
  nullify(times)
  
  ! set flag for final output
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final) then
      cur_waypoint%print_snap_output = &
        surf_realization%output_option%print_final_snap
      exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  ! use final time in conditional below
  if (associated(cur_waypoint)) then
    final_time = cur_waypoint%time
  else
    option%io_buffer = 'Final time not found in RealizSurfAddWaypointsToList'
    call printErrMsg(option)
  endif

  ! add update of flow conditions
  cur_flow_condition => surf_realization%surf_flow_conditions%first
  do
    if (.not.associated(cur_flow_condition)) exit
    if (cur_flow_condition%sync_time_with_update) then
      do isub_condition = 1, cur_flow_condition%num_sub_conditions
        sub_condition => cur_flow_condition%sub_condition_ptr(isub_condition)%ptr
        !TODO(geh): check if this updated more than simply the flow_dataset (i.e. datum and gradient)
        !geh: followup - no, datum/gradient are not considered.  Should they be considered?
        call TimeStorageGetTimes(sub_condition%dataset%time_storage, option, &
                                final_time, times)        
        if (associated(times)) then
          if (size(times) > 1000) then
            option%io_buffer = 'For flow condition "' // &
              trim(cur_flow_condition%name) // &
              '" dataset "' // trim(sub_condition%name) // &
              '", the number of times is excessive for synchronization ' // &
              'with waypoints.'
            call printErrMsg(option)
          endif
          do itime = 1, size(times)
            waypoint => WaypointCreate()
            waypoint%time = times(itime)
            waypoint%update_conditions = PETSC_TRUE
            call WaypointInsertInList(waypoint,waypoint_list)
          enddo
          deallocate(times)
          nullify(times)
        endif
      enddo
    endif
    cur_flow_condition => cur_flow_condition%next
  enddo
     
end subroutine RealizSurfAddWaypointsToList

end module Realization_Surface_class
