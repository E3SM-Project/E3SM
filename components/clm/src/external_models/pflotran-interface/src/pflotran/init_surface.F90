module Init_Surface_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none


  public :: SurfaceInitReadRequiredCards, &
            InitSurfaceSetupRealization, &
            InitSurfaceSetupSolvers
contains

! ************************************************************************** !

subroutine SurfaceInitReadRequiredCards(surf_realization,input)
  ! 
  ! This routine reads the required input file cards related to surface flows
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/18/12
  ! 

  use Option_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module

  use Realization_Surface_class
  use Surface_Auxiliary_module

  implicit none

  class(realization_surface_type) :: surf_realization
  type(input_type), pointer :: input

  type(discretization_type), pointer :: discretization
  character(len=MAXSTRINGLENGTH) :: string
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  patch          => surf_realization%patch
  option         => surf_realization%option
  discretization => surf_realization%discretization
  
! Read in select required cards
!.........................................................................
 
  ! GRID information
!  string = "GRID"
!  call InputFindStringInFile(input,option,string)
!  call InputFindStringErrorMsg(input,option,string)

  ! SURFACE_FLOW information
  string = "SURFACE_FLOW"
  call InputFindStringInFile(input,option,string)
  if (InputError(input)) return
  option%surf_flow_on = PETSC_TRUE
  option%nsurfflowdof = 1
  
  string = "SURF_GRID"
  call InputFindStringInFile(input,option,string)
!  call SurfaceFlowReadRequiredCardsFromInput(surf_realization,input,option)
  call SurfaceInit(surf_realization,input,option)

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      patch%surf_or_subsurf_flag = SURFACE
      if (.not.associated(surf_realization%patch_list)) then
        surf_realization%patch_list => PatchCreateList()
      endif
      call PatchAddToList(patch,surf_realization%patch_list)
      surf_realization%patch => patch
  end select
    
end subroutine SurfaceInitReadRequiredCards

! ************************************************************************** !

subroutine SurfaceInit(surf_realization,input,option)
  ! 
  ! This routine reads required surface flow data from the input file
  ! grids.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/09/12
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Surface_Material_module
  use Realization_Surface_class
  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Grid_Unstructured_Aux_module
  use Discretization_module
  use Region_module
  use Condition_module
  use Grid_Unstructured_Aux_module

  implicit none

  class(realization_surface_type) :: surf_realization
  type(discretization_type),pointer :: discretization
  type(grid_type), pointer :: grid
  type(input_type), pointer :: input
  type(option_type) :: option
  type(grid_unstructured_type), pointer :: un_str_sfgrid
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: unstructured_grid_ctype
  PetscInt :: unstructured_grid_itype

  discretization => surf_realization%discretization

  input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  call InputReadPflotranString(input,option)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'keyword','SURFACE_FLOW')
  call StringToUpper(word)
    
  select case(trim(word))
    case ('TYPE')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'keyword','TYPE')
      call StringToUpper(word)

      select case(trim(word))
        case ('UNSTRUCTURED')
          unstructured_grid_itype = IMPLICIT_UNSTRUCTURED_GRID
          unstructured_grid_ctype = 'implicit unstructured'
          discretization%itype = UNSTRUCTURED_GRID
          call InputReadNChars(input,option, &
                               discretization%filename, &
                               MAXSTRINGLENGTH, &
                               PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','filename')

          grid => GridCreate()
          un_str_sfgrid => UGridCreate()
          un_str_sfgrid%grid_type = TWO_DIM_GRID
          if (index(discretization%filename,'.h5') > 0) then
#if defined(PETSC_HAVE_HDF5)
            call UGridReadHDF5SurfGrid( un_str_sfgrid, &
                                        discretization%filename, &
                                        option)
#endif
          else
            call UGridReadSurfGrid(un_str_sfgrid, &
                                   surf_realization%subsurf_filename, &
                                   discretization%filename, &
                                   option)
          endif
          grid%unstructured_grid => un_str_sfgrid
          discretization%grid => grid
          grid%itype = unstructured_grid_itype
          grid%ctype = unstructured_grid_ctype

        case default
          option%io_buffer = 'Surface-flow supports only unstructured grid'
          call printErrMsg(option)
      end select
  end select

end subroutine SurfaceInit

! ************************************************************************** !

subroutine InitSurfaceSetupRealization(surf_realization,subsurf_realization, &
                                       waypoint_list)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Surface_Flow_module
  use Realization_Surface_class
  use Surface_TH_module
  use Surface_Global_module
  use Timestepper_Base_class
  use Realization_Subsurface_class
  
  use Option_module
  use Waypoint_module
  use Condition_Control_module
  use EOS_module
  
  implicit none
  
  class(realization_surface_type), pointer :: surf_realization
  class(realization_subsurface_type), pointer :: subsurf_realization
  type(waypoint_list_type) :: waypoint_list
  
  type(option_type), pointer :: option
  PetscErrorCode :: ierr
  
  option => surf_realization%option

  ! set reference densities if not specified in input file.
  call EOSReferenceDensity(option)

  call RealizSurfCreateDiscretization(surf_realization)

  ! Check if surface-flow is compatible with the given flowmode
  select case(option%iflowmode)
    case(TH_MODE)
    case default
      option%io_buffer = 'For surface-flow only TH mode implemented'
      call printErrMsgByRank(option)
  end select

  call SurfaceInitReadRegionFiles(surf_realization)
  call RealizSurfMapSurfSubsurfGrids(subsurf_realization,surf_realization)
  call RealizSurfLocalizeRegions(surf_realization)
  call RealizSurfPassFieldPtrToPatches(surf_realization)
  call RealizSurfProcessMatProp(surf_realization)
  call RealizSurfProcessCouplers(surf_realization)
  call RealizSurfProcessConditions(surf_realization)
  !call RealProcessFluidProperties(surf_realization)
  call SurfaceInitMatPropToRegions(surf_realization)
  call RealizSurfInitAllCouplerAuxVars(surf_realization)
  !call SurfaceRealizationPrintCouplers(surf_realization)

  ! add waypoints associated with boundary conditions, source/sinks etc. to list
  call RealizSurfAddWaypointsToList(surf_realization,waypoint_list)

  select case(option%iflowmode)
    case(TH_MODE)
      call SurfaceTHSetup(surf_realization)
  end select

  call SurfaceGlobalSetup(surf_realization)
  ! initialize FLOW
  ! set up auxillary variable arrays

  ! assign initial conditionsRealizAssignFlowInitCond
  call CondControlAssignFlowInitCondSurface(surf_realization)

  ! override initial conditions if they are to be read from a file
  if (len_trim(option%surf_initialize_flow_filename) > 1) then
    option%io_buffer = 'For surface-flow initial conditions cannot be read from file'
    call printErrMsgByRank(option)
  endif
  
  select case(option%iflowmode)
    case(TH_MODE)
      call SurfaceTHUpdateAuxVars(surf_realization)
    case default
      option%io_buffer = 'For surface-flow only RICHARDS and TH mode implemented'
      call printErrMsgByRank(option)
  end select
  
end subroutine InitSurfaceSetupRealization

! ************************************************************************** !

subroutine InitSurfaceSetupSolvers(surf_realization,solver,final_time)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
#include "petsc/finclude/petscts.h"
  use petscts
  use Realization_Surface_class
  use Option_module
  
  use Solver_module
  use Convergence_module
  use Discretization_module
  use Surface_Flow_module
  use Surface_TH_module
  
  implicit none
  
  class(realization_surface_type) :: surf_realization
  type(solver_type), pointer :: solver
  PetscReal :: final_time
  
  type(option_type), pointer :: option
  type(convergence_context_type), pointer :: convergence_context
  SNESLineSearch :: linesearch
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  
  option => surf_realization%option
  
  call printMsg(option,"  Beginning setup of FLOW SNES ")

  ! Setup PETSc TS for explicit surface flow solution
  call printMsg(option,"  Beginning setup of SURF FLOW TS ")

  call SolverCreateTS(solver,option%mycomm)
  call TSSetProblemType(solver%ts,TS_NONLINEAR, &
                        ierr);CHKERRQ(ierr)
  call TSSetDuration(solver%ts,ONE_INTEGER,final_time,ierr);CHKERRQ(ierr)
  
end subroutine InitSurfaceSetupSolvers

! ************************************************************************** !

subroutine SurfaceInitMatPropToRegions(surf_realization)
  ! 
  ! This routine assigns surface material properties to associated regions in
  ! the model (similar to assignMaterialPropToRegions)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/13/12
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Surface_class
  use Discretization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Surface_Field_module
  use Surface_Material_module
  
  use HDF5_module

  implicit none

  class(realization_surface_type) :: surf_realization
  
  PetscReal, pointer :: man0_p(:)
  PetscReal, pointer :: vec_p(:)
  
  PetscInt :: icell, local_id, ghosted_id, natural_id, surf_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(surface_field_type), pointer :: surf_field
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: patch  
  type(patch_type), pointer :: cur_patch

  type(surface_material_property_type), pointer :: surf_material_property
  type(surface_material_property_type), pointer :: null_surf_material_property
  type(region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  
  option => surf_realization%option
  discretization => surf_realization%discretization
  surf_field => surf_realization%surf_field

  ! loop over all patches and allocation material id arrays
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    if (.not.associated(cur_patch%imat)) then
      allocate(cur_patch%imat(cur_patch%grid%ngmax))
      ! initialize to "unset"
      cur_patch%imat = UNINITIALIZED_INTEGER
      ! also allocate saturation function id
      allocate(cur_patch%sat_func_id(cur_patch%grid%ngmax))
      cur_patch%sat_func_id = UNINITIALIZED_INTEGER
    endif
    cur_patch => cur_patch%next
  enddo

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    grid => cur_patch%grid
    strata => cur_patch%strata_list%first
    do
      if (.not.associated(strata)) exit
      ! Read in cell by cell material ids if they exist
      if (.not.associated(strata%region) .and. strata%active) then
        option%io_buffer = 'Reading of material prop from file for' // &
          ' surface flow is not implemented.'
        call printErrMsgByRank(option)
        !call readMaterialsFromFile(realization,strata%realization_dependent, &
        !                           strata%material_property_filename)
      ! Otherwise, set based on region
      else if (strata%active) then
        update_ghosted_material_ids = PETSC_TRUE
        region => strata%region
        surf_material_property => strata%surf_material_property
        if (associated(region)) then
          istart = 1
          iend = region%num_cells
        else
          istart = 1
          iend = grid%nlmax
        endif
        do icell=istart, iend
          if (associated(region)) then
            local_id = region%cell_ids(icell)
          else
            local_id = icell
          endif
          ghosted_id = grid%nL2G(local_id)
          cur_patch%imat(ghosted_id) = surf_material_property%internal_id
        enddo
      endif
      strata => strata%next
    enddo
    cur_patch => cur_patch%next
  enddo

  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call RealizSurfLocalToLocalWithArray(surf_realization,MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_surf_material_property => SurfaceMaterialPropertyCreate()
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    call VecGetArrayF90(surf_field%mannings0,man0_p,ierr);CHKERRQ(ierr)

    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      surf_material_id = cur_patch%imat(ghosted_id)
      if (surf_material_id == 0) then ! accomodate inactive cells
        surf_material_property = null_surf_material_property
      else if ( surf_material_id > 0 .and. &
                surf_material_id <= &
                size(surf_realization%surf_material_property_array)) then
        surf_material_property => &
          surf_realization%surf_material_property_array(surf_material_id)%ptr
        if (.not.associated(surf_material_property)) then
          write(dataset_name,*) surf_material_id
          option%io_buffer = 'No material property for surface material id ' // &
                              trim(adjustl(dataset_name)) &
                              //  ' defined in input file.'
          call printErrMsgByRank(option)
        endif
      else if (Uninitialized(surf_material_id)) then 
        write(dataset_name,*) grid%nG2A(ghosted_id)
        option%io_buffer = 'Uninitialized surface material id in patch at cell ' // &
                            trim(adjustl(dataset_name))
        call printErrMsgByRank(option)
      else if (surf_material_id > size(surf_realization%surf_material_property_array)) then
        write(option%io_buffer,*) surf_material_id
        option%io_buffer = 'Unmatched surface material id in patch:' // &
          adjustl(trim(option%io_buffer))
        call printErrMsgByRank(option)
      else
        option%io_buffer = 'Something messed up with surface material ids. ' // &
          ' Possibly material ids not assigned to all grid cells. ' // &
          ' Contact Glenn!'
        call printErrMsgByRank(option)
      endif
      man0_p(local_id) = surf_material_property%mannings
    enddo ! local_id - loop

    call VecRestoreArrayF90(surf_field%mannings0,man0_p,ierr);CHKERRQ(ierr)
      
    cur_patch => cur_patch%next
  enddo ! looping over patches
  
  call SurfaceMaterialPropertyDestroy(null_surf_material_property)
  nullify(null_surf_material_property)

  call DiscretizationGlobalToLocal(discretization,surf_field%mannings0, &
                                   surf_field%mannings_loc,ONEDOF)

end subroutine SurfaceInitMatPropToRegions

! ************************************************************************** !

subroutine SurfaceInitReadRegionFiles(surf_realization)
  ! 
  ! This routine reads surface region files
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/20/12
  ! 

  use Realization_Surface_class
  use Region_module
  use HDF5_module
  use Grid_module
  use Option_module

  implicit none

  class(realization_surface_type) :: surf_realization
  
  type(option_type), pointer :: option
  type(region_type), pointer :: surf_region
  PetscBool :: cell_ids_exists
  PetscBool :: face_ids_exists
  PetscBool :: vert_ids_exists

  option => surf_realization%option
  surf_region => surf_realization%surf_regions%first
  do 
    if (.not.associated(surf_region)) exit
    if (len_trim(surf_region%filename) > 1) then
      if (index(surf_region%filename,'.h5') > 0) then
#if defined(PETSC_HAVE_HDF5)
        call HDF5QueryRegionDefinition(surf_region, surf_region%filename, &
                                       surf_realization%option, &
                                       cell_ids_exists, face_ids_exists, &
                                       vert_ids_exists)
        if ( (.not. cell_ids_exists) .and. &
             (.not. face_ids_exists) .and. &
             (.not. vert_ids_exists)) then
          option%io_buffer = '"Regions/' // trim(surf_region%name) // &
              ' is not defined by "Cell Ids" or "Face Ids" or "Vertex Ids".'
          call printErrMsg(option)
        end if
        if (cell_ids_exists .or. face_ids_exists) then
          call HDF5ReadRegionFromFile(surf_realization%patch%grid, &
                                   surf_region, surf_region%filename, option)
        else
          call HDF5ReadRegionDefinedByVertex(option,surf_region, &
                                             surf_region%filename)
        endif
#endif      
      else if (index(surf_region%filename,'.ss') > 0) then
        surf_region%sideset => RegionCreateSideset()
        call RegionReadFromFile(surf_region%sideset,surf_region%filename, &
                                surf_realization%option)
      else
        call RegionReadFromFile(surf_region,surf_realization%option, &
                                surf_region%filename)
      endif
    endif
    surf_region => surf_region%next
  enddo

end subroutine SurfaceInitReadRegionFiles


end module Init_Surface_module
