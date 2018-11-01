module Factory_Surf_Subsurf_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Surf_Subsurf_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: SurfSubsurfaceInitialize, &
            SurfSubsurfaceReadFlowPM

contains

! ************************************************************************** !

subroutine SurfSubsurfaceInitialize(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 
  
  implicit none
  
  class(simulation_surfsubsurface_type) :: simulation

  ! NOTE: PETSc must already have been initialized here!
  call SurfSubsurfaceInitializePostPETSc(simulation)
  
end subroutine SurfSubsurfaceInitialize

! ************************************************************************** !

subroutine SurfSubsurfaceInitializePostPETSc(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Simulation_Surface_class
  use Simulation_Subsurface_class
  use Factory_Surface_module
  use Factory_Subsurface_module
  use Option_module
  use Init_Common_module
  use Init_Surface_module
  use Surface_Flow_module
  use Surface_TH_module
  use Simulation_Aux_module
  use PMC_Base_class
  use PMC_Surface_class
  use PFLOTRAN_Constants_module
  use PM_Base_class
  use PM_Base_Pointer_module
  use PM_Surface_class
  use PM_Surface_Flow_class
  use PM_Surface_TH_class
  use Input_Aux_module
  use Realization_Subsurface_class
  use String_module
  use Waypoint_module
  use Realization_Surface_class
  use Timestepper_Surface_class
  use Logging_module
  use Output_Aux_module
  
  implicit none

  class(simulation_surfsubsurface_type) :: simulation
  
  type(option_type), pointer :: option
  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_surface_type), pointer :: surf_realization
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pm_surface_flow_type), pointer :: pm_surface_flow
  class(pm_surface_th_type), pointer :: pm_surface_th
  class(pm_base_type), pointer :: cur_pm, prev_pm
  class(pmc_surface_type), pointer :: pmc_surface
  class(timestepper_surface_type), pointer :: timestepper
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: init_status
  VecScatter :: vscat_surf_to_subsurf
  VecScatter :: vscat_subsurf_to_surf
  Vec :: vec_subsurf_pres
  Vec :: vec_subsurf_pres_top_bc
  Vec :: vec_surf_head
  type(waypoint_type), pointer :: waypoint
  PetscErrorCode :: ierr
  
  option => simulation%option
  ! process command line arguments specific to subsurface
  !call SurfSubsurfInitCommandLineSettings(option)
  
  ! we need to remove the surface pm from the list while leaving the
  ! the subsurface pm linkage intact
  nullify(prev_pm)
  nullify(pm_surface_flow)
  nullify(pm_surface_th)

  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_surface_th_type)
        pm_surface_th => cur_pm
        if (associated(prev_pm)) then
          prev_pm%next => cur_pm%next
        else
          simulation%process_model_list => cur_pm%next
        endif
      class is(pm_surface_flow_type)
        pm_surface_flow => cur_pm
        if (associated(prev_pm)) then
          prev_pm%next => cur_pm%next
        else
          simulation%process_model_list => cur_pm%next
        endif
      class default
    end select
    prev_pm => cur_pm
    cur_pm => cur_pm%next
  enddo
  call SubsurfaceInitializePostPetsc(simulation)
  ! in SubsurfaceInitializePostPetsc, the first pmc in the list is set as
  ! the master, we need to negate this setting
  simulation%process_model_coupler_list%is_master = PETSC_FALSE

  if (associated(pm_surface_flow) .or. associated(pm_surface_th)) then
    simulation%surf_realization => RealizSurfCreate(option)
    surf_realization => simulation%surf_realization
    surf_realization%output_option => OutputOptionDuplicate(simulation%output_option)
    nullify(surf_realization%output_option%output_snap_variable_list)
    nullify(surf_realization%output_option%output_obs_variable_list)
    surf_realization%output_option%output_snap_variable_list => OutputVariableListCreate()
    surf_realization%output_option%output_obs_variable_list => OutputVariableListCreate()
    subsurf_realization => simulation%realization
    input => InputCreate(IN_UNIT,option%input_filename,option)
    surf_realization%subsurf_filename = &
      subsurf_realization%discretization%filename
    call SurfaceInitReadRequiredCards(simulation%surf_realization,input)
  
    call setSurfaceFlowMode(option)
  
    if (associated(pm_surface_flow)) then
      pm_surface_flow%output_option => simulation%output_option
      pmc_surface => PMCSurfaceCreate()
      pmc_surface%name = 'PMCSurface'
      simulation%surf_flow_process_model_coupler => pmc_surface
      pmc_surface%option => option
      pmc_surface%checkpoint_option => simulation%checkpoint_option
      pmc_surface%waypoint_list => simulation%waypoint_list_surfsubsurface
      pmc_surface%pm_list => pm_surface_flow
      pmc_surface%pm_ptr%pm => pm_surface_flow
      pmc_surface%surf_realization => simulation%surf_realization
      pmc_surface%subsurf_realization => simulation%realization
      timestepper => TimestepperSurfaceCreate()
      pmc_surface%timestepper => timestepper
      ! set up logging stage
      string = trim(pm_surface_flow%name) // 'Surface'
      call LoggingCreateStage(string,pmc_surface%stage)
    endif

    if (associated(pm_surface_th)) then
      pm_surface_th%output_option => simulation%output_option
      pmc_surface => PMCSurfaceCreate()
      pmc_surface%name = 'PMCSurface'
      simulation%surf_flow_process_model_coupler => pmc_surface
      pmc_surface%option => option
      pmc_surface%checkpoint_option => simulation%checkpoint_option
      pmc_surface%waypoint_list => simulation%waypoint_list_surfsubsurface
      pmc_surface%pm_list => pm_surface_th
      pmc_surface%pm_ptr%pm => pm_surface_th
      pmc_surface%surf_realization => simulation%surf_realization
      pmc_surface%subsurf_realization => simulation%realization
      timestepper => TimestepperSurfaceCreate()
      pmc_surface%timestepper => timestepper
      ! set up logging stage
      string = trim(pm_surface_th%name) // 'Surface'
      call LoggingCreateStage(string,pmc_surface%stage)
    endif
    
    input => InputCreate(IN_UNIT,option%input_filename,option)    
    string = 'SURFACE_FLOW'
    call InputFindStringInFile(input,option,string)
    call InputFindStringErrorMsg(input,option,string)  
    call SurfaceReadInput(surf_realization,timestepper%solver, &
                          simulation%waypoint_list_surfsubsurface,input)
    ! Add first waypoint
    waypoint => WaypointCreate()
    waypoint%time = 0.d0
    call WaypointInsertInList(waypoint,simulation%waypoint_list_surfsubsurface)
  
    ! Add final_time waypoint to surface_realization
    waypoint => WaypointCreate()
    waypoint%final = PETSC_TRUE
    waypoint%time = simulation%waypoint_list_subsurface%last%time
    waypoint%print_snap_output = PETSC_TRUE
    call WaypointInsertInList(waypoint,simulation%waypoint_list_surfsubsurface)   
    ! merge in outer waypoints (e.g. checkpoint times)
    call WaypointListCopyAndMerge(simulation%waypoint_list_surfsubsurface, &
                                  simulation%waypoint_list_subsurface,option)
    call WaypointListCopyAndMerge(simulation%waypoint_list_surfsubsurface, &
                                  simulation%waypoint_list_outer,option)
    call InitSurfaceSetupRealization(surf_realization,subsurf_realization, &
                                     simulation%waypoint_list_surfsubsurface)
    call InitCommonAddOutputWaypoints(option,simulation%output_option, &
                                      simulation%waypoint_list_surfsubsurface)
    ! fill in holes in waypoint data
    call WaypointListFillIn(simulation%waypoint_list_surfsubsurface,option)
    call WaypointListRemoveExtraWaypnts(simulation% &
                                            waypoint_list_surfsubsurface,option)
      
    if (associated(simulation%surf_flow_process_model_coupler)) then
      if (associated(simulation%surf_flow_process_model_coupler% &
                     timestepper)) then
        simulation%surf_flow_process_model_coupler%timestepper%cur_waypoint => &
          simulation%waypoint_list_surfsubsurface%first
      endif
    endif
    call InitSurfaceSetupSolvers(surf_realization,timestepper%solver, &
                             simulation%waypoint_list_surfsubsurface%last%time)

    pmc_surface%timestepper%dt_init = surf_realization%dt_init
    pmc_surface%timestepper%dt_max = surf_realization%dt_max
    option%surf_subsurf_coupling_flow_dt = surf_realization%dt_coupling
    option%surf_flow_dt=pmc_surface%timestepper%dt_init

    if (associated(pm_surface_flow)) then
      call pm_surface_flow%PMSurfaceSetRealization(surf_realization)
      call pm_surface_flow%Setup()
      call TSSetRHSFunction(timestepper%solver%ts, &
                            pm_surface_flow%residual_vec, &
                            PMRHSFunction, &
                            pmc_surface%pm_ptr, &
                            ierr);CHKERRQ(ierr)
    endif

    if (associated(pm_surface_th)) then
      call pm_surface_th%PMSurfaceSetRealization(surf_realization)
      call pm_surface_th%Setup()
      call TSSetRHSFunction(timestepper%solver%ts, &
                            pm_surface_th%residual_vec, &
                            PMRHSFunction, &
                            pmc_surface%pm_ptr, &
                            ierr);CHKERRQ(ierr)
    endif

    timestepper%dt = option%surf_flow_dt
    
    nullify(simulation%process_model_coupler_list)
  endif

  ! sim_aux: Create PETSc Vectors and VectorScatters
  if (option%surf_flow_on .and. option%subsurf_surf_coupling /= DECOUPLED) then

    call SurfSubsurfCreateSubsurfVecs(subsurf_realization, option, &
                                      vec_subsurf_pres, vec_subsurf_pres_top_bc)
    call SimAuxCopySubsurfVec(simulation%sim_aux, vec_subsurf_pres)
    call SimAuxCopySubsurfTopBCVec(simulation%sim_aux, vec_subsurf_pres_top_bc)
    call VecDestroy(vec_subsurf_pres, ierr);CHKERRQ(ierr)
    call VecDestroy(vec_subsurf_pres_top_bc, ierr);CHKERRQ(ierr)

    call SurfSubsurfCreateSurfVecs(surf_realization, option, &
                                   vec_surf_head)
    call SimAuxCopySurfVec(simulation%sim_aux, vec_surf_head)
    call VecDestroy(vec_surf_head, ierr);CHKERRQ(ierr)

    call SurfSubsurfCreateSurfSubSurfVScats(subsurf_realization, &
          surf_realization, vscat_surf_to_subsurf, vscat_subsurf_to_surf)
    call SimAuxCopyVecScatter(simulation%sim_aux, vscat_surf_to_subsurf, SURF_TO_SUBSURF)
    call SimAuxCopyVecScatter(simulation%sim_aux, vscat_subsurf_to_surf, SUBSURF_TO_SURF)
    call VecScatterDestroy(vscat_surf_to_subsurf, ierr);CHKERRQ(ierr)
    call VecScatterDestroy(vscat_surf_to_subsurf, ierr);CHKERRQ(ierr)
  endif

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if (associated(simulation%rt_process_model_coupler)) &
    simulation%rt_process_model_coupler%sim_aux => simulation%sim_aux
  if (option%surf_flow_on .and. &
     associated(simulation%surf_flow_process_model_coupler)) &
    simulation%surf_flow_process_model_coupler%sim_aux => simulation%sim_aux

  ! set surface flow as master
  simulation%surf_flow_process_model_coupler%is_master = PETSC_TRUE
  ! link surface flow and master
  simulation%process_model_coupler_list => &
    simulation%surf_flow_process_model_coupler
  ! link subsurface flow as peer
  simulation%process_model_coupler_list%peer => &
    simulation%flow_process_model_coupler

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  if (associated(cur_process_model_coupler%peer)) then
    cur_process_model_coupler => cur_process_model_coupler%peer
    call cur_process_model_coupler%SetAuxData()
  endif

end subroutine SurfSubsurfaceInitializePostPETSc

! ************************************************************************** !

subroutine SurfSubsurfaceReadFlowPM(input, option, pm)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 11/29/15
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PMC_Base_class
  use PM_Base_class
  use PM_Surface_Flow_class
  use PM_Surface_TH_class
  use Init_Common_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,SURFACE_SUBSURFACE'

  option%surf_flow_on = PETSC_TRUE

  nullify(pm)
  word = ''
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('MODE')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)
        select case(word)
          case('RICHARDS')
            pm => PMSurfaceFlowCreate()
          case('TH')
            pm => PMSurfaceTHCreate()
          case default
            error_string = trim(error_string) // ',MODE'
            call InputKeywordUnrecognized(word,error_string,option)
        end select
        pm%option => option
        exit
      case default
        error_string = trim(error_string) // ',SURFACE_FLOW'
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo

  if (.not.associated(pm)) then
    option%io_buffer = 'A flow MODE (card) must be included in the ' // &
      'SURFACE_SUBSURFACE block in ' // trim(error_string) // '.'
    call printErrMsg(option)
  endif

end subroutine SurfSubsurfaceReadFlowPM

! ************************************************************************** !

subroutine SurfSubsurfCreateSurfSubSurfVScats(realization, surf_realization, &
                                              surf_to_subsurf, subsurf_to_surf)
  ! 
  ! This routine creates VecScatter between surface-subsurface grids.
  ! Algorithm:
  ! - It uses a similar logic of Matrix-Vector multiplication used in
  ! UGridMapSideSet() subroutine. The algorithm here is extended to use
  ! Matrix-Matrix mulitplication
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
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
  use Realization_Surface_class

  implicit none

  class(realization_subsurface_type),pointer :: realization
  class(realization_surface_type),pointer :: surf_realization
  VecScatter :: surf_to_subsurf
  VecScatter :: subsurf_to_surf

  type(option_type),pointer :: option
  type(grid_unstructured_type),pointer :: subsurf_grid
  type(grid_unstructured_type),pointer :: surf_grid
  type(patch_type),pointer :: cur_patch
  type(region_type),pointer :: cur_region,top_region
  type(region_type),pointer :: patch_region

  Mat :: Mat_vert_to_face_subsurf
  Mat :: Mat_vert_to_face_subsurf_transp
  Mat :: Mat_vert_to_face_surf
  Mat :: Mat_vert_to_face_surf_transp
  Mat :: prod
  Vec :: subsurf_petsc_ids,surf_petsc_ids

  PetscViewer :: viewer

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer ::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4)
  PetscInt :: nvertices
  PetscInt :: iface
  PetscInt :: local_id,ii,jj
  PetscInt :: cell_type
  PetscInt :: ivertex,vertex_id_local
  PetscReal :: real_array4(4)
  PetscReal,pointer :: vec_ptr(:)

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
  call MatZeroEntries(Mat_vert_to_face_subsurf,ierr);CHKERRQ(ierr)
  real_array4 = 1.d0

  call VecGetArrayF90(subsurf_petsc_ids,vec_ptr,ierr);CHKERRQ(ierr)

  offset=0
  call MPI_Exscan(top_region%num_cells,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  do ii = 1,top_region%num_cells
    local_id = top_region%cell_ids(ii)
    vec_ptr(ii) = subsurf_grid%cell_ids_petsc(local_id)
    iface    = top_region%faces(ii)
    cell_type = subsurf_grid%cell_type(local_id)
    !nfaces = UCellGetNFaces(cell_type,option)

    call UCellGetNFaceVertsandVerts(option,cell_type,iface,nvertices, &
                                    int_array4)
    ! For this matrix:
    !   irow = local face id
    !   icol = natural (global) vertex id
    do ivertex = 1,nvertices
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
  offset=0
  call MPI_Exscan(surf_grid%nlmax,offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  call VecGetArrayF90(surf_petsc_ids,vec_ptr,ierr);CHKERRQ(ierr)

  do local_id = 1,surf_grid%nlmax
    cell_type = surf_grid%cell_type(local_id)
    vec_ptr(local_id) = surf_grid%cell_ids_petsc(local_id)

    int_array4_0 = 0
    nvertices = surf_grid%cell_vertices(0,local_id)
    do ivertex = 1,nvertices
      vertex_id_local = surf_grid%cell_vertices(ivertex,local_id)
      int_array4_0(ivertex) = &
        surf_grid%vertex_ids_natural(vertex_id_local)-1
    enddo
    call MatSetValues(Mat_vert_to_face_surf, &
                      1,local_id-1+offset, &
                      nvertices,int_array4_0, &
                      real_array4, &
                      INSERT_VALUES,ierr);CHKERRQ(ierr)
    call MatSetValues(Mat_vert_to_face_surf_transp, &
                      nvertices,int_array4_0, &
                      1,local_id-1+offset, &
                      real_array4, &
                      INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo

  call VecRestoreArrayF90(surf_petsc_ids,vec_ptr,ierr);CHKERRQ(ierr)

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

  call SurfSubsurfCreateSurfSubSurfVScat(realization,surf_realization,prod, &
                                     surf_petsc_ids,surf_to_subsurf)
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
  call SurfSubsurfCreateSurfSubSurfVScat(realization,surf_realization,prod, &
                                        subsurf_petsc_ids,subsurf_to_surf)

  call MatDestroy(prod,ierr);CHKERRQ(ierr)

  call MatDestroy(Mat_vert_to_face_subsurf,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_vert_to_face_subsurf_transp,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_vert_to_face_surf,ierr);CHKERRQ(ierr)
  call MatDestroy(Mat_vert_to_face_surf_transp,ierr);CHKERRQ(ierr)

  call VecDestroy(subsurf_petsc_ids,ierr);CHKERRQ(ierr)
  call VecDestroy(surf_petsc_ids,ierr);CHKERRQ(ierr)

end subroutine SurfSubsurfCreateSurfSubSurfVScats

! ************************************************************************** !

subroutine SurfSubsurfCreateSurfSubSurfVScat( &
              realization, &       !<
              surf_realization, &  !<
              prod_mat, &          !< Mat-Mat-Mult matrix
              source_petsc_ids, &   !< MPI-Vector containing cell ids in PETSc order
              scatter &
              )
  ! 
  ! This subroutine creates a single vector scatter context
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
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
  use Realization_Surface_class

  implicit none

  class(realization_subsurface_type),pointer :: realization
  class(realization_surface_type),pointer :: surf_realization
  Mat :: prod_mat
  Vec :: source_petsc_ids
  VecScatter :: scatter

  Mat :: prod_loc_mat
  Vec :: source_loc_vec
  Vec :: corr_dest_ids_vec
  Vec :: corr_dest_ids_vec_ndof
  Vec :: source_petsc_ids_ndof
  IS :: is_tmp1,is_tmp2
  IS :: is_tmp3,is_tmp4
  PetscInt,pointer :: corr_v2_ids(:)
  VecScatter :: scatter_ndof

  PetscViewer :: viewer

  type(option_type),pointer :: option
  type(field_type),pointer :: field
  type(surface_field_type),pointer :: surf_field

  type(dm_ptr_type),pointer :: dm_ptr
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt,pointer ::int_array(:)
  PetscInt :: offset
  PetscInt :: int_array4(4)
  PetscInt :: int_array4_0(4)
  PetscReal :: real_array4(4)
  PetscInt :: ii,jj
  PetscReal,pointer :: vec_ptr(:)
  PetscInt :: ivertex,cell_id,vertex_id_local
  PetscReal :: max_value

  PetscInt,pointer :: ia_p(:),ja_p(:)
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
    call MatGetRowIJF90(prod_loc_mat,ONE_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                        nrow,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArrayF90(prod_loc_mat,aa_v,ierr);CHKERRQ(ierr)
  else
    ! Get i and j indices of the local-matrix
    call MatGetRowIJF90(prod_mat,ONE_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                        nrow,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
    ! Get values stored in the local-matrix
    call MatSeqAIJGetArrayF90(prod_mat,aa_v,ierr);CHKERRQ(ierr)
  endif

  ! For each row of the local-matrix,find the column with the largest value
  allocate(corr_v2_ids(nrow))
  row = 1
  col = 0
  do ii = 1,nrow
    max_value = 0.d0
    do jj = ia_p(ii),ia_p(ii + 1) - 1
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
  do ii = 1,nrow
    int_array(ii) = (ii-1)+offset
  enddo
  call ISCreateGeneral(option%mycomm,nrow, &
                       int_array,PETSC_COPY_VALUES,is_tmp1,ierr);CHKERRQ(ierr)
  call ISCreateBlock(option%mycomm,option%nflowdof,nrow, &
                     int_array,PETSC_COPY_VALUES,is_tmp3,ierr);CHKERRQ(ierr)

  do ii = 1,nrow
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

  call VecDestroy(corr_dest_ids_vec,ierr);CHKERRQ(ierr)
  if (option%mycommsize>1) then
    call MatSeqAIJRestoreArrayF90(prod_loc_mat,aa_v,ierr);CHKERRQ(ierr)
    call MatDestroy(prod_loc_mat,ierr);CHKERRQ(ierr)
  else
    call MatSeqAIJRestoreArrayF90(prod_mat,aa_v,ierr);CHKERRQ(ierr)
  endif

end subroutine SurfSubsurfCreateSurfSubSurfVScat

! ************************************************************************** !

subroutine SurfSubsurfCreateSubsurfVecs(subsurf_realization, option, &
                                        subsurf_pres, subsurf_pres_top_bc)
  ! 
  ! This routine creates subsurface vectors.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Coupler_module
  use Option_module
  use String_module

  implicit none

  class(realization_subsurface_type),pointer :: subsurf_realization
  type(option_type),pointer :: option
  Vec :: subsurf_pres
  Vec :: subsurf_pres_top_bc

  type(coupler_list_type),pointer :: coupler_list
  type(coupler_type),pointer :: coupler

  PetscInt :: num_conn
  PetscInt :: found
  PetscInt :: found_global
  PetscErrorCode :: ierr

  call VecCreate(option%mycomm,subsurf_pres,ierr);CHKERRQ(ierr)
  call VecSetSizes(subsurf_pres, &
                   subsurf_realization%discretization%grid%nlmax, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(subsurf_pres,ierr);CHKERRQ(ierr)
  call VecSet(subsurf_pres,0.d0,ierr);CHKERRQ(ierr)

#if 1
  found = 0
  num_conn = 0
  coupler_list => subsurf_realization%patch%source_sink_list
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit
    if (StringCompare(coupler%name,'from_surface_ss')) then
      num_conn = coupler%connection_set%num_connections
      found = 1
    endif
    coupler => coupler%next
  enddo

  call MPI_AllReduce(found,found_global,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                     option%mycomm,ierr)
#endif

  !found_global = 0
  if (found_global == 0) then
    coupler_list => subsurf_realization%patch%boundary_condition_list
    coupler => coupler_list%first
    do
      if (.not.associated(coupler)) exit
      if (StringCompare(coupler%name,'from_surface_bc')) then
        num_conn = coupler%connection_set%num_connections
        found = 1
      endif
      coupler => coupler%next
    enddo
    call MPI_AllReduce(found,found_global,ONE_INTEGER_MPI,MPI_INTEGER,MPI_MAX, &
                       option%mycomm,ierr)
  endif

  if (found_global > 0) then
    call VecCreate(option%mycomm,subsurf_pres_top_bc,ierr);CHKERRQ(ierr)
    call VecSetSizes(subsurf_pres_top_bc,num_conn,PETSC_DECIDE, &
                     ierr);CHKERRQ(ierr)
    call VecSetFromOptions(subsurf_pres_top_bc,ierr);CHKERRQ(ierr)
    call VecSet(subsurf_pres_top_bc,0.d0,ierr);CHKERRQ(ierr)
  endif

end subroutine SurfSubsurfCreateSubsurfVecs

! ************************************************************************** !

subroutine SurfSubsurfCreateSurfVecs(surf_realization,option,surf_head)
  ! 
  ! This routine creates surface vectors.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Surface_class
  use Option_module

  implicit none

  class(realization_surface_type),pointer :: surf_realization
  type(option_type),pointer :: option
  Vec :: surf_head

  PetscErrorCode :: ierr

  call VecCreate(option%mycomm,surf_head,ierr);CHKERRQ(ierr)
  call VecSetSizes(surf_head,surf_realization%discretization%grid%nlmax, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetFromOptions(surf_head,ierr);CHKERRQ(ierr)
  call VecSet(surf_head,0.d0,ierr);CHKERRQ(ierr)

end subroutine SurfSubsurfCreateSurfVecs

! ************************************************************************** !

subroutine SurfSubsurfInitCommandLineSettings(option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/01/13
  ! 

  use Option_module
  use Input_Aux_module
  
  implicit none
  
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscBool :: bool_flag
  
  string = '-multisimulation'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = MULTISIMULATION_SIM_TYPE
  endif

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) then
    option%subsurface_simulation_type = STOCHASTIC_SIM_TYPE
  endif
  
end subroutine SurfSubsurfInitCommandLineSettings

end module Factory_Surf_Subsurf_module
