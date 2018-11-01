#ifdef CLM_PFLOTRAN

module pflotran_clm_setmapping_module

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petsclog.h"
#include "petsc/finclude/petscviewer.h"
  use petscsys
  use petscvec

  use Option_module, only : option_type
  use Simulation_Base_class, only : simulation_base_type
  use Realization_Base_class, only : realization_base_type
  use PFLOTRAN_Constants_module
  use Utility_module, only : where_checkerr

  use Mapping_module

  implicit none

  private

  public ::                                  &
       ! mesh-mapping
       pflotranModelSetupMappingFiles,       &
       pflotranModelInitMapping,             &
       pflotranModelGetTopFaceArea

  private ::                                 &
       pflotranModelInitMappingSub2Sub,      &
       pflotranModelInitMapFaceToFace

  ! mesh ids
  PetscInt, parameter :: CLM_3DSUB_MESH = 1
  PetscInt, parameter :: PF_3DSUB_MESH  = 2
  PetscInt, parameter :: CLM_FACE_MESH  = 3
  PetscInt, parameter :: PF_FACE_MESH   = 4

  !PetscInt, parameter :: CLM_2DTOP_MESH = 5
  !PetscInt, parameter :: PF_2DTOP_MESH  = 6

  character(len=MAXWORDLENGTH) :: subname = ""

contains

! ************************************************************************** !

! ************************************************************************** !

  subroutine pflotranModelSetupMappingFiles(model)
  ! 
  ! pflotranModelSetupMappingFiles
  ! create the mapping objects, reopen the input file and read the file names
  ! before model integration is performed by the call to StepperRun()
  ! routine
  ! NOTE(bja, 2013-07) this really needs to be moved out of pflotran
  ! CLM should be responsible for passing data in the correct
  ! format. That may require pflotran to provide a call back function
  ! for grid info.
  !
  ! Author: Gautam Bisht
  ! Date: 9/10/2010
  ! Revised by Fengming YUAN

    use String_module
    use Option_module
    use Input_Aux_module
    use pflotran_clm_main_module, only : pflotran_model_type

    !use Simulation_Subsurface_class, only : subsurface_simulation_type
    !use Realization_class, only : realization_type

    implicit none

    type(pflotran_model_type), pointer, intent(inout) :: model
    type(input_type), pointer :: input

    PetscBool :: clm2pf_3dsub_file
    PetscBool :: pf2clm_3dsub_file

    PetscBool :: clm2pf_bctop_file
    PetscBool :: pf2clm_bctop_file

    PetscBool :: clm2pf_bcbot_file
    PetscBool :: pf2clm_bcbot_file

    character(len=MAXSTRINGLENGTH) :: string
    character(len=MAXWORDLENGTH) :: word

    !nullify(model%pf_cells)
    !nullify(model%clm_cells)
    nullify(model%map_clm_sub_to_pf_sub)
    nullify(model%map_clm_2dtop_to_pf_2dtop)
    nullify(model%map_pf_sub_to_clm_sub)
    !
    nullify(model%map_clm_2dbot_to_pf_2dbot)
    nullify(model%map_pf_2dbot_to_clm_2dbot)
    nullify(model%map_pf_2dtop_to_clm_2dtop)

    model%nlclm = -1
    model%ngclm = -1

    input => InputCreate(IUNIT_TEMP, &
                    model%option%input_filename, model%option)

    ! Read names of mapping file
    clm2pf_3dsub_file=PETSC_FALSE
    pf2clm_3dsub_file=PETSC_FALSE
    clm2pf_bctop_file=PETSC_FALSE
    pf2clm_bctop_file=PETSC_FALSE
    clm2pf_bcbot_file=PETSC_FALSE
    pf2clm_bcbot_file=PETSC_FALSE
    
    string = "MAPPING_FILES"
    call InputFindStringInFile(input,model%option,string)

    do
      call InputReadPflotranString(input, model%option)
      if (InputCheckExit(input, model%option)) exit
      if (input%ierr /= 0) exit

      call InputReadWord(input, model%option, word, PETSC_TRUE)
      call InputErrorMsg(input, model%option, 'keyword', 'MAPPING_FILES')
      call StringToUpper(word)

      select case(trim(word))
        case('CLM2PF_SUB_FILE')
          model%map_clm_sub_to_pf_sub => MappingCreate()
          call InputReadNChars(input, model%option, model%map_clm_sub_to_pf_sub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm_sub_to_pf_sub%filename = &
            trim(model%map_clm_sub_to_pf_sub%filename)//CHAR(0)
          model%map_clm_sub_to_pf_sub%id = CLM_3DSUB_TO_PF_3DSUB
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')   
          clm2pf_3dsub_file=PETSC_TRUE
        case('PF2CLM_SUB_FILE')
          model%map_pf_sub_to_clm_sub => MappingCreate()
          call InputReadNChars(input, model%option, model%map_pf_sub_to_clm_sub%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_pf_sub_to_clm_sub%filename = &
            trim(model%map_pf_sub_to_clm_sub%filename)//CHAR(0)
          model%map_pf_sub_to_clm_sub%id = PF_3DSUB_TO_CLM_3DSUB
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_3dsub_file=PETSC_TRUE

        case('CLM2PF_BCTOP_FILE')
          model%map_clm_2dtop_to_pf_2dtop => MappingCreate()
          call InputReadNChars(input, model%option, model%map_clm_2dtop_to_pf_2dtop%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm_2dtop_to_pf_2dtop%filename = &
            trim(model%map_clm_2dtop_to_pf_2dtop%filename)//CHAR(0)
          model%map_clm_2dtop_to_pf_2dtop%id = CLM_2DTOP_TO_PF_2DTOP
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_bctop_file=PETSC_TRUE
        case('PF2CLM_BCTOP_FILE')
          model%map_pf_2dtop_to_clm_2dtop => MappingCreate()
          call InputReadNChars(input, model%option, model%map_pf_2dtop_to_clm_2dtop%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_pf_2dtop_to_clm_2dtop%filename = &
              trim(model%map_pf_2dtop_to_clm_2dtop%filename)//CHAR(0)
          model%map_pf_2dtop_to_clm_2dtop%id = PF_2DTOP_TO_CLM_2DTOP
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_bctop_file=PETSC_TRUE
        case('CLM2PF_BCBOT_FILE')
          model%map_clm_2dbot_to_pf_2dbot => MappingCreate()
          call InputReadNChars(input, model%option, model%map_clm_2dbot_to_pf_2dbot%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_clm_2dbot_to_pf_2dbot%filename = &
            trim(model%map_clm_2dbot_to_pf_2dbot%filename)//CHAR(0)
          model%map_clm_2dbot_to_pf_2dbot%id = CLM_2DBOT_TO_PF_2DBOT
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          clm2pf_bcbot_file=PETSC_TRUE
        case('PF2CLM_BCBOT_FILE')
          model%map_pf_2dbot_to_clm_2dbot => MappingCreate()
          call InputReadNChars(input, model%option, model%map_pf_2dbot_to_clm_2dbot%filename, &
               MAXSTRINGLENGTH, PETSC_TRUE)
          model%map_pf_2dbot_to_clm_2dbot%filename = &
            trim(model%map_pf_2dbot_to_clm_2dbot%filename)//CHAR(0)
          model%map_pf_2dbot_to_clm_2dbot%id = PF_2DBOT_TO_CLM_2DBOT
          call InputErrorMsg(input, &
                             model%option, 'type', 'MAPPING_FILES')
          pf2clm_bcbot_file=PETSC_TRUE
        case default
          model%option%io_buffer='Keyword ' // trim(word) // &
            ' in input file not recognized and ignored'
          call printMsg(model%option)
      end select

    enddo
    call InputDestroy(input)

    if ((.not. clm2pf_3dsub_file) .and. &
        (.not. pf2clm_3dsub_file) ) then
      model%option%io_buffer='One of the 3D soil-mesh mapping files not found - So,' // &
      ' CLM grids conversion to PF structured-cartesian grid USED!'
      call printMsg(model%option)

      if(.not.associated(model%map_clm_sub_to_pf_sub)) then
        model%map_clm_sub_to_pf_sub => MappingCreate()
        model%map_clm_sub_to_pf_sub%id = CLM_3DSUB_TO_PF_3DSUB
      endif

      if(.not.associated(model%map_pf_sub_to_clm_sub)) then
        model%map_pf_sub_to_clm_sub => MappingCreate()
        model%map_pf_sub_to_clm_sub%id = PF_3DSUB_TO_CLM_3DSUB
      endif

    endif
    
    if(model%option%iflowmode==TH_MODE) then
      if(.not.clm2pf_bctop_file .and. .not.pf2clm_bctop_file) then
         model%option%io_buffer='Running in TH_MODE ' // &
          ' without pair of top-cell mesh mapping files - SO, ' // &
          ' CLM grids conversion to PF structured-cartesian grid USED!'
        call printMsg(model%option)

        if(.not.associated(model%map_clm_2dtop_to_pf_2dtop)) then
          model%map_clm_2dtop_to_pf_2dtop => MappingCreate()
          model%map_clm_2dtop_to_pf_2dtop%id = CLM_2DTOP_TO_PF_2DTOP
        endif

        if(.not.associated(model%map_pf_2dtop_to_clm_2dtop)) then
          model%map_pf_2dtop_to_clm_2dtop => MappingCreate()
          model%map_pf_2dtop_to_clm_2dtop%id = PF_2DTOP_TO_CLM_2DTOP
        endif

      endif

      if(.not.clm2pf_bcbot_file .and. .not.pf2clm_bcbot_file) then
        model%option%io_buffer='Running in TH_MODE/Richards_MODE ' // &
          ' without pair of bottom-cell mesh mapping files - SO, ' // &
          ' CLM grids conversion to PF structured-cartesian grid USED!'
        call printMsg(model%option)

        if(.not.associated(model%map_clm_2dbot_to_pf_2dbot)) then
          model%map_clm_2dbot_to_pf_2dbot => MappingCreate()
          model%map_clm_2dbot_to_pf_2dbot%id = CLM_2DBOT_TO_PF_2DBOT
        endif

        if(.not.associated(model%map_pf_2dbot_to_clm_2dbot)) then
          model%map_pf_2dbot_to_clm_2dbot => MappingCreate()
          model%map_pf_2dbot_to_clm_2dbot%id = PF_2DBOT_TO_CLM_2DBOT
        endif

      endif

    endif

  end subroutine pflotranModelSetupMappingFiles

! ************************************************************************** !

  subroutine pflotranModelInitMapping(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
  ! 
  ! Initialize mapping between the two model grid
  ! (CLM and PFLTORAN)
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/24/2011
  ! Revised by Fengming YUAN

    use Option_module
    use String_module

    use pflotran_clm_main_module, only : pflotran_model_type

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id
    
    select case (map_id)
      case (CLM_3DSUB_TO_PF_3DSUB, PF_3DSUB_TO_CLM_3DSUB)
        call pflotranModelInitMappingSub2Sub(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
      ! A more generalized Mapping for Faces (sidesets) is now implemented (F.-M. Yuan)
      !case (CLM_2DTOP_TO_PF_2DTOP, PF_2DTOP_TO_CLM_2DTOP)
      !  call pflotranModelInitMapTopTo2DSub(pflotran_model,  &
      !                                      grid_clm_cell_ids_nindex, &
      !                                      grid_clm_npts_local, &
      !                                      map_id)
      case (CLM_2DTOP_TO_PF_2DTOP, PF_2DTOP_TO_CLM_2DTOP, &
            CLM_2DBOT_TO_PF_2DBOT, PF_2DBOT_TO_CLM_2DBOT)
        call pflotranModelInitMapFaceToFace(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)

      case default
        pflotran_model%option%io_buffer = 'Invalid map_id argument to ' // &
          'pflotranModelInitMapping'
        call printErrMsg(pflotran_model%option)
    end select

  end subroutine pflotranModelInitMapping

! ************************************************************************** !

  subroutine pflotranModelInitMappingSub2Sub(pflotran_model,  &
                                      grid_clm_cell_ids_nindex, &
                                      grid_clm_npts_local, &
                                      map_id)
  ! 
  ! Initialize mapping between 3D subsurface
  ! CLM grid and 3D subsurface PFLOTRAN grid.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/09/2013
  ! Revised by Fengming YUAN

    use Option_module
    use Grid_module
    use Patch_module
    use String_module
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_subsurface_class, only : realization_subsurface_type

    use pflotran_clm_main_module, only : pflotran_model_type
    use clm_pflotran_interface_data

    implicit none

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id

    ! local
    PetscInt                           :: local_id, ghosted_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_subsurface_type), pointer   :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch

    !-----------------------------------------------------------------------------------
    option          => pflotran_model%option
    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         option%io_buffer = "subroutine: " // trim(subname) //&
          " only works on subsurface simulations."
         call printErrMsg(option)
    end select
    patch           => realization%patch
    grid            => patch%grid

    !
    ! Mapping to/from entire PFLOTRAN 3D subsurface domain
    !

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM_3DSUB_TO_PF_3DSUB)
        map => pflotran_model%map_clm_sub_to_pf_sub
        source_mesh_id = CLM_3DSUB_MESH
        dest_mesh_id = PF_3DSUB_MESH
      case(PF_3DSUB_TO_CLM_3DSUB)
        map => pflotran_model%map_pf_sub_to_clm_sub
        source_mesh_id = PF_3DSUB_MESH
        dest_mesh_id = CLM_3DSUB_MESH
      case default
        option%io_buffer = 'Invalid map_id argument to pflotranModelInitMapping'
        call printErrMsg(option)
    end select

    if (len(trim(map%filename)) > 0) then
      ! Read mapping file
      if (index(map%filename, '.h5') > 0) then
        call MappingReadHDF5(map, map%filename, option)
      else
        call MappingReadTxtFile(map, map%filename, option)
      endif

      ! checking if the CLM-PF has same number of soil layers for mapping
      if (map%pflotran_nlev /= clm_pf_idata%nzclm_mapped) then
         option%io_buffer = 'Invalid mapping soil layers between CLM and PFLOTRAN!'
        call printErrMsg(option)
      end if

    elseif(.not.option%mapping_files) then
      ! directly mapping between CLM and PF meshes, if no user-defined mapping file
      map%pflotran_nlev_mapped = grid%structured_grid%nz
      map%clm_nlev_mapped = clm_pf_idata%nzclm_mapped
      if (dest_mesh_id == PF_3DSUB_MESH) then
        call MappingFromCLMGrids(map, grid, PETSC_TRUE, option)
      elseif(dest_mesh_id == CLM_3DSUB_MESH) then
        call MappingFromCLMGrids(map, grid, PETSC_FALSE, option)
      endif

    else

        option%io_buffer = 'MUST provide mapping files between CLM and PFLOTRAN! '
        call printErrMsg(option)

    endif

    grid_clm_npts_ghost=0
    ! Allocate memory to identify if CLM cells are local or ghosted.
    ! Note: Presently all CLM cells are local
    allocate(grid_clm_local_nindex(grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_nindex(local_id) = local_id     ! LOCAL ID
    enddo

    ! Find cell IDs for PFLOTRAN grid
    grid_pf_npts_local = grid%nlmax
    grid_pf_npts_ghost = grid%ngmax - grid%nlmax


    allocate(grid_pf_cell_ids_nindex(grid%ngmax))
    allocate(grid_pf_local_nindex(grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      local_id = grid%nG2L(ghosted_id)
      grid_pf_cell_ids_nindex(ghosted_id) = grid%nG2A(ghosted_id)-1
      if (local_id <= 0) then
        grid_pf_local_nindex(ghosted_id) = 0        ! GHOST
      else
        grid_pf_local_nindex(ghosted_id) = local_id        ! LOCAL ID
      endif
    enddo

    select case(source_mesh_id)
      case(CLM_3DSUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                                      grid_clm_npts_ghost, &
                                                      grid_clm_cell_ids_nindex, &
                                                      grid_clm_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                                           grid_pf_npts_ghost, &
                                                           grid_pf_cell_ids_nindex, &
                                                           grid_pf_local_nindex)
      case(PF_3DSUB_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                                      grid_pf_npts_ghost, &
                                                      grid_pf_cell_ids_nindex, &
                                                      grid_pf_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                                           grid_clm_npts_ghost, &
                                                           grid_clm_cell_ids_nindex, &
                                                           grid_clm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to pflotranModelInitMapping'
        call printErrMsg(option)
    end select

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)
    call MappingFreeNotNeeded(map)

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_clm_local_nindex)

    ! Setting the number of cells constituting the 3D
    ! subsurface domain for each model.
    ! NOTE: no need for CLM's cell numbers, which have already set in clm_interface.
    select case(map_id)
      case(CLM_3DSUB_TO_PF_3DSUB)
        ! none
      case(PF_3DSUB_TO_CLM_3DSUB)
        clm_pf_idata%nlpf_sub  = grid_pf_npts_local
        clm_pf_idata%ngpf_sub  = grid_pf_npts_ghost+grid_pf_npts_local
      case default
        option%io_buffer = 'map_id argument NOT yet supported in ' // &
                        'pflotranModelInitMappingSubToSub'
        call printErrMsg(option)
    end select

  end subroutine pflotranModelInitMappingSub2Sub

! ************************************************************************** !

  subroutine pflotranModelGetTopFaceArea(pflotran_model)
  ! 
  ! This subroutine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 6/10/2013
  ! 

    use Option_module
    use Discretization_module
    use Patch_module
    use Grid_module
    use Grid_Unstructured_Aux_module
    use Grid_Unstructured_Cell_module
    use Grid_Unstructured_module
    use Utility_module, only : DotProduct, CrossProduct
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use pflotran_clm_main_module, only : pflotran_model_type
     use clm_pflotran_interface_data

    implicit none

    type(pflotran_model_type), pointer :: pflotran_model

    type(option_type), pointer         :: option
    class(realization_subsurface_type), pointer   :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(discretization_type), pointer :: discretization

    PetscInt :: local_id
    PetscInt :: ghosted_id
    PetscInt :: iface
    PetscInt :: face_id
    PetscInt :: cell_type
    PetscInt :: vertex_ids(4)

    PetscReal :: area1

    PetscScalar, pointer :: area_p(:)
    PetscErrorCode :: ierr

    subname = "ModelGetTopFaceArea"
    !-----------------------------------------------------------------------------------

    option          => pflotran_model%option
    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         option%io_buffer = "subroutine: " // trim(subname) //&
          " only works on subsurface simulations."
         call printErrMsg(option)
    end select
    discretization  => realization%discretization
    patch           => realization%patch
    grid            => patch%grid
    !------------

    call VecGetArrayF90(clm_pf_idata%area_top_face_pfp, area_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    if(grid%itype == STRUCTURED_GRID) then
      ! Structured grid
      do ghosted_id=1,grid%ngmax
        local_id = grid%nG2L(ghosted_id)
        if(local_id>0) then
          area1 = grid%structured_grid%dx(ghosted_id)* &
                  grid%structured_grid%dy(ghosted_id)
          area_p(local_id) = area1
        endif
      enddo
    else if (grid%itype == UNSTRUCTURED_GRID) then
      ! Unstructured grid
      do local_id = 1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        cell_type = grid%unstructured_grid%cell_type(local_id)

        ! Find iface
        if (cell_type == HEX_TYPE) then
          iface = 6
        else if (cell_type == WEDGE_TYPE) then
          iface = 5
        else
          call printErrMsg(pflotran_model%option, &
            'Only hex and wedge cell_type supported in CLM-PFLOTRAN')
        endif

        ! Get face-id
        face_id = grid%unstructured_grid%cell_to_face_ghosted(iface, ghosted_id)

        ! Save face area
        area_p(local_id) = grid%unstructured_grid%face_area(face_id)
      enddo
    endif
    call VecRestoreArrayF90(clm_pf_idata%area_top_face_pfp, area_p, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call MappingSourceToDestination(pflotran_model%map_pf_sub_to_clm_sub, &
                                    pflotran_model%option, &
                                    clm_pf_idata%area_top_face_pfp, &
                                    clm_pf_idata%area_top_face_clms)

  end subroutine pflotranModelGetTopFaceArea

! ************************************************************************** !

  subroutine pflotranModelInitMapFaceToFace(pflotran_model,  &
                                            grid_clm_cell_ids_nindex, &
                                            grid_clm_npts_local, &
                                            map_id)
  !
  ! This routine maps CLM grids/columns structure onto BC faces
  ! (TOP, BOTTOM, EAST, WEST, NORTH, or, SOUTH, which type depends on BC condition-name in PF input cards)
  ! of PFLOTRAN 3D Domain grid,
  ! by extending GB's code - Fengming Yuan, ORNL
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/09/13
  !
  ! 02/14/2014 - TOP/BOTTOM faces, from CLM => PF, finished
  ! 05/12/2014 - TOP/BOTTOM faces, from PF => CLM, finished
  !
  ! NOTE: for TOP face, BC condition name: 'clm_gflux_bc') ('g' for ground);
  !       for BOTTOM face, BC condition name: 'clm_bflux_bc') ('b' for bottom);

    use Option_module
    use Grid_module
    use Patch_module
    use Region_module
    use String_module
    use Simulation_Subsurface_class, only : simulation_subsurface_type
    use Realization_Subsurface_class, only : realization_subsurface_type

    use pflotran_clm_main_module, only : pflotran_model_type
    use clm_pflotran_interface_data

    implicit none

#include "petsc/finclude/petscviewer.h"

    type(pflotran_model_type), intent(inout), pointer :: pflotran_model
    PetscInt, intent(in), pointer                     :: grid_clm_cell_ids_nindex(:)
    PetscInt, intent(in)                              :: grid_clm_npts_local
    PetscInt, intent(in)                              :: map_id

    ! local
    PetscInt                           :: local_id, grid_pf_npts_local, grid_pf_npts_ghost
    PetscInt                           :: grid_clm_npts_ghost, source_mesh_id
    PetscInt                           :: dest_mesh_id
    PetscInt, pointer                  :: grid_pf_cell_ids_nindex(:)
    PetscInt, pointer                  :: grid_pf_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_local_nindex(:)
    PetscInt, pointer                  :: grid_clm_cell_ids_nindex_copy(:)
    PetscInt                           :: count
    PetscInt                           :: ghosted_id, natural_id, icell
    PetscInt                           :: iconn
    PetscInt                           :: istart
    PetscInt, pointer                  :: int_array(:)
    PetscBool                          :: found
    PetscScalar,pointer                :: v_loc(:)
    PetscErrorCode                     :: ierr

    Vec                                :: face_ids
    Vec                                :: face_ids_loc
    IS                                 :: is_from
    IS                                 :: is_to
    VecScatter                         :: vec_scat

    type(mapping_type), pointer        :: map
    type(option_type), pointer         :: option
    class(realization_subsurface_type), pointer   :: realization
    type(grid_type), pointer           :: grid
    type(patch_type), pointer          :: patch
    type(region_type), pointer         :: cur_region
    character(len=MAXSTRINGLENGTH)     :: region_name

    subname = "ModelInitMapFaceToFace"
    !-----------------------------------------------------------------------------------

    option          => pflotran_model%option
    select type (simulation => pflotran_model%simulation)
      class is (simulation_subsurface_type)
         realization => simulation%realization
      class default
         option%io_buffer = "subroutine: " // trim(subname) //&
          " only works on subsurface simulations."
         call printErrMsg(option)
    end select
    patch           => realization%patch
    grid            => patch%grid
    !------------

    allocate(grid_clm_cell_ids_nindex_copy(grid_clm_npts_local))
    grid_clm_cell_ids_nindex_copy = grid_clm_cell_ids_nindex

    ! Choose the appriopriate map
    select case(map_id)
      case(CLM_2DTOP_TO_PF_2DTOP)
        map => pflotran_model%map_clm_2dtop_to_pf_2dtop
        source_mesh_id = CLM_FACE_MESH
        dest_mesh_id = PF_FACE_MESH
        region_name = 'top'

      case(PF_2DTOP_TO_CLM_2DTOP)
        map => pflotran_model%map_pf_2dtop_to_clm_2dtop
        source_mesh_id = PF_FACE_MESH
        dest_mesh_id = CLM_FACE_MESH
        region_name = 'top'

      case(CLM_2DBOT_TO_PF_2DBOT)
        map => pflotran_model%map_clm_2dbot_to_pf_2dbot
        source_mesh_id = CLM_FACE_MESH
        dest_mesh_id = PF_FACE_MESH
        region_name = 'bottom'

      case(PF_2DBOT_TO_CLM_2DBOT)
        map => pflotran_model%map_pf_2dbot_to_clm_2dbot
        source_mesh_id = PF_FACE_MESH
        dest_mesh_id = CLM_FACE_MESH
        region_name = 'bottom'

      case default
        option%io_buffer = 'map_id argument NOT yet support ' // &
          'pflotranModelInitMappingFaceToFace'
        call printErrMsg(option)
    end select

    if (len(trim(map%filename)) > 0) then
      ! Read mapping file
      if (index(map%filename, '.h5') > 0) then
        call MappingReadHDF5(map, map%filename, option)
      else
        call MappingReadTxtFile(map, map%filename, option)
      endif

    elseif(.not.option%mapping_files) then
      ! directly mapping between CLM and PF meshes, if no user-defined mapping file
      map%pflotran_nlev_mapped = grid%structured_grid%nz
      map%clm_nlev_mapped = clm_pf_idata%nzclm_mapped
      if (dest_mesh_id == PF_FACE_MESH) then
        call MappingFromCLMGrids(map, grid, PETSC_TRUE, option)
      elseif(dest_mesh_id == CLM_FACE_MESH) then
        call MappingFromCLMGrids(map, grid, PETSC_FALSE, option)
      endif

    else

        option%io_buffer = 'MUST provide mapping files between CLM and PFLOTRAN! '
        call printErrMsg(option)

    endif

    grid_clm_npts_ghost=0
    ! Allocate memory to identify if CLM cells are local or ghosted.
    ! Note: Presently all CLM cells are local
    allocate(grid_clm_local_nindex(grid_clm_npts_local))
    do local_id = 1, grid_clm_npts_local
      grid_clm_local_nindex(local_id) = local_id ! LOCAL
    enddo


    found=PETSC_FALSE
    grid_pf_npts_local = 0
    grid_pf_npts_ghost = 0
    ! find the specified face of PFLOTRAN domain, by region name
    if (source_mesh_id == PF_FACE_MESH .or. &
      dest_mesh_id == PF_FACE_MESH) then

      cur_region => patch%region_list%first
      do
        if (.not.associated(cur_region)) exit
        ! the top/bottom cells of CLM soil domain (3-D) with a REGION name 'top'/'bottom' in PF input card
        if (StringCompareIgnoreCase(trim(cur_region%name), trim(region_name))) then

          found=PETSC_TRUE
          ! Allocate memory to save cell ids and flag for local cells
          allocate(grid_pf_cell_ids_nindex(cur_region%num_cells))
          allocate(grid_pf_local_nindex(cur_region%num_cells))
          grid_pf_npts_local = cur_region%num_cells

          do icell = 1, grid_pf_npts_local
            local_id   = cur_region%cell_ids(icell)  ! 'cur_region%cell_ids' are local_ids?
            ghosted_id = grid%nL2G(local_id)
            natural_id = grid%nG2A(ghosted_id)

            grid_pf_cell_ids_nindex(icell) = natural_id - 1
            grid_pf_local_nindex(icell)    = local_id        ! LOCAL ID

          end do

          if(found) exit

        endif

        cur_region => cur_region%next
      end do

      if(.not.found) then
        pflotran_model%option%io_buffer = 'region name ' // trim(region_name) // &
          ' not found in PF REGIONs'
        call printErrMsg(pflotran_model%option)
      endif

    else
        option%io_buffer='Unknown PFLOTRAN face mesh id'
        call printErrMsg(option)

    endif

    !
    ! Step-1: Find face cells-ids of PFLOTRAN subsurface domain
    !
    call VecCreateMPI(option%mycomm, grid%nlmax, PETSC_DETERMINE, face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecSet(face_ids, -1.d0, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! Set 1.0 to all cells that make up a face of PFLOTRAN subsurface domain
    allocate(v_loc(grid_pf_npts_local))
    v_loc = 1.d0
    call VecSetValues(face_ids, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                      v_loc, INSERT_VALUES, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    deallocate(v_loc)

    call VecAssemblyBegin(face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecAssemblyEnd(face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(face_ids, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    count = 0
    do local_id=1,grid%nlmax
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif

    enddo
    call VecRestoreArrayF90(face_ids, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    allocate(int_array(grid_pf_npts_local))
    do iconn = 1, grid_pf_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_pf_npts_local, grid_pf_cell_ids_nindex, &
                         PETSC_COPY_VALUES, is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! create scatter context
    call VecCreateSeq(PETSC_COMM_SELF, grid_pf_npts_local, face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecScatterDestroy(vec_scat, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    count = 0
    do iconn = 1, grid_pf_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_pf_cell_ids_nindex(count) = INT(v_loc(iconn))
      endif

    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDestroy(face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    ! Step-2: Recompute 'map%s2d_i/jscr' for pf mesh
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                     PETSC_COPY_VALUES, is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    do iconn = 1, map%s2d_nwts
      if (dest_mesh_id == PF_FACE_MESH) then
         int_array(iconn) = map%s2d_icsr(iconn)
      elseif (source_mesh_id == PF_FACE_MESH) then
         int_array(iconn) = map%s2d_jcsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                 PETSC_COPY_VALUES, is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecScatterDestroy(vec_scat, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (dest_mesh_id == PF_FACE_MESH) then
          map%s2d_icsr(count) = INT(v_loc(iconn))
        elseif (source_mesh_id == PF_FACE_MESH) then
          map%s2d_jcsr(count) = INT(v_loc(iconn))
        endif
      endif

    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDestroy(face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of face cells in mapping dataset (  ' // map%filename // &
        ') does not match face cells on which BC is applied - PFLOTRAN. '
      call printErrMsg(option)
    endif
    call VecDestroy(face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    ! Step-3: Find face cells-ids of CLM soil/below-ground domain
    !
    allocate(v_loc(grid_clm_npts_local))
    v_loc = 1.d0
    call VecCreateSeq(PETSC_COMM_SELF, grid_clm_npts_local, face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecCreateMPI(option%mycomm, clm_pf_idata%nlclm_sub, PETSC_DECIDE, face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecSet(face_ids, -1.d0, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    ! Set 1.0 to all cells that make up surface of CLM subsurface domain
    call VecSetValues(face_ids, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                      v_loc, INSERT_VALUES, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    deallocate(v_loc)
    call VecAssemblyBegin(face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecAssemblyEnd(face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(face_ids, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) count = count + 1
    enddo

    istart = 0
    call MPI_Exscan(count, istart, ONE_INTEGER_MPI, MPIU_INTEGER, MPI_SUM, &
                    option%mycomm, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    count = 0
    do local_id=1,clm_pf_idata%nlclm_sub
      if(v_loc(local_id) == 1.d0) then
        v_loc(local_id) = istart + count
        count = count + 1
      endif
    enddo
    call VecRestoreArrayF90(face_ids, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    allocate(int_array(grid_clm_npts_local))
    do iconn = 1, grid_clm_npts_local
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    deallocate(int_array)

    call ISCreateGeneral(option%mycomm, grid_clm_npts_local, grid_clm_cell_ids_nindex_copy, &
                         PETSC_COPY_VALUES, is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)


    ! create scatter context
    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, &
                          ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecScatterDestroy(vec_scat, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    count = 0
    do iconn = 1, grid_clm_npts_local
      if (v_loc(iconn)>-1) then
        count = count + 1
        grid_clm_cell_ids_nindex_copy(count) = INT(v_loc(iconn))
      endif
    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDestroy(face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    ! Step-4: Recompute 'map%s2d_i/jscr' for clm mesh
    !
    call VecCreateSeq(PETSC_COMM_SELF, map%s2d_nwts, face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    allocate(int_array(map%s2d_nwts))
    do iconn = 1, map%s2d_nwts
      int_array(iconn) = iconn - 1
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)


    do iconn = 1, map%s2d_nwts
      if (source_mesh_id == CLM_FACE_MESH) then
         int_array(iconn) = map%s2d_jcsr(iconn)
      elseif (dest_mesh_id == CLM_FACE_MESH) then
         int_array(iconn) = map%s2d_icsr(iconn)
      endif
    enddo
    call ISCreateGeneral(option%mycomm, map%s2d_nwts, int_array, &
                         PETSC_COPY_VALUES, is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    deallocate(int_array)

    ! create scatter context
    call VecScatterCreate(face_ids, is_from, face_ids_loc, is_to, vec_scat, &
                          ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_from, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call ISDestroy(is_to, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecScatterBegin(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecScatterEnd(vec_scat, face_ids, face_ids_loc, INSERT_VALUES, &
                        SCATTER_FORWARD, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecScatterDestroy(vec_scat, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    call VecGetArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    count = 0
    do iconn = 1, map%s2d_nwts
      if (v_loc(iconn)>-1) then
        count = count + 1
        if (source_mesh_id == CLM_FACE_MESH) then
           map%s2d_jcsr(count) = INT(v_loc(iconn))
        elseif (dest_mesh_id == CLM_FACE_MESH) then
           map%s2d_icsr(count) = INT(v_loc(iconn))
        endif
      endif
    enddo
    call VecRestoreArrayF90(face_ids_loc, v_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)
    call VecDestroy(face_ids_loc, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    if(count /= map%s2d_nwts) then
      option%io_buffer='No. of face cells in mapping dataset does not ' // &
        'match face cells on which BC is applied - CLM.'
      call printErrMsg(option)
    endif

    call VecDestroy(face_ids, ierr)
    call where_checkerr(ierr, subname, __FILE__, __LINE__)

    !
    select case(source_mesh_id)
      case(CLM_FACE_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_clm_npts_local, &
                                                      grid_clm_npts_ghost, &
                                                      grid_clm_cell_ids_nindex_copy, &
                                                      grid_clm_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_pf_npts_local, &
                                                           grid_pf_npts_ghost, &
                                                           grid_pf_cell_ids_nindex, &
                                                           grid_pf_local_nindex)
      case(PF_FACE_MESH)
        call MappingSetSourceMeshCellIds(map, option, grid_pf_npts_local, &
                                                      grid_pf_npts_ghost, &
                                                      grid_pf_cell_ids_nindex, &
                                                      grid_pf_local_nindex)
        call MappingSetDestinationMeshCellIds(map, option, grid_clm_npts_local, &
                                                           grid_clm_npts_ghost, &
                                                           grid_clm_cell_ids_nindex_copy, &
                                                           grid_clm_local_nindex)
      case default
        option%io_buffer = 'Invalid argument source_mesh_id passed to ' // &
          'pflotranModelInitMappingFaceToFace'
        call printErrMsg(option)
    end select

    deallocate(grid_pf_cell_ids_nindex)
    deallocate(grid_pf_local_nindex)
    deallocate(grid_clm_cell_ids_nindex_copy)
    deallocate(grid_clm_local_nindex)

    call MappingDecompose(map, option)
    call MappingFindDistinctSourceMeshCellIds(map, option)
    call MappingCreateWeightMatrix(map, option)
    call MappingCreateScatterOfSourceMesh(map, option)
    call MappingFreeNotNeeded(map)

    ! Setting the number of cells constituting the face of the 3D
    ! subsurface domain for each model.
    ! NOTE: no need for CLM's cell numbers, which have already set in clm_interface.
    select case(map_id)
      case(CLM_2DTOP_TO_PF_2DTOP, CLM_2DBOT_TO_PF_2DBOT)
        ! none
      case(PF_2DTOP_TO_CLM_2DTOP)
        clm_pf_idata%nlpf_2dtop  = grid_pf_npts_local
        clm_pf_idata%ngpf_2dtop  = grid_pf_npts_ghost+grid_pf_npts_local
      case(PF_2DBOT_TO_CLM_2DBOT)
        clm_pf_idata%nlpf_2dbot  = grid_pf_npts_local
        clm_pf_idata%ngpf_2dbot  = grid_pf_npts_ghost+grid_pf_npts_local
      case default
        option%io_buffer = 'map_id argument NOT yet supported in ' // &
                        'pflotranModelInitMappingFaceToFace'
        call printErrMsg(option)
    end select

  end subroutine pflotranModelInitMapFaceToFace

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  
end module pflotran_clm_setmapping_module

#endif

