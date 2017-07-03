module Output_Geomechanics_module

  use Output_Aux_module  
  use Output_Tecplot_module
  use Output_Common_module
  use Output_HDF5_module
  use PFLOTRAN_Constants_module
  
  implicit none

  private
  
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscdm.h"
#include "petsc/finclude/petscdm.h90"
#include "petsc/finclude/petsclog.h"

  PetscInt, save, public :: max_local_node_size_saved = -1
  PetscBool :: geomech_hdf5_first

  public :: OutputGeomechanics, &
            OutputGeomechInit, &
            OutputGeomechGetVarFromArray
      
contains

! ************************************************************************** !

subroutine OutputGeomechInit(num_steps)
  ! 
  ! Initializes module variables for geomechanics variables
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/2/13
  ! 

  use Option_module

  implicit none
  
  PetscInt :: num_steps

  if (num_steps == 0) then
    geomech_hdf5_first = PETSC_TRUE
  else
    geomech_hdf5_first = PETSC_FALSE
  endif
  
end subroutine OutputGeomechInit

! ************************************************************************** !

subroutine OutputGeomechanics(geomech_realization,snapshot_plot_flag, &
                              observation_plot_flag,massbal_plot_flag)
  ! 
  ! Main driver for all geomechanics output
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/2/13
  ! 

  use Geomechanics_Realization_class
  use Logging_module
  use Option_module, only : OptionCheckTouch, option_type, &
                            printMsg, printErrMsg

  implicit none

  type(realization_geomech_type) :: geomech_realization
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend
  type(option_type), pointer :: option

  option => geomech_realization%option

  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)

  ! check for plot request from active directory
  if (.not.snapshot_plot_flag) then

    if (option%use_touch_options) then
      string = 'plot'
      if (OptionCheckTouch(option,string)) then
        geomech_realization%output_option%plot_name = 'plot'
        snapshot_plot_flag = PETSC_TRUE
      endif
    endif
  endif

!.....................................
  if (snapshot_plot_flag) then
  
    if (geomech_realization%output_option%print_hdf5) then
       call OutputHDF5UGridXDMFGeomech(geomech_realization, &
                                       INSTANTANEOUS_VARS)
    endif
   
    if (geomech_realization%output_option%print_tecplot) then
      call PetscTime(tstart,ierr);CHKERRQ(ierr)
      call OutputTecplotGeomechanics(geomech_realization)
      call PetscTime(tend,ierr);CHKERRQ(ierr)
    endif

  endif

!......................................
  if (observation_plot_flag) then
  endif

!......................................
  if (massbal_plot_flag) then
  endif

  ! Increment the plot number
  if (snapshot_plot_flag) then
    geomech_realization%output_option%plot_number = &
      geomech_realization%output_option%plot_number + 1
  endif

  call PetscLogStagePop(ierr);CHKERRQ(ierr)

end subroutine OutputGeomechanics

! ************************************************************************** !

subroutine OutputTecplotGeomechanics(geomech_realization)
  ! 
  ! Tecplot output for geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/2/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module
  use Option_module
  
  implicit none

  type(realization_geomech_type) :: geomech_realization
  
  PetscInt, parameter :: icolumn = -1
  character(len=MAXSTRINGLENGTH) :: filename, string, string2
  character(len=MAXSTRINGLENGTH) :: tmp_global_prefix
  character(len=MAXWORDLENGTH) :: word
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: geomech_field
  type(geomech_patch_type), pointer :: patch 
  type(output_option_type), pointer :: output_option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  Vec :: global_cconn_vec
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: ivar, isubvar, var_type
  PetscErrorCode :: ierr  
  
  type(gmdm_type), pointer :: gmdm_element
  
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_field => geomech_realization%geomech_field
  output_option => geomech_realization%output_option

  tmp_global_prefix = option%global_prefix 
  option%global_prefix = trim(tmp_global_prefix) // '-geomech'
  filename = OutputFilename(output_option,option,'tec','')
  option%global_prefix = tmp_global_prefix
    
  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write tecplot geomech output file: ' // &
        trim(filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=filename,action="write")
    call OutputTecplotHeader(OUTPUT_UNIT,geomech_realization,icolumn)
  endif

  ! write blocks
  ! write out data sets
  call GeomechDiscretizationCreateVector(geomech_discretization,ONEDOF, &
                                         global_vec,GLOBAL,option)
  call GeomechDiscretizationCreateVector(geomech_discretization,ONEDOF, &
                                         natural_vec,NATURAL,option)

  ! write out coordinates
  call WriteTecplotGeomechGridVertices(OUTPUT_UNIT,geomech_realization)

  ! loop over variables and write to file
  cur_variable => output_option%output_snap_variable_list%first
  do
    if (.not.associated(cur_variable)) exit
    call OutputGeomechGetVarFromArray(geomech_realization,global_vec, &
                                      cur_variable%ivar, &
                                      cur_variable%isubvar)
    call GeomechDiscretizationGlobalToNatural(geomech_discretization, &
                                              global_vec, &
                                              natural_vec,ONEDOF)
    if (cur_variable%iformat == 0) then
      call WriteTecplotDataSetGeomechFromVec(OUTPUT_UNIT,geomech_realization, &
                                             natural_vec, &
                                             TECPLOT_REAL)
    else
      call WriteTecplotDataSetGeomechFromVec(OUTPUT_UNIT,geomech_realization, &
                                             natural_vec, &
                                             TECPLOT_INTEGER)
    endif
    cur_variable => cur_variable%next
  enddo

  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)

  ! write vertices
  call WriteTecplotGeomechGridElements(OUTPUT_UNIT,geomech_realization)

  if (option%myrank == option%io_rank) close(OUTPUT_UNIT)

end subroutine OutputTecplotGeomechanics

! ************************************************************************** !

subroutine WriteTecplotGeomechGridElements(fid,geomech_realization)
  ! 
  ! This subroutine writes unstructured grid elements
  ! for geomechanics grid
  ! 
  ! Author: Satish Karra
  ! Date: 07/03/2013
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Patch_module
  
  implicit none

  PetscInt :: fid
  type(realization_geomech_type) :: geomech_realization

  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch 
  Vec :: global_cconn_vec
  type(gmdm_type), pointer :: gmdm_element
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  Vec :: global_vec
  Vec :: natural_vec

  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  
  call GMCreateGMDM(grid,gmdm_element,EIGHT_INTEGER,option)
  call GMGridDMCreateVectorElem(grid,gmdm_element,global_vec, &
                            GLOBAL,option) 
  call GMGridDMCreateVectorElem(grid,gmdm_element,natural_vec, &
                            NATURAL,option) 
  call OutputGetCellVerticesGeomech(grid,global_vec)
  call VecScatterBegin(gmdm_element%scatter_gton_elem,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(gmdm_element%scatter_gton_elem,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call WriteTecplotDataSetNumPerLineGeomech(fid,geomech_realization,vec_ptr, &
                                            TECPLOT_INTEGER, &
                                            grid%nlmax_elem*8, &
                                            EIGHT_INTEGER)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call GMDMDestroy(gmdm_element)

end subroutine WriteTecplotGeomechGridElements

! ************************************************************************** !

subroutine OutputGetCellVerticesGeomech(grid,vec)
  ! 
  ! This routine returns a vector containing vertex ids
  ! in natural order of local cells for geomech grid
  ! 
  ! Author: Satish Karra
  ! Date: 07/03/2013
  ! 

  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  
  implicit none
  
  type(geomech_grid_type) :: grid
  Vec :: vec
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr  
    
  call GeomechGridVecGetArrayF90(grid,vec,vec_ptr,ierr)

  ! initialize
  vec_ptr = UNINITIALIZED_DOUBLE
  do local_id = 1, grid%nlmax_elem
    ghosted_id = local_id
    select case(grid%elem_type(ghosted_id))
      case(HEX_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 8
          vec_ptr(offset + ivertex) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(ivertex,local_id))
        enddo
      case(WEDGE_TYPE)
        offset = (local_id-1)*8
        vec_ptr(offset + 1) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(1,local_id))
        vec_ptr(offset + 2) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(1,local_id))
        vec_ptr(offset + 3) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(4,local_id))
        vec_ptr(offset + 4) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(4,local_id))
        vec_ptr(offset + 5) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(3,local_id))
        vec_ptr(offset + 6) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(2,local_id))
        vec_ptr(offset + 7) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(5,local_id))
        vec_ptr(offset + 8) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(6,local_id))
      case (PYR_TYPE)
        offset = (local_id-1)*8
        ! from Tecplot 360 Data Format Guide
        ! n1=vert1,n2=vert2,n3=vert3,n4=vert4,n5=n6=n7=n8=vert5
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(ivertex,local_id))
        enddo
        do ivertex = 5, 8
          vec_ptr(offset + ivertex) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(5,local_id))
        enddo
      case (TET_TYPE)
        offset = (local_id-1)*8
        ! from Tecplot 360 Data Format Guide
        ! n1=vert1,n2=vert2,n3=n4=vert3,n5=vert5=n6=n7=n8=vert4
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(ivertex,local_id))
        enddo
        vec_ptr(offset + 4) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(3,local_id))
        do ivertex = 5, 8
          vec_ptr(offset + ivertex) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(4,local_id))
        enddo
      case (QUAD_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(ivertex,local_id))
        enddo
      case (TRI_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
            grid%node_ids_ghosted_natural(grid%elem_nodes(ivertex,local_id))
        enddo
        ivertex = 4
        vec_ptr(offset + ivertex) = &
          grid%node_ids_ghosted_natural(grid%elem_nodes(3,local_id))
    end select
  enddo

  call GeomechGridVecRestoreArrayF90(grid,vec,vec_ptr,ierr)

end subroutine OutputGetCellVerticesGeomech

! ************************************************************************** !

subroutine OutputTecplotHeader(fid,geomech_realization,icolumn)
  ! 
  ! Prints Tecplot header for geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/2/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Patch_module

  implicit none

  PetscInt :: fid
  type(realization_geomech_type) :: geomech_realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  PetscInt :: variable_count
  
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  output_option => geomech_realization%output_option

  ! write header
  ! write title
  write(fid,'(''TITLE = "'',1es13.5," [",a1,'']"'')') &
                option%geomech_time/output_option%tconv,output_option%tunit

  ! initial portion of header
  string = 'VARIABLES=' // &
           '"X [m]",' // &
           '"Y [m]",' // &
           '"Z [m]"'
  write(fid,'(a)',advance="no") trim(string)

  call OutputWriteVariableListToHeader(fid, &
                                      output_option%output_snap_variable_list, &
                                       '',icolumn,PETSC_TRUE,variable_count)
 ! need to terminate line
  write(fid,'(a)') ''
  ! add x, y, z variables to count
  variable_count = variable_count + 3

  !geh: due to pgi bug, cannot embed functions with calls to write() within
  !     write statement
  call OutputWriteTecplotZoneHeader(fid,geomech_realization,variable_count, &
                                    output_option%tecplot_format)

end subroutine OutputTecplotHeader

! ************************************************************************** !

subroutine OutputWriteTecplotZoneHeader(fid,geomech_realization, &
                                        variable_count,tecplot_format)
  ! 
  ! Prints zone header to a tecplot file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/2/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use String_module
  
  implicit none

  PetscInt :: fid
  type(realization_geomech_type) :: geomech_realization
  PetscInt :: variable_count
  PetscInt :: tecplot_format
  
  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  
  grid => geomech_realization%geomech_patch%geomech_grid
  option => geomech_realization%option
  output_option => geomech_realization%output_option
  
  
  string = 'ZONE T="' // &
           trim(StringFormatDouble(option%time/output_option%tconv)) // &
           '"'
  string2 = ''
  select case(tecplot_format)
    case (TECPLOT_POINT_FORMAT)
      string2 = 'POINT format not supported for geomechanics ' // &
                'unstructured'
      string2 = trim(string2) // &
              ', DATAPACKING=POINT'
    case default !(TECPLOT_BLOCK_FORMAT,TECPLOT_FEBRICK_FORMAT)
      string2 = ', N=' // &
                trim(StringFormatInt(grid%nmax_node)) // &
                ', ELEMENTS=' // &
                trim(StringFormatInt(grid%nmax_elem))
      string2 = trim(string2) // ', ZONETYPE=FEBRICK' 
      
      string3 = ', VARLOCATION=(NODAL)'

      string2 = trim(string2) // trim(string3) // ', DATAPACKING=BLOCK'
  end select
  
  write(fid,'(a)') trim(string) // trim(string2)

end subroutine OutputWriteTecplotZoneHeader

! ************************************************************************** !

subroutine WriteTecplotGeomechGridVertices(fid,geomech_realization)
  ! 
  ! Prints zone header to a tecplot file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/2/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Option_module
  use Geomechanics_Patch_module
  use Variables_module

  implicit none

  PetscInt :: fid
  type(realization_geomech_type) :: geomech_realization 
  
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch 
  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vertex_vec
  PetscInt :: local_size
  PetscErrorCode :: ierr
  
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%nmax_node, &
                    global_vertex_vec,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call OutputGetVertexCoordinatesGeomech(grid,global_vertex_vec,X_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call WriteTecplotDataSetGeomech(fid,geomech_realization,vec_ptr, &
                                  TECPLOT_REAL,local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr);CHKERRQ(ierr)

  call OutputGetVertexCoordinatesGeomech(grid,global_vertex_vec,Y_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call WriteTecplotDataSetGeomech(fid,geomech_realization,vec_ptr, &
                                  TECPLOT_REAL,local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr);CHKERRQ(ierr)

  call OutputGetVertexCoordinatesGeomech(grid,global_vertex_vec,Z_COORDINATE,option)
  call VecGetArrayF90(global_vertex_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call WriteTecplotDataSetGeomech(fid,geomech_realization,vec_ptr, &
                                  TECPLOT_REAL,local_size)
  call VecRestoreArrayF90(global_vertex_vec,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(global_vertex_vec,ierr);CHKERRQ(ierr)

end subroutine WriteTecplotGeomechGridVertices

! ************************************************************************** !

subroutine OutputGetVertexCoordinatesGeomech(grid,vec,direction,option)
  ! 
  ! Extracts vertex coordinates of cells into
  ! a PetscVec
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/02/2013
  ! 

  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Variables_module, only : X_COORDINATE, Y_COORDINATE, Z_COORDINATE
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  type(geomech_grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  type(option_type) :: option
  
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscInt, allocatable :: indices(:)
  PetscReal, allocatable :: values(:)
  PetscErrorCode :: ierr
  
  if (option%mycommsize == 1) then
    call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          vec_ptr(ivertex) = grid%nodes(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          vec_ptr(ivertex) = grid%nodes(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          vec_ptr(ivertex) = grid%nodes(ivertex)%z
        enddo
    end select
    call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  else
    ! initialize to UNINITIALIZED_DOUBLE to catch bugs
    call VecSet(vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    allocate(values(grid%nlmax_node))
    allocate(indices(grid%nlmax_node))
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          values(ivertex) = grid%nodes(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          values(ivertex) = grid%nodes(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%nlmax_node
          values(ivertex) = grid%nodes(ivertex)%z
        enddo
    end select
    indices(:) = grid%node_ids_local_natural(:)-1
    call VecSetValues(vec,grid%nlmax_node, &
                      indices,values,INSERT_VALUES,ierr);CHKERRQ(ierr)
    call VecAssemblyBegin(vec,ierr);CHKERRQ(ierr)
    deallocate(values)
    deallocate(indices)
    call VecAssemblyEnd(vec,ierr);CHKERRQ(ierr)
  endif
  
end subroutine OutputGetVertexCoordinatesGeomech

! ************************************************************************** !

subroutine OutputGeomechGetVarFromArray(geomech_realization,vec,ivar,isubvar, &
                                        isubvar1)
  ! 
  ! Gets variables from an array
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/3/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Field_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petsclog.h"

  class(realization_geomech_type) :: geomech_realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubvar1
  
  PetscErrorCode :: ierr

  call GeomechRealizGetDataset(geomech_realization,vec,ivar,isubvar,isubvar1)

end subroutine OutputGeomechGetVarFromArray

! ************************************************************************** !

subroutine WriteTecplotDataSetGeomechFromVec(fid,geomech_realization,vec, &
                                             datatype)
  ! 
  ! Writes data from a Petsc Vec within a block
  ! of a Tecplot file
  ! 
  ! Author: Satish Karra
  ! Date: 07/03//13
  ! 

  use Geomechanics_Realization_class

  implicit none

  PetscInt :: fid
  type(realization_geomech_type) :: geomech_realization
  Vec :: vec
  PetscInt :: datatype
  PetscErrorCode :: ierr  
  
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  call WriteTecplotDataSetGeomech(fid,geomech_realization,vec_ptr, &
                                  datatype,ZERO_INTEGER) 
  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine WriteTecplotDataSetGeomechFromVec

! ************************************************************************** !

subroutine WriteTecplotDataSetGeomech(fid,geomech_realization,array,datatype, &
                                      size_flag)
  ! 
  ! Writes data from an array within a block
  ! of a Tecplot file
  ! 
  ! Author: Satish Karra
  ! Date: 07/02//13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Patch_module

  implicit none

  PetscInt :: fid
  type(realization_geomech_type) :: geomech_realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size

  PetscInt, parameter :: num_per_line = 10

  call WriteTecplotDataSetNumPerLineGeomech(fid,geomech_realization,array, &
                                            datatype,size_flag,num_per_line) 
  
end subroutine WriteTecplotDataSetGeomech

! ************************************************************************** !

subroutine WriteTecplotDataSetNumPerLineGeomech(fid,geomech_realization, &
                                                array,datatype, &
                                                size_flag,num_per_line)
  ! 
  ! WriteTecplotDataSetNumPerLine: Writes data from an array within a block
  ! of a Tecplot file with a specified number
  ! of values per line
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07, 12/02/11, Satish Karra 07/02/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Geomechanics_Patch_module

  implicit none
  
  PetscInt :: fid
  type(realization_geomech_type) :: geomech_realization
  PetscReal :: array(:)
  PetscInt :: datatype
  PetscInt :: size_flag ! if size_flag /= 0, use size_flag as the local size
  PetscInt :: num_per_line
  
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch  
  PetscInt :: i
  PetscInt :: max_proc, max_proc_prefetch
  PetscMPIInt :: iproc_mpi, recv_size_mpi
  PetscInt :: max_local_size
  PetscMPIInt :: local_size_mpi
  PetscInt :: istart, iend, num_in_array
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscInt, allocatable :: integer_data(:), integer_data_recv(:)
  PetscReal, allocatable :: real_data(:), real_data_recv(:)
  PetscErrorCode :: ierr  
  
1000 format(100(i2,1x))
1001 format(100(i4,1x))
1002 format(100(i6,1x))
1003 format(100(i8,1x))
1004 format(100(i10,1x))
1010 format(100(es13.6,1x))

  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option

  ! if num_per_line exceeds 100, need to change the format statement below
  if (num_per_line > 100) then
    option%io_buffer = 'Number of values to be written to line in ' // &
      'WriteTecplotDataSetNumPerLine() exceeds 100.  ' // &
      'Must fix format statements.'
    call printErrMsg(option)
  endif

  ! maximum number of initial messages  
#define HANDSHAKE  
  max_proc = option%io_handshake_buffer_size
  max_proc_prefetch = option%io_handshake_buffer_size / 10

  if (size_flag /= 0) then
    call MPI_Allreduce(size_flag,max_local_size,ONE_INTEGER_MPI,MPIU_INTEGER, &
                       MPI_MAX,option%mycomm,ierr)
    local_size_mpi = size_flag
  else 
  ! if first time, determine the maximum size of any local array across 
  ! all procs
    if (max_local_node_size_saved < 0) then
      call MPI_Allreduce(grid%nlmax_node,max_local_size,ONE_INTEGER_MPI, &
                         MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
      max_local_node_size_saved = max_local_size
      write(option%io_buffer,'("max_local_node_size_saved: ",i9)') &
          max_local_size
      call printMsg(option)
    endif
    max_local_size = max_local_node_size_saved
    local_size_mpi = grid%nlmax_node
  endif
  
  ! transfer the data to an integer or real array
  if (datatype == TECPLOT_INTEGER) then
    allocate(integer_data(max_local_size+10))
    allocate(integer_data_recv(max_local_size))
    do i=1,local_size_mpi
      integer_data(i) = int(array(i))
    enddo
  else
    allocate(real_data(max_local_size+10))
    allocate(real_data_recv(max_local_size))
    do i=1,local_size_mpi
      real_data(i) = array(i)
    enddo
  endif
  
  ! communicate data to processor 0, round robin style
  if (option%myrank == option%io_rank) then
    if (datatype == TECPLOT_INTEGER) then
      ! This approach makes output files identical, regardless of processor
      ! distribution.  It is necessary when diffing files.
      iend = 0
      do
        istart = iend+1
        if (iend+num_per_line > local_size_mpi) exit
        iend = istart+(num_per_line-1)
        i = abs(maxval(integer_data(istart:iend)))
        if (i < 10) then
          write(fid,1000) integer_data(istart:iend)
        else if (i < 1000) then
          write(fid,1001) integer_data(istart:iend)
        else if (i < 100000) then
          write(fid,1002) integer_data(istart:iend)
        else if (i < 10000000) then
          write(fid,1003) integer_data(istart:iend)
        else
          write(fid,1004) integer_data(istart:iend)
        endif
      enddo
      ! shift remaining data to front of array
      integer_data(1:local_size_mpi-iend) = integer_data(iend+1:local_size_mpi)
      num_in_array = local_size_mpi-iend
    else
      iend = 0
      do
        istart = iend+1
        if (iend+num_per_line > local_size_mpi) exit
        iend = istart+(num_per_line-1)
        ! if num_per_line exceeds 100, need to change the format statement below
        write(fid,1010) real_data(istart:iend)
      enddo
      ! shift remaining data to front of array
      real_data(1:local_size_mpi-iend) = real_data(iend+1:local_size_mpi)
      num_in_array = local_size_mpi-iend
    endif
    do iproc_mpi=1,option%mycommsize-1
#ifdef HANDSHAKE    
      if (option%io_handshake_buffer_size > 0 .and. &
          iproc_mpi+max_proc_prefetch >= max_proc) then
        max_proc = max_proc + option%io_handshake_buffer_size
        call MPI_Bcast(max_proc,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                       option%mycomm,ierr)
      endif
#endif      
      call MPI_Probe(iproc_mpi,MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
      recv_size_mpi = status_mpi(MPI_TAG)
      if (datatype == TECPLOT_INTEGER) then
        call MPI_Recv(integer_data_recv,recv_size_mpi,MPIU_INTEGER,iproc_mpi, &
                      MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
        if (recv_size_mpi > 0) then
          integer_data(num_in_array+1:num_in_array+recv_size_mpi) = &
                                             integer_data_recv(1:recv_size_mpi)
          num_in_array = num_in_array+recv_size_mpi
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+num_per_line > num_in_array) exit
          iend = istart+(num_per_line-1)
          i = abs(maxval(integer_data(istart:iend)))
          if (i < 10) then
            write(fid,1000) integer_data(istart:iend)
          else if (i < 1000) then
            write(fid,1001) integer_data(istart:iend)
          else if (i < 100000) then
            write(fid,1002) integer_data(istart:iend)
          else if (i < 10000000) then
            write(fid,1003) integer_data(istart:iend)
          else
            write(fid,1004) integer_data(istart:iend)
          endif
        enddo
        if (iend > 0) then
          integer_data(1:num_in_array-iend) = integer_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      else
        call MPI_Recv(real_data_recv,recv_size_mpi, &
                      MPI_DOUBLE_PRECISION,iproc_mpi, &
                      MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
        if (recv_size_mpi > 0) then
          real_data(num_in_array+1:num_in_array+recv_size_mpi) = &
                                             real_data_recv(1:recv_size_mpi)
          num_in_array = num_in_array+recv_size_mpi
        endif
        iend = 0
        do
          istart = iend+1
          if (iend+num_per_line > num_in_array) exit
          iend = istart+(num_per_line-1)
          ! if num_per_line exceeds 100, need to change the format statement below
          write(fid,1010) real_data(istart:iend)
        enddo
        if (iend > 0) then
          real_data(1:num_in_array-iend) = real_data(iend+1:num_in_array)
          num_in_array = num_in_array-iend
        endif
      endif
    enddo
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      max_proc = -1
      call MPI_Bcast(max_proc,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                     option%mycomm,ierr)
    endif
#endif      
    ! Print the remaining values, if they exist
    if (datatype == TECPLOT_INTEGER) then
      if (num_in_array > 0) then
        i = abs(maxval(integer_data(1:num_in_array)))
        if (i < 10) then
          write(fid,1000) integer_data(1:num_in_array)
        else if (i < 1000) then
          write(fid,1001) integer_data(1:num_in_array)
        else if (i < 100000) then
          write(fid,1002) integer_data(1:num_in_array)
        else if (i < 10000000) then
          write(fid,1003) integer_data(1:num_in_array)
        else
          write(fid,1004) integer_data(1:num_in_array)
        endif
      endif
    else
      if (num_in_array > 0) &
        write(fid,1010) real_data(1:num_in_array)
    endif
  else
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        if (option%myrank < max_proc) exit
        call MPI_Bcast(max_proc,1,MPIU_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
      enddo
    endif
#endif    
    if (datatype == TECPLOT_INTEGER) then
      call MPI_Send(integer_data,local_size_mpi,MPIU_INTEGER,option%io_rank, &
                    local_size_mpi,option%mycomm,ierr)
    else
      call MPI_Send(real_data,local_size_mpi, &
                    MPI_DOUBLE_PRECISION,option%io_rank, &
                    local_size_mpi,option%mycomm,ierr)
    endif
#ifdef HANDSHAKE    
    if (option%io_handshake_buffer_size > 0) then
      do
        call MPI_Bcast(max_proc,1,MPIU_INTEGER,option%io_rank,option%mycomm, &
                       ierr)
        if (max_proc < 0) exit
      enddo
    endif
#endif
#undef HANDSHAKE
  endif
      
  if (datatype == TECPLOT_INTEGER) then
    deallocate(integer_data)
  else
    deallocate(real_data)
  endif

end subroutine WriteTecplotDataSetNumPerLineGeomech

! ************************************************************************** !

subroutine OutputXMFHeaderGeomech(fid,time,nmax,xmf_vert_len,ngvert,filename)
  ! 
  ! This subroutine writes header to a .xmf file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/3/13
  ! 

  implicit none

  PetscInt :: fid, vert_count
  PetscReal :: time
  PetscInt :: nmax,xmf_vert_len,ngvert
  character(len=MAXSTRINGLENGTH) :: filename

  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  
  string="<?xml version=""1.0"" ?>"
  write(fid,'(a)') trim(string)
  
  string="<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>"
  write(fid,'(a)') trim(string)

  string="<Xdmf>"
  write(fid,'(a)') trim(string)

  string="  <Domain>"
  write(fid,'(a)') trim(string)

  string="    <Grid Name=""Mesh"">"
  write(fid,'(a)') trim(string)

  write(string2,'(es13.5)') time
  string="      <Time Value = """ // trim(adjustl(string2)) // """ />"
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string="      <Topology Type=""Mixed"" NumberOfElements=""" // &
    trim(adjustl(string2)) // """ >"
  write(fid,'(a)') trim(string)

  write(string2,*) xmf_vert_len
  string="        <DataItem Format=""HDF"" DataType=""Int"" Dimensions=""" // &
    trim(adjustl(string2)) // """>"
  write(fid,'(a)') trim(string)

  string="          "//trim(filename) //":/Domain/Cells"
  write(fid,'(a)') trim(string)

  string="        </DataItem>"
  write(fid,'(a)') trim(string)

  string="      </Topology>"
  write(fid,'(a)') trim(string)

  string="      <Geometry GeometryType=""XYZ"">"
  write(fid,'(a)') trim(string)

  write(string2,*) ngvert
  string="        <DataItem Format=""HDF"" Dimensions=""" // &
      trim(adjustl(string2)) // " 3"">"
  write(fid,'(a)') trim(string)

  string="          "//trim(filename) //":/Domain/Vertices"
  write(fid,'(a)') trim(string)

  string="        </DataItem>"
  write(fid,'(a)') trim(string)

  string="      </Geometry>"
  write(fid,'(a)') trim(string)
 
#if 0
  string="      <Attribute Name=""X"" AttributeType=""Scalar""  Center=""Node"">"
  write(fid,'(a)') trim(string)

  write(string2,*) ngvert
  string="        <DataItem Dimensions=""" // &
      trim(adjustl(string2)) // " 1"" Format=""HDF""> "
  write(fid,'(a)') trim(string)

  string="        " // trim(filename) //":/Domain/X"
  write(fid,'(a)') trim(string)

  string="        </DataItem> " 
  write(fid,'(a)') trim(string)

  string="      </Attribute>"
  write(fid,'(a)') trim(string)  
#endif

end subroutine OutputXMFHeaderGeomech

! ************************************************************************** !

subroutine OutputXMFFooterGeomech(fid)
  ! 
  ! This subroutine writes footer to a .xmf file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/3/13
  ! 

  implicit none

  PetscInt :: fid

  character(len=MAXSTRINGLENGTH) :: string

  string="    </Grid>"
  write(fid,'(a)') trim(string)

  string="  </Domain>"
  write(fid,'(a)') trim(string)

  string="</Xdmf>"
  write(fid,'(a)') trim(string)

end subroutine OutputXMFFooterGeomech

! ************************************************************************** !

subroutine OutputXMFAttributeGeomech(fid,nmax,attname,att_datasetname)
  ! 
  ! This subroutine writes an attribute to a .xmf file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/3/13
  ! 

  implicit none

  PetscInt :: fid,nmax
  
  character(len=MAXSTRINGLENGTH) :: attname, att_datasetname
  character(len=MAXSTRINGLENGTH) :: string,string2
  string="      <Attribute Name=""" // trim(attname) // &
    """ AttributeType=""Scalar""  Center=""Node"">"
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string="        <DataItem Dimensions=""" // &
      trim(adjustl(string2)) // " 1"" Format=""HDF""> "
  write(fid,'(a)') trim(string)

  string="        " // trim(att_datasetname)
  write(fid,'(a)') trim(string)

  string="        </DataItem> " 
  write(fid,'(a)') trim(string)

  string="      </Attribute>"
  write(fid,'(a)') trim(string)

end subroutine OutputXMFAttributeGeomech

! ************************************************************************** !

subroutine OutputHDF5UGridXDMFGeomech(geomech_realization,var_list_type)
  ! 
  ! This routine writes unstructured grid data
  ! in HDF5 XDMF format
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/3/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Option_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module

#if  !defined(PETSC_HAVE_HDF5)

  implicit none
  
  type(realization_geomech_type) :: geomech_realization
  PetscInt :: var_list_type

  call printMsg(geomech_realization%option,'')
  write(geomech_realization%option%io_buffer, &
        '("PFLOTRAN must be compiled with HDF5 to &
        &write HDF5 formatted structured grids Darn.")')
  call printErrMsg(geomech_realization%option)

#else

! 64-bit stuff
#ifdef PETSC_USE_64BIT_INDICES
!#define HDF_NATIVE_INTEGER H5T_STD_I64LE
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#else
#define HDF_NATIVE_INTEGER H5T_NATIVE_INTEGER
#endif

  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use HDF5_Aux_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petsclog.h"

  type(realization_geomech_type) :: geomech_realization
  PetscInt :: var_list_type

#if defined(SCORPIO_WRITE)
  integer :: file_id
  integer :: data_type
  integer :: grp_id
  integer :: file_space_id
  integer :: memory_space_id
  integer :: data_set_id
  integer :: realization_set_id
  integer :: prop_id
  PetscMPIInt :: rank
  integer :: rank_mpi,file_space_rank_mpi
  integer :: dims(3)
  integer :: start(3), length(3), stride(3)
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: rank
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
#endif

  type(geomech_grid_type), pointer :: grid
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: field
  type(geomech_patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(option_type), pointer :: option
  type(output_variable_type), pointer :: cur_variable

  Vec :: global_vec
  Vec :: natural_vec
  PetscReal, pointer :: v_ptr

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: xmf_filename, att_datasetname, group_name
  character(len=MAXSTRINGLENGTH) :: string, string2,string3
  character(len=MAXWORDLENGTH) :: word
  character(len=2) :: free_mol_char, tot_mol_char, sec_mol_char
  PetscReal, pointer :: array(:)
  PetscInt :: istart
  PetscInt :: i
  PetscInt :: nviz_flow, nviz_tran, nviz_dof
  PetscInt :: current_component
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscFortranAddr :: app_ptr
  PetscMPIInt :: hdf5_err  
  PetscBool :: first
  PetscInt :: ivar, isubvar, var_type
  PetscInt :: vert_count
  PetscErrorCode :: ierr

  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  option => geomech_realization%option
  field => geomech_realization%geomech_field
  output_option => geomech_realization%output_option

  select case (var_list_type)
    case (INSTANTANEOUS_VARS)
      string2=''
      write(string3,'(i4)') output_option%plot_number
      xmf_filename = OutputFilename(output_option,option,'xmf','geomech')
    case (AVERAGED_VARS)
      string2='-aveg'
      write(string3,'(i4)') &
          int(option%time/output_option%periodic_snap_output_time_incr)
      xmf_filename = OutputFilename(output_option,option,'xmf','geomech_aveg')
  end select
  if (output_option%print_single_h5_file) then
    first = geomech_hdf5_first
    filename = trim(option%global_prefix) // trim(string2) // &
               trim(option%group_prefix) // '-geomech.h5'
  else
    string = OutputHDF5FilenameID(output_option,option,var_list_type)
    select case (var_list_type)
      case (INSTANTANEOUS_VARS)
        if (mod(output_option%plot_number, &
              output_option%times_per_h5_file)==0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
      case (AVERAGED_VARS)
        if (mod((option%time-output_option%periodic_snap_output_time_incr)/ &
                output_option%periodic_snap_output_time_incr, &
                dble(output_option%times_per_h5_file))==0) then
          first = PETSC_TRUE
        else
          first = PETSC_FALSE
        endif
    end select

    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               trim(string2) // '-' // trim(string) // '-geomech.h5'
  endif

  grid => patch%geomech_grid

#ifdef SCORPIO_WRITE
   option%io_buffer='OutputHDF5UGridXDMF not supported with SCORPIO_WRITE'
   call printErrMsg(option)
#endif

    ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif

  if (.not.first) then
    call h5eset_auto_f(OFF,hdf5_err)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
    if (hdf5_err /= 0) first = PETSC_TRUE
    call h5eset_auto_f(ON,hdf5_err)
  endif
  if (first) then
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                     H5P_DEFAULT_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (first) then
    option%io_buffer = '--> creating hdf5 geomech output file: ' // &
                       trim(filename)
  else
    option%io_buffer = '--> appending to hdf5 geomech output file: ' // &
                       trim(filename)
  endif
  call printMsg(option)

  if (first) then
    ! create a group for the coordinates data set
    string = "Domain"
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
    call WriteHDF5CoordinatesXDMFGeomech(geomech_realization,option,grp_id)
    call h5gclose_f(grp_id,hdf5_err)
  endif 

  if (option%myrank == option%io_rank) then
    option%io_buffer = '--> write xmf geomech output file: ' // &
                       trim(xmf_filename)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=xmf_filename,action="write")
    call OutputXMFHeaderGeomech(OUTPUT_UNIT, &
                         option%time/output_option%tconv, &
                         grid%nmax_elem, &
                         geomech_realization%output_option%xmf_vert_len, &
                         grid%nmax_node,filename)
  endif

  ! create a group for the data set
  write(string,'(''Time'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  if (len_trim(output_option%plot_name) > 2) then
    string = trim(string) // ' ' // output_option%plot_name
  endif
  string = trim(string3) // ' ' // trim(string)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  group_name=string
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  ! write out data sets 
  call GeomechDiscretizationCreateVector(geomech_discretization,ONEDOF, &
                                         global_vec, &
                                         GLOBAL,option)
  call GeomechDiscretizationCreateVector(geomech_discretization,ONEDOF, &
                                         natural_vec, &
                                         NATURAL,option)

  select case (var_list_type)

    case (INSTANTANEOUS_VARS)
      ! loop over variables and write to file
      cur_variable => output_option%output_snap_variable_list%first
      do
        if (.not.associated(cur_variable)) exit
        call OutputGeomechGetVarFromArray(geomech_realization,global_vec, &
                                          cur_variable%ivar, &
                                          cur_variable%isubvar)
        call GeomechDiscretizationGlobalToNatural(geomech_discretization, &
                                                  global_vec, &
                                                  natural_vec,ONEDOF)
        string = cur_variable%name
        if (len_trim(cur_variable%units) > 0) then
          word = cur_variable%units
          call HDF5MakeStringCompatible(word)
          string = trim(string) // ' [' // trim(word) // ']'
        endif
        if (cur_variable%iformat == 0) then
          call HDF5WriteDataSetFromVec(string,option, &
                                          natural_vec,grp_id,H5T_NATIVE_DOUBLE)
        else
          call HDF5WriteDataSetFromVec(string,option, &
                                          natural_vec,grp_id,H5T_NATIVE_INTEGER)
        endif
        att_datasetname = trim(filename) // ":/" // &
            trim(group_name) // "/" // trim(string)
        if (option%myrank == option%io_rank) then
          call OutputXMFAttributeGeomech(OUTPUT_UNIT,grid%nmax_node,string, &
                                         att_datasetname)
        endif
        cur_variable => cur_variable%next
      enddo

#if 0
    case (AVERAGED_VARS)
      if (associated(output_option%aveg_output_variable_list%first)) then
        cur_variable => output_option%aveg_output_variable_list%first
        do ivar = 1,output_option%aveg_output_variable_list%nvars
          string = 'Aveg. ' // cur_variable%name
          if (len_trim(cur_variable%units) > 0) then
            word = cur_variable%units
            call HDF5MakeStringCompatible(word)
            string = trim(string) // ' [' // trim(word) // ']'
          endif

          call GeomechDiscretizationGlobalToNatural(geomech_discretization, &
                                            field%avg_vars_vec(ivar), &
                                            natural_vec,ONEDOF)
          call HDF5WriteDataSetFromVec(string,option, &
                                          natural_vec,grp_id,H5T_NATIVE_DOUBLE)
          att_datasetname = trim(filename) // ":/" // &
              trim(group_name) // "/" // trim(string)
          if (option%myrank == option%io_rank) then
            call OutputXMFAttributeGeomech(OUTPUT_UNIT,grid%nlmax_node,string, &
                                           att_datasetname)
          endif
          cur_variable => cur_variable%next
        enddo
      endif
#endif
! 0

  end select

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call h5gclose_f(grp_id,hdf5_err)

  call h5fclose_f(file_id,hdf5_err)
  call h5close_f(hdf5_err)

  if (option%myrank == option%io_rank) then
    call OutputXMFFooterGeomech(OUTPUT_UNIT)
    close(OUTPUT_UNIT)
  endif

  geomech_hdf5_first = PETSC_FALSE
  
#endif
! !defined(PETSC_HAVE_HDF5)

end subroutine OutputHDF5UGridXDMFGeomech

#if defined(PETSC_HAVE_HDF5)

! ************************************************************************** !

subroutine WriteHDF5CoordinatesXDMFGeomech(geomech_realization, &
                                                option,file_id)
  ! 
  ! Writes the geomech coordinates in HDF5 file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/3/13
  ! 

  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_module
  use Variables_module
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petsclog.h"

  type(realization_geomech_type) :: geomech_realization
  type(option_type), pointer :: option

#if defined(SCORPIO_WRITE)
  integer :: file_id
  integer :: data_type
  integer :: grp_id
  integer :: file_space_id
  integer :: memory_space_id
  integer :: data_set_id
  integer :: realization_set_id
  integer :: prop_id
  integer :: dims(3)
  integer :: start(3), length(3), stride(3)
  integer :: rank_mpi,file_space_rank_mpi
  integer :: hdf5_flag
  integer, parameter :: ON=1, OFF=0
#else
  integer(HID_T) :: file_id
  integer(HID_T) :: data_type
  integer(HID_T) :: grp_id
  integer(HID_T) :: file_space_id
  integer(HID_T) :: realization_set_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: data_set_id
  integer(HID_T) :: prop_id
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  PetscMPIInt :: rank_mpi,file_space_rank_mpi
  PetscMPIInt :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0
#endif

  PetscInt :: istart
  type(geomech_grid_type), pointer :: grid
  character(len=MAXSTRINGLENGTH) :: string
  PetscMPIInt :: hdf5_err  

  PetscInt :: local_size,vert_count,nverts
  PetscInt :: i,j
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, pointer :: double_array(:)
  Vec :: global_x_vertex_vec,global_y_vertex_vec,global_z_vertex_vec
  Vec :: natural_x_vertex_vec,natural_y_vertex_vec,natural_z_vertex_vec

  PetscReal, pointer :: vec_ptr(:)
  Vec :: global_vec, natural_vec
  PetscInt, pointer :: int_array(:)
  type(gmdm_type),pointer :: gmdm_element
  PetscErrorCode :: ierr

  PetscInt :: TET_ID_XDMF = 6
  PetscInt :: PYR_ID_XDMF = 7
  PetscInt :: WED_ID_XDMF = 8
  PetscInt :: HEX_ID_XDMF = 9

  grid => geomech_realization%geomech_patch%geomech_grid

  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%nmax_node, &
                    global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%nmax_node, &
                    global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, &
                    grid%nmax_node, &
                    global_z_vertex_vec,ierr);CHKERRQ(ierr)

  call VecGetLocalSize(global_x_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_y_vertex_vec,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_z_vertex_vec,local_size,ierr);CHKERRQ(ierr)

  call OutputGetVertexCoordinatesGeomech(grid,global_x_vertex_vec, &
                                   X_COORDINATE,option)
  call OutputGetVertexCoordinatesGeomech(grid,global_y_vertex_vec, &
                                   Y_COORDINATE,option)
  call OutputGetVertexCoordinatesGeomech(grid,global_z_vertex_vec, &
                                   Z_COORDINATE,option)

  call VecGetArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)

#if defined(SCORPIO_WRITE)
  write(*,*) 'SCORPIO_WRITE'
  option%io_buffer = 'WriteHDF5CoordinatesUGrid not supported for SCORPIO_WRITE'
  call printErrMsg(option)
#else

  !
  !        not(SCORPIO_WRITE)
  !

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims = 0
  dims(1) = local_size * 3
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
   
  ! file space which is a 2D block
  rank_mpi = 2
  dims = 0
  dims(2) = grid%nmax_node
  dims(1) = 3
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Vertices" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(local_size, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(2) = istart
  start(1) = 0
  
  length(2) = local_size
  length(1) = 3

  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(double_array(local_size*3))
  do i=1,local_size
    double_array((i-1)*3+1) = vec_x_ptr(i)
    double_array((i-1)*3+2) = vec_y_ptr(i)
    double_array((i-1)*3+3) = vec_z_ptr(i)
  enddo

  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,double_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)

  deallocate(double_array)
  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(global_x_vertex_vec,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_y_vertex_vec,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_z_vertex_vec,vec_z_ptr,ierr);CHKERRQ(ierr)

  ! Vertex X/Y/Z
  ! X coord
  string = "X" // CHAR(0)
  
  call HDF5WriteDataSetFromVec(string,option, &
                                           global_x_vertex_vec, &
                                           file_id,H5T_NATIVE_DOUBLE)

  ! Y coord
  string = "Y" // CHAR(0)
  
  call HDF5WriteDataSetFromVec(string,option, &
                                           global_y_vertex_vec, &
                                           file_id,H5T_NATIVE_DOUBLE)
                                           
  ! Z coord
  string = "Z" // CHAR(0)
  
  call HDF5WriteDataSetFromVec(string,option, &
                                           global_z_vertex_vec, &
                                           file_id,H5T_NATIVE_DOUBLE) 


  call VecDestroy(global_x_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_y_vertex_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_z_vertex_vec,ierr);CHKERRQ(ierr)

  !
  !  Write elements
  !
  
  call GMCreateGMDM(grid,gmdm_element,EIGHT_INTEGER,option)
  call GMGridDMCreateVectorElem(grid,gmdm_element,global_vec, &
                            GLOBAL,option) 
  call GMGridDMCreateVectorElem(grid,gmdm_element,natural_vec, &
                            NATURAL,option) 
  call OutputGetCellVerticesGeomech(grid,global_vec)
  call VecScatterBegin(gmdm_element%scatter_gton_elem,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(gmdm_element%scatter_gton_elem,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  
  local_size = grid%nlmax_elem

  vert_count=0
  do i=1,local_size*EIGHT_INTEGER
    if (int(vec_ptr(i)) >0 ) vert_count=vert_count+1
  enddo
  vert_count=vert_count+grid%nlmax_elem

  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = vert_count
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)

  call MPI_Allreduce(vert_count,dims(1),ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)
  geomech_realization%output_option%xmf_vert_len=int(dims(1))

  ! file space which is a 2D block
  rank_mpi = 1
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)

  string = "Cells" // CHAR(0)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(file_id,string,data_set_id,hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON,hdf5_err)
  if (hdf5_flag < 0) then
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(file_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif

  call h5pclose_f(prop_id,hdf5_err)

  istart = 0
  call MPI_Exscan(vert_count, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  start(1) = istart
  length(1) = vert_count
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)

    ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif

  allocate(int_array(vert_count))

  vert_count=0
  do i=1,local_size
    nverts=0
    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) nverts=nverts+1
    enddo
    vert_count=vert_count+1
    select case (nverts)
      case (4) ! Tetrahedron
        int_array(vert_count) = TET_ID_XDMF
      case (5) ! Pyramid
        int_array(vert_count) = PYR_ID_XDMF
      case (6) ! Wedge
        int_array(vert_count) = WED_ID_XDMF
      case (8) ! Hexahedron
        int_array(vert_count) = HEX_ID_XDMF
    end select

    do j=1,8
      if (vec_ptr((i-1)*8+j)>0) then
        vert_count=vert_count+1
        int_array(vert_count) = INT(vec_ptr((i-1)*8+j))-1
      endif
    enddo
  enddo

  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,int_array,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)

  deallocate(int_array)
  call h5pclose_f(prop_id,hdf5_err)

  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  call GMDMDestroy(gmdm_element)
                                  
#endif
!if defined(SCORPIO_WRITE)

end subroutine WriteHDF5CoordinatesXDMFGeomech
#endif
! defined(PETSC_HAVE_HDF5)

end module Output_Geomechanics_module
