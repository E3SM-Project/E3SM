module Output_Common_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Logging_module 
  use Output_Aux_module
  
  !note: only realization_base_type can be used throughout this module.
  use Realization_Base_class, only : realization_base_type

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, save, public :: max_local_size_saved = -1

  !geh: would prefer that this be local to Output_Tecplot_module, but needed
  !     in Output_Surface_module
  PetscInt, parameter, public :: TECPLOT_INTEGER = 0
  PetscInt, parameter, public :: TECPLOT_REAL = 1  
  
  public :: OutputCommonInit, &
            OutputGetVariableArray, &
            OutputGetVariableAtCell, &
            OutputGetVariableAtCoord, &
            OutputGetCellCenteredVelocities, &
            OutputConvertArrayToNatural, &
            OutputGetCellCoordinates, &
            OutputGetVertexCoordinates, &
            OutputFilenameID, &
            OutputFilename, &
            OutputGetCellVertices, &
            OutputXMFHeader, &
            OutputXMFAttribute, &
            OutputXMFFooter, &
            OutputGetFaceVelUGrid, &
            OutputGetFaceFlowrateUGrid, &
            OutputGetExplicitFlowrates, &
            OutputGetCellVerticesExplicit, &
!            OutputXMFHeaderExplicit, &
!            OutputXMFAttributeExplicit, &
            OutputGetExplicitIDsFlowrates, &
            OutputGetExplicitAuxVars, &
            OutputGetExplicitCellInfo, &
            OutputCollectVelocityOrFlux
              
contains

! ************************************************************************** !

subroutine OutputCommonInit()
  ! 
  ! Initializes module variables for common formats
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/13
  ! 

  use Option_module

  implicit none
  
  ! set size to -1 in order to re-initialize parallel communication blocks
  max_local_size_saved = -1

end subroutine OutputCommonInit

! ************************************************************************** !

function OutputFilenameID(output_option,option)
  ! 
  ! Creates an ID for filename
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/12
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option

  character(len=MAXWORDLENGTH) :: OutputFilenameID
  
  if (output_option%plot_number < 10) then
    write(OutputFilenameID,'("00",i1)') output_option%plot_number
  else if (output_option%plot_number < 100) then
    write(OutputFilenameID,'("0",i2)') output_option%plot_number  
  else if (output_option%plot_number < 1000) then
    write(OutputFilenameID,'(i3)') output_option%plot_number  
  else if (output_option%plot_number < 10000) then
    write(OutputFilenameID,'(i4)') output_option%plot_number  
  else if (output_option%plot_number < 100000) then
    write(OutputFilenameID,'(i5)') output_option%plot_number  
  else
    option%io_buffer = 'Plot number exceeds current maximum of 10^5. &
      &Email pflotran-dev@googlegroups.com and ask for a higher maximum.'
    call printErrMsg(option)
  endif 
  
  OutputFilenameID = adjustl(OutputFilenameID)

end function OutputFilenameID

! ************************************************************************** !

function OutputFilename(output_option,option,suffix,optional_string)
  ! 
  ! Creates a filename for a Tecplot file
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/12
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option
  character(len=*) :: suffix
  character(len=*) :: optional_string
  
  character(len=MAXSTRINGLENGTH) :: OutputFilename

  character(len=MAXWORDLENGTH) :: final_suffix
  character(len=MAXSTRINGLENGTH) :: final_optional_string


  if (len_trim(optional_string) > 0) then
    final_optional_string = '-' // optional_string
  else
    final_optional_string = ''
  endif
  final_suffix = '.' // suffix
  
  ! open file
  if (len_trim(output_option%plot_name) > 2) then
    OutputFilename = trim(output_option%plot_name) // &
            trim(final_optional_string) // &
            final_suffix
  else  
    OutputFilename = trim(option%global_prefix) // &
            trim(option%group_prefix) // &
            trim(final_optional_string) // &
            '-' // &
            trim(OutputFilenameID(output_option,option)) // &
            final_suffix
  endif
  
end function OutputFilename

! ************************************************************************** !

subroutine OutputGetVariableArray(realization_base,vec,variable)
  ! 
  ! Extracts variables indexed by ivar from a multivar array
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : RealizationGetVariable

  implicit none

  class(realization_base_type) :: realization_base
  Vec :: vec
  type(output_variable_type) :: variable
  
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_output_get_var_from_array, &
                          ierr);CHKERRQ(ierr)
                        
  call RealizationGetVariable(realization_base,vec,variable%ivar, &
                              variable%isubvar,variable%isubsubvar)

  call PetscLogEventEnd(logging%event_output_get_var_from_array, &
                        ierr);CHKERRQ(ierr)
  
end subroutine OutputGetVariableArray

! ************************************************************************** !

subroutine OutputConvertArrayToNatural(indices,array,local_size,global_size,option)
  ! 
  ! Converts an array  to natural ordering
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none

  PetscInt :: local_size, global_size
  PetscInt :: indices(:)
  PetscReal, pointer :: array(:)
  type(option_type) :: option
  
  Vec :: natural_vec
  PetscInt, allocatable :: indices_zero_based(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  call VecCreate(option%mycomm,natural_vec,ierr);CHKERRQ(ierr)
  call VecSetSizes(natural_vec,PETSC_DECIDE,global_size,ierr);CHKERRQ(ierr)
  call VecSetType(natural_vec,VECMPI,ierr);CHKERRQ(ierr)

  allocate(indices_zero_based(local_size))
  indices_zero_based(1:local_size) = indices(1:local_size)-1

  call VecSetValues(natural_vec,local_size,indices_zero_based, &
                    array,INSERT_VALUES,ierr);CHKERRQ(ierr)

  call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)

  call VecGetLocalSize(natural_vec,local_size,ierr);CHKERRQ(ierr)
  deallocate(array)
  allocate(array(local_size))
  
  call VecGetArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)
  array(1:local_size) = vec_ptr(1:local_size)
  call VecRestoreArrayF90(natural_vec,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  
end subroutine OutputConvertArrayToNatural

! ************************************************************************** !

function OutputGetVariableAtCell(realization_base,ghosted_id,variable)
  ! 
  ! Extracts variables indexed by ivar from a multivar array
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use Realization_Base_class, only : RealizGetVariableValueAtCell
  use Grid_module
  use Option_module

  implicit none
  
  PetscReal :: OutputGetVariableAtCell
  class(realization_base_type) :: realization_base
  PetscInt :: ghosted_id
  type(output_variable_type) :: variable
 
  OutputGetVariableAtCell = &
    RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                 variable%ivar,variable%isubvar, &
                                 variable%isubsubvar)

end function OutputGetVariableAtCell

! ************************************************************************** !

function OutputGetVariableAtCoord(realization_base,variable,x,y,z, &
                                  num_cells,ghosted_ids)
  ! 
  ! Extracts variables indexed by ivar from a multivar array
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use Realization_Base_class, only : RealizGetVariableValueAtCell
  use Grid_module
  use Option_module

  implicit none
  
  PetscReal :: OutputGetVariableAtCoord
  class(realization_base_type) :: realization_base
  type(output_variable_type) :: variable
  PetscReal :: x,y,z
  PetscInt :: num_cells
  PetscInt :: ghosted_ids(num_cells)

  type(grid_type), pointer :: grid
  PetscInt :: icell
  PetscInt :: ghosted_id
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: isubsubvar
  PetscReal :: dx, dy, dz
  PetscReal :: value, sum_value
  PetscReal :: weight, sum_weight, sum_root
  
  sum_value = 0.d0
  sum_weight = 0.d0
  
  grid => realization_base%patch%grid
  
  ivar = variable%ivar
  isubvar = variable%isubvar
  isubsubvar = variable%isubsubvar

  do icell=1, num_cells
    ghosted_id = ghosted_ids(icell)
    dx = x-grid%x(ghosted_id)
    dy = y-grid%y(ghosted_id)
    dz = z-grid%z(ghosted_id)
    sum_root = sqrt(dx*dx+dy*dy+dz*dz)
    value = 0.d0
    value = RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                         ivar,isubvar,isubsubvar)
    if (sum_root < 1.d-40) then ! bail because it is right on this coordinate
      sum_weight = 1.d0
      sum_value = value
      exit
    endif
    weight = 1.d0/sum_root
    sum_weight = sum_weight + weight
    sum_value = sum_value + weight * value
  enddo
  
  OutputGetVariableAtCoord = sum_value/sum_weight

end function OutputGetVariableAtCoord

! ************************************************************************** !

subroutine OutputGetCellCenteredVelocities(realization_base,vec_x,vec_y, &
                                           vec_z, iphase)
  ! 
  ! Computes the cell-centered velocity component
  ! as an averages of cell face velocities
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07; refactored 01/31/14
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Logging_module
  use Patch_module
  use Grid_module

  implicit none

  class(realization_base_type) :: realization_base
  Vec :: vec_x,vec_y,vec_z
  PetscInt :: direction
  PetscInt :: iphase
  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_x_ptr(:),vec_y_ptr(:),vec_z_ptr(:)
  PetscReal, allocatable :: velocities(:,:)
  
  call PetscLogEventBegin(logging%event_output_get_cell_vel, &
                          ierr);CHKERRQ(ierr)
                            
  allocate(velocities(3,realization_base%patch%grid%nlmax))
  call PatchGetCellCenteredVelocities(realization_base%patch,iphase,velocities)

  call VecGetArrayF90(vec_x,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(vec_y,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(vec_z,vec_z_ptr,ierr);CHKERRQ(ierr)

  vec_x_ptr(:) = velocities(X_DIRECTION,:)*realization_base%output_option%tconv
  vec_y_ptr(:) = velocities(Y_DIRECTION,:)*realization_base%output_option%tconv
  vec_z_ptr(:) = velocities(Z_DIRECTION,:)*realization_base%output_option%tconv

  call VecRestoreArrayF90(vec_x,vec_x_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(vec_y,vec_y_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(vec_z,vec_z_ptr,ierr);CHKERRQ(ierr)

  deallocate(velocities)
  
  call PetscLogEventEnd(logging%event_output_get_cell_vel,ierr);CHKERRQ(ierr)

end subroutine OutputGetCellCenteredVelocities

! ************************************************************************** !

subroutine OutputGetCellCoordinates(grid,vec,direction)
  ! 
  ! Extracts coordinates of cells into a PetscVec
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Variables_module
  
  implicit none

  type(grid_type) :: grid
  Vec :: vec
  PetscInt :: direction
  PetscErrorCode :: ierr
  
  PetscInt :: i
  PetscReal, pointer :: vec_ptr(:)
  
  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  
  if (direction == X_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%x(grid%nL2G(i))
    enddo
  else if (direction == Y_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%y(grid%nL2G(i))
    enddo
  else if (direction == Z_COORDINATE) then
    do i = 1,grid%nlmax
      vec_ptr(i) = grid%z(grid%nL2G(i))
    enddo
  endif
  
  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  
end subroutine OutputGetCellCoordinates

! ************************************************************************** !

subroutine OutputGetVertexCoordinates(grid,vec,direction,option)
  ! 
  ! Extracts vertex coordinates of cells into a PetscVec
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/01/2011
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Option_module
  use Variables_module, only : X_COORDINATE, Y_COORDINATE, Z_COORDINATE
  
  implicit none

  type(grid_type) :: grid
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
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          vec_ptr(ivertex) = grid%unstructured_grid%vertices(ivertex)%z
        enddo
    end select
    call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  else
    ! initialize to UNINITIALIZED_INTEGER to catch bugs
    call VecSet(vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    allocate(values(grid%unstructured_grid%num_vertices_local))
    allocate(indices(grid%unstructured_grid%num_vertices_local))
    select case(direction)
      case(X_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%x
        enddo
      case(Y_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%y
        enddo
      case(Z_COORDINATE)
        do ivertex = 1,grid%unstructured_grid%num_vertices_local
          values(ivertex) = grid%unstructured_grid%vertices(ivertex)%z
        enddo
    end select
    indices(:) = grid%unstructured_grid%vertex_ids_natural(:)-1
    call VecSetValues(vec,grid%unstructured_grid%num_vertices_local, &
                      indices,values,INSERT_VALUES,ierr);CHKERRQ(ierr)
    call VecAssemblyBegin(vec,ierr);CHKERRQ(ierr)
    deallocate(values)
    deallocate(indices)
    call VecAssemblyEnd(vec,ierr);CHKERRQ(ierr)
  endif
  
end subroutine OutputGetVertexCoordinates

! ************************************************************************** !

subroutine OutputGetCellVertices(grid, vec)
  ! 
  ! This routine returns a vector containing vertex ids in natural order of
  ! local cells for unstructured grid.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 05/31/12
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module

  implicit none

  type(grid_type) :: grid
  type(grid_unstructured_type),pointer :: ugrid
  Vec :: vec
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: ivertex
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  ugrid => grid%unstructured_grid
  
  call VecGetArrayF90( vec, vec_ptr, ierr);CHKERRQ(ierr)

  
  ! initialize
  vec_ptr = UNINITIALIZED_DOUBLE
  do local_id=1, ugrid%nlmax
    ghosted_id = local_id
    select case(ugrid%cell_type(ghosted_id))
      case(HEX_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 8
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case(WEDGE_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 6
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        vec_ptr(offset + 7) = 0
        vec_ptr(offset + 8) = 0
      case (PYR_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 5
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        do ivertex = 6, 8
          vec_ptr(offset + ivertex) = 0
        enddo
      case (TET_TYPE)
        offset = (local_id-1)*8
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        do ivertex = 5, 8
          vec_ptr(offset + ivertex) = 0
        enddo
      case (QUAD_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 4
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
      case (TRI_TYPE)
        offset = (local_id-1)*4
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
            ugrid%vertex_ids_natural(ugrid%cell_vertices(ivertex,local_id))
        enddo
        ivertex = 4
        vec_ptr(offset + 4) = 0
    end select
  enddo

  call VecRestoreArrayF90( vec, vec_ptr, ierr);CHKERRQ(ierr)

end subroutine OutputGetCellVertices

! ************************************************************************** !

subroutine OutputGetCellVerticesExplicit(grid, vec)
  ! 
  ! returns a vector containing vertex ids in natural order of
  ! local cells for unstructured grid of explicit type
  ! 
  ! Author: Satish Karra
  ! Date: 07/16/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module

  implicit none

  type(grid_type) :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(unstructured_explicit_type), pointer :: explicit_grid
  Vec :: vec
  PetscInt :: offset
  PetscInt :: ivertex, icell
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  
  ugrid => grid%unstructured_grid
  explicit_grid => ugrid%explicit_grid
  
  call VecGetArrayF90( vec, vec_ptr, ierr);CHKERRQ(ierr)

  ! initialize
  vec_ptr = UNINITIALIZED_DOUBLE
  do icell = 1, explicit_grid%num_elems
    select case(explicit_grid%cell_vertices(0,icell))
      case(8)
        offset = (icell-1)*8
        do ivertex = 1, 8
          vec_ptr(offset + ivertex) = &
            explicit_grid%cell_vertices(ivertex,icell)
        enddo
      case(6)
        offset = (icell-1)*8
        do ivertex = 1, 6
          vec_ptr(offset + ivertex) = &
            explicit_grid%cell_vertices(ivertex,icell)
        enddo
        vec_ptr(offset + 7) = 0
        vec_ptr(offset + 8) = 0
      case (5)
        offset = (icell-1)*8
        do ivertex = 1, 5
          vec_ptr(offset + ivertex) = &
            explicit_grid%cell_vertices(ivertex,icell)
        enddo
        do ivertex = 6, 8
          vec_ptr(offset + ivertex) = 0
        enddo
      case (4)
        if (grid%unstructured_grid%grid_type /= TWO_DIM_GRID) then
          offset = (icell-1)*8
          do ivertex = 1, 4
            vec_ptr(offset + ivertex) = &
              explicit_grid%cell_vertices(ivertex,icell)
          enddo
          do ivertex = 5, 8
            vec_ptr(offset + ivertex) = 0
          enddo
        else
          offset = (icell-1)*8
          do ivertex = 1, 4
            vec_ptr(offset + ivertex) = &
              explicit_grid%cell_vertices(ivertex,icell)
          enddo
          do ivertex = 5, 8
            vec_ptr(offset + ivertex) = 0
          enddo          
        endif
      case (3)
        offset = (icell-1)*8
        do ivertex = 1, 3
          vec_ptr(offset + ivertex) = &
           explicit_grid%cell_vertices(ivertex,icell)
        enddo
        do ivertex = 4, 8
          vec_ptr(offset + ivertex) = 0
        enddo
    end select
  enddo

  call VecRestoreArrayF90( vec, vec_ptr, ierr);CHKERRQ(ierr)

end subroutine OutputGetCellVerticesExplicit

! ************************************************************************** !

subroutine OutputXMFHeader(fid,time,nmax,xmf_vert_len,ngvert,filename, &
                           include_cell_centers)
  ! 
  ! This subroutine writes header to a .xmf file
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/29/12
  ! 

  implicit none

  PetscReal :: time
  PetscInt :: fid, nmax, xmf_vert_len, ngvert
  character(len=MAXSTRINGLENGTH) :: filename
  PetscBool :: include_cell_centers

  character(len=MAXSTRINGLENGTH) :: string, string2
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: i
  
  string='<?xml version="1.0" ?>'
  write(fid,'(a)') trim(string)
  
  string='<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(fid,'(a)') trim(string)

  string='<Xdmf>'
  write(fid,'(a)') trim(string)

  string='  <Domain>'
  write(fid,'(a)') trim(string)

  string='    <Grid Name="Mesh">'
  write(fid,'(a)') trim(string)

  write(string2,'(es13.5)') time
  string='      <Time Value = "' // trim(adjustl(string2)) // '" />'
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string='      <Topology Type="Mixed" NumberOfElements="' // &
    trim(adjustl(string2)) // '">'
  write(fid,'(a)') trim(string)

  write(string2,*) xmf_vert_len
  string='        <DataItem Format="HDF" DataType="Int" Dimensions="' // &
    trim(adjustl(string2)) // '">'
  write(fid,'(a)') trim(string)

  string='          ' // trim(filename) // ':/Domain/Cells'
  write(fid,'(a)') trim(string)

  string='        </DataItem>'
  write(fid,'(a)') trim(string)

  string='      </Topology>'
  write(fid,'(a)') trim(string)

  string='      <Geometry GeometryType="XYZ">'
  write(fid,'(a)') trim(string)

  write(string2,*) ngvert
  string='        <DataItem Format="HDF" Dimensions="' // &
         trim(adjustl(string2)) // ' 3">'
  write(fid,'(a)') trim(string)

  string='          ' // trim(filename) // ':/Domain/Vertices'
  write(fid,'(a)') trim(string)

  string='        </DataItem>'
  write(fid,'(a)') trim(string)

  string="      </Geometry>"
  write(fid,'(a)') trim(string)

  if (include_cell_centers) then

    do i = 1, 3
      select case(i)
        case(1)
          word = 'XC'
        case(2)
          word = 'YC'
        case(3)
          word = 'ZC'
      end select
 
      string='      <Attribute Name="' // trim(word) // &
             '" AttributeType="Scalar"  Center="Cell">'
      write(fid,'(a)') trim(string)
    
      write(string2,*) nmax
      string='        <DataItem Dimensions="' // trim(adjustl(string2)) // &
             ' 1" Format="HDF"> '
      write(fid,'(a)') trim(string)
    
      string='          ' // trim(filename) // ':/Domain/' // trim(word)
      write(fid,'(a)') trim(string)
    
      string='        </DataItem> ' 
      write(fid,'(a)') trim(string)
    
      string='      </Attribute>'
      write(fid,'(a)') trim(string)
    enddo

  endif
    
end subroutine OutputXMFHeader

! ************************************************************************** !

subroutine OutputXMFFooter(fid)
  ! 
  ! This subroutine writes footer to a .xmf file
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 10/29/12
  ! 

  implicit none

  PetscInt :: fid

  character(len=MAXSTRINGLENGTH) :: string

  string='    </Grid>'
  write(fid,'(a)') trim(string)

  string='  </Domain>'
  write(fid,'(a)') trim(string)

  string='</Xdmf>'
  write(fid,'(a)') trim(string)

end subroutine OutputXMFFooter

! ************************************************************************** !

subroutine OutputXMFAttribute(fid,nmax,attname,att_datasetname,mesh_type)
  ! 
  ! Header for xdmf attribute with explicit grid
  ! 
  ! Author: Gautam Bisht, Satish Karra, Glenn Hammond
  ! Date: 10/29/12, 07/17/13, 03/06/17
  ! 
  implicit none

  PetscInt :: fid, nmax, mesh_type
  character(len=MAXSTRINGLENGTH) :: attname, att_datasetname
  
  character(len=MAXSTRINGLENGTH) :: string,string2
  character(len=MAXWORDLENGTH) :: mesh_type_word

  if (mesh_type == VERTEX_CENTERED_OUTPUT_MESH) then
    mesh_type_word = 'Node'
  else if (mesh_type == CELL_CENTERED_OUTPUT_MESH) then
    mesh_type_word = 'Cell'
  end if

  string='      <Attribute Name="' // trim(attname) // &
         '" AttributeType="Scalar"  Center="' // trim(mesh_type_word) // '">'
  write(fid,'(a)') trim(string)

  write(string2,*) nmax
  string='        <DataItem Dimensions="' // trim(adjustl(string2)) // &
         ' 1" Format="HDF"> '
  write(fid,'(a)') trim(string)

  string='        ' // trim(att_datasetname)
  write(fid,'(a)') trim(string)

  string='        </DataItem> ' 
  write(fid,'(a)') trim(string)

  string='      </Attribute>'
  write(fid,'(a)') trim(string)

end subroutine OutputXMFAttribute

! ************************************************************************** !

subroutine OutputGetFaceVelUGrid(realization_base)
  ! 
  ! This subroutine saves:
  !  - Face elocities at x/y/z directions, or
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/15/2016
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use HDF5_module
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Variables_module
  use Connection_module
  use Coupler_module
  use HDF5_Aux_module
  use Output_Aux_module
  use Field_module
  
  implicit none

  class(realization_base_type) :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition
  type(ugdm_type),pointer :: ugdm
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: idual
  PetscInt :: iconn
  PetscInt :: face_id
  PetscInt :: local_id_up,local_id_dn
  PetscInt :: ghosted_id_up,ghosted_id_dn
  PetscInt :: iface_up,iface_dn
  PetscInt :: dof
  PetscInt :: sum_connection
  PetscInt :: offset
  PetscInt :: cell_type
  PetscInt :: local_size
  PetscInt :: i
  PetscInt :: iface
  PetscInt :: ndof
  PetscInt :: idx

  PetscReal, pointer :: flowrates(:,:,:)
  PetscReal, pointer :: vx(:,:,:)
  PetscReal, pointer :: vy(:,:,:)
  PetscReal, pointer :: vz(:,:,:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal, pointer :: vec_ptr3(:)
  PetscReal, pointer :: vx_ptr(:)
  PetscReal, pointer :: vy_ptr(:)
  PetscReal, pointer :: vz_ptr(:)
  PetscReal, pointer :: double_array(:)
  PetscReal :: vel_vector(3)
  PetscReal :: dtime

  Vec :: natural_flowrates_vec
  Vec :: natural_vx_vec
  Vec :: natural_vy_vec
  Vec :: natural_vz_vec

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: unit_string

  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  output_option =>realization_base%output_option
  option => realization_base%option
  field => realization_base%field

  ! Create UGDM for
  call UGridCreateUGDM(grid%unstructured_grid,ugdm, &
                       (option%nflowspec*MAX_FACE_PER_CELL + 1),option)

  ! Create vectors in natural order for velocity in x/y/z direction
  call UGridDMCreateVector(grid%unstructured_grid,ugdm,natural_vx_vec, &
                           NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm,natural_vy_vec, &
                           NATURAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm,natural_vz_vec, &
                           NATURAL,option)

  allocate(vx(option%nflowspec,MAX_FACE_PER_CELL,ugrid%nlmax))
  allocate(vy(option%nflowspec,MAX_FACE_PER_CELL,ugrid%nlmax))
  allocate(vz(option%nflowspec,MAX_FACE_PER_CELL,ugrid%nlmax))

  vx = 0.d0
  vy = 0.d0
  vz = 0.d0

  call VecGetArrayF90(field%vx_face_inst,vx_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%vy_face_inst,vy_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%vz_face_inst,vz_ptr,ierr);CHKERRQ(ierr)

  vx_ptr = 0.d0
  vy_ptr = 0.d0
  vz_ptr = 0.d0

  offset = 1 + option%nflowspec*MAX_FACE_PER_CELL

  ! Save the number of faces of all cell
  do local_id = 1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    cell_type = ugrid%cell_type(ghosted_id)

    vx_ptr((local_id-1)*offset+1) = UCellGetNFaces(cell_type,option)
    vy_ptr((local_id-1)*offset+1) = UCellGetNFaces(cell_type,option)
    vz_ptr((local_id-1)*offset+1) = UCellGetNFaces(cell_type,option)
  enddo

  ! Interior Flowrates Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      face_id = cur_connection_set%face_id(iconn)
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      do iface_up = 1,MAX_FACE_PER_CELL
        if (face_id==ugrid%cell_to_face_ghosted(iface_up,local_id_up)) exit
      enddo

      iface_dn=-1
      if (local_id_dn>0) then
        do iface_dn = 1,MAX_FACE_PER_CELL
          if (face_id==ugrid%cell_to_face_ghosted(iface_dn,local_id_dn)) exit
        enddo
      endif

      do dof=1,option%nflowspec

        ! Save velocity for iface_up of local_id_up cell using flowrate up-->dn
        vel_vector = cur_connection_set%dist(1:3,iconn)* &
                     patch%internal_velocities(dof,sum_connection)

        vx(dof,iface_up,local_id_up) = vel_vector(1)
        vy(dof,iface_up,local_id_up) = vel_vector(2)
        vz(dof,iface_up,local_id_up) = vel_vector(3)

        idx = (local_id_up-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_up + 1

        vx_ptr(idx) = vel_vector(1)
        vy_ptr(idx) = vel_vector(2)
        vz_ptr(idx) = vel_vector(3)

        if (iface_dn>0) then

          ! Save velocity for iface_dn of local_id_dn cell using -ve flowrate up-->dn

          idx = (local_id_dn-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_dn + 1

          vx(dof,iface_dn,local_id_dn) = -vel_vector(1)
          vy(dof,iface_dn,local_id_dn) = -vel_vector(2)
          vz(dof,iface_dn,local_id_dn) = -vel_vector(3)

          vx_ptr(idx) = -vel_vector(1)
          vx_ptr(idx) = -vel_vector(2)
          vx_ptr(idx) = -vel_vector(3)

        endif

      enddo ! dof-loop

    enddo ! iconn-loop

    cur_connection_set => cur_connection_set%next

  enddo

  ! Boundary Flowrates Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do 
    if (.not.associated(boundary_condition)) exit

    cur_connection_set => boundary_condition%connection_set

    do iconn = 1, cur_connection_set%num_connections

      sum_connection = sum_connection + 1
      face_id = cur_connection_set%face_id(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_dn = grid%nG2L(ghosted_id_dn)

      do iface_dn = 1,MAX_FACE_PER_CELL
        if (face_id==ugrid%cell_to_face_ghosted(iface_dn,local_id_dn)) exit
      enddo

      do dof=1,option%nflowspec

        ! Save velocity for iface_dn of local_id_dn cell using -ve flowrate up-->dn
        vel_vector = cur_connection_set%dist(1:3,iconn)* &
                     patch%boundary_velocities(dof,sum_connection)

        idx = (local_id_dn-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_dn + 1

        vx(dof,iface_dn,local_id_dn) = -vel_vector(1)
        vy(dof,iface_dn,local_id_dn) = -vel_vector(2)
        vz(dof,iface_dn,local_id_dn) = -vel_vector(3)

        vx_ptr(idx) = -vel_vector(1)
        vx_ptr(idx) = -vel_vector(2)
        vx_ptr(idx) = -vel_vector(3)

      enddo ! dof-loop

    enddo ! iconn-loop

    boundary_condition => boundary_condition%next

  enddo

  deallocate(vx)
  deallocate(vy)
  deallocate(vz)

  call VecRestoreArrayF90(field%vx_face_inst,vx_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%vy_face_inst,vy_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%vz_face_inst,vz_ptr,ierr);CHKERRQ(ierr)

  ! Scatter flowrate from Global --> Natural order
  call VecScatterBegin(ugdm%scatter_gton,field%vx_face_inst,natural_vx_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm%scatter_gton,field%vx_face_inst,natural_vx_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  call VecScatterBegin(ugdm%scatter_gton,field%vy_face_inst,natural_vy_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm%scatter_gton,field%vy_face_inst,natural_vy_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  call VecScatterBegin(ugdm%scatter_gton,field%vz_face_inst,natural_vz_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm%scatter_gton,field%vz_face_inst,natural_vz_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  ! X-direction
  call VecGetArrayF90(natural_vx_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%vx_face_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  ! Copy the vectors
  vec_ptr2 = vec_ptr
  call VecRestoreArrayF90(natural_vx_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%vx_face_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  ! Y-direction
  call VecGetArrayF90(natural_vy_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%vy_face_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  ! Copy the vectors
  vec_ptr2 = vec_ptr
  call VecRestoreArrayF90(natural_vy_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%vy_face_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  ! Z-direction
  call VecGetArrayF90(natural_vz_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%vz_face_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  ! Copy the vectors
  vec_ptr2 = vec_ptr
  call VecRestoreArrayF90(natural_vz_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%vz_face_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  call VecDestroy(natural_vx_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vy_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vz_vec,ierr);CHKERRQ(ierr)

  call UGridDMDestroy(ugdm)
  
end subroutine OutputGetFaceVelUGrid

! ************************************************************************** !

subroutine OutputGetFaceFlowrateUGrid(realization_base)
  ! 
  ! This subroutine saves:
  !  - Mass/energy flowrate at all faces of a control volume
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/15/2016
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use HDF5_module
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Variables_module
  use Connection_module
  use Coupler_module
  use HDF5_Aux_module
  use Output_Aux_module
  use Field_module
  
  implicit none

  class(realization_base_type) :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition
  type(ugdm_type),pointer :: ugdm
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: idual
  PetscInt :: iconn
  PetscInt :: face_id
  PetscInt :: local_id_up,local_id_dn
  PetscInt :: ghosted_id_up,ghosted_id_dn
  PetscInt :: iface_up,iface_dn
  PetscInt :: dof
  PetscInt :: sum_connection
  PetscInt :: offset
  PetscInt :: cell_type
  PetscInt :: local_size
  PetscInt :: i
  PetscInt :: iface
  PetscInt :: ndof
  PetscInt :: idx

  PetscReal, pointer :: flowrates(:,:,:)
  PetscReal, pointer :: vx(:,:,:)
  PetscReal, pointer :: vy(:,:,:)
  PetscReal, pointer :: vz(:,:,:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal, pointer :: vec_ptr3(:)
  PetscReal, pointer :: vx_ptr(:)
  PetscReal, pointer :: vy_ptr(:)
  PetscReal, pointer :: vz_ptr(:)
  PetscReal, pointer :: double_array(:)
  PetscReal :: vel_vector(3)
  PetscReal :: dtime

  Vec :: natural_flowrates_vec
  Vec :: natural_vx_vec
  Vec :: natural_vy_vec
  Vec :: natural_vz_vec

  PetscMPIInt :: hdf5_err
  PetscErrorCode :: ierr
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: unit_string

  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  output_option =>realization_base%output_option
  option => realization_base%option
  field => realization_base%field

  ! Create UGDM for
  call UGridCreateUGDM(grid%unstructured_grid,ugdm, &
                       (option%nflowdof*MAX_FACE_PER_CELL + 1),option)

  ! Create a flowrate vector in natural order
  call UGridDMCreateVector(grid%unstructured_grid,ugdm,natural_flowrates_vec, &
                           NATURAL,option)

  allocate(flowrates(option%nflowdof,MAX_FACE_PER_CELL,ugrid%nlmax))
  flowrates = 0.d0

  call VecGetArrayF90(field%flowrate_inst,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr = 0.d0

  offset = 1 + option%nflowdof*MAX_FACE_PER_CELL
  ! Save the number of faces of all cell
  do local_id = 1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    cell_type = ugrid%cell_type(ghosted_id)
    vec_ptr((local_id-1)*offset+1) = UCellGetNFaces(cell_type,option)
  enddo

  ! Interior Flowrates Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      face_id = cur_connection_set%face_id(iconn)
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)

      do iface_up = 1,MAX_FACE_PER_CELL
        if (face_id==ugrid%cell_to_face_ghosted(iface_up,local_id_up)) exit
      enddo

      iface_dn=-1
      if (local_id_dn>0) then
        do iface_dn = 1,MAX_FACE_PER_CELL
          if (face_id==ugrid%cell_to_face_ghosted(iface_dn,local_id_dn)) exit
        enddo
      endif

      do dof=1,option%nflowdof

        ! Save flowrate for iface_up of local_id_up cell using flowrate up-->dn
        flowrates(dof,iface_up,local_id_up) = &
          patch%internal_flow_fluxes(dof,sum_connection)

        idx = (local_id_up-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_up + 1
        vec_ptr(idx) = patch%internal_flow_fluxes(dof,sum_connection)

        if (iface_dn>0) then
          ! Save flowrate for iface_dn of local_id_dn cell using -ve flowrate up-->dn
          flowrates(dof,iface_dn,local_id_dn) = -patch%internal_flow_fluxes(dof,sum_connection)

          idx = (local_id_dn-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_dn + 1
          vec_ptr(idx) = -patch%internal_flow_fluxes(dof,sum_connection)
        endif

      enddo ! dof-loop

    enddo ! iconn-loop

    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flowrates Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do 
    if (.not.associated(boundary_condition)) exit

    cur_connection_set => boundary_condition%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      face_id = cur_connection_set%face_id(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      do iface_dn = 1,MAX_FACE_PER_CELL
        if (face_id==ugrid%cell_to_face_ghosted(iface_dn,local_id_dn)) exit
      enddo

      do dof=1,option%nflowdof

        ! Save flowrate for iface_dn of local_id_dn cell using -ve flowrate up-->dn
        idx = (local_id_dn-1)*offset + (dof-1)*MAX_FACE_PER_CELL + iface_dn + 1
        flowrates(dof,iface_dn,local_id_dn) = &
          -patch%boundary_flow_fluxes(dof,sum_connection)
        vec_ptr(idx) = &
          -patch%boundary_flow_fluxes(dof,sum_connection)

      enddo ! dof-loop

    enddo ! iconn-loop

    boundary_condition => boundary_condition%next
  enddo

  deallocate(flowrates)
  call VecRestoreArrayF90(field%flowrate_inst,vec_ptr,ierr);CHKERRQ(ierr)

  ! Scatter flowrate from Global --> Natural order
  call VecScatterBegin(ugdm%scatter_gton,field%flowrate_inst,natural_flowrates_vec, &
                        INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm%scatter_gton,field%flowrate_inst,natural_flowrates_vec, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(natural_flowrates_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flowrate_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  ! Copy the vectors
  vec_ptr2 = vec_ptr

  if (output_option%print_hdf5_aveg_mass_flowrate.or. &
    output_option%print_hdf5_aveg_energy_flowrate) then

    dtime = option%time-output_option%aveg_var_time
    call VecGetArrayF90(field%flowrate_aveg,vec_ptr3,ierr);CHKERRQ(ierr)
    vec_ptr3 = vec_ptr3 + vec_ptr2/dtime
    call VecRestoreArrayF90(field%flowrate_aveg,vec_ptr3,ierr);CHKERRQ(ierr)
  endif

  call VecRestoreArrayF90(natural_flowrates_vec,vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flowrate_inst,vec_ptr2,ierr);CHKERRQ(ierr)

  call VecDestroy(natural_flowrates_vec,ierr);CHKERRQ(ierr)

  call UGridDMDestroy(ugdm)
  
end subroutine OutputGetFaceFlowrateUGrid

! ************************************************************************** !

subroutine OutputGetExplicitIDsFlowrates(realization_base,count,vec_proc, &
                                         ids_up,ids_dn)
  ! 
  ! Calculates the ids of the nodes of a
  ! connection for flow rats output
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 04/24/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Field_module
  use Connection_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(field_type), pointer :: field
  type(ugdm_type), pointer :: ugdm  
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal, pointer :: vec_proc_ptr(:)
  PetscInt, pointer :: ids_up(:),ids_dn(:)
  PetscInt :: offset
  PetscInt :: istart,iend
  PetscInt :: iconn
  PetscErrorCode :: ierr
  PetscReal :: val
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscReal :: proc_up, proc_dn, conn_proc
  PetscInt :: sum_connection, count
  Vec :: global_vec
  Vec :: local_vec
  Vec :: vec_proc
  PetscInt :: idof
  
  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  option => realization_base%option
  field => realization_base%field
  
  call VecCreateMPI(option%mycomm, &
                    size(grid%unstructured_grid%explicit_grid%connections,2), &
                    PETSC_DETERMINE,vec_proc,ierr);CHKERRQ(ierr)
  call VecSet(vec_proc,0.d0,ierr);CHKERRQ(ierr)
  
  call UGridCreateUGDM(grid%unstructured_grid,ugdm,ONE_INTEGER,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm,global_vec, &
                           GLOBAL,option)
  call UGridDMCreateVector(grid%unstructured_grid,ugdm,local_vec, &
                           LOCAL,option)
  call VecGetArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr = option%myrank
  call VecRestoreArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
  
  call VecScatterBegin(ugdm%scatter_gtol,global_vec,local_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(ugdm%scatter_gtol,global_vec,local_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
                     
  call VecGetArrayF90(local_vec,vec_ptr2,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)
  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  count = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      proc_up = vec_ptr2(ghosted_id_up)
      proc_dn = vec_ptr2(ghosted_id_dn)
      proc_up = min(option%myrank,int(proc_up))
      proc_dn = min(option%myrank,int(proc_dn))
      conn_proc = min(proc_up,proc_dn)
      vec_proc_ptr(sum_connection) = conn_proc
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  call VecRestoreArrayF90(local_vec,vec_ptr2,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)

  call VecAssemblyBegin(vec_proc,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(vec_proc,ierr);CHKERRQ(ierr)
  
  ! Count the number of connections on a local process
  call VecGetArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  count = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      if (option%myrank == int(vec_proc_ptr(sum_connection))) &
        count = count + 1
    enddo
    cur_connection_set => cur_connection_set%next
  enddo  
  call VecRestoreArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)
  

  ! Count the number of connections on a local process
  allocate(ids_up(count))
  allocate(ids_dn(count))
  call VecGetArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  count = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)      
      if (option%myrank == int(vec_proc_ptr(sum_connection))) then
        count = count + 1
        ids_up(count) = grid%nG2A(ghosted_id_up)
        ids_dn(count) = grid%nG2A(ghosted_id_dn)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo  
  call VecRestoreArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)

end subroutine OutputGetExplicitIDsFlowrates

! ************************************************************************** !

subroutine OutputGetExplicitFlowrates(realization_base,count,vec_proc, &
                                      flowrates,darcy,area)
  ! 
  ! Forms a vector of magnitude of flowrates
  ! which will be printed out to file for particle tracking.
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 04/24/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Field_module
  use Connection_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(field_type), pointer :: field
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  
  PetscReal, pointer :: vec_proc_ptr(:)
  PetscReal, pointer :: flowrates(:,:)
  PetscReal, pointer :: darcy(:), area(:)
  PetscInt :: offset
  PetscInt :: iconn
  PetscErrorCode :: ierr
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscReal :: proc_up, proc_dn, conn_proc
  PetscInt :: sum_connection, count
  Vec :: vec_proc
  PetscInt :: idof
  
  patch => realization_base%patch
  grid => patch%grid
  ugrid => grid%unstructured_grid
  option => realization_base%option
  field => realization_base%field

  ! Count the number of connections on a local process
  allocate(flowrates(count,option%nflowdof))
  allocate(darcy(count))
  allocate(area(count))
  call VecGetArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  count = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)      
      if (option%myrank == int(vec_proc_ptr(sum_connection))) then
        count = count + 1
        do idof = 1,option%nflowdof
          flowrates(count,option%nflowdof) = &
            patch%internal_flow_fluxes(idof,sum_connection)
        enddo
        darcy(count) = patch%internal_velocities(1,sum_connection)
        area(count) = cur_connection_set%area(iconn)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo  
  call VecRestoreArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)

end subroutine OutputGetExplicitFlowrates

! ************************************************************************** !

subroutine OutputGetExplicitAuxVars(realization_base,count,vec_proc,density)
  ! 
  ! Calculates density at the face
  ! between a connection
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/17/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Field_module
  use Connection_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(field_type), pointer :: field
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(global_auxvar_type), pointer :: global_auxvar(:)
  type(material_parameter_type), pointer :: material_parameter


  PetscReal, pointer :: vec_proc_ptr(:)
  PetscReal, pointer :: flowrates(:,:)
  PetscReal, pointer :: darcy(:)
  PetscReal, pointer :: density(:)
  PetscInt :: offset
  PetscInt :: iconn
  PetscErrorCode :: ierr
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscReal :: proc_up, proc_dn, conn_proc
  PetscInt :: sum_connection, count
  Vec :: vec_proc
  PetscInt :: idof
  PetscInt :: icap_up, icap_dn
  PetscReal :: sir_up, sir_dn
  PetscReal, parameter :: eps = 1.D-8
  PetscReal :: upweight

  
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  grid => patch%grid
  global_auxvar => patch%aux%Global%auxvars
  material_parameter => patch%aux%Material%material_parameter
 
  allocate(density(count))
  call VecGetArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  count = 0 
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn) 
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
      if (option%myrank == int(vec_proc_ptr(sum_connection))) then
        count = count + 1
        sir_up = material_parameter%soil_residual_saturation(1,icap_up)
        sir_dn = material_parameter%soil_residual_saturation(1,icap_dn)

        if (global_auxvar(ghosted_id_up)%sat(1) > sir_up .or. &
            global_auxvar(ghosted_id_dn)%sat(1) > sir_dn) then
          if (global_auxvar(ghosted_id_up)%sat(1) <eps) then 
            upweight = 0.d0
          else if (global_auxvar(ghosted_id_dn)%sat(1) <eps) then 
            upweight = 1.d0
          endif
    
          density(count) = upweight*global_auxvar(ghosted_id_up)%den(1)+ &
                  (1.D0 - upweight)*global_auxvar(ghosted_id_dn)%den(1)
        endif
      endif  
          
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
      
  call VecRestoreArrayF90(vec_proc,vec_proc_ptr,ierr);CHKERRQ(ierr)


end subroutine OutputGetExplicitAuxVars

! ************************************************************************** !

subroutine OutputGetExplicitCellInfo(realization_base,num_cells,ids,sat,por, &
                                     density,pressure)
  ! 
  ! Calculates porosity, saturation, density
  ! and pressure in a cell (explicit)
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 08/21/13
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Grid_Unstructured_Aux_module
  use Field_module
  use Connection_module
  use Global_Aux_module

  implicit none

  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(grid_unstructured_type),pointer :: ugrid
  type(field_type), pointer :: field
  type(global_auxvar_type), pointer :: global_auxvar(:)


  PetscErrorCode :: ierr
  PetscInt :: num_cells
  PetscReal, pointer :: sat(:)
  PetscReal, pointer :: por(:)
  PetscReal, pointer :: density(:)
  PetscReal, pointer :: pressure(:)
  PetscInt, pointer :: ids(:)
  PetscInt :: local_id, ghosted_id
  
  patch => realization_base%patch
  option => realization_base%option
  field => realization_base%field
  grid => patch%grid
  global_auxvar => patch%aux%Global%auxvars
  
  num_cells = grid%nlmax
  allocate(sat(num_cells))
  allocate(por(num_cells))
  allocate(ids(num_cells))
  allocate(density(num_cells))
  allocate(pressure(num_cells))
  
  do local_id = 1, num_cells
    ghosted_id = grid%nL2G(local_id)
    ids(local_id) = grid%nG2A(ghosted_id)
    sat(local_id) = global_auxvar(ghosted_id)%sat(1)
    por(local_id) = patch%aux%Material%auxvars(ghosted_id)%porosity
    density(local_id) = global_auxvar(ghosted_id)%den(1)
    pressure(local_id) = global_auxvar(ghosted_id)%pres(1)
  enddo

end subroutine OutputGetExplicitCellInfo

! ************************************************************************** !

subroutine OutputCollectVelocityOrFlux(realization_base, iphase, direction, &
                                       output_flux, array)
  ! 
  ! Accumulates fluxes or velocities for a structured grid into a 1D array.
  ! This routine is called for HDF5 and Tecplot flux/velocity output.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/26/18
  ! 
  use Realization_Base_class, only : realization_base_type
  use Discretization_module
  use Patch_module
  use Grid_module
  use Grid_Structured_module
  use Option_module
  use Field_module
  use Connection_module
  use Coupler_module
  use DM_Kludge_module
  
  implicit none

  class(realization_base_type) :: realization_base
  PetscInt :: iphase
  PetscInt :: direction
  PetscBool :: output_flux
  PetscReal :: array(*)

  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(grid_structured_type), pointer :: structured_grid
  type(option_type), pointer :: option
  type(dm_ptr_type), pointer :: dm_ptr

  PetscInt :: local_id, ghosted_id
  PetscInt :: local_size
  PetscInt :: nx_local, ny_local, nz_local
  PetscInt :: i, j, k
  PetscReal :: scale, value, dist(-1:3)
  PetscReal :: face_location_pert
  PetscReal :: max_global
  PetscReal :: min_global
  PetscReal, pointer :: coord_ptr(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, parameter :: perturbation = 1.d-6
  PetscInt :: count, iconn, sum_connection

  Vec :: local_vec
  Vec :: global_vec
  PetscErrorCode :: ierr

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition

  discretization => realization_base%discretization
  patch => realization_base%patch
  grid => patch%grid
  structured_grid => grid%structured_grid
  option => realization_base%option
  field => realization_base%field

  local_size = grid%nlmax
!GEH - Structured Grid Dependence - Begin
  nx_local = structured_grid%nlx
  ny_local = structured_grid%nly
  nz_local = structured_grid%nlz
  select case(direction)
    case(X_DIRECTION)
      if (structured_grid%gxe-structured_grid%lxe == 0) then
        local_size = grid%nlmax-structured_grid%nlyz
        nx_local = structured_grid%nlx-1
      endif
    case(Y_DIRECTION)
      if (structured_grid%gye-structured_grid%lye == 0) then
        local_size = grid%nlmax-structured_grid%nlxz
        ny_local = structured_grid%nly-1
      endif
    case(Z_DIRECTION)
      if (structured_grid%gze-structured_grid%lze == 0) then
        local_size = grid%nlmax-structured_grid%nlxy
        nz_local = structured_grid%nlz-1
      endif
  end select


  ! must use a local vec so that potential boundary values can be 
  ! accumulated across ghosted cells
  call DiscretizationCreateVector(discretization,ONEDOF,local_vec,LOCAL, &
                                  option) 
  call VecZeroEntries(local_vec,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(local_vec,vec_ptr,ierr);CHKERRQ(ierr)
  
  ! place interior velocities in a vector
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id = cur_connection_set%id_up(iconn)
      local_id = grid%nG2L(ghosted_id) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      if (local_id <= 0 .or. &
          dabs(cur_connection_set%dist(direction,iconn)) < 0.99d0) cycle
      if (output_flux) then
        ! iphase here is really the dof
        vec_ptr(ghosted_id) = patch%internal_flow_fluxes(iphase,sum_connection)
      else
        vec_ptr(ghosted_id) = patch%internal_velocities(iphase,sum_connection)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! add contribution of boundary velocities
  select case(direction)
    case(X_DIRECTION)
      coord_ptr => grid%x
      max_global = grid%x_max_global
      min_global = grid%x_min_global
    case(Y_DIRECTION)
      coord_ptr => grid%y
      max_global = grid%y_max_global
      min_global = grid%y_min_global
    case(Z_DIRECTION)
      coord_ptr => grid%z
      max_global = grid%z_max_global
      min_global = grid%z_min_global
  end select
  boundary_condition => patch%boundary_condition_list%first 
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      dist = cur_connection_set%dist(:,iconn)
      if (dabs(dist(direction)) < 0.99d0) cycle
      scale = 1.d0
      if (dist(direction) < 0.d0) scale = -1.d0
      ! if the connection is on the domain boundary, we need to skip it.
      ! use a small perturbation to determine 
      face_location_pert = coord_ptr(ghosted_id) - &
                         (1.d0+perturbation)*dist(0)*dist(direction)
      if (face_location_pert >= max_global .or. &
          face_location_pert <= min_global) then
        cycle
      endif
      ! velocities are stored as the downwind face of the upwind cell.
      ! if the direction is positive, then the value needs to be assigned
      ! to the downwind face of the next cell upwind in the specified 
      ! direction.
      select case(direction)
        case(X_DIRECTION)
          if (scale > 0.d0) ghosted_id = ghosted_id - 1
        case(Y_DIRECTION)
          if (scale > 0.d0) ghosted_id = ghosted_id - structured_grid%ngx
        case(Z_DIRECTION)
          if (scale > 0.d0) ghosted_id = ghosted_id - structured_grid%ngxy
      end select
      if (ghosted_id <= 0) then
        option%io_buffer = 'Negative ghosted id in OutputFluxVelocities&
          &TecplotBlk while adding boundary values. Please contact &
          &pflotran-dev@googlegroups.com with this message.'
        call printErrMsgByRank(option)
      endif
      ! I don't know why one would do this, but it is possible that a 
      ! boundary condition could be applied to an interior face shared
      ! by two active cells. Thus, we must sum.
      if (output_flux) then
        value = patch%boundary_flow_fluxes(iphase,sum_connection)
      else
        value = patch%boundary_velocities(iphase,sum_connection)
      endif
      vec_ptr(ghosted_id) = vec_ptr(ghosted_id) + scale*value
    enddo
    boundary_condition => boundary_condition%next 
  enddo
  call VecRestoreArrayF90(local_vec,vec_ptr,ierr);CHKERRQ(ierr)

  ! sum values across processes
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,ONEDOF)
  ! for a given cell, ghosted values for that cell may only be summed
  ! using DMLocalToGlobalBegin/End with ADD_VALUES. LocalToLocal does not
  ! work
  call DiscretizationCreateVector(discretization,ONEDOF,global_vec,GLOBAL, &
                                  option) 
  call VecZeroEntries(global_vec,ierr);CHKERRQ(ierr)
  call DMLocalToGlobalBegin(dm_ptr%dm,local_vec,ADD_VALUES,global_vec, &
                            ierr);CHKERRQ(ierr)
  call DMLocalToGlobalEnd(dm_ptr%dm,local_vec,ADD_VALUES,global_vec, &
                          ierr);CHKERRQ(ierr)

  call VecGetArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
  ! write out data set 
  count = 0 
  do k=1,nz_local 
    do j=1,ny_local 
      do i=1,nx_local 
        count = count + 1 
        local_id = i+(j-1)*structured_grid%nlx+ &
                   (k-1)*structured_grid%nlxy 
        array(count) = vec_ptr(local_id)
      enddo 
    enddo 
  enddo 
  call VecRestoreArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
   
  call VecDestroy(local_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)

end subroutine OutputCollectVelocityOrFlux

end module Output_Common_module
