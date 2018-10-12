module Matrix_Buffer_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Grid_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: matrix_buffer_type
    PetscInt :: matrix_type
    type(grid_type), pointer :: grid
    PetscInt :: i_width
    PetscInt :: j_width
    PetscInt :: k_width
    PetscInt, pointer :: icol(:,:)
    PetscReal, pointer :: values(:,:)
  end type matrix_buffer_type

  PetscInt, parameter, public :: AIJ = 1
  PetscInt, parameter, public :: HYPRESTRUCT = 2

  public :: MatrixBufferCreate, &
            MatrixBufferInit, &
            MatrixBufferDestroy, &
            MatrixBufferZero, &
            MatrixBufferAdd, &
            MatrixBufferZeroRows, &
            MatrixBufferSetValues
 
contains

! ************************************************************************** !

function MatrixBufferCreate()
  ! 
  ! Creates a matrix object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/13/09
  ! 
  
  implicit none

  type(matrix_buffer_type), pointer :: MatrixBufferCreate
  
  type(matrix_buffer_type), pointer :: matrix_buffer
  
  allocate(matrix_buffer)
  matrix_buffer%matrix_type = 0
  nullify(matrix_buffer%grid)
  matrix_buffer%i_width = 0
  matrix_buffer%j_width = 0
  matrix_buffer%k_width = 0
  nullify(matrix_buffer%values)
  nullify(matrix_buffer%icol)
  MatrixBufferCreate => matrix_buffer

end function MatrixBufferCreate

! ************************************************************************** !

subroutine MatrixBufferInit(A,matrix_buffer,grid)
  ! 
  ! Initializes matrix buffer object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/13/09
  ! 

  use Grid_Structured_module

  implicit none

  Mat :: A
  type(matrix_buffer_type), pointer :: matrix_buffer  
  type(grid_type), target :: grid

  MatType :: mat_type
  PetscInt :: local_id, ghosted_id
  PetscInt :: i, j, k
  type(grid_structured_type), pointer :: structured_grid
  PetscErrorCode :: ierr

  matrix_buffer%grid => grid
  structured_grid => grid%structured_grid

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)

  select case(mat_type)
    case(MATMPIAIJ,MATSEQAIJ,MATMPIBAIJ,MATSEQBAIJ)
      matrix_buffer%matrix_type = AIJ
      matrix_buffer%i_width = 1
      matrix_buffer%j_width = structured_grid%ngx 
      matrix_buffer%k_width = structured_grid%ngxy
      allocate(matrix_buffer%icol(7,structured_grid%ngmax))
      allocate(matrix_buffer%values(7,structured_grid%ngmax))
      ! compute fill indices
      ! initialize all to zero
      matrix_buffer%icol = 0 
      local_id = 0
      do k=structured_grid%kstart,structured_grid%kend
        do j=structured_grid%jstart,structured_grid%jend
          do i=structured_grid%istart,structured_grid%iend
            local_id = local_id + 1
            ghosted_id = grid%nL2G(local_id)
            if (k > 0) then
              matrix_buffer%icol(1,ghosted_id) = ghosted_id - structured_grid%ngxy
            endif
            if (j > 0) then
              matrix_buffer%icol(2,ghosted_id) = ghosted_id - structured_grid%ngx
            endif
            if (i > 0) then
              matrix_buffer%icol(3,ghosted_id) = ghosted_id - 1
            endif
            matrix_buffer%icol(4,ghosted_id) = ghosted_id
            if (i < structured_grid%iend .or. &
                structured_grid%gxe-structured_grid%lxe > 0) then
              matrix_buffer%icol(5,ghosted_id) = ghosted_id + 1
            endif
            if (j < structured_grid%jend .or. &
                structured_grid%gye-structured_grid%lye > 0) then
              matrix_buffer%icol(6,ghosted_id) = ghosted_id + structured_grid%ngx
            endif
            if (k < structured_grid%kend .or. &
                structured_grid%gze-structured_grid%lze > 0) then
              matrix_buffer%icol(7,ghosted_id) = ghosted_id + structured_grid%ngxy
            endif
          enddo
        enddo
      enddo
      ! convert to zero-based
      ! inactive entries will become -1, which is constistent with PETSc
      matrix_buffer%icol = matrix_buffer%icol - 1 
    case(MATHYPRESTRUCT)
      matrix_buffer%matrix_type = HYPRESTRUCT
      matrix_buffer%i_width = 1
      matrix_buffer%j_width = structured_grid%ngx 
      matrix_buffer%k_width = structured_grid%ngxy
      allocate(matrix_buffer%values(7,structured_grid%ngmax))
  end select
  matrix_buffer%values = 0.d0

end subroutine MatrixBufferInit

! ************************************************************************** !

subroutine MatrixBufferZero(matrix_buffer)
  ! 
  ! Zeros matrix buffer values
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/13/09
  ! 

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer  

  matrix_buffer%values = 0.d0

end subroutine MatrixBufferZero

! ************************************************************************** !

subroutine MatrixBufferAdd(matrix_buffer,irow,icol,value)
  ! 
  ! Adds values to matrix buffer object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/13/09
  ! 

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer  
  PetscInt :: irow, icol
  PetscReal :: value
  
  PetscInt :: index

  if (icol < 0) then
    index = 4
  else
    index = irow-icol
    if (index == 0) then
      index = 4
    else if (index == -matrix_buffer%i_width) then
      index = 5
    else if (index == matrix_buffer%i_width) then
      index = 3
    else if (index == -matrix_buffer%j_width) then
      index = 6
    else if (index == matrix_buffer%j_width) then
      index = 2
    else if (index == -matrix_buffer%k_width) then
      index = 7
    else if (index == matrix_buffer%k_width) then
      index = 1
    endif
  endif
  matrix_buffer%values(index,irow) = matrix_buffer%values(index,irow) + value

end subroutine MatrixBufferAdd

! ************************************************************************** !

subroutine MatrixBufferZeroRows(matrix_buffer,nrows,rows)
  ! 
  ! Zeros rows in a matrix buffer object, placing 1
  ! on the diagonal
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/14/09
  ! 

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer  
  PetscInt :: nrows
  PetscInt, pointer :: rows(:)

  PetscInt :: irow, irow_local
  
  select case(matrix_buffer%matrix_type)
    case(AIJ)
      do irow = 1, nrows
        irow_local = matrix_buffer%grid%nG2L(rows(irow))
        matrix_buffer%values(:,irow_local) = 0.d0
        matrix_buffer%values(4,irow_local) = 1.d0
      enddo
    case(HYPRESTRUCT)
      do irow = 1, nrows
        matrix_buffer%values(:,rows(irow)) = 0.d0
        matrix_buffer%values(4,rows(irow)) = 1.d0
      enddo
  end select

end subroutine MatrixBufferZeroRows

! ************************************************************************** !

subroutine MatrixBufferSetValues(A,matrix_buffer)
  ! 
  ! Sets values in PETSc matrix
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/14/09
  ! 

  implicit none

  Mat :: A
  type(matrix_buffer_type), pointer :: matrix_buffer  

  select case(matrix_buffer%matrix_type)
    case(AIJ)
      call MatrixBufferSetValuesAij(A,matrix_buffer)
    case(HYPRESTRUCT)
      call MatrixBufferSetValuesHypre(A,matrix_buffer)
  end select

end subroutine MatrixBufferSetValues

! ************************************************************************** !

subroutine MatrixBufferSetValuesHypre(A,matrix_buffer)
  ! 
  ! Sets values in PETSc Hypre matrix
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/13/09
  ! 

  implicit none

  Mat :: A
  type(matrix_buffer_type), pointer :: matrix_buffer  

  PetscInt :: icol
  PetscErrorCode :: ierr

  do icol = 1, 7
    call MatSetValuesLocal(A,1,0,1,icol-1, &
                           matrix_buffer%values(icol,1),INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

end subroutine MatrixBufferSetValuesHypre

! ************************************************************************** !

subroutine MatrixBufferSetValuesAij(A,matrix_buffer)
  ! 
  ! MatrixBufferSetValues: Sets values in PETSc Aij matrix
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/13/09
  ! 

  implicit none

  Mat :: A
  type(matrix_buffer_type), pointer :: matrix_buffer  

  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr

  do local_id = 1, matrix_buffer%grid%nlmax
    ghosted_id = matrix_buffer%grid%nL2G(local_id)
    call MatSetValuesLocal(A,1,matrix_buffer%icol(4,ghosted_id), &
                           size(matrix_buffer%icol,1), &
                           matrix_buffer%icol(:,ghosted_id), &
                           matrix_buffer%values(:,ghosted_id),INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

end subroutine MatrixBufferSetValuesAij

! ************************************************************************** !

subroutine MatrixBufferDestroy(matrix_buffer)
  ! 
  ! Destroys a matrix buffer object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/13/09
  ! 

  implicit none
  
  type(matrix_buffer_type), pointer :: matrix_buffer
  
  if (.not.associated(matrix_buffer)) return
  
  if (associated(matrix_buffer%values)) deallocate(matrix_buffer%values)
  nullify(matrix_buffer%values)
  if (associated(matrix_buffer%icol)) deallocate(matrix_buffer%icol)
  nullify(matrix_buffer%icol)
  
  deallocate(matrix_buffer)
  nullify(matrix_buffer)
  
end subroutine MatrixBufferDestroy

end module Matrix_Buffer_module
