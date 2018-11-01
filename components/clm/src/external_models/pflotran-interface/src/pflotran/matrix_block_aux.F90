module Matrix_Block_Aux_module

  ! this module cannot depend on any other modules beside Option_module

#include "petsc/finclude/petscsys.h"
   use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private 


  type, public :: matrix_block_auxvar_type

    PetscReal, pointer :: dtotal(:,:,:)
    
  end type matrix_block_auxvar_type
  
  type, public :: matrix_block_info_type
    PetscInt :: dim1
    PetscInt :: dim2
    PetscInt :: dim3
  end type matrix_block_info_type
  
  interface MatrixBlockAuxVarInit
    module procedure MatrixBlockAuxVarInit1
    module procedure MatrixBlockAuxVarInit2
  end interface MatrixBlockAuxVarInit
  
  public :: MatrixBlockAuxVarCreate, &
            MatrixBlockAuxVarInit, &
            MatrixBlockAuxVarCopy, &
            MatrixBlockAuxVarDestroy, &
            MatrixBlockInfoCreate, &
            MatrixBlockInfoDestroy
            
contains

! ************************************************************************** !

function MatrixBlockAuxVarCreate(option)
  ! 
  ! MatrixBlockAuxCreate: Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/2010
  ! 

  use Option_module

  implicit none
  
  type(option_type) :: option
  type(matrix_block_auxvar_type), pointer :: MatrixBlockAuxVarCreate
  
  type(matrix_block_auxvar_type), pointer :: aux

  allocate(aux)  
  nullify(aux%dtotal)

  MatrixBlockAuxVarCreate => aux
  
end function MatrixBlockAuxVarCreate

! ************************************************************************** !

subroutine MatrixBlockAuxVarInit1(auxvar,dim1,dim2,dim3,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/2010
  ! 

  use Option_module

  implicit none
  
  type(matrix_block_auxvar_type) :: auxvar
  type(matrix_block_info_type) :: matrix_info
  PetscInt :: dim1
  PetscInt :: dim2
  PetscInt :: dim3
  type(option_type) :: option  
  
  allocate(auxvar%dtotal(dim1,dim2,dim3))
  auxvar%dtotal = 0.d0
  
end subroutine MatrixBlockAuxVarInit1

! ************************************************************************** !

subroutine MatrixBlockAuxVarInit2(auxvar,matrix_info,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/2010
  ! 

  use Option_module

  implicit none
  
  type(matrix_block_auxvar_type) :: auxvar
  type(matrix_block_info_type) :: matrix_info
  type(option_type) :: option  
  
  allocate(auxvar%dtotal(matrix_info%dim1,matrix_info%dim2,matrix_info%dim3))
  auxvar%dtotal = 0.d0
  
end subroutine MatrixBlockAuxVarInit2

! ************************************************************************** !

subroutine MatrixBlockAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copys an auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/2010
  ! 

  use Option_module

  implicit none
  
  type(matrix_block_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option  
  
  auxvar%dtotal = auxvar2%dtotal
  
end subroutine MatrixBlockAuxVarCopy

! ************************************************************************** !

subroutine MatrixBlockAuxVarDestroy(auxvar)
  ! 
  ! Deallocates a matrix block auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/2010
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(matrix_block_auxvar_type), pointer :: auxvar
  
  if (.not.associated(auxvar)) return
  
  call DeallocateArray(auxvar%dtotal)
  
  deallocate(auxvar)
  nullify(auxvar)
  
end subroutine MatrixBlockAuxVarDestroy

! ************************************************************************** !

function MatrixBlockInfoCreate(dim1,dim2,dim3,option)
  ! 
  ! Allocate and initialize matrix block info object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/09/2010
  ! 

  use Option_module

  implicit none
  
  PetscInt :: dim1
  PetscInt :: dim2
  PetscInt :: dim3
  type(option_type) :: option
  
  type(matrix_block_info_type), pointer :: MatrixBlockInfoCreate
  
  type(matrix_block_info_type), pointer :: info

  allocate(info) 
  info%dim1 = dim1 
  info%dim2 = dim2 
  info%dim3 = dim3 

  MatrixBlockInfoCreate => info
  
end function MatrixBlockInfoCreate

! ************************************************************************** !

subroutine MatrixBlockInfoDestroy(matrix_info)
  ! 
  ! Deallocates a matrix block info object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/08/2010
  ! 

  implicit none

  type(matrix_block_info_type), pointer :: matrix_info
  
  if (.not.associated(matrix_info)) return
  
  if (associated(matrix_info))deallocate(matrix_info)
  nullify(matrix_info)
  
end subroutine MatrixBlockInfoDestroy

end module Matrix_Block_Aux_module
