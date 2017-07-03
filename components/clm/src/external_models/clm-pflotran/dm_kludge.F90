module DM_Kludge_module
  
  use Grid_Unstructured_Aux_module, only : ugdm_type
  
  use PFLOTRAN_Constants_module

  implicit none

  private
 
#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscdm.h"
#include "petsc/finclude/petscdm.h90"
#include "petsc/finclude/petscdmda.h"
#include "petsc/finclude/petscdmshell.h90"

  type, public :: dm_ptr_type
    DM :: dm  ! PETSc DM
    type(ugdm_type), pointer :: ugdm
      ! Unstructured grid "private" dm.  This gets wrapped in a PETSc DM via 
      ! DMShell routines.
  end type dm_ptr_type

end module DM_Kludge_module
