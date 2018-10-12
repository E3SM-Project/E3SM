module DM_Kludge_module
#include "petsc/finclude/petscdm.h"
  use petscdm
  use Grid_Unstructured_Aux_module, only : ugdm_type
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: dm_ptr_type
    DM :: dm  ! PETSc DM
    type(ugdm_type), pointer :: ugdm
      ! Unstructured grid "private" dm.  This gets wrapped in a PETSc DM via 
      ! DMShell routines.
  end type dm_ptr_type

end module DM_Kludge_module
