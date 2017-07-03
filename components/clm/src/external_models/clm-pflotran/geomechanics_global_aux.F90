module Geomechanics_Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: geomech_global_auxvar_type
    PetscReal, pointer :: disp_vector(:)       ! [m]
    PetscReal, pointer :: rel_disp_vector(:)   ! [m]
    PetscReal, pointer :: strain(:)            ! dimensionless -- xx, yy, zz, xy, yz, zx
    PetscReal, pointer :: stress(:)            ! [Pa]
    PetscInt :: count                      ! Number of elements shared by a vertex
    ! The count above will be used for averaging the strains and stresses
    ! over the elements
  end type geomech_global_auxvar_type
  
  type, public :: geomech_global_type
    PetscInt :: num_aux
    type(geomech_global_auxvar_type), pointer :: aux_vars(:)
  end type geomech_global_type
  
  interface GeomechGlobalAuxVarDestroy
    module procedure GeomechGlobalAuxVarSingleDestroy
    module procedure GeomechGlobalAuxVarArrayDestroy
  end interface GeomechGlobalAuxVarDestroy

  public :: GeomechGlobalAuxCreate, &
            GeomechGlobalAuxDestroy, &
            GeomechGlobalAuxVarInit, &
            GeomechGlobalAuxVarCopy, &
            GeomechGlobalAuxVarDestroy, &
            GeomechGlobalAuxVarStrip

contains

! ************************************************************************** !

function GeomechGlobalAuxCreate()
  ! 
  ! Creates a geomech global aux
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  implicit none
  
  type(geomech_global_type), pointer :: GeomechGlobalAuxCreate
  
  type(geomech_global_type), pointer :: aux

  allocate(aux) 
  aux%num_aux = 0
  nullify(aux%aux_vars)

  GeomechGlobalAuxCreate => aux
  
end function GeomechGlobalAuxCreate

! ************************************************************************** !

subroutine GeomechGlobalAuxVarInit(aux_var,option)
  ! 
  ! Initializes a geomech global aux
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  use Option_module

  implicit none
  
  type(geomech_global_auxvar_type) :: aux_var
  type(option_type) :: option
  
  allocate(aux_var%disp_vector(option%ngeomechdof))
  allocate(aux_var%rel_disp_vector(option%ngeomechdof))
  allocate(aux_var%strain(SIX_INTEGER))
  allocate(aux_var%stress(SIX_INTEGER))
  aux_var%disp_vector = 0.d0
  aux_var%rel_disp_vector = 0.d0
  aux_var%strain = 0.d0
  aux_var%stress = 0.d0
  
end subroutine GeomechGlobalAuxVarInit

! ************************************************************************** !

subroutine GeomechGlobalAuxVarCopy(aux_var,aux_var2,option)
  ! 
  ! Copies a geomech global aux to another
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  use Option_module

  implicit none
  
  type(geomech_global_auxvar_type) :: aux_var, aux_var2
  type(option_type) :: option

  aux_var%disp_vector = aux_var2%disp_vector
  aux_var%rel_disp_vector = aux_var2%rel_disp_vector
  aux_var%strain = aux_var2%strain
  aux_var%stress = aux_var2%stress
  
end subroutine GeomechGlobalAuxVarCopy

! ************************************************************************** !

subroutine GeomechGlobalAuxVarSingleDestroy(aux_var)
  ! 
  ! Destroys a geomech global aux
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  implicit none

  type(geomech_global_auxvar_type), pointer :: aux_var
  
  if (associated(aux_var)) then
    call GeomechGlobalAuxVarStrip(aux_var)
    deallocate(aux_var)
  endif
  nullify(aux_var)

end subroutine GeomechGlobalAuxVarSingleDestroy

! ************************************************************************** !

subroutine GeomechGlobalAuxVarArrayDestroy(aux_vars)
  ! 
  ! Destroys an array of geomech global auxvar
  ! type
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  implicit none

  type(geomech_global_auxvar_type), pointer :: aux_vars(:)
  
  PetscInt :: iaux
  
  if (associated(aux_vars)) then
    do iaux = 1, size(aux_vars)
      call GeomechGlobalAuxVarStrip(aux_vars(iaux))
    enddo  
    deallocate(aux_vars)
  endif
  nullify(aux_vars)

end subroutine GeomechGlobalAuxVarArrayDestroy

! ************************************************************************** !

subroutine GeomechGlobalAuxVarStrip(aux_var)
  ! 
  ! Strips a geomech global auxvar
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(geomech_global_auxvar_type) :: aux_var
  
  call DeallocateArray(aux_var%disp_vector)
  call DeallocateArray(aux_var%rel_disp_vector)
  call DeallocateArray(aux_var%strain)
  call DeallocateArray(aux_var%stress)

end subroutine GeomechGlobalAuxVarStrip

! ************************************************************************** !

subroutine GeomechGlobalAuxDestroy(aux)
  ! 
  ! Destroys a geomech global type
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13
  ! 

  implicit none

  type(geomech_global_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call GeomechGlobalAuxVarDestroy(aux%aux_vars)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GeomechGlobalAuxDestroy

end module Geomechanics_Global_Aux_module
