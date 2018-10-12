module Surface_Auxiliary_module

  use Surface_Global_Aux_module
!  use Surface_Flow_Aux_module
  use Surface_TH_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: surface_auxiliary_type
    type(surface_global_type), pointer :: SurfaceGlobal
    type(surface_th_type), pointer :: SurfaceTH
  end type surface_auxiliary_type
  
  public :: SurfaceAuxInit, &
            SurfaceAuxDestroy

contains

! ************************************************************************** !

subroutine SurfaceAuxInit(surf_aux)
  ! 
  ! This routine initializes a surface-auxiliary object
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  implicit none
  
  type(surface_auxiliary_type) :: surf_aux
  
  nullify(surf_aux%SurfaceGlobal)
  nullify(surf_aux%SurfaceTH)
  
end subroutine SurfaceAuxInit

! ************************************************************************** !

subroutine SurfaceAuxDestroy(surf_aux)
  ! 
  ! This routine deallocates pointers in a surface-auxiliary object
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  implicit none
  
  type(surface_auxiliary_type) :: surf_aux
  
  call SurfaceGlobalAuxDestroy(surf_aux%SurfaceGlobal)
  call SurfaceTHAuxDestroy(surf_aux%SurfaceTH)

  nullify(surf_aux%SurfaceGlobal)
  nullify(surf_aux%SurfaceTH)

end subroutine SurfaceAuxDestroy

end module Surface_Auxiliary_module
