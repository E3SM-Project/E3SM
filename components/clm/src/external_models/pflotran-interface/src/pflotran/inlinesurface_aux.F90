module InlineSurface_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none

  private


  type, public :: inlinesurface_auxvar_type

    PetscReal :: half_cell_height
    PetscReal :: surface_water_depth
    PetscReal :: density
    PetscReal :: Mannings_coeff

  end type inlinesurface_auxvar_type

  type, public :: inlinesurface_type

    PetscBool :: auxvars_up_to_date
    PetscInt  :: num_aux, num_aux_bc, num_aux_ss
    type(inlinesurface_auxvar_type), pointer :: auxvars(:)
    type(inlinesurface_auxvar_type), pointer :: auxvars_bc(:)
    type(inlinesurface_auxvar_type), pointer :: auxvars_ss(:)

  end type inlinesurface_type

  public ::                        &
       InlineSurfaceAuxCreate,     &
       InlineSurfaceAuxDestroy,    &
       InlineSurfaceAuxVarCompute, &
       InlineSurfaceAuxVarInit,    &
       InlineSurfaceAuxVarCopy

contains

  !************************************************************************** !

  subroutine InlineSurfaceAuxVarInit(auxvar,option)
    ! 
    ! Initialize auxiliary variable
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    ! 
    use Option_module
    implicit none
    type(inlinesurface_auxvar_type) :: auxvar
    type(option_type)               :: option
    PetscReal      :: NaN

    NaN = 0.d0
    NaN = 1.d0/NaN
    NaN = 0.d0*NaN
    auxvar%half_cell_height    = NaN
    auxvar%surface_water_depth = NaN
    auxvar%Mannings_coeff      = NaN
    auxvar%density             = NaN
  end subroutine InlineSurfaceAuxVarInit

  !************************************************************************** !

  subroutine InlineSurfaceAuxVarCopy(source,copy)
    ! 
    ! Copies an auxiliary variable
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    ! 
    implicit none
    type(inlinesurface_auxvar_type) :: source, copy

    copy%half_cell_height    = source%half_cell_height
    copy%surface_water_depth = source%surface_water_depth
    copy%density             = source%density
    copy%Mannings_coeff      = source%Mannings_coeff

  end subroutine InlineSurfaceAuxVarCopy

  !************************************************************************** !

  subroutine InlineSurfaceAuxVarCompute(auxvar,global_auxvar,option)
    ! 
    ! Evaluates the auxiliary variable, assumes global_auxvar
    ! has been updated previously
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    ! 
    use Option_module
    use Global_Aux_module
    implicit none
    type(inlinesurface_auxvar_type) :: auxvar
    type(global_auxvar_type)        :: global_auxvar
    type(option_type)               :: option
    PetscReal :: Pref,Pl,g,rho

    ! We need to ensure that surface density stays consistent with the subsurface
    auxvar%density = global_auxvar%den(1)

    Pref = option%reference_pressure
    Pl   = global_auxvar%pres(1)
    rho  = global_auxvar%den_kg(1)
    g    = ABS(option%gravity(3))
    auxvar%surface_water_depth = MAX(0.0d0,(Pl-Pref)/(rho*g)-auxvar%half_cell_height)

  end subroutine InlineSurfaceAuxVarCompute

  !************************************************************************** !

  function InlineSurfaceAuxCreate() result(aux)
    ! 
    ! Allocates a inlinesurface auxiliary object
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    !
    implicit none
    type(inlinesurface_type), pointer :: aux

    allocate(aux)
    aux%auxvars_up_to_date = PETSC_FALSE
    aux%num_aux    = 0
    aux%num_aux_bc = 0
    aux%num_aux_ss = 0
    nullify(aux%auxvars)
    nullify(aux%auxvars_bc)
    nullify(aux%auxvars_ss)

  end function InlineSurfaceAuxCreate

  !************************************************************************** !

  subroutine InlineSurfaceAuxDestroy(aux)
    ! 
    ! Deallocates a inlinesurface auxiliary object
    ! 
    ! Author: Nathan Collier
    ! Date: 09/2015
    !
    implicit none
    type(inlinesurface_type), pointer :: aux

    if (.not.associated(aux)) return
    if (associated(aux%auxvars)) then
      deallocate(aux%auxvars)
    endif
    nullify(aux%auxvars)
    if (associated(aux%auxvars_bc)) then
      deallocate(aux%auxvars_bc)
    endif
    nullify(aux%auxvars_bc)
    if (associated(aux%auxvars_ss)) then
      deallocate(aux%auxvars_ss)
    endif
    nullify(aux%auxvars_ss)
    deallocate(aux)
    nullify(aux)

  end subroutine InlineSurfaceAuxDestroy

end module InlineSurface_Aux_module
