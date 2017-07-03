module Geomechanics_Global_module

  use Geomechanics_Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

  public :: GeomechGlobalSetup, &
            GeomechGlobalSetAuxVarScalar, &
            GeomechGlobalSetAuxVarVecLoc, &
            GeomechGlobalUpdateAuxVars

contains

! ************************************************************************** !

subroutine GeomechGlobalSetup(geomech_realization)
  ! 
  ! Set up global aux vars in a realization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  
  implicit none

  class(realization_geomech_type) :: geomech_realization
  
  ! There is only one patch in each realization
  call GeomechGlobalSetupPatch(geomech_realization)
  
end subroutine GeomechGlobalSetup

! ************************************************************************** !

subroutine GeomechGlobalSetupPatch(geomech_realization)
  ! 
  ! Strips a geomech global auxvar
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Option_module
  use Geomechanics_Coupler_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
 
  implicit none
  
  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(geomech_patch_type),pointer :: patch
  type(geomech_grid_type), pointer :: grid
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id
  type(geomech_global_auxvar_type), pointer :: aux_vars(:)
  PetscInt :: ivertex
  
  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid

  patch%geomech_aux%GeomechGlobal => GeomechGlobalAuxCreate()
  
  allocate(aux_vars(grid%ngmax_node))
  do ghosted_id = 1, grid%ngmax_node
    call GeomechGlobalAuxVarInit(aux_vars(ghosted_id),option)
  enddo
  patch%geomech_aux%GeomechGlobal%aux_vars => aux_vars
  patch%geomech_aux%GeomechGlobal%num_aux = grid%ngmax_node
   
end subroutine GeomechGlobalSetupPatch

! ************************************************************************** !

subroutine GeomechGlobalSetAuxVarScalar(geomech_realization,value,ivar)
  ! 
  ! Strips a geomech global auxvar
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  PetscReal :: value
  PetscInt :: ivar
  
  type(geomech_patch_type), pointer :: cur_patch
  
  cur_patch => geomech_realization%geomech_patch
  call GeomechGlobalSetAuxVarScalarPatch(geomech_realization,value,ivar)
 
end subroutine GeomechGlobalSetAuxVarScalar

! ************************************************************************** !

subroutine GeomechGlobalSetAuxVarScalarPatch(geomech_realization,value,ivar)
  ! 
  ! Strips a geomech global auxvar
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Option_module
  use Geomechanics_Patch_module
  use Variables_module, only : GEOMECH_DISP_X, &
                               GEOMECH_DISP_Y, &
                               GEOMECH_DISP_Z
  
  implicit none

  class(realization_geomech_type) :: geomech_realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
    
  PetscInt :: i
  
  patch => geomech_realization%geomech_patch
  option => geomech_realization%option  
  
  select case(ivar)
    case(GEOMECH_DISP_X)
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux
        patch%geomech_aux%GeomechGlobal%aux_vars(i)%disp_vector(&
          GEOMECH_DISP_X_DOF) = value
      enddo
    case(GEOMECH_DISP_Y)
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux
        patch%geomech_aux%GeomechGlobal%aux_vars(i)%disp_vector(&
          GEOMECH_DISP_Y_DOF) = value
      enddo
    case(GEOMECH_DISP_Z)
      do i=1, patch%geomech_aux%GeomechGlobal%num_aux
        patch%geomech_aux%GeomechGlobal%aux_vars(i)%disp_vector(&
          GEOMECH_DISP_Z_DOF) = value
      enddo
  end select
  
end subroutine GeomechGlobalSetAuxVarScalarPatch

! ************************************************************************** !

subroutine GeomechGlobalSetAuxVarVecLoc(geomech_realization,vec_loc,ivar, &
                                        isubvar)
  ! 
  ! Strips a geomech global auxvar
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module

  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"  

  class(realization_geomech_type) :: geomech_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar
  
  type(geomech_patch_type), pointer :: cur_patch

  cur_patch => geomech_realization%geomech_patch
  call GeomechGlobalSetAuxVarVecLocPatch(geomech_realization,vec_loc,ivar,isubvar)

end subroutine GeomechGlobalSetAuxVarVecLoc

! ************************************************************************** !

subroutine GeomechGlobalSetAuxVarVecLocPatch(geomech_realization,vec_loc,ivar,&
                                             isubvar)
  ! 
  ! Strips a geomech global auxvar
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Option_module
  use Variables_module, only : GEOMECH_DISP_X, &
                               GEOMECH_DISP_Y, &
                               GEOMECH_DISP_Z
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_geomech_type) :: geomech_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
  type(geomech_grid_type), pointer :: grid
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  
  call GeomechGridVecGetArrayF90(grid,vec_loc,vec_loc_p,ierr)
  
  select case(ivar)
    case(GEOMECH_DISP_X)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax_node
            patch%geomech_aux%GeomechGlobal%aux_vars(&
              ghosted_id)%disp_vector(GEOMECH_DISP_X_DOF) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GEOMECH_DISP_Y)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax_node
            patch%geomech_aux%GeomechGlobal%aux_vars(&
              ghosted_id)%disp_vector(GEOMECH_DISP_Y_DOF) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GEOMECH_DISP_Z)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax_node
            patch%geomech_aux%GeomechGlobal%aux_vars(&
              ghosted_id)%disp_vector(GEOMECH_DISP_Z_DOF) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
  end select

  call GeomechGridVecRestoreArrayF90(grid,vec_loc,vec_loc_p,ierr)

end subroutine GeomechGlobalSetAuxVarVecLocPatch

! ************************************************************************** !

subroutine GeomechGlobalUpdateAuxVars(geomech_realization,time_level)
  ! 
  ! Strips a geomech global auxvar
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  ! 

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Option_module
  use Geomechanics_Discretization_module
  use Variables_module, only : GEOMECH_DISP_X, &
                               GEOMECH_DISP_Y, &
                               GEOMECH_DISP_Z
  
  class(realization_geomech_type) :: geomech_realization
  PetscInt :: time_level
  
  type(geomech_field_type), pointer :: geomech_field
  type(option_type), pointer :: option
  
  option => geomech_realization%option
  geomech_field => geomech_realization%geomech_field
  
  ! x displacement
  call GeomechRealizGetDataset(geomech_realization,geomech_field%work, &
                               GEOMECH_DISP_X,ZERO_INTEGER)
  call GeomechDiscretizationGlobalToLocal(&
                              geomech_realization%geomech_discretization, &
                              geomech_field%work,geomech_field%work_loc,ONEDOF)
  call GeomechGlobalSetAuxVarVecLoc(geomech_realization,&
                                    geomech_field%work_loc, &
                                    GEOMECH_DISP_X,time_level)
                                  
  ! y displacement
  call GeomechRealizGetDataset(geomech_realization,geomech_field%work, &
                               GEOMECH_DISP_Y,ZERO_INTEGER)
  call GeomechDiscretizationGlobalToLocal(&
                              geomech_realization%geomech_discretization, &
                              geomech_field%work,geomech_field%work_loc,ONEDOF)
  call GeomechGlobalSetAuxVarVecLoc(geomech_realization, &
                                    geomech_field%work_loc, &
                                    GEOMECH_DISP_Y,time_level)

  ! z displacement
  call GeomechRealizGetDataset(geomech_realization,geomech_field%work, &
                               GEOMECH_DISP_Z,ZERO_INTEGER)
  call GeomechDiscretizationGlobalToLocal(&
                              geomech_realization%geomech_discretization, &
                              geomech_field%work,geomech_field%work_loc,ONEDOF)
  call GeomechGlobalSetAuxVarVecLoc(geomech_realization, &
                                    geomech_field%work_loc, &
                                    GEOMECH_DISP_Z,time_level)


end subroutine GeomechGlobalUpdateAuxVars

end module Geomechanics_Global_module
