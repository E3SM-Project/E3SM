module Surface_Global_module

  use Surface_Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

  public SurfaceGlobalSetup, &
         SurfaceGlobalSetAuxVarScalar, &
         SurfaceGlobalSetAuxVarVecLoc, &
         SurfaceGlobalUpdateAuxVars

contains

! ************************************************************************** !

subroutine SurfaceGlobalSetup(surf_realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Realization_Surface_class
  use Patch_module
  
  implicit none

  class(realization_surface_type) :: surf_realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    surf_realization%patch => cur_patch
    call SurfaceGlobalSetupPatch(surf_realization)
    cur_patch => cur_patch%next
  enddo

end subroutine SurfaceGlobalSetup

! ************************************************************************** !

subroutine SurfaceGlobalSetupPatch(surf_realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Realization_Surface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
 
  implicit none
  
  class(realization_surface_type) :: surf_realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id, iconn, sum_connection
  type(surface_global_auxvar_type), pointer :: auxvars(:)
  type(surface_global_auxvar_type), pointer :: auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: auxvars_ss(:)
  
  option => surf_realization%option
  patch => surf_realization%patch
  grid => patch%grid

  patch%surf_aux%SurfaceGlobal => SurfaceGlobalAuxCreate()
  
  ! allocate auxvar data structures for all grid cells  
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else  
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call SurfaceGlobalAuxVarInit(auxvars(ghosted_id),option)
  enddo
  patch%surf_aux%SurfaceGlobal%auxvars => auxvars
  patch%surf_aux%SurfaceGlobal%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! auxvar data structures for them  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceGlobalAuxVarInit(auxvars_bc(iconn),option)
    enddo
    patch%surf_aux%SurfaceGlobal%auxvars_bc => auxvars_bc
  endif
  patch%surf_aux%SurfaceGlobal%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  source_sink => patch%source_sink_list%first
  sum_connection = 0    
  do 
    if (.not.associated(source_sink)) exit
    sum_connection = sum_connection + &
                     source_sink%connection_set%num_connections
    source_sink => source_sink%next
  enddo

  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call SurfaceGlobalAuxVarInit(auxvars_ss(iconn),option)
    enddo
    patch%surf_aux%SurfaceGlobal%auxvars_ss => auxvars_ss
  endif
  patch%surf_aux%SurfaceGlobal%num_aux_ss = sum_connection

  option%iflag = 0
  
end subroutine SurfaceGlobalSetupPatch

! ************************************************************************** !

subroutine SurfaceGlobalSetAuxVarScalar(surf_realization,value,ivar)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Realization_Surface_class
  use Patch_module

  implicit none

  class(realization_surface_type) :: surf_realization
  PetscReal :: value
  PetscInt :: ivar
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    surf_realization%patch => cur_patch
    call SurfaceGlobalSetAuxVarScalarPatch(surf_realization,value,ivar)
    cur_patch => cur_patch%next
  enddo
  
end subroutine SurfaceGlobalSetAuxVarScalar

! ************************************************************************** !

subroutine SurfaceGlobalSetAuxVarScalarPatch(surf_realization,value,ivar)

  use Realization_Surface_class
  use Option_module
  use Patch_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               SURFACE_LIQUID_DENSITY
  
  implicit none

  class(realization_surface_type) :: surf_realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
    
  PetscInt :: i
  
  patch => surf_realization%patch
  option => surf_realization%option  
  
  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux
        patch%surf_aux%SurfaceGlobal%auxvars(i)%head = value
      enddo
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux_bc
        patch%surf_aux%SurfaceGlobal%auxvars_bc(i)%head = value
      enddo
    case(SURFACE_LIQUID_TEMPERATURE)
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux
        patch%surf_aux%SurfaceGlobal%auxvars(i)%temp = value
      enddo
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux_bc
        patch%surf_aux%SurfaceGlobal%auxvars_bc(i)%temp = value
      enddo
    case(SURFACE_LIQUID_DENSITY)
      do i=1, patch%surf_aux%SurfaceGlobal%num_aux
        patch%surf_aux%SurfaceGlobal%auxvars(i)%den_kg(option%liquid_phase) = value
      enddo
      do i=1, surf_realization%patch%surf_aux%SurfaceGlobal%num_aux_bc
        patch%surf_aux%SurfaceGlobal%auxvars_bc(i)%den_kg(option%liquid_phase) = value
      enddo
  end select
  
end subroutine SurfaceGlobalSetAuxVarScalarPatch

! ************************************************************************** !

subroutine SurfaceGlobalSetAuxVarVecLoc(surf_realization,vec_loc,ivar,isubvar)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Surface_class
  use Patch_module

  implicit none

  class(realization_surface_type) :: surf_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => surf_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    surf_realization%patch => cur_patch
    call SurfaceGlobalSetAuxVarVecLocPatch(surf_realization,vec_loc,ivar,isubvar)
    cur_patch => cur_patch%next
  enddo

end subroutine SurfaceGlobalSetAuxVarVecLoc

! ************************************************************************** !

subroutine SurfaceGlobalSetAuxVarVecLocPatch(surf_realization,vec_loc,ivar,isubvar)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Surface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               SURFACE_LIQUID_DENSITY
  
  implicit none

  class(realization_surface_type) :: surf_realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => surf_realization%patch
  grid => patch%grid
  option => surf_realization%option
  
  call VecGetArrayF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(SURFACE_LIQUID_HEAD)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%surf_aux%SurfaceGlobal%auxvars(ghosted_id)%head(option%liquid_phase) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
    case(SURFACE_LIQUID_TEMPERATURE)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%surf_aux%SurfaceGlobal%auxvars(ghosted_id)%temp = vec_loc_p(ghosted_id)
          enddo
      end select
    case(SURFACE_LIQUID_DENSITY)
      select case(isubvar)
        case default
          do ghosted_id=1, grid%ngmax
            patch%surf_aux%SurfaceGlobal%auxvars(ghosted_id)%den_kg(option%liquid_phase) &
              = vec_loc_p(ghosted_id)
          enddo
      end select
  end select

  call VecRestoreArrayF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine SurfaceGlobalSetAuxVarVecLocPatch

! ************************************************************************** !

subroutine SurfaceGlobalUpdateAuxVars(surf_realization,time_level)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/07/13
  ! 

  use Realization_Surface_class
  use Surface_Field_module
  use Option_module
  use Discretization_module
  use Variables_module, only : SURFACE_LIQUID_HEAD, &
                               SURFACE_LIQUID_TEMPERATURE, &
                               LIQUID_DENSITY
  
  class(realization_surface_type) :: surf_realization
  PetscInt :: time_level
  
  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  
  option => surf_realization%option
  surf_field => surf_realization%surf_field
  
  ! liquid density
  call RealizSurfGetVariable(surf_realization,surf_field%work,LIQUID_DENSITY, &
                             ZERO_INTEGER)
  call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%work,surf_field%work_loc,ONEDOF)
  call SurfaceGlobalSetAuxVarVecLoc(surf_realization,surf_field%work_loc, &
                                    LIQUID_DENSITY,time_level)

  select case(option%iflowmode)
    case(TH_MODE)
      ! head
      call RealizSurfGetVariable(surf_realization,surf_field%work, &
              SURFACE_LIQUID_HEAD,ZERO_INTEGER)
      call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                  surf_field%work,surf_field%work_loc,ONEDOF)
      call SurfaceGlobalSetAuxVarVecLoc(surf_realization,surf_field%work_loc, &
              SURFACE_LIQUID_HEAD,time_level)
 
      ! temperature
      call RealizSurfGetVariable(surf_realization,surf_field%work, &
              SURFACE_LIQUID_TEMPERATURE, ZERO_INTEGER)
      call DiscretizationGlobalToLocal(surf_realization%discretization, &
                                   surf_field%work,surf_field%work_loc,ONEDOF)
      call SurfaceGlobalSetAuxVarVecLoc(surf_realization,surf_field%work_loc, &
              SURFACE_LIQUID_TEMPERATURE,time_level)
      

  end select

end subroutine SurfaceGlobalUpdateAuxVars

end module Surface_Global_module
