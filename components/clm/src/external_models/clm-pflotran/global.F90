module Global_module

  use Global_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"
  
  public GlobalSetup, &
         GlobalSetAuxVarScalar, &
         GlobalSetAuxVarVecLoc, &
         GlobalGetAuxVarVecLoc, &
         GlobalWeightAuxVars, &
         GlobalUpdateState, &
         GlobalUpdateAuxVars

contains

! ************************************************************************** !

subroutine GlobalSetup(realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
 
  implicit none
  
  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink

  PetscInt :: ghosted_id, iconn, sum_connection
  type(global_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: auxvars_bc(:)
  type(global_auxvar_type), pointer :: auxvars_ss(:)
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Global => GlobalAuxCreate()
  
  ! allocate auxvar data structures for all grid cells  
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else  
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call GlobalAuxVarInit(auxvars(ghosted_id),option)
  enddo
  patch%aux%Global%auxvars => auxvars
  patch%aux%Global%num_aux = grid%ngmax
  
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
      call GlobalAuxVarInit(auxvars_bc(iconn),option)
    enddo
    patch%aux%Global%auxvars_bc => auxvars_bc
  endif
  patch%aux%Global%num_aux_bc = sum_connection

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
      call GlobalAuxVarInit(auxvars_ss(iconn),option)
    enddo
    patch%aux%Global%auxvars_ss => auxvars_ss
  endif
  patch%aux%Global%num_aux_ss = sum_connection

  option%iflag = 0
  
end subroutine GlobalSetup

! ************************************************************************** !

subroutine GlobalSetAuxVarScalar(realization,value,ivar)
  ! 
  ! Sets values of auxvar data using a scalar value.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/19/08
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               LIQUID_DENSITY, GAS_PRESSURE, &
                               GAS_DENSITY, GAS_SATURATION, &
                               TEMPERATURE, LIQUID_DENSITY_MOL
  
  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: value
  PetscInt :: ivar

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
    
  PetscInt :: i
  
  patch => realization%patch
  option => realization%option  
  
  select case(ivar)
    case(LIQUID_PRESSURE)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%auxvars(i)%pres = value
      enddo
      do i=1, patch%aux%Global%num_aux_bc
        patch%aux%Global%auxvars_bc(i)%pres = value
      enddo
    case(TEMPERATURE)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%auxvars(i)%temp = value
      enddo
      do i=1, patch%aux%Global%num_aux_bc
        patch%aux%Global%auxvars_bc(i)%temp = value
      enddo
    case(LIQUID_DENSITY)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%auxvars(i)%den_kg(option%liquid_phase) = value
      enddo
      do i=1, realization%patch%aux%Global%num_aux_bc
        patch%aux%Global%auxvars_bc(i)%den_kg(option%liquid_phase) = value
      enddo
    case(LIQUID_DENSITY_MOL)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%auxvars(i)%den(option%liquid_phase) = value/FMWH2O
      enddo
      do i=1, realization%patch%aux%Global%num_aux_bc
        patch%aux%Global%auxvars_bc(i)%den(option%liquid_phase) = value/FMWH2O
      enddo
    case(LIQUID_SATURATION)
      do i=1, patch%aux%Global%num_aux
        patch%aux%Global%auxvars(i)%sat(option%liquid_phase) = value
      enddo
      do i=1, patch%aux%Global%num_aux_bc
        patch%aux%Global%auxvars_bc(i)%sat(option%liquid_phase) = value
      enddo
  end select
  
end subroutine GlobalSetAuxVarScalar

! ************************************************************************** !

subroutine GlobalSetAuxVarVecLoc(realization,vec_loc,ivar,isubvar)
  ! 
  ! Sets values of auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/19/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               LIQUID_DENSITY, GAS_PRESSURE, &
                               GAS_DENSITY, GAS_SATURATION, &
                               TEMPERATURE, SC_FUGA_COEFF, GAS_DENSITY_MOL, &
                               STATE, OIL_PRESSURE, OIL_SATURATION, &
                               OIL_DENSITY, OIL_DENSITY_MOL
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_subsurface_type) :: realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(LIQUID_PRESSURE)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres_store(option%liquid_phase,TIME_T) = &
                vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres_store(option%liquid_phase,TIME_TpDT) = &
                vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres(option%liquid_phase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GAS_PRESSURE)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres_store(option%gas_phase,TIME_T) = &
                vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres_store(option%gas_phase,TIME_TpDT) = &
                vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres(option%gas_phase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(OIL_PRESSURE)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres_store(option%oil_phase,TIME_T) = &
                vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres_store(option%oil_phase,TIME_TpDT) = &
                vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%pres(option%oil_phase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(TEMPERATURE)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%temp_store(TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%temp_store(TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%temp = vec_loc_p(ghosted_id)
          enddo
      end select
    case(LIQUID_DENSITY)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg_store(option%liquid_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg_store(option%liquid_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg(option%liquid_phase) = vec_loc_p(ghosted_id)
            patch%aux%Global%auxvars(ghosted_id)%den(option%liquid_phase) = &
              vec_loc_p(ghosted_id)/FMWH2O
          enddo
      end select
    case(GAS_SATURATION)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat_store(option%gas_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat_store(option%gas_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat(option%gas_phase) = &
              vec_loc_p(ghosted_id)
          enddo
      end select
    case(OIL_SATURATION)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat_store(option%oil_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat_store(option%oil_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat(option%oil_phase) = &
              vec_loc_p(ghosted_id)
          enddo
      end select

    case(GAS_DENSITY)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg_store(option%gas_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg_store(option%gas_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg(option%gas_phase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(GAS_DENSITY_MOL)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_store(option%gas_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_store(option%gas_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den(option%gas_phase) = vec_loc_p(ghosted_id)
          enddo
      end select

    case(OIL_DENSITY)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg_store(option%oil_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg_store(option%oil_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_kg(option%oil_phase) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(OIL_DENSITY_MOL)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_store(option%oil_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den_store(option%oil_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%den(option%oil_phase) = vec_loc_p(ghosted_id)
          enddo
      end select

    case(LIQUID_SATURATION)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat_store(option%liquid_phase,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat_store(option%liquid_phase,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%sat(option%liquid_phase) = &
              vec_loc_p(ghosted_id)
          enddo
      end select
    case(SC_FUGA_COEFF)
      select case(isubvar)
        case(TIME_T)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%fugacoeff_store(1,TIME_T) = &
              vec_loc_p(ghosted_id)
          enddo
        case(TIME_TpDT)
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%fugacoeff_store(1,TIME_TpDT) = &
              vec_loc_p(ghosted_id)
          enddo
        case default
          do ghosted_id=1, grid%ngmax
            patch%aux%Global%auxvars(ghosted_id)%fugacoeff(1) = vec_loc_p(ghosted_id)
          enddo
      end select
    case(STATE)
      do ghosted_id=1, grid%ngmax
        patch%aux%Global%auxvars(ghosted_id)%istate = &
          int(vec_loc_p(ghosted_id)+1.d-10)
      enddo
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine GlobalSetAuxVarVecLoc

! ************************************************************************** !

subroutine GlobalGetAuxVarVecLoc(realization,vec_loc,ivar,isubvar)
  ! 
  ! Sets values of auxvar data using a vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/19/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Variables_module, only : STATE
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(realization_subsurface_type) :: realization
  Vec :: vec_loc
  PetscInt :: ivar
  PetscInt :: isubvar  
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  
  PetscInt :: ghosted_id
  PetscReal, pointer :: vec_loc_p(:)
  PetscErrorCode :: ierr
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  
  call VecGetArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(STATE)
      do ghosted_id=1, grid%ngmax
        vec_loc_p(ghosted_id) = &
          dble(patch%aux%Global%auxvars(ghosted_id)%istate)
      enddo
    case default
      option%io_buffer = 'Variable unrecognized in GlobalGetAuxVarVecLoc.'
      call printErrMsg(option)
  end select

  call VecRestoreArrayReadF90(vec_loc,vec_loc_p,ierr);CHKERRQ(ierr)

end subroutine GlobalGetAuxVarVecLoc

! ************************************************************************** !

subroutine GlobalWeightAuxVars(realization,weight)
  ! 
  ! Updates the densities and saturations in auxiliary
  ! variables associated with reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/03/08
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Material_module, only : MaterialWeightAuxVars
  
  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: weight
  
  type(option_type), pointer :: option
  type(global_auxvar_type), pointer :: auxvars(:)
  PetscInt :: ghosted_id
  
  option => realization%option
  auxvars => realization%patch%aux%Global%auxvars
  
  do ghosted_id = 1, realization%patch%aux%Global%num_aux
    ! interpolate density and saturation based on weight
    auxvars(ghosted_id)%den_kg(:) = &
      (weight*auxvars(ghosted_id)%den_kg_store(:,TIME_TpDT)+ &
       (1.d0-weight)*auxvars(ghosted_id)%den_kg_store(:,TIME_T))
    auxvars(ghosted_id)%sat(:) = &
      (weight*auxvars(ghosted_id)%sat_store(:,TIME_TpDT)+ &
       (1.d0-weight)*auxvars(ghosted_id)%sat_store(:,TIME_T))
  enddo
  
  select case(option%iflowmode) 
    case(G_MODE)
      do ghosted_id = 1, realization%patch%aux%Global%num_aux
        auxvars(ghosted_id)%pres(:) = &
          (weight*auxvars(ghosted_id)%pres_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%pres_store(:,TIME_T))
        auxvars(ghosted_id)%temp = &
          (weight*auxvars(ghosted_id)%temp_store(TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%temp_store(TIME_T))
      enddo  
    case(MPH_MODE,FLASH2_MODE)
      ! need future implementation for ims_mode too    
      do ghosted_id = 1, realization%patch%aux%Global%num_aux
        auxvars(ghosted_id)%pres(:) = &
          (weight*auxvars(ghosted_id)%pres_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%pres_store(:,TIME_T))
        auxvars(ghosted_id)%temp = &
          (weight*auxvars(ghosted_id)%temp_store(TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%temp_store(TIME_T))
        auxvars(ghosted_id)%fugacoeff(:) = &
          (weight*auxvars(ghosted_id)%fugacoeff_store(:,TIME_TpDT)+ &
           (1.d0-weight)*auxvars(ghosted_id)%fugacoeff_store(:,TIME_T))
        if (weight<1D-12) auxvars(ghosted_id)%reaction_rate(:)=0D0
  !      auxvars(ghosted_id)%den(:) = &
  !        (weight*auxvars(ghosted_id)%den_store(:,TIME_TpDT)+ &
  !         (1.d0-weight)*auxvars(ghosted_id)%den_store(:,TIME_T))
      enddo  
  end select
  
end subroutine GlobalWeightAuxVars

! ************************************************************************** !

subroutine GlobalUpdateState(realization)
  ! 
  ! Updates global aux var variables for use in
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/09
  ! 

  use Realization_Subsurface_class
  use Realization_Base_class, only : RealizationGetVariable
  use Communicator_Base_module
  use Variables_module, only : STATE
  
  class(realization_subsurface_type) :: realization
  
  call RealizationGetVariable(realization,realization%field%work,STATE, &
                              ZERO_INTEGER)
  call realization%comm1%GlobalToLocal(realization%field%work, &
                                       realization%field%work_loc)
  call GlobalSetAuxVarVecLoc(realization,realization%field%work_loc,STATE, &
                             ZERO_INTEGER)
  
end subroutine GlobalUpdateState

! ************************************************************************** !

subroutine GlobalUpdateAuxVars(realization,time_level,time)
  ! 
  ! Updates global aux var variables for use in
  ! reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/09
  ! 

  use Realization_Subsurface_class
  use Realization_Base_class, only : RealizationGetVariable
  use Field_module
  use Option_module
  use Communicator_Base_module
  use Material_module, only : MaterialUpdateAuxVars
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION, &
                               LIQUID_DENSITY, GAS_PRESSURE, &
                               GAS_DENSITY, GAS_SATURATION, &
                               TEMPERATURE, SC_FUGA_COEFF, GAS_DENSITY_MOL
  
  class(realization_subsurface_type) :: realization
  PetscReal :: time
  PetscInt :: time_level
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  
  option => realization%option
  field => realization%field
  
  select case(time_level)
    case(TIME_T)
      realization%patch%aux%Global%time_t = time
    case(TIME_TpDT)
      realization%patch%aux%Global%time_tpdt = time
  end select  
  
  ! liquid density
  call RealizationGetVariable(realization,field%work,LIQUID_DENSITY, &
                             ZERO_INTEGER)
  call realization%comm1%GlobalToLocal(field%work,field%work_loc)
  call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_DENSITY, &
                             time_level)

  ! liquid saturation
  call RealizationGetVariable(realization,field%work,LIQUID_SATURATION, &
                             ZERO_INTEGER)
  call realization%comm1%GlobalToLocal(field%work,field%work_loc)
  call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_SATURATION, &
                             time_level)
  select case(option%iflowmode)
    case(MPH_MODE,FLASH2_MODE)
      ! Gas density
      call RealizationGetVariable(realization,field%work,GAS_DENSITY, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY, &
                                 time_level)
      call RealizationGetVariable(realization,field%work,GAS_DENSITY_MOL, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY_MOL, &
                                 time_level)
 
 
      ! Gas saturation
      call RealizationGetVariable(realization,field%work,GAS_SATURATION, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_SATURATION, &
                                 time_level)                         
   
      ! liquid pressure
      call RealizationGetVariable(realization,field%work,LIQUID_PRESSURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_PRESSURE, &
                                 time_level)                      
 
      ! gas pressure
      call RealizationGetVariable(realization,field%work,GAS_PRESSURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_PRESSURE, &
                                 time_level)
 
      ! temperature
      call RealizationGetVariable(realization,field%work,TEMPERATURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,TEMPERATURE, &
                                 time_level)
      
      ! fugacity coeff
      call RealizationGetVariable(realization,field%work,SC_FUGA_COEFF, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,SC_FUGA_COEFF, &
                                 time_level)
    case(TH_MODE)
      ! pressure
      call RealizationGetVariable(realization,field%work,LIQUID_PRESSURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_PRESSURE, &
                                 time_level)
 
      ! temperature
      call RealizationGetVariable(realization,field%work,TEMPERATURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,TEMPERATURE, &
                                 time_level)
    case(G_MODE)
      ! pressure
      call RealizationGetVariable(realization,field%work,LIQUID_PRESSURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_PRESSURE, &
                                 time_level)
      ! temperature
      call RealizationGetVariable(realization,field%work,TEMPERATURE, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,TEMPERATURE, &
                                 time_level)
      ! Gas density
      call RealizationGetVariable(realization,field%work,GAS_DENSITY, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY, &
                                 time_level)
      ! Gas saturation
      call RealizationGetVariable(realization,field%work,GAS_SATURATION, &
                                  ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_SATURATION, &
                                 time_level)                         
    case(IMS_MODE)
      ! Gas density
      call RealizationGetVariable(realization,field%work,GAS_DENSITY, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY, &
                                 time_level)
      call RealizationGetVariable(realization,field%work,GAS_DENSITY_MOL, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY_MOL, &
                                 time_level)
 
 
      ! Gas saturation
      call RealizationGetVariable(realization,field%work,GAS_SATURATION, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_SATURATION, &
                                 time_level)
   
      ! pressure
      call RealizationGetVariable(realization,field%work,LIQUID_PRESSURE, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_PRESSURE, &
                                 time_level)
 
      ! temperature
      call RealizationGetVariable(realization,field%work,TEMPERATURE, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,TEMPERATURE, &
                                 time_level)

    case(MIS_MODE)
      ! Gas density
      call RealizationGetVariable(realization,field%work,GAS_DENSITY, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY, &
                                 time_level)
      call RealizationGetVariable(realization,field%work,GAS_DENSITY_MOL, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_DENSITY_MOL, &
                                 time_level)
 
 
      ! Gas saturation
      call RealizationGetVariable(realization,field%work,GAS_SATURATION, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,GAS_SATURATION, &
                                 time_level)                  
   
      ! pressure
      call RealizationGetVariable(realization,field%work,LIQUID_PRESSURE, &
                                 ZERO_INTEGER)
      call realization%comm1%GlobalToLocal(field%work,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,LIQUID_PRESSURE, &
                                 time_level)
 
      ! temperature
      call RealizationGetVariable(realization,field%work,TEMPERATURE, &
                                 ZERO_INTEGER)
      call realization%comm1%LocalToLocal(field%work_loc,field%work_loc)
      call GlobalSetAuxVarVecLoc(realization,field%work_loc,TEMPERATURE, &
                                 time_level)
  end select

end subroutine GlobalUpdateAuxVars

end module Global_module
