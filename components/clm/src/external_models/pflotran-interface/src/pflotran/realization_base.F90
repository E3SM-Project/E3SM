module Realization_Base_class

  use Patch_module

  use Discretization_module
  use Option_module
  use Input_Aux_module
  use Debug_module
  use Output_Aux_module
  use Field_module
  use Reaction_Aux_module
  use Data_Mediator_Base_class
  use Communicator_Base_module
  use Waypoint_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
  type, public :: realization_base_type

    PetscInt :: id
    type(discretization_type), pointer :: discretization
    class(communicator_type), pointer :: comm1
    type(patch_list_type), pointer :: patch_list
    type(patch_type), pointer :: patch

    type(option_type), pointer :: option
    type(field_type), pointer :: field
    type(debug_type), pointer :: debug
    type(output_option_type), pointer :: output_option
    class(data_mediator_base_type), pointer :: flow_data_mediator_list
    class(data_mediator_base_type), pointer :: tran_data_mediator_list
    
    type(reaction_type), pointer :: reaction
    
  end type realization_base_type
  
  public :: RealizationBaseInit, &
            RealizationGetVariable, &
            RealizGetVariableValueAtCell, &
            RealizationSetVariable, &
            RealizCreateTranMassTransferVec, &
            RealizCreateFlowMassTransferVec, &
            RealizationBaseStrip

contains

! ************************************************************************** !

subroutine RealizationBaseInit(realization_base,option)
  ! 
  ! Initializes variables/objects in base realization class
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/13
  ! 

  implicit none
  
  class(realization_base_type) :: realization_base
  type(option_type), pointer :: option
  
  realization_base%id = 0
  if (associated(option)) then
    realization_base%option => option
  else
    realization_base%option => OptionCreate()
  endif
  realization_base%discretization => DiscretizationCreate()
  nullify(realization_base%comm1)  
  realization_base%field => FieldCreate()
  realization_base%debug => DebugCreate()
  nullify(realization_base%output_option)

  realization_base%patch_list => PatchCreateList()

  nullify(realization_base%reaction)

  nullify(realization_base%patch)
  nullify(realization_base%flow_data_mediator_list)
  nullify(realization_base%tran_data_mediator_list)

end subroutine RealizationBaseInit

! ************************************************************************** !

subroutine RealizationGetVariable(realization_base,vec,ivar,isubvar, &
                                  isubsubvar)
  ! 
  ! Extracts variables indexed by ivar and isubvar from a
  ! realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module

  implicit none

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubsubvar
  
  PetscInt :: isubsubvar_temp
  
  isubsubvar_temp = 0
  if (present(isubsubvar)) isubsubvar_temp = isubsubvar
  
  call PatchGetVariable(realization_base%patch,realization_base%field, &
                       realization_base%reaction,realization_base%option, &
                       realization_base%output_option,vec,ivar,isubvar, &
                       isubsubvar_temp)

end subroutine RealizationGetVariable

! ************************************************************************** !

function RealizGetVariableValueAtCell(realization_base,ghosted_id, &
                                      ivar,isubvar,isubsubvar)
  ! 
  ! Extracts variables indexed by ivar and isubvar
  ! from a realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  ! 

  use Option_module

  implicit none
  
  PetscReal :: RealizGetVariableValueAtCell
  class(realization_base_type) :: realization_base
  PetscInt :: ghosted_id
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt, optional :: isubsubvar
  
  PetscReal :: value
  PetscInt :: isubsubvar_temp
  
  isubsubvar_temp = 0
  if (present(isubsubvar)) isubsubvar_temp = isubsubvar
  
  value = PatchGetVariableValueAtCell(realization_base%patch, &
                                      realization_base%field, &
                                      realization_base%reaction, &
                                      realization_base%option, &
                                      realization_base%output_option, &
                                      ghosted_id,ivar,isubvar,isubsubvar_temp)
  RealizGetVariableValueAtCell = value

end function RealizGetVariableValueAtCell

! ************************************************************************** !

subroutine RealizationSetVariable(realization_base,vec,vec_format,ivar,isubvar)
  ! 
  ! Sets variables indexed by ivar and isubvar in a
  ! realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module

  implicit none

  class(realization_base_type) :: realization_base
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar

  call PatchSetVariable(realization_base%patch,realization_base%field, &
                       realization_base%option, &
                       vec,vec_format,ivar,isubvar)

end subroutine RealizationSetVariable

! ************************************************************************** !

subroutine RealizCreateFlowMassTransferVec(this)
  ! 
  ! Creates the Vec where mass transfer is summed prior to being added to
  ! the reactive transport residual.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/15
  !
#include <petsc/finclude/petscvec.h>
  use petscvec
  implicit none
  
  class(realization_base_type) :: this
  
  PetscInt :: ierr
  
  if (this%field%flow_mass_transfer == PETSC_NULL_VEC) then
    call VecDuplicate(this%field%flow_xx,this%field%flow_mass_transfer, &
                      ierr);CHKERRQ(ierr)
  endif

end subroutine RealizCreateFlowMassTransferVec

! ************************************************************************** !

subroutine RealizCreateTranMassTransferVec(this)
  ! 
  ! Creates the Vec where mass transfer is summed prior to being added to
  ! the reactive transport residual.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/15
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none
  
  class(realization_base_type) :: this
  
  PetscInt :: ierr
  
  if (this%field%tran_mass_transfer == PETSC_NULL_VEC) then
    call VecDuplicate(this%field%tran_xx,this%field%tran_mass_transfer, &
                      ierr);CHKERRQ(ierr)
  endif

end subroutine RealizCreateTranMassTransferVec

! ************************************************************************** !

subroutine RealizationBaseStrip(this)
  ! 
  ! Deallocates members of base realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  ! 
  use Data_Mediator_module
  
  implicit none
  
  class(realization_base_type) :: this
  
  call FieldDestroy(this%field)

  nullify(this%output_option)
  
  call DiscretizationDestroy(this%discretization)
  
  if (associated(this%comm1)) then
    call this%comm1%Destroy()
    deallocate(this%comm1)
  endif
  nullify(this%comm1)
  
  call PatchDestroyList(this%patch_list)
  nullify(this%patch)

  call DebugDestroy(this%debug)
  
  call DataMediatorDestroy(this%flow_data_mediator_list)
  call DataMediatorDestroy(this%tran_data_mediator_list)

end subroutine RealizationBaseStrip

end module Realization_Base_class
