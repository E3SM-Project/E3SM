module PM_Base_Pointer_module

  use PM_Base_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  ! Since the context (ctx) for procedures passed to PETSc must be declared 
  ! as a "type" instead of a "class", object is a workaround for passing the 
  ! process model as context of a procedure where one can pass the
  ! pm_base_pointer_type to a procedure, declaring it as e.g.
  !
  ! type(pm_base_pointer_type) :: pm_ptr
  !
  ! and use the ptr:
  !
  ! pm_ptr%this%Residual
  !  
  type, public :: pm_base_pointer_type
    class(pm_base_type), pointer :: pm
  end type pm_base_pointer_type

  public :: PMResidual, &
            PMJacobian, &
            PMCheckUpdatePre, &
            PMCheckUpdatePost, &
            PMRHSFunction, &
            PMResidualPtr, &
            PMJacobianPtr, &
            PMCheckUpdatePrePtr, &
            PMCheckUpdatePostPtr, &
            PMRHSFunctionPtr

contains

! ************************************************************************** !

subroutine PMResidual(snes,xx,r,this,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module
  use Realization_Subsurface_class
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscsnes.h"

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(pm_base_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMResidual()'
#endif

  call this%Residual(snes,xx,r,ierr)

end subroutine PMResidual

! ************************************************************************** !

subroutine PMResidualPtr(snes,xx,r,this,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module
  use Realization_Subsurface_class
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscsnes.h"

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(pm_base_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMResidualPtr()'
#endif

  call this%pm%Residual(snes,xx,r,ierr)

end subroutine PMResidualPtr

! ************************************************************************** !

subroutine PMJacobian(snes,xx,A,B,this,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscsnes.h"

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  class(pm_base_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMJacobian()'
#endif

  call this%Jacobian(snes,xx,A,B,ierr)
    
end subroutine PMJacobian

! ************************************************************************** !

subroutine PMJacobianPtr(snes,xx,A,B,this,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscsnes.h"

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(pm_base_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMJacobianPtr()'
#endif

  call this%pm%Jacobian(snes,xx,A,B,ierr)
    
end subroutine PMJacobianPtr

! ************************************************************************** !

subroutine PMRHSFunction(ts,time,xx,ff,this,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/12/13
  ! 

  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscts.h"

  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  class(pm_base_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMRHSFunction()'
#endif

  call this%RHSFunction(ts,time,xx,ff,ierr)

end subroutine PMRHSFunction

! ************************************************************************** !

subroutine PMRHSFunctionPtr(ts,time,xx,ff,this,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/12/13
  ! 

  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscts.h"

  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  type(pm_base_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMRHSFunctionPtr()'
#endif

  call this%pm%RHSFunction(ts,time,xx,ff,ierr)

end subroutine PMRHSFunctionPtr

! ************************************************************************** !

subroutine PMCheckUpdatePre(line_search,X,dX,changed,this,ierr)
  ! 
  ! Wrapper for native call to XXXCheckUpdatePre
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  ! 
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscsnes.h"

  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  class(pm_base_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMCheckUpdatePre()'
#endif

  call this%CheckUpdatePre(line_search,X,dX,changed,ierr)
    
end subroutine PMCheckUpdatePre

! ************************************************************************** !

subroutine PMCheckUpdatePrePtr(line_search,X,dX,changed,this,ierr)
  ! 
  ! Wrapper for native call to XXXCheckUpdatePre
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  ! 
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscsnes.h"

  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  type(pm_base_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMCheckUpdatePrePtr()'
#endif

  call this%pm%CheckUpdatePre(line_search,X,dX,changed,ierr)
    
end subroutine PMCheckUpdatePrePtr

! ************************************************************************** !

subroutine PMCheckUpdatePost(line_search,X0,dX,X1,dX_changed,X1_changed,this, &
                             ierr)
  ! 
  ! Wrapper for native call to XXXCheckUpdatePost
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  ! 
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscsnes.h"

  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  class(pm_base_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMCheckUpdatePost()'
#endif

  call this%CheckUpdatePost(line_search,X0,dX,X1,dX_changed,X1_changed,ierr)
    
end subroutine PMCheckUpdatePost

! ************************************************************************** !

subroutine PMCheckUpdatePostPtr(line_search,X0,dX,X1,dX_changed,X1_changed, &
                                this,ierr)
  ! 
  ! Wrapper for native call to XXXCheckUpdatePost
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  ! 
  
  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscsnes.h"

  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  type(pm_base_pointer_type) :: this
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  print *, 'PMCheckUpdatePostPtr()'
#endif

  call this%pm%CheckUpdatePost(line_search,X0,dX,X1,dX_changed,X1_changed,ierr)
    
end subroutine PMCheckUpdatePostPtr

end module PM_Base_Pointer_module
